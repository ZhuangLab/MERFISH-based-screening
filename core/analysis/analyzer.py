#George Emanuel
#emanuega0@gmail.com
#Copyright Presidents and Fellows of Harvard College, 2017.

import numpy as np

from core.analysis import analysis
from core.analysis import phenotype

from scipy.optimize import least_squares
from scipy.optimize import curve_fit

class Analyzer(analysis.AbstractAnalysisTask):

    def __init__(
            self, measuredBarcodes, bcToAA, cellIntensities, 
            phenotypes, phenotypeParameters, analysisParameters):

        self.measuredBarcodes = measuredBarcodes
        self.bcToAA = bcToAA
        self.cellIntensities = cellIntensities
        self.phenotypes = phenotypes
        self.phenotypeParameters = phenotypeParameters
        self.analysisParameters = analysisParameters
    
        self.analysisResults = None
        self.complete = False

    def run(self):
        if self.measuredBarcodes is not None and self.bcToAA is not None:
            self.barcodeAggregator = phenotype.BarcodeAggregator(
                    self.measuredBarcodes, self.bcToAA)

        self.analysisResults = []
        for p in self.analysisParameters:
            self.analysisResults.append(
                    self._run_analysis(p, self.analysisResults))

        self.complete = True

    def _run_analysis(self, analysisParameters, previousResults):
        result = None


        if analysisParameters['type'] == 'Background' and \
                analysisParameters['method'] == 'Correlation':

            measurements1, groupedMeasurements1 = \
                    self._extract_measurements(
                            analysisParameters['intensity 1'], 
                            analysisParameters['source 1'] \
                                    if 'source 1' in analysisParameters \
                                    else None,
                            previousResults)

            measurements2, groupedMeasurements2 = \
                    self._extract_measurements(
                            analysisParameters['intensity 2'], 
                            analysisParameters['source 2'] \
                                    if 'source 2' in analysisParameters \
                                    else None,
                            previousResults)

            correlationCoefficients = {k: \
                    np.corrcoef([x \
                            for x in groupedMeasurements1[k].values()],
                        [x \
                            for x in groupedMeasurements2[k].values()]
                        )[0][1]\
                    for k in groupedMeasurements1.keys() \
                        if len(groupedMeasurements1[k]) > \
                            analysisParameters['count threshold']}

            indexBelowThreshold = [k for k,v \
                    in correlationCoefficients.items() \
                    if v < analysisParameters['correlation threshold']]
            mediansBelowThreshold = [np.median(
                [x for x in groupedMeasurements1[k].values()]) \
                        for k in indexBelowThreshold]

            result = np.median(mediansBelowThreshold)

        elif analysisParameters['type'] == 'Background' and \
                analysisParameters['method'] == 'Predetermined':

            result = analysisParameters['value']

        elif analysisParameters['type'] == 'Relative Intensity':
            measurements1 = \
                    self._extract_measurements(
                            analysisParameters['intensity 1'], 
                            analysisParameters['source 1'] \
                                    if 'source 1' in analysisParameters \
                                    else None,
                            previousResults, group = False)

            measurements2 = \
                    self._extract_measurements(
                            analysisParameters['intensity 2'], 
                            analysisParameters['source 2'] \
                                    if 'source 2' in analysisParameters \
                                    else None,
                            previousResults, group = False)
            
            background = self._extract_background(
                    analysisParameters, previousResults)

            result = {k: (measurements1[k]-background)/measurements2[k] \
                    for k in measurements1.keys()}

        elif analysisParameters['type'] == 'Double Exponential Fit':
            background = self._extract_background(
                    analysisParameters, previousResults)

            result = {k: self.extract_double_exponential(intensities,
                analysisParameters, background) for k,intensities \
                        in self.cellIntensities.items()}

        elif analysisParameters['type'] == 'Fit Parameter':
            fitResults = self._extract_measurements(
                    analysisParameters['fit name'], 'analysis',
                    previousResults, group=False)
            fitIndex = analysisParameters['parameter index']
            result = {k: x[fitIndex] for k, x in fitResults.items()}

        elif analysisParameters['type'] == 'Relative Difference':
            intensity1 = self._extract_measurements(
                    analysisParameters['intensity 1'],
                    analysisParameters['source 1'] \
                            if 'source 1' in analysisParameters \
                            else None,
                    previousResults, group=False)
            intensity2 = self._extract_measurements(
                    analysisParameters['intensity 2'],
                    analysisParameters['source 2'] \
                            if 'source 2' in analysisParameters \
                            else None,
                    previousResults, group=False)
            intensity3 = self._extract_measurements(
                    analysisParameters['intensity 3'],
                    analysisParameters['source 3'] \
                            if 'source 1' in analysisParameters \
                            else None,
                    previousResults, group=False)

            result = {k: (intensity1[k] - intensity2[k])/intensity3[k] \
                    for k in intensity1.keys()}


        return result

    def extract_double_exponential(self, intensities, analysisParameters, 
            background):
        startFrame = analysisParameters['start frame']
        endFrame = analysisParameters['end frame']
        imagingRound = analysisParameters['round']
        inputRate = analysisParameters['fast rate']
        normalizationFrame = analysisParameters['normalization frame']

        bleachCurve = \
                (intensities[imagingRound][startFrame:endFrame] \
                    -background)/ \
                    intensities[imagingRound][normalizationFrame]

        if 'fast rate type' not in analysisParameters \
                or analysisParameters['fast rate type'] == 'fixed':
            def double_exp_constrained(
                    u, fastAmplitude, slowAmplitude, slowRate):
                return fastAmplitude*np.exp(inputRate*u)+ \
                        slowAmplitude*np.exp(slowRate*u) 

            u = np.array(range(len(bleachCurve)))
            try:
                fitparams = curve_fit(double_exp_constrained, u, bleachCurve,
                        p0 = [bleachCurve[0]-bleachCurve[1], bleachCurve[2], 
                            -0.01]) 
            except RuntimeError:
                fitparams = [[np.nan, np.nan, np.nan]]

        elif analysisParameters['fast rate type'] == 'constrained':
            def double_exp(
                    u, fastAmplitude, slowAmplitude, slowRate, fastRate):
                return fastAmplitude*np.exp(fastRate*u)+ \
                        slowAmplitude*np.exp(slowRate*u)

            u = np.array(range(len(bleachCurve)))
            try: 
                fitparams = curve_fit(double_exp, u, bleachCurve,
                        p0 = [1, 1, -0.01, -2.3],
                        bounds = ((-np.inf, -np.inf, -np.inf, -np.inf),
                            (np.inf, np.inf, np.inf, inputRate)))
            except RuntimeError:
                fitparams = [[np.nan, np.nan, np.nan, np.nan]]
        else: 
            print('Double exponential fit type not recognized')

        return fitparams[0]

    def extract_amplitude_ratio(self, intensities, analysisParameters, 
            background):
        startFrame = analysisParameters['start frame']
        endFrame = analysisParameters['end frame']
        imagingRound = analysisParameters['round']
        fastRate = analysisParameters['fast rate']

        bleachCurve = intensities[imagingRound][startFrame:endFrame]

        def double_exp_constrained(x, u, y):
            return background+x[0]*np.exp(fastRate*u)+x[1]*np.exp(x[2]*u) - y

        u = np.array(range(len(bleachCurve)))
        fitparams = least_squares(double_exp_constrained, [1, 1, -0.01], 
                args = (u, bleachCurve))['x']

        return fitparams[0]/(fitparams[1]+fitparams[0])

    def _extract_background(self, analysisParameters, analysisResults):
        if 'background estimate' in analysisParameters:
            index = [i for i,v in enumerate(self.analysisParameters) \
                    if v['name'] == \
                        analysisParameters['background estimate']][0]
            background = analysisResults[index]
        else:
            background = 0

        return background

    def _extract_measurements(self, index, source = None, 
            previousResults = None, group = True):
        if source is None or source == 'phenotypes':
            measurements = \
                    self._extract_phenotype_measurements(index, group=group)
        elif source == 'analysis':
            measurements = \
                    self._extract_analysis_measurements(index,
                            previousResults, group=group)

        return measurements

    def _extract_phenotype_measurements(self, phenotypeName, group):
        index = [i for i,v in enumerate(
            self.phenotypeParameters)\
                if v['name'] == phenotypeName][0]

        measurements1 = {k: x[index] for k,x in self.phenotypes.items()}
        if group:
            groupedMeasurements1 = self.barcodeAggregator.partition_phenotypes(
                    self.phenotypes, index = index)
            return measurements1, groupedMeasurements1
        else:
            return measurements1

    def _extract_analysis_measurements(self, analysisName, analysisResults,
            group):
        index = [i for i,v in enumerate(
            self.analysisParameters) \
                if v['name'] == analysisName][0]

        if group:
            measurements1 = self.barcodeAggregator.partition_phenotypes(
                    analysisResults[index])
            return analysisResults[index], measurements1
        else:
            return analysisResults[index]



    def result(self):
        if self.complete:
            return self.analysisResults

    def to_string():
        return None
