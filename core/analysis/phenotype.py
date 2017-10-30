#George Emanuel
#emanuega0@gmail.com

import multiprocessing
from functools import reduce

import numpy as np
from scipy.optimize import least_squares
from scipy.optimize import minimize_scalar

from core.analysis import analysis
from core.data import sequence

class BarcodeAggregator(object):

    def __init__(self, matchedBarcodes, bcToAATable):
        self.matchedBarcodes = matchedBarcodes
        self.bcToAATable = bcToAATable

        uniqueSequencedBarcodes = set(
                self.bcToAATable.get_unique_barcodes_as_ints())
        matchingUnique = {k: b for k,b in self.matchedBarcodes.items() \
                if uniqueSequencedBarcodes.__contains__(b)}

        uniqueBarcodes, uniqueInverse = np.unique(
                list(matchingUnique.values()), return_inverse = True)
        keys = np.array(list(matchingUnique.keys()))

        self.groupedCells = {uniqueBarcodes[i]: \
                [j for j in keys[uniqueInverse==i]]\
                for i in range(len(uniqueBarcodes))}

    def get_cells_by_barcode(self):
        return self.groupedCells 

    def reverse(self, x):
        return x[::-1]

    def get_barcodes_abundance_sorted(self, abundanceThreshold = 0):
        groupedCells = self.get_cells_by_barcode()

        cellCounts = {i: len(x) for i,x in groupedCells.items()}
        sortedCells = [(key, value) for key, value \
                in sorted(cellCounts.items(), key=self.reverse, reverse=True)]
        return [k for k,v in sortedCells if v > abundanceThreshold]

    def partition_phenotypes(self, phenotypeDictionary, abundanceThreshold=0,
            index = None):
        groupedCells = self.get_cells_by_barcode()

        partitionedPhenotypes = None
        if index is None:
            partitionedPhenotypes = {k: {j: phenotypeDictionary[j] \
                    for j in v} for k,v  in groupedCells.items() \
                    if len(v) > abundanceThreshold}
        else:
            partitionedPhenotypes = {k: {j: phenotypeDictionary[j][index] \
                    for j in v} for k,v  in groupedCells.items() \
                    if len(v) > abundanceThreshold}

        return partitionedPhenotypes

    def median_partition_phenotypes(
            self, phenotypeDictionary, abundanceThreshold=0, index=None):

        partitionedPhenotypes = self.partition_phenotypes(
                phenotypeDictionary, abundanceThreshold, index)
        return {k: np.median([x for x in v.values()], axis=0) \
                for k,v in partitionedPhenotypes.items() \
                if len(v) > abundanceThreshold}

    def std_partition_phenotypes(
            self, phenotypeDictionary, abundanceThreshold=0, index=None):

        partitionedPhenotypes = self.partition_phenotypes(
                phenotypeDictionary, abundanceThreshold, index)
        return {k: np.std([x for x in v.values()], axis=0) \
                for k,v in partitionedPhenotypes.items() \
                if len(v) > abundanceThreshold}

    def partition_by_mutant(
            self, phenotypeDictionary, abundanceThreshold=0, index = None):

        partitionedByBC = self.partition_phenotypes(
                phenotypeDictionary, abundanceThreshold = abundanceThreshold,
                index = index)
        mediansByBC = self.median_partition_phenotypes(
                phenotypeDictionary, abundanceThreshold = abundanceThreshold,
                index = index)

        byMutant = self._sort_measurements_by_mutant(partitionedByBC,
                {k: len(v) for k,v in partitionedByBC.items()})

        return byMutant

    def median_partition_by_mutant(
            self, phenotypeDictionary, abundanceThreshold=0, index = None):

        byMutant = self.partition_by_mutant(phenotypeDictionary,
                abundanceThreshold = 0, index = index)


        valuesByMutant = [reduce((lambda x,y: x+y),
                [list(x[0].values()) for x in currentMutant]) \
                        for currentMutant in byMutant]

        mediansByMutant = [(np.median(v), [x[3] for x in y]) \
                for y,v in zip(byMutant, valuesByMutant) \
                if len(v)>abundanceThreshold]

        return mediansByMutant


    def _sort_measurements_by_mutant(self, measurements, counts):

        barcodes = self.bcToAATable.get_barcode_list_as_ints()
        sequences = self.bcToAATable.get_sequence_list()
        uniqueSequences, uniqueI, uniqueCounts = np.unique(
                sequences, return_inverse = True, return_counts = True)

        medians = [[(measurements[barcodes[j]], counts[barcodes[j]], j, \
                self.bcToAATable.bcToAA[j])\
                for j in np.where(uniqueI==i)[0]\
                if barcodes[j] in measurements.keys()]\
                for i in range(len(uniqueSequences))]
    
        medians = [x for x in medians if len(x) > 0]
        return medians


class PhenotypeExtractor(analysis.AbstractAnalysisTask):

    def __init__(self, cellIntensities, phenotypeParameters):
        self.cellIntensities = cellIntensities
        self.phenotypeParameters = phenotypeParameters

        self.complete = False
        self.coreCount = 1

    def run(self):
        if self.coreCount == 1:
            self.phenotypes = {id: self._cell_intensities_to_phenotypes(
                        intensities) \
                    for id, intensities in self.cellIntensities.items()}
        else:
            pool = multiprocessing.Pool(processes = self.coreCount)
            self.phenotypes = dict(pool.map(
                self._intensities_to_id_phenotypes, 
                self.cellIntensities.items()))

        self.complete = True

    def set_cpu_count(self, cpuCount):
        self.coreCount = cpuCount

    def _cell_intensities_to_phenotypes(self, intensities):
        return [self._extract_phenotype(intensities, p) for p \
                in self.phenotypeParameters]

    def _intensities_to_id_phenotypes(self, intensityId):
        return (intensityId[0], self._cell_intensities_to_phenotypes(
            intensityId[1]))

    def _extract_phenotype(self, intensities, phenotypeProperties):
        if phenotypeProperties['type'] == 'Intensity':
            measuredIntensity = self._intensity_at_index(
                    intensities, phenotypeProperties['index'])
            if 'normalization' in phenotypeProperties:
                normalization = self._intensity_at_index(
                        intensities, phenotypeProperties['normalization'])
                measuredIntensity /= normalization
            return measuredIntensity

        if phenotypeProperties['type'] == 'Intensity Difference':
            intensity1 = self._intensity_at_index(
                    intensities, phenotypeProperties['intensity_1'])
            intensity2 = self._intensity_at_index(
                    intensities, phenotypeProperties['intensity_2'])
            difference = intensity1 - intensity2
            if 'normalization' in phenotypeProperties:
                normalization = self._intensity_at_index(
                        intensities, phenotypeProperties['normalization'])
                difference /= normalization
            return difference

        elif phenotypeProperties['type'] == 'Recovery':
            initialIntensity = self._intensity_at_index(
                    intensities, phenotypeProperties['initial'])
            finalIntensity = self._intensity_at_index(
                    intensities, phenotypeProperties['final'])
            recoveredIntensity = self._intensity_at_index(
                    intensities, phenotypeProperties['recovered'])
            return (recoveredIntensity-finalIntensity)/ \
                    (initialIntensity-finalIntensity)

        return np.NaN

    def _intensity_at_index(self, intensities, index):
        return intensities[index[0]][index[1]]

    def result(self):
        if self.complete:
            return self.phenotypes
        
    def to_string(self):
        pass




