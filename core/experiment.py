#George Emanuel
#emanuega0@gmail.com
#Copyright Presidents and Fellows of Harvard College, 2017.

import multiprocessing
import gc
import json
import pickle
import os
import json
import logging

import numpy as np

from core.analysis import normalization 
from core.data import imagenormalize
from core.data.image import DaxImagingDataSet
from core.data import data
from core.data import alignment
from core.data import sequence
from core.data.cells import CellCollection
from core.analysis import cellfinder 
from core.analysis import aligner
from core.analysis import cellfinder
from core.analysis import phenotype
from core.analysis import analyzer
from core.analysis import bitcalling
from core.analysis import stackanalyzer
from core.visual import summaryplots



class SequenceFunctionExperiment(object):

    def __init__(self, parametersFile = None, imagingExperiment = None, 
            sequencingExperiment = None, sequencingLibrary = None, 
            positionFile = None, cores=-1, overwriteParameters = False):
        if parametersFile is None and imagingExperiment is None:
            logging.error('Unable to load experiment if neither ' \
                    + 'experiment or parameters are specified')
            return

        if parametersFile is None:
            self.experimentName = imagingExperiment
            self.parametersName = self._list_analysis_for_experiment(
                    imagingExperiment)[0]
            self._load_parameters()
        else:
            self.parametersName = parametersFile.rsplit('.', 1)[0]
            self.parametersFile = parametersFile

            oldParameters = None
            if imagingExperiment is not None and \
                    self._list_analysis_for_experiment(imagingExperiment) \
                    is not None: 
                self.experimentName = imagingExperiment
                analysisList = self._list_analysis_for_experiment(
                        imagingExperiment)
                if self.parametersName in analysisList:
                    self._load_parameters()
                    oldParameters = self.parameters
                elif any([self.parametersName in x for \
                        x in analysisList]):
                    parametersIndex = np.where([
                            self.parametersName in x for \
                                    x in analysisList]) 
                    self.parametersName = analysisList[parametersIndex[0][0]]
                    self._load_parameters()
                    oldParameters = self.parameters


            if oldParameters is None or overwriteParameters:
                with open(self.parametersFile) as pFile:
                    self.parameters = json.load(pFile) 

            if oldParameters is not None:
                for k,v in oldParameters.items():
                    if k not in self.parameters:
                        self.parameters[k] = v

        if imagingExperiment is not None:
            self.parameters['imaging_experiment'] = imagingExperiment
        if sequencingExperiment is not None:
            self.parameters['sequencing_experiment'] = sequencingExperiment
        if sequencingLibrary is not None:
            self.parameters['sequencing_library'] = sequencingLibrary
        if positionFile is not None:
            self.parameters['position_file'] = positionFile


        self.experimentName = self.parameters['imaging_experiment']
        self._save_parameters()


        self.coreCount = multiprocessing.cpu_count() if cores is -1 else cores 

        self.positions = np.loadtxt(
                data.__POSITIONSHOME__ + '/' + self.parameters['position_file'],
                delimiter = ',')

        self.imageData = DaxImagingDataSet.from_experiment(
                self.experimentName, 
                self.parameters['images_per_location'])

        self.normalizerClass = 'MinMeanImageNormalization'
        self.alignmentSetClass = 'ImageSetAlignment'
        self.cellListClass = 'CellCollection'
        self.bitMeasurementsSaveName = 'BitMeasurements'
        self.normalizedBitMeasurementSaveName = 'NormalizedBitMeasurements'
        self.thresholdNormBitMeasurementSaveName = \
                'ThresholdNormalizedBitMeasurements'
        self.phenotypeSaveName = 'PhenotypeMeasurements'
        self.normalizedPhenotypeSaveName = 'NormalizedPhenotypeMeasurements'
        self.intensitiesSaveName = 'IntensityDict'
        self.bitAssignmentsSaveName = 'CellAssignments'
        self.bitCallerParametersSaveName = 'BitCallParameters'
        self.analysisSaveName = 'AnalysisResults'

        self.normalizer = None
        self.cellList = None
        self.alignments = None
        self.bitMeasurements = None
        self.normalizedBitMeasurements = None
        self.thresholdNormalizedBitMeasurements = None
        self.phenotypes = None
        self.normalizedPhenotypes = None
        self.intensities = None
        self.measuredBarcodes = None
        self.analysisResults = None

        logging.basicConfig(
                filename = os.sep.join([self._cache_home(), 'log.txt']),
                format = '%(levelname)s:%(asctime)s-%(message)s',
                level = logging.INFO)
        logging.info('---------------------------------------------')
        logging.info('Beginning analysis with parameters: ' 
                + self.parametersName)
        if imagingExperiment is not None:
            logging.info('Imaging experiment: ' + imagingExperiment)
        if sequencingExperiment is not None:
            logging.info('Sequencing experiment: ' + sequencingExperiment)
        if sequencingLibrary is not None:
            logging.info('Sequencing library: ' + sequencingLibrary)
        if positionFile is not None:
            logging.info('Positions: ' + positionFile)
        logging.info('Cache home: ' + self._cache_home())    
        logging.info('Figure home: ' + self._figure_home())

    def get_phenotypes(self, reanalyze = False):
        logging.info('Getting phenotypes')
        if self.phenotypes is None:
            self.phenotypes = self._load(
                    self._generate_load_path(self.phenotypeSaveName))

        if self.phenotypes is None or reanalyze:
            logging.info('Unable to load phenotypes. Generating phenotypes')

            cells = self.get_cells()
            cellIntensities = self.get_cell_intensities()

            phenotypeExtractor = phenotype.PhenotypeExtractor( 
                    cellIntensities, self.parameters['phenotypes'])
            phenotypeExtractor.set_cpu_count(self.coreCount)
            phenotypeExtractor.run()

            self.phenotypes = phenotypeExtractor.result()
            self._save(self.phenotypes, self.phenotypeSaveName)
            logging.info('Phenotypes extracted and saved')

        return self.phenotypes

    def get_analysis(self, reanalyze = False):
        logging.info('Generating analysis')
        if self.analysisResults is None:
            self.analysisResults = self._load(
                    self._generate_load_path(self.analysisSaveName))

        if self.analysisResults is None or reanalyze:
            logging.info('Unable to load analysis results. Beginning ' + 
                    'to generate analysis results.')
            a = analyzer.Analyzer(self.get_matched_barcodes(),
                    self.get_barcode_to_AA_table(),
                    self.get_cell_intensities(),
                    self.get_phenotypes(reanalyze=reanalyze),
                    self.parameters['phenotypes'],
                    self.parameters['analysis'])
            a.run()

            self.analysisResults = a.result()

            self._save(self.analysisResults, 
                    self.analysisSaveName)
            logging.info('Successfully generated analysis results')

        return self.analysisResults


    def set_normalizer(self, normalization):
        self.normalization = normalization
        self._save(self.normalization)


    def get_normalizer(self):
        logging.info('Getting normalizer')
        if self.normalizer is None:
            self.normalization = self._load(
                    self._generate_load_path(self.normalizerClass))
        
        if self.normalization is None:
            logging.info('Unable to load normlizer. Beginning to '
                    +  'generate normalizer')
            normGenerator = normalization.MinMeanNormalizationGenerator(
                    self.imageData, 
                    self.parameters['normalization_calculate_index'])
            normGenerator.run()
            self.normalization = normGenerator.result()
            summaryplots.normalizer_images(
                    self.normalization, self._figure_home())
            self._save(self.normalization)
            logging.info('Successfully generated normalizer')

        normalizationHints = self.parameters['normalization_hints'] \
                if 'normalization_hints' in  self.parameters \
                else None
        self.normalizer = normalization.MinMeanNormalizer(self.normalization,
                normalizationHints)

        return self.normalizer

    def get_image_alignments(self):
        logging.info('Getting alignments')
        if self.alignments is None:
            self.alignments = self._load(self._generate_load_path(
                self.alignmentSetClass))

        if self.alignments is None:
            logging.info('Unable to load alignments. Begening to generate ' +
                    'alignments')
            self._analyze_image_stacks()
            logging.info('Successfully generated alignments')

        return self.alignments

    def get_normalized_aligned_images_for_region(self, regionIndex):
        alignments = self.get_image_alignments()
        normalizer = self.get_normalizer()

        imageStack = self.imageData.images_at_position(regionIndex)
        imageStack = normalizer.normalize_image_stack(imageStack)

        imageAligner = aligner.ImageAligner(
                imageStack, 
                alignments.transformations_for_position(regionIndex))
        imageAligner.run()
        alignedStack = imageAligner.result()

        return alignedStack

    def get_cells(self, reanalyze=False):
        logging.info('Getting cells')
        if self.cellList is None and not reanalyze:
            self.cellList = self._load(self._generate_load_path(
                self.cellListClass))

        if self.cellList is None or reanalyze:
            logging.info('Unable to load cells. Beginning to generate ' + 
                'cells')
            self._analyze_image_stacks()
            logging.info('Successfully generated cells')

        return self.cellList

    def get_cell_intensities(self, reanalyze=False):
        logging.info('Getting cell intensities')
        if self.intensities is None and not reanalyze:
            self.intensities = self._load(self._generate_load_path(
                self.intensitiesSaveName))

        if self.intensities is None or reanalyze:
            logging.info('Unable to load cell intensities. Beginning to ' 
                + 'generate cell intensities')
            self._analyze_image_stacks()

        return self.intensities

    def _analyze_image_stacks(self):
        logging.info('Extracting alignments, cells, and cell intensities '
                + ' from image stacks')
        normalizer = self.get_normalizer()

        stackAnalyzer = stackanalyzer.ImageStackAnalyzer(
                self.imageData, normalizer, self.positions, self.parameters)
        stackAnalyzer.set_core_count(self.coreCount)
        stackAnalyzer.run()

        self.alignments = stackAnalyzer.get_alignments()
        self.cellList = stackAnalyzer.get_cells()
        self.intensities = stackAnalyzer.get_intensities()

        self._save(self.alignments)
        self._save(self.cellList)
        self._save(self.intensities)

        summaryplots.alignment_offset_histogram(self.alignments,
                self._figure_home())

    def get_bit_measurements(self):
        logging.info('Getting bit measurements')
        if (self.bitMeasurements is None):
            self.bitMeasurements = self._load(
                    self._generate_load_path(
                        self.bitMeasurementsSaveName))

        if (self.bitMeasurements is None):
            logging.info('Unable to load bit measurements. Beginning to ' 
                + 'generate bit measurements')
            intensities = self.get_cell_intensities()

            readoutIndexes = self.parameters['readouts']
            readoutRounds = [x[1] for x in readoutIndexes]
            readoutRoundCount = max(readoutRounds) - min(readoutRounds) + 1 
            readoutStart = \
                    self.imageData.imaging_round_count() - readoutRoundCount 

            
            bitIndexes = self.parameters['bits']
            self.bitMeasurements = {\
                    id: np.array(self._bits_from_measurements(
                        i, bitIndexes, readoutIndexes, readoutStart))\
                    for id, i in intensities.items()}

            self._save(self.bitMeasurements, self.bitMeasurementsSaveName)
            logging.info('Finished extracting bit measurements')

        return self.bitMeasurements


    def get_normalized_bit_measurements(self):
        logging.info('Getting normalized bit measurements')
        if self.normalizedBitMeasurements is None:
            self.normalizedBitMeasurements = self._load(
                    self._generate_load_path(
                        self.normalizedBitMeasurementSaveName))

        if self.normalizedBitMeasurements is None:
            logging.info('Unable to load normalized bit measurements. ' 
                + 'Beginning to generate normalized bit measurements')
            bitMeasurements = self.get_bit_measurements()
            allMeasurements = np.array(list(bitMeasurements.values()))

            medians = np.median(allMeasurements, axis=0)
            mins = np.min(allMeasurements, axis=0)

            self.normalizedBitMeasurements = \
                {k: (v-mins)/(medians-mins) for k,v in bitMeasurements.items()}

            self._save(
                    self.normalizedBitMeasurements, 
                    self.normalizedBitMeasurementSaveName)
            logging.info('Finished extracting normalizeb bit measurements')

        return self.normalizedBitMeasurements

    def get_thresholded_normalized_bit_measurements(self):
        logging.info('Getting thresholded bit measurements')
        if self.thresholdNormalizedBitMeasurements is None:
            self.thresholdNormalizedBitMeasurements = self._load(
                self._generate_load_path(
                    self.thresholdNormBitMeasurementSaveName))

        if self.thresholdNormalizedBitMeasurements is None:
            logging.info('Unable to load thresholded bit measurements. '  
                    + 'Beginning to generate thresholded bit measurements')
            normBitMeasurements = self.get_normalized_bit_measurements()

            threshold = self.parameters['bit_threshold']
            self.thresholdNormalizedBitMeasurements = \
                    {k: v for k,v in normBitMeasurements.items() if \
                    bitcalling.bit_measurement_passes_threshold(v, threshold)}

            summaryplots.cumulative_bit_intensity_plot(
                    self.thresholdNormalizedBitMeasurements, 
                    self._figure_home())
            summaryplots.bit_measurement_histograms2d(
                    self.thresholdNormalizedBitMeasurements,
                    self._figure_home())

            bitCaller = bitcalling.BitCaller(
                    self.thresholdNormalizedBitMeasurements)
            measuredBarcodeInts = self.barcodes_as_ints(
                    bitCaller._call_barcodes())
            summaryplots.barcode_abundance_plot(
                    measuredBarcodeInts, nameSuffix = 'Unoptimized',
                    savePath = self._figure_home())


            self._save(
                    self.thresholdNormalizedBitMeasurements,
                    self.thresholdNormBitMeasurementSaveName)
            logging.info('Finished extracting thersholded bit measurements')

        return self.thresholdNormalizedBitMeasurements



    def get_measured_barcodes(self):
        logging.info('Getting measured barcodes')
        if self.measuredBarcodes is None:
            self.measuredBarcodes = self._load(
                    self._generate_load_path(self.bitAssignmentsSaveName))

        if self.measuredBarcodes is None:
            logging.info('Unable to load measured barcodes. ' +
                    'Beginning to generate mesaured barcodes')
            barcodeToAA = self.get_barcode_to_AA_table()
            thresholdNormalizedMeasurements = \
                    self.get_thresholded_normalized_bit_measurements()

            
            bitCaller = bitcalling.BitCaller(
                    thresholdNormalizedMeasurements,
                    knownBarcodes = barcodeToAA.get_barcode_list())
            logging.info('Initially ' + str(bitCaller.fraction_matching()) \
                    + ' percent match')

            bitCaller.run()
            self.measuredBarcodes = bitCaller.result()

            summaryplots.bit_measurement_histograms2d_with_threshold(
                    thresholdNormalizedMeasurements,
                    bitCaller.get_slopes(), bitCaller.get_intercepts(),
                    self._figure_home())

            summaryplots.call_rate_vs_threshold_plot(
                    bitCaller, 
                    self.get_normalized_bit_measurements(),
                    self._figure_home())

            barcodeInts = self.barcodes_as_ints(self.measuredBarcodes)
            summaryplots.barcode_abundance_plot(
                    barcodeInts, 'All', self._figure_home())

            sequencedBarcodeSet = set(barcodeToAA.get_barcode_list_as_ints())
            matchedBarcodes = {k: b for k,b in barcodeInts.items() \
                    if sequencedBarcodeSet.__contains__(b)}
            summaryplots.barcode_abundance_plot(
                    matchedBarcodes, 'Matched', self._figure_home())

            self._save(bitCaller.get_optimal_parameters(), 
                    self.bitCallerParametersSaveName)
            self._save(bitCaller.result(),
                    self.bitAssignmentsSaveName)
            logging.info('Finished extracting measured barcodes')

        return self.measuredBarcodes


    def get_matched_barcodes(self):
        if 'sequencing_experiment' not in self.parameters:
            return None

        barcodeInts = self.barcodes_as_ints(self.get_measured_barcodes())
        bcToAA = self.get_barcode_to_AA_table()
        sequencedBarcodeSet = set(bcToAA.get_unique_barcodes_as_ints())
        matchedBarcodes = {k: b for k,b in barcodeInts.items() \
                if sequencedBarcodeSet.__contains__(b)}

        return matchedBarcodes


    def get_cells_with_barcode(self, barcode):
        bcToAA = self.get_barcode_to_AA_table()
        barcodeInt = bcToAA.get_barcode_as_int(barcode)
        measuredBarcodeInts = self.barcodes_as_ints(
                self.get_measured_barcodes())
        matchedCells = [k for k,b in measuredBarcodeInts.items() \
                if b == barcodeInt]
        return matchedCells


    def get_cells_with_barcode_int(self, barcodeInt):
        measuredBarcodeInts = self.barcodes_as_ints(
                self.get_measured_barcodes())
        matchedCells = [k for k,b in measuredBarcodeInts.items() \
                if b == barcodeInt]
        return matchedCells


    def get_measurements_for_mutant(self, mutantSequence, analysisIndex=1):
        bcToAA = self.get_barcode_to_AA_table()
        ba = phenotype.BarcodeAggregator(self.get_matched_barcodes(), bcToAA)
        pm = ba.partition_by_mutant(self.get_analysis()[analysisIndex])

        mutantIndex = np.where([x[0][3]['aaSequence'] == mutantSequence \
                for x in pm])[0]

        if len(mutantIndex) == 0:
            return None

        return pm[mutantIndex[0]]

    def get_analysis_by_mutant(self, analysisIndex, normalizer = 1):
        analysis = self.get_analysis()[analysisIndex]
        if normalizer != 1:
            normalizedAnalysis = {k: v/normalizer for k,v in analysis.items()}
        else:
            normalizedAnalysis = analysis

        ba = phenotype.BarcodeAggregator(self.get_matched_barcodes(),
                self.get_barcode_to_AA_table())
        pm = ba.partition_by_mutant(normalizedAnalysis)

        return pm

    def get_phenotypes_by_mutant(self, phenotypeIndex):
        phenotypes = self.get_phenotypes()
        desiredPhenotype = {k: x[phenotypeIndex] for k,x in \
                phenotypes.items()}

        ba = phenotype.BarcodeAggregator(self.get_matched_barcodes(),
                self.get_barcode_to_AA_table())
        pm = ba.partition_by_mutant(desiredPhenotype)

        return pm

    def barcodes_as_ints(self, barcodes):
        return {k: sequence.bin_array_to_int(b) for k,b in barcodes.items()}


    def get_barcode_to_AA_table(self):
        if 'sequencing_experiment' not in self.parameters:
            return None

        return sequence.MatlabBarcodeToAA(
                self.parameters['sequencing_experiment'],
                self.parameters['sequencing_library'],
                len(self.parameters['bits']))

    def _bits_from_measurements(
            self, measurements, bitIndexes, readoutIndexes, readoutOffset):

        return [self._bit_measurements(
            measurements, b, readoutIndexes, readoutOffset) \
                    for b in bitIndexes]

    def _bit_measurements(
            self, measurements, bitIndexes, readoutIndexes, readoutOffset):
        
        readout0 = readoutIndexes[bitIndexes[0]-1]
        bit0 = measurements[readoutOffset + readout0[1]-1][readout0[0]-1]

        readout1 = readoutIndexes[bitIndexes[1]-1]
        bit1 = measurements[readoutOffset + readout1[1]-1][readout1[0]-1]

        return [bit0, bit1]

    def _load(self, loadPath):
        if not os.path.exists(loadPath):
            return None

        loadedObject = None
        with open(loadPath, 'rb') as readFile:
            loadedObject = np.load(readFile).tolist()
        logging.info('Loaded ' + loadPath)

        return loadedObject

    def _save(self, saveObject, saveName = None):

        savePath = self._generate_save_path(saveObject, saveName)

        if not os.path.exists(os.path.dirname(savePath)):
            try:
                os.makedirs(os.path.dirname(savePath))
            except OSError as exc:
                if exc.errno != errno.EEXIST:
                    raise

        logging.info('Saving ' + saveObject.__class__.__name__ + ' to ' 
                + savePath)
        with open(savePath, 'wb') as writeFile:
            np.save(writeFile, saveObject)

    def _figure_home(self):
        figureHome = os.sep.join([self._cache_home(), 'figures'])
        self._create_path_if_needed(figureHome)

        return figureHome

    def _create_path_if_needed(self, path):
        if not os.path.exists(path):
            logging.debug('Path ' + path + ' doesn\'t exist, so it is being '
                    + 'created')
            try:
                os.makedirs(path)
            except OSError as exc:
                if exc.errno != errno.EEXIST:
                    raise

    def _cache_home(self):
        cacheHome = os.sep.join([data.__CACHEHOME__, self.experimentName,
            self.parametersName])
        self._create_path_if_needed(cacheHome)

        return cacheHome
        
    def _generate_load_path(self, loadObjectClassName):
        return os.sep.join([self._cache_home(), loadObjectClassName + '.npy'])

    def _generate_save_path(self, saveObject, saveName = None):
        if saveName is None:
            saveName = type(saveObject).__name__
        if saveName.find('.') != -1:
            return os.sep.join([self._cache_home(), saveName])
        else:
            return os.sep.join([self._cache_home(), saveName + '.npy'])

    def _write_fasta(self, fastaName, fastaList):
        savePath = self._generate_save_path(
                fastaList, fastaName + '.fasta')

        if not os.path.exists(os.path.dirname(savePath)):
            try:
                os.makedirs(os.path.dirname(savePath))
            except OSError as exc:
                if exc.errno != errno.EEXIST:
                    raise

        logging.info('Saving ' + fastaName + ' to ' + savePath)

        with open(savePath, 'w') as fastaFile:
            for fastaEntry in fastaList:
                fastaFile.write('>' + fastaEntry['name'] + '\n')
                fastaFile.write(fastaEntry['sequence'] + '\n')


    def _save_parameters(self):
        parametersSave = os.sep.join([self._cache_home(), 'parameters.json'])
        with open(parametersSave, 'w') as f:
            json.dump(self.parameters, f)

    def _load_parameters(self):
        parametersPath = os.sep.join([self._cache_home(), 'parameters.json'])
        with open(parametersPath) as f:
            self.parameters = json.load(f) 

    def _list_analysis_for_experiment(self, imagingExperimentName):
        experimentCachePath = os.sep.join(
                [data.__CACHEHOME__, imagingExperimentName])

        if not os.path.exists(experimentCachePath):
            return None

        return next(os.walk(experimentCachePath))[1]


         

