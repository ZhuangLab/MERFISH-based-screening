#George Emanuel
#emanuega0@gmail.com
#Copyright Presidents and Fellows of Harvard College, 2017.

import json
import os
from enum import Enum

import scipy.io as sio
import numpy as np

from core.data import data


aaList = ['A','R','N','D','C','E','Q','G','H','I','L','K','M',\
        'F','P','S','T','W','Y','V','*']

class MutationType(Enum):
    unknown = 0
    none = 1
    substitution = 2
    deletion = 3
    insertion = 4
    multisubstitution = 5


def bin_array_to_int(binaryArray):
    return binaryArray.dot(1 << np.arange(binaryArray.shape[-1] -1, -1, -1))


class MatlabBarcodeToAA(data.AbstractData):

    def __init__(self, experiment, sample, bitCount):
        matPath = os.sep.join([data.__SEQUENCEHOME__, experiment, sample,
            'BCtoAA.mat'])
        self.bcToAAMatFull = sio.loadmat(matPath)
        self.bcToAAMat = sio.loadmat(matPath)['bcToAA'][0]
        self.bcToAA = [self._mat_BC_to_AA(x) \
                for x in enumerate(self.bcToAAMat) \
                if len(x[1][1]) > 0]
        self.bitCount = bitCount

        if len(self.bcToAA) > 0:
            self.sequenceLength = len(self.bcToAA[0]['aaSequence'])
    
    def get_aa_list(self):
        return aaList

    def get_entry_for_UMI(self, umi):
        for i in range(len(self.bcToAA)):
            if self.bcToAA[i]['umi'] == umi:
                return self.bcToAA[i]
        
        return None

    def get_entries_for_sequence(self, sequence):
        return [x for x in self.bcToAA \
                if x['aaSequence'] == sequence]

    def get_entries_for_barcode(self, barcodeInt):
        return [x for x in self.bcToAA \
                if self.get_barcode_as_int(x['barcode']) == barcodeInt]

    def get_wt_sequence(self):
        wtIndexes = self.get_wt_indexes()
        if len(wtIndexes) is 0:
            return ''

        return self.bcToAA[wtIndexes[0]]['aaSequence']
    
    def get_sequence_list(self):
        return [x['aaSequence'] for x in self.bcToAA]

    def get_barcode_list(self):
        if self.bitCount is None:
            return [x['barcode'] for x in self.bcToAA]
        else:
            return [x['barcode'][0:self.bitCount] for x in self.bcToAA]

    def get_barcode_as_int(self, barcode):
        if self.bitCount is None:
            return bin_array_to_int(barcode)
        else:
            return bin_array_to_int(barcode[0:self.bitCount]) 

    def get_barcode_list_as_ints(self):
        return [self.get_barcode_as_int(x['barcode']) for x in self.bcToAA]

    def get_unique_barcodes_as_ints(self):
        barcodeList = self.get_barcode_list_as_ints()

        uniqueBarcodes, counts = np.unique(barcodeList, return_counts=True)

        return [x for x,c in zip(uniqueBarcodes, counts) if c == 1]

    def get_barcode_frequencies(self):
        frequencies = np.zeros(2**len(self.bcToAA[0]['barcode']))
        for bcInt in self.get_barcode_list_as_ints():
            frequencies[bcInt] += 1
        return frequencies

    def get_wt_indexes(self):
        return [i for i,x in enumerate(self.bcToAA) if self.is_wt(x)]

    def _mat_BC_to_AA(self, matEntry):
        i = matEntry[0]
        matEntry = matEntry[1]
        bcToAAElement = {}
        bcToAAElement['index'] = i
        bcToAAElement['umi'] = matEntry[0][0]
        bcToAAElement['aaSequence'] = str(matEntry[1][0].tostring(),'utf-8')
        bcToAAElement['alignment'] = matEntry[2][:]
        bcToAAElement['barcode'] = matEntry[3][0]
         
        if len(matEntry) > 4:
            bcToAAElement['aaConfidence'] = matEntry[4][0][0]
        
        if len(matEntry) > 5:
            bcToAAElement['bcConfidence'] = matEntry[5][0][0]

        return bcToAAElement

    def entries_for_type(self, mutationType):
        return [x for x in self.bcToAA \
                if self._classify_mutation(x) is mutationType]

    def _classify_mutation(self, bcToAAEntry):
        entryAlignment = bcToAAEntry['alignment']
        sequenceLength = len(entryAlignment[0])
        containsGap = entryAlignment[0].count('-') == 1
        mismatchCount = len(entryAlignment[1]) - entryAlignment[1].count('|')
        if (mismatchCount == 1 and not containsGap):
            return MutationType.substitution
        elif mismatchCount > 1 and not containsGap:
            return MutationType.multisubstitution
        elif (containsGap and mismatchCount == 2):
            gapPositionTop = entryAlignment[0].find('-')
            gapPositionBottom = entryAlignment[2].find('-')
            if (gapPositionTop == 0 or gapPositionTop == sequenceLength - 1):
                return MutationType.deletion
            elif (gapPositionBottom == 0 or 
                    gapPositionBottom == sequenceLength - 1):
                return MutationType.insertion
            else:
                return MutationType.unknown
        elif (mismatchCount == 0):
            return MutationType.none
        return MutationType.unknown


    def is_wt(self, bcToAAEntry):
        return self._classify_mutation(bcToAAEntry) is MutationType.none

    def is_single_AA_change(self, bcToAAEntry):
        return self._classify_mutation(bcToAAEntry) is MutationType.substitution

    def single_AA_change_position(self, bcToAAEntry):
        return max(bcToAAEntry['alignment'][1].find(' '),
                bcToAAEntry['alignment'][1].find(':'))

    def insertion_AA_change_position(self, bcToAAEntry):
        gapPositionTop = bcToAAEntry['alignment'][0].find('-')
        gapPositionBottom = bcToAAEntry['alignment'][2].find('-')
        sequenceLength = len(bcToAAEntry['alignment'][0])
        return gapPositionTop


    def locate_AA_change(self, bcToAAEntry):
        changeLocation = self.single_AA_change_position(bcToAAEntry)
        changeAA = bcToAAEntry['alignment'][2][changeLocation]
        return changeLocation, changeAA

    def locate_AA_insertion(self, bcToAAEntry):
        insertLocation = self.insertion_AA_change_position(bcToAAEntry)
        insertAA = bcToAAEntry['alignment'][2][insertLocation]
        return insertLocation, insertAA


    def get_sequence_length(self):
        return self.sequenceLength


    def single_AA_for_barcode(self, barcodeInt):
        matchedEntry = next(
                (x for x in self.bcToAA \
                        if bin_array_to_int(x['barcode']) == barcodeInt))

        if not self.is_single_AA_change(matchedEntry):
            return None

        changeLocation = self.locate_AA_change(matchedEntry)
        return changeLocation[0], aaList.index(changeLocation[1])


    def wt_normalized_measurement_medians_by_AA(self, measurementsUniqueSorted,
            measurementIndex):

        measurementMedians = self.measurement_medians_by_AA(
                measurementsUniqueSorted, measurementIndex)

        wtMeasurements = [x for x in measurementsUniqueSorted \
                if len(x) > 0 and self.is_wt(x[0][3])]
        wtMedian = np.median([y[0][measurementIndex] \
                for x in wtMeasurements for y in x])

        normalizedMeasurements = [[x/wtMedian for x in y] \
                for y in measurementMedians]

        return normalizedMeasurements

    def hamming_distance(self, barcode1, barcode2):
        return np.sum(barcode1 - barcode2 != 0)

    def get_min_hamming_distance_distribution(self):
        barcodes = self.get_barcode_list()

        minDistances = np.zeros(len(barcodes))
        for i in range(len(barcodes)):
            currentDistances = [self.hamming_distance(barcodes[i], y) \
                    for y in barcodes]
            currentDistances[i] = self.bitCount

            minDistances[i] = np.min(currentDistances)

        return minDistances

    def find_closest_entry_to_barcode(self, testBarcode):
        barcodes = self.get_barcode_list()

        return self.bcToAA[np.argmin(
            [self.hamming_distance(testBarcode, x) for x in barcodes])]

    def measurement_medians_by_AA(self, measurementsUniqueSorted, 
            measurementIndex):
        measurementsByAA = self.sort_measurements_by_AA(
                measurementsUniqueSorted)

        '''
        measurementMedians = \
                [[np.median([z[0] for z in x]) \
                if len(x)>0 else float('nan') for x in y] \
                for y in measurementsByAA]
        '''


        return measurementsByAA

    def sort_measurements_by_AA(self, measurementsUniqueSorted):
        measurementsByAA = [[np.nan for x in range(len(aaList))] \
                for y in range(self.sequenceLength)]

        wtMeasurements = 0 
        wtSequence = self.get_wt_sequence()

        for i, currentMeasurement in enumerate(measurementsUniqueSorted):
            currentEntry = currentMeasurement[1][0]
            if self.is_single_AA_change(currentEntry):
                aaChange = self.locate_AA_change(currentEntry)
                measurementsByAA[aaChange[0]][aaList.index(aaChange[1])] = \
                        currentMeasurement[0]
            elif self.is_wt(currentEntry):
                wtMeasurements = currentMeasurement[0]

        if wtMeasurements != 0:
            for i, wtAA in enumerate(wtSequence):
                measurementsByAA[i][aaList.index(wtAA)] = wtMeasurements

        return measurementsByAA


    def get_single_AA_coverage(self):
        aaChange = [self.locate_AA_change(x) for x in self.bcToAA\
                if self.is_single_AA_change(x)]

        countArray = np.zeros((self.sequenceLength, len(aaList)))
        
        for currentChange in aaChange:
            if currentChange[1] in aaList:
                countArray[currentChange[0], aaList.index(currentChange[1])] \
                        += 1
        
        return countArray
        
    def _show_alignment(bcToAAEntry):
        print(bcToAAEntry['alignment'][0])
        print(bcToAAEntry['alignment'][1])
        print(bcToAAEntry['alignment'][2])

    def to_string(self):
        pass

class SequencingDataSet(data.AbstractData):
    pass


'''
class IlluminaSequencingDataSet(SequencingDataSet):

    def __init__(self, read1File, read2File):
        self.read1File = read1File
        self.read2File = read2File

        self.loaded = False

    @classmethod
    def from_experiment(cls, experiment, sampleName, sampleNumber):
        read1Name = '_'.join(
                [sampleName, 'S' + str(sampleNumber), 'L001_R1_001.fastq'])
        read1 = '/'.join([data.__DATAHOME__, experiment, read1Name])

        read2Name = '_'.join(
                [sampleName, 'S' + str(sampleNumber), 'L001_R2_001.fastq'])
        read2 = '/'.join([data.__DATAHOME__, experiment, read2Name])

        return cls(read1, read2)

    def read_iterators(self):
        read1Iterator = SeqIO.parse(self.read1File, 'fastq')
        read2Iterator = SeqIO.parse(self.read2File, 'fastq')
        return read1Iterator, read2Iterator

    def to_string(self):
        pass
'''
