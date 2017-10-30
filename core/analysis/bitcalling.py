#George Emanuel
#emanuega0@gmail.com
#Copyright Presidents and Fellows of Harvard College, 2017.

from functools import reduce
import random

import numpy as np
from scipy import optimize


def bit_measurement_passes_threshold(measurement, threshold):
    return np.min(np.max(measurement, axis=1)) > threshold


def bin_array_to_int(binaryArray):
    return binaryArray.dot(1 << np.arange(binaryArray.shape[-1] -1, -1, -1))


class BitCaller(object):

    def __init__(self, bitMeasurements, knownBarcodes = None):
        self.bitMeasurements = bitMeasurements

        if knownBarcodes is not None:
            self.knownBarcodes = np.array(knownBarcodes)
            self.knownInts = bin_array_to_int(self.knownBarcodes)

        self.barcodesToOptimize = None

        self.bitCount = len(next(iter(bitMeasurements.values())))
        self.slopes = [1]*self.bitCount
        self.intercepts = [0]*self.bitCount

        self.measuredBarcodes = None

        self.complete = False

    def run(self):
        self.optimize_slopes()
        self.measuredBarcodes = self._call_barcodes()
        self.complete = True

    def result(self):
        return self.measuredBarcodes

    def get_optimal_parameters(self):
        return [self.slopes, self.intercepts]

    def get_slopes(self):
        return self.slopes

    def get_intercepts(self):
        return self.intercepts

    def optimize_slopes(self):
        if self.knownBarcodes is None:
            return

        self.barcodesToOptimize = self.knownInts;
        if len(self.barcodesToOptimize) > 150:
            self.barcodesToOptimize = random.sample(
                    list(self.barcodesToOptimize), 150)

        optimizedValues = optimize.fmin(
                func = self.fraction_matching, 
                x0 = np.concatenate((self.slopes, self.intercepts)))

        self.slopes = optimizedValues[0:self.bitCount]
        self.intercepts = optimizedValues[self.bitCount:]

        self.measuredBarcodes = None
        self.barcodesToOptimize = None

    def fraction_matching(self, slopeIntercepts = None, measurements = None): 
        if measurements is None: measurements = self.bitMeasurements
        matchedInts = self.barcodesToOptimize \
                if self.barcodesToOptimize is not None \
                else self.knownInts

        if slopeIntercepts is None:
            slopes = self.slopes
            intercepts = self.intercepts
        else:
            slopes = slopeIntercepts[0:self.bitCount]
            intercepts = slopeIntercepts[self.bitCount:]

        calledBarcodes = self._call_barcodes(
                slopes, intercepts, measurements=measurements)
        if len(calledBarcodes) == 0:
            return 0

        matchedCount = reduce(
                np.add, 
                (1 for b in calledBarcodes.values() if \
                        bin_array_to_int(b) in matchedInts), 0) 

        return -matchedCount/len(calledBarcodes)

    def _call_barcodes(self, slopes = None, intercepts = None,
            measurements = None):
        if slopes is None: slopes = self.slopes
        if intercepts is None: intercepts = self.intercepts
        if measurements is None: measurements = self.bitMeasurements

        return {index:  self._bits_to_barcode(
                        bits, slopes = slopes, intercepts = intercepts) \
                    for index,bits in measurements.items() }

    def _bits_to_barcode(self, bits, slopes = None, intercepts = None):
        if slopes is None:
            slopes = self.slopes
        if intercepts is None:
            intercepts = self.intercepts

        return bits[:,0] > slopes*bits[:,1] 

    def get_barcodes(self):
        if self.measuredBarcodes is None:
            self.measuredBarcodes = self._call_barcodes()
                        
        return self.measuredBarcodes

    def get_matched_barcodes(self):
        pass




