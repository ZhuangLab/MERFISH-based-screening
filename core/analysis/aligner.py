#George Emanuel
#emanuega0@gmail.com
#Copyright Presidents and Fellows of Harvard College, 2017.

from functools import reduce

import numpy as np
from scipy.signal import correlate2d
from skimage.feature import register_translation
import skimage.transform as tf

from core.data.alignment import ImageAlignment
from core.analysis import analysis


class ImageAligner(analysis.AbstractAnalysisTask):

    def __init__(self, imageSet, imageAlignment):
        self.imageSet = imageSet
        self.imageAlignment = imageAlignment

        self.complete = False

    def run(self):
        alignedImages = [np.array([tf.warp(
                    ci, self.imageAlignment.transformation_for_round(i)) \
                            for ci in image])
                        for i, image in enumerate(self.imageSet)]
        self.alignedImages = alignedImages

        goodPixels = reduce(
                np.bitwise_and, map(lambda x: x[0] != 0, alignedImages))

        xAny = list(map(any, goodPixels))
        x1 = next(i for (i,x) in enumerate(xAny) if x)
        x2 = next(len(xAny)-i-1 for (i,x) in enumerate(reversed(xAny)) if x)

        yAny = list(map(any, np.transpose(goodPixels)))
        y1 = next(i for (i,y) in enumerate(yAny) if y)
        y2 = next(len(yAny)-i-1 for (i,y) in enumerate(reversed(yAny)) if y)

        self.transformedImages = [i[:,x1:x2,y1:y2] for i in alignedImages]

        self.complete = True
        
        return self.transformedImages

    def result(self):
        if self.complete:
            return self.transformedImages

    def to_string(self):
        pass


class ImageAlignmentFinder(analysis.AbstractAnalysisTask):

    def __init__(self, imageSet, registrationIndexes):
        self.imageSet = imageSet
        self.registrationIndexes = registrationIndexes \
                if hasattr(registrationIndexes, '__len__') \
                else [registrationIndexes]*len(imageSet)

        self.complete = False

    def run(self):
        fixedImage = self.imageSet[0][self.registrationIndexes[0]]
        offsetErrors = [register_translation(fixedImage, x[ri], 10) \
                    for ri, x \
                    in zip(self.registrationIndexes[1:], self.imageSet[1:])]

        offsets = [tf.SimilarityTransform()] \
                + [tf.SimilarityTransform(
                        translation=[-offset[1], -offset[0]]) \
                    for offset, error, nothing in offsetErrors]
        errors = [0] + [error for offset, error, nothing in offsetErrors]

        self.alignments = ImageAlignment(offsets, errors)

        self.complete = True

        return self.alignments

    def result(self):
        if self.complete:
            return self.alignments

    def to_string(self):
        pass

