#George Emanuel
#emanuega0@gmail.com
#Copyright Presidents and Fellows of Harvard College, 2017.

from functools import reduce

import numpy as np
import skimage as sk

from core.data import imagenormalize
from core.analysis import analysis


class MinMeanNormalizer(analysis.AbstractAnalysisTask):

    def __init__(self, normalization, normalizationHints = None):
        self.normalization = normalization
        self.hintDict = None

        if normalizationHints is not None:
            self.set_normalization_hints(normalizationHints)

    def run():
        pass

    def result():
        pass

    def set_normalization_hints(self, normalizationHints):
        self.hintDict = {x['imaging_round']: x['hints'] \
                for x in normalizationHints}

    def normalize_image_stack(self, imageStackIn):
        return [self.normalize_image_sequence(stack, i) \
                for i,stack in enumerate(imageStackIn)]

    def normalize_image_sequence(self, imageSequence, imageRound):
        sequenceLength = imageSequence.shape[0]
        normalizationIndexes = np.array(range(sequenceLength))
        if self.hintDict is not None:
            if imageRound in self.hintDict.keys():
                for currentHint in self.hintDict[imageRound]:
                    normalizationIndexes[\
                            currentHint['start_frame']:\
                            currentHint['end_frame']]\
                            = currentHint['normalization_index']

        return np.array([self.normalize_image(image, i) for (i, image) in \
                zip(normalizationIndexes, imageSequence)])

    def normalize_image(self, imageIn, normalizationIndex):
        return (sk.img_as_float(imageIn) \
                - self.normalization.get_min()[normalizationIndex]) \
                / \
                (self.normalization.get_mean_image()[normalizationIndex] \
                - self.normalization.get_min()[normalizationIndex])

    def to_string(self):
        pass
                        

class MinMeanNormalizationGenerator(analysis.AbstractAnalysisTask):

    def __init__(self, dataSet, normalizationRound):
        self.dataSet = dataSet
        self.normalizationRound = normalizationRound

        self.complete = False

    def run(self):

        if not isinstance(self.normalizationRound, list):
            imageList = self.dataSet.images_in_round(self.normalizationRound)

            minImage = sk.img_as_float(
                    reduce(lambda x,y: np.minimum(x, y), imageList))
            sumImage = sk.img_as_float(
                    reduce(lambda x,y: np.add(
                        sk.img_as_float(x), sk.img_as_float(y)), 
                        imageList))

        else:
            imageDimensions = self.dataSet.image_dimensions()
            minImage = sk.img_as_float(np.zeros(
                    (len(self.normalizationRound), 
                        imageDimensions[0], imageDimensions[1])))
            sumImage = sk.img_as_float(np.zeros(
                    (len(self.normalizationRound),
                        imageDimensions[0], imageDimensions[1])))

            for i in range(len(self.normalizationRound)):
                for j in range(self.dataSet.position_count()):
                    imageList = self.dataSet.images_at_position_in_round(
                            j, self.normalizationRound[i])

                    minImage[i] = sk.img_as_float(
                            reduce(lambda x,y: np.minimum(x, y), 
                                [minImage[i], imageList[i]]))

                    sumImage[i] = sk.img_as_float(
                            reduce(lambda x,y: np.add(
                                sk.img_as_float(x), sk.img_as_float(y)), 
                                [sumImage[i], imageList[i]]))



        self.normalizer = imagenormalize.MinMeanImageNormalization(
                sumImage/self.dataSet.position_count(), minImage)

        self.complete = True

    def result(self):
        if self.complete:
            return self.normalizer

    def to_string(self):
        pass






