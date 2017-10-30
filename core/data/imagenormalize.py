#George Emanuel
#emanuega0@gmail.com
#Copyright Presidents and Fellows of Harvard College, 2017.

from abc import abstractmethod

import numpy as np

from core.data import data

class ImageNormalization(data.AbstractData):

    @abstractmethod
    def normalization_count(self, i):
        pass


class MinMeanImageNormalization(ImageNormalization):

    def __init__(self, meanImage, minImage):
        self.meanImage = meanImage
        self.minImage = minImage
        self.mins = np.amin(minImage, axis=(1,2))

    def get_mean_image(self):
        return self.meanImage

    def get_min_image(self):
        return self.minImage

    def get_min(self):
        return self.mins

    def normalization_count(self):
        return len(self.meanImage)

    def load(self):
        pass

    def save(self):
        pass

    def to_string(self):
        pass

