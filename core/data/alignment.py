#George Emanuel
#emanuega0@gmail.com
#Copyright Presidents and Fellows of Harvard College, 2017.

import numpy as np

from core.data.data import AbstractData

class ImageAlignment(AbstractData):

    def __init__(self, transformationList, errorList):
        self.transformationList = transformationList
        self.errorList = errorList

    def max_offset(self):
        xOffsets = [t.translation[0] for t in self.transformationList]
        yOffsets = [t.translation[1] for t in self.transformationList]

        return np.max(xOffsets), np.max(yOffsets)

    def min_offset(self):
        xOffsets = [t.translation[0] for t in self.transformationList]
        yOffsets = [t.translation[1] for t in self.transformationList]

        return np.min(xOffsets), np.min(yOffsets)

    def transformation_for_round(self, roundIndex):
        return self.transformationList[roundIndex]

    def round_count(self):
        return len(self.transformationList)

    def load(self):
        pass

    def save(self):
        pass

    def to_string(self):
        pass

    def __iter__(self):
        return self.transformationList.__iter__()

class ImageSetAlignment(AbstractData):

    def __init__(self, transformationList):
        self.transformationList = transformationList

    def transformations_for_position(self, positionIndex):
        return self.transformationList[positionIndex]

    def positions_count(self):
        return len(self.transformationList)

    def __iter__(self):
        return self.transformationList.values().__iter__()

    def load(self):
        pass

    def save(self):
        pass

    def to_string(self):
        pass

