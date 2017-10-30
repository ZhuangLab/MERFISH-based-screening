#George Emanuel
#emanuega0@gmail.com
#Copyright Presidents and Fellows of Harvard College, 2017.

from abc import abstractmethod
from os import listdir
import logging
import sys

import numpy as np

from core.data import data


class ImagingDataSet(data.AbstractRawData):

    @abstractmethod
    def imaging_round_count(self):
        pass

    @abstractmethod
    def images_per_location(self):
        pass


class DaxImagingDataSet(ImagingDataSet):

    def __init__(self, daxPaths, framesPerLocation):
        self.daxFiles = [DaxFile(f) for f in sorted(daxPaths)]
        self.framesPerLocation = \
                framesPerLocation \
                if hasattr(framesPerLocation, '__len__') \
                else [framesPerLocation]*len(self.daxFiles)

    @classmethod
    def from_experiment(cls, experiment, framesPerLocation):
        basePath = '/'.join([data.__DATAHOME__, experiment])
        daxFiles = sorted([f for f in listdir(basePath) if f.endswith('dax')])

        return cls(
                ['/'.join([basePath, f]) for f in daxFiles], 
                framesPerLocation)

    def _start_frame(self, positionIndex, roundNumber):
        return self.framesPerLocation[roundNumber]*positionIndex

    def _end_frame(self, positionIndex, roundNumber):
        return self.framesPerLocation[roundNumber]*(positionIndex+1)

    def imaging_round_count(self):
        return len(self.daxFiles)

    def images_at_position_in_round(self, positionIndex, roundNumber):
        return self.daxFiles[roundNumber].image_data(
                startFrame = self._start_frame(positionIndex, roundNumber),
                endFrame = self._end_frame(positionIndex, roundNumber))

    def images_at_position(self, positionIndex):
        return [np.array(self.images_at_position_in_round(positionIndex, i)) 
                for  i in range(len(self.daxFiles))]
        
    def images_in_round(self, roundNumber):
        fullData = np.array(self.daxFiles[roundNumber].image_data())
        positionCount = int(len(fullData)/self.framesPerLocation[roundNumber])
        return np.reshape(fullData, 
                (positionCount, self.framesPerLocation[roundNumber],
                    fullData.shape[1], fullData.shape[2]))

    def image_dimensions(self):
        return self.daxFiles[0].dimensions()

    def images_per_location(self):
        return self.framesPerLocation

    def position_count(self):
        return int(
                self.daxFiles[0].frame_count()/self.images_per_location()[0])

    def to_string(self):
        pass
        

class AlignedImageStack(ImagingDataSet):

    def __init__(self, imageData):
        self.imageData = imageData

    def imaging_round_count(self):
        return len(imageData)

    def images_per_round(self):
        pass



class ImageFile(data.AbstractRawData):

    @abstractmethod
    def frame_count(self):
        pass

    @abstractmethod
    def dimensions(self):
        pass

    @abstractmethod
    def image_data(self):
        pass


class DaxFile(ImageFile):

    def __init__(self, daxPath):
        self.daxPath = daxPath
        self.infoPath = daxPath.rsplit('.', 1)[0] + '.inf'

        with open(self.infoPath, 'r') as infoFile:
            self.infoDictionary = dict(
                [line.split(' = ') for line in infoFile \
                        if len(line.split(' = ')) is 2])

    @classmethod
    def from_experiment(cls, experiment, daxIndex):
        basePath = '/'.join([data.__DATAHOME__, experiment])
        daxFiles = sorted([f for f in listdir(basePath) if f.endswith('.dax')])

        return cls('/'.join([basePath, daxFiles[daxIndex]]))

    def frame_count(self):
        return int(self.infoDictionary['number of frames'])

    def dimensions(self):
        return list(map(
            int, self.infoDictionary['frame dimensions'].split(' x ')))

    def frame_size(self):
        return int(self.infoDictionary['frame size'])

    def image_data(self, startFrame = 0, endFrame = -1):
        if (endFrame == -1):
            endFrame = self.frame_count()

        frameSize = self.frame_size()
        framesToRead = endFrame - startFrame

        with open(self.daxPath, 'rb') as daxFile:
            daxFile.seek(2*startFrame*frameSize)
            imageData = [np.fromfile(daxFile, dtype='>u2', 
                count=frameSize) \
                    .reshape(
                        self.dimensions()) for i in range(framesToRead)]

        return imageData

    def to_string(self):
        pass


        



            





