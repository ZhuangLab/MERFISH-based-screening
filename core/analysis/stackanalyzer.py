#George Emanuel
#emanuega0@gmail.com
#Copyright Presidents and Fellows of Harvard College, 2017.

import itertools
from multiprocessing import Process
import multiprocessing
from multiprocessing import Value
from multiprocessing import Lock
from multiprocessing import Queue
import logging

import numpy as np

from core.data import alignment
from core.data import intensities
from core.data.cells import CellCollection
from core.analysis import analysis
from core.analysis import aligner
from core.analysis import focus
from core.analysis import cellfinder

class Counter(object):

    def __init__(self, manager = None):
        if manager is None:
            self.idValue = multiprocessing.Value('i', 0)
            self.lock = multiprocessing.Lock()
        else:
            self.idValue = manager.Value('i', 0)
            self.lock = manager.Lock()

    def getAndIncrement(self):
        with self.lock:
            returnValue = self.idValue.value
            self.idValue.value += 1
        return returnValue


class ImageStackAnalyzer(analysis.AbstractMulticoreAnalysisTask):

    def __init__(self, imageData, normalizer, positions, parameters):
        self.imageData = imageData
        self.normalizer = normalizer
        self.positionCount = len(positions)
        self.parameters = parameters
        self.coreCount = 1
        manager = multiprocessing.Manager()
        self.idCounter = Counter(manager)

        self.complete = False

    def set_core_count(self, coreCount):
        self.coreCount = coreCount

    def run(self):

        pool = multiprocessing.Pool(self.coreCount) 
        mapResults = pool.map(
                self._analyze_images_by_position, range(self.positionCount))

        self.alignments = alignment.ImageSetAlignment(
                {i: x[0] for i,x in enumerate(mapResults) \
                        if x is not None and x[0] is not None})

        cells = list(
                itertools.chain.from_iterable(
                    (x[1] for x in mapResults \
                            if x is not None and x[1] is not None)))
        self.cellList = CellCollection(cells)

        self.intensities = intensities.IntensityDict(
            self.imageData.images_per_location())
        for i in iter(mapResults):
            if i is not None and i[2] is not None: 
                self.intensities.update(i[2])

        self.complete = True

        logging.info('Finished analyzing image stacks. ' + 
                str(self.cellList.cell_count()) + ' cells found. ' + 
                str(sum([x is None or x[0] is None for x in mapResults])) \
                        + ' of ' + 
                str(len(mapResults)) + ' regions were rejected for bad focus.')

    def get_alignments(self):
        if self.complete:
            return self.alignments

    def get_cells(self):
        if self.complete:
            return self.cellList

    def get_intensities(self):
        if self.complete:
            return self.intensities
    
    def _analyze_images_by_position(self, positionIndex):

        stackParameters = None
        if 'stack_hints' in self.parameters:
            stackParameters = self.parameters['stack_hints']

        imageStack = self.normalizer.normalize_image_stack(
                self.imageData.images_at_position(positionIndex))

        alignmentIndexes = self.parameters['alignment_indexes']
        if not hasattr(alignmentIndexes, '__len__'):
            alignmentIndexes = [alignmentIndexes] * len(imageStack)

        if not all(focus.is_focused(image[index]) for image,index \
                in zip(imageStack, alignmentIndexes)):
            return None, None, None

        alignFinder = aligner.ImageAlignmentFinder(
                imageStack, self.parameters['alignment_indexes'])
        alignFinder.run()

        imageAligner = aligner.ImageAligner(
                imageStack, alignFinder.result())
        imageAligner.run()

        alignedImages = imageAligner.result()

        if alignedImages[0].shape[2] ==  1 or alignedImages[0].shape[1] ==  1:
            return None, None, None


        traceIndex = self.parameters['cell_trace_image']
        offset = np.abs(alignFinder.result().min_offset())
        cellFinder = cellfinder.CellFinder(
                alignedImages[traceIndex[0]][traceIndex[1]], positionIndex,
                (offset[1], offset[0]), self.idCounter)
        cellFinder.run()



        currentCells = cellFinder.result()
        intensityExtractor = cellfinder.CellIntensityExtractor(
                currentCells, alignedImages, hints = stackParameters)
        intensityExtractor.run()

        return alignFinder.result(), cellFinder.result(), \
                intensityExtractor.result()

    def result(self):
        pass

    def to_string(self):
        pass
    
