#George Emanuel
#emanuega0@gmail.com
#Copyright Presidents and Fellows of Harvard College, 2017.

import numpy as np
from scipy import ndimage as ndi
from skimage import feature
from skimage import measure
from skimage import morphology

from core.data import cells
from core.analysis import analysis

class CellFinder(analysis.AbstractAnalysisTask):

    def __init__(self, cellImage, regionIndex, offset, idCounter,
            method = 'canny', sigma = 2):
        self.cellImage = cellImage
        self.regionIndex = regionIndex
        self.offset = offset
        self.idCounter = idCounter

        self.method =  method
        self.sigma = sigma

        self.complete = False

    def _extract_from_binary(self, binaryImage):
        labeledEdges = measure.label(binaryImage)
        labelProperties = measure.regionprops(labeledEdges)

        goodProperties = [x for x in labelProperties \
                if x.solidity > 0.9 and x.filled_area > 20]

        return [cells.Cell(self.regionIndex, x.coords, 
                    id=self.idCounter.getAndIncrement()) \
                for x in goodProperties]

    def run(self): 
        if 0 in self.cellImage.shape:
            self.cells = []
            self.complete = True
            return 

        filledEdges = ndi.binary_fill_holes(feature.canny(
                self.cellImage, sigma=self.sigma))

        self.cells = self._extract_from_binary(filledEdges)

        filledCells = np.zeros(self.cellImage.shape, dtype=bool)

        for currentCell in self.cells:
            cellPixels = currentCell.get_pixels()
            filledCells[cellPixels[:,0], cellPixels[:,1]] = True

        filledEdges2 = morphology.binary_erosion(
                ndi.binary_fill_holes(
                    morphology.binary_dilation(
                        filledEdges & np.invert(filledCells))))
        
        self.cells += self._extract_from_binary(filledEdges2)

        [cell.set_raw_offset(self.offset) for cell in self.cells]

        self.complete = True

    def result(self):
        if self.complete:
            return self.cells

    def to_string(self):
        pass


class CellIntensityExtractor(analysis.AbstractAnalysisTask):

    def __init__(self, cells, imageStack, hints = None):
        self.cells = cells
        self.imageStack = imageStack 

        self.complete  = False
        self.hints = hints

    def run(self):

        metricFunction = np.mean
        if self.hints is not None:
            if self.hints['intensity_metric'] == 'median':
                metricFunction = np.median
            elif self.hints['intensity_metric'] == 'mean':
                metricFunction = np.mean
            elif self.hints['intensity_metric'] == 'max':
                metricFunction = np.max

        self.intensities = {currentCell.get_id(): \
            [np.array(
                [metricFunction(list(map(
                lambda x: i[x[0],x[1]],
                currentCell.get_pixels()))) for i in iStack]) \
            for iStack in self.imageStack] \
            for currentCell in self.cells}

        self.complete = True

    def result(self):
        if self.complete:
            return self.intensities

    def to_string(self):
        pass



