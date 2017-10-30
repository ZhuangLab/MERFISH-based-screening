#George Emanuel
#emanuega0@gmail.com
#Copyright Presidents and Fellows of Harvard College, 2017.

from abc import abstractmethod

import numpy as np
import scipy.io as  sio

from core.data import data

class Cell(data.AbstractData):

    def __init__(self, region, pixels, intensities = None, id = None): 
        self.region = region
        self.pixels = np.array(pixels, np.int16)
        self.id = id 

        self.rawLocation = self.pixels[0]

    def set_raw_offset(self, offset):
        self.rawLocation = [self.pixels[0][0] + offset[0], \
                self.pixels[0][1] + offset[1]]

    def set_intensities(self, intensities):
        self.intensities = intensities

    def get_region(self):
        return self.region

    def get_pixels(self):
        return self.pixels

    def get_raw_offset(self):
        return self.rawLocation

    def get_id(self):
        return self.id

    def to_string(self):
        pass


class CellCollection(data.AbstractData):

    def __init__(self, cellList):
        self.cells = cellList

    def cell_count(self):
        return len(self.cells)

    def cell(self, index):
        return self.cells[index]

    def get_cells(self):
        return self.cells

    def cell_with_id(self, id):
        return next(x for x in self.cells if x.get_id() == id)

    def to_string(self):
        pass







