#George Emanuel
#emanuega0@gmail.com
#Copyright Presidents and Fellows of Harvard College, 2017.

import itertools

import numpy as np

class IntensityDict(dict):

    def __init__(self, imagesPerLocation):
        super(IntensityDict, self).__init__()
        self.imagesCumulative = np.cumsum(imagesPerLocation)[:-1]

    def __getitem__(self, key):
        rawItem = dict.__getitem__(self, key)
        return self._convert_value(rawItem)

    def __setitem__(self, key, value):
        if type(value).__module__ == np.__name__:
            dict.__setitem__(self, key, value)
        else:
            dict.__setitem__(
                    self, key, np.array(
                        list(itertools.chain.from_iterable(value))))
    
    def get(self, key):
        return self.__getitem__(key)

    def update(self, *args, **kwargs):
        if args:
            if len(args) > 1:
                raise TypeError('update expected at most 1 arguments, '
                        'got %d' % len(args))
            other = dict(args[0])
            for key in other:
                self[key] = other[key]
        for key in kwargs:
            self[key] = kwargs[key]

    def setdefault(self, key, value=None):
        if key not in self:
            self[key] = value
        return self[key]

    def values(self):
        return IntensityDictValueView(self)

    def items(self):
        return IntensityDictItemView(self)

    def _convert_value(self, value):
        return np.array_split(value, self.imagesCumulative)

class IntensityDictItemView(object):

    def __init__(self, intensityDict):
        self.iDict = intensityDict
        self.keyIterator  = iter(intensityDict.keys())

    def __iter__(self):
        return self

    def __next__(self):
        nextKey =  next(self.keyIterator)
        return nextKey, self.iDict[nextKey]

class IntensityDictValueView(object):

    def __init__(self, intensityDict):
        self.iDict = intensityDict
        self.keyIterator = iter(intensityDict.keys())

    def  __iter__(self):
        return self

    def __next__(self):
        return self.iDict[next(self.keyIterator)]


    
