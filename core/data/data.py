#George Emanuel
#emanuega0@gmail.com
#Copyright Presidents and Fellows of Harvard College, 2017.

from abc import ABC, abstractmethod

global __DATAHOME__ 
__DATAHOME__ = "//data_home"
global  __CACHEHOME__
__CACHEHOME__ = "//sequencing_home/"
global __POSITIONSHOME__
__POSITIONSHOME__ = "//10.245.74.90/analysis/htseq_fun/stage_positions/flow chamber"
global __SEQUENCEHOME__
__SEQUENCEHOME__ = "//10.245.74.90/analysis/htseq_fun/cache"


class AbstractRawData(ABC):

    @abstractmethod
    def to_string(self):
        pass


class AbstractData(ABC):

    @abstractmethod
    def to_string(self):
        pass
