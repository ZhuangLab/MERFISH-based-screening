#George Emanuel
#emanuega0@gmail.com
#Copyright Presidents and Fellows of Harvard College, 2017.

from abc import ABC, abstractmethod


class AbstractAnalysisTask(ABC):

    @abstractmethod
    def run(self):
        pass

    @abstractmethod
    def result(self):
        pass

    @abstractmethod
    def to_string(self):
        pass


class AbstractMulticoreAnalysisTask(AbstractAnalysisTask):

    @abstractmethod
    def set_core_count(self, coreCount):
        pass

