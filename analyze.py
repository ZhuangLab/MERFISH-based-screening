#George Emanuel
#emanuega0@gmail.com
#Copyright Presidents and Fellows of Harvard College, 2017.

#Helper file for analyzing imaging data from the command line

import getopt
import sys

import numpy as np
import matplotlib.pyplot as plt

from core.data import data
import core.data.image as image
import core.analysis.aligner as aligner
from core.analysis import cellfinder
from core.analysis import normalization
from core import experiment


data.__DATAHOME__ = #path to data
data.__CACHEHOME__ = #path to sequencing results

coreCount = -1
imagingExperiment = None
parametersFile = None
positionFile = None
sequencingExperiment = None
sequencingLibrary = None

opts, args = getopt.getopt(
        sys.argv[1:],
        'c:i:s:l:p:',
        ['cores=','imaging_experiment=','sequencing_experiment=',
            'sequencing_library=', 'position_file=', 'parameters_file='])
for opt, arg in opts:
    if opt in ('-c', '--cores'):
        coreCount = int(arg)
    elif opt in ('-i', '--imaging_experiment'):
        imagingExperiment = arg
    elif opt in ('-s', '--sequencing_experiment'):
        sequencingExperiment = arg
    elif opt in ('-l', '--sequencing_library'):
        sequencingLibrary = arg
    elif opt in ('-p', '--parameters_file'):
        parametersFile = arg
    elif opt in ('--position_file'):
        positionFile = arg
    else:
        print(arg)



print("Imaging experiment " + imagingExperiment)

outputExperiment = experiment.SequenceFunctionExperiment(
        imagingExperiment = imagingExperiment,
        sequencingExperiment = sequencingExperiment,
        sequencingLibrary = sequencingLibrary,
        parametersFile = parametersFile)

outputExperiment.coreCount = coreCount
analysis = outputExperiment.get_analysis()
