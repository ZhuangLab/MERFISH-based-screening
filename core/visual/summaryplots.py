#George Emanuel
#emanuega0@gmail.com
#Copyright Presidents and Fellows of Harvard College, 2017.

import math
import os
import copy

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from core.analysis import bitcalling
from core.analysis import aligner

def normalizer_images(normalizer, savePath=None):
    meanImages = normalizer.get_mean_image()
    imageCount = len(meanImages)

    plt.figure(facecolor='white', figsize=(12,8))
    grid = gridspec.GridSpec(2, (imageCount+1)//2)
    grid.update(wspace=0.01, hspace=0.01)
    for i in range(imageCount): 
        plt.subplot(grid[i])
        plt.imshow(meanImages[i], interpolation='nearest', cmap='hot')
        plt.axis('off')
    plt.subplots_adjust(left=0, right=1, bottom=0, top=1)

    saveName = 'NormalizerImages'
    save_figure(saveName, savePath)

def alignment_offset_histogram(alignments, savePath=None):
    xOff = np.ravel(
            np.array(
                [[y.translation[0] for y in x] for x in alignments]))
    yOff = np.ravel(
            np.array(
                [[y.translation[1] for y in x] for x in alignments]))
    
    plt.figure(facecolor="white", figsize=(4,6))
    plt.subplot(2,1,1)
    plt.hist(xOff, range=(-10,10), bins=30)
    xFraction = sum(np.logical_and(xOff < 5, xOff > -5))/len(xOff)
    plt.title('X offsets (%.3f in (-5,5))' % xFraction)
    plt.ylabel('Counts')
    plt.xlabel('Offset (pixels)')
    plt.subplot(2,1,2)
    plt.hist(yOff, range=(-10,10), bins=30)
    yFraction = sum(np.logical_and(yOff < 5, yOff > -5))/len(yOff)
    plt.title('Y offsets (%.3f in (-5,5))' % yFraction)
    plt.ylabel('Counts')
    plt.xlabel('Offset (pixels)')

    plt.tight_layout()

    saveName = 'AlignmentOffsetHistogram'
    save_figure(saveName, savePath)

def bit_measurement_histograms2d(
        bitMeasurements, savePath=None):

    measurementCount = len(next(iter(bitMeasurements.values())))
    rowCount = min(4, measurementCount)
    columnCount = math.ceil(measurementCount/rowCount)

    plt.figure(facecolor='white')
    for i in range(measurementCount):
        plt.subplot(rowCount, columnCount, i+1)
        bit0 = np.array([x[i][0] for x in bitMeasurements.values()]) 
        bit1 = np.array([x[i][1] for x in bitMeasurements.values()])
        plt.hist2d(np.log(bit0+0.00001), np.log(bit1+0.00001), bins=50) 
        plt.title('Bit ' + str(i+1))
        plt.xlim([-1, 3])
        plt.ylim([-1, 3])
        plt.xticks([0,1,2,3])
        plt.yticks([0,1,2,3])
    plt.tight_layout()    

    saveName = 'BitMeasurementHistogram'
    save_figure(saveName, savePath)

def bit_measurement_histograms2d_with_threshold(
        bitMeasurements, slopes, intercepts, savePath=None):

    measurementCount = len(next(iter(bitMeasurements.values())))
    rowCount = 4
    columnCount = math.ceil(measurementCount/rowCount)

    plt.figure(facecolor='white')
    for i in range(measurementCount):
        plt.subplot(rowCount, columnCount, i+1)
        bit0 = np.array([x[i][0] for x in bitMeasurements.values()]) 
        bit1 = np.array([x[i][1] for x in bitMeasurements.values()])
        plt.hist2d(np.log(bit0+0.00001), np.log(bit1+0.00001), bins=50) 
        plt.title('Bit ' + str(i+1))
        plt.plot([0+intercepts[i], 0+intercepts[i]+4*slopes[i]], [0, 4], 'w')
        plt.xlim([-1, 3])
        plt.ylim([-1, 3])
        plt.xticks([0,1,2,3])
        plt.yticks([0,1,2,3])
    plt.tight_layout()    

    saveName = 'BitMeasurementHistogramWThreshold'
    save_figure(saveName, savePath)

def call_rate_vs_threshold_plot(bitCaller, bitMeasurements, savePath=None):

    thresholdLevels = np.arange(0, 5, 0.1)
    cellCounts = np.zeros(len(thresholdLevels))
    callingRates = np.zeros(len(thresholdLevels))
    for (i,t) in enumerate(thresholdLevels):
        thresholdedBits = {k: v for k,v in bitMeasurements.items() \
                if bitcalling.bit_measurement_passes_threshold(v, t)}
        cellCounts[i] = len(thresholdedBits)
        callingRates[i] = -bitCaller.fraction_matching(
                measurements = thresholdedBits)
    
    fig = plt.figure(facecolor='white', figsize=(5,5))
    rateAxis = fig.add_subplot(111)
    rateLine = rateAxis.plot(thresholdLevels, callingRates, 'b')
    plt.ylabel('Fraction matched', color = 'b', fontname='Arial',
            fontsize=16)
    for tl in rateAxis.get_yticklabels():
        tl.set_color('b')


    countAxis = fig.add_subplot(111, sharex=rateAxis, frameon=False)
    countLine = countAxis.semilogy(thresholdLevels, cellCounts, 'r')
    countAxis.yaxis.tick_right()
    countAxis.yaxis.set_label_position('right')
    plt.ylabel('Cells above threshold', color = 'r',
            fontname='Arial', fontsize=16)
    for tl in countAxis.get_yticklabels():
        tl.set_color('r')

    plt.xlabel('Threshold', fontname='Arial', fontsize=16)

    plt.yticks(fontname='Arial', fontsize=14)
    plt.xticks(fontname='Arial', fontsize=14)

    plt.tight_layout()
    saveName  = 'CallingRatevsThreshold'
    save_figure(saveName, savePath)

    return cellCounts, callingRates

def barcode_abundance_plot(barcodes, nameSuffix = '', savePath=None):
    uniqueBarcodes, counts = np.unique(
            list(barcodes.values()), return_counts=True)

    fig = plt.figure(facecolor='white')
    plt.loglog(range(len(counts)), sorted(counts, reverse=True))
    plt.ylabel('Counts')
    plt.xlabel('Barcode (sorted)')

    plt.tight_layout()
    saveName = 'BarcodeAbundances' + nameSuffix
    save_figure(saveName,savePath)


def cumulative_bit_intensity_plot(bitMeasurements, savePath=None):
    allMeasurements = np.array(list(bitMeasurements.values()))

    reshapedMeasurements = np.reshape(
            allMeasurements,
            (allMeasurements.shape[0], 
                allMeasurements.shape[1]*allMeasurements.shape[2]))
    
    sortedMeasurements = np.sort(reshapedMeasurements, axis=0)        

    plt.figure()
    plt.semilogx(sortedMeasurements, range(sortedMeasurements.shape[0]))
    saveName = 'CumulativeBitIntensities'
    save_figure(saveName, savePath)


def save_figure(saveName, savePath):
    if savePath is not None:
        fullPath = os.sep.join([savePath, saveName])
        plt.savefig(fullPath + '.png', pad_inches=0)
        plt.savefig(fullPath + '.pdf', transparent=True, pad_inches=0)

