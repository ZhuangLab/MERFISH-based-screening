# High-Throughput Sequence-Function Analysis Code

Please check github.com/ZhuangLab for the latest version of this code.

Example data can be found at zhuang.harvard.edu/merfish.html.

-----Python image analysis-----

This analysis software is used for segmenting images into cells, extracting barcodes from each cell, and analyzing cell phenotypes. This analysis software is written in an object-oriented manner and found in the core directory. Analysis is controlled through core\experiment.py, which takes the following parameters:

imagingExperiment - the path to the imaging data in .dax format
sequencingExperiment - the path to the root sequence analysis folder
sequencingLibrary - the directory within the folder sequencingExperiment where the BCtoAA.mat output from sequence analysis is present
parametersFile - a file in json format containing the following parameters:

position_file - the name of a txt file containing a list of microscope stage positions where images were acquired for each imaging round
cell_trace_image - the index (imaging round and frame number) to be used for cell segmentation
images_per_location - a list containing the number of images at each position for each imaging round
alignment_indexes - a list containing the frame number to use to align images from different imaging rounds
phenotype_rounds - the number of imaging rounds that were used to measure the phenotype before barcode readout began
phenotypes - a dictionary of phenotype analysis parameters, see core\analysis\phenotype.py for phenotype analysis options
analysis -  a dictionary of analysis parameters, see core\analysis\analyzer.py for analysis options
bit_threshold - the threshold for median barcode readout intensity for removing barcodes that are too dim
readouts - a list specifying the frame where each bit is imaged. Each element is a list containing [imaging round, frame number in imaging round]
bits - a list specifying the readouts corresponding to each bit. Each element is a list containing [readout 1 index, readout 0 index]

Before analysis, the path to data should be specified in core\data\data.py with the following variables:

__DATAHOME__ - The path to the root imaging data directory
__CACHEHOME__ - The path to save analysis output 
__SEQUENCEHOME__ - The path to the root sequencing data directory
__POSITIONSHOME__ - The path to the root position file directory


------Matlab sequence analysis------

The Matlab sequence analysis for building the barcode to genetic variant lookup table is comprised of three scripts:

UMItoBC.m - Analyze sequencing reads to construct a UMI to barcode lookup table
UMItoAA.m - Analyze sequencing reads to construct a UMI to genetic variant lookup table
AAtoBC.m - Merge the UMI to barcode lookup table with the UMI to genetic variant lookup table by matching UMIs to create the final barcode to genetic variant lookup table

The accessory file AdjacencyMatrix.java is used to help pair UMIs from each sequencing read by determining the adjacency matrix of the UMI graph.

## Code authors
George Emanuel [emanuega0@gmail.com]
