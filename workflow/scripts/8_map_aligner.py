#MapAlignerPoseClustering is used to perform a linear retention time alignment, basically correct for linear shifts in retention time.

from pyopenms import *

def consensus(filename):
    fmap = FeatureMap()
    FeatureXMLFile().load(filename, fmap)
    MapAlignmentAlgorithmPoseClustering.align(fmap)
    FeatureXMLFile().store(snakemake.output[0], fmap)

consensus(snakemake.input[0])