from pyopenms import *
import os
import glob

for complete_map in sorted(glob.glob(snakemake.input[0])):
    for requant_map in sorted(glob.glob(snakemake.input[1])):
        if os.path.basename(complete_map)[9:] == os.path.basename(requant_map)[6:]:
            fm_ffm = FeatureMap()
            FeatureXMLFile().load(complete_map, fm_ffm)
            fm_ffmid = FeatureMap()
            FeatureXMLFile().load(requant_map, fm_ffmid)
            for f in fm_ffmid:
                fm_ffm.push_back(f)
            fm_ffm.setUniqueIds()
            FeatureXMLFile().store(snakemake.output[0], fm_ffm)