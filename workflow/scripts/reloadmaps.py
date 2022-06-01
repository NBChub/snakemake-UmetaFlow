from pyopenms import *
import glob

# first load feature files in an OpenMS format 
featurexml_files= glob.glob(snakemake.input[0])
feature_maps = []
for featurexml_file in featurexml_files:
    fmap = FeatureMap()
    FeatureXMLFile().load(featurexml_file, fmap)
    feature_maps.append(fmap)

consensus_map = ConsensusMap()
ConsensusXMLFile().load(snakemake.input[1], consensus_map)
to_keep_ids = [item for sublist in [[feature.getUniqueId() for feature in cf.getFeatureList()] for cf in consensus_map] for item in sublist]

for fm in feature_maps:
    fm_filterd = FeatureMap(fm)
    fm_filterd.clear(False)
    for f in fm:
        if f.getUniqueId() in to_keep_ids:
            fm_filterd.push_back(f)
    FeatureXMLFile().store(snakemake.output[0], fm_filterd)