from pyopenms import *
import glob

# first load feature files in an OpenMS format 
feature_maps = []
for featurexml_file in glob.glob(snakemake.input[0]):
    fmap = FeatureMap()
    FeatureXMLFile().load(featurexml_file, fmap)
    feature_maps.append(fmap)

# reconstruct complete FeatureMaps
to_keep_ids = [item for sublist in [[feature.getUniqueId() for feature in cf.getFeatureList()] for cf in complete] for item in sublist]
for fm in feature_maps:
    fm_filterd = FeatureMap(fm)
    fm_filterd.clear(False)
    for f in fm:
        if f.getUniqueId() in to_keep_ids:
            fm_filterd.push_back(f)
    FeatureXMLFile().store(snakemake.output[0], fm_filterd)