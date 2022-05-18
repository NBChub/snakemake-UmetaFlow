from pyopenms import *
#remove convex hulls (otherwise Sirius cannot run)
feature_files = sorted(glob.glob(snakemake.input[0]))

feature_maps = []
for file in feature_files:
    fmap = FeatureMap()
    FeatureXMLFile().load(file, fmap)
    feature_maps.append(fmap)

for fmap in feature_maps:
    fm_no_sub = FeatureMap(fmap)
    fm_no_sub.clear(False)
    for f in fmap:
        f.setConvexHulls([])
        f.setSubordinates([])
        fm_no_sub.push_back(f)
    FeatureXMLFile().store(snakemake.output[0], fm_no_sub)