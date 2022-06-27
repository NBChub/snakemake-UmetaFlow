import pandas as pd
import numpy as np
from pyopenms import *

consensus_map = ConsensusMap()
ConsensusXMLFile().load(snakemake.input[0], consensus_map)
df = consensus_map.get_df()
for cf in consensus_map:
    if cf.metaValueExists("best ion"):
        df["adduct"] = [cf.getMetaValue("best ion") for cf in consensus_map]
        break
df["feature_ids"] = [[handle.getUniqueId() for handle in cf.getFeatureList()] for cf in consensus_map]
df= df.reset_index()
df= df.drop(columns= ["sequence"])

df.to_csv(snakemake.output[0], sep="\t", index = False)