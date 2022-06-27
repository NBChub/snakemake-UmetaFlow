from pyopenms import *

# Split the ConsensusMap into features that have no missing values, 
# and features that have at least one missing value; 
# requantify only the missing values. 

# split ConsensusMap
consensus_map = ConsensusMap()
ConsensusXMLFile().load(snakemake.input[0], consensus_map)

headers = consensus_map.getColumnHeaders()

complete = ConsensusMap(consensus_map)
complete.clear(False)
missing = ConsensusMap(consensus_map)
missing.clear(False)

for cf in consensus_map:
    if len(cf.getFeatureList()) < len(headers): #missing values
        missing.push_back(cf)
    else:
        complete.push_back(cf) #no missing values

ConsensusXMLFile().store(snakemake.output[0], complete)
ConsensusXMLFile().store(snakemake.output[1], missing)