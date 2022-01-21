import pandas as pd
consensus = "results/GNPSexport/interim/consensus.csv"
with open(consensus, 'r') as file:
    for i,line in enumerate(file):
        if '#CONSENSUS' in line:
            header = line.split('\t')
            break

HEADERS = ["rt_cf", "mz_cf", "charge_cf"]
positions = [i for i,col in enumerate(header) if col in HEADERS]

def thin():
    with open(consensus, 'r') as file:
        for i,line in enumerate(file):
            if '#CONSENSUS' in line:
                header = line
            if '#' in line:
                continue
            if 'MAP' in line:
                continue
            if 'RUN' in line:
                continue
            row = line.split('\t')
            row = [row[i] for i in positions]
            yield row #a generator that won't keep the data matrix (lines) in memory but will provide them when needed

DF_features = pd.DataFrame(thin(), columns=HEADERS)
DF_features = DF_features[["rt_cf","mz_cf", "charge_cf"]]
DF_features["charge_cf"] = pd.to_numeric(DF_features["charge_cf"], downcast="integer")
DF_features["mz_cf"] = pd.to_numeric(DF_features["mz_cf"], downcast="float")
DF_features["rt_cf"] = pd.to_numeric(DF_features["rt_cf"], downcast="float")
DF_features= DF_features.rename(columns={"charge_cf": "Charge", "mz_cf": "Mass", "rt_cf": "RetentionTime"})

for ind in DF_features.index:
    if DF_features["Charge"][ind] == 0:
        DF_features.loc[ind, "Mass"]= DF_features.loc[ind,"Mass"]
    if DF_features["Charge"][ind] == 1:
        DF_features.loc[ind, "Mass"]= DF_features.loc[ind,"Mass"]- 1.007825
    if DF_features["Charge"][ind] == 2:
        DF_features.loc[ind, "Mass"]= (DF_features.loc[ind,"Mass"]*2)- 2.015650
    if DF_features["Charge"][ind] == 3:
        DF_features.loc[ind, "Mass"]= (DF_features.loc[ind,"Mass"]*3)- 3.023475

DF_features["Charge"]= DF_features["Charge"].astype(str)
for ind in DF_features.index:
    if DF_features["Charge"][ind] == "1":
        DF_features.loc[ind, "Charge"]= "+" + DF_features.loc[ind,"Charge"]
    if DF_features["Charge"][ind] == "2":
        DF_features.loc[ind, "Charge"]= "+" + DF_features.loc[ind,"Charge"]
    if DF_features["Charge"][ind] == "3":
        DF_features.loc[ind, "Charge"]= "+" + DF_features.loc[ind,"Charge"]

import numpy as np
DF_features['CompoundName'] = np.arange(len(DF_features))
DF_features['CompoundName'] = "feature_" + DF_features['CompoundName'].astype(str)
DF_features["SumFormula"] = " "
DF_features["RetentionTimeRange"]= "0"
DF_features["IsoDistribution"]= "0"
DF_features= DF_features[["CompoundName","SumFormula", "Mass","Charge","RetentionTime","RetentionTimeRange", "IsoDistribution"]]
DF_features.to_csv("resources/MetaboliteIdentification.tsv", sep="\t", index= None)
