import pandas as pd

DF_features=pd.read_csv(snakemake.input[0], sep="\t")

DF_features= DF_features.rename(columns={ "charge":"Charge", "mz": "Mass", "RT": "RetentionTime"})

for ind in DF_features.index:
    if DF_features["Charge"][ind] == 0:
        DF_features.loc[ind, "Mass"]= DF_features.loc[ind,"Mass"]- 1.007825
    if DF_features["Charge"][ind] == 1:
        DF_features.loc[ind, "Mass"]= DF_features.loc[ind,"Mass"]- 1.007825
    if DF_features["Charge"][ind] == 2:
        DF_features.loc[ind, "Mass"]= (DF_features.loc[ind,"Mass"]*2)- 2.015650
    if DF_features["Charge"][ind] == 3:
        DF_features.loc[ind, "Mass"]= (DF_features.loc[ind,"Mass"]*3)- 3.023475

DF_features["Charge"]= DF_features["Charge"].astype(str)
for ind in DF_features.index:
    if DF_features["Charge"][ind] == "0":
        DF_features.loc[ind, "Charge"]= "+1"
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
DF_features.to_csv(snakemake.output[0], sep="\t", index= None)
