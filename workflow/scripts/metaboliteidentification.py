import pandas as pd
import numpy as np 

FeatureMatrix= snakemake.input[0]
<<<<<<< HEAD
=======


with open(FeatureMatrix, 'r') as file:
    for i,line in enumerate(file):
        if '#CONSENSUS' in line:
            header = line.split('\t')
            break

positions = [i for i,col in enumerate(header)]

def thin():
    with open(FeatureMatrix, 'r') as file:
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
            yield row

with open(FeatureMatrix, 'r') as file:
    for i,line in enumerate(file):
        if '#MAP' in line:
            maps= line.split('\t')
            break

map_IDs = [i for i, col in enumerate(maps)]

def map():
    with open(FeatureMatrix, 'r') as file:
        for i,line in enumerate(file):
            if '#CONSENSUS' in line:
                continue
            if '#' in line:
                continue
            if '#MAP' in line:
                maps=line
            if 'RUN' in line:
                continue
            if 'CONSENSUS' in line:
                continue
            files = line.split('\t')
            files= [files[i] for i in map_IDs]
            yield files

DF_features = pd.DataFrame(thin(), columns=header)
DF_features= DF_features.rename(columns={"rt_cf":"RT", "mz_cf":"mz", "intensity_cf": "intensity", "charge_cf":"charge", "quality_cf":"quality"})
Features_flt=DF_features.filter(regex=fr'(rt_\d+|mz_\d+|quality_\d+|width_\d+|charge_\d+)')
cols= Features_flt.columns
DF_features= DF_features.drop(columns=["#CONSENSUS", "width_cf"])
DF_features= DF_features.drop(columns=cols)
DF_features.columns = DF_features.columns.str.replace(r'intensity_', 'MAP_')

DF_maps= pd.DataFrame(map(), columns=maps)
DF_maps= DF_maps.drop(columns=["label", "size\n"])
DF_maps["#MAP"]= DF_maps["#MAP"]+ "_"+ DF_maps["id"]
DF_maps= DF_maps.drop(columns="id")
DF_maps= DF_maps.set_index("#MAP")
DF_maps["filename"]= DF_maps["filename"].str.replace("results/Interim/mzML/PCpeak_", "")
DF_maps["filename"]= DF_maps["filename"].str.replace(".mzML", "")

DF_features= DF_features.set_index(["RT", "mz"])
cols= DF_features.columns
for col in cols:
    for i, idx in enumerate(DF_maps.index):
        if idx== col:
            col_name= DF_maps["filename"].iloc[i]
            DF_features.rename(columns={col: col_name}, inplace=True)
>>>>>>> 36aac52cfd4423bb92386ddfb3059a9655b26263


with open(FeatureMatrix, 'r') as file:
    for i,line in enumerate(file):
        if '#CONSENSUS' in line:
            header = line.split('\t')
            break

positions = [i for i,col in enumerate(header)]

def thin():
    with open(FeatureMatrix, 'r') as file:
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
            yield row


DF_features = pd.DataFrame(thin(), columns=header)
DF_features = DF_features[["rt_cf","mz_cf", "charge_cf"]]
DF_features= DF_features.rename(columns={ "charge_cf":"Charge", "mz_cf": "Mass", "rt_cf": "RetentionTime"})
DF_features["Charge"]= DF_features["Charge"].astype(str)
DF_features["Mass"]= DF_features["Mass"].astype(float)

for ind in DF_features.index:
    if DF_features["Charge"][ind] == "0":
        DF_features.loc[ind, "Mass"]= DF_features.loc[ind,"Mass"]- 1.007825
    if DF_features["Charge"][ind] == "1":
        DF_features.loc[ind, "Mass"]= DF_features.loc[ind,"Mass"]- 1.007825
    if DF_features["Charge"][ind] == "2":
        DF_features.loc[ind, "Mass"]= (DF_features.loc[ind,"Mass"]*2)- 2.015650
    if DF_features["Charge"][ind] == "3":
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

DF_features["CompoundName"] = np.arange(len(DF_features))
DF_features["CompoundName"] = "feature_" + DF_features["CompoundName"].astype(str)
DF_features["SumFormula"] = " "
DF_features["RetentionTimeRange"]= "0"
DF_features["IsoDistribution"]= "0"
DF_features= DF_features[["CompoundName","SumFormula", "Mass","Charge","RetentionTime","RetentionTimeRange", "IsoDistribution"]]
DF_features.to_csv(snakemake.output[0], sep="\t", index= None)
