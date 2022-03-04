import pandas as pd
import numpy as np

FeatureMatrix= snakemake.input[0]


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
DF_features= DF_features.rename(columns={"rt_cf":"RT", "mz_cf":"mz", "intensity_cf": "intensity", "charge_cf":"charge"})
Features_flt=DF_features.filter(regex=fr'(rt_\d+|mz_\d+|quality_\d+|width_\d+|charge_\d+)')
cols= Features_flt.columns
DF_features= DF_features.drop(columns=["#CONSENSUS", "width_cf", "quality_cf"])
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

DF_features= DF_features.reset_index()            
DF_features.to_csv(snakemake.output[0],  sep='\t', index = False)