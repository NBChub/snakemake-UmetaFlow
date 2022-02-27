import pandas as pd

consensus = snakemake.input[0]
with open(consensus, 'r') as file:
    for i,line in enumerate(file):
        if '#CONSENSUS' in line:
            header = line.split('\t')
            break

positions = [i for i,col in enumerate(header)]

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

DF_features = pd.DataFrame(thin(), columns= header)
DF_features= DF_features.drop(columns="#CONSENSUS")
df_width=DF_features.filter(regex=fr"width_")
df_qual=DF_features.filter(regex=fr"quality_")
df_rt=DF_features.filter(regex=fr"rt_\d")
df_charge=DF_features.filter(regex=fr"charge_\d")
df_mz=DF_features.filter(regex=fr"mz_\d")

cols_width= df_width.columns
cols_qual= df_qual.columns
cols_rt= df_rt.columns
cols_charge= df_charge.columns
cols_mz= df_mz.columns

DF_features=DF_features.drop(columns="intensity_cf")
DF_features=DF_features.drop(columns=cols_width)
DF_features=DF_features.drop(columns=cols_qual)
DF_features=DF_features.drop(columns=cols_rt)
DF_features=DF_features.drop(columns=cols_charge)
DF_features=DF_features.drop(columns=cols_mz)
DF_features=DF_features.rename(columns={"rt_cf": "RT", "mz_cf": "mz", "charge_cf": "charge"})

DF_features.to_csv(snakemake.output[0], sep="\t", index= None)