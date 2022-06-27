import glob
import pandas as pd
import os

DF_CSI = pd.DataFrame()
list_of_df=[]
for csv in input_formulas:
    df= pd.read_csv(csv, sep=",", index_col="Unnamed: 0")
    s= df["opt_global_rank"]
    pd.to_numeric(s)
    df= df.loc[df["opt_global_rank"]==1]
    df= df.rename(columns={"opt_global_featureId":"featureId"})
    df_score=df.filter(regex=fr"Score")
    df_opt=df.filter(regex=fr"opt")
    cols_score= df_score.columns
    cols_opt= df_opt.columns
    df= df.drop(columns=cols_score)
    df= df.drop(columns= cols_opt)
    df=df.reset_index()
    list_of_df.append(df)
DF_CSI= pd.concat(list_of_df,ignore_index=True)
DF_CSI= DF_CSI.drop(columns="index")
DF_CSI= DF_CSI.rename(columns= {"chemical_formula": "formulas", "exp_mass_to_charge": "mz", "retention_time": "RT"})
DF_CSI["featureId"]= DF_CSI["featureId"].str.replace(r"id_", "")
DF_CSI.to_csv(snakemake.output[0], sep="\t", index= None)