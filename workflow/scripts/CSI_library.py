import glob
import pandas as pd
import numpy as np

input_structures= glob.glob(snakemake.input[0])
DF_CSI= []
for i, formulas in enumerate(input_structures):
    df= pd.read_csv(formulas, index_col="Unnamed: 0")
    df= df.loc[df["opt_global_rank"]==1]
    df_rank= df.filter(regex=fr"opt_global_rank")
    df_score=df.filter(regex=fr"best_search_engine_score")
    df_opt=df.filter(regex=fr"opt")
    cols_score= df_score.columns
    cols_opt= df_opt.columns
    for i, row in df_rank.iterrows():
        if row.sum()>=2:
            df.loc[[i],row.index] = np.nan
    df= df.dropna()
    df= df.drop(columns=cols_score)
    df= df.drop(columns= cols_opt)
    df= df.drop(columns= "identifier")
    df=df.reset_index()
    df= df.drop(columns="index")
    DF_CSI.append(df)


df_structures= pd.concat(DF_CSI, axis=0).sort_values("chemical_formula")
df_structures = df_structures.drop_duplicates(subset=['inchi_key'], keep='first')
df_structures= df_structures.drop(columns=["inchi_key"]) #leave smiles for visualisationdf_structures= df_structures.rename(columns={"chemical_formula": "formulas", "exp_mass_to_charge": "mz", "retention_time": "RT"})
df_structures= df_structures.rename(columns={"chemical_formula":"formulas"})
df_structures_helper= df_structures.copy(deep=True)
df_structures= df_structures.set_index("formulas")
df_singletons=df_structures.reset_index().drop_duplicates(subset="formulas", keep=False)
df_singletons= df_singletons.set_index("formulas")
idx= df_singletons.index
df_CSI= df_structures.drop(labels=idx, axis=0)
new_df= pd.DataFrame()
df= pd.DataFrame()
idx= df_CSI.index
for i, index in enumerate(idx):
    new_index= new_df.index
    if index not in new_index:
        s= df_CSI.iloc[i]
        new_df= new_df.append(s)
    else:
        #print(index)
        mz_0= df_CSI["exp_mass_to_charge"][i]
        mz_1= new_df["exp_mass_to_charge"][index]
        time_0= df_CSI["retention_time"][i]
        time_1= new_df["retention_time"][index]
        #(print(mz_0, time_0, mz_1, time_1))
        mass_delta = (abs(mz_0 - mz_1)/mz_0)*1000000
        maxdeltaRT = time_0 + 30.0
        mindeltaRT = time_0 - 30.0
        if (mindeltaRT<= time_1 <= maxdeltaRT) & (mass_delta<= 20.0):
            pass
        else:
            m= df_CSI.iloc[i]
            df= df.append(m)


DF_CSI= pd.concat([new_df, df], axis=0)
DF_CSI_final= pd.concat([DF_CSI, df_singletons], axis=0)
DF_CSI_final= DF_CSI_final.reset_index()
DF_CSI_final= DF_CSI_final.rename(columns={"index":"formulas"})
DF_CSI_final.to_csv(snakemake.output[0], sep="\t", index= None)