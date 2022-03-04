import glob
import pandas as pd
import numpy as np

input_formulas= glob.glob(snakemake.input[0])

DF_SIRIUS= []
for i, formulas in enumerate(input_formulas):
    df= pd.read_csv(formulas, index_col="Unnamed: 0")
    df= df.loc[df["opt_global_rank"]==1]
    df_rank= df.filter(regex=fr"opt_global_rank")
    df_score=df.filter(regex=fr"Score")
    df_opt=df.filter(regex=fr"opt")
    cols_score= df_score.columns
    cols_opt= df_opt.columns
    df= df.drop(columns=cols_score)
    df= df.drop(columns= cols_opt)
    df=df.reset_index()
    df= df.drop(columns="index")
    DF_SIRIUS.append(df)
    
df_formulas= pd.concat(DF_SIRIUS, axis=0).sort_values("chemical_formula")
df_formulas= df_formulas.rename(columns= {"chemical_formula": "formulas", "exp_mass_to_charge": "mz", "retention_time": "RT"})
df_formulas_helper= df_formulas.copy(deep=True)
df_formulas= df_formulas.set_index("formulas")

df_singletons=df_formulas.reset_index().drop_duplicates(subset="formulas", keep=False)
df_singletons= df_singletons.set_index("formulas")
idx= df_singletons.index
df_sirius= df_formulas.drop(labels=idx, axis=0)

new_df= pd.DataFrame()
df= pd.DataFrame()
idx= df_sirius.index
for i, index in enumerate(idx):
    new_index= new_df.index
    if index not in new_index:
        s= df_sirius.iloc[i]
        new_df= new_df.append(s)
    else:
        #print(index)
        mz_0= df_sirius["mz"][i]
        mz_1= new_df["mz"][index]
        time_0= df_sirius["RT"][i]
        time_1= new_df["RT"][index]
        #(print(mz_0, time_0, mz_1, time_1))
        mass_delta = (abs(mz_0 - mz_1)/mz_0)*1000000
        maxdeltaRT = time_0 + 30.0
        mindeltaRT = time_0 - 30.0
        if (mindeltaRT<= time_1 <= maxdeltaRT) & (mass_delta<= 20.0):
            pass
        else:
            m= df_sirius.iloc[i]
            df= df.append(m)


DF_SIRIUS= pd.concat([new_df, df], axis=0)
DF_SIRIUS_final= pd.concat([DF_SIRIUS, df_singletons], axis=0)
DF_SIRIUS_final= DF_SIRIUS_final.reset_index()
DF_SIRIUS_final= DF_SIRIUS_final.rename(columns={"index":"formulas"})
DF_SIRIUS_final.to_csv(snakemake.output[0], sep="\t")
