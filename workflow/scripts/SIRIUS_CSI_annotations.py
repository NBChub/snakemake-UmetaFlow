import pandas as pd
import glob
import os

<<<<<<< HEAD
input_formulas= glob.glob(os.path.join("results", "SiriusCSI", "formulas_*.tsv"))
=======
input_formulas= glob.glob(os.path.join("results", "SIRIUS", "formulas_*.tsv"))
>>>>>>> 8c5ccedc94fc856850dde76e921e22a705428575
DF_SIRIUS = pd.DataFrame()
list_of_df=[]
for csv in input_formulas:
    df= pd.read_csv(csv, sep="\t", index_col="Unnamed: 0")
    s= df["opt_global_rank"]
    pd.to_numeric(s)
    df= df.loc[df["opt_global_rank"]==1]
    df= df.rename(columns={"opt_global_featureId":"featureId"})
    df= df.drop(columns=df.filter(regex=fr"Score").columns)
    df= df.drop(columns= df.filter(regex=fr"opt").columns)
    df=df.reset_index()
    list_of_df.append(df)
DF_SIRIUS= pd.concat(list_of_df,ignore_index=True)
DF_SIRIUS= DF_SIRIUS.drop(columns="index")
DF_SIRIUS= DF_SIRIUS.rename(columns= {"chemical_formula": "formulas", "exp_mass_to_charge": "mz", "retention_time": "RT"})
DF_SIRIUS["featureId"]= DF_SIRIUS["featureId"].str.replace(r"id_", "")
for i, rows in DF_SIRIUS.iterrows():
    DF_SIRIUS["featureId"][i]= DF_SIRIUS["featureId"][i].split(",")

<<<<<<< HEAD
input_structures= glob.glob(os.path.join("results", "SiriusCSI", "structures_*.tsv"))
=======
input_structures= glob.glob(os.path.join("results", "CSI", "structures_*.tsv"))
>>>>>>> 8c5ccedc94fc856850dde76e921e22a705428575
DF_CSI = pd.DataFrame()
list_of_df=[]
for csv in input_structures:
    df= pd.read_csv(csv, sep="\t", index_col="Unnamed: 0")
    s= df["opt_global_rank"]
    pd.to_numeric(s)
    df= df.loc[df["opt_global_rank"]==1]
    df= df.rename(columns={"opt_global_featureId":"featureId"})
    df= df.drop(columns=df.filter(regex=fr"Score").columns)
    df= df.drop(columns= df.filter(regex=fr"opt").columns)
    df=df.reset_index()
    list_of_df.append(df)
DF_CSI= pd.concat(list_of_df,ignore_index=True)
DF_CSI= DF_CSI.drop(columns="index")
DF_CSI= DF_CSI.rename(columns= {"chemical_formula": "formulas", "exp_mass_to_charge": "mz", "retention_time": "RT", "description":"name"})
DF_CSI["featureId"]= DF_CSI["featureId"].str.replace(r"id_", "")
for i, rows in DF_CSI.iterrows():
    DF_CSI["featureId"][i]= DF_CSI["featureId"][i].split(",")

DF_features= pd.read_csv(snakemake.input[0], sep="\t")
DF_features=DF_features.set_index(["mz", "RT"])
DF_features= DF_features.drop(columns=["charge", "quality", "id"])
DF_features= DF_features.fillna(0)
DF_features["feature_ids"]= [ids[1:-1].split(",") for ids in DF_features["feature_ids"]]
DF_features= DF_features.reset_index()

DF_features.insert(0, "CSI_predictions_name", "")
DF_features.insert(0, "CSI_predictions_formula", "")
DF_features.insert(0, "CSI_predictions_smiles", "")

for i, id in zip(DF_features.index, DF_features["feature_ids"]):
    hits1 = []
    hits2= []
    hits3=[]
    for name, formula, smiles, Pred_id in zip(DF_CSI["name"], DF_CSI["formulas"], DF_CSI["smiles"], DF_CSI["featureId"]): 
        for x,y in zip(id,Pred_id):
            if x==y:
                hit1 = f"{name}"
                hit2 = f"{formula}"
                hit3= f"{smiles}"
                if hit1 not in hits1:
                    hits1.append(hit1)
                    hits2.append(hit2)
                    hits3.append(hit3)
    DF_features["CSI_predictions_name"][i] = " ## ".join(hits1)
    DF_features["CSI_predictions_formula"][i] = " ## ".join(hits2)
    DF_features["CSI_predictions_smiles"][i] = " ## ".join(hits3)

DF_features.insert(0, "SIRIUS_predictions", "")

for i, id in zip(DF_features.index, DF_features["feature_ids"]):
    hits = []
    for name, Pred_id in zip(DF_SIRIUS["formulas"], DF_SIRIUS["featureId"]): 
        for x,y in zip(id,Pred_id):
            if x==y:
                hit = f"{name}"
                if hit not in hits:
                    hits.append(hit)
    DF_features["SIRIUS_predictions"][i] = " ## ".join(hits)

DF_features.to_csv(snakemake.output[0], sep="\t", index= None)