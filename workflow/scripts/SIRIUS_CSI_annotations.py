import pandas as pd

DF_SIRIUS= pd.read_csv(snakemake.input[0], sep="\t")
DF_CSI= pd.read_csv(snakemake.input[1], sep="\t")
DF_features= pd.read_csv(snakemake.input[2], sep="\t")

DF_features=DF_features.set_index(["mz", "RT"])
DF_features= DF_features.drop(columns=["charge", "quality", "id"])
DF_features= DF_features.fillna(0)
DF_features["feature_ids"]= [ids[1:-1].split(",") for ids in DF_features["feature_ids"]]
DF_features= DF_features.reset_index()

DF_features.insert(0, "CSI_predictions_name", "")
DF_features.insert(0, "CSI_predictions_formula", "")
DF_features.insert(0, "CSI_predictions_smiles", "")


for i, id in zip(DF_features.index, DF_features["id_list"]):
    hits1 = []
    hits2= []
    hits3=[]
    for i, id in zip(DF_features.index, DF_features["feature_ids"]):
        for name, formula, smiles, Pred_id in zip(DF_CSI["name"], DF_CSI["formulas"], DF_CSI["smiles"], DF_CSI["featureId"]): 
            for feature in id:
                if feature in Pred_id:
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
        for feature in id:
            if feature in Pred_id:
                hit = f"{name}"
                if hit not in hits:
                    hits.append(hit)
        DF_features["SIRIUS_predictions"][i] = " ## ".join(hits)

DF_features.to_csv(snakemake.output[0], sep="\t", index= None)