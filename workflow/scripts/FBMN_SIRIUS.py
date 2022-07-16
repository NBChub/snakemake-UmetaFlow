# Integrating into Graphml
import requests
import pandas as pd
import networkx as nx
import glob
import os
import sys

def integration(input_graphml, output_graphml):
    G = nx.read_graphml(input_graphml)
    # Adding sirius information
    input_formulas= glob.glob(os.path.join("results", "SiriusCSI", "formulas_*.tsv"))
    DF_SIRIUS = pd.DataFrame()
    list_of_df=[]
    for tsv in input_formulas:
        df= pd.read_csv(tsv, sep="\t", index_col="Unnamed: 0")
        s= df["opt_global_rank"]
        pd.to_numeric(s)
        df= df.loc[df["opt_global_rank"]==1]
        df=df.reset_index()
        list_of_df.append(df)
    DF_SIRIUS= pd.concat(list_of_df,ignore_index=True)
    DF_SIRIUS= DF_SIRIUS.drop(columns="index")
    DF_SIRIUS["opt_global_featureId"]= DF_SIRIUS["opt_global_featureId"].str.replace(r"id_", "")
   
    for result in DF_SIRIUS.to_dict(orient="records"):
        scan = str(result["opt_global_compoundScanNumber"])
        if scan in G:
            G.nodes[scan]["sirius:molecularFormula"] = result["chemical_formula"]
            G.nodes[scan]["sirius:adduct"] = result["opt_global_adduct"]
            G.nodes[scan]["sirius:TreeScore"] = result["TreeScore"]
            G.nodes[scan]["sirius:IsotopeScore"] = result["IsotopeScore"]
            G.nodes[scan]["sirius:explainedPeaks"] = result["opt_global_explainedPeaks"]
            G.nodes[scan]["sirius:explainedIntensity"] = result["opt_global_explainedIntensity"]
            G.nodes[scan]["sirius:explainedPeaks"] = result["opt_global_explainedPeaks"]

    input_structures= glob.glob(os.path.join("results", "SiriusCSI", "structures_*.tsv"))
    DF_CSI = pd.DataFrame()
    list_of_df=[]
    for tsv in input_structures:
        df= pd.read_csv(tsv, sep="\t", index_col="Unnamed: 0")
        s= df["opt_global_rank"]
        pd.to_numeric(s)
        df= df.loc[df["opt_global_rank"]==1]
        df=df.reset_index()
        list_of_df.append(df)
    DF_CSI= pd.concat(list_of_df,ignore_index=True)
    DF_CSI= DF_CSI.drop(columns="index")
    DF_CSI["opt_global_featureId"]= DF_CSI["opt_global_featureId"].str.replace(r"id_", "")

    # Adding CSI:FingerID information
    for result in DF_CSI.to_dict(orient="records"):
        scan = str(result["opt_global_compoundScanNumber"])
        if scan in G:
            G.nodes[scan]["csifingerid:smiles"] = result["smiles"]
            G.nodes[scan]["csifingerid:Confidence_Score"] = result["best_search_engine_score[1]"]
            G.nodes[scan]["csifingerid:dbflags"] = result["opt_global_dbflags"]

    nx.write_graphml(G, output_graphml)

if __name__ == "__main__":
    integration(sys.argv[1], sys.argv[2])