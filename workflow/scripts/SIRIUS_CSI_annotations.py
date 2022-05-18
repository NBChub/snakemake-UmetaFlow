import pandas as pd

DF_SIRIUS= pd.read_csv(snakemake.input[0], sep="\t")
DF_CSI= pd.read_csv(snakemake.input[1], sep="\t")
DF_features= pd.read_csv(snakemake.input[2], sep="\t")

DF_features.insert(0, 'CSI_predictions_name', '')
DF_features.insert(0, 'CSI_predictions_formula', '')
DF_features.insert(0, 'CSI_predictions_smiles', '')


for i, mz, rt in zip(DF_features.index, DF_features['mz'], DF_features['RT']):
    hits1 = []
    hits2= []
    hits3=[]
    for name, smiles, formula, Pred_mz, Pred_rt, in zip(DF_CSI['description'], DF_CSI['smiles'], DF_CSI['formulas'], DF_CSI['exp_mass_to_charge'], DF_CSI['retention_time']):
        mass_delta = (abs(Pred_mz-mz)/Pred_mz)*1000000.0 if Pred_mz != 0 else 0
        if (Pred_rt >= rt-30.0) & (Pred_rt <= rt+30.0) & (mass_delta<= 10.0):
            hit1 = f'{name}'
            hit2 = f'{formula}'
            hit3= f'{smiles}'
            if hit1 not in hits1:
                hits1.append(hit1)
                hits2.append(hit2)
                hits3.append(hit3)
    DF_features['CSI_predictions_name'][i] = ' ## '.join(hits1)
    DF_features['CSI_predictions_formula'][i] = ' ## '.join(hits2)
    DF_features['CSI_predictions_smiles'][i] = ' ## '.join(hits3)

DF_features.insert(0, 'SIRIUS_predictions', '')

for i, mz, rt in zip(DF_features.index, DF_features['mz'], DF_features['RT']):
    hits = []
    for name, Pred_mz, Pred_rt, in zip(DF_SIRIUS['formulas'], DF_SIRIUS['mz'], DF_SIRIUS['RT']):
        mass_delta = (abs(Pred_mz-mz)/Pred_mz)*1000000.0 if Pred_mz != 0 else 0
        if (Pred_rt >= rt-30.0) & (Pred_rt <= rt+30.0) & (mass_delta<= 10.0):
            hit = f'{name}'
            if hit not in hits:
                hits.append(hit)
    DF_features['SIRIUS_predictions'][i] = ' ## '.join(hits)

DF_features.to_csv(snakemake.output[0], sep="\t", index= None)