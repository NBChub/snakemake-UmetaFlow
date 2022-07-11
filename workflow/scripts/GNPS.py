import pandas as pd
import numpy as np
import sys

def GNPS_annotations(lib, featurematrix, gnps):
    df= pd.read_csv(lib, sep='\t', encoding='latin-1')
    df.drop(df.index[df['IonMode'] == "negative"], inplace=True)
    df.drop(df.index[df['MZErrorPPM'] > 10.0], inplace=True)
    GNPS=df.filter(["Compound_Name", "RT_Query", "Precursor_MZ"])
    GNPS=GNPS.rename(columns= {"RT_Query": "RetentionTime"})
    GNPS=GNPS.drop_duplicates(subset="Compound_Name", keep='first')

    DF_features= pd.read_csv(featurematrix, sep="\t")


    DF_features.insert(0, 'GNPS_IDs', '')

    for i, mz, rt in zip(DF_features.index, DF_features['mz'], DF_features['RT']):
        hits = []
        for name, GNPS_mz, GNPS_rt, in zip(GNPS['Compound_Name'], GNPS['Precursor_MZ'], GNPS['RetentionTime']):
            mass_delta = (abs(GNPS_mz-mz)/GNPS_mz)*1000000.0 if GNPS_mz != 0 else np.nan
            if (GNPS_rt >= rt-30.0) & (GNPS_rt <= rt+30.0) & (mass_delta<= 10.0):
                hit = f'{name}'
                if hit not in hits:
                    hits.append(hit)
        DF_features['GNPS_IDs'][i] = ' ## '.join(hits)

    DF_features.to_csv(gnps, sep="\t", index = False)
    return DF_features


if __name__ == "__main__":
    GNPS_annotations(sys.argv[1], sys.argv[2], sys.argv[3])