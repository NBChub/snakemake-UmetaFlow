#Explanation of columns
#mz= mass-to-charge ratio (m/z)
#RT= retention time (sec)
#intensity = intensity of the feature (AU-arbitrary units)
#FWHM= Full Width of the peak at Half its Maximum height
#num_of_masstraces = number of mass traces detected (single mass traces are excluded). This is relevant to the isotopic pattern
#isotope_distances = distance in mz between the isotopes (jumps of app. 1 is important to confirm that this is a real feature)

from pyopenms import *
from pandas import DataFrame
import pandas as pd
import pyteomics
from pyteomics.openms import featurexml
import numpy as np
import sys
from pyteomics import mztab

with featurexml.read(snakemake.input[0]) as f:
    features_list = [FXML for FXML in f]
    
df = pd.DataFrame() 
for feat in features_list:
    idx = feat['id']
    for key in feat.keys():
        if key == 'id':
           pass
        elif key == 'position':
            pos_list = feat['position']
            for pos in pos_list:
                if pos['dim'] == '0':
                    df.loc[idx, 'position_0'] = pos['position']
                elif pos['dim'] == '1':
                    df.loc[idx, 'position_1'] = pos['position']
        elif key == 'quality':
            qual_list = feat['quality']
            for qual in qual_list:
                if qual['dim'] == '0':
                    df.loc[idx, 'quality_0'] = qual['quality']
                elif qual['dim'] == '1':
                    df.loc[idx, 'quality_1'] = qual['quality']
        else:
            df.loc[idx, key] = feat[key]
df_tidy = df.rename(columns = {'position_0': 'RT', 'position_1': 'mz'}, inplace = False)
df_tidy=df_tidy.drop(columns= ["quality_0", "quality_1", "overallquality", "label", "legal_isotope_pattern"])
df_tidy.reset_index(drop=True, inplace=True) 
df_tidy.to_csv(snakemake.output[0])