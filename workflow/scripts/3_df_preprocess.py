#Explanation of columns
#mz= mass-to-charge ratio (m/z)
#RT= retention time (sec)
#intensity = intensity of the feature (AU-arbitrary units)
#FWHM= Full Width of the peak at Half its Maximum height
#num_of_masstraces = number of mass traces detected (single mass traces are excluded). This is relevant to the isotopic pattern
#isotope_distances = distance in mz between the isotopes (jumps of app. 1 is important to confirm that this is a real feature)

from pyopenms import *
import pandas as pd
from pyteomics.openms import featurexml


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
                    df.loc[idx, 'RT'] = pos['position']
                elif pos['dim'] == '1':
                    df.loc[idx, 'mz'] = pos['position']
        elif key == 'intensity':
            df.loc[idx, 'intensity'] = feat['intensity']
        elif key == 'charge':
            df.loc[idx, 'charge'] = feat['charge']
        elif key == 'FWHM':
            df.loc[idx, 'FWHM'] = feat['FWHM']
        elif key == 'max_height':
            df.loc[idx, 'max_height'] = feat['max_height']
        elif key == 'num_of_masstraces':
            df.loc[idx, 'num_of_masstraces'] = feat['num_of_masstraces']
        elif key == 'masstrace_intesity':
            df.loc[idx, 'masstrace_intesity'] = feat['masstrace_intesity']
        elif key == 'masstrace_centroid_rt':
            df.loc[idx, 'masstrace_centroid_rt'] = feat['masstrace_centroid_rt']
        elif key == 'masstrace_centroid_mz':
            df.loc[idx, 'masstrace_centroid_mz'] = feat['masstrace_centroid_mz']
        elif key == 'isotope_distances':
            df.loc[idx, 'isotope_distances'] = feat['isotope_distances']
        elif key == "adducts":
            if 'adducts' not in df.columns:
                df['adducts'] = ''
                df.loc[idx, "adducts"] = feat["adducts"]
            else:
                df.loc[idx, "adducts"] = feat["adducts"]
df= df[df["num_of_masstraces"]>=2.0]
df=df.reset_index()
df= df.rename(columns={"index": "FeatureID"})
df= df.set_index("FeatureID")
df.to_csv(snakemake.output[0])