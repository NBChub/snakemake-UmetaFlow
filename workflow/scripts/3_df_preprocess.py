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
            int_list = feat['intensity']
            df.loc[idx, 'intensity'] = feat['intensity']
        elif key == 'charge':
            charge_list = feat['charge']
            df.loc[idx, 'charge'] = feat['charge']
        elif key == 'FWHM':
            FWHM_list = feat['FWHM']
            df.loc[idx, 'FWHM'] = feat['FWHM']
        elif key == 'max_height':
            max_height_list = feat['max_height']
            df.loc[idx, 'max_height'] = feat['max_height']
        elif key == 'num_of_masstraces':
            num_of_masstraces_list = feat['num_of_masstraces']
            df.loc[idx, 'num_of_masstraces'] = feat['num_of_masstraces']
        elif key == 'masstrace_intesity':
            masstrace_intesity_list = feat['masstrace_intesity']
            df.loc[idx, 'masstrace_intesity'] = feat['masstrace_intesity']
        elif key == 'masstrace_centroid_rt':
            masstrace_intesity_list = feat['masstrace_centroid_rt']
            df.loc[idx, 'masstrace_centroid_rt'] = feat['masstrace_centroid_rt']
        elif key == 'masstrace_centroid_mz':
            masstrace_intesity_list = feat['masstrace_centroid_mz']
            df.loc[idx, 'masstrace_centroid_mz'] = feat['masstrace_centroid_mz']
        elif key == 'isotope_distances':
            masstrace_intesity_list = feat['isotope_distances']
            df.loc[idx, 'isotope_distances'] = feat['isotope_distances']
        elif key == 'dc_charge_adducts':
            dc_charge_adducts_list = feat['dc_charge_adducts']
            df.loc[idx, 'dc_charge_adducts'] = feat['dc_charge_adducts']

df.to_csv(snakemake.output[0])