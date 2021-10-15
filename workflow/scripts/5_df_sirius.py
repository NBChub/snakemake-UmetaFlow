from pyopenms import *
from pandas import DataFrame
import pandas as pd
import pyteomics
from pyteomics import mztab

def df_sirius(filename):
    sirius=  pyteomics.mztab.MzTab(filename, encoding='UTF8', table_format='df')
    sirius.metadata
    df= sirius.small_molecule_table
    SIRIUS_DF= df.drop(columns= ["identifier", "smiles", "inchi_key", "description", "database","calc_mass_to_charge", "charge", "taxid", "species", "database_version", "spectra_ref", "search_engine", "modifications"])
#    SIRIUS_DF=SIRIUS_DF[SIRIUS_DF["opt_global_explainedIntensity"] >= 0.6]
    SIRIUS_DF= SIRIUS_DF.sort_values(by= "exp_mass_to_charge")
    SIRIUS_DF.to_csv(snakemake.output[0])


df_sirius(snakemake.input[0])