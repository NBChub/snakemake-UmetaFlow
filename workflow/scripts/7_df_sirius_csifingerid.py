from pyopenms import *
from pandas import DataFrame
import pandas as pd
import pyteomics
from pyteomics.openms import featurexml
import numpy as np
import sys
from pyteomics import mztab

def df_sirius_csi(filename):
    sirius=  pyteomics.mztab.MzTab(filename, encoding='UTF8', table_format='df')
    sirius.metadata
    df= sirius.small_molecule_table
    SIRIUS_DF= df.drop(columns= ["identifier", "smiles", "inchi_key", "description", "calc_mass_to_charge", "charge", "taxid", "species","database", "database_version", "spectra_ref", "search_engine", "modifications"])
    SIRIUS_DF.to_csv(snakemake.output[0])

    
    CSI=  pyteomics.mztab.MzTab(snakemake.input[1], encoding='UTF8', table_format='df')
    CSI.metadata
    df= CSI.small_molecule_table
    csifingerID= df.drop(columns= ["calc_mass_to_charge", "charge", "taxid", "species","database", "database_version", "spectra_ref", "search_engine", "modifications"])
    csifingerID.to_csv(snakemake.output[1])

df_sirius_csi(snakemake.input[0])