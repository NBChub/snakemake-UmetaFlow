import os
import fnmatch
import pandas as pd
from snakemake.utils import validate
from snakemake.utils import min_version
import glob

min_version("5.18.0")

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
#singularity: "docker://continuumio/miniconda3"

##### load config and sample sheets #####

configfile: "config/config.yaml"
validate(config, schema="../schemas/config.schema.yaml")

# set up sample
samples = pd.read_csv(config["samples"], sep="\t").set_index("sample_name", drop=False)
samples.index.names = ["samples"]


##### Wildcard constraints #####
wildcard_constraints:
    sample="|".join(samples.index),

##### Helper functions #####

def find_files(base, pattern):
    '''Return list of files matching pattern in base folder.'''
    return [n for n in fnmatch.filter(os.listdir(base), pattern) if
        os.path.isfile(os.path.join(base, n))]

SAMPLES = samples.sample_name.to_list()

##### 7. Customize final output based on config["rule"] values #####
def get_final_output():
    """
    Generate final output for rule all given a TRUE value in config["rules"]
    """
    # dictionary of rules and its output files
    rule_dict = {"fileconversion" : expand("data/mzML/{samples}.mzML", samples=SAMPLES),
                "preprocessing" : [expand("results/Interim/mzML/PCpeak_{samples}.mzML", samples=SAMPLES),
        expand("results/Interim/Preprocessed/FFM_{samples}.featureXML", samples=SAMPLES),
        expand("results/Interim/mzML/PCfeature_{samples}.mzML", samples=SAMPLES),
        expand(["results/Interim/Preprocessed/MapAligned_{samples}.featureXML", "results/Interim/Preprocessed/MapAligned_{samples}.trafoXML"], samples=SAMPLES),
        expand("results/GNPSexport/mzML/Aligned_{samples}.mzML", samples=SAMPLES),
        expand("results/Interim/Preprocessed/MFD_{samples}.featureXML", samples=SAMPLES),
        expand("results/Interim/Preprocessed/Preprocessed.consensusXML"),
        expand("results/Preprocessed/FeatureMatrix.tsv")],
                "requantification" : [expand(["results/Interim/Requantified/Complete.consensusXML", "results/Interim/Requantified/Missing.consensusXML", "results/Interim/Requantified/Complete_{samples}.featureXML"], samples=SAMPLES),
        expand("results/Interim/Requantified/MetaboliteIdentification.tsv"),
        expand("results/Interim/Requantified/FFMID_{samples}.featureXML", samples=SAMPLES),
        expand("results/Interim/Requantified/Merged_{samples}.featureXML", samples=SAMPLES),
        expand("results/Interim/Requantified/MFD_{samples}.featureXML", samples=SAMPLES),
        expand("results/Interim/Requantified/IDMapper_{samples}.featureXML", samples=SAMPLES),
        expand("results/Interim/Requantified/Requantified.consensusXML"),
        expand("results/Requantified/FeatureMatrix.tsv")],
                "GNPSexport" : [expand("results/Interim/GNPSexport/filtered.consensusXML"),
        expand("results/GNPSexport/MSMS.mgf"),
        expand("results/GNPSexport/FeatureQuantificationTable.txt"),
        expand("results/GNPSexport/SuppPairs.csv"),
        expand("results/GNPSexport/metadata.tsv")],
                "sirius_csi" : [expand(["results/Interim/sirius/formulas_{samples}.mzTab", "results/Interim/sirius/structures_{samples}.mzTab"], samples=SAMPLES),
        expand(["results/SIRIUS/formulas_{samples}.tsv", "results/CSI/structures_{samples}.tsv"], samples=SAMPLES)],
                "sirius" : [expand("results/Interim/sirius/formulas_{samples}.mzTab", samples=SAMPLES),
        expand("results/SIRIUS/formulas_{samples}.tsv", samples=SAMPLES)],
                "annotate" : [expand("results/annotations/SIRIUS_CSI_annotated_FeatureTable.tsv"),
        expand("results/annotations/GNPS_annotated_FeatureTable.tsv"),
        ]
                }
    
    # get keys from config
    opt_rules = config["rules"].keys()

    # if values are TRUE add output files to rule all
    final_output = [rule_dict[r] for r in opt_rules if config["rules"][r]]

    return final_output