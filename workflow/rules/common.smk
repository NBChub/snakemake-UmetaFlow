import os
import pandas as pd
from snakemake.utils import validate
from snakemake.utils import min_version

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


SAMPLES = samples.sample_name.to_list()

##### 7. Customize final output based on config["rule"] values #####
def get_final_output():
    """
    Generate final output for rule all given a TRUE value in config["rules"]
    """
    # dictionary of rules and its output files
    rule_dict = {"fileconversion" : expand("data/mzML/{samples}.mzML", samples=SAMPLES),
                "preprocessing" : [expand("results/Interim/mzML/PCpeak_{samples}.mzML", samples=SAMPLES),
        expand("results/Interim/preprocessed/FFM_{samples}.featureXML", samples=SAMPLES),
        expand("results/Interim/preprocessed/MFD_{samples}.featureXML", samples=SAMPLES),
        expand("results/Interim/feature_tables/features_{samples}.csv", samples=SAMPLES),
        expand("results/Interim/mzML/PCfeature_{samples}.mzML", samples=SAMPLES),
        expand(["results/Interim/preprocessed/MapAlignerPoseClustering_{samples}.featureXML", "results/Interim/preprocessed/MapAlignerPoseClustering_{samples}.trafoXML"], samples=SAMPLES),
        expand("results/Interim/preprocessed/IDMapper_{samples}.featureXML", samples=SAMPLES),
        expand("results/Interim/preprocessed/FeatureLinkerUnlabeledKD.consensusXML"),
        expand("results/Interim/preprocessed/FeatureQuantificationTable.txt"),
        expand("results/Preprocessed/FeatureQuantificationTable.csv")],
                "GNPSexport" : [expand("results/Interim/GNPSexport/filtered.consensusXML"),
        expand("results/GNPSexport/MSMS.mgf"),
        expand("results/GNPSexport/FeatureQuantificationTable.txt")],
                "requantification" : [expand("results/Interim/Requantified/MetaboliteIdentification.tsv"),
        expand("results/Interim/Requantified/Aligned_{samples}.mzML", samples=SAMPLES),
        expand("results/Interim/Requantified/FFMID_{samples}.featureXML", samples=SAMPLES),
        expand("results/Interim/Requantified/Requantified.consensusXML"),
        expand("results/Interim/Requantified/consensus.tsv"),
        expand("results/Requantified/FeatureMatrix_requant.tsv")],
                "sirius" : [expand("results/Interim/sirius/FFM_nch_{samples}.featureXML", samples=SAMPLES),
        expand("results/Interim/sirius/MFD_nch_{samples}.featureXML", samples=SAMPLES),
        expand("results/Interim/sirius/PCfeature_nch_{samples}.mzML", samples=SAMPLES),
        expand(["results/Interim/sirius/formulas_{samples}.mzTab", "results/Interim/sirius/structures_{samples}.mzTab"], samples=SAMPLES),
        expand(["results/SIRIUS/formulas_{samples}.csv", "results/CSI/structures_{samples}.csv"], samples=SAMPLES)],
                "annotate" : [expand("results/SIRIUS/SIRIUS_library.csv"),
        expand("results/CSI/CSI_library.csv"),
        expand("results/annotations/SIRIUS_CSI_annotated_FeatureTable.csv"),
        expand("results/annotations/SIRIUS_CSI_annotated_FeatureTable_Requant.csv"),
        expand("results/annotations/GNPS_annotated_FeatureTable.csv"),
        expand("results/annotations/GNPS_annotated_FeatureTable_Requant.csv")]
                }
    
    # get keys from config
    opt_rules = config["rules"].keys()

    # if values are TRUE add output files to rule all
    final_output = [rule_dict[r] for r in opt_rules if config["rules"][r]]

    return final_output