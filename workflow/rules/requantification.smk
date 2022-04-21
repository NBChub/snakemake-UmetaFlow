# Re-quantify the features in all data (missing value correction)
# This rule is currently not used - in progress
import glob
from os.path import join

# 1) Build a library of features detected in all files from the consensus feature (using pandas)

rule build_library:
    input:
        "results/Preprocessed/FeatureQuantificationTable.csv"
    output:
        "results/Interim/Requantified/MetaboliteIdentification.tsv"
    conda:
        "../envs/python.yaml"   
    script:
        "../scripts/metaboliteidentification.py"

# 2) MapAlignerPoseClustering is used to perform a linear retention time alignment, to correct for linear shifts in retention time between different runs.

rule aligner:
    input:
        var1= "results/Interim/mzML/PCpeak_{samples}.mzML",
        var2= "results/Interim/preprocessed/MapAlignerPoseClustering_{samples}.trafoXML"
    output:
        "results/Interim/Requantified/Aligned_{samples}.mzML"
    params:
        bin_path= glob.glob(join('.snakemake','conda','*','bin','MapRTTransformer'))[0]
    conda:
        "../envs/python.yaml"
    shell:
        """
        {params.bin_path} -in {input.var1} -trafo_in {input.var2} -out {output}
        """ 

# 3) Re-quantify all the raw files to cover missing values (missing value imputation can be avoided with that step)

rule metaboident:
    input:
        var1= "results/Interim/Requantified/MetaboliteIdentification.tsv",
        var2= "results/Interim/Requantified/Aligned_{samples}.mzML"
    output:
        "results/Interim/Requantified/FFMID_{samples}.featureXML"
    params:
        bin_path= glob.glob(join('.snakemake','conda','*','bin','FeatureFinderMetaboIdent'))[0]
    conda:
        "../envs/python.yaml"
    shell:
        """
        {params.bin_path} -id {input.var1} -in {input.var2} -out {output} -extract:mz_window 5.0 -detect:peak_width 40.0
        """


# 4) The FeatureLinkerUnlabeledKD is used to aggregate the feature information (from single files) into a ConsensusFeature, linking features from different sfiles together, which have a smiliar m/z and rt (MS1 level).

rule FeatureLinker:
    input:
        expand("results/Interim/Requantified/FFMID_{samples}.featureXML", samples=SAMPLES)
    output:
        "results/Interim/Requantified/Requantified.consensusXML"
    params:
        bin_path= glob.glob(join('.snakemake','conda','*','bin','FeatureLinkerUnlabeledKD'))[0]
    conda:
        "../envs/python.yaml"
    shell:
        """
        {params.bin_path} -in {input} -out {output}
        """

# 5) export the consensusXML file to a csv file to produce a single matrix for PCA

rule matrix:
    input:
        "results/Interim/Requantified/Requantified.consensusXML"
    output:
        "results/Interim/Requantified/consensus.tsv" 
    params:
        bin_path= glob.glob(join('.snakemake','conda','*','bin','TextExporter'))[0]
    conda:
        "../envs/python.yaml"
    shell:
        """
        {params.bin_path} -in {input} -out {output}
        """
        
# 6) Convert the table to an easily readable format:

rule cleanup:
    input:
        "results/Interim/Requantified/consensus.tsv" 
    output:
        "results/Requantified/FeatureMatrix_requant.tsv"
    conda:
        "../envs/python.yaml"   
    script:
        "../scripts/cleanup.py"