# Re-quantify the features in all data (missing value correction)
# This rule is currently not used - in progress

# 1) rename tsv to csv for pandas processing

rule rename:
    input:
        "results/Interim/GNPSexport/consensus.tsv"
    output:
        "results/Interim/Requantified/consensus.csv"
    shell:
        """
        mv {input} {output}
        """

# 2) Build a library of features detected in all files from the consensus feature (using pandas)

rule build_library:
    input:
        "results/Interim/Requantified/consensus.csv"
    output:
        "results/Interim/Requantified/MetaboliteIdentification.tsv"
    conda:
        "../envs/python.yaml"   
    script:
        "../scripts/metaboliteidentification.py"

# 3) MapAlignerPoseClustering is used to perform a linear retention time alignment, to correct for linear shifts in retention time between different runs.

rule aligner:
    input:
        var1= "results/Interim/mzML/PCpeak_{samples}.mzML",
        var2= "results/Interim/preprocessed/MapAlignerPoseClustering_{samples}.trafoXML"
    output:
        "results/Interim/Requantified/Aligned_{samples}.mzML"
    shell:
        """
        resources/OpenMS-2.7.0/bin/MapRTTransformer -in {input.var1} -trafo_in {input.var2} -out {output}
        """ 

# 4) Re-quantify all the raw files to cover missing values (missing value imputation can be avoided with that step)

rule metaboident:
    input:
        var1= "results/Interim/Requantified/MetaboliteIdentification.tsv",
        var2= "results/Interim/Requantified/Aligned_{samples}.mzML"
    output:
        "results/Interim/Requantified/FFMID_{samples}.featureXML"
    shell:
        """
        resources/OpenMS-2.7.0/bin/FeatureFinderMetaboIdent -id {input.var1} -in {input.var2} -out {output} -extract:mz_window 5.0 -detect:peak_width 20.0 -no_progress
        """


# 5) The FeatureLinkerUnlabeledKD is used to aggregate the feature information (from single files) into a ConsensusFeature, linking features from different files together, which have a smiliar m/z and rt (MS1 level).

rule FeatureLinker:
    input:
        expand("results/Interim/Requantified/FFMID_{samples}.featureXML", samples=SAMPLES)
    output:
        "results/Interim/Requantified/Requantified.consensusXML"
    shell:
        """
        resources/OpenMS-2.7.0/bin/FeatureLinkerUnlabeledKD -in {input} -out {output} -algorithm:warp:enabled false
        """

# 6) export the consensusXML file to a csv file to produce a single matrix for PCA

rule matrix:
    input:
        "results/Interim/Requantified/Requantified.consensusXML"
    output:
        "results/Interim/Requantified/consensus.tsv" 
    shell:
        """
        resources/OpenMS-2.7.0/bin/TextExporter -in {input} -out {output}
        """
        
# 7) Clean up consensus.tsv to a proper FeatureMatrix table:

rule cleanup:
    input:
        "results/Interim/Requantified/consensus.tsv" 
    output:
        "results/Requantified/FeatureMatrix.tsv"
    conda:
        "../envs/python.yaml"   
    script:
        "../scripts/cleanup.py"