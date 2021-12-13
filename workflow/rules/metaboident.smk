# Re-quantify the features in all data (missing value correction)
# This rule is currently not used - in progress

# 1) rename tsv to csv for pandas processing

rule rename:
    input:
        "results/GNPSexport/interim/consensus.tsv"
    output:
        "results/GNPSexport/interim/consensus.csv"
    shell:
        """
        mv {input} {output}
        """

# 2) Build a library of features detected in all files from the consensus feature (using pandas)

rule build_library:
    input:
        "results/GNPSexport/interim/consensus.csv"
    output:
        "resources/MetaboliteIdentification.tsv"
    conda:
        "../envs/file_conversion.yaml"   
    script:
        "../scripts/metaboliteidentification.py"

# 3) MapAlignerPoseClustering is used to perform a linear retention time alignment, to correct for linear shifts in retention time between different runs.

rule aligner:
    input:
        expand("results/{samples}/interim/PCpeak_{samples}.mzML", samples=SAMPLES)
    output:
        expand("results/Requant/interim/MapAlignerPoseClustering_{samples}.mzML", samples=SAMPLES)
    shell:
        """
        /nfs/wsi/abi/scratch/alka/openms/openms_build/bin/MapAlignerPoseClustering -in {input} -out {output} -algorithm:max_num_peaks_considered 3000 
        """ 

# 4) Re-quantify all the raw files to cover missing values (missing value imputation can be avoided with that step)

rule metaboident:
    input:
        var1= "resources/MetaboliteIdentification.tsv",
        var2= "results/Requant/interim/MapAlignerPoseClustering_{samples}.mzML"
    output:
        "results/Requant/interim/FFMID_{samples}.featureXML"
    shell:
        """
        /nfs/wsi/abi/scratch/alka/openms/openms_build/bin/FeatureFinderMetaboIdent -id {input.var1} -in {input.var2} -out {output} -extract:mz_window 5.0 -detect:peak_width 20.0
        """

# 5) Export the consensusXML file to a csv file 

rule FFMI_df:
    input:
        "results/Requant/interim/FFMID_{samples}.featureXML"
    output:
        "results/Requant/FFMID_{samples}.csv" 
    shell:
        """
        /nfs/wsi/abi/scratch/alka/openms/openms_build/bin/TextExporter -in {input} -out {output}
        """

# 6) The FeatureLinkerUnlabeledKD is used to aggregate the feature information (from single files) into a ConsensusFeature, linking features from different files together, which have a smiliar m/z and rt (MS1 level).

rule FeatureLinker:
    input:
        expand("results/Requant/interim/FFMID_{samples}.featureXML", samples=SAMPLES)
    output:
        "results/Requant/interim/Requant.consensusXML"
    shell:
        """
        /nfs/wsi/abi/scratch/alka/openms/openms_build/bin/FeatureLinkerUnlabeledKD -in {input} -out {output} 
        """

# 7) export the consensusXML file to a csv file to produce a single matrix for PCA

rule matrix:
    input:
        "results/Requant/interim/Requant.consensusXML"
    output:
        "results/Requant/consensus.tsv" 
    shell:
        """
        /nfs/wsi/abi/scratch/alka/openms/openms_build/bin/TextExporter -in {input} -out {output}
        """