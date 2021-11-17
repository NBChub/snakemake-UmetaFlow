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
        "../envs/file_convertion.yaml"   
    script:
        "../scripts/metaboliteidentidication.py"

# 3) Re-quantify all the raw files to cover missing values (missing value imputation can be avoided with that step)

rule metaboident:
    input:
        "resources/MetaboliteIdentification.tsv",
        "results/{samples}/interim/{samples}.mzML"
    output:
        "results/Requant/interim/FFMID_{samples}.featureXML"
    shell:
        """
        resources/OpenMS-2.7.0/bin/FeatureFinderMetaboIdent -id {input[0]} -extract:mz_window 5.0 -in {input[1]} -out {output}
        """

# 4) Export the consensusXML file to a csv file 

rule FFMI_df:
    input:
        "results/Requant/interim/FFMID_{samples}.featureXML"
    output:
        "results/Requant/FFMID_{samples}.csv" 
    shell:
        """
        resources/OpenMS-2.7.0/bin/TextExporter -in {input} -out {output}
        """

# 5) Introduce the features to a protein identification file (idXML)- the only way to create an aggregated ConsensusXML file currently (run FeatureLinkerUnlabeledKD)  

rule IDMapper_FFMID:
    input:
        "resources/emptyfile.idXML",
        "results/Requant/interim/FFMID_{samples}.featureXML",
        "results/{samples}/interim/{samples}.mzML"
    output:
        "results/Requant/interim/IDMapper_FFMID{samples}.featureXML"
    shell:
        """
        resources/OpenMS-2.7.0/bin/IDMapper -id {input[0]} -in {input[1]}  -spectra:in {input[2]} -out {output} 
        """

# 6) The FeatureLinkerUnlabeledKD is used to aggregate the feature information (from single files) into a ConsensusFeature, linking features from different files together, which have a smiliar m/z and rt (MS1 level).

rule FeatureLinkerUnlabeledKD_requant:
    input:
        expand("results/Requant/interim/IDMapper_FFMID{samples}.featureXML", samples=SAMPLES)
    output:
        "results/Requant/interim/FeatureLinkerUnlabeledKD.consensusXML"
    shell:
        """
        resources/OpenMS-2.7.0/bin/FeatureLinkerUnlabeledKD -in {input} -out {output} 
        """

# 7) export the consensusXML file to a csv file to produce a single matrix for PCA

rule matrix:
    input:
        "results/Requant/interim/FeatureLinkerUnlabeledKD.consensusXML"
    output:
        "results/Requant/consensus.tsv" 
    shell:
        """
        resources/OpenMS-2.7.0/bin/TextExporter -in {input} -out {output}
        """