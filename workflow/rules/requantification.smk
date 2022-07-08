# Re-quantify the features in all data (missing values only)

# 1) Split the consensus map to features with no missing values (complete) and features with missing values (missing) and re-load the complete consensus to individual feature maps

rule split_consensus:
    input:
        "results/Interim/preprocessed/preprocessed.consensusXML",
    output:
        "results/Interim/Requantified/Complete.consensusXML",
        "results/Interim/Requantified/Missing.consensusXML"
    threads: 4
    conda:
        "../envs/pyopenms.yaml"
    script:
        "../scripts/split.py"

rule reload_maps:
    input:
        "results/Interim/preprocessed/MapAligned_{samples}.featureXML",
        "results/Interim/Requantified/Complete.consensusXML"
    output:
        "results/Interim/Requantified/Complete_{samples}.featureXML"
    conda:
        "../envs/pyopenms.yaml"
    script:
        "../scripts/reloadmaps.py"

# 2) Build a library of features from the consensus with missing values

rule text_export:
    input:
        "results/Interim/Requantified/Missing.consensusXML"
    output:
        "results/Interim/Requantified/FeatureQuantificationTable.txt" 
    shell:
        """
        /Users/eeko/openms-develop/openms_build/bin/TextExporter -in {input} -out {output}
        """

rule build_library:
    input:
        "results/Interim/Requantified/FeatureQuantificationTable.txt" 
    output:
        "results/Interim/Requantified/MetaboliteIdentification.tsv"
    threads: 4
    conda:
        "../envs/pyopenms.yaml"
    script:
        "../scripts/metaboliteidentification.py"    

# 3) Re-quantify all the raw files to cover missing values (missing value imputation can be avoided with that step)

rule requantify:
    input:
        var1= "results/Interim/Requantified/MetaboliteIdentification.tsv",
        var2= "results/GNPSexport/mzML/Aligned_{samples}.mzML"
    output:
        "results/Interim/Requantified/FFMID_{samples}.featureXML"
    threads: 4
    shell:
        """
        /Users/eeko/openms-develop/openms_build/bin/FeatureFinderMetaboIdent -id {input.var1} -in {input.var2} -out {output} -extract:mz_window 10.0 -threads {threads}
        """

# 4) Merge the re-quantified with the complete feature files

rule merge:
    input:
        "results/Interim/Requantified/Complete_{samples}.featureXML",
        "results/Interim/Requantified/FFMID_{samples}.featureXML"
    output:
        "results/Interim/Requantified/Merged_{samples}.featureXML"
    threads: 4
    conda:
        "../envs/pyopenms.yaml"
    script:
        "../scripts/merge.py"    


# 5) Decharger: Decharging algorithm for adduct assignment

rule adduct_annotations_FFMident:
    input:
        "results/Interim/Requantified/Merged_{samples}.featureXML"
    output:
        "results/Interim/Requantified/MFD_{samples}.featureXML"
    shell:
        """
        /Users/eeko/openms-develop/openms_build/bin/MetaboliteAdductDecharger -in {input} -out_fm {output} -algorithm:MetaboliteFeatureDeconvolution:potential_adducts "H:+:0.6" "Na:+:0.1" "NH4:+:0.1" "H-1O-1:+:0.1" "H-3O-2:+:0.1" -algorithm:MetaboliteFeatureDeconvolution:charge_max "1" -algorithm:MetaboliteFeatureDeconvolution:charge_span_max "1"  -algorithm:MetaboliteFeatureDeconvolution:max_neutrals "1"
        """
# 6) Introduce the features to a protein identification file (idXML)- the only way to annotate MS2 spectra for GNPS FBMN  

rule IDMapper:
    input:
        var1= "resources/emptyfile.idXML",
        var2= "results/Interim/Requantified/MFD_{samples}.featureXML",
        var3= "results/Interim/mzML/PCfeature_{samples}.mzML"
    output:
        "results/Interim/Requantified/IDMapper_{samples}.featureXML"
    shell:
        """
        /Users/eeko/openms-develop/openms_build/bin/IDMapper -id {input.var1} -in {input.var2}  -spectra:in {input.var3} -out {output} 
        """

# 7) The FeatureLinkerUnlabeledKD is used to aggregate the feature information (from single files) into a ConsensusFeature, linking features from different sfiles together, which have a smiliar m/z and rt (MS1 level).

rule FeatureLinker:
    input:
        expand("results/Interim/Requantified/IDMapper_{samples}.featureXML", samples=SAMPLES)
    output:
        "results/Interim/Requantified/Requantified.consensusXML"
    threads: 4
    shell:
        """
        /Users/eeko/openms-develop/openms_build/bin/FeatureLinkerUnlabeledKD -in {input} -out {output} -threads {threads}
        """

# 8) export the consensusXML file to a tsv file to produce a single matrix for PCA

rule FFMident_matrix:
    input:
        "results/Interim/Requantified/Requantified.consensusXML"
    output:
        "results/Requantified/FeatureMatrix.tsv"
    conda:
        "../envs/pyopenms.yaml"
    script:
        "../scripts/cleanup.py"


