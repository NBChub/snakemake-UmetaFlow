# Re-quantify the features in all data (missing values only)

# 1) Split the consensus map to features with no missing values (complete) and features with missing values (missing) and re-load the complete consensus to individual feature maps
<<<<<<< HEAD
rule split_consensus:
    input:
        "results/Interim/preprocessed/preprocessed.consensusXML",
    output:
        "results/Interim/Requantified/Complete.consensusXML",
        "results/Interim/Requantified/Missing.consensusXML"
    threads: 4
    conda:
        "../envs/exe.yaml"
    script:
        "../scripts/split.py"

rule reload_maps:
    input:
        "results/Interim/preprocessed/MapAligned_{samples}.featureXML",
        "results/Interim/Requantified/Complete.consensusXML"
    output:
        "results/Interim/Requantified/Complete_{samples}.featureXML"
    conda:
        "../envs/exe.yaml"
    script:
        "../scripts/reloadmaps.py"

=======

rule split_consensus:
    input:
        "results/Interim/preprocessed/preprocessed.consensusXML",
    output:
        "results/Interim/Requantified/Complete.consensusXML",
        "results/Interim/Requantified/Missing.consensusXML",
    threads: 4
    conda:
        "../envs/python.yaml"   
    script:
        "../scripts/split.py"

rule reload_maps:
    input:
        expand("results/Interim/preprocessed/MapAligned_{samples}.featureXML", samples=SAMPLES)
    output:
        "results/Interim/Requantified/Complete_{samples}.featureXML"
    conda:
        "../envs/python.yaml"   
    script:
        "../scripts/reloadmaps.py"

>>>>>>> 36aac52cfd4423bb92386ddfb3059a9655b26263
# 2) Build a library of features from the consensus with missing values

rule txt_export_2:
    input:
        "results/Interim/Requantified/Missing.consensusXML"
    output:
        "results/Interim/Requantified/FeatureQuantificationTable.txt" 
    shell:
        """
<<<<<<< HEAD
        /Users/eeko/openms-develop/openms_build/bin/TextExporter -in {input} -out {output}
=======
        TextExporter -in {input} -out {output}
>>>>>>> 36aac52cfd4423bb92386ddfb3059a9655b26263
        """

rule build_library:
    input:
        "results/Interim/Requantified/FeatureQuantificationTable.txt" 
    output:
        "results/Interim/Requantified/MetaboliteIdentification.tsv"
    threads: 4
    conda:
<<<<<<< HEAD
        "../envs/exe.yaml"
=======
        "../envs/python.yaml"
>>>>>>> 36aac52cfd4423bb92386ddfb3059a9655b26263
    script:
        "../scripts/metaboliteidentification.py"    

# 3) MapRTTransformer is used to perform a linear retention time alignment, to correct for linear shifts in retention time between different runs using the transformation files from the reprocessing rule MapAlignerPoseClustering (faster computationally)

rule aligner:
    input:
        var1= "results/Interim/mzML/PCfeature_{samples}.mzML",
        var2= "results/Interim/preprocessed/MapAligned_{samples}.trafoXML"
    output:
        "results/Interim/Requantified/Aligned_{samples}.mzML"
    shell:
        """
        /Users/eeko/openms-develop/openms_build/bin/MapRTTransformer -in {input.var1} -trafo_in {input.var2} -out {output}
        """ 

# 4) Re-quantify all the raw files to cover missing values (missing value imputation can be avoided with that step)

rule requant:
    input:
        var1= "results/Interim/Requantified/MetaboliteIdentification.tsv",
        var2= "results/Interim/Requantified/Aligned_{samples}.mzML"
    output:
        "results/Interim/Requantified/FFMID_{samples}.featureXML"
    threads: 4
<<<<<<< HEAD
    shell:
        """
        /Users/eeko/openms-develop/openms_build/bin/FeatureFinderMetaboIdent -id {input.var1} -in {input.var2} -out {output} -extract:mz_window 10.0 -extract:rt_window 30.0 -detect:peak_width 60.0 -threads {threads}
        """

# 5) Merge the re-quantified with the complete feature files

rule merge:
    input:
        "results/Interim/Requantified/Complete_{samples}.featureXML",
        "results/Interim/Requantified/FFMID_{samples}.featureXML"
    output:
        "results/Interim/Requantified/Merged_{samples}.featureXML"
    threads: 4
    conda:
        "../envs/exe.yaml"
    script:
        "../scripts/merge.py"    

# 6) Introduce the features to a protein identification file (idXML)- the only way to annotate MS2 spectra for GNPS FBMN  

rule IDMapper:
    input:
        var1= "resources/emptyfile.idXML",
        var2= "results/Interim/Requantified/Merged_{samples}.featureXML",
        var3= "results/Interim/mzML/PCfeature_{samples}.mzML"
    output:
        "results/Interim/Requantified/IDMapper_{samples}.featureXML"
    shell:
        """
        /Users/eeko/openms-develop/openms_build/bin/IDMapper -id {input.var1} -in {input.var2}  -spectra:in {input.var3} -out {output} 
        """

# 7) Decharger: Decharging algorithm for adduct assignment

rule decharge:
    input:
        "results/Interim/Requantified/IDMapper_{samples}.featureXML"
    output:
        "results/Interim/Requantified/MFD_{samples}.featureXML"
    shell:
        """
        /Users/eeko/openms-develop/openms_build/bin/MetaboliteAdductDecharger -in {input} -out_fm {output} -algorithm:MetaboliteFeatureDeconvolution:potential_adducts "H:+:0.6" "Na:+:0.1" "NH4:+:0.1" "H-1O-1:+:0.1" "H-3O-2:+:0.1" -algorithm:MetaboliteFeatureDeconvolution:charge_max "1" -algorithm:MetaboliteFeatureDeconvolution:charge_span_max "1"  -algorithm:MetaboliteFeatureDeconvolution:max_neutrals "1"
        """

=======
    shell:
        """
        FeatureFinderMetaboIdent -id {input.var1} -in {input.var2} -out {output} -extract:mz_window 10.0 -extract:rt_window 30.0 -detect:peak_width 60.0 -threads {threads}
        """

# 5) Merge the re-quantified with the complete feature files

rule merge:
    input:
        "results/Interim/Requantified/Complete_{samples}.featureXML",
        "results/Interim/Requantified/FFMID_{samples}.featureXML"
    output:
        "results/Interim/Requantified/Merged_{samples}.featureXML"
    threads: 4
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/merge.py"    

# 6) Introduce the features to a protein identification file (idXML)- the only way to annotate MS2 spectra for GNPS FBMN  

rule IDMapper:
    input:
        var1= "resources/emptyfile.idXML",
        var2= "results/Interim/Requantified/Merged_{samples}.featureXML",
        var3= "results/Interim/mzML/PCfeature_{samples}.mzML"
    output:
        "results/Interim/Requantified/IDMapper_{samples}.featureXML"
    shell:
        """
        IDMapper -id {input.var1} -in {input.var2}  -spectra:in {input.var3} -out {output} 
        """

# 7) Decharger: Decharging algorithm for adduct assignment

rule decharge:
    input:
        "results/Interim/Requantified/IDMapper_{samples}.featureXML"
    output:
        "results/Interim/Requantified/MFD_{samples}.featureXML"
    shell:
        """
        MetaboliteAdductDecharger -in {input} -out_fm {output} -algorithm:MetaboliteFeatureDeconvolution:potential_adducts "H:+:0.6" "Na:+:0.1" "NH4:+:0.1" "H-1O-1:+:0.1" "H-3O-2:+:0.1" -algorithm:MetaboliteFeatureDeconvolution:charge_max "1" -algorithm:MetaboliteFeatureDeconvolution:charge_span_max "1"  -algorithm:MetaboliteFeatureDeconvolution:max_neutrals "1"
        """

>>>>>>> 36aac52cfd4423bb92386ddfb3059a9655b26263
# 8) The FeatureLinkerUnlabeledKD is used to aggregate the feature information (from single files) into a ConsensusFeature, linking features from different sfiles together, which have a smiliar m/z and rt (MS1 level).

rule FeatureLinker:
    input:
        expand("results/Interim/Requantified/MFD_{samples}.featureXML", samples=SAMPLES)
    output:
        "results/Interim/Requantified/Requantified.consensusXML"
    threads: 4
    shell:
        """
<<<<<<< HEAD
        /Users/eeko/openms-develop/openms_build/bin/FeatureLinkerUnlabeledKD -in {input} -out {output} -threads {threads}
=======
        FeatureLinkerUnlabeledKD -in {input} -out {output} -threads {threads}
>>>>>>> 36aac52cfd4423bb92386ddfb3059a9655b26263
        """

# 9) export the consensusXML file to a tsv file to produce a single matrix for PCA

rule matrix:
    input:
        "results/Interim/Requantified/Requantified.consensusXML"
    output:
        "results/Interim/Requantified/consensus.tsv" 
    shell:
        """
        /Users/eeko/openms-develop/openms_build/bin/TextExporter -in {input} -out {output}
        """
        
# 10) Convert the table to an easily readable format:

rule cleanup:
    input:
        "results/Interim/Requantified/consensus.tsv" 
    output:
        "results/Requantified/FeatureMatrix.tsv"
    conda:
        "../envs/exe.yaml"
    script:
        "../scripts/cleanup.py"


