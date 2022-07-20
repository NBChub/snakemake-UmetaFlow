# Re-quantify the features in all data (missing values only)

# 1) Split the consensus map to features with no missing values (complete) and features with missing values (missing) and re-load the complete consensus to individual feature maps

rule split_consensus:
    input:
        in_cmap= "results/Interim/Preprocessed/Preprocessed.consensusXML",
    output:
        out_complete= "results/Interim/Requantified/Complete.consensusXML",
        out_missing= "results/Interim/Requantified/Missing.consensusXML"
    log: "workflow/report/logs/requantification/split_consensus.log"
    threads: 4
    conda:
        "../envs/pyopenms.yaml"
    shell:    
        """
        python workflow/scripts/split.py {input.in_cmap} {output.out_complete} {output.out_missing} 2>> {log}
        """

rule reload_maps:
    input:
        in_aligned= "results/Interim/Preprocessed/MapAligned_{samples}.featureXML",
        in_complete= "results/Interim/Requantified/Complete.consensusXML"
    output:
        out_complete= "results/Interim/Requantified/Complete_{samples}.featureXML"
    log: "workflow/report/logs/requantification/reload_maps_{samples}.log"
    threads: 4
    conda:
        "../envs/pyopenms.yaml"
    shell:    
        """
        python workflow/scripts/reloadmaps.py {input.in_aligned} {input.in_complete} {output.out_complete} 2>> {log}
        """

# 2) Build a library of features from the consensus with missing values

rule text_export:
    input:
        "results/Interim/Requantified/Missing.consensusXML"
    output:
        "results/Interim/Requantified/FeatureQuantificationTable.txt" 
    log: "workflow/report/logs/requantification/text_export.log"
    conda:
        "../envs/openms.yaml"
    shell:
        """
        TextExporter -in {input} -out {output} -log {log} 2>> {log}
        """

rule build_library:
    input:
        matrix= "results/Interim/Requantified/FeatureQuantificationTable.txt" 
    output:
        lib= "results/Interim/Requantified/MetaboliteIdentification.tsv"
    log: "workflow/report/logs/requantification/build_library.log"
    threads: 4
    conda:
        "../envs/pyopenms.yaml"
    shell:    
        """
        python workflow/scripts/metaboliteidentification.py {input.matrix} {output.lib} 2>> {log}   
        """

# 3) Re-quantify all the raw files to cover missing values (missing value imputation can be avoided with that step)

rule requantify:
    input:
        var1= "results/Interim/Requantified/MetaboliteIdentification.tsv",
        var2= "results/GNPSexport/mzML/Aligned_{samples}.mzML"
    output:
        "results/Interim/Requantified/FFMID_{samples}.featureXML"
    log: "workflow/report/logs/requantification/requantify_{samples}.log"
    conda:
        "../envs/openms.yaml"
    threads: 4
    shell:
        """
        FeatureFinderMetaboIdent -id {input.var1} -in {input.var2} -out {output} -extract:mz_window 10.0 -threads {threads} -log {log} 2>> {log}
        """

# 4) Merge the re-quantified with the complete feature files

rule merge:
    input:
        in_complete= "results/Interim/Requantified/Complete_{samples}.featureXML",
        in_requant= "results/Interim/Requantified/FFMID_{samples}.featureXML"
    output:
        out_merged= "results/Interim/Requantified/Merged_{samples}.featureXML"
    log: "workflow/report/logs/requantification/merge_{samples}.log"
    threads: 4
    conda:
        "../envs/pyopenms.yaml"
    shell:    
        """
        python workflow/scripts/merge.py {input.in_complete} {input.in_requant} {output.out_merged} 2>> {log}
        """


# 5) Decharger: Decharging algorithm for adduct assignment

rule adduct_annotations_FFMident:
    input:
        "results/Interim/Requantified/Merged_{samples}.featureXML"
    output:
        "results/Interim/Requantified/MFD_{samples}.featureXML"
    log: "workflow/report/logs/requantification/adduct_annotations_FFMident_{samples}.log"
    threads: 4
    conda:
        "../envs/openms.yaml"
    shell:
        """
        MetaboliteAdductDecharger -in {input} -out_fm {output} -algorithm:MetaboliteFeatureDeconvolution:potential_adducts "H:+:0.6" "Na:+:0.1" "NH4:+:0.1" "H-1O-1:+:0.1" "H-3O-2:+:0.1" -algorithm:MetaboliteFeatureDeconvolution:charge_max "1" -algorithm:MetaboliteFeatureDeconvolution:charge_span_max "1"  -algorithm:MetaboliteFeatureDeconvolution:max_neutrals "1" -threads {threads} -log {log} 2>> {log}
        """
# 6) Introduce the features to a protein identification file (idXML)- the only way to annotate MS2 spectra for GNPS FBMN  

rule IDMapper_FFMident:
    input:
        var1= "resources/emptyfile.idXML",
        var2= "results/Interim/Requantified/MFD_{samples}.featureXML",
        var3= "results/Interim/mzML/PCfeature_{samples}.mzML"
    output:
        "results/Interim/Requantified/IDMapper_{samples}.featureXML"
    log: "workflow/report/logs/requantification/IDMapper_FFMident_{samples}.log"
    threads: 4
    conda:
        "../envs/openms.yaml"
    shell:
        """
        IDMapper -id {input.var1} -in {input.var2} -spectra:in {input.var3} -out {output} -threads {threads} -log {log} 2>> {log}
        """

# 7) The FeatureLinkerUnlabeledKD is used to aggregate the feature information (from single files) into a ConsensusFeature, linking features from different sfiles together, which have a smiliar m/z and rt (MS1 level).

rule FeatureLinker_FFMident:
    input:
        expand("results/Interim/Requantified/IDMapper_{samples}.featureXML", samples=SAMPLES)
    output:
        "results/Interim/Requantified/Requantified.consensusXML"
    log: "workflow/report/logs/requantification/FeatureLinker_FFMident.log"
    conda:
        "../envs/openms.yaml"
    threads: 4
    shell:
        """
        FeatureLinkerUnlabeledKD -in {input} -out {output} -algorithm:warp:enabled false -algorithm:link:rt_tol 30.0 -algorithm:link:mz_tol 8.0 -threads {threads} -log {log} 2>> {log} 
        """

# 8) export the consensusXML file to a tsv file to produce a single matrix for PCA

rule FFMident_matrix:
    input:
        input_cmap= "results/Interim/Requantified/Requantified.consensusXML"
    output:
        output_tsv= "results/Requantified/FeatureMatrix.tsv"
    log: "workflow/report/logs/requantification/FFMident_matrix.log"
    conda:
        "../envs/pyopenms.yaml"
    shell:
        """
        python workflow/scripts/cleanup.py {input.input_cmap} {output.output_tsv} 2>> {log}
        """


