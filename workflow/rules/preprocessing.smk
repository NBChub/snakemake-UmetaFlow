# 1) Correct the MS2 precursor on a peak level (To the "highest intensity MS1 peak")

rule precursorcorrection_peak:
    input:
        "data/mzML/{samples}.mzML"
    output:
        "results/Interim/mzML/PCpeak_{samples}.mzML"
    shell:
        """
        /Users/eeko/openms-develop/openms_build/bin/HighResPrecursorMassCorrector -in {input} -out {output} -highest_intensity_peak:mz_tolerance "100.0"
        """

# 2) Preprocessing: Feature finding algorithm that detects peaks 
    
rule preprocess:
    input:
        "results/Interim/mzML/PCpeak_{samples}.mzML"
    output:
        "results/Interim/preprocessed/nofilter_{samples}.featureXML"
    shell:
        """
        /Users/eeko/openms-develop/openms_build/bin/FeatureFinderMetabo -in {input} -out {output} -algorithm:common:noise_threshold_int "1.0e04" -algorithm:mtd:mass_error_ppm "10.0" -algorithm:epd:width_filtering "fixed" -algorithm:ffm:isotope_filtering_model "none" -algorithm:ffm:remove_single_traces "true" -algorithm:ffm:report_convex_hulls "true"
        """

# 3) Quality filter: Remove features that are of quality lower than 0.0005: 
    
rule quality:
    input:
        "results/Interim/preprocessed/nofilter_{samples}.featureXML"
    output:
        "results/Interim/preprocessed/FFM_{samples}.featureXML"
    shell:
        """
        /Users/eeko/openms-develop/openms_build/bin/FileFilter -in {input} -out {output} -feature:q 0.0005:
        """

# 4) Tables of features (individual)

rule individual_feature_tables:
    input:
        "results/Interim/preprocessed/MFD_{samples}.featureXML"
    output:
        "results/Interim/feature_tables/features_{samples}.csv"
    shell:
        """
        /Users/eeko/openms-develop/openms_build/bin/TextExporter -in {input} -out {output}
        """        

# 5) Correct the MS2 precursor in a feature level (for GNPS FBMN and SIRIUS).        

rule precursorcorrection_feature:
    input:
        var1= "results/Interim/mzML/PCpeak_{samples}.mzML",
        var2= "results/Interim/preprocessed/MFD_{samples}.featureXML"
    output:
        "results/Interim/mzML/PCfeature_{samples}.mzML"
    shell:
        """
        /Users/eeko/openms-develop/openms_build/bin/HighResPrecursorMassCorrector -in {input.var1} -feature:in {input.var2} -out {output}  -nearest_peak:mz_tolerance "100.0"
        """ 

# 6) MapAlignerPoseClustering is used to perform a linear retention time alignment, to correct for linear shifts in retention time between different runs.

rule MapAlignerPoseClustering:
    input:
        expand("results/Interim/preprocessed/MFD_{samples}.featureXML", samples=SAMPLES)
    output:
        var1= expand("results/Interim/preprocessed/MapAlignerPoseClustering_{samples}.featureXML", samples=SAMPLES),
        var2= expand("results/Interim/preprocessed/MapAlignerPoseClustering_{samples}.trafoXML", samples=SAMPLES)
    shell:
        """
        /Users/eeko/openms-develop/openms_build/bin/MapAlignerPoseClustering -algorithm:max_num_peaks_considered -1 -algorithm:superimposer:mz_pair_max_distance 0.05 -algorithm:pairfinder:distance_MZ:max_difference 10.0 -algorithm:pairfinder:distance_MZ:unit ppm -in {input} -out {output.var1} -trafo_out {output.var2}
        """ 

# 7) The FeatureLinkerUnlabeledKD is used to aggregate the feature information (from single files) into a ConsensusFeature, linking features from different files together, which have a similar m/z and rt (MS1 level).

rule FeatureLinkerUnlabeledKD:
    input:
        expand("results/Interim/preprocessed/IDMapper_{samples}.featureXML", samples=SAMPLES)
    output:
        "results/Interim/preprocessed/FeatureLinkerUnlabeledKD.consensusXML"
    shell:
        """
        /Users/eeko/openms-develop/openms_build/bin/FeatureLinkerUnlabeledKD -in {input} -out {output} 
        """

# 8) export the consensusXML file to a txt file

rule txt_export:
    input:
        "results/Interim/preprocessed/FeatureLinkerUnlabeledKD.consensusXML"
    output:
        "results/Interim/preprocessed/FeatureQuantificationTable.txt" 
    shell:
        """
        /Users/eeko/openms-develop/openms_build/bin/TextExporter -in {input} -out {output}
        """

# 9) Convert the table to an easily readable format

rule feature_table:
    input:
        "results/Interim/preprocessed/FeatureQuantificationTable.txt" 
    output:
        "results/Preprocessed/FeatureQuantificationTable.csv"
    threads: 4
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/cleanup.py"    

# Re-quantify the features in all data (missing value correction)
# This rule is currently not used - in progress


# 10) Build a library of features detected in all files from the consensus feature (using pandas)

rule build_library:
    input:
        "results/Preprocessed/FeatureQuantificationTable.csv"
    output:
        "results/Interim/Requantified/MetaboliteIdentification.tsv"
    conda:
        "../envs/python.yaml"   
    script:
        "../scripts/metaboliteidentification.py"

# 11) MapAlignerPoseClustering is used to perform a linear retention time alignment, to correct for linear shifts in retention time between different runs.

rule aligner:
    input:
        var1= "results/Interim/mzML/PCpeak_{samples}.mzML",
        var2= "results/Interim/preprocessed/MapAlignerPoseClustering_{samples}.trafoXML"
    output:
        "results/Interim/Requantified/Aligned_{samples}.mzML"
    shell:
        """
        /Users/eeko/openms-develop/openms_build/bin/MapRTTransformer -in {input.var1} -trafo_in {input.var2} -out {output}
        """ 

# 12) Re-quantify all the raw files to cover missing values (missing value imputation can be avoided with that step)

rule metaboident:
    input:
        var1= "results/Interim/Requantified/MetaboliteIdentification.tsv",
        var2= "results/Interim/Requantified/Aligned_{samples}.mzML"
    output:
        "results/Interim/Requantified/FFMID_{samples}.featureXML"
    shell:
        """
        /Users/eeko/openms-develop/openms_build/bin/FeatureFinderMetaboIdent -id {input.var1} -in {input.var2} -out {output} -extract:mz_window 5.0 -detect:peak_width 20.0
        """

# 13) Decharger: Decharging algorithm for adduct assignment

rule decharge:
    input:
        "results/Interim/Requantified/FFMID_{samples}.featureXML"
    output:
        "results/Interim/Requantified/MFD_{samples}.featureXML"
    shell:
        """
        /Users/eeko/openms-develop/openms_build/bin/MetaboliteAdductDecharger -in {input} -out_fm {output} -algorithm:MetaboliteFeatureDeconvolution:potential_adducts "H:+:0.6" "Na:+:0.1" "NH4:+:0.1" "H-1O-1:+:0.1" "H-3O-2:+:0.1" -algorithm:MetaboliteFeatureDeconvolution:charge_max "1" -algorithm:MetaboliteFeatureDeconvolution:charge_span_max "1"  -algorithm:MetaboliteFeatureDeconvolution:max_neutrals "1"
        """

# 14) The FeatureLinkerUnlabeledKD is used to aggregate the feature information (from single files) into a ConsensusFeature, linking features from different sfiles together, which have a smiliar m/z and rt (MS1 level).

rule FeatureLinker:
    input:
        expand("results/Interim/Requantified/MFD_{samples}.featureXML", samples=SAMPLES)
    output:
        "results/Interim/Requantified/Requantified.consensusXML"
    shell:
        """
        /Users/eeko/openms-develop/openms_build/bin/FeatureLinkerUnlabeledKD -in {input} -out {output}
        """

# 15) export the consensusXML file to a csv file to produce a single matrix for PCA

rule matrix:
    input:
        "results/Interim/Requantified/Requantified.consensusXML"
    output:
        "results/Interim/Requantified/consensus.tsv" 
    shell:
        """
        /Users/eeko/openms-develop/openms_build/bin/TextExporter -in {input} -out {output}
        """
        
# 16) Convert the table to an easily readable format:

rule cleanup:
    input:
        "results/Interim/Requantified/consensus.tsv" 
    output:
        "results/Requantified/FeatureMatrix_requant.tsv"
    conda:
        "../envs/python.yaml"   
    script:
        "../scripts/cleanup.py"