# 1) Correct the MS2 precursor on a peak level (To the "highest intensity MS1 peak")

rule precursorcorrection_peak:
    input:
        "data/mzML/{samples}.mzML"
    output:
        "results/Interim/mzML/PCpeak_{samples}.mzML"
    shell:
        """
        ./resources/OpenMS-2.7.0/bin/HighResPrecursorMassCorrector -in {input} -out {output} -highest_intensity_peak:mz_tolerance "100.0"
        """

# 2) Preprocessing: Feature finding algorithm that detects peaks 
    
rule preprocess:
    input:
        "results/Interim/mzML/PCpeak_{samples}.mzML"
    output:
        "results/Interim/preprocessed/FFM_{samples}.featureXML"
    shell:
        """
        ./resources/OpenMS-2.7.0/bin/FeatureFinderMetabo -in {input} -out {output} -algorithm:common:noise_threshold_int "1.0e04" -algorithm:mtd:mass_error_ppm "10.0" -algorithm:epd:width_filtering "fixed" -algorithm:ffm:isotope_filtering_model "none" -algorithm:ffm:remove_single_traces "true" -algorithm:ffm:report_convex_hulls "true"
        """

# 3) Decharger: Decharging algorithm for adduct assignment

rule decharge:
    input:
        "results/Interim/preprocessed/FFM_{samples}.featureXML"
    output:
        "results/Interim/preprocessed/MFD_{samples}.featureXML"
    shell:
        """
        ./resources/OpenMS-2.7.0/bin/MetaboliteAdductDecharger -in {input} -out_fm {output} -algorithm:MetaboliteFeatureDeconvolution:potential_adducts "H:+:0.4" "Na:+:0.2" "NH4:+:0.2" "H-1O-1:+:0.1" "H-3O-2:+:0.1" -algorithm:MetaboliteFeatureDeconvolution:charge_max "1" -algorithm:MetaboliteFeatureDeconvolution:charge_span_max "1"  -algorithm:MetaboliteFeatureDeconvolution:max_neutrals "1"
        """

# 4) Tables of features (individual)

rule individual_feature_tables:
    input:
        "results/Interim/preprocessed/MFD_{samples}.featureXML"
    output:
        "results/Interim/feature_tables/features_{samples}.csv"
    shell:
        """
        ./resources/OpenMS-2.7.0/bin/TextExporter -in {input} -out {output}
        """


# 4) Correct the MS2 precursor in a feature level (for GNPS FBMN).        

rule precursorcorrection_feature:
    input:
        var1= "results/Interim/mzML/PCpeak_{samples}.mzML",
        var2= "results/Interim/preprocessed/MFD_{samples}.featureXML"
    output:
        "results/Interim/mzML/PCfeature_{samples}.mzML"
    shell:
        """
        ./resources/OpenMS-2.7.0/bin/HighResPrecursorMassCorrector -in {input.var1} -feature:in {input.var2} -out {output}  -nearest_peak:mz_tolerance "100.0"
        """ 

# 5) MapAlignerPoseClustering is used to perform a linear retention time alignment, to correct for linear shifts in retention time between different runs.

rule MapAlignerPoseClustering:
    input:
        expand("results/Interim/preprocessed/MFD_{samples}.featureXML", samples=SAMPLES)
    output:
        var1= expand("results/Interim/preprocessed/MapAlignerPoseClustering_{samples}.featureXML", samples=SAMPLES),
        var2= expand("results/Interim/preprocessed/MapAlignerPoseClustering_{samples}.trafoXML", samples=SAMPLES)
    shell:
        """
        ./resources/OpenMS-2.7.0/bin/MapAlignerPoseClustering -algorithm:max_num_peaks_considered -1 -algorithm:superimposer:mz_pair_max_distance 0.05 -algorithm:pairfinder:distance_MZ:max_difference 10.0 -algorithm:pairfinder:distance_MZ:unit ppm -in {input} -out {output.var1} -trafo_out {output.var2}
        """ 

# 6) Introduce the features to a protein identification file (idXML)- the only way to annotate MS2 spectra for GNPS FBMN  

rule IDMapper:
    input:
        var1= "resources/emptyfile.idXML",
        var2= "results/Interim/preprocessed/MapAlignerPoseClustering_{samples}.featureXML",
        var3= "results/Interim/mzML/PCfeature_{samples}.mzML"
    output:
        "results/Interim/preprocessed/IDMapper_{samples}.featureXML"
    shell:
        """
        ./resources/OpenMS-2.7.0/bin/IDMapper -id {input.var1} -in {input.var2}  -spectra:in {input.var3} -out {output} 
        """

# 7) The FeatureLinkerUnlabeledKD is used to aggregate the feature information (from single files) into a ConsensusFeature, linking features from different files together, which have a similar m/z and rt (MS1 level).

rule FeatureLinkerUnlabeledKD:
    input:
        expand("results/Interim/preprocessed/IDMapper_{samples}.featureXML", samples=SAMPLES)
    output:
        "results/Interim/preprocessed/FeatureLinkerUnlabeledKD.consensusXML"
    shell:
        """
        ./resources/OpenMS-2.7.0/bin/FeatureLinkerUnlabeledKD -in {input} -out {output} 
        """

# 8) export the consensusXML file to a txt file

rule txt_export:
    input:
        "results/Interim/preprocessed/FeatureLinkerUnlabeledKD.consensusXML"
    output:
        "results/Interim/preprocessed/FeatureQuantificationTable.txt" 
    shell:
        """
        ./resources/OpenMS-2.7.0/bin/TextExporter -in {input} -out {output}
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