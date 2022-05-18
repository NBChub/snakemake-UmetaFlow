# 1) Correct the MS2 precursor on a peak level (To the "highest intensity MS1 peak")

rule precursorcorrection_peak:
    input:
        "data/mzML/{samples}.mzML"
    output:
        "results/Interim/mzML/PCpeak_{samples}.mzML"
    shell:
        """
        HighResPrecursorMassCorrector -in {input} -out {output} -highest_intensity_peak:mz_tolerance "100.0"
        """

# 2) Preprocessing: Feature finding algorithm that detects peaks 
    
rule preprocess:
    input:
        "results/Interim/mzML/PCpeak_{samples}.mzML"
    output:
        "results/Interim/preprocessed/FFM_{samples}.featureXML"
    threads: 4
    shell:
        """
        FeatureFinderMetabo -in {input} -out {output} -algorithm:common:noise_threshold_int "1.0e04" -algorithm:mtd:mass_error_ppm "10.0" -algorithm:epd:width_filtering "fixed" -algorithm:ffm:isotope_filtering_model "none" -algorithm:ffm:remove_single_traces "true" -algorithm:ffm:report_convex_hulls "true" -threads {threads}
        """

# 3) Correct the MS2 precursor in a feature level (for GNPS FBMN).        

rule precursorcorrection_feature:
    input:
        var1= "results/Interim/mzML/PCpeak_{samples}.mzML",
        var2= "results/Interim/preprocessed/FFM_{samples}.featureXML"
    output:
        "results/Interim/mzML/PCfeature_{samples}.mzML"
    shell:
        """
        HighResPrecursorMassCorrector -in {input.var1} -feature:in {input.var2} -out {output}  -nearest_peak:mz_tolerance "100.0"
        """ 

# 4) MapAlignerPoseClustering is used to perform a linear retention time alignment, to correct for linear shifts in retention time between different runs.

rule MapAlignerPoseClustering:
    input:
        expand("results/Interim/preprocessed/FFM_{samples}.featureXML", samples=SAMPLES)
    output:
        var1= expand("results/Interim/preprocessed/MapAligned_{samples}.featureXML", samples=SAMPLES),
        var2= expand("results/Interim/preprocessed/MapAligned_{samples}.trafoXML", samples=SAMPLES)
    shell:
        """
        MapAlignerPoseClustering -algorithm:max_num_peaks_considered -1 -algorithm:superimposer:mz_pair_max_distance 0.05 -algorithm:pairfinder:distance_MZ:max_difference 10.0 -algorithm:pairfinder:distance_MZ:unit ppm -in {input} -out {output.var1} -trafo_out {output.var2}
        """ 

# 5) The FeatureLinkerUnlabeledKD is used to aggregate the feature information (from single files) into a ConsensusFeature, linking features from different files together, which have a similar m/z and rt (MS1 level).

rule FeatureLinkerUnlabeledKD:
    input:
        expand("results/Interim/preprocessed/MapAligned_{samples}.featureXML", samples=SAMPLES)
    output:
        "results/Interim/preprocessed/preprocessed.consensusXML"
    threads: 4
    shell:
        """
        FeatureLinkerUnlabeledKD -in {input} -out {output} -threads {threads}
        """

