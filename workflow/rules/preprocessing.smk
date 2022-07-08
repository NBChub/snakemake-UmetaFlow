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
        "results/Interim/Preprocessed/FFM_{samples}.featureXML"
    threads: 4
    shell:
        """
        /Users/eeko/openms-develop/openms_build/bin/FeatureFinderMetabo -in {input} -out {output} -algorithm:common:noise_threshold_int "1.0e04" -algorithm:mtd:mass_error_ppm "10.0" -algorithm:epd:width_filtering "fixed" -algorithm:ffm:isotope_filtering_model "none" -algorithm:ffm:remove_single_traces "true" -algorithm:ffm:report_convex_hulls "true" -threads {threads}
        """

# 3) Correct the MS2 precursor in a feature level (for GNPS FBMN).        

rule precursorcorrection_feature:
    input:
        var1= "results/Interim/mzML/PCpeak_{samples}.mzML",
        var2= "results/Interim/Preprocessed/FFM_{samples}.featureXML"
    output:
        "results/Interim/mzML/PCfeature_{samples}.mzML"
    shell:
        """
        /Users/eeko/openms-develop/openms_build/bin/HighResPrecursorMassCorrector -in {input.var1} -feature:in {input.var2} -out {output}  -nearest_peak:mz_tolerance "100.0"
        """ 

# 4) (i) MapAlignerPoseClustering is used to perform a linear retention time alignment, to correct for linear shifts in retention time between different runs.

rule MapAligner:
    input:
        expand("results/Interim/Preprocessed/FFM_{samples}.featureXML", samples=SAMPLES)
    output:
        var1= expand("results/Interim/Preprocessed/MapAligned_{samples}.featureXML", samples=SAMPLES),
        var2= expand("results/Interim/Preprocessed/MapAligned_{samples}.trafoXML", samples=SAMPLES)
    shell:
        """
        /Users/eeko/openms-develop/openms_build/bin/MapAlignerPoseClustering -algorithm:max_num_peaks_considered -1 -algorithm:superimposer:mz_pair_max_distance 0.05 -algorithm:pairfinder:distance_MZ:max_difference 10.0 -algorithm:pairfinder:distance_MZ:unit ppm -in {input} -out {output.var1} -trafo_out {output.var2}
        """ 

# 4) (ii) MapRTTransformer is used to perform a linear retention time alignment, to correct for linear shifts in retention time between different runs using the transformation files from the reprocessing rule MapAlignerPoseClustering (faster computationally)

rule mzMLaligner:
    input:
        var1= "results/Interim/mzML/PCfeature_{samples}.mzML",
        var2= "results/Interim/Preprocessed/MapAligned_{samples}.trafoXML"
    output:
        "results/GNPSexport/mzML/Aligned_{samples}.mzML"
    shell:
        """
        /Users/eeko/openms-develop/openms_build/bin/MapRTTransformer -in {input.var1} -trafo_in {input.var2} -out {output}
        """ 

# 5) Decharger: Decharging algorithm for adduct assignment

rule adduct_annotations_FFM:
    input:
        "results/Interim/Preprocessed/MapAligned_{samples}.featureXML"
    output:
        "results/Interim/Preprocessed/MFD_{samples}.featureXML" 
    shell:
        """
        /Users/eeko/openms-develop/openms_build/bin/MetaboliteAdductDecharger -in {input} -out_fm {output} -algorithm:MetaboliteFeatureDeconvolution:potential_adducts "H:+:0.6" "Na:+:0.1" "NH4:+:0.1" "H-1O-1:+:0.1" "H-3O-2:+:0.1" -algorithm:MetaboliteFeatureDeconvolution:charge_max "1" -algorithm:MetaboliteFeatureDeconvolution:charge_span_max "1"  -algorithm:MetaboliteFeatureDeconvolution:max_neutrals "1"
        """    

<<<<<<< HEAD
# 6) Introduce the features to a protein identification file (idXML)- the only way to annotate MS2 spectra for GNPS FBMN  

rule IDMapper_FFM:
    input:
        var1= "resources/emptyfile.idXML",
        var2= "results/Interim/Preprocessed/MFD_{samples}.featureXML",
        var3= "results/GNPSexport/mzML/Aligned_{samples}.mzML"
    output:
        "results/Interim/Preprocessed/IDMapper_{samples}.featureXML"
    shell:
        """
        /Users/eeko/openms-develop/openms_build/bin/IDMapper -id {input.var1} -in {input.var2}  -spectra:in {input.var3} -out {output} 
        """

# 7) The FeatureLinkerUnlabeledKD is used to aggregate the feature information (from single files) into a ConsensusFeature, linking features from different files together, which have a similar m/z and rt (MS1 level).

rule FeatureLinkerUnlabeledKD:
    input:
        expand("results/Interim/Preprocessed/IDMapper_{samples}.featureXML", samples=SAMPLES)
    output:
        "results/Interim/Preprocessed/Preprocessed.consensusXML"
    threads: 4
    shell:
        """
        /Users/eeko/openms-develop/openms_build/bin/FeatureLinkerUnlabeledKD -in {input} -out {output} -threads {threads}
        """

# 8) export the consensusXML file to a tsv file to produce a single matrix for PCA
=======
# 6) The FeatureLinkerUnlabeledKD is used to aggregate the feature information (from single files) into a ConsensusFeature, linking features from different files together, which have a similar m/z and rt (MS1 level).

rule FeatureLinkerUnlabeledKD:
    input:
        expand("results/Interim/Preprocessed/MapAligned_{samples}.featureXML", samples=SAMPLES)
    output:
        "results/Interim/Preprocessed/Preprocessed.consensusXML"
    threads: 4
    shell:
        """
        /Users/eeko/openms-develop/openms_build/bin/FeatureLinkerUnlabeledKD -in {input} -out {output} -threads {threads}
        """

# 7) export the consensusXML file to a tsv file to produce a single matrix for PCA
>>>>>>> 8c5ccedc94fc856850dde76e921e22a705428575

rule FFM_matrix:
    input:
        "results/Interim/Preprocessed/Preprocessed.consensusXML"
    output:
        "results/Preprocessed/FeatureMatrix.tsv"
    conda:
        "../envs/pyopenms.yaml"
    script:
        "../scripts/cleanup.py"