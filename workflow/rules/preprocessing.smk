#1) PrecursorCorrection (To the "highest intensity MS1 peak")

rule precursorcorrection:
    input:
        "results/{samples}/interim/{samples}.mzML"
    output:
        "results/{samples}/interim/precursorcorrected_{samples}.mzML",
        "results/{samples}/interim/precursorcorrected_{samples}.csv" 
    shell:
        """
        resources/OpenMS-2.7.0/bin/HighResPrecursorMassCorrector -in {input} -out {output[0]} -highest_intensity_peak:mz_tolerance "100.0"
        """

# 2) Preprocessing: Feature finding algorithm that detects peaks 
    
rule preprocess:
    input:
        "results/{samples}/interim/precursorcorrected_{samples}.mzML"
    output:
        "results/{samples}/interim/preprocessed/FFM_{samples}.featureXML"
    shell:
        """
        resources/OpenMS-2.7.0/bin/FeatureFinderMetabo -in {input} -out {output} -algorithm:common:noise_threshold_int "1.0e04" -algorithm:mtd:mass_error_ppm "10.0" -algorithm:epd:width_filtering "fixed" -algorithm:ffm:isotope_filtering_model "none" -algorithm:ffm:remove_single_traces "true" -algorithm:ffm:report_convex_hulls "true"
        """

# 3) Deacharger: Decharging algorithm for adduct assignment

rule decharge:
    input:
        "results/{samples}/interim/preprocessed/FFM_{samples}.featureXML"
    output:
        "results/{samples}/interim/preprocessed/MFD_{samples}.featureXML"
    shell:
        """
        resources/OpenMS-2.7.0/bin/MetaboliteAdductDecharger -in {input} -out_fm {output} -algorithm:MetaboliteFeatureDeconvolution:potential_adducts "H:+:0.6" "Na:+:0.2" "NH4:+:0.1" "H2O:-:0.1"  -algorithm:MetaboliteFeatureDeconvolution:charge_max "1" -algorithm:MetaboliteFeatureDeconvolution:charge_span_max "1"  -algorithm:MetaboliteFeatureDeconvolution:max_neutrals "1"
        """

rule df_preprocess:
    input: 
        "results/{samples}/interim/preprocessed/MFD_{samples}.featureXML"
    output:
        "results/{samples}/features_{samples}.csv"
    shell:
        """
        resources/OpenMS-2.7.0/bin/TextExporter -in {input} -out {output}
        """

