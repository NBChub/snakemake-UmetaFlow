# 1) Correct the MS2 precursor on a peak level (To the "highest intensity MS1 peak")

rule precursorcorrection_peak:
    input:
        "results/{samples}/interim/{samples}.mzML"
    output:
        "results/{samples}/interim/PCpeak_{samples}.mzML"
    shell:
        """
        HighResPrecursorMassCorrector -in {input} -out {output} -highest_intensity_peak:mz_tolerance "100.0"
        """

# 2) Preprocessing: Feature finding algorithm that detects peaks 
    
rule preprocess:
    input:
        "results/{samples}/interim/PCpeak_{samples}.mzML"
    output:
        "results/{samples}/interim/preprocessed/FFM_{samples}.featureXML"
    shell:
        """
        FeatureFinderMetabo -in {input} -out {output} -algorithm:common:noise_threshold_int "1.0e04" -algorithm:mtd:mass_error_ppm "10.0" -algorithm:epd:width_filtering "fixed" -algorithm:ffm:isotope_filtering_model "none" -algorithm:ffm:remove_single_traces "true" -algorithm:ffm:report_convex_hulls "true"
        """

# 3) Decharger: Decharging algorithm for adduct assignment

rule decharge:
    input:
        "results/{samples}/interim/preprocessed/FFM_{samples}.featureXML"
    output:
        "results/{samples}/interim/preprocessed/MFD_{samples}.featureXML"
    shell:
        """
        MetaboliteAdductDecharger -in {input} -out_fm {output} -algorithm:MetaboliteFeatureDeconvolution:potential_adducts "H:+:0.6" "Na:+:0.2" "NH4:+:0.1" "H2O:-:0.1"  -algorithm:MetaboliteFeatureDeconvolution:charge_max "1" -algorithm:MetaboliteFeatureDeconvolution:charge_span_max "1"  -algorithm:MetaboliteFeatureDeconvolution:max_neutrals "1"
        """

# 4) Correct the MS2 precursor in a feature level (for GNPS FBMN).        

rule precursorcorrection_feature:
    input:
        var1= "results/{samples}/interim/PCpeak_{samples}.mzML",
        var2= "results/{samples}/interim/preprocessed/MFD_{samples}.featureXML"
    output:
        "results/{samples}/interim/PCfeature_{samples}.mzML"
    shell:
        """
        HighResPrecursorMassCorrector -in {input.var1} -feature:in {input.var2} -out {output}  -nearest_peak:mz_tolerance "100.0"
        """ 

# 5) Convert the featureXML to a csv file

rule df_preprocess:
    input: 
        "results/{samples}/interim/preprocessed/MFD_{samples}.featureXML"
    output:
        "results/{samples}/features_{samples}.csv"
    shell:
        """
        TextExporter -in {input} -out {output}
        """

