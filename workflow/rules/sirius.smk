# 1) Correction of wrong assignment of the mono-isotopic mass for precursors (MS2 level)

rule preprocess_noconvexhulls_sirius:
    input:
        "results/{samples}/interim/PCpeak_{samples}.mzML"
    output:
        "results/{samples}/interim/sirius/FFM_nch_{samples}.featureXML"
    shell:
        """
        OpenMS/OpenMS-build/bin/FeatureFinderMetabo -in {input} -out {output} -algorithm:common:noise_threshold_int "1.0e04" -algorithm:mtd:mass_error_ppm "10.0" -algorithm:epd:width_filtering "fixed" -algorithm:ffm:isotope_filtering_model "none" -algorithm:ffm:remove_single_traces "true"
        """

# 2) Decharger: Decharging algorithm for adduct assignment

rule sirius_decharge:
    input:
        "results/{samples}/interim/sirius/FFM_nch_{samples}.featureXML"
    output:
        "results/{samples}/interim/sirius/MFD_nch_{samples}.featureXML"
    shell:
        """
        OpenMS/OpenMS-build/bin/MetaboliteAdductDecharger -in {input} -out_fm {output} -algorithm:MetaboliteFeatureDeconvolution:potential_adducts "H:+:0.6" "Na:+:0.2" "NH4:+:0.2" "H-2O-1:0:0.05" -algorithm:MetaboliteFeatureDeconvolution:charge_max "1" -algorithm:MetaboliteFeatureDeconvolution:charge_span_max "1"  -algorithm:MetaboliteFeatureDeconvolution:max_neutrals "1"
        """

# 3) Correct the MS2 precursor in a feature level.        

rule precursorcorrection_feature_sirius:
    input:
        var1= "results/{samples}/interim/PCpeak_{samples}.mzML",
        var2= "results/{samples}/interim/sirius/MFD_nch_{samples}.featureXML"
    output:
        "results/{samples}/interim/sirius/PCfeature_nch_{samples}.mzML"
    shell:
        """
        OpenMS/OpenMS-build/bin/HighResPrecursorMassCorrector -in {input.var1} -feature:in {input.var2} -out {output}  -nearest_peak:mz_tolerance "100.0"
        """ 

# 4) SIRIUS generates formula prediction from scores calculated from 1) MS2 fragmentation scores (ppm error + intensity) and 2) MS1 isotopic pattern scores.        

rule sirius:
    input: 
        var1= "resources/Sirius/sirius/bin/sirius",
        var2= "results/{samples}/interim/sirius/PCfeature_nch_{samples}.mzML", 
        var3= "results/{samples}/interim/sirius/MFD_nch_{samples}.featureXML"        
    output:
        "results/{samples}/interim/sirius/formulas_{samples}.mzTab"
    shell:
        """
        OpenMS/OpenMS-build/bin/SiriusAdapter -executable {input.var1} -in {input.var2} -in_featureinfo {input.var3} -out_sirius {output} -preprocessing:filter_by_num_masstraces 2 -preprocessing:feature_only -sirius:profile orbitrap -sirius:db none -sirius:ions_considered "[M+H]+, [M-H2O+H]+, [M+Na]+, [M+NH4]+" -sirius:elements_enforced CHN[15]O[57]S[4] -project:processors 2 -debug 3 
        """

# 5) Convert the mzTab to a csv file

rule df_sirius:
    input: 
        "results/{samples}/interim/sirius/formulas_{samples}.mzTab"
    output:
        "results/{samples}/formulas_{samples}.csv"
    conda:
        "../envs/file_conversion.yaml"
    script:
        "/datadrive/snakemake/snakemake-metabolomics/workflow/scripts/df_sirius.py"

