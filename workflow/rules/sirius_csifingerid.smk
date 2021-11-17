# 1) Correction of wrong assignment of the mono-isotopic mass for precursors (MS2 level)

rule preprocess_noconvexhulls:
    input:
        "results/{samples}/interim/{samples}.mzML"
    output:
        "results/{samples}/interim/sirius_csifingerid/FFM_noconvexhulls_{samples}.featureXML"
    shell:
        """
        resources/OpenMS-2.7.0/bin/FeatureFinderMetabo -in {input} -out {output} -algorithm:common:noise_threshold_int "1.0e04" -algorithm:mtd:mass_error_ppm "10.0" -algorithm:epd:width_filtering "fixed" -algorithm:ffm:isotope_filtering_model "none" -algorithm:ffm:remove_single_traces "true"
        """

# 2) Decharger: Decharging algorithm for adduct assignment

rule sirius_decharge:
    input:
        "results/{samples}/interim/sirius_csifingerid/FFM_noconvexhulls_{samples}.featureXML"
    output:
        "results/{samples}/interim/sirius_csifingerid/MFD_noconvexhulls_{samples}.featureXML"
    shell:
        """
        resources/OpenMS-2.7.0/bin/MetaboliteAdductDecharger -in {input} -out_fm {output} -algorithm:MetaboliteFeatureDeconvolution:potential_adducts "H:+:0.6" "Na:+:0.2" "NH4:+:0.1" "H2O:-:0.1"  -algorithm:MetaboliteFeatureDeconvolution:charge_max "1" -algorithm:MetaboliteFeatureDeconvolution:charge_span_max "1"  -algorithm:MetaboliteFeatureDeconvolution:max_neutrals "1"
        """

# 3) Correct the MS2 precursor in a feature level.        

rule precursorcorrection_feature:
    input:
        "results/{samples}/interim/{samples}.mzML",
        "results/{samples}/interim/sirius_csifingerid/MFD_noconvexhulls_{samples}.featureXML"
    output:
        "results/{samples}/interim/sirius_csifingerid/PrecursorMassCorrector_{samples}.mzML"
    shell:
        """
        resources/OpenMS-2.7.0/bin/HighResPrecursorMassCorrector -in {input[0]} -feature:in {input[1]} -out {output} -nearest_peak:mz_tolerance "100.0"
        """ 

# 4) SIRIUS generates formula prediction from scores calculated from 1) MS2 fragmentation scores (ppm error + intensity) and 2) MS1 isotopic pattern scores.        
#    The CSI_fingerID function is another algorithm from the Boecher lab, just like SIRIUS adapter and is using the formula predictions from SIRIUS, to search in structural libraries and predict the structure of each formula.

rule sirius_csifingerid:
    input: 
        "resources/Sirius/sirius.app/Contents/MacOS/sirius",
        "results/{samples}/interim/sirius_csifingerid/PrecursorMassCorrector_{samples}.mzML", 
        "results/{samples}/interim/sirius_csifingerid/MFD_noconvexhulls_{samples}.featureXML" 
    output:
        "results/{samples}/interim/sirius_csifingerid/formulas_{samples}.mzTab",
        "results/{samples}/interim/sirius_csifingerid/structures_{samples}.mzTab"
    shell:
        """
        resources/OpenMS-2.7.0/bin/SiriusAdapter -executable {input[0]} -in {input[1]} -in_featureinfo {input[2]} -out_sirius {output[0]} -out_fingerid {output[1]} -preprocessing:filter_by_num_masstraces 2 -preprocessing:feature_only -sirius:profile orbitrap -sirius:db none -sirius:ions_considered "[M+H]+, [M-H2O+H]+, [M+Na]+, [M+NH4]+" -sirius:elements_enforced CHN[15]OS[4] -fingerid:db BIO -project:processors 2 -debug 3 -fingerid:candidates 10
        """

# 5) Convert the mzTab to a csv file

rule df_sirius_csifingerid:
    input: 
        "results/{samples}/interim/sirius_csifingerid/formulas_{samples}.mzTab",
        "results/{samples}/interim/sirius_csifingerid/structures_{samples}.mzTab"
    output:
        "results/{samples}/sirius_csi_formulas_{samples}.csv",
        "results/{samples}/sirius_csi_structures_{samples}.csv"
    conda:
        "../envs/file_conversion.yaml"
    script:
        "../scripts/df_sirius_csifingerid.py"
    