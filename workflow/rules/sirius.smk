# 1) Correction of wrong assignment of the mono-isotopic mass for precursors (MS2 level)

rule preprocess_noconvexhulls:
    input:
        "results/{samples}/interim/precursorcorrected_{samples}.mzML"
    output:
        "results/{samples}/interim/sirius/FFM_noconvexhulls_{samples}.featureXML",
        "results/{samples}/interim/sirius/MFD_noconvexhulls_{samples}.featureXML"
    shell:
        """
        resources/OpenMS-2.7.0/bin/FeatureFinderMetabo -in {input} -out {output} -algorithm:common:noise_threshold_int "1.0e04" -algorithm:mtd:mass_error_ppm "10.0" -algorithm:epd:width_filtering "fixed" -algorithm:ffm:isotope_filtering_model "none" -algorithm:ffm:remove_single_traces "true"
        """

# 2) SIRIUS generates formula prediction from scores calculated from 1) MS2 fragmentation scores (ppm error + intensity) and 2) MS1 isotopic pattern scores.        

rule sirius:
    input: 
        "resources/Sirius/sirius.app/Contents/MacOS/sirius",
        "results/{samples}/interim/precursorcorrected_{samples}.mzML", 
        "results/{samples}/interim/sirius/MFD_noconvexhulls_{samples}.featureXML"        
    output:
        "results/{samples}/interim/sirius/formulas_{samples}.mzTab"
    shell:
        """
        resources/OpenMS-2.7.0/bin/SiriusAdapter -executable {input[0]} -in {input[1]} -in_featureinfo {input[2]} -out_sirius {output[0]} -preprocessing:filter_by_num_masstraces 2 -preprocessing:feature_only -sirius:profile orbitrap -sirius:db none -sirius:ions_considered "[M+H]+, [M-H2O+H]+, [M+Na]+, [M+NH4]+" -sirius:elements_enforced CHN[15]OS[4] -project:processors 2 -debug 3 
        """

rule df_sirius:
    input: 
        "results/{samples}/interim/sirius/formulas_{samples}.mzTab"
    output:
        "results/{samples}/formulas_{samples}.csv"
    conda:
        "../envs/pyopenms.yaml"
    script:
        "../scripts/5_df_sirius.py"

