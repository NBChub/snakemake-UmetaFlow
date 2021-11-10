# 1) Correction of wrong assignment of the mono-isotopic mass for precursors (MS2 level)

rule preprocess_noconvexhulls:
    input:
        "results/{samples}/interim/precursorcorrected_{samples}.mzML"
    output:
        "results/{samples}/interim/sirius/FFM_noconvexhulls_{samples}.featureXML",
        "results/{samples}/interim/sirius/MFD_noconvexhulls_{samples}.featureXML"
    conda:
        "../envs/pyopenms.yaml"
    script:
        "../scripts/2_noconvexhulls_preprocessing.py"

# 2) SIRIUS generates formula prediction from scores calculated from 1) MS2 fragmentation scores (ppm error + intensity) and 2) MS1 isotopic pattern scores.        

rule sirius:
    input: 
        "results/{samples}/interim/sirius/MFD_noconvexhulls_{samples}.featureXML",
        "results/{samples}/interim/precursorcorrected_{samples}.mzML", 
        "resources/Sirius/sirius/bin/sirius",
    output:
        "results/{samples}/interim/sirius/formulas_{samples}.mzTab"
    conda:
        "../envs/pyopenms.yaml"
    script:
        "../scripts/4_sirius.py"

rule df_sirius:
    input: 
        "results/{samples}/interim/sirius/formulas_{samples}.mzTab"
    output:
        "results/{samples}/formulas_{samples}.csv"
    conda:
        "../envs/pyopenms.yaml"
    script:
        "../scripts/5_df_sirius.py"

