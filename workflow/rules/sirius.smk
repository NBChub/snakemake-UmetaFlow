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
        
rule sirius:
    input: 
        "results/{samples}/interim/sirius/MFD_noconvexhulls_{samples}.featureXML",
        "results/{samples}/interim/precursorcorrected_{samples}.mzML", 
        "resources/Sirius/sirius.app/Contents/MacOS/sirius",
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

