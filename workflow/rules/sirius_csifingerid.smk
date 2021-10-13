rule preprocess_noconvexhulls:
    input:
        "results/{samples}/interim/precursorcorrected_{samples}.mzML"
    output:
        "results/{samples}/interim/MFD_no_convex_hulls_{samples}.featureXML"
    conda:
        "../envs/pyopenms.yaml"
    script:
        "../scripts/2_sirius_preprocessing.py"

rule sirius_csifingerid:
    input: 
        "results/{samples}/interim/MFD_no_convex_hulls_{samples}.featureXML",
        "results/{samples}/interim/precursorcorrected_{samples}.mzML", 
        "resources/Sirius/sirius.app/Contents/MacOS/sirius"
    output:
        "results/{samples}/interim/formulas_{samples}.mzTab",
        "results/{samples}/interim/structures_{samples}.mzTab"
    conda:
        "../envs/pyopenms.yaml"
    script:
        "../scripts/6_sirius_csifingerid.py"
        
rule df_sirius_csifingerid:
    input: 
        "results/{samples}/interim/formulas_{samples}.mzTab",
        "results/{samples}/interim/structures_{samples}.mzTab"
    output:
        "results/{samples}/formulas_{samples}.csv",
        "results/{samples}/structures_{samples}.csv"
    conda:
        "../envs/pyopenms.yaml"
    script:
        "../scripts/7_df_sirius_csifingerid.py"
    