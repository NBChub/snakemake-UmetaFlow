rule preprocess_noconvexhulls:
    input:
        "results/{samples}/interim/precursorcorrected_{samples}.mzML"
    output:
        "results/{samples}/interim/sirius_csifingerid/FFM_noconvexhulls_{samples}.featureXML",
        "results/{samples}/interim/sirius_csifingerid/MFD_noconvexhulls_{samples}.featureXML"
    conda:
        "../envs/pyopenms.yaml"
    script:
        "../scripts/2_noconvexhulls_preprocessing.py"

rule sirius_csifingerid:
    input: 
        "results/{samples}/interim/sirius_csifingerid/MFD_noconvexhulls_{samples}.featureXML",
        "results/{samples}/interim/precursorcorrected_{samples}.mzML", 
        "resources/Sirius/sirius.app/Contents/MacOS/sirius"
    output:
        "results/{samples}/interim/sirius_csifingerid/formulas_{samples}.mzTab",
        "results/{samples}/interim/sirius_csifingerid/structures_{samples}.mzTab"
    conda:
        "../envs/pyopenms.yaml"
    script:
        "../scripts/6_sirius_csifingerid.py"
        
rule df_sirius_csifingerid:
    input: 
        "results/{samples}/interim/sirius_csifingerid/formulas_{samples}.mzTab",
        "results/{samples}/interim/sirius_csifingerid/structures_{samples}.mzTab"
    output:
        "results/{samples}/sirius_csi_formulas_{samples}.csv",
        "results/{samples}/sirius_csi_structures_{samples}.csv"
    conda:
        "../envs/pyopenms.yaml"
    script:
        "../scripts/7_df_sirius_csifingerid.py"
    