# 1) Correction of wrong assignment of the mono-isotopic mass for precursors (MS2 level)

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

# 2) SIRIUS generates formula prediction from scores calculated from 1) MS2 fragmentation scores (ppm error + intensity) and 2) MS1 isotopic pattern scores.        
#    The CSI_fingerID function is another algorithm from the Boecher lab, just like SIRIUS adapter and is using the formula predictions from SIRIUS, to search in structural libraries and predict the structure of each formula.

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
    