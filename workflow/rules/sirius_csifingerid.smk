rule sirius_csifingerid:
    input: 
        "results/{samples}/interim/MFD_{samples}.featureXML",
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
    