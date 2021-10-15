rule sirius:
    input: 
        "results/{samples}/interim/MFD_{samples}.featureXML",
        "results/{samples}/interim/precursorcorrected_{samples}.mzML", 
        "resources/Sirius/sirius.app/Contents/MacOS/sirius",
    output:
        "results/{samples}/interim/formulas_{samples}.mzTab"
    conda:
        "../envs/pyopenms.yaml"
    script:
        "../scripts/4_sirius.py"

rule df_sirius:
    input: 
        "results/{samples}/interim/formulas_{samples}.mzTab"
    output:
        "results/{samples}/formulas_{samples}.csv"
    conda:
        "../envs/pyopenms.yaml"
    script:
        "../scripts/5_df_sirius.py"

