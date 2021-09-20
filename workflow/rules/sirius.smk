rule sirius:
    input: 
        "results/{samples}/interim/preprocessed_{samples}.featureXML",
        "results/{samples}/interim/precursorcorrected_{samples}.mzML", 
        "resources/Sirius/sirius/bin/sirius",
        "results/{samples}/interim/{samples}.mzML"
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

