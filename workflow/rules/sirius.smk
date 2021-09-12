rule sirius:
    input: 
        "results/preprocessed/{samples}.featureXML",
        "results/precursorcorrection/{samples}.mzML", 
        "./resources/Contents/MacOS/sirius",
        "data/mzML/{samples}.mzML"
    output:
        "results/SIRIUS/{samples}.mzTab"
    conda:
        "../envs/pyopenms.yaml"
    script:
        "../scripts/4_sirius.py"

rule df_sirius:
    input: 
        "results/SIRIUS/{samples}.mzTab"
    output:
        "results/dataframes/sirius_output/sirius_{samples}.csv"
    conda:
        "../envs/pyopenms.yaml"
    script:
        "../scripts/5_df_sirius.py"

