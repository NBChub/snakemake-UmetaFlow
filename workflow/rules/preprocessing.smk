rule precursorcorrection:
    input:
        "data/mzML/{samples}.mzML"
    output:
        "results/precursorcorrection/{samples}.mzML" 
    conda:
        "../envs/pyopenms.yaml"   
    script:
        "../scripts/1_precursorcorrection.py"

rule preprocess:
    input:
        "results/precursorcorrection/{samples}.mzML"
    output:
        "results/preprocessed/{samples}.featureXML" 
    conda:
        "../envs/pyopenms.yaml"
    script:
        "../scripts/2_preprocessing.py"

rule df_preprocess:
    input: 
        "results/preprocessed/{samples}.featureXML",
    output:
        "results/dataframes/features_tables/{samples}.csv",
    conda:
        "../envs/pyopenms.yaml"
    script:
        "../scripts/3_df_preprocess.py"

