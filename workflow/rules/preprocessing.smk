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
        exp="results/precursorcorrection/{samples}.mzML"
    output:
        "results/preprocessed/{samples}.featureXML" 
    conda:
        "../envs/pyopenms.yaml"
    script:
        "../scripts/2_preprocessing.py"

rule df_preprocess:
    input: 
        featureinfo= "results/preprocessed/{samples}.featureXML",
    output:
        "results/dataframes/features_tables/{features}.csv",
    conda:
        "../envs/pyopenms.yaml"
    script:
        "../scripts/3_df_preprocess.py"

