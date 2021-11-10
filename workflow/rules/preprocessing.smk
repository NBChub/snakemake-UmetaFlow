#1) PrecursorCorrection (To the "highest intensity MS1 peak")

rule precursorcorrection:
    input:
        "results/{samples}/interim/{samples}.mzML"
    output:
        "results/{samples}/interim/precursorcorrected_{samples}.mzML" 
    conda:
        "../envs/pyopenms.yaml"   
    script:
        "../scripts/1_precursorcorrection.py"

# 2) Preprocessing: Feature finding algorithm that detects peaks 
    
rule preprocess:
    input:
        "results/{samples}/interim/precursorcorrected_{samples}.mzML"
    output:
        "results/{samples}/interim/preprocessed/FFM_{samples}.featureXML",
        "results/{samples}/interim/preprocessed/MFD_{samples}.featureXML"
    conda:
        "../envs/pyopenms.yaml"
    script:
        "../scripts/2_preprocessing.py"

rule df_preprocess:
    input: 
        "results/{samples}/interim/preprocessed/MFD_{samples}.featureXML"
    output:
        "results/{samples}/features_{samples}.csv"
    conda:
        "../envs/pyopenms.yaml"
    script:
        "../scripts/3_df_preprocess.py"

