rule sirius_csifingerid:
    input: 
        "results/preprocessed/{samples}.featureXML",
        "results/precursorcorrection/{samples}.mzML", 
        "resources/Sirius/sirius/bin/sirius",
        "data/mzML/{samples}.mzML"
    output:
        "results/SIRIUS/{samples}.mzTab",
        "results/CSIFingerID/{samples}.mzTab"
    conda:
        "../envs/pyopenms.yaml"
    script:
        "../scripts/6_sirius_csifingerid.py"
        
rule df_sirius_csifingerid:
    input: 
        "results/SIRIUS/{samples}.mzTab",
        "results/CSIFingerID/{samples}.mzTab"
    output:
        "results/dataframes/sirius_output/sirius_{samples}.csv",
        "results/dataframes/CSIFingerID_output/csi_{samples}.csv"
    conda:
        "../envs/pyopenms.yaml"
    script:
        "../scripts/7_df_sirius_csifingerid.py"
    