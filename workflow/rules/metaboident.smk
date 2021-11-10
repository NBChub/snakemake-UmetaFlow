# Re-quantify the features in all data (missing value correction)
# This rule is currently not used - in progress

rule metaboident:
    input:
        "resources/MetaboliteIdentification.tsv",
        "results/{samples}/interim/precursorcorrected_{samples}.mzML"
    output:
        "results/MetaboliteAnnotation/interim/FFMID_{samples}.featureXML"
    shell:
        """
        resources/OpenMS-2.7.0/bin/FeatureFinderMetaboIdent -id {input[0]} -extract:mz_window 5.0 -in {input[1]} -out {output}
        """

rule df_metaboident:
    input:
        "results/MetaboliteAnnotation/interim/FFMID_{samples}.featureXML"
    output:
        "results/MetaboliteAnnotation/FFMID_{samples}.csv"
    conda:
        "../envs/pyopenms.yaml"   
    script:
        "../scripts/8_df_FFMI.py"