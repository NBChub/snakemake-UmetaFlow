
rule metaboident:
    input:
        "resources/MetaboliteIdentification.tsv",
        "results/{samples}/interim/precursorcorrected_{samples}.mzML"
    output:
        "results/MetaboIdentification/FFMID_{samples}.featureXML"
    shell:
        """
        resources/OpenMS-2.7.0/bin/FeatureFinderMetaboIdent -id {input[0]} -in {input[1]} -out {output}
        """

rule df_metaboident:
    input:
        "results/MetaboIdentification/FFMID_{samples}.featureXML"
    output:
        "results/MetaboIdentification/FFMID_{samples}.csv"
    conda:
        "../envs/pyopenms.yaml"   
    script:
        "../scripts/8_df_metaboident.py"