# Re-quantify the features in all data (missing value correction)
# This rule is currently not used - in progress

rule rename:
    input:
        "results/GNPSexport/interim/consensus.tsv"
    output:
        "results/GNPSexport/interim/consensus.csv"
    shell:
        """
        mv {input} {output}
        """

rule build_library:
    input:
        "results/GNPSexport/interim/consensus.csv"
    output:
        "resources/MetaboliteIdentification.tsv"
    conda:
        "../envs/pyopenms.yaml"   
    script:
        "../scripts/8_metaboliteidentidication.py"

rule metaboident:
    input:
        "resources/MetaboliteIdentification.tsv",
        "results/{samples}/interim/precursorcorrected_{samples}.mzML"
    output:
        "results/Requant/interim/FFMID_{samples}.featureXML"
    shell:
        """
        resources/OpenMS-2.7.0/bin/FeatureFinderMetaboIdent -id {input[0]} -extract:mz_window 5.0 -in {input[1]} -out {output}
        """
# Export the consensusXML file to a txt file for GNPS
rule FFMI_df:
    input:
        "results/Requant/interim/FFMID_{samples}.featureXML"
    output:
        "results/Requant/FFMID_{samples}.csv" 
    shell:
        """
        resources/OpenMS-2.7.0/bin/TextExporter -in {input} -out {output}
        """