# 1) Filter out the features that do not have an MS2 pattern (no protein ID annotations)

rule FileFilter:
    input:
        "results/Interim/preprocessed/FeatureLinkerUnlabeledKD.consensusXML"
    output:
        "results/Interim/GNPSexport/filtered.consensusXML"
    shell:
        """
        FileFilter -id:remove_unannotated_features -in {input} -out {output} 
        """

# 1) GNPS_export creates an mgf file with only the MS2 information of all files (introduce mzml files with spaces between them)

rule GNPS_export:
    input:
        var1= "results/Interim/GNPSexport/filtered.consensusXML",
        var2= expand("results/Interim/mzML/PCfeature_{samples}.mzML", samples=SAMPLES)
    output:
        "results/GNPSexport/MSMS.mgf" 
    shell:
        """
        GNPSExport -ini resources/GNPSExport.ini -in_cm {input.var1} -in_mzml {input.var2} -out {output} 
        """

# 3) export the consensusXML file to a txt file for GNPS

rule GNPS_txt_export:
    input:
        "results/Interim/GNPSexport/filtered.consensusXML"
    output:
        "results/GNPSexport/FeatureQuantificationTable.txt" 
    shell:
        """
        TextExporter -in {input} -out {output}
        """
