# 1) copy all the original mzml files (precursor-corrected ones) in the GNPSExport folder for easier use
rule FileCopy:
    input:
        "results/Interim/mzML/PCfeature_{samples}.mzML"
    output:
        "results/GNPSexport/{samples}.mzML"
    shell:
        """
        cp {input} {output}
        """ 


# 2) Filter out the features that do not have an MS2 pattern (no protein ID annotations)

rule FileFilter:
    input:
        "results/Interim/preprocessed/FeatureLinkerUnlabeledKD.consensusXML"
    output:
        "results/Interim/GNPSexport/filtered.consensusXML"
    shell:
        """
        /Users/eeko/openms-develop/openms_build/bin/FileFilter -id:remove_unannotated_features -in {input} -out {output} 
        """

# 3) GNPS_export creates an mgf file with only the MS2 information of all files (introduce mzml files with spaces between them)

rule GNPS_export:
    input:
        var1= "results/Interim/GNPSexport/filtered.consensusXML",
        var2= expand("results/Interim/mzML/PCfeature_{samples}.mzML", samples=SAMPLES)
    output:
        "results/GNPSexport/MSMS.mgf" 
    shell:
        """
        /Users/eeko/openms-develop/openms_build/bin/GNPSExport -ini resources/GNPSExport.ini -in_cm {input.var1} -in_mzml {input.var2} -out {output} 
        """

# 4) export the consensusXML file to a txt file for GNPS

rule GNPS_txt_export:
    input:
        "results/Interim/GNPSexport/filtered.consensusXML"
    output:
        "results/GNPSexport/FeatureQuantificationTable.txt" 
    shell:
        """
        /Users/eeko/openms-develop/openms_build/bin/TextExporter -in {input} -out {output}
        """
