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

# 2) GNPS_export creates an mgf file with only the MS2 information of all files (introduce mzml files with spaces between them)

rule GNPS_export:
    input:
        var1= "results/Interim/GNPSexport/filtered.consensusXML",
        var2= expand("results/Interim/mzML/PCfeature_{samples}.mzML", samples=SAMPLES)
    output:
        out1= "results/GNPSexport/MSMS.mgf",
        out2= "results/GNPSexport/FeatureQuantificationTable.txt", 
        out3= "results/GNPSexport/SuppPairs.csv"
    shell:
        """
        GNPSExport -in_cm {input.var1} -in_mzml {input.var2} -out {output.out1} -out_quantification {output.out2} -out_pairs {output.out3}
        """