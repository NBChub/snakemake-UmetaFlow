# 1) Filter out the features that do not have an MS2 pattern (no protein ID annotations)

<<<<<<< HEAD
if config["rules"]["requantification"]==True:
    rule FileFilter:
        input:
            "results/Interim/Requantified/Requantified.consensusXML"
        output:
            "results/Interim/GNPSexport/filtered.consensusXML"
        shell:
            """
            /Users/eeko/openms-develop/openms_build/bin/FileFilter -id:remove_unannotated_features -in {input} -out {output} 
            """
else:            
    rule FileFilter:
        input:
            "results/Interim/Preprocessed/Preprocessed.consensusXML"
        output:
            "results/Interim/GNPSexport/filtered.consensusXML"
        shell:
            """
            /Users/eeko/openms-develop/openms_build/bin/FileFilter -id:remove_unannotated_features -in {input} -out {output} 
            """        
=======
rule FileFilter:
    input:
        "results/Interim/Requantified/Requantified.consensusXML"
    output:
        "results/Interim/GNPSexport/filtered.consensusXML"
    shell:
        """
        /Users/eeko/openms-develop/openms_build/bin/FileFilter -id:remove_unannotated_features -in {input} -out {output} 
        """
>>>>>>> 8c5ccedc94fc856850dde76e921e22a705428575

# 2) GNPS_export creates an mgf file with only the MS2 information of all files (introduce mzml files with spaces between them)

rule GNPS_export:
    input:
        var1= "results/Interim/GNPSexport/filtered.consensusXML",
        var2= expand("results/GNPSexport/mzML/Aligned_{samples}.mzML", samples=SAMPLES)
    output:
        out1= "results/GNPSexport/MSMS.mgf",
        out2= "results/GNPSexport/FeatureQuantificationTable.txt", 
        out3= "results/GNPSexport/SuppPairs.csv",
        out4= "results/GNPSexport/metadata.tsv"
    shell:
        """
        /Users/eeko/openms-develop/openms_build/bin/GNPSExport -in_cm {input.var1} -in_mzml {input.var2} -out {output.out1} -out_quantification {output.out2} -out_pairs {output.out3} -out_meta_values {output.out4}
        """