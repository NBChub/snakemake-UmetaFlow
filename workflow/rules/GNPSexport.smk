# 1) Filter out the features that do not have an MS2 pattern (no protein ID annotations)

if config["rules"]["requantification"]==True:
    rule FileFilter:
        input:
            "results/Interim/Requantified/Requantified.consensusXML"
        output:
            "results/Interim/GNPSexport/filtered.consensusXML"
        log: "workflow/report/logs/GNPSexport/FileFilter.log"
        conda:
            "../envs/openms.yaml"
        shell:
            """
            FileFilter -id:remove_unannotated_features -in {input} -out {output} -log {log} 2>> {log}
            """
else:            
    rule FileFilter:
        input:
            "results/Interim/Preprocessed/Preprocessed.consensusXML"
        output:
            "results/Interim/GNPSexport/filtered.consensusXML"
        log: "workflow/report/logs/GNPSexport/FileFilter.log"
        conda:
            "../envs/openms.yaml"
        shell:
            """
            FileFilter -id:remove_unannotated_features -in {input} -out {output} -log {log} 2>> {log}
            """        

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
    log: "workflow/report/logs/GNPSexport/GNPS_export.log"
    conda:
        "../envs/openms.yaml"
    shell:
        """
        GNPSExport -in_cm {input.var1} -in_mzml {input.var2} -out {output.out1} -out_quantification {output.out2} -out_pairs {output.out3} -out_meta_values {output.out4} -log {log} 2>> {log}
        """