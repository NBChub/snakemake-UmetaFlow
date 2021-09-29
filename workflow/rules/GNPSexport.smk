rule MapAlignerPoseClustering:
    input:
        lambda w: expand("results/{samples}/interim/MFD_{samples}.featureXML", samples=SAMPLES)
    output:
        expand("results/GNPSexport/interim/MapAlignerPoseClustering_{samples}.featureXML", samples=SAMPLES)
    shell:
        """
        resources/OpenMS-2.7.0/bin/MapAlignerPoseClustering -in {input} -out {output} 
        """ 
       
rule IDMapper:
    input:
        "resources/emptyfile.idXML",
        "results/GNPSexport/interim/MapAlignerPoseClustering_{samples}.featureXML",
        "results/{samples}/interim/precursorcorrected_{samples}.mzML"
    output:
        "results/GNPSexport/interim/IDMapper_{samples}.featureXML"
    shell:
        """
        resources/OpenMS-2.7.0/bin/IDMapper -id {input[0]} -in {input[1]}  -spectra:in {input[2]} -out {output} 
        """

rule FeatureLinkerUnlabeledKD:
    input:
        lambda w: expand("results/GNPSexport/interim/IDMapper_{samples}.featureXML", samples=SAMPLES)
    output:
        "results/GNPSexport/interim/FeatureLinkerUnlabeledKD.consensusXML"
    shell:
        """
        resources/OpenMS-2.7.0/bin/FeatureLinkerUnlabeledKD -in {input} -out {output} 
        """

rule FileFilter:
    input:
        "results/GNPSexport/interim/FeatureLinkerUnlabeledKD.consensusXML"
    output:
        "results/GNPSexport/interim/filtered.consensusXML"
    shell:
        """
        resources/OpenMS-2.7.0/bin/FileFilter -id:remove_unannotated_features -in {input} -out {output} 
        """

rule GNPS_export:
    input:
        "results/GNPSexport/interim/filtered.consensusXML",
        lambda w: expand("results/{samples}/interim/precursorcorrected_{samples}.mzML", samples=SAMPLES),
    output:
        "results/GNPSexport/MSMS.mgf" 
    shell:
        """
        resources/OpenMS-2.7.0/bin/GNPSExport -in_cm {input[0]} -in_mzml {input[1]} -out {output} 
        """


rule txt_export:
    input:
        "results/GNPSexport/interim/filtered.consensusXML"
    output:
        "results/GNPSexport/FeatureQuantificationTable.txt" 
    shell:
        """
        resources/OpenMS-2.7.0/bin/TextExporter -in {input} -out {output}
        """