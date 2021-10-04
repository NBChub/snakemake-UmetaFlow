#Create a metadata csv file for GNPS from the samples.tsv file 
import pandas
df= pandas.read_csv("config/samples.tsv", sep= "\t", header= 0)
metadata= df.rename(columns= {"sample_name": "filename", "comment": "ATTRIBUTE_comment"})
metadata["filename"]= "precursorcorrected_" + metadata["filename"].astype(str) +".mzml"
metadata.to_csv("results/GNPSexport/metadata.tsv", sep='\t')

#MapAlignerPoseClustering is used to perform a linear retention time alignment, basically correct for linear shifts in retention time.
rule MapAlignerPoseClustering:
    input:
        expand("results/{samples}/interim/MFD_{samples}.featureXML", samples=SAMPLES)
    output:
        expand("results/GNPSexport/interim/MapAlignerPoseClustering_{samples}.featureXML", samples=SAMPLES)
    shell:
        """
        resources/OpenMS-2.7.0/bin/MapAlignerPoseClustering -in {input} -out {output}
        """ 

#Introduce the features to a protein identification file (idXML)       
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

#The FeatureLinkerUnlabeledKD is used to aggregate the feature information (from single files) into a ConsensusFeature, linking features from different files together, which have a smiliar m/z and rt (no MS2 data).
rule FeatureLinkerUnlabeledKD:
    input:
        expand("results/GNPSexport/interim/IDMapper_{samples}.featureXML", samples=SAMPLES)
    output:
        "results/GNPSexport/interim/FeatureLinkerUnlabeledKD.consensusXML"
    shell:
        """
        resources/OpenMS-2.7.0/bin/FeatureLinkerUnlabeledKD -in {input} -out {output} 
        """

#Filter out the features that do not have an MS2 pattern
rule FileFilter:
    input:
        "results/GNPSexport/interim/FeatureLinkerUnlabeledKD.consensusXML"
    output:
        "results/GNPSexport/interim/filtered.consensusXML"
    shell:
        """
        resources/OpenMS-2.7.0/bin/FileFilter -id:remove_unannotated_features -in {input} -out {output} 
        """

#GNPS_export creates an mgf file with only the MS2 information of all files (introduce mzml files with spaces between them)
rule GNPS_export:
    input:
        var1= "results/GNPSexport/interim/filtered.consensusXML",
        var2= expand("results/{samples}/interim/precursorcorrected_{samples}.mzML", samples=SAMPLES)
    output:
        "results/GNPSexport/MSMS.mgf" 
    shell:
        """
        resources/OpenMS-2.7.0/bin/GNPSExport -in_cm {input.var1} -in_mzml {input.var2} -out {output} 
        """

#export the consensusXML file to a txt file for GNPS
rule txt_export:
    input:
        "results/GNPSexport/interim/filtered.consensusXML"
    output:
        "results/GNPSexport/FeatureQuantificationTable.txt" 
    shell:
        """
        resources/OpenMS-2.7.0/bin/TextExporter -in {input} -out {output}
        """
