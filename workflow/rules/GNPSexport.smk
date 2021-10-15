#Create a metadata csv file for GNPS from the samples.tsv file 
import pandas as pd
import numpy as np 
import os

path= "results/GNPSexport/"
isExist= os.path.exists(path)
if not isExist:
    os.mkdir("results/GNPSexport/")
df= pd.read_csv("config/samples.tsv", sep= "\t", header= 0)
metadata= df.rename(columns= {"sample_name": "filename", "comment": "ATTRIBUTE_comment", "MAPnumber": "ATTRIBUTE_MAPnumber"})
metadata["filename"]= metadata["filename"].astype(str) +".mzml"
metadata['ATTRIBUTE_MAPnumber'] = np.arange(len(metadata))
metadata["ATTRIBUTE_comment"]= metadata["ATTRIBUTE_comment"].astype(str) + "_MAP" + metadata["ATTRIBUTE_MAPnumber"].astype(str)
metadata= metadata.drop(columns= "ATTRIBUTE_MAPnumber")
metadata['ATTRIBUTE_genomeID']=metadata['filename'].str.extract(r'(NBC_?\w{5})')
metadata['ATTRIBUTE_genomeIDMDNA']=metadata['filename'].str.extract(r'(MDNA_?\w{5})')
metadata['ATTRIBUTE_genomeID']=metadata['ATTRIBUTE_genomeID'].fillna(metadata['ATTRIBUTE_genomeIDMDNA'])
metadata=metadata.drop(columns="ATTRIBUTE_genomeIDMDNA")
#metadata['ATTRIBUTE_genomeID']= metadata['ATTRIBUTE_genomeID'].replace(to_replace= r'NBC', value= 'NBC_', regex= True)
#metadata['ATTRIBUTE_genomeID']= metadata['ATTRIBUTE_genomeID'].replace(to_replace= r'MDNAWGS', value= 'MDNA_WGS_', regex= True)
metadata.to_csv("results/GNPSexport/metadata.tsv", sep='\t')

#copy all the original mzml files (precursor corrected ones) in the GNPSExport folder for easier used
rule FileCopy:
    input:
        "results/{samples}/interim/precursorcorrected_{samples}.mzML"
    output:
        "results/GNPSexport/{samples}.mzML"
    shell:
        """
        cp {input} {output}
        """ 

#MapAlignerPoseClustering is used to perform a linear retention time alignment, basically correct for linear shifts in retention time.
#add trafoXML files as an output also (TransformationXMLFile()) for transforming also the MS2 spectra later on

rule MapAlignerPoseClustering:
    input:
        expand("results/{samples}/interim/MFD_{samples}.featureXML", samples=SAMPLES)
    output:
        var1= expand("results/Consensus/interim/MapAlignerPoseClustering_{samples}.featureXML", samples=SAMPLES),
        var2= expand("results/Consensus/interim/MapAlignerPoseClustering_{samples}.trafoXML", samples=SAMPLES)
    shell:
        """
        resources/OpenMS-2.7.0/bin/MapAlignerPoseClustering -in {input} -out {output.var1} -trafo_out {output.var2}
        """ 

#Introduce the features to a protein identification file (idXML)- the only way to create a ConsensusXML file currently (run FeatureLinkerUnlabeledKD)       
rule IDMapper:
    input:
        "resources/emptyfile.idXML",
        "results/Consensus/interim/MapAlignerPoseClustering_{samples}.featureXML",
        "results/{samples}/interim/precursorcorrected_{samples}.mzML"
    output:
        "results/Consensus/interim/IDMapper_{samples}.featureXML"
    shell:
        """
        resources/OpenMS-2.7.0/bin/IDMapper -id {input[0]} -in {input[1]}  -spectra:in {input[2]} -out {output} 
        """

#The FeatureLinkerUnlabeledKD is used to aggregate the feature information (from single files) into a ConsensusFeature, linking features from different files together, which have a smiliar m/z and rt (no MS2 data).
rule FeatureLinkerUnlabeledKD:
    input:
        expand("results/Consensus/interim/IDMapper_{samples}.featureXML", samples=SAMPLES)
    output:
        "results/Consensus/interim/FeatureLinkerUnlabeledKD.consensusXML"
    shell:
        """
        resources/OpenMS-2.7.0/bin/FeatureLinkerUnlabeledKD -in {input} -out {output} 
        """

#Filter out the features that do not have an MS2 pattern
rule FileFilter:
    input:
        "results/Consensus/interim/FeatureLinkerUnlabeledKD.consensusXML"
    output:
        "results/Consensus/filtered.consensusXML"
    shell:
        """
        resources/OpenMS-2.7.0/bin/FileFilter -id:remove_unannotated_features -in {input} -out {output} 
        """

#GNPS_export creates an mgf file with only the MS2 information of all files (introduce mzml files with spaces between them)
rule GNPS_export:
    input:
        var1= "results/Consensus/filtered.consensusXML",
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
        "results/Consensus/filtered.consensusXML"
    output:
        "results/GNPSexport/FeatureQuantificationTable.txt" 
    shell:
        """
        resources/OpenMS-2.7.0/bin/TextExporter -in {input} -out {output}
        """
