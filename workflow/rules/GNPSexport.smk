# Create a directory with all files necessary for FBMN

# 1) Create a metadata csv file for GNPS from the samples.tsv file 

import pandas as pd
import numpy as np 
path= "results/GNPSexport/interim"
isExist= os.path.exists(path)
if not isExist:
    os.mkdir(path)

df= pd.read_csv("config/samples.tsv", sep= "\t", index_col= "Unnamed: 0")
df["sample_name"]=df["sample_name"].replace(to_replace= r'MDNAWGS', value= 'MDNA_WGS_', regex= True)
metadata= df.rename(columns= {"sample_name": "filename", "comment": "ATTRIBUTE_comment", "MAPnumber": "ATTRIBUTE_MAPnumber"})
metadata["filename"]= metadata["filename"].astype(str) +".mzml"
metadata['ATTRIBUTE_MAPnumber'] = np.arange(len(metadata))
metadata["ATTRIBUTE_MAP_ID"]= "MAP" + metadata["ATTRIBUTE_MAPnumber"].astype(str)
metadata= metadata.drop(columns= "ATTRIBUTE_MAPnumber")
metadata['ATTRIBUTE_genomeID']=metadata['filename'].str.extract(r'(NBC_?\d*)')
metadata['ATTRIBUTE_genomeIDMDNA']=metadata['filename'].str.extract(r'(MDNAWGS?\d*|MDNA_WGS_?\d*)')
metadata['ATTRIBUTE_genomeID']=metadata['ATTRIBUTE_genomeID'].fillna(metadata['ATTRIBUTE_genomeIDMDNA'])
metadata["ATTRIBUTE_media"]= metadata['filename'].str.extract(r'(ISP2|DNPM|FPY12\d*)')
metadata["ATTRIBUTE_comment"]= metadata['ATTRIBUTE_genomeID'].astype(str) +"_" + metadata["ATTRIBUTE_media"].astype(str)
metadata=metadata.drop(columns="ATTRIBUTE_genomeIDMDNA")
#metadata['ATTRIBUTE_genomeID']= metadata['ATTRIBUTE_genomeID'].replace(to_replace= r'NBC', value= 'NBC_', regex= True)
#metadata['ATTRIBUTE_genomeID']= metadata['ATTRIBUTE_genomeID'].replace(to_replace= r'MDNAWGS', value= 'MDNA_WGS_', regex= True)
metadata.to_csv("results/GNPSexport/metadata.tsv", sep='\t')


# 2) copy all the original mzml files (precursor-corrected ones) in the GNPSExport folder for easier use
rule FileCopy:
    input:
        "results/{samples}/interim/PCfeature_{samples}.mzML"
    output:
        "results/GNPSexport/{samples}.mzML"
    shell:
        """
        cp {input} {output}
        """ 

# 3) MapAlignerPoseClustering is used to perform a linear retention time alignment, to correct for linear shifts in retention time between different runs.

rule MapAlignerPoseClustering:
    input:
        expand("results/{samples}/interim/preprocessed/MFD_{samples}.featureXML", samples=SAMPLES)
    output:
        var1= expand("results/GNPSexport/interim/MapAlignerPoseClustering_{samples}.featureXML", samples=SAMPLES),
        var2= expand("results/GNPSexport/interim/MapAlignerPoseClustering_{samples}.trafoXML", samples=SAMPLES)
    shell:
        """
        OpenMS/OpenMS-build/bin/MapAlignerPoseClustering -algorithm:max_num_peaks_considered -1 -algorithm:superimposer:mz_pair_max_distance 0.05 -algorithm:pairfinder:distance_MZ:max_difference 10.0 -algorithm:pairfinder:distance_MZ:unit ppm -in {input} -out {output.var1} -trafo_out {output.var2}
        """ 

# 4) Introduce the features to a protein identification file (idXML)- the only way to create an aggregated ConsensusXML file currently (run FeatureLinkerUnlabeledKD)  

rule IDMapper:
    input:
        var1= "resources/emptyfile.idXML",
        var2= "results/GNPSexport/interim/MapAlignerPoseClustering_{samples}.featureXML",
        var3= "results/{samples}/interim/PCfeature_{samples}.mzML"
    output:
        "results/GNPSexport/interim/IDMapper_{samples}.featureXML"
    shell:
        """
        OpenMS/OpenMS-build/bin/IDMapper -id {input.var1} -in {input.var2}  -spectra:in {input.var3} -out {output} 
        """

# 5) The FeatureLinkerUnlabeledKD is used to aggregate the feature information (from single files) into a ConsensusFeature, linking features from different files together, which have a smiliar m/z and rt (MS1 level).

rule FeatureLinkerUnlabeledKD:
    input:
        expand("results/GNPSexport/interim/IDMapper_{samples}.featureXML", samples=SAMPLES)
    output:
        "results/GNPSexport/interim/FeatureLinkerUnlabeledKD.consensusXML"
    shell:
        """
        OpenMS/OpenMS-build/bin/FeatureLinkerUnlabeledKD -in {input} -out {output} 
        """

# 6) export the consensusXML file to a csv file for FFMI (later)

rule csv_export:
    input:
        "results/GNPSexport/interim/FeatureLinkerUnlabeledKD.consensusXML"
    output:
        "results/GNPSexport/interim/consensus.tsv" 
    shell:
        """
        OpenMS/OpenMS-build/bin/TextExporter -in {input} -out {output}
        """

# 7) Filter out the features that do not have an MS2 pattern (no protein ID annotations)

rule FileFilter:
    input:
        "results/GNPSexport/interim/FeatureLinkerUnlabeledKD.consensusXML"
    output:
        "results/GNPSexport/interim/filtered.consensusXML"
    shell:
        """
        OpenMS/OpenMS-build/bin/FileFilter -id:remove_unannotated_features -in {input} -out {output} 
        """

# 8) GNPS_export creates an mgf file with only the MS2 information of all files (introduce mzml files with spaces between them)

rule GNPS_export:
    input:
        var1= "results/GNPSexport/interim/filtered.consensusXML",
        var2= expand("results/{samples}/interim/PCfeature_{samples}.mzML", samples=SAMPLES)
    output:
        "results/GNPSexport/MSMS.mgf" 
    shell:
        """
        OpenMS/OpenMS-build/bin/GNPSExport -ini resources/GNPSExport.ini -in_cm {input.var1} -in_mzml {input.var2} -out {output} 
        """

# 9) export the consensusXML file to a txt file for GNPS

rule txt_export:
    input:
        "results/GNPSexport/interim/filtered.consensusXML"
    output:
        "results/GNPSexport/FeatureQuantificationTable.txt" 
    shell:
        """
        OpenMS/OpenMS-build/bin/TextExporter -in {input} -out {output}
        """
