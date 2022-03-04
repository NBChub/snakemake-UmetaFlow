# Create a directory with all files necessary for FBMN

# 1) Create a metadata csv file for GNPS from the samples.tsv file 

import pandas as pd
import numpy as np 
path= "results/GNPSexport/"
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
        "results/Interim/mzML/PCfeature_{samples}.mzML"
    output:
        "results/GNPSexport/{samples}.mzML"
    shell:
        """
        cp {input} {output}
        """ 


# 3) Filter out the features that do not have an MS2 pattern (no protein ID annotations)

rule FileFilter:
    input:
        "results/Interim/preprocessed/FeatureLinkerUnlabeledKD.consensusXML"
    output:
        "results/Interim/GNPSexport/filtered.consensusXML"
    shell:
        """
        FileFilter -id:remove_unannotated_features -in {input} -out {output} 
        """

# 4) GNPS_export creates an mgf file with only the MS2 information of all files (introduce mzml files with spaces between them)

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

# 5) export the consensusXML file to a txt file for GNPS

rule GNPS_txt_export:
    input:
        "results/Interim/GNPSexport/filtered.consensusXML"
    output:
        "results/GNPSexport/FeatureQuantificationTable.txt" 
    shell:
        """
        TextExporter -in {input} -out {output}
        """
