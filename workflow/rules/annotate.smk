# 1) Create a sirius library from all the tables with formula predictions by only taking into acount the rank #1 predictions for simplicity. Mind that there are cases where SIRIUS predicts the correct formula ranked as >1. 

rule sirius_library:
    input:
        expand("results/SIRIUS/formulas_{samples}.csv", samples=SAMPLES)
    output:
        "results/SIRIUS/SIRIUS_library.csv"
    threads: 4
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/SIRIUS_library.py"    
        
# 2) Create a CSI library from all the tables with formula predictions by only taking into acount the rank #1 predictions for simplicity. Mind that there are cases where SIRIUS:CSI predicts the correct formula ranked as >1. 

rule csi_library:
    input:
        expand("results/CSI/structures_{samples}.csv", samples=SAMPLES)
    output:
        "results/CSI/CSI_library.csv"
    threads: 4
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/CSI_library.py"    

# 3) Annotate the Feature Matrix with structural and formula predictions. I recomment that you use both SIRIUS and CSI , even if there are repeating formula values occassionaly, since CSI would sometimes not predict a formula (due to size?)       
    # You can chose between the preprocessed or the final requantified matrix

rule SIRIUS_annotations:
    input:
        "results/SIRIUS/SIRIUS_library.csv",
        "results/CSI/CSI_library.csv",
        "results/Preprocessed/FeatureQuantificationTable.csv"
    output:
        "results/annotations/SIRIUS_CSI_annotated_FeatureTable.csv"
    threads: 4
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/SIRIUS_CSI_annotations.py"    

# 4) Annotate the Requantified Feature Matrix with structural and formula predictions. I recomment that you use both SIRIUS and CSI , even if there are repeating formula values occassionaly, since CSI would sometimes not predict a formula (due to size?)       
    # You can chose between the preprocessed or the final requantified matrix

rule SIRIUS_annotations_requant:
    input:
        "results/SIRIUS/SIRIUS_library.csv",
        "results/CSI/CSI_library.csv",
        "results/Requantified/FeatureMatrix_requant.tsv"
    output:
        "results/annotations/SIRIUS_CSI_annotated_FeatureTable_Requant.csv"
    threads: 4
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/SIRIUS_CSI_annotations.py"    


# 5) Run your raw mzml files directly to GNPS for MS/MS library matching to generate a csv table of metabolites. 
# Filter out the ones that have a mass error > 20.0 ppm and also metabolites that originate from libraries such as HMDB when your samples are generated from bacteria.
# Annotate compounds in FeatureMatrix

rule GNPS_library:
    input:
        "resources/MS2_LIBRARYSEARCH_all_identifications.tsv",
        "results/Preprocessed/FeatureQuantificationTable.csv"
    output:
        "results/annotations/GNPS_annotated_FeatureTable.csv"
    threads: 4
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/GNPS.py"     

# 6) Run your raw mzml files directly to GNPS for MS/MS library matching to generate a csv table of metabolites. 
# Filter out the ones that have a mass error > 20.0 ppm and also metabolites that originate from libraries such as HMDB when your samples are generated from bacteria.
# Annotate compounds in FeatureMatrix

#rule GNPS_library_requant:
#    input:
#        "resources/MS2_LIBRARYSEARCH_all_identifications.tsv",
#        "results/Requantified/FeatureMatrix_requant.tsv"
#    output:
#        "results/annotations/GNPS_annotated_FeatureTable_Requant.csv"
#    threads: 4
#    conda:
#        "../envs/python.yaml"
#    script:
#        "../scripts/GNPS.py"     