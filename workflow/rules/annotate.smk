# 1) Create a sirius library from all the tables with formula predictions by only taking into acount the rank #1 predictions for simplicity. Mind that there are cases where SIRIUS predicts the correct formula ranked as >1. 

rule sirius_library:
    input:
        glob.glob(os.path.join("results", "SIRIUS", "formulas_*.tsv"))
    output:
        "results/SIRIUS/SIRIUS_library.tsv"
    threads: 4
    conda:
        "../envs/openms.yaml"
    script:
        "../scripts/SIRIUS_library.py"    
        
# 2) Create a CSI library from all the tables with formula predictions by only taking into acount the rank #1 predictions for simplicity. Mind that there are cases where SIRIUS:CSI predicts the correct formula ranked as >1. 

rule csi_library:
    input:
        glob.glob(os.path.join("results", "CSI", "structures_*.tsv"))
    output:
        "results/CSI/CSI_library.tsv"
    threads: 4
    conda:
        "../envs/openms.yaml"
    script:
        "../scripts/CSI_library.py"    

# 3) Annotate the Feature Matrix with structural and formula predictions. I recomment that you use both SIRIUS and CSI , even if there are repeating formula values occassionaly, since CSI would sometimes not predict a formula (due to size?)       

rule SIRIUS_annotations:
    input:
        "results/SIRIUS/SIRIUS_library.tsv",
        "results/CSI/CSI_library.tsv",
        "results/Requantified/FeatureMatrix.tsv"
    output:
        "results/annotations/SIRIUS_CSI_annotated_FeatureTable.tsv"
    threads: 4
    conda:
        "../envs/openms.yaml"
    script:
        "../scripts/SIRIUS_CSI_annotations.py"    


# 4) Run your raw mzml files directly to GNPS for MS/MS library matching to generate a tsv table of metabolites. 
# Filter out the ones that have a mass error > 20.0 ppm and also metabolites that originate from libraries such as HMDB when your samples are generated from bacteria.
# Annotate compounds in FeatureMatrix

rule GNPS_library:
    input:
        "resources/MS2_LIBRARYSEARCH_all_identifications.tsv",
        "results/annotations/SIRIUS_CSI_annotated_FeatureTable.tsv"
    output:
        "results/annotations/GNPS_annotated_FeatureTable.tsv"
    threads: 4
    conda:
        "../envs/openms.yaml"
    script:
        "../scripts/GNPS.py"       