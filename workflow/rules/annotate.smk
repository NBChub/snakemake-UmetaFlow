import os
import glob
import os.path 

# 1) Create a sirius library from all the tables with formula predictions by only taking into acount the rank #1 predictions for simplicity. Mind that there are cases where SIRIUS predicts the correct formula ranked as >1. 

if config["rules"]["sirius_csi"]==True:
    if config["rules"]["requantification"]==True:
        rule sirius_annotations:
            input:
                "results/Requantified/FeatureMatrix.tsv",
<<<<<<< HEAD
                expand("results/SiriusCSI/formulas_{samples}.tsv", samples=SAMPLES),
                expand("results/SiriusCSI/structures_{samples}.tsv", samples=SAMPLES)
            output:
                "results/annotations/annotated_FeatureTable.tsv"
=======
                expand("results/SIRIUS/formulas_{samples}.tsv", samples=SAMPLES),
                expand("results/CSI/structures_{samples}.tsv", samples=SAMPLES)
            output:
                "results/annotations/SIRIUS_CSI_annotated_FeatureTable.tsv"
>>>>>>> 8c5ccedc94fc856850dde76e921e22a705428575
            threads: 4
            conda:
                "../envs/openms.yaml"
            script:
                "../scripts/SIRIUS_CSI_annotations.py"    
    else:
        rule sirius_annotations:
            input:
                "results/Preprocessed/FeatureMatrix.tsv",
<<<<<<< HEAD
                expand("results/SiriusCSI/formulas_{samples}.tsv", samples=SAMPLES),
                expand("results/SiriusCSI/structures_{samples}.tsv", samples=SAMPLES)
            output:
                "results/annotations/annotated_FeatureTable.tsv"
=======
                expand("results/SIRIUS/formulas_{samples}.tsv", samples=SAMPLES),
                expand("results/CSI/structures_{samples}.tsv", samples=SAMPLES)
            output:
                "results/annotations/SIRIUS_CSI_annotated_FeatureTable.tsv"
>>>>>>> 8c5ccedc94fc856850dde76e921e22a705428575
            threads: 4
            conda:
                "../envs/openms.yaml"
            script:
                "../scripts/SIRIUS_CSI_annotations.py"   
else:
    if config["rules"]["requantification"]==True:
        rule sirius_annotations:
            input:
                "results/Requantified/FeatureMatrix.tsv",
<<<<<<< HEAD
                expand("results/Sirius/formulas_{samples}.tsv", samples=SAMPLES)
            output:
                "results/annotations/annotated_FeatureTable.tsv"
=======
                expand("results/SIRIUS/formulas_{samples}.tsv", samples=SAMPLES)
            output:
                "results/annotations/SIRIUS_annotated_FeatureTable.tsv"
>>>>>>> 8c5ccedc94fc856850dde76e921e22a705428575
            threads: 4
            conda:
                "../envs/openms.yaml"
            script:
                "../scripts/SIRIUS_annotations.py"    
    else:
        rule sirius_annotations:
            input:
                "results/Preprocessed/FeatureMatrix.tsv",
<<<<<<< HEAD
                expand("results/Sirius/formulas_{samples}.tsv", samples=SAMPLES)
            output:
                "results/annotations/annotated_FeatureTable.tsv"
=======
                expand("results/SIRIUS/formulas_{samples}.tsv", samples=SAMPLES)
            output:
                "results/annotations/SIRIUS_annotated_FeatureTable.tsv"
>>>>>>> 8c5ccedc94fc856850dde76e921e22a705428575
            threads: 4
            conda:
                "../envs/openms.yaml"
            script:
                "../scripts/SIRIUS_annotations.py" 

# 2) Run your aligned mzml files (from directory GNPSexport/mzML) directly to GNPS for MS/MS library matching to generate a tsv table of metabolites. 
# Filter out the ones that have a mass error > 10.0 ppm and also metabolites that originate from libraries such as HMDB when your samples are generated from bacteria.
# Annotate compounds in FeatureMatrix. The annotation can be considered a bit ambiguous because it is done with mz and RT information.

GNPS_library = find_files("resources", "*.tsv")
if GNPS_library:
    rule GNPS_annotations:
        input:
            glob.glob(os.path.join("resources", "*.tsv")),
<<<<<<< HEAD
            "results/annotations/annotated_FeatureTable.tsv"
=======
            "results/annotations/SIRIUS_CSI_annotated_FeatureTable.tsv"
>>>>>>> 8c5ccedc94fc856850dde76e921e22a705428575
        output:
            "results/annotations/GNPS_annotated_FeatureTable.tsv"
        threads: 4
        conda:
            "../envs/openms.yaml"
        script:
            "../scripts/GNPS.py" 
else:
    print("no file found")
    rule GNPS_annotations:
        input:
            "results/annotations/SIRIUS_CSI_annotated_FeatureTable.tsv"
        output:
            "results/annotations/GNPS_annotated_FeatureTable.tsv"
        shell:
            """ 
            echo "No GNPS metabolite identification file was found" > {output}
            """