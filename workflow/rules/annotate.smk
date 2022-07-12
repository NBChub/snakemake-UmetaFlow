import os
import glob
import os.path 

# 1) Create a sirius library from all the tables with formula predictions by only taking into acount the rank #1 predictions for simplicity. Mind that there are cases where SIRIUS predicts the correct formula ranked as >1. 

if config["rules"]["sirius_csi"]==True:
    if config["rules"]["requantification"]==True:
        rule sirius_annotations:
            input:
                matrix= "results/Requantified/FeatureMatrix.tsv",
                sirius= expand("results/SiriusCSI/formulas_{samples}.tsv", samples=SAMPLES),
                csi= expand("results/SiriusCSI/structures_{samples}.tsv", samples=SAMPLES)
            output:
                annotated= "results/annotations/annotated_FeatureTable.tsv"
            log: "workflow/report/logs/annotate/sirius_annotations.log"
            threads: 4
            conda:
                "../envs/openms.yaml"
            shell:
                """
                python workflow/scripts/SIRIUS_CSI_annotations.py {input.matrix} {output.annotated} 2>> {log}
                """
    else:
        rule sirius_annotations:
            input:
                matrix="results/Preprocessed/FeatureMatrix.tsv",
                sirius= expand("results/SiriusCSI/formulas_{samples}.tsv", samples=SAMPLES),
                csi= expand("results/SiriusCSI/structures_{samples}.tsv", samples=SAMPLES)
            output:
                annotated= "results/annotations/annotated_FeatureTable.tsv"
            log: "workflow/report/logs/annotate/sirius_annotations.log"
            threads: 4
            conda:
                "../envs/openms.yaml"
            shell:
                """
                python workflow/scripts/SIRIUS_CSI_annotations.py {input.matrix} {output.annotated} 2>> {log}
                """ 
else:
    if config["rules"]["requantification"]==True:
        rule sirius_annotations:
            input:
                matrix= "results/Requantified/FeatureMatrix.tsv",
                sirius= expand("results/Sirius/formulas_{samples}.tsv", samples=SAMPLES)
            output:
                annotated= "results/annotations/annotated_FeatureTable.tsv"
            log: "workflow/report/logs/annotate/sirius_annotations.log"
            threads: 4
            conda:
                "../envs/openms.yaml"
            shell:
                """
                python workflow/scripts/SIRIUS_annotations.py {input.matrix} {output.annotated} 2>> {log}    
                """
    else:
        rule sirius_annotations:
            input:
                matrix= "results/Preprocessed/FeatureMatrix.tsv",
                sirius= expand("results/Sirius/formulas_{samples}.tsv", samples=SAMPLES)
            output:
                annotated= "results/annotations/annotated_FeatureTable.tsv"
            log: "workflow/report/logs/annotate/sirius_annotations.log"
            threads: 4
            conda:
                "../envs/openms.yaml"
            shell:
                """
                python workflow/scripts/SIRIUS_annotations.py {input.matrix} {output.annotated} 2>> {log}    
                """

# 2) Run your aligned mzml files (from directory GNPSexport/mzML) directly to GNPS for MS/MS library matching to generate a tsv table of metabolites. 
# Filter out the ones that have a mass error > 10.0 ppm and also metabolites that originate from libraries such as HMDB when your samples are generated from bacteria.
# Annotate compounds in FeatureMatrix. The annotation can be considered a bit ambiguous because it is done with mz and RT information.

GNPS_library = find_files("resources", "*.tsv")
if GNPS_library:
    rule GNPS_annotations:
        input:
            lib= glob.glob(os.path.join("resources", "*.tsv")),
            featurematrix= "results/annotations/annotated_FeatureTable.tsv"
        output:
            gnps= "results/annotations/GNPS_annotated_FeatureTable.tsv"
        log: "workflow/report/logs/annotate/GNPS_annotations.log"
        threads: 4
        conda:
            "../envs/openms.yaml"
        shell:
            """
            python workflow/scripts/GNPS.py {input.lib} {input.featurematrix} {output.gnps} 2>> {log}
            """
else:
    print("no file found")
    rule GNPS_annotations:
        input:
            "results/annotations/SIRIUS_CSI_annotated_FeatureTable.tsv"
        output:
            "results/annotations/GNPS_annotated_FeatureTable.tsv"
        log: "workflow/report/logs/annotate/GNPS_annotations.log"
        conda:
            "../envs/openms.yaml"
        shell:
            """ 
            echo "No GNPS metabolite identification file was found" > {output} 2>> {log}
            """