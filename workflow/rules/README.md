# Snakemake rules 

The pipeline consists of three separate workflows that are interconnected, and one data analysis guide:

### `1) File conversion:`

Convert raw files from Thermo to open community-driven format mzML centroid (see documentation [here](https://github.com/compomics/ThermoRawFileParser))

### `2) Pre-processing:`

Converting raw data to a feature table with a series of OpenMS algorithms (see documentation [here](https://abibuilder.informatik.uni-tuebingen.de/archive/openms/Documentation/nightly/html/index.html)). Additinaly, the feature table is introduced to SIRIUS and CSI:FingerID for formula and structural predictions (see documentation [here](https://boecker-lab.github.io/docs.sirius.github.io/)).
CSI:FingerID is using external Web servers (from the Boecher lab in Jena) for the structural library seach and all computations for the structural predictions. The disadvantage in this case is that the workflow is dependent on the functionality of their servers, queued jobs, etc. 

CSI_FingeID can be easily removed from the workflow by replacing rule sirius and df_ sirius with the following script:

```
# 4) SIRIUS 

rule sirius:
    input: 
        var1= "resources/Sirius/sirius.app/Contents/MacOS/sirius",
        var2= "results/{samples}/interim/sirius/PCfeature_nch_{samples}.mzML", 
        var3= "results/{samples}/interim/sirius/MFD_nch_{samples}.featureXML"        
    output:
        "results/{samples}/interim/sirius/formulas_{samples}.mzTab"
    threads: 2
    shell:
        """
        resources/OpenMS-2.7.0/bin/SiriusAdapter -executable {input.var1} -in {input.var2} -in_featureinfo {input.var3} -out_sirius {output[0]} -preprocessing:filter_by_num_masstraces 2 -preprocessing:feature_only -sirius:profile orbitrap -sirius:db none -sirius:ions_considered "[M+H]+, [M-H2O+H]+, [M+Na]+, [M+NH4]+" -sirius:elements_enforced CHN[15]O[40]S[4] -debug 3 -project:processors {threads} -threads {threads}
        """
# 5) Convert the mzTab to a csv file

rule df_sirius:
    input: 
        "results/{samples}/interim/sirius/formulas_{samples}.mzTab"
    output:
        "results/{samples}/formulas_{samples}.csv"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/df_sirius.py"

```

and removing the output files also from [Snakefile](workflow/Snakefile) by replacing lines 14 and 15 with the following:
```
        expand("results/{samples}/interim/sirius/formulas_{samples}.mzTab", samples=SAMPLES),
        expand("results/{samples}/formulas_{samples}.csv", samples=SAMPLES),
```



![dag](/images/Preprocessing+SIRIUS_CSI_FingerID.svg)

### `3) GNPSexport:` 

Generate all the files necessary to create a FBMN job at GNPS (see documentation [here](https://ccms-ucsd.github.io/GNPSDocumentation/featurebasedmolecularnetworking-with-openms/)). 

![dag](/images/GNPSExport.svg) 

### `4) Re-quantification:` 

Re-quantify all raw files to avoid missing values resulted by the pre-processing workflow for statistical analysis and data exploration. Generate a FeatureMatrix for further statistical analysis.

![dag](/images/Re-quantification.svg) 
