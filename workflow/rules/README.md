# Snakemake rules 

The pipeline consists of three separate workflows that are interconnected, and one data analysis guide:

### `1) File conversion:`

Conversion of raw files from Thermo to open community-driven format mzML centroid (see documentation [here](https://github.com/compomics/ThermoRawFileParser)).

If you have Agilent or Bruker files, skip that step (write "FALSE" for rule fileconversion in the [config.yaml](config/config.yaml) file, convert the files independently using proteowizard (see https://proteowizard.sourceforge.io/) and add them to the data/mzML/ directory.

### `2) Pre-processing:`

Converting raw data to a feature table with a series of OpenMS algorithms (see documentation [here](https://abibuilder.informatik.uni-tuebingen.de/archive/openms/Documentation/nightly/html/index.html)). 

![dag](/images/Preprocessing.svg) 

### `3) Re-quantification:` 

Re-quantify all raw files to avoid missing values resulted by the pre-processing steps for statistical analysis and data exploration. Generate a FeatureMatrix for further statistical analysis.

![dag](/images/Re-quantification.svg) 

### `4) SIRIUS and CSI:FingerID:`

The pre-processed feature tables are then introduced to SIRIUS and CSI:FingerID for formula and structural predictions (see documentation [here](https://boecker-lab.github.io/docs.sirius.github.io/)).

CSI:FingerID is using external Web servers (from the Boecher lab in Jena) for the structural library seach and all computations for the structural predictions. The disadvantage in this case is that the workflow is dependent on the functionality of their servers, queued jobs, etc. 

CSI_FingeID can be easily removed from the workflow by replacing #4) rule sirius and #5) df_ sirius from [sirius.smk](workflow/rules/sirius.smk) with the following script:

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

and also removing the output files also from [Snakefile](workflow/Snakefile) by replacing lines 19 and 20 with the following:

```
        expand("results/{samples}/interim/sirius/formulas_{samples}.mzTab", samples=SAMPLES),
        expand("results/{samples}/formulas_{samples}.csv", samples=SAMPLES),
```

![dag](/images/SIRIUS_CSI_FingerID.svg)

### `5) GNPSexport:` 

Generate all the files necessary to create a FBMN job at GNPS (see documentation [here](https://ccms-ucsd.github.io/GNPSDocumentation/featurebasedmolecularnetworking-with-openms/)). 

Create a metadata csv file for GNPS from the samples.tsv file using the Jupyter notebook [file](Create_sampletsv_file.ipynb)

![dag](/images/GNPSExport.svg) 

