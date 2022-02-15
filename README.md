# Metabolomics workflow for Linux systems

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.7.0-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Build Status](https://travis-ci.org/snakemake-workflows/snakemake-bgc-analytics.svg?branch=master)](https://travis-ci.org/snakemake-workflows/snakemake-bgc-analytics)

This is a snakemake implementation of the [Metabolomics OpenMS workflow](./OpenMS_workflow.ipynb) tailored by [Eftychia Eva Kontou](https://github.com/eeko-kon)

## Workflow overview

The pipeline consists of three separate workflows that are interconnected, and one data analysis guide:

1) Pre-processing: converting raw data to a feature table with a series of algorithms 

2) GNPSexport: generate all the files necessary to create a FBMN job at GNPS. (see https://ccms-ucsd.github.io/GNPSDocumentation/featurebasedmolecularnetworking-with-openms/) 

3) Re-quantification: Re-quantify all raw files to avoid missing values resulted by the pre-processing workflow for statistical analysis and data exploration.

![dag](/images/MetabolomicsFlow.svg)

## Usage

### Step 1: Clone the workflow

[Clone](https://help.github.com/en/articles/cloning-a-repository) this repository to your local system, into the place where you want to perform the data analysis.

(Make sure to have the right access / SSH Key. If **not**, follow the steps:
Step 1: https://docs.github.com/en/github/authenticating-to-github/connecting-to-github-with-ssh/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent

Step 2: https://docs.github.com/en/github/authenticating-to-github/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account)


    git clone git@github.com:NBChub/snakemake-metabolomics.git

### Step 2: Configure workflow

Configure the workflow according to your needs via editing the files in the `config/` folder. Adjust `config.yaml` to configure the workflow execution, and `samples.tsv` to specify the samples (files) that will be processed + analysed. 

**Suggestion: Use the Jupyter notebook [Create_sampletsv_file](./Create_sampletsv_file.ipynb) after you add all your files in the data/raw/ directory**

`samples.tsv` example:

|  sample_name |       comment                |
|-------------:|-----------------------------:|
| NBC_00162    | pyracrimicin                 |
| MDNA_WGS_14  | epemicins_A_B                |


Further formatting rules can be defined in the `workflow/schemas/` folder.


### Step 3: Create a conda environment& install snakemake

Installing conda
Get the latest version of conda with the following command:

    wget -O - https://www.anaconda.com/distribution/ 2>/dev/null | sed -ne 's@.*\(https:\/\/repo\.anaconda\.com\/archive\/Anaconda3-.*-Linux-x86_64\.sh\)\">64-Bit (x86) Installer.*@\1@p' | xargs wget


Installing Snakemake using [Mamba](https://github.com/mamba-org/mamba) is advised. In case you don’t use [Mambaforge](https://github.com/conda-forge/miniforge#mambaforge) you can always install [Mamba](https://github.com/mamba-org/mamba) into any other Conda-based Python distribution with:

    conda install -n base -c conda-forge mamba

Then install Snakemake with:

    mamba create -c conda-forge -c bioconda -n snakemake snakemake

For installation details, see the [instructions in the Snakemake documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).


### Step 4: Execute workflow

Activate the conda environment:

    conda activate snakemake

Install mono with sudo (https://www.mono-project.com/download/stable/#download-lin):

    sudo apt install mono-devel

If sudo cannot find the package, then follow the directions in the above link for the Ubuntu version that you work with.

Press enter (RETURN) to continue and install wget if it is not already installed in your Linux machine:
    
    brew install wget


#### Get example input data (only for testing the workflow with the example dataset)

    (cd data && wget https://zenodo.org/record/5511115/files/raw.zip && unzip *.zip -d raw)


#### Get the necessary executables (OpenMS binaries, ThermoRawFileParser & sirius executables):

    (cd resources && wget https://abibuilder.informatik.uni-tuebingen.de/archive/openms/OpenMSInstaller/release/latest/OpenMS-2.7.0-Debian-Linux-x86_64.deb && sudo apt-get install sudo apt-get install gdebi
    && sudo gdebi OpenMS-2.7.0-Debian-Linux-x86_64.deb)

If you encounter errors installing OpenMS, follow the steps indicated here: https://abibuilder.informatik.uni-tuebingen.de/archive/openms/Documentation/nightly/html/install_linux_bin.html
    
    (cd resources/ThermoRawFileParser && wget https://github.com/compomics/ThermoRawFileParser/releases/download/v1.3.4/ThermoRawFileParser.zip && unzip ThermoRawFileParser.zip)
    
    (cd resources/Sirius/ && wget https://github.com/boecker-lab/sirius/releases/download/v4.9.9/sirius-4.9.9-linux64-headless.zip  && unzip *.zip)

Test your configuration by performing a dry-run via

    snakemake --use-conda -n

Execute the workflow locally via

    snakemake --use-conda --cores $N

See the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executable.html) for further details.


### Step 5: Investigate results

All the results are in a csv format and can be opened simply with excel or using pandas dataframes. 

## Developer Notes
### Config & Schemas

* [Config & schemas](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html) define the input formatting and are important to generate `wildcards`. The idea of using `samples` and `units` came from [here](https://github.com/snakemake-workflows/dna-seq-gatk-variant-calling). I think we should use `units.txt` as a central metadata of the runs which are regulary updated (and should be the same for all use case). Then, `samples.txt` can be used to decide which strains need to be assembled per use case. 

### Rules

* [Snakefile](workflow/Snakefile): the main entry of the pipeline which tells the final output to be generated and the rules being used
* [common.smk](workflow/rules/common.smk): a rule that generate the variable used (strain names) & other helper scripts
* [The main rules (*.smk)](workflow/rules/): is the bash code that has been chopped into modular units, with defined input & output. Snakemake then chain this rules together to generate required jobs. This should be intuitive and makes things easier for adding / changing steps in the pipeline.

### Environments

* Conda environments are defined as .yaml file in `workflow/envs`
* Note that not all dependencies are compatible/available as conda libraries. Once installed, the virtual environment are stored in `.snakemake/conda` with unique hashes. The ALE and pilon are example where environment needs to be modified / dependencies need to be installed.
* It might be better to utilise containers / dockers and cloud execution for "hard to install" dependencies
* Custom dependencies and databases are stored in the `resources/` folder.
* Snakemake dependencies with conda packages is one of the drawbacks and why [Nextflow](https://www.nextflow.io/) might be more preferable. Nevertheless, the pythonic language of snakemake enable newcomers to learn and develop their own pipeline faster.

### Test Data (only for testing the workflow with the example dataset)

* Current test data are built from known metabolite producer strains or standard samples that have been analysed with a Thermo IDX mass spectrometer. The presence of the metabolites and their fragmentation patterns has been manually confirmed using TOPPView.

### Citations

Röst, H.L., Sachsenberg, T., Aiche, S., Bielow, C., Weisser, H., Aicheler, F., Andreotti, S., Ehrlich, H.-C., Gutenbrunner, P., Kenar, E., Liang, X., Nahnsen, S., Nilse, L., Pfeuffer, J., Rosenberger, G., Rurik, M., Schmitt, U., Veit, J., Walzer, M., Wojnar, D., Wolski, W.E.,Schilling, O., Choudhary, J.S., Malmström, L., Aebersold, R., Reinert, K., Kohlbacher, O. OpenMS: A flexible open-source software platform for mass spectrometry data analysis. Nature Methods, vol. 13, 2016. doi:10.1038/nmeth.3959

Kai Dührkop, Markus Fleischauer, Marcus Ludwig, Alexander A. Aksenov, Alexey V. Melnik, Marvin Meusel, Pieter C. Dorrestein, Juho Rousu, and Sebastian Böcker, SIRIUS 4: Turning tandem mass spectra into metabolite structure information. Nature Methods 16, 299–302, 2019 doi:10.1038/s41592-019-0344-8

Kai Dührkop, Huibin Shen, Marvin Meusel, Juho Rousu, and Sebastian Böcker, Searching molecular structure databases with tandem mass spectra using CSI:FingerID, PNAS October 13, 2015 112 (41) 12580-12585, doi:10.1073/pnas.1509788112

Nothias, L.-F., Petras, D., Schmid, R. et al. Feature-based molecular networking in the GNPS analysis environment. Nat. Methods 17, 905–908 (2020).
