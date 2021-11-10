# Metabolomics workflow for Linux systems
[![Snakemake](https://img.shields.io/badge/snakemake-≥6.7.0-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Build Status](https://travis-ci.org/snakemake-workflows/snakemake-bgc-analytics.svg?branch=master)](https://travis-ci.org/snakemake-workflows/snakemake-bgc-analytics)

This is a snakemake implementation of the [Metabolomics OpenMS workflow](snakemake-metabolomics/OpenMS_workflow.ipynb) tailored by [Eftychia Eva Kontou](https://github.com/eeko-kon)
## Workflow overview
![dag](/images/pyOpenMS_workflow.svg)

View the workflow with interactive comments in lucid: https://lucid.app/lucidchart/4dc81d37-bca3-4b2d-8253-33341ac79ab4/edit?viewport_loc=-71%2C53%2C2422%2C1416%2C0_0&invitationId=inv_5c1c0383-052a-4905-8146-dd842ee528fb 

## Usage
### Step 1: Clone the workflow

[Clone](https://help.github.com/en/articles/cloning-a-repository) this repository to your local system, into the place where you want to perform the data analysis.

(Make sure to have the right access / SSH Key. If **not**, follow the steps:
Step 1: https://docs.github.com/en/github/authenticating-to-github/connecting-to-github-with-ssh/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent

Step 2: https://docs.github.com/en/github/authenticating-to-github/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account)


    git clone git@github.com:NBChub/snakemake-metabolomics.git

### Step 2: Configure workflow

Configure the workflow according to your needs via editing the files in the `config/` folder. Adjust `config.yaml` to configure the workflow execution, and `samples.tsv` to specify the samples (files) that will be processed + analysed. 

**Suggestion: Always add the genome ID in the filename for easier connection later to BGC_analytics.**

`samples.tsv` example:

|  sample_name |       comment                |
|-------------:|-----------------------------:|
| NBC_00162    | pyracrimicin                 |
| MDNA_WGS_14  | epemicins_A_B                |


Further formatting rules will be defined in the `workflow/schemas/` folder.


### Step 3: Create a conda environment& install snakemake

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
  
Press enter (RETURN) to continue 
    
    brew install wget

#### Get input data (only for testing the workflow with the example dataset)

    (cd data && wget https://zenodo.org/record/5511115/files/raw.zip && unzip *.zip -d raw)

#### Get the necessary executables (ThermoRawFileParser & sirius):
    
    (cd resources/ThermoRawFileParser && wget https://github.com/compomics/ThermoRawFileParser/releases/download/v1.3.4/ThermoRawFileParser.zip && unzip ThermoRawFileParser.zip)
    
    (cd resources/Sirius/ && wget https://github.com/boecker-lab/sirius/releases/download/v4.9.3/sirius-4.9.3-linux64-headless.zip  && unzip *.zip)

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
* Current test data are built from real runs of known metabolite producer strains or standard samples that have been already alanysed with the GUI Software Freestyle and confirmed the presence of fragmentation patterns for the specific metabolites


