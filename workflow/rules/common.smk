import pandas as pd
from snakemake.utils import validate
from snakemake.utils import min_version

min_version("5.18.0")

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
#singularity: "docker://continuumio/miniconda3"

##### load config and sample sheets #####

configfile: "config/config.yaml"
validate(config, schema="../schemas/config.schema.yaml")

# set up sample
samples = pd.read_csv(config["samples"], sep="\t").set_index("sample_name", drop=False)
samples.index.names = ["samples"]

## watermark snakemake version to samples
#for i in samples.index:
#    samples.loc[i, "sample_name"] = str(samples.loc[i, "sample_name"]) + __version__

validate(samples, schema="../schemas/samples.schema.yaml")


##### Wildcard constraints #####
wildcard_constraints:
    sample="|".join(samples.index),

##### Helper functions #####


SAMPLES = samples.sample_name.to_list()
