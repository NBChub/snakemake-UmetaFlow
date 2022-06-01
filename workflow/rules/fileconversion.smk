# Convert files from Thermo
# Convert profile MS1 data from Thermo raw files to centroid (MS1 and MS2) mzml

import glob
from os.path import join

rule mzml_conversion:
    input:
        "data/raw/{samples}.raw"
    output:
        "data/mzML/{samples}.mzML" 
    params:
        exec_path = glob.glob(join('.snakemake','conda','*','bin','ThermoRawFileParser.exe'))[0]
    shell:
        """
        mono {params.exec_path} -i={input} -b={output}
        """
