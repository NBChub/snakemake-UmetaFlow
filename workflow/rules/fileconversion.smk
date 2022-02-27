# Convert files from Thermo
# Convert profile MS1 data from Thermo raw files to centroid (MS1 and MS2) mzml

rule mzml_conversion:
    input:
        "data/raw/{samples}.raw"
    output:
        "data/mzML/{samples}.mzML" 
    shell:
        """
        mono resources/ThermoRawFileParser/ThermoRawFileParser.exe -i={input} -b={output}
        """
