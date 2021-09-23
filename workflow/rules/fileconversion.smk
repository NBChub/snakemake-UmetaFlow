rule mzml_conversion:
    input:
        "data/raw/{samples}.raw"
    output:
        "results/{samples}/interim/{samples}.mzML" 
    shell:
        """
        mono resources/ThermoRawFileParser/ThermoRawFileParser.exe -i={input} -b={output}
        """
