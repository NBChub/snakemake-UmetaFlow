rule mgf_conversion:
    input:
        "data/raw/{samples}.raw"
    output:
        "data/mgf/{samples}.mgf"
    shell:
        """
        mono resources/ThermoRawFileParser/ThermoRawFileParser.exe -i={input} -b={output} -f=0 -g 
        """

rule mzml_conversion:
    input:
        "data/raw/{samples}.raw"
    output:
        "data/mzML/{samples}.mzML" 
    shell:
        """
        mono resources/ThermoRawFileParser/ThermoRawFileParser.exe -i={input} -b={output}
        """
