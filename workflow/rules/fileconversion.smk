rule mgf_conversion:
    input:
        "data/raw/{samples}.raw"
    output:
        "results/{samples}/interim/{samples}.mgf.gzip"
    shell:
        """
        mono resources/ThermoRawFileParser/ThermoRawFileParser.exe -i={input} -b={output} -f=0 -g 
        """

rule mzml_conversion:
    input:
        "data/raw/{samples}.raw"
    output:
        "results/{samples}/interim/{samples}.mzML" 
    shell:
        """
        mono resources/ThermoRawFileParser/ThermoRawFileParser.exe -i={input} -b={output}
        """
