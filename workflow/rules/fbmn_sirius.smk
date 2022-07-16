# This rule is for integrating the FBMN results with the SIRIUS annotations.
# After the FBMN jon is done, download the graphml file under the directory workflow/GNPSexport and run the following rule:
rule graphml:
    input:
        input_graphml= glob.glob("results/GNPSexport/*.graphml")
    output:    
        output_graphml= "results/GNPSexport/fbmn_network_sirius.graphml"
    log: "workflow/report/logs/GNPSexport/fbmn_sirius.log"
    threads: 4
    conda:
        "../envs/pyopenms.yaml"
    shell:
        """
        python workflow/scripts/FBMN_SIRIUS.py {input.input_graphml} {output.output_graphml} 2>> {log}
        """