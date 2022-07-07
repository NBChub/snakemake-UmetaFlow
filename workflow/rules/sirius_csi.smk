import glob
from os.path import join

# 1) SIRIUS generates formula predictions from scores calculated from 1) MS2 fragmentation scores (ppm error + intensity) and 2) MS1 isotopic pattern scores.        
#    The CSI_fingerID function is another algorithm from the Boecher lab, just like SIRIUS adapter and is using the formula predictions from SIRIUS, to search in structural libraries and predict the structure of each formula.
# "CSI:FingerID identifies the structure of a compound by searching in a molecular structure database. “Structure” refers to the identity and connectivity (with bond multiplicities) of the atoms, but no stereochemistry information. Elucidation of stereochemistry is currently beyond the power of automated search engines."
# "CANOPUS predicts compound classes from the molecular fingerprint predicted by CSI:FingerID without any database search involved. Hence, it provides structural information for compounds for which neither spectral nor structural reference data are available."

if config["rules"]["requantification"]==True:
    rule sirius:
        input: 
            var1= "results/GNPSexport/mzML/Aligned_{samples}.mzML",
            var2= "results/Interim/Requantification/MFD_{samples}.featureXML" 
        output:
            out1= "results/Interim/sirius/formulas_{samples}.mzTab",
            out2= "results/Interim/sirius/structures_{samples}.mzTab"
        params:
            exec_path = glob.glob(join('.snakemake','conda','exe', 'bin', 'sirius'))
        threads: 4
        shell:
            """
            /Users/eeko/openms-develop/openms_build/bin/SiriusAdapter -sirius_executable {params.exec_path} -in {input.var1} -in_featureinfo {input.var2} -out_sirius {output.out1} -out_fingerid {output.out2} -preprocessing:filter_by_num_masstraces 2 -preprocessing:feature_only -sirius:profile orbitrap -sirius:db none -sirius:ions_considered "[M+H]+, [M-H2O+H]+, [M+Na]+, [M+NH4]+" -sirius:elements_enforced CHN[15]OS[4]Cl[2]P[2] -debug 3 -fingerid:candidates 5 -project:processors {threads} -threads {threads}
            """
else:
    rule sirius:
        input: 
            var1= "results/GNPSexport/mzML/Aligned_{samples}.mzML",
            var2= "results/Interim/preprocessed/MFD_{samples}.featureXML" 
        output:
            out1= "results/Interim/sirius/formulas_{samples}.mzTab",
            out2= "results/Interim/sirius/structures_{samples}.mzTab"
        params:
            exec_path = glob.glob(join('.snakemake','conda','exe', 'bin', 'sirius'))
        threads: 4
        shell:
            """
            /Users/eeko/openms-develop/openms_build/bin/SiriusAdapter -sirius_executable {params.exec_path} -in {input.var1} -in_featureinfo {input.var2} -out_sirius {output.out1} -out_fingerid {output.out2} -preprocessing:filter_by_num_masstraces 2 -preprocessing:feature_only -sirius:profile orbitrap -sirius:db none -sirius:ions_considered "[M+H]+, [M-H2O+H]+, [M+Na]+, [M+NH4]+" -sirius:elements_enforced CHN[15]OS[4]Cl[2]P[2] -debug 3 -fingerid:candidates 5 -project:processors {threads} -threads {threads}
            """

# 2) Convert the mzTab to a tsv file

rule df_sirius:
    input: 
        "results/Interim/sirius/formulas_{samples}.mzTab",
        "results/Interim/sirius/structures_{samples}.mzTab"
    output:
        "results/SIRIUS/formulas_{samples}.tsv",
        "results/CSI/structures_{samples}.tsv"
    conda:
        "../envs/openms.yaml"
    script:
        "../scripts/df_SIRIUS_CSI.py"

