import glob
from os.path import join

# 1) Remove convex hulls from the feature files

rule noconvexhulls:
    input:
        "results/Interim/preprocessed/MFD_{samples}.featureXML"
    output:
        "results/Interim/sirius/MFD_nch_{samples}.featureXML"
    threads: 4
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/nch.py"  


# 3) SIRIUS generates formula predictions from scores calculated from 1) MS2 fragmentation scores (ppm error + intensity) and 2) MS1 isotopic pattern scores.        
#    The CSI_fingerID function is another algorithm from the Boecher lab, just like SIRIUS adapter and is using the formula predictions from SIRIUS, to search in structural libraries and predict the structure of each formula.
# "CSI:FingerID identifies the structure of a compound by searching in a molecular structure database. “Structure” refers to the identity and connectivity (with bond multiplicities) of the atoms, but no stereochemistry information. Elucidation of stereochemistry is currently beyond the power of automated search engines."
# "CANOPUS predicts compound classes from the molecular fingerprint predicted by CSI:FingerID without any database search involved. Hence, it provides structural information for compounds for which neither spectral nor structural reference data are available."

rule sirius:
    input: 
        var1= "results/Interim/mzML/PCfeature_{samples}.mzML", 
        var2= "results/Interim/sirius/MFD_nch_{samples}.featureXML"        
    output:
        "results/Interim/sirius/formulas_{samples}.mzTab",
        "results/Interim/sirius/structures_{samples}.mzTab"
    params:
        exec_path = glob.glob(join('.snakemake','conda','*','bin','sirius'))[0],
    conda:
        "../envs/python.yaml"
    threads: 4
    shell:
        """
        SiriusAdapter -sirius_executable {params.exec_path} -in {input.var1} -in_featureinfo {input.var2} -out_sirius {output[0]} -out_fingerid {output[1]} -preprocessing:filter_by_num_masstraces 2 -preprocessing:feature_only -sirius:profile orbitrap -sirius:db none -sirius:ions_considered "[M+H]+, [M-H2O+H]+, [M+Na]+, [M+NH4]+" -sirius:elements_enforced CHN[15]O[40]S[4]P[3] -debug 3 -fingerid:candidates 5 -project:processors {threads} -threads {threads}
        """

# 5) Convert the mzTab to a csv file

rule df_sirius:
    input: 
        "results/Interim/sirius/formulas_{samples}.mzTab",
        "results/Interim/sirius/structures_{samples}.mzTab"
    output:
        "results/SIRIUS/formulas_{samples}.csv",
        "results/CSI/structures_{samples}.csv"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/df_SIRIUS_CSI.py"

