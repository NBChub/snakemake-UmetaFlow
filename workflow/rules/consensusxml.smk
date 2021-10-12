#MapAlignerPoseClustering is used to perform a linear retention time alignment, basically correct for linear shifts in retention time.
#add trafoXML files as an output also (TransformationXMLFile()) for transforming also the MS2 spectra later on
rule MapAlignerPoseClustering:
    input:
        expand("results/{samples}/interim/MFD_{samples}.featureXML", samples=SAMPLES)
    output:
        expand("results/Consensus/interim/MapAlignerPoseClustering_{samples}.featureXML", samples=SAMPLES)
    conda:
        "../envs/pyopenms.yaml"   
    script:
        "../scripts/8_map_aligner.py"

#Introduce the features to a protein identification file (idXML)- the only way to create a ConsensusXML file currently (run FeatureLinkerUnlabeledKD)       
rule IDMapper:
    input:
        "resources/emptyfile.idXML",
        "results/Consensus/interim/MapAlignerPoseClustering_{samples}.featureXML",
        "results/{samples}/interim/precursorcorrected_{samples}.mzML"
    output:
        "results/Consensus/interim/IDMapper_{samples}.featureXML"
    conda:
        "../envs/pyopenms.yaml"   
    script:
        "../scripts/9_idmapper.py"

#The FeatureLinkerUnlabeledKD is used to aggregate the feature information (from single files) into a ConsensusFeature, linking features from different files together, which have a smiliar m/z and rt (no MS2 data).
rule FeatureLinkerUnlabeledKD:
    input:
        expand("results/Consensus/interim/IDMapper_{samples}.featureXML", samples=SAMPLES)
    output:
        "results/Consensus/interim/FeatureLinkerUnlabeledKD.consensusXML"
    conda:
        "../envs/pyopenms.yaml"   
    script:
        "../scripts/10_featurelinker.py"

#Filter out the features that do not have an MS2 pattern
rule FileFilter:
    input:
        "results/Consensus/interim/FeatureLinkerUnlabeledKD.consensusXML"
    output:
        "results/Consensus/filtered.consensusXML"
    conda:
        "../envs/pyopenms.yaml"   
    script:
        "../scripts/11_filefilter.py"

