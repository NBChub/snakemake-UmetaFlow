
######## snakemake preamble start (automatically inserted, do not edit) ########
import sys; sys.path.extend(['/Users/eeko/mambaforge/envs/snakemake/lib/python3.9/site-packages', '/Users/eeko/Desktop/Metabolomics/Snakemake/workflow/scripts']); import pickle; snakemake = pickle.loads(b'\x80\x04\x95\xdd\x04\x00\x00\x00\x00\x00\x00\x8c\x10snakemake.script\x94\x8c\tSnakemake\x94\x93\x94)\x81\x94}\x94(\x8c\x05input\x94\x8c\x0csnakemake.io\x94\x8c\nInputFiles\x94\x93\x94)\x81\x94(\x8c)results/preprocessed/Epemicins.featureXML\x94\x8c*results/precursorcorrection/Epemicins.mzML\x94\x8c!./resources/Contents/MacOS/sirius\x94\x8c\x18data/mzML/Epemicins.mzML\x94e}\x94(\x8c\x06_names\x94}\x94\x8c\x12_allowed_overrides\x94]\x94(\x8c\x05index\x94\x8c\x04sort\x94eh\x13\x8c\tfunctools\x94\x8c\x07partial\x94\x93\x94h\x06\x8c\x19Namedlist._used_attribute\x94\x93\x94\x85\x94R\x94(h\x19)}\x94\x8c\x05_name\x94h\x13sNt\x94bh\x14h\x17h\x19\x85\x94R\x94(h\x19)}\x94h\x1dh\x14sNt\x94bub\x8c\x06output\x94h\x06\x8c\x0bOutputFiles\x94\x93\x94)\x81\x94(\x8c\x1eresults/SIRIUS/Epemicins.mzTab\x94\x8c#results/CSIFingerID/Epemicins.mzTab\x94e}\x94(h\x0f}\x94h\x11]\x94(h\x13h\x14eh\x13h\x17h\x19\x85\x94R\x94(h\x19)}\x94h\x1dh\x13sNt\x94bh\x14h\x17h\x19\x85\x94R\x94(h\x19)}\x94h\x1dh\x14sNt\x94bub\x8c\x06params\x94h\x06\x8c\x06Params\x94\x93\x94)\x81\x94}\x94(h\x0f}\x94h\x11]\x94(h\x13h\x14eh\x13h\x17h\x19\x85\x94R\x94(h\x19)}\x94h\x1dh\x13sNt\x94bh\x14h\x17h\x19\x85\x94R\x94(h\x19)}\x94h\x1dh\x14sNt\x94bub\x8c\twildcards\x94h\x06\x8c\tWildcards\x94\x93\x94)\x81\x94\x8c\tEpemicins\x94a}\x94(h\x0f}\x94\x8c\x07samples\x94K\x00N\x86\x94sh\x11]\x94(h\x13h\x14eh\x13h\x17h\x19\x85\x94R\x94(h\x19)}\x94h\x1dh\x13sNt\x94bh\x14h\x17h\x19\x85\x94R\x94(h\x19)}\x94h\x1dh\x14sNt\x94bhJhGub\x8c\x07threads\x94K\x01\x8c\tresources\x94h\x06\x8c\tResources\x94\x93\x94)\x81\x94(K\x01K\x01\x8c0/var/folders/c_/ysz9v_bd1yb7h3ymmkn6m199jbv7x7/T\x94e}\x94(h\x0f}\x94(\x8c\x06_cores\x94K\x00N\x86\x94\x8c\x06_nodes\x94K\x01N\x86\x94\x8c\x06tmpdir\x94K\x02N\x86\x94uh\x11]\x94(h\x13h\x14eh\x13h\x17h\x19\x85\x94R\x94(h\x19)}\x94h\x1dh\x13sNt\x94bh\x14h\x17h\x19\x85\x94R\x94(h\x19)}\x94h\x1dh\x14sNt\x94bh]K\x01h_K\x01hahZub\x8c\x03log\x94h\x06\x8c\x03Log\x94\x93\x94)\x81\x94}\x94(h\x0f}\x94h\x11]\x94(h\x13h\x14eh\x13h\x17h\x19\x85\x94R\x94(h\x19)}\x94h\x1dh\x13sNt\x94bh\x14h\x17h\x19\x85\x94R\x94(h\x19)}\x94h\x1dh\x14sNt\x94bub\x8c\x06config\x94}\x94\x8c\x04rule\x94\x8c\x12sirius_csifingerid\x94\x8c\x0fbench_iteration\x94N\x8c\tscriptdir\x94\x8c;/Users/eeko/Desktop/Metabolomics/Snakemake/workflow/scripts\x94ub.'); from snakemake.logging import logger; logger.printshellcmds = False; __real_file__ = __file__; __file__ = '/Users/eeko/Desktop/Metabolomics/Snakemake/workflow/scripts/6_sirius_csifingerid.py';
######## snakemake preamble end #########
from pyopenms import *

def sirius_csifingerid(filename):
    fmap = FeatureMap()
    FeatureXMLFile().load(filename, fmap)
    exp = MSExperiment()
    MzMLFile().load(snakemake.input[1], exp)

    PrecursorCorrection.correctToNearestFeature(fmap, exp, 0.0, 100.0, True, False, False, False, 3, 0)

    sirius_algo = SiriusAdapterAlgorithm()

    sirius_algo_par = sirius_algo.getDefaults()

    sirius_algo_par.setValue("preprocessing:filter_by_num_masstraces", 2) 
    sirius_algo_par.setValue("preprocessing:precursor_mz_tolerance", 10.0) #default
    sirius_algo_par.setValue("preprocessing:precursor_mz_tolerance_unit", "ppm")
    sirius_algo_par.setValue("preprocessing:precursor_rt_tolerance", 5.0) #default
    sirius_algo_par.setValue("preprocessing:feature_only", "true")
    sirius_algo_par.setValue("sirius:profile", "orbitrap")
    sirius_algo_par.setValue("sirius:db", "none")
    sirius_algo_par.setValue("sirius:ions_considered", "[M+H]+, [M-H2O+H]+, [M+Na]+, [M+NH4]+")
    sirius_algo_par.setValue("sirius:candidates", 10)
    sirius_algo_par.setValue("sirius:elements_enforced", "CHNOS") 
    sirius_algo_par.setValue("project:processors", 2)
    sirius_algo_par.setValue("fingerid:db", "BIO")
    sirius_algo.setParameters(sirius_algo_par)
    
    fm_info = FeatureMapping_FeatureMappingInfo()
    feature_mapping = FeatureMapping_FeatureToMs2Indices() 
    sirius_algo.preprocessingSirius(filename,
                                    exp,
                                    fm_info,
                                    feature_mapping)

    sirius_algo.logFeatureSpectraNumber(filename, 
                                        feature_mapping,
                                        exp)


    msfile = SiriusMSFile()
    debug_level = 3
    sirius_tmp = SiriusTemporaryFileSystemObjects(debug_level)
    siriusstring= String(sirius_tmp.getTmpMsFile())

    feature_only = sirius_algo.isFeatureOnly()
    isotope_pattern_iterations = sirius_algo.getIsotopePatternIterations()
    no_mt_info = sirius_algo.isNoMasstraceInfoIsotopePattern()
    compound_info = []

    msfile.store(exp, 
                String(siriusstring),
                feature_mapping, 
                feature_only,
                isotope_pattern_iterations, 
                no_mt_info, 
                compound_info)

    executable= snakemake.input[2]
    out_csifingerid = snakemake.output[1]
    subdirs = sirius_algo.callSiriusQProcess((siriusstring),
                                            String(sirius_tmp.getTmpOutDir()),
                                            String(executable),
                                            String(out_csifingerid),
                                            False) # debug option not yet implemented

    candidates = sirius_algo.getNumberOfSiriusCandidates()
    sirius_result = MzTab()
    siriusfile = MzTabFile()
    SiriusMzTabWriter.read(subdirs,
                            snakemake.input[3],
                            candidates,
                            sirius_result)
    siriusfile.store(snakemake.output[0], sirius_result)

    top_hits= 5
    csi_result=MzTab()
    csi_file=MzTabFile()
    CsiFingerIdMzTabWriter.read(subdirs,
                        snakemake.input[3],
                        top_hits,
                        csi_result)

    csi_file.store(snakemake.output[1], csi_result)

sirius_csifingerid(snakemake.input[0])




