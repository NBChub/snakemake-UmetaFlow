#  PrecursorCorrection (To the "nearest feature”): This algorithm is used after feature detection& adduct grouping. It basically allows the precursor correction on MS2 level.
#Which means that if there are MS2 spectra in my feature space which have been measured in isotope traces, it “corrects” the MS2 spectrum annotation to the monoisotopic trace. That is why you have a high mass deviation 100 pm, but 0.0 rt tolerance. So it basically corrects the MS2 to the feature centroid that can be found/mapped by SIRIUS::preprocessing later on.
#Wrong assignment of the mono-isotopic mass for precursors are assumed: if precursor_mz matches the mz of a non-monoisotopic feature mass trace and
#in the case that believe_charge is true: if feature_charge matches the precursor_charge In the case of wrong mono-isotopic assignment several options for correction are available: keep_original will create a copy of the precursor and tandem spectrum for the new mono-isotopic mass trace and retain the original one. all_matching_features does this not for only the closest feature but all features in a question.

# SIRIUS adapter: The algorithm generates formula prediction from scores calculated from 1) MS2 fragmentation scores (ppm error + intensity) and 2) MS1 isotopic pattern scores.
# It can only compute feautures that are singly charged. There is also a timeout for compounds (compound timeout so that it doesn't compute for longer than 100 seconds per feature, which normally happens with larger molecules).
#-sirius:compound_timeout
#Maximal computation time in seconds for a single compound. 0 for an infinite amount of time. (default: '100' min: '0')
#This algorith can help in data dereplication and analysis for direct library search.


from pyopenms import *
def sirius(filename):
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
    sirius_algo_par.setValue("sirius:candidates", 5)
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
    out_csifingerid = "" 
    subdirs = sirius_algo.callSiriusQProcess((siriusstring),
                                            String(sirius_tmp.getTmpOutDir()),
                                            String(executable),
                                            String(out_csifingerid),
                                            False) # debug option not yet implemented

    candidates = sirius_algo.getNumberOfSiriusCandidates()
    sirius_result = MzTab()
    siriusfile = MzTabFile()
    SiriusMzTabWriter.read(subdirs,
                            snakemake.input[1],
                            candidates,
                            sirius_result)
    siriusfile.store(snakemake.output[0], sirius_result)


sirius(snakemake.input[0])