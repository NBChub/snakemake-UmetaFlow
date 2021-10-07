from pyopenms import *
from pandas import DataFrame
import pandas as pd
import pyteomics
from pyteomics.openms import featurexml
import numpy as np
import sys
from pyteomics import mztab

example_filename= "Epemicins.mzML"
exp = MSExperiment()
MzMLFile().load(example_filename, exp)
exp.sortSpectra(True)
    
delta_mzs= []
mzs = []
rts= []
PrecursorCorrection.correctToHighestIntensityMS1Peak(exp, 100.0, True, delta_mzs, mzs, rts)

# Run mass trace detection
mass_traces = []
mtd = MassTraceDetection()
mtd_par = mtd.getDefaults()
mtd_par.setValue("mass_error_ppm", 10.0) 
mtd_par.setValue("noise_threshold_int", 1.0e04)
mtd_par.setValue("chrom_peak_snr", 3.0)
mtd_par.setValue("chrom_fwhm", 1.5)
mtd_par.setValue("min_trace_length", 3.0)
mtd_par.setValue("max_trace_length", 60.0)
mtd.setParameters(mtd_par)
mtd.run(exp, mass_traces, 0)
# Run elution peak detection
mass_traces_split = []
mass_traces_final = []
epd = ElutionPeakDetection()
epd_par = epd.getDefaults()
epd_par.setValue("width_filtering", "auto")
epd_par.setValue("min_fwhm", 1.0)
epd_par.setValue("max_fwhm", 30.0)
epd.setParameters(epd_par)
epd.detectPeaks(mass_traces, mass_traces_split)
    
if (epd.getParameters().getValue("width_filtering") == "fixed"):
     epd.filterByPeakWidth(mass_traces_split, mass_traces_final)
else:
    mass_traces_final = mass_traces_split

# Run feature detection
feature_map_FFM = FeatureMap()
feat_chrom = []
ffm = FeatureFindingMetabo()
ffm_par = ffm.getDefaults() 
ffm_par.setValue("isotope_filtering_model", "none")
ffm_par.setValue("remove_single_traces", "false")
ffm_par.setValue("mz_scoring_by_elements", "false")
ffm_par.setValue("report_convex_hulls", "true")
ffm.setParameters(ffm_par)
ffm.run(mass_traces_final, feature_map_FFM, feat_chrom)
feature_map_FFM.setUniqueIds()
fh = FeatureXMLFile()
fh.store('./EpemicinsFeatureFindingMetabo.featureXML', feature_map_FFM)

# Run metabolite adduct decharging detection
# With SIRIUS you are only able to use singly charged adducts
mfd = MetaboliteFeatureDeconvolution()
mdf_par = mfd.getDefaults()
mdf_par.setValue("potential_adducts",  [b"H:+:0.6",b"Na:+:0.2",b"NH4:+:0.1", b"H2O:-:0.1"])
mdf_par.setValue("charge_min", 1, "Minimal possible charge")
mdf_par.setValue("charge_max", 1, "Maximal possible charge")
mdf_par.setValue("charge_span_max", 1)
mdf_par.setValue("max_neutrals", 1)
mfd.setParameters(mdf_par)
    
feature_map_DEC = FeatureMap()
cons_map0 = ConsensusMap()
cons_map1 = ConsensusMap()
mfd.compute(feature_map_FFM, feature_map_DEC, cons_map0, cons_map1)
fxml = FeatureXMLFile()
fh.store('./EpemicinsMFD.featureXML', feature_map_FFM)

with featurexml.read('./EpemicinsMFD.featureXML') as f:
    features_list = [FXML for FXML in f]
    
df = pd.DataFrame() 
for feat in features_list:
    idx = feat['id']
    for key in feat.keys():
        if key == 'id':
           pass
        elif key == 'position':
            pos_list = feat['position']
            for pos in pos_list:
                if pos['dim'] == '0':
                    df.loc[idx, 'position_0'] = pos['position']
                elif pos['dim'] == '1':
                    df.loc[idx, 'position_1'] = pos['position']
        elif key == 'quality':
            qual_list = feat['quality']
            for qual in qual_list:
                if qual['dim'] == '0':
                    df.loc[idx, 'quality_0'] = qual['quality']
                elif qual['dim'] == '1':
                    df.loc[idx, 'quality_1'] = qual['quality']
        else:
            df.loc[idx, key] = feat[key]
df_tidy = df.rename(columns = {'position_0': 'RT', 'position_1': 'mz'}, inplace = False)
df_tidy=df_tidy.drop(columns= ["quality_0", "quality_1", "overallquality", "label", "legal_isotope_pattern"])
df_tidy.reset_index(drop=True, inplace=True) 
df_tidy.to_csv("EpemicinsFFM.csv")

# Precursor corrector

PrecursorCorrection.correctToNearestFeature(feature_map_DEC, exp, 0.0, 100.0, True, False, False, False, 3, 0)

# SIRIUS & CSI_FingerID

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
    
featureinfo = "./EpemicinsMFD.featureXML"
fm_info = FeatureMapping_FeatureMappingInfo()
feature_mapping = FeatureMapping_FeatureToMs2Indices() 
sirius_algo.preprocessingSirius(featureinfo,
                                exp,
                                fm_info,
                                feature_mapping)
sirius_algo.logFeatureSpectraNumber(featureinfo, 
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
             String(sirius_tmp.getTmpMsFile()),
             feature_mapping, 
             feature_only,
             isotope_pattern_iterations, 
             no_mt_info, 
             compound_info)
out_csifingerid = "./Epemicins_csifingerID.mzTab" 
executable= "/Users/eeko/Desktop/software/Contents/MacOS/sirius"
subdirs = sirius_algo.callSiriusQProcess(String(sirius_tmp.getTmpMsFile()),
                                         String(sirius_tmp.getTmpOutDir()),
                                         String(executable),
                                         String(out_csifingerid),
                                         False)
candidates = sirius_algo.getNumberOfSiriusCandidates()
sirius_result = MzTab()
siriusfile = MzTabFile()
SiriusMzTabWriter.read(subdirs,
                       filename,
                       candidates,
                       sirius_result)
siriusfile.store("./Epemicins_sirius.mzTab", sirius_result) 
siriusfile= "./Epemicins_sirius.mzTab"
sirius=  pyteomics.mztab.MzTab(siriusfile, encoding='UTF8', table_format='df')
sirius.metadata
df= sirius.small_molecule_table
SIRIUS_DF= df.drop(columns= ["identifier", "smiles", "inchi_key", "description", "calc_mass_to_charge", "charge", "taxid", "species","database", "database_version", "spectra_ref", "search_engine", "modifications"])
SIRIUS_DF.to_csv("./Epemicins_SIRIUS_DF.csv")