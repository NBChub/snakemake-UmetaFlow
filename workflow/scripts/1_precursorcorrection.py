#Selection of the peak with the highest intensity as corrected precursor mass in a given mass range (e.g. precursor mass +/- 0.2 Da)
# For each MS2 spectrum the corresponding MS1 spectrum is determined by using the rt information of the precursor. 
# In the MS1, the peak with the highest intensity in a given mass range to the uncorrected precursor m/z is selected and used as corrected precursor m/z.
# Using our thermo orbitrap instruments, the fragmentation is by default done to the highest intensity peaks of the MS1 spectrum, but the instrument could still assign the wrone procursor. 

from pyopenms import *
def precursorcorrect(filename):
    exp = MSExperiment()
    MzMLFile().load(filename, exp)
    exp.sortSpectra(True)
    delta_mzs= []
    mzs = []
    rts= []
    PrecursorCorrection.correctToHighestIntensityMS1Peak(exp, 100.0, True, delta_mzs, mzs, rts)
    MzMLFile().store(snakemake.output[0], exp)

precursorcorrect(snakemake.input[0])