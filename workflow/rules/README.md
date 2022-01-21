# Snakemake rules 

The pipeline consists of three separate workflows that are interconnected, and one data analysis guide:

1) File conversion: convert raw files from Thermo to open community-driven format mzML centroid (see documentation [here](https://github.com/compomics/ThermoRawFileParser))

2) Pre-processing: converting raw data to a feature table with a series of OpenMS algorithms (see documentation [here](https://abibuilder.informatik.uni-tuebingen.de/archive/openms/Documentation/nightly/html/index.html)). Additinaly, the feature table is introduced to SIRIUS and CSI:FingerID for formula and structural predictions (see documentation [here](https://boecker-lab.github.io/docs.sirius.github.io/)).

![dag](/images/Preprocessing+SIRIUS_CSI_FingerID.svg) 

3) GNPSexport: generate all the files necessary to create a FBMN job at GNPS (see documentation [here](https://ccms-ucsd.github.io/GNPSDocumentation/featurebasedmolecularnetworking-with-openms/)). 

![dag](/images/GNPSExport.svg) 

4) Re-quantification: missing value imputation resulted by the pre-processing workflow using OpenMS algorithms. Generate a FeatureMatrix.

![dag](/images/Re-quantification.svg) 

5) Data analysis: Annotate the FeatureMatrix with GNPS metabolites, formula and structural predictions, and perform outlier analysis to find features detected less frequently in the matrix.

![dag](/images/Data_analysis.svg) 