#Explanation of columns
#mz= mass-to-charge ratio (m/z)
#RT= retention time (sec)
#intensity = intensity of the feature (AU-arbitrary units)
#FWHM= Full Width of the peak at Half its Maximum height
#num_of_masstraces = number of mass traces detected (single mass traces are excluded). This is relevant to the isotopic pattern
#isotope_distances = distance in mz between the isotopes (jumps of app. 1 is important to confirm that this is a real feature)

from pyopenms import *

from collections import defaultdict
from functools import reduce
from pathlib import Path
from time import perf_counter
import sys

from IPython.core.display import display
from pandas import CategoricalDtype
import numpy as np
from pyopenms import *
import pandas as pd
import os

class FeatureMapDF(FeatureMap):      
    def __init__(self):
        super().__init__()

    def get_df(self):
        def gen(fmap: FeatureMap, fun):
            for f in fmap:
                yield from fun(f)

        def extractMetaData(f: Feature):
            metavals = []
            types = []
            for k in metavals:
                if k == b"num_of_masstraces":
                    types.append('f')
                elif k== b"FWHM":
                    types.append('f')
                elif k== b"max_height":
                    types.append('i4')
                elif k== b"masstrace_intensity":
                    types.append("f")
            # subfeatures = f.getFeatureList()  # type: list[FeatureHandle]
            pep = f.getPeptideIdentifications()  # type: list[PeptideIdentification]
            bb = f.getConvexHull().getBoundingBox2D()
            for u in metavals:
                f.getMetaValue(u)
            if len(pep) != 0:
                hits = pep[0].getHits()
                if len(hits) != 0:
                    besthit = hits[0]  # type: PeptideHit
                    yield f.getUniqueId(), besthit.getSequence().toString(), f.getCharge(), f.getRT(), f.getMZ(), bb[0][0], bb[1][0], f.getMetaValue("num_of_masstraces"), f.getOverallQuality(), f.getIntensity()
                else:
                    yield f.getUniqueId(), None, f.getCharge(), f.getRT(), f.getMZ(), bb[0][0], bb[1][0], f.getMetaValue("num_of_masstraces"), f.getOverallQuality(), f.getIntensity()
            else:
                yield f.getUniqueId(), None, f.getCharge(), f.getRT(), f.getMZ(), bb[0][0], bb[1][0], f.getMetaValue("num_of_masstraces"), f.getOverallQuality(), f.getIntensity()

        cnt = self.size()

        mddtypes = [('id', np.dtype('uint64')), ('sequence', 'U200'), ('charge', 'i4'), ('RT', 'f'), ('mz', 'f'),
                    ('RTstart', 'f'), ('RTend', 'f'), ("num_of_masstraces", 'f'),
                    ('quality', 'f'), ('intensity', 'f')]
        mdarr = np.fromiter(iter=gen(self, extractMetaData), dtype=mddtypes, count=cnt)
        df= pd.DataFrame(mdarr).set_index('id').sort_values("mz").drop(columns= "sequence")
        df= df[df["num_of_masstraces"]>=2]
        return df

def preprocessDF(filename):
    fmap = FeatureMapDF()
    FeatureXMLFile().load(filename, fmap)
    DF= fmap.get_df()
    DF.to_csv(snakemake.output[0])

filename=snakemake.input[0]