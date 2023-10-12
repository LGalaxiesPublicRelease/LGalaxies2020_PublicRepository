# -*- coding: utf-8 -*-
"""
Created on Tue Nov  9 17:52:26 2021

@author: robyates
"""

"""
make_lgals_dictionary.py
  ;Makes dictionaries for the sample data needed to make plots.
  ;
  ;Rob Yates 09-11-2021
  ;
  ;08-11-22: This version was adapted for use at the L-Galaxies workshop 2022
  ;
"""
import numpy as np

def make_lgals_dictionary(LABEL, G_samp, PlotDir, COSMOLOGY, SIMULATION, FILE_TYPE, MODEL, \
                          VERSION, MRI_cutoff, NumGals, MRI_gals, MRII_gals, SAMPLE_TYPE, \
                          char_lower_redshift, char_upper_redshift, Volume_MRI, Volume_MRII, \
                          FullSnapnumList_MRI, FullSnapnumList_MRII, CALC_RINGS_AND_SFH_INFO, \
                          RNUM=None, RingRadii=None, RingCenRadii=None, RReRings=None, \
                          RReSFHRings=None, RingArea=None, SFH_bins_num=None, SFH_bins_lbt_allGals=None) : 
    Volumes = np.full(len(G_samp),0.0)
    Volumes[0:MRI_gals] = Volume_MRI
    Volumes[MRI_gals:] = Volume_MRII
    
    if CALC_RINGS_AND_SFH_INFO == 1 :
        Samp = {
                "G_samp": G_samp,
                "Cosmology": COSMOLOGY,
                "Simulation": SIMULATION,
                "File_type" : FILE_TYPE,
                "Model" : MODEL,
                "Version" : VERSION,
                "MRI_cutoff" : MRI_cutoff,
                "Sample_type" : SAMPLE_TYPE,
                "NumGals" : NumGals,
                "MRI_gals" : MRI_gals,
                "MRII_gals" : MRII_gals,
                "z_lower" : char_lower_redshift,
                "z_upper" : char_upper_redshift,
                "Volumes" : Volumes,
                "FullSnapnumList_MRI" : FullSnapnumList_MRI,
                "FullSnapnumList_MRII" : FullSnapnumList_MRII,
                "RNUM" : RNUM,
                "RingRadii" : RingRadii,
                "RingCenRadii" : RingCenRadii,
                "RReRings" : RReRings,
                "RReSFHRings" : RReSFHRings,
                "RingArea" : RingArea,
                "SFH_bins_num" : SFH_bins_num,
                "SFH_bins_lbt_allGals" : SFH_bins_lbt_allGals,
                "Label" : LABEL
                }
    else :
        Samp = {
                "G_samp": G_samp,
                "Cosmology": COSMOLOGY,
                "Simulation": SIMULATION,
                "File_type" : FILE_TYPE,
                "Model" : MODEL,
                "Version" : VERSION,
                "MRI_cutoff" : MRI_cutoff,
                "Sample_type" : SAMPLE_TYPE,
                "NumGals" : NumGals,
                "MRI_gals" : MRI_gals,
                "MRII_gals" : MRII_gals,
                "z_lower" : char_lower_redshift,
                "z_upper" : char_upper_redshift,
                "Volumes" : Volumes,
                "FullSnapnumList_MRI" : FullSnapnumList_MRI,
                "FullSnapnumList_MRII" : FullSnapnumList_MRII,
                "Label" : LABEL
                }
    return Samp
