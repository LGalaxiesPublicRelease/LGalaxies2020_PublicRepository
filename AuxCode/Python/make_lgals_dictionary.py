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

def make_lgals_dictionary(LABEL, G_samp, PlotDir, COSMOLOGY, SIMULATION, FILE_TYPE, MODEL, \
                          VERSION, NumGals, SAMPLE_TYPE, char_lower_redshift, char_upper_redshift, \
                          RNUM, RingRadii, RReRings, RReSFHRings, RingArea, \
                          SFH_bins_num, SFH_bins_lbt_allGals) : 
    Samp = {
            "G_samp": G_samp,
            "Cosmology": COSMOLOGY,
            "Simulation": SIMULATION,
            "File_type" : FILE_TYPE,
            "Model" : MODEL,
            "Version" : VERSION,
            "Sample_type" : SAMPLE_TYPE,
            "NumGals" : NumGals,
            "z_lower" : char_lower_redshift,
            "z_upper" : char_upper_redshift,
            "RNUM" : RNUM,
            "RingRadii" : RingRadii,
            "RReRings" : RReRings,
            "RReSFHRings" : RReSFHRings,
            "RingArea" : RingArea,
            "SFH_bins_num" : SFH_bins_num,
            "SFH_bins_lbt_allGals" : SFH_bins_lbt_allGals,
            #"mark" : mark,
            "Label" : LABEL
            }
    return Samp

