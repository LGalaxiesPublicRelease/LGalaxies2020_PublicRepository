# -*- coding: utf-8 -*-
"""
Created on Thu Aug 29 16:41:24 2024

@author: ry22aas
"""

import numpy as np

def read_SFH_Bins(FileDir, Simulation) :

    if Simulation == 'Mil-I' :
        FileName = "SFH_Bins_MilI"
    elif Simulation == 'Mil-II' :
        FileName = "SFH_Bins_MilII"
    else :
        print("***** ERROR: read_SFH_Bins.py: SIMULATION unknown, so correct SFH_BIN file is unknown. *****") 
    
    SFH_Struct = np.dtype([
    ('snapnum',np.int32), # snapshot number
    ('bin',np.int32), # index of current bin
    ('lookbacktime',np.float64), # yr # lookback time to centre of current bin # Already corrected for Hubble_h and UnitTime_in_years inside L-Galaxies
    ('dt',np.float64), # yr # width of the current bin # Already corrected for Hubble_h and UnitTime_in_years inside L-Galaxies
    ('nbins',np.int32) # number of highest resolution bins used to create current bin
    ])
    
    f = open(FileDir+FileName,"rb") # open file
    tot_nbins = np.fromfile(f, np.int32, count=1)[0] # Total number of rows in the file
    SFH_bins = np.fromfile(f, SFH_Struct, tot_nbins) # Load data into numpy structure
    f.close() # close file
    
    SFH_bins["bin"] = SFH_bins["bin"]+1 # SFH bin numbers count from 0 in these binary files, so increment them all by 1 so that the final SFH bin number for each snapshot gives the total number of SFH bins in that snapshot.
    # SFH_bins["lookbacktime"] = SFH_bins["lookbacktime"]/1.e9 # Gyr # Convert lookback times from yr to Gyr
    # SFH_bins["dt"] = SFH_bins["dt"]/1.e9 # Gyr # Convert bin widths from yr to Gyr
    
    return SFH_bins

    # SFH_bins contains:
    #
    # SFH_bins["snapnum"]
    # SFH_bins["bin"]
    # SFH_bins["lookbacktime"]
    # SFH_bins["dt"]
    # SFH_bins["nbins"]
