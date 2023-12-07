# -*- coding: utf-8 -*-
"""
Created on Fri Sep 16 12:04:18 2022

@author: robyates
"""

"""
robs_snapnum_from_redshift.py
  ;Function to return a numpy array of snapshot numbers corresponding to the inputted numpy array of redshifts
  ;
  ;Rob Yates 16-09-2022
  ;
  ;UPDATES: 
  ;08-11-22: This version was adapted for use at the L-Galaxies workshop 2022
  ;06-12-23: Changed the match between the redshifts listed below and the inputted redshifts so that they are 
  ;          compared after being rounded to the *same* number of decimal places.
  ;
"""

#Basic packages:
import numpy as np


##########
def robs_snapnum_from_redshift(sim, cos, redshifts) :
    if (sim == 'Mil-I') :
        if (cos == 'Planck') :
            All_redshifts = [56.418465,49.24166,44.754864,26.812078,17.76479,16.264746,14.901791,13.662391,
                            12.534529,11.507145,10.570324,9.715587,8.934729,8.220924,7.567806,6.969652,
                            6.421481,5.918548,5.456737,5.032382,4.642067,4.282843,3.951917,3.646797,
                            3.365268,3.105276,2.865033,2.642801,2.437088,2.246648,2.070046,1.906216,
                            1.754237,1.612945,1.481728,1.359643,1.245949,1.140232,1.041466,0.949732,
                            0.863829,0.783975,0.709089,0.639476,0.574025,0.513287,0.456079,0.40291,
                            0.352951,0.306211,0.262623,0.221449,0.183387,0.147548,0.11378,0.082661,
                            0.053316,0.025612,2.64E-4,-0.023571,-0.046002,-0.066833,-0.085948,-0.103822]
    elif (sim == 'Mil-II') :
        if (cos == 'Planck') :
            All_redshifts = [56.418465,49.24166,37.83254,29.547876,26.951857,24.604391,22.4798,20.555162,\
                            17.76479,16.264746,14.901791,13.662391,12.534529,11.507145,10.570324,9.715587,\
                            8.934729,8.220924,7.567806,6.969652,6.421481,5.918548,5.456737,5.032382,4.642067,\
                            4.282843,3.951917,3.646797,3.365268,3.105276,2.865033,2.642801,2.437088,2.246648,\
                            2.070046,1.906216,1.754237,1.612945,1.481728,1.359643,1.245949,1.140232,1.041466,\
                            0.949732,0.863829,0.783975,0.709089,0.639476,0.574025,0.513287,0.456079,0.40291,\
                            0.352951,0.306211,0.262623,0.221449,0.183387,0.147548,0.11378,0.082661,0.053316,\
                            0.025612,2.64E-4,-0.023571,-0.046002,-0.066833,-0.085948,-0.103822]

    All_redshifts_dp = np.around(All_redshifts, decimals=decimal_places)
    redshifts_dp = np.around(redshifts, decimals=decimal_places)
    snapnums = np.empty(len(redshifts), dtype=int)
    for iz in range(0,len(redshifts)) :
        #snap = np.where(All_redshifts_dp == redshifts[iz])
        snap = np.where(All_redshifts_dp == redshifts_dp[iz])
        snapnums[iz] = snap[0] #.astype(int)
    
    return snapnums
    
