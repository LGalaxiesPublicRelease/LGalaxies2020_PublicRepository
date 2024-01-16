# -*- coding: utf-8 -*-
"""
Created on Wed May  3 19:07:29 2023

@author: ry22aas
"""
# Adapted from robs_lumdist_from_redshift.pro IDL script.
#Basic packages:
import numpy as np
from scipy.integrate import simps

def robs_lumdist_from_redshift(redshift, OmegaM=0.315, OmegaL=0.683, OmegaK=0.0, OmegaR=0.0, Hubble_h=0.673, \
                               res=1000) :
    #Function to calculate luminosity distance from the redshift.
    #
    #EXAMPLE:
    #To print the luminosity distance (in Mpc) for an object at redshift 0.5:
    #print(robs_lumdist_from_redshift(0.5))
    #
    #NON-STANDARD REQUIREMENTS:
    #
    #--------------- 
    sol = 299792.458 #speed of light [km/s] 
    if isinstance(redshift, float): redshift = np.array([redshift]) #Convert float to single-element float array, so len() can be used.
    lum_dist = np.full(len(redshift),np.nan)
    for ii in range(len(redshift)) :
        z_dash = np.linspace(0.0, redshift[ii], res)
        Efunc = np.sqrt(OmegaL + OmegaK*(1.+z_dash)**2 + OmegaM*(1.+z_dash)**3 + OmegaR*(1.+z_dash)**4)
        lum_dist[ii] = (((1.+redshift[ii])*sol)/(Hubble_h*100.)) * simps(1./Efunc, z_dash)
    return lum_dist #in Mpc
    
  