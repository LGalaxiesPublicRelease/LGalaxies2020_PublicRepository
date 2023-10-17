# -*- coding: utf-8 -*-
"""
Created on Wed May  3 19:08:57 2023

@author: ry22aas
"""
#Basic packages:
import numpy as np

#Local packages:
import sys
sys.path.append('/vol/ph/astro_data2/rmyates/python/python_libraries/Robs_Routines/')
from robs_lumdist_from_redshift import robs_lumdist_from_redshift

# Adapted from robs_arcsec_to_kpc.pro IDL script:
def robs_arcsec_to_kpc(r_arcsec, redshift, dist=None, OmegaM=0.315, OmegaL=0.683, Hubble_h=0.673, \
                       res=1000, comoving=None, verbose=None) :
    #Function to calculate a sky-plane PROPER or COMOVING distance (e.g. projected diameter of an extended object) in kpc, 
    #given the angular distance and redshift. (Gives same result as arcsec2parsec.joseonorbe.com to within at least 
    #3 sig. figs. [~4 sig. figs. for nearby objects, and ~6 sig. figs. for large & nearby objects]).
    #
    #INPUTS:
    #r_arcsec: arcsec: The radial distance to be converted to kpc
    #redshift: The redshift of object.
    #dist: kpc: LoS distance to the object in question. If not known, the luminosity distance will be calculated from the redshift
    #OmegaM: Assumed matter density in the universe. [Default: Planck value of 0.315]
    #OmegaL: Assumed dark-energy density in the universe. [Default: Planck value of 0.683]
    #Hubble_h: Dimensionless Hubble parameter at z=0. [Default: Planck value of 0.673]
    #comoving: Will cause the outputted distance to be in comoving units. Default is to return the radial distance in proper units.
    #
    #OUTPUTS:
    #r_kpc: kpc: The PROPER or COMOVING radial distance.
    #
    #NON-STANDARD REQUIREMENTS:
    #robs_lumdist_from_redshift()
    #
    #---------------   
    r_radians = r_arcsec * np.pi/(180.*60.*60.) #*4.84814e-6 #
    scaleFactor = 1./(1.+redshift)
    if dist is None :
        dist = robs_lumdist_from_redshift(redshift, OmegaM=OmegaM, OmegaL=OmegaL, Hubble_h=Hubble_h, res=res)*1.e3 #in kpc
    r_kpc = scaleFactor**2 * r_radians * dist
    if verbose :
        print("robs_arcsec_to_kpc:")
        print("-------------------")
        print("r_arcsec = ", r_arcsec)
        print("redshift = ", redshift)
        print("dist (kpc) =  ", dist)
        print("r_kpc =  ", r_kpc)
        print(" ")
    if comoving :
        return r_kpc*(1.+redshift) 
    else :
        return r_kpc
