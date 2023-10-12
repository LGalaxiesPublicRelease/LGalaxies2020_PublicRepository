# -*- coding: utf-8 -*-
"""
Created on Thu Nov  4 15:22:51 2021

@author: robyates
"""

"""
make_lgals_sample.py
  ;Script to make a sample of galaxies from the read-in output data.
  ;
  ;Rob Yates 04-11-2021
  ;
  ;08-11-22: Adapted for use at the L-Galaxies workshop 2022
  ;12-10-23: Adapted for use with the Yates+23 version of L-Galaxies
  ;
"""

import numpy as np
from robs_redshift_from_snapnum import robs_redshift_from_snapnum

def make_lgals_sample(G_lgal, FILE_TYPE, SampType, Simulation, Cosmology, Hubble_h, DMParticleMass, ParticleMassRes, StellarMassRes, \
                      FullRedshiftList, RedshiftsToRead) : 
    if FILE_TYPE == 'snapshots' : 
        #Sample selection criteria:
        tH_z0 = ((1./(100.*Hubble_h))*3.086e+19)*3.17098e-8 #Hubble time at z=0 [in years]
        if SampType == 'All' :
            G_samp = G_lgal[(G_lgal['Mvir']*Hubble_h >= DMParticleMass*ParticleMassRes) \
                          & (G_lgal['StellarMass'] >= StellarMassRes)]              
        elif SampType == 'Discs' :            
            G_samp = G_lgal[(G_lgal['Mvir']*Hubble_h >= DMParticleMass*ParticleMassRes) \
                          & (G_lgal['StellarMass'] >= StellarMassRes) \
                          & (G_lgal['Type'] != 2) \
                          & (G_lgal['BulgeMass']/G_lgal['StellarMass'] < 0.3) \
                          & (G_lgal['Sfr']/G_lgal['StellarMass'] >= ((2.*(1.+robs_redshift_from_snapnum(Simulation, Cosmology, G_lgal['SnapNum']))**2)/tH_z0)/10.)]
        elif SampType == 'ETGs' :
            G_samp = G_lgal[(G_lgal['Mvir']*Hubble_h >= DMParticleMass*ParticleMassRes) \
                          & (G_lgal['StellarMass'] >= StellarMassRes) \
                          & (G_lgal['Type'] != 2) \
                          & (G_lgal['BulgeMass']/G_lgal['StellarMass'] > 0.7) \
                          & (G_lgal['Sfr']/G_lgal['StellarMass'] < ((2.*(1.+robs_redshift_from_snapnum(Simulation, Cosmology, G_lgal['SnapNum']))**2)/tH_z0)/10.)]     
        elif SampType == 'Dwarfs' :
            G_samp = G_lgal[(G_lgal['Mvir']*Hubble_h >= DMParticleMass*150.) \
                          & (G_lgal['StellarMass'] >= StellarMassRes) \
                          & (G_lgal['StellarMass'] <= 1.e9) \
                          & (G_lgal['StellarMass'] > 0.0) \
                          & (G_lgal['Type'] < 2) \
                          & (G_lgal['MassWeightAge'] >= 0.0)]
        else : print("***** ERROR: Sample type not chosen. Please enter a valid string for SampType *****") 
    elif FILE_TYPE == 'galtree' :
        yy = 0.0 #Not tested yet! (04-11-21)
    else : print("***** ERROR: File type not chosen. Please enter either 'snapshots' or 'galtree' for FILE_TYPE *****") 
    
    print('Number of galaxies selected: ', len(G_samp), '\n')
    
    return G_samp
