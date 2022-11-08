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
  ;08-11-22: This version was adapted for use at the L-Galaxies workshop 2022
  ;
"""

import numpy as np

def make_lgals_sample(G_lgal, FILE_TYPE, SampType, Hubble_h, DMParticleMass, ParticleMassRes, StellarMassRes, \
                      FullRedshiftList, RedshiftsToRead) : 
    if FILE_TYPE == 'snapshots' : 
        #Sample selection criteria:
        if SampType == 'All' :
            G_samp = G_lgal[(G_lgal['Mvir']*Hubble_h >= DMParticleMass*ParticleMassRes) \
                          & (G_lgal['StellarMass'] >= StellarMassRes)]
                          #& (G_lgal['SnapNum'] == snap_z0)]                     
        elif SampType == 'Discs' :
            tH_z0 = ((1./(100.*Hubble_h))*3.086e+19)*3.17098e-8 #Hubble time at z=0 [in years]
            sSFR_MS_H20 = (2.*(1.+FullRedshiftList[0])**2)/tH_z0            
            G_samp = G_lgal[(G_lgal['Mvir']*Hubble_h >= DMParticleMass*ParticleMassRes) \
                          & (G_lgal['StellarMass'] >= 1.e9) \
                          & (G_lgal['Type'] != 2) \
                          & (G_lgal['BulgeMass']/G_lgal['StellarMass'] < 0.3) \
                          & (G_lgal['Sfr']/G_lgal['StellarMass'] >= 10**(np.log10(sSFR_MS_H20)-1.))]
                          #& (G_lgal['Sfr']/G_lgal['StellarMass'] >= 1.259e-11)]        
        elif SampType == 'Dwarfs' :
            G_samp = G_lgal[(G_lgal['Mvir']*Hubble_h >= DMParticleMass*150.) \
                          & (G_lgal['StellarMass'] <= 1.e9) \
                          & (G_lgal['StellarMass'] > 0.0) \
                          & (G_lgal['Type'] < 2) \
                          & (G_lgal['MassWeightAge'] >= 0.0) \
                          & (G_lgal['ColdGas_elements'][:,0] > 0.0) \
                          & (G_lgal['ColdGas_elements'][:,4] > 0.0) \
                          & (G_lgal['DiskMass_elements'][:,0] > 0.0) \
                          & (G_lgal['DiskMass_elements'][:,4] > 0.0) \
                          & (G_lgal['DiskMass_elements'][:,10] > 0.0) \
                          & (G_lgal['BulgeMass_elements'][:,0] >= 0.0) \
                          & (G_lgal['BulgeMass_elements'][:,4] >= 0.0) \
                          & (G_lgal['BulgeMass_elements'][:,10] >= 0.0)]            
        else : print("***** ERROR: Sample type not chosen. Please enter a valid string for SampType *****") 
    elif FILE_TYPE == 'galtree' :
        yy = 0.0 #Not tested yet! (04-11-21)
    else : print("***** ERROR: File type not chosen. Please enter either 'snapshots' or 'galtree' for FILE_TYPE *****") 
    
    print('Number of galaxies selected: ', len(G_samp), '\n')
    
    return G_samp
