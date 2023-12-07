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
  ;07-12-23: Adapted to enable reading of GALAXYTREE outputs too
  ;
"""

import numpy as np
from robs_redshift_from_snapnum import robs_redshift_from_snapnum

def make_lgals_sample(G_lgal, FILE_TYPE, SampType, Simulation, Cosmology, Hubble_h, DMParticleMass, ParticleMassRes, StellarMassRes, \
                      FullRedshiftList, RedshiftsToRead, select_main_progenitors=0) : 
    tH_z0 = ((1./(100.*Hubble_h))*3.086e+19)*3.17098e-8 #Hubble time at z=0 [in years]
    if (FILE_TYPE == 'snapshots') | ((FILE_TYPE == 'galtree') & (select_main_progenitors == 0)) : 
        #Sample selection criteria:      
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
    elif ((FILE_TYPE == 'galtree') & (select_main_progenitors == 1)) :  
        #Sample selection criteria:
        Gz = G_lgal[np.around(G_lgal['Redshift'], decimals=2) == FullRedshiftList[0]]
        if SampType == 'All' :
            G_sampz = Gz[(Gz['Mvir']*Hubble_h >= DMParticleMass*ParticleMassRes) \
                          & (Gz['StellarMass'] >= StellarMassRes)]
        elif SampType == 'Discs' :
            G_sampz = Gz[(Gz['Mvir']*Hubble_h >= DMParticleMass*ParticleMassRes) \
                          & (Gz['StellarMass'] >= StellarMassRes) \
                          & (Gz['Type'] != 2) \
                          & (Gz['BulgeMass']/Gz['StellarMass'] < 0.3) \
                          & (Gz['Sfr']/Gz['StellarMass'] >= ((2.*(1.+robs_redshift_from_snapnum(Simulation, Cosmology, Gz['SnapNum']))**2)/tH_z0)/10.)]
        elif SampType == 'ETGs' :
            G_sampz = Gz[(Gz['Mvir']*Hubble_h >= DMParticleMass*ParticleMassRes) \
                          & (Gz['StellarMass'] >= StellarMassRes) \
                          & (Gz['Type'] != 2) \
                          & (Gz['BulgeMass']/Gz['StellarMass'] > 0.7) \
                          & (Gz['Sfr']/Gz['StellarMass'] < ((2.*(1.+robs_redshift_from_snapnum(Simulation, Cosmology, Gz['SnapNum']))**2)/tH_z0)/10.)]
        elif SampType == 'Dwarfs' :
            G_sampz = Gz[(Gz['Mvir']*Hubble_h >= DMParticleMass*150.) \
                          & (Gz['StellarMass'] <= 1.e9) \
                          & (Gz['StellarMass'] > 0.0) \
                          & (Gz['Type'] < 2) \
                          & (Gz['MassWeightAge'] >= 0.0)]
        elif SampType == 'MWAs' :
            G_sampz = Gz[(np.log10(Gz['StellarMass']) >= 10.2) \
                          & (np.log10(Gz['StellarMass']) <= 10.8) \
                          & (Gz['BulgeMass']/Gz['StellarMass'] < 0.3)
                          & (Gz['Type'] == 0) \
                          & (Gz['Sfr'] >= 1.0) \
                          & (Gz['Sfr'] <= 5.0) \
                          & (Gz['MassWeightAge'] >= 0.0)]
        else : print("***** ERROR: Sample type not chosen. Please enter a valid string for SampType *****")                  
        print("Selecting all main progenitor branch for "+str(len(G_sampz))+" ["+SampType+"] galaxies (may take a few minutes)...")
        indices = np.where(G_lgal['MainLeafId'] == G_sampz['MainLeafId'][0])
        for iselec in range(1,len(G_sampz)) :
            indices = np.append(indices, np.where(G_lgal['MainLeafId'] == G_sampz['MainLeafId'][iselec]))
        G_samp = G_lgal[indices]
        print(str(len(G_samp))+" main progenitors selected.")
    else : print("***** ERROR: File type not chosen. Please enter either 'snapshots' or 'galtree' for FILE_TYPE *****") 
    
    print('Number of galaxies selected: ', len(G_samp), '\n')
    
    return G_samp
