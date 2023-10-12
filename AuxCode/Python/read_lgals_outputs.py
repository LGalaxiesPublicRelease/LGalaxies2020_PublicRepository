# -*- coding: utf-8 -*-
"""
Created on Thu Nov  4 14:24:44 2021

@author: robyates
"""

"""
read_lgals_outputs.py
  ;Script to read-in L-Galaxies binary outputs and convert some units.
  ;Also reads-in SFH bin info.
  ;Calls the read_snap() or read_tree() functions from procedures.py.
  ;
  ;Rob Yates 04-11-2021
  ;
  ;08-11-22: Adapted for use at the L-Galaxies workshop 2022
  ;12-10-23: Adapted for use with the Yates+23 version of L-Galaxies
  ;
"""

#Base packages:
import numpy as np
from importlib import reload
import astropy
from astropy.io import fits

#Local packages
import procedures
reload (procedures)
from procedures import read_snap, read_tree
# import sys
# sys.path.append('../awk/output/python/')




#################
def read_lgals_outputs(BaseDir, OutputDir, Hubble_h, SIMULATION, FILE_TYPE, STRUCT_TYPE, MODEL, VERSION, \
                       FirstFile, LastFile, FullRedshiftList, RedshiftsToRead) :         
    if VERSION == '' :
        if MODEL == 'default' : 
            model_suffix = 'DM'
        elif MODEL == 'modified' : 
            model_suffix = 'MM'
        else : print("***** ERROR: Model unknown. Please choose from: default, modified *****")    
    else :
        if MODEL == 'default' : 
            model_suffix = 'DM'+"_"+VERSION
        elif MODEL == 'modified' : 
            model_suffix = 'MM'+"_"+VERSION
        else : print("***** ERROR: Model unknown. Please choose from: default, modified *****")  
    print('\n-------------')
    
    if FILE_TYPE == 'snapshots' :
        if STRUCT_TYPE == 'liteOutput' :
            from LGalaxy_snapshots_liteOutput import LGalaxiesStruct
            from LGalaxy_snapshots_liteOutput import properties_used
        elif STRUCT_TYPE == 'normal' :
            from LGalaxy_snapshots_normal import LGalaxiesStruct
            from LGalaxy_snapshots_normal import properties_used
        elif STRUCT_TYPE == 'ringSFHs' :
            from LGalaxy_snapshots_ringSFHs import LGalaxiesStruct
            from LGalaxy_snapshots_ringSFHs import properties_used
        # elif STRUCT_TYPE == '[ADD_YOUR_OWN_STRUCT_TYPE_HERE]' :
        #     from LGalaxy_snapshots_new import LGalaxiesStruct
        #     from LGalaxy_snapshots_new import properties_used
        else : print("***** ERROR: Ouput structure type unknown. Please choose from: LG2020, plusBinaries *****")  
        if SIMULATION == 'Mil-I' : 
            (G_lgal, SnapshotList) = read_snap(OutputDir, FirstFile, LastFile, \
                                               properties_used, LGalaxiesStruct, \
                                               RedshiftsToRead, FullRedshiftList, model_suffix)
        elif SIMULATION == 'Mil-II' :
            (G_lgal, SnapshotList) = read_snap(OutputDir+'MRII/', FirstFile, LastFile, \
                                               properties_used, LGalaxiesStruct, \
                                               RedshiftsToRead, FullRedshiftList, model_suffix)
    elif FILE_TYPE == 'galtree' :  #<--WARNING: Untested! (21-04-21) 
        if STRUCT_TYPE == 'liteOutput' :
            from LGalaxy_galtree_liteOutput import LGalaxiesStruct
            from LGalaxy_galtree_liteOutput import properties_used
        elif STRUCT_TYPE == 'normal' :
            from LGalaxy_galtree_normal import LGalaxiesStruct
            from LGalaxy_galtree_normal import properties_used
        elif STRUCT_TYPE == 'ringSFHs' :
            from LGalaxy_galtree_ringSFHs import LGalaxiesStruct
            from LGalaxy_galtree_ringSFHs import properties_used
        # elif STRUCT_TYPE == '[ADD_YOUR_OWN_STRUCT_TYPE_HERE]' :
        #     from LGalaxy_galtree_new import LGalaxiesStruct
        #     from LGalaxy_galtree_new import properties_used
        else : print("***** ERROR: Output structure type unknown. Please choose from: LG2020, plusBinaries *****") 
        (G_lgal) = read_tree(OutputDir,FirstFile,LastFile,properties_used,LGalaxiesStruct)    
        SnapshotList = np.zeros(len(FullRedshiftList),dtype=np.int32)
        for ii in range(0,len(FullRedshiftList)):                  
            G0=G_lgal[ np.rint(G_lgal['Redshift']*100.) == FullRedshiftList[ii]*100. ]             
            SnapshotList[ii]=G0['SnapNum'][0]
    else : print("***** ERROR: File type unknown. Please choose from: snapshots, galtree *****")
    
    print('\nReading done')
    
    
    
    
    
    #################
    #Convert properties to common units:
    #Masses [Msun]:
    G_lgal['Mvir'] = (G_lgal['Mvir']*1.e10)/Hubble_h 
    G_lgal['CentralMvir'] = (G_lgal['CentralMvir']*1.e10)/Hubble_h 
    G_lgal['HaloM_Crit200'] = (G_lgal['HaloM_Crit200']*1.e10)/Hubble_h 
    if not ("liteOutput" in (STRUCT_TYPE)) :
        G_lgal['HaloM_TopHat'] = (G_lgal['HaloM_TopHat']*1.e10)/Hubble_h 
        G_lgal['MassFromInSitu'] = (G_lgal['MassFromInSitu']*1.e10)/Hubble_h 
        G_lgal['MassFromMergers'] = (G_lgal['MassFromMergers']*1.e10)/Hubble_h 
        G_lgal['MassFromBursts'] = (G_lgal['MassFromBursts']*1.e10)/Hubble_h 
    G_lgal['ColdGas'] = (G_lgal['ColdGas']*1.e10)/Hubble_h 
    G_lgal['HotGas'] = (G_lgal['HotGas']*1.e10)/Hubble_h 
    G_lgal['StellarMass'] = (G_lgal['StellarMass']*1.e10)/Hubble_h 
    G_lgal['DiskMass'] = (G_lgal['DiskMass']*1.e10)/Hubble_h 
    G_lgal['BulgeMass'] = (G_lgal['BulgeMass']*1.e10)/Hubble_h 
    G_lgal['EjectedMass'] = (G_lgal['EjectedMass']*1.e10)/Hubble_h 
    G_lgal['BlackHoleMass'] = (G_lgal['BlackHoleMass']*1.e10)/Hubble_h 
    G_lgal['ICM'] = (G_lgal['ICM']*1.e10)/Hubble_h #N.B. this is actually the mass of the stellar halo (called 'ICM' for legacy reasons)
    G_lgal['ColdGasRings'] = (G_lgal['ColdGasRings']*1.e10)/Hubble_h 
    G_lgal['DiskMassRings'] = (G_lgal['DiskMassRings']*1.e10)/Hubble_h 
    G_lgal['BulgeMassRings'] = (G_lgal['BulgeMassRings']*1.e10)/Hubble_h 
    G_lgal['MetalsColdGas'] = (G_lgal['MetalsColdGas']*1.e10)/Hubble_h 
    G_lgal['MetalsHotGas'] = (G_lgal['MetalsHotGas']*1.e10)/Hubble_h 
    G_lgal['MetalsDiskMass'] = (G_lgal['MetalsDiskMass']*1.e10)/Hubble_h 
    G_lgal['MetalsBulgeMass'] = (G_lgal['MetalsBulgeMass']*1.e10)/Hubble_h 
    G_lgal['MetalsEjectedMass'] = (G_lgal['MetalsEjectedMass']*1.e10)/Hubble_h 
    G_lgal['MetalsICM'] = (G_lgal['MetalsICM']*1.e10)/Hubble_h #N.B. this is actually the mass of metals in the stellar halo (called 'MetalsICM' for legacy reasons)
    G_lgal['MetalsColdGasRings'] = (G_lgal['MetalsColdGasRings']*1.e10)/Hubble_h 
    G_lgal['MetalsDiskMassRings'] = (G_lgal['MetalsDiskMassRings']*1.e10)/Hubble_h 
    G_lgal['MetalsBulgeMassRings'] = (G_lgal['MetalsBulgeMassRings']*1.e10)/Hubble_h 
    if not ("liteOutput" in (STRUCT_TYPE)) :
        G_lgal['sfh_DiskMass'] = (G_lgal['sfh_DiskMass']*1.e10)/Hubble_h 
        G_lgal['sfh_BulgeMass'] = (G_lgal['sfh_BulgeMass']*1.e10)/Hubble_h 
        G_lgal['sfh_ICM'] = (G_lgal['sfh_ICM']*1.e10)/Hubble_h #N.B. this is actually the SFH for the stellar halo (called 'sfh_ICM' for legacy reasons)
        G_lgal['sfh_MetalsDiskMass'] = (G_lgal['sfh_MetalsDiskMass']*1.e10)/Hubble_h 
        G_lgal['sfh_MetalsBulgeMass'] = (G_lgal['sfh_MetalsBulgeMass']*1.e10)/Hubble_h 
        G_lgal['sfh_MetalsICM'] = (G_lgal['sfh_MetalsICM']*1.e10)/Hubble_h #N.B. this is actually the metal enrichment history for the stellar halo (called 'sfh_MetalsICM' for legacy reasons)
        if ("ringSFHs" in (STRUCT_TYPE)) :
            G_lgal['sfh_MetalsDiskMassRings'] = (G_lgal['sfh_MetalsDiskMassRings']*1.e10)/Hubble_h
            G_lgal['sfh_MetalsBulgeMassRings'] = (G_lgal['sfh_MetalsBulgeMassRings']*1.e10)/Hubble_h
    #N.B. Chemical element and dust masses are actually outputted in units of Msun (i.e. with Hubble_h already factored in, according to the COSMOLOGY L-Galaxies was run with).    
    
    #Lengths & positions [Mpc]:
    G_lgal['Rvir'] = G_lgal['Rvir']/Hubble_h
    G_lgal['Pos'] = G_lgal['Pos']/Hubble_h #positions are in comoving Mpc, so divide by (1+z) when comparing galaxy separations from their positions (delta_pos) to e.g. Rvir.
    G_lgal['DiskRadius'] = G_lgal['DiskRadius']/Hubble_h    
    G_lgal['ColdGasRadius'] = G_lgal['ColdGasRadius']/Hubble_h
    G_lgal['StellarHalfMassRadius'] = G_lgal['StellarHalfMassRadius']/Hubble_h
    if not ("liteOutput" in (STRUCT_TYPE)) :
        G_lgal['StellarHalfLightRadius'] = G_lgal['StellarHalfLightRadius']/Hubble_h
      
    #Spins [Mpc/h km/s] (unmodified)
    
    #Mass rates [Msun/yr] (unmodified)
    
    #Velocities [km/s] (unmodified)
    
    print('Unit conversions done')
    
    
    
    
    #################
    #Read-in SFH bin info:
    SFH_bins = astropy.io.fits.open(BaseDir+'AuxCode/Python/'+'Database_SFH_table.fits')
    SFH_bins = SFH_bins[1].data
    
    print('SFH bins read')
    print('-------------\n')
    
    
    
    
    return G_lgal, SFH_bins
    