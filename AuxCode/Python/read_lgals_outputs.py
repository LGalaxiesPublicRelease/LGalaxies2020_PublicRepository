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
  ;07-12-23: Adapted to enable reading of GALAXYTREE outputs too
  ;15-07-26: Adapted to enable the use of the automatically-generated struct
  ;15-07-26: Adapted to enable reading HDF5 output files
  ;
"""

#Base packages:
import numpy as np
from importlib import reload
# import astropy
# from astropy.io import fits
import h5py
import sys
sys.path.append('./Robs_python_routines/')

#Local packages
import procedures
reload (procedures)
from procedures import read_snap, read_tree
from robs_snapnum_from_redshift import robs_snapnum_from_redshift


#################
def read_lgals_outputs(OutputDir, Hubble_h, SIMULATION, FILE_TYPE, STRUCT_TYPE, MODEL, VERSION, \
                       COSMOLOGY, HDF5_OUTPUT, MULTI_REDSHIFT_PLOTS, GENERAL_PLOTS, PAPER_PLOTS, GALTREE_REDSHIFT_TO_PLOT, \
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

    # Read-in HDF5 output file(s):
    if HDF5_OUTPUT == 1:
        datalist = []
        if (FILE_TYPE == "snapshots"):
            the_snapnums = robs_snapnum_from_redshift(SIMULATION, COSMOLOGY, FullRedshiftList)
            for iredshift in range(0,len(FullRedshiftList)):
                if RedshiftsToRead[iredshift]: 
                    char_redshift="%0.2f" % FullRedshiftList[iredshift]
                    the_snapnum = str(the_snapnums[iredshift]) 
                    for ifile in range(FirstFile, LastFile + 1):
                        # Collect data from each treefile of each snapshot into a single list:
                        filename = 'SA_' + model_suffix +'_z' + char_redshift + "_%d" % ifile
                        print("Loading HDF5 snapshot file: "+filename+'\n')
                        with h5py.File(OutputDir + filename + ".h5", "r") as f:
                            datalist.append(f[the_snapnum][:])
        elif (FILE_TYPE == "galtree"):
            for ifile in range(FirstFile, LastFile + 1):
                # Collect data from each treefile of each snapshot into a single list:
                filename = 'SA_' + model_suffix +'_%d' % ifile
                print("Loading HDF5 galtree file: "+filename+'\n')
                with h5py.File(OutputDir + filename + ".h5", "r") as f:
                    # Only select galaxies of a certain redshift (GALTREE_REDSHIFT_TO_PLOT), if required:
                    if (MULTI_REDSHIFT_PLOTS != 1) & ((GENERAL_PLOTS == 1) | (PAPER_PLOTS == 1)) :
                        print("\nSelecting only galaxies at redshift = "+str(GALTREE_REDSHIFT_TO_PLOT)+"\n")
                        snap_to_select = int(robs_snapnum_from_redshift(SIMULATION, COSMOLOGY, [GALTREE_REDSHIFT_TO_PLOT]))
                        the_snapnum = f['GalTree']['SnapNum'][:]
                        the_samp = f['GalTree'][the_snapnum == snap_to_select]
                        datalist.append(the_samp)
                    else:
                        datalist.append(f['GalTree'][:])

        # Concatenate data list into one structured numpy array:        
        G_lgal = np.concatenate(datalist) 
        print("HDF5 file structure:\n")
        for name, (dtype, offset) in G_lgal.dtype.fields.items():
            print(f"{name:30s} {dtype}")
        print("\n")

    # Read-in binary output file(s):  
    else:
        if FILE_TYPE == 'snapshots' :
            if STRUCT_TYPE == 'auto' :
                from auto_LGalaxy_struct import LGalaxiesStruct
                from auto_LGalaxy_struct import properties_used
            elif STRUCT_TYPE == 'liteOutput' :
                from LGalaxy_snapshots_liteOutput import LGalaxiesStruct
                from LGalaxy_snapshots_liteOutput import properties_used
            elif STRUCT_TYPE == 'normal' :
                from LGalaxy_snapshots_normal import LGalaxiesStruct
                from LGalaxy_snapshots_normal import properties_used
            elif STRUCT_TYPE == 'ringSFHs' :
                from LGalaxy_snapshots_ringSFHs import LGalaxiesStruct
                from LGalaxy_snapshots_ringSFHs import properties_used
            elif STRUCT_TYPE == 'liteOutput_noDust' :
                from LGalaxy_snapshots_liteOutput_noDust import LGalaxiesStruct
                from LGalaxy_snapshots_liteOutput_noDust import properties_used
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
        elif FILE_TYPE == 'galtree' :  
            if STRUCT_TYPE == 'auto' :
                from auto_LGalaxy_struct import LGalaxiesStruct
                from auto_LGalaxy_struct import properties_used
            elif STRUCT_TYPE == 'liteOutput' :
                from LGalaxy_galtree_liteOutput import LGalaxiesStruct
                from LGalaxy_galtree_liteOutput import properties_used
            elif STRUCT_TYPE == 'normal' :
                from LGalaxy_galtree_normal import LGalaxiesStruct
                from LGalaxy_galtree_normal import properties_used
            elif STRUCT_TYPE == 'ringSFHs' :
                from LGalaxy_galtree_ringSFHs import LGalaxiesStruct
                from LGalaxy_galtree_ringSFHs import properties_used
            elif STRUCT_TYPE == 'liteOutput_noDust' :
                from LGalaxy_galtree_liteOutput_noDust import LGalaxiesStruct
                from LGalaxy_galtree_liteOutput_noDust import properties_used
            # elif STRUCT_TYPE == '[ADD_YOUR_OWN_STRUCT_TYPE_HERE]' :
            #     from LGalaxy_galtree_new import LGalaxiesStruct
            #     from LGalaxy_galtree_new import properties_used
            else : print("***** ERROR: Output structure type unknown. Please choose from: LG2020, plusBinaries *****") 
            (G_lgal) = read_tree(OutputDir, FirstFile, LastFile, properties_used, LGalaxiesStruct, model_suffix) 
        else : print("***** ERROR: File type unknown. Please choose from: snapshots, galtree *****")
    
    print('\nReading done')
    
        
    #################
    #Convert properties to common units:
    #Masses [Msun]:
    mass_props = ['Mvir', 'CentralMvir', 'HaloM_Crit200', 'HaloM_TopHat', 'MassFromInSitu', 'MassFromMergers', \
                  'MassFromBursts', 'ColdGas', 'HotGas', 'StellarMass', 'DiskMass', 'BulgeMass', 'EjectedMass', \
                  'BlackHoleMass', 'ICM', 'ColdGasRings', 'DiskMassRings', 'BulgeMassRings', 'MetalsColdGas', \
                  'MetalsHotGas', 'MetalsDiskMass', 'MetalsBulgeMass', 'MetalsEjectedMass', 'MetalsICM', \
                  'MetalsColdGasRings', 'MetalsDiskMassRings', 'MetalsBulgeMassRings', 'sfh_DiskMass', \
                  'sfh_BulgeMass', 'sfh_ICM', 'sfh_MetalsDiskMass', 'sfh_MetalsBulgeMass', 'sfh_MetalsICM',\
                  'sfh_DiskMassRings', 'sfh_BulgeMassRings', 'sfh_MetalsDiskMassRings', 'sfh_MetalsBulgeMassRings']
    for prop in mass_props:
        if prop in G_lgal.dtype.names:
            G_lgal[prop] = (G_lgal[prop]*1.e10)/Hubble_h  
            
    #Lengths & positions [Mpc]:
    len_pos_props = ['Rvir', 'Pos', 'DiskRadius', 'ColdGasRadius', 'StellarHalfMassRadius', 'StellarHalfLightRadius']
    for prop in len_pos_props:
        if prop in G_lgal.dtype.names:
            G_lgal[prop] = G_lgal[prop]/Hubble_h 
    
    print('Unit conversions done')
    
    
    # #################
    # #Read-in SFH bin info:
    # SFH_bins = astropy.io.fits.open(BaseDir+'AuxCode/Python/'+'Database_SFH_table.fits')
    # SFH_bins = SFH_bins[1].data
    
    # print('SFH bins read')
    # print('-------------\n')
    
    
    #################
    return G_lgal #, SFH_bins
    