# -*- coding: utf-8 -*-
"""
Created on Wed Jan 19 11:49:20 2022

@author: robyates
"""

"""
mass_checks.py
  ;Script to check:
  ;a) how many NaNs and negatives there are in key mass properties,
  ;b) if sub-component masses add up to component masses properly (e.g. SUM(ColdGas_elements) = ColdGas, etc)
  ;
  ;Rob Yates 19-01-2022
  ;
  ;08-11-22: This version was adapted for use at the L-Galaxies workshop 2022
  ;
"""

#Base packages:
import numpy as np
import matplotlib.pyplot as plt




#################
def mass_checks(G_lgal, COSMOLOGY, SIMULATION, FILE_TYPE, STRUCT_TYPE, MODEL, VERSION, SAMPLE_TYPE, plots=None) :  
    #Properties to check:
    prop_name = ['ColdGas','BulgeMass','DiskMass', \
                 'MetalsBulgeMass_SNII','MetalsBulgeMass_SNIa','MetalsBulgeMass_AGB', \
                 'MetalsDiskMass_SNII','MetalsDiskMass_SNIa','MetalsDiskMass_AGB', \
                 'MetalsDiskMassRings_SNII','MetalsDiskMassRings_SNIa','MetalsDiskMassRings_AGB', \
                 'ColdGas_elements']    
    prop = [G_lgal['ColdGas'],G_lgal['BulgeMass'],G_lgal['DiskMass'], \
            G_lgal['MetalsBulgeMass'][:,0],G_lgal['MetalsBulgeMass'][:,1],G_lgal['MetalsBulgeMass'][:,2], \
            G_lgal['MetalsDiskMass'][:,0],G_lgal['MetalsDiskMass'][:,1],G_lgal['MetalsDiskMass'][:,2], \
            G_lgal['MetalsDiskMassRings'][:,:,0],G_lgal['MetalsDiskMassRings'][:,:,1],G_lgal['MetalsDiskMassRings'][:,:,2], \
            G_lgal['ColdGas_elements']]   
    
    ##############
    #Check for NaNs:
    print("\n------------------")
    print("Check for NaNs:")
    print("------------------")
    for ii in range(len(prop)) :     
        nans = np.count_nonzero(np.isnan(prop[ii]))
        if nans == 0 :
            print(prop_name[ii], ":", nans, "/", prop[ii].size)
        else :
            print("*************************************")
            print("WARNING: NaNs in", prop_name[ii], "=", nans, "/", prop[ii].size)
            print("*************************************")
            
    ##############
    #Check for negatives:
    print("\n------------------")
    print("Check for < 0:")
    print("------------------")
    for ii in range(len(prop)) :     
        negs = np.sum(prop[ii] < 0.)
        if negs == 0 :
            print(prop_name[ii], ":", negs, "/", prop[ii].size)
        else :
            var = prop[ii]
            varneg = var[var < 0.]            
            print("*************************************")
            print("WARNING: Negatives in", prop_name[ii], ":", negs, "/", prop[ii].size)
            print("Min =", np.min(varneg), " | ", "Max =", np.max(varneg))
            print("*************************************")
    
    ##############
    #Check for sums:
    print("\n------------------")
    print("Check sub-component masses:")
    print("------------------")
    #Mass sum properties:
    prop1_name = ['ColdGas','DiskMass','MetalsColdGas','MetalsDiskMass','MetalsBulgeMass', \
                  'ColdGas','DiskMass','BulgeMass','HotGas','HaloStars', \
                  'MetalsColdGas','MetalsDiskMass','MetalsBulgeMass', \
                  'MetalsHotGas','MetalsHaloStars', \
                  'Gal number', 'Gal number', 'Gal number', \
                  'Gal number', 'Gal number', 'Gal number', \
                  'Gal number', 'Gal number', 'Gal number']    
    prop1 = [G_lgal['ColdGas'],G_lgal['DiskMass'],np.nansum(G_lgal['MetalsColdGas'], axis=1),np.nansum(G_lgal['MetalsDiskMass'], axis=1),np.nansum(G_lgal['MetalsBulgeMass'], axis=1), \
             G_lgal['ColdGas'],G_lgal['DiskMass'],G_lgal['BulgeMass'],G_lgal['HotGas'],G_lgal['ICM'], \
             np.nansum(G_lgal['MetalsColdGas'], axis=1),np.nansum(G_lgal['MetalsDiskMass'], axis=1),np.nansum(G_lgal['MetalsBulgeMass'], axis=1), \
             np.nansum(G_lgal['MetalsHotGas'], axis=1),np.nansum(G_lgal['MetalsICM'], axis=1), \
             np.arange(len(G_lgal['DiskMass'])), np.arange(len(G_lgal['DiskMass'])), np.arange(len(G_lgal['DiskMass'])), \
             np.arange(len(G_lgal['BulgeMass'])), np.arange(len(G_lgal['BulgeMass'])), np.arange(len(G_lgal['BulgeMass'])), \
             np.arange(len(G_lgal['ICM'])), np.arange(len(G_lgal['ICM'])), np.arange(len(G_lgal['ICM']))]
    prop2_name = ['ColdGasRings','DiskMassRings','MetalsColdGasRings','MetalsDiskMassRings', 'MetalsBulgeMassRings', \
                  'ColdGasElements','DiskMassElements','BulgeMassElements','HotGasElements','HaloStarElements', \
                  'MetalsColdGasElements','MetalsDiskMassElements','MetalsBulgeMassElements', \
                  'MetalsHotGasElements','MetalsHaloStarsElements', \
                  'MetalsDiskMassSNII','MetalsDiskMassSNIa','MetalsDiskMassAGB', \
                  'MetalsBulgeMassSNII','MetalsBulgeMassSNIa','MetalsBulgeMassAGB', \
                  'MetalsHaloStarsSNII','MetalsHaloStarsSNIa','MetalsHaloStarsAGB']
    prop2 = [G_lgal['ColdGasRings'],G_lgal['DiskMassRings'],np.nansum(G_lgal['MetalsColdGasRings'], axis=2),np.nansum(G_lgal['MetalsDiskMassRings'], axis=2),np.nansum(G_lgal['MetalsBulgeMassRings'], axis=2), \
             G_lgal['ColdGas_elements'],G_lgal['DiskMass_elements'],G_lgal['BulgeMass_elements'],G_lgal['HotGas_elements'],G_lgal['ICM_elements'],\
             G_lgal['ColdGas_elements'][:,2:],G_lgal['DiskMass_elements'][:,2:],G_lgal['BulgeMass_elements'][:,2:], \
             G_lgal['HotGas_elements'][:,2:],G_lgal['ICM_elements'][:,2:], \
             G_lgal['MetalsDiskMass'][:,0],G_lgal['MetalsDiskMass'][:,1],G_lgal['MetalsDiskMass'][:,2], \
             G_lgal['MetalsBulgeMass'][:,0],G_lgal['MetalsBulgeMass'][:,1],G_lgal['MetalsBulgeMass'][:,2], \
             G_lgal['MetalsICM'][:,0],G_lgal['MetalsICM'][:,1],G_lgal['MetalsICM'][:,2]]
    
    if plots :
        fig = plt.figure(figsize=(15,38)) #15,25
        plt.rc('font', family='serif', size='14.0')
        rowno = 9 #6
        colno = 3
    
    for ii in range(len(prop1)) :       
        theprop1 = prop1[ii]
        if (prop1_name[ii] != 'Gal number')         :
            theprop2 = np.nansum(prop2[ii], axis=1) 
            if (ii == 7) | (ii == 8) : #If prop1 or prop2 include a multiple of G_lgal['H2fraction']  
                ww = np.where((theprop1 > 100.0) & (theprop2 > 100.0) & (G_lgal['H2fraction'] > 1.e-04))
                theprop3 = G_lgal['H2fraction'][ww]
            else :
                ww = np.where((theprop1 > 100.0) & (theprop2 > 100.0)) #Only consider galaxies with whole-component mass of >100 Msun
            theprop1 = theprop1[ww]
            theprop2 = theprop2[ww]

            if np.allclose(theprop1, theprop2) :
                print("SUM(", prop2_name[ii], ") !=", prop1_name[ii], ": 0 /", prop1[ii].size)#
            else :
                num = np.count_nonzero(theprop1 != theprop2)
                vardiff = theprop1 - theprop2
                vardiff = abs(vardiff[vardiff != 0.0])
                vardifffrac = (theprop1 - theprop2) / theprop1
                vardifffrac = abs(vardifffrac[vardifffrac != 0.0])
                print("*************************************")
                print("WARNING: SUM(", prop2_name[ii], ") !=", prop1_name[ii], ":", num, "/", theprop1.size)
                print("Min diff =", np.min(vardiff), " | ", "Mean diff =", np.mean(vardiff), " | ", "Max diff =", np.max(vardiff))
                print("Min diff frac =", np.min(vardifffrac), " | ", "Mean diff frac =", np.mean(vardifffrac), " | ", "Max diff frac =", np.max(vardifffrac))
                print("*************************************")
        else :
            theprop2 = prop2[ii]
            
        if plots :
            #Set-up plot:
            plt.subplot(rowno, colno, ii+1) 
            if (prop1_name[ii] != 'Gal number') :
                xlimits = np.array([0.0,np.max([np.max(theprop1),np.max(theprop2)])])
                ylimits = np.array([0.0,np.max([np.max(theprop1),np.max(theprop2)])])
                plt.plot(xlimits, ylimits, linestyle='--', color='grey')  
            else :
                xlimits = np.array([np.min(theprop1),np.max(theprop1)])
                ylimits = np.array([np.min(theprop2),np.max(theprop2)])
                plt.plot(xlimits, np.array([0.0,0.0]), linestyle='--', color='grey')  
            #fig = plt.figure(figsize=(5,5))
            plt.xlabel(prop1_name[ii], fontsize=12)
            plt.ylabel(prop2_name[ii], fontsize=12)
            plt.xlim(xlimits)
            plt.ylim(ylimits)
            plt.tick_params(direction="in", top=True, bottom=True, left=True, right=True)
            plt.scatter(theprop1, theprop2, marker='o')
            if ii == 0 :
                plt.text(0.25, 0.92, ['Mass checks', 'z=0', str(len(G_lgal))+' gals', COSMOLOGY, SIMULATION, FILE_TYPE, MODEL, VERSION, SAMPLE_TYPE], transform=fig.transFigure, horizontalalignment='left', color='blue')         
   
    print('Mass checks done')
    