# -*- coding: utf-8 -*-
"""
Created on Thu Nov  4 14:06:32 2021

@author: robyates
"""

"""
main_lgals.py
  ;Main script for analysing L-Galaxies 2020 models. It calls read_lgals_outputs.py
  ;to read-in the L-Galaxies output binaries, then calls make_lgals_samples.py to
  ;select a galaxy sample, then calls various plotting functions to make plots.
  ;
  ;Rob Yates 04-11-2021
  ;
  ;08-11-22: This version was adapted for use at the L-Galaxies workshop 2022
  ;
"""

#Basic packages:
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import astropy
from astropy.io import fits

#Local packages
import sys
sys.path.append('./Robs_python_routines/')
from robs_snapnum_from_redshift import robs_snapnum_from_redshift
from read_lgals_outputs import read_lgals_outputs
from mass_checks import mass_checks
from make_lgals_sample import make_lgals_sample
from make_lgals_dictionary import make_lgals_dictionary
from plot_lgals import *

#Directories:
BaseDir = '../../' #project root directory
OutputDir = BaseDir+'output/' #where the L-Galaxies output files are kept
PlotDir = BaseDir+'figures/' #where plots created by this script are sent




#################
#################
# USER DEFINED INPUTS:
#################
#################
#Switches:
LOAD_SAMPLE = 0 #If on, a pre-made sample of galaxies is loaded from the pickeled save file (quicker). If off, L-Galaxies outputs are read in and a sample is made (slower).
MASS_CHECKS = 0 #If on (and LOAD_SAMPLE is off), key mass properties will be checked for Nans, negatives, and whether sub-components add up to component masses properly.
GENERAL_PLOTS = 1 #If on, general plots are plotted, such as the SMF, cosmic SFRD, etc.
MULTI_REDSHIFT_PLOTS = 0 #OVERRIDES ALL OTHER PLOT SWITCHES ABOVE. If on, plots requiring multiple redshift outputs will be calculated and made (inc. e.g. SFRD evo).

#################
#Output file info:
COSMOLOGY = 'Planck' #Choose from: 'WMAP1', 'Planck'
SIMULATION  = 'Mil-I' #Choose from: 'Mil-I', 'Mil-II' 
FILE_TYPE = 'snapshots' #Choose from: 'snapshots', 'galtree'
STRUCT_TYPE = 'snapshots' #Choose from: 'snapshots', 'galtree'
MODEL = 'default' #Choose from: 'default', 'modified'
VERSION = 'test1' #Pick a memorable name
LABEL = MODEL+' model '+VERSION #NEEDS TO BE IN SIMPLE ASCII (so it can be used in a filename ok in Linux, and read by latex). White space is ok [removed later]). A label that wll be added to plots to denote this model if MULTIPLE_MODELS is on.
    
#################       
#Files to load:
if (SIMULATION == 'Mil-I') :
    FirstFile = 5
    LastFile = 5
elif (SIMULATION == 'Mil-II') :
    FirstFile = 40
    LastFile = 79
else : print("***** ERROR: SIMULATION unknown. Cannot determine FirstFile, LastFile... *****")
TreeFilesUsed = LastFile-FirstFile+1
TotTreeFiles = 512

################# 
#Select sample type:
SAMPLE_TYPE = 'All' #Select from: 'All', 'Discs', 'Dwarfs'

################# 
#Redshifts to load: (in snapshot mode, only works for single redshifts currently 21-04-21)
FullRedshiftList = [0.00] #[0.00,1.04,2.07,3.11,3.95,5.03,5.92,6.97,8.22,8.93]
RedshiftsToRead = [True] #[True,True,True,True,True,True,True,True,True,True]

#################  
#Plot suffix:
mark = 'mk1'




#################
#################
# KEY PARAMETERS:
#################
#################
char_z_low="%0.2f" % FullRedshiftList[0]
char_z_high="%0.2f" % FullRedshiftList[-1]
if (MULTI_REDSHIFT_PLOTS == 1) :    
    FullSnapnumList = robs_snapnum_from_redshift(SIMULATION, COSMOLOGY, FullRedshiftList)

#################
#Treefiles to read:
if SIMULATION == 'Mil-I' :  
    if COSMOLOGY == 'Planck' : snap_z0 = 58
    elif COSMOLOGY == 'WMAP1' : snap_z0 = 63
    else : print("***** ERROR: Cosmology unknown. Please choose from: WMAP1, PLANCK *****") 
elif SIMULATION == 'Mil-II' :
    if COSMOLOGY == 'Planck' : snap_z0 = 62
    elif COSMOLOGY == 'WMAP1' : snap_z0 = 67
    else : print("***** ERROR: Cosmology unknown. Please choose from: WMAP1, PLANCK *****") 
else : print("***** ERROR: Simulation unknown. Please choose from: Mil-I, Mil-II *****")
    
#################
#Set cosmology and mass resolution:
if COSMOLOGY == 'WMAP1' :
    Hubble_h      = 0.73
    Omega_M       = 0.25 
    Omega_Lambda  = 0.75  
elif COSMOLOGY == 'Planck' :
    Hubble_h      = 0.673
    Omega_M       = 0.315 
    Omega_Lambda  = 0.683  
else : print("***** ERROR: Cosmology unknown. Please choose from: WMAP1, PLANCK *****") 

if SIMULATION == 'Mil-I' :
    BoxSideLength = 500.*0.960558 #Mpc/h
    DMParticleMass = 0.0961104*1e10 #Msun/h 
    StellarMassRes = 1.e9 #Rough recommended minimum stellar mass for well-resolved galaxies for Mil-I
elif SIMULATION == 'Mil-II' :
    BoxSideLength = 100.*0.960558 #Mpc/h
    DMParticleMass = 0.000768884*1e10 #Msun/h
    StellarMassRes = 1.e8 #Rough recommended minimum stellar mass for well-resolved galaxies for Mil-II
else : print("***** ERROR: Simulation unknown. Please choose from: Mil-I, Mil-II *****")
ParticleMassRes = 20. #times the particle mass #N.B. DMParticleMass*150. = 1.14e11 Msun/h (for Mil-I) or 1.16e9 Msun/h (for Mil-II), which is the advised minimum virial mass for "well-resolved" subhaloes.
Volume = (BoxSideLength**3.0) * TreeFilesUsed / TotTreeFiles #(Mpc/h)^3 #Note, this volume needs to be divided by h^3 to get the cosmology-dependent value in Mpc^3
    


    
#################
#################
#LOAD SAMPLE:
#################
#################
print('\n***************') 
print('L-GALAXIES 2020')
print('***************')  
print('MODEL:'+COSMOLOGY+"_"+SIMULATION+"_"+FILE_TYPE+"_"+MODEL+"_"+VERSION+"_"+'z'+char_z_low+"-"+char_z_high+"_"+SAMPLE_TYPE+'\n') 
 
if LOAD_SAMPLE == 0 :
    G_lgal, SFH_bins = read_lgals_outputs(BaseDir, OutputDir, Hubble_h, SIMULATION, FILE_TYPE, STRUCT_TYPE, MODEL, VERSION, \
                                          FirstFile, LastFile, FullRedshiftList, RedshiftsToRead)

    if MASS_CHECKS == 1 :
        mass_checks(G_lgal, COSMOLOGY, SIMULATION, FILE_TYPE, STRUCT_TYPE, MODEL, VERSION, SAMPLE_TYPE, plots=1)

    #################                                                            
    #Make and pickle sample data (NOTE: "mark" is not included in the sample filename, so these will be overwritten if the same model is run again):       
    G_samp1 = make_lgals_sample(G_lgal, FILE_TYPE, SAMPLE_TYPE, Hubble_h, DMParticleMass, ParticleMassRes, StellarMassRes, \
                           FullRedshiftList, RedshiftsToRead)
    np.save(OutputDir+COSMOLOGY+"_"+SIMULATION+"_"+FILE_TYPE+"_"+MODEL+"_"+VERSION+"_"+'z'+char_z_low+"-"+char_z_high+"_"+SAMPLE_TYPE, G_samp1, allow_pickle=True)
    NumGals = len(G_samp1) #Number of galaxies in the sample selected in make_lgals_sample (not the total number of objects in the treefiles used).
    print('Sample data pickled\n')

elif LOAD_SAMPLE == 1 :
    G_samp1 = np.load(OutputDir+COSMOLOGY+"_"+SIMULATION+"_"+FILE_TYPE+"_"+MODEL+"_"+VERSION+"_"+'z'+char_z_low+"-"+char_z_high+"_"+SAMPLE_TYPE+".npy")
    NumGals = len(G_samp1)
    print('Pickled sample data loaded\n')
else : print("***** ERROR: No sample calculated or pre-loaded. Please set LOAD_SAMPLE parameter to 0 or 1. *****")


#################     
#Calculate ring info:
RNUM = 12
RBINS = 20 #Max number of SFH bins allowed in L-Galaxies.
RingArray = np.arange(1,RNUM+1)
RingRadii = 0.01*np.power(2.,RingArray)/Hubble_h #in kpc
RingArea = np.empty(RNUM) #in kpc^2
RingCenRadii = np.empty(RNUM)
RReRings = np.empty((NumGals,RNUM))    
RReSFHRings  = np.empty((NumGals,RNUM,RBINS))

#Calculate standard ring areas and central radii:
for ii in range(RNUM) :
    if ii == 0 :
        RingArea[ii] = np.pi * RingRadii[ii]**2
        RingCenRadii[ii] = 0.75*RingRadii[ii]
    else :
        RingArea[ii] = np.pi * (RingRadii[ii]**2 - RingRadii[ii-1]**2)
        RingCenRadii[ii] = RingRadii[ii-1] + ((RingRadii[ii] - RingRadii[ii-1])/2.) 
    RReRings[:,ii] = RingCenRadii[ii]/(1.e+3*G_samp1['StellarHalfLightRadius'][:])
    for jj in range(RBINS) :
        RReSFHRings[:,ii,jj] = RingCenRadii[ii]/(1.e+3*G_samp1['StellarHalfLightRadius'][:])
    
#################
#Calculate SFH properties:
#N.B. If making new samples, SFH_bins IS ALREADY STORED IN THE DICTIONARY BELOW
SFH_bins = astropy.io.fits.open(BaseDir+'AuxCode/Python/'+'Database_SFH_table.fits')
SFH_bins = SFH_bins[1].data    
if FILE_TYPE == 'snapshots' :
    theSnap = G_samp1['SnapNum'][0] #The snapshot number corresponding to this snapshot file
elif FILE_TYPE == 'galtree' :
    theSnap = snap_z0
SFH_bins_lbt = SFH_bins.LOOKBACKTIME[SFH_bins.SNAPNUM == theSnap]/1.e9 #Gyr #Lookback times from z=0 to centre of each z=theSnap SFH bin
SFH_bins_dt = SFH_bins.DT[SFH_bins.SNAPNUM == theSnap] #yr #Width of bins from the z=theSnap SFH
SFH_bins_num = SFH_bins.BIN[SFH_bins.SNAPNUM == theSnap][-1] #Number of bins active in the z=theSnap SFH
SFH_bins_lbt_allGals = np.empty([NumGals,RNUM,RBINS])
SFH_bins_lbt_allGals[:] = np.nan
for ii in range(NumGals) :
    for jj in range(RNUM) :
        SFH_bins_lbt_allGals[ii][jj][0:SFH_bins_num] = SFH_bins_lbt
            
#################
#Make dictionaries for each sample:
Samp1 = make_lgals_dictionary(LABEL, G_samp1, PlotDir, COSMOLOGY, SIMULATION, FILE_TYPE, MODEL, \
                              VERSION, NumGals, SAMPLE_TYPE, char_z_low, char_z_high, \
                              RNUM, RingRadii, RReRings, RReSFHRings, RingArea, \
                              SFH_bins_num, SFH_bins_lbt_allGals)
struct1 = STRUCT_TYPE



        
#################
#################
#MAKE PLOTS:
#################
#################
pdf = PdfPages(PlotDir+Samp1['Cosmology']+"_"+Samp1['Simulation']+"_"+Samp1['File_type']+"_"+Samp1['Model']+"_"+Samp1['Version']+"_"+'z'+char_z_low+"-"+char_z_high+"_"+Samp1['Sample_type']+"_"+mark+".pdf")

MassBins = np.array([[9.0,10.2],[10.2,11.4]]) 
AgeBins = np.array([[0.0,5.0],[5.0,9.0],[9.0,14.0]])

if MULTI_REDSHIFT_PLOTS == 1:
    plot_cosmic_SFRD(Volume, Hubble_h, FullRedshiftList, FullSnapnumList, RedshiftsToRead, Samp1, struct1, char_z_low, pdf=pdf)
else :
    if GENERAL_PLOTS == 1 :
        plot_dists(Samp1, struct1, char_z_low, pdf=pdf)
        plot_smf(Volume, Hubble_h, Samp1, char_z_low, pdf=pdf)
        plot_himf(Volume, Hubble_h, Samp1, char_z_low, pdf=pdf)
        plot_mssfr(Volume, Hubble_h, Samp1, char_z_low, pdf=pdf)
        plot_mzgr(Samp1, struct1, char_z_low, pdf=pdf)
        plot_mzsr(Samp1, char_z_low, pdf=pdf)
        plot_sfhs(Samp1, SFH_bins, snap_z0, char_z_low, pdf=pdf)
        plot_sfrd_prof(Samp1, struct1, char_z_low, MassBins, pdf=pdf)

pdf.close()
print("DONE!")

