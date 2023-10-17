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
#from astropy.io import fits

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
LOAD_SAMPLE = 1 #If on, a pre-made sample of galaxies is loaded from the pickeled save file (quicker). If off, L-Galaxies outputs are read in and a sample is made (slower).
MASS_CHECKS = 0 #If on (and LOAD_SAMPLE is off), key mass properties will be checked for Nans, negatives, and whether sub-components add up to component masses properly.
COMBINE_MODELS = 0 #If on, MR-I and MR-II versions of the same model are combined (if the other version has a .npy sample already saved).
STELLAR_MASS_CUT = 1 #If on, only galaxies above the mass resolution thresholds of log(M*/Msun) >= 8.0 for Millennium-I and log(M*/Msun) >= 7.0 for Millennium-II will be selected.
CALC_SFH_INFO = 0 #If on, SFH info will be calculated and added to the sample dictionaries
GENERAL_PLOTS = 1 #If on, general plots are produced, such as the SMF, MZR, etc.
PAPER_PLOTS = 0 #If on, the plots presented in Yates+23 are produced.
MULTI_REDSHIFT_PLOTS = 0 #If on, in combination with either GENERAL_PLOTS or PAPER_PLOTS, plots requiring multiple redshift outputs will be calculated and made.

#################
#Sample info:
COSMOLOGY = 'Planck' #Choose from: 'WMAP1', 'Planck'
SIMULATION  = 'Mil-I' #Choose from: 'Mil-I', 'Mil-II' 
FILE_TYPE = 'snapshots' #Choose from: 'snapshots', 'galtree'
STRUCT_TYPE = 'liteOutput' #Choose from: 'normal', 'liteOutput', 'ringSFHs'
MODEL = 'modified' #Choose from: 'default', 'modified'
VERSION = 'test1' #Pick a memorable name
LABEL = MODEL+' model '+VERSION #NEEDS TO BE IN SIMPLE ASCII (so it can be used in a filename ok in Linux, and read by latex). White space is ok [removed later]). A label that wll be added to plots to denote this model if MULTIPLE_MODELS is on.

SOLAR_ABUNDANCE_SET = 'A09' #'GAS07' #'AG89_phot' #'AG89_mete' #Sets which solar abundances are assumed when normalising abundances and enhancements inplots
SAMPLE_TYPE = 'All' #Select from: 'All', 'Discs', 'ETGs', 'Dwarfs'

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
FullSnapnumList_MRI = robs_snapnum_from_redshift('Mil-I', COSMOLOGY, FullRedshiftList)
FullSnapnumList_MRII = robs_snapnum_from_redshift('Mil-II', COSMOLOGY, FullRedshiftList)
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

MRI_cutoff = 1.e9 #Msun
MRI_BoxSideLength = 500.*0.960558 #Mpc/h
MRII_BoxSideLength = 100.*0.960558 #Mpc/h
if SIMULATION == 'Mil-I' :
    BoxSideLength = MRI_BoxSideLength #Mpc/h
    DMParticleMass = 0.0961104*1e10 #Msun/h 
    if STELLAR_MASS_CUT == 1 :
        StellarMassRes = 1.e8 #Rough recommended minimum stellar mass for well-resolved galaxies for Mil-I
    else :
        StellarMassRes = 0.0
elif SIMULATION == 'Mil-II' :
    BoxSideLength = MRII_BoxSideLength #Mpc/h
    DMParticleMass = 0.000768884*1e10 #Msun/h
    if STELLAR_MASS_CUT == 1 :
        StellarMassRes = 1.e7 #Rough recommended minimum stellar mass for well-resolved galaxies for Mil-II
    else :
        StellarMassRes = 0.0
else : print("***** ERROR: Simulation unknown. Please choose from: Mil-I, Mil-II *****")
ParticleMassRes = 20. #times the particle mass #N.B. DMParticleMass*150. = 1.14e11 Msun/h (for Mil-I) or 1.16e9 Msun/h (for Mil-II), which is the advised minimum virial mass for "well-resolved" subhaloes.
Volume = (BoxSideLength**3.0) * TreeFilesUsed / TotTreeFiles #(Mpc/h)^3 #Note, this volume needs to be divided by h^3 to get the cosmology-dependent value in Mpc^3
    
    
#################
#################
#LOAD SAMPLE:
#################
#################
print('\n**************************') 
print('L-GALAXIES 2020 (Yates+23)')
print('**************************')  
print('MODEL:'+COSMOLOGY+"_"+SIMULATION+"_"+FILE_TYPE+"_"+MODEL+"_"+VERSION+"_"+'z'+char_z_low+"-"+char_z_high+"_"+SAMPLE_TYPE+'\n')  

if LOAD_SAMPLE == 0 :
    G_lgal, SFH_bins = read_lgals_outputs(BaseDir, OutputDir, Hubble_h, SIMULATION, FILE_TYPE, STRUCT_TYPE, MODEL, VERSION, \
                                          FirstFile, LastFile, FullRedshiftList, RedshiftsToRead)

    if MASS_CHECKS == 1 :
        mass_checks(G_lgal, COSMOLOGY, SIMULATION, FILE_TYPE, STRUCT_TYPE, MODEL, VERSION, SAMPLE_TYPE, plots=1)

    #################                                                            
    #Make and pickle sample data (NOTE: "mark" is not included in the sample filename, so these will be overwritten if the same model is run again):          
    G_samp1 = make_lgals_sample(G_lgal, FILE_TYPE, SAMPLE_TYPE, SIMULATION, COSMOLOGY, Hubble_h, DMParticleMass, ParticleMassRes, StellarMassRes, \
                           FullRedshiftList, RedshiftsToRead)
    np.save(OutputDir+'samples/'+COSMOLOGY+"_"+SIMULATION+"_"+FILE_TYPE+"_"+MODEL+"_"+VERSION+"_"+'z'+char_z_low+"-"+char_z_high+"_"+SAMPLE_TYPE, G_samp1, allow_pickle=True)
    NumGals = len(G_samp1) #Number of galaxies in the sample selected in make_lgals_sample (not the total number of objects in the treefiles used).
    #Save treefile info:
    np.save(OutputDir+'samples/'+COSMOLOGY+"_"+SIMULATION+"_"+FILE_TYPE+"_"+MODEL+"_"+VERSION+"_"+'z'+char_z_low+"-"+char_z_high+"_"+SAMPLE_TYPE+"_"+"treefiles", [FirstFile,LastFile], allow_pickle=True)
    print('Sample data pickled\n')

elif LOAD_SAMPLE == 1 :
    G_samp1 = np.load(OutputDir+"/samples/"+COSMOLOGY+"_"+SIMULATION+"_"+FILE_TYPE+"_"+MODEL+"_"+VERSION+"_"+'z'+char_z_low+"-"+char_z_high+"_"+SAMPLE_TYPE+".npy")
    NumGals = len(G_samp1)
    print('Pickled sample data loaded\n')
    
    #Check treefile info:
    treefiles = np.load(OutputDir+'samples/'+COSMOLOGY+"_"+SIMULATION+"_"+FILE_TYPE+"_"+MODEL+"_"+VERSION+"_"+'z'+char_z_low+"-"+char_z_high+"_"+SAMPLE_TYPE+"_"+"treefiles"+".npy")
    if ((treefiles[0] != FirstFile) | (treefiles[1] != LastFile)) :
        errString = "***** ERROR: The treefile range "+str(treefiles)+" from the loaded sample doesn't match the FirstFile ["+str(FirstFile)+"] and LastFile ["+str(LastFile)+"] requested. *****"
        sys.exit(errString)
else : print("***** ERROR: No sample calculated or pre-loaded. Please set LOAD_SAMPLE parameter to 0 or 1. *****")



#################
#################
#COMBINE MIL-I/II:
#################
#################
if SIMULATION == 'Mil-I' :
    MRI_gals = NumGals
    MRII_gals = 0
    Volume_MRI = Volume #(Mpc/h)^3 #Note, this volume needs to be divided by h^3 to get the cosmology-dependent value in Mpc^3
    Volume_MRII = 0.0
elif SIMULATION == 'Mil-II' :
    MRI_gals = 0
    MRII_gals = NumGals
    Volume_MRI = 0.0
    Volume_MRII = Volume #(Mpc/h)^3 #Note, this volume needs to be divided by h^3 to get the cosmology-dependent value in Mpc^3

#Combine MRI and MRII data, with MRI galaxies always first:
if COMBINE_MODELS == 1 :
    if SIMULATION == 'Mil-I' :
        G_sampb = np.load(OutputDir+'samples/'+COSMOLOGY+"_"+"Mil-II"+"_"+FILE_TYPE+"_"+MODEL+"_"+"MRII_"+VERSION+"_"+'z'+char_z_low+"-"+char_z_high+"_"+SAMPLE_TYPE+".npy")
        treefilesb = np.load(OutputDir+'samples/'+COSMOLOGY+"_"+"Mil-II"+"_"+FILE_TYPE+"_"+MODEL+"_"+"MRII_"+VERSION+"_"+'z'+char_z_low+"-"+char_z_high+"_"+SAMPLE_TYPE+"_"+"treefiles"+".npy")
        MRII_gals = len(G_sampb)
        Volume_MRII = (MRII_BoxSideLength**3.0) * (treefilesb[1]-treefilesb[0]+1) / TotTreeFiles #(Mpc/h)^3 #Note, this volume needs to be divided by h^3 to get the cosmology-dependent value in Mpc^3
        G_samp1 = np.hstack((G_samp1,G_sampb))
        NumGals = len(G_samp1)
    elif SIMULATION == "Mil-II" :   
        G_sampb = np.load(OutputDir+'samples/'+COSMOLOGY+"_"+"Mil-I"+"_"+FILE_TYPE+"_"+MODEL+"_"+VERSION.replace("MRII_","")+"_"+'z'+char_z_low+"-"+char_z_high+"_"+SAMPLE_TYPE+".npy")
        treefilesb = np.load(OutputDir+'samples/'+COSMOLOGY+"_"+"Mil-I"+"_"+FILE_TYPE+"_"+MODEL+"_"+VERSION.replace("MRII_","")+"_"+'z'+char_z_low+"-"+char_z_high+"_"+SAMPLE_TYPE+"_"+"treefiles"+".npy")
        MRI_gals = len(G_sampb)
        Volume_MRI = (MRI_BoxSideLength**3.0) * (treefilesb[1]-treefilesb[0]+1) / TotTreeFiles #(Mpc/h)^3 #Note, this volume needs to be divided by h^3 to get the cosmology-dependent value in Mpc^3
        G_samp1 = np.hstack((G_sampb,G_samp1))
        NumGals = len(G_samp1)
    NumGalsb = len(G_sampb)
    print("Mil-I and Mil-II models combined.")

print("NumGals = ", NumGals)
print("MRI_gals = ", MRI_gals)
print("MRII_gals = ", MRII_gals, "\n")


#################
#################
#CALCULATE RING AND SFH INFO:
#################
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
if ("liteOutput" in (STRUCT_TYPE)) :
    print("NOTE: Using StellarHalfMassRadius as effective radius (Re)\n")
else :
    print("NOTE: Using StellarHalfLightRadius as effective radius (Re)\n")
for ii in range(RNUM) :
    if ii == 0 :
        RingArea[ii] = np.pi * RingRadii[ii]**2
        RingCenRadii[ii] = 0.75*RingRadii[ii]
    else :
        RingArea[ii] = np.pi * (RingRadii[ii]**2 - RingRadii[ii-1]**2)
        RingCenRadii[ii] = RingRadii[ii-1] + ((RingRadii[ii] - RingRadii[ii-1])/2.)      
#Calculate R/Re for every ring of every galaxy:
    if ("liteOutput" in (STRUCT_TYPE)) :
        RReRings[:,ii] = RingCenRadii[ii]/(1.e+3*G_samp1['StellarHalfMassRadius'][:])
        for jj in range(RBINS) :
            RReSFHRings[:,ii,jj] = RingCenRadii[ii]/(1.e+3*G_samp1['StellarHalfMassRadius'][:])      
    else :           
        RReRings[:,ii] = RingCenRadii[ii]/(1.e+3*G_samp1['StellarHalfLightRadius'][:])
        for jj in range(RBINS) :
            RReSFHRings[:,ii,jj] = RingCenRadii[ii]/(1.e+3*G_samp1['StellarHalfLightRadius'][:])
                
if CALC_SFH_INFO == 1 :        
    #Calculate SFH properties:
    #N.B. If making new samples, SFH_bins IS ALREADY STORED IN THE DICTIONARY BELOW
    SFH_bins = astropy.io.fits.open(BaseDir+'AuxCode/Python/'+'Database_SFH_table.fits')
    SFH_bins = SFH_bins[1].data    
    if FILE_TYPE == 'snapshots' :
        theSnap = G_samp1['SnapNum'][0] #The snapshot number corresponding to this snapshot file
    elif FILE_TYPE == 'galtree' :
        theSnap = snap_z0
    SFH_bins_lbt = SFH_bins.LOOKBACKTIME[SFH_bins.SNAPNUM == theSnap]/1.e9 #Gyr #Lookback times from z=theSnap to centre of each z=theSnap SFH bin
    SFH_bins_dt = SFH_bins.DT[SFH_bins.SNAPNUM == theSnap] #yr #Width of bins from the z=theSnap SFH
    SFH_bins_num = SFH_bins.BIN[SFH_bins.SNAPNUM == theSnap][-1] #Number of bins active in the z=theSnap SFH
    SFH_bins_lbt_allGals = np.empty([NumGals,RNUM,RBINS])
    SFH_bins_lbt_allGals[:] = np.nan
    for ii in range(NumGals) :
        for jj in range(RNUM) :
            SFH_bins_lbt_allGals[ii][jj][0:SFH_bins_num] = SFH_bins_lbt
          
            
#################
#################
#MAKE SAMPLE DICTIONARY:
#################
#################
if CALC_SFH_INFO == 1 :
    Samp1 = make_lgals_dictionary(LABEL, G_samp1, PlotDir, COSMOLOGY, SIMULATION, FILE_TYPE, MODEL, \
                                  VERSION, MRI_cutoff, NumGals, MRI_gals, MRII_gals, SAMPLE_TYPE, char_z_low, char_z_high, \
                                  Volume_MRI, Volume_MRII, FullSnapnumList_MRI, FullSnapnumList_MRII, CALC_SFH_INFO, \
                                  RNUM, RingRadii, RingCenRadii, RReRings, RReSFHRings, RingArea, SFH_bins_num, SFH_bins_lbt_allGals)
else :
    Samp1 = make_lgals_dictionary(LABEL, G_samp1, PlotDir, COSMOLOGY, SIMULATION, FILE_TYPE, MODEL, \
                                  VERSION, MRI_cutoff, NumGals, MRI_gals, MRII_gals, SAMPLE_TYPE, char_z_low, char_z_high, \
                                  Volume_MRI, Volume_MRII, FullSnapnumList_MRI, FullSnapnumList_MRII, CALC_SFH_INFO, \
                                  RNUM, RingRadii, RingCenRadii, RReRings, RReSFHRings, RingArea)
struct1 = STRUCT_TYPE

        
#################
#################
#MAKE PLOTS:
#################
#################
if GENERAL_PLOTS == 1 :
    pdf = PdfPages(PlotDir+Samp1['Cosmology']+"_"+Samp1['Simulation']+"_"+Samp1['File_type']+"_" \
                   +Samp1['Model']+"_"+Samp1['Version']+"_"+'z'+char_z_low+"-"+char_z_high+"_"+Samp1['Sample_type']+"_"+mark+".pdf") 
    MassBins = np.array([[9.0,10.2],[10.2,11.4]])
    AgeBins = np.array([[0.0,5.0],[5.0,9.0],[9.0,14.0]])
    if MULTI_REDSHIFT_PLOTS == 1:
        plot_cosmic_SFRD(Volume, Hubble_h, FullRedshiftList, FullSnapnumList, RedshiftsToRead, Samp1, struct1, char_z_low, pdf=pdf)
    else :
        plot_dists(Samp1, struct1, char_z_low, pdf=pdf)
        plot_smf(Volume, Hubble_h, Samp1, char_z_low, pdf=pdf)
        plot_himf(Volume, Hubble_h, Samp1, char_z_low, pdf=pdf)
        plot_mssfr(Volume, Hubble_h, Samp1, char_z_low, pdf=pdf)
        plot_mzgr(Samp1, struct1, char_z_low, pdf=pdf)
        plot_mzsr(Samp1, char_z_low, pdf=pdf)
        if CALC_RINGS_AND_SFH_INFO == 1 :
            plot_sfrd_prof(Samp1, struct1, char_z_low, MassBins, pdf=pdf)
            if not ("liteOutput" in (STRUCT_TYPE)) :
                plot_sfhs(Samp1, SFH_bins, snap_z0, char_z_low, pdf=pdf)
    pdf.close()
        
if PAPER_PLOTS == 1 :
    MassBins = np.array([[9.0,9.75],[9.75,10.5],[10.5,11.25]])
    AgeBins = np.array([[0.0,5.0],[5.0,9.0],[9.0,14.0]])
    if MULTI_REDSHIFT_PLOTS == 1 :
        plot_dust_scaling_relations_evo(Hubble_h, FullRedshiftList, FullSnapnumList_MRI, FullSnapnumList_MRII, RedshiftsToRead, Samp1, struct1, \
                                        props=['DTM','DTG'], xprop='OH', SolarNorm=SOLAR_ABUNDANCE_SET, aveOnly=1, incHotDust=None, SFRWeighted=None)
        plot_dust_scaling_relations_evo(Hubble_h, FullRedshiftList, FullSnapnumList_MRI, FullSnapnumList_MRII, RedshiftsToRead, Samp1, struct1, \
                                        props=['Mdust'], xprop='stellarMass', SolarNorm=SOLAR_ABUNDANCE_SET, aveOnly=1, incHotDust=None, SFRWeighted=None)
        if COMBINE_MODELS != 1 :
            plot_cosmic_dust_evos(Volume, Hubble_h, Omega_M, Omega_b, FullRedshiftList, FullSnapnumList, RedshiftsToRead, \
                                  char_z_low, char_z_high, Samp1, obs=1, plotTotal=None)
    else :
        plot_timescales(Samp1, struct1, REDSHIFT, xprop='OH', yprops='All', contourLines='graded', \
                        outlierFrac=0.0, SolarNorm=SOLAR_ABUNDANCE_SET) #xprop: Choose from: ['OH','logZ','H2D']
        plot_smf(Hubble_h, Samp1, REDSHIFT, Add_Edd_bias=None)
        plot_himf(Hubble_h, Samp1, struct1, REDSHIFT, prop="H")
        plot_mssfr(Hubble_h, Samp1, REDSHIFT, contourLines='graded')
        plot_snrates(Samp1, REDSHIFT)
        plot_mzgr(Hubble_h, Samp1, struct1, REDSHIFT, dustCorrec=1, contourLines='graded')
        plot_mzsr(Hubble_h, Omega_M, Omega_Lambda, Samp1, REDSHIFT, ApertureCorrec=1, contourLines='graded')
        plot_dust_scaling_relations(Hubble_h, Samp1, struct1, REDSHIFT, props='All', xprop='OH', \
                                    SolarNorm=SOLAR_ABUNDANCE_SET, levels='3sig', contourOnly=1, \
                                    contourLines="graded", SFRWeighted=1, incHotDust=None) #xprop: Choose from: ['OH','MH','stellarMass']
print("DONE!")

