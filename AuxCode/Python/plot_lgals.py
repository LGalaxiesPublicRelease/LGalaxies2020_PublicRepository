# -*- coding: utf-8 -*-
"""
Created on Thu Nov  4 17:36:17 2021

@author: ry0005
"""

"""
plot_lgals.py
  ;Functions to make plots for L-Galaxies
  ;
  ;Rob Yates 04-11-2021
  ;
  ;
"""

#Basic packages:
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.collections as mcoll
import astropy
from astropy.io import fits

#Local packages:
import sys
sys.path.append('./Robs_python_routines/')
from chemistry import *
from robs_contour_plot import robs_contour_plot
from robs_plot_average import robs_plot_average
from robs_plot_profile import robs_plot_profile
from robs_solarnorm_finder import robs_solarnorm_finder
from robs_element_list_finder import robs_element_list_finder
from robs_plot_text import robs_plot_text
from robs_plot_panels import robs_plot_panels


#################
#Default parameters:
plt.rcParams['text.usetex'] = True
plt.rcParams['font.size'] = 14
plt.rcParams['font.serif'] = 'Times'
defaultOutlierFrac = 0.25


#################
#Directories:
BaseDir = '../../' #project root directory
OutputDir = BaseDir+'output/' #where the L-Galaxies output files are kept
PlotDir = BaseDir+'figures/' #where plots created by this script are sent


#################
def plot_dists(Samp1, struct1, char_z_low, pdf=None) :  
    #Set-up plot:
    fig = plt.figure(figsize=(12,35))
    plt.rc('font', family='serif', size='14.0')
    binno = 40
    rowno = 6 #9
    colno = 2
    
    p1 = plt.subplot(rowno, colno, 1)
    plt.xlabel(r'log($M_{\rm vir}$/M$_{\odot}$)')
    plt.ylabel(r'Number')
    xlimits = np.array([10.0,14.5])
    plt.xlim(np.array(xlimits))
    plt.tick_params(direction="in", top=True, bottom=True, left=True, right=True)
    plt.hist(np.log10(Samp1['G_samp']['Mvir']), bins=binno, range=xlimits) #, density=False
 
    p2 = plt.subplot(rowno, colno, 2)
    plt.xlabel(r'log($M_{*}$/M$_{\odot}$)')
    plt.ylabel(r'Number')
    xlimits = np.array([9.0,12.0])
    plt.xlim(np.array(xlimits))
    plt.tick_params(direction="in", top=True, bottom=True, left=True, right=True)
    plt.hist(np.log10(Samp1['G_samp']['StellarMass']), bins=binno, range=xlimits)

    p3 = plt.subplot(rowno, colno, 3)
    plt.xlabel(r'Age$_{\rm mw}$/Gyr')
    plt.ylabel(r'Number')
    xlimits = np.array([0.0,14.0])
    plt.xlim(np.array(xlimits))
    plt.tick_params(direction="in", top=True, bottom=True, left=True, right=True)
    plt.hist(Samp1['G_samp']['MassWeightAge'], bins=binno, range=xlimits)

    p4 = plt.subplot(rowno, colno, 4)
    plt.xlabel(r'log(sSFR/yr)')
    plt.ylabel(r'Number')
    xlimits = np.array([-16.0,-8.0])
    plt.xlim(np.array(xlimits))
    plt.tick_params(direction="in", top=True, bottom=True, left=True, right=True)
    plt.hist(np.log10(Samp1['G_samp']['Sfr'][Samp1['G_samp']['Sfr']>0.0]/Samp1['G_samp']['StellarMass'][Samp1['G_samp']['Sfr']>0.0]), bins=binno, range=xlimits) #, range=np.array([-17.,-8.])

    p5 = plt.subplot(rowno, colno, 5)
    plt.xlabel(r'log(M$_{\rm HI}$/M$_{\odot}$)')
    plt.ylabel(r'Number')
    xlimits = np.array([4.0,12.0])
    plt.xlim(np.array(xlimits))
    plt.tick_params(direction="in", top=True, bottom=True, left=True, right=True)
    plt.hist(np.log10(Samp1['G_samp']['ColdGas'][Samp1['G_samp']['ColdGas']>0.0]*(1.-Samp1['G_samp']['H2fraction'][Samp1['G_samp']['ColdGas']>0.0])), bins=binno, range=xlimits) #, range=np.array([4.,10.])   

    p6 = plt.subplot(rowno, colno, 6)
    plt.xlabel(r'log(M$_{\rm H2}$/M$_{\odot}$)')
    plt.ylabel(r'Number')
    xlimits = np.array([1.0,11.0])
    plt.xlim(np.array(xlimits))
    plt.tick_params(direction="in", top=True, bottom=True, left=True, right=True)
    plt.hist(np.log10(Samp1['G_samp']['ColdGas'][Samp1['G_samp']['H2fraction']*Samp1['G_samp']['ColdGas']>0.0]*Samp1['G_samp']['H2fraction'][Samp1['G_samp']['H2fraction']*Samp1['G_samp']['ColdGas']>0.0]), bins=binno, range=xlimits) #, range=np.array([4.,10.])

    p7 = plt.subplot(rowno, colno, 7)
    plt.xlabel(r'12+log(O/H)')
    plt.ylabel(r'Number')
    xlimits = np.array([5.5,10.5])
    plt.xlim(np.array(xlimits))
    plt.tick_params(direction="in", top=True, bottom=True, left=True, right=True)
    Z_OH = 12. + np.log10((Samp1['G_samp']['ColdGas_elements'][:,4] / Samp1['G_samp']['ColdGas_elements'][:,0]) * (H_aw/O_aw))
    plt.hist(Z_OH[Z_OH > 0.0], bins=binno, range=xlimits)

    p8 = plt.subplot(rowno, colno, 8)
    plt.xlabel(r'[Fe/H]$_{\rm bulge+disc}$')
    plt.ylabel(r'Number')
    xlimits = np.array([-3.0,0.5])
    plt.xlim(np.array(xlimits))
    plt.tick_params(direction="in", top=True, bottom=True, left=True, right=True)
    FeH = np.log10(((Samp1['G_samp']['DiskMass_elements'][:,10]+Samp1['G_samp']['BulgeMass_elements'][:,10]) / (Samp1['G_samp']['DiskMass_elements'][:,0]+Samp1['G_samp']['BulgeMass_elements'][:,0]))) - np.log10(FeH_mf_A09)
    plt.hist(FeH[FeH > -10.0], bins=binno, range=xlimits)
  
    p9 = plt.subplot(rowno, colno, 9)
    plt.xlabel(r'log(M$_{\rm O,gas}$/M$_{\odot}$)')
    plt.ylabel(r'Number')
    xlimits = np.array([2.0,10.0])
    plt.xlim(np.array(xlimits))
    plt.tick_params(direction="in", top=True, bottom=True, left=True, right=True)
    Omass = np.log10(Samp1['G_samp']['ColdGas_elements'][:,4])
    plt.hist(Omass[Omass > 0.0], bins=binno, range=xlimits)

    p10 = plt.subplot(rowno, colno, 10)
    plt.xlabel(r'log(M$_{\rm H, gas}$/M$_{\odot}$)')
    plt.ylabel(r'Number')
    xlimits = np.array([4.0,12.0])
    plt.xlim(np.array(xlimits))
    plt.tick_params(direction="in", top=True, bottom=True, left=True, right=True)
    Hmass = np.log10(Samp1['G_samp']['ColdGas_elements'][:,0])
    plt.hist(Hmass[Hmass > 0.0], bins=binno, range=xlimits, label=Samp1['Label'])

    p11 = plt.subplot(rowno, colno, 11)
    plt.xlabel(r'log$(R_{\rm SNII} / 100\textnormal{ yr})$')
    plt.ylabel(r'Number')
    SNIIRate = Samp1['G_samp']['DiskSNIIRate'] + Samp1['G_samp']['BulgeSNIIRate'] + Samp1['G_samp']['ICMSNIIRate'] #[1/yr]
    logSNIIRate_per_cen = np.log10(SNIIRate * 100.)
    xlimits = [-7.0,3.0]
    plt.xlim(np.array(xlimits))
    plt.tick_params(direction="in", top=True, bottom=True, left=True, right=True)
    plt.hist(logSNIIRate_per_cen[SNIIRate > 0.0], bins=binno, range=xlimits, label=Samp1['Label'])

    #title:
    plt.text(.5, 0.91, ['Prop dists', 'z='+char_z_low, str(Samp1['NumGals'])+' gals', Samp1['Cosmology'], Samp1['Simulation'], Samp1['File_type'], Samp1['Model'], Samp1['Version'], Samp1['Sample_type']], transform=fig.transFigure, horizontalalignment='center', color='blue')  
    
    #print figure to pdf:
    pdf.savefig()  
    print("Property distributions plotted")




#################
def plot_smf(Volume, Hubble_h, Samp1, char_z_low, pdf=None) :  
    #Set axis limits:
    xlimits = np.array([9.0,12.0])
    ylimits = np.array([-5.9,-0.5])
    
    #Set-up plot:
    fig = plt.figure(figsize=(6,5))
    plt.xlabel(r'log$(M_{*} / \textnormal{M}_{\odot})$')
    plt.ylabel(r'log$(\phi / \textnormal{Mpc}^{-3} \textnormal{ dex}^{-1})$')
    plt.xlim(xlimits)
    plt.ylim(ylimits)
    plt.tick_params(direction="in", top=True, bottom=True, left=True, right=True)
    
    #Plot stellar mass function:
    if ((Samp1["MRI_gals"] > 0) & (Samp1["MRII_gals"] > 0)) :
        prop_MRI = np.log10(Samp1['G_samp']['StellarMass'][0:Samp1['MRI_gals']])
        prop_MRII = np.log10(Samp1['G_samp']['StellarMass'][Samp1['MRI_gals']:])
        bin_arr_MRI = np.arange(np.log10(Samp1["MRI_cutoff"]),xlimits[1]+binwidth,binwidth)
        bin_arr_MRII = np.arange(xlimits[0],np.log10(Samp1["MRI_cutoff"])+(2.*binwidth),binwidth)
        hist_MRI = np.histogram(prop_MRI, bins=bin_arr_MRI, range=(np.log10(Samp1["MRI_cutoff"]),xlimits[1])) 
        hist_MRII = np.histogram(prop_MRII, bins=bin_arr_MRII, range=(xlimits[0],np.log10(Samp1["MRI_cutoff"]))) 
        plt.plot(hist_MRI[1][0:len(hist_MRI[1][:])-1]+binwidth/2., np.log10((hist_MRI[0][:]/(Samp1["Volumes"][0]*binwidth))*Hubble_h**3), label=Samp1['Label'])
        plt.plot(hist_MRII[1][0:len(hist_MRII[1][:])-1]+binwidth/2., np.log10((hist_MRII[0][:]/(Samp1["Volumes"][Samp1['MRI_gals']]*binwidth))*Hubble_h**3), linestyle='dotted')
    else :
        prop = np.log10(Samp1['G_samp']['StellarMass'])
        hist = np.histogram(prop, bins=bin_arr, range=(xlimits[0],xlimits[1])) 
        plt.plot(hist[1][0:len(hist[1][:])-1]+binwidth/2., np.log10((hist[0][:]/(Samp1['Volumes'][0]*binwidth))*Hubble_h**3), label=Samp1['Label'])

    ##########
    #title:
    plt.text(.5, 0.91, ['SMF', 'z='+char_z_low, str(Samp1['NumGals'])+' gals', Samp1['Cosmology'], Samp1['Simulation'], Samp1['File_type'], Samp1['Model'], Samp1['Version'], Samp1['Sample_type']], transform=fig.transFigure, horizontalalignment='center', color='blue')  
    
    #print figure to pdf:
    pdf.savefig()    
    print("SMF plotted")




#################
def plot_himf(Volume, Hubble_h, Samp1, char_z_low, pdf=None) :  
    #Set axis limits:
    xlimits = np.array([7.0,11.0])
    ylimits = np.array([-5.9,-0.5])
    
    #Set-up plot:
    fig = plt.figure(figsize=(6,5))
    plt.xlabel(r'log$(M_{\textnormal{HI}} / \textnormal{M}_{\odot})$')
    plt.ylabel(r'log$(\phi / \textnormal{Mpc}^{-3} \textnormal{ dex}^{-1})$')
    plt.xlim(xlimits)
    plt.ylim(ylimits)
    plt.tick_params(direction="in", top=True, bottom=True, left=True, right=True)
    
    #Plot HI mass function:
    binwidth = 0.1
    bin_arr=np.arange(xlimits[0],xlimits[1]+binwidth,binwidth)
    if ((Samp1["MRI_gals"] > 0) & (Samp1["MRII_gals"] > 0)) :
        MRI_HI_cutoff = 10**(9.5) #Samp1["MRI_cutoff"]
        prop_MRI = np.log10(Samp1['G_samp']['ColdGas'][0:Samp1['MRI_gals']]*(1.-Samp1['G_samp']['H2fraction'][0:Samp1['MRI_gals']]))
        prop_MRII = np.log10(Samp1['G_samp']['ColdGas'][Samp1['MRI_gals']:]*(1.-Samp1['G_samp']['H2fraction'][Samp1['MRI_gals']:]))
        bin_arr_MRI = np.arange(np.log10(MRI_HI_cutoff),xlimits[1]+binwidth,binwidth)
        bin_arr_MRII = np.arange(xlimits[0],np.log10(MRI_HI_cutoff)+(2.*binwidth),binwidth)
        hist_MRI = np.histogram(prop_MRI, bins=bin_arr_MRI, range=(np.log10(MRI_HI_cutoff),xlimits[1])) 
        hist_MRII = np.histogram(prop_MRII, bins=bin_arr_MRII, range=(xlimits[0],np.log10(MRI_HI_cutoff))) 
        plt.plot(hist_MRI[1][0:len(hist_MRI[1][:])-1]+binwidth/2., np.log10((hist_MRI[0][:]/(Samp1["Volumes"][0]*binwidth))*Hubble_h**3), label=Samp1['Label']) 
        plt.plot(hist_MRII[1][0:len(hist_MRII[1][:])-1]+binwidth/2., np.log10((hist_MRII[0][:]/(Samp1["Volumes"][Samp1['MRI_gals']]*binwidth))*Hubble_h**3), linestyle='dotted') 
    else :
        prop = np.log10(Samp1['G_samp']['ColdGas']*(1.-Samp1['G_samp']['H2fraction']))
        hist = np.histogram(prop, bins=bin_arr, range=(xlimits[0],xlimits[1])) 
        plt.plot(hist[1][0:len(hist[1][:])-1]+binwidth/2., np.log10((hist[0][:]/(Samp1['Volumes'][0]*binwidth))*Hubble_h**3), label=Samp1['Label']) 

    ##########
    #title:
    plt.text(.5, 0.91, ['HIMF', 'z='+char_z_low, str(Samp1['NumGals'])+' gals', Samp1['Cosmology'], Samp1['Simulation'], Samp1['File_type'], Samp1['Model'], Samp1['Version'], Samp1['Sample_type']], transform=fig.transFigure, horizontalalignment='center', color='blue')  
        
    #print figure to pdf:
    pdf.savefig()    
    print("HIMF plotted")
    
    
    

#################
def plot_mssfr(Volume, Hubble_h, Samp1, char_z_low, pdf=None) :  
    #Set axis limits:
    xlimits = np.array([9.0,12.0])
    ylimits = np.array([-14.0,-7.5])
    
    #Set-up plot:
    fig = plt.figure(figsize=(6,5))
    plt.xlabel(r'log$(M_{*} / \textnormal{M}_{\odot})$')
    plt.ylabel(r'log$(\textnormal{sSFR} / \textnormal{yr}^{-1})$')
    plt.xlim(xlimits)
    plt.ylim(ylimits)
    plt.tick_params(direction="in", top=True, bottom=True, left=True, right=True)
    
    #Plot M* - sSFR relation:
    binno = 50
    theX = np.log10(Samp1['G_samp']['StellarMass'][Samp1['G_samp']['Sfr']>0.0])
    theY = np.log10(Samp1['G_samp']['Sfr'][Samp1['G_samp']['Sfr']>0.0]/Samp1['G_samp']['StellarMass'][Samp1['G_samp']['Sfr']>0.0])
    #robs_contour_plot(theX, theY, binno, '2sig', noOutliers=1)
    robs_contour_plot(theX, theY, binno, '2sig', noOutliers=1, theRange=[xlimits,ylimits], \
                      weights=1./Samp1['Volumes'][Samp1['G_samp']['Sfr']>0.0])
    #robs_plot_average(theX, theY, binno=10, aveType='Median', minInBin=50, linewidth=4., inputCentres=1)

    ##########
    #title:
    plt.text(.5, 0.92, ['M*-sSFR', 'z='+char_z_low, str(Samp1['NumGals'])+' gals', Samp1['Cosmology'], Samp1['Simulation'], Samp1['File_type'], Samp1['Model'], Samp1['Version'], Samp1['Sample_type']], transform=fig.transFigure, horizontalalignment='center', color='blue')          

    #print figure to pdf:
    pdf.savefig()    
    print("M*-sSFR plotted")
    
    
    
    
#################
def plot_mzgr(Samp1, struct1, char_z_low, pdf=None) :  
    ##########    
    # Plot ColdGas metallicity as mass-weighted M_Z/M_tot, using metals arrays, normalised to Sun:
    ##########     
    #Set axis limits:
    xlimits = np.array([9.0,12.0])
    ylimits = np.array([-1.5,0.5])
    
    #Set-up plot:
    fig = plt.figure(figsize=(6,5))
    plt.xlabel(r'log$(M_{*} / \textnormal{M}_{\odot})$')
    plt.ylabel(r'$Z_{\rm gas} / \textnormal{Z}_{\odot}$')
    plt.xlim(xlimits)
    plt.ylim(ylimits)
    plt.tick_params(direction="in", top=True, bottom=True, left=True, right=True)
    
    #Plot:
    binno = 50
    TotMetalsCold = np.nansum(Samp1['G_samp']['MetalsColdGas'], axis=1)
    Zs = np.log10(TotMetalsCold / Samp1['G_samp']['ColdGas']) - np.log10(Z_mf_bulk_A09)
    theX = np.log10(Samp1['G_samp']['StellarMass'][(TotMetalsCold > 0.0) & (Samp1['G_samp']['ColdGas'] > 0.0)])
    theY = Zs[(TotMetalsCold > 0.0) & (Samp1['G_samp']['ColdGas'] > 0.0)]
    robs_contour_plot(theX, theY, binno, '2sig', noOutliers=1)  
    robs_plot_average(theX, theY, binno=10, aveType='Median', minInBin=50, linewidth=4., inputCentres=1)

    #title:
    plt.text(.5, 0.92, ['M*-Zg (global, mw, metals arrays)', 'z='+char_z_low, str(Samp1['NumGals'])+' gals', Samp1['Cosmology'], Samp1['Simulation'], Samp1['File_type'], Samp1['Model'], Samp1['Version'], Samp1['Sample_type']], transform=fig.transFigure, horizontalalignment='center', color='blue')  
    
    #print figure to pdf:
    pdf.savefig()  
    
    ##########    
    # Plot ColdGas metallicity as number-weighted 12+log(O/H):
    ########## 
    #Set axis limits:
    xlimits = np.array([9.0,12.0])
    ylimits = np.array([7.0,9.5])
    
    #Set-up plot:
    fig = plt.figure(figsize=(6,5))
    plt.xlabel(r'log$(M_{*} / \textnormal{M}_{\odot})$')
    plt.ylabel(r'12+log(O/H)')
    plt.xlim(xlimits)
    plt.ylim(ylimits)
    plt.tick_params(direction="in", top=True, bottom=True, left=True, right=True)
    
    #Plot M* - Zg:
    binno = 50
    Zg = 12. + np.log10((Samp1['G_samp']['ColdGas_elements'][:,4] / Samp1['G_samp']['ColdGas_elements'][:,0]) * (H_aw/O_aw))
    ww = (Zg > 0.0) & (np.isfinite(Zg)) & (~np.isnan(Zg))      
    theX = np.log10(Samp1['G_samp']['StellarMass'][ww])
    theY = Zg[ww]
    robs_contour_plot(theX, theY, binno, '2sig', noOutliers=1)
    robs_plot_average(theX, theY, binno=10, aveType='Median', minInBin=50, linewidth=4., inputCentres=1)

    #title:
    plt.text(.5, 0.92, ['M*-Zg (global, number weighted)', 'z='+char_z_low, str(Samp1['NumGals'])+' gals', Samp1['Cosmology'], Samp1['Simulation'], Samp1['File_type'], Samp1['Model'], Samp1['Version'], Samp1['Sample_type']], transform=fig.transFigure, horizontalalignment='center', color='blue')  

    #print figure to pdf:
    pdf.savefig()
    
    print("M*-Zg plotted")
    
    
    
    
#################
def plot_mzsr(Samp1, char_z_low, pdf=None) : 
    ##########    
    # Plot Stellar metallicity as mass-weighted M_Z/M_tot, using metals arrays, normalised to Sun:
    ########## 
    #Set axis limits:
    xlimits = np.array([9.0,12.0])
    ylimits = np.array([-1.5,0.5])
    
    #Set-up plot:
    fig = plt.figure(figsize=(6,5))
    plt.xlabel(r'log$(M_{*} / \textnormal{M}_{\odot})$')
    plt.ylabel(r'$Z_{*} / \textnormal{Z}_{\odot}$')
    plt.xlim(xlimits)
    plt.ylim(ylimits)
    plt.tick_params(direction="in", top=True, bottom=True, left=True, right=True)
    
    #Plot:
    binno = 50
    TotMetalsStars = np.nansum(Samp1['G_samp']['MetalsDiskMass'], axis=1) + np.nansum(Samp1['G_samp']['MetalsBulgeMass'], axis=1)
    Zs = np.log10(TotMetalsStars / Samp1['G_samp']['StellarMass']) - np.log10(Z_mf_bulk_A09)
    theX = np.log10(Samp1['G_samp']['StellarMass'][TotMetalsStars > 0.0])
    theY = Zs[TotMetalsStars > 0.0]
    theX1_ma = theX.copy()
    theY1_ma = theY.copy()
    robs_contour_plot(theX, theY, binno, '2sig', noOutliers=1)  
    robs_plot_average(theX, theY, binno=10, aveType='Median', minInBin=50, linewidth=4., inputCentres=1)

    ##########
    #title:
    plt.text(.5, 0.92, ['M*-Zs (global, mw, metals arrays)', 'z='+char_z_low, str(Samp1['NumGals'])+' gals', Samp1['Cosmology'], Samp1['Simulation'], Samp1['File_type'], Samp1['Model'], Samp1['Version'], Samp1['Sample_type']], transform=fig.transFigure, horizontalalignment='center', color='blue')  

    #print figure to pdf:
    pdf.savefig()    
    
    
    ##########
    # Plot Stellar metallicity as mass-weighted [Fe/H], normalised to Sun:
    ##########           
    #Set axis limits:
    xlimits = np.array([9.0,12.0])
    ylimits = np.array([-1.5,0.5])
    
    #Set-up plot:
    fig = plt.figure(figsize=(6,5))
    plt.xlabel(r'log$(M_{*} / \textnormal{M}_{\odot})$')
    plt.ylabel(r'[Fe/H]')
    plt.xlim(xlimits)
    plt.ylim(ylimits)
    plt.tick_params(direction="in", top=True, bottom=True, left=True, right=True)
    
    #Plot:
    binno = 50
    M_stars = np.nansum(Samp1['G_samp']['DiskMass_elements'][:,:], axis=1) +  np.nansum(Samp1['G_samp']['BulgeMass_elements'][:,:], axis=1)
    MH_stars = Samp1['G_samp']['DiskMass_elements'][:,0] + Samp1['G_samp']['BulgeMass_elements'][:,0]
    MZ_stars = Samp1['G_samp']['DiskMass_elements'][:,10] + Samp1['G_samp']['BulgeMass_elements'][:,10]
    Zs = np.log10(MZ_stars / MH_stars) - np.log10(FeH_mf_A09)
    theX = np.log10(Samp1['G_samp']['StellarMass'][(M_stars > 0.0) & (MZ_stars > 0.0)])
    theY = Zs[(M_stars > 0.0) & (MZ_stars > 0.0)]
    robs_contour_plot(theX, theY, binno, '2sig', noOutliers=1)  
    robs_plot_average(theX, theY, binno=10, aveType='Median', minInBin=50, linewidth=4., inputCentres=1)

    #title:
    plt.text(.5, 0.92, ['M*-Zs (global, mw, elements arrays)', 'z='+char_z_low, str(Samp1['NumGals'])+' gals', Samp1['Cosmology'], Samp1['Simulation'], Samp1['File_type'], Samp1['Model'], Samp1['Version'], Samp1['Sample_type']], transform=fig.transFigure, horizontalalignment='center', color='blue')  
    
    #print figure to pdf:
    pdf.savefig()
               
    print("M*-Zs plotted")






#################
def plot_sfhs(Samp1, SFH_bins, snap_z0, char_z_low, pdf=None) :   
    #Get SFH bin info:
    if Samp1['File_type'] == 'snapshots' :
        theSnap = Samp1['G_samp']['SnapNum'][0] #The snapshot number corresponding to this snapshot file
    elif Samp1['File_type'] == 'galtree' :
        theSnap = snap_z0
    SFH_bins_lbt_z0 = SFH_bins.LOOKBACKTIME[SFH_bins.SNAPNUM == theSnap]/1.e9 #Gyr #Lookback times from z=0 to centre of each z=0 SFH bin
    SFH_bins_dt_z0 = SFH_bins.DT[SFH_bins.SNAPNUM == theSnap] #yr #Width of bins from the z=0 SFH
    SFH_bins_num = SFH_bins.BIN[SFH_bins.SNAPNUM == theSnap][-1] #Number of bins active in the z=0 SFH
    
    #Set-up plot:
    fig = plt.figure(figsize=(15,12))
    rowno = 2
    colno = 2
     
    p1 = plt.subplot(rowno, colno, 1)
    plt.xlabel(r'$t_{\rm lookback}$ / Gyr')
    plt.ylabel(r'SFR / M$_{\odot}$/yr')
    plt.tick_params(direction="in", top=True, bottom=True, left=True, right=True)
    plt.yscale('log')
    #################
    numSFHs = Samp1['NumGals'] #the number of galaxy SFHs to plot here
    #Plot individual SFHs:
    allSFHs = (Samp1['G_samp']['sfh_DiskMass'][0:numSFHs,0:SFH_bins_num]+Samp1['G_samp']['sfh_BulgeMass'][0:numSFHs,0:SFH_bins_num])/SFH_bins_dt_z0
    NumToPlot = 100
    pick = np.random.choice(numSFHs, NumToPlot, replace=False)
    SFHs = (Samp1['G_samp']['sfh_DiskMass'][pick,0:SFH_bins_num]+Samp1['G_samp']['sfh_BulgeMass'][pick,0:SFH_bins_num])/SFH_bins_dt_z0
    lbts = np.repeat([SFH_bins_lbt_z0], NumToPlot, axis=0)
    SFHindp, = p1.plot(lbts[0], SFHs[0], '-', color='blue', label='Random subset', alpha=0.1)
    segments = np.array(list(zip(lbts,SFHs))).swapaxes(1,2)
    line_segments = mcoll.LineCollection(segments, colors='blue', alpha=0.1)
    p1.add_collection(line_segments)
    p1.set(xlim=(0.0, 13.5), ylim=(0.01,500.))

    #Plot mean SFR in bins of lookbacktime:
    SFHmeanp, = p1.plot(SFH_bins_lbt_z0, np.nanmean(allSFHs[0:numSFHs], axis=0), 'o', color='blue', markeredgecolor='blue', label='Mean SFR')
    SFHmedianp, = p1.plot(SFH_bins_lbt_z0, np.nanmedian(allSFHs[0:numSFHs], axis=0), 'o', color='white', markeredgecolor='blue',  markeredgewidth=2., label='Median SFR')
    p1.errorbar(SFH_bins_lbt_z0, np.nanmedian(allSFHs[0:numSFHs], axis=0), \
                yerr=[np.nanmedian(allSFHs[0:numSFHs], axis=0)-np.nanpercentile(allSFHs[0:numSFHs], 16, axis=0), np.nanpercentile(allSFHs[0:numSFHs], 84, axis=0)-np.nanmedian(allSFHs[0:numSFHs], axis=0)], \
                color='white', marker='.', ecolor='blue', linestyle='', linewidth=1.49, capsize=3)
    legend1 = plt.legend(handles=[SFHmedianp,SFHmeanp,SFHindp], \
                         title=Samp1['Version'], loc='upper left', fontsize='small', frameon=True, labelspacing=0.3)
    
    #title:
    plt.text(.0, 0.91, ['SFHs', 'z='+char_z_low, str(Samp1['NumGals'])+' gals', Samp1['Cosmology'], Samp1['Simulation'], Samp1['File_type'], Samp1['Model'], Samp1['Version'], Samp1['Sample_type']], transform=fig.transFigure, horizontalalignment='left', color='blue')  
    
    pdf.savefig() 
    
    print("SFHs plotted")



    
#################
def plot_sfrd_prof(Samp1, struct1, char_z_low, MassBins = np.array([[9.0,9.6],[9.6,10.2],[10.2,10.8],[10.8,11.4]]), \
                   pdf=None) :  
    #Set plot parameters:
    plt.figure(figsize=(12,10))    
    xlimits = np.array([0.0,5.0])
    ylimits = np.array([-5.,1.])
    xlabel=r'$R/\textnormal{R}_{\textnormal{e}}$'
    ylabel=r'log($\Sigma_{\textnormal{SFR}} / \textnormal{M}_{\odot}\textnormal{yr}^{-1}\textnormal{kpc}^{-2}$)'
    binwidth = 0.2 #R/Re

    #Calc properties to plot:
    theX = Samp1['RReRings'][:,:]               
    theY = np.log10(Samp1['G_samp']['SfrRings'][:,:]/Samp1['RingArea']) #to keep Msun/yr/kpc^2
    panelProp = np.log10(Samp1['G_samp']['StellarMass'])
    xrange = [theX[np.isfinite(theX)].min(), theX[np.isfinite(theX)].max()]
              
    #Plot profiles:
    robs_plot_profile(theX, theY, panelProp, xlimits, ylimits, xlabel, ylabel, \
                      (xrange[1]-xrange[0])/binwidth, MassBins, aveType='All', \
                      minInBin=1, plot1sig=1, plot2sig=1, inputCentres=1, panelLabel=r'log$(M_{*}/\textnormal{M}_{\odot}) =$ ')               
    #Title:
    ypad = (ylimits[1]-ylimits[0])/15.
    plt.text(xlimits[1], ylimits[1]+ypad, ['SFR density profiles', 'z='+char_z_low, str(Samp1['NumGals'])+' gals', Samp1['Cosmology'], Samp1['Simulation'], Samp1['File_type'], Samp1['Model'], Samp1['Version'], Samp1['Sample_type']], horizontalalignment='right', color='blue')  

    #print figure to pdf:
    pdf.savefig()    
    print("SFR density profiles plotted") 




def plot_cosmic_SFRD(Volume, Hubble_h, FullRedshiftList, FullSnapnumList, RedshiftsToRead, Samp1, struct1, \
                     char_z_low, Samp2=None, struct2=None, pdf=None) :  
    
    #Set-up plot:
    xlimits = np.array([0.0,9.2])
    ylimits = np.array([-2.5,-0.5])
    fig = plt.figure(figsize=(9,5))
    plt.xlabel(r'Redshift')
    plt.ylabel(r'SFRD')
    plt.xlim(xlimits)
    plt.ylim(ylimits)
    plt.tick_params(direction="in", top=True, bottom=True, left=True, right=True)
    
    #Madau & Dickinson 2014 (eqn. 15):
    MD14_redshift = np.arange(xlimits[0], xlimits[1], (xlimits[1]-xlimits[0])/100.)
    MD14_h = 0.7
    SFRD_MD14 = 0.015*(((1.+MD14_redshift)**(2.7)) / (1.+((1.+MD14_redshift)/2.9)**(5.6)))
    logSFRD_MD14_cor = np.log10(SFRD_MD14) + np.log10(0.63) #MD14 assumed a Salpter IMF, so a factor of 0.63 is required to convert to a Chabrier IMF.
    logSFRD_MD14_cor = logSFRD_MD14_cor + np.log10(MD14_h**3/Hubble_h**3) #Adding a extra factor of ~0.038 dex to account for differences in Hubble h in the volume. (20-12-20)

    #Driver+18 (from Bruno [21-08-18]):
    D18_h = 0.7
    D18_redshift = np.array([0.05,0.10,0.17,0.24,0.32,0.40,0.50,0.62,0.75,0.91,1.10,1.32,1.60,1.98,2.40,2.93,3.50,4.00,4.63])
    D18_logSFRD = np.array([-1.95,-1.82,-1.90,-1.77,-1.75,-1.79,-1.73,-1.56,-1.42,-1.29,-1.31,-1.27,-1.17,-1.30,-1.29,-1.28,-1.33,-1.42,-1.45])
    D18_cor_EddBias   = np.array([0.03, 0.03, 0.02, 0.01, 0.01, 0.01, 0.04, 0.05, 0.06, 0.05, 0.04, 0.03, 0.02, 0.04, 0.04, 0.04, 0.03, 0.04, 0.03])
    D18_err_poisson   = np.array([0.00, 0.01, 0.00, 0.00, 0.00, 0.01, 0.01, 0.00, 0.01, 0.00, 0.00, 0.00, 0.00, 0.01, 0.01, 0.01, 0.01, 0.04, 0.04])
    D18_err_cv        = np.array([0.07, 0.05, 0.04, 0.05, 0.06, 0.06, 0.09, 0.07, 0.06, 0.07, 0.05, 0.06, 0.06, 0.07, 0.04, 0.04, 0.03, 0.05, 0.04])
    D18_err_agn       = np.array([0.00, 0.01, 0.00, 0.00, 0.01, 0.01, 0.03, 0.02, 0.04, 0.01, 0.01, 0.02, 0.03, 0.06, 0.09, 0.11, 0.08, 0.02, 0.04])
    D18_logSFRD_err = np.sqrt((D18_err_poisson)**2 + (D18_err_cv)**2 + (D18_err_agn)**2) #Eddington bias correction has already been *subtracted* from the D18_logSFRD.
    D18_logSFRD_cor = D18_logSFRD + np.log10(D18_h**3/Hubble_h**3)
    D18_logSFRD_cor_err = D18_logSFRD_err + np.log10(D18_h**3/Hubble_h**3)
    
    #Plot obs data:
    MD14p, = plt.plot(MD14_redshift, logSFRD_MD14_cor, 'orange', label='Madau \& Dickinson (2014)') 
    D18p, = plt.plot(D18_redshift, D18_logSFRD_cor, color='orange', marker='o', linestyle='', label='Driver et al. (2018)') #Adding a extra factor of ~0.038 dex to account for differences in Hubble h in the volume. (20-12-20)
    plt.errorbar(D18_redshift, D18_logSFRD_cor, yerr=D18_logSFRD_cor_err, color='orange', marker='o', linestyle='', linewidth=1.49, capsize=3)

    #Plot model data:
    logSFRD = np.empty(len(FullRedshiftList))
    for iz in range(0,len(FullRedshiftList)) : #loop through redshifts
            if RedshiftsToRead[iz] : 
                wclean = np.where((Samp1['G_samp']['SnapNum'] == FullSnapnumList[iz]) & (Samp1['G_samp']['Sfr'] > 0.0))
                SFR_clean = np.nansum(Samp1['G_samp']['Sfr'][wclean])
                logSFRD[iz] = np.log10(SFR_clean / (Volume/Hubble_h**3))
    FRL = np.array(FullRedshiftList)
    RTR = np.array(RedshiftsToRead)
    plt.scatter(FRL[RTR], logSFRD, marker='o', color='green')
 
    #title:
    plt.text(.0, 0.91, ['cosmic SFRD', Samp1['Cosmology'], Samp1['Simulation'], Samp1['File_type'], Samp1['Model'], Samp1['Version'], Samp1['Sample_type']], transform=fig.transFigure, horizontalalignment='left', color='blue')  
    
    pdf.savefig() 
    
    print("SFHs plotted")
    

#################
#Plot tau_acc, tau_shock, and tau_sput,CGM, and tau_sput,Ejec:
def plot_timescales(Samp1, struct1, redshift, xprop='OH', \
                    yprops='All', contourLines='solid', outlierFrac=defaultOutlierFrac, SolarNorm=None, label=1) : 
    pdf = PdfPages(PlotDir+"dust_timescales"+"_"+xprop+"_"+Samp1['Model']+"_"+Samp1['Sample_type']+"_z"+str(redshift)+".pdf")
    SolarAbunds = robs_solarnorm_finder(SolarNorm)
    El1 = robs_element_list_finder(struct1)
    if yprops == 'All' :
        yprops = ['tau_acc','tau_shock','tau_sput_hotGas','tau_sput_ejecta']
    
    #Set-up multi-panel plot:
    fig = plt.figure(figsize=(4.5*len(yprops),4))
    panelNum = 0
    rowno = 1
    colno = len(yprops)
    binno = 35
    VolRings = np.repeat(Samp1["Volumes"][:,np.newaxis], Samp1["RNUM"], axis=1) #Extend Volumes array to an extra dimension
    Samp1Linestyle = contourLines
    Samp1Alpha = 1.0
    Samp1Fill = 1 #0
    if (outlierFrac == 0.0) | (outlierFrac == None) :
        Samp1NoOutliers = 1
    else :
        Samp1NoOutliers = 0
    Samp1Linecolour = 'black'
    Samp1zorder = 2
        
    #####
    #tau_acc for rings: 
    if 'tau_acc' in yprops :
        panelNum += 1
        ax = plt.subplot(rowno, colno, panelNum)
        if xprop == 'OH' :
            xlimits = np.array([6.5,9.75])
            plt.xlabel(r'12+log(O/H)')
            ylimits = np.array([4.5,12.0])
            plt.ylabel(r'log$(\tau_{\rm acc}/\textnormal{yr})$')
        if xprop == 'logZ' :
            xlimits = np.array([-2.5,1.0])
            plt.xlabel(r'log($Z/Z_{\odot}$)')
            ylimits = np.array([-2.0,7.0])
            plt.ylabel(r'log$(\tau_{\rm acc}/\textnormal{Myr})$')
        elif xprop == 'H2D' :
            xlimits = np.array([-2.0,4.0])
            plt.xlabel(r'log($\Sigma_{\rm H2} / \textnormal{M}_{\odot}\textnormal{ pc}^{-2}$)')
            ylimits = np.array([-2.0,7.0])
            plt.ylabel(r'log$(\tau_{\rm acc}/\textnormal{Myr})$')
        plt.xlim(xlimits)
        plt.ylim(ylimits)
       
        t_acc_rings = Samp1['G_samp']['t_acc'][:,:]
        if xprop == 'OH' :
            Zg = 12. + np.log10((Samp1['G_samp']['ColdGasRings_elements'][:,:,El1["O_NUM"]] / Samp1['G_samp']['ColdGasRings_elements'][:,:,El1["H_NUM"]]) * (H_aw/O_aw))        
            ww = np.where((t_acc_rings > 0.0) & (np.isfinite(t_acc_rings)) & (~np.isnan(t_acc_rings)) & (Zg > 0.0) & (np.isfinite(Zg)) & (~np.isnan(Zg)))  
            theX = Zg[ww]
            theY = np.log10(t_acc_rings[ww]) #log(t_acc/yr)
        elif xprop == 'logZ' :
            Zg = np.log10((np.nansum(Samp1['G_samp']['MetalsColdGasRings'], axis=2) / Samp1['G_samp']['ColdGasRings']) / SolarAbunds['Z_mf_sun'])
            ww = np.where((t_acc_rings > 0.0) & (np.isfinite(t_acc_rings)) & (~np.isnan(t_acc_rings)) & (np.isfinite(Zg)) & (~np.isnan(Zg)))  
            theX = Zg[ww]
            theY = np.log10(t_acc_rings[ww]/1.e6) #log(t_acc/Myr)
        elif xprop == 'H2D' :
            H2D = np.log10((Samp1['G_samp']['ColdGasRings_elements'][:,:,El1["H_NUM"]] * (Samp1['G_samp']['H2fractionRings']))/Samp1['RingArea']/(1000.**2)) #to get Msun/pc^2 
            ww = np.where((t_acc_rings > 0.0) & (np.isfinite(t_acc_rings)) & (~np.isnan(t_acc_rings)) & (np.isfinite(H2D)) & (~np.isnan(H2D)))  
            theX = H2D[ww]
            theY = np.log10(t_acc_rings[ww]/1.e6) #log(t_acc/Myr)
        theWeights = 1./VolRings[ww]
        robs_contour_plot(theX, theY, binno, '4sig', noOutliers=Samp1NoOutliers, fill=Samp1Fill, alpha=Samp1Alpha, theRange=[xlimits,ylimits], weights=theWeights, \
                          linestyle=Samp1Linestyle, linecolour=Samp1Linecolour, outlierFrac=outlierFrac, zorder=Samp1zorder)
        
        robs_plot_text(ax, r'z = '+str(redshift), vpos='bottom', hpos='left')
        if label :
            robs_plot_text(ax, [r'MM+\textsc{binary\_c}', r'MM+singleStars'], vpos='top', hpos='right', colour=[model1avecol,model2col], padder=0.06)
           
    #####
    #tau_shock for rings: 
    if 'tau_shock' in yprops :
        panelNum += 1
        xlimits = np.array([6.5,9.75]) 
        ylimits = np.array([7.0,19.0])
        plt.subplot(rowno, colno, panelNum)
        plt.xlabel(r'12+log(O/H)')
        plt.ylabel(r'log$(\tau_{\rm shock}/\textnormal{yr})$')
        plt.xlim(xlimits)
        plt.ylim(ylimits)
     
        Zg = 12. + np.log10((Samp1['G_samp']['ColdGasRings_elements'][:,:,El1["O_NUM"]] / Samp1['G_samp']['ColdGasRings_elements'][:,:,El1["H_NUM"]]) * (H_aw/O_aw))        
        t_des_rings = Samp1['G_samp']['t_des'][:,:]
        ww = np.where((t_des_rings > 0.0) & (np.isfinite(t_des_rings)) & (~np.isnan(t_des_rings)) & (Zg > 0.0) & (np.isfinite(Zg)) & (~np.isnan(Zg)))  
        theX = Zg[ww]
        theY = np.log10(t_des_rings[ww])
        theWeights = 1./VolRings[ww]
        robs_contour_plot(theX, theY, binno, '4sig', noOutliers=Samp1NoOutliers, fill=Samp1Fill, alpha=Samp1Alpha, theRange=[xlimits,ylimits], weights=theWeights, \
                          linestyle=Samp1Linestyle, linecolour=Samp1Linecolour, outlierFrac=outlierFrac, zorder=Samp1zorder)

    #####
    #tau_sput for HotGas:  
    if 'tau_sput_hotGas' in yprops :
        panelNum += 1
        xlimits = np.array([8.5,12.0])
        plt.subplot(rowno, colno, panelNum)
        plt.xlabel(r'log$(M_{*} / \textnormal{M}_{\odot})$')
        plt.ylabel(r'log$(\tau_{\rm sput,CGM}/\textnormal{yr})$')
        plt.xlim(xlimits)
        plt.ylim(ylimits)
      
        tau_sput = Samp1['G_samp']['t_sput_HotGas']*1.e9 #convert from Gyr to yr      
        ww = np.where((tau_sput > 0.0) & (np.isfinite(tau_sput)) & (~np.isnan(tau_sput)))  
        theX = np.log10(Samp1['G_samp']['StellarMass'][ww])
        theY = np.log10(tau_sput[ww])
        theWeights = 1./Samp1["Volumes"][ww]
        robs_contour_plot(theX, theY, binno, '3sig', noOutliers=Samp1NoOutliers, fill=Samp1Fill, alpha=Samp1Alpha, theRange=[xlimits,ylimits], weights=theWeights, \
                          linestyle=Samp1Linestyle, linecolour=Samp1Linecolour, outlierFrac=outlierFrac, zorder=Samp1zorder)

    #####
    #tau_sput for EjectedMass:  
    if 'tau_sput_ejecta' in yprops :
        panelNum += 1
        xlimits = np.array([8.5,12.0])
        plt.subplot(rowno, colno, panelNum)
        plt.xlabel(r'log$(M_{*} / \textnormal{M}_{\odot})$')
        plt.ylabel(r'log$(\tau_{\rm sput,Ejecta}/\textnormal{yr})$')
        plt.xlim(xlimits)
        plt.ylim(ylimits)
          
        tau_sput = Samp1['G_samp']['t_sput_EjectedMass']*1.e9 #convert from Gyr to yr     
        ww = np.where((tau_sput > 0.0) & (np.isfinite(tau_sput)) & (~np.isnan(tau_sput)))  
        theX = np.log10(Samp1['G_samp']['StellarMass'][ww])
        theY = np.log10(tau_sput[ww])
        theWeights = 1./Samp1["Volumes"][ww]
        robs_contour_plot(theX, theY, binno, '3sig', noOutliers=Samp1NoOutliers, fill=Samp1Fill, alpha=Samp1Alpha, theRange=[xlimits,ylimits], weights=theWeights, \
                          linestyle=Samp1Linestyle, linecolour=Samp1Linecolour, outlierFrac=outlierFrac, zorder=Samp1zorder)

    pdf.savefig(bbox_inches='tight') 
    pdf.close()
    
    
#################
def plot_smf(Hubble_h, Samp1, redshift, Add_Edd_bias=None) :  
    pdf = PdfPages(PlotDir+"smf"+"_"+Samp1['Model']+"_"+Samp1['Sample_type']+"_z"+str(redshift)+".pdf")
    
    xlimits = np.array([8.0,12.0])
    ylimits = np.array([-5.9,-0.5])
    
    #Set-up plot:
    fig, ax = plt.subplots(figsize=(5,5))
    plt.xlabel(r'log$(M_{*} / \textnormal{M}_{\odot})$')
    plt.ylabel(r'log$(\phi / \textnormal{Mpc}^{-3} \textnormal{ dex}^{-1})$')
    plt.xlim(xlimits)
    plt.ylim(ylimits)
    binwidth = 0.1
    bin_arr = np.arange(xlimits[0],xlimits[1]+binwidth,binwidth)
    
    if ((Samp1["MRI_gals"] > 0) & (Samp1["MRII_gals"] > 0)) :
        prop_MRI = np.log10(Samp1['G_samp']['StellarMass'][0:Samp1['MRI_gals']])
        prop_MRII = np.log10(Samp1['G_samp']['StellarMass'][Samp1['MRI_gals']:])
        if Add_Edd_bias :
            prop_MRI += np.random.randn(len(Samp1['G_samp']['StellarMass'][0:Samp1['MRI_gals']]))*0.08*(1+float(Samp1['z_lower'])) #Add Eddington bias to more fa
            prop_MRII += np.random.randn(len(Samp1['G_samp']['StellarMass'][Samp1['MRI_gals']:]))*0.08*(1+float(Samp1['z_lower']))
        bin_arr_MRI = np.arange(np.log10(Samp1["MRI_cutoff"]),xlimits[1]+binwidth,binwidth)
        bin_arr_MRII = np.arange(xlimits[0],np.log10(Samp1["MRI_cutoff"])+(2.*binwidth),binwidth)
        hist_MRI = np.histogram(prop_MRI, bins=bin_arr_MRI, range=(np.log10(Samp1["MRI_cutoff"]),xlimits[1])) 
        hist_MRII = np.histogram(prop_MRII, bins=bin_arr_MRII, range=(xlimits[0],np.log10(Samp1["MRI_cutoff"]))) 
        plt.plot(hist_MRI[1][0:len(hist_MRI[1][:])-1]+binwidth/2., np.log10((hist_MRI[0][:]/(Samp1["Volumes"][0]*binwidth))*Hubble_h**3), \
                 linewidth=defaultLinewidth, color=model1avecol)#, label=Samp1['Label']) 
        plt.plot(hist_MRII[1][0:len(hist_MRII[1][:])-1]+binwidth/2., np.log10((hist_MRII[0][:]/(Samp1["Volumes"][Samp1['MRI_gals']]*binwidth))*Hubble_h**3), \
                 linewidth=defaultLinewidth, linestyle='dotted', color=model1avecol) 
    else :
        prop = np.log10(Samp1['G_samp']['StellarMass'])
        if Add_Edd_bias :
            prop += np.random.randn(len(Samp1['G_samp']['StellarMass']))*0.08*(1+float(Samp1['z_lower']))
        hist = np.histogram(prop, bins=bin_arr, range=(xlimits[0],xlimits[1])) 
        plt.plot(hist[1][0:len(hist[1][:])-1]+binwidth/2., np.log10((hist[0][:]/(Samp1['Volumes'][0]*binwidth))*Hubble_h**3), \
                 linewidth=defaultLinewidth, color=model1avecol, label=Samp1['Label'])   
   
    #Labels:
    robs_plot_text(ax, r'a)', hpos='right', vpos='top')
    
    pdf.savefig(bbox_inches='tight')     
    pdf.close()

    
#################
#Plot HI mass function (HIMF):
def plot_himf(Hubble_h, Samp1, struct1, redshift, prop="ColdGas") :  
    pdf = PdfPages(PlotDir+"himf"+"_"+Samp1['Model']+"_"+Samp1['Sample_type']+"_z"+str(redshift)+".pdf")
    if prop == "H" :
        El1 = robs_element_list_finder(struct1)
    xlimits = np.array([7.0,11.0])
    ylimits = np.array([-5.9,0.0])
    
    #Set-up plot:
    fig, ax = plt.subplots(figsize=(5,5))
    plt.xlabel(r'log$(M_{\rm HI} / \textnormal{M}_{\odot})$')
    plt.ylabel(r'log$(\phi / \textnormal{Mpc}^{-3} \textnormal{ dex}^{-1})$')
    plt.xlim(xlimits)
    plt.ylim(ylimits)
    
    #Plot HIMF:
    binwidth = 0.1
    bin_arr = np.arange(xlimits[0],xlimits[1]+binwidth,binwidth)
    if ((Samp1["MRI_gals"] > 0) & (Samp1["MRII_gals"] > 0)) :
        MRI_HI_cutoff = 10**(9.5)
        if prop == "ColdGas" :
            prop_MRI = np.log10(Samp1['G_samp']['ColdGas'][0:Samp1['MRI_gals']]*(1.-Samp1['G_samp']['H2fraction'][0:Samp1['MRI_gals']]))
            prop_MRII = np.log10(Samp1['G_samp']['ColdGas'][Samp1['MRI_gals']:]*(1.-Samp1['G_samp']['H2fraction'][Samp1['MRI_gals']:]))
        elif prop == "H" :
            prop_MRI = np.log10(Samp1['G_samp']['ColdGas_elements'][0:Samp1['MRI_gals'],El1["H_NUM"]]*(1.-Samp1['G_samp']['H2fraction'][0:Samp1['MRI_gals']]))
            prop_MRII = np.log10(Samp1['G_samp']['ColdGas_elements'][Samp1['MRI_gals']:,El1["H_NUM"]]*(1.-Samp1['G_samp']['H2fraction'][Samp1['MRI_gals']:]))
        bin_arr_MRI = np.arange(np.log10(MRI_HI_cutoff),xlimits[1]+binwidth,binwidth)
        bin_arr_MRII = np.arange(xlimits[0],np.log10(MRI_HI_cutoff)+(2.*binwidth),binwidth)
        hist_MRI = np.histogram(prop_MRI, bins=bin_arr_MRI, range=(np.log10(MRI_HI_cutoff),xlimits[1])) 
        hist_MRII = np.histogram(prop_MRII, bins=bin_arr_MRII, range=(xlimits[0],np.log10(MRI_HI_cutoff))) 
        plt.plot(hist_MRI[1][0:len(hist_MRI[1][:])-1]+binwidth/2., np.log10((hist_MRI[0][:]/(Samp1["Volumes"][0]*binwidth))*Hubble_h**3), \
                  linewidth=defaultLinewidth, color=model1avecol)#, label=Samp1['Label']) 
        plt.plot(hist_MRII[1][0:len(hist_MRII[1][:])-1]+binwidth/2., np.log10((hist_MRII[0][:]/(Samp1["Volumes"][Samp1['MRI_gals']]*binwidth))*Hubble_h**3), \
                  linewidth=defaultLinewidth, linestyle='dotted', color=model1avecol) 
    else :
        if prop == "ColdGas" :
            propAll = np.log10(Samp1['G_samp']['ColdGas']*(1.-Samp1['G_samp']['H2fraction']))
        elif prop == "H" :
            propAll = np.log10(Samp1['G_samp']['ColdGas_elements'][:,El1["H_NUM"]]*(1.-Samp1['G_samp']['H2fraction']))
        hist = np.histogram(propAll, bins=bin_arr, range=(xlimits[0],xlimits[1])) 
        plt.plot(hist[1][0:len(hist[1][:])-1]+binwidth/2., np.log10((hist[0][:]/(Samp1['Volumes'][0]*binwidth))*Hubble_h**3), \
             linewidth=defaultLinewidth, color=model1avecol)  
 
    #Plot panel label:
    robs_plot_text(ax, r'b)', hpos='right', vpos='top')

    pdf.savefig(bbox_inches='tight')     
    pdf.close()
    

#################
#Plot M* - sSFR relation:
def plot_mssfr(Hubble_h, Samp1, redshift, contourLines=None) :  
    pdf = PdfPages(PlotDir+"mstar-ssfr"+"_"+Samp1['Model']+"_"+Samp1['Sample_type']+"_z"+str(redshift)+".pdf")
    
    if Samp1["MRII_gals"] > 0 :
        xlimits = np.array([8.5,12.0])
    else :
        xlimits = np.array([9.0,12.0])
    ylimits = np.array([-14.5,-8.0])
    
    if contourLines is None :
        ls1 = 'solid'
    else :
        ls1 = contourLines
    
    #Set-up plot:
    fig, ax = plt.subplots(figsize=(5,5))
    plt.xlabel(r'log$(M_{*} / \textnormal{M}_{\odot})$')
    plt.ylabel(r'log$(\textnormal{sSFR} / \textnormal{yr}^{-1})$')
    plt.xlim(xlimits)
    plt.ylim(ylimits)
    
    binno = 25
    theX = np.log10(Samp1['G_samp']['StellarMass'][Samp1['G_samp']['Sfr']>0.0])
    theY = np.log10(Samp1['G_samp']['Sfr'][Samp1['G_samp']['Sfr']>0.0]/Samp1['G_samp']['StellarMass'][Samp1['G_samp']['Sfr']>0.0])
    robs_contour_plot(theX, theY, binno, '3sig', noOutliers=0, fill=1, alpha=1.0, theRange=[xlimits,ylimits], \
                      colour=model1col, linestyle=ls1, weights=1./Samp1['Volumes'][Samp1['G_samp']['Sfr']>0.0], outlierFrac=defaultOutlierFrac)

    #Plot panel label:
    robs_plot_text(ax, r'c)', hpos='right', vpos='top')
    
    pdf.savefig(bbox_inches='tight')     
    pdf.close()  
    

#################
def plot_snrates(Samp1, redshift) :     
    pdf = PdfPages(PlotDir+"sSFR-SNRates"+"_"+Samp1['Model']+"_"+Samp1['Sample_type']+"_z"+str(redshift)+".pdf")
    xlimits = np.array([-13.0,-9.0])
    xticks = [-13.0,-12.0,-11.0,-10.0,-9.0]
    ylimits = np.array([[-4.5,3.0],[-5.5,1.5]])
    theMinInBin = 1
    
    #Set-up plot:
    rows=2
    columns=1
    numPanels = rows*columns
    xlabs = [r'log(sSFR / yr$^{-1}$)',r'log(sSFR / yr$^{-1}$)']
    ylabs = [r'log($R_{\rm SNII}$ / SNuM)',r'log($R_{\rm SNIa}$ / SNuM)']
    xlims = [xlimits,xlimits]
    ylims = ylimits
    yticks = [[-4.0,-2.0,0.0]]*rows
    fig, ax = plt.subplots(figsize=(5,5))
    for ii in range(numPanels) :
        panel = robs_plot_panels(ii, rows=rows, columns=columns, xlimits=xlims, ylimits=ylims, \
                                 xlab=xlabs, ylab=ylabs, xticks=xticks, yticks=yticks, SecondXAxis=None)
        if ii == 0 :
            #Plot SN-II rate:
            binno = defaultBinno #50
            SNIIRate = Samp1['G_samp']['DiskSNIIRate'] + Samp1['G_samp']['BulgeSNIIRate'] + Samp1['G_samp']['ICMSNIIRate'] #[1/yr]
            SNIIRate_SNuM = SNIIRate * 1.e2 / (Samp1['G_samp']['StellarMass']/1.e10) #[1/100yr * 1/10^10Msun]
            theX = np.log10(Samp1['G_samp']['Sfr'][(SNIIRate_SNuM > 0.0) & (Samp1['G_samp']['Sfr'] > 0.0)] / Samp1['G_samp']['StellarMass'][(SNIIRate_SNuM > 0.0) & (Samp1['G_samp']['Sfr'] > 0.0)])
            theY = np.log10(SNIIRate_SNuM[(SNIIRate_SNuM > 0.0) & (Samp1['G_samp']['Sfr'] > 0.0)])
            theWeights = 1./Samp1["Volumes"][(SNIIRate_SNuM > 0.0) & (Samp1['G_samp']['Sfr'] > 0.0)]
            robs_plot_average(theX, theY, binno=10, aveType='Median', minInBin=theMinInBin, colour=model1col, \
                              linewidth=defaultLinewidth, plot1sig=1, plot2sig=1, inputCentres=1, invertLinestyle=1, \
                              weights=theWeights, plotlims=[xlimits,ylimits[ii]])     
        elif ii == 1 :
            #Plot SN-Ia rate:
            SNIaRate = Samp1['G_samp']['DiskSNIaRate'] + Samp1['G_samp']['BulgeSNIaRate'] + Samp1['G_samp']['ICMSNIaRate'] #[1/yr]
            SNIaRate_SNuM = SNIaRate * 1.e2 / (Samp1['G_samp']['StellarMass']/1.e10) #[1/100yr * 1/10^10Msun]
            theX = np.log10(Samp1['G_samp']['Sfr'][(SNIaRate_SNuM > 0.0) & (Samp1['G_samp']['Sfr'] > 0.0)] / Samp1['G_samp']['StellarMass'][(SNIaRate_SNuM > 0.0) & (Samp1['G_samp']['Sfr'] > 0.0)])
            theY = np.log10(SNIaRate_SNuM[(SNIaRate_SNuM > 0.0) & (Samp1['G_samp']['Sfr'] > 0.0)])
            theWeights = 1./Samp1["Volumes"][(SNIaRate_SNuM > 0.0) & (Samp1['G_samp']['Sfr'] > 0.0)]
            robs_plot_average(theX, theY, binno=10, aveType='Median', minInBin=theMinInBin, colour=model1col, \
                              linewidth=defaultLinewidth, plot1sig=1, plot2sig=1, inputCentres=1, invertLinestyle=1, \
                                  weights=theWeights, plotlims=[xlimits,ylimits[ii]]) #linestyle='--',         
    
        #Plot panel label:
        if ii == 0 :   
            robs_plot_text(ax, r'd)', hpos='right', vpos='top')
    
    pdf.savefig(bbox_inches='tight')     
    pdf.close() 
    
    
#################
#Plot M* - Zg relation: 12+log(O/H), rings-based, SFR weighted:
def plot_mzgr(Hubble_h, Samp1, struct1, redshift, dustCorrec=None, \
              contourLines=None, plotAve=None, plotRedshift=None) :     
    pdf = PdfPages(PlotDir+"mzgr"+"_"+Samp1['Model']+"_"+Samp1['Sample_type']+"_z"+str(redshift)+".pdf")   
    El1 = robs_element_list_finder(struct1)
    
    if Samp1["MRII_gals"] > 0 :
        xlimits = np.array([8.5,12.0])
    else :
        xlimits = np.array([9.0,12.0])
    if redshift > 6.0 :
        xlimits = np.array([7.0,10.6])
        ylimits = np.array([6.0,8.75])
    else :
        ylimits = np.array([7.0,9.75])
    
    if contourLines is None :
        ls1 = 'solid'
    else :
        ls1 = contourLines
    
    #Set-up plot:
    fig, ax = plt.subplots(figsize=(5,5))
    plt.xlabel(r'log$(M_{*} / \textnormal{M}_{\odot})$')
    plt.ylabel(r'12+log(O/H)')
    plt.xlim(xlimits)
    plt.ylim(ylimits)
    binno = 25
    SFRmax = np.max(Samp1['G_samp']['SfrRings'], axis=1) #Max SFR from all rings, for every galaxy
    SFRfrac = Samp1['G_samp']['SfrRings']/SFRmax[:,None] #Fraction of max SFR in each ring, for every galaxy
    Mstar = Samp1['G_samp']['StellarMass']
    Vol = Samp1['Volumes']
    SFRtot = np.nansum(SFRfrac, axis=1) #Sum of the SFR fractions from all rings, for evey galaxy. (Minimum must be 1.0)
    if dustCorrec :    
        H_rings = Samp1['G_samp']['ColdGasCloudsRings_elements'][:,:,El1["H_NUM"]] + Samp1['G_samp']['ColdGasDiffRings_elements'][:,:,El1["H_NUM"]] \
                - Samp1['G_samp']['DustColdGasCloudsRings_elements'][:,:,El1["H_NUM"]] - Samp1['G_samp']['DustColdGasDiffRings_elements'][:,:,El1["H_NUM"]]
        O_rings = Samp1['G_samp']['ColdGasCloudsRings_elements'][:,:,El1["O_NUM"]] + Samp1['G_samp']['ColdGasDiffRings_elements'][:,:,El1["O_NUM"]] \
                - Samp1['G_samp']['DustColdGasCloudsRings_elements'][:,:,El1["O_NUM"]] - Samp1['G_samp']['DustColdGasDiffRings_elements'][:,:,El1["O_NUM"]]
    else :
        H_rings = Samp1['G_samp']['ColdGasRings_elements'][:,:,El1["H_NUM"]]
        O_rings = Samp1['G_samp']['ColdGasRings_elements'][:,:,El1["O_NUM"]]
    Zg_rings = (O_rings / H_rings) * (H_aw/O_aw) * SFRfrac
    Zg_tot = np.nansum(Zg_rings, axis=1)
    Zg = 12. + np.log10(Zg_tot[Zg_tot > 0.0] / SFRtot[Zg_tot > 0.0])
    Mstar = Mstar[Zg_tot > 0.0] 
    Vol = Vol[Zg_tot > 0.0] 
    theX = np.log10(Mstar[Zg > 0.0])
    theY = Zg[Zg > 0.0]
    robs_contour_plot(theX, theY, binno, '3sig', noOutliers=0, fill=1, alpha=1.0, theRange=[xlimits,ylimits], \
                      colour=model1col, linestyle=ls1, weights=1./Vol[Zg > 0.0], outlierFrac=defaultOutlierFrac)
    if plotAve :
        Samp1_ave = robs_plot_average(theX, theY, binno=10, aveType='Median', minInBin=50, \
                          linewidth=4., inputCentres=1, returnValues=1) 
            
    #Plot panel label:
    if plotRedshift :
        robs_plot_text(ax, "z="+char_z_low, hpos='right', vpos='top')
    else :
        robs_plot_text(ax, r'e)', hpos='right', vpos='top')
    
    pdf.savefig(bbox_inches='tight')     
    pdf.close() 


#################
#Plot M* - Zs relation: Z*/Zsun, rings-based, mass weighted, within 3 arcsec (i.e. first 7 rings, for z=0.037):
def plot_mzsr(Hubble_h, Omega_M, Omega_Lambda, Samp1, redshift, \
              SolarNorm="A09_bulk", ApertureCorrec=None, contourLines=None) :  
    pdf = PdfPages(PlotDir+"mzsr"+"_"+Samp1['Model']+"_"+Samp1['Sample_type']+"_z"+str(redshift)+".pdf")
    SolarAbunds = robs_solarnorm_finder(SolarNorm)
    #N.B. SolarNorm default is "A09_bulk" to match the observations that are plotted for comparison in Yates+23.

    if Samp1["MRII_gals"] > 0 :
        xlimits = np.array([minMassPlotMRII,12.0])
    else :
        xlimits = np.array([9.0,12.0])
    ylimits = np.array([-1.5,0.5])
    
    if contourLines is None :
        ls1 = 'solid'
    else :
        ls1 = contourLines
    
    #Set-up plot:
    fig, ax = plt.subplots(figsize=(5,5))
    plt.xlabel(r'log$(M_{*} / \textnormal{M}_{\odot})$')
    plt.ylabel(r'$Z_{*\textnormal{,}<3\textnormal{arcsec}} / \textnormal{Z}_{\odot}$')
    plt.xlim(xlimits)
    plt.ylim(ylimits)
    plt.tick_params(direction="in", top=True, bottom=True, left=True, right=True)
    
    binno = 25
    if ApertureCorrec is None :
        max_ring = Samp1['RNUM']
    else :
        if redshift == 0.0 :
            scale_redshift = 0.037 #Average redshift of MaNGA
        else :
            scale_redshift = redshift
        aperture_diameter_arcsec = 3. #Diameter of an SDSS fibre
        aperture_diameter_kpc = robs_arcsec_to_kpc(aperture_diameter_arcsec, scale_redshift, \
                                                   OmegaM=Omega_M, OmegaL=Omega_Lambda, Hubble_h=Hubble_h)
        for ii in range(Samp1['RNUM']) :
            if Samp1['RingRadii'][ii]*2. < aperture_diameter_kpc :
                max_ring = ii+2 #Has +2 to get the ring that is just beyond the required aperture (rather than just within)
    TotMetalsStarsRings = np.nansum(Samp1['G_samp']['MetalsDiskMassRings'][:,0:max_ring,:], axis=2) + np.nansum(Samp1['G_samp']['MetalsBulgeMassRings'][:,0:max_ring,:], axis=2)
    TotMassRings = Samp1['G_samp']['DiskMassRings'][:,0:max_ring] + Samp1['G_samp']['BulgeMassRings'][:,0:max_ring]
    TotMetalsStars = np.nansum(TotMetalsStarsRings, axis=1)
    ZsRings = np.log10(TotMetalsStarsRings / TotMassRings) - np.log10(SolarAbunds['Z_mf_sun']) #np.log10(Z_mf_bulk_A09)
    Zs = np.nansum(ZsRings*TotMassRings, axis=1) / np.nansum(TotMassRings, axis=1) #Mass-weighted mean ZsRings
    theX = np.log10(Samp1['G_samp']['StellarMass'][TotMetalsStars > 0.0])
    theY = Zs[TotMetalsStars > 0.0]
    robs_contour_plot(theX, theY, binno, '3sig', noOutliers=0, fill=1, alpha=1.0, theRange=[xlimits,ylimits], \
                      colour=model1col, linestyle=ls1, weights=1./Samp1['Volumes'][TotMetalsStars > 0.0], outlierFrac=defaultOutlierFrac)
    #robs_plot_average(theX, theY, binno=10, aveType='Median', minInBin=50, linewidth=4., inputCentres=1)
    
    if Samp2 :
        TotMetalsStarsRings = np.nansum(Samp2['G_samp']['MetalsDiskMassRings'][:,0:max_ring,:], axis=2) + np.nansum(Samp2['G_samp']['MetalsBulgeMassRings'][:,0:max_ring,:], axis=2)
        TotMassRings = Samp2['G_samp']['DiskMassRings'][:,0:max_ring] + Samp2['G_samp']['BulgeMassRings'][:,0:max_ring]
        TotMetalsStars = np.nansum(TotMetalsStarsRings, axis=1)
        ZsRings = np.log10(TotMetalsStarsRings / TotMassRings) - np.log10(SolarAbunds['Z_mf_sun'])
        Zs = np.nansum(ZsRings*TotMassRings, axis=1) / np.nansum(TotMassRings, axis=1) #Mass-weighted mean ZsRings
        theX = np.log10(Samp2['G_samp']['StellarMass'][TotMetalsStars > 0.0])
        theY = Zs[TotMetalsStars > 0.0]
        robs_contour_plot(theX, theY, binno, '3sig', noOutliers=1, fill=0, alpha=1.0, theRange=[xlimits,ylimits], \
                          linestyle=ls2, colour=model2col, weights=1./Samp2['Volumes'][TotMetalsStars > 0.0])
        #robs_plot_average(theX, theY, binno=10, aveType='Median', minInBin=50, colour=samp2col, linestyle='--', linewidth=4., inputCentres=1, invertLinestyle=1)

    if obs :
        Z17 = Z17_MZsR_obs_data()
        Y21a = Y21a_MZgR_obs_data()
        Z17p = plt.errorbar(Z17['logMstar'], Z17['logZs_mw'], \
                yerr=[Z17['logZs_mw_lower_err'],Z17['logZs_mw_upper_err']], \
                color=obscol, marker='s', linestyle='', linewidth=1.6, capsize=3, markeredgecolor='black', \
                #color=obscol, marker='s', markersize=obsMarkerSize, linestyle='', linewidth=1., capsize=3, markeredgecolor='black', \
                markersize=obsMarkerSize, label=r'Zahid+17')
            
        Y21a = Y21a_MZsR_obs_data()
        Y21ap = plt.errorbar(Y21a['logMstar'], Y21a['Zs_mean'], \
                yerr=[Y21a['Zs_mean']-Y21a['Zs_16p'],Y21a['Zs_84p']-Y21a['Zs_mean']], \
                color=obscol, marker='o', linestyle='', linewidth=1.6, capsize=3, markeredgecolor='black', \
                #color=obscol, marker='o', markersize=obsMarkerSize, linestyle='', linewidth=1., capsize=3, markeredgecolor='black', \
                markersize=obsMarkerSize, label=r'Yates+21a')
    
    #Plot panel label:
    robs_plot_text(ax, r'f)', hpos='right', vpos='top')
    
    pdf.savefig(bbox_inches='tight')     
    pdf.close() 
    

#################
#Plot DTM, DTG, and Mdust versus 12+log(O/H) or [M/H] or M*:
def plot_dust_scaling_relations(Hubble_h, Samp1, struct1, redshift, props='All', xprop='OH', \
                                SolarNorm='A09', levels='3sig', contourOnly=None, contourLines=None, \
                                SFRWeighted=None, incHotDust=None, outlierFrac=defaultOutlierFrac, SFing_only=None) :  
    pdf = PdfPages(PlotDir+xprop+'_vs_'+"dust"+"_"+Samp1['Model']+"_"+Samp1['Sample_type']+"_z"+str(redshift)+".pdf")
    if props == 'All' :
        props=['DTM','DTG','Mdust']
    SolarAbunds = robs_solarnorm_finder(SolarNorm)
    El1 = robs_element_list_finder(struct1)
  
    if contourLines is None :
        ls1 = 'solid'
    else :
        ls1 = contourLines
    
    #Set-up panels:
    figsizescale = 5.
    plt.figure(figsize=(5,len(props)*figsizescale))
    if xprop == 'OH' :
        xlimits = np.array([6.0,9.6])
        xlabel = r'12+log(O/H)'
        xticks = [6.0,7.0,8.0,9.0]
    elif xprop == 'stellarMass' :
        xlimits = np.array([8.0,12.0]) 
        xlabel = r'log$(M_{*} / \textnormal{M}_{\odot})$'
        xticks = [8.0,9.0,10.0,11.0]
    elif xprop == 'MH' :
        xlimits = np.array([-3.0,1.0])
        xlabel = r'[M/H]$_{\rm tot}$'
        xticks = [-3.,-2.,-1.,0.,1.]
    binno = 40
    minInBin = 50
    rowno = len(props)
    colno = 1
    numPanels = rowno*colno
    xlabs = np.full(colno,xlabel)
    xlims = np.full((colno,len(xlimits)),xlimits)
    ylims = np.full((rowno,len(xlimits)),-99.)
    yticks = [None]*rowno
    ylabs = ['y']*rowno
    
    #Calc x-axis values:
    if xprop == 'OH' :
        if SFRWeighted :
            SFRmax = np.max(Samp1['G_samp']['SfrRings'], axis=1) #Max SFR from all rings, for every galaxy
            SFRfrac = Samp1['G_samp']['SfrRings']/SFRmax[:,None] #Fraction of max SFR in each ring, for every galaxy
            SFRtot = np.nansum(SFRfrac, axis=1) #Sum of the SFR fractions from all rings, for evey galaxy. (Minimum must be 1.0)   
            H_rings = Samp1['G_samp']['ColdGasCloudsRings_elements'][:,:,El1["H_NUM"]] + Samp1['G_samp']['ColdGasDiffRings_elements'][:,:,El1["H_NUM"]] \
                    - Samp1['G_samp']['DustColdGasCloudsRings_elements'][:,:,El1["H_NUM"]] - Samp1['G_samp']['DustColdGasDiffRings_elements'][:,:,El1["H_NUM"]]
            O_rings = Samp1['G_samp']['ColdGasCloudsRings_elements'][:,:,El1["O_NUM"]] + Samp1['G_samp']['ColdGasDiffRings_elements'][:,:,El1["O_NUM"]] \
                    - Samp1['G_samp']['DustColdGasCloudsRings_elements'][:,:,El1["O_NUM"]] - Samp1['G_samp']['DustColdGasDiffRings_elements'][:,:,El1["O_NUM"]]
            Zg_rings = (O_rings / H_rings) * (H_aw/O_aw) * SFRfrac
            Zg_tot = np.nansum(Zg_rings, axis=1)
            Zg = 12. + np.log10(Zg_tot / SFRtot)
            theX1 = Zg
        else :
            MH_gas = Samp1['G_samp']['ColdGasClouds_elements'][:,El1["H_NUM"]] + Samp1['G_samp']['ColdGasDiff_elements'][:,El1["H_NUM"]] \
                    - Samp1['G_samp']['DustColdGasClouds_elements'][:,El1["H_NUM"]] - Samp1['G_samp']['DustColdGasDiff_elements'][:,El1["H_NUM"]]
            MO_gas = Samp1['G_samp']['ColdGasClouds_elements'][:,El1["O_NUM"]] + Samp1['G_samp']['ColdGasDiff_elements'][:,El1["O_NUM"]] \
                   - Samp1['G_samp']['DustColdGasClouds_elements'][:,El1["O_NUM"]] - Samp1['G_samp']['DustColdGasDiff_elements'][:,El1["O_NUM"]]
            theX1 = 12. + np.log10((MO_gas/MH_gas) * (H_aw/O_aw))
    elif xprop == 'stellarMass' :
        theX1 = np.log10(Samp1['G_samp']['StellarMass'])
    elif xprop == 'MH' :
        H_mass = Samp1['G_samp']['ColdGas_elements'][:,El1["H_NUM"]]
        Met_mass = np.nansum(Samp1['G_samp']['MetalsColdGas'], axis=1)
        theX1 = np.log10(Met_mass/H_mass) - np.log10(SolarAbunds['ZH_mf_sun'])
    theY1 = np.full((rowno,Samp1["NumGals"]),-99.)
    theWeights1 = 1./Samp1['Volumes']
    
    #Calc y-axis values:
    if incHotDust :
        TotMdust1 = np.nansum(Samp1['G_samp']['DustColdGasClouds_elements'], axis=1) + np.nansum(Samp1['G_samp']['DustColdGasDiff_elements'], axis=1) \
                  + np.nansum(Samp1['G_samp']['DustHotGas_elements'], axis=1) + np.nansum(Samp1['G_samp']['DustEjectedMass_elements'], axis=1)
    else :
        TotMdust1 = np.nansum(Samp1['G_samp']['DustColdGasClouds_elements'], axis=1) + np.nansum(Samp1['G_samp']['DustColdGasDiff_elements'], axis=1)
    for ii in range(len(props)) :
        if props[ii] == 'DTM' :
            ylims[ii] = np.array([-3.4,0.8])
            yticks[ii] = [-3.0,-2.0,-1.0,0.0]
            ylabs[ii] = r'log$(M_{\rm dust} / M_{\rm metal,tot})$'
            TotMmetal1 = np.nansum(Samp1['G_samp']['MetalsColdGas'], axis=1)
            theY1[ii] = np.log10(TotMdust1/TotMmetal1)
            if Samp2 :
                TotMmetal2 = np.nansum(Samp2['G_samp']['MetalsColdGas'], axis=1)
                theY2[ii] = np.log10(TotMdust2/TotMmetal2)
        elif props[ii] == 'DTG' :
            ylims[ii] = np.array([-5.9,-0.5])
            yticks[ii] = [-5.0,-4.0,-3.0,-2.0,-1.0]
            ylabs[ii] = r'log$(M_{\rm dust} /M_{\rm HI+H2})$'
            H_mass1 = Samp1['G_samp']['ColdGas_elements'][:,El1["H_NUM"]]
            theY1[ii] = np.log10(TotMdust1/H_mass1)
            if Samp2 :
                H_mass2 = Samp2['G_samp']['ColdGas_elements'][:,El2["H_NUM"]]
                theY2[ii] = np.log10(TotMdust2/H_mass2)
        elif props[ii] == 'Mdust' :
            ylims[ii] = np.array([1.5,9.5])
            yticks[ii] = [2,3,4,5,6,7,8,9]
            ylabs[ii] = r'log$(M_{\rm dust} / \textnormal{M}_{\odot})$'
            theY1[ii] = np.log10(TotMdust1)
            if Samp2 :
                theY2[ii] = np.log10(TotMdust2)
        
    for kk in range(numPanels) :
        panel = robs_plot_panels(kk, rows=rowno, columns=colno, xlimits=xlims, ylimits=ylims, \
                                 xlab=xlabs, ylab=ylabs, xticks=xticks, yticks=yticks, SecondXAxis=None)
        # Clean model samples:
        if SFing_only :
            tH_z0 = ((1./(100.*Hubble_h))*3.086e+19)*3.17098e-8 #Hubble time at z=0 [in years]
            w_clean1 = np.where((np.isfinite(theX1)) & (~np.isnan(theX1)) & \
                              (np.isfinite(theY1[kk])) & (~np.isnan(theY1[kk])) & \
                              (np.isfinite(theWeights1)) & (~np.isnan(theWeights1)) & \
                              (Samp1['G_samp']['Sfr']/Samp1['G_samp']['StellarMass'] >= ((2.*(1.+robs_redshift_from_snapnum(Samp1['Simulation'], Samp1['Cosmology'], Samp1['G_samp']['SnapNum']))**2)/tH_z0)/10.))
            theX1_clean = theX1[w_clean1[0]]
            theY1_clean = theY1[kk][w_clean1[0]]
            theWeights1_clean = theWeights1[w_clean1[0]]
        else :
            w_clean1 = np.where((np.isfinite(theX1)) & (~np.isnan(theX1)) & \
                              (np.isfinite(theY1[kk])) & (~np.isnan(theY1[kk])) & \
                              (np.isfinite(theWeights1)) & (~np.isnan(theWeights1)))
            theX1_clean = theX1[w_clean1[0]]
            theY1_clean = theY1[kk][w_clean1[0]]
            theWeights1_clean = theWeights1[w_clean1[0]]  
        
        # Plot models:
        robs_contour_plot(theX1_clean, theY1_clean, binno, levels=levels, noOutliers=0, fill=1, alpha=1.0, theRange=[xlimits,ylims[kk]], \
                          linestyle=ls1, colour=model1col, outlierFrac=outlierFrac, weights=theWeights1_clean)
        if not contourOnly :
            robs_plot_average(theX1_clean, theY1_clean, binno=binno, aveType='Median', minInBin=minInBin, linewidth=defaultLinewidth*scaleFactor, \
                              noShade=0, colour=model1col, inputCentres=1, weights=theWeights1_clean, zorder=10)      
    
        #Labels:
        if kk == 0 :
            robs_plot_text(panel, r'z = '+char_z_low, vpos='top', hpos='left')
    
    pdf.savefig(bbox_inches='tight')     
    pdf.close()  
    

#################
#Plot various dust scaling relations back to high redshift:
def plot_dust_scaling_relations_evo(Hubble_h, FullRedshiftList, FullSnapnumList_MRI, FullSnapnumList_MRII, RedshiftsToRead, Samp1, struct1, \
                                    props='All', xprop='OH', SolarNorm='A09', aveOnly=None, contourOnly=None, incHotDust=None, SFRWeighted=None) :  
    pdf = PdfPages(PlotDir+xprop+'_vs_'+"dust_evo"+"_"+Samp1['Model']+"_"+Samp1['Sample_type']+"_z"+Samp1['z_lower']+"-"+Samp1['z_upper']+".pdf")    
    El1 = robs_element_list_finder(struct1)

    #Set-up panels:
    if props == 'All' :
        props = ['DTM','DTG','Mdust','ColdGas','OH','sSFR']
    rowno = len(props)
    colno = len(RedshiftsToRead)
    figsizescale = 3.
    plt.figure(figsize=(4*figsizescale,len(props)*figsizescale)) #plt.figure(figsize=(12,len(props)*figsizescale)) #plt.figure(figsize=(12,12))
    binno = 50
    thelinewidth = defaultLinewidth
    minInBin = 20
    numPanels = rowno*colno
    if xprop == 'OH' :
        xlimits = np.array([6.0,9.6])
        xlabel = r'12+log(O/H)'
        xticks = [6.0,7.0,8.0,9.0]
    elif xprop == 'stellarMass' :
        xlimits = np.array([8.0,12.0]) 
        xlabel = r'log$(M_{*} / \textnormal{M}_{\odot})$'
        xticks = [8.0,9.0,10.0,11.0]
    xlabs = np.full(colno,xlabel)
    xlims = np.full((colno,len(xlimits)),xlimits)
    
    #Set-up arrays:
    ylabs = ['y']*rowno
    ylims = np.full((rowno,len(xlimits)),-99.)
    yticks = [None]*rowno
    if xprop == 'OH' :
        if SFRWeighted :
            SFRmax = np.max(Samp1['G_samp']['SfrRings'], axis=1) #Max SFR from all rings, for every galaxy
            SFRfrac = Samp1['G_samp']['SfrRings']/SFRmax[:,None] #Fraction of max SFR in each ring, for every galaxy
            SFRtot = np.nansum(SFRfrac, axis=1) #Sum of the SFR fractions from all rings, for evey galaxy. (Minimum must be 1.0)   
            H_rings = Samp1['G_samp']['ColdGasCloudsRings_elements'][:,:,El1["H_NUM"]] + Samp1['G_samp']['ColdGasDiffRings_elements'][:,:,El1["H_NUM"]] \
                    - Samp1['G_samp']['DustColdGasCloudsRings_elements'][:,:,El1["H_NUM"]] - Samp1['G_samp']['DustColdGasDiffRings_elements'][:,:,El1["H_NUM"]]
            O_rings = Samp1['G_samp']['ColdGasCloudsRings_elements'][:,:,El1["O_NUM"]] + Samp1['G_samp']['ColdGasDiffRings_elements'][:,:,El1["O_NUM"]] \
                    - Samp1['G_samp']['DustColdGasCloudsRings_elements'][:,:,El1["O_NUM"]] - Samp1['G_samp']['DustColdGasDiffRings_elements'][:,:,El1["O_NUM"]]
            Zg_rings = (O_rings / H_rings) * (H_aw/O_aw) * SFRfrac
            Zg_tot = np.nansum(Zg_rings, axis=1)
            Zg = 12. + np.log10(Zg_tot / SFRtot)
            theX1 = Zg
        else :
            MH_gas = Samp1['G_samp']['ColdGasClouds_elements'][:,El1["H_NUM"]] + Samp1['G_samp']['ColdGasDiff_elements'][:,El1["H_NUM"]] \
                    - Samp1['G_samp']['DustColdGasClouds_elements'][:,El1["H_NUM"]] - Samp1['G_samp']['DustColdGasDiff_elements'][:,El1["H_NUM"]]
            MO_gas = Samp1['G_samp']['ColdGasClouds_elements'][:,El1["O_NUM"]] + Samp1['G_samp']['ColdGasDiff_elements'][:,El1["O_NUM"]] \
                   - Samp1['G_samp']['DustColdGasClouds_elements'][:,El1["O_NUM"]] - Samp1['G_samp']['DustColdGasDiff_elements'][:,El1["O_NUM"]]
            theX1 = 12. + np.log10((MO_gas/MH_gas) * (H_aw/O_aw))
    elif xprop == 'stellarMass' :
        theX1 = np.log10(Samp1['G_samp']['StellarMass'])
    theY1 = np.full((rowno,Samp1["NumGals"]),-99.)
    theWeights1 = 1./Samp1['Volumes']
    
    #Get y-axis values:
    if incHotDust :
        TotMdust1 = np.nansum(Samp1['G_samp']['DustColdGasClouds_elements'], axis=1) + np.nansum(Samp1['G_samp']['DustColdGasDiff_elements'], axis=1) \
                  + np.nansum(Samp1['G_samp']['DustHotGas_elements'], axis=1) + np.nansum(Samp1['G_samp']['DustEjectedMass_elements'], axis=1)
    else :
        TotMdust1 = np.nansum(Samp1['G_samp']['DustColdGasClouds_elements'], axis=1) + np.nansum(Samp1['G_samp']['DustColdGasDiff_elements'], axis=1)
    for ii in range(len(props)) :
        if props[ii] == 'DTM' :
            ylims[ii] = np.array([-3.4,0.8])
            yticks[ii] = [-3.0,-2.0,-1.0,0.0]
            ylabs[ii] = r'log$(M_{\rm dust} / M_{\rm metal,tot})$'
            TotMmetal1 = np.nansum(Samp1['G_samp']['MetalsColdGas'], axis=1)
            theY1[ii] = np.log10(TotMdust1/TotMmetal1)
        elif props[ii] == 'DTG' : 
            ylims[ii] = np.array([-6.4,0.4])
            yticks[ii] = [-6.0,-5.0,-4.0,-3.0,-2.0,-1.0,0.0]
            ylabs[ii] = r'log$(M_{\rm dust} /M_{\rm HI+H2})$'
            H_mass1 = Samp1['G_samp']['ColdGas_elements'][:,El1["H_NUM"]]
            theY1[ii] = np.log10(TotMdust1/H_mass1)
        elif props[ii] == 'Mdust' :
            ylims[ii] = np.array([1.5,9.9])
            yticks[ii] = [2,3,4,5,6,7,8,9]
            ylabs[ii] = r'log$(M_{\rm dust} / \textnormal{M}_{\odot})$'
            theY1[ii] = np.log10(TotMdust1)
        elif props[ii] == 'ColdGas' :
            ylims[ii] = np.array([1.8,12.4])
            yticks[ii] = [2,3,4,5,6,7,8,9,10,11,12]
            ylabs[ii] = r'log$(M_{\rm cold} / \textnormal{M}_{\odot})$'
            theY1[ii] = np.log10(Samp1['G_samp']['ColdGas'])
        elif props[ii] == 'OH' :
            ylims[ii] = np.array([6.0,9.6])
            yticks[ii] = [6.0,7.0,8.0,9.0]
            ylabs[ii] = r'12+log(O/H)'
            if SFRWeighted :
                SFRmax = np.max(Samp1['G_samp']['SfrRings'], axis=1) #Max SFR from all rings, for every galaxy
                SFRfrac = Samp1['G_samp']['SfrRings']/SFRmax[:,None] #Fraction of max SFR in each ring, for every galaxy
                SFRtot = np.nansum(SFRfrac, axis=1) #Sum of the SFR fractions from all rings, for evey galaxy. (Minimum must be 1.0)   
                H_rings = Samp1['G_samp']['ColdGasCloudsRings_elements'][:,:,El1["H_NUM"]] + Samp1['G_samp']['ColdGasDiffRings_elements'][:,:,El1["H_NUM"]] \
                        - Samp1['G_samp']['DustColdGasCloudsRings_elements'][:,:,El1["H_NUM"]] - Samp1['G_samp']['DustColdGasDiffRings_elements'][:,:,El1["H_NUM"]]
                O_rings = Samp1['G_samp']['ColdGasCloudsRings_elements'][:,:,El1["O_NUM"]] + Samp1['G_samp']['ColdGasDiffRings_elements'][:,:,El1["O_NUM"]] \
                        - Samp1['G_samp']['DustColdGasCloudsRings_elements'][:,:,El1["O_NUM"]] - Samp1['G_samp']['DustColdGasDiffRings_elements'][:,:,El1["O_NUM"]]
                Zg_rings = (O_rings / H_rings) * (H_aw/O_aw) * SFRfrac
                Zg_tot = np.nansum(Zg_rings, axis=1)
                theY1[ii] = 12. + np.log10(Zg_tot / SFRtot)
            else :
                MH_gas = Samp1['G_samp']['ColdGasClouds_elements'][:,El1["H_NUM"]] + Samp1['G_samp']['ColdGasDiff_elements'][:,El1["H_NUM"]] \
                        - Samp1['G_samp']['DustColdGasClouds_elements'][:,El1["H_NUM"]] - Samp1['G_samp']['DustColdGasDiff_elements'][:,El1["H_NUM"]]
                MO_gas = Samp1['G_samp']['ColdGasClouds_elements'][:,El1["O_NUM"]] + Samp1['G_samp']['ColdGasDiff_elements'][:,El1["O_NUM"]] \
                       - Samp1['G_samp']['DustColdGasClouds_elements'][:,El1["O_NUM"]] - Samp1['G_samp']['DustColdGasDiff_elements'][:,El1["O_NUM"]]
                theY1[ii] = 12. + np.log10((MO_gas/MH_gas) * (H_aw/O_aw))
        elif props[ii] == 'sSFR' :
            ylims[ii] = np.array([-15.0,-6.5])
            yticks[ii] = [-15.,-14.,-13.,-12.,-11.,-10.,-9.,-8.,-7.]
            ylabs[ii] = r'log$(\textnormal{sSFR} / \textnormal{yr}^{-1})$'
            theY1[ii] = np.log10(Samp1['G_samp']['Sfr'] / Samp1['G_samp']['StellarMass'])
                
    #Plot panels:
    for kk in range(numPanels) :
        panel = robs_plot_panels(kk, rows=rowno, columns=colno, xlimits=xlims, ylimits=ylims, \
                                 xlab=xlabs, ylab=ylabs, xticks=xticks, yticks=yticks, SecondXAxis=None)
        theRowNum = int(kk/colno)
        theColNum = kk%colno
        if (Samp1['MRI_gals'] > 0) & (Samp1['MRII_gals'] > 0) :
            wz1 = np.where(((Samp1['Volumes'] == Samp1['Volumes'][0]) & (Samp1['G_samp']['SnapNum'] == FullSnapnumList_MRI[theColNum]) & (np.isfinite(theY1[theRowNum])) & (~np.isnan(theY1[theRowNum]))) | \
                           ((Samp1['Volumes'] == Samp1['Volumes'][Samp1['MRI_gals']]) & (Samp1['G_samp']['SnapNum'] == FullSnapnumList_MRII[theColNum]) & (np.isfinite(theY1[theRowNum])) & (~np.isnan(theY1[theRowNum]))))
        elif Samp1['Simulation'] == 'Mil-I' :
            wz1 = np.where((Samp1['G_samp']['SnapNum'] == FullSnapnumList_MRI[theColNum]) & (np.isfinite(theY1[theRowNum])) & (~np.isnan(theY1[theRowNum])))
        elif Samp1['Simulation'] == 'Mil-II' :
            wz1 = np.where((Samp1['G_samp']['SnapNum'] == FullSnapnumList_MRII[theColNum]) & (np.isfinite(theY1[theRowNum])) & (~np.isnan(theY1[theRowNum])))

        if aveOnly :
            robs_plot_average(theX1[wz1], theY1[theRowNum][wz1], binno=binno, aveType='Median', minInBin=minInBin, linewidth=defaultLinewidth*scaleFactor, \
                              noShade=0, colour=model1col, inputCentres=1, levels='2sig', weights=theWeights1[wz1])#, plotlims=[xlimits,ylims[theRowNum]])
            
        else :
            robs_contour_plot(theX1[wz1], theY1[theRowNum][wz1], binno, levels='4sig', noOutliers=0, fill=1, alpha=0.3, theRange=[xlimits,ylims[theRowNum]], \
                              colour=model1col, outlierFrac=defaultOutlierFrac, weights=theWeights1[wz1])
            if not contourOnly :
                robs_plot_average(theX1[wz1], theY1[theRowNum][wz1], binno=binno, aveType='Median', minInBin=minInBin, linewidth=defaultLinewidth*scaleFactor, \
                                  noShade=0, colour=model1col, inputCentres=1, levels='2sig', weights=theWeights1[wz1])
          
        #Plot SFing cut-off line: N.B. This cut-off seems to work nicely for Mil-I and Mil-II right back to z~10 (19-06-23):
        if props[theRowNum] == 'sSFR' :
            tH_z0 = ((1./(100.*Hubble_h))*3.086e+19)*3.17098e-8 #Hubble time at z=0 [in years]
            logsSFR_cutoff = np.log10((2. * (1. + FullRedshiftList[theColNum])**2/tH_z0)/10.)
            plt.plot(xlimits, [logsSFR_cutoff,logsSFR_cutoff], '--', linewidth=thelinewidth, color='red') #, label=phase[ii])            
   
        #Labels:
        textPadder = 15. #20.
        xpadder = (xlimits[1]-xlimits[0])/textPadder
        ypadder = (ylims[theRowNum][1]-ylims[theRowNum][0])/textPadder
        if theRowNum == 0 :
            if xprop == 'OH' :
                plt.text(xlimits[0]+xpadder, ylims[theRowNum][1]-ypadder, r'z = '+str(FullRedshiftList[theColNum]), horizontalalignment='left', \
                         verticalalignment='top', color='black') #r'z = '+Samp1['z_lower']
            elif xprop == 'stellarMass' :
                plt.text(xlimits[1]-(1.25*xpadder), ylims[theRowNum][0]+ypadder, r'z = '+str(FullRedshiftList[theColNum]), horizontalalignment='right', \
                         verticalalignment='bottom', color='black')
        
    pdf.savefig(bbox_inches='tight')     
    pdf.close()  
    
    
#################
#Plot evolution of dust production/destruction rate densities and dust mass densities:
def plot_cosmic_dust_evos(Volume, Hubble_h, Omega_M, Omega_b, FullRedshiftList, FullSnapnumList, \
                          RedshiftsToRead, Samp1, plotTotal=None) :
    pdf = PdfPages(PlotDir+"cosmic_dust_evos"+"_"+Samp1['Model']+"_"+Samp1['Sample_type']+"_z"+Samp1['z_lower']+"-"+Samp1['z_upper']+".pdf")
    NumSamps = 1
    
    # Set-up axes:
    xlimits = [Samp1['z_lower']-0.2, Samp1['z_upper']+0.2] #[0.0,8.0]
    ylimits = [[-9.4,-2.2],[0.4,6.2]]
    rowno=2
    colno=1
    xlabs = [r'Redshift']
    ylabs = [r'log$(\dot{\rho}_{\rm dust}\ /\ \textnormal{M}_{\odot}\ \textnormal{yr}^{-1}\ \textnormal{cMpc}^{-3})$', \
             r'log$(\rho_{\rm dust}\ /\ \textnormal{M}_{\odot}\ \textnormal{cMpc}^{-3})$']
    xlims = xlimits
    ylims = ylimits
    xticks = [0,1,2,3,4,5,6,7,8]
    yticks = [[-9,-8,-7,-6,-5,-4,-3],[1.,2.,3.,4.,5.,6.]]
    fig = plt.figure(figsize=(5,10))
    
    # Set cols and processes:
    thecols1 = ['blue','orange','green','red','cyan','purple','pink']
    theLabs1 = ['AGB (total)              ','SN-II (total)','SN-Ia (total)','Grain growth','SN shocks + astration (ISM)','Sputtering (CGM)','Sputtering (Ejecta)']
    thecols2 = ['darkblue','cornflowerblue','green','pink']
    theLabs2 = ['ISM (mol. clouds)         ', 'ISM (diffuse gas)', 'CGM', 'Ejecta']  
    zorders2 = [10,-1,-1,-1]
    allLabs = theLabs1 + theLabs2
    
    # X-axis values:
    FRL = np.array(FullRedshiftList)
    RTR = np.array(RedshiftsToRead)
    
    # Calc cosmic densities:
    logDRate = np.empty((NumSamps,len(FullRedshiftList),len(theLabs1))) #Each dust prod/dest rate at each redshift
    logTotProdDRate = np.empty((NumSamps,len(FullRedshiftList)))
    logTotDestDRate = np.empty((NumSamps,len(FullRedshiftList)))
    tot_prod_rates = np.zeros((NumSamps,len(FullRedshiftList)))
    tot_dest_rates = np.zeros((NumSamps,len(FullRedshiftList)))
    logDMass = np.empty((NumSamps,len(FullRedshiftList),len(theLabs2)))
    logTotDMass = np.empty((NumSamps,len(FullRedshiftList)))
    tot_dust_mass = np.zeros((NumSamps,len(FullRedshiftList)))
    for iz in range(len(FullRedshiftList)) : #loop over redshifts
        if RedshiftsToRead[iz] : 
            wz1 = np.where(Samp1['G_samp']['SnapNum'] == FullSnapnumList[iz])
            # Dust rates:
            for ii in range(len(theLabs1)) : #loop over rates
                if (ii <= 2) : #AGB, SN-II, SN-Ia
                    Rate = np.nansum(Samp1['G_samp']['DustColdGasRates'][wz1,ii] + Samp1['G_samp']['DustHotGasRates'][wz1,ii])
                    tot_prod_rates[0,iz] += Rate
                elif (ii == 3) : #Grain growth      
                    Rate = np.nansum(Samp1['G_samp']['DustColdGasRates'][wz1,ii])
                    tot_prod_rates[0,iz] += Rate
                elif (ii == 4) : #SN dest      
                    Rate = np.nansum(Samp1['G_samp']['DustColdGasRates'][wz1,ii])
                    tot_dest_rates[0,iz] += Rate
                elif (ii == 5) :  #Sputtering dest (HotGas)
                    Rate = np.nansum(Samp1['G_samp']['DustHotGasRates'][wz1,ii-2])
                    tot_dest_rates[0,iz] += Rate
                elif (ii == 6) : #Sputtering dest (EjectedMass)
                    Rate = np.nansum(Samp1['G_samp']['DustEjectedMassRates'][wz1])
                    tot_dest_rates[0,iz] += Rate
                logDRate[0,iz,ii] = np.log10(Rate / (Volume/Hubble_h**3))
            # Dust masses:
            for ii in range(len(theLabs2)) : #loop over rates
                if (ii == 0) : #Clouds (ISM)
                    Masses = Samp1['G_samp']['DustColdGasClouds_elements'][wz1,:]
                elif (ii == 1) : #Diffuse (ISM)
                    Masses = Samp1['G_samp']['DustColdGasDiff_elements'][wz1,:]
                elif (ii == 2) : #CGM     
                    Masses = Samp1['G_samp']['DustHotGas_elements'][wz1,:]
                elif (ii == 3) : #Ejecta     
                    Masses = Samp1['G_samp']['DustEjectedMass_elements'][wz1,:]
                Masses[np.where(Masses < 0.0)] = 0.0
                Mass = np.nansum(Masses)
                tot_dust_mass[0,iz] += Mass
                logDMass[0,iz,ii] = np.log10(Mass / (Volume/Hubble_h**3)) 
    logTotProdDRate = np.log10(tot_prod_rates / (Volume/Hubble_h**3))
    logTotDestDRate = np.log10(tot_dest_rates / (Volume/Hubble_h**3))
    logTotDMass = np.log10(tot_dust_mass / (Volume/Hubble_h**3))
    
    # Plot panels:
    for ii in range(rowno*colno) :
        panel = robs_plot_panels(ii, rows=rowno, columns=colno, xlimits=xlims, ylimits=ylims, \
                                 xlab=xlabs, ylab=ylabs, xticks=xticks, yticks=yticks, SecondXAxis=None)
        #Add second x axis:
        thezs = np.array([0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0])
        thelbts = np.array([0.0,4.0,8.0,10.0,11.0,12.0,13.0])
        thelbts_zs = robs_get_lookbacktimes(Hubble_h, Omega_M, Omega_b, thezs=thezs, thelbts=thelbts)
        p1b = panel.twiny()
        p1b.set_xlim(xlimits[0],xlimits[1])
        p1b.set_xticks(thelbts_zs, minor=False) 
        p1b.tick_params(direction="in") 
        p1b.set_xticks([0], minor=True)
        
        #####################
        # Cosmic dust rates:
        if ii == 0 :
            # Finish-off second x axis:
            p1b.set_xticklabels(['{:g}'.format(lbt) for lbt in thelbts])
            p1b.set_xlabel(r'Lookback time / Gyr')
            plt.setp(panel.get_xticklabels(), visible=False)
                  
            #Plot cosmic rate densities:
            for ii in range(len(theLabs1)) : #loop over rates
                coeffs = np.polyfit(FRL[RTR], logDRate[0,:,ii], 5)
                poly = np.poly1d(coeffs)
                Samp1_x = np.linspace(FRL[RTR][0], FRL[RTR][-1])
                Samp1_y = poly(Samp1_x)
                plt.plot(Samp1_x, Samp1_y, linewidth=defaultLinewidth, linestyle='solid', color=thecols1[ii], zorder=-1)
                  
            # Plot tot prod and dest rates:
            if plotTotal :
                linstyle = ['solid','dashed']
                for samp in range(NumSamps) :
                    coeffs = np.polyfit(FRL[RTR], logTotProdDRate[samp,:], 5)
                    poly = np.poly1d(coeffs)
                    Samp1_x = np.linspace(FRL[RTR][0], FRL[RTR][-1])
                    Samp1_y = poly(Samp1_x)
                    plt.plot(Samp1_x, Samp1_y, linewidth=2.*defaultLinewidth, linestyle=linstyle[samp], color='black', zorder=-2)
                    coeffs = np.polyfit(FRL[RTR], logTotDestDRate[samp,:], 5)
                    poly = np.poly1d(coeffs)
                    Samp1_x = np.linspace(FRL[RTR][0], FRL[RTR][-1])
                    Samp1_y = poly(Samp1_x)
                    plt.plot(Samp1_x, Samp1_y, linewidth=2.*defaultLinewidth, linestyle=linstyle[samp], color='grey', zorder=-2)
   
            # Labels:
            numLeft = 4
            padder=0.05
            # Production labs:
            theProdLabs = theLabs1[0:numLeft]
            theProdCols = thecols1[0:numLeft]
            theProdLabs = [r'\textbf{Dust production/growth:}'] + theProdLabs
            theProdCols =  ['black'] + theProdCols
            if plotTotal :
                theProdLabs += ['Total']
                theProdCols += ['black']
            robs_plot_text(panel, theProdLabs, vpos='top', hpos='right', padder=padder, colour=theProdCols)
            # Destruction labs:
            theDestLabs = list(reversed(theLabs1[numLeft:]))
            theDestCols = list(reversed(thecols1[numLeft:]))
            theDestLabs = [r'\textbf{Dust destruction:}'] + theDestLabs
            theDestCols =  ['black'] + theDestCols
            if plotTotal :
                theDestLabs += ['Total']
                theDestCols += ['grey']
            robs_plot_text(panel, theDestLabs, vpos='bottom', hpos='left', padder=padder, colour=theDestCols)
                
        #####################
        # Cosmic dust masses:
        elif ii == 1 :
            # Finish-off second x axis:
            p1b.set_xticklabels([])
            
            #Plot cosmic mass densities:
            for ii in range(len(theLabs2)) : #loop over rates
                coeffs = np.polyfit(FRL[RTR], logDMass[0,:,ii], 5)
                poly = np.poly1d(coeffs)
                Samp1_x = np.linspace(FRL[RTR][0], FRL[RTR][-1])
                Samp1_y = poly(Samp1_x)
                plt.plot(Samp1_x, Samp1_y, linewidth=defaultLinewidth, linestyle='solid', color=thecols2[ii], zorder=zorders2[ii])
            
            # Plot MCs + diffuse:
            coeffs = np.polyfit(FRL[RTR], np.log10(10**(logDMass[0,:,0]) + 10**(logDMass[0,:,1])), 5)
            poly = np.poly1d(coeffs)
            Samp1_x = np.linspace(FRL[RTR][0], FRL[RTR][-1])
            Samp1_y = poly(Samp1_x)
            plt.plot(Samp1_x, Samp1_y, linewidth=defaultLinewidth, linestyle='solid', color='red', zorder=10)
   
            # Plot tot prod and dest masses:
            if plotTotal :
                linstyle = ['solid','dashed']
                for samp in range(NumSamps) :
                    coeffs = np.polyfit(FRL[RTR], logTotDMass[samp,:], 5)
                    poly = np.poly1d(coeffs)
                    Samp1_x = np.linspace(FRL[RTR][0], FRL[RTR][-1])
                    Samp1_y = poly(Samp1_x)
                    plt.plot(Samp1_x, Samp1_y, linewidth=2.*defaultLinewidth, linestyle=linstyle[samp], color='black', zorder=-2)
            
            # Labels:
            theLabs2 = list(reversed(['ISM (total)'] + theLabs2))
            thecols2 = list(reversed(['red'] + thecols2))
            theLabs2 = [r'\textbf{Dust mass:}'] + theLabs2
            thecols2 =  ['black'] + thecols2
            if plotTotal :
                theLabs2 += ['Total']
                thecols2 += ['black']
            padder=0.05
            robs_plot_text(panel, theLabs2, vpos='bottom', hpos='left', padder=padder, colour=thecols2)    
            robs_plot_text(panel, r'$H_{0} =$ '+str(round(Hubble_h*100.,1)), vpos='top', hpos='right', padder=3.3*padder)
            
            #Add second y-axis:
            PH20_h = 0.7
            PH20_rho_crit = 1.36e11 * ((100.*Hubble_h)**2/(100.*PH20_h)**2) #Msun/cMpc^3 Critical mass density of the Universe (PH20, section 2.1.1), corrected for Hubble parameter. N.B. rho_crit,0 = 3*H0^2 / 8*pi*G
            p2b = panel.twinx()
            p2b.set_ylim(panel.get_ylim() - np.log10(PH20_rho_crit))
            p2b.set_ylabel(r'log($\Omega_{\rm dust}$)')
            p2b.tick_params(direction="in")

    pdf.savefig(bbox_inches='tight')     
    pdf.close()
    
    
#################
def plot_profs_multiplot(Samp1, struct1, Hubble_h, Omega_M, Omega_Lambda, \
                         MassBins, props='All', SolarNorm=None, stellarComp='Disc+Bulge') :
    pdf = PdfPages(PlotDir+"Profs_multiplot"+"_"+Samp1['Model']+"_"+Samp1['Sample_type']+"_z"+Samp1['z_lower']+".pdf")
    SolarAbunds = robs_solarnorm_finder(SolarNorm)
    El1 = robs_element_list_finder(struct1)
    if props == 'All' :
        props = ['OH_cold','Z_stars','FeH_stars','OFE_cold','OFE_stars','DTM','CvsD_DTM','DTG','DMD', \
                 'DCloudsD','DDiffD','CloudsD','DiffD','MmetD','MHD','MHID','MH2D','SFRD','MstarD']

    #Set-up panels:
    figsizescale = 3.
    figsizescale_x = 4.
    plt.figure(figsize=(len(MassBins)*figsizescale_x,len(props)*figsizescale))
    xlimits = np.array([0.0,5.0])
    xlabel = r'$R/\textnormal{R}_{\textnormal{e}}$'
    binwidth = 0.2
    binno = 50
    thelinewidth = defaultLinewidth
    rowno = len(props)
    colno = len(MassBins)
    numPanels = rowno*colno
    xlabs = np.full(colno,xlabel)
    xlims = np.full((colno,len(xlimits)),xlimits)
    theX = Samp1['RReRings'][:,:] 
    ylims = np.full((rowno,len(xlimits)),-99.)
    yticks = [None]*rowno
    ylabs = ['y']*rowno
    theY = np.full((rowno,Samp1["NumGals"],Samp1["RNUM"]),-99.)
    theY1a = np.full((rowno,Samp1["NumGals"],Samp1["RNUM"]),-99.)
    theY1b = np.full((rowno,Samp1["NumGals"],Samp1["RNUM"]),-99.)
    
    #Get y-axis values:
    for ii in range(len(props)) :
        if props[ii] == 'OH_cold' :
            ylims[ii] = np.array([7.51,9.6])
            yticks[ii] = [8.0,8.5,9.0,9.5]
            ylabs[ii] = r'12+log(O/H)$_{\rm gas}$'
            MH_gas = Samp1['G_samp']['ColdGasCloudsRings_elements'][:,:,El1["H_NUM"]] + Samp1['G_samp']['ColdGasDiffRings_elements'][:,:,El1["H_NUM"]] \
                   - Samp1['G_samp']['DustColdGasCloudsRings_elements'][:,:,El1["H_NUM"]] - Samp1['G_samp']['DustColdGasDiffRings_elements'][:,:,El1["H_NUM"]]
            MO_gas = Samp1['G_samp']['ColdGasCloudsRings_elements'][:,:,El1["O_NUM"]] + Samp1['G_samp']['ColdGasDiffRings_elements'][:,:,El1["O_NUM"]] \
                   - Samp1['G_samp']['DustColdGasCloudsRings_elements'][:,:,El1["O_NUM"]] - Samp1['G_samp']['DustColdGasDiffRings_elements'][:,:,El1["O_NUM"]]
            theY[ii] = 12. + np.log10((MO_gas / MH_gas) * (H_aw/O_aw))
        elif props[ii] == 'Z_stars' :
            ylims[ii] = np.array([-1.7,0.7])
            yticks[ii] = [-1.5,-1.0,-0.5,0.0,0.5]
            if stellarComp == 'Disc' :
                ylabs[ii] = r'log($Z_{*\rm ,disc}/\textnormal{Z}_{\odot}$)'
                TotMetRings = np.nansum(Samp1['G_samp']['MetalsDiskMassRings'][:,:,:], axis=2)
                theY[ii] = np.log10(TotMetRings / Samp1['G_samp']['DiskMassRings'][:,:]) - np.log10(SolarAbunds['Z_mf_sun'])
            elif (stellarComp == 'Disc+Bulge') | (stellarComp == 'Bulge+Disc') :
                ylabs[ii] = r'log($Z_{*}/\textnormal{Z}_{\odot}$)'
                TotMetRings = np.nansum(Samp1['G_samp']['MetalsDiskMassRings'][:,:,:], axis=2) + np.nansum(Samp1['G_samp']['MetalsBulgeMassRings'][:,:,:], axis=2)
                theY[ii] = np.log10(TotMetRings / (Samp1['G_samp']['DiskMassRings'][:,:] + Samp1['G_samp']['BulgeMassRings'][:,:])) - np.log10(SolarAbunds['Z_mf_sun'])
            elif stellarComp == 'Bulge' :
                ylabs[ii] = r'log($Z_{*\rm ,bulge}/\textnormal{Z}_{\odot}$)'
                TotMetRings = np.nansum(Samp1['G_samp']['MetalsBulgeMassRings'][:,:,:], axis=2)
                theY[ii] = np.log10(TotMetRings / Samp1['G_samp']['BulgeMassRings'][:,:]) - np.log10(SolarAbunds['Z_mf_sun'])
        elif props[ii] == 'FeH_stars' :
            ylims[ii] = np.array([-1.8,0.4])
            yticks[ii] = [-1.5,-1.0,-0.5,0.0]
            if stellarComp == 'Disc' :
                ylabs[ii] = r'[Fe/H]$_{*\rm ,disc}$'
                HRings = Samp1['G_samp']['DiskMassRings_elements'][:,:,El1["H_NUM"]]
                FeRings = Samp1['G_samp']['DiskMassRings_elements'][:,:,El1["Fe_NUM"]]
            elif (stellarComp == 'Disc+Bulge') | (stellarComp == 'Bulge+Disc') :
                ylabs[ii] = r'[Fe/H]$_{*}$'
                HRings = Samp1['G_samp']['DiskMassRings_elements'][:,:,El1["H_NUM"]] + Samp1['G_samp']['BulgeMassRings_elements'][:,:,El1["H_NUM"]]
                FeRings = Samp1['G_samp']['DiskMassRings_elements'][:,:,El1["Fe_NUM"]] + Samp1['G_samp']['BulgeMassRings_elements'][:,:,El1["Fe_NUM"]]
            elif stellarComp == 'Bulge' :
                ylabs[ii] = r'[Fe/H]$_{*\rm ,bulge}$'
                HRings = Samp1['G_samp']['BulgeMassRings_elements'][:,:,El1["H_NUM"]]
                FeRings = Samp1['G_samp']['BulgeMassRings_elements'][:,:,El1["Fe_NUM"]]
            theY[ii] = np.log10(FeRings/HRings) - np.log10(SolarAbunds['FeH_mf_sun'])
        elif props[ii] == 'OFE_cold' :
            ylims[ii] = np.array([0.08,0.46])
            yticks[ii] = [0.1,0.2,0.3,0.4]
            ylabs[ii] = r'[O/Fe]$_{\rm ISM,gas+dust}$'
            ORings = Samp1['G_samp']['ColdGasRings_elements'][:,:,El1["O_NUM"]]
            FeRings = Samp1['G_samp']['ColdGasRings_elements'][:,:,El1["Fe_NUM"]]
            theY[ii] = np.log10(ORings/FeRings) - np.log10(SolarAbunds['OH_mf_sun']/SolarAbunds['FeH_mf_sun'])
        elif props[ii] == 'OFE_stars' :
            ylims[ii] = np.array([0.15,0.51])
            yticks[ii] = [0.2,0.3,0.4,0.5]
            if (stellarComp == 'Disc+Bulge') | (stellarComp == 'Bulge+Disc') :
                ylabs[ii] = r'[O/Fe]$_{*\rm ,disc+bulge}$'
                ORings = Samp1['G_samp']['DiskMassRings_elements'][:,:,El1["O_NUM"]] + Samp1['G_samp']['BulgeMassRings_elements'][:,:,El1["O_NUM"]]
                FeRings = Samp1['G_samp']['DiskMassRings_elements'][:,:,El1["Fe_NUM"]] + Samp1['G_samp']['BulgeMassRings_elements'][:,:,El1["Fe_NUM"]]
            elif stellarComp == 'Disc' :
                ylabs[ii] = r'[O/Fe]$_{*\rm ,disc}$'
                ORings = Samp1['G_samp']['DiskMassRings_elements'][:,:,El1["O_NUM"]]
                FeRings = Samp1['G_samp']['DiskMassRings_elements'][:,:,El1["Fe_NUM"]]
            elif stellarComp == 'Bulge' :
                ylabs[ii] = r'[O/Fe]$_{*\rm ,bulge}$'
                ORings = Samp1['G_samp']['BulgeMassRings_elements'][:,:,El1["O_NUM"]]
                FeRings = Samp1['G_samp']['BulgeMassRings_elements'][:,:,El1["Fe_NUM"]]
            theY[ii] = np.log10(ORings/FeRings) - np.log10(SolarAbunds['OH_mf_sun']/SolarAbunds['FeH_mf_sun'])  
        elif props[ii] == 'DTM' :
            ylims[ii] = np.array([-2.7,0.0])
            yticks[ii] = [-2.5,-2.0,-1.5,-1.0,-0.5]
            ylabs[ii] = r'log$(M_{\rm dust} /M_{\rm metal,tot})$'
            TotMdustRings = np.nansum(Samp1['G_samp']['DustColdGasCloudsRings_elements'], axis=2) + np.nansum(Samp1['G_samp']['DustColdGasDiffRings_elements'], axis=2)
            TotMmetalRings = np.nansum(Samp1['G_samp']['MetalsColdGasRings'], axis=2)
            theY[ii] = np.log10(TotMdustRings/TotMmetalRings)
        elif props[ii] == 'DTG' :
            ylims[ii] = np.array([-5.2,-0.8])
            yticks[ii] = [-5.,-4.,-3.,-2.,-1.]
            ylabs[ii] = r'log$(M_{\rm dust} /M_{\rm HI+H2})$'
            TotMdustRings = np.nansum(Samp1['G_samp']['DustColdGasCloudsRings_elements'], axis=2) + np.nansum(Samp1['G_samp']['DustColdGasDiffRings_elements'], axis=2)
            TotHMassRings = Samp1['G_samp']['ColdGasRings_elements'][:,:,El1["H_NUM"]]
            theY[ii] = np.log10(TotMdustRings/TotHMassRings)             
        elif props[ii] == 'DMD' :
            ylims[ii] = np.array([-5.5,2.5])
            yticks[ii] = [-5.,-4.,-3.,-2.,-1.,0.,1.,2.]
            ylabs[ii] = r'log($\Sigma_{\rm dust} / \textnormal{M}_{\odot}\textnormal{ pc}^{-2}$)'
            TotMdustRings = np.nansum(Samp1['G_samp']['DustColdGasCloudsRings_elements'], axis=2) + np.nansum(Samp1['G_samp']['DustColdGasDiffRings_elements'], axis=2)
            theY[ii] = np.log10(TotMdustRings/Samp1['RingArea']/(1000.**2)) #to get Msun/pc^2  
        elif props[ii] == 'DCloudsD' :
            ylims[ii] = np.array([-5.5,2.5])
            yticks[ii] = [-5.,-4.,-3.,-2.,-1.,0.,1.,2.]
            ylabs[ii] = r'log($\Sigma_{\rm dust,clouds} / \textnormal{M}_{\odot}\textnormal{ pc}^{-2}$)'
            theY[ii] = np.log10(np.nansum(Samp1['G_samp']['DustColdGasCloudsRings_elements'], axis=2)/Samp1['RingArea']/(1000.**2)) #to get Msun/pc^2  
        elif props[ii] == 'DDiffD' :
            ylims[ii] = np.array([-5.5,2.5])
            yticks[ii] = [-5.,-4.,-3.,-2.,-1.,0.,1.,2.]
            ylabs[ii] = r'log($\Sigma_{\rm dust,diffuse} / \textnormal{M}_{\odot}\textnormal{ pc}^{-2}$)'
            theY[ii] = np.log10(np.nansum(Samp1['G_samp']['DustColdGasDiffRings_elements'], axis=2)/Samp1['RingArea']/(1000.**2)) #to get Msun/pc^2   
        elif props[ii] == 'CloudsD' :
            ylims[ii] = np.array([-3.5,2.5])
            yticks[ii] = [-5.,-4.,-3.,-2.,-1.,0.,1.,2.]
            ylabs[ii] = r'log($\Sigma_{\rm clouds} / \textnormal{M}_{\odot}\textnormal{ pc}^{-2}$)'
            theY[ii] = np.log10(np.nansum(Samp1['G_samp']['ColdGasCloudsRings_elements'], axis=2)/Samp1['RingArea']/(1000.**2)) #to get Msun/pc^2  
        elif props[ii] == 'DiffD' :
            ylims[ii] = np.array([-1.5,2.5])
            yticks[ii] = [-5.,-4.,-3.,-2.,-1.,0.,1.,2.]
            ylabs[ii] = r'log($\Sigma_{\rm diffuse} / \textnormal{M}_{\odot}\textnormal{ pc}^{-2}$)'
            theY[ii] = np.log10(np.nansum(Samp1['G_samp']['ColdGasDiffRings_elements'], axis=2)/Samp1['RingArea']/(1000.**2)) #to get Msun/pc^2  
        elif props[ii] == 'CvsD_DTM' :
            ylims[ii] = np.array([-2.7,0.0])
            yticks[ii] = [-2.5,-2.0,-1.5,-1.0,-0.5]
            ylabs[ii] = r'log$(M_{\rm dust} /M_{\rm metal,tot})$'
            DustMassClouds = np.nansum(Samp1['G_samp']['DustColdGasCloudsRings_elements'], axis=2)
            MetalMassClouds = np.nansum(Samp1['G_samp']['ColdGasCloudsRings_elements'][:,:,El1["He_NUM"]+1:], axis=2) #2:11
            theY1a[ii] = np.log10(DustMassClouds/MetalMassClouds)
            DustMassDiff = np.nansum(Samp1['G_samp']['DustColdGasDiffRings_elements'], axis=2)
            MetalMassDiff = np.nansum(Samp1['G_samp']['ColdGasDiffRings_elements'][:,:,El1["He_NUM"]+1:], axis=2)
            theY1b[ii] = np.log10(DustMassDiff/MetalMassDiff)
        elif props[ii] == 'MmetD' :
            ylims[ii] = np.array([-3.5,1.5])
            yticks[ii] = [-3.,-2.,-1.,0.,1.]
            ylabs[ii] = r'log($\Sigma_{\rm met,tot} / \textnormal{M}_{\odot}\textnormal{ pc}^{-2}$)'
            theY[ii] = np.log10(np.nansum(Samp1['G_samp']['MetalsColdGasRings'], axis=2)/Samp1['RingArea']/(1000.**2)) #to get Msun/pc^2  
        elif props[ii] == 'MHD' :
            ylims[ii] = np.array([-1.5,3.5])
            yticks[ii] = [-1.,0.,1.,2.,3.]
            ylabs[ii] = r'log($\Sigma_{\rm H,tot} / \textnormal{M}_{\odot}\textnormal{ pc}^{-2}$)'
            theY[ii] = np.log10(Samp1['G_samp']['ColdGasRings_elements'][:,:,El1["H_NUM"]]/Samp1['RingArea']/(1000.**2)) #to get Msun/pc^2  
        elif props[ii] == 'MHID' :
            ylims[ii] = np.array([-0.9,1.7])
            yticks[ii] = [-0.5,0.,0.5,1.,1.5]
            ylabs[ii] = r'log($\Sigma_{\rm HI} / \textnormal{M}_{\odot}\textnormal{ pc}^{-2}$)'
            theY[ii] = np.log10((Samp1['G_samp']['ColdGasRings_elements'][:,:,El1["H_NUM"]]*(1. - Samp1['G_samp']['H2fractionRings']))/Samp1['RingArea']/(1000.**2)) #to get Msun/pc^2  
        elif props[ii] == 'MH2D' :
            ylims[ii] = np.array([-2.4,3.5])
            yticks[ii] = [-2.,-1.,0.,1.,2.,3.]
            ylabs[ii] = r'log($\Sigma_{\rm H2} / \textnormal{M}_{\odot}\textnormal{ pc}^{-2}$)'
            theY[ii] = np.log10((Samp1['G_samp']['ColdGasRings_elements'][:,:,El1["H_NUM"]]*Samp1['G_samp']['H2fractionRings'])/Samp1['RingArea']/(1000.**2)) #to get Msun/pc^2  
        elif props[ii] == 'SFRD' :
            ylims[ii] = np.array([-6.5,0.5])
            yticks[ii] = [-6.,-5.,-4.,-3.,-2.,-1.,0.]
            ylabs[ii] = r'log($\Sigma_{\rm SFR} / \textnormal{ M}_{\odot}\textnormal{ yr}^{-1}\textnormal{ kpc}^{-2}$)'
            theY[ii] = np.log10(Samp1['G_samp']['SfrRings']/Samp1['RingArea']) #in Msun/kpc^2    /(1000.**2)) #to get Msun/pc^2                       
        elif props[ii] == 'MstarD' :
            ylims[ii] = np.array([-2.4,3.9])
            yticks[ii] = [-2.,-1.,0.,1.,2.,3.]
            ylabs[ii] = r'log($\Sigma_{*} / \textnormal{M}_{\odot}\textnormal{ pc}^{-2}$)'
            theY[ii] = np.log10((Samp1['G_samp']['DiskMassRings'] + Samp1['G_samp']['BulgeMassRings'])/Samp1['RingArea']/(1000.**2)) #to get Msun/pc^2   
        
    for kk in range(numPanels) :
        panel = robs_plot_panels(kk, rows=rowno, columns=colno, xlimits=xlims, ylimits=ylims, \
                                 xlab=xlabs, ylab=ylabs, xticks=[0,1,2,3,4], yticks=yticks, SecondXAxis=None)
        if kk >= numPanels-colno :
            thePanelLabel = r'log$(M_{*}/\textnormal{M}_{\odot}) =$ '
        else :
            thePanelLabel = None
        theRowNum = int(kk/colno)
        theColNum = kk%colno
        panelProp = np.log10(Samp1['G_samp']['StellarMass'])
        xrange = [theX[np.isfinite(theX)].min(), theX[np.isfinite(theX)].max()]
        if props[theRowNum] == 'CvsD_DTM' :
            robs_plot_profile(theX, theY1a[theRowNum], panelProp, xlims[theColNum], ylims[theRowNum], binno=(xrange[1]-xrange[0])/binwidth, \
                                                panelBins=[MassBins[theColNum]], rowno=rowno, colno=colno, aveType='Median', \
                                                linewidth=defaultLinewidth*scaleFactor, textPadder=12., colour=model1col, minInBin=1, \
                                                newSubplot=None, plot1sig=1, plot2sig=1, inputCentres=1, panelLabel=thePanelLabel, panelLabelDecs=2, zorder=-1)
            robs_plot_profile(theX, theY1b[theRowNum], panelProp, xlims[theColNum], ylims[theRowNum], binno=(xrange[1]-xrange[0])/binwidth, \
                                                panelBins=[MassBins[theColNum]], rowno=rowno, colno=colno, aveType='Median', linestyle='dashed', \
                                                linewidth=defaultLinewidth*scaleFactor, textPadder=12., colour=model1col, minInBin=1, \
                                                newSubplot=None, plot1sig=1, plot2sig=1, inputCentres=1, zorder=-2) #colour=model1bcol
            cloudsp, = plt.plot(np.nan, np.nan, linewidth=defaultLinewidth, linestyle='solid', color="black", zorder=-1, label=r'Molecular clouds')
            diffp, = plt.plot(np.nan, np.nan, linewidth=defaultLinewidth, linestyle='dashed', color="black", zorder=-1, label=r'Diffuse gas')             
        else :
            robs_plot_profile(theX, theY[theRowNum], panelProp, xlims[theColNum], ylims[theRowNum], binno=(xrange[1]-xrange[0])/binwidth, \
                                            panelBins=[MassBins[theColNum]], rowno=rowno, colno=colno, aveType='Median', linewidth=defaultLinewidth*scaleFactor, \
                                            textPadder=12., colour=model1col, minInBin=1, panelLabelDecs=2, \
                                            newSubplot=None, plot1sig=1, plot2sig=1, inputCentres=1, panelLabel=thePanelLabel, zorder=-1)
        #legend:
        if ((props[theRowNum] == 'CvsD_DTM') & (theColNum == 0)) :               
            legend0 = plt.legend(handles=[cloudsp, diffp], \
                  loc='lower left', fontsize='small' , frameon=False, labelspacing=0.2)
            legend0._legend_box.align = "left"              
            plt.gca().add_artist(legend0)
        if ((props[theRowNum] == 'Z_stars') | (props[theRowNum] == 'FeH_stars')) & (theColNum == 0) :
            robs_plot_text(panel, SolarNorm, vpos='top', hpos='left')
    
    pdf.savefig(bbox_inches='tight')     
    pdf.close()
    
    
#################
def plot_profs_together(Hubble_h, Omega_M, Omega_Lambda, Samp1, struct1, MassBins=None, props='All') :      
    pdf = PdfPages(PlotDir+"various_profs"+"_"+comp+"_"+Samp1['Model']+"_"+Samp1['Sample_type']+"_z"+Samp1['z_lower']+".pdf")
    El1 = robs_element_list_finder(struct1)
        
    #Set-up plot:
    xlimits = np.array([0.0,3.0])
    ylimits = np.array([2.0,11.99])
    binno = defaultBinno #50 
    theMinInBin = 2
    binwidth = 0.2
    thelinewidth = 2.
    SFRScaleFactor = 1.e+08
    DTGScaleFactor = 1.e+06
    title = None
    xlab = [r'$R/\textnormal{R}_{\textnormal{e}}$']
    rowno = 2
    colno = 1
    
    #For comparison to the HERACLES sample (Abdurro'uf+22a):
    MassBins=np.array([[9.7,11.2]])
    title = 'HERACLES sample'
    kpcBins = None

    if props == 'All' :
        props = [['MstarD','SFRD','McoldD','OISMD'],['HID','H2D','MdustD','DTG']]
        
    #Set-up panels:
    numPanels = rowno*colno
    ylab = [r'log($y$)']
    xlabs = np.full(colno,xlab)
    ylabs = np.full(rowno,ylab)
    xlims = [xlimits,xlimits]
    ylims = [ylimits,ylimits]
    fig = plt.figure(figsize=(4,3.5*rowno))
        
    for kk in range(numPanels) :
        panel = robs_plot_panels(kk, rows=rowno, columns=colno, xlimits=xlims, ylimits=ylims, \
                                     xlab=xlabs, ylab=ylabs, SecondXAxis=None, title=title)
        labs = props[kk]
        for ii in range(len(labs)) :      
            #Calc properties to plot:   
            theX = Samp1['RReRings'][:,:] 
            if labs[ii] == 'MstarD' :  
                ylegend = r'$\Sigma_{*} / \textnormal{M}_{\odot}\textnormal{ kpc}^{-2}$'
                linecol = "orange"  
                avelinecol = robs_get_colour(linecol, shade=0.5)
                MstarDp, = plt.plot([np.nan,np.nan], [np.nan,np.nan], '-', linewidth=thelinewidth, color=avelinecol, label=ylegend)
                PropRings = Samp1['G_samp']['DiskMassRings'] + Samp1['G_samp']['BulgeMassRings']
            elif labs[ii] == 'SFRD' :
                ylegend = r'$10^{8}\,\Sigma_{\rm SFR} / \textnormal{ M}_{\odot}\textnormal{ yr}^{-1}\textnormal{ kpc}^{-2}$'
                linecol = "blue"
                avelinecol = robs_get_colour(linecol, shade=0.5)
                SFRDp, = plt.plot([np.nan,np.nan], [np.nan,np.nan], '-', linewidth=thelinewidth, color=avelinecol, label=ylegend)
                PropRings = Samp1['G_samp']['SfrRings'][:,:] * SFRScaleFactor
            elif labs[ii] == 'HID' :  
                ylegend = r'$\Sigma_{\rm HI} / \textnormal{M}_{\odot}\textnormal{ kpc}^{-2}$'
                linecol = "green"
                avelinecol = robs_get_colour(linecol, shade=0.5)
                HIDp, = plt.plot([np.nan,np.nan], [np.nan,np.nan], '-', linewidth=thelinewidth, color=avelinecol, label=ylegend)
                PropRings = Samp1['G_samp']['ColdGasDiffRings_elements'][:,:,El1["H_NUM"]]
            elif labs[ii] == 'H2D' :
                ylegend = r'$\Sigma_{\rm H2} / \textnormal{M}_{\odot}\textnormal{ kpc}^{-2}$'
                linecol = "black"
                avelinecol = robs_get_colour(linecol, shade=0.5)
                H2Dp, = plt.plot([np.nan,np.nan], [np.nan,np.nan], '-', linewidth=thelinewidth, color=avelinecol, label=ylegend)
                PropRings = Samp1['G_samp']['ColdGasCloudsRings_elements'][:,:,El1["H_NUM"]]
            elif labs[ii] == 'MdustD' :
                ylegend = r'$\Sigma_{\rm dust} / \textnormal{M}_{\odot}\textnormal{ kpc}^{-2}$'
                linecol = "red"
                avelinecol = robs_get_colour(linecol, shade=0.5)
                MdustDp, = plt.plot([np.nan,np.nan], [np.nan,np.nan], '-', linewidth=thelinewidth, color=avelinecol, label=ylegend)
                PropRings = np.nansum(Samp1['G_samp']['DustColdGasCloudsRings_elements'], axis=2) + np.nansum(Samp1['G_samp']['DustColdGasDiffRings_elements'], axis=2)
            elif labs[ii] == 'McoldD' :
                ylegend = r'$\Sigma_{\rm gas} / \textnormal{M}_{\odot}\textnormal{ kpc}^{-2}$'
                linecol = "purple"
                avelinecol = robs_get_colour(linecol, shade=0.9)
                McoldDp, = plt.plot([np.nan,np.nan], [np.nan,np.nan], '-', linewidth=thelinewidth, color=avelinecol, label=ylegend)
                PropRings = Samp1['G_samp']['ColdGasRings'][:,:] - np.nansum(Samp1['G_samp']['DustColdGasCloudsRings_elements'], axis=2) - np.nansum(Samp1['G_samp']['DustColdGasDiffRings_elements'], axis=2)
            elif labs[ii] == 'OISMD' :
                ylegend = r'$\Sigma_{\rm O,ISM} / \textnormal{M}_{\odot}\textnormal{ kpc}^{-2}$'
                linecol = "pink"
                avelinecol = robs_get_colour(linecol, shade=0.5)
                OISMDp, = plt.plot([np.nan,np.nan], [np.nan,np.nan], '-', linewidth=thelinewidth, color=avelinecol, label=ylegend)
                PropRings = Samp1['G_samp']['ColdGasRings_elements'][:,:,El1["O_NUM"]]# * OISMScaleFactor #Oxygen in gas and dust
            elif labs[ii] == 'OStarsD' :
                ylegend = r'$\Sigma_{\rm O,*} / \textnormal{M}_{\odot}\textnormal{ kpc}^{-2}$'
                linecol = "grey"
                avelinecol = robs_get_colour(linecol, shade=0.5)
                OStarsDp, = plt.plot([np.nan,np.nan], [np.nan,np.nan], '-', linewidth=thelinewidth, color=avelinecol, label=ylegend)
                PropRings = Samp1['G_samp']['DiskMassRings_elements'][:,:,El1["O_NUM"]]# * OISMScaleFactor #Oxygen in stellar disc
            elif labs[ii] == 'DTG' :
                ylegend = r'$10^{6}\,M_{\rm dust} /M_{\rm H}$'
                linecol = "cyan"
                avelinecol = robs_get_colour(linecol, shade=0.5)
                DTGp, = plt.plot([np.nan,np.nan], [np.nan,np.nan], '-', linewidth=thelinewidth, color=avelinecol, label=ylegend)
                TotMdustRings = np.nansum(Samp1['G_samp']['DustColdGasCloudsRings_elements'], axis=2) + np.nansum(Samp1['G_samp']['DustColdGasDiffRings_elements'], axis=2) 
                TotHMassRings = Samp1['G_samp']['ColdGasRings_elements'][:,:,El1["H_NUM"]]
                PropRings = (TotMdustRings/TotHMassRings) * DTGScaleFactor * Samp1['RingArea'] #Multiply by Ring areas here to cancel out the dividion by ring area to get the Y below (DTD is not a density property).      
            theY = np.log10(PropRings/Samp1['RingArea']) #in Msun/kpc^2 (or 1.e8*Msun/yr/kpc^2 for the SFRD) (or 1.e6*Msun/yr/kpc^2 for the DTG)
            panelProp = np.log10(Samp1['G_samp']['StellarMass'])
            xrange = [theX[np.isfinite(theX)].min(), theX[np.isfinite(theX)].max()]
                            
            noShade = 0
            if kk == numPanels-1 :
                thePanelLabel = r'log$(M_{*}/\textnormal{M}_{\odot}) =$ '
                basicLabel = None
                basicLabelPos = None
            elif kk == 0 :
                thePanelLabel = None
                basicLabel = None
                basicLabelPos = None
            robs_plot_profile(theX, theY, panelProp, xlimits, ylimits, xlab, ylab, \
                              (xrange[1]-xrange[0])/binwidth, MassBins, aveType='Median', rowno=rowno, colno=colno, \
                               colour=linecol, noShade=noShade, minInBin=1, plot1sig=1, plot2sig=1, linewidth=thelinewidth, inputCentres=1, \
                               newSubplot=None, textPadder=15., panelLabel=thePanelLabel, basicLabel=basicLabel, basicLabelPos=basicLabelPos,
                               discreteBins=kpcBins)
            
        #Legend:
        if kk == 0 :
            theHandles = [MstarDp,SFRDp,McoldDp]
            legend0 = plt.legend(handles=theHandles, loc='upper right', fontsize='small', frameon=False, labelspacing=0.2)
        else :
            theHandles = [HIDp,H2Dp,MdustDp]
            legend0 = plt.legend(handles=theHandles, loc='upper right', fontsize='small', frameon=False, labelspacing=0.2)
        legend0._legend_box.align = "left"
        plt.gca().add_artist(legend0)
    
    pdf.savefig(bbox_inches='tight')     
    pdf.close()
    