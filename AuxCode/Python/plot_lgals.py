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

plt.rcParams['text.usetex'] = True
plt.rcParams['font.size'] = 14
plt.rcParams['font.serif'] = 'Times'




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
    binwidth = 0.1
    bin_arr=np.arange(xlimits[0],xlimits[1]+binwidth,binwidth)
    hist=np.histogram(np.log10(Samp1['G_samp']['StellarMass']), bins=bin_arr, range=(xlimits[0],xlimits[1]))     
    plt.plot(hist[1][0:len(hist[1][:])-1]+binwidth/2., np.log10((hist[0][:]/(Volume*binwidth))*Hubble_h**3), label=Samp1['Label'])    

    ##########
    #title:
    plt.text(.5, 0.91, ['SMF', 'z='+char_z_low, str(Samp1['NumGals'])+' gals', Samp1['Cosmology'], Samp1['Simulation'], Samp1['File_type'], Samp1['Model'], Samp1['Version'], Samp1['Sample_type']], transform=fig.transFigure, horizontalalignment='center', color='blue')  
    
    #print figure to pdf:
    pdf.savefig()    
    print("SMF plotted")




#################
def plot_himf(Volume, Hubble_h, Samp1, char_z_low, Samp2=None, pdf=None) :  
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
    hist=np.histogram(np.log10(Samp1['G_samp']['ColdGas']*(1.-Samp1['G_samp']['H2fraction'])), bins=bin_arr, range=(xlimits[0],xlimits[1]))     
    plt.plot(hist[1][0:len(hist[1][:])-1]+binwidth/2., np.log10((hist[0][:]/(Volume*binwidth))*Hubble_h**3), label=Samp1['Label'])    

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
    robs_contour_plot(theX, theY, binno, '2sig', noOutliers=1)
    robs_plot_average(theX, theY, binno=10, aveType='Median', minInBin=50, linewidth=4., inputCentres=1)

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
    # Plot ColdGas metallicity as mass-weighted M_Z/M_tot, using elements arrays, normalised to Sun:
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
    M_cold = np.nansum(Samp1['G_samp']['ColdGas_elements'][:,:], axis=1)
    MZ_cold = np.nansum(Samp1['G_samp']['ColdGas_elements'][:,2:], axis=1)
    Zs = np.log10(MZ_cold / M_cold) - np.log10(Z_mf_bulk_A09)
    theX = np.log10(Samp1['G_samp']['StellarMass'][(M_cold > 0.0) & (MZ_cold > 0.0)])
    theY = Zs[(M_cold > 0.0) & (MZ_cold > 0.0)]
    robs_contour_plot(theX, theY, binno, '2sig', noOutliers=1)  
    robs_plot_average(theX, theY, binno=10, aveType='Median', minInBin=50, linewidth=4., inputCentres=1)

    #title:
    plt.text(.5, 0.92, ['M*-Zg (global, mw, elements arrays)', 'z='+char_z_low, str(Samp1['NumGals'])+' gals', Samp1['Cosmology'], Samp1['Simulation'], Samp1['File_type'], Samp1['Model'], Samp1['Version'], Samp1['Sample_type']], transform=fig.transFigure, horizontalalignment='center', color='blue')  
    
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
    # Plot Stellar metallicity as mass-weighted M_Z/M_tot, using elements arrays, normalised to Sun:
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
    M_stars = np.nansum(Samp1['G_samp']['DiskMass_elements'][:,:], axis=1) +  np.nansum(Samp1['G_samp']['BulgeMass_elements'][:,:], axis=1)
    MZ_stars = np.nansum(Samp1['G_samp']['DiskMass_elements'][:,2:], axis=1) + np.nansum(Samp1['G_samp']['BulgeMass_elements'][:,2:], axis=1)
    Zs = np.log10(MZ_stars / M_stars) - np.log10(Z_mf_bulk_A09)
    theX = np.log10(Samp1['G_samp']['StellarMass'][(M_stars > 0.0) & (MZ_stars > 0.0)])
    theY = Zs[(M_stars > 0.0) & (MZ_stars > 0.0)]
    robs_contour_plot(theX, theY, binno, '2sig', noOutliers=1)  
    robs_plot_average(theX, theY, binno=10, aveType='Median', minInBin=50, linewidth=4., inputCentres=1)

    #title:
    plt.text(.5, 0.92, ['M*-Zs (global, mw, elements arrays)', 'z='+char_z_low, str(Samp1['NumGals'])+' gals', Samp1['Cosmology'], Samp1['Simulation'], Samp1['File_type'], Samp1['Model'], Samp1['Version'], Samp1['Sample_type']], transform=fig.transFigure, horizontalalignment='center', color='blue')  
    
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
    


 