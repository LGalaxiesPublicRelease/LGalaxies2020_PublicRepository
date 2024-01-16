# -*- coding: utf-8 -*-
"""
Created on Thu Nov 11 16:18:02 2021

@author: ry0005
"""

"""
robs_plot_average.py
  ;Function to plot the running average (mean/median/mode) of a 2D relation.
  ;
  ;Rob Yates 11-11-2021
  ;
  ;OPTIONAL ARGUMENTS:
  ;discreteBins: None or 1: If on, it is assumed that bins are equal to the elements in the theX array (useful e.g. when plotting radial profiles in kpc, where the x-axis is the central radius of the L-Galaxies rings for ALL galaxies)
  ;
  ;UPDATES:
  ;14-02-22: Changed the names of xlimits and ylimits to xlims and ylims, and made them a copy of plotlims (rather than = to), so that the xlimits and ylimits variables in the function from which this is called aren't edited too.
  ;07-09-23: Added the cleaning of the input theX and theY arrays, to remove all NaNs and Infinities.
  ;
"""
#Basic packages:
import numpy as np
import matplotlib.pyplot as plt
import copy
from scipy import stats

#Local packages:
import sys
sys.path.append('/vol/ph/astro_data2/rmyates/python/python_libraries/Robs_Routines/')
from robs_truncate_colormap import *
from robs_get_colour_map import robs_get_colour_map

def robs_plot_average(theX, theY, binno, aveType='Mean', minInBin=0, colour='blue', shade=None, linewidth=2., linestyle='solid', \
                      colTrunc=[0.2,0.8], inputCentres=None, plotlims=None, plot1sig=None, plot2sig=None, levels=None, noRanges=None, \
                      invertLinestyle=None, label=None, noShade=1, plotPoints=None, pointsOnly=None, markeredgecolor=None, \
                      markersize=None, weights=None, plotMode=None, themarker='o', printxy=None, zorder=None, innerLineScaling=2., \
                      discreteBins=None, alpha=0.3, backgroundLineColour='white', zprop=None, returnValues=None) : 
    #Clean inputs:
    clean = np.where((np.isfinite(theX)) & (~np.isnan(theX)) & (np.isfinite(theY)) & (~np.isnan(theY)))
    theX = theX[clean[0]]
    theY = theY[clean[0]]
    if weights is not None :
        weights = weights[clean[0]]
    
    if plotlims is not None :
        xlims = np.copy(plotlims[0])
        ylims = np.copy(plotlims[1])
    else :         
        # # xlims = [np.amin(theX),np.amax(theX)]
        # # ylims = [np.amin(theY),np.amax(theY)]
        # xlims = [np.nanmin(theX),np.nanmax(theX)]
        # ylims = [np.nanmin(theY),np.nanmax(theY)]
        xlims = [np.nanmin(theX[np.isfinite(theX)]),np.nanmax(theX[np.isfinite(theX)])]
        ylims = [np.nanmin(theY[np.isfinite(theY)]),np.nanmax(theY[np.isfinite(theY)])]

    if discreteBins :
        binno = len(theX[0,:])
    else :
        if inputCentres :                        
            binwidth = (xlims[1]-xlims[0])/(binno-1)
            xlims[0] = xlims[0]-(0.5*binwidth)
            xlims[1] = xlims[1]+(0.5*binwidth)
        else : 
            binwidth = (xlims[1]-xlims[0])/binno
    
    #Check weights:
    if weights is None :
        #weights = np.ones(len(theY))
        weights = np.ones(theY.shape)
    
    #0th bin:
    if discreteBins :
        theY1D = theY[:,0] #theY[:,0] #theY[theX == theX[0]]
        weights1D = weights[:,0]
        if zprop is not None :
            zprop1D = zprop[:,0]
        if plotlims :
            theGalsY = theY1D[(theY1D >= ylims[0]) & (theY1D < ylims[1]) & (weights1D > 0.0) & (np.isfinite(weights1D))]
            theWeights = weights1D[(theY1D >= ylims[0]) & (theY1D < ylims[1]) & (weights1D > 0.0) & (np.isfinite(weights1D))]
            if zprop is not None :
                theZprop = zprop1D[(theY1D >= ylims[0]) & (theY1D < ylims[1])]
        else :
            theGalsY = theY1D[(weights1D > 0.0) & (np.isfinite(weights1D))]
            theWeights = weights1D[(weights1D > 0.0) & (np.isfinite(weights1D))]
            if zprop is not None :
                theZprop = zprop1D
        binCentre = theX[0,0]
    else :
        if plotlims : 
            theGalsY = theY[(theX >= xlims[0]) & (theX < xlims[0]+binwidth) & (theY >= ylims[0]) & (theY < ylims[1]) & (weights > 0.0) & (np.isfinite(weights))]
            theWeights = weights[(theX >= xlims[0]) & (theX < xlims[0]+binwidth) & (theY >= ylims[0]) & (theY < ylims[1]) & (weights > 0.0) & (np.isfinite(weights))]
            if zprop is not None :
                theZprop = zprop[(theX >= xlims[0]) & (theX < xlims[0]+binwidth) & (theY >= ylims[0]) & (theY < ylims[1])]
        else : 
            theGalsY = theY[(theX >= xlims[0]) & (theX < xlims[0]+binwidth) & (weights > 0.0) & (np.isfinite(weights))] 
            theWeights = weights[(theX >= xlims[0]) & (theX < xlims[0]+binwidth) & (weights > 0.0) & (np.isfinite(weights))] 
            if zprop is not None :
                theZprop = zprop[(theX >= xlims[0]) & (theX < xlims[0]+binwidth)] 
        binCentre = [xlims[0]+(0.5*binwidth)]
    
    if zprop is None :
        zpropCond = 1
    elif len(theZprop) == 0 :
        zpropCond = 1
    else : #Insert your own custom condition for each radial bin here:
        # (07-0-523): To match the Erroz-Ferrer+19 radial bin selection of >=10% of spaxels are HII-region dominated (rather than DIG dominated):
        if ((len(theZprop[theZprop == 1.]) / len(theZprop)) >= 0.1) :
            zpropCond = 1
        else :
            zpropCond = 0
        yy = 0.0
        
    #if len(theGalsY) >= minInBin :   
    if (len(theGalsY) >= minInBin) & (zpropCond == 1) :   
        #theMean = [np.mean(theGalsY)]
        theMean = [np.average(theGalsY, weights=theWeights)]
        theMedian = [np.median(theGalsY)] #[np.percentile(theGalsY, 50)] #
        theMode = [stats.mode(theGalsY).mode]
        the03p = [np.percentile(theGalsY, 0.3)]
        the2p = [np.percentile(theGalsY, 2)]
        the16p = [np.percentile(theGalsY, 16)]
        the84p = [np.percentile(theGalsY, 84)]
        the98p = [np.percentile(theGalsY, 98)]
        the997p = [np.percentile(theGalsY, 99.7)] 
        theTotal = [np.nansum(theGalsY)]
        theCount = len(theGalsY)
    else :
        theMean = [np.nan]
        theMedian = [np.nan]
        theMode = [np.nan]
        the03p = [np.nan]
        the2p = [np.nan]
        the16p = [np.nan]
        the84p = [np.nan]
        the98p = [np.nan]
        the997p = [np.nan]
        theTotal = [np.nan]
        theCount = [np.nan]
    
    #Other bins:
    for theBin in range(1,binno) : 
        if discreteBins :
            theY1D = theY[:,theBin]
            weights1D = weights[:,theBin]
            if zprop is not None :
                zprop1D = zprop[:,theBin]
            if plotlims :
                theGalsY = theY1D[(theY1D >= ylims[0]) & (theY1D < ylims[1]) & (weights1D > 0.0) & (np.isfinite(weights1D))]
                theWeights = weights1D[(theY1D >= ylims[0]) & (theY1D < ylims[1]) & (weights1D > 0.0) & (np.isfinite(weights1D))]
                if zprop is not None :
                    theZprop = zprop1D[(theY1D >= ylims[0]) & (theY1D < ylims[1])]
            else :
                theGalsY = theY1D[(weights1D > 0.0) & (np.isfinite(weights1D))]
                theWeights = weights1D[(weights1D > 0.0) & (np.isfinite(weights1D))]
                if zprop is not None :
                    theZprop = zprop1D
            binCentre = np.append(binCentre, theX[0,theBin])
        else :
            if plotlims :         
                theGalsY = theY[(theX >= xlims[0]+(theBin*binwidth)) & (theX < xlims[0]+((theBin+1.)*binwidth)) & (theY >= ylims[0]) & (theY < ylims[1]) & \
                                (weights > 0.0) & (np.isfinite(weights))]
                theWeights = weights[(theX >= xlims[0]+(theBin*binwidth)) & (theX < xlims[0]+((theBin+1.)*binwidth)) & (theY >= ylims[0]) & (theY < ylims[1]) \
                                     & (weights > 0.0) & (np.isfinite(weights))]
                if zprop is not None :
                    theZprop = zprop[(theX >= xlims[0]+(theBin*binwidth)) & (theX < xlims[0]+((theBin+1.)*binwidth)) & (theY >= ylims[0]) & (theY < ylims[1])]
            else :
                theGalsY = theY[(theX >= xlims[0]+(theBin*binwidth)) & (theX < xlims[0]+((theBin+1.)*binwidth)) & (weights > 0.0) & (np.isfinite(weights))]
                theWeights = weights[(theX >= xlims[0]+(theBin*binwidth)) & (theX < xlims[0]+((theBin+1.)*binwidth)) & (weights > 0.0) & (np.isfinite(weights))]
                if zprop is not None :
                    theZprop = zprop[(theX >= xlims[0]+(theBin*binwidth)) & (theX < xlims[0]+((theBin+1.)*binwidth))]
            binCentre = np.append(binCentre, (xlims[0]+(theBin*binwidth))+(0.5*binwidth))
        
        if zprop is None :
            zpropCond = 1
        elif len(theZprop) == 0 :
            zpropCond = 1
        else : #Insert your own custom condition for each radial bin here:
            # (07-0-523): To match the Erroz-Ferrer+19 radial bin selection of >=10% of spaxels are HII-region dominated (rather than DIG dominated):
            if ((len(theZprop[theZprop == 1.]) / len(theZprop)) >= 0.1) :
                zpropCond = 1
            else :
                zpropCond = 0
            yy = 0.0
        
        #if len(theGalsY) >= minInBin : 
        if (len(theGalsY) >= minInBin) & (zpropCond == 1) :
            #theMean = np.append(theMean, np.mean(theGalsY))
            theMean = np.append(theMean, np.average(theGalsY, weights=theWeights))
            theMedian = np.append(theMedian, np.median(theGalsY))
            theMode = np.append(theMode, stats.mode(theGalsY).mode)
            the03p = np.append(the03p, np.percentile(theGalsY, 0.3))
            the2p = np.append(the2p, np.percentile(theGalsY, 2))
            the16p = np.append(the16p, np.percentile(theGalsY, 16))
            the84p = np.append(the84p, np.percentile(theGalsY, 84))
            the98p = np.append(the98p, np.percentile(theGalsY, 98))
            the997p = np.append(the997p, np.percentile(theGalsY, 99.7))
            theTotal = np.append(theTotal, np.nansum(theGalsY))
            theCount = np.append(theCount, len(theGalsY))
        else :
            theMean = np.append(theMean, np.nan)
            theMedian = np.append(theMedian, np.nan)
            theMode = np.append(theMode, np.nan)
            the03p = np.append(the03p, np.nan)
            the2p = np.append(the2p, np.nan)
            the16p = np.append(the16p, np.nan)
            the84p = np.append(the84p, np.nan)
            the98p = np.append(the98p, np.nan)
            the997p = np.append(the997p, np.nan)
            theTotal = np.append(theTotal, np.nan)
            theCount = np.append(theCount, np.nan)
    
    # Return binned properties:
    if returnValues :
        Dict = {
            "theX": binCentre,
            "theMean": theMean,
            "theMedian": theMedian,
            "theMode": theMode,
            "the03p": the03p,
            "the2p": the2p,
            "the16p": the16p,
            "the84p": the84p,
            "the98p": the98p,
            "the997p": the997p,
            "theTotal": theTotal,
            "theCount": theCount
        }
    
    #Select colours:
    if colTrunc :
        new_cmap = robs_get_colour_map(colour, colTrunc=[colTrunc[0],colTrunc[1]])
        if colour == 'rainbow' :
            new_cmap_ranges = robs_get_colour_map('grey', colTrunc=[colTrunc[0],colTrunc[1]])
        else :
            new_cmap_ranges = new_cmap
        
    #Set zorder:
    if zorder is None :
        zorder = 0

    #Plot percentile ranges:
    if pointsOnly :
        thelinestyle = 'None'  
        if not noRanges :
            if levels is not None :    
                if levels == '1sig' :    
                    plt.errorbar(binCentre, theMean, yerr=[theMean-the16p, the84p-theMean], \
                                color=new_cmap_ranges(0.25), marker=None, linestyle=thelinestyle, linewidth=1.49, capsize=3, zorder=-1)
                elif levels == '2sig' :    
                    plt.errorbar(binCentre, theMean, yerr=[theMean-the16p, the84p-theMean], \
                                color=new_cmap_ranges(0.25), marker=None, linestyle=thelinestyle, linewidth=1.49, capsize=3, zorder=-1)
                    plt.errorbar(binCentre, theMean, yerr=[theMean-the2p, the98p-theMean], \
                                color=new_cmap_ranges(0.25), marker=None, linestyle=thelinestyle, linewidth=1.49, capsize=3, zorder=-1)
                elif levels == '3sig' :
                    plt.errorbar(binCentre, theMean, yerr=[theMean-the16p, the84p-theMean], \
                                color=new_cmap_ranges(0.25), marker=None, linestyle=thelinestyle, linewidth=1.49, capsize=3, zorder=-1)
                    plt.errorbar(binCentre, theMean, yerr=[theMean-the2p, the98p-theMean], \
                                color=new_cmap_ranges(0.25), marker=None, linestyle=thelinestyle, linewidth=1.49, capsize=3, zorder=-1)
                    plt.errorbar(binCentre, theMean, yerr=[theMean-the03p, the997p-theMean], \
                                color=new_cmap_ranges(0.25), marker=None, linestyle=thelinestyle, linewidth=1.49, capsize=3, zorder=-1)           
            else :
                if plot2sig : 
                    plt.errorbar(binCentre, theMean, yerr=[theMean-the2p, the98p-theMean], \
                                color=new_cmap_ranges(0.25), marker=None, linestyle=thelinestyle, linewidth=1.49, capsize=3, zorder=-1)
                if plot1sig :        
                    plt.errorbar(binCentre, theMean, yerr=[theMean-the16p, the84p-theMean], \
                                color=new_cmap_ranges(0.75), marker=None, linestyle=thelinestyle, linewidth=1.49, capsize=3, zorder=-1)
    else :
        thelinestyle = copy.copy(linestyle)
        if not noRanges :
            if levels is not None :  
                if levels == '1sig' :    
                    plt.fill_between(binCentre, the16p, the84p, color=new_cmap_ranges(0.75), alpha=alpha, zorder=zorder-1)   
                elif levels == '2sig' :    
                    plt.fill_between(binCentre, the2p, the98p, color=new_cmap_ranges(0.25), alpha=alpha, zorder=zorder-1)
                    plt.fill_between(binCentre, the16p, the84p, color=new_cmap_ranges(0.75), alpha=alpha, zorder=zorder-1)  
                elif levels == '3sig' :    
                    plt.fill_between(binCentre, the03p, the997p, color=new_cmap_ranges(0.1), alpha=alpha, zorder=zorder-1)
                    plt.fill_between(binCentre, the2p, the98p, color=new_cmap_ranges(0.25), alpha=alpha, zorder=zorder-1)
                    plt.fill_between(binCentre, the16p, the84p, color=new_cmap_ranges(0.75), alpha=alpha, zorder=zorder-1)
            else :
                if plot2sig :     
                    plt.fill_between(binCentre, the2p, the98p, color=new_cmap_ranges(0.25), alpha=alpha, zorder=zorder-1)    
                if plot1sig :
                    plt.fill_between(binCentre, the16p, the84p, color=new_cmap_ranges(0.75), alpha=alpha, zorder=zorder-1)       
       
    #Plot average:
    if aveType != None :
        if (noShade == 1) :
            avecol = colour
            if (plotPoints == 1) | (pointsOnly == 1) :
                themarkeredgecolor = colour
            else :
                themarker = None
                themarkeredgecolor = None
        else :
            if shade is None :
                if (colour == 'black') : 
                    avecol = new_cmap(1.0)
                else : 
                    avecol = new_cmap(0.5)
            else :
                avecol = new_cmap(shade)
            if (plotPoints == 1) | (pointsOnly == 1) :
                if markeredgecolor == None :
                    themarkeredgecolor = new_cmap(1.0) #'black'
                else :
                    themarkeredgecolor = markeredgecolor
            else : 
                themarker = None
                themarkeredgecolor = None
                
        if invertLinestyle :    
            if aveType == 'Mean' :    
                plt.plot(binCentre, theMean, color=backgroundLineColour, linewidth=linewidth, linestyle=thelinestyle, marker=themarker, markeredgecolor=themarkeredgecolor, markersize=markersize, zorder=zorder) #zorder=1
                plt.plot(binCentre, theMean, color=avecol, linewidth=linewidth/innerLineScaling, linestyle=thelinestyle, marker=themarker, markeredgecolor=themarkeredgecolor, markersize=markersize, zorder=zorder) #color=colour
            elif aveType == 'Total' :    
                plt.plot(binCentre, theTotal, color=backgroundLineColour, linewidth=linewidth, linestyle=thelinestyle, marker=themarker, markeredgecolor=themarkeredgecolor, zorder=zorder) #zorder=1
                plt.plot(binCentre, theTotal, color=avecol, linewidth=linewidth/innerLineScaling, linestyle=thelinestyle, marker=themarker, markeredgecolor=themarkeredgecolor, markersize=markersize, zorder=zorder) #color=colour
            elif aveType == 'Median' :    
                plt.plot(binCentre, theMedian, color=avecol, linewidth=linewidth, linestyle=thelinestyle, marker=themarker, markeredgecolor=themarkeredgecolor, markersize=markersize, zorder=zorder)            
    #            plt.plot(binCentre, theMedian, color=backgroundLineColour, linewidth=linewidth, linestyle=thelinestyle)
    #            plt.plot(binCentre, theMedian, color=avecol, linewidth=1., linestyle=thelinestyle, marker=themarker, markeredgecolor=themarkeredgecolor)
            elif aveType == 'All' :    
                plt.plot(binCentre, theMean, color=avecol, linewidth=linewidth, linestyle=thelinestyle, marker=themarker, markeredgecolor=themarkeredgecolor, markersize=markersize, zorder=zorder)       
                plt.plot(binCentre, theMedian, color=backgroundLineColour, linewidth=linewidth, linestyle=thelinestyle, zorder=zorder)
                plt.plot(binCentre, theMedian, color=avecol, linewidth=linewidth/innerLineScaling, linestyle=thelinestyle, marker=themarker, markeredgecolor=themarkeredgecolor, markersize=markersize, zorder=zorder)
                if plotMode == 1 : plt.plot(binCentre, theMode, color=avecol, linewidth=linewidth, linestyle='None', marker="s", markeredgecolor='white', zorder=zorder)            
            else : print("****robs_contour_plot: WARNING: Type of average not provided. Please choose between 'Median', 'Mean', 'All'. (Default: 'Mean') *****") 
        else :
            if aveType == 'Mean' :    
                plt.plot(binCentre, theMean, color=avecol, linewidth=linewidth, linestyle=thelinestyle, marker=themarker, markeredgecolor=themarkeredgecolor, markersize=markersize, zorder=zorder)
                #plt.plot(binCentre, theMean, color=backgroundLineColour, linewidth=1., linestyle=thelinestyle, marker=themarker, markeredgecolor=themarkeredgecolor)
            elif aveType == 'Total' :    
                plt.plot(binCentre, theTotal, color=avecol, linewidth=linewidth, linestyle=thelinestyle, marker=themarker, markeredgecolor=themarkeredgecolor, markersize=markersize, zorder=zorder)
            elif aveType == 'Median' :    
                # plt.plot(binCentre, theMedian, color=avecol, linewidth=linewidth, linestyle=thelinestyle, marker=themarker, markeredgecolor=themarkeredgecolor, zorder=zorder)
                # plt.plot(binCentre, theMedian, color=backgroundLineColour, linewidth=linewidth/3., linestyle=thelinestyle, marker=themarker, markeredgecolor=themarkeredgecolor, zorder=zorder)
                plt.plot(binCentre, theMedian, color=backgroundLineColour, linewidth=linewidth, linestyle=thelinestyle, marker=themarker, markeredgecolor=themarkeredgecolor, markersize=markersize, zorder=zorder)
                plt.plot(binCentre, theMedian, color=avecol, linewidth=linewidth/innerLineScaling, linestyle=thelinestyle, marker=themarker, markeredgecolor=themarkeredgecolor, markersize=markersize, zorder=zorder)
            elif aveType == 'All' :    
                plt.plot(binCentre, theMean, color=avecol, linewidth=linewidth, linestyle=thelinestyle, marker=themarker, markeredgecolor=themarkeredgecolor, markersize=markersize, zorder=zorder)       
                # plt.plot(binCentre, theMedian, color=avecol, linewidth=linewidth, linestyle=thelinestyle, zorder=zorder)
                # plt.plot(binCentre, theMedian, color=backgroundLineColour, linewidth=linewidth/3., linestyle=thelinestyle, marker=themarker, markeredgecolor=themarkeredgecolor, zorder=zorder)
                plt.plot(binCentre, theMedian, color=backgroundLineColour, linewidth=linewidth, linestyle=thelinestyle, zorder=zorder)
                plt.plot(binCentre, theMedian, color=avecol, linewidth=linewidth/innerLineScaling, linestyle=thelinestyle, marker=themarker, markeredgecolor=themarkeredgecolor, markersize=markersize, zorder=zorder)
                if plotMode == 1 : plt.plot(binCentre, theMode, color=backgroundLineColour, linewidth=linewidth, linestyle='None', marker="s", markeredgecolor=new_cmap(1.0), zorder=zorder)
            else : print("****robs_contour_plot: WARNING: Type of average not provided. Please choose between 'Median', 'Mean', 'All'. (Default: 'Mean') *****") 
    
    if printxy :
        ww = np.arange(len(binCentre)) #write full arrays  
        #ww = np.where((binCentre > 2.0) & (binCentre < 2.2))        
        print("\n**********\nrobs_plot_average():")        
        print("x:")
        print(binCentre[ww])
        print("dx:")
        print(binwidth)
        print("\nMean:")
        print(theMean[ww])
        print("\nMedian:")
        print(theMedian[ww])
        print("\n98th percentile:")
        print(the98p[ww])
        print("\n84th percentile:")
        print(the84p[ww])
        print("\n16th percentile:")
        print(the16p[ww])
        print("\n2nd percentile:")
        print(the2p[ww])
        print("\nTotal:")
        print(theTotal[ww])
        print("\nCount:")
        print(theCount[ww])
        print("**********\n")
        
    if returnValues :
        return Dict
    
        
