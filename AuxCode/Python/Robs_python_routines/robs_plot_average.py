# -*- coding: utf-8 -*-
"""
Created on Thu Nov 11 16:18:02 2021

@author: robyates
"""

"""
robs_plot_average.py
  ;Function to plot the running average (mean/meduan/mode) of a 2D relation.
  ;
  ;Rob Yates 11-11-2021
  ;
  ;UPDATES:
  ;14-02-22: Changed the names of xlimits and ylimits to xlims and ylims, and made them a copy of plotlims (rather than = to), so that the xlimits and ylimits variables in the function from which this is called aren't edited too.
  ;08-11-22: This version was adapted for use at the L-Galaxies workshop 2022
  ;
"""

#Basic packages:
import numpy as np
import matplotlib.pyplot as plt
import copy
from scipy import stats

#Local packages:
# import sys
# sys.path.append('./Robs_python_routines/')
from robs_truncate_colormap import *

def robs_plot_average(theX, theY, binno, aveType='Mean', minInBin=0, colour='blue', linewidth=2., linestyle='solid', \
                      colTrunc=[0.2,0.8], inputCentres=None, plotlims=None, plot1sig=None, plot2sig=None, \
                      invertLinestyle=None, label=None, noShade=1, plotPoints=None, pointsOnly=None, weights=None, \
                      plotMode=None, themarker='o', printxy=None) : 
    if plotlims :
        xlims = np.copy(plotlims[0])
        ylims = np.copy(plotlims[1])
    else :         
        xlims = [np.amin(theX),np.amax(theX)]
        ylims = [np.amin(theY),np.amax(theY)]
        
    if inputCentres :                        
        binwidth = (xlims[1]-xlims[0])/(binno-1)
        xlims[0] = xlims[0]-(0.5*binwidth)
        xlims[1] = xlims[1]+(0.5*binwidth)
    else : 
        binwidth = (xlims[1]-xlims[0])/binno
    
    #Check weights:
    if weights == None :
        weights = np.ones(len(theY))
    
    #0th bin:
    if plotlims : 
        theGalsY = theY[(theX >= xlims[0]) & (theX < xlims[0]+binwidth) & (theY >= ylims[0]) & (theY < ylims[1]) & (weights > 0.0)]
        theWeights = weights[(theX >= xlims[0]) & (theX < xlims[0]+binwidth) & (theY >= ylims[0]) & (theY < ylims[1]) & (weights > 0.0)]
    else : 
        theGalsY = theY[(theX >= xlims[0]) & (theX < xlims[0]+binwidth) & (weights > 0.0)] 
        theWeights = weights[(theX >= xlims[0]) & (theX < xlims[0]+binwidth) & (weights > 0.0)] 
    binCentre = [xlims[0]+(0.5*binwidth)]
    if len(theGalsY) >= minInBin :    
        theMean = [np.average(theGalsY, weights=theWeights)]
        theMedian = [np.median(theGalsY)] #[np.percentile(theGalsY, 50)] #
        theMode = [stats.mode(theGalsY).mode]
        the2p = [np.percentile(theGalsY, 2)]
        the16p = [np.percentile(theGalsY, 16)]
        the84p = [np.percentile(theGalsY, 84)]
        the98p = [np.percentile(theGalsY, 98)]  
    else :
        theMean = [np.nan]
        theMedian = [np.nan]
        theMode = [np.nan]
        the2p = [np.nan]
        the16p = [np.nan]
        the84p = [np.nan]
        the98p = [np.nan]
    
    #Other bins:
    for theBin in range(1,binno) : 
        if plotlims :         
            theGalsY = theY[(theX >= xlims[0]+(theBin*binwidth)) & (theX < xlims[0]+((theBin+1.)*binwidth)) & (theY >= ylims[0]) & (theY < ylims[1]) & (weights > 0.0)]
            theWeights = weights[(theX >= xlims[0]+(theBin*binwidth)) & (theX < xlims[0]+((theBin+1.)*binwidth)) & (theY >= ylims[0]) & (theY < ylims[1]) & (weights > 0.0)]
        else : 
            theGalsY = theY[(theX >= xlims[0]+(theBin*binwidth)) & (theX < xlims[0]+((theBin+1.)*binwidth)) & (weights > 0.0)]
            theWeights = weights[(theX >= xlims[0]+(theBin*binwidth)) & (theX < xlims[0]+((theBin+1.)*binwidth)) & (weights > 0.0)]
        binCentre = np.append(binCentre, (xlims[0]+(theBin*binwidth))+(0.5*binwidth))
        if len(theGalsY) >= minInBin :  
            theMean = np.append(theMean, np.average(theGalsY, weights=theWeights))
            theMedian = np.append(theMedian, np.median(theGalsY))
            theMode = np.append(theMode, stats.mode(theGalsY).mode)
            the2p = np.append(the2p, np.percentile(theGalsY, 2))
            the16p = np.append(the16p, np.percentile(theGalsY, 16))
            the84p = np.append(the84p, np.percentile(theGalsY, 84))
            the98p = np.append(the98p, np.percentile(theGalsY, 98))
        else :
            theMean = np.append(theMean, np.nan)
            theMedian = np.append(theMedian, np.nan)
            theMode = np.append(theMode, np.nan)
            the2p = np.append(the2p, np.nan)
            the16p = np.append(the16p, np.nan)
            the84p = np.append(the84p, np.nan)
            the98p = np.append(the98p, np.nan)
    
    #Select colours:
    if colour == 'blue' :        
        new_cmap = robs_truncate_colormap(plt.get_cmap('Blues'), colTrunc[0], colTrunc[1])             
    elif colour == 'green' :                   
        new_cmap = robs_truncate_colormap(plt.get_cmap('Greens'), colTrunc[0], colTrunc[1])
    elif colour == 'red' :        
        new_cmap = robs_truncate_colormap(plt.get_cmap('Reds'), colTrunc[0], colTrunc[1])  
    elif colour == 'grey' :    
        new_cmap = robs_truncate_colormap(plt.get_cmap('Greys'), colTrunc[0], colTrunc[1])  
    elif colour == 'black' :        
        new_cmap = robs_truncate_colormap(plt.get_cmap('binary'), colTrunc[0], 1.0)  
    elif colour == 'purple' :        
        new_cmap = robs_truncate_colormap(plt.get_cmap('Purples'), colTrunc[0], colTrunc[1])  
    elif (colour == 'yellow') | (colour == 'orange') :        
        new_cmap = robs_truncate_colormap(plt.get_cmap('Wistia'), colTrunc[0], colTrunc[1])
    else : print("****robs_contour_plot: WARNING: Incorrect colour provided. Please use either: 'blue', 'green', 'black', 'red', 'grey', 'purple', 'yellow', 'orange'. (Default: 'blue') *****") 
    
    #Plot percentile ranges:
    if pointsOnly :
        thelinestyle = 'None'        
        if plot2sig : 
            plt.errorbar(binCentre, theMedian, yerr=[theMedian-the2p, the98p-theMedian], \
                        color=new_cmap(0.25), marker=None, linestyle=thelinestyle, linewidth=1.49, capsize=3)
        if plot1sig :        
            plt.errorbar(binCentre, theMedian, yerr=[theMedian-the16p, the84p-theMedian], \
                        color=new_cmap(0.75), marker=None, linestyle=thelinestyle, linewidth=1.49, capsize=3)
    else :
        thelinestyle = copy.copy(linestyle)
        if plot2sig :     
            plt.fill_between(binCentre, the2p, the98p, color=new_cmap(0.25), alpha=0.3)    
        if plot1sig :
            plt.fill_between(binCentre, the16p, the84p, color=new_cmap(0.75), alpha=0.3)
       
    #Plot average:
    if aveType != None :
        if (colour == 'black') : 
            avecol = new_cmap(1.0)
        else : avecol = new_cmap(0.5)
        if (noShade == 1) :
            avecol = colour
        if (plotPoints == 1) | (pointsOnly == 1) :
            themarkeredgecolor = new_cmap(1.0) #'black'
        else : 
            themarker = None
            themarkeredgecolor = None
        if invertLinestyle :    
            if aveType == 'Mean' :    
                plt.plot(binCentre, theMean, color='white', linewidth=linewidth, linestyle=thelinestyle, marker=themarker, markeredgecolor=themarkeredgecolor) #zorder=1
                plt.plot(binCentre, theMean, color=avecol, linewidth=1., linestyle=thelinestyle, marker=themarker, markeredgecolor=themarkeredgecolor) #color=colour
            elif aveType == 'Median' :    
                plt.plot(binCentre, theMedian, color=avecol, linewidth=linewidth, linestyle=thelinestyle, marker=themarker, markeredgecolor=themarkeredgecolor)            
            elif aveType == 'All' :    
                plt.plot(binCentre, theMean, color=avecol, linewidth=linewidth, linestyle=thelinestyle, marker=themarker, markeredgecolor=themarkeredgecolor)       
                plt.plot(binCentre, theMedian, color='white', linewidth=linewidth, linestyle=thelinestyle)
                plt.plot(binCentre, theMedian, color=avecol, linewidth=1., linestyle=thelinestyle, marker=themarker, markeredgecolor=themarkeredgecolor)
                if plotMode == 1 : plt.plot(binCentre, theMode, color=avecol, linewidth=linewidth, linestyle='None', marker="s", markeredgecolor='white')            
            else : print("****robs_contour_plot: WARNING: Type of average not provided. Please choose between 'Median', 'Mean', 'All'. (Default: 'Mean') *****") 
        else :
            if aveType == 'Mean' :    
                plt.plot(binCentre, theMean, color=avecol, linewidth=linewidth, linestyle=thelinestyle, marker=themarker, markeredgecolor=themarkeredgecolor)
            elif aveType == 'Median' :    
                plt.plot(binCentre, theMedian, color=avecol, linewidth=linewidth, linestyle=thelinestyle, marker=themarker, markeredgecolor=themarkeredgecolor)
                plt.plot(binCentre, theMedian, color='white', linewidth=1., linestyle=thelinestyle, marker=themarker, markeredgecolor=themarkeredgecolor)
            elif aveType == 'All' :    
                plt.plot(binCentre, theMean, color=avecol, linewidth=linewidth, linestyle=thelinestyle, marker=themarker, markeredgecolor=themarkeredgecolor)       
                plt.plot(binCentre, theMedian, color=avecol, linewidth=linewidth, linestyle=thelinestyle)
                plt.plot(binCentre, theMedian, color='white', linewidth=1., linestyle=thelinestyle, marker=themarker, markeredgecolor=themarkeredgecolor)
                if plotMode == 1 : plt.plot(binCentre, theMode, color='white', linewidth=linewidth, linestyle='None', marker="s", markeredgecolor=new_cmap(1.0))
            else : print("****robs_contour_plot: WARNING: Type of average not provided. Please choose between 'Median', 'Mean', 'All'. (Default: 'Mean') *****") 
    
    if printxy :
        ww = np.arange(len(binCentre)) #write full arrays  
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
        print("**********\n")
        
