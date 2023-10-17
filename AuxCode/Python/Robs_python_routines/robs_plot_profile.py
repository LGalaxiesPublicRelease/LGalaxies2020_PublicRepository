# -*- coding: utf-8 -*-
"""
Created on Mon Feb 14 13:13:37 2022

@author: robyates
"""

"""
robs_plot_profile.py
  ;Function to plot profiles, for data binned into panels, using robs_plot_average().
  ;Requires the figure to be instantiated before the call to robs_plot_profile (see example below).
  ;
  ;Rob Yates 14-02-2021
  ;
  ;ARGUMENTS:
  ;
  ;OPTIONAL ARGUMENTS:
  ;panelPropDim: string: If 'one' [default], then it is assumed the property used to bin items is 1-domensional and a property of the first dimension of the theX and theY arrays (e.g. a whole-galaxy property). If 'same', then the panelProp is assumed to have the same dimensionality as theX and theY (e.g. a ring or SFH property).
  ;samePanel: int: [Default: None]: If on, all bins will be plotted on the same panel, using different colours. If off, each bin is plotted in a separate panel, using the same user-defined colour.
  ;
  ;NON-STANDARD REQUIREMENTS:
  ;robs_plot_average()
  ;robs_truncate_colormap()
  ;
  ;EXAMPLE CALL:
  ;>> plt.figure(figsize=(12,10))
  ;>> robs_plot_profile(theX, theY, xlimits, ylimits, xlabel, ylabel, binno, panelProp, panelBins)
  ;
  ;UPDATES:
  ;24-02-22: Added panelPropDim, which allows the used to define the dimensionality of the property used to assign items into each bin (see above).
  ;24-02-22: Added samePanel, which will plot all bins on the same panel if on. This will prob only work if there are 5 or fewer bins, as I've only set five colours so far (see below).
  ;10-06-22: Replaced samePanel with binnedProp, binnedBins, and binLabel. So now, the property used to separate out the panels is separate to that used to bin data in each panel.
  ;08-11-22: This version was adapted for use at the L-Galaxies workshop 2022
  ;17-10-23: This version was adapted for use with the publicly-available code for the Yates+23 model
  ;
"""

#Basic packages:
import numpy as np
import matplotlib.pyplot as plt

#Local packages:
# import sys
# sys.path.append('./Robs_python_routines/')
from robs_plot_average import *

def robs_plot_profile(theX, theY, panelProp=None, xlims=None, ylims=None, xlabel='x', ylabel='y', binno=10, \
                      panelBins=None, panelPropDim='one', xran=None, yran=None, aveType='Mean', minInBin=1, rowno=None, colno=None, \
                      colour='blue', linewidth=2., linestyle='solid', colTrunc=[0.2,0.8], shade=None, noShade=0, panelLabel=None, \
                      inputCentres=1, plot1sig=1, plot2sig=1, invertLinestyle=None, plotPoints=None, pointsOnly=None, discreteBins=None, \
                      binnedProp=None, binnedBins=None, binLabel=None, weights=None, marker='o', markeredgecolor=None, markersize=None, \
                      basicLabel=None, basicLabelPos='bottom right', printxy=None, returnProps=None, textPadder=20., forceMinX=None, \
                      innerLineScaling=2., panelLabelDecs=1, newSubplot=1, labFontSize='medium', zorder=None, zprop=None) :

    #Set any unset parameters:
    if xlims is None :
        xlims = np.array([theX[np.isfinite(theX)].min(),theX[np.isfinite(theX)].max()]) #[np.nanmin(theX),np.nanmax(theX)]
    if ylims is None :
        ylims = np.array([theY[np.isfinite(theY)].min(),theY[np.isfinite(theY)].max()])
    if (rowno is None) & (colno is None) :
        rowno = 2
        colno = 2
        
    if panelProp is None :
        panelBins = [[0]]
    elif (panelBins is None) :
        panelBins = np.array([[panelProp[np.isfinite(panelProp)].min(),panelProp[np.isfinite(panelProp)].max()]])
    #Define range within which to consider data (not the same as the limits of the plot):
    if xran is None :    
        xran = np.array([theX[np.isfinite(theX)].min(),theX[np.isfinite(theX)].max()])
    if yran is None :    
        yran = np.array([theY[np.isfinite(theY)].min(),theY[np.isfinite(theY)].max()])
        
    #Force minimum radius to be 0.0:
    if forceMinX :
        xlims[xlims < 0.0] = 0.0
        xran[xran < 0.0] = 0.0
    
    #Set-up multi-panel plot:
    plt.rc('font', family='serif', size='14.0')
    binno = int(binno) #binno.astype(int)
    xpadder = (xlims[1]-xlims[0])/textPadder #0.25
    ypadder = (ylims[1]-ylims[0])/textPadder #0.1
    
    if binnedProp is not None :
        colours = ['blue','green','red','grey','purple']
        
    for ii in range(len(panelBins)) : 
        #Set-up panel:
        if newSubplot :
            plt.subplot(rowno, colno, ii+1)
            plt.xlabel(xlabel)
            plt.ylabel(ylabel)
            plt.xlim(xlims)
            plt.ylim(ylims)
            plt.tick_params(direction="in", top=True, bottom=True, left=True, right=True)
            
        if panelProp is None :
            theXBins = theX
            theYBins = theY
            theWeightsBins = weights
            if zprop is not None :
               theZpropBins = zprop
        else :
            wBins = (panelProp >= panelBins[ii][0]) & (panelProp < panelBins[ii][1])
            if panelPropDim == 'same' :
                theXBins = theX[wBins]
                theYBins = theY[wBins]
                if weights is not None :
                    theWeightsBins = weights[wBins]
                if zprop is not None :
                    theZpropBins = zprop[wBins]
            elif panelPropDim == 'one' :          
                theXBins = theX[wBins,:]
                theYBins = theY[wBins,:]
                if weights is not None :
                    theWeightsBins = weights[wBins,:]
                if zprop is not None :
                    theZpropBins = zprop[wBins,:]
            else : print("***** ERROR: robs_plot_profile: panelPropDim unknown. Choose from: 'one' [default], 'same'. *****") 
        if not discreteBins : #for discreteBins, keeps theXBins and theYBins in their original 2D shape (needed for binning along the second dimension in robs_plot_average())
            ww = (np.isfinite(theYBins)) & (~np.isnan(theYBins))  
            theXBins = theXBins[ww]
            theYBins = theYBins[ww]
            if weights is not None :
                theWeightsBins = theWeightsBins[ww]
            if zprop is not None :
                theZpropBins = theZpropBins[ww]
        if weights is None :
            theWeightsBins = None
        if zprop is None :
            theZpropBins = None
        if binnedProp is not None :
            binnedPropBins = binnedProp[ww]            
            for jj in range(len(binnedBins)) : 
                wBins2 = (binnedPropBins >= binnedBins[jj][0]) & (binnedPropBins < binnedBins[jj][1])
                if weights is None :
                    theWeights = None
                else :
                    theWeights = theWeightsBins[wBins2]
                if zprop is None :
                    thezprop = None
                else :
                    thezprop = theZpropBins[wBins2]
                robs_plot_average(theXBins[wBins2], theYBins[wBins2], binno=binno, aveType=aveType, plotlims=[xran,yran], \
                                  minInBin=minInBin, plot1sig=plot1sig, plot2sig=plot2sig, inputCentres=inputCentres, linewidth=linewidth, linestyle=linestyle, \
                                  invertLinestyle=invertLinestyle, shade=shade, noShade=noShade, plotPoints=plotPoints, pointsOnly=pointsOnly, innerLineScaling=innerLineScaling, \
                                  markeredgecolor=markeredgecolor, markersize=markersize, weights=theWeights, colour=colours[jj], themarker=marker, printxy=printxy, \
                                  discreteBins=discreteBins, zorder=zorder, zprop=thezprop)       
                if binLabel is not None : 
                    plt.text(xlims[1]-xpadder, ylims[1]-(1.5*ypadder*(jj+1)), str(binnedBins[jj][0])+r' $<$ '+binLabel+r' $\leq{}$ '+str(binnedBins[jj][1]), \
                             horizontalalignment='right', color=colours[jj], fontsize=labFontSize) 
            #if panelLabel is not None :
            if panelProp is not None :
                if panelLabelDecs == 0 :
                    fullPanelLabel = panelLabel+'['+str(int(np.rint(panelBins[ii][0])))+' - '+str(int(np.rint(panelBins[ii][1])))+']'
                else :
                    fullPanelLabel = panelLabel+'['+str(np.round(panelBins[ii][0],panelLabelDecs))+' - '+str(np.round(panelBins[ii][1],panelLabelDecs))+']'
                plt.text(xlims[0]+xpadder, ylims[0]+ypadder, fullPanelLabel, horizontalalignment='left', color='black', fontsize=labFontSize)
        else :           
            robs_plot_average(theXBins, theYBins, binno=binno, aveType=aveType, plotlims=[xran,yran], colour=colour, \
                              minInBin=minInBin, plot1sig=plot1sig, plot2sig=plot2sig, inputCentres=inputCentres, linewidth=linewidth, linestyle=linestyle, \
                              invertLinestyle=invertLinestyle, shade=shade, noShade=noShade, plotPoints=plotPoints, pointsOnly=pointsOnly, innerLineScaling=innerLineScaling,  \
                              markeredgecolor=markeredgecolor, markersize=markersize, weights=theWeightsBins, themarker=marker, printxy=printxy, \
                              discreteBins=discreteBins, zorder=zorder, zprop=theZpropBins)       
            if (panelProp is not None) & (panelLabel is not None) :
                if panelLabelDecs == 0 :
                    fullPanelLabel = panelLabel+'['+str(int(np.rint(panelBins[ii][0])))+' - '+str(int(np.rint(panelBins[ii][1])))+']'
                else :
                    fullPanelLabel = panelLabel+'['+str(np.round(panelBins[ii][0],panelLabelDecs))+' - '+str(np.round(panelBins[ii][1],panelLabelDecs))+']'
                plt.text(xlims[0]+xpadder, ylims[0]+ypadder, fullPanelLabel, horizontalalignment='left', color='black', fontsize=labFontSize)
        if basicLabel is not None :
            if basicLabelPos == 'bottom right' :            
                plt.text(xlims[1]-xpadder, ylims[0]+ypadder, basicLabel, horizontalalignment='right', color='black', fontsize=labFontSize) 
            elif basicLabelPos == 'bottom left' :
                plt.text(xlims[0]+xpadder, ylims[0]+ypadder, basicLabel, horizontalalignment='left', color='black', fontsize=labFontSize)
            elif basicLabelPos == 'top right' :
                plt.text(xlims[1]-xpadder, ylims[1]-ypadder, basicLabel, horizontalalignment='right', verticalalignment='top', color='black', fontsize=labFontSize)
            elif basicLabelPos == 'top left' :            
                plt.text(xlims[0]+xpadder, ylims[1]-ypadder, basicLabel, horizontalalignment='left', verticalalignment='top', color='black', fontsize=labFontSize)
                
    if returnProps :
        return rowno, colno
            