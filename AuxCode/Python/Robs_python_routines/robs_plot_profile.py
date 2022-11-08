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
                      colour='blue', linewidth=2., linestyle='solid', colTrunc=[0.2,0.8], panelLabel=None, \
                      inputCentres=1, plot1sig=1, plot2sig=1, invertLinestyle=None, plotPoints=None, pointsOnly=None, \
                      binnedProp=None, binnedBins=None, binLabel=None, weights=None, marker='o', \
                      basicLabel=None, basicLabelPos='bottom right', printxy=None) :
    #Set any unset parameters:
    if xlims is None :
        xlims = [theX[np.isfinite(theX)].min(),theX[np.isfinite(theX)].max()]
    if ylims is None :
        ylims = [theY[np.isfinite(theY)].min(),theY[np.isfinite(theY)].max()]
    if (rowno is None) & (colno is None) :
        rowno = 2 #int(len(panelBins)/2)
        colno = 2
    if panelProp is None :
        panelBins = [[0]]
    elif (panelBins is None) :
        panelBins = [[panelProp[np.isfinite(panelProp)].min(),panelProp[np.isfinite(panelProp)].max()]]
    
    #Define range within which to consider data (not the same as the limits of the plot):
    if xran is None :    
        xran = [theX[np.isfinite(theX)].min(),theX[np.isfinite(theX)].max()]
    if yran is None :    
        yran = [theY[np.isfinite(theY)].min(),theY[np.isfinite(theY)].max()] 
    
    #Set-up multi-panel plot:
    plt.rc('font', family='serif', size='14.0')
    binno = int(binno) #binno.astype(int)
    if (rowno is None) & (colno is None) :
        rowno = 2 #int(len(panelBins)/2)
        colno = 2
    xpadder = (xlims[1]-xlims[0])/25. #0.25
    ypadder = (ylims[1]-ylims[0])/25. #0.1
    
    if binnedProp is not None :
        colours = ['blue','green','red','grey','purple']
    
    for ii in range(len(panelBins)) : 
        #Set-up panel:
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
        else :
            wBins = (panelProp >= panelBins[ii][0]) & (panelProp < panelBins[ii][1])
            if panelPropDim == 'same' :
                theXBins = theX[wBins]
                theYBins = theY[wBins]
                if weights != None : theWeightsBins = weights[wBins]
            elif panelPropDim == 'one' :          
                theXBins = theX[wBins,:]
                theYBins = theY[wBins,:]
                if weights != None : theWeightsBins = weights[wBins,:]
            else : print("***** ERROR: robs_plot_profile: panelPropDim unknown. Choose from: 'one' [default], 'same'. *****") 
        ww = (np.isfinite(theYBins)) & (~np.isnan(theYBins))  
        theXBins = theXBins[ww]
        theYBins = theYBins[ww]
        if weights != None : theWeightsBins = theWeightsBins[ww]
        else : theWeightsBins = None
        if binnedProp is not None :
            binnedPropBins = binnedProp[ww]            
            for jj in range(len(binnedBins)) : 
                wBins2 = (binnedPropBins >= binnedBins[jj][0]) & (binnedPropBins < binnedBins[jj][1])
                if weights == None :
                    theWeights = None
                else :
                    theWeights = theWeightsBins[wBins2]
                robs_plot_average(theXBins[wBins2], theYBins[wBins2], binno=binno, aveType=aveType, plotlims=[xran,yran], \
                                  minInBin=minInBin, plot1sig=plot1sig, plot2sig=plot2sig, inputCentres=inputCentres, \
                                  invertLinestyle=invertLinestyle, noShade=0, plotPoints=plotPoints, pointsOnly=pointsOnly, \
                                  weights=theWeights, colour=colours[jj], themarker=marker, printxy=printxy)       
                if binLabel is not None : 
                    plt.text(xlims[1]-xpadder, ylims[1]-(1.5*ypadder*(jj+1)), str(binnedBins[jj][0])+r' $<$ '+binLabel+r' $\leq{}$ '+str(binnedBins[jj][1]), \
                             horizontalalignment='right', color=colours[jj]) 
            if panelProp is not None :
                plt.text(xlims[0]+xpadder, ylims[0]+ypadder, panelLabel+'['+str(panelBins[ii][0])+' - '+str(panelBins[ii][1])+']', \
                         horizontalalignment='left', color='black')
        else :           
            robs_plot_average(theXBins, theYBins, binno=binno, aveType=aveType, plotlims=[xran,yran], colour=colour, \
                              minInBin=minInBin, plot1sig=plot1sig, plot2sig=plot2sig, inputCentres=inputCentres, \
                              invertLinestyle=invertLinestyle, noShade=0, plotPoints=plotPoints, pointsOnly=pointsOnly, \
                              weights=theWeightsBins, themarker=marker, printxy=printxy)       
            if (panelProp is not None) & (panelLabel is not None) :
                plt.text(xlims[0]+xpadder, ylims[0]+ypadder, panelLabel+'['+str(panelBins[ii][0])+' - '+str(panelBins[ii][1])+']', \
                         horizontalalignment='left', color='black') 
        if basicLabel is not None :
            if basicLabelPos == 'bottom right' :            
                plt.text(xlims[1]-xpadder, ylims[0]+ypadder, basicLabel, horizontalalignment='right', color='black') 
            elif basicLabelPos == 'bottom left' :
                plt.text(xlims[0]+xpadder, ylims[0]+ypadder, basicLabel, horizontalalignment='left', color='black')
            elif basicLabelPos == 'top right' :
                plt.text(xlims[1]-xpadder, ylims[1]-ypadder, basicLabel, horizontalalignment='right', color='black')
            elif basicLabelPos == 'top left' :            
                plt.text(xlims[0]+xpadder, ylims[1]-ypadder, basicLabel, horizontalalignment='left', color='black')
            
