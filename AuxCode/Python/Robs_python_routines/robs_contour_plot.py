# -*- coding: utf-8 -*-
"""
Created on Fri Nov  5 11:57:06 2021

@author: robyates
"""

"""
robs_contour_plot.py
  ;Function to make a contour plot from sample X and Y values (and weights, if required).
  ;Outlying points can also be plotted, if required.
  ;Only works for certain pre-defined levels and colours currently (see below).
  ;
  ;Rob Yates 05-11-2021
  ;
  ;ARGUMENTS:
  ;theX: 1D float/double array: Array of X values to bin and plot.
  ;theY: 1D float/double array: Array of Y values to bin and plot.
  ;binno: integer: Number of bins to consider for both the X and Y values.
  ;
  ;OPTIONAL ARGUMENTS:
  ;levels: string: Density levels from which to draw the contours. Currently, only 2-4 sigma allowed. (Default: '2sig')
  ;colour: string: If fill=None: colour of the contour lines. If fill!=None: colour "scheme" used to fill the contour levels. (Default:'blue')
  ;linewidth: float/double: Contour line width. (Default: 1.0)
  ;linestyle: string: Contour line style. (Default: 'solid')
  ;weights: 1D float/double array: Values used to weight the bins from which contours are made. I.e. the 'height' of the contouts at a given X and Y position. (Default: None)
  ;theRange = float/double array of form [2,2]: The lower/upper limits for the X values (1st element) and Y values (2nd element) to consider when binning. (Default: None)
  ;label: string: Label for plot [Untested] (Default: None)
  ;fill: switch: If fill=None: contours will be plotted as lines with no coloured area between them. If fill!=None: contours will be plotted as black lines with the area between them coloured. (Default: None)
  ;noOutliers: switch: If noOutliers=None: Data points beyond the lowest level provided by 'levels' will be plotted as individual points. If noOutliers!=None: these outlying points won't be plotted at all. (Default: None)
  ;outlierFrac = float (0.0 - 1.0): [default = 1.0] The fraction of outlier points to actually plot. These are randomly selected using np.random.choice(). [Useful for decreasing plot file size when the sample is large and there are many outliers.]
  ;colTrunc: 2-element float/double array: Truncates the edges of the chosen colour map using robs_truncate_colormap(), so that a smaller colour range is used for plotting filled contours, for the same range of actual X and Y values considered. (Default: np.array([0.2,0.8]))
  ;alpha = float: For translucent filled contours.
  ;
  ;NON-STANDARD REQUIREMENTS:
  ;robs_truncate_colormap()
  ;
  ;EXAMPLE CALL:
  ;>> robs_contour_plot(theX, theY, binno, '2sig', fill=1, theRange=[xlimits,ylimits])
  ;
  ;UPDATES:
  ;07-02-22: Added the colTrunc optional argument.
  ;07-03-22: Added the alpha optional argument, for translucent filled contours.
  ;22-03-22: Added the levels='None' option, which will just print every point, with no contours.
  ;27-04-23: Added the ability to give custom contour levels as a numpy array.
  ;27-04-23: Added outlierFrac optional argument, allowing user to scale down the number of outliers actually plotted (e.g. if reduced file size is important)
  ;07-09-23: Added the cleaning of the input theX and theY arrays, to remove all NaNs and Infinities.
  ;17-10-23: This version was adapted for use with the publicly-available code for the Yates+23 model
  ;
"""
#Basic packages:
import numpy as np
import matplotlib.pyplot as plt
import sys

#Local packages:
# sys.path.append('/vol/ph/astro_data2/rmyates/python/python_libraries/Robs_Routines/')
from robs_truncate_colormap import *
from robs_get_colour_map import robs_get_colour_map
from robs_get_colour import robs_get_colour

def robs_contour_plot(theX, theY, binno, levels='3sig', colour='blue', linewidth=1.0, linestyle='solid', linecolour='black', \
                      weights=None, theRange=None, label=None, fill=None, noOutliers=None, alpha=None, colTrunc=[0.2,0.8], \
                      outlierFrac=1.0, zorder=None) : 
    #Clean inputs:
    clean = np.where((np.isfinite(theX)) & (~np.isnan(theX)) & (np.isfinite(theY)) & (~np.isnan(theY)))
    theX = theX[clean[0]]
    theY = theY[clean[0]]
    if weights is not None :
        weights = weights[clean[0]]
    
    #Set sigmas:
    Sig1 = 68.
    Sig2 = 95.
    Sig3 = 99.7
    Sig4 = 99.99
    Sig5 = 99.9999
    
    #Find contours:    
    bincount, xbins, ybins = np.histogram2d(theX, theY, bins=binno, range=theRange, weights=weights)  
    if levels == 'None' :    
        myLevels = [((100.-Sig1)/100.)*np.amax(bincount)]
        #minVal = 1.1*np.amax(bincount) #((100.-67.8)/100.)*np.amax(bincount) #intentionally set a little lower, to make sure all the "outlier" points are actually plotted
    elif levels == '2sig' :    
        myLevels = [((100.-Sig2)/100.)*np.amax(bincount), ((100.-Sig1)/100.)*np.amax(bincount)]
        #minVal = ((100.-94.8)/100.)*np.amax(bincount) #intentionally set a little lower, to make sure all the "outlier" points are actually plotted
    elif levels == '3sig' :
        myLevels = [((100.-Sig3)/100.)*np.amax(bincount), ((100.-Sig2)/100.)*np.amax(bincount), ((100.-Sig1)/100.)*np.amax(bincount)]
        #minVal = ((100.-99.5)/100.)*np.amax(bincount) #intentionally set a little lower, to make sure all the "outlier" points are actually plotted
    elif levels == '4sig' :
        myLevels = [((100.-Sig4)/100.)*np.amax(bincount), ((100.-Sig3)/100.)*np.amax(bincount), ((100.-Sig2)/100.)*np.amax(bincount), ((100.-Sig1)/100.)*np.amax(bincount)]
        #minVal = ((100.-99.95)/100.)*np.amax(bincount) #intentionally set a little lower, to make sure all the "outlier" points are actually plotted
    elif levels == '5sig' :
        myLevels = [((100.-Sig5)/100.)*np.amax(bincount), ((100.-Sig4)/100.)*np.amax(bincount), ((100.-Sig3)/100.)*np.amax(bincount), ((100.-Sig2)/100.)*np.amax(bincount), ((100.-Sig1)/100.)*np.amax(bincount)]
        #minVal = ((100.-99.999)/100.)*np.amax(bincount) #intentionally set a little lower, to make sure all the "outlier" points are actually plotted
    else : #For when exact percentages are inputted
        levels = sorted(levels, reverse=True) #Make sure the levels list is ordered from outermost level to innermost level
        myLevels = [((100.-levels[0])/100.)*np.amax(bincount)]
        for ii in range(1,len(levels)) :
            myLevels = np.append(myLevels,((100.-levels[ii])/100.)*np.amax(bincount))
        #minVal = ((100.-(1.0*levels[0]))/100.)*np.amax(bincount) #intentionally set a little lower, to make sure all the "outlier" points are actually plotted
    minVal = myLevels[0]
    
    #Find outliers:
    if (noOutliers != 1) | (levels == 'None') :    
        outliers = np.full(len(theX), 1, dtype=int)
        for xx in range(binno) :
            for yy in range(binno) :
                points = (theX >= xbins[xx]) & (theX < xbins[xx+1]) & (theY >= ybins[yy]) & (theY < ybins[yy+1])  
                if np.any(weights) : 
                    val = np.sum(weights[points])
                else : 
                    val = np.count_nonzero(points)
                if val >= minVal : outliers[points] = 0
                
    #Select colours:
    new_cmap = robs_get_colour_map(colour, colTrunc=[colTrunc[0],colTrunc[1]])
 
    #Set linestyle and colour:
    if linestyle == 'graded' :
        if len(myLevels) == 1 :
            linestyle = 'solid'
        elif len(myLevels) == 2 :
            linestyle = ['dashed','solid']
        elif len(myLevels) == 3 :
            linestyle = ['dotted','dashed','solid']
        elif len(myLevels) == 4 :
            linestyle = [(0,(1,4)),'dotted','dashed','solid']
        elif len(myLevels) == 5 :
            linestyle = [(0,(1,7)),(0,(1,4)),'dotted','dashed','solid']    
    if linecolour is not 'black' :
        linecolour = robs_get_colour(colour, shade=1.0)
               
    #Plot contours and outliers:        
    if (noOutliers != 1) | (levels == 'None') :
        cnorm = plt.Normalize(vmin=myLevels[0],vmax=myLevels[-1])  
        if ((outlierFrac < 0.0) | (outlierFrac > 1.0)) :
            print("****robs_contour_plot: ERROR: outlierFrac = "+str(outlierFrac)+" is beyond acceptable range of 0.0 - 1.0.*****\nExiting...") 
            sys.exit()
        outlierIndices = np.where(outliers == 1)
        outliersToPlot = np.random.choice(outlierIndices[0], int(outlierFrac*len(outlierIndices[0])), replace=False)
    if fill :
        if (noOutliers != 1) | (levels == 'None') :
            plt.scatter(theX[outliersToPlot], theY[outliersToPlot], marker='o', color=new_cmap(cnorm(myLevels[0])), alpha=alpha, zorder=zorder)
        if levels != 'None' :
            plt.contourf(bincount.transpose(), myLevels, extent=[xbins.min(),xbins.max(), ybins.min(),ybins.max()], \
                         cmap=new_cmap, label=label, extend="max", alpha=alpha, zorder=zorder)  
                            
            plt.contour(bincount.transpose(), myLevels, extent=[xbins.min(),xbins.max(), ybins.min(),ybins.max()], \
                        linewidths=linewidth, colors=linecolour, linestyles=linestyle, label=label)#, zorder=zorder) #colors='black''
    else :       
        if (noOutliers != 1) | (levels == 'None') : 
            plt.scatter(theX[outliersToPlot], theY[outliersToPlot], marker='o', color=new_cmap(cnorm(myLevels[0])), alpha=alpha, zorder=zorder)
        if levels != 'None' :
            plt.contour(bincount.transpose(), myLevels, extent=[xbins.min(),xbins.max(), ybins.min(),ybins.max()], \
                        linewidths=linewidth, colors=colour, linestyles=linestyle, label=label)#, zorder=zorder) #colors=new_cmap(1.0)
                