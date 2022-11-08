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
  ;08-11-22: This version was adapted for use at the L-Galaxies workshop 2022
  ;
  ;
"""

#Basic packages:
import numpy as np
import matplotlib.pyplot as plt

#Local packages:
# import sys
# sys.path.append('./Robs_python_routines/')
from robs_truncate_colormap import *

def robs_contour_plot(theX, theY, binno, levels='3sig', colour='blue', linewidth=1.0, linestyle='solid', \
                      weights=None, theRange=None, label=None, fill=None, noOutliers=None, alpha=None, \
                      colTrunc=[0.2,0.8]) : 
    
    #Find contours:    
    bincount, xbins, ybins = np.histogram2d(theX, theY, bins=binno, range=theRange, weights=weights)  
    if levels == 'None' :    
        myLevels = [((100.-68.)/100.)*np.amax(bincount)]
        minVal = 1.1*np.amax(bincount) 
    elif levels == '2sig' :    
        myLevels = [((100.-95.)/100.)*np.amax(bincount), ((100.-68.)/100.)*np.amax(bincount)]
        minVal = ((100.-94.8)/100.)*np.amax(bincount) #intentionally set a little lower, to make sure all the "outlier" points are actually plotted
    elif levels == '3sig' :
        myLevels = [((100.-99.7)/100.)*np.amax(bincount), ((100.-95.)/100.)*np.amax(bincount), ((100.-68.)/100.)*np.amax(bincount)]
        minVal = ((100.-99.5)/100.)*np.amax(bincount) #intentionally set a little lower, to make sure all the "outlier" points are actually plotted
    elif levels == '4sig' :
        myLevels = [((100.-99.9)/100.)*np.amax(bincount), ((100.-99.7)/100.)*np.amax(bincount), ((100.-95.)/100.)*np.amax(bincount), ((100.-68.)/100.)*np.amax(bincount)]
        minVal = ((100.-99.7)/100.)*np.amax(bincount) #intentionally set a little lower, to make sure all the "outlier" points are actually plotted
    else : print("****robs_contour_plot: WARNING: Incorrect levels provided. Please use either: '2sig', '3sig', '4sig'. (Default: '3sig') *****")

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
    if colour == 'blue' :        
        new_cmap = robs_truncate_colormap(plt.get_cmap('Blues'), colTrunc[0], colTrunc[1])             
    elif colour == 'green' :                   
        new_cmap = robs_truncate_colormap(plt.get_cmap('Greens'), colTrunc[0], colTrunc[1])
    elif colour == 'red' :        
        new_cmap = robs_truncate_colormap(plt.get_cmap('Reds'), colTrunc[0], colTrunc[1])  
    elif (colour == 'grey') | (colour == 'black') :        
        new_cmap = robs_truncate_colormap(plt.get_cmap('Greys'), colTrunc[0], colTrunc[1])  
    elif colour == 'purple' :        
        new_cmap = robs_truncate_colormap(plt.get_cmap('Purples'), colTrunc[0], colTrunc[1])  
    elif colour == 'yellow' :        
        new_cmap = robs_truncate_colormap(plt.get_cmap('Wistia'), colTrunc[0], colTrunc[1])
    else : print("****robs_contour_plot: WARNING: Incorrect colour provided. Please use either: 'blue', 'green', 'red', 'grey', 'black', 'purple', 'yellow'. (Default: 'blue') *****") 
                        
    #Plot contours and outliers:        
    if (noOutliers != 1) | (levels == 'None') : cnorm = plt.Normalize(vmin=myLevels[0],vmax=myLevels[-1])      
    if fill :
        if (noOutliers != 1) | (levels == 'None') :
                plt.scatter(theX[outliers == 1], theY[outliers == 1], marker='o', color=new_cmap(cnorm(myLevels[0])), alpha=alpha)                      
        if levels != 'None' :
            plt.contourf(bincount.transpose(), myLevels, extent=[xbins.min(),xbins.max(), ybins.min(),ybins.max()], \
                         cmap=new_cmap, label=label, extend="max", alpha=alpha)  
                            
            plt.contour(bincount.transpose(), myLevels, extent=[xbins.min(),xbins.max(), ybins.min(),ybins.max()], \
                        linewidths=linewidth, colors="black", linestyles=linestyle, label=label)
    else :       
        if (noOutliers != 1) | (levels == 'None') : 
            plt.scatter(theX[outliers == 1], theY[outliers == 1], marker='o', color=new_cmap(cnorm(myLevels[0])), alpha=alpha)        
        if levels != 'None' :
            plt.contour(bincount.transpose(), myLevels, extent=[xbins.min(),xbins.max(), ybins.min(),ybins.max()], \
                        linewidths=linewidth, colors=colour, linestyles=linestyle, label=label)
                