#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 23 12:48:28 2021

@author: ry0005
"""

"""
robs_plot_panels.py
  ; Simple script to set-up axes for a NxM plot.
  ; 2nd y-axis not included, so if needed should be added after calling this function.
  ; Based on the more restricted plot2x2.py routine used for Yates+21b.
  ;
  ;EXAMPLE USE:
  xlimits = [0.0,1.0]
  ylimits = [0.0,1.0]
  rows=2
  columns=1
  xlabs = [r'x label']*columns
  ylabs = [r'y label 1', r'y label 2']
  xlims = xlimits #OR: np.array([[0.0,1.0],[2.5,3.4],...])
  ylims = ylimits #OR: np.array([[0.0,1.0],[2.5,3.4],...])
  xticks = [0.0,0.5,1.0]
  yticks = [0.0,0.5,1.0]
  title = r'title'
  fig = plt.figure(figsize=(4,7))
  for ii in range(rows*columns) :
      panel = robs_plot_panels(ii, rows=rows, columns=columns, xlimits=xlims, ylimits=ylims, \
                               xlab=xlabs, ylab=ylabs, xticks=xticks, yticks=yticks, title=title)
      
  ;Rob Yates 04-04-2023
  ;
"""

import numpy as np
import matplotlib.pyplot as plt

def robs_plot_panels(ii, rows, columns, xlimits, ylimits, xlab='x', ylab='y', xticks=None, yticks=None, \
                    SecondXAxis=None, xlimits2=None, xnos2=None, xlab2=None, title=None, trueNumPanels=None) :  
    ##########
    # Set params:
    numPanels = rows*columns
    theRowNum = int(ii/columns) #int((ii/columns)%columns)
    theColNum = ii%columns
    if trueNumPanels == None :
        trueNumPanels = numPanels
    
    #If a single array for various axis attributes is provided for use in all panels:
    if np.array(xlimits[0]).shape == () :
        xlimits = np.array([xlimits]*columns)
    if xticks is not None :
        if np.array(xticks[0]).shape == () :
            xticks = np.array([xticks]*columns)
    if np.array(ylimits[0]).shape == () :
        ylimits = np.array([ylimits]*rows)
    if yticks :
        if np.array(yticks[0]).shape == () :
            yticks = np.array([yticks]*rows)
 
    #Create subplot:
    panel = plt.subplot(rows, columns, ii+1)
    
    #All panels in 1st column:
    if ii % columns == 0 :
        plt.ylabel(ylab[theRowNum])
        if xticks is not None :
            panel.set_xticks(xticks[theColNum])     
        if yticks : #[theRowNum] :
            panel.set_yticks(yticks[theRowNum])    
            panel.set_yticklabels(['{:g}'.format(yt) for yt in yticks[theRowNum]])
        if ii == numPanels-columns :
            plt.xlabel(xlab[theColNum])
            if xticks is not None :
                panel.set_xticklabels(['{:g}'.format(xt) for xt in xticks[theColNum]])
        else :
            plt.setp(panel.get_xticklabels(), visible=False)
            
    #All panels in last row (except the first, which is covered above):
    elif ii >= numPanels-columns :
        plt.xlabel(xlab[theColNum])
        if xticks is not None :
            panel.set_xticks(xticks[theColNum])
            panel.set_xticklabels(['{:g}'.format(xt) for xt in xticks[theColNum]])
        if yticks : #[theRowNum] :
            panel.set_yticks(yticks[theRowNum])  
        if ii != numPanels-columns :
            plt.setp(panel.get_yticklabels(), visible=False)
    elif ii == trueNumPanels-columns :
        plt.xlabel(xlab[theColNum])
        if xticks is not None :
            xticksThisCol = xticks[theColNum][1:]
            panel.set_xticks(xticksThisCol)
            panel.set_xticklabels(['{:g}'.format(xt) for xt in xticksThisCol])
        if yticks : #[theRowNum] :
            panel.set_yticks(yticks[theRowNum])  
        if ii != numPanels-columns :
            plt.setp(panel.get_yticklabels(), visible=False)
            
    #All remaining panels:
    else :
        plt.setp(panel.get_xticklabels(), visible=False)
        plt.setp(panel.get_yticklabels(), visible=False)
        if xticks is not None :
            panel.set_xticks(xticks[theColNum])
        if yticks : #[theRowNum] :
            panel.set_yticks(yticks[theRowNum])    
    
    #Plot title:
    if title :
        if ((theRowNum == 0) & (columns%2 == 0)) | ((theRowNum == 0) & (theColNum == int(np.floor(columns/2)))) :
            plt.title(title)
    
    #In all cases:
    plt.xlim(xlimits[theColNum])
    plt.ylim(ylimits[theRowNum])
    
    # #Add 2nd x axis:
    # if SecondXAxis :
    #     panelb = panel.twiny()
    #     panelb.set_xlim(xlimits[0],xlimits[1])
    #     panelb.set_xticks(xlimits2)
    #     panelb.set_xlabel(xlab2)
    #     panelb.tick_params(direction="in")

    # remove gaps between subplots
    plt.subplots_adjust(hspace=.0)
    plt.subplots_adjust(wspace=.0)
    
    return panel
