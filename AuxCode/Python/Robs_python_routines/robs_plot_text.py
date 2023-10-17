# -*- coding: utf-8 -*-
"""
Created on Fri Jun  9 13:01:58 2023

@author: ry22aas
"""

"""
robs_plot_text.py
  Simple script to add text to a plot efficiently.
  
  ARGUMENTS:
  ax = axis or subplot onto which the text should be plotted
  text = string or list of strings containing the text to be plotted
  
  EXAMPLE USAGE:
  >> robs_plot_text(['Model 1', 'Model 2', 'Model 3'], vpos='top', hpos='left', \
                    colour=['red', 'green', 'blue'])
  
  Rob Yates 09-06-2023
"""

import numpy as np
import matplotlib.pyplot as plt
import sys

def robs_plot_text(ax, text, xlimits=np.array([0.0,1.0]), ylimits=np.array([0.0,1.0]), vpos='top', hpos='left', \
                   padder=0.05, xpadder=0.05, colour='black', dataCoords=0) :  
    # Set coord system (default: axis):
    if dataCoords == 0 :
        xlimits=np.array([0.0,1.0])
        ylimits=np.array([0.0,1.0])
        transform = ax.transAxes
    else :
        transform = None
    
    #Convert single text string to 1-element string list:
    if isinstance(text, str) :
        text = [text]
        colour = [colour]
    
    #Convert colour from single string to list of strings, iff text has multiple elements
    #(i.e. this is if you want all strings in the text list to be the same colour):
    if (len(text) > 1) & isinstance(colour, str) :
        colour = [colour] * len(text)
    
    #Define padded space between axes and text:
    #xpadder = padder*(xlimits[1]-xlimits[0])
    ypadder = padder*(ylimits[1]-ylimits[0])

    #Determine text x position:
    if hpos == 'left' :
        xpos = xlimits[0]+xpadder
    elif hpos == 'right' :
        xpos = xlimits[1]-xpadder
    elif hpos == 'middle' :
        xpos = xlimits[0]+0.5*(xlimits[1]-xlimits[0]) 
    else :
        #sys.exit("ERROR: robs_plot_text.py: hpos unrecognised. Please choose from: 'left',' right', 'middle'.")
        xpos = hpos
        hpos = 'left'

    #Determine text y position:
    if vpos == 'top' :
        ypos = ylimits[1]-ypadder
        ysign = -1.
    elif vpos == 'bottom' :
        ypos = ylimits[0]+ypadder
        ysign = 1.
        text.reverse()
        colour.reverse()
    elif vpos == 'middle' :
        ypos = ylimits[0]+0.5*(ylimits[1]-ylimits[0])
        ysign = -1.
    else :
        #sys.exit("ERROR: robs_plot_text.py: vpos unrecognised. Please choose from: 'top',' bottom', 'middle'.")
        ypos = vpos
        vpos = 'bottom'
        ysign = 1.

    #Plot text:
    for ii in range(len(text)) :
        plt.text(xpos, ypos+(ysign*ii*ypadder), text[ii], horizontalalignment=hpos, verticalalignment=vpos, \
                 color=colour[ii], transform=transform)
