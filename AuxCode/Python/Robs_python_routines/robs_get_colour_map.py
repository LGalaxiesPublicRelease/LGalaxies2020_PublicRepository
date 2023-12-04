# -*- coding: utf-8 -*-
"""
Created on Mon Apr 10 13:31:20 2023

@author: ry22aas
"""
"""
robs_get_colour_map.py:
    Routine to return a colour map from an inputted colour.
"""

import matplotlib.pyplot as plt

#Local packages:
import sys
sys.path.append('/vol/ph/astro_data2/rmyates/python/python_libraries/Robs_Routines/')
from robs_truncate_colormap import *

def robs_get_colour_map(colour, colTrunc=[0.2,0.8]):
    if colour.lower() == 'blue' :        
        new_cmap = robs_truncate_colormap(plt.get_cmap('Blues'), colTrunc[0], colTrunc[1])             
    elif colour.lower() == 'green' :                   
        new_cmap = robs_truncate_colormap(plt.get_cmap('Greens'), colTrunc[0], colTrunc[1])
    elif colour.lower() == 'lightgreen' :                   
        new_cmap = robs_truncate_colormap(plt.get_cmap('Greens'), colTrunc[0], 0.7)
    elif colour.lower() == 'darkgreen' :                   
        new_cmap = robs_truncate_colormap(plt.get_cmap('Greens'), 0.7, colTrunc[1])
    elif colour.lower() == 'red' :        
        new_cmap = robs_truncate_colormap(plt.get_cmap('Reds'), colTrunc[0], colTrunc[1])  
    elif colour.lower() == 'grey' :    
        new_cmap = robs_truncate_colormap(plt.get_cmap('Greys'), colTrunc[0], colTrunc[1])  
    elif colour.lower() == 'black' :        
        new_cmap = robs_truncate_colormap(plt.get_cmap('binary'), colTrunc[0], 1.0)  
    elif colour.lower() == 'purple' :        
        new_cmap = robs_truncate_colormap(plt.get_cmap('Purples'), colTrunc[0], colTrunc[1])  
    elif (colour.lower() == 'yellow') | (colour.lower() == 'orange') :        
        new_cmap = robs_truncate_colormap(plt.get_cmap('Wistia'), colTrunc[0], colTrunc[1])
    elif (colour.lower() == 'cyan') :        
        new_cmap = robs_truncate_colormap(plt.get_cmap('GnBu'), colTrunc[0], colTrunc[1])
    elif colour.lower() == 'pink' :        
        new_cmap = robs_truncate_colormap(plt.get_cmap('spring'), colTrunc[0], 0.5)  
    elif colour.lower() == 'rainbow' :                   
        new_cmap = robs_truncate_colormap(plt.get_cmap('tab20b'), 0.0, 1.0)
    else : print("****robs_contour_plot: WARNING: Incorrect colour provided. Please use either: 'blue', 'green', 'black', 'red', 'grey', 'purple', 'yellow', 'orange', 'cyan', 'pink'. (Default: 'blue') *****") 

    return new_cmap
