# -*- coding: utf-8 -*-
"""
Created on Tue Jan 25 15:55:28 2022

@author: robyates

;UPDATES:
;08-11-22: This version was adapted for use at the L-Galaxies workshop 2022
"""

import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np

def robs_truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap