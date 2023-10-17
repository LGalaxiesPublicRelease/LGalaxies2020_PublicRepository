# -*- coding: utf-8 -*-
"""
Created on Tue Jun 14 16:23:25 2022

@author: ry0005
"""

"""
robs_solarnorm_finder.py
  ; Function to define the correct number and ordering of chemical elements in the L-Galaxies outputs, 
  ; for use with the _elements arrays.
  ;
  ;Rob Yates 18-05-2023
  ;
  ;UPDATES: 
  ;
"""
#Basic packages:
import numpy as np

#Local packages:
import sys
sys.path.append('/vol/ph/astro_data2/rmyates/python/python_libraries/Robs_Routines/')
from chemistry import *

# NOTE!: At some point, the following should be automated by e.g. reading the H_NUM, etc definitions from h_metals.h in the L-Galaxies code:
def robs_element_list_finder(struct_type) :
    #if (struct_type == 'plusBinaries32_hotDust_minimal') | (struct_type == 'plusBinaries32_hotDust_extraSFHs') :     
    if 'plusBinaries32' in struct_type :
        Dict = {
            "H_NUM" : 0,
            "He_NUM" : 1,
            "C_NUM" : 2, #NOTE: We are calling this C_NUM here (and therefore in all python loading/plotting scripts), rather than the Cb_NUM used inside L-Galaxies. This allows better matching to user-inputted element names.
            "N_NUM" : 3,
            "O_NUM" : 4,
            "Ne_NUM" : 5,
            "Mg_NUM" : 6,
            "Al_NUM" : 7,
            "Si_NUM" : 8,
            "S_NUM" : 9,
            "Ca_NUM" : 10,
            "Mn_NUM" : 11,
            "Fe_NUM" : 12
        }
    else :
        Dict = {
            "H_NUM" : 0,
            "He_NUM" : 1,
            "C_NUM" : 2, #NOTE: We are calling this C_NUM here (and therefore in all python loading/plotting scripts), rather than the Cb_NUM used inside L-Galaxies. This allows better matching to user-inputted element names.
            "N_NUM" : 3,
            "O_NUM" : 4,
            "Ne_NUM" : 5,
            "Mg_NUM" : 6,
            "Si_NUM" : 7,
            "S_NUM" : 8,
            "Ca_NUM" : 9,
            "Fe_NUM" : 10
        }
    return Dict
    # else : 
    #     print("***** ERROR: robs_element_list_finder.py: No element list obtainable. *****")
    #     return sys.exit(0)        
