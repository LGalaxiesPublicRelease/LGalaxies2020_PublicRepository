# -*- coding: utf-8 -*-
"""
Created on Tue Jun 14 16:23:25 2022

@author: ry0005
"""

"""
robs_solarnorm_finder.py
  ;Function to define the correct solar abundances for use when normalising chemical abundances
  ;
  ;Rob Yates 14-06-2022
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

def robs_solarnorm_finder(SolarNorm) :
    if SolarNorm == "A09" :     
        SolarDict = {
            "HeH_sun": HeH_A09,
            "CH_sun": CH_A09,
            "NH_sun": NH_A09,
            "OH_sun": OH_A09,
            "NeH_sun": NeH_A09,
            "MgH_sun": MgH_A09,
            "AlH_sun": AlH_A09,
            "SiH_sun": SiH_A09,
            "SH_sun": SH_A09,
            "CaH_sun": CaH_A09,
            "MnH_sun": MnH_A09,
            "FeH_sun": FeH_A09,
            "HeH_mf_sun": HeH_mf_A09,
            "CH_mf_sun": CH_mf_A09,
            "NH_mf_sun": NH_mf_A09,
            "OH_mf_sun": OH_mf_A09,
            "NeH_mf_sun": NeH_mf_A09,
            "MgH_mf_sun": MgH_mf_A09,
            "AlH_mf_sun": AlH_mf_A09,
            "SiH_mf_sun": SiH_mf_A09,
            "SH_mf_sun": SH_mf_A09,
            "CaH_mf_sun": CaH_mf_A09,
            "MnH_mf_sun": MnH_mf_A09,
            "FeH_mf_sun": FeH_mf_A09,
            "X_mf_sun": X_mf_A09,
            "Y_mf_sun": Y_mf_A09,
            "Z_mf_sun": Z_mf_A09,
            "ZH_mf_sun": ZH_mf_A09
        }
        return SolarDict

    elif SolarNorm == "A09_bulk" :     
        SolarDict = {
            "HeH_sun": HeH_A09,
            "CH_sun": CH_A09,
            "NH_sun": NH_A09,
            "OH_sun": OH_A09,
            "NeH_sun": NeH_A09,
            "MgH_sun": MgH_A09,
            "AlH_sun": AlH_A09,
            "SiH_sun": SiH_A09,
            "SH_sun": SH_A09,
            "CaH_sun": CaH_A09,
            "MnH_sun": MnH_A09,
            "FeH_sun": FeH_A09,
            "HeH_mf_sun": HeH_mf_A09,
            "CH_mf_sun": CH_mf_A09,
            "NH_mf_sun": NH_mf_A09,
            "OH_mf_sun": OH_mf_A09,
            "NeH_mf_sun": NeH_mf_A09,
            "MgH_mf_sun": MgH_mf_A09,
            "AlH_mf_sun": AlH_mf_A09,
            "SiH_mf_sun": SiH_mf_A09,
            "SH_mf_sun": SH_mf_A09,
            "CaH_mf_sun": CaH_mf_A09,
            "MnH_mf_sun": MnH_mf_A09,
            "FeH_mf_sun": FeH_mf_A09,
            "X_mf_sun": X_mf_bulk_A09,
            "Y_mf_sun": Y_mf_bulk_A09,
            "Z_mf_sun": Z_mf_bulk_A09,
            "ZH_mf_sun": Z_mf_bulk_A09/X_mf_bulk_A09
        }
        return SolarDict
    
    elif SolarNorm == "AG89_phot" :
        SolarDict = {
            "HeH_sun": HeH_AG89,
            "CH_sun": CH_AG89,
            "NH_sun": NH_AG89,
            "OH_sun": OH_AG89,
            "NeH_sun": NeH_AG89,
            "MgH_sun": MgH_AG89,
            "AlH_sun": AlH_AG89,
            "SiH_sun": SiH_AG89,
            "SH_sun": SH_phot_AG89,
            "CaH_sun": CaH_phot_AG89,
            "MnH_sun": MnH_phot_AG89,
            "FeH_sun": FeH_phot_AG89,
            "HeH_mf_sun": HeH_mf_AG89,
            "CH_mf_sun": CH_mf_AG89,
            "NH_mf_sun": NH_mf_AG89,
            "OH_mf_sun": OH_mf_AG89,
            "NeH_mf_sun": NeH_mf_AG89,
            "MgH_mf_sun": MgH_mf_AG89,
            "AlH_mf_sun": AlH_mf_AG89,
            "SiH_mf_sun": SiH_mf_AG89,
            "SH_mf_sun": SH_phot_mf_AG89,
            "CaH_mf_sun": CaH_phot_mf_AG89,
            "MnH_mf_sun": MnH_phot_mf_AG89,
            "FeH_mf_sun": FeH_phot_mf_AG89,
            "X_mf_sun": X_mf_AG89,
            "Y_mf_sun": Y_mf_AG89,
            "Z_mf_sun": Z_mete_mf_AG89,
            "ZH_mf_sun": Z_mf_phot_AG89/X_mf_AG89
        }
        return SolarDict

    elif SolarNorm == "AG89_mete" :
        SolarDict = {
            "HeH_sun": HeH_AG89,
            "CH_sun": CH_AG89,
            "NH_sun": NH_AG89,
            "OH_sun": OH_AG89,
            "NeH_sun": NeH_AG89,
            "MgH_sun": MgH_AG89,
            "AlH_sun": AlH_AG89,
            "SiH_sun": SiH_AG89,
            "SH_sun": SH_mete_AG89,
            "CaH_sun": CaH_mete_AG89,
            "MnH_sun": MnH_mete_AG89,
            "FeH_sun": FeH_mete_AG89,
            "HeH_mf_sun": HeH_mf_AG89,
            "CH_mf_sun": CH_mf_AG89,
            "NH_mf_sun": NH_mf_AG89,
            "OH_mf_sun": OH_mf_AG89,
            "NeH_mf_sun": NeH_mf_AG89,
            "MgH_mf_sun": MgH_mf_AG89,
            "AlH_mf_sun": AlH_mf_AG89,
            "SiH_mf_sun": SiH_mf_AG89,
            "SH_mf_sun": SH_mete_mf_AG89,
            "CaH_mf_sun": CaH_mete_mf_AG89,
            "MnH_mf_sun": MnH_mete_mf_AG89,
            "FeH_mf_sun": FeH_mete_mf_AG89,
            "X_mf_sun": X_mf_AG89,
            "Y_mf_sun": Y_mf_AG89,
            "Z_mf_sun": Z_mete_mf_AG89,
            "ZH_mf_sun": Z_mf_mete_AG89/X_mf_AG89
        }
        return SolarDict
    elif SolarNorm == "GAS07" :     
        SolarDict = {
            "HeH_sun": HeH_GAS07,
            "CH_sun": CH_GAS07,
            "NH_sun": NH_GAS07,
            "OH_sun": OH_GAS07,
            "NeH_sun": NeH_GAS07,
            "MgH_sun": MgH_GAS07,
            "AlH_sun": AlH_GAS07,
            "SiH_sun": SiH_GAS07,
            "SH_sun": SH_GAS07,
            "CaH_sun": CaH_GAS07,
            "MnH_sun": MnH_GAS07,
            "FeH_sun": FeH_GAS07,
            "HeH_mf_sun": HeH_mf_GAS07,
            "CH_mf_sun": CH_mf_GAS07,
            "NH_mf_sun": NH_mf_GAS07,
            "OH_mf_sun": OH_mf_GAS07,
            "NeH_mf_sun": NeH_mf_GAS07,
            "MgH_mf_sun": MgH_mf_GAS07,
            "AlH_mf_sun": AlH_mf_GAS07,
            "SiH_mf_sun": SiH_mf_GAS07,
            "SH_mf_sun": SH_mf_GAS07,
            "CaH_mf_sun": CaH_mf_GAS07,
            "MnH_mf_sun": MnH_mf_GAS07,
            "FeH_mf_sun": FeH_mf_GAS07,
            "X_mf_sun": X_mf_GAS07,
            "Y_mf_sun": Y_mf_GAS07,
            "Z_mf_sun": Z_mf_GAS07,
            "ZH_mf_sun": ZH_mf_GAS07
        }
        return SolarDict
    else : 
        print("***** ERROR: No Solar normalisation provided. *****")
        return sys.exit(0)        
