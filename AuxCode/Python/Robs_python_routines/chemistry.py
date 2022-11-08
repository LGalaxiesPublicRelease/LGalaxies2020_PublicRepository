#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 22 13:42:16 2021

@author: robyates
"""

"""
chemistry.py
  ;Provides all the chemical constants required, e.g. atomic masses, solar 
  ;abundances, etc, as global variables.
  ;
  ;Rob Yates 22-04-2021
  ;
  ;08-11-22: This version was adapted for use at the L-Galaxies workshop 2022
  ;
"""

#Atomic weights: 
H_aw = 1.008
He_aw = 4.003
C_aw = 12.01
N_aw = 14.01
O_aw = 16.0
Ne_aw = 20.18
Mg_aw = 24.31
Si_aw = 28.09
S_aw = 32.07
Ca_aw = 40.08
Fe_aw = 55.84

#Solar mass fractions (from Wiersma et al. 2009a):
H_mf = 0.7065
He_mf = 0.2806
C_mf = 2.07E-3
N_mf = 8.36E-4
O_mf = 5.49E-3
Ne_mf = 1.41E-3
Mg_mf = 5.91E-4
Si_mf = 6.83E-4
S_mf = 4.09E-4
Ca_mf = 6.44E-5
Fe_mf = 1.1E-3

#Solar n_i/n_H ratios (from Wiersma et al. 2009a):
H_nr = 1
He_nr = 0.1
C_nr = 2.46E-4
N_nr = 8.51E-5
O_nr = 4.9E-4
Ne_nr = 1.0E-4
Mg_nr = 3.47E-5
Si_nr = 3.47E-5 #Same as Mg
S_nr = 1.86E-5
Ca_nr = 2.29E-6
Fe_nr = 2.82E-5

##########  
#Assumed solar metallicity for many old models:
Z_mf_mod = 0.02


##########
#ASPLUND ET AL. (2009):
#Solar mass ratios:
X_mf_A09 = 0.7381 #This is M_hydrogen/M_baryons in the solar photosphere (Asplund+09, section 3.12).
Y_mf_A09 = 0.2485 #This is M_helium/M_baryons in the solar photosphere (Asplund+09, section 3.12).
Z_mf_A09 = 0.0134 #This is M_metals/M_baryons in the solar photosphere, whereas M_metals/M_hydrogen = 0.0181 from Asplund+09.
ZH_mf_A09 = 0.0181 # = Z_mf_A09 / X_mf_A09. This is M_metals/M_hydrogen in the solar photosphere
X_mf_bulk_A09 = 0.7154 #This is M_hydrogen/M_baryons in the whole bulk of the Sun (Asplund+09, section 3.12).
Y_mf_bulk_A09 = 0.2703 #This is M_helium/M_baryons in the whole bulk of the Sun (Asplund+09, section 3.12).
Z_mf_bulk_A09 = 0.0142 #This is M_metals/M_baryons in the whole bulk of the Sun (Asplund+09, section 3.12).

#Solar abundances (in 12+log(x/H)):
HeH_A09 = 10.93
CH_A09 = 8.43
NH_A09 = 7.83
OH_A09 = 8.69
NeH_A09 = 7.93
MgH_A09 = 7.60
SiH_A09 = 7.51
SH_A09 = 7.12
CaH_A09 = 6.34
FeH_A09 = 7.50

#Solar abundances (mass fractions for x/H):
HeH_mf_A09 = (He_aw/H_aw)*10.0**(HeH_A09-12.)
CH_mf_A09 = (C_aw/H_aw)*10.0**(CH_A09-12.)
NH_mf_A09 = (N_aw/H_aw)*10.0**(NH_A09-12.)
OH_mf_A09 = (O_aw/H_aw)*10.0**(OH_A09-12.)
NeH_mf_A09 = (Ne_aw/H_aw)*10.0**(NeH_A09-12.)
MgH_mf_A09 = (Mg_aw/H_aw)*10.0**(MgH_A09-12.)
SiH_mf_A09 = (Si_aw/H_aw)*10.0**(SiH_A09-12.)
SH_mf_A09 = (S_aw/H_aw)*10.0**(SH_A09-12.)
CaH_mf_A09 = (Ca_aw/H_aw)*10.0**(CaH_A09-12.)
FeH_mf_A09 = (Fe_aw/H_aw)*10.0**(FeH_A09-12.)

#Solar abundances (mass fractions for x/baryons):
H_mf_bary_A09 = X_mf_A09
He_mf_bary_A09 = HeH_mf_A09 * X_mf_A09
C_mf_bary_A09 = CH_mf_A09 * X_mf_A09
N_mf_bary_A09 = NH_mf_A09 * X_mf_A09
O_mf_bary_A09 = OH_mf_A09 * X_mf_A09
Ne_mf_bary_A09 = NeH_mf_A09 * X_mf_A09
Mg_mf_bary_A09 = MgH_mf_A09 * X_mf_A09
Si_mf_bary_A09 = SiH_mf_A09 * X_mf_A09
S_mf_bary_A09 = SH_mf_A09 * X_mf_A09
Ca_mf_bary_A09 = CaH_mf_A09 * X_mf_A09
Fe_mf_bary_A09 = FeH_mf_A09 * X_mf_A09

#Solar abundances (number fractions):
HeH_nf_A09 = 10.0**(HeH_A09-12.)
CH_nf_A09 = 10.0**(CH_A09-12.)
NH_nf_A09 = 10.0**(NH_A09-12.)
OH_nf_A09 = 10.0**(OH_A09-12.)
NeH_nf_A09 = 10.0**(NeH_A09-12.)
MgH_nf_A09 = 10.0**(MgH_A09-12.)
SiH_nf_A09 = 10.0**(SiH_A09-12.)
SH_nf_A09 = 10.0**(SH_A09-12.)
CaH_nf_A09 = 10.0**(CaH_A09-12.)
FeH_nf_A09 = 10.0**(FeH_A09-12.)




##########
#ANDERS & GREVESSE (1989):
#Solar mass ratios:
X_mf_AG89 = 0.70683
Z_mf_phot_AG89 = 0.01941
Z_mf_mete_AG89 = 0.01886
Y_mf_phot_AG89 = 1.0 - Z_mf_phot_AG89 - X_mf_AG89
Y_mf_mete_AG89 = 1.0 - Z_mf_mete_AG89 - X_mf_AG89

#Solar abundances (in 12+log(x/H)):
HeH_AG89 = 10.99
CH_AG89 = 8.56
NH_AG89 = 8.05
OH_AG89 = 8.93
NeH_AG89 = 8.09
MgH_AG89 = 7.58
SiH_AG89 = 7.55
SH_phot_AG89 = 7.21
SH_mete_AG89 = 7.27
CaH_phot_AG89 = 6.36
CaH_mete_AG89 = 6.34
FeH_phot_AG89 = 7.67
FeH_mete_AG89 = 7.51
  
#Solar abundances (mass fractions):
HeH_mf_AG89 = (He_aw/H_aw)*10.0**(HeH_AG89-12.)
CH_mf_AG89 = (C_aw/H_aw)*10.0**(CH_AG89-12.)
NH_mf_AG89 = (N_aw/H_aw)*10.0**(NH_AG89-12.)
OH_mf_AG89 = (O_aw/H_aw)*10.0**(OH_AG89-12.)
NeH_mf_AG89 = (Ne_aw/H_aw)*10.0**(NeH_AG89-12.)
MgH_mf_AG89 = (Mg_aw/H_aw)*10.0**(MgH_AG89-12.)
SiH_mf_AG89 = (Si_aw/H_aw)*10.0**(SiH_AG89-12.)
SH_mf_phot_AG89 = (S_aw/H_aw)*10.0**(SH_phot_AG89-12.)
SH_mf_mete_AG89 = (S_aw/H_aw)*10.0**(SH_mete_AG89-12.)
CaH_mf_phot_AG89 = (Ca_aw/H_aw)*10.0**(CaH_phot_AG89-12.)
CaH_mf_mete_AG89 = (Ca_aw/H_aw)*10.0**(CaH_mete_AG89-12.)
FeH_mf_phot_AG89 = (Fe_aw/H_aw)*10.0**(FeH_phot_AG89-12.)  
FeH_mf_mete_AG89 = (Fe_aw/H_aw)*10.0**(FeH_mete_AG89-12.) 




##########
#GREVESSE, ASPLUND, & SAUVAL (2007):
#Solar mass ratios:
#X_mf_GAS07 = 0.70683
Y_mf_GAS07 = 0.2485 #This is M_helium/M_baryons in the solar photosphere
Z_mf_GAS07 = 0.0122 #This is M_metals/M_baryons in the solar photosphere
X_mf_GAS07 = 1.0 - Y_mf_GAS07 - Z_mf_GAS07 #This is M_hydrogen/M_baryons in the solar photosphere
Y_mf_proto_GAS07 = 0.2735 #This is an estimate of the protosolar M_helium/M_baryons
Z_mf_proto_GAS07 = 0.0132 #This is an estimate of the protosolar M_metals/M_baryons
X_mf_proto_GAS07 = 1.0 - Y_mf_proto_GAS07 - Z_mf_proto_GAS07 #This is an estimate of the protosolar M_hydrogen/M_baryons

#Solar PHOTOSPHERIC abundances (in 12+log(x/H)):
HeH_GAS07 = 10.93
CH_GAS07 = 8.39
NH_GAS07 = 7.78
OH_GAS07 = 8.66
NeH_GAS07 = 7.84 #Indirect estimate
MgH_GAS07 = 7.53
SiH_GAS07 = 7.51
SH_GAS07 = 7.14
CaH_GAS07 = 6.31
FeH_GAS07 = 7.45
  
#Solar abundances (mass fractions):
HeH_mf_GAS07 = (He_aw/H_aw)*10.0**(HeH_GAS07-12.)
CH_mf_GAS07 = (C_aw/H_aw)*10.0**(CH_GAS07-12.)
NH_mf_GAS07 = (N_aw/H_aw)*10.0**(NH_GAS07-12.)
OH_mf_GAS07 = (O_aw/H_aw)*10.0**(OH_GAS07-12.)
NeH_mf_GAS07 = (Ne_aw/H_aw)*10.0**(NeH_GAS07-12.)
MgH_mf_GAS07 = (Mg_aw/H_aw)*10.0**(MgH_GAS07-12.)
SiH_mf_GAS07 = (Si_aw/H_aw)*10.0**(SiH_GAS07-12.)
SH_mf_GAS07 = (S_aw/H_aw)*10.0**(SH_GAS07-12.)
CaH_mf_GAS07 = (Ca_aw/H_aw)*10.0**(CaH_GAS07-12.)
FeH_mf_GAS07 = (Fe_aw/H_aw)*10.0**(FeH_GAS07-12.)  
