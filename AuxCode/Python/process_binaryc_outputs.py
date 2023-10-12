# -*- coding: utf-8 -*-
"""
Created on Thu Oct 12 09:28:57 2023

@author: R. Yates & D. Hendricks
"""

"""
Simple script to read-in binary_c ensemble outputs, process them into look-up
tables for L-Galaxies, and plot diagnostic plots.
"""

import json
import pandas as pd
import numpy as np
from collections import OrderedDict, Counter
from ensemble_functions.ensemble_functions import inflate_ensemble_with_lists, find_columnames_recursively
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from datetime import datetime
now = datetime.now()
import sqlite3
import filecmp

#Local packages:
import sys
sys.path.append(r'C:\Users\ry22aas\robyates\Astro\python_libraries\Robs_Routines\\')
from robs_plot_panels import robs_plot_panels
from robs_plot_text import robs_plot_text

##########
# Switches:
##########
LOAD_DATA = 1 #Turn on if you want to load and/or plot any data.
SINGLE_METALLICITY = 0 #If on, only one metallicity (determiend by singleZ below) is considered when loading/writing/plotting.
CONVERT_DATA_TO_LINEAR_TIME_BINS = 0 #Turn off when generating yields for input into L-Galaxies (the log-to-linear conversion is done there).
CONVERT_DATA_TO_PER_SOLARMASS = 0 #TURN OFF ALWAYS. IT'S NOT NEEDED! Apparently, this is already done before the raw binary_c yields are tabulated.
CORRECT_SN_RATE_UNITS = 1 #If on, SN "rates" from the binary_c ensemble files will be corrected from a SN number to a SN rate per logt, using the log timestep width from the file (i.e. bc_logdt, calculated below).
ORIGINAL_ENSEMBLES = 1 #If on, this code will run with the ensembles based on the original parameter settings (e.g. Hurley+02 fallback fractions, pre-Temmink SNe-Ia, Kroupa IMF, etc)
SINGLE_STARS = 0 #If on, the single star stellar populations from binary_c will be loaded, rather than the binary star ones. (30-05-23)
RESCALED_SINGLE_STAR_POP = 0 #If on, the "single star" ensembles providing only the single star component of a mixed ensemble are loaded (good to assess impact of binaries in a given population). If off, the 100% single-star ensembles are loaded (good for comp to e.g. oldYields). ****Only valid if SINGLE_STARS is onD****
READ_SIM_INFO = 0 #Only works if ORIGINAL_ENSEMBLES is off. If on, the sim_info_file will be read and some key meta data printed. (N.B> This info file is not available for the older binary_c ensembles).
IMF = 'kroupa2001' #Only works if ORIGINAL_ENSEMBLES is off. Choose from: 'kroupa2001', 'chabrier2003' #Only relevant for single-star SPs, as binary-star SPs only come for one IMF option (old SPs: 'kroupa2001', new SPs: 'chabrier2003')
SN_IA_CHAND_COAL_FACTOR = 1. #10000. #This is the fudge factor used to arbitrarily increase the number of (and yield from) DD-scenario SNe-Ia (i.e. SN_IA_CHAND_coal) in order to match that predicted by Claeys+14 using an older version of binary_c. (07-07-23)
WRITE_BINARYC_YIELDS = 0
WRITE_SN_RATES = 0
WRITE_HEADER_FILE = 0 #If on, this script will generate a new h_metals.h file for L-Galaxies, which includes (a) the correct number of chemical elements to consider, and (b) their order.
AL26 = 0 #If on, calcium will be replaced by Al26 when extracting yields/metals/amsses from the binary_c ensemble files.
NEW_LINE_PER_ELEMENT = 0 #Turn off when generating yields for input into L-Galaxies (they may not be read-in correctly otherwise).
ROBS_PLOT_LAYOUT = 1 #If on, my style of plots (times new roman font, etc) is used. If off, David's is used.
PLOT_BINARYC_YIELDS = 1
PLOT_ONE_ELEMENT = 1 #If on, the normalised mass ejection rate for just one chemical element (given by ele_to_plot) will be plotted from the binary_c yields. If off, the sum of all the metal elements is plotted. (25-08-23)
PLOT_METAL_ELEMENTS = 0 #If on (and PLOT_ONE_ELEMENTS is off), only the elements heavier than He are summed when plotting the metal ejection rate from binary_c.
CONVERT_PLOTS_TO_LINEAR_TIME_BINS = 1 #When on, this will convert the binary_c yields from dM/dlogt * 1/Msun to dM/dt * 1/Msun for plots only.
CONVERT_PLOTS_TO_PER_SOLARMASS = 0 #TURN OFF ALWAYS. IT'S NOT NEEDED! Apparently, this is already done before the raw binary_c yields are tabulated.
PLOT_LGALAXIES_YIELDS = 0
PLOT_YIELD_COMPARISONS = 0
PLOT_SN_RATES = 1 #If on, the SN rates will be plotted in a panel underneath the total mass and metal mass ejection rates.
PLOT_SINGLE_AND_BINARY_STARS = 1 #If on, both the singleStar and binaryStar ensembles from binary_c will be plotted together in the PLOT_YIELD_COMPARISONS plots.
PLOT_ONE_COLUMN = 1 #If on, the individual element evolutions will be plotted all in one colum (for Supp. Material). If off, they will be split into 2 columns (for main paper).
PLOT_PIECHARTS = 1 #If on, pie charts showing the relatvie masses of each element ejected by the binary_c and Yates+13 yields is plotted next to the burst yield evolutions.
MAXSNIIMASS25_YIELDS = 0 # If on, data from L-Galaxies when assuming the max mass for SN-II progenitors is 25Msun (rather than the defulat 120Msun) is used.
READ_SQL_FILE = 0
PLOT_UNPROCESSED_COMP = 0 #If on, plots of the assumed unprocessed component of an SPs elected metals are made. (This is the metal mass subtracted from the stars in L-GALAXIES).

switches_dict = {
    'LOAD_DATA' : LOAD_DATA,
    'SINGLE_METALLICITY' : SINGLE_METALLICITY,
    'CONVERT_DATA_TO_LINEAR_TIME_BINS'  : CONVERT_DATA_TO_LINEAR_TIME_BINS,
    'CONVERT_DATA_TO_PER_SOLARMASS' : CONVERT_DATA_TO_PER_SOLARMASS,
    'CORRECT_SN_RATE_UNITS' : CORRECT_SN_RATE_UNITS,
    'ORIGINAL_ENSEMBLES' : ORIGINAL_ENSEMBLES,
    'SINGLE_STARS' : SINGLE_STARS,
    'RESCALED_SINGLE_STAR_POP' : RESCALED_SINGLE_STAR_POP,
    'READ_SIM_INFO' : READ_SIM_INFO,
    'IMF' : IMF,
    'SN_IA_CHAND_COAL_FACTOR' : SN_IA_CHAND_COAL_FACTOR,
    'WRITE_BINARYC_YIELDS' : WRITE_BINARYC_YIELDS,
    'WRITE_SN_RATES' : WRITE_SN_RATES,
    'WRITE_HEADER_FILE' : WRITE_HEADER_FILE,
    'AL26' : AL26,
    'NEW_LINE_PER_ELEMENT'  : NEW_LINE_PER_ELEMENT,
    'ROBS_PLOT_LAYOUT'  : ROBS_PLOT_LAYOUT,
    'PLOT_BINARYC_YIELDS' : PLOT_BINARYC_YIELDS,
    'PLOT_ONE_ELEMENT' : PLOT_ONE_ELEMENT,
    'PLOT_METAL_ELEMENTS' : PLOT_METAL_ELEMENTS,
    'CONVERT_PLOTS_TO_LINEAR_TIME_BINS'  : CONVERT_PLOTS_TO_LINEAR_TIME_BINS,
    'CONVERT_PLOTS_TO_PER_SOLARMASS' : CONVERT_PLOTS_TO_PER_SOLARMASS,
    'PLOT_LGALAXIES_YIELDS' : PLOT_LGALAXIES_YIELDS,
    'PLOT_YIELD_COMPARISONS' : PLOT_YIELD_COMPARISONS,
    'PLOT_SN_RATES': PLOT_SN_RATES,
    'PLOT_SINGLE_AND_BINARY_STARS' : PLOT_SINGLE_AND_BINARY_STARS,
    'PLOT_ONE_COLUMN' : PLOT_ONE_COLUMN,
    'PLOT_PIECHARTS' : PLOT_PIECHARTS,
    'MAXSNIIMASS25_YIELDS' : MAXSNIIMASS25_YIELDS,
    'READ_SQL_FILE' : READ_SQL_FILE,
    'PLOT_UNPROCESSED_COMP' : PLOT_UNPROCESSED_COMP
}

##########
# Files & Directories:
##########
#basedir = r"C:\Users\ry0005\robyates\Astro\L-Galaxies\LG2020_plusBinaries\misc\David_loading_data_into_dataframe_and_plotting"
basedir = r"C:\Users\ry22aas\robyates\Astro\L-Galaxies\L-Galaxies2020_plusBinaries\misc\David_loading_data_into_dataframe_and_plotting\\"
##########
# ORIGINAL ENSEMBLES:
##########
if ORIGINAL_ENSEMBLES == 1 :
    ensembledir = basedir+"ensemble_suite/original_ensembles\\"
    ensemble_file_metallicities = np.array([0.0001,0.001,0.004,0.008,0.01,0.03])
    if SINGLE_STARS == 1 :
        ensembledir = ensembledir+'/singleStars\\'
        if RESCALED_SINGLE_STAR_POP == 1 : #Load "single star" ensembles providing only the single star component of a mixed ensemble
            ensembledir = ensembledir+'/rescaled\\'
            if IMF == 'kroupa2001' :
                IMFdir = "kroupa2001\\"
                ensemble_file_codes = ['dc08d370','d5383a16','b55c6746','05ebbd6d','373ff2e4','ae1cf746'] #['05ebbd6d','373ff2e4','ae1cf746','b55c6746','d5383a16','dc08d370']
                bcYields_tag = 'binaryc_origEnsembles_singleStars_rescaled_kroupa2001'
                ensembledir2 = basedir+'ensemble_suite/original_ensembles\\'+'/binaryStars\\'
                ensemble_file_codes2 = ['6d984eac','dab8bea4','3083cc1f','47cd834d','d723bdf4','a28ea3be']
                if PLOT_SINGLE_AND_BINARY_STARS :
                    bcYields_tag2 = 'binaryc_origEnsembles_binaryStars_kroupa2001'
            else :
                errString = "***** ERROR: IMF not recognised. For ORIGINAL_ENSEMBLES, please choose from: 'kroupa2001'. *****"
                sys.exit(errString) 
        else : #100% single-star ensembles
            ensembledir = ensembledir+'/100pc_singleStars\\'
            if IMF == 'kroupa2001' :
                IMFdir = "kroupa2001\\"
                ensemble_file_codes = ['319d1c05','bcc22d43','da1bfc0c','4fb2a217','b3e037d4','1e923d7c'] #['1e923d7c','4fb2a217','319d1c05','b3e037d4','bcc22d43','da1bfc0c']
                bcYields_tag = 'binaryc_origEnsembles_singleStars_100pcSingleStars_kroupa2001'
                if PLOT_SINGLE_AND_BINARY_STARS :
                    bcYields_tag2 = 'binaryc_origEnsembles_binaryStars_kroupa2001'
            else :
                errString = "***** ERROR: IMF not recognised. For ORIGINAL_ENSEMBLES, please choose from: 'kroupa2001'. *****"
                sys.exit(errString)         
    else :
        ensembledir = ensembledir+'/binaryStars\\'
        if IMF == 'kroupa2001' :
            IMFdir = "kroupa2001\\"
            ensemble_file_codes = ['6d984eac','dab8bea4','3083cc1f','47cd834d','d723bdf4','a28ea3be'] #N.B. Run READ_SQL_FILE to get these.
            bcYields_tag = 'binaryc_origEnsembles_binaryStars_kroupa2001'
            if PLOT_SINGLE_AND_BINARY_STARS :
                bcYields_tag2 = 'binaryc_origEnsembles_singleStars_kroupa2001'
        else :
            errString = "***** ERROR: IMF not recognised. For ORIGINAL_ENSEMBLES, please choose from: 'kroupa2001'. *****"
            sys.exit(errString) 
##########
# NEW ENSEMBLES:
##########
else :
    ensembledir = basedir+"ensemble_suite/new_ensembles\\"
    ensemble_file_metallicities = np.array([0.0001,0.001,0.004,0.008,0.01,0.02,0.03])
    if SINGLE_STARS == 1 :
        ensemble_file_metallicities = np.array([0.0001,0.001,0.004,0.008,0.01,0.02,0.03])
        simulation_info_file = r'simulation_20230529_115118_settings'
        ensembledir = ensembledir+'/singleStars\\'
        if IMF == 'kroupa2001' :
            IMFdir = "kroupa2001\\"
            ensemble_file_codes = ['463427c5','f762a7ef','651b4026','76b48666','edfb6989','e38affd0','56a1621a']
            bcYields_tag = 'binaryc_newEnsembles_singleStars_kroupa2001'
        elif IMF == 'chabrier2003' :
            IMFdir = "chabrier2003\\"
            ensemble_file_codes = ['ebfb2eeb','88e4a468','bb1e98e1','af89a73a','a6b81372','18b80cc1','7d46969a']
            bcYields_tag = 'binaryc_newEnsembles_singleStars_chabrier2003'
        else :
            errString = "***** ERROR: IMF not recognised. For NEW_ENSEMBLES with single stars, please choose from: 'kroupa2001', 'chabrier2003'. *****"
            sys.exit(errString) 
    else : #BINARY_SPs:
        ensemble_file_metallicities = np.array([0.0001])
        simulation_info_file = r'simulation_20230601_074244_settings'
        ensembledir = ensembledir+'/binaryStars\\'
        if IMF == 'chabrier2003' :
            IMFdir = "chabrier2003\\"
            ensemble_file_codes = ['cfd7bbf6']
            bcYields_tag = 'binaryc_newEnsembles_binaryStars_chabrier2003'
        else :
            errString = "***** ERROR: IMF not recognised. For NEW_ENSEMBLES with binary stars, please choose from: 'chabrier2003'. *****"
            sys.exit(errString) 

ensemble_datafile = 'ensemble_output' #'ensemble_output-100' #basedir+'/other_ensembles/'+'ensemble_output-f7b882ec' #
#Single metallicity options:
singleZ =0.004 # 0.0001 #
#Written-out files:
bc_yieldsdir = basedir+"processed_binaryc_yields\\"
bc_yieldsfile = ensemble_datafile
LGalaxies_timestepsfile = 'timestep_times'
LGalaxies_dtsfile = 'dt_Myr'
#Burst yields:
burstYields_folder = r"C:\Users\ry22aas\robyates\Astro\L-Galaxies\L-Galaxies2020_plusBinaries\L-Galaxies2020_plusBinaries_version\YieldTables\BurstYields\\" #basedir+"BurstYields\\" #"/Burst_yields_LGals/"
LGalaxies_metallicities = ['0.0004','0.0040','0.0080','0.0200','0.0500','1.0000']
LGalaxies_burstyieldsfile = ['AGB_M01_burst_yields','SNII_P98_burst_yields','SNIa_T03_burst_yields'] #['AGB_burst_yields_Z0.0200','SNII_burst_yields_Z0.0200','SNIa_burst_yields_Z0.0200']
BinaryC_burstyieldsfile = ['AGB_'+bcYields_tag+'_burst_yields','SNII_'+bcYields_tag+'_burst_yields','SNIa_'+bcYields_tag+'_burst_yields']
if PLOT_SINGLE_AND_BINARY_STARS :
    BinaryC_burstyieldsfile2 = ['AGB_'+bcYields_tag2+'_burst_yields','SNII_'+bcYields_tag2+'_burst_yields','SNIa_'+bcYields_tag2+'_burst_yields']
#Total yields & SN numbers:
ensemble_Zs_string = ["{0:0.4f}".format(i) for i in ensemble_file_metallicities] #ensemble_file_metallicities.astype(str) #['0.0001','0.0010','0.0040','0.0080','0.0100','0.0300']
original_comp_Z = ['0.0004','0.0004','0.0040','0.0080','0.0200','0.0500']
orig_yield_ref = ["Marigo (2001)","Portinari et al. (1998)","Thielemann et al. (2003)"] #["Marigo+01","Portinari+98","Thielemann+03"]
orig_GCE_ref = "Yates et al. (2013)" #(Yates+13) #(Henriques+20)
totalYields_folder = r"C:\Users\ry22aas\robyates\Astro\L-Galaxies\L-Galaxies2020_plusBinaries\L-Galaxies2020_plusBinaries_version\YieldTables\TotalYields\\" #basedir+"/TotalYields/" #"/comparing_old_and_new_yields/"
LGals_timestepsfile = 'lookbackTimes_Myr' #'LGals_lookbackTimes_Myr'
OriginalYields_EjectedMassesFile = ['AGB_M01_TotalEjectedMass','SNII_P98_TotalEjectedMass','SNIa_T03_TotalEjectedMass']
OriginalYields_TotalMetalsFile = ['AGB_M01_TotalMetals','SNII_P98_TotalMetals','SNIa_T03_TotalMetals']
OriginalYields_SNNumFile = ['AGB_M01_TotalSNNum','SNII_P98_TotalSNNum','SNIa_T03_TotalSNNum']
bcYields_EjectedMassesFile = ['AGB_'+bcYields_tag+'_TotalEjectedMass','SNII_'+bcYields_tag+'_TotalEjectedMass','SNIa_'+bcYields_tag+'_TotalEjectedMass'] #['LGals_AGBMassEjec_'+bcYields_tag+'','LGals_SNIIMassEjec_'+bcYields_tag+'','LGals_SNIaMassEjec_'+bcYields_tag+'']
bcYields_TotalMetalsFile = ['AGB_'+bcYields_tag+'_TotalMetals','SNII_'+bcYields_tag+'_TotalMetals','SNIa_'+bcYields_tag+'_TotalMetals']
bcYields_SNNumFile = ['AGB_'+bcYields_tag+'_TotalSNNum','SNII_'+bcYields_tag+'_TotalSNNum','SNIa_'+bcYields_tag+'_TotalSNNum']
if PLOT_SINGLE_AND_BINARY_STARS :
    bcYields_EjectedMassesFile2 = ['AGB_'+bcYields_tag2+'_TotalEjectedMass','SNII_'+bcYields_tag2+'_TotalEjectedMass','SNIa_'+bcYields_tag2+'_TotalEjectedMass'] #['LGals_AGBMassEjec_'+bcYields_tag2+'','LGals_SNIIMassEjec_'+bcYields_tag2+'','LGals_SNIaMassEjec_'+bcYields_tag2+'']
    bcYields_TotalMetalsFile2 = ['AGB_'+bcYields_tag2+'_TotalMetals','SNII_'+bcYields_tag2+'_TotalMetals','SNIa_'+bcYields_tag2+'_TotalMetals']
    bcYields_SNNumFile2 = ['AGB_'+bcYields_tag2+'_TotalSNNum','SNII_'+bcYields_tag2+'_TotalSNNum','SNIa_'+bcYields_tag2+'_TotalSNNum']
#Ensemble suite:
SQLdir = ensembledir #basedir+"ensemble_suite\\"
SQL_file = "ensembles.sql"
#Unprocessed comp:
Original_metallicities = np.array([0.0004,0.0040,0.0080,0.0200,0.0500,1.0000]) 
#Version:
mark = "_mk126" #"" #"_mk105_DDonlyx10000" #


##########
#Input variables:
##########
if AL26 == 1 :
    #element_names = ['Hydrogen','Helium','Carbon','Nitrogen','Oxygen','Neon','Magnesium','Silicon','Sulphur','Al26','Iron'] #['Iron'] #
    elements = ['H','He','C','N','O','Ne','Mg','Si','S','Al26','Fe']
else :
    #element_names = ['Hydrogen','Helium','Carbon','Nitrogen','Oxygen','Neon','Magnesium','Silicon','Sulphur','Calcium','Iron'] #['Iron'] #
    #elements = ['H1|H2','He3|He4','C1','N1|N2','O1','Ne','Mg','Si','S3','Ca','Fe']
    elements = ['H','He','C','N','O','Ne','Mg','Si','S','Ca','Fe'] #['H','He','C','N','O','Ne','Mg','Al','Si','S','Ca','Mn','Fe'] #
ele_to_plot = 'C12' #'elements[2] # 'N14' #If PLOT_ONE_ELEMENT is on, ele_to_plot defines which element it is that's plotted from the binary_c yields.
selec_channels = {
    'AGB'  : ['AGB', 'RLOF', 'Wind', 'TZ', 'GB', 'Comenv'],
    'SNII' : ['SNII', 'SNIIa', 'SNIbc', 'WR', 'SNeCAP', 'SNAIC'],
    'SNIa' : ['SNIa_ELD', 'SNIa_CHAND_Coal', 'SNIa_CHAND', 'Novae']
}
selec_channel_types = list(selec_channels.keys())
channel_names = ['Wind group', 'SN-II group', 'SN-Ia group']
all_channels = np.sort(np.array(['Wind', 'WR', 'Comenv', 'SNIbc', 'TZ', 'AGB', 'SNII', 'GB', 'SNAIC', 'SNeCAP', 'Novae', 'SNIa_CHAND', 'SNIa_ELD', 'SNIIa', 'SNIa_CHAND_Coal', 'RLOF']))
if ORIGINAL_ENSEMBLES == 1 :
    #NOTE!: Some channel names are witten differently in the SN scalars compared to the yield sources! E.g. 'SNIa_CHAND_Coal' vs. 'SN_IA_CHAND_Coal'
    sn_channels = { #As listed in ensemble_data['scalars']
        # Old names:
        'AGB'  : ['SN_TZ'], #'TPAGB'
        'SNII' : ['SN_II', 'SN_IIa', 'SN_IBC', 'SN_AIC', 'SN_AIC_BH'], #Rob I says to remove these: 'SN_BH_NS', 'SN_BH_BH' 
        'SNIa' : ['SN_IA_CHAND', 'SN_IA_CHAND_Coal', 'SN_IA_ELD'] #, 'SN_IA_He', 'SN_IA_He_Coal', 'SN_IA_Hybrid_HeCOWD', 'SN_IA_Hybrid_HeCOWD_subluminous'] #, 'WDWD_DD' 
    }
else :
    sn_channels = { #As listed in ensemble_data['scalars']
        # New names:
        'AGB'  : ['SN_TZ'],
        'SNII' : ['SN_II', 'SN_IIa', 'SN_IBC', 'SN_AIC', 'SN_AIC_BH', 'SN_II_BRIGHT'], #Is SN_II_BRIGHT a sub-set of 'SN+II'?... #Rob I says to remove these: 'SN_BH_NS', 'SN_BH_BH' 
        'SNIa' : ['SN_ChandIa', 'SN_CHAND_Coal', 'SN_ELDIa'] #, 'WDWD_DD'] #, 'SN_IA_He', 'SN_IA_He_Coal', 'SN_IA_Hybrid_HeCOWD', 'SN_IA_Hybrid_HeCOWD_subluminous']
    }
#sn_channel_types = list(sn_channels.keys())


##########
# Read binary_c SQL file describing the suite of ensembles available:
##########
if READ_SQL_FILE == 1 :
    # Load file
    db = sqlite3.connect(SQLdir+SQL_file)
    
    # Check which tables exist
    #tables = pd.read_sql_query("SELECT name FROM sqlite_master WHERE type='table'", db)
    
    # Check the ensemble table
    ensemble_table = pd.read_sql("SELECT * FROM ensembles", db)
    
    # Find those that are finished
    finished = ensemble_table[ensemble_table.status=='finished']
    
    # Dictionary containing the values that we want.
    value_dict = {
        'b_use_periastron_Roche_radius': False,
        'b_post_ce_adaptive_menv': False,
        'b_no_thermohaline_mixing': '1',
        'b_alpha_ce': '1',
        'b_lambda_ce': '0.5',
        'b_RLOF_method': '0',
        'b_E2_prescription': '0',
        'b_gb_reimers_eta': '0.4',
        'b_sn_kick_dispersion_II': '190',
        'b_WRLOF_method': '0',
        'b_wind_mass_loss': '3',
        'b_tpagbwind': '0',
        'b_sn_kick_dispersion_IBC': '190',
        'b_gbwind': '0',
        'b_CRAP_parameter': '0',
        'b_rotationally_enhanced_mass_loss': '0',
        'b_BH_prescription': '0',
        #'b_rotationally_enhanced_mass_loss': '0',
        'b_lambda_ionisation': '0.5',
        'b_mass_accretion_for_eld': '0.15',
        'b_type_Ia_MCh_supernova_algorithm': '0',
        'b_tidal_strength_factor': '1',
        'b_mc13_pocket_multiplier': '1',
        'b_lambda_min': '0',
        'b_delta_mcmin': '0',
        'e_binaries': True
    }
    
    # Filter the dataframe
    filtered_df = finished.iloc[:]
    for column in value_dict:
        filtered_df = filtered_df.query("{} == '{}'".format(column, value_dict[column]))
    
    # Convert to list
    sorted_and_selected = filtered_df.sort_values(by='b_metallicity')[['b_metallicity', 'e_uuid', 'e_outdir']]
    print(sorted_and_selected)
    list_form = sorted_and_selected.values.tolist()
    print(list_form)
    

if LOAD_DATA == 1 : 
    print("\n#################")
    print("BINARY_C")
    print("#################\n")
    print("Original ensembles: ", ORIGINAL_ENSEMBLES)
    # Read simulation info:
    if READ_SIM_INFO == 1:
        with open(ensembledir+simulation_info_file+'.json', 'r', encoding="utf-8") as f:
            sim_info = json.loads(f.read())
        el_dict = sim_info['binary_c_version_info']['elements'] #Load element info
        el_df = pd.DataFrame.from_dict(el_dict) #Convert dictionary to data frame
        el_names = el_df.loc['name',:] #Get element names
        el_isoNum = el_df.loc['amt_isotopes',:] #Get isotope numbers
        el_names.index = el_names.index.astype(int) #Convert indices (i.e. atomic numbers) from strings to integers
        #pd.set_option('display.max_rows', el_df.shape[1]+1) #Set to print all rows (not just the first and last N rows)
        #Print info:
        print("Binary_c version:", sim_info['binary_c_version_info']['miscellaneous']['version'])
        print("Binary_c build:", sim_info['binary_c_version_info']['miscellaneous']['build'], "\n")
        if SINGLE_STARS == 1 :
            print("SP type: Single stars\n")
        else :
            print("SP type: Binary stars\n")
        print("Elements included:")
        print(np.array(el_names.sort_index())) #Print elements (in order of atomic number) as an array
        print("Total elements: "+str(len(el_names)-1))
        print("Total isotopes: "+str(el_isoNum.sum()-1))
    
    ##########
    # LOOP OVER THE BINARY_C METALLICITIES AVAILABLE:
    ##########
    loopa = 0 #Is 0 for first loop over metalicities, and 1 thereafter
    the_metallicities = np.array([]) #Array to stor the metallicity files actually loaded (for printing to metadata).
    if SINGLE_METALLICITY == 1 :
        zz_start = int(np.where(ensemble_file_metallicities == singleZ)[0])
        zz_stop = int(np.where(ensemble_file_metallicities == singleZ)[0])+1
        #for zz in range(np.where(ensemble_file_metallicities == singleZ),np.where(ensemble_file_metallicities == singleZ)+1) :
    else :      
        zz_start = 0
        zz_stop = len(ensemble_file_metallicities)
        #for zz in range(len(ensemble_file_metallicities)) :            
    for zz in range(zz_start,zz_stop) :
        ##########
        # Load metadata:
        ##########
        #with open(ensemble_datafile+'-100'+'.json', 'r', encoding="utf-8") as f : meta_data = json.loads(f.read())['metadata']
        with open(ensembledir+IMFdir+ensemble_datafile+'-'+ensemble_file_codes[zz]+'.json', 'r', encoding="utf-8") as f:
            meta_data = json.loads(f.read())['metadata']
        
          
        ##########
        # Load ensemble data:
        ##########
        with open(ensembledir+IMFdir+ensemble_datafile+'-'+ensemble_file_codes[zz]+'.json', 'r', encoding="utf-8") as f:
            ensemble_data = json.loads(f.read())['ensemble']
        
        #Get key parameters:
        bc_metallicity = str(meta_data['settings']['population_settings']['bse_options']['metallicity']) #N.B. might have to remove the starting "0." from the metallicity, for L-Galaxies to read in the filename properly? No, it seems fine (02-05-22)
        the_metallicities = np.append(the_metallicities,bc_metallicity)
        #N.B. If that metallicity cannot be found in 'bse_options', it may be because this ensemble has the default metallicity, which can be found in: meta_data['settings']['binary_c_defaults']['metallicity'] (16-05-22)
        bc_total_probability_weighted_mass = meta_data['total_probability_weighted_mass'] #This is M_av,sp, the "arithmetic mean mass of stellar systems in the population" to which the yields are normalised.
        if RESCALED_SINGLE_STAR_POP == 1 :
            with open(ensembledir2+IMFdir+ensemble_datafile+'-'+ensemble_file_codes2[zz]+'.json', 'r', encoding="utf-8") as f2:
                bc_total_probability_weighted_mass2 = json.loads(f2.read())['metadata']['total_probability_weighted_mass']
        bc_timsteps = np.sort(np.array(list(ensemble_data['Xyield']['time'].keys()), dtype=float))
        bc_dlogt = (bc_timsteps[-1]-bc_timsteps[0])/(len(bc_timsteps)-1) #This is assuming you're reading in an ensemble with log timesteps
        print("\n***********")
        print("Z = "+str(bc_metallicity)+":")
        print("***********")
        print("Ensemble file: "+ensemble_file_codes[zz])
        print("log timesteps?: "+str(bool(meta_data['settings']['population_settings']['bse_options']['ensemble_logtimes']))+" ("+str(meta_data['settings']['population_settings']['bse_options']['ensemble_logtimes'])+str(")"))
        print("log(dt/Myr) = "+str(meta_data['settings']['population_settings']['bse_options']['ensemble_logdt']))
        print("Max evo time (Myr) = "+str(meta_data['settings']['population_settings']['bse_options']['max_evolution_time']))
        print("Multiplicity = "+str(meta_data['settings']['population_settings']['bse_options']['multiplicity'])+"  (1 = single stars, 2 = binaries, 3 = triples, 4 = quadruples)")
        print("Binaries = "+str(meta_data['settings']['population_settings']['custom_options']['binaries']))
        print("Mmax (Msun) = "+str(meta_data['settings']['population_settings']['custom_options']['mmax']))
        print("SNIa_MCh algorithm = "+str(meta_data['settings']['binary_c_defaults']['type_Ia_MCh_supernova_algorithm'])+" (0 = DD7 Iwamoto 1999, 1 = Seitenzahl 2013 3D hydro yields)")
        if ORIGINAL_ENSEMBLES != 1 :
            print("IMF = "+str(meta_data['settings']['population_settings']['population_options']['Moe2017_options']['IMF_distribution']))
            print("Initial abundance mix = "+str(meta_data['settings']['population_settings']['bse_options']['initial_abundance_mix']))
            print("BH prescription = "+str(meta_data['settings']['population_settings']['bse_options']['BH_prescription']))
            print("qcrit_nuclear_burning = "+str(meta_data['settings']['population_settings']['bse_options']['qcrit_nuclear_burning']))
        print("\n")
        # print([(key, value) for key, value in meta_data['settings']['binary_c_defaults'].items() if key.startswith("qcrit")])
        # print([(key, value) for key, value in meta_data['settings']['binary_c_defaults'].items() if "qcrit" in key])     
        
        # print(ensemble_data.keys())
        yield_data = ensemble_data['Xyield']
        if (WRITE_SN_RATES) | (PLOT_BINARYC_YIELDS) :
            sn_data = ensemble_data['scalars']
            sn_rates = {}
            # print([i for i in sn_data.keys() if i == 'SN_IIa']) #Check if a certain SN type is in the ensemble.
            # print([i for i in sn_data.keys() if 'IIa' in i]) #Check if a certain SN type is in the ensemble.
    
        # First I inflate it to a list of lists
        data_list = inflate_ensemble_with_lists(yield_data)
        
        # Load it as an array and rotate
        data_array = np.array(data_list).T
        
        # I fetch the names of the columns like so:
        columnames = find_columnames_recursively(yield_data)
        
        # Which I update with a final column name
        columnames = columnames + ['yield_per_solarmass'] 
        
        # Then I load it into a dataframe
        df = pd.DataFrame(data_array, columns=columnames)
        
        # Make sure we set the correct type in the columns
        # (because we turn the list into an array and transpose it, the whole array is filled with strings rather than a mix of both)
        df = df.astype({"time": float, "yield_per_solarmass": float})

    
    #########################################################
        if WRITE_BINARYC_YIELDS :
            df_blank = df.groupby(by='time')['yield_per_solarmass'].sum()*0.0 #A dataframe of all times and a set of 0.0s for yield_per_solarmass (for use later to ensure times with no element ejection are not missed-out)
            #sn_df_blank = df.groupby(by='time')['yield_per_solarmass'].sum() / df.groupby(by='time')['yield_per_solarmass'].sum()
            
            # Boost all DD scenario SNIa yields by SN_IA_CHAND_COAL_FACTOR:
            #df[df['source']=='SNIa_CHAND_Coal']['yield_per_solarmass'] = df[df['source']=='SNIa_CHAND_Coal']['yield_per_solarmass'] * SN_IA_CHAND_COAL_FACTOR
            #print("1) ", df[df['source']=='SNIa_CHAND_Coal']['yield_per_solarmass'])
            dfA = df['yield_per_solarmass']
            dfB = dfA[df['source']=='SNIa_CHAND_Coal']
            dfB *= SN_IA_CHAND_COAL_FACTOR
            for idx in dfB.keys() :
                df['yield_per_solarmass'].loc[idx] = dfB[idx]
            #print("2) ", df[df['source']=='SNIa_CHAND_Coal']['yield_per_solarmass'])
            
            for group in selec_channel_types :
                source_df = df[df['source'].isin(selec_channels[group])]
                print(group, " group:")
                print("Yield sources requested: ", np.sort(selec_channels[group]))
                print("Yield sources available: ", np.sort(source_df['source'].unique()))
                if WRITE_SN_RATES == 1 :
                    print("SN sources requested: ", np.sort(sn_channels[group]))
                    sn_avail = []
                    for sn_channel in sn_channels[group] :
                        try :
                            sn_avail.append([i for i in sn_data.keys() if sn_channel == i][0])
                        except :
                            yy = 0.0
                    print("SN sources available: ", np.sort(sn_avail), "\n")
                #####
                if AL26 == 1 :
                    Al26df = source_df[source_df['isotope'].str.contains('Al26')]
                    grouped_by_source = Al26df.groupby(by='source')['yield_per_solarmass'].sum()    
                    print(group+' group:', 'Al26', grouped_by_source) #ensemble_file_metallicities[zz],    
                #####
                #(a) element yields for each channel (AGB, SNe-II, SNe_Ia):
                for element in elements : 
                    ##ele = (df["isotope"].str.contains(element)) & (((df["isotope"].str.replace(element,'')).str.isnumeric()) | (df["isotope"].str.replace(element,'') == '')) #(16-05-23): This compound condition selects only those entries in the dataframe related to the element in question. In detail, it selects entries where the isotope name contains only the given element name (e.g. H, He, etc) plus numbers (or contains only the given isotope name [Al26, H34, etc])
                    ##element_df = source_df[(df["isotope"].str.contains(element)) & (((df["isotope"].str.replace(element,'')).str.isnumeric()) | (df["isotope"].str.replace(element,'') == ''))] 
                    # New way:
                    ele = (source_df["isotope"].str.contains(element)) \
                        & (((source_df["isotope"].str.replace(element,'')).str.isnumeric()) | (source_df["isotope"].str.replace(element,'') == '')) #(16-05-23): This compound condition selects only those entries in the dataframe related to the element in question. In detail, it selects entries where the isotope name contains only the given element name (e.g. H, He, etc) plus numbers (or contains only the given isotope name [Al26, H34, etc])
                    element_df = source_df[ele]
                    # Old way:
                    #element_df = source_df[source_df['isotope'].str.contains(element)]
                    #print(element, element_df["isotope"].shape)
                    grouped_by_time = element_df.groupby(by='time')['yield_per_solarmass'].sum()           
                    if CONVERT_DATA_TO_LINEAR_TIME_BINS == 1 :
                        grouped_by_time = grouped_by_time / (10**(np.array(grouped_by_time.index, dtype=float)) * np.log(10))
                    if CONVERT_DATA_TO_PER_SOLARMASS == 1 :
                        grouped_by_time = grouped_by_time / bc_total_probability_weighted_mass
                    if RESCALED_SINGLE_STAR_POP == 1 :
                        grouped_by_time = grouped_by_time * bc_total_probability_weighted_mass/bc_total_probability_weighted_mass2 #Accounting for different scaling of binary+single star ensembles and single-star-component ensembles
                    the_yield_array = np.nan_to_num(np.array(df_blank+grouped_by_time))
                    the_yield_list = list(np.array(the_yield_array, dtype=str))
                    if (element == elements[0]) : 
                        writemode = 'w'
                        tot_ej_masses_array = the_yield_array
                        #print("1: "+element)
                    else : 
                        writemode = 'a'
                        tot_ej_masses_array += the_yield_array
                        #print("2: "+element)
                        if (element == elements[2]) :
                            tot_metal_masses_array = the_yield_array              
                        elif (element != elements[1]) :
                            tot_metal_masses_array += the_yield_array                        
                    with open(bc_yieldsdir+bc_yieldsfile+'_'+group+'_'+'Z'+bc_metallicity+'_'+'Yields'+mark+'.txt', writemode) as f:
                        f.write(' '.join(the_yield_list)+' ')
                        if NEW_LINE_PER_ELEMENT == 1 :
                            f.write('\n')
                #(b) total ejected masses for each channel (AGB, SNe-II, SNe_Ia):
                the_tot_ej_masses_list = list(np.array(tot_ej_masses_array, dtype=str))
                with open(bc_yieldsdir+bc_yieldsfile+'_'+group+'_'+'Z'+bc_metallicity+'_'+'EjectedMasses'+mark+'.txt', 'w') as f:
                     f.write(' '.join(the_tot_ej_masses_list))
                #(c) total ejected metal masses for each channel (AGB, SNe-II, SNe_Ia):
                the_tot_metal_masses_list = list(np.array(tot_metal_masses_array, dtype=str))
                with open(bc_yieldsdir+bc_yieldsfile+'_'+group+'_'+'Z'+bc_metallicity+'_'+'TotalMetals'+mark+'.txt', 'w') as f:
                     f.write(' '.join(the_tot_metal_masses_list))
                if WRITE_SN_RATES == 1 :
                    # (d) SN rates for each channel (AGB, SNe-II, SNe_Ia):
                    #bc_timesteps = list(np.array(df_blank.index, dtype=str))
                    bc_timesteps = np.array(df_blank.index, dtype=str)
                    #Convert whole-number floats to integers in bc_timesteps, to match sn_source keys:
                    for bc_timestep in bc_timesteps :
                        if (np.mod(float(bc_timestep),1) == 0.0) :
                            bc_timesteps[np.where(bc_timesteps == bc_timestep)] = str(int(float(bc_timestep)))
                    bc_timesteps = list(bc_timesteps)
                    sn_rates = dict.fromkeys(bc_timesteps, 0.0)
                    for sn_source in sn_channels[group] :
                        if sn_source in sn_data :
                            for bc_timestep in bc_timesteps :
                                if bc_timestep in sn_data[sn_source] :
                                    #sn_rates[bc_timestep] += sn_data[sn_source][bc_timestep]
                                    rate = sn_data[sn_source][bc_timestep]
                                    if CORRECT_SN_RATE_UNITS == 1 :
                                        #sn_rates[bc_timestep] /= bc_dlogt
                                        rate /= bc_dlogt
                                    if CONVERT_DATA_TO_LINEAR_TIME_BINS == 1 :
                                        #sn_rates[bc_timestep] /= 10**(float(bc_timestep)) * np.log(10)
                                        rate /= 10**(float(bc_timestep)) * np.log(10)
                                    if CONVERT_DATA_TO_PER_SOLARMASS == 1 :
                                        #sn_rates[bc_timestep] /= bc_total_probability_weighted_mass
                                        rate /= bc_total_probability_weighted_mass
                                    if RESCALED_SINGLE_STAR_POP == 1 :
                                        rate *= bc_total_probability_weighted_mass/bc_total_probability_weighted_mass2 #Accounting for different scaling of binary+single star ensembles and single-star-component ensembles
                                    if sn_source == 'SN_IA_CHAND_Coal' :
                                        #print("SN_IA_CHAND_Coal: ", rate, rate * SN_IA_CHAND_COAL_FACTOR)
                                        rate *= SN_IA_CHAND_COAL_FACTOR #Multiply DD scenario SN-Ia rate by arbirtrary factor to match Claeys+14 rates.
                                    sn_rates[bc_timestep] += rate
                                    # if group == 'SNIa' :
                                    #     print(str(sn_source)+" | "+str(bc_timestep)+" | "+str(rate*bc_dlogt)+" | "+str(rate)+" | "+str(sn_rates[bc_timestep]))
                        else : print('N.B. SN type '+sn_source+' is not in '+ensemble_datafile+'-'+ensemble_file_codes[zz])
                    with open(bc_yieldsdir+bc_yieldsfile+'_'+group+'_'+'Z'+bc_metallicity+'_'+'Rates'+mark+'.txt', 'w') as f:
                        f.write(' '.join(np.array(list(sn_rates.values()), dtype=str)))
                    # if group == 'SNII' : 
                        # print("Data: bc_dlogt: ", bc_dlogt)
                        # #print("Data: SNII: ", sn_rates.values())
                        # print("Data: SNII: ", sn_rates)
                        # print("Data: Sum SNII: ", sum(sn_rates.values()))
            
            #(e) timesteps:
            with open(bc_yieldsdir+'ind_Z_timestep_files/'+bc_yieldsfile+'_'+'Z'+bc_metallicity+'_'+'logt_timesteps'+mark+'.txt', 'w') as f:
                f.write(' '.join(list(np.array(df_blank.index, dtype=str))))
            # if (zz > 0) & (SINGLE_METALLICITY != 1) :
            #     newfile = bc_yieldsdir+'ind_Z_timestep_files/'+bc_yieldsfile+'_'+'Z'+str(ensemble_file_metallicities[zz])+'_'+'logt_timesteps'+mark+'.txt'
            #     prevfile = bc_yieldsdir+'ind_Z_timestep_files/'+bc_yieldsfile+'_'+'Z'+str(ensemble_file_metallicities[zz-1])+'_'+'logt_timesteps'+mark+'.txt'
            #     if filecmp.cmp(newfile, prevfile, shallow=False) == True :
            #         with open(bc_yieldsdir+bc_yieldsfile+'_'+'logt_timesteps'+mark+'.txt', 'w') as f:
            #             f.write(' '.join(list(np.array(df_blank.index, dtype=str))))
            #     else : 
            #         print('N.B. Timesteps for '+ensemble_datafile+'-'+ensemble_file_codes[zz-1]+' (Z='+str(ensemble_file_metallicities[zz-1])+') do not match '+ensemble_datafile+'-'+ensemble_file_codes[zz]+' (Z='+str(ensemble_file_metallicities[zz])+')')
            if (zz > 0) & (SINGLE_METALLICITY != 1) :
                newfile = bc_yieldsdir+'ind_Z_timestep_files/'+bc_yieldsfile+'_'+'Z'+str(the_metallicities[zz])+'_'+'logt_timesteps'+mark+'.txt'
                prevfile = bc_yieldsdir+'ind_Z_timestep_files/'+bc_yieldsfile+'_'+'Z'+str(the_metallicities[zz-1])+'_'+'logt_timesteps'+mark+'.txt'
                if filecmp.cmp(newfile, prevfile, shallow=False) == True :
                    with open(bc_yieldsdir+bc_yieldsfile+'_'+'logt_timesteps'+mark+'.txt', 'w') as f:
                        f.write(' '.join(list(np.array(df_blank.index, dtype=str))))
                else : 
                    print('N.B. Timesteps for '+ensemble_datafile+'-'+ensemble_file_codes[zz-1]+' (Z='+str(the_metallicities[zz-1])+') do not match '+ensemble_datafile+'-'+ensemble_file_codes[zz]+' (Z='+str(the_metallicities[zz])+')')
            
            #(f) metallicities:
            with open(bc_yieldsdir+bc_yieldsfile+'_'+'metallicities'+mark+'.txt', 'w') as f:
                f.write(' '.join(list(np.array(ensemble_file_metallicities, dtype=str))))    
            
            # #(g) elements considered:
            # if (loopa == 0) : #Only need to write this once, as the chosen elements are the same for all metallicities
            #     with open(bc_yieldsdir+bc_yieldsfile+'_'+'elements'+mark+'.txt', 'w') as f:
            #         f.write(' '.join(list(np.array(elements, dtype=str))))
        
        ##########
        #PLOTS:
        ##########
        if ROBS_PLOT_LAYOUT :
            from robs_custom_mpl_settings import robs_load_mpl_rc, robs_load_mpl_rc_paper
            robs_load_mpl_rc_paper()
        default_cols = plt.rcParams['axes.prop_cycle'].by_key()['color'] # These are the default colours python uses when plotting, in order.
        #########################################################
        if PLOT_BINARYC_YIELDS :
            ##########
            #Print binary_c plots to pdf:
            ##########
            pdf = PdfPages(basedir+"/plots/"+ensemble_datafile+'_Z'+bc_metallicity+mark+".pdf") #, bbox_inches='tight'
            xlimits = [-1.0, 4.4]
                   
            ##########
            #Plotting the total of all elements ejected by each channel as a function of time:
            ##########
            if PLOT_METAL_ELEMENTS == 1 :
                ylimits = [-12.,-2.] #[-13.,-2.]
            else :
                ylimits = [-12.,-1.] #[1.e-12, 1.0]
            #fig = plt.figure(figsize=(15,15))
            fig = plt.figure() #figsize=(6,6)
            counta=0
            
            #Plot total line
            #grouped_by_time = df.groupby(by='time')['yield_per_solarmass'].sum()
            if PLOT_ONE_ELEMENT == 1 :
                ele = (df["isotope"].str.contains(ele_to_plot)) \
                    & (((df["isotope"].str.replace(ele_to_plot,'')).str.isnumeric()) | (df["isotope"].str.replace(ele_to_plot,'') == '')) #(16-05-23): This compound condition selects only those entries in the dataframe related to the element in question. In detail, it selects entries where the isotope name contains only the given element name (e.g. H, He, etc) plus numbers (or contains only the given isotope name [Al26, H34, etc])
                total_df = df[ele]  
                grouped_by_time = total_df.groupby(by='time')['yield_per_solarmass'].sum()
            else :
                if PLOT_METAL_ELEMENTS == 1 :
                    eles = elements[2:]
                else :
                    eles = elements
                for element in eles : 
                    ele = (df["isotope"].str.contains(element)) \
                        & (((df["isotope"].str.replace(element,'')).str.isnumeric()) | (df["isotope"].str.replace(element,'') == '')) #(16-05-23): This compound condition selects only those entries in the dataframe related to the element in question. In detail, it selects entries where the isotope name contains only the given element name (e.g. H, He, etc) plus numbers (or contains only the given isotope name [Al26, H34, etc])
                    total_df = df[ele]
                    if element == eles[0] :
                        grouped_by_time = total_df.groupby(by='time')['yield_per_solarmass'].sum()
                    else :
                        grouped_by_time2 = total_df.groupby(by='time')['yield_per_solarmass'].sum()
                        grouped_by_time_all = pd.concat((grouped_by_time,grouped_by_time2))
                        grouped_by_time = grouped_by_time_all.groupby(by='time').sum()             
            if CONVERT_PLOTS_TO_LINEAR_TIME_BINS == 1 :
                grouped_by_time = grouped_by_time / (10**(np.array(grouped_by_time.index, dtype=float)) * np.log(10))
            if CONVERT_PLOTS_TO_PER_SOLARMASS == 1 :
                grouped_by_time = grouped_by_time / bc_total_probability_weighted_mass  
            if RESCALED_SINGLE_STAR_POP == 1 :
                grouped_by_time = grouped_by_time * bc_total_probability_weighted_mass/bc_total_probability_weighted_mass2 #Accounting for different scaling of binary+single star ensembles and single-star-component ensembles
            plt.plot(np.array(grouped_by_time.index, dtype=float), np.log10(np.array(grouped_by_time)), \
                      label='Total', linestyle='--', color='black', linewidth=3., zorder=10)
            # plt.plot(np.array(grouped_by_time.index, dtype=float), np.array(grouped_by_time), \
            #           label='Total', linestyle='--', color='black', linewidth=5.)
            
            #Plot individual sources:
            #Get default colour cycle:
            prop_cycle = plt.rcParams['axes.prop_cycle']
            colours = np.array(np.append(prop_cycle.by_key()['color'],prop_cycle.by_key()['color'])) #np.array(prop_cycle.by_key()['color'])
            #cols = enumerate(colours)
            for source in np.sort(df['source'].unique()) :
            #for source in df['source'].unique():
            #for source in selec_channels['SNIa']:
                #print(source)
                source_df_all = df[df.source==source]
                col = colours[np.where(all_channels == source)]
                #grouped_by_time = source_df.groupby(by='time')['yield_per_solarmass'].sum()
                if PLOT_ONE_ELEMENT == 1 :
                    ele = (source_df_all["isotope"].str.contains(ele_to_plot)) \
                        & (((source_df_all["isotope"].str.replace(ele_to_plot,'')).str.isnumeric()) | (source_df_all["isotope"].str.replace(ele_to_plot,'') == '')) #(16-05-23): This compound condition selects only those entries in the dataframe related to the element in question. In detail, it selects entries where the isotope name contains only the given element name (e.g. H, He, etc) plus numbers (or contains only the given isotope name [Al26, H34, etc])
                    source_df = source_df_all[ele]  
                    grouped_by_time = source_df.groupby(by='time')['yield_per_solarmass'].sum()
                else :
                    if PLOT_METAL_ELEMENTS == 1 :
                        eles = elements[2:]
                    else :
                        eles = elements
                    for element in eles :
                        ele = (source_df_all["isotope"].str.contains(element)) \
                            & (((source_df_all["isotope"].str.replace(element,'')).str.isnumeric()) | (source_df_all["isotope"].str.replace(element,'') == '')) #(16-05-23): This compound condition selects only those entries in the dataframe related to the element in question. In detail, it selects entries where the isotope name contains only the given element name (e.g. H, He, etc) plus numbers (or contains only the given isotope name [Al26, H34, etc])
                        source_df = source_df_all[ele]
                        if element == eles[0] :
                            grouped_by_time = source_df.groupby(by='time')['yield_per_solarmass'].sum()
                        else :
                            #pd.concat((source_df,source_df_all[ele])).groupby('yield_per_solarmass',as_index=False).sum()
                            grouped_by_time2 = source_df.groupby(by='time')['yield_per_solarmass'].sum()
                            #grouped_by_time += grouped_by_time2
                            #pd.concat((grouped_by_time,grouped_by_time2)).groupby(by='time').sum()
                            grouped_by_time_all = pd.concat((grouped_by_time,grouped_by_time2))
                            grouped_by_time = grouped_by_time_all.groupby(by='time').sum()
                if CONVERT_PLOTS_TO_LINEAR_TIME_BINS == 1 :
                    grouped_by_time = grouped_by_time / (10**(np.array(grouped_by_time.index, dtype=float)) * np.log(10))
                if CONVERT_PLOTS_TO_PER_SOLARMASS == 1 :
                    grouped_by_time = grouped_by_time / bc_total_probability_weighted_mass
                if RESCALED_SINGLE_STAR_POP == 1 :
                    grouped_by_time = grouped_by_time * bc_total_probability_weighted_mass/bc_total_probability_weighted_mass2 #Accounting for different scaling of binary+single star ensembles and single-star-component ensembles
                # if counta < 10 :
                #     plt.plot(np.array(grouped_by_time.index, dtype=float), np.log10(np.array(grouped_by_time)), \
                #              label=source.replace("_", "-"), marker='o', color=col[0])
                # else :
                #     plt.plot(np.array(grouped_by_time.index, dtype=float), np.log10(np.array(grouped_by_time)), \
                #              label=source.replace("_", "-"), marker='o', color=col[0], linestyle=':')
                if source in all_channels[0:10] :
                    plt.plot(np.array(grouped_by_time.index, dtype=float), np.log10(np.array(grouped_by_time)), \
                             label=source.replace("_", "-"), marker='o', color=col[0])
                else :
                    plt.plot(np.array(grouped_by_time.index, dtype=float), np.log10(np.array(grouped_by_time)), \
                             label=source.replace("_", "-"), marker='o', color=col[0], linestyle=':')
                counta+=1
                        
            #plt.yscale('log')
            plt.title(r"$Z$ = "+bc_metallicity)
            plt.xlabel(r'log$(t/\textnormal{Myr})$')
            if CONVERT_PLOTS_TO_LINEAR_TIME_BINS == 1 :
                if PLOT_ONE_ELEMENT :
                    plt.ylabel(r'$\textnormal{log}(\dot{M}_{\rm '+ele_to_plot+r',norm,SP}/\textnormal{Myr}^{-1})$')
                elif PLOT_METAL_ELEMENTS == 1 :
                    plt.ylabel(r'$\textnormal{log}(\dot{M}_{\rm Z,norm,SP}/\textnormal{Myr}^{-1})$')
                else :
                    plt.ylabel(r'$\textnormal{log}(\dot{M}_{\rm norm,SP}/\textnormal{Myr}^{-1})$')
            else :
                plt.ylabel(r'Metal yield $\ [\textnormal{log}(\textnormal{d}M/\textnormal{dlog}(t)\ 1/\textnormal{M}_{\odot})]$')
            plt.xlim(xlimits)
            plt.ylim(ylimits)
            #if zz == zz_stop-1 :
            if PLOT_METAL_ELEMENTS == 1 :
                plt.legend(labelspacing=0.15)
            else :
                plt.legend()
            plt.show()
            pdf.savefig(fig, bbox_inches='tight')  
            pdf.close()
    
    
    #########################################################
        if PLOT_LGALAXIES_YIELDS :
            ##########
            #Read L-Galaxies yield data:
            ##########
            #(np.abs(ensemble_file_metallicities - LGalaxies_metallicities[3])).argmin()
            LGals_dts = np.loadtxt(burstYields_folder+LGalaxies_dtsfile+".txt", delimiter=",") #These are the timestep widths of each LGals timestep in Myr
            LGals_times = (0.5*LGals_dts[0]) + np.loadtxt(burstYields_folder+LGalaxies_timestepsfile+".txt", delimiter=",") #These are the (linear) cosmic times from the start of the L-Galaxies simulation to the midlle of each timestep in Myr
            
            ##########
            #Print binary_c plots to pdf:
            ##########
            pdf = PdfPages(basedir+"/plots/"+"L-Galaxies_yields"+'_Z'+LGalaxies_metallicities[zz]+mark+".pdf")
            fig = plt.figure(figsize=(40,60)) #figsize=(40,60)
            rowno = 6
            colno = 2
        
        
            ##########
            #Plotting the total element ejected by 'AGB', 'SNII', and 'SNIa' as a function of time:
            ##########
            LGals_yields_comb = np.loadtxt(burstYields_folder+LGalaxies_burstyieldsfile[0]+'_Z'+LGalaxies_metallicities[zz]+".txt") + \
                                np.loadtxt(burstYields_folder+LGalaxies_burstyieldsfile[1]+'_Z'+LGalaxies_metallicities[zz]+".txt") + \
                                np.loadtxt(burstYields_folder+LGalaxies_burstyieldsfile[2]+'_Z'+LGalaxies_metallicities[zz]+".txt")
            
            for element in range(len(elements)) :
                plt.subplot(rowno, colno, element+1)
                for group in range(len(selec_channel_types)) :
                    LGals_yields = np.loadtxt(burstYields_folder+LGalaxies_burstyieldsfile[group]+'_Z'+LGalaxies_metallicities[zz]+".txt")
                    plt.plot(np.log10(LGals_times), LGals_yields[element]/LGals_dts, label=selec_channel_types[group])
                    
                # Add total line (selected channels):
                plt.plot(np.log10(LGals_times), LGals_yields_comb[element,:]/LGals_dts, label='Selected total', linestyle=':', color='black', linewidth=5.)
                  
                if (element == 0) | (element == 1) : plt.title("Z = "+LGalaxies_metallicities[zz])
                plt.yscale('log')
                plt.xlabel(r'log(time/Myr)')
                #plt.ylabel(element_names[element]+r' yield')
                plt.ylabel(elements[element]+r' yield')
                plt.xlim(xlimits)
                plt.ylim(ylimits)
                plt.legend()
                #plt.show()  
                
            pdf.savefig(fig, bbox_inches='tight')  
            pdf.close()
        
    
    #########################################################
        if PLOT_YIELD_COMPARISONS == 1 : 
            LGals_timesteps = np.loadtxt(totalYields_folder+LGals_timestepsfile+".txt") #, delimiter=','
            LGals_dts = np.loadtxt(burstYields_folder+LGalaxies_dtsfile+".txt", delimiter=",") #These are the timestep widths of each LGals timestep in Myr
            age_of_Universe_z0 = LGals_timesteps[0] #13800.0 #Myr
            logt = np.log10(age_of_Universe_z0-LGals_timesteps)
            #SFHBin = 'All' #0 #5 #0 #13
            #Zbin = str(zz) #'Z0.0080' #'Z0.0200' #3   
            
            #Set-up multi-panel plot:
            xlimits = [-1.0, 4.4]
            ylimits = [[-6.,0.],[-6.,0.],[-9.5,-1.]] #[[-6.,0.],[-6.,0.],[-6.,-1.]]
            yticks = [[-5.,-4.,-3.,-2.,-1.],[-5.,-4.,-3.,-2.,-1.],[-9.,-8.,-7.,-6.,-5.,-4.,-3.,-2.]]
            
            if PLOT_SN_RATES == 1 :
                rows=3
            else :
                rows=2
            columns=1
            xlabs = [r'log$(t/\textnormal{Myr})$']*columns
            ylabs = [r'$\textnormal{log}(\dot{M}_{\rm tot}/\textnormal{M}_{\odot}\,\textnormal{yr}^{-1})$', \
                     r'$\textnormal{log}(\dot{M}_{\rm Z}/\textnormal{M}_{\odot}\,\textnormal{yr}^{-1})$', \
                     r'$\textnormal{log}(R_{\rm SN}/\textnormal{yr}^{-1})$']
            # xticks = [0.0,0.5,1.0]
            title = r'$Z$ = '+ensemble_Zs_string[zz]
            pdf = PdfPages(basedir+"/plots/"+"yields_integrals_comp"+'_Zorig'+original_comp_Z[zz]+'_Zbc'+ensemble_Zs_string[zz]+"_"+"multiPlot"+mark+".pdf")
            fig = plt.figure(figsize=(6,6*rows)) #figsize=(5,16)
            for ii in range(rows*columns) :
                panel = robs_plot_panels(ii, rows=rows, columns=columns, xlimits=xlimits, ylimits=ylimits, \
                                         xlab=xlabs, ylab=ylabs, yticks=yticks, title=title) #, xticks=xticks)
                if SINGLE_STARS == 1 :
                    bclabel = ["binary\_c (single stars only)"]
                    if PLOT_SINGLE_AND_BINARY_STARS == 1 :
                        bclabel = ["binary\_c (single stars only)", "binary\_c (mixed)"]
                else :
                    bclabel = ["binary\_c (mixed)"]
                    if PLOT_SINGLE_AND_BINARY_STARS == 1 :
                        bclabel = ["binary\_c (mixed)", "binary\_c (single stars only)"]
                    
                ##########
                # Total Ejected Masses per Msun formed:
                ##########
                if ii == 0 :
                    for group in range(len(selec_channel_types)) :
                        if MAXSNIIMASS25_YIELDS == 1 :
                            Original_yields = np.loadtxt(totalYields_folder+OriginalYields_EjectedMassesFile[group]+"_Z"+original_comp_Z[zz]+"_MaxSNIIMass25"+".txt")
                            plt.plot(logt, np.log10(Original_yields), linestyle=':')#, label=orig_yield_ref[group]+" (MaxSNIIMass25)") #*(LGals_dts*1.e6)
                        else :
                            Original_yields = np.loadtxt(totalYields_folder+OriginalYields_EjectedMassesFile[group]+"_Z"+original_comp_Z[zz]+".txt")
                            plt.plot(logt, np.log10(Original_yields), linestyle=':')#, label=orig_yield_ref[group]) #*(LGals_dts*1.e6)
                    plt.gca().set_prop_cycle(None)
                    for group in range(len(selec_channel_types)) :
                        bc_yields = np.loadtxt(totalYields_folder+bcYields_EjectedMassesFile[group]+"_Z"+ensemble_Zs_string[zz]+".txt")
                        plt.plot(logt, np.log10(bc_yields))#, label=bclabel[0]) #*(LGals_dts*1.e6)
                    if PLOT_SINGLE_AND_BINARY_STARS == 1 :
                        plt.gca().set_prop_cycle(None)
                        for group in range(len(selec_channel_types)) :
                            bc_yields = np.loadtxt(totalYields_folder+bcYields_EjectedMassesFile2[group]+"_Z"+ensemble_Zs_string[zz]+".txt")
                            plt.plot(logt, np.log10(bc_yields), linestyle='dashed')#, label=bclabel[1])
                    #Prep legend:
                    # plt.plot(np.nan, np.nan, linestyle=':', color='black', label=orig_GCE_ref)
                    # if PLOT_SINGLE_AND_BINARY_STARS == 1 :
                    #     plt.plot(np.nan, np.nan, linestyle='dashed', color='black', label=bclabel[1])
                    # plt.plot(np.nan, np.nan, linestyle='solid', color='black', label=bclabel[0])
                    # plt.legend()
                    # robs_plot_text(panel, [selec_channel_types[0]+' group', selec_channel_types[1]+' group', selec_channel_types[2]+' group'], \
                    #                vpos='top', hpos='left', colour=default_cols[0:len(selec_channels)])
                    robs_plot_text(panel, ['Wind group', 'SN-II group', 'SN-Ia group'], \
                                   vpos='top', hpos='left', colour=default_cols[0:len(selec_channels)])
                ##########
                # Total Metal Masses per Msun formed:
                ##########
                elif ii == 1 :
                    for group in range(len(selec_channel_types)) :
                        if MAXSNIIMASS25_YIELDS == 1 :
                            Original_yields = np.loadtxt(totalYields_folder+OriginalYields_TotalMetalsFile[group]+"_Z"+original_comp_Z[zz]+"_MaxSNIIMass25"+".txt")
                            plt.plot(logt, np.log10(Original_yields), linestyle=':')#, label=selec_channel_types[group]+" "+orig_yield_ref[group]+" (MaxSNIIMass25)")
                        else :
                            Original_yields = np.loadtxt(totalYields_folder+OriginalYields_TotalMetalsFile[group]+"_Z"+original_comp_Z[zz]+".txt")
                            plt.plot(logt, np.log10(Original_yields), linestyle=':')#, label=selec_channel_types[group]+" "+orig_yield_ref[group]) #*(LGals_dts*1.e6)
                    plt.gca().set_prop_cycle(None)
                    for group in range(len(selec_channel_types)) :
                        bc_yields = np.loadtxt(totalYields_folder+bcYields_TotalMetalsFile[group]+"_Z"+ensemble_Zs_string[zz]+".txt")
                        plt.plot(logt, np.log10(bc_yields))#, label=selec_channel_types[group]+bclabel[0]) #*(LGals_dts*1.e6)
                    if PLOT_SINGLE_AND_BINARY_STARS == 1 :
                        plt.gca().set_prop_cycle(None)
                        for group in range(len(selec_channel_types)) :
                            bc_yields = np.loadtxt(totalYields_folder+bcYields_TotalMetalsFile2[group]+"_Z"+ensemble_Zs_string[zz]+".txt")
                            plt.plot(logt, np.log10(bc_yields), linestyle='dashed')#, label=selec_channel_types[group]+bclabel[1])
                    #Prep legend:
                    plt.plot(np.nan, np.nan, linestyle=':', color='black', label=orig_GCE_ref)
                    if PLOT_SINGLE_AND_BINARY_STARS == 1 :
                        plt.plot(np.nan, np.nan, linestyle='dashed', color='black', label=bclabel[1])
                    plt.plot(np.nan, np.nan, linestyle='solid', color='black', label=bclabel[0])
                    plt.legend()
                ##########
                # Total SN number per Msun formed:
                ##########
                elif ii == 2 :
                    if PLOT_SN_RATES == 1 :
                        for group in range(len(selec_channel_types)) :
                            if group != 0 : #'AGB'
                                if MAXSNIIMASS25_YIELDS == 1 :  
                                    Original_SNNum = np.loadtxt(totalYields_folder+OriginalYields_SNNumFile[group]+"_Z"+original_comp_Z[zz]+"_MaxSNIIMass25"+".txt")
                                    plt.plot(logt, np.log10(Original_SNNum), linestyle=':')#, label=selec_channel_types[group]+" "+orig_GCE_ref+" (MaxSNIIMass25)")
                                else :
                                    Original_SNNum = np.loadtxt(totalYields_folder+OriginalYields_SNNumFile[group]+"_Z"+original_comp_Z[zz]+".txt")
                                    plt.plot(logt, np.log10(Original_SNNum), linestyle=':')#, label=selec_channel_types[group]+" "+orig_GCE_ref)
                            else : plt.plot(-99.,-99.)
                        plt.gca().set_prop_cycle(None)
                        for group in range(len(selec_channel_types)) :
                            if group != 0 :
                                bc_SNNum = np.loadtxt(totalYields_folder+bcYields_SNNumFile[group]+"_Z"+ensemble_Zs_string[zz]+".txt")
                                plt.plot(logt, np.log10(bc_SNNum))#, label=selec_channel_types[group]+bclabel[0])                    
                            else : plt.plot(-99.,-99.)
                        if PLOT_SINGLE_AND_BINARY_STARS == 1 :
                            plt.gca().set_prop_cycle(None)
                            for group in range(len(selec_channel_types)) :
                                if group != 0 :
                                    bc_SNNum = np.loadtxt(totalYields_folder+bcYields_SNNumFile2[group]+"_Z"+ensemble_Zs_string[zz]+".txt")
                                    plt.plot(logt, np.log10(bc_SNNum), linestyle='dashed')#, label=selec_channel_types[group]+bclabel[1])
                                else : plt.plot(-99.,-99.)
                        # plt.legend()                    
            pdf.savefig(fig, bbox_inches='tight')
            pdf.close()
                        
            ##########
            # Compare element yields generated by yields_integrals.c between old yields and binary_c yields:
            ##########
            ylimits = [-13.,-1.]
            pdf = PdfPages(basedir+"/plots/"+"yields_integrals_comp"+'_Zorig'+original_comp_Z[zz]+'_Zbc'+ensemble_Zs_string[zz]+"_"+"yields_alt"+mark+".pdf")
            LGals_dts = np.loadtxt(burstYields_folder+LGalaxies_dtsfile+".txt", delimiter=",") #These are the timestep widths of each LGals timestep in Myr
            LGals_times = (0.5*LGals_dts[0]) + np.loadtxt(burstYields_folder+LGalaxies_timestepsfile+".txt", delimiter=",") #These are the (linear) cosmic times from the start of the L-Galaxies simulation to the midlle of each timestep in Myr
            
            #fig = plt.figure(figsize=(16,18)) #figsize=(16,22) figsize=(16,24) figsize=(40,60)
            #fig, ax = plt.subplots(figsize=(16,18))
            if PLOT_ONE_COLUMN == 1:
                rowno = 11
                if PLOT_PIECHARTS == 1 :
                    fig = plt.figure(figsize=(14,32)) #figsize=(21,32)
                    colno = 2 #len(selec_channel_types)+1 #
                    #fig, ax = plt.subplots(rowno, colno, figsize=(14,32), gridspec_kw={'width_ratios': [3, 1]})
                else :
                    fig = plt.figure(figsize=(7,32)) #figsize=(16,22) figsize=(16,24) figsize=(40,60)
                    colno = 1
            else :
                fig = plt.figure(figsize=(16,18)) #figsize=(16,22) figsize=(16,24) figsize=(40,60)
                rowno = 6
                colno = 2
            
            if MAXSNIIMASS25_YIELDS == 1 :  
                LGals_yields_comb = np.loadtxt(burstYields_folder+LGalaxies_burstyieldsfile[0]+'_Z'+original_comp_Z[zz]+"_MaxSNIIMass25"+".txt") + \
                                    np.loadtxt(burstYields_folder+LGalaxies_burstyieldsfile[1]+'_Z'+original_comp_Z[zz]+"_MaxSNIIMass25"+".txt") + \
                                    np.loadtxt(burstYields_folder+LGalaxies_burstyieldsfile[2]+'_Z'+original_comp_Z[zz]+"_MaxSNIIMass25"+".txt")
            else :
                LGals_yields_comb = np.loadtxt(burstYields_folder+LGalaxies_burstyieldsfile[0]+'_Z'+original_comp_Z[zz]+".txt") + \
                                    np.loadtxt(burstYields_folder+LGalaxies_burstyieldsfile[1]+'_Z'+original_comp_Z[zz]+".txt") + \
                                    np.loadtxt(burstYields_folder+LGalaxies_burstyieldsfile[2]+'_Z'+original_comp_Z[zz]+".txt")
            bc_yields_comb = np.loadtxt(burstYields_folder+BinaryC_burstyieldsfile[0]+'_Z'+ensemble_Zs_string[zz]+".txt") + \
                                np.loadtxt(burstYields_folder+BinaryC_burstyieldsfile[1]+'_Z'+ensemble_Zs_string[zz]+".txt") + \
                                np.loadtxt(burstYields_folder+BinaryC_burstyieldsfile[2]+'_Z'+ensemble_Zs_string[zz]+".txt")
            
            normMass_lgal = np.zeros(3)
            totNormMass_lgal = 0.0
            totNormHHeMass_lgal = 0.0
            normMass_bc = np.zeros(3)
            totNormMass_bc = 0.0
            totNormHHeMass_bc = 0.0
            if PLOT_SINGLE_AND_BINARY_STARS == 1 :
                normMass_bc_2 = np.zeros(3)
                totNormMass_bc_2 = 0.0
                totNormHHeMass_bc_2 = 0.0
            for element in range(len(elements)) :
                panel = plt.subplot(rowno, colno, (element*colno)+1)
                #panel = plt.subplots(rowno, colno)
                for group in range(len(selec_channel_types)) :
                    if MAXSNIIMASS25_YIELDS == 1 :  
                        LGals_yields = np.loadtxt(burstYields_folder+LGalaxies_burstyieldsfile[group]+'_Z'+original_comp_Z[zz]+"_MaxSNIIMass25"+".txt")
                        plt.plot(np.log10(LGals_times), np.log10(LGals_yields[element]/LGals_dts), linestyle=':')#, label=selec_channel_types[group]+" "+orig_yield_ref[group]+" (MaxSNIIMass25)") #np.log10(LGals_times)
                        normMass_lgal[group] = np.trapz(LGals_yields[element]/LGals_dts, x=LGals_times) #Calculate total mass ejected of this element from this channel
                        totNormMass_lgal += normMass_lgal[group] #Cumulative total (normalised) mass of all elements ejected by all channels
                        if element <= 1 :
                            totNormHHeMass_lgal += normMass_lgal[group] #Cumulative total (normalised) mass of H and He ejected by all channels
                    else :
                        LGals_yields = np.loadtxt(burstYields_folder+LGalaxies_burstyieldsfile[group]+'_Z'+original_comp_Z[zz]+".txt")
                        plt.plot(np.log10(LGals_times), np.log10(LGals_yields[element]/LGals_dts), linestyle=':')#, label=selec_channel_types[group]+" "+orig_yield_ref[group])
                        normMass_lgal[group] = np.trapz(LGals_yields[element]/LGals_dts, x=LGals_times) #Calculate total mass ejected of this element from this channel
                        totNormMass_lgal += normMass_lgal[group] #Cumulative total (normalised) mass of all elements ejected by all channels
                        if element <= 1 :
                            totNormHHeMass_lgal += normMass_lgal[group] #Cumulative total (normalised) mass of H and He ejected by all channels
                        #print(elements[element], channel_names[group], normMass_lgal[group], totNormMass_lgal)
                plt.gca().set_prop_cycle(None)
                for group in range(len(selec_channel_types)) :
                    bc_yields = np.loadtxt(burstYields_folder+BinaryC_burstyieldsfile[group]+'_Z'+ensemble_Zs_string[zz]+".txt")
                    plt.plot(np.log10(LGals_times), np.log10(bc_yields[element]/LGals_dts))#, label=selec_channel_types[group]+" (binary\_c)")
                    normMass_bc[group] = np.trapz(bc_yields[element]/LGals_dts, x=LGals_times) #Calculate total mass ejected of this element from this channel
                    totNormMass_bc += normMass_bc[group] #Cumulative total (normalised) mass of all elements ejected by all channels
                    if element <= 1 :
                        totNormHHeMass_bc += normMass_bc[group] #Cumulative total (normalised) mass of H and He ejected by all channels
                if PLOT_SINGLE_AND_BINARY_STARS == 1 :
                    plt.gca().set_prop_cycle(None)
                    for group in range(len(selec_channel_types)) :
                        bc_yields = np.loadtxt(burstYields_folder+BinaryC_burstyieldsfile2[group]+'_Z'+ensemble_Zs_string[zz]+".txt")
                        plt.plot(np.log10(LGals_times), np.log10(bc_yields[element]/LGals_dts), linestyle='dashed')
                        normMass_bc_2[group] = np.trapz(bc_yields[element]/LGals_dts, x=LGals_times) #Calculate total mass ejected of this element from this channel
                        totNormMass_bc_2 += normMass_bc_2[group] #Cumulative total (normalised) mass of all elements ejected by all channels
                        if element <= 1 :
                            totNormHHeMass_bc_2 += normMass_bc_2[group] #Cumulative total (normalised) mass of H and He ejected by all channels
                    
                # # Add total line (selected channels):
                # plt.plot(np.log10(LGals_times), np.log10(LGals_yields_comb[element,:]/LGals_dts), label='Original total', linestyle=':', color='grey', linewidth=4.)
                # plt.plot(np.log10(LGals_times), np.log10(bc_yields_comb[element,:]/LGals_dts), label='binary_c total', color='black', linewidth=4.)
                    
                if (element == 0) | ((PLOT_ONE_COLUMN != 1) & (element == 1)) : plt.title(r'$Z$ = '+ensemble_Zs_string[zz]) #plt.title("Original Z: "+original_comp_Z[zz]+" $|$ binary\_c Z: "+ensemble_Zs_string[zz])
                if (element >= len(elements)-2) : plt.xlabel(r'log$(t/\textnormal{Myr})$')
                plt.ylabel(r'$\textnormal{log}(\dot{M}_{\rm '+elements[element]+r', burst}/\textnormal{M}_{\odot}\,\textnormal{yr}^{-1})$')
                plt.xlim(xlimits)
                plt.ylim(ylimits)
                if element < len(elements)-2 :
                    plt.setp(panel.get_xticklabels(), visible=False)
                if (element == 0) : 
                    robs_plot_text(panel, channel_names, \
                                    vpos='bottom', hpos='left', colour=default_cols[0:len(selec_channels)], padder=0.08, xpadder=0.03) #['Wind group', 'SN-II group', 'SN-Ia group']
                if (element == 1) : 
                    #Prep legend:
                    plt.plot(np.nan, np.nan, linestyle=':', color='black', label=orig_GCE_ref)
                    if PLOT_SINGLE_AND_BINARY_STARS == 1 :
                        plt.plot(np.nan, np.nan, linestyle='dashed', color='black', label=bclabel[1])
                    plt.plot(np.nan, np.nan, linestyle='solid', color='black', label=bclabel[0])
                    plt.legend()
                #plt.show() 
                plt.subplots_adjust(hspace=.0)
                
                if PLOT_PIECHARTS == 1 :
                    # for group in range(len(selec_channel_types)) :
                    #     panel = plt.subplot(rowno, colno, (element*(len(selec_channel_types)+1))+group+2)
                    #     #plt.gca().set_prop_cycle(None)
                    #     if PLOT_SINGLE_AND_BINARY_STARS == 1 :
                    #         pieLabels = [orig_GCE_ref, bclabel[0], bclabel[1]]
                    #         segments = np.array([normMass_lgal[group],normMass_bc[group],normMass_bc_2[group]])
                    #     else :
                    #         pieLabels = [orig_GCE_ref, bclabel[0]]
                    #         segments = np.array([normMass_lgal[group],normMass_bc[group]])
                    #     plt.pie(segments, colors=pieCols) #, labels=pieLabels, labeldistance=.6) #, autopct='%1.1f%%'
                    panel = plt.subplot(rowno, colno, (element*colno)+2)
                    rad = 1.0
                    size_outer = 0.15
                    size_inner = 0.7
                    if PLOT_SINGLE_AND_BINARY_STARS == 1 :
                        pieLabels = [orig_GCE_ref, r'binary\_c (single)', bclabel[0]] #[orig_GCE_ref, bclabel[0], bclabel[1]]
                        vals = np.array([normMass_lgal, normMass_bc_2, normMass_bc])
                    else :
                        pieLabels = [orig_GCE_ref, bclabel[0]]
                        vals = np.array([normMass_lgal, normMass_bc])
                    cmap = plt.colormaps["tab20c"]
                    # outer_colors = cmap(np.arange(3)*4)
                    # inner_colors = cmap([1, 2, 3, 5, 6, 7, 9, 10, 11])
                    outer_colors = ['lightgrey','darkgrey','dimgrey']
                    inner_colors = cmap(np.arange(3)*4)
                    lab1 = [s + "\n" for s in pieLabels]
                    #lab2 = [str(f) for f in np.round(np.sum(vals,axis=1),4)]
                    #lab2 = [s + r"$\textnormal{M}_{\odot}$" for s in lab2]
                    lab2 = [str(f) for f in np.round((np.sum(vals,axis=1)/np.sum(vals[0]))*100.,0).astype(int)]
                    lab2 = [s + r"\%" for s in lab2]
                    pcs = np.round((vals.T/np.sum(vals,axis=1))*100.,0) #np.round((vals.T/np.sum(vals,axis=1))*100.,1)
                    pcs_lab = pcs.astype(int).astype(str).T.flatten()
                    pcs_lab = [s + r"\%" for s in pcs_lab]
                    for ii in range(len(pcs_lab)) :
                        if pcs.T.flatten()[ii] < 6. :
                            pcs_lab[ii] = '' # Remove labels for sources which contribute <6%
                    if element == 0 :
                        plt.pie(vals.sum(axis=1), radius=rad, colors=outer_colors, wedgeprops=dict(width=size_outer, edgecolor='w'), \
                                labels=list(map(str.__add__, lab1, lab2))) #, explode=np.full((3),0.1)
                        # patches, texts = plt.pie(vals.flatten(), radius=1-size_outer, colors=inner_colors, wedgeprops=dict(width=size_inner, edgecolor='w'), \
                        #         labels=pcs_lab, labeldistance=.6, textprops=dict(color="w"))
                        # for text in texts:
                        #     text.set_horizontalalignment('center')
                        #     text.set_fontsize(10)
                    else :
                        plt.pie(vals.sum(axis=1), radius=rad, colors=outer_colors, wedgeprops=dict(width=size_outer, edgecolor='w'), \
                                labels=lab2)
                        #plt.pie(vals.flatten(), radius=1-size_outer, colors=inner_colors, wedgeprops=dict(width=size_inner, edgecolor='w'))
                    patches, texts = plt.pie(vals.flatten(), radius=rad-size_outer, colors=inner_colors, wedgeprops=dict(width=size_inner, edgecolor='w'), \
                            labels=pcs_lab, labeldistance=.7, textprops=dict(color="black"))
                    for text in texts:
                        text.set_horizontalalignment('center')
                        text.set_fontsize(10)
                        
                    
            pdf.savefig(fig, bbox_inches='tight') 
            pdf.close()
            
            # Get total metal yields:
            totMetalMass_lgal = totNormMass_lgal - totNormHHeMass_lgal
            #print("Total mass ejected:", totNormMass_lgal, "\nTotal H+He mass ejected:", totNormHHeMass_lgal, "\nTotal metal mass ejected:", totMetalMass_lgal)
            totMetalMass_bc = totNormMass_bc - totNormHHeMass_bc
            if PLOT_SINGLE_AND_BINARY_STARS == 1 :
                totMetalMass_bc_2 = totNormMass_bc_2 - totNormHHeMass_bc_2
            
        loopa = 1 #Is 0 for first loop over metlalicities, and 1 thereafter.
        
if PLOT_UNPROCESSED_COMP == 1 :
    Z_in = Original_metallicities
    # SFR = np.array([0.01,0.1,1.0,10.]) #Msun/yr
    # dt_SFH = 10.0 #Myr
    # M_formed = SFR*dt_SFH*1.e6 #Msun
    M_formed = np.array([0.01,0.1,1.0,10.,100.]) #Msun
    M_out = 0.6*M_formed #Mun
    M_unpro = np.outer(Z_in, M_out)
    
    pdf = PdfPages(basedir+"/plots/"+"unprocessed_comp"+mark+".pdf")
    fig = plt.figure(figsize=(15,15))
    for theZ in range(len(Z_in)) :
        plt.plot(M_formed, M_unpro[theZ]/M_formed, label="Z = "+str(Z_in[theZ]))
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r'M_formed / Msun') #r'$M_{\textnormal formed}$/Msun'
    plt.ylabel(r'M_Z,unpro,out / M_formed')
    #plt.title("Original Z: "+original_comp_Z[zz]+" | binary_c Z: "+ensemble_Zs_string[zz])
    #plt.ylim(ylimits)
    plt.legend()
    pdf.savefig(fig) 
 
    
##########
# WRITE L-GALAXIES HEADER FILE:
##########
if WRITE_HEADER_FILE :
    with open(bc_yieldsdir+'h_metals'+mark+'.h', 'w') as f:
        # f.write('THIS IS A TEST\n' \
        #         'THIS IS ANOTHER TEST'+str(len(elements))+'.')
        f.write("//NOTE: This header file was written by the binary_c yields generator on "+now.strftime("%d/%m/%Y")+':\n')
        f.write( \
                "\n#ifdef DETAILED_METALS_AND_MASS_RETURN\
                \n    #define NUM_METAL_CHANNELS 3 //[SNe-II][SNe-Ia][AGBs]\
                \n#ifdef INDIVIDUAL_ELEMENTS\
                \n#ifndef MAINELEMENTS\
                \n#ifdef BINARYC\
                \n    #define NUM_ELEMENTS "+str(len(elements)))
        for ele in elements :
            if ele == 'C' :
                f.write("\n    #define "+ele+"b_NUM "+str(elements.index(ele)))
            else :
                f.write("\n    #define "+ele+"_NUM "+str(elements.index(ele)))
        f.write( \
                "\n#else //BINARYC\
                \n    #define NUM_ELEMENTS 11\
                \n    #define H_NUM 0\
                \n    #define He_NUM 1\
                \n    #define Cb_NUM 2\
                \n    #define N_NUM 3\
                \n    #define O_NUM 4\
                \n    #define Ne_NUM 5\
                \n    #define Mg_NUM 6\
                \n    #define Si_NUM 7\
                \n    #define S_NUM 8\
                \n    #define Ca_NUM 9\
                \n    #define Fe_NUM 10\
                \n#endif //BINARYC\
                \n#else //MAINELEMENTS\
                \n    #define NUM_ELEMENTS 5\
                \n    #define H_NUM 0\
                \n    #define He_NUM 1\
                \n    #define O_NUM 2\
                \n    #define Mg_NUM 3\
                \n    #define Fe_NUM 4\
                \n#endif //MAINELEMENTS\
                \n#endif //INDIVIDUAL_ELEMENTS\n\
                \n#else //DETAILED_METALS_AND_MASS_RETURN\
                \n    #define NUM_METAL_CHANNELS 1\
                \n#endif //DETAILED_METALS_AND_MASS_RETURN"
            )


##########
# WRITE METADATA:
##########
if LOAD_DATA :
    with open(bc_yieldsdir+bc_yieldsfile+'_'+'info'+mark+'.txt', 'w') as f:
        f.write('INFO FOR '+mark+' FILES:'+'\n')
        f.write('--------------------\n\n')
        f.write('Date/time produced:\n')
        f.write(now.strftime("%d/%m/%Y %H:%M:%S")+'\n')   
        f.write('\nbinary_c files location:\n')
        f.write(ensembledir+IMFdir+'\n\n')
        f.write('binary_c files:\n')
        for zz in range(zz_start,zz_stop) :
            f.write(ensemble_datafile+'-'+ensemble_file_codes[zz]+'.json'+'\n')         
        f.write('\nWritten files location:\n')
        f.write(bc_yieldsdir+'\n\n') 
        f.write('Written files:\n')
        f.write(bc_yieldsfile+'_'+'logt_timesteps'+mark+'.txt\n')
        f.write(bc_yieldsfile+'_'+'metallicities'+mark+'.txt\n')
        if WRITE_HEADER_FILE :
            f.write('h_metals'+'_'+mark+'.h\n')
        for mets in the_metallicities :
            for group in selec_channels :
                f.write(bc_yieldsfile+'_'+group+'_'+'Z'+mets+'_'+'Yields'+mark+'.txt\n')
                f.write(bc_yieldsfile+'_'+group+'_'+'Z'+mets+'_'+'EjectedMasses'+mark+'.txt\n')
                f.write(bc_yieldsfile+'_'+group+'_'+'Z'+mets+'_'+'TotalMetals'+mark+'.txt\n')
                if WRITE_SN_RATES :
                    f.write(bc_yieldsfile+'_'+group+'_'+'Z'+mets+'_'+'Rates'+mark+'.txt\n')   
        f.write('\nSwitches:\n')
        f.write(json.dumps(switches_dict)+'\n\n')
        f.write('elements:\n')
        f.write(json.dumps(elements)+'\n\n')
        f.write('selec_channels:\n')
        f.write(json.dumps(selec_channels)+'\n\n')
        f.write('sn_channels:\n')
        f.write(json.dumps(sn_channels)+'\n\n')
        f.write('Outputted yield units:\n')
        if CONVERT_DATA_TO_LINEAR_TIME_BINS == 1 :
            f.write('d(M/Msun)/d(t/Myr) * 1/Msun'+'\n\n')
        else :
            f.write('d(M/Msun)/dlog(t/Myr) * 1/Msun'+'\n\n')
        f.write("END")
    

################
print("\nDONE!")
