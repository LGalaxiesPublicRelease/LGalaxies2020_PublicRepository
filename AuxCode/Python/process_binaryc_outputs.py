# -*- coding: utf-8 -*-
"""
Created on Thu Oct 12 09:28:57 2023

@author: R. Yates & D. Hendricks
"""

"""
Simple script to read-in binary_c ensemble outputs, process them into look-up
tables for L-Galaxies.
"""

import json
import pandas as pd
import numpy as np
import os
import shutil
from datetime import datetime
now = datetime.now()

#Local packages:
import sys
sys.path.append('./Robs_python_routines/')
from ensemble_functions import inflate_ensemble_with_lists, find_columnames_recursively


##########
# Switches:
##########
USE_NEW_ENSEMBLES = 0 #If on, new binary_c ensembles stored in /YieldTables/binary_c_yields/new/ensembles/ will be used. If off, the default ensembles in /YieldTables/binary_c_yields/default/ensembles/ will be used.
SINGLE_STARS = 0 #If on, the single star stellar populations from binary_c will be loaded, rather than the binary star ones. (30-05-23)
WRITE_HEADER_FILE = 1 #If on, this script will generate a new h_metals.h file for L-Galaxies, which includes (a) the correct number of chemical elements to consider, and (b) their order.

switches_dict = {
    'USE_NEW_ENSEMBLES' : USE_NEW_ENSEMBLES,
    'SINGLE_STARS' : SINGLE_STARS,
    'WRITE_HEADER_FILE' : WRITE_HEADER_FILE
}


##############
# DIRECTORIES:
##############
basedir = '../../'
if USE_NEW_ENSEMBLES == 1 :
    ensembledir = basedir+'YieldTables/binary_c_yields/new/ensembles/' #Where the binary_c ensembles to be processed will be read from
else : #default ensembles
    ensembledir = basedir+'YieldTables/binary_c_yields/default/ensembles/' #Where the binary_c ensembles to be processed will be read from
outputdir = basedir+'YieldTables/binary_c_yields/new/' #Where the yield tables for L-Galaxies will be written to
if SINGLE_STARS == 1 :
    ensembledir = ensembledir+'singleStars/'
    outputdir = outputdir+'singleStars/'
else : #binary stars      
    ensembledir = ensembledir+'binaryStars/'
    outputdir = outputdir+'binaryStars/'
    
    
##########
# ENSEMBLES:
##########
ensemble_file_metallicities = np.array([0.0001,0.001,0.004,0.008,0.01,0.03])
if SINGLE_STARS == 1 :
    ensemble_file_codes = ['319d1c05','bcc22d43','da1bfc0c','4fb2a217','b3e037d4','1e923d7c']
    bcYields_tag = 'binaryc_origEnsembles_singleStars_kroupa2001' #can lose these tags, and be cleverer about how we call files below.
    bcYields_tag2 = 'binaryc_origEnsembles_binaryStars_kroupa2001'
else : #binary stars      
    ensemble_file_codes = ['6d984eac','dab8bea4','3083cc1f','47cd834d','d723bdf4','a28ea3be']
    bcYields_tag = 'binaryc_origEnsembles_binaryStars_kroupa2001'
    bcYields_tag2 = 'binaryc_origEnsembles_singleStars_kroupa2001'


#################
#INPUT VARIABLES:
#################
elements = ['H','He','C','N','O','Ne','Mg','Si','S','Ca','Fe'] #Put the elements/isotopes you want to track in L-Galaxies here
selec_channels = {
    'AGB'  : ['AGB', 'RLOF', 'Wind', 'TZ', 'GB', 'Comenv'],
    'SNII' : ['SNII', 'SNIIa', 'SNIbc', 'WR', 'SNeCAP', 'SNAIC'],
    'SNIa' : ['SNIa_ELD', 'SNIa_CHAND_Coal', 'SNIa_CHAND', 'Novae']
}
selec_channel_types = list(selec_channels.keys())
channel_names = ['Wind group', 'SN-II group', 'SN-Ia group']
all_channels = np.sort(np.array(['Wind', 'WR', 'Comenv', 'SNIbc', 'TZ', 'AGB', 'SNII', 'GB', 'SNAIC', 'SNeCAP', 'Novae', 'SNIa_CHAND', 'SNIa_ELD', 'SNIIa', 'SNIa_CHAND_Coal', 'RLOF']))
sn_channels = { #As listed in ensemble_data['scalars']
    'AGB'  : ['SN_TZ'], #'TPAGB'
    'SNII' : ['SN_II', 'SN_IIa', 'SN_IBC', 'SN_AIC', 'SN_AIC_BH'],
    'SNIa' : ['SN_IA_CHAND', 'SN_IA_CHAND_Coal', 'SN_IA_ELD']
}


###########
#LOAD DATA:
###########
print("\n########")
print("BINARY_C")
print("########")
the_metallicities = np.array([]) #Array to stor the metallicity files actually loaded (for printing to metadata).         
for zz in range(0,len(ensemble_file_metallicities)) :
    with open(ensembledir+'ensemble_output-'+ensemble_file_codes[zz]+'.json', 'r', encoding="utf-8") as f:
        meta_data = json.loads(f.read())['metadata']
    with open(ensembledir+'ensemble_output-'+ensemble_file_codes[zz]+'.json', 'r', encoding="utf-8") as f:    
        ensemble_data = json.loads(f.read())['ensemble']

    #Get key parameters:
    bc_metallicity = str(meta_data['settings']['population_settings']['bse_options']['metallicity'])
    the_metallicities = np.append(the_metallicities,bc_metallicity)
    bc_total_probability_weighted_mass = meta_data['total_probability_weighted_mass'] #This is M_av,sp, the "arithmetic mean mass of stellar systems in the population" to which the yields are normalised.
    bc_timsteps = np.sort(np.array(list(ensemble_data['Xyield']['time'].keys()), dtype=float))
    bc_dlogt = (bc_timsteps[-1]-bc_timsteps[0])/(len(bc_timsteps)-1) #This is assuming you're reading in an ensemble with log timesteps
    if SINGLE_STARS != 1 :
        cor_factor = 2 #To correct printing error in binary_c ensemble outputs which doubles some factors
    else :
        cor_factor = 1
    print("\n***********")
    print("Z = "+str(bc_metallicity)+":")
    print("***********")
    print("Ensemble file: "+ensemble_file_codes[zz])
    print("log timesteps?: "+str(bool(meta_data['settings']['population_settings']['bse_options']['ensemble_logtimes'])))
    print("log(dt/Myr) = "+str(meta_data['settings']['population_settings']['bse_options']['ensemble_logdt']/cor_factor))
    print("Max evo time (Myr) = "+str(meta_data['settings']['population_settings']['bse_options']['max_evolution_time']/cor_factor))
    print("Multiplicity = "+str(int(meta_data['settings']['population_settings']['bse_options']['multiplicity']/cor_factor))+"  (1 = single stars, 2 = binaries, 3 = triples, 4 = quadruples)")
    print("Binaries = "+str(meta_data['settings']['population_settings']['custom_options']['binaries']))
    print("Mmax (Msun) = "+str(meta_data['settings']['population_settings']['custom_options']['mmax']))

    # Inflate yield data to lists, fetch column names, and load into a pandas dataframe:
    yield_data = ensemble_data['Xyield']
    data_list = inflate_ensemble_with_lists(yield_data)
    data_array = np.array(data_list).T
    columnames = find_columnames_recursively(yield_data)
    columnames = columnames + ['yield_per_solarmass'] 
    df = pd.DataFrame(data_array, columns=columnames)
    df = df.astype({"time": float, "yield_per_solarmass": float}) #Make sure we set the correct type in the columns
    
    sn_data = ensemble_data['scalars']
    sn_rates = {}


##############
#WRITE YIELDS:
##############
    df_blank = df.groupby(by='time')['yield_per_solarmass'].sum()*0.0   
    for group in selec_channel_types :
        source_df = df[df['source'].isin(selec_channels[group])]
        sn_avail = []
        for sn_channel in sn_channels[group] :
            try :
                sn_avail.append([i for i in sn_data.keys() if sn_channel == i][0])
            except :
                y = 0.0

        #(a) element yields for each channel (AGB, SNe-II, SNe_Ia):
        for element in elements : 
            ele = (source_df["isotope"].str.contains(element)) \
                & (((source_df["isotope"].str.replace(element,'')).str.isnumeric()) | (source_df["isotope"].str.replace(element,'') == '')) #Selects entries where the isotope name contains only the given element name (e.g. H, He, etc) plus numbers, or contains only the given isotope name (e.g. Al26, H34, etc).
            element_df = source_df[ele]
            grouped_by_time = element_df.groupby(by='time')['yield_per_solarmass'].sum()           
            the_yield_array = np.nan_to_num(np.array(df_blank+grouped_by_time))
            the_yield_list = list(np.array(the_yield_array, dtype=str))
            if (element == elements[0]) : 
                writemode = 'w'
                tot_ej_masses_array = the_yield_array
            else : 
                writemode = 'a'
                tot_ej_masses_array += the_yield_array
                if (element == elements[2]) :
                    tot_metal_masses_array = the_yield_array              
                elif (element != elements[1]) :
                    tot_metal_masses_array += the_yield_array                        
            with open(outputdir+'ensemble_output_'+group+'_'+'Z'+bc_metallicity+'_'+'Yields.txt', writemode) as f:
                f.write(' '.join(the_yield_list)+' ')
        
        #(b) total ejected masses for each channel (AGB, SNe-II, SNe_Ia):
        the_tot_ej_masses_list = list(np.array(tot_ej_masses_array, dtype=str))
        with open(outputdir+'ensemble_output_'+group+'_'+'Z'+bc_metallicity+'_'+'EjectedMasses.txt', 'w') as f:
              f.write(' '.join(the_tot_ej_masses_list))
        
        #(c) total ejected metal masses for each channel (AGB, SNe-II, SNe_Ia):
        the_tot_metal_masses_list = list(np.array(tot_metal_masses_array, dtype=str))
        with open(outputdir+'ensemble_output_'+group+'_'+'Z'+bc_metallicity+'_'+'TotalMetals.txt', 'w') as f:
              f.write(' '.join(the_tot_metal_masses_list))

        # (d) SN rates for each channel (AGB, SNe-II, SNe_Ia):
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
                        rate = sn_data[sn_source][bc_timestep]
                        rate /= bc_dlogt #correct from a SN number to a SN rate per logt, using the log timestep width from the file
                        sn_rates[bc_timestep] += rate
            else : print('SN type '+sn_source+' is not in ensemble_output-'+ensemble_file_codes[zz])
        with open(outputdir+'ensemble_output_'+group+'_'+'Z'+bc_metallicity+'_'+'Rates.txt', 'w') as f:
            f.write(' '.join(np.array(list(sn_rates.values()), dtype=str)))

    #(e) timesteps:
    if zz == 0 :
        with open(outputdir+'ensemble_output_'+'logt_timesteps.txt', 'w') as f:
            logt_timesteps = ' '.join(list(np.array(df_blank.index, dtype=str)))
            f.write(logt_timesteps)
    #Check that the binary_c timestep structure is the same for all ensembles:
    if zz > 0 :
        logt_timesteps_new = ' '.join(list(np.array(df_blank.index, dtype=str)))
        if logt_timesteps_new != logt_timesteps :
            print('WARNING: Timesteps for ensemble_output-'+ensemble_file_codes[0]+' (Z='+str(the_metallicities[0])+') do not match with ensemble_output-'+ensemble_file_codes[zz]+' (Z='+str(the_metallicities[zz])+')')

    #(f) metallicities:
    with open(outputdir+'ensemble_output_'+'metallicities.txt', 'w') as f:
        f.write(' '.join(list(np.array(ensemble_file_metallicities, dtype=str)))) 
                
                
###################
# WRITE h_metals.h:
###################
if WRITE_HEADER_FILE :
    #Rename old header file:
    os.replace(basedir+'code/'+'h_metals.h', basedir+'code/'+'h_metals_old.h')
    #Write new header file:
    with open(basedir+'code/'+'h_metals.h', 'w') as f:
        f.write("//NOTE: This header file was written by process_binaryc_outputs.py on "+now.strftime("%m/%d/%Y at %H:%M:%S")+'.\n')
        f.write( \
                "\n#ifdef DETAILED_METALS_AND_MASS_RETURN\
                \n    #define NUM_METAL_CHANNELS 3 //[SN-II group][SN-Ia group][Wind group]\
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
    #Make a copy of the new header file in the yields folder:
    shutil.copy2(basedir+'code/'+'h_metals.h', outputdir+'h_metals.h')
    
    print("\nHeader file written.")


##########
# WRITE METADATA:
##########
with open(outputdir+'ensemble_output_'+'info.txt', 'w') as f:
    f.write('INFO FOR BINARY_C YIELD FILES:'+'\n')
    f.write('------------------------------\n\n')
    f.write('Date/time produced:\n')
    f.write(now.strftime("%d/%m/%Y at %H:%M:%S")+'\n')   
    f.write('\nbinary_c files location:\n')
    f.write(ensembledir+'\n\n')
    f.write('binary_c files:\n')
    for zz in range(0,len(ensemble_file_metallicities)) :
        f.write('ensemble_output-'+ensemble_file_codes[zz]+'.json'+'\n')         
    f.write('\nWritten files location:\n')
    f.write(outputdir+'\n\n') 
    f.write('Written files:\n')
    f.write('ensemble_output_'+'logt_timesteps.txt\n')
    f.write('ensemble_output_'+'metallicities.txt\n')
    if WRITE_HEADER_FILE :
        f.write('h_metals.h\n')
    for mets in the_metallicities :
        for group in selec_channels :
            f.write('ensemble_output_'+group+'_'+'Z'+mets+'_'+'Yields.txt\n')
            f.write('ensemble_output_'+group+'_'+'Z'+mets+'_'+'EjectedMasses.txt\n')
            f.write('ensemble_output_'+group+'_'+'Z'+mets+'_'+'TotalMetals.txt\n')
            f.write('ensemble_output_'+group+'_'+'Z'+mets+'_'+'Rates.txt\n')   
    f.write('\nSwitches:\n')
    f.write(json.dumps(switches_dict)+'\n\n')
    f.write('elements:\n')
    f.write(json.dumps(elements)+'\n\n')
    f.write('selec_channels:\n')
    f.write(json.dumps(selec_channels)+'\n\n')
    f.write('sn_channels:\n')
    f.write(json.dumps(sn_channels)+'\n\n')
    f.write('Outputted yield units:\n')
    f.write('d(M/Msun)/dlog(t/Myr) * 1/Msun'+'\n\n')
    f.write("END")

################
print("\nDONE!")
