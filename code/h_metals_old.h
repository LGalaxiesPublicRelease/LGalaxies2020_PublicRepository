//NOTE: This header file was written by process_binaryc_outputs.py on 10/12/2023 at 12:59:15:

#ifdef DETAILED_METALS_AND_MASS_RETURN                
    #define NUM_METAL_CHANNELS 3 //[SN-II group][SN-Ia group][Wind group]                
#ifdef INDIVIDUAL_ELEMENTS                
#ifndef MAINELEMENTS                
#ifdef BINARYC                
    #define NUM_ELEMENTS 11
    #define H_NUM 0
    #define He_NUM 1
    #define Cb_NUM 2
    #define N_NUM 3
    #define O_NUM 4
    #define Ne_NUM 5
    #define Mg_NUM 6
    #define Si_NUM 7
    #define S_NUM 8
    #define Ca_NUM 9
    #define Fe_NUM 10
#else //BINARYC                
    #define NUM_ELEMENTS 11                
    #define H_NUM 0                
    #define He_NUM 1                
    #define Cb_NUM 2                
    #define N_NUM 3                
    #define O_NUM 4                
    #define Ne_NUM 5                
    #define Mg_NUM 6                
    #define Si_NUM 7                
    #define S_NUM 8                
    #define Ca_NUM 9                
    #define Fe_NUM 10                
#endif //BINARYC                
#else //MAINELEMENTS                
    #define NUM_ELEMENTS 5                
    #define H_NUM 0                
    #define He_NUM 1                
    #define O_NUM 2                
    #define Mg_NUM 3                
    #define Fe_NUM 4                
#endif //MAINELEMENTS                
#endif //INDIVIDUAL_ELEMENTS
                
#else //DETAILED_METALS_AND_MASS_RETURN                
    #define NUM_METAL_CHANNELS 1                
#endif //DETAILED_METALS_AND_MASS_RETURN