#ifdef DETAILED_METALS_AND_MASS_RETURN
#define NUM_METAL_CHANNELS 3

#ifdef INDIVIDUAL_ELEMENTS
//Number of chemical elements tracked:
#ifndef MAINELEMENTS
#define NUM_ELEMENTS 11 //All: [H][He][C][N][O][Ne][Mg][Si][S][Ca][Fe]
#else
#define NUM_ELEMENTS 5 //Only [H][He][O][Mg][Fe]
#endif
#endif //INDIVIDUAL_ELEMENTS

#else //DETAILED_METALS_AND_MASS_RETURN
#define NUM_METAL_CHANNELS 1

#endif //DETAILED_METALS_AND_MASS_RETURN
