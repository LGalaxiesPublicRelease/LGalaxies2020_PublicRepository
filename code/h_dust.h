//ROB: 05-11-21: I have now created this dust header file for the NDUST and DustRates definitions. Before, they were in h_metals.h.
//ROB: 27-01-22: As, one by one, all the things defined here were removed/changed in the code, this file is now redundant.
//				 The only remaining parameter, NUM_DUST_RATES, was moved to h_params.h.

#ifdef DETAILED_DUST
//Dust arrays contain 9 array elements, one for each metal element tracked (i.e. all excluding H and He).
//Of these 9, only 6 (C, O, Mg, Si, Ca, Fe) are ever actually updated, as the remaining 3 (N, Ne, S) are considered volatile.
//ROB: THose arrays that had dimension NUM_DUST_ELEMENTS have now been changed to have dimension NUM_ELEMENTS. So this definition is no longer required. (27-01-22)
//#define NUM_DUST_ELEMENTS 9

/*struct elements
{
  float H;
  float He;
#ifndef MAINELEMENTS
  float Cb; //NOTE: Carbon (C) is stored as Cb here
  float N;
#endif
  float O;
#ifndef MAINELEMENTS
  float Ne;
#endif
  float Mg;
#ifndef MAINELEMENTS
  float Si;
  float S;
  float Ca;
#endif
  float Fe;
};
*/

/*#ifdef FULL_DUST_RATES
#define NUM_DUST_RATES 5 // Needed for HDF5 table creation - should match struct below
struct DustRates
{
  float AGB;
  float SNII;
  float SNIA;
  float GROW;
  float DEST;
};
#endif //FULL_DUST_RATES
*/
#endif //DETAILED_DUST
