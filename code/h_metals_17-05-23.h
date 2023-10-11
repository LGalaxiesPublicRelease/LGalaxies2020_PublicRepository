#ifdef DETAILED_METALS_AND_MASS_RETURN
/*struct elements_str
{
    float H;
    float He;
    float O;
    float Mg;
    float Fe;
};
#endif  //MAINELEMENTS*/
/* Define two views into the same element structure/array.
   For example define:
      union elements DiskMassElements;
   Then access as either
      DiskMassElements.str.H
   or
      DiskMassElements.arr[0]
*/
/*
union elements
{
    struct elements_str str;
    float arr[NUM_ELEMENTS];
};*/
#define NUM_METAL_CHANNELS 3 //[SNe-II][SNe-Ia][AGBs] //N.B. This ordering is different to that used in Dust-Master, where a metals struct is used. (17-01-22)

#ifdef INDIVIDUAL_ELEMENTS
//Number of chemical elements tracked:
#ifndef MAINELEMENTS
#define NUM_ELEMENTS 11 //All: [H][He][C][N][O][Ne][Mg][Si][S][Ca][Fe] //For BINARYC, NUM_ELEMENTS is now defined by reading in the number of elements from ensemble_output_elements.txt in read_parameters.c. (17-05-23)
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
#ifdef AL26
#define Al26_NUM 9 //When the Al26 switch is on, we assume that the binary_c yields actually contain Al26 in array element 9, rather than Calcium. (09-07-22)
#endif //AL26

#else //MAINELEMENTS
#define NUM_ELEMENTS 5 //Only [H][He][O][Mg][Fe] //For BINARYC, NUM_ELEMENTS is now defined by reading in the number of elements from ensemble_output_elements.txt in read_parameters.c. (17-05-23)
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

//ROB: These dust definitions are now in a new file called h_dust.h (05-11-21):
/*#ifdef DETAILED_DUST 			//Scott 15/11/2015
#ifdef FULL_DUST_RATES
#define NDUST 5 // Needed for HDF5 table creation - should match struct below
struct DustRates
{
  float AGB;
  float SNII;
  float SNIA;
  float GROW;
  float DEST;
};
#endif //FULL_DUST_RATES
#endif //DETAILED_DUST*/
