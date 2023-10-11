/*
 * yield_integrals.c
 *
 * Pre-calculates the normalised ejecta rates at every timestep, assuming 1 Msun populations.
 * Multiply by SFR from SFH bins (and interpolate between default metallicities) to obtain
 * true ejecta rates (done in model_yields.c).
 *
 *  Created on: 10.05.2012
 *      Author: robyates
 *
 * Updates:
 * 17-11-21: Cleaned up to remove the counters, comments, etc, that aren't frequently needed. An uncleaned version called yields_integrals_17-11-21.c has been saved.
 * 04-04-22: Revised/simplified coding to make things clearer. Should return exactly the same result.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "allvars.h"
#include "proto.h"

#ifdef DETAILED_METALS_AND_MASS_RETURN

void init_integrated_yields()
{
  int ii, jj, kk, ll;

  for(ii=0;ii<STEPS*LastDarkMatterSnapShot;ii++)
    for(jj=0;jj<SFH_NBIN;jj++)
      for(kk=0;kk<LIFETIME_Z_NUM;kk++)
	{
	  NormSNIIMassEjecRate[ii][jj][kk]=0.;
	  NormSNIIMetalEjecRate[ii][jj][kk]=0.;
#ifdef INDIVIDUAL_ELEMENTS
	  for(ll=0;ll<NUM_ELEMENTS;ll++)
	    NormSNIIYieldRate[ii][jj][kk][ll]=0.;
#endif
	  NormAGBMassEjecRate[ii][jj][kk]=0.;
	  NormAGBMetalEjecRate[ii][jj][kk]=0.;
#ifdef INDIVIDUAL_ELEMENTS
	  for(ll=0;ll<NUM_ELEMENTS;ll++)
	    NormAGBYieldRate[ii][jj][kk][ll]=0.;
#endif
	  NormSNIaMassEjecRate[ii][jj][kk]=0.;
	  NormSNIaMetalEjecRate[ii][jj][kk]=0.;
#ifdef INDIVIDUAL_ELEMENTS
	  for(ll=0;ll<NUM_ELEMENTS;ll++)
	    NormSNIaYieldRate[ii][jj][kk][ll]=0.;
#endif
	  NormSNIINum[ii][jj][kk]=0.;
	  NormSNIaNum[ii][jj][kk]=0.;
	  NormAGBNum[ii][jj][kk]=0.;
	}
}

void integrate_yields()
{
  double previoustime, newtime, deltaT;
  int snap, step, i, mb, Zi;
  double timet;
  int mbmax;

  int Mi_lower, Mi_upper, t_lower_lifetime, t_upper_lifetime;
  int Mi_lower_lt_SNII, Mi_upper_lt_SNII, Mi_lower_SNII, Mi_upper_SNII, Zi_SNII; //, Zi_correc;
  int Mi_lower_lt_AGB, Mi_upper_lt_AGB, Mi_lower_AGB, Mi_upper_AGB, Zi_AGB;
  double dt, t_lower, t_upper, DTD_lower, DTD_upper, Mi_lower_actual, Mi_upper_actual;
  double Mi_lower_actual_SNII, Mi_upper_actual_SNII, SNIIEjectedMasses_lower_actual, SNIIEjectedMasses_upper_actual, SNIITotalMetals_lower_actual, SNIITotalMetals_upper_actual;
#ifdef INDIVIDUAL_ELEMENTS
  double SNIIYields_lower_actual[NUM_ELEMENTS], SNIIYields_upper_actual[NUM_ELEMENTS];
#endif
  double Mi_lower_actual_AGB, Mi_upper_actual_AGB, AGBEjectedMasses_lower_actual, AGBEjectedMasses_upper_actual, AGBTotalMetals_lower_actual, AGBTotalMetals_upper_actual;
#ifdef INDIVIDUAL_ELEMENTS
  double AGBYields_lower_actual[NUM_ELEMENTS], AGBYields_upper_actual[NUM_ELEMENTS];
#endif
  double Tot_SNII, Tot_SNIa, Tot_SNe, Tot_NormSNIINum; //Tot_SNIIRate, Tot_SNIaRate, Tot_SNe, Tot_SNII_SP, lifetime_lower_actual, lifetime_upper_actual,
  double NormFactor, ML, MU, EjectedMasses_lower, EjectedMasses_upper, TotalMetals_lower, TotalMetals_upper, Yields_lower[NUM_ELEMENTS], Yields_upper[NUM_ELEMENTS];
  double lifetimes_lower, lifetimes_upper;
  //Tot_SNIIRate = 0.0;
  //Tot_SNIaRate = 0.0;
  Tot_SNII = 0.0;
  Tot_SNIa = 0.0;
  Tot_SNe = 0.0;
  //Tot_SNII_SP = 0.0;
  Tot_NormSNIINum = 0.0;
  //lifetime_lower_actual = 0.0;
  //lifetime_upper_actual = 0.0;
  ML = 0.0;
  MU = 0.0;
  /*int lifetime_lower, lifetime_upper;
  lifetime_lower = 0;
  lifetime_upper = 0;*/
  /*double Tot_N_SNII, Tot_N_SNIa, Tot_N_AGB, Tot_O_SNII, Tot_O_SNIa, Tot_O_AGB, Tot_Fe_SNII, Tot_Fe_SNIa, Tot_Fe_AGB;
  Tot_N_SNII = 0.0;
  Tot_N_SNIa = 0.0;
  Tot_N_AGB = 0.0;
  Tot_O_SNII = 0.0;
  Tot_O_SNIa = 0.0;
  Tot_O_AGB = 0.0;
  Tot_Fe_SNII = 0.0;
  Tot_Fe_SNIa = 0.0;
  Tot_Fe_AGB = 0.0;
  double Tot_N_SNII_burst, Tot_N_SNIa_burst, Tot_N_AGB_burst, Tot_O_SNII_burst, Tot_O_SNIa_burst, Tot_O_AGB_burst, Tot_Fe_SNII_burst, Tot_Fe_SNIa_burst, Tot_Fe_AGB_burst;
  double Tot_EjecMass_SNII_burst, Tot_EjecMass_SNIa_burst, Tot_EjecMass_AGB_burst;
  Tot_N_SNII_burst = 0.0;
  Tot_N_SNIa_burst = 0.0;
  Tot_N_AGB_burst = 0.0;
  Tot_O_SNII_burst = 0.0;
  Tot_O_SNIa_burst = 0.0;
  Tot_O_AGB_burst = 0.0;
  Tot_Fe_SNII_burst = 0.0;
  Tot_Fe_SNIa_burst = 0.0;
  Tot_Fe_AGB_burst = 0.0;
  Tot_EjecMass_SNII_burst = 0.0;
  Tot_EjecMass_SNIa_burst = 0.0;
  Tot_EjecMass_AGB_burst = 0.0;*/

  //FRACCOUNTA = 0;

  int L1a,L2a,L3a,L4a,L4aa,L5a,L1b,L2b,L3b,L4b;
  L1a=0;L2a=0;L3a=0;L4a=0;L4aa=0;L5a=0;L1b=0;L2b=0;L3b=0;L4b=0;
  int L1AGB,L2AGB,L3AGB,L4AGB,L5AGB,L6AGB,L7AGB;
  L1AGB=0;L2AGB=0;L3AGB=0;L4AGB=0;L5AGB=0;L6AGB=0;L7AGB=0;

  double First_SFH_bin_width, Tot_NormSNIIMetalEjecRate, Tot_NormSNIIMassEjecRate, Tot_SNII_MetEjecMass; //, Tot_SFH;
  First_SFH_bin_width = 0.0;
  Tot_NormSNIIMetalEjecRate = 0.0;
  Tot_NormSNIIMassEjecRate = 0.0;
  Tot_SNII_MetEjecMass = 0.0;
  int Zi_pick, i_pick;
  Zi_pick = 1; //3 //Choose which of the discrete lifetime metallicities to assume when calculating the SN-II yields: [0.0004, 0.004, 0.008, 0.02, 0.05, 1]
  i_pick = 0; //0 //Choose which SFH bin you want to print out SNII numbers in timestep "step" for. (22-05-20)

  /*int counta;
  TheSFH[0] = 1.0; ///(tau_dt[0]*UnitTime_in_years/Hubble_h);
  for(counta=1;counta<SFH_NBIN;counta++)
    {
      TheSFH[counta] = 0.0; ///(tau_dt[counta]*UnitTime_in_years/Hubble_h);
    }*/

  //Set KALPHA //Integral of the Chabrier IMF (by number) from 0.1 - MAX Msun:
  double KALPHA;
  double F316;
  if (IMF_MAX_MASS == 120.0)
    {
      if (IMF_SLOPE == 2.0)
	{
	  KALPHA = 1.00065;
	  F316 = 0.0417281;
	}
      else if (IMF_SLOPE == 2.15)
	{
	  KALPHA = 1.24409;
	  F316 = 0.0410135;
	}
      else if (IMF_SLOPE == 2.3) //i.e. Normal Chabrier IMF
	{
	  KALPHA = 1.47717;
	  F316 = 0.0385147;
	}
      else if (IMF_SLOPE == 2.6)
	{
	  KALPHA = 1.86960;
	  F316 = 0.0306260;
	}
      else
	{
	  printf("****************\n");
	  printf("In yield_integrals.c:\n");
	  printf("No valid IMF_SLOPE selected. Default x=2.3 (Chabrier IMF) chosen.\n");
	  printf("****************\n");
	  KALPHA = 1.47717; //1.4765;
	  F316 = 0.0385147;
	}
    }
  else if (IMF_MAX_MASS == 100.0)
    {
      KALPHA = 1.49272;
      F316 = 0.0389188;
    }
  else if (IMF_MAX_MASS == 70.0)
    {
      KALPHA = 1.52679;
      F316 = 0.0398185;
    }
  else if (IMF_MAX_MASS == 60.0)
    {
      KALPHA = 1.54319;
      F316 = 0.0402507;
    }
  else if (IMF_MAX_MASS == 50.0)
    {
      KALPHA = 1.56408;
      F316 = 0.0408018;
    }
  else if (IMF_MAX_MASS == 40.0)
    {
      KALPHA = 1.59203;
      F316 = 0.0415416;
    }
  else if (IMF_MAX_MASS == 30.0)
    {
      KALPHA = 1.63252;
      F316 = 0.0426178;
    }
  else if (IMF_MAX_MASS == 25.0)
      {
        KALPHA = 1.66121;
        F316 = 0.0433495;
      }
  else
    {
      KALPHA = 1.49272;
      F316 = 0.0389188;
      printf("****************\n");
      printf("In yield_integrals.c:\n");
      printf("\nIMF_MAX_MASS neither 25, 30, 40, 50, 60, 70, 100, nor 120 Msun\n");
      printf("KALPHA set to 1.49272\n");
      printf("F316 set to 0.0389188\n");
      printf("(These are the values for IMF_MAX_MASS = 100.0)\n");
      printf("****************\n\n");
    }

  //***** LOOP OVER SNAPSHOTS *****
  for(snap=0;snap<LastDarkMatterSnapShot;snap++) {
      previoustime = NumToTime(snap); //Time to z=0 from start of current snapshot [in code units]
      newtime = NumToTime(snap+1); //Time to z=0 from end of current snapshot [in code units]
      deltaT = previoustime - newtime; //Width of current snapshot [in code units]

      //***** LOOP OVER TIMESTEPS *****
      for(step=0;step<STEPS;step++) {
		  dt = deltaT/STEPS;  //Time-width of a timestep in current snapshot [in code units]
		  timet = previoustime - ((step + 0.5) * dt); //Time from middle of the current timestep to z=0 [in code units]
		  First_SFH_bin_width = SFH_dt[snap][step][0]*(UnitTime_in_years/Hubble_h); //Time width of the first (oldest) SFH bin for this timestep [IN YEARS, so it can cancel out with TheSFH which is in years.]

		  //***** LOOP OVER SFH BINS *****
		  for (i=0;i<=SFH_ibin[snap][step];i++) {
			  mbmax = (int)(ceilf(SFH_dt[snap][step][i]/dt)); //New width of SFH bin in number of current timesteps [in code units] //Rounds up to an extra minibin if the SFH bin width isn't exactly divisible by the current timestep width. (17-05-20)
			  if(mbmax < 1) {
				  printf("****************\n");
				  printf("In yield_integrals.c:\n");
				  printf("\nWARNING: SFH bin width is less than the width of one current timestep: mbmax = %i\n", mbmax);
				  printf("****************\n\n");
			  }

			  //***** LOOP OVER MINI BINS *****
			  for (mb=1;mb<=mbmax;mb++) {
				  //From lower/upper edge of mini-bin to middle of current timestep:
				  //N.B. The lowest-z SFH bin actually encompasses the current timestep, so t_lower = -0.5dt and t_upper = +0.5dt for the final minibin of that SFH bin. (26-04-22)
				  //This is accounted for immediately after in find_initial_mass(), which sets Mi_upper to the maximum mass bin in the lifetime arrays when t_lower <= 0.0. (26-04-22)
				  t_lower = (SFH_t[snap][step][i] + (SFH_dt[snap][step][i]) - (mb*(SFH_dt[snap][step][i]/mbmax)) - timet) * UnitTime_in_years/Hubble_h; //IN YEARS //Time from low-z (lower) edge of SFH mini-bin mb to middle of current timestep
				  t_upper = (SFH_t[snap][step][i] + (SFH_dt[snap][step][i]) - (mb*(SFH_dt[snap][step][i]/mbmax)) + (SFH_dt[snap][step][i]/mbmax) - timet) * UnitTime_in_years/Hubble_h; //IN YEARS //Time from high-z (upper) edge of SFH mini-bin mb to middle of current timestep
				  //if (t_lower <= 0.0) printf("%i %i %i %i %i | mbmax = %i | t_lower = %f | t_upper = %f | SFH_t[snap][step][i] = %f | SFH_dt[snap][step][i] = %f | timet = %f\n", snap, step, i, mb, Zi, mbmax, t_lower, t_upper, SFH_t[snap][step][i], SFH_dt[snap][step][i], timet);

				  //***** LOOP OVER INITIAL METALLICITIES *****
				  for (Zi=0;Zi<LIFETIME_Z_NUM;Zi++) {
					  /*
					  Mi_lower = find_initial_mass(t_upper, Zi); //Mass bin (lifetime arrays) corresp. to lowest mass of star to 'die' in current timestep, from SFH bin i.
					  Mi_upper = find_initial_mass(t_lower, Zi); //Mass bin (lifetime arrays) corresp. to highest mass of star to 'die' in current timestep, from SFH bin i.

					  Mi_lower_actual = lifetimeMasses[Mi_lower] + ((lifetimeMasses[Mi_lower+1]-lifetimeMasses[Mi_lower]) * ((t_upper-lifetimes[Zi][Mi_lower])/(lifetimes[Zi][Mi_lower+1]-lifetimes[Zi][Mi_lower]))); //IN MSUN  //Lowest mass of star to 'die' in current timestep from SFH bin i.
					  Mi_upper_actual = lifetimeMasses[Mi_upper] + ((lifetimeMasses[Mi_upper+1]-lifetimeMasses[Mi_upper]) * ((t_lower-lifetimes[Zi][Mi_upper])/(lifetimes[Zi][Mi_upper+1]-lifetimes[Zi][Mi_upper]))); //IN MSUN  //Highest mass of star to 'die' in current timestep from SFH bin i.

					  if (Mi_upper_actual <= 0.0 || Mi_upper_actual > SNII_MAX_MASS || Mi_upper == LIFETIME_MASS_NUM-1 || lifetimeMasses[Mi_upper] >= SNII_MAX_MASS) //No stars of mass above max. SN-II progenitor assumed to contribute chemically.
						  Mi_upper_actual = SNII_MAX_MASS; //Mi_upper_actual could be < 0.0 for the highest Zi, as lifetimes[5][149+1] = 0.0 (i.e. doesn't exist!)

					  if (Mi_lower_actual < AGB_MIN_MASS)
						  Mi_lower_actual = AGB_MIN_MASS; //No stars below AGB_MIN_MASS = 0.85 Msun contribute to chemical enrichment.
					  */
					  //**********
					  //Calculate overall mass limits:
					  //**********
					  Mi_lower = find_initial_mass(t_upper, Zi); //Mass bin (lifetime arrays) corresp. to lowest mass of star to 'die' in current timestep, from SFH bin i.
					  Mi_upper = find_initial_mass(t_lower, Zi); //Mass bin (lifetime arrays) corresp. to highest mass of star to 'die' in current timestep, from SFH bin i.
					  Mi_lower_actual = lifetimeMasses[Mi_lower] + ((lifetimeMasses[Mi_lower+1]-lifetimeMasses[Mi_lower]) * ((t_upper-lifetimes[Zi][Mi_lower])/(lifetimes[Zi][Mi_lower+1]-lifetimes[Zi][Mi_lower]))); //IN MSUN  //Lowest mass of star to 'die' in current timestep from SFH bin i.
					  Mi_upper_actual = lifetimeMasses[Mi_upper] + ((lifetimeMasses[Mi_upper+1]-lifetimeMasses[Mi_upper]) * ((t_lower-lifetimes[Zi][Mi_upper])/(lifetimes[Zi][Mi_upper+1]-lifetimes[Zi][Mi_upper]))); //IN MSUN  //Highest mass of star to 'die' in current timestep from SFH bin i.

					  if (Mi_lower_actual < AGB_MIN_MASS)
						  Mi_lower_actual = AGB_MIN_MASS; //No stars below AGB_MIN_MASS = 0.85 Msun contribute to chemical enrichment.
					  /*This correction isn't needed, because:
					  //1) Setting a lower limit to an even lower value is a bad idea in general, and
					  //2) Mi_lower_actual is only > SNII_MAX_MASS when t_upper < the shortest possible lifetime for chemically-enriching stars,
					  //   and these situations are filtered out below anyway with: if (t_upper >= lifetimes[Zi][LIFETIME_MASS_NUM-1])
					  if (Mi_lower_actual > SNII_MAX_MASS) {
						  printf("%i %i %i %i %i | mbmax = %i | Mi_lower_actual = %f | Mi_upper_actual = %f | Mi_lower = %i | Mi_upper = %i | t_lower = %f | t_upper = %f | lifetimes[Zi][LIFETIME_MASS_NUM-1] = %f\n", snap, step, i, mb, Zi, mbmax, Mi_lower_actual, Mi_upper_actual, Mi_lower, Mi_upper, t_lower, t_upper, lifetimes[Zi][LIFETIME_MASS_NUM-1]); //if (t_upper >= lifetimes[Zi][LIFETIME_MASS_NUM-1])
						  Mi_lower_actual = SNII_MAX_MASS;
					  }*/
					  //N.B. Only the first condition is ever actually met here (in the current version of the code). But the others are possible eventualities which could arise with e.g. different compilers, new bugs added, etc. (26-04-22)
					  if (Mi_upper_actual > SNII_MAX_MASS || Mi_upper_actual <= 0.0 || Mi_upper == LIFETIME_MASS_NUM-1 || lifetimeMasses[Mi_upper] >= SNII_MAX_MASS)
						  Mi_upper_actual = SNII_MAX_MASS;

					  if (mb == 1) { //i.e. if this is the first minibin of this SFH bin
						  Mi_lower_lastTS[i][Zi][0] = Mi_lower_lastTS[i][Zi][1]; //Shift lowest mass threshold from previous timestep to the 0th array element (for use below).
						  Mi_lower_lastTS[i][Zi][1] = Mi_lower_actual; //Write the current timestep's lowest mass threshold to the 1st array element (ready for the subsequent timestep).
					  }
					  //TURN ON FOR DEFAULT MODE (overwrites Mi_lower_lastTS[i][Zi][0] = Mi_lower_lastTS[i][Zi][1] above). TURN OFF WHEN WANTING TO CALCULATE THE SN-II RATE FROM A 1Msun SP BORN IN THE FIRST TIMESTEP (maintaining effective 1-timestep resolution by preventing any mass range over/under-shooting when doing the integrations). THIS ONLY WORKS FOR THE 0th SFH BIN (or at full SFH resolution, probably) (20-05-20):
					  Mi_lower_lastTS[i][Zi][0] = Mi_upper_actual; //DEFAULT MODE: Uncomment this to go back to the default mode, with mass over/underlaps.

					  if (t_upper >= lifetimes[Zi][LIFETIME_MASS_NUM-1]) { //If the longest time from SFH bin i to current timestep is shorter than the shortest possible lifetime, there's no enrichment, so skip calculations.
						  //*****************************************
						  //SNe-II (Disc and Bulge):
						  //*****************************************
						  Zi_SNII = find_initial_metallicity_comp(Zi, i, 2); //Metallicity bin (SNe-II arrays) corresp. to metallicity Zi. The find_initial_metallicity_comp() function automatically accounts for the mismatch in number of metallicities in the lifetime and SNII metallicity arrays (i.e. accounts for the difference between Zi and Zi_correc made below for other variables) (06-04-22)

						  /*//Here, Zi_lt_SNII (formerly Zi_correc) is used to force ejected material at metallicities above the max Z of the SN-II tables to equal the ejected material AT the upper SN-II Z limit. i.e. We choose to ignore the Z=1.0 values from the Portinari+98 input tables. (09-06-16)
						  //For Portinari+98 SN-II yields, Zi_lt_SNII happens to be identical to Zi_SNII
						  if (Zi==5) Zi_lt_SNII = 4; //Metallicity bin (lifetime arrays) corresp. to metallicity required for SN-II integrations.
						  else Zi_lt_SNII = Zi;
						  //printf("Zi = %i | Zi_SNII = %i | Zi_lt_SNII = %i\n", Zi, Zi_SNII, Zi_lt_SNII);*/

						  //Check if mass range is within range for SN-II progenitor stars:
						  Mi_lower_lt_SNII = max_Mi_lower(Mi_lower,2); //Mass bin (lifetime arrays) corresp. to lowest mass of star to 'die' in current timestep from SFH bin i.
						  Mi_upper_lt_SNII = min_Mi_upper(Mi_upper,2); //Mass bin (lifetime arrays) corresp. to highest mass of star to 'die' in current timestep from SFH bin i.
						  if (Mi_lower_lt_SNII <= Mi_upper_lt_SNII) {
							  Mi_lower_SNII = find_SNII_mass_bin(lifetimeMasses[Mi_lower_lt_SNII]); //Mass bin (SNe-II arrays) corresp. to lowest mass of star to 'die' in current timestep from SFH bin i.
							  Mi_upper_SNII = find_SNII_mass_bin(lifetimeMasses[Mi_upper_lt_SNII]); //Mass bin (SNe-II arrays) corresp. to highest mass of star to 'die' in current timestep from SFH bin i.


#ifdef INDIVIDUAL_ELEMENTS
							  //Find true yields at the true upper and lower masses, given by 'Mi_upper_actual' and 'Mi_lower_actual':
							  find_actual_ejecta_limits(2, Mi_lower_actual, Mi_upper_actual, Mi_lower_SNII, Mi_upper_SNII, Zi_SNII,
														&SNIIEjectedMasses_lower_actual, &SNIIEjectedMasses_upper_actual, &SNIITotalMetals_lower_actual, &SNIITotalMetals_upper_actual,
														SNIIYields_lower_actual, SNIIYields_upper_actual);
#else
							  find_actual_ejecta_limits(2, Mi_lower_actual, Mi_upper_actual, Mi_lower_SNII, Mi_upper_SNII, Zi_SNII,
														&SNIIEjectedMasses_lower_actual, &SNIIEjectedMasses_upper_actual, &SNIITotalMetals_lower_actual, &SNIITotalMetals_upper_actual);
#endif

							  //NUMERICALLY INTEGRATE OVER THE MASS RANGE APPLICABLE FOR SNe-II:
							  int j, kk;
							  for (j=Mi_lower_SNII;j<=Mi_upper_SNII;j++) {
								  //**********
								  //Calculate mass limits over which to numerically integrate for SNe-II:
								  //**********
								  Mi_upper_actual_SNII = Mi_upper_actual; //Cannot be > SNII_MAX_MASS, due to checks above.
								  if (Mi_lower_actual < SNII_MIN_MASS)
									  Mi_lower_actual_SNII = SNII_MIN_MASS;
								  else
									  Mi_lower_actual_SNII = Mi_lower_actual;

#ifdef INDIVIDUAL_ELEMENTS
								  //double Masses[], double EjectedMasses[AGB_Z_NUM][AGB_MASS_NUM], double TotalMetals[AGB_Z_NUM][AGB_MASS_NUM], double Yields[AGB_Z_NUM][NUM_ELEMENTS][AGB_MASS_NUM],
								  /*calc_mass_limits_SNII(j, Zi_SNII, Mi_lower_SNII, Mi_upper_SNII, Mi_lower_actual_SNII, Mi_upper_actual_SNII,
								  					SNIIEjectedMasses_lower_actual, SNIIEjectedMasses_upper_actual, SNIITotalMetals_lower_actual, SNIITotalMetals_upper_actual, SNIIYields_lower_actual, SNIIYields_upper_actual,
								  					&ML, &MU, &EjectedMasses_lower, &EjectedMasses_upper, &TotalMetals_lower, &TotalMetals_upper, &Yields_lower, &Yields_upper);*/
								  calc_mass_limits(j, Zi_SNII, Mi_lower_SNII, Mi_upper_SNII, Mi_lower_actual_SNII, Mi_upper_actual_SNII,
												   SNII_MASS_NUM, SNII_Z_NUM, SNIIMasses, SNIIEjectedMasses, SNIITotalMetals, SNIIYields,
										  	  	   SNIIEjectedMasses_lower_actual, SNIIEjectedMasses_upper_actual, SNIITotalMetals_lower_actual, SNIITotalMetals_upper_actual, SNIIYields_lower_actual, SNIIYields_upper_actual,
												   &ML, &MU, &EjectedMasses_lower, &EjectedMasses_upper, &TotalMetals_lower, &TotalMetals_upper, &Yields_lower, &Yields_upper);
#else
								  calc_mass_limits(j, Zi_SNII, Mi_lower_SNII, Mi_upper_SNII, Mi_lower_actual_SNII, Mi_upper_actual_SNII,
										  	  	   SNII_MASS_NUM, SNII_Z_NUM, SNIIMasses, SNIIEjectedMasses, SNIITotalMetals,
												   SNIIEjectedMasses_lower_actual, SNIIEjectedMasses_upper_actual, SNIITotalMetals_lower_actual, SNIITotalMetals_upper_actual,
												   &ML, &MU, &EjectedMasses_lower, &EjectedMasses_upper, &TotalMetals_lower, &TotalMetals_upper);
#endif //INDIVIDUAL_ELEMENTS
								  if (mb == mbmax && MU == Mi_upper_actual && MU < Mi_lower_lastTS[i][Zi][0])
									  MU = Mi_lower_lastTS[i][Zi][0]; //only for last mass bin integrated over
								  if (MU > Mi_lower_lastTS[i][Zi][0] && Mi_lower_lastTS[i][Zi][0] > 0.0)
									  MU = Mi_lower_lastTS[i][Zi][0];
								  if (MU < ML)
									  ML = MU;
								  //**********
								  //Calculate rates:
								  //**********
								  if (SNIIMasses[j] <= SNIA_MAX_MASS)
									  NormFactor = 1.0-A_FACTOR; //For mass range where both SN-II and SN-Ia progenitors are possible
								  else
									  NormFactor = 1.0; //For mass range where only SN-II progenitors are possible
								  //SNII_Rate[(STEPS*snap)+step][Zi] += NormFactor * (MU-ML) * ((Chabrier_IMF(ML)*TheSFH[i]/First_SFH_bin_width) + (Chabrier_IMF(MU)*TheSFH[i]/First_SFH_bin_width))/2.0;
								  SNII_Rate[(STEPS*snap)+step][Zi] += NormFactor * (MU-ML) * ((Chabrier_IMF(ML) + Chabrier_IMF(MU))*TheSFH[i]/First_SFH_bin_width)/2.0;
								  NormSNIINum[(STEPS*snap)+step][i][Zi] += NormFactor * (MU-ML) * (Chabrier_IMF(ML) + Chabrier_IMF(MU))/2.0; //Number of SNe-II exploding in this timestep from minibin mb of SFHBin i [units: # / Msun]
								  NormSNIIMassEjecRate[(STEPS*snap)+step][i][Zi] += NormFactor * (MU-ML) * ((EjectedMasses_lower + EjectedMasses_upper)/2.0); //This is a unitless quantity, which must be multiplied by SFR (i.e. M_starsFormedInSFHBin*dt_SFH, in model_yields.c) to get units of [Msun/yr], and then by dt (i.e. timestep width, also in model_yields.c) to get units of [Msun]. (22-05-20) //NO! --> //IN [MSun/yr] //The rate of ejection of mass from a 1 Msun stellar population.
								  NormSNIIMetalEjecRate[(STEPS*snap)+step][i][Zi] += NormFactor * (MU-ML) * ((TotalMetals_lower + TotalMetals_upper)/2.0); //This is a unitless quantity, which must be multiplied by SFR (in model_yields.c) to get units of [Msun/yr], and then by dt (also in model_yields.c) to get units of [Msun]. (22-05-20) //NO! --> //IN [MSun/yr] //The rate of ejection of 'newly synthesised' metals from a 1 Msun stellar population.
#ifdef INDIVIDUAL_ELEMENTS
								  for (kk=0;kk<NUM_ELEMENTS;kk++) {
#ifndef MAINELEMENTS
									  NormSNIIYieldRate[(STEPS*snap)+step][i][Zi][kk] += NormFactor * (MU-ML) * ((Yields_lower[kk] + Yields_upper[kk])/2.0); //IN [Msun/yr] //The rate of ejection of individual elements from a 1 Msun stellar population.
#else
									  switch(kk){case 0: kk=0; break; case 1: kk=1; break; case 2: kk=4; break; case 3: kk=6; break; case 4: kk=10; break;}
									  NormSNIIYieldRate[(STEPS*snap)+step][i][Zi][kk] += NormFactor * (MU-ML) * ((Yields_lower[kk] + Yields_upper[kk])/2.0); //IN [Msun/yr] //The rate of ejection of individual elements from a 1 Msun stellar population.
#endif
								  }
#endif //INDIVIDUAL_ELEMENTS
							  } //for (j=Mi_lower_SNII;j<=Mi_upper_SNII;j++)
						  } //if (Mi_lower_lt_SNII <= Mi_upper_lt_SNII)


						  //*****************************************
						  //SNe-Ia (Disc and Bulge):
						  //*****************************************
#ifdef DTD
						  t_lower_lifetime = Mi_upper; //Lifetime bin (lifetime arrays) corresp. to longest lifetime (lowest mass) of star to 'die' in current timestep, from SFH bin i.
						  t_upper_lifetime = Mi_lower; //Lifetime bin (lifetime arrays) corresp. to shortest lifetime (highest mass) of star to 'die' in current timestep, from SFH bin i.

						  if ((lifetimes[Zi][Mi_lower] > SNIA_MIN_TIME) && (lifetimes[Zi][Mi_upper] < SNIA_MAX_TIME) && (t_upper > SNIA_MIN_TIME) && (t_lower < SNIA_MAX_TIME)) {
							  int j;
							  for (j=t_lower_lifetime; j>=t_upper_lifetime; j--) { //This is the corrected loop definition, without the -1. (13-06-16)
								  //**********
								  //Calculate time limits over which to numerically integrate:
								  //**********
								  int stat = 0; //Just an integer to store which of the if/else statements in calc_time_limits_SNIa() ends up being entered, for analysis later. Not needed for the functioning of the code. (07-04-22)
								  calc_time_limits_SNIa(j, Zi, t_lower_lifetime, t_upper_lifetime, t_lower, t_upper,
										  	  	  	    &lifetimes_lower, &lifetimes_upper, &stat);
								  DTD_lower = DTDcalc(lifetimes_lower) * 1.0e-9; //NOTE: DTDcalc returns SNe/Gyr, not SNe/yr, hence the multiple (1.0e-9).
								  DTD_upper = DTDcalc(lifetimes_upper) * 1.0e-9;
								  //**********
								  //Calculate rates:
								  //**********
								  SNIa_Rate[(STEPS*snap)+step][Zi] += A_FACTOR * F316 * KALPHA * (lifetimes_upper-lifetimes_lower) * ((DTD_lower + DTD_upper) * TheSFH[i]/First_SFH_bin_width)/2.0;
								  NormSNIaNum[(STEPS*snap)+step][i][Zi] += A_FACTOR * F316 * KALPHA * (lifetimes_upper-lifetimes_lower) * (DTD_lower + DTD_upper)/2.0; //Number of SNe-Ia exploding in this timestep from minibin mb of SFHBin i [units: # / Msun]
								  NormSNIaMassEjecRate[(STEPS*snap)+step][i][Zi] += A_FACTOR * F316 * KALPHA * (lifetimes_upper-lifetimes_lower) * ((DTD_lower + DTD_upper)*SNIAEJECMASS)/2.0; //IN [Msun/yr]
								  NormSNIaMetalEjecRate[(STEPS*snap)+step][i][Zi] = NormSNIaMassEjecRate[(STEPS*snap)+step][i][Zi];
								  //if (snap == 2 && step < 10) printf("%i %i %i %i %i %i | %i | NormSNIaMassEjecRate[%i][%i][%i] = %f\n", snap, step, i, mb, Zi, j, stat, (STEPS*snap)+step, i, Zi, NormSNIaMassEjecRate[(STEPS*snap)+step][i][Zi]);
#ifdef INDIVIDUAL_ELEMENTS
								  int k;
								  for (k=0;k<NUM_ELEMENTS;k++) {
									  int kk; //Iterator for the SNIa input yield arrays [0 to 41]
#ifndef MAINELEMENTS
									  switch(k){case 0: kk=0; break; case 1: kk=1; break; case 2: kk=5; break; case 3: kk=6; break; case 4: kk=7; break; case 5: kk=9; break; case 6: kk=11; break; case 7: kk=13; break; case 8: kk=15; break; case 9: kk=19; break; case 10: kk=25; break;}
									  NormSNIaYieldRate[(STEPS*snap)+step][i][Zi][k] += A_FACTOR * F316 * KALPHA * (lifetimes_upper-lifetimes_lower) * ((DTD_lower + DTD_upper) * SNIaYields[kk])/2.0; //IN [Msun/yr]
#else //MAINELEMENTS
									  switch(k){case 0: kk=0; break; case 1: kk=1; break; case 2: kk=7; break; case 3: kk=11; break; case 4: kk=25; break;}
									  NormSNIaYieldRate[(STEPS*snap)+step][i][Zi][k] += A_FACTOR * F316 * KALPHA * (lifetimes_upper-lifetimes_lower) * ((DTD_lower + DTD_upper) * SNIaYields[kk])/2.0; //IN [Msun/yr]
#endif //MAINELEMENTS
								  }
#endif //INDIVIDUAL_ELEMENTS
								  //if (NormSNIaMassEjecRate[(STEPS*snap)+step][i][Zi] < 0.0)
									  //printf("%i %i %i %i %i %i | %i | NormSNIaMassEjecRate[%i][%i][%i] = %.7e | lifetimes_lower = %.7e | lifetimes_upper = %.7e\n", snap, step, i, mb, Zi, j, stat, (STEPS*snap)+step, i, Zi, NormSNIaMassEjecRate[(STEPS*snap)+step][i][Zi], lifetimes_lower, lifetimes_upper);
							  } //for (j=t_lower_lifetime; j>=t_upper_lifetime; j--)
						  } //if ((lifetimes[Zi][Mi_lower] > SNIA_MIN_TIME) && (lifetimes[Zi][Mi_upper] < SNIA_MAX_TIME) && (t_upper > SNIA_MIN_TIME) && (t_lower < SNIA_MAX_TIME))
#endif //DTD


						  //*****************************************
						  //AGB Winds (Disc and Bulge):
						  //*****************************************
						  Zi_AGB = find_initial_metallicity_comp(Zi, i, 4);

						  /*//Here, Zi_correc is used to force ejected material at metallicities below/above the Z limits of the AGB tables to equal the ejected material AT the lower/upper AGB Z limits.
						  if (Zi==0) Zi_correc = 1; if (Zi==1) Zi_correc = 1; if (Zi==2) Zi_correc = 2; if (Zi==3) Zi_correc = 3; if (Zi==4) Zi_correc = 3; if (Zi==5) Zi_correc = 3;
						  //if (Zi_AGB != Zi_correc) printf("Zi_AGB = %i | Zi_correc = %i | Zi = %i\n", Zi_AGB, Zi_correc, Zi);*/

						  /*Mi_lower = find_initial_mass(t_upper, Zi_correc); //Mass bin (lifetime arrays) corresp. to lowest mass of star to 'die' in current timestep, from SFH bin i.
						  Mi_upper = find_initial_mass(t_lower, Zi_correc); //Mass bin (lifetime arrays) corresp. to highest mass of star to 'die' in current timestep, from SFH bin i.
						  Mi_lower_actual = lifetimeMasses[Mi_lower] + ((lifetimeMasses[Mi_lower+1]-lifetimeMasses[Mi_lower]) * ((t_upper-lifetimes[Zi_correc][Mi_lower])/(lifetimes[Zi_correc][Mi_lower+1]-lifetimes[Zi_correc][Mi_lower]))); //IN MSUN  //Lowest mass of star to 'die' in current timestep from SFH bin i.
						  Mi_upper_actual = lifetimeMasses[Mi_upper] + ((lifetimeMasses[Mi_upper+1]-lifetimeMasses[Mi_upper]) * ((t_lower-lifetimes[Zi_correc][Mi_upper])/(lifetimes[Zi_correc][Mi_upper+1]-lifetimes[Zi_correc][Mi_upper]))); //IN MSUN  //Highest mass of star to 'die' in current timestep from SFH bin i.
						  if (Mi_upper_actual <= 0.0 || Mi_upper_actual > AGB_MAX_MASS || Mi_upper == LIFETIME_MASS_NUM-1) Mi_upper_actual = AGB_MAX_MASS; //ROB (18-06-21)*/

						  //Check if mass range is within range for AGB stars:
						  Mi_lower_lt_AGB = max_Mi_lower(Mi_lower,4); //Mass bin (lifetime arrays) corresp. to lowest mass of star to 'die' in current timestep from SFH bin i.
						  Mi_upper_lt_AGB = min_Mi_upper(Mi_upper,4); //Mass bin (lifetime arrays) corresp. to highest mass of star to 'die' in current timestep from SFH bin i.
						  //if (Mi_lower_actual > AGB_MAX_MASS) printf("%i %i %i %i %i | Mi_lower=%i | Mi_upper=%i | Mi_lower_lt_AGB=%i | Mi_upper_lt_AGB=%i | Mi_lower_actual=%f | Mi_upper_actual=%f\n", snap, step, i, mb, Zi, Mi_lower, Mi_upper, Mi_lower_lt_AGB, Mi_upper_lt_AGB, Mi_lower_actual, Mi_upper_actual); //if (snap == 2 && step < 10)
						  if (Mi_lower_lt_AGB <= Mi_upper_lt_AGB && Mi_lower_actual <= AGB_MAX_MASS) {
							  Mi_lower_AGB = find_agb_mass_bin(lifetimeMasses[Mi_lower_lt_AGB]); //Mass bin (AGB arrays) corresp. to lowest mass of star to 'die' in current timestep from SFH bin i.
							  Mi_upper_AGB = find_agb_mass_bin(lifetimeMasses[Mi_upper_lt_AGB]); //Mass bin (AGB arrays) corresp. to highest mass of star to 'die' in current timestep from SFH bin i.

#ifdef INDIVIDUAL_ELEMENTS
							  //Find true yields at the true upper and lower masses, given by 'Mi_upper_actual' and 'Mi_lower_actual':
							  find_actual_ejecta_limits(4, Mi_lower_actual, Mi_upper_actual, Mi_lower_AGB, Mi_upper_AGB, Zi_AGB,
										&AGBEjectedMasses_lower_actual, &AGBEjectedMasses_upper_actual, &AGBTotalMetals_lower_actual, &AGBTotalMetals_upper_actual,
										AGBYields_lower_actual, AGBYields_upper_actual);
#else
							  find_actual_ejecta_limits(4, Mi_lower_actual, Mi_upper_actual, Mi_lower_AGB, Mi_upper_AGB, Zi_AGB,
										&AGBEjectedMasses_lower_actual, &AGBEjectedMasses_upper_actual, &AGBTotalMetals_lower_actual, &AGBTotalMetals_upper_actual);
#endif
							  //if (Mi_upper_actual >= 7.0) printf("%i %i %i %i %i | Mi_lower_actual = %f | Mi_upper_actual = %f | AGBEjectedMasses_lower_actual = %f | AGBEjectedMasses_upper_actual = %f\n", snap, step, i, mb, Zi, Mi_lower_actual, Mi_upper_actual, AGBEjectedMasses_lower_actual, AGBEjectedMasses_upper_actual);

							  //NUMERICALLY INTEGRATE OVER THE MASS RANGE APPLICABLE FOR SNe-II:
							  int j, kk;
							  for (j=Mi_lower_AGB;j<=Mi_upper_AGB;j++) {
								  //**********
								  //Calculate mass limits over which to numerically integrate for SNe-II:
								  //**********
								  Mi_lower_actual_AGB = Mi_lower_actual; //Cannot be < AGB_MIN_MASS, due to checks above.
								  if (Mi_upper_actual > AGB_MAX_MASS)
									  Mi_upper_actual_AGB = AGB_MAX_MASS;
								  else
									  Mi_upper_actual_AGB = Mi_upper_actual;
#ifdef INDIVIDUAL_ELEMENTS
								  /*calc_mass_limits_AGB(j, Zi_AGB, Mi_lower_AGB, Mi_upper_AGB, Mi_lower_actual_AGB, Mi_upper_actual_AGB,
								  						AGBEjectedMasses_lower_actual, AGBEjectedMasses_upper_actual, AGBTotalMetals_lower_actual, AGBTotalMetals_upper_actual, AGBYields_lower_actual, AGBYields_upper_actual,
								  						&ML, &MU, &EjectedMasses_lower, &EjectedMasses_upper, &TotalMetals_lower, &TotalMetals_upper, &Yields_lower, &Yields_upper);*/
								  calc_mass_limits(j, Zi_AGB, Mi_lower_AGB, Mi_upper_AGB, Mi_lower_actual_AGB, Mi_upper_actual_AGB,
										  	  	   AGB_MASS_NUM, AGB_Z_NUM, AGBMasses, AGBEjectedMasses, AGBTotalMetals, AGBYields,
												   AGBEjectedMasses_lower_actual, AGBEjectedMasses_upper_actual, AGBTotalMetals_lower_actual, AGBTotalMetals_upper_actual, AGBYields_lower_actual, AGBYields_upper_actual,
												   &ML, &MU, &EjectedMasses_lower, &EjectedMasses_upper, &TotalMetals_lower, &TotalMetals_upper, &Yields_lower, &Yields_upper);
#else
								  calc_mass_limits(j, Zi_AGB, Mi_lower_AGB, Mi_upper_AGB, Mi_lower_actual_AGB, Mi_upper_actual_AGB,
										  	  	   AGB_MASS_NUM, AGB_Z_NUM, AGBMasses, AGBEjectedMasses, AGBTotalMetals,
												   AGBEjectedMasses_lower_actual, AGBEjectedMasses_upper_actual, AGBTotalMetals_lower_actual, AGBTotalMetals_upper_actual,
												   &ML, &MU, &EjectedMasses_lower, &EjectedMasses_upper, &TotalMetals_lower, &TotalMetals_upper);
#endif //INDIVIDUAL_ELEMENTS
								  /*if (mb == mbmax && MU == Mi_upper_actual && MU < Mi_lower_lastTS[i][Zi][0])
									  MU = Mi_lower_lastTS[i][Zi][0]; //only for last mass bin integrated over
								  if (MU > Mi_lower_lastTS[i][Zi][0] && Mi_lower_lastTS[i][Zi][0] > 0.0)
									  MU = Mi_lower_lastTS[i][Zi][0];
								  if (MU < ML)
									  ML = MU;*/
								  //**********
								  //Calculate rates:
								  //**********
								  AGB_Rate[(STEPS*snap)+step][Zi] += (MU-ML) * ((Chabrier_IMF(ML) + Chabrier_IMF(MU))*TheSFH[i]/First_SFH_bin_width)/2.0;
								  NormAGBNum[(STEPS*snap)+step][i][Zi] += (MU-ML) * (Chabrier_IMF(ML) + Chabrier_IMF(MU))/2.0; //Number of SNe-II exploding in this timestep from minibin mb of SFHBin i [units: # / Msun]
								  NormAGBMassEjecRate[(STEPS*snap)+step][i][Zi] += (MU-ML) * ((EjectedMasses_lower + EjectedMasses_upper)/2.0); //This is a unitless quantity, which must be multiplied by SFR (i.e. M_starsFormedInSFHBin*dt_SFH, in model_yields.c) to get units of [Msun/yr], and then by dt (i.e. timestep width, also in model_yields.c) to get units of [Msun]. (22-05-20) //NO! --> //IN [MSun/yr] //The rate of ejection of mass from a 1 Msun stellar population.
								  NormAGBMetalEjecRate[(STEPS*snap)+step][i][Zi] += (MU-ML) * ((TotalMetals_lower + TotalMetals_upper)/2.0); //This is a unitless quantity, which must be multiplied by SFR (in model_yields.c) to get units of [Msun/yr], and then by dt (also in model_yields.c) to get units of [Msun]. (22-05-20) //NO! --> //IN [MSun/yr] //The rate of ejection of 'newly synthesised' metals from a 1 Msun stellar population.
#ifdef INDIVIDUAL_ELEMENTS
								  for (kk=0;kk<NUM_ELEMENTS;kk++) {
#ifndef MAINELEMENTS
									  NormAGBYieldRate[(STEPS*snap)+step][i][Zi][kk] += (MU-ML) * ((Yields_lower[kk] + Yields_upper[kk])/2.0); //IN [Msun/yr] //The rate of ejection of individual elements from a 1 Msun stellar population.
#else
									  switch(kk){case 0: kk=0; break; case 1: kk=1; break; case 2: kk=4; break; case 3: kk=6; break; case 4: kk=10; break;}
									  NormAGBYieldRate[(STEPS*snap)+step][i][Zi][kk] += (MU-ML) * ((Yields_lower[kk] + Yields_upper[kk])/2.0); //IN [Msun/yr] //The rate of ejection of individual elements from a 1 Msun stellar population.
#endif
								  }
#endif //INDIVIDUAL_ELEMENTS
							  } //for (j=Mi_lower_AGB;j<=Mi_upper_AGB;j++)
						  } //if (Mi_lower_AGB <= Mi_upper_AGB)
						  //*****************************************
					  } //if (t_upper >= lifetimes[Zi][Mi_lower+1])
				  } //for (Zi=0;Zi<LIFETIME_Z_NUM;Zi++)
			  } //for (mb=1;mb<=mbmax;mb++) //MINI_BINS
		  } //for (i=0;i<=SFH_ibin_structure[(SFH_NBIN*snap)+step];i++)
		  Tot_SNII += SNII_Rate[(STEPS*snap)+step][Zi_pick]*(dt*(UnitTime_in_years/Hubble_h));
		  Tot_SNIa += SNIa_Rate[(STEPS*snap)+step][Zi_pick]*(dt*(UnitTime_in_years/Hubble_h));
		  Tot_NormSNIINum += NormSNIINum[(STEPS*snap)+step][i_pick][Zi_pick]; //Summing up all the SNe that exploded in every timestep from all the minibins of SFHBin i_pick at metallicity Zi_pick. //*(dt*(UnitTime_in_years/Hubble_h))/First_SFH_bin_width;
		  Tot_SNe = Tot_SNII+Tot_SNIa;
		  Tot_NormSNIIMetalEjecRate += NormSNIIMetalEjecRate[(STEPS*snap)+step][0][Zi_pick];
		  Tot_NormSNIIMassEjecRate += NormSNIIMassEjecRate[(STEPS*snap)+step][0][Zi_pick];
		  Tot_SNII_MetEjecMass += ((dt*(UnitTime_in_years/Hubble_h))/First_SFH_bin_width) * (NormSNIIMetalEjecRate[(STEPS*snap)+step][0][0] + lifetimeMetallicities[Zi_pick]*NormSNIIMetalEjecRate[(STEPS*snap)+step][0][Zi_pick]);
      } //for(step=0;step<STEPS;step++)
   } //for(snap=0;snap<LastDarkMatterSnapShot;snap++)
#ifdef PARALLEL
  	if(ThisTask == 0)
#endif
    printf("Yield integrals calculated.\n");

  int ii, jj, kk;
  for(snap=0;snap<LastDarkMatterSnapShot;snap++) //LOOP OVER SNAPSHOTS
    for(step=0;step<STEPS;step++) //LOOP O
      for(jj=0;jj<SFH_NBIN;jj++)
	for(kk=0;kk<LIFETIME_Z_NUM;kk++)
	  {

	    ii=(STEPS*snap)+step;

	    if(NormAGBMassEjecRate[ii][jj][kk]>0.3)
	      {
		printf("ii=%d snap=%d step=%d jj=%d kk=%d AGB_Rate=%0.2f\n",ii,snap, step,jj,kk, NormAGBMassEjecRate[ii][jj][kk]);
		terminate("AGB rate too high");
	      }
	    if(NormSNIaMassEjecRate[ii][jj][kk]>0.7)
	      {
		printf("ii=%d jj=%d kk=%d SNIa_Rate=%0.2f\n",ii,jj,kk, NormSNIaMassEjecRate[ii][jj][kk]);
		terminate("SNIa rate too high");
	      }

	    if(NormSNIIMassEjecRate[ii][jj][kk]<0.)
	      NormSNIIMassEjecRate[ii][jj][kk]=0.;
	    if(NormSNIaMassEjecRate[ii][jj][kk]<0.)
	      NormSNIaMassEjecRate[ii][jj][kk]=0.;
	    if(NormAGBMassEjecRate[ii][jj][kk]<0.)
	      NormAGBMassEjecRate[ii][jj][kk]=0.;
	    if(NormSNIIMetalEjecRate[ii][jj][kk]<0.)
	      NormSNIIMetalEjecRate[ii][jj][kk]=0.;
	    if(NormSNIaMetalEjecRate[ii][jj][kk]<0.)
	      NormSNIaMetalEjecRate[ii][jj][kk]=0.;
	    if(NormAGBMetalEjecRate[ii][jj][kk]<0.)
	      NormAGBMetalEjecRate[ii][jj][kk]=0.;
	  }

}

int find_initial_metallicity_comp(int Zi, int sfh_bin, int table_type)
{
  int i, Zi_bin;
  double Z_in;

  Zi_bin = -1;
  i = 0;
  Z_in = lifetimeMetallicities[Zi];

  switch (table_type)
  {
    case 1: //Lifetime metallicity table
      while (Zi_bin == -1)
	{
	  if (lifetimeMetallicities[i] < Z_in)
	    {
	      i++;
	      if (i == LIFETIME_Z_NUM) Zi_bin = i-1; //If galaxy's Z is higher than max Z from table, then just take max Z from table
	    }
	  else Zi_bin = i;
	}
      break;
    case 2: //SN-II metallicity table
      while (Zi_bin == -1)
	{
	  if (SNIIMetallicities[i] < Z_in)
	    {
	      i++;
	      if (i == SNII_Z_NUM) Zi_bin = i-1;
	    }
	  else Zi_bin = i;
	}
      break;
    //case 3 //SNIa yields are NOT metallicity dependent
    case 4: //AGB metallicity table
      while (Zi_bin == -1)
	{
	  if (AGBMetallicities[i] < Z_in)
	    {
	      i++;
	      if (i == AGB_Z_NUM) Zi_bin = i-1;
	    }
	  else Zi_bin = i;
	}
      break;
  }
  return Zi_bin;
}

int find_lifetime(double mass) //Function to find the bin in mass dimension of lifetimes[Zi][M] which corresponds to the "mass" given as an input.
{
	int Li_bin;
	Li_bin = -1;
	int i;
	i = 0;
	while (Li_bin == -1)
	{
		if (lifetimeMasses[i] < mass)
		{
			i++;
			if (i == LIFETIME_MASS_NUM) Li_bin = i-1;
		}
		else Li_bin = i-1;
	}
	return Li_bin; //This returns element number i for the lifetimes[Zi][i] array BELOW the true lifetime corresponding to "mass".
}

int find_initial_mass(double lifetime, int Zi_bin)
{
  if (lifetime <= 0.0) return LIFETIME_MASS_NUM-2; //If the bin 'touches now', then return max mass (ie: star of shortest lifetime) ie: bin for 120Msun
  else if (lifetime > lifetimes[Zi_bin][0]) return 0; //If true lifetime is longer than max lifetime in table (shouldn't be), then return element number 0
  else
    {
      int Mi_bin;

      Mi_bin = -1;
      int i;
      i = 0;
      while (Mi_bin == -1)
	{
	  if (lifetimes[Zi_bin][i] > lifetime)
	    {
	      i++;
	      if (i == LIFETIME_MASS_NUM-1) Mi_bin = i; //If lifetime is shorter than min lifetime from table, then just return max mass (120 Msun)
	    }
	  else Mi_bin = i;
	}
      return Mi_bin-1; //This returns element number i for lifetimeMasses[Zi][i] array BELOW the true initial mass corresponding to t_lower or t_upper.
    }
}

int max_Mi_lower(int Mi_lower, int channel_type)
//This function checks whether the mass corresponding to Mi_lower (i.e. the mass bin in the lifetime arrays
//corresponding to lowest mass of star born in SFH bin i to 'die' in current timestep) is above the minimum
//allowed mass for a particular channel (i.e. SNe-II, SNe-Ia, or AGBs).
//If it is, Mi_lower is returned. If it isn't, the mass bin corresponding to the minimum allowed mass is
//returned.
{
  switch (channel_type)
  {
    case 2: //SNII mass limits
      if (lifetimeMasses[Mi_lower] > SNII_MIN_MASS) return Mi_lower;
      else {//i.e. If lower mass is below SNII_MIN_MASS, return bin number for mass=SNII_MIN_MASS
		  int i;
		  i = 0;
		  do { i++; }
		  while (lifetimeMasses[i] < SNII_MIN_MASS);
		  return i;
      }
      break;
#ifndef DTD
    case 3: //SNIa mass limits
      if (lifetimeMasses[Mi_lower] > 0.85) return Mi_lower; //NB: Lifetimes of SNe-Ia binaries depend on M2, not Mb. (ie: 0.85<=M2<=8.0)
      else
	{
	  int i;
	  i = 0;
	  do { i++; }
	  while (lifetimeMasses[i] < 0.85); //NB: Lifetimes of SNe-Ia binaries depend on M2, not Mb. (ie: 0.85<=M2<=8.0)
	  return i;
	}
      break;
#endif
    case 4: //AGB mass limits
      if (lifetimeMasses[Mi_lower] > AGB_MIN_MASS) return Mi_lower;
      else
	{
	  int i;
	  i = 0;
	  do { i++; }
	  while (lifetimeMasses[i] < AGB_MIN_MASS);
	  return i;
	}
      break;
    default: printf("Wrong ejection mode chosen in max_Mi_lower: Use either 2 (SNe-II), 3 (SNe-Ia) or 4 (AGB winds)"); exit(1);
  }
}

int min_Mi_upper(int Mi_upper, int channel_type)
//This function checks whether the mass corresponding to Mi_upper (i.e. the mass bin in the lifetime arrays
//corresponding to highest mass of star born in SFH bin i to 'die' in current timestep) is below the maximum
//allowed mass for a particular channel (i.e. SNe-II, SNe-Ia, or AGBs).
//If it is, Mi_upper is returned. If it isn't, the mass bin corresponding to the maximum allowed mass is
//returned.
{
  switch (channel_type)
  {
    case 2: //SNII mass limits
      if (lifetimeMasses[Mi_upper] < SNII_MAX_MASS) return Mi_upper;
      //else return LIFETIME_MASS_NUM-1;
      else
	{
	  int i;
	  //i = LIFETIME_MASS_NUM-1;
	  i = LIFETIME_MASS_NUM;
	  do { i--; }
	  while (lifetimeMasses[i] > SNII_MAX_MASS);
	  return i;
	}
      break;
#ifndef DTD
    case 3: //SNIa mass limits
      if (lifetimeMasses[Mi_upper] < 8.0) return Mi_upper; //NB: Lifetimes of SNe-Ia binaries depends on M2, not Mb. (ie: 0.85<=M2<=8.0)
      else
	{
	  int i;
	  i = LIFETIME_MASS_NUM-1;
	  do { i--; }
	  while (lifetimeMasses[i] > 8.0); //NB: Lifetimes of SNe-Ia binaries depends on M2, not Mb. (ie: 0.85<=M2<=8.0)
	  return i;
	}
      break;
#endif
    case 4: //AGB mass limits
      if (lifetimeMasses[Mi_upper] < AGB_MAX_MASS) return Mi_upper;
      else
	{
	  int i;
	  //i = LIFETIME_MASS_NUM-1;
	  i = LIFETIME_MASS_NUM;
	  do { i--; }
	  while (lifetimeMasses[i] > AGB_MAX_MASS);
	  return i;
	}
      break;
    default: printf("Wrong ejection mode chosen in min_Mi_upper: Use either 2 (SNe-II), 3 (SNe-Ia) or 4 (AGB winds)"); exit(1);
  }
}

int find_SNII_mass_bin(double masslimit)
{
  if (masslimit == SNII_MAX_MASS) return SNII_MASS_NUM-1;
  else
    {
      int Mi_bin;

      Mi_bin = -1;
      int i;
      i = 0;
      while (Mi_bin == -1)
	{
	  if (SNIIMasses[i] < masslimit)
	    {
	      i++;
	      if (i == SNII_MASS_NUM) Mi_bin = i-1; //If mass is greater than max mass for SNe-II (shouldn't be), then just return max mass (120.0 Msun)
	    }
	  else Mi_bin = i;
	}
      return Mi_bin;
    }
}

int find_agb_mass_bin(double masslimit)
{
  if (masslimit == AGB_MAX_MASS) return AGB_MASS_NUM-1;
  else
    {
      int Mi_bin;

      Mi_bin = -1;
      int i;
      i = 0;
      while (Mi_bin == -1)
	{
	  if (AGBMasses[i] < masslimit)
	    {
	      i++;
	      if (i == AGB_MASS_NUM) Mi_bin = i-1; //If mass is greater than max mass for AGB winds (shouldn't be), then just return max mass (5.0 Msun)
	    }
	  else Mi_bin = i;
	}
      return Mi_bin;
    }
}

#ifdef DTD
double DTDcalc (double timevalue)
{
	if (timevalue == 0.0)
		return 0.0;
	else {
#ifdef BIMODALDTD
  double timevalueM, DTDvalueM;
  timevalueM = log10(timevalue); //IN [log(YEARS)]
  if (timevalueM < 7.93) //Characteristic time == 10^7.93 yrs (NB: timevalue is in log10(yrs) here)
    {
      DTDvalueM = (1.4 - 50.0*(timevalueM-7.7)*(timevalueM-7.7)) / DTD_NORM; //When using Mannucci bi-modal DTD
      //DTDvalueM = (0.74 - 50.0*(timevalueM-7.7)*(timevalueM-7.7)) / DTD_NORM; //When using CUSTOM bi-modal DTD, with 20% in prompt component
    }
  else
    {
      DTDvalueM = (-0.8 - 0.9*(timevalueM-8.7)*(timevalueM-8.7)) / DTD_NORM; //Same eqn for Custom bi-modal DTD too
    }
  return pow(10.0,DTDvalueM); //IN [SNe/Gyr]
#endif

#ifdef CUSTOMDTD
  double timevalueM, DTDvalueM;
  timevalueM = log10(timevalue); //IN [log(YEARS)]
  if (timevalueM < 7.93) //Characteristic time == 10^7.93 yrs (NB: timevalue is in log10(yrs) here)
    {
      //DTDvalueM = (0.74 - 50.0*(timevalueM-7.7)*(timevalueM-7.7)) / DTD_NORM; //When using CUSTOM bi-modal DTD, with 20% in prompt component (and tmin = 26Myrs)
      DTDvalueM = (1.0 - 50.0*(timevalueM-7.7)*(timevalueM-7.7)) / DTD_NORM; //When using CUSTOM bi-modal DTD, with 20% in prompt component (and tmin = 26Myrs)
    }
  else
    {
      DTDvalueM = (-0.8 - 0.9*(timevalueM-8.7)*(timevalueM-8.7)) / DTD_NORM;
    }
  return pow(10.0,DTDvalueM); //IN [SNe/Gyr]
#endif

#ifdef GAUSSIANDTD
  double timevalueG, pivalue, DTDvalueG; //, tauCharac, sigmatd
  timevalueG = timevalue/1.0e9; //IN [Gyrs]
  pivalue = 3.14159;
  DTDvalueG = ((1./sqrt(2.*pivalue*SIGMA_TD*SIGMA_TD)) * exp(-(pow((timevalueG-TAUCHARAC),2.))/(2.*SIGMA_TD*SIGMA_TD)));// / DTD_NORM;
  return DTDvalueG; //IN [SNe/Gyr]
#endif

#ifdef POWERLAWDTD
  double timevalueP, DTDvalueP;
  timevalueP = timevalue/1.0e9; //IN [Gyrs]
  DTDvalueP = pow(timevalueP, DTD_SLOPE) / DTD_NORM;
  return DTDvalueP; //IN [SNe/Gyr]
#endif

#ifdef RUITERDTD
  double timevalueR, pivalue, DTDvalueR;
  timevalueR = timevalue/1.0e9; //IN [Gyrs]
  if (timevalueR <= 1.0) //Time between Gaussian and power-law components = 1.0 Gyrs
    {
      pivalue = 3.14159;
      DTDvalueR = (0.143 * (1./sqrt(2.*pivalue*SIGMA_TD*SIGMA_TD)) * exp(-(pow((timevalueR-TAUCHARAC),2.))/(2.*SIGMA_TD*SIGMA_TD))) / DTD_NORM; //0.143 factor ensures Gaussian component is 13% of total
    }
  if (timevalueR > 1.0)
    {
      DTDvalueR = pow(timevalueR, DTD_SLOPE) / DTD_NORM;
    }
  return DTDvalueR; //IN [SNe/Gyr]
#endif
	}
}
#endif //DTD


#ifdef INDIVIDUAL_ELEMENTS
void find_actual_ejecta_limits(int channel_type, double Mi_lower_actual, double Mi_upper_actual, int Mi_lower, int Mi_upper, int Zi,
			       double* EjectedMasses_lower_actual, double* EjectedMasses_upper_actual, double* TotalMetals_lower_actual, double* TotalMetals_upper_actual,
			       double* Yields_lower_actual, double* Yields_upper_actual)
#else
void find_actual_ejecta_limits(int channel_type, double Mi_lower_actual, double Mi_upper_actual, int Mi_lower, int Mi_upper, int Zi,
			       double* EjectedMasses_lower_actual, double* EjectedMasses_upper_actual, double* TotalMetals_lower_actual, double* TotalMetals_upper_actual)
#endif
{
  switch (channel_type)
  {
    case 2: //SNII
      if (Mi_lower_actual <= SNII_MIN_MASS) {
		  //printf("find_actual_ejecta_limits(): Zi = %i | Mi_lower = %i | SNIIEjectedMasses[Zi][Mi_lower] = %f\n", Zi, Mi_lower, SNIIEjectedMasses[Zi][Mi_lower]);
		  *EjectedMasses_lower_actual = SNIIEjectedMasses[Zi][Mi_lower]; //If Mi_lower_actual was less than or equal to Mi_lower (e.g. 7 Msun bin for SNe-II), then EjectedMasses_lower_actual is set to ejected mass from an e.g. 7 Msun star.
		  *TotalMetals_lower_actual = SNIITotalMetals[Zi][Mi_lower]; //If Mi_lower_actual was less than or equal to Mi_lower (e.g. 7 Msun bin for SNe-II), then TotalMetals_lower_actual is set to ejected mass in metals from an e.g. 7 Msun star.
#ifdef INDIVIDUAL_ELEMENTS
		  int k;
		  for (k=0;k<NUM_ELEMENTS;k++) {
			  int kk;
#ifndef MAINELEMENTS
			  kk=k;
			  Yields_lower_actual[k] = SNIIYields[Zi][kk][Mi_lower]; //Note that Zi is actually Zi_SNII here (see the variables where find_actual_ejecta_limits() is called)
#else
			  switch(k){case 0: kk=0; break; case 1: kk=1; break; case 2: kk=4; break; case 3: kk=6; break; case 4: kk=10; break;}
			  Yields_lower_actual[k] = SNIIYields[Zi][kk][Mi_lower];
#endif
		  }
#endif //INDIVIDUAL_ELEMENTS
      }
      else {
		  *EjectedMasses_lower_actual = SNIIEjectedMasses[Zi][Mi_lower] + ((SNIIEjectedMasses[Zi][Mi_lower+1]-SNIIEjectedMasses[Zi][Mi_lower]) * ((Mi_lower_actual-SNIIMasses[Mi_lower])/(SNIIMasses[Mi_lower+1]-SNIIMasses[Mi_lower])));
		  *TotalMetals_lower_actual = SNIITotalMetals[Zi][Mi_lower] + ((SNIITotalMetals[Zi][Mi_lower+1]-SNIITotalMetals[Zi][Mi_lower]) * ((Mi_lower_actual-SNIIMasses[Mi_lower])/(SNIIMasses[Mi_lower+1]-SNIIMasses[Mi_lower])));
#ifdef INDIVIDUAL_ELEMENTS
		  int k;
		  for (k=0;k<NUM_ELEMENTS;k++) {
			  int kk;
#ifndef MAINELEMENTS
			  kk=k;
			  Yields_lower_actual[k] = SNIIYields[Zi][kk][Mi_lower] + ((SNIIYields[Zi][kk][Mi_lower+1]-SNIIYields[Zi][kk][Mi_lower]) * ((Mi_lower_actual-SNIIMasses[Mi_lower])/(SNIIMasses[Mi_lower+1]-SNIIMasses[Mi_lower])));
#else
			  switch(k){case 0: kk=0; break; case 1: kk=1; break; case 2: kk=4; break; case 3: kk=6; break; case 4: kk=10; break;}
			  Yields_lower_actual[k] = SNIIYields[Zi][kk][Mi_lower] + ((SNIIYields[Zi][kk][Mi_lower+1]-SNIIYields[Zi][kk][Mi_lower]) * ((Mi_lower_actual-SNIIMasses[Mi_lower])/(SNIIMasses[Mi_lower+1]-SNIIMasses[Mi_lower])));
#endif
		  }
#endif //INDIVIDUAL_ELEMENTS
      }

      if (Mi_upper == SNII_MASS_NUM-1)
	{
	  *EjectedMasses_upper_actual = SNIIEjectedMasses[Zi][SNII_MASS_NUM-1]; //If Mi_upper_actual was more than or equal to 120 Msun, then EjectedMasses_upper_actual is set to ejected mass from 120 Msun star.
	  *TotalMetals_upper_actual = SNIITotalMetals[Zi][SNII_MASS_NUM-1]; //If Mi_upper_actual was more than or equal to 120 Msun, then TotalMetals_upper_actual is set to ejected mass in metals from 120 Msun star.
#ifdef INDIVIDUAL_ELEMENTS
	  int k;
	  for (k=0;k<NUM_ELEMENTS;k++)
	    {
	      int kk;
#ifndef MAINELEMENTS
	      kk=k;
	      Yields_upper_actual[k] = SNIIYields[Zi][kk][SNII_MASS_NUM-1];
#else
	      switch(k){case 0: kk=0; break; case 1: kk=1; break; case 2: kk=4; break; case 3: kk=6; break; case 4: kk=10; break;}
	      Yields_upper_actual[k] = SNIIYields[Zi][kk][SNII_MASS_NUM-1];
#endif
	    }
#endif //INDIVIDUAL_ELEMENTS
	}
      else
	{
	  *EjectedMasses_upper_actual = SNIIEjectedMasses[Zi][Mi_upper] + ((SNIIEjectedMasses[Zi][Mi_upper+1]-SNIIEjectedMasses[Zi][Mi_upper]) * ((Mi_upper_actual-SNIIMasses[Mi_upper])/(SNIIMasses[Mi_upper+1]-SNIIMasses[Mi_upper])));
	  *TotalMetals_upper_actual = SNIITotalMetals[Zi][Mi_upper] + ((SNIITotalMetals[Zi][Mi_upper+1]-SNIITotalMetals[Zi][Mi_upper]) * ((Mi_upper_actual-SNIIMasses[Mi_upper])/(SNIIMasses[Mi_upper+1]-SNIIMasses[Mi_upper])));
#ifdef INDIVIDUAL_ELEMENTS
	  int k;
	  for (k=0;k<NUM_ELEMENTS;k++)
	    {
	      int kk;
#ifndef MAINELEMENTS
	      kk=k;
	      Yields_upper_actual[k] = SNIIYields[Zi][kk][Mi_upper] + ((SNIIYields[Zi][kk][Mi_upper+1]-SNIIYields[Zi][kk][Mi_upper]) * ((Mi_upper_actual-SNIIMasses[Mi_upper])/(SNIIMasses[Mi_upper+1]-SNIIMasses[Mi_upper])));
#else
	      switch(k){case 0: kk=0; break; case 1: kk=1; break; case 2: kk=4; break; case 3: kk=6; break; case 4: kk=10; break;}
	      Yields_upper_actual[k] = SNIIYields[Zi][kk][Mi_upper] + ((SNIIYields[Zi][kk][Mi_upper+1]-SNIIYields[Zi][kk][Mi_upper]) * ((Mi_upper_actual-SNIIMasses[Mi_upper])/(SNIIMasses[Mi_upper+1]-SNIIMasses[Mi_upper])));
#endif
	    }
#endif //INDIVIDUAL_ELEMENTS
	}
      break;

    case 4: //AGB
      if (Mi_lower == 0)
	{
	  *EjectedMasses_lower_actual = AGBEjectedMasses[Zi][0]; //If Mi_lower_actual was less than or equal to Mi_lower (e.g. 0.85 Msun bin for AGB), then EjectedMasses_lower_actual is set to ejected mass from an e.g. 0.85 Msun star.
	  *TotalMetals_lower_actual = AGBTotalMetals[Zi][0]; //If Mi_lower_actual was less than or equal to Mi_lower (e.g. 0.85 Msun bin for AGB), then TotalMetals_lower_actual is set to ejected mass in metals from an e.g. 0.85 Msun star.
#ifdef INDIVIDUAL_ELEMENTS
	  int k;
	  for (k=0;k<NUM_ELEMENTS;k++)
	    {
	      int kk;
#ifndef MAINELEMENTS
	      kk=k;
	      Yields_lower_actual[k] = AGBYields[Zi][kk][0];
#else
	      switch(k){case 0: kk=0; break; case 1: kk=1; break; case 2: kk=4; break; case 3: kk=6; break; case 4: kk=10; break;}
	      Yields_lower_actual[k] = AGBYields[Zi][kk][0];
#endif
	    }
	  if(Yields_lower_actual[10] > 0.0) {printf("YLA1: Yields_lower_actual[10] = %f\n", Yields_lower_actual[10]);}
#endif //INDIVIDUAL_ELEMENTS
	}
      else if (Mi_lower+1 == AGB_MASS_NUM) //This option added to prevent numerical integration beyond the maximum AGB mass bin with AGBMasses[Mi_lower+1]. ROB (18-06-21)
	{
	  *EjectedMasses_lower_actual = AGBEjectedMasses[Zi][AGB_MASS_NUM-1];
	  *TotalMetals_lower_actual = AGBTotalMetals[Zi][AGB_MASS_NUM-1];
#ifdef INDIVIDUAL_ELEMENTS
	  int k;
	  for (k=0;k<NUM_ELEMENTS;k++)
	    {
	      int kk;
#ifndef MAINELEMENTS
	      kk=k;
	      Yields_lower_actual[k] = AGBYields[Zi][kk][AGB_MASS_NUM-1];
#else
	      switch(k){case 0: kk=0; break; case 1: kk=1; break; case 2: kk=4; break; case 3: kk=6; break; case 4: kk=10; break;}
	      Yields_lower_actual[k] = AGBYields[Zi][kk][AGB_MASS_NUM-1];
#endif
	    }
#endif //INDIVIDUAL_ELEMENTS
	}
      else
      	{
      	  *EjectedMasses_lower_actual = AGBEjectedMasses[Zi][Mi_lower] + ((AGBEjectedMasses[Zi][Mi_lower+1]-AGBEjectedMasses[Zi][Mi_lower]) * ((Mi_lower_actual-AGBMasses[Mi_lower])/(AGBMasses[Mi_lower+1]-AGBMasses[Mi_lower])));
      	  *TotalMetals_lower_actual = AGBTotalMetals[Zi][Mi_lower] + ((AGBTotalMetals[Zi][Mi_lower+1]-AGBTotalMetals[Zi][Mi_lower]) * ((Mi_lower_actual-AGBMasses[Mi_lower])/(AGBMasses[Mi_lower+1]-AGBMasses[Mi_lower])));
      #ifdef INDIVIDUAL_ELEMENTS
      	  int k;
      	  for (k=0;k<NUM_ELEMENTS;k++)
      	    {
      	      int kk;
      #ifndef MAINELEMENTS
      	      kk=k;
      	      Yields_lower_actual[k] = AGBYields[Zi][kk][Mi_lower] + ((AGBYields[Zi][kk][Mi_lower+1]-AGBYields[Zi][kk][Mi_lower]) * ((Mi_lower_actual-AGBMasses[Mi_lower])/(AGBMasses[Mi_lower+1]-AGBMasses[Mi_lower])));
      #else
      	      switch(k){case 0: kk=0; break; case 1: kk=1; break; case 2: kk=4; break; case 3: kk=6; break; case 4: kk=10; break;}
      	      Yields_lower_actual[k] = AGBYields[Zi][kk][Mi_lower] + ((AGBYields[Zi][kk][Mi_lower+1]-AGBYields[Zi][kk][Mi_lower]) * ((Mi_lower_actual-AGBMasses[Mi_lower])/(AGBMasses[Mi_lower+1]-AGBMasses[Mi_lower])));
      #endif
      	    }
      #endif //INDIVIDUAL_ELEMENTS
      	}

      if (Mi_upper == AGB_MASS_NUM-1)
	{
	  *EjectedMasses_upper_actual = AGBEjectedMasses[Zi][AGB_MASS_NUM-1]; //If Mi_upper_actual was more than or equal to 5 Msun, then EjectedMasses_upper_actual is set to ejected mass from 5 Msun star.
	  *TotalMetals_upper_actual = AGBTotalMetals[Zi][AGB_MASS_NUM-1]; //If Mi_upper_actual was more than or equal to 5 Msun, then TotalMetals_upper_actual is set to ejected mass in metals from 5 Msun star.
#ifdef INDIVIDUAL_ELEMENTS
	  int k;
	  for (k=0;k<NUM_ELEMENTS;k++)
	    {
	      int kk;
#ifndef MAINELEMENTS
	      kk=k;
	      Yields_upper_actual[k] = AGBYields[Zi][kk][AGB_MASS_NUM-1];
#else
	      switch(k){case 0: kk=0; break; case 1: kk=1; break; case 2: kk=4; break; case 3: kk=6; break; case 4: kk=10; break;}
	      Yields_upper_actual[k] = AGBYields[Zi][kk][AGB_MASS_NUM-1];
#endif
	    }
#endif //INDIVIDUAL_ELEMENTS
	}
      else
	{
	  *EjectedMasses_upper_actual = AGBEjectedMasses[Zi][Mi_upper] + ((AGBEjectedMasses[Zi][Mi_upper+1]-AGBEjectedMasses[Zi][Mi_upper]) * ((Mi_upper_actual-AGBMasses[Mi_upper])/(AGBMasses[Mi_upper+1]-AGBMasses[Mi_upper])));
	  *TotalMetals_upper_actual = AGBTotalMetals[Zi][Mi_upper] + ((AGBTotalMetals[Zi][Mi_upper+1]-AGBTotalMetals[Zi][Mi_upper]) * ((Mi_upper_actual-AGBMasses[Mi_upper])/(AGBMasses[Mi_upper+1]-AGBMasses[Mi_upper])));
#ifdef INDIVIDUAL_ELEMENTS
	  int k;
	  for (k=0;k<NUM_ELEMENTS;k++)
	    {
	      int kk;
#ifndef MAINELEMENTS
	      kk=k;
	      Yields_upper_actual[k] = AGBYields[Zi][kk][Mi_upper] + ((AGBYields[Zi][kk][Mi_upper+1]-AGBYields[Zi][kk][Mi_upper]) * ((Mi_upper_actual-AGBMasses[Mi_upper])/(AGBMasses[Mi_upper+1]-AGBMasses[Mi_upper])));
#else
	      switch(k){case 0: kk=0; break; case 1: kk=1; break; case 2: kk=4; break; case 3: kk=6; break; case 4: kk=10; break;}
	      Yields_upper_actual[k] = AGBYields[Zi][kk][Mi_upper] + ((AGBYields[Zi][kk][Mi_upper+1]-AGBYields[Zi][kk][Mi_upper]) * ((Mi_upper_actual-AGBMasses[Mi_upper])/(AGBMasses[Mi_upper+1]-AGBMasses[Mi_upper])));
#endif
	    }
#endif //INDIVIDUAL_ELEMENTS
	}
      break;
  }
}


#ifdef INDIVIDUAL_ELEMENTS
//double Masses[], double EjectedMasses[AGB_Z_NUM][AGB_MASS_NUM], double TotalMetals[AGB_Z_NUM][AGB_MASS_NUM], double Yields[AGB_Z_NUM][NUM_ELEMENTS][AGB_MASS_NUM],
void calc_mass_limits(int j, int Zi_chan, int Mi_lower_chan, int Mi_upper_chan, double Mi_lower_actual, double Mi_upper_actual,
					  int MASS_NUM, int Z_NUM, double Masses[MASS_NUM], double EjectedMasses[Z_NUM][MASS_NUM], double TotalMetals[Z_NUM][MASS_NUM], double Yields[Z_NUM][NUM_ELEMENTS][MASS_NUM],
					  double EjectedMasses_lower_actual, double EjectedMasses_upper_actual, double TotalMetals_lower_actual, double TotalMetals_upper_actual,
					  double Yields_lower_actual[], double Yields_upper_actual[],
					  double* ML, double* MU, double* EjectedMasses_lower, double* EjectedMasses_upper, double* TotalMetals_lower, double* TotalMetals_upper,
					  double* Yields_lower, double* Yields_upper)
#else
void calc_mass_limits(int j, int Zi_chan, int Mi_lower_chan, int Mi_upper_chan, double Mi_lower_actual, double Mi_upper_actual,
					  int MASS_NUM, int Z_NUM, double Masses[MASS_NUM], double EjectedMasses[Z_NUM][MASS_NUM], double TotalMetals[Z_NUM][MASS_NUM],
					  double EjectedMasses_lower_actual, double EjectedMasses_upper_actual, double TotalMetals_lower_actual, double TotalMetals_upper_actual,
					  double* ML, double* MU, double* EjectedMasses_lower, double* EjectedMasses_upper, double* TotalMetals_lower, double* TotalMetals_upper)
#endif //INDIVIDUAL_ELEMENTS
{
	int kk;
	//1) If mass bin j is NEITHER the lowest NOR highest mass bin to be integrated over:
	if (j != Mi_lower_chan && j != Mi_upper_chan) {
		*ML = Masses[j];
		*MU = Masses[j+1];
		*EjectedMasses_lower = EjectedMasses[Zi_chan][j];
		*EjectedMasses_upper = EjectedMasses[Zi_chan][j+1];
		*TotalMetals_lower = TotalMetals[Zi_chan][j];
		*TotalMetals_upper = TotalMetals[Zi_chan][j+1];
#ifdef INDIVIDUAL_ELEMENTS
		for (kk=0; kk<NUM_ELEMENTS; kk++) {
			Yields_lower[kk] = Yields[Zi_chan][kk][j];
			Yields_upper[kk] = Yields[Zi_chan][kk][j+1];
		}
#endif //INDIVIDUAL_ELEMENTS
	}
	//2) If mass bin j is BOTH the lowest AND highest mass bin to be integrated over:
	else if (j == Mi_lower_chan && j == Mi_upper_chan) {
		*ML = Mi_lower_actual;
		*MU = Mi_upper_actual;
		*EjectedMasses_lower = EjectedMasses_lower_actual;
		*EjectedMasses_upper = EjectedMasses_upper_actual;
		*TotalMetals_lower = TotalMetals_lower_actual;
		*TotalMetals_upper = TotalMetals_upper_actual;
#ifdef INDIVIDUAL_ELEMENTS
		for (kk=0; kk<NUM_ELEMENTS; kk++) {
			Yields_lower[kk] = Yields_lower_actual[kk];
			Yields_upper[kk] = Yields_upper_actual[kk];
		}
#endif //INDIVIDUAL_ELEMENTS
	}
	//3) If mass bin j IS the lowest but IS NOT the highest mass bin to be integrated over:
	else if (j == Mi_lower_chan && j != Mi_upper_chan) {
		*ML = Mi_lower_actual;
		*MU = Masses[j+1];
		*EjectedMasses_lower = EjectedMasses_lower_actual;
		*EjectedMasses_upper = EjectedMasses[Zi_chan][j+1];
		*TotalMetals_lower = TotalMetals_lower_actual;
		*TotalMetals_upper = TotalMetals[Zi_chan][j+1];
#ifdef INDIVIDUAL_ELEMENTS
		for (kk=0; kk<NUM_ELEMENTS; kk++) {
			Yields_lower[kk] = Yields_lower_actual[kk];
			Yields_upper[kk] = Yields[Zi_chan][kk][j+1];
		}
#endif //INDIVIDUAL_ELEMENTS
	}
	//4) If mass bin j IS NOT the lowest but IS the highest mass bin to be integrated over:
	else if (j != Mi_lower_chan && j == Mi_upper_chan) {
		*ML = Masses[j];
		*MU = Mi_upper_actual;
		*EjectedMasses_lower = EjectedMasses[Zi_chan][j];
		*EjectedMasses_upper = EjectedMasses_upper_actual;
		*TotalMetals_lower = TotalMetals[Zi_chan][j];
		*TotalMetals_upper = TotalMetals_upper_actual;
#ifdef INDIVIDUAL_ELEMENTS
		for (kk=0; kk<NUM_ELEMENTS; kk++) {
			Yields_lower[kk] = Yields[Zi_chan][kk][j];
			Yields_upper[kk] = Yields_upper_actual[kk];
		}
#endif //INDIVIDUAL_ELEMENTS
	}
	else printf("yield_integrals.c: calc_mass_limits(): No option entered!\n");
}


#ifdef INDIVIDUAL_ELEMENTS
void calc_mass_limits_SNII(int j, int Zi_SNII, int Mi_lower_SNII, int Mi_upper_SNII, double Mi_lower_actual, double Mi_upper_actual,
					  double SNIIEjectedMasses_lower_actual, double SNIIEjectedMasses_upper_actual, double SNIITotalMetals_lower_actual, double SNIITotalMetals_upper_actual,
					  double SNIIYields_lower_actual[], double SNIIYields_upper_actual[],
					  double* ML, double* MU, double* EjectedMasses_lower, double* EjectedMasses_upper, double* TotalMetals_lower, double* TotalMetals_upper,
					  double* Yields_lower, double* Yields_upper)
#else
void calc_mass_limits_SNII(int j, int kk, int Zi_SNII, int Mi_lower_SNII, int Mi_upper_SNII, double Mi_lower_actual, double Mi_upper_actual,
					  double SNIIEjectedMasses_lower_actual, double SNIIEjectedMasses_upper_actual, double SNIITotalMetals_lower_actual, double SNIITotalMetals_upper_actual,
					  double* ML, double* MU, double* EjectedMasses_lower, double* EjectedMasses_upper, double* TotalMetals_lower, double* TotalMetals_upper)
#endif //INDIVIDUAL_ELEMENTS
{
	int kk;
	//1) If mass bin j is NEITHER the lowest NOR highest mass bin to be integrated over:
	//This will work for either (SNIIMasses[j] <= SNIA_MAX_MASS) or (SNIIMasses[j] > SNIA_MAX_MASS)
	if (j != Mi_lower_SNII && j != Mi_upper_SNII) {
		*ML = SNIIMasses[j];
		*MU = SNIIMasses[j+1];
		*EjectedMasses_lower = SNIIEjectedMasses[Zi_SNII][j];
		*EjectedMasses_upper = SNIIEjectedMasses[Zi_SNII][j+1];
		*TotalMetals_lower = SNIITotalMetals[Zi_SNII][j];
		*TotalMetals_upper = SNIITotalMetals[Zi_SNII][j+1];
#ifdef INDIVIDUAL_ELEMENTS
		for (kk=0; kk<NUM_ELEMENTS; kk++) {
			Yields_lower[kk] = SNIIYields[Zi_SNII][kk][j];
			Yields_upper[kk] = SNIIYields[Zi_SNII][kk][j+1];
		}
#endif //INDIVIDUAL_ELEMENTS
	}
	//2) If mass bin j is BOTH the lowest AND highest mass bin to be integrated over:
	//This will work for either (SNIIMasses[j] <= SNIA_MAX_MASS) or (SNIIMasses[j] > SNIA_MAX_MASS)
	else if (j == Mi_lower_SNII && j == Mi_upper_SNII) { // && Mi_lower_actual >= SNII_MIN_MASS
		*ML = Mi_lower_actual;
		*MU = Mi_upper_actual;
		*EjectedMasses_lower = SNIIEjectedMasses_lower_actual;
		*EjectedMasses_upper = SNIIEjectedMasses_upper_actual;
		*TotalMetals_lower = SNIITotalMetals_lower_actual;
		*TotalMetals_upper = SNIITotalMetals_upper_actual;
#ifdef INDIVIDUAL_ELEMENTS
		for (kk=0; kk<NUM_ELEMENTS; kk++) {
			Yields_lower[kk] = SNIIYields_lower_actual[kk];
			Yields_upper[kk] = SNIIYields_upper_actual[kk];
		}
#endif //INDIVIDUAL_ELEMENTS
	}
	/*No longer required, as Mi_lower_actual >= SNII_MIN_MASS is now already enforced above:
	//3) If mass bin j is BOTH the lowest AND highest mass bin to be integrated over, AND 'Mi_lower_actual' is below min. mass for SNe-II. (Still only counts stars from SNII_MIN_MASS to Mi_upper_actual):
	//This will only be entered for (SNIIMasses[j] <= SNIA_MAX_MASS)
	else if (j == Mi_lower_SNII && j == Mi_upper_SNII && Mi_lower_actual < SNII_MIN_MASS) {
		*ML = SNII_MIN_MASS;
		*MU = Mi_upper_actual;
		*EjectedMasses_lower = SNIIEjectedMasses_lower_actual; //Is this ok? Did find_actual_ejecta_limits() account for Mi_lower_actual < SNII_MIN_MASS when calculating SNIIEjectedMasses_lower_actual? <-- Yes, I think it does, os this is ok. (26-04-22)
		*EjectedMasses_upper = SNIIEjectedMasses_upper_actual;
		*TotalMetals_lower = SNIITotalMetals_lower_actual;
		*TotalMetals_upper = SNIITotalMetals_upper_actual;
#ifdef INDIVIDUAL_ELEMENTS
		for (kk=0; kk<NUM_ELEMENTS; kk++) {
			Yields_lower[kk] = SNIIYields_lower_actual[kk];
			Yields_upper[kk] = SNIIYields_upper_actual[kk];
		}
#endif //INDIVIDUAL_ELEMENTS
	}*/
	//3) If mass bin j IS the lowest but IS NOT the highest mass bin to be integrated over:
	//This will work for either (SNIIMasses[j] <= SNIA_MAX_MASS) or (SNIIMasses[j] > SNIA_MAX_MASS)
	else if (j == Mi_lower_SNII && j != Mi_upper_SNII) { // && Mi_lower_actual >= SNII_MIN_MASS
		*ML = Mi_lower_actual;
		*MU = SNIIMasses[j+1];
		*EjectedMasses_lower = SNIIEjectedMasses_lower_actual;
		*EjectedMasses_upper = SNIIEjectedMasses[Zi_SNII][j+1];
		*TotalMetals_lower = SNIITotalMetals_lower_actual;
		*TotalMetals_upper = SNIITotalMetals[Zi_SNII][j+1];
#ifdef INDIVIDUAL_ELEMENTS
		for (kk=0; kk<NUM_ELEMENTS; kk++) {
			Yields_lower[kk] = SNIIYields_lower_actual[kk];
			Yields_upper[kk] = SNIIYields[Zi_SNII][kk][j+1];
		}
#endif //INDIVIDUAL_ELEMENTS
	}
	/*No longer required, as Mi_lower_actual >= SNII_MIN_MASS is now already enforced above:
	//5) If mass bin j IS the lowest but IS NOT the highest mass bin to be integrated over, and Mi_lower_actual < SNII_MIN_MASS:
	//This will only be entered for (SNIIMasses[j] <= SNIA_MAX_MASS)
	else if (j == Mi_lower_SNII && j != Mi_upper_SNII && Mi_lower_actual < SNII_MIN_MASS) {
		*ML = SNII_MIN_MASS;
		*MU = SNIIMasses[j+1];
		*EjectedMasses_lower = SNIIEjectedMasses_lower_actual;
		*EjectedMasses_upper = SNIIEjectedMasses[Zi_SNII][j+1];
		*TotalMetals_lower = SNIITotalMetals_lower_actual;
		*TotalMetals_upper = SNIITotalMetals[Zi_SNII][j+1];
#ifdef INDIVIDUAL_ELEMENTS
		for (kk=0; kk<NUM_ELEMENTS; kk++) {
			Yields_lower[kk] = SNIIYields_lower_actual[kk];
			Yields_upper[kk] = SNIIYields[Zi_SNII][kk][j+1];
		}
#endif //INDIVIDUAL_ELEMENTS
	}*/
	//4) If mass bin j IS NOT the lowest but IS the highest mass bin to be integrated over:
	//This will work for either (SNIIMasses[j] <= SNIA_MAX_MASS) or (SNIIMasses[j] > SNIA_MAX_MASS)
	else if (j != Mi_lower_SNII && j == Mi_upper_SNII) {
		*ML = SNIIMasses[j];
		*MU = Mi_upper_actual;
		*EjectedMasses_lower = SNIIEjectedMasses[Zi_SNII][j];
		*EjectedMasses_upper = SNIIEjectedMasses_upper_actual;
		*TotalMetals_lower = SNIITotalMetals[Zi_SNII][j];
		*TotalMetals_upper = SNIITotalMetals_upper_actual;
#ifdef INDIVIDUAL_ELEMENTS
		for (kk=0; kk<NUM_ELEMENTS; kk++) {
			Yields_lower[kk] = SNIIYields[Zi_SNII][kk][j];
			Yields_upper[kk] = SNIIYields_upper_actual[kk];
		}
#endif //INDIVIDUAL_ELEMENTS
	}
	else printf("yield_integrals.c: calc_mass_limits_SNII(): No loop entered!\n");
}


#ifdef INDIVIDUAL_ELEMENTS
void calc_mass_limits_AGB(int j, int Zi_AGB, int Mi_lower_AGB, int Mi_upper_AGB, double Mi_lower_actual, double Mi_upper_actual,
					  double AGBEjectedMasses_lower_actual, double AGBEjectedMasses_upper_actual, double AGBTotalMetals_lower_actual, double AGBTotalMetals_upper_actual,
					  double AGBYields_lower_actual[], double AGBYields_upper_actual[],
					  double* ML, double* MU, double* EjectedMasses_lower, double* EjectedMasses_upper, double* TotalMetals_lower, double* TotalMetals_upper,
					  double* Yields_lower, double* Yields_upper)
#else
void calc_mass_limits_AGB(int j, int kk, int Zi_AGB, int Mi_lower_AGB, int Mi_upper_AGB, double Mi_lower_actual, double Mi_upper_actual,
					  double AGBEjectedMasses_lower_actual, double AGBEjectedMasses_upper_actual, double AGBTotalMetals_lower_actual, double AGBTotalMetals_upper_actual,
					  double* ML, double* MU, double* EjectedMasses_lower, double* EjectedMasses_upper, double* TotalMetals_lower, double* TotalMetals_upper)
#endif //INDIVIDUAL_ELEMENTS
{
	int kk;
	//1) If mass bin j is NEITHER the lowest NOR highest mass bin to be integrated over:
	if (j != Mi_lower_AGB && j != Mi_upper_AGB) {
		*ML = AGBMasses[j];
		*MU = AGBMasses[j+1];
		*EjectedMasses_lower = AGBEjectedMasses[Zi_AGB][j];
		*EjectedMasses_upper = AGBEjectedMasses[Zi_AGB][j+1];
		*TotalMetals_lower = AGBTotalMetals[Zi_AGB][j];
		*TotalMetals_upper = AGBTotalMetals[Zi_AGB][j+1];
#ifdef INDIVIDUAL_ELEMENTS
		for (kk=0; kk<NUM_ELEMENTS; kk++) {
			Yields_lower[kk] = AGBYields[Zi_AGB][kk][j];
			Yields_upper[kk] = AGBYields[Zi_AGB][kk][j+1];
		}
#endif //INDIVIDUAL_ELEMENTS
	}
	//2) If mass bin j is BOTH the lowest AND highest mass bin to be integrated over:
	else if (j == Mi_lower_AGB && j == Mi_upper_AGB) { // && Mi_lower_actual >= AGB_MIN_MASS && Mi_lower_actual <= AGB_MAX_MASS //&& Mi_upper_actual <= AGB_MAX_MASS
		*ML = Mi_lower_actual;
		*MU = Mi_upper_actual;
		*EjectedMasses_lower = AGBEjectedMasses_lower_actual;
		*EjectedMasses_upper = AGBEjectedMasses_upper_actual;
		*TotalMetals_lower = AGBTotalMetals_lower_actual;
		*TotalMetals_upper = AGBTotalMetals_upper_actual;
#ifdef INDIVIDUAL_ELEMENTS
		for (kk=0; kk<NUM_ELEMENTS; kk++) {
			Yields_lower[kk] = AGBYields_lower_actual[kk];
			Yields_upper[kk] = AGBYields_upper_actual[kk];
		}
#endif //INDIVIDUAL_ELEMENTS
	}
	/*//3) If full mass range is above limit for AGB stars: (No longer needed, as Mi_lower_actual <= AGB_MAX_MASS is already checked before entering the AGB part above) (26-04-22)
	else if (Mi_lower_actual > AGB_MAX_MASS) { //(j == Mi_lower_AGB && j == Mi_upper_AGB && Mi_lower_actual > AGB_MAX_MASS && Mi_upper_actual <= AGB_MAX_MASS)
		printf("TEST1: Loop entered: Mi_lower_actual = %f | Mi_upper_actual = %f | AGB_MAX_MASS = %f\n", Mi_lower_actual, Mi_upper_actual, AGB_MAX_MASS);
		*ML = 0.0;
		*MU = 0.0;
		*EjectedMasses_lower = 0.0;
		*EjectedMasses_upper = 0.0;
		*TotalMetals_lower = 0.0;
		*TotalMetals_upper = 0.0;
#ifdef INDIVIDUAL_ELEMENTS
		for (kk=0; kk<NUM_ELEMENTS; kk++) {
			Yields_lower[kk] = 0.0;
			Yields_upper[kk] = 0.0;
		}
#endif //INDIVIDUAL_ELEMENTS
	}*/
	/*No longer required, as Mi_upper_actual <= AGB_MAX_MASS is now already enforced above:
	//3) If mass bin j is BOTH the lowest AND highest mass bin to be integrated over, AND Mi_upper_actual > AGB_MAX_MASS:
	else if (j == Mi_lower_AGB && j == Mi_upper_AGB && Mi_upper_actual > AGB_MAX_MASS) {
		*ML = Mi_lower_actual;
		*MU = AGB_MAX_MASS;
		*EjectedMasses_lower = AGBEjectedMasses_lower_actual;
		*EjectedMasses_upper = AGBEjectedMasses_upper_actual;
		*TotalMetals_lower = AGBTotalMetals_lower_actual;
		*TotalMetals_upper = AGBTotalMetals_upper_actual;
#ifdef INDIVIDUAL_ELEMENTS
		for (kk=0; kk<NUM_ELEMENTS; kk++) {
			Yields_lower[kk] = AGBYields_lower_actual[kk];
			Yields_upper[kk] = AGBYields_upper_actual[kk];
		}
#endif //INDIVIDUAL_ELEMENTS
	}*/
	//3) If mass bin j IS the lowest but IS NOT the highest mass bin to be integrated over:
	else if (j == Mi_lower_AGB && j != Mi_upper_AGB) { //&& Mi_lower_actual >= AGB_MIN_MASS //&& Mi_lower_actual <= AGB_MAX_MASS
		*ML = Mi_lower_actual;
		*MU = AGBMasses[j+1];
		*EjectedMasses_lower = AGBEjectedMasses_lower_actual;
		*EjectedMasses_upper = AGBEjectedMasses[Zi_AGB][j+1];
		*TotalMetals_lower = AGBTotalMetals_lower_actual;
		*TotalMetals_upper = AGBTotalMetals[Zi_AGB][j+1];
#ifdef INDIVIDUAL_ELEMENTS
		for (kk=0; kk<NUM_ELEMENTS; kk++) {
			Yields_lower[kk] = AGBYields_lower_actual[kk];
			Yields_upper[kk] = AGBYields[Zi_AGB][kk][j+1];
		}
#endif //INDIVIDUAL_ELEMENTS
	}
	//4) If mass bin j IS NOT the lowest but IS the highest mass bin to be integrated over:
	else if (j != Mi_lower_AGB && j == Mi_upper_AGB) { // && Mi_upper_actual <= AGB_MAX_MASS
		*ML = AGBMasses[j];
		*MU = Mi_upper_actual;
		*EjectedMasses_lower = AGBEjectedMasses[Zi_AGB][j];
		*EjectedMasses_upper = AGBEjectedMasses_upper_actual;
		*TotalMetals_lower = AGBTotalMetals[Zi_AGB][j];
		*TotalMetals_upper = AGBTotalMetals_upper_actual;
#ifdef INDIVIDUAL_ELEMENTS
		for (kk=0; kk<NUM_ELEMENTS; kk++) {
			Yields_lower[kk] = AGBYields[Zi_AGB][kk][j];
			Yields_upper[kk] = AGBYields_upper_actual[kk];
		}
#endif //INDIVIDUAL_ELEMENTS
	}
	/*No longer required, as Mi_upper_actual <= AGB_MAX_MASS is now already enforced above:
	//6) If mass bin j IS NOT the lowest but IS the highest mass bin to be integrated over, AND Mi_upper_actual > AGB_MAX_MASS:
	else if (j != Mi_lower_AGB && j == Mi_upper_AGB && Mi_upper_actual > AGB_MAX_MASS) {
		*ML = AGBMasses[j];
		*MU = AGB_MAX_MASS;
		*EjectedMasses_lower = AGBEjectedMasses[Zi_AGB][j];
		*EjectedMasses_upper = AGBEjectedMasses_upper_actual;
		*TotalMetals_lower = AGBTotalMetals[Zi_AGB][j];
		*TotalMetals_upper = AGBTotalMetals_upper_actual;
#ifdef INDIVIDUAL_ELEMENTS
		for (kk=0; kk<NUM_ELEMENTS; kk++) {
			Yields_lower[kk] = AGBYields[Zi_AGB][kk][j];
			Yields_upper[kk] = AGBYields_upper_actual[kk];
		}
#endif //INDIVIDUAL_ELEMENTS
	}*/
	else printf("yield_integrals.c: calc_mass_limits_AGB(): No loop entered!\n");
}


void calc_time_limits_SNIa(int j, int Zi, int t_lower_lifetime, int t_upper_lifetime, double t_lower, double t_upper,
						   double* lifetimes_lower, double* lifetimes_upper, int* stat)
{
	//1) If lifetime bin j is NEITHER the lowest NOR the highest bin to be integrated over:
	if (j != t_lower_lifetime && j != t_upper_lifetime && lifetimes[Zi][j+1] >= SNIA_MIN_TIME) {
		*lifetimes_lower = lifetimes[Zi][j+1]; //i.e. the shorter of the two lifetime limits (corresp. to the upper mass limit)
		*lifetimes_upper = lifetimes[Zi][j]; //i.e. the longer of the two lifetime limits (corresp. to the lower mass limit)
		*stat = 1;
	}
	//2) If lifetime bin j is NEITHER the lowest NOR the highest bin to be integrated over:
	else if (j != t_lower_lifetime && j != t_upper_lifetime && lifetimes[Zi][j+1] < SNIA_MIN_TIME && lifetimes[Zi][j] >= SNIA_MIN_TIME) {
		*lifetimes_lower = SNIA_MIN_TIME;
		*lifetimes_upper = lifetimes[Zi][j];
		*stat = 2;
	}
	//3) Special case where SNIA_MIN_TIME is greater than lower edge of middle lifetime bin (i.e. lifetimes[Zi][j]):
	else if (j != t_lower_lifetime && j != t_upper_lifetime && lifetimes[Zi][j+1] < SNIA_MIN_TIME && lifetimes[Zi][j] < SNIA_MIN_TIME) {
		*lifetimes_lower = 0.0;
		*lifetimes_upper = 0.0;
		*stat = 3;
	}
	//4) If lifetime bin j is BOTH the lowest AND the highest bin to be integrated over:
	else if (j == t_lower_lifetime && j == t_upper_lifetime && t_lower >= SNIA_MIN_TIME && t_upper <= SNIA_MAX_TIME) {
		*lifetimes_lower = t_lower;
		*lifetimes_upper = t_upper;
		*stat = 4;
	}
	//5) If lifetime bin j is BOTH the lowest AND the highest bin to be integrated over, AND t_lower is less than SNIA_MIN_TIME:
	else if (j == t_lower_lifetime && j == t_upper_lifetime && t_lower < SNIA_MIN_TIME) {
		*lifetimes_lower = SNIA_MIN_TIME;
		*lifetimes_upper = t_upper;
		*stat = 5;
	}
	//6) If lifetime bin j is BOTH the lowest AND the highest bin to be integrated over, AND t_upper is greater than SNIA_MAX_TIME:
	else if (j == t_lower_lifetime && j == t_upper_lifetime && t_upper > SNIA_MAX_TIME) {
		*lifetimes_lower = t_lower;
		*lifetimes_upper = SNIA_MAX_TIME;
		*stat = 6;
	}
	//7) If lifetime bin j IS the lowest bin, but IS NOT the highest bin to be integrated over:
	else if (j == t_lower_lifetime && j != t_upper_lifetime && t_lower >= SNIA_MIN_TIME) {
		*lifetimes_lower = t_lower;
		*lifetimes_upper = lifetimes[Zi][j];
		*stat = 7;
	}
	//8) If lifetime bin j IS the lowest bin, but IS NOT the highest bin to be integrated over, and SNIA_MIN_TIME is less than/equal to lower edge of middle lifetime bin (i.e. lifetimes[Zi][j]):
	else if (j == t_lower_lifetime && j != t_upper_lifetime && t_lower < SNIA_MIN_TIME && lifetimes[Zi][j] >= SNIA_MIN_TIME) {
		*lifetimes_lower = SNIA_MIN_TIME;
		*lifetimes_upper = lifetimes[Zi][j];
		*stat = 8;
	}
	//9) Special case where SNIA_MIN_TIME is greater than lower edge of last lifetime bin (i.e. lifetimes[Zi][j]):
	else if (j == t_lower_lifetime && j != t_upper_lifetime && t_lower < SNIA_MIN_TIME && lifetimes[Zi][j] < SNIA_MIN_TIME) {
		*lifetimes_lower = 0.0;
		*lifetimes_upper = 0.0;
		*stat = 9;
	}
	//10) If lifetime bin j IS NOT the lowest bin, but IS the highest bin to be integrated over:
	else if (j != t_lower_lifetime && j == t_upper_lifetime && t_upper <= SNIA_MAX_TIME && lifetimes[Zi][j+1] >= SNIA_MIN_TIME) {
		*lifetimes_lower = lifetimes[Zi][j+1];
		*lifetimes_upper = t_upper;
		*stat = 10;
	}
	//11) If lifetime bin j IS NOT the lowest bin, but IS the highest bin to be integrated over, and SNIA_MIN_TIME is greater than upper edge of middle lifetime bin (i.e. lifetimes[Zi][j+1])::
	else if (j != t_lower_lifetime && j == t_upper_lifetime && t_upper <= SNIA_MAX_TIME && lifetimes[Zi][j+1] < SNIA_MIN_TIME) {
		*lifetimes_lower = SNIA_MIN_TIME;
		*lifetimes_upper = t_upper;
		*stat = 11;
	}
	//12) If lifetime bin j IS NOT the lowest bin, but IS the highest bin to be integrated over, and t_upper is above SNIA_MAX_TIME:
	else if (j != t_lower_lifetime && j == t_upper_lifetime && t_upper > SNIA_MAX_TIME) {
		*lifetimes_lower = lifetimes[Zi][j+1];
		*lifetimes_upper = SNIA_MAX_TIME;
		*stat = 12;
	}
	else printf("yield_integrals.c: calc_time_limits_SNIa(): No loop entered!\n");
}

#endif //DETAILED_ENRICHMENT

