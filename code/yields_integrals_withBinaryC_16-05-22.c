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
 * 05-05-22: Modified to account for binary_c yields, which need to be integrated over time for AGBs, SNe-II, and SNe-Ia.
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

#ifdef WRITE_YIELD_DATA
  for(ii=0;ii<STEPS*LastDarkMatterSnapShot;ii++) {
	  for(kk=0;kk<LIFETIME_Z_NUM;kk++) {
		  TotNormSNIIEjMass[ii][kk] = 0.0;
		  TotNormSNIaEjMass[ii][kk] = 0.0;
		  TotNormAGBEjMass[ii][kk] = 0.0;
		  TotNormSNIINum[ii][kk] = 0.0;
		  TotNormSNIaNum[ii][kk] = 0.0;
		  TotNormAGBNum[ii][kk] = 0.0;
		  for(ll=0;ll<NUM_ELEMENTS;ll++) {
			  NormSNIIYieldRate_burst[ii][kk][ll] = 0.0;
			  NormSNIaYieldRate_burst[ii][kk][ll] = 0.0;
			  NormAGBYieldRate_burst[ii][kk][ll] = 0.0;
		  }
	  }
  }
#endif //WRITE_YIELD_DATA
}

void integrate_yields()
{
  double previoustime, newtime, deltaT;
  int snap, step, i, mb, Zi, j, kk;
  double timet;
  int mbmax;
  double dt, t_lower, t_upper;

#ifndef BINARYC
  int Mi_lower, Mi_upper, t_lower_lifetime, t_upper_lifetime;
  double Mi_lower_actual, Mi_upper_actual;
  int Mi_lower_lt_SNII, Mi_upper_lt_SNII, Mi_lower_SNII, Mi_upper_SNII, Zi_SNII; //, Zi_correc;
  double Mi_lower_actual_SNII, Mi_upper_actual_SNII, NormFactor;
  int Mi_lower_lt_AGB, Mi_upper_lt_AGB, Mi_lower_AGB, Mi_upper_AGB, Zi_AGB;
  double Mi_lower_actual_AGB, Mi_upper_actual_AGB;
  double DTD_lower, DTD_upper, lifetimes_lower, lifetimes_upper;
#endif //BINARYC

  double SNIIEjectedMasses_lower_actual, SNIIEjectedMasses_upper_actual, SNIITotalMetals_lower_actual, SNIITotalMetals_upper_actual;
  double AGBEjectedMasses_lower_actual, AGBEjectedMasses_upper_actual, AGBTotalMetals_lower_actual, AGBTotalMetals_upper_actual;
#ifdef INDIVIDUAL_ELEMENTS
  double SNIIYields_lower_actual[NUM_ELEMENTS], SNIIYields_upper_actual[NUM_ELEMENTS];
  double AGBYields_lower_actual[NUM_ELEMENTS], AGBYields_upper_actual[NUM_ELEMENTS];
#endif
  double Tot_SNII, Tot_SNIa, Tot_SNe, Tot_NormSNIINum; //Tot_SNIIRate, Tot_SNIaRate, Tot_SNe, Tot_SNII_SP, lifetime_lower_actual, lifetime_upper_actual,
  double EjectedMasses_lower, EjectedMasses_upper, TotalMetals_lower, TotalMetals_upper, Yields_lower[NUM_ELEMENTS], Yields_upper[NUM_ELEMENTS];
  double IL, IU;

#ifdef BINARYC
  int ti_lower, ti_upper, Zi_BC;
  double ti_lower_actual, ti_upper_actual;
  double SNIaEjectedMasses_lower_actual, SNIaEjectedMasses_upper_actual, SNIaTotalMetals_lower_actual, SNIaTotalMetals_upper_actual;
  double SNIIRates_lower_actual, SNIIRates_upper_actual;
  double SNIaRates_lower_actual, SNIaRates_upper_actual;
  double AGBRates_lower_actual, AGBRates_upper_actual;
  double Rates_lower, Rates_upper;
#ifdef INDIVIDUAL_ELEMENTS
  double SNIaYields_lower_actual[NUM_ELEMENTS], SNIaYields_upper_actual[NUM_ELEMENTS];
#endif //INDIVIDUAL_ELEMENTS
#endif //BINARYC

  Tot_SNII = 0.0;
  Tot_SNIa = 0.0;
  Tot_SNe = 0.0;
  Tot_NormSNIINum = 0.0;
  IL = 0.0;
  IU = 0.0;

#ifdef WRITE_YIELD_DATA
  double Solar_mass_ratios[NUM_ELEMENTS];
  int ee2=0.;
  for (int ee=0;ee<11;ee++) {
#ifndef MAINELEMENTS
	  Solar_mass_ratios[ee] = solarAbundMassRatios[ee]; //Taken from Asplund+09 for [H][He][C][N][O][Ne][Mg][Si][S][Ca][Fe]
#else
	  if (ee==0 || ee==1 | ee==4 | ee==6 | ee==10) {
		  Solar_mass_ratios[ee2] = solarAbundMassRatios[ee];
		  ee2++;
	  }
#endif
  }
#endif //WRITE_YIELD_DATA

  double First_SFH_bin_width, Tot_NormSNIIMetalEjecRate, Tot_NormSNIIMassEjecRate, Tot_SNII_MetEjecMass; //, Tot_SFH;
  First_SFH_bin_width = 0.0;
  Tot_NormSNIIMetalEjecRate = 0.0;
  Tot_NormSNIIMassEjecRate = 0.0;
  Tot_SNII_MetEjecMass = 0.0;
  int Zi_pick, i_pick;
  Zi_pick = 3; //3 //Choose which of the discrete lifetime metallicities to assume when calculating the SN-II yields: [0.0004, 0.004, 0.008, 0.02, 0.05, 1]
  i_pick = 0; //13; //5; //Choose which SFH bin you want to print out SNII numbers in timestep "step" for. (22-05-20)

  /*int counta;
  TheSFH[0] = 1.0; ///(tau_dt[0]*UnitTime_in_years/Hubble_h);
  for(counta=1;counta<SFH_NBIN;counta++)
    {
      TheSFH[counta] = 0.0; ///(tau_dt[counta]*UnitTime_in_years/Hubble_h);
    }*/

#ifndef BINARYC
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
#endif //BINARYC

  //***** LOOP OVER SNAPSHOTS *****
  for(snap=0;snap<LastDarkMatterSnapShot;snap++) {
      previoustime = NumToTime(snap); //Time to z=0 from start of current snapshot [in code units]
      newtime = NumToTime(snap+1); //Time to z=0 from end of current snapshot [in code units]
      deltaT = previoustime - newtime; //Width of current snapshot [in code units]

      //***** LOOP OVER TIMESTEPS *****
      for(step=0;step<STEPS;step++) {
		  dt = deltaT/STEPS;  //Time-width of a timestep in current snapshot [in code units]
		  //if (snap < 5) printf("snap %i | step %i | dt = %f Myr\n", snap, step, dt*(UnitTime_in_years/Hubble_h)/1.e6);
		  //printf("%f,", dt*(UnitTime_in_years/Hubble_h)/1.e6);
		  timet = previoustime - ((step + 0.5) * dt); //Time from middle of the current timestep to z=0 [in code units]
		  //if (snap >= 60) printf("snap %i | step %i | t_cosmic = %f Myr\n", snap, step, (13780.169634+0.439881)-(timet*(UnitTime_in_years/Hubble_h)/1.e6));
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
				  //printf("snap=%i step=%i i=%i mb=%i | mbmax = %i (%f) | dt = %e [yr] | mbwidth = %e [yr] | t_lower = %e [yr] | t_upper = %e [yr]\n", snap, step, i, mb, mbmax, SFH_dt[snap][step][i]/dt, dt*UnitTime_in_years/Hubble_h, (SFH_dt[snap][step][i]/mbmax)*UnitTime_in_years/Hubble_h, t_lower, t_upper);

				  //***** LOOP OVER INITIAL METALLICITIES *****
#ifdef BINARYC
				  for (Zi=0;Zi<BC_Z_NUM;Zi++) {
					  //printf("\nsnap %i | step %i | SFH bin %i | minibin %i | Zi %i\n", snap, step, i, mb, Zi);
					  ti_lower_actual = t_lower; //Renaming the actual lower time [in YEARS] to something analogous to Mi_lower_actual for other yield sets.
					  ti_upper_actual = t_upper; //Renaming the actual higher time [in YEARS] to something analogous to Mi_upper_actual for other yield sets.
					  //N.B. What do we do about the negative t_lower for the lowest-z minibin? Maybe set: if(ti_lower_actual < 0.0) ti_lower_actual = 0.0;? Or just leave it, as the returned yield will be 0.0 for t<0.0 anyway? (06-05-22)
					  //Find lower/upper time bins:
					  //ti_lower = find_bc_timebin(ti_upper_actual/1.e6); //N.B. the larger lookback time (i.e. t_upper) corresponds to the shorter cosmic time in the binaryc tables (i.e. ti_lower).
					  ti_lower = find_bc_timebin(ti_lower_actual); //The lookback time, ti_lower_actual, corresponds to the cosmic time, ti_lower, from the birth of the binary_c SP.
					  ti_upper = find_bc_timebin(ti_upper_actual);
					  //printf("find_bc_timebin(): ti_lower_actual = %f [Myr] | ti_upper_actual = %f [Myr] | ti_lower = %i | ti_upper = %i\n"
					  //		 "bcTimes[ti_lower] = %f [Myr] | bcTimes[ti_upper] = %f [Myr]\n", ti_lower_actual/1.e6, ti_upper_actual/1.e6, ti_lower, ti_upper, bcTimes[ti_lower]/1.e6, bcTimes[ti_upper]/1.e6);

					  //Find binary_c metallicity bin:
					  Zi_BC = find_initial_metallicity_comp(Zi, 5);
					  //printf("find_initial_metallicity_comp(): Zi = %i | Zi_BC = %i | lifetimeMetallicities[Zi] = %f | bcMetallicities[Zi_BC] = %f\n", Zi, Zi_BC, lifetimeMetallicities[Zi], bcMetallicities[Zi_BC]);
					  //********
					  //SNe-II:
					  //********
					  //Calc EjectedMass/TotalMetal/yield/rate at lower/upper times:
#ifdef INDIVIDUAL_ELEMENTS
					  calc_ejecta_limits(ti_lower_actual, ti_upper_actual, ti_lower, ti_upper, Zi_BC,
										 BC_TIME_NUM, BC_Z_NUM, bcTimes, bcSNIIEjectedMasses, bcSNIITotalMetals, bcSNIIYields, bcSNIIRates,
										 &SNIIEjectedMasses_lower_actual, &SNIIEjectedMasses_upper_actual, &SNIITotalMetals_lower_actual, &SNIITotalMetals_upper_actual,
										 SNIIYields_lower_actual, SNIIYields_upper_actual, &SNIIRates_lower_actual, &SNIIRates_upper_actual);
#else //INDIVIDUAL_ELEMENTS
					  calc_ejecta_limits(ti_lower_actual, ti_upper_actual, ti_lower, ti_upper, Zi_BC,
								 	 	 BC_TIME_NUM, BC_Z_NUM, bcTimes, bcSNIIEjectedMasses, bcSNIITotalMetals, bcSNIIRates,
										 &SNIIEjectedMasses_lower_actual, &SNIIEjectedMasses_upper_actual, &SNIITotalMetals_lower_actual, &SNIITotalMetals_upper_actual,
										 &SNIIRates_lower_actual, &SNIIRates_upper_actual);
#endif //INDIVIDUAL_ELEMENTS
					  //printf("SNeII: calc_ejecta_limits(): SNIIEjectedMasses_lower_actual = %e | SNIIEjectedMasses_upper_actual = %e\n", SNIIEjectedMasses_lower_actual, SNIIEjectedMasses_upper_actual);
					  for (j=ti_lower;j<=ti_upper;j++) { //loop over binary_c timebins of relevance
						  //Calc time integral limits for each time bin:
#ifdef INDIVIDUAL_ELEMENTS
						  calc_integral_limits(j, Zi_BC, ti_lower, ti_upper, ti_lower_actual, ti_upper_actual,
										   	   BC_TIME_NUM, BC_Z_NUM, bcTimes, bcSNIIEjectedMasses, bcSNIITotalMetals, bcSNIIYields, bcSNIIRates,
											   SNIIEjectedMasses_lower_actual, SNIIEjectedMasses_upper_actual, SNIITotalMetals_lower_actual, SNIITotalMetals_upper_actual,
											   SNIIYields_lower_actual, SNIIYields_upper_actual, SNIIRates_lower_actual, SNIIRates_upper_actual,
											   &IL, &IU, &EjectedMasses_lower, &EjectedMasses_upper, &TotalMetals_lower, &TotalMetals_upper,
											   Yields_lower, Yields_upper, &Rates_lower, &Rates_upper);
#else //INDIVIDUAL_ELEMENTS
						  calc_integral_limits(j, Zi_BC, ti_lower, ti_upper, ti_lower_actual, ti_upper_actual,
								  	  	  	   BC_TIME_NUM, BC_Z_NUM, bcTimes, bcSNIIEjectedMasses, bcSNIITotalMetals, bcSNIIRates,
											   SNIIEjectedMasses_lower_actual, SNIIEjectedMasses_upper_actual, SNIITotalMetals_lower_actual, SNIITotalMetals_upper_actual,
											   SNIIRates_lower_actual, SNIIRates_upper_actual,
											   &IL, &IU, &EjectedMasses_lower, &EjectedMasses_upper, &TotalMetals_lower, &TotalMetals_upper, &Rates_lower, &Rates_upper);
#endif //INDIVIDUAL_ELEMENTS
						  //printf("SNeII: calc_integral_limits(): Bin %i | IL = %e | IU = %e | EjectedMasses_lower = %e | EjectedMasses_upper = %e\n", j, IL, IU, EjectedMasses_lower, EjectedMasses_upper);
						  /*if (mb == mbmax && IU == Mi_upper_actual && IU < Mi_lower_lastTS[i][Zi][0])
							  IU = Mi_lower_lastTS[i][Zi][0]; //only for last mass bin integrated over
						  if (IU > Mi_lower_lastTS[i][Zi][0] && Mi_lower_lastTS[i][Zi][0] > 0.0)
							  IU = Mi_lower_lastTS[i][Zi][0];
						  if (IU < IL)
							  IL = IU;*/

						  //Calc yield ejected in timestep step from stars born in minibin mb of SFH bin i:
						  //NormSNIINum[(STEPS*snap)+step][i][Zi] += 0.0;  //NOTE: Set to 0.0 for now. Should get from binary_c and add as read-in table. (06-05-22) //(IU-IL) * (Chabrier_IMF(IL) + Chabrier_IMF(IU))/2.0; //Number of SNe-II exploding in this timestep from minibin mb of SFHBin i [units: # / Msun]
						  NormSNIINum[(STEPS*snap)+step][i][Zi] += (IU-IL) * ((Rates_lower + Rates_upper)/2.0); //This quantity is in Num/Msun, so should be multiplied by SFR (i.e. M_starsFormedInSFHBin*dt_SFH, in model_yields.c) to get units of [Num/yr]
						  NormSNIIMassEjecRate[(STEPS*snap)+step][i][Zi] += (IU-IL) * ((EjectedMasses_lower + EjectedMasses_upper)/2.0); //This is a unitless quantity, which must be multiplied by SFR (i.e. M_starsFormedInSFHBin*dt_SFH, in model_yields.c) to get units of [Msun/yr], and then by dt (i.e. timestep width, also in model_yields.c) to get units of [Msun]. (22-05-20)
						  NormSNIIMetalEjecRate[(STEPS*snap)+step][i][Zi] += (IU-IL) * ((TotalMetals_lower + TotalMetals_upper)/2.0); //This is a unitless quantity, which must be multiplied by SFR (in model_yields.c) to get units of [Msun/yr], and then by dt (also in model_yields.c) to get units of [Msun]. (22-05-20) //NO! --> //IN [MSun/yr] //The rate of ejection of 'newly synthesised' metals from a 1 Msun stellar population.
						  //printf("SNeII: NormSNIIMassEjecRate[%i][%i][%i] += %e\n", (STEPS*snap)+step, i, Zi, NormSNIIMassEjecRate[(STEPS*snap)+step][i][Zi]);
#ifdef INDIVIDUAL_ELEMENTS
						  for (kk=0;kk<NUM_ELEMENTS;kk++) {
#ifndef MAINELEMENTS
							  NormSNIIYieldRate[(STEPS*snap)+step][i][Zi][kk] += (IU-IL) * ((Yields_lower[kk] + Yields_upper[kk])/2.0); //IN [Msun/yr] //The rate of ejection of individual elements from a 1 Msun stellar population.
#else //MAINELEMENTS
							  switch(kk){case 0: kk=0; break; case 1: kk=1; break; case 2: kk=4; break; case 3: kk=6; break; case 4: kk=10; break;}
							  NormSNIIYieldRate[(STEPS*snap)+step][i][Zi][kk] += (IU-IL) * ((Yields_lower[kk] + Yields_upper[kk])/2.0); //IN [Msun/yr] //The rate of ejection of individual elements from a 1 Msun stellar population.
#endif //MAINELEMENTS
						  }
#endif //INDIVIDUAL_ELEMENTS
					  }
					  //********
					  //SNe-Ia:
					  //********
					  //Calc yield at lower/upper times:
#ifdef INDIVIDUAL_ELEMENTS
					  calc_ejecta_limits(ti_lower_actual, ti_upper_actual, ti_lower, ti_upper, Zi_BC,
										 BC_TIME_NUM, BC_Z_NUM, bcTimes, bcSNIaEjectedMasses, bcSNIaTotalMetals, bcSNIaYields, bcSNIaRates,
										 &SNIaEjectedMasses_lower_actual, &SNIaEjectedMasses_upper_actual, &SNIaTotalMetals_lower_actual, &SNIaTotalMetals_upper_actual,
										 SNIaYields_lower_actual, SNIaYields_upper_actual, &SNIaRates_lower_actual, &SNIaRates_upper_actual);
#else //INDIVIDUAL_ELEMENTS
					  calc_ejecta_limits(ti_lower_actual, ti_upper_actual, ti_lower, ti_upper, Zi_BC,
								 	 	 BC_TIME_NUM, BC_Z_NUM, bcTimes, bcSNIaEjectedMasses, bcSNIaTotalMetals, bcSNIaRates,
										 &SNIaEjectedMasses_lower_actual, &SNIaEjectedMasses_upper_actual, &SNIaTotalMetals_lower_actual, &SNIaTotalMetals_upper_actual,
										 &SNIaRates_lower_actual, &SNIaRates_upper_actual);
#endif //INDIVIDUAL_ELEMENTS
					  //printf("SNeIa: calc_ejecta_limits(): SNIaEjectedMasses_lower_actual = %e | SNIaEjectedMasses_upper_actual = %e\n", SNIaEjectedMasses_lower_actual, SNIaEjectedMasses_upper_actual);
					  for (j=ti_lower;j<=ti_upper;j++) { //loop over binary_c timebins of relevance
						  //Calc time integral limits for each time bin:
#ifdef INDIVIDUAL_ELEMENTS
						  calc_integral_limits(j, Zi_BC, ti_lower, ti_upper, ti_lower_actual, ti_upper_actual,
											   BC_TIME_NUM, BC_Z_NUM, bcTimes, bcSNIaEjectedMasses, bcSNIaTotalMetals, bcSNIaYields, bcSNIaRates,
											   SNIaEjectedMasses_lower_actual, SNIaEjectedMasses_upper_actual, SNIaTotalMetals_lower_actual, SNIaTotalMetals_upper_actual,
											   SNIaYields_lower_actual, SNIaYields_upper_actual, SNIaRates_lower_actual, SNIaRates_upper_actual,
											   &IL, &IU, &EjectedMasses_lower, &EjectedMasses_upper, &TotalMetals_lower, &TotalMetals_upper,
											   Yields_lower, Yields_upper, &Rates_lower, &Rates_upper);
#else //INDIVIDUAL_ELEMENTS
						  calc_integral_limits(j, Zi_BC, ti_lower, ti_upper, ti_lower_actual, ti_upper_actual,
											   BC_TIME_NUM, BC_Z_NUM, bcTimes, bcSNIaEjectedMasses, bcSNIaTotalMetals, bcSNIaRates,
											   SNIaEjectedMasses_lower_actual, SNIaEjectedMasses_upper_actual, SNIaTotalMetals_lower_actual, SNIaTotalMetals_upper_actual,
											   SNIaRates_lower_actual, SNIaRates_upper_actual,
											   &IL, &IU, &EjectedMasses_lower, &EjectedMasses_upper, &TotalMetals_lower, &TotalMetals_upper, &Rates_lower, &Rates_upper);
#endif //INDIVIDUAL_ELEMENTS
						  //printf("SNeIa: calc_integral_limits(): Bin %i | IL = %e | IU = %e | EjectedMasses_lower = %e | EjectedMasses_upper = %e\n", j, IL, IU, EjectedMasses_lower, EjectedMasses_upper);
						  //Calc yield ejected in timestep step from stars born in minibin mb of SFH bin i:
						  NormSNIaNum[(STEPS*snap)+step][i][Zi] += (IU-IL) * ((Rates_lower + Rates_upper)/2.0); //This quantity is in Num/Msun, so should be multiplied by SFR (i.e. M_starsFormedInSFHBin*dt_SFH, in model_yields.c) to get units of [Num/yr]
						  NormSNIaMassEjecRate[(STEPS*snap)+step][i][Zi] += (IU-IL) * ((EjectedMasses_lower + EjectedMasses_upper)/2.0); //This is a unitless quantity, which must be multiplied by SFR (i.e. M_starsFormedInSFHBin*dt_SFH, in model_yields.c) to get units of [Msun/yr], and then by dt (i.e. timestep width, also in model_yields.c) to get units of [Msun]. (22-05-20)
						  NormSNIaMetalEjecRate[(STEPS*snap)+step][i][Zi] += (IU-IL) * ((TotalMetals_lower + TotalMetals_upper)/2.0); //This is a unitless quantity, which must be multiplied by SFR (in model_yields.c) to get units of [Msun/yr], and then by dt (also in model_yields.c) to get units of [Msun]. (22-05-20) //NO! --> //IN [MSun/yr] //The rate of ejection of 'newly synthesised' metals from a 1 Msun stellar population.
						  //printf("SNeIa: NormSNIaMassEjecRate[%i][%i][%i] += %e\n", (STEPS*snap)+step, i, Zi, NormSNIaMassEjecRate[(STEPS*snap)+step][i][Zi]);
#ifdef INDIVIDUAL_ELEMENTS
						  for (kk=0;kk<NUM_ELEMENTS;kk++) {
#ifndef MAINELEMENTS
							  NormSNIaYieldRate[(STEPS*snap)+step][i][Zi][kk] += (IU-IL) * ((Yields_lower[kk] + Yields_upper[kk])/2.0); //IN [Msun/yr] //The rate of ejection of individual elements from a 1 Msun stellar population.
#else //MAINELEMENTS
							  switch(kk){case 0: kk=0; break; case 1: kk=1; break; case 2: kk=4; break; case 3: kk=6; break; case 4: kk=10; break;}
							  NormSNIaYieldRate[(STEPS*snap)+step][i][Zi][kk] += (IU-IL) * ((Yields_lower[kk] + Yields_upper[kk])/2.0); //IN [Msun/yr] //The rate of ejection of individual elements from a 1 Msun stellar population.
#endif //MAINELEMENTS
						  }
#endif //INDIVIDUAL_ELEMENTS
					  }
					  //********
					  //AGBs:
					  //********
					  //Calc yield at lower/upper times:
#ifdef INDIVIDUAL_ELEMENTS
					  calc_ejecta_limits(ti_lower_actual, ti_upper_actual, ti_lower, ti_upper, Zi_BC,
										 BC_TIME_NUM, BC_Z_NUM, bcTimes, bcAGBEjectedMasses, bcAGBTotalMetals, bcAGBYields, bcAGBRates,
										 &AGBEjectedMasses_lower_actual, &AGBEjectedMasses_upper_actual, &AGBTotalMetals_lower_actual, &AGBTotalMetals_upper_actual,
										 AGBYields_lower_actual, AGBYields_upper_actual, &AGBRates_lower_actual, &AGBRates_upper_actual);
#else //INDIVIDUAL_ELEMENTS
					  calc_ejecta_limits(ti_lower_actual, ti_upper_actual, ti_lower, ti_upper, Zi_BC,
								 	 	 BC_TIME_NUM, BC_Z_NUM, bcTimes, bcAGBEjectedMasses, bcAGBTotalMetals, bcAGBRates,
										 &AGBEjectedMasses_lower_actual, &AGBEjectedMasses_upper_actual, &AGBTotalMetals_lower_actual, &AGBTotalMetals_upper_actual,
										 &AGBRates_lower_actual, &AGBRates_upper_actual);
#endif //INDIVIDUAL_ELEMENTS
					  //printf("AGBs:  calc_ejecta_limits(): AGBEjectedMasses_lower_actual = %e | AGBEjectedMasses_upper_actual = %e\n", AGBEjectedMasses_lower_actual, AGBEjectedMasses_upper_actual);
					  for (j=ti_lower;j<=ti_upper;j++) { //loop over binary_c timebins of relevance
						  //Calc time integral limits for each time bin:
#ifdef INDIVIDUAL_ELEMENTS
						  calc_integral_limits(j, Zi_BC, ti_lower, ti_upper, ti_lower_actual, ti_upper_actual,
											   BC_TIME_NUM, BC_Z_NUM, bcTimes, bcAGBEjectedMasses, bcAGBTotalMetals, bcAGBYields, bcAGBRates,
											   AGBEjectedMasses_lower_actual, AGBEjectedMasses_upper_actual, AGBTotalMetals_lower_actual, AGBTotalMetals_upper_actual,
											   AGBYields_lower_actual, AGBYields_upper_actual, AGBRates_lower_actual, AGBRates_upper_actual,
											   &IL, &IU, &EjectedMasses_lower, &EjectedMasses_upper, &TotalMetals_lower, &TotalMetals_upper,
											   Yields_lower, Yields_upper, &Rates_lower, &Rates_upper);
#else //INDIVIDUAL_ELEMENTS
						  calc_integral_limits(j, Zi_BC, ti_lower, ti_upper, ti_lower_actual, ti_upper_actual,
											   BC_TIME_NUM, BC_Z_NUM, bcTimes, bcAGBEjectedMasses, bcAGBTotalMetals, bcAGBRates,
											   AGBEjectedMasses_lower_actual, AGBEjectedMasses_upper_actual, AGBTotalMetals_lower_actual, AGBTotalMetals_upper_actual,
											   AGBRates_lower_actual, AGBRates_upper_actual,
											   &IL, &IU, &EjectedMasses_lower, &EjectedMasses_upper, &TotalMetals_lower, &TotalMetals_upper, &Rates_lower, &Rates_upper);
#endif //INDIVIDUAL_ELEMENTS
						  //printf("AGB:  calc_integral_limits(): Bin %i | IL = %e | IU = %e | EjectedMasses_lower = %e | EjectedMasses_upper = %e\n", j, IL, IU, EjectedMasses_lower, EjectedMasses_upper);
						  //Calc yield ejected in timestep step from stars born in minibin mb of SFH bin i:
						  NormAGBNum[(STEPS*snap)+step][i][Zi] += (IU-IL) * ((Rates_lower + Rates_upper)/2.0); //This quantity is in Num/Msun, so should be multiplied by SFR (i.e. M_starsFormedInSFHBin*dt_SFH, in model_yields.c) to get units of [Num/yr]
						  NormAGBMassEjecRate[(STEPS*snap)+step][i][Zi] += (IU-IL) * ((EjectedMasses_lower + EjectedMasses_upper)/2.0); //This is a unitless quantity, which must be multiplied by SFR (i.e. M_starsFormedInSFHBin*dt_SFH, in model_yields.c) to get units of [Msun/yr], and then by dt (i.e. timestep width, also in model_yields.c) to get units of [Msun]. (22-05-20)
						  NormAGBMetalEjecRate[(STEPS*snap)+step][i][Zi] += (IU-IL) * ((TotalMetals_lower + TotalMetals_upper)/2.0); //This is a unitless quantity, which must be multiplied by SFR (in model_yields.c) to get units of [Msun/yr], and then by dt (also in model_yields.c) to get units of [Msun]. (22-05-20) //NO! --> //IN [MSun/yr] //The rate of ejection of 'newly synthesised' metals from a 1 Msun stellar population.
						  //printf("AGB:  NormAGBMassEjecRate[%i][%i][%i] += %e\n", (STEPS*snap)+step, i, Zi, NormAGBMassEjecRate[(STEPS*snap)+step][i][Zi]);
#ifdef INDIVIDUAL_ELEMENTS
						  for (kk=0;kk<NUM_ELEMENTS;kk++) {
#ifndef MAINELEMENTS
							  NormAGBYieldRate[(STEPS*snap)+step][i][Zi][kk] += (IU-IL) * ((Yields_lower[kk] + Yields_upper[kk])/2.0); //IN [Msun/yr] //The rate of ejection of individual elements from a 1 Msun stellar population.
#else //MAINELEMENTS
							  switch(kk){case 0: kk=0; break; case 1: kk=1; break; case 2: kk=4; break; case 3: kk=6; break; case 4: kk=10; break;}
							  NormAGBYieldRate[(STEPS*snap)+step][i][Zi][kk] += (IU-IL) * ((Yields_lower[kk] + Yields_upper[kk])/2.0); //IN [Msun/yr] //The rate of ejection of individual elements from a 1 Msun stellar population.
#endif //MAINELEMENTS
						  }
#endif //INDIVIDUAL_ELEMENTS
					  } //for (j=ti_lower;j<=ti_upper;j++)

#else //BINARYC
					for (Zi=0;Zi<LIFETIME_Z_NUM;Zi++) {
					  //**********
					  //Calculate overall mass limits:
					  //**********
					  Mi_lower = find_initial_mass(t_upper, Zi); //Mass bin (lifetime arrays) corresp. to lowest mass of star to 'die' in current timestep, from SFH bin i.
					  Mi_upper = find_initial_mass(t_lower, Zi); //Mass bin (lifetime arrays) corresp. to highest mass of star to 'die' in current timestep, from SFH bin i.
					  Mi_lower_actual = lifetimeMasses[Mi_lower] + ((lifetimeMasses[Mi_lower+1]-lifetimeMasses[Mi_lower]) * ((t_upper-lifetimes[Zi][Mi_lower])/(lifetimes[Zi][Mi_lower+1]-lifetimes[Zi][Mi_lower]))); //IN MSUN  //Lowest mass of star to 'die' in current timestep from minimbin mb in SFH bin i.
					  Mi_upper_actual = lifetimeMasses[Mi_upper] + ((lifetimeMasses[Mi_upper+1]-lifetimeMasses[Mi_upper]) * ((t_lower-lifetimes[Zi][Mi_upper])/(lifetimes[Zi][Mi_upper+1]-lifetimes[Zi][Mi_upper]))); //IN MSUN  //Highest mass of star to 'die' in current timestep from minibin mb in SFH bin i.

					  if (Mi_lower_actual < AGB_MIN_MASS)
						  Mi_lower_actual = AGB_MIN_MASS; //No stars below AGB_MIN_MASS = 0.85 Msun contribute to chemical enrichment.
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
						  Zi_SNII = find_initial_metallicity_comp(Zi, 2); //Metallicity bin (SNe-II arrays) corresp. to metallicity Zi. The find_initial_metallicity_comp() function automatically accounts for the mismatch in number of metallicities in the lifetime and SNII metallicity arrays (i.e. accounts for the difference between Zi and Zi_correc made below for other variables) (06-04-22)

						  //Check if mass range is within range for SN-II progenitor stars:
						  Mi_lower_lt_SNII = max_Mi_lower(Mi_lower,2); //Mass bin (lifetime arrays) corresp. to lowest mass of star to 'die' in current timestep from SFH bin i.
						  Mi_upper_lt_SNII = min_Mi_upper(Mi_upper,2); //Mass bin (lifetime arrays) corresp. to highest mass of star to 'die' in current timestep from SFH bin i.
						  if (Mi_lower_lt_SNII <= Mi_upper_lt_SNII) {
							  Mi_lower_SNII = find_SNII_mass_bin(lifetimeMasses[Mi_lower_lt_SNII]); //Mass bin (SNe-II arrays) corresp. to lowest mass of star to 'die' in current timestep from SFH bin i.
							  Mi_upper_SNII = find_SNII_mass_bin(lifetimeMasses[Mi_upper_lt_SNII]); //Mass bin (SNe-II arrays) corresp. to highest mass of star to 'die' in current timestep from SFH bin i.
							  //Redefine the lower/upper allowed masses, given the mass range for SNe-II:
							  Mi_upper_actual_SNII = Mi_upper_actual; //Cannot be > SNII_MAX_MASS, due to checks above.
							  if (Mi_lower_actual < SNII_MIN_MASS)
								  Mi_lower_actual_SNII = SNII_MIN_MASS;
							  else
								  Mi_lower_actual_SNII = Mi_lower_actual;
#ifdef INDIVIDUAL_ELEMENTS
							  //Find true yields at the true upper and lower mass limits, given by 'Mi_upper_actual_SNII' and 'Mi_lower_actual_SNII':
							  calc_ejecta_limits(Mi_lower_actual_SNII, Mi_upper_actual_SNII, Mi_lower_SNII, Mi_upper_SNII, Zi_SNII,
									  	  	  	 SNII_MASS_NUM, SNII_Z_NUM, SNIIMasses, SNIIEjectedMasses, SNIITotalMetals, SNIIYields,
												 &SNIIEjectedMasses_lower_actual, &SNIIEjectedMasses_upper_actual, &SNIITotalMetals_lower_actual, &SNIITotalMetals_upper_actual,
												 SNIIYields_lower_actual, SNIIYields_upper_actual);
#else
							  calc_ejecta_limits(Mi_lower_actual_SNII, Mi_upper_actual_SNII, Mi_lower_SNII, Mi_upper_SNII, Zi_SNII,
												 SNII_MASS_NUM, SNII_Z_NUM, SNIIMasses, SNIIEjectedMasses, SNIITotalMetals,
												 &SNIIEjectedMasses_lower_actual, &SNIIEjectedMasses_upper_actual, &SNIITotalMetals_lower_actual, &SNIITotalMetals_upper_actual);
#endif

							  //NUMERICALLY INTEGRATE OVER THE MASS RANGE APPLICABLE FOR SNe-II:
							  int j, kk;
							  for (j=Mi_lower_SNII;j<=Mi_upper_SNII;j++) {
								  //**********
								  //Calculate mass limits over which to numerically integrate for SNe-II:
								  //**********
#ifdef INDIVIDUAL_ELEMENTS
								  calc_integral_limits(j, Zi_SNII, Mi_lower_SNII, Mi_upper_SNII, Mi_lower_actual_SNII, Mi_upper_actual_SNII,
												   	   SNII_MASS_NUM, SNII_Z_NUM, SNIIMasses, SNIIEjectedMasses, SNIITotalMetals, SNIIYields,
													   SNIIEjectedMasses_lower_actual, SNIIEjectedMasses_upper_actual, SNIITotalMetals_lower_actual, SNIITotalMetals_upper_actual, SNIIYields_lower_actual, SNIIYields_upper_actual,
													   &IL, &IU, &EjectedMasses_lower, &EjectedMasses_upper, &TotalMetals_lower, &TotalMetals_upper, &Yields_lower, &Yields_upper);
#else
								  calc_integral_limits(j, Zi_SNII, Mi_lower_SNII, Mi_upper_SNII, Mi_lower_actual_SNII, Mi_upper_actual_SNII,
										  	  	   	   SNII_MASS_NUM, SNII_Z_NUM, SNIIMasses, SNIIEjectedMasses, SNIITotalMetals,
													   SNIIEjectedMasses_lower_actual, SNIIEjectedMasses_upper_actual, SNIITotalMetals_lower_actual, SNIITotalMetals_upper_actual,
													   &IL, &IU, &EjectedMasses_lower, &EjectedMasses_upper, &TotalMetals_lower, &TotalMetals_upper);
#endif //INDIVIDUAL_ELEMENTS
								  if (mb == mbmax && IU == Mi_upper_actual && IU < Mi_lower_lastTS[i][Zi][0])
									  IU = Mi_lower_lastTS[i][Zi][0]; //only for last mass bin integrated over
								  if (IU > Mi_lower_lastTS[i][Zi][0] && Mi_lower_lastTS[i][Zi][0] > 0.0)
									  IU = Mi_lower_lastTS[i][Zi][0];
								  if (IU < IL)
									  IL = IU;
								  //**********
								  //Calculate rates:
								  //**********
								  if (SNIIMasses[j] <= SNIA_MAX_MASS)
									  NormFactor = 1.0-A_FACTOR; //For mass range where both SN-II and SN-Ia progenitors are possible
								  else
									  NormFactor = 1.0; //For mass range where only SN-II progenitors are possible
								  //SNII_Rate[(STEPS*snap)+step][Zi] += NormFactor * (IU-IL) * ((Chabrier_IMF(IL)*TheSFH[i]/First_SFH_bin_width) + (Chabrier_IMF(IU)*TheSFH[i]/First_SFH_bin_width))/2.0;
								  SNII_Rate[(STEPS*snap)+step][Zi] += NormFactor * (IU-IL) * ((Chabrier_IMF(IL) + Chabrier_IMF(IU))*TheSFH[i]/First_SFH_bin_width)/2.0;
								  NormSNIINum[(STEPS*snap)+step][i][Zi] += NormFactor * (IU-IL) * (Chabrier_IMF(IL) + Chabrier_IMF(IU))/2.0; //Number of SNe-II exploding in this timestep from minibin mb of SFHBin i [units: # / Msun]
								  NormSNIIMassEjecRate[(STEPS*snap)+step][i][Zi] += NormFactor * (IU-IL) * ((EjectedMasses_lower + EjectedMasses_upper)/2.0); //This is a unitless quantity, which must be multiplied by SFR (i.e. M_starsFormedInSFHBin*dt_SFH, in model_yields.c) to get units of [Msun/yr], and then by dt (i.e. timestep width, also in model_yields.c) to get units of [Msun]. (22-05-20)
								  NormSNIIMetalEjecRate[(STEPS*snap)+step][i][Zi] += NormFactor * (IU-IL) * ((TotalMetals_lower + TotalMetals_upper)/2.0); //This is a unitless quantity, which must be multiplied by SFR (in model_yields.c) to get units of [Msun/yr], and then by dt (also in model_yields.c) to get units of [Msun]. (22-05-20) //NO! --> //IN [MSun/yr] //The rate of ejection of 'newly synthesised' metals from a 1 Msun stellar population.
#ifdef INDIVIDUAL_ELEMENTS
								  for (kk=0;kk<NUM_ELEMENTS;kk++) {
#ifndef MAINELEMENTS
									  NormSNIIYieldRate[(STEPS*snap)+step][i][Zi][kk] += NormFactor * (IU-IL) * ((Yields_lower[kk] + Yields_upper[kk])/2.0); //IN [Msun/yr] //The rate of ejection of individual elements from a 1 Msun stellar population.
#else
									  switch(kk){case 0: kk=0; break; case 1: kk=1; break; case 2: kk=4; break; case 3: kk=6; break; case 4: kk=10; break;}
									  NormSNIIYieldRate[(STEPS*snap)+step][i][Zi][kk] += NormFactor * (IU-IL) * ((Yields_lower[kk] + Yields_upper[kk])/2.0); //IN [Msun/yr] //The rate of ejection of individual elements from a 1 Msun stellar population.
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
						  Zi_AGB = find_initial_metallicity_comp(Zi, 4);

						  //Check if mass range is within range for AGB stars:
						  Mi_lower_lt_AGB = max_Mi_lower(Mi_lower,4); //Mass bin (lifetime arrays) corresp. to lowest mass of star to 'die' in current timestep from SFH bin i.
						  Mi_upper_lt_AGB = min_Mi_upper(Mi_upper,4); //Mass bin (lifetime arrays) corresp. to highest mass of star to 'die' in current timestep from SFH bin i.
						  //if (Mi_lower_actual > AGB_MAX_MASS) printf("%i %i %i %i %i | Mi_lower=%i | Mi_upper=%i | Mi_lower_lt_AGB=%i | Mi_upper_lt_AGB=%i | Mi_lower_actual=%f | Mi_upper_actual=%f\n", snap, step, i, mb, Zi, Mi_lower, Mi_upper, Mi_lower_lt_AGB, Mi_upper_lt_AGB, Mi_lower_actual, Mi_upper_actual); //if (snap == 2 && step < 10)
						  if (Mi_lower_lt_AGB <= Mi_upper_lt_AGB && Mi_lower_actual <= AGB_MAX_MASS) {
							  Mi_lower_AGB = find_agb_mass_bin(lifetimeMasses[Mi_lower_lt_AGB]); //Mass bin (AGB arrays) corresp. to lowest mass of star to 'die' in current timestep from SFH bin i.
							  Mi_upper_AGB = find_agb_mass_bin(lifetimeMasses[Mi_upper_lt_AGB]); //Mass bin (AGB arrays) corresp. to highest mass of star to 'die' in current timestep from SFH bin i.
							  //Redefine the lower/upper allowed masses, given the mass range for AGBs:
							  Mi_lower_actual_AGB = Mi_lower_actual; //Cannot be < AGB_MIN_MASS, due to checks above.
							  if (Mi_upper_actual > AGB_MAX_MASS)
								  Mi_upper_actual_AGB = AGB_MAX_MASS;
							  else
								  Mi_upper_actual_AGB = Mi_upper_actual;
#ifdef INDIVIDUAL_ELEMENTS
							  //Find true yields at the true upper and lower mass limits, given by 'Mi_upper_actual_AGB' and 'Mi_lower_actual_AGB':
							  calc_ejecta_limits(Mi_lower_actual_AGB, Mi_upper_actual_AGB, Mi_lower_AGB, Mi_upper_AGB, Zi_AGB,
												 AGB_MASS_NUM, AGB_Z_NUM, AGBMasses, AGBEjectedMasses, AGBTotalMetals, AGBYields,
												 &AGBEjectedMasses_lower_actual, &AGBEjectedMasses_upper_actual, &AGBTotalMetals_lower_actual, &AGBTotalMetals_upper_actual,
												 AGBYields_lower_actual, AGBYields_upper_actual);
#else
							  calc_ejecta_limits(Mi_lower_actual_AGB, Mi_upper_actual_AGB, Mi_lower_AGB, Mi_upper_AGB, Zi_AGB,
												 AGB_MASS_NUM, AGB_Z_NUM, AGBMasses, AGBEjectedMasses, AGBTotalMetals,
												 &AGBEjectedMasses_lower_actual, &AGBEjectedMasses_upper_actual, &AGBTotalMetals_lower_actual, &AGBTotalMetals_upper_actual);
#endif
							  //if (Mi_upper_actual >= 7.0) printf("%i %i %i %i %i | Mi_lower_actual = %f | Mi_upper_actual = %f | AGBEjectedMasses_lower_actual = %f | AGBEjectedMasses_upper_actual = %f\n", snap, step, i, mb, Zi, Mi_lower_actual, Mi_upper_actual, AGBEjectedMasses_lower_actual, AGBEjectedMasses_upper_actual);

							  //NUMERICALLY INTEGRATE OVER THE MASS RANGE APPLICABLE FOR SNe-II:
							  int j, kk;
							  for (j=Mi_lower_AGB;j<=Mi_upper_AGB;j++) {
								  //**********
								  //Calculate mass limits over which to numerically integrate for SNe-II:
								  //**********
#ifdef INDIVIDUAL_ELEMENTS
								  calc_integral_limits(j, Zi_AGB, Mi_lower_AGB, Mi_upper_AGB, Mi_lower_actual_AGB, Mi_upper_actual_AGB,
										  	  	   	   AGB_MASS_NUM, AGB_Z_NUM, AGBMasses, AGBEjectedMasses, AGBTotalMetals, AGBYields,
													   AGBEjectedMasses_lower_actual, AGBEjectedMasses_upper_actual, AGBTotalMetals_lower_actual, AGBTotalMetals_upper_actual, AGBYields_lower_actual, AGBYields_upper_actual,
													   &IL, &IU, &EjectedMasses_lower, &EjectedMasses_upper, &TotalMetals_lower, &TotalMetals_upper, &Yields_lower, &Yields_upper);
#else
								  calc_integral_limits(j, Zi_AGB, Mi_lower_AGB, Mi_upper_AGB, Mi_lower_actual_AGB, Mi_upper_actual_AGB,
										  	  	   	   AGB_MASS_NUM, AGB_Z_NUM, AGBMasses, AGBEjectedMasses, AGBTotalMetals,
													   AGBEjectedMasses_lower_actual, AGBEjectedMasses_upper_actual, AGBTotalMetals_lower_actual, AGBTotalMetals_upper_actual,
													   &IL, &IU, &EjectedMasses_lower, &EjectedMasses_upper, &TotalMetals_lower, &TotalMetals_upper);
#endif //INDIVIDUAL_ELEMENTS
								  /*if (mb == mbmax && IU == Mi_upper_actual && IU < Mi_lower_lastTS[i][Zi][0])
									  IU = Mi_lower_lastTS[i][Zi][0]; //only for last mass bin integrated over
								  if (IU > Mi_lower_lastTS[i][Zi][0] && Mi_lower_lastTS[i][Zi][0] > 0.0)
									  IU = Mi_lower_lastTS[i][Zi][0];
								  if (IU < IL)
									  IL = IU;*/
								  //**********
								  //Calculate rates:
								  //**********
								  AGB_Rate[(STEPS*snap)+step][Zi] += (IU-IL) * ((Chabrier_IMF(IL) + Chabrier_IMF(IU))*TheSFH[i]/First_SFH_bin_width)/2.0;
								  NormAGBNum[(STEPS*snap)+step][i][Zi] += (IU-IL) * (Chabrier_IMF(IL) + Chabrier_IMF(IU))/2.0; //Number of SNe-II exploding in this timestep from minibin mb of SFHBin i [units: # / Msun]
								  NormAGBMassEjecRate[(STEPS*snap)+step][i][Zi] += (IU-IL) * ((EjectedMasses_lower + EjectedMasses_upper)/2.0); //This is a unitless quantity, which must be multiplied by SFR (i.e. M_starsFormedInSFHBin*dt_SFH, in model_yields.c) to get units of [Msun/yr], and then by dt (i.e. timestep width, also in model_yields.c) to get units of [Msun]. (22-05-20) //NO! --> //IN [MSun/yr] //The rate of ejection of mass from a 1 Msun stellar population.
								  NormAGBMetalEjecRate[(STEPS*snap)+step][i][Zi] += (IU-IL) * ((TotalMetals_lower + TotalMetals_upper)/2.0); //This is a unitless quantity, which must be multiplied by SFR (in model_yields.c) to get units of [Msun/yr], and then by dt (also in model_yields.c) to get units of [Msun]. (22-05-20) //NO! --> //IN [MSun/yr] //The rate of ejection of 'newly synthesised' metals from a 1 Msun stellar population.
#ifdef INDIVIDUAL_ELEMENTS
								  for (kk=0;kk<NUM_ELEMENTS;kk++) {
#ifndef MAINELEMENTS
									  NormAGBYieldRate[(STEPS*snap)+step][i][Zi][kk] += (IU-IL) * ((Yields_lower[kk] + Yields_upper[kk])/2.0); //IN [Msun/yr] //The rate of ejection of individual elements from a 1 Msun stellar population.
#else
									  switch(kk){case 0: kk=0; break; case 1: kk=1; break; case 2: kk=4; break; case 3: kk=6; break; case 4: kk=10; break;}
									  NormAGBYieldRate[(STEPS*snap)+step][i][Zi][kk] += (IU-IL) * ((Yields_lower[kk] + Yields_upper[kk])/2.0); //IN [Msun/yr] //The rate of ejection of individual elements from a 1 Msun stellar population.
#endif
								  }
#endif //INDIVIDUAL_ELEMENTS
							  } //for (j=Mi_lower_AGB;j<=Mi_upper_AGB;j++)
						  } //if (Mi_lower_AGB <= Mi_upper_AGB)
						  //*****************************************
					  } //if (t_upper >= lifetimes[Zi][Mi_lower+1])
#endif //BINARYC

#ifdef WRITE_YIELD_DATA
				      if (i==0 && mb==1) {
						  for (int ee=0;ee<NUM_ELEMENTS;ee++) {
#if defined(CHIEFFI) || defined(BINARYC)
							  NormSNIIYieldRate_burst[(STEPS*snap)+step][Zi][ee] = NormSNIIYieldRate[(STEPS*snap)+step][i][Zi][ee];
							  NormSNIaYieldRate_burst[(STEPS*snap)+step][Zi][ee] = NormSNIaYieldRate[(STEPS*snap)+step][i][Zi][ee];
							  NormAGBYieldRate_burst[(STEPS*snap)+step][Zi][ee] = NormAGBYieldRate[(STEPS*snap)+step][i][Zi][ee];
#elif defined(PORTINARI)
							  NormSNIIYieldRate_burst[(STEPS*snap)+step][Zi][ee] = NormSNIIYieldRate[(STEPS*snap)+step][i][Zi][ee] + (Solar_mass_ratios[ee] * NormSNIIMassEjecRate[(STEPS*snap)+step][i][Zi]);
							  NormSNIaYieldRate_burst[(STEPS*snap)+step][Zi][ee] = NormSNIaYieldRate[(STEPS*snap)+step][i][Zi][ee];
							  NormAGBYieldRate_burst[(STEPS*snap)+step][Zi][ee] = NormAGBYieldRate[(STEPS*snap)+step][i][Zi][ee] + (Solar_mass_ratios[ee] * NormAGBMassEjecRate[(STEPS*snap)+step][i][Zi]);
#endif
						  }
					  }
#endif //WRITE_YIELD_DATA
				  } //for (Zi=0;Zi<BC_Z_NUM;Zi++) -OR- for (Zi=0;Zi<LIFETIME_Z_NUM;Zi++) //METALLICITIES
			  } //for (mb=1;mb<=mbmax;mb++) //MINI_BINS
#ifdef WRITE_YIELD_DATA
			  //Calculate total mass ejecta rates in each timestep (assuming 1Msun of SF in every SFH bin):
			  for (Zi=0;Zi<LIFETIME_Z_NUM;Zi++) {
				  TotNormSNIIEjMass[(STEPS*snap)+step][Zi] += NormSNIIMassEjecRate[(STEPS*snap)+step][i][Zi];
				  TotNormSNIaEjMass[(STEPS*snap)+step][Zi] += NormSNIaMassEjecRate[(STEPS*snap)+step][i][Zi];
				  TotNormAGBEjMass[(STEPS*snap)+step][Zi] += NormAGBMassEjecRate[(STEPS*snap)+step][i][Zi];
				  TotNormSNIINum[(STEPS*snap)+step][Zi] += NormSNIINum[(STEPS*snap)+step][i][Zi];
				  TotNormSNIaNum[(STEPS*snap)+step][Zi] += NormSNIaNum[(STEPS*snap)+step][i][Zi];
				  TotNormAGBNum[(STEPS*snap)+step][Zi] += NormAGBNum[(STEPS*snap)+step][i][Zi];
			  }
#endif //WRITE_YIELD_DATA
		  } //for (i=0;i<=SFH_ibin[snap][step];i++)
		  Tot_SNII += SNII_Rate[(STEPS*snap)+step][Zi_pick]*(dt*(UnitTime_in_years/Hubble_h));
		  Tot_SNIa += SNIa_Rate[(STEPS*snap)+step][Zi_pick]*(dt*(UnitTime_in_years/Hubble_h));
		  Tot_NormSNIINum += NormSNIINum[(STEPS*snap)+step][i_pick][Zi_pick]; //Summing up all the SNe that exploded in every timestep from all the minibins of SFHBin i_pick at metallicity Zi_pick. //*(dt*(UnitTime_in_years/Hubble_h))/First_SFH_bin_width;
		  Tot_SNe = Tot_SNII+Tot_SNIa;
		  Tot_NormSNIIMetalEjecRate += NormSNIIMetalEjecRate[(STEPS*snap)+step][0][Zi_pick];
		  Tot_NormSNIIMassEjecRate += NormSNIIMassEjecRate[(STEPS*snap)+step][0][Zi_pick];
		  Tot_SNII_MetEjecMass += ((dt*(UnitTime_in_years/Hubble_h))/First_SFH_bin_width) * (NormSNIIMetalEjecRate[(STEPS*snap)+step][0][0] + lifetimeMetallicities[Zi_pick]*NormSNIIMetalEjecRate[(STEPS*snap)+step][0][Zi_pick]);
		  /*if (snap >= 45 && snap < 50) {
			  printf("\nSnap %i | Step %i | Lookback time = %e [Myr]:\n", snap, step, timet*(UnitTime_in_years/Hubble_h)/1.e6);
			  printf("SNeII: NormSNIIMassEjecRate[%i][%i][%i] += %e\n", (STEPS*snap)+step, i_pick, Zi_pick, NormSNIIMassEjecRate[(STEPS*snap)+step][i_pick][Zi_pick]);
			  printf("SNeIa: NormSNIaMassEjecRate[%i][%i][%i] += %e\n", (STEPS*snap)+step, i_pick, Zi_pick, NormSNIaMassEjecRate[(STEPS*snap)+step][i_pick][Zi_pick]);
			  printf("AGB:   NormAGBMassEjecRate[%i][%i][%i]  += %e\n", (STEPS*snap)+step, i_pick, Zi_pick, NormAGBMassEjecRate[(STEPS*snap)+step][i_pick][Zi_pick]);
		  }*/
		  //Print timestep lookback times to z=0 in Myr:
		  //printf("%e,", timet*(UnitTime_in_years/Hubble_h)/1.e6);
		  //Print (unitless) total ejected mass yields for every timestep from SFHbin i_pick at metallicity Zi_pick:
		  //printf("%.4e,", NormSNIIMassEjecRate[(STEPS*snap)+step][i_pick][Zi_pick]);
		  //printf("%.4e,", NormSNIaMassEjecRate[(STEPS*snap)+step][i_pick][Zi_pick]);
		  //printf("%e,", NormAGBMassEjecRate[(STEPS*snap)+step][i_pick][Zi_pick]);
		  //Print total unitless mass ejected from all SFH bins of each timestep:
		  /*double EjRate=0.;
		  for (int iSFH=0; iSFH<=SFH_ibin[snap][step]; iSFH++) {
			  //EjRate += NormSNIIMassEjecRate[(STEPS*snap)+step][iSFH][Zi_pick];
			  //EjRate += NormSNIaMassEjecRate[(STEPS*snap)+step][iSFH][Zi_pick];
			  EjRate += NormAGBMassEjecRate[(STEPS*snap)+step][iSFH][Zi_pick];
		  }
		  printf("%e,", EjRate);
		  printf("%e\n\n", TotNormAGBEjMass[(STEPS*snap)+step][Zi_pick]);*/
      } //for(step=0;step<STEPS;step++)
   } //for(snap=0;snap<LastDarkMatterSnapShot;snap++)
#ifdef PARALLEL
  	if(ThisTask == 0)
#endif
    printf("Yield integrals calculated.\n");

  //Check that final yields don't exceed some pre-defined minimum and maximum:
  int ii, jj, k;
  for(snap=0;snap<LastDarkMatterSnapShot;snap++) //LOOP OVER SNAPSHOTS
    for(step=0;step<STEPS;step++)
      for(jj=0;jj<SFH_NBIN;jj++)
    	  for(k=0;k<LIFETIME_Z_NUM;k++)
	  {

	    ii=(STEPS*snap)+step;

	    /*if(NormAGBMassEjecRate[ii][jj][k]>0.3)
	      {
		printf("ii=%d snap=%d step=%d jj=%d k=%d AGB_Rate=%0.2f\n",ii,snap, step,jj,k, NormAGBMassEjecRate[ii][jj][k]);
		terminate("AGB rate too high");
	      }
	    if(NormSNIaMassEjecRate[ii][jj][k]>0.7)
	      {
		printf("ii=%d jj=%d k=%d SNIa_Rate=%0.2f\n",ii,jj,k, NormSNIaMassEjecRate[ii][jj][k]);
		terminate("SNIa rate too high");
	      }*/

	    if(NormSNIIMassEjecRate[ii][jj][k]<0.)
	      NormSNIIMassEjecRate[ii][jj][k]=0.;
	    if(NormSNIaMassEjecRate[ii][jj][k]<0.)
	      NormSNIaMassEjecRate[ii][jj][k]=0.;
	    if(NormAGBMassEjecRate[ii][jj][k]<0.)
	      NormAGBMassEjecRate[ii][jj][k]=0.;
	    if(NormSNIIMetalEjecRate[ii][jj][k]<0.)
	      NormSNIIMetalEjecRate[ii][jj][k]=0.;
	    if(NormSNIaMetalEjecRate[ii][jj][k]<0.)
	      NormSNIaMetalEjecRate[ii][jj][k]=0.;
	    if(NormAGBMetalEjecRate[ii][jj][k]<0.)
	      NormAGBMetalEjecRate[ii][jj][k]=0.;
	  }

#ifdef WRITE_YIELD_DATA
    //Write yields from a 1Msun burst (for each metallicity) at the start of the simulation to file:
	//SNeII:
	for(kk=0;kk<LIFETIME_Z_NUM;kk++) { //LOOP OVER METALLICITIES
	  //SNeII:
	  FILE *SNIIyf;
	  char filename[100];
	  //sprintf(filename, "./YieldTables/BurstYields/SNII_burst_yields_Z%.4f.txt", lifetimeMetallicities[kk]);
#if defined(BINARYC)
	  char yieldSet[] = "binaryc";
#elif defined(PORTINARI)
	  char yieldSet[] = "P98";
#elif defined(CHIEFFI)
	  char yieldSet[] = "CL01";
#endif
	  sprintf(filename, "./YieldTables/BurstYields/SNII_%s_burst_yields_Z%.4f.txt", yieldSet, lifetimeMetallicities[kk]);
	  SNIIyf = fopen(filename, "w");
	  if (SNIIyf == NULL) {
		  printf("Error opening file: %s!\n", filename);
		  exit(1);
	  }
	  int ee;
	  for(ee=0;ee<NUM_ELEMENTS;ee++) { //LOOP OVER ELEMENTS
		  for(snap=0;snap<LastDarkMatterSnapShot;snap++) { //LOOP OVER SNAPSHOTS
			  for(step=0;step<STEPS;step++) { //LOOP OVER STEPS
					  ii=(STEPS*snap)+step;
					  fprintf(SNIIyf, "%e ", NormSNIIYieldRate_burst[ii][kk][ee]);
			  }
		  }
		  fprintf(SNIIyf, "\n");
	  }
	  fclose(SNIIyf);
	}

	//SNeIa:
	for(kk=0;kk<LIFETIME_Z_NUM;kk++) { //LOOP OVER METALLICITIES
	  FILE *SNIayf;
	  char filename[100];
	  //sprintf(filename, "./YieldTables/BurstYields/SNIa_burst_yields_Z%.4f.txt", lifetimeMetallicities[kk]);
#ifdef BINARYC
	  char yieldSet[] = "binaryc";
#else
	  char yieldSet[] = "T03";
#endif
	  sprintf(filename, "./YieldTables/BurstYields/SNIa_%s_burst_yields_Z%.4f.txt", yieldSet, lifetimeMetallicities[kk]);
	  SNIayf = fopen(filename, "w");
	  if (SNIayf == NULL) {
		  printf("Error opening file: %s!\n", filename);
		  exit(1);
	  }
	  int ee;
	  for(ee=0;ee<NUM_ELEMENTS;ee++) { //LOOP OVER ELEMENTS
		  for(snap=0;snap<LastDarkMatterSnapShot;snap++) { //LOOP OVER SNAPSHOTS
			  for(step=0;step<STEPS;step++) { //LOOP OVER STEPS
					  ii=(STEPS*snap)+step;
					  fprintf(SNIayf, "%e ", NormSNIaYieldRate_burst[ii][kk][ee]);
			  }
		  }
		  fprintf(SNIayf, "\n");
	  }
	  fclose(SNIayf);
	}

	//AGB:
	for(kk=0;kk<LIFETIME_Z_NUM;kk++) { //LOOP OVER METALLICITIES
	  FILE *AGByf;
	  char filename[100];
	  //sprintf(filename, "./YieldTables/BurstYields/AGB_burst_yields_Z%.4f.txt", lifetimeMetallicities[kk]);
#ifdef BINARYC
	  char yieldSet[] = "binaryc";
#else
	  char yieldSet[] = "M01";
#endif
	  sprintf(filename, "./YieldTables/BurstYields/AGB_%s_burst_yields_Z%.4f.txt", yieldSet, lifetimeMetallicities[kk]);
	  AGByf = fopen(filename, "w");
	  if (AGByf == NULL) {
		  printf("Error opening file: %s!\n", filename);
		  exit(1);
	  }
	  int ee;
	  for(ee=0;ee<NUM_ELEMENTS;ee++) { //LOOP OVER ELEMENTS
		  for(snap=0;snap<LastDarkMatterSnapShot;snap++) { //LOOP OVER SNAPSHOTS
			  for(step=0;step<STEPS;step++) { //LOOP OVER STEPS
					  ii=(STEPS*snap)+step;
					  fprintf(AGByf, "%e ", NormAGBYieldRate_burst[ii][kk][ee]);
			  }
		  }
		  fprintf(AGByf, "\n");
	  }
	  fclose(AGByf);
	} //for(kk=0;kk<LIFETIME_Z_NUM;kk++)
	printf("Burst yield files written.\n");

	//--------------------
	//WRITE TIMESTEPS IN MYR:
	FILE *Tsf;
	char filename0[100];
	sprintf(filename0, "./YieldTables/TotalYields/lookbackTimes_Myr.txt");
	Tsf = fopen(filename0, "w");
	if (Tsf == NULL) {
	printf("Error opening file: %s!\n", filename0);
	exit(1);
	}
	for(snap=0;snap<LastDarkMatterSnapShot;snap++) {//LOOP OVER SNAPSHOTS
		previoustime = NumToTime(snap); //Time to z=0 from start of current snapshot [in code units]
		newtime = NumToTime(snap+1); //Time to z=0 from end of current snapshot [in code units]
		deltaT = previoustime - newtime; //Width of current snapshot [in code units]
		for(step=0;step<STEPS;step++) { //LOOP OVER TIMESTEPS
			dt = deltaT/STEPS;  //Time-width of a timestep in current snapshot [in code units]
			timet = previoustime - ((step + 0.5) * dt); //Time from middle of the current timestep to z=0 [in code units]
			fprintf(Tsf, "%e ", timet*(UnitTime_in_years/Hubble_h)/1.e6);
		}
	}
	fclose(Tsf);

	//WRITE TOTAL EJECTED MASSES AND SN NUMBERS:
	const char *channels[NUM_CHANNELS] = {"SNII","SNIa","AGB"};
#if defined(PORTINARI)
	const char *Orig_yieldSets[NUM_CHANNELS] = {"P98","T03","M01"};
#elif defined(CHIEFFI)
	const char *Orig_yieldSets[NUM_CHANNELS] = {"CL04","T03","M01"};
#endif
	for(int ll=0;ll<NUM_CHANNELS;ll++) { //LOOP OVER CHANNELS (SNe-II, SNe-Ia, AGBs)
#ifdef BINARYC
		char yieldSet[] = "binaryc";
#else
		char yieldSet[4];
		strcpy(yieldSet, Orig_yieldSets[ll]);
#endif //BINARYC
		for(kk=0;kk<LIFETIME_Z_NUM;kk++) { //LOOP OVER METALLICITIES
		  FILE *Ejf, *Raf;
		  char filename1[100], filename2[100];
		  sprintf(filename1, "./YieldTables/TotalYields/%s_%s_TotalEjectedMass_Z%.4f.txt", channels[ll], yieldSet, lifetimeMetallicities[kk]);
		  Ejf = fopen(filename1, "w");
		  if (Ejf == NULL) {
			printf("Error opening file: %s!\n", filename1);
			exit(1);
		  }
		  sprintf(filename2, "./YieldTables/TotalYields/%s_%s_TotalSNNum_Z%.4f.txt", channels[ll], yieldSet, lifetimeMetallicities[kk]);
		  Raf = fopen(filename2, "w");
		  if (Raf == NULL) {
			printf("Error opening file: %s!\n", filename2);
			exit(1);
		  }
		  for(snap=0;snap<LastDarkMatterSnapShot;snap++) {//LOOP OVER SNAPSHOTS
			for(step=0;step<STEPS;step++) { //LOOP OVER TIMESTEPS
				if (channels[ll]=="SNII") {
					fprintf(Ejf, "%e ", TotNormSNIIEjMass[(STEPS*snap)+step][kk]);
					fprintf(Raf, "%e ", TotNormSNIINum[(STEPS*snap)+step][kk]);
				}
				else if(channels[ll]=="SNIa") {
					fprintf(Ejf, "%e ", TotNormSNIaEjMass[(STEPS*snap)+step][kk]);
					fprintf(Raf, "%e ", TotNormSNIaNum[(STEPS*snap)+step][kk]);
				}
				else if(channels[ll]=="AGB") {
					fprintf(Ejf, "%e ", TotNormAGBEjMass[(STEPS*snap)+step][kk]);
					fprintf(Raf, "%e ", TotNormAGBNum[(STEPS*snap)+step][kk]);
				}
			}
		  }
		  fclose(Ejf);
		  fclose(Raf);
		} //for(kk=0;kk<LIFETIME_Z_NUM;kk++)
	} //for(int ll=0;ll<NUM_CHANNELS;ll++)


	printf("Total yield files written.\n");
#endif //WRITE_YIELD_DATA

}


int find_initial_metallicity_comp(int Zi, int table_type)
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
#ifndef BINARYC
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
#else //BINARYC
    case 5: //Binary_c metallicity table
		while (Zi_bin == -1) {
			if (bcMetallicities[i] < Z_in) {
			  i++;
			  if (i == BC_Z_NUM) Zi_bin = i-1;
			}
			else Zi_bin = i;
		}
    break;
#endif //BINARYC
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

#ifdef BINARYC
int find_bc_timebin(double actual_time) //actual_time should be in yr for this function
{
  int i=-1;
  if (actual_time <= 0.0) {
	  i=1;
	  //printf("find_bc_timebin(): Loop 1: actual_time = %f | bcTimes[%i] = %f\n", actual_time, i-1, bcTimes[i-1]);
	  return i-1; //If the bin 'touches now', then return the 0th time bin
  }
  else {
	  /*do {
		  i++;
		  bcBinwidth = BC_LOGT_BINWIDTH*bcTimes[i]*log(10.); //Width of binary_c timebin i in (linear) yr
	  }
	  while (bcTimes[i]+(0.5*bcBinwidth) < actual_time);
	  return i-1; //Returns the binary_c dampen in which LGals_time lies, given that bcTimes[] contains the mid-time in each bin*/
	  do {
		  i++;
	  }
	  while (bcTimes[i] < actual_time);
	  //printf("find_bc_timebin(): Loop 2: actual_time = %f [yr] | bcTimes[%i] = %f [yr]\n", actual_time, i-1, bcTimes[i-1]);
	  return i-1; //Returns the binary_c timebin (which are cosmic times) directly BELOW the actual_time given (mimicking what is done in find_initial_mass())
  }

}
#endif //BINARYC

#ifndef BINARYC
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
#endif //BINARYC

#ifdef INDIVIDUAL_ELEMENTS
#ifdef BINARYC
void calc_ejecta_limits(double x_lower_actual, double x_upper_actual, int x_lower, int x_upper, int Zi,
						int X_NUM, int Z_NUM, double x_bins[X_NUM], double EjectedMasses[Z_NUM][X_NUM], double TotalMetals[Z_NUM][X_NUM],
						double Yields[Z_NUM][NUM_ELEMENTS][X_NUM], double Rates[Z_NUM][X_NUM],
						double* EjectedMasses_lower_actual, double* EjectedMasses_upper_actual, double* TotalMetals_lower_actual, double* TotalMetals_upper_actual,
						double* Yields_lower_actual, double* Yields_upper_actual, double* Rates_lower_actual, double* Rates_upper_actual)
#else //BINARYC
void calc_ejecta_limits(double x_lower_actual, double x_upper_actual, int x_lower, int x_upper, int Zi,
						int X_NUM, int Z_NUM, double x_bins[X_NUM], double EjectedMasses[Z_NUM][X_NUM], double TotalMetals[Z_NUM][X_NUM], double Yields[Z_NUM][NUM_ELEMENTS][X_NUM],
						double* EjectedMasses_lower_actual, double* EjectedMasses_upper_actual, double* TotalMetals_lower_actual, double* TotalMetals_upper_actual,
						double* Yields_lower_actual, double* Yields_upper_actual)
#endif //BINARYC
#else //INDIVIDUAL_ELEMENTS
#ifdef BINARYC
void calc_ejecta_limits(double x_lower_actual, double x_upper_actual, int x_lower, int x_upper, int Zi,
						int X_NUM, int Z_NUM, double x_bins[X_NUM], double EjectedMasses[Z_NUM][X_NUM], double TotalMetals[Z_NUM][X_NUM], double Rates[Z_NUM][X_NUM],
						double* EjectedMasses_lower_actual, double* EjectedMasses_upper_actual, double* TotalMetals_lower_actual, double* TotalMetals_upper_actual,
						double* Rates_lower, double* Rates_upper)
#else //BINARYC
void calc_ejecta_limits(double x_lower_actual, double x_upper_actual, int x_lower, int x_upper, int Zi,
						int X_NUM, int Z_NUM, double x_bins[X_NUM], double EjectedMasses[Z_NUM][X_NUM], double TotalMetals[Z_NUM][X_NUM],
						double* EjectedMasses_lower_actual, double* EjectedMasses_upper_actual, double* TotalMetals_lower_actual, double* TotalMetals_upper_actual)
#endif //BINARYC
#endif //INDIVIDUAL_ELEMENTS
{
  int k, kk;
  //Lower limits:
  *EjectedMasses_lower_actual = EjectedMasses[Zi][x_lower] + ((EjectedMasses[Zi][x_lower+1]-EjectedMasses[Zi][x_lower]) * ((x_lower_actual-x_bins[x_lower])/(x_bins[x_lower+1]-x_bins[x_lower])));
//  printf("calc_ejecta_limits(): x_lower = %i | x_lower+1 = %i | x_lower_actual = %f | x_bins[x_lower] = %f | x_bins[x_lower+1] = %f\n",
//		  x_lower, x_lower+1, x_lower_actual, x_bins[x_lower], x_bins[x_lower+1]);
//  printf("calc_ejecta_limits(): %e = %e + ((%e-%e) * ((%f-%f)/(%f-%f)))\n",
//		  *EjectedMasses_lower_actual, EjectedMasses[Zi][x_lower], EjectedMasses[Zi][x_lower+1], EjectedMasses[Zi][x_lower], x_lower_actual, x_bins[x_lower], x_bins[x_lower+1], x_bins[x_lower]);
  *TotalMetals_lower_actual = TotalMetals[Zi][x_lower] + ((TotalMetals[Zi][x_lower+1]-TotalMetals[Zi][x_lower]) * ((x_lower_actual-x_bins[x_lower])/(x_bins[x_lower+1]-x_bins[x_lower])));
#ifdef INDIVIDUAL_ELEMENTS
  for (k=0;k<NUM_ELEMENTS;k++) {
#ifndef MAINELEMENTS
	  kk=k;
	  Yields_lower_actual[k] = Yields[Zi][kk][x_lower] + ((Yields[Zi][kk][x_lower+1]-Yields[Zi][kk][x_lower]) * ((x_lower_actual-x_bins[x_lower])/(x_bins[x_lower+1]-x_bins[x_lower])));
#else
	  switch(k){case 0: kk=0; break; case 1: kk=1; break; case 2: kk=4; break; case 3: kk=6; break; case 4: kk=10; break;}
	  Yields_lower_actual[k] = Yields[Zi][kk][x_lower] + ((Yields[Zi][kk][x_lower+1]-Yields[Zi][kk][x_lower]) * ((x_lower_actual-x_bins[x_lower])/(x_bins[x_lower+1]-x_bins[x_lower])));
#endif
  }
#endif //INDIVIDUAL_ELEMENTS
#ifdef BINARYC
  *Rates_lower_actual = Rates[Zi][x_lower] + ((Rates[Zi][x_lower+1]-Rates[Zi][x_lower]) * ((x_lower_actual-x_bins[x_lower])/(x_bins[x_lower+1]-x_bins[x_lower])));
#endif //BINARYC

  //Upper limits:
  *EjectedMasses_upper_actual = EjectedMasses[Zi][x_upper] + ((EjectedMasses[Zi][x_upper+1]-EjectedMasses[Zi][x_upper]) * ((x_upper_actual-x_bins[x_upper])/(x_bins[x_upper+1]-x_bins[x_upper])));
  *TotalMetals_upper_actual = TotalMetals[Zi][x_upper] + ((TotalMetals[Zi][x_upper+1]-TotalMetals[Zi][x_upper]) * ((x_upper_actual-x_bins[x_upper])/(x_bins[x_upper+1]-x_bins[x_upper])));
#ifdef INDIVIDUAL_ELEMENTS
  for (k=0;k<NUM_ELEMENTS;k++) {
#ifndef MAINELEMENTS
	  kk=k;
	  Yields_upper_actual[k] = Yields[Zi][kk][x_upper] + ((Yields[Zi][kk][x_upper+1]-Yields[Zi][kk][x_upper]) * ((x_upper_actual-x_bins[x_upper])/(x_bins[x_upper+1]-x_bins[x_upper])));
#else
	  switch(k){case 0: kk=0; break; case 1: kk=1; break; case 2: kk=4; break; case 3: kk=6; break; case 4: kk=10; break;}
	  Yields_upper_actual[k] = Yields[Zi][kk][x_upper] + ((Yields[Zi][kk][x_upper+1]-Yields[Zi][kk][x_upper]) * ((x_upper_actual-x_bins[x_upper])/(x_bins[x_upper+1]-x_bins[x_upper])));
#endif
  }
#endif //INDIVIDUAL_ELEMENTS
#ifdef BINARYC
  *Rates_upper_actual = Rates[Zi][x_upper] + ((Rates[Zi][x_upper+1]-Rates[Zi][x_upper]) * ((x_upper_actual-x_bins[x_upper])/(x_bins[x_upper+1]-x_bins[x_upper])));
#endif //BINARYC
}


#ifdef INDIVIDUAL_ELEMENTS
#ifdef BINARYC
void calc_integral_limits(int j, int Zi_chan, int x_lower_chan, int x_upper_chan, double x_lower_actual, double x_upper_actual,
						  int X_NUM, int Z_NUM, double x_bins[X_NUM], double EjectedMasses[Z_NUM][X_NUM], double TotalMetals[Z_NUM][X_NUM],
						  double Yields[Z_NUM][NUM_ELEMENTS][X_NUM], double Rates[Z_NUM][X_NUM],
						  double EjectedMasses_lower_actual, double EjectedMasses_upper_actual, double TotalMetals_lower_actual, double TotalMetals_upper_actual,
						  double Yields_lower_actual[], double Yields_upper_actual[], double Rates_lower_actual, double Rates_upper_actual,
						  double* IL, double* IU, double* EjectedMasses_lower, double* EjectedMasses_upper, double* TotalMetals_lower, double* TotalMetals_upper,
						  double* Yields_lower, double* Yields_upper, double* Rates_lower, double* Rates_upper)
#else //BINARYC
void calc_integral_limits(int j, int Zi_chan, int x_lower_chan, int x_upper_chan, double x_lower_actual, double x_upper_actual,
						  int X_NUM, int Z_NUM, double x_bins[X_NUM], double EjectedMasses[Z_NUM][X_NUM], double TotalMetals[Z_NUM][X_NUM], double Yields[Z_NUM][NUM_ELEMENTS][X_NUM],
						  double EjectedMasses_lower_actual, double EjectedMasses_upper_actual, double TotalMetals_lower_actual, double TotalMetals_upper_actual,
						  double Yields_lower_actual[], double Yields_upper_actual[],
						  double* IL, double* IU, double* EjectedMasses_lower, double* EjectedMasses_upper, double* TotalMetals_lower, double* TotalMetals_upper,
						  double* Yields_lower, double* Yields_upper)
#endif //BINARYC
#else //INDIVIDUAL_ELEMENTS
#ifdef BINARYC
void calc_integral_limits(int j, int Zi_chan, int x_lower_chan, int x_upper_chan, double x_lower_actual, double x_upper_actual,
						  int X_NUM, int Z_NUM, double x_bins[X_NUM], double EjectedMasses[Z_NUM][X_NUM], double TotalMetals[Z_NUM][X_NUM], double Rates[Z_NUM][X_NUM],
						  double EjectedMasses_lower_actual, double EjectedMasses_upper_actual, double TotalMetals_lower_actual, double TotalMetals_upper_actual,
						  double Rates_lower_actual, double Rates_upper_actual,
						  double* IL, double* IU, double* EjectedMasses_lower, double* EjectedMasses_upper, double* TotalMetals_lower, double* TotalMetals_upper,
						  double* Rates_lower, double* Rates_upper);
#else //BINARYC
void calc_integral_limits(int j, int Zi_chan, int x_lower_chan, int x_upper_chan, double x_lower_actual, double x_upper_actual,
						  int X_NUM, int Z_NUM, double x_bins[X_NUM], double EjectedMasses[Z_NUM][X_NUM], double TotalMetals[Z_NUM][X_NUM],
						  double EjectedMasses_lower_actual, double EjectedMasses_upper_actual, double TotalMetals_lower_actual, double TotalMetals_upper_actual,
						  double* IL, double* IU, double* EjectedMasses_lower, double* EjectedMasses_upper, double* TotalMetals_lower, double* TotalMetals_upper)
#endif //BINARYC
#endif //INDIVIDUAL_ELEMENTS
{
	int kk;
	//1) If mass bin j is NEITHER the lowest NOR highest mass bin to be integrated over:
	if (j != x_lower_chan && j != x_upper_chan) {
		*IL = x_bins[j];
		*IU = x_bins[j+1];
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
#ifdef BINARYC
		*Rates_lower = Rates[Zi_chan][j];
		*Rates_upper = Rates[Zi_chan][j+1];
#endif //BINARYC
	}
	//2) If mass bin j is BOTH the lowest AND highest mass bin to be integrated over:
	else if (j == x_lower_chan && j == x_upper_chan) {
		*IL = x_lower_actual;
		*IU = x_upper_actual;
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
#ifdef BINARYC
		*Rates_lower = Rates_lower_actual;
		*Rates_upper = Rates_upper_actual;
#endif //BINARYC
	}
	//3) If mass bin j IS the lowest but IS NOT the highest mass bin to be integrated over:
	else if (j == x_lower_chan && j != x_upper_chan) {
		*IL = x_lower_actual;
		*IU = x_bins[j+1];
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
#ifdef BINARYC
		*Rates_lower = Rates_lower_actual;
		*Rates_upper = Rates[Zi_chan][j+1];
#endif //BINARYC
	}
	//4) If mass bin j IS NOT the lowest but IS the highest mass bin to be integrated over:
	else if (j != x_lower_chan && j == x_upper_chan) {
		*IL = x_bins[j];
		*IU = x_upper_actual;
		*EjectedMasses_lower = EjectedMasses[Zi_chan][j];
		*EjectedMasses_upper = EjectedMasses_upper_actual;
		*TotalMetals_lower = TotalMetals[Zi_chan][j];
		*TotalMetals_upper = TotalMetals_upper_actual;
#ifdef INDIVIDUAL_ELEMENTS
		for (kk=0; kk<NUM_ELEMENTS; kk++) {
			Yields_lower[kk] = Yields[Zi_chan][kk][j];
			Yields_upper[kk] = Yields_upper_actual[kk];
		}
#endif //INDIVIDUAL_ELEMENTS#
#ifdef BINARYC
		*Rates_lower = Rates[Zi_chan][j];
		*Rates_upper = Rates_upper_actual;
#endif //BINARYC
	}
	else printf("yield_integrals.c: calc_integral_limits(): No option entered!\n");
}


#ifndef BINARYC
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
#endif //BINARYC

#endif //DETAILED_ENRICHMENT

