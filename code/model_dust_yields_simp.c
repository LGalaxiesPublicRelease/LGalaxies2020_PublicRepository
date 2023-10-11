/*
 * model_dustyields.c
 *
 *  Created on: Oct2016
 *  Last modified: Nov 2017
 *      Author: scottclay
 * 
 *  Adds a model of dust production (via AGB stars, SNe remnants and grain growth
 *  in molecular clouds) and dust destruction (via SNe shock waves).
 *
 *  NOTE: This is a simplified version of model_dust_yields.c, where it is
 *  implicitly assumed that H2_AND_RINGS is on and MAINELEMENTS is off. (07-03-22)
 *
 *  Edited from: Nov2021
 *  Author: robyates
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "allvars.h"
#include "proto.h"

/*

Notes:
* Dust elements are a sub-component of diffuse/cloud elements, which are a sub-component of ColdGas elements. (ROB: Correct?)
* ColdGasDiff + ColdGasClouds = ColdGas
* fI = fraction of metals locked up in dust in the diffuse medium
     = DustColdGasDiff/MetalsColdGasDiff
* fC = fraction of metals locked up in dust in cold clouds
     = DustColdGasClouds/MetalsColdGasClouds
* Any dust production mechanism except for growth in molecular clouds enriches the diffuse medium. 
* When SNe destroys dust, a fraction of dust in both the diffuse and the molecular medium is destroyed.
* There is also a slow sputtering of dust (by cosmic rays and hot gas) that occurs in the diffuse medium.
* Also we do not include Nitrogen, Neon and Sulphur in our dust growth model.
*/

void update_fractions(float dt_bin, float t_acc, float t_exch, float mu, int p, int n, int jj) {
    float fc, acc_fac, exch_fac, fmean, fdiff;

    acc_fac=1.-exp(-dt_bin/t_acc);
    exch_fac=exp(-dt_bin/((1.-mu)*t_exch));
    fc = Gal[p].f_c[jj][n] + (Gal[p].f_cmax[jj][n] - Gal[p].f_c[jj][n])*acc_fac;
    Gal[p].f_c[jj][n] = fc;
    fmean = mu*Gal[p].f_c[jj][n]+(1. - mu)*Gal[p].f_i[jj][n];
    fdiff = (Gal[p].f_c[jj][n] - Gal[p].f_i[jj][n])*exch_fac;
    Gal[p].f_c[jj][n] = fmean + (1. - mu)*fdiff;
    Gal[p].f_i[jj][n] = fmean - mu*fdiff;
    return;
}

void drop_fnan(int jj, int n, int p) {
	if (isnan(Gal[p].f_i[jj][n])) {Gal[p].f_i[jj][n] = 0.;}
    if (isnan(Gal[p].f_c[jj][n])) {Gal[p].f_c[jj][n] = 0.;}
    return;
}

void update_dust_mass(int p, int centralgal, double dt, int nstep, int halonr) {
	int ee;
	int i,k;
	int j;
	double Clouds_tot[RNUM], DustClouds_init[RNUM], DustDiff_init[RNUM], DustClouds_tot[RNUM], DustDiff_new[RNUM], DustClouds_new[RNUM];
	for(j=0;j<RNUM;j++) {
		Clouds_tot[j] = 0.0;
		DustClouds_init[j] = 0.0;
		DustDiff_init[j] = 0.0;
		DustClouds_tot[j] = 0.0;
		DustDiff_new[j] = 0.0;
		DustClouds_new[j] = 0.0;
	} //H2_AND_RINGS

    double DustDiff_Growth, DustClouds_Growth;
#ifdef DUST_DESTRUCTION
	double total_dust_before=0., total_dust_after=0.;
#endif
	double timestep_width; //Width of current timestep in CODE UNITS
	int TimeBin; //Bin in Yield arrays corresponding to current timestep
	double NormAGBDustYieldRate_actual[AGB_DUST_TYPE_NUM];
	double DiskSFR, DiskSFR_physical_units, step_width_times_DiskSFR_physical_units;
	double Dust_Forsterite, Dust_Fayalite, Dust_Enstatite, Dust_Ferrosilite, Dust_Quartz, Dust_Iron, Dust_SiC, Dust_Carbon;
	double H2frac;
	double ColdGasClouds_avail, ColdGasDiff_avail, AGBAllElementsDiff_avail, AGBAllElementsClouds_avail, SNIIAllElementsDiff_avail, SNIIAllElementsClouds_avail;
	double New_dust_diff, New_dust_clouds;
	double num_CO, O_clouds; //Cb_clouds
	double tacc;
	double previoustime, newtime, deltaT;
	previoustime = NumToTime(Gal[p].SnapNum);
	newtime = NumToTime(Gal[p].SnapNum+1);
	deltaT = previoustime - newtime; //ROB (26-11-21): Is this not just the same as dt*STEPS?

	mass_checks(p,"Start of model_dust_yields.c",__LINE__);
	
	timestep_width = dt; //Width of current timestep in CODE UNITS (units cancel out when dividing by SFH bin width, sfh_dt) (12-04-12)
	TimeBin = (STEPS*(Halo[Gal[p].HaloNr].SnapNum-1.0))+nstep; //TimeBin = (STEPS*Gal[p].SnapNum)+nstep; //Bin in Yield tables corresponding to current timestep //TEST!: BRUNO: Snapnum would be +1 too low for a 'jumping' galaxy (14-11-13)


//**************************************************
//Dust grain destruction from supernova shock waves
//**************************************************

#ifdef DUST_DESTRUCTION
	for (j=0;j<RNUM;j++) {
		if ((Gal[p].MetalsColdGasRings[j][0]+Gal[p].MetalsColdGasRings[j][1]+Gal[p].MetalsColdGasRings[j][2]) > 0.0) {
			//For dust destruction we follow the prescription of McKee1989.
			double tdes, survive_frac, R_SN=0.0;

			//Calculate destruction timescale and destruction fraction: ROB: N.B. dust in ColdGas assumed to only be destroyed by SNe in disc (not bulge or ICL).
			/* ROB: (27-11-21): DiskSNIIRate_current_ts and DiskSNIaRate_current_ts only give the SNII rate in the DISC in the LAST TIMESTEP of that snapshot
			 * Gal[p].DiskSNIIRate, on the other hand, outputs the average SNII rate from all stellar components combined (disc+bulge+ICL) from the whole snapshot (similarly to SFR).
			*/
			R_SN = (DiskSNIIRate_current_ts[j] + DiskSNIaRate_current_ts[j]) / (UnitTime_in_s / SEC_PER_YEAR); //To convert from 1/code_time_units to 1/yr. NOTE: I have put the denominator there in brackets, so the correct conversion is UnitTime_in_s / SEC_PER_YEAR, not UnitTime_in_s * SEC_PER_YEAR. (09-02-22)
			tdes = (Gal[p].ColdGasRings[j]*(1.0e10/Hubble_h))/(M_CLEARED * F_SN * R_SN);
			survive_frac = exp(-dt*UnitTime_in_years/tdes);

			//We assume that the SNR will destroy equal amounts of dust in cold clouds and
			//the diffuse medium, but all those will end up as diffuse gas, I guess.
			//Then some will be reaccreted onto cold clouds.
			//Simplest approximation is just to destroy the same fraction in each.
#ifdef FULL_DUST_RATES
			for (ee=0; ee<NUM_ELEMENTS; ee++) {
				total_dust_before += Gal[p].DustColdGasDiffRings_elements[j][ee] + Gal[p].DustColdGasCloudsRings_elements[j][ee];
			}
#endif //FULL_DUST_RATES

			//Update dust elements due to SN dust destruction:
			for (ee=0; ee<NUM_ELEMENTS; ee++) {
				Gal[p].DustColdGasDiff_elements[ee] += Gal[p].DustColdGasDiffRings_elements[j][ee]*(survive_frac-1.); //This will subtract the original DustColdGasDiffRings_elements[j][ee] and add the new "survival" fraction. (08-02-22)
				Gal[p].DustColdGasClouds_elements[ee] += Gal[p].DustColdGasCloudsRings_elements[j][ee]*(survive_frac-1.);
				Gal[p].DustColdGasDiffRings_elements[j][ee] = Gal[p].DustColdGasDiffRings_elements[j][ee]*survive_frac;
				Gal[p].DustColdGasCloudsRings_elements[j][ee] = Gal[p].DustColdGasCloudsRings_elements[j][ee]*survive_frac;
			}

#ifdef FULL_DUST_RATES
			for (ee=0; ee<NUM_ELEMENTS; ee++) {
				total_dust_after += Gal[p].DustColdGasDiffRings_elements[j][ee] + Gal[p].DustColdGasCloudsRings_elements[j][ee];
			}
			mass_checks(p,"Dust from destruction: model_dust_yields.c",__LINE__);
			Gal[p].DustColdGasRates[4] += (total_dust_before - total_dust_after)/(deltaT * UnitTime_in_years);
#endif
		}
	} //for (j=0;j<RNUM;j++)
#endif //DUST_DESTRUCTION

  //***************************************
  partition_gas_and_dust_elements(p);
  mass_checks(p,"After partition_gas_and_dust_elements: model_dust_yields.c",__LINE__);

  for (i=0;i<=Gal[p].sfh_ibin;i++) { //LOOP OVER SFH BINS
	  for(j=0;j<RNUM;j++) { //LOOP OVER RINGS

//*****************************************
//DUST ENRICHMENT FROM AGB DISK STARS INTO COLD GAS:
//*****************************************
#ifdef DUST_AGB
	  if ((Gal[p].sfh_DiskMassRings[j][i] > 0.0) && (Gal[p].MetalsColdGasRings[j][2] > 0.0)) {
		  DiskSFR = Gal[p].sfh_DiskMassRings[j][i]/Gal[p].sfh_dt[i];

		  //pre-calculations to speed up the code
		  DiskSFR_physical_units = DiskSFR * (1.0e10/Hubble_h); //Note: This is NOT in physical units (i.e. NOT in Msun/yr, but in Msun/[code_time_units]). But this is ok, as code_time_units cancel out when multiplying by timestep_width to get 'step_width_times_DiskSFR_physical_units' on the line below ('DiskSFR_physical_units' is never used itself).
		  step_width_times_DiskSFR_physical_units = timestep_width * DiskSFR_physical_units; //ROB: This is the same as DiskSFRxStep_Phys in moedl_yields.c. (17-01-22)
    	
    	//interpolates yields from lookup tables we produced in dust_yield_integrals.c
	    for (k=0;k<AGB_DUST_TYPE_NUM;k++) {
	    	NormAGBDustYieldRate_actual[k] = NormAGBDustYieldRate[TimeBin][i][Zi_saved[j][i]][k]
										   + ((NormAGBDustYieldRate[TimeBin][i][Zi_saved[j][i]+1][k] - NormAGBDustYieldRate[TimeBin][i][Zi_saved[j][i]][k])*Zi_disp_saved[j][i]);
#ifdef FULL_DUST_RATES
	    	Gal[p].DustColdGasRates[0] +=  NormAGBDustYieldRate_actual[k] * DiskSFR_physical_units*dt / (deltaT*UnitTime_in_years);
#endif
	    }

		//Calculate the amount of dust CREATED ----------------------------------------------------------------------
		//These are calculated based on pre-code calculations in dustyield_integrals.c and then multiplied
		//by the SFR here to get the amount of dust created for each specific type(quartz, iron, carbon etc.)
		//and for 3 types of star (M,C,S). 
		//N.B. There is no "unprocessed" component to dust yields, so just the newly-synthesised masses are calculated here:
		Dust_Forsterite = (1.0-fwind_AGB) * step_width_times_DiskSFR_physical_units * NormAGBDustYieldRate_actual[0]; //M_forsterite
		Dust_Fayalite = (1.0-fwind_AGB) * step_width_times_DiskSFR_physical_units * NormAGBDustYieldRate_actual[1]; //M_fayalite
		Dust_Enstatite = (1.0-fwind_AGB) * step_width_times_DiskSFR_physical_units * NormAGBDustYieldRate_actual[2]; //M_enstatite
		Dust_Ferrosilite = (1.0-fwind_AGB) * step_width_times_DiskSFR_physical_units * NormAGBDustYieldRate_actual[3]; //M_ferrosilite
		Dust_Quartz = (1.0-fwind_AGB) * step_width_times_DiskSFR_physical_units * NormAGBDustYieldRate_actual[4]; //M_quartz
		Dust_Quartz += (1.0-fwind_AGB) * step_width_times_DiskSFR_physical_units * NormAGBDustYieldRate_actual[6]; //S_quartz
		Dust_Iron = (1.0-fwind_AGB) * step_width_times_DiskSFR_physical_units * NormAGBDustYieldRate_actual[5]; //M_iron
		Dust_Iron += (1.0-fwind_AGB) * step_width_times_DiskSFR_physical_units * NormAGBDustYieldRate_actual[7]; //S_iron
		Dust_Iron += (1.0-fwind_AGB) * step_width_times_DiskSFR_physical_units * NormAGBDustYieldRate_actual[9]; //C_iron
		Dust_SiC = (1.0-fwind_AGB) * step_width_times_DiskSFR_physical_units * NormAGBDustYieldRate_actual[8]; //C_SiC
		Dust_Carbon = (1.0-fwind_AGB) * step_width_times_DiskSFR_physical_units * NormAGBDustYieldRate_actual[10]; //C_carbon

		for (k=0;k<AGB_DUST_TYPE_NUM;k++)
			if ((1.0-fwind_AGB) * step_width_times_DiskSFR_physical_units * NormAGBDustYieldRate_actual[k] < 0.0)
				printf("***** WARNING: Negative AGB dust yield: SFRxStep_Phys * NormAGBDustYieldRate_actual[%i] = %f\n", k, step_width_times_DiskSFR_physical_units * NormAGBDustYieldRate_actual[k]);

		//Element Conversion -----------------------------------------------------------------------------------
		//Conversion of dust species (i.e. Ferrosilite) into Actual elements to store
		//in correct arrays (i.e. Forsterite -> Mg/Si/O)
		//All the following conversions are done by mass fraction
		//Corrections added at places to avoid dust masses going beyond the mass in gas-phase elements left to form dust from
		//(given that the newly-formed elements have already been added to the gas in model_yields.c).
		for (ee=0; ee<NUM_ELEMENTS; ee++) {
			H2frac = Gal[p].H2fractionRings[j];

			if (ee == Cb_NUM) { //Cb
				New_dust_diff = ((Dust_SiC * SILICONCARBIDE_Cb_FRAC) + (Dust_Carbon * 1.0)) * (1.0 - H2frac);
				New_dust_clouds = ((Dust_SiC * SILICONCARBIDE_Cb_FRAC) + (Dust_Carbon * 1.0)) * H2frac;
			}
			else if (ee == Si_NUM) { //Si
				New_dust_diff = ((Dust_Forsterite * FORSTERITE_Si_FRAC) + (Dust_Fayalite * FAYALITE_Si_FRAC) + (Dust_Enstatite * ENSTATITE_Si_FRAC)
							  + (Dust_Ferrosilite * FERROSILITE_Si_FRAC) + (Dust_Quartz * QUARTZ_Si_FRAC) + (Dust_SiC * SILICONCARBIDE_Si_FRAC)) * (1.0 - H2frac);
				New_dust_clouds = ((Dust_Forsterite * FORSTERITE_Si_FRAC) + (Dust_Fayalite * FAYALITE_Si_FRAC) + (Dust_Enstatite * ENSTATITE_Si_FRAC)
								+ (Dust_Ferrosilite * FERROSILITE_Si_FRAC) + (Dust_Quartz * QUARTZ_Si_FRAC) + (Dust_SiC * SILICONCARBIDE_Si_FRAC)) * H2frac;
			}
			else if (ee == O_NUM) { //O
				New_dust_diff = ((Dust_Forsterite * FORSTERITE_O_FRAC) + (Dust_Fayalite * FAYALITE_O_FRAC) + (Dust_Enstatite * ENSTATITE_O_FRAC)
							  + (Dust_Ferrosilite * FERROSILITE_O_FRAC) + (Dust_Quartz * QUARTZ_O_FRAC)) * (1.0 - H2frac);
				New_dust_clouds = ((Dust_Forsterite * FORSTERITE_O_FRAC) + (Dust_Fayalite * FAYALITE_O_FRAC) + (Dust_Enstatite * ENSTATITE_O_FRAC)
								+ (Dust_Ferrosilite * FERROSILITE_O_FRAC) + (Dust_Quartz * QUARTZ_O_FRAC)) * H2frac;
			}
			else if (ee == Mg_NUM) { //Mg
				New_dust_diff = ((Dust_Forsterite * FORSTERITE_Mg_FRAC) + (Dust_Enstatite * ENSTATITE_Mg_FRAC)) * (1.0 - H2frac);
				New_dust_clouds = ((Dust_Forsterite * FORSTERITE_Mg_FRAC) + (Dust_Enstatite * ENSTATITE_Mg_FRAC)) * H2frac;
			}
			else if (ee == Fe_NUM) { //Fe
				New_dust_diff = ((Dust_Fayalite * FAYALITE_Fe_FRAC) + (Dust_Ferrosilite * FERROSILITE_Fe_FRAC) + (Dust_Iron * 1.0)) * (1.0 - H2frac);
				New_dust_clouds = ((Dust_Fayalite * FAYALITE_Fe_FRAC) + (Dust_Ferrosilite * FERROSILITE_Fe_FRAC) + (Dust_Iron * 1.0)) * H2frac;
			}
			else {
				New_dust_diff = 0.0;
				New_dust_clouds = 0.0;
			}

			//Check how much gas there is actually available to form dust:
			ColdGasDiff_avail = Gal[p].ColdGasDiffRings_elements[j][ee] - Gal[p].DustColdGasDiffRings_elements[j][ee]; //Total diffuse gas available to form dust
			ColdGasClouds_avail = Gal[p].ColdGasCloudsRings_elements[j][ee] - Gal[p].DustColdGasCloudsRings_elements[j][ee]; //Total cloud gas available to form dust
			AGBAllElementsDiff_avail = (1.0-fwind_AGB) * DiskAGBAllElements_ts[j][ee] * (1.0-Gal[p].H2fractionRings[j]); //Newly-ejected element mass into diffuse gas available to form dust
			AGBAllElementsClouds_avail = (1.0-fwind_AGB) * DiskAGBAllElements_ts[j][ee] * Gal[p].H2fractionRings[j]; //Newly-ejected element mass into clouds available to form dust
			//Add newly-formed dust to ColdGas:
			Gal[p].DustColdGasDiffRings_elements[j][ee]  += min(New_dust_diff, min(ColdGasDiff_avail, AGBAllElementsDiff_avail));
			Gal[p].DustColdGasCloudsRings_elements[j][ee]  += min(New_dust_clouds, min(ColdGasClouds_avail, AGBAllElementsClouds_avail));
			//Adding the same small increments to the global quantities within this loop over rings should work here: (02-02-22):
			Gal[p].DustColdGasDiff_elements[ee]  += min(New_dust_diff, min(ColdGasDiff_avail, AGBAllElementsDiff_avail));
			Gal[p].DustColdGasClouds_elements[ee]  += min(New_dust_clouds, min(ColdGasClouds_avail, AGBAllElementsClouds_avail));
		}
	} //if ( (Gal[p].sfh_DiskMass[i] > 0.0) && (Gal[p].MetalsColdGas[2] > 0.0) )
    
    mass_checks(p,"Dust from AGBs: model_dust_yields.c",__LINE__);
#endif //DUST_AGB


//*****************************************
//DUST ENRICHMENT FROM SNII FROM DISK STARS INTO COLD GAS:
//*****************************************

#ifdef DUST_SNII
    if ((Gal[p].sfh_DiskMassRings[j][i] > 0.0) && (Gal[p].MetalsColdGasRings[j][0] > 0.0)) {
#ifdef FULL_DUST_RATES
    	//This is estimating the dust in various compounds (say silicates) using the amount of the particular element (say silicon) produced in a process
    	Gal[p].DustColdGasRates[1] += (SNII_prevstep_Cold_Si[j][i] * eta_SNII_Sil * A_Sil_dust/A_Si)/(deltaT * UnitTime_in_years);
		Gal[p].DustColdGasRates[1] += (SNII_prevstep_Cold_Si[j][i] * eta_SNII_SiC * A_SiC_dust/A_Si)/(deltaT * UnitTime_in_years);
		Gal[p].DustColdGasRates[1] += (SNII_prevstep_Cold_Cb[j][i] * eta_SNII_Cb  * A_Cb_dust/A_Cb) /(deltaT * UnitTime_in_years);
		Gal[p].DustColdGasRates[1] += (SNII_prevstep_Cold_Fe[j][i] * eta_SNII_Fe  * A_Fe_dust/A_Fe )/(deltaT * UnitTime_in_years);
#endif //FULL_DUST_RATES

		//Create dust (based on the prescription of Zhukovska2008)---------------------------
		//SNII_prevstep_x is calculated in model_yields.c
		//It is the amount of a specific metal (i.e. Si) returned to gas phase from SFH_BIN i. //produced in the last timestep
		double Dust_Silicates = SNII_prevstep_Cold_Si[j][i] * eta_SNII_Sil * A_Sil_dust/A_Si;
		double Dust_Iron      = SNII_prevstep_Cold_Fe[j][i] * eta_SNII_Fe  * A_Fe_dust/A_Fe;
		double Dust_SiC	      = SNII_prevstep_Cold_Si[j][i] * eta_SNII_SiC * A_SiC_dust/A_Si;
		double Dust_Carbon    = SNII_prevstep_Cold_Cb[j][i] * eta_SNII_Cb  * A_Cb_dust/A_Cb;
            
		//Element conversion -----------------------------------------------------------------
		//Conversion of dust species (i.e. Silicates) into Actual elements to store
		//in correct arrays (i.e. Silicates -> Mg/Si/Fe/O)
		//All the following conversions are done by mass fraction
		//SNII Silicates -------------------
		//ROB: (04-01-22) The mass fractions used here don't make sense to me.
		//Assuming that the silicate mass is just olivines + pyroxenes (and f_ol = 0.32), I get A_Sil_dust = 121.62, close to the quoted value used here (in h_params.h).
		//But, e.g. A_Si / A_Sil_dust = 0.2310, whereas 0.210432 is used here...
		for (ee=0; ee<NUM_ELEMENTS; ee++) {
			H2frac = Gal[p].H2fractionRings[j];
				if (ee == Cb_NUM) { //Cb
					New_dust_diff = ((Dust_SiC * SILICONCARBIDE_Cb_FRAC) + (Dust_Carbon * 1.0)) * (1.0 - H2frac);
					New_dust_clouds = ((Dust_SiC * SILICONCARBIDE_Cb_FRAC) + (Dust_Carbon * 1.0)) * H2frac;
				}
				else if (ee == Si_NUM) { //Si
					New_dust_diff = ((Dust_Silicates * SILICATES_Si_FRAC) + (Dust_SiC * SILICONCARBIDE_Si_FRAC)) * (1.0 - H2frac);
					New_dust_clouds = ((Dust_Silicates * SILICATES_Si_FRAC) + (Dust_SiC * SILICONCARBIDE_Si_FRAC)) * H2frac;
				}
				else if (ee == O_NUM) { //O
					New_dust_diff = (Dust_Silicates * SILICATES_O_FRAC) * (1.0 - H2frac);
					New_dust_clouds = (Dust_Silicates * SILICATES_O_FRAC) * H2frac;
				}
				else if (ee == Mg_NUM) { //Mg
					New_dust_diff = (Dust_Silicates * SILICATES_Mg_FRAC) * (1.0 - H2frac);
					New_dust_clouds = (Dust_Silicates * SILICATES_Mg_FRAC) * H2frac;
				}
				else if (ee == Fe_NUM) { //Fe
					New_dust_diff = ((Dust_Silicates * SILICATES_Fe_FRAC) + (Dust_Iron * 1.0)) * (1.0 - H2frac);
					New_dust_clouds = ((Dust_Silicates * SILICATES_Fe_FRAC) + (Dust_Iron * 1.0)) * H2frac;
				}
				else {
					New_dust_diff = 0.0;
					New_dust_clouds = 0.0;
				}

				//Check how much gas there is actually available to form dust:
				ColdGasDiff_avail = Gal[p].ColdGasDiffRings_elements[j][ee] - Gal[p].DustColdGasDiffRings_elements[j][ee]; //Total diffuse gas available to form dust
				ColdGasClouds_avail = Gal[p].ColdGasCloudsRings_elements[j][ee] - Gal[p].DustColdGasCloudsRings_elements[j][ee]; //Total cloud gas available to form dust
				SNIIAllElementsDiff_avail = (1.0-fwind_SNII) * DiskSNIIAllElements_ts[j][ee] * (1.0-Gal[p].H2fractionRings[j]); //Newly-ejected element mass into diffuse gas available to form dust
				SNIIAllElementsClouds_avail = (1.0-fwind_SNII) * DiskSNIIAllElements_ts[j][ee] * Gal[p].H2fractionRings[j]; //Newly-ejected element mass into clouds available to form dust
				//Add newly-formed dust to ColdGas:
				Gal[p].DustColdGasDiffRings_elements[j][ee]  += min(New_dust_diff, min(ColdGasDiff_avail, SNIIAllElementsDiff_avail));
				Gal[p].DustColdGasCloudsRings_elements[j][ee]  += min(New_dust_clouds, min(ColdGasClouds_avail, SNIIAllElementsClouds_avail));
				//Adding the same small increments to the global quantities within this loop over rings should work here: (02-02-22):
				Gal[p].DustColdGasDiff_elements[ee]  += min(New_dust_diff, min(ColdGasDiff_avail, SNIIAllElementsDiff_avail));
				Gal[p].DustColdGasClouds_elements[ee]  += min(New_dust_clouds, min(ColdGasClouds_avail, SNIIAllElementsClouds_avail));
			} //for (ee=0; ee<NUM_ELEMENTS; ee++)

			mass_checks(p,"Dust from SNe-II: model_dust_yields.c",__LINE__);
	}//if ((Gal[p].sfh_DiskMass[i] > 0.0) && (Gal[p].MetalsColdGas[0] > 0.0))
#endif //DUST_SNII
	
	
//*****************************************
//DUST ENRICHMENT FROM SNIA FROM DISK STARS INTO COLD GAS:
//*****************************************
#ifdef DUST_SNIA
    if ((Gal[p].sfh_DiskMassRings[j][i] > 0.0) && (Gal[p].MetalsColdGasRings[j][1] > 0.0)) {
		double Dust_Iron = SNIa_prevstep_Cold_Fe[j][i] * eta_SNIa_Fe  * A_Fe_dust/A_Fe;
        #ifdef FULL_DUST_RATES
		    Gal[p].DustColdGasRates[2]  += (SNIa_prevstep_Cold_Fe[j][i] * eta_SNIa_Fe  * A_Fe_dust/A_Fe)/(deltaT * UnitTime_in_years);
        #endif
		
		//ROB: I have added an extra condition here - that the newly-added dust doesn't exceed the newly-added metals from this enrichment channel: (28-01-22):
		Gal[p].DustColdGasDiffRings_elements[j][Fe_NUM] += min(Dust_Iron * 1.0 * (1.0-Gal[p].H2fractionRings[j]),
												   min(Gal[p].ColdGasDiffRings_elements[j][Fe_NUM] - Gal[p].DustColdGasDiffRings_elements[j][ee],
													   (1.0-fwind_SNIa) * DiskSNIaAllElements_ts[j][Fe_NUM] * (1.0-Gal[p].H2fractionRings[j])));
		Gal[p].DustColdGasCloudsRings_elements[j][Fe_NUM] += min(Dust_Iron * 1.0 * Gal[p].H2fractionRings[j],
													 min(Gal[p].ColdGasCloudsRings_elements[j][Fe_NUM] - Gal[p].DustColdGasCloudsRings_elements[j][Fe_NUM],
														 (1.0-fwind_SNIa) * DiskSNIaAllElements_ts[j][Fe_NUM] * Gal[p].H2fractionRings[j]));
		Gal[p].DustColdGasDiff_elements[Fe_NUM] += min(Dust_Iron * 1.0 * (1.0-Gal[p].H2fractionRings[j]),
												   min(Gal[p].ColdGasDiffRings_elements[j][Fe_NUM] - Gal[p].DustColdGasDiffRings_elements[j][Fe_NUM],
													   (1.0-fwind_SNIa) * DiskSNIaAllElements_ts[j][Fe_NUM] * (1.0-Gal[p].H2fractionRings[j])));
		Gal[p].DustColdGasClouds_elements[Fe_NUM] += min(Dust_Iron * 1.0 * Gal[p].H2fractionRings[j],
													 min(Gal[p].ColdGasCloudsRings_elements[j][Fe_NUM] - Gal[p].DustColdGasCloudsRings_elements[j][Fe_NUM],
														 (1.0-fwind_SNIa) * DiskSNIaAllElements_ts[j][Fe_NUM] * Gal[p].H2fractionRings[j]));

		mass_checks(p,"Dust from SNe-Ia: model_dust_yields.c",__LINE__);
		
	}
#endif //DUST_SNIA
  } //for(j=0;j<RNUM;j++)
} //loop over SFH bins


//*****************************************
//Dust grain growth inside molecular clouds 
//*****************************************
#ifdef DUST_GROWTH
for(j=0;j<RNUM;j++) {
   	if (((Gal[p].MetalsColdGasRings[j][0]+Gal[p].MetalsColdGasRings[j][1]+Gal[p].MetalsColdGasRings[j][2])>0.0) && (Gal[p].ColdGasRings[j] > 0.0)) {

        mass_checks(p,"Dust from growth: model_dust_yields.c",__LINE__);
        
        //The number of CO molecules that can be produced from available Carbon and Oxygen
        //in the clouds. By default, assuming only 30% C is locked up as CO in clouds.
        //ROB (10-01-22): This sets the number of CO molecules in the gas phase of the ISM to the minimum of either:
        //		(a) number of carbon atoms in gas phase
        //		(b) Cmax_CO (= 0.3): Fraction of total carbon atoms in gas and dust (this is the max amount of carbon allowed in CO, set in input file)
        //		(c) number of oxygen atoms in gas phase
        num_CO = min((Gal[p].ColdGasCloudsRings_elements[j][Cb_NUM]-Gal[p].DustColdGasCloudsRings_elements[j][Cb_NUM])/A_Cb,
                     min(Gal[p].ColdGasCloudsRings_elements[j][Cb_NUM]/A_Cb*Cmax_CO, (Gal[p].ColdGasCloudsRings_elements[j][O_NUM]-Gal[p].DustColdGasCloudsRings_elements[j][O_NUM])/A_O));

        //Mass of Carbon and Oxygen available in the gas+dust of clouds for grain growth:
        //double Cb_clouds = Gal[p].ColdGasClouds_elements[Cb_NUM] - num_CO*A_Cb; //ROB: Note, Cb_clouds isn't used anywhere... (17-01-22)
        O_clouds = Gal[p].ColdGasCloudsRings_elements[j][O_NUM] - num_CO*A_O;
        
        //*********************************************
        //Silicates:
        //*********************************************
        /*
        We assume that oxygen is volatile such that it can exist in the dust phase only in the 
        form of compounds. We adopt from Zhukovska et al.(2008) that the dominant silicate species
        are olivine and pyroxene in the ratio 32:68. 
        Olivine: [Mg_x Fe_{1-x}]_2 Si O_4
        Pyroxene: Mg_x Fe_{1-x} Si O_3 
        x is assumed to be 0.8 in their paper, this value doesn't seem to affect the dust masses 
        that much. It is because this doesn't significantly affect the amount of oxygen left in 
        the cold gas.
        */
        /* ROB (10-01-22): Olivine ([Mg,Fe]2SiO4) describes the range of magnesium iron silicates between Forsterite (Mg2SiO4) and Fayalite (Fe2SiO4).
         * Forsterite contains no Fe, and Fayalite contains no Mg. All Olivines in between contain some mixture of Mg and Fe.
         * Here, we follow the simplification of Zhukovska+08 by only considering a typical Olivine mix of Mg:Fe = 0.8:0.2 = number ratio of Mg to Fe.
         * The same typical mix is assumed for Pyroxene ([Mg,Fe]SiO3).
         * Also, it is assumed that the number ratio of olivine to pyroxene is 32:68.
         * This gives rise to the parameters x =  0.8 and f_ol = 0.32, which I have now made parameters in h_params.h.
         * I have also re-written these "number of atoms per typical olivine/pyroxene silicate molecule" for O, Mg, Si, and Fe calculations accordingly.
         */
        float N_O = OLIVINE_NUMFRAC*4. + (1.-OLIVINE_NUMFRAC)*3.;
		float N_Mg = OLIVINE_NUMFRAC*2.*OLIVINE_Mg_NUMFRAC + (1.-OLIVINE_NUMFRAC)*PYROXENE_Mg_NUMFRAC;
		float N_Si = OLIVINE_NUMFRAC*1. + (1.-OLIVINE_NUMFRAC)*1.;
		float N_Fe = OLIVINE_NUMFRAC*2.*(1-OLIVINE_Mg_NUMFRAC) + (1.-OLIVINE_NUMFRAC)*(1.-PYROXENE_Mg_NUMFRAC);

		//ROB: This is a nested min() function to find number of olivine/pyroxene silicate molecules formed, given the mass available of the least-abundant constituent element in the gas phase of the molecular clouds:
		float num_Silicates = max(0., min((O_clouds-Gal[p].DustColdGasCloudsRings_elements[j][O_NUM])/(N_O*A_O),
        								 min((Gal[p].ColdGasCloudsRings_elements[j][Mg_NUM]-Gal[p].DustColdGasCloudsRings_elements[j][Mg_NUM])/(N_Mg*A_Mg),
        									min((Gal[p].ColdGasCloudsRings_elements[j][Si_NUM]-Gal[p].DustColdGasCloudsRings_elements[j][Si_NUM])/(N_Si*A_Si),
        									   (Gal[p].ColdGasCloudsRings_elements[j][Fe_NUM]-Gal[p].DustColdGasCloudsRings_elements[j][Fe_NUM])/(N_Fe*A_Fe)))));

        //*********************************************
        //Iron Oxide:
        //*********************************************
        float N_O_IronOxide = (HEMATITE_NUMFRAC*4. + (1.-HEMATITE_NUMFRAC)*3.);
        float N_Fe_IronOxide = (HEMATITE_NUMFRAC*3. + (1.-HEMATITE_NUMFRAC)*2.);

        //*********************************************
        //Calculate the f_max values for C and O (i.e. the maximum condensation fractions):
        float num_iron_oxide = max(0., min((O_clouds-num_Silicates*N_O*A_O-Gal[p].DustColdGasCloudsRings_elements[j][O_NUM])/(N_O_IronOxide*A_O),
										  (Gal[p].ColdGasCloudsRings_elements[j][Fe_NUM]-num_Silicates*N_Fe*A_Fe-Gal[p].DustColdGasCloudsRings_elements[j][Fe_NUM])/(N_Fe_IronOxide*A_Fe)));
        float f_O_max =  (num_Silicates*N_O*A_O + num_iron_oxide*N_O_IronOxide*A_O + Gal[p].DustColdGasCloudsRings_elements[j][O_NUM])/Gal[p].ColdGasCloudsRings_elements[j][O_NUM]; //Mass fraction of all oxygen in clouds that is in dust
        float f_Cb_max =  1.0 - num_CO*A_Cb/Gal[p].ColdGasCloudsRings_elements[j][Cb_NUM]; //Mass fraction of all carbon in clouds not in CO

        //*********************************************
        //Calculate total mass in clouds and dust mass in clouds for updating Gal[p].t_acc, and dust mass in diffuse gas for total gas growth rates:
		for (ee=0; ee<NUM_ELEMENTS; ee++) {
			Clouds_tot[j] += Gal[p].ColdGasCloudsRings_elements[j][ee];
			DustClouds_init[j] += Gal[p].DustColdGasCloudsRings_elements[j][ee];
			DustDiff_init[j] += Gal[p].DustColdGasDiffRings_elements[j][ee];
		}

		if (DustClouds_init[j] > 0.0) {
			tacc = Dust_tAcc0*(Clouds_tot[j]/DustClouds_init[j]);
		}
		else tacc = 1e15;

		//*********************************************
		//Calculate the new masses of each element that are in dust in the diffuse gas and clouds:

        //Here we use an approximation for calculating the fraction of dust produced in
        //clouds and the diffused ISM which is valid for a constant accretion timescale
        //across the timestep. Since that is not the case as it depends on the amount of 
        //dust in the clouds, we divide the calculation in constant bins of the time step
        //in between snaps. The dust mass in clouds is calculated again to update t_acc.
        
        //Note: Dust molecules injected into the ISM has elements bonded with itself and 
        //only certain other elements. 
        // * Cb: Cb, Si
        // * O: Mg, Si, Fe
        // * Mg: Si, O
        // * Si: with all
        // * Fe: Si, Fe 
        // * Maybe bring this into play like implemented by Zhukovska
        
        //*****
		//Check that f_Cb_max and f_O_max are between 0. and 1., and assign their values to Gal[p].f_cmax (for use in update_fractions()):
        Gal[p].f_cmax[j][Cb_NUM] = max(0., min(1., f_Cb_max));
		if (isnan(Gal[p].f_cmax[j][Cb_NUM])) {Gal[p].f_cmax[j][Cb_NUM] = 1.;}
		Gal[p].f_cmax[j][O_NUM] = max(0., min(1., f_O_max));
		if (isnan(Gal[p].f_cmax[j][O_NUM])) {Gal[p].f_cmax[j][O_NUM] = 1.;}
        
		//*****
		//Update dust element fractions:
		float dt_bins = (dt*UnitTime_in_years); //Width of one timestep in years
        for(ee=0;ee<NUM_ELEMENTS;ee++) {
			Gal[p].f_i[j][ee] = min(1., Gal[p].DustColdGasDiffRings_elements[j][ee]/Gal[p].ColdGasDiffRings_elements[j][ee]); //initial fraction of element locked-up in dust in the diffuse medium
			Gal[p].f_c[j][ee] = min(1., Gal[p].DustColdGasCloudsRings_elements[j][ee]/Gal[p].ColdGasCloudsRings_elements[j][ee]);
			drop_fnan(j, ee, p);
			update_fractions(dt_bins, tacc, Dust_tExch, Gal[p].H2fractionRings[j], p, ee, j); //updates the fraction of element locked-up in clouds and diffuse gas, using t_acc, t_exch, f_cmax, initial f_c and f_i, and the H2 fraction.
			drop_fnan(j, ee, p);
			//Updating the amount of dust after this mini-step:
			Gal[p].DustColdGasCloudsRings_elements[j][ee] = Gal[p].f_c[j][ee]*Gal[p].ColdGasCloudsRings_elements[j][ee];
			//Gal[p].DustColdGasClouds_elements[ee] += Gal[p].f_c[j][ee]*Gal[p].ColdGasCloudsRings_elements[j][ee];
		}
            
        //*****
		//Recalculate properties again, after updating dust fractions:
        num_CO = min((Gal[p].ColdGasCloudsRings_elements[j][Cb_NUM]-Gal[p].DustColdGasCloudsRings_elements[j][Cb_NUM])/A_Cb,
                     min(Gal[p].ColdGasCloudsRings_elements[j][Cb_NUM]/A_Cb*Cmax_CO, (Gal[p].ColdGasCloudsRings_elements[j][O_NUM]-Gal[p].DustColdGasCloudsRings_elements[j][O_NUM])/A_O));
        O_clouds = Gal[p].ColdGasCloudsRings_elements[j][O_NUM] - num_CO*A_O;
        num_Silicates = max(0., min((O_clouds-Gal[p].DustColdGasCloudsRings_elements[j][O_NUM])/(N_O*A_O),
									 min((Gal[p].ColdGasCloudsRings_elements[j][Mg_NUM]-Gal[p].DustColdGasCloudsRings_elements[j][Mg_NUM])/(N_Mg*A_Mg),
										min((Gal[p].ColdGasCloudsRings_elements[j][Si_NUM]-Gal[p].DustColdGasCloudsRings_elements[j][Si_NUM])/(N_Si*A_Si),
										   (Gal[p].ColdGasCloudsRings_elements[j][Fe_NUM]-Gal[p].DustColdGasCloudsRings_elements[j][Fe_NUM])/(N_Fe*A_Fe)))));

        num_iron_oxide = max(0., min((O_clouds-num_Silicates*N_O*A_O-Gal[p].DustColdGasCloudsRings_elements[j][O_NUM])/(N_O_IronOxide*A_O),
        							 (Gal[p].ColdGasCloudsRings_elements[j][Fe_NUM]-num_Silicates*N_Fe*A_Fe-Gal[p].DustColdGasCloudsRings_elements[j][Fe_NUM])/(N_Fe_IronOxide*A_Fe)));
        f_O_max = (num_Silicates*N_O*A_O + num_iron_oxide*N_O_IronOxide*A_O + Gal[p].DustColdGasCloudsRings_elements[j][O_NUM])/Gal[p].ColdGasCloudsRings_elements[j][O_NUM];

        f_Cb_max =  1.0 - num_CO*A_Cb/Gal[p].ColdGasCloudsRings_elements[j][Cb_NUM];
        Gal[p].f_cmax[j][Cb_NUM] = max(0., min(1., f_Cb_max));
        if (isnan(Gal[p].f_cmax[j][Cb_NUM])) {Gal[p].f_cmax[j][Cb_NUM] = 1.;}

		Gal[p].f_cmax[j][O_NUM] = max(0., min(1., f_O_max));
		if (isnan(Gal[p].f_cmax[j][O_NUM])) {Gal[p].f_cmax[j][O_NUM] = 1.;}

		for (ee=0; ee<NUM_ELEMENTS; ee++)
			DustClouds_tot[j] += Gal[p].DustColdGasCloudsRings_elements[j][ee];
		tacc = Dust_tAcc0*(Clouds_tot[j]/DustClouds_tot[j]);

		for(ee=0;ee<NUM_ELEMENTS;ee++) {
			Gal[p].DustColdGasDiffRings_elements[j][ee] = Gal[p].f_i[j][ee] * Gal[p].ColdGasDiffRings_elements[j][ee];
			Gal[p].DustColdGasCloudsRings_elements[j][ee] = Gal[p].f_c[j][ee] * Gal[p].ColdGasCloudsRings_elements[j][ee];
			//Gal[p].DustColdGasDiff_elements[ee] += Gal[p].f_i[j][ee] * Gal[p].ColdGasDiffRings_elements[j][ee];
			//Gal[p].DustColdGasClouds_elements[ee] += Gal[p].f_c[j][ee] * Gal[p].ColdGasCloudsRings_elements[j][ee];
			DustDiff_new[j] += Gal[p].ColdGasCloudsRings_elements[j][ee];
			DustClouds_new[j] += Gal[p].DustColdGasCloudsRings_elements[j][ee];
		}
		DustDiff_Growth = DustDiff_new[j] - DustDiff_init[j];
		DustClouds_Growth = DustClouds_new[j] - DustClouds_init[j];

#ifdef FULL_DUST_RATES
		Gal[p].DustColdGasRates[3] += (DustDiff_Growth + DustClouds_Growth)/(deltaT * UnitTime_in_years);
#endif
        
		mass_checks(p,"End of dust from growth: model_dust_yields.c",__LINE__);
        
  } //if (((Gal[p].MetalsColdGas[0]+Gal[p].MetalsColdGas[1]+Gal[p].MetalsColdGas[2])>0.0) && (Gal[p].ColdGas > 0.0))
} //for(j=0;j<RNUM;j++)

//Calculate final global dust elements masses from all rings:
for(ee=0;ee<NUM_ELEMENTS;ee++) {
   	double DustColdGasDiff_elements_allRings=0., DustColdGasClouds_elements_allRings=0.;
	for(j=0;j<RNUM;j++) {
		DustColdGasDiff_elements_allRings += Gal[p].DustColdGasDiffRings_elements[j][ee];
		DustColdGasClouds_elements_allRings += Gal[p].DustColdGasCloudsRings_elements[j][ee];
	}
	Gal[p].DustColdGasDiff_elements[ee] = DustColdGasDiff_elements_allRings;
	Gal[p].DustColdGasClouds_elements[ee] = DustColdGasClouds_elements_allRings;
}

//Calculate final global t_des from all rings:
double tot_SNRate=0., tot_ColdGasRings=0.;
for (int jjj=0;jjj<RNUM;jjj++) {
	tot_SNRate += DiskSNIIRate_current_ts[jjj] + DiskSNIaRate_current_ts[jjj];
	tot_ColdGasRings += Gal[p].ColdGasRings[jjj];
}
Gal[p].t_des = (tot_ColdGasRings*(1.0e10/Hubble_h))/(M_CLEARED * F_SN * (tot_SNRate/(UnitTime_in_s / SEC_PER_YEAR))); //Note, this will be outputted as the (global) destruction timescale of the last timestep of each snapshot, rather than the average across the whole snapshot. (09-02-22)

//Calculate final global t_acc from all rings:
double Clouds_tot_allRings=0., DustClouds_tot_allRings=0.;
for(j=0;j<RNUM;j++) {
	Clouds_tot_allRings += Clouds_tot[j];
	DustClouds_tot_allRings += DustClouds_tot[j];
}
Gal[p].t_acc = Dust_tAcc0*(Clouds_tot_allRings/DustClouds_tot_allRings);
#endif //DUST_GROWTH
    
}//update_dust_mass()

