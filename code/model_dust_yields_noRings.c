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

void update_fractions(float dt_bin, float t_acc, float t_exch, float mu, int p, int n) {
    float fc, acc_fac, exch_fac, fmean, fdiff;
    
    acc_fac=1.-exp(-dt_bin/t_acc);
    exch_fac=exp(-dt_bin/((1.-mu)*t_exch));
    fc = Gal[p].f_c[n] + (Gal[p].f_cmax[n] - Gal[p].f_c[n])*acc_fac;
    Gal[p].f_c[n] = fc;
    fmean = mu*Gal[p].f_c[n]+(1. - mu)*Gal[p].f_i[n];
    fdiff = (Gal[p].f_c[n] - Gal[p].f_i[n])*exch_fac;
    Gal[p].f_c[n] = fmean + (1. - mu)*fdiff;
    Gal[p].f_i[n] = fmean - mu*fdiff;
    
    return;
}

void drop_fnan(int n, int p) {
	if (isnan(Gal[p].f_i[n])) {Gal[p].f_i[n] = 0.;}
    if (isnan(Gal[p].f_c[n])) {Gal[p].f_c[n] = 0.;}

    return;
}

void update_dust_mass(int p, int centralgal, double dt, int nstep, int halonr) {
	int ee;
	double timestep_width; //Width of current timestep in CODE UNITS
	int TimeBin; //Bin in Yield arrays corresponding to current timestep
	double NormAGBDustYieldRate_actual[AGB_DUST_TYPE_NUM];
	double DiskSFR, DiskSFR_physical_units, step_width_times_DiskSFR_physical_units;
	double Dust_Forsterite, Dust_Fayalite, Dust_Enstatite, Dust_Ferrosilite, Dust_Quartz, Dust_Iron, Dust_SiC, Dust_Carbon;
	double ColdGasClouds_avail, ColdGasDiff_avail, AGBAllElementsDiff_avail, AGBAllElementsClouds_avail, SNIIAllElementsDiff_avail, SNIIAllElementsClouds_avail;
	double New_dust_diff, New_dust_clouds;
	double previoustime, newtime, deltaT;
	previoustime = NumToTime(Gal[p].SnapNum);
	newtime = NumToTime(Gal[p].SnapNum+1);
	deltaT = previoustime - newtime; //ROB (26-11-21): Is this not just the same as dt*STEPS?

	mass_checks(p,"Start of model_dust_yields.c",__LINE__);
	
	timestep_width = dt; //Width of current timestep in CODE UNITS (units cancel out when dividing by SFH bin width, sfh_dt) (12-04-12)
	TimeBin = (STEPS*(Halo[Gal[p].HaloNr].SnapNum-1.0))+nstep; //TimeBin = (STEPS*Gal[p].SnapNum)+nstep; //Bin in Yield tables corresponding to current timestep //TEST!: BRUNO: Snapnum would be +1 too low for a 'jumping' galaxy (14-11-13)

//**************************************************
//Dust grain destruction from supernovae shock waves
//**************************************************

#ifdef DUST_DESTRUCTION
	if ((Gal[p].MetalsColdGas[0]+Gal[p].MetalsColdGas[1]+Gal[p].MetalsColdGas[2])>0.0) {
		//For dust destruction we follow the prescription of McKee1989.
		float tdes, survive_frac, R_SN=0.0; //, des_frac
        
        //Calculate destruction timescale and destruction fraction: ROB: N.B. dust in ColdGas assumed to only be destroyed by SNe in disc (not bulge or ICL).
		//ROB: (27-11-21): DiskSNIIRate_current_ts and DiskSNIaRate_current_ts only give the SNII rate in the DISC in the LAST TIMESTEP of that snapshot
		//Gal[p].DiskSNIIRate, on the other hand, outputs the average SNII rate from all stellar components combined (disc+bulge+ICL) from the whole snapshot (similarly to SFR).
		R_SN = (DiskSNIIRate_current_ts + DiskSNIaRate_current_ts) / (UnitTime_in_s / SEC_PER_YEAR); //To convert from 1/code_time_units to 1/yr
		tdes = (Gal[p].ColdGas*(1.0e10/Hubble_h))/(M_CLEARED * F_SN * R_SN);
		Gal[p].t_des = tdes;
		survive_frac = exp(-dt*UnitTime_in_years/tdes);
        
		//We assume that the SNR will destroy equal amounts of dust in cold clouds and 
		//the diffuse medium, but all those will end up as diffuse gas, I guess.  
		//Then some will be reaccreted onto cold clouds.
	    //Simplest approximation is just to destroy the same fraction in each.
#ifdef FULL_DUST_RATES
		double total_dust_before=0.;
		for (ee=0; ee<NUM_ELEMENTS; ee++)
			total_dust_before += Gal[p].DustColdGasDiff_elements[ee] + Gal[p].DustColdGasClouds_elements[ee];
#endif
		
		//Update dust elements due to SN dust destruction:
		for (ee=0; ee<NUM_ELEMENTS; ee++) {
	    	Gal[p].DustColdGasDiff_elements[ee] = Gal[p].DustColdGasDiff_elements[ee]*survive_frac;
	    	Gal[p].DustColdGasClouds_elements[ee] = Gal[p].DustColdGasClouds_elements[ee]*survive_frac;
	    }

		printf("Snapnum %i | tdes = %e | survive_frac = %f\n", Gal[p].SnapNum, tdes, survive_frac);

	    /*
		#ifdef DCR_Dest
		    //****************************
		    //Destruction of dust from miscellaneous processes, only done in the diffuse medium
		    //****************************
		    //Timescale of destruction chosen to be 1Gyr (arbitrary)
		    float CR_timescale = 1e9;
	        survive_frac = exp(-(dt*UnitTime_in_years/CR_timescale));
		    
		    Gal[p].DustColdGasDiff_elements=elements_add(elements_init(),Gal[p].DustColdGasDiff_elements,survive_frac);

	    #endif    
		
		#ifdef DDestHIIregion
		    //****************************
		    //Dust destruction in the diffused medium by evaporation in HII regions
		    //****************************
		    //This is a rough implementation based on the mass of HII regions created by OB stars
		    //as presented in Tielens' book on ISM. Need to get better values of the photon flux for
		    //the star of average mass produced in the particular time step. Current implementation
		    //has only negligible effect on the dust mass.
		    float MetClouds_tot = elements_total(Gal[p].ColdGasClouds_elements);
		    
	        float MetDiff_tot = elements_total(Gal[p].ColdGasDiff_elements);
		    
		    double sfr = Gal[p].Sfr * UNITMASS_IN_G / UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS;
		    double gas_cleared = (0.0496173*sfr*dt*UnitTime_in_years/13.7054)*2.24;
		    double frac_diff = gas_cleared/MetDiff_tot;
		    double frac_clouds = 0.;
		    
		    if (frac_diff > 1.) {
		    
		        frac_diff = 1.;
		        frac_clouds = (gas_cleared-MetDiff_tot)/MetClouds_tot;
	        }
	        
            Gal[p].DustColdGasDiff_elements=elements_add(Gal[p].DustColdGasDiff_elements,Gal[p].DustColdGasDiff_elements,-frac_diff);
		    Gal[p].DustColdGasClouds_elements=elements_add(Gal[p].DustColdGasDiff_elements,Gal[p].DustColdGasDiff_elements,-frac_clouds);
		#endif    
		*/
		
#ifdef FULL_DUST_RATES
	    double total_dust_after=0.;
	    for (ee=0; ee<NUM_ELEMENTS; ee++)
	    	total_dust_after += Gal[p].DustColdGasDiff_elements[ee] + Gal[p].DustColdGasClouds_elements[ee];
	    Gal[p].DustColdGasRates[4] += (total_dust_before - total_dust_after)/(deltaT * UnitTime_in_years);
#endif
	    mass_checks(p,"Dust from destruction: model_dust_yields.c",__LINE__);
		}	
#endif //DUST_DESTRUCTION

  //(28-01-22): partition_gas_and_dust_elements() was moved here from just above grain growth. Didn't make any difference to the outputs, but is more consistent when checking the diff/cloud mass fractions for AGB/SN dust production.
  //(02-02-22): Note that partition_gas_and_dust_elements() now has it's own internal loop over rings, so make sure it is only called outside of a loop over rings in this file (i.e. here at the top).
  partition_gas_and_dust_elements(p);

  int i,k; //j
  for (i=0;i<=Gal[p].sfh_ibin;i++) { //LOOP OVER SFH BINS

	  /*ROB: (02-02-22): When incorporating rings, should probably add a big loop over rings starting here. Something like this (from model_yields.c):
	   * #ifdef H2_AND_RINGS
	   * 	for(jj=0;jj<RNUM;jj++)
	   * 		if (Gal[p].DiskMassRings[jj] > 0.0) //ROB: Discs can be destroyed (i.e. converted in to bulges). So only calculate enrichment from stars born in the disc if there is still a disc
	   * 			if (Gal[p].sfh_DiskMassRings[jj][ii] > 0.0) {
	   * #else
	   * 	if (Gal[p].DiskMass > 0.0) //ROB: Discs can be destroyed (i.e. converted in to bulges). So only calculate enrichment from stars born in the disc if there is still a disc at the current timestep. This if statement has scope until the end of the following if statement.
	   * 		if (Gal[p].sfh_DiskMass[ii] > 0.0) {
	   * #endif
	   */

//*****************************************
//DUST ENRICHMENT FROM AGB DISK STARS INTO COLD GAS:
//*****************************************

#ifdef DUST_AGB		
    if ( (Gal[p].sfh_DiskMass[i] > 0.0) && (Gal[p].MetalsColdGas[2] > 0.0) ) {
     	//pre-calculations to speed up the code
    	DiskSFR = Gal[p].sfh_DiskMass[i]/Gal[p].sfh_dt[i];
    	DiskSFR_physical_units = DiskSFR * (1.0e10/Hubble_h); //Note: This is NOT in physical units (i.e. NOT in Msun/yr, but in Msun/[code_time_units]). But this is ok, as code_time_units cancel out when multiplying by timestep_width to get 'step_width_times_DiskSFR_physical_units' on the line below ('DiskSFR_physical_units' is never used itself).
    	step_width_times_DiskSFR_physical_units = timestep_width * DiskSFR_physical_units; //ROB: This is the same as DiskSFRxStep_Phys in moedl_yields.c. (17-01-22)

    	/* ROB: The Zi and Zi_disp variables are not needed in model_dust_yields.c, as they calculated in model_yields.c and stored in the global variables Zi_saved and Zi_disp_saved instead. (04-01-22)
    	 * 	// Note: This approach is ok, but only if update_dust_mass() is called after update_yields_and_return_mass() in main.c, and both are looking at the same galaxy.
    	 * 	// Alternative approach would be to run find_initial_metallicity_dust() here, to re-calculate Zi and Zi_disp again, independently of model_yields.c.
    	 */
    	
    	//interpolates yields from lookup tables we produced in dust_yield_integrals.c
	    for (k=0;k<AGB_DUST_TYPE_NUM;k++) {
	    	NormAGBDustYieldRate_actual[k] = NormAGBDustYieldRate[TimeBin][i][Zi_saved[i]][k] + ((NormAGBDustYieldRate[TimeBin][i][Zi_saved[i]+1][k] - NormAGBDustYieldRate[TimeBin][i][Zi_saved[i]][k])*Zi_disp_saved[i]);
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
#ifndef MAINELEMENTS
			if (ee == 2) { //Cb
				New_dust_diff = ((Dust_SiC * SILICONCARBIDE_Cb_FRAC) + (Dust_Carbon * 1.0)) * (1.0 - Gal[p].H2fraction);
				New_dust_clouds = ((Dust_SiC * SILICONCARBIDE_Cb_FRAC) + (Dust_Carbon * 1.0)) * Gal[p].H2fraction;
			}
			else if (ee == 4) { //O
				New_dust_diff = ((Dust_Forsterite * FORSTERITE_O_FRAC) + (Dust_Fayalite * FAYALITE_O_FRAC) + (Dust_Enstatite * ENSTATITE_O_FRAC)
							  + (Dust_Ferrosilite * FERROSILITE_O_FRAC) + (Dust_Quartz * QUARTZ_O_FRAC)) * (1.0 - Gal[p].H2fraction);
				New_dust_clouds = ((Dust_Forsterite * FORSTERITE_O_FRAC) + (Dust_Fayalite * FAYALITE_O_FRAC) + (Dust_Enstatite * ENSTATITE_O_FRAC)
						  	    + (Dust_Ferrosilite * FERROSILITE_O_FRAC) + (Dust_Quartz * QUARTZ_O_FRAC)) * Gal[p].H2fraction;
			}
			else if (ee == 6) { //Mg
				New_dust_diff = ((Dust_Forsterite * FORSTERITE_Mg_FRAC) + (Dust_Enstatite * ENSTATITE_Mg_FRAC)) * (1.0 - Gal[p].H2fraction);
				New_dust_clouds = ((Dust_Forsterite * FORSTERITE_Mg_FRAC) + (Dust_Enstatite * ENSTATITE_Mg_FRAC)) * Gal[p].H2fraction;
			}
			else if (ee == 7) { //Si
				New_dust_diff = ((Dust_Forsterite * FORSTERITE_Si_FRAC) + (Dust_Fayalite * FAYALITE_Si_FRAC) + (Dust_Enstatite * ENSTATITE_Si_FRAC)
							  + (Dust_Ferrosilite * FERROSILITE_Si_FRAC) + (Dust_Quartz * QUARTZ_Si_FRAC) + (Dust_SiC * SILICONCARBIDE_Si_FRAC)) * (1.0 - Gal[p].H2fraction);
				New_dust_clouds = ((Dust_Forsterite * FORSTERITE_Si_FRAC) + (Dust_Fayalite * FAYALITE_Si_FRAC) + (Dust_Enstatite * ENSTATITE_Si_FRAC)
								+ (Dust_Ferrosilite * FERROSILITE_Si_FRAC) + (Dust_Quartz * QUARTZ_Si_FRAC) + (Dust_SiC * SILICONCARBIDE_Si_FRAC)) * Gal[p].H2fraction;
			}
			else if (ee == 10) { //Fe
				New_dust_diff = ((Dust_Fayalite * FAYALITE_Fe_FRAC) + (Dust_Ferrosilite * FERROSILITE_Fe_FRAC) + (Dust_Iron * 1.0)) * (1.0 - Gal[p].H2fraction);
				New_dust_clouds = ((Dust_Fayalite * FAYALITE_Fe_FRAC) + (Dust_Ferrosilite * FERROSILITE_Fe_FRAC) + (Dust_Iron * 1.0)) * Gal[p].H2fraction;
			}
#else
			if (ee == 2) { //O
				New_dust_diff = ((Dust_Forsterite * FORSTERITE_O_FRAC) + (Dust_Fayalite * FAYALITE_O_FRAC) + (Dust_Enstatite * ENSTATITE_O_FRAC)
							  + (Dust_Ferrosilite * FERROSILITE_O_FRAC) + (Dust_Quartz * QUARTZ_O_FRAC)) * (1.0 - Gal[p].H2fraction);
				New_dust_clouds = ((Dust_Forsterite * FORSTERITE_O_FRAC) + (Dust_Fayalite * FAYALITE_O_FRAC) + (Dust_Enstatite * ENSTATITE_O_FRAC)
								+ (Dust_Ferrosilite * FERROSILITE_O_FRAC) + (Dust_Quartz * QUARTZ_O_FRAC)) * Gal[p].H2fraction;
			}
			else if (ee == 3) { //Mg
				New_dust_diff = ((Dust_Forsterite * FORSTERITE_Mg_FRAC) + (Dust_Enstatite * ENSTATITE_Mg_FRAC)) * (1.0 - Gal[p].H2fraction);
				New_dust_clouds = ((Dust_Forsterite * FORSTERITE_Mg_FRAC) + (Dust_Enstatite * ENSTATITE_Mg_FRAC)) * Gal[p].H2fraction;
			}
			else if (ee == 4) { //Fe
				New_dust_diff = ((Dust_Fayalite * FAYALITE_Fe_FRAC) + (Dust_Ferrosilite * FERROSILITE_Fe_FRAC) + (Dust_Iron * 1.0)) * (1.0 - Gal[p].H2fraction);
				New_dust_clouds = ((Dust_Fayalite * FAYALITE_Fe_FRAC) + (Dust_Ferrosilite * FERROSILITE_Fe_FRAC) + (Dust_Iron * 1.0)) * Gal[p].H2fraction;
			}
#endif
			else {
				New_dust_diff = 0.0;
				New_dust_clouds = 0.0;
			}
			//Check how much gas there is actually available to form dust:
			ColdGasDiff_avail = Gal[p].ColdGasDiff_elements[ee] - Gal[p].DustColdGasDiff_elements[ee]; //Total diffuse gas available to form dust
			ColdGasClouds_avail = Gal[p].ColdGasClouds_elements[ee] - Gal[p].DustColdGasClouds_elements[ee]; //Total cloud gas available to form dust
			AGBAllElementsDiff_avail = (1.0-fwind_AGB) * DiskAGBAllElements_ts[ee] * (1.0-Gal[p].H2fraction); //Newly-ejected element mass into diffuse gas available to form dust
			AGBAllElementsClouds_avail = (1.0-fwind_AGB) * DiskAGBAllElements_ts[ee] * Gal[p].H2fraction; //Newly-ejected element mass into clouds available to form dust
			//Add newly-formed dust to ColdGas:
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
	if ((Gal[p].sfh_DiskMass[i] > 0.0) && (Gal[p].MetalsColdGas[0] > 0.0)) {
	
	    #ifdef FULL_DUST_RATES
	        //This is estimating the dust in various compounds (say silicates) using the amount of the particular element (say silicon) produced in a process
			Gal[p].DustColdGasRates[1] += (SNII_prevstep_Cold_Si[i] * eta_SNII_Sil * A_Sil_dust/A_Si)/(deltaT * UnitTime_in_years);
			Gal[p].DustColdGasRates[1] += (SNII_prevstep_Cold_Fe[i] * eta_SNII_Fe  * A_Fe_dust/A_Fe )/(deltaT * UnitTime_in_years);
			Gal[p].DustColdGasRates[1] += (SNII_prevstep_Cold_Si[i] * eta_SNII_SiC * A_SiC_dust/A_Si)/(deltaT * UnitTime_in_years);
			Gal[p].DustColdGasRates[1] += (SNII_prevstep_Cold_Cb[i] * eta_SNII_Cb  * A_Cb_dust/A_Cb) /(deltaT * UnitTime_in_years);
	    #endif

			//Create dust (based on the prescription of Zhukovska2008)---------------------------
			//SNII_prevstep_x is calculated in model_yields.c 
			//It is the amount of a specific metal (i.e. Si) returned to gas phase from SFH_BIN i. //produced in the last timestep

			//ROB: Here, the mass of newly-formed dust is determined by the mass of the "key element" that has just been returned by SNe.
			//So, e.g. if the key element makes up 25% of the total dust mass, then the mass of dust formed is 1/0.25 = 4 times greater than the mass of the key element formed * the conversion efficiency (eta).
			//This is why the A_dust/A_key_element terms are included here.
			double Dust_Silicates = SNII_prevstep_Cold_Si[i] * eta_SNII_Sil * A_Sil_dust/A_Si;
			double Dust_Iron      = SNII_prevstep_Cold_Fe[i] * eta_SNII_Fe  * A_Fe_dust/A_Fe;
			double Dust_SiC	      = SNII_prevstep_Cold_Si[i] * eta_SNII_SiC * A_SiC_dust/A_Si;
			double Dust_Carbon    = SNII_prevstep_Cold_Cb[i] * eta_SNII_Cb  * A_Cb_dust/A_Cb;	
            
			//Element conversion -----------------------------------------------------------------
			//Conversion of dust species (i.e. Silicates) into Actual elements to store
			//in correct arrays (i.e. Silicates -> Mg/Si/Fe/O)
			//All the following conversions are done by mass fraction
			//SNII Silicates -------------------
			//ROB: (04-01-22) The mass fractions used here don't make sense to me.
			//Assuming that the silicate mass is just olivines + pyroxenes (and f_ol = 0.32), I get A_Sil_dust = 121.62, close to the quoted value used here (in h_params.h).
			//But, e.g. A_Si / A_Sil_dust = 0.2310, whereas 0.210432 is used here...
			for (ee=0; ee<NUM_ELEMENTS; ee++) {
#ifndef MAINELEMENTS
				if (ee == 2) { //Cb
					New_dust_diff = ((Dust_SiC * SILICONCARBIDE_Cb_FRAC) + (Dust_Carbon * 1.0)) * (1.0 - Gal[p].H2fraction);
					New_dust_clouds = ((Dust_SiC * SILICONCARBIDE_Cb_FRAC) + (Dust_Carbon * 1.0)) * Gal[p].H2fraction;
				}
				else if (ee == 4) { //O
					New_dust_diff = (Dust_Silicates * SILICATES_O_FRAC) * (1.0 - Gal[p].H2fraction);
					New_dust_clouds = (Dust_Silicates * SILICATES_O_FRAC) * Gal[p].H2fraction;
				}
				else if (ee == 6) { //Mg
					New_dust_diff = (Dust_Silicates * SILICATES_Mg_FRAC) * (1.0 - Gal[p].H2fraction);
					New_dust_clouds = (Dust_Silicates * SILICATES_Mg_FRAC) * Gal[p].H2fraction;
				}
				else if (ee == 7) { //Si
					New_dust_diff = ((Dust_Silicates * SILICATES_Si_FRAC) + (Dust_SiC * SILICONCARBIDE_Si_FRAC)) * (1.0 - Gal[p].H2fraction);
					New_dust_clouds = ((Dust_Silicates * SILICATES_Si_FRAC) + (Dust_SiC * SILICONCARBIDE_Si_FRAC)) * Gal[p].H2fraction;
				}
				else if (ee == 10) { //Fe
					New_dust_diff = ((Dust_Silicates * SILICATES_Fe_FRAC) + (Dust_Iron * 1.0)) * (1.0 - Gal[p].H2fraction);
					New_dust_clouds = ((Dust_Silicates * SILICATES_Fe_FRAC) + (Dust_Iron * 1.0)) * Gal[p].H2fraction;
				}
#else
				if (ee == 2) { //O
					New_dust_diff = (Dust_Silicates * SILICATES_O_FRAC) * (1.0 - Gal[p].H2fraction);
					New_dust_clouds = (Dust_Silicates * SILICATES_O_FRAC) * Gal[p].H2fraction;
				}
				else if (ee == 3) { //Mg
					New_dust_diff = (Dust_Silicates * SILICATES_Mg_FRAC) * (1.0 - Gal[p].H2fraction);
					New_dust_clouds = (Dust_Silicates * SILICATES_Mg_FRAC) * Gal[p].H2fraction;
				}
				else if (ee == 4) { //Fe
					New_dust_diff = ((Dust_Silicates * SILICATES_Fe_FRAC) + (Dust_Iron * 1.0)) * (1.0 - Gal[p].H2fraction);
					New_dust_clouds = ((Dust_Silicates * SILICATES_Fe_FRAC) + (Dust_Iron * 1.0)) * Gal[p].H2fraction;
				}
#endif
				else {
					New_dust_diff = 0.0;
					New_dust_clouds = 0.0;
				}
				//Check how much gas there is actually available to form dust:
				ColdGasDiff_avail = Gal[p].ColdGasDiff_elements[ee] - Gal[p].DustColdGasDiff_elements[ee]; //Total diffuse gas available to form dust
				ColdGasClouds_avail = Gal[p].ColdGasClouds_elements[ee] - Gal[p].DustColdGasClouds_elements[ee]; //Total cloud gas available to form dust
				SNIIAllElementsDiff_avail = (1.0-fwind_SNII) * DiskSNIIAllElements_ts[ee] * (1.0-Gal[p].H2fraction); //Newly-ejected element mass into diffuse gas available to form dust
				SNIIAllElementsClouds_avail = (1.0-fwind_SNII) * DiskSNIIAllElements_ts[ee] * Gal[p].H2fraction; //Newly-ejected element mass into clouds available to form dust
				//Add newly-formed dust to ColdGas:
				Gal[p].DustColdGasDiff_elements[ee]  += min(New_dust_diff, min(ColdGasDiff_avail, SNIIAllElementsDiff_avail));
				Gal[p].DustColdGasClouds_elements[ee]  += min(New_dust_clouds, min(ColdGasClouds_avail, SNIIAllElementsClouds_avail));
			}

			mass_checks(p,"Dust from SNe-II: model_dust_yields.c",__LINE__);
	}//if ((Gal[p].sfh_DiskMass[i] > 0.0) && (Gal[p].MetalsColdGas[0] > 0.0))
#endif //DUST_SNII
	
	
//*****************************************
//DUST ENRICHMENT FROM SNIA FROM DISK STARS INTO COLD GAS:
//*****************************************

#ifdef DUST_SNIA		
	if ((Gal[p].sfh_DiskMass[i] > 0.0) && (Gal[p].MetalsColdGas[1] >0.0)) {
		double Dust_Iron = SNIa_prevstep_Cold_Fe[i] * eta_SNIa_Fe  * A_Fe_dust/A_Fe;
	
        #ifdef FULL_DUST_RATES		
		    Gal[p].DustColdGasRates[2]  += (SNIa_prevstep_Cold_Fe[i] * eta_SNIa_Fe  * A_Fe_dust/A_Fe)/(deltaT * UnitTime_in_years);
        #endif	
		
		//ROB: I have added an extra condition here - that the newly-added dust doesn't exceed the newly-added metals from this enrichment channel: (28-01-22):
#ifndef MAINELEMENTS
		Gal[p].DustColdGasDiff_elements[10] += min(Dust_Iron * 1.0 * (1.0-Gal[p].H2fraction),
												   min(Gal[p].ColdGasDiff_elements[10] - Gal[p].DustColdGasDiff_elements[10],
													   (1.0-fwind_SNIa) * DiskSNIaAllElements_ts[10] * (1.0-Gal[p].H2fraction)));
		Gal[p].DustColdGasClouds_elements[10] += min(Dust_Iron * 1.0 * Gal[p].H2fraction,
												 	 min(Gal[p].ColdGasClouds_elements[10] - Gal[p].DustColdGasClouds_elements[10],
												 		 (1.0-fwind_SNIa) * DiskSNIaAllElements_ts[10] * Gal[p].H2fraction));
#else
		Gal[p].DustColdGasDiff_elements[4] += min(Dust_Iron * 1.0 * (1.0-Gal[p].H2fraction),
												   min(Gal[p].ColdGasDiff_elements[4] - Gal[p].DustColdGasDiff_elements[4],
													   (1.0-fwind_SNIa) * DiskSNIaAllElements_ts[4] * (1.0-Gal[p].H2fraction)));
		Gal[p].DustColdGasClouds_elements[4] += min(Dust_Iron * 1.0 * Gal[p].H2fraction,
													 min(Gal[p].ColdGasClouds_elements[4] - Gal[p].DustColdGasClouds_elements[4],
														 (1.0-fwind_SNIa) * DiskSNIaAllElements_ts[4] * Gal[p].H2fraction));
#endif

		mass_checks(p,"Dust from SNe-Ia: model_dust_yields.c",__LINE__);
		
	}
#endif //DUST_SNIA

} //loop over SFH bins


//*****************************************
//Dust grain growth inside molecular clouds 
//*****************************************

#ifdef DUST_GROWTH
    if (((Gal[p].MetalsColdGas[0]+Gal[p].MetalsColdGas[1]+Gal[p].MetalsColdGas[2])>0.0) && (Gal[p].ColdGas > 0.0)) {
			
        //Growth implementation requires the fraction of elements in both the media
        //to be distributed according to the method of molecular hydrogen computation 
        //used, since H2_fraction (formerly mu_gas) is passed into the dust growth equations.

        //partition_gas_and_dust_elements(p); //ROB: Now done at the top of this file. (28-01-22)
        //ROB: Now that Gal[p].H2fraction is used, is this function redundant? Is the partition of gas/dust into clouds and diff gas already correct at this point? (19-01-22)
        //(25-01-22): No. When partition_gas_and_dust_elements() is commented-out, Cloud and Diff gas masses don't match ColdGas*H2fraction any more.
        mass_checks(p,"Dust from growth: model_dust_yields.c",__LINE__);
        
        //The number of CO molecules that can be produced from available Carbon and Oxygen
        //in the clouds. xxAssuming only 30% C is locked up as CO in cloudsxx
        //ROB (10-01-22): This sets the number of CO molecules in the gas phase of the ISM to the minimum of either:
        //		(a) number of carbon atoms in gas phase
        //		(b) Cmax_CO (= 0.3): Fraction of total carbon atoms in gas and dust (this is the max amount of carbon allowed in CO, set in input file)
        //		(c) number of oxygen atoms in gas phase
        float num_CO = min((Gal[p].ColdGasClouds_elements[2]-Gal[p].DustColdGasClouds_elements[2])/A_Cb,
                     min(Gal[p].ColdGasClouds_elements[2]/A_Cb*Cmax_CO, (Gal[p].ColdGasClouds_elements[4]-Gal[p].DustColdGasClouds_elements[4])/A_O));

        //Mass of Carbon and Oxygen available in the gas+dust of clouds for grain growth:
        double Cb_clouds = Gal[p].ColdGasClouds_elements[2] - num_CO*A_Cb; //ROB: Note, Cb_clouds isn't used anywhere... (17-01-22)
        double O_clouds = Gal[p].ColdGasClouds_elements[4] - num_CO*A_O;
        
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
        float num_Silicates = max(0., min((O_clouds-Gal[p].DustColdGasClouds_elements[4])/(N_O*A_O),
        								 min((Gal[p].ColdGasClouds_elements[6]-Gal[p].DustColdGasClouds_elements[6])/(N_Mg*A_Mg),
        									min((Gal[p].ColdGasClouds_elements[7]-Gal[p].DustColdGasClouds_elements[7])/(N_Si*A_Si),
        									   (Gal[p].ColdGasClouds_elements[10]-Gal[p].DustColdGasClouds_elements[10])/(N_Fe*A_Fe)))));
        
        //*********************************************
        //Iron Oxide:
        //*********************************************
        /* ROB (10-01-22): I don't like the way the number of iron oxide molecules is calculated here.
         * There is a hidden normalisation constant of N_Fe_IronOxide floating around, because N_Oxide is actually the number ratio of N_O_IronOxide/N_Fe_IronOxide.
         * Therefore, I have re-written this bit to be exactly the same as for silicates above, with variables like N_O_IronOxide and num_iron_oxide properly representing the unnormalised "number" of molecules.
         */

        float N_O_IronOxide = (HEMATITE_NUMFRAC*4. + (1.-HEMATITE_NUMFRAC)*3.);
        float N_Fe_IronOxide = (HEMATITE_NUMFRAC*3. + (1.-HEMATITE_NUMFRAC)*2.);
		float num_iron_oxide = max(0., min((O_clouds-num_Silicates*N_O*A_O-Gal[p].DustColdGasClouds_elements[4])/(N_O_IronOxide*A_O),
										  (Gal[p].ColdGasClouds_elements[10]-num_Silicates*N_Fe*A_Fe-Gal[p].DustColdGasClouds_elements[10])/(N_Fe_IronOxide*A_Fe)));

        //*********************************************
        //Calculate the f_max values for C and O (i.e. the maximum condensation fractions):

        float f_Cb_max =  1.0 - num_CO*A_Cb/Gal[p].ColdGasClouds_elements[2]; //Mass fraction of all carbon in clouds not in CO
        float f_O_max =  (num_Silicates*N_O*A_O + num_iron_oxide*N_O_IronOxide*A_O + Gal[p].DustColdGasClouds_elements[4])/Gal[p].ColdGasClouds_elements[4]; //Mass fraction of all oxygen in clouds that is in dust

        //*********************************************
        //Calculate total mass in clouds and dust mass in clouds for updating Gal[p].t_acc, and dust mass in diffuse gas for total gas growth rates:
        
        float Clouds_tot=0., DustClouds_init=0., DustDiff_init=0.;
		for (ee=0; ee<NUM_ELEMENTS; ee++) {
			Clouds_tot += Gal[p].ColdGasClouds_elements[ee];
			DustClouds_init += Gal[p].DustColdGasClouds_elements[ee];
			DustDiff_init += Gal[p].DustColdGasDiff_elements[ee];
		}

		//*********************************************
		//Calculate the new masses of each element that are in dust in the diffuse gas and clouds:
        if (DustClouds_init == 0) {
            Gal[p].t_acc = 1e15;
        }
        else {
            Gal[p].t_acc = Dust_tAcc0*(Clouds_tot/DustClouds_init);
        }
               
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
        
        //ROB: Check that f_Cb_max and f_O_max are between 0. and 1., and assign their values to Gal[p].f_cmax (for use in update_fractions()):
#ifndef MAINELEMENTS
        Gal[p].f_cmax[2] = max(0., min(1., f_Cb_max));
		if (isnan(Gal[p].f_cmax[2])) {Gal[p].f_cmax[2] = 1.;}
		Gal[p].f_cmax[4] = max(0., min(1., f_O_max));
		if (isnan(Gal[p].f_cmax[4])) {Gal[p].f_cmax[4] = 1.;}
#else
		Gal[p].f_cmax[2] = max(0., min(1., f_O_max));
		if (isnan(Gal[p].f_cmax[2])) {Gal[p].f_cmax[2] = 1.;}
#endif
        
		//*****
		float DustClouds_tot = 0.0;
		float dt_bins = (dt*UnitTime_in_years); //Width of one timestep in years
        for(ee=0;ee<NUM_ELEMENTS;ee++) {
			Gal[p].f_i[ee] = min(1., Gal[p].DustColdGasDiff_elements[ee]/Gal[p].ColdGasDiff_elements[ee]); //initial fraction of element locked-up in dust in the diffuse medium
			Gal[p].f_c[ee] = min(1., Gal[p].DustColdGasClouds_elements[ee]/Gal[p].ColdGasClouds_elements[ee]);
			drop_fnan(ee, p);
			update_fractions(dt_bins, Gal[p].t_acc, Dust_tExch, Gal[p].H2fraction, p, ee); //updates the fraction of element locked-up in clouds and diffuse gas, using t_acc, t_exch, f_cmax, initial f_c and f_i, and the H2 fraction.
			drop_fnan(ee, p);
			//Updating the amount of dust after this mini-step:
			Gal[p].DustColdGasClouds_elements[ee] = Gal[p].f_c[ee]*Gal[p].ColdGasClouds_elements[ee];
		}
            
		/* ROB: This re-update seems to be redundant. The update is already done inside the loop above: (27-01-22)
		 *
		 * for(ee=0;ee<NUM_ELEMENTS;ee++)
		 * 		Gal[p].DustColdGasClouds_elements[ee] = Gal[p].f_c[ee]*Gal[p].ColdGasClouds_elements[ee]; //ROB: This should work fine, as f_c for H, He, etc should be 0.0. (27-01-22)
		 */
            
		//Recalculate properties again, after updating dust fractions:
		num_CO = min((Gal[p].ColdGasClouds_elements[2]-Gal[p].DustColdGasClouds_elements[2])/A_Cb,
					 min(Gal[p].ColdGasClouds_elements[2]/A_Cb*Cmax_CO, (Gal[p].ColdGasClouds_elements[4]-Gal[p].DustColdGasClouds_elements[4])/A_O));

		O_clouds = Gal[p].ColdGasClouds_elements[4] - num_CO*A_O;

		f_Cb_max =  1.0 - num_CO*A_Cb/Gal[p].ColdGasClouds_elements[2];

		num_Silicates = max(0., min((O_clouds-Gal[p].DustColdGasClouds_elements[4])/(N_O*A_O),
									min((Gal[p].ColdGasClouds_elements[6]-Gal[p].DustColdGasClouds_elements[6])/(N_Mg*A_Mg),
											min((Gal[p].ColdGasClouds_elements[7]-Gal[p].DustColdGasClouds_elements[7])/(N_Si*A_Si),
													(Gal[p].ColdGasClouds_elements[10]-Gal[p].DustColdGasClouds_elements[10])/(N_Fe*A_Fe)))));

		num_iron_oxide = max(0., min((O_clouds-num_Silicates*N_O*A_O-Gal[p].DustColdGasClouds_elements[4])/(N_O_IronOxide*A_O),
									(Gal[p].ColdGasClouds_elements[10]-num_Silicates*N_Fe*A_Fe-Gal[p].DustColdGasClouds_elements[10])/(N_Fe_IronOxide*A_Fe)));

		f_O_max = (num_Silicates*N_O*A_O + num_iron_oxide*N_O_IronOxide*A_O + Gal[p].DustColdGasClouds_elements[4])/Gal[p].ColdGasClouds_elements[4]; //Mass fraction of all oxygen in clouds that is in dust


#ifndef MAINELEMENTS
        Gal[p].f_cmax[2] = max(0., min(1., f_Cb_max));
		if (isnan(Gal[p].f_cmax[2])) {Gal[p].f_cmax[2] = 1.;}
		Gal[p].f_cmax[4] = max(0., min(1., f_O_max));
		if (isnan(Gal[p].f_cmax[4])) {Gal[p].f_cmax[4] = 1.;}
#else
		Gal[p].f_cmax[2] = max(0., min(1., f_O_max));
		if (isnan(Gal[p].f_cmax[2])) {Gal[p].f_cmax[2] = 1.;}
#endif

		for (ee=0; ee<NUM_ELEMENTS; ee++)
			DustClouds_tot += Gal[p].DustColdGasClouds_elements[ee]; //Not sure this is what we want if we change to N > 1 ... (ROB 11-01-22)

		Gal[p].t_acc = Dust_tAcc0*(Clouds_tot/DustClouds_tot);


		//*****
        for(ee=0;ee<NUM_ELEMENTS;ee++) {
			Gal[p].DustColdGasDiff_elements[ee] = Gal[p].f_i[ee] * Gal[p].ColdGasDiff_elements[ee];
			Gal[p].DustColdGasClouds_elements[ee] = Gal[p].f_c[ee] * Gal[p].ColdGasClouds_elements[ee];
		}
		
		float Dust_Diff_new=0., Dust_Clouds_new=0.;
		for (ee=0; ee<NUM_ELEMENTS; ee++) {
			Dust_Diff_new += Gal[p].ColdGasClouds_elements[ee];
			Dust_Clouds_new += Gal[p].DustColdGasClouds_elements[ee];
		}
		float Dust_DiffGrowth = Dust_Diff_new - DustDiff_init;
		float Dust_CloudsGrowth = Dust_Clouds_new - DustClouds_init;
		
        #ifdef FULL_DUST_RATES
		Gal[p].DustColdGasRates[3] += (Dust_DiffGrowth + Dust_CloudsGrowth)/(deltaT * UnitTime_in_years);
        #endif
        
		mass_checks(p,"End of dust from growth: model_dust_yields.c",__LINE__);
        
}
#endif //DUST_GROWTH
    
}//update_dust_mass()

