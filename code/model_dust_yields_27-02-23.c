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

#ifdef H2_AND_RINGS
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
#else
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
#endif //H2_AND_RINGS

#ifdef H2_AND_RINGS
void drop_fnan(int jj, int n, int p) {
	if (isnan(Gal[p].f_i[jj][n])) {Gal[p].f_i[jj][n] = 0.;}
    if (isnan(Gal[p].f_c[jj][n])) {Gal[p].f_c[jj][n] = 0.;}
    return;
}
#else
void drop_fnan(int n, int p) {
	if (isnan(Gal[p].f_i[n])) {Gal[p].f_i[n] = 0.;}
    if (isnan(Gal[p].f_c[n])) {Gal[p].f_c[n] = 0.;}
    return;
}
#endif //H2_AND_RINGS

double calc_cleared_mass(double frac_Cb) {
	double m_cleared; //_clouds, m_cleared_diff;
	/*
	//Option 1) Default way (use fixed value of m_cleared used by Vijayan+19):
	m_cleared = M_CLEARED; //[Msun]
	*/

	/*
	//Option 2) Harder way (requires an estimate of n_H and fitting formulae from Hu+19 (eqn.29)):
	double log_nH, m_cleared_C_th, m_cleared_Si_th, m_cleared_C_nth, m_cleared_Si_nth;
	log_nH = log10(((1.-Gal[p].H2fraction)*Gal[p].ColdGas_elements[0])/ringVol); //[1/cm^3] The volumetric number density of atomic H in the Cold Gas is used here, to mimic the "hydrogen number density of gas" used in Hu+19
	m_cleared_C_th   = pow(10, 2.65 + 0.061*log_nH - 0.100*(log_nH*log_nH));
	m_cleared_C_nth  = pow(10, 2.74 - 0.360*log_nH + 0.023*(log_nH*log_nH) + 0.049*(log_nH*log_nH*log_nH) - 0.028*(log_nH*log_nH*log_nH*log_nH));
	m_cleared_Si_th  = pow(10, 2.70 + 0.018*log_nH - 0.100*(log_nH*log_nH));
	m_cleared_Si_nth = pow(10, 2.93 - 0.370*log_nH + 0.010*(log_nH*log_nH) + 0.040*(log_nH*log_nH*log_nH) - 0.024*(log_nH*log_nH*log_nH*log_nH));
	//...
	*/

	//Option 3) Easier way (assume the "typical SN environment density of n_H = 0.3 cm^-3" quoted by Hu+19, and the C/(C+Si) ratio to weight the overall m_cleared):
	m_cleared = (frac_Cb*M_CLEARED_Cb) + ((1.-frac_Cb)*M_CLEARED_Si);

	return m_cleared;
}

void update_dust_mass(int p, int centralgal, double dt, int nstep, int halonr) {
	int ee;
	int i,k;
#ifdef H2_AND_RINGS
	int j;
	double Clouds_tot[RNUM], DustClouds_init[RNUM], DustDiff_init[RNUM], DustClouds_tot[RNUM], DustDiff_new[RNUM], DustClouds_new[RNUM];
	double DustColdGasCloudsRings_elements_test1, DustColdGasCloudsRings_elements_test2, DustColdGasCloudsRings_elements_test2a, DustColdGasCloudsRings_elements_test2b, DustColdGasCloudsRings_elements_test2c, DustColdGasCloudsRings_elements_test3, DustColdGasCloudsRings_elements_test4;
	for(j=0;j<RNUM;j++) {
		Clouds_tot[j] = 0.0;
		DustClouds_init[j] = 0.0;
		DustDiff_init[j] = 0.0;
		DustClouds_tot[j] = 0.0;
		DustDiff_new[j] = 0.0;
		DustClouds_new[j] = 0.0;
	} //H2_AND_RINGS
#else
    double Clouds_tot=0., DustClouds_init=0., DustDiff_init=0., DustClouds_tot=0., DustDiff_new=0., DustClouds_new=0.;
#endif //H2_AND_RINGS

    double DustDiff_Growth, DustClouds_Growth;
#ifdef DUST_DESTRUCTION
	double total_dust_before=0., total_dust_after=0.;
#ifdef DUST_HOTGAS
	double frac_destroyed_HotGas;
	double Tvir, HotGasDensity, Delta_grainRadius, tau_sput, HotGasDust, Delta_Msput;
#endif //DUST_HOTGAS
#endif //DUST_DESTRUCTION
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
#ifdef DUST_DESTRUCTION
	double frac_Cb, m_cleared, tdes, survive_frac, R_SN=0.0;
	//double frac_Cb_clouds, frac_Cb_diff, m_cleared_clouds, m_cleared_diff, tdes_clouds, tdes_diff, survive_frac_clouds, survive_frac_diff, R_SN=0.0;
#endif //DUST_DESTRUCTION
	double previoustime, newtime, deltaT;
	previoustime = NumToTime(Gal[p].SnapNum);
	newtime = NumToTime(Gal[p].SnapNum+1);
	deltaT = previoustime - newtime; //ROB (26-11-21): Is this not just the same as dt*STEPS?

	mass_checks(p,"Start of model_dust_yields.c",__LINE__);
	
	timestep_width = dt; //Width of current timestep in CODE UNITS (units cancel out when dividing by SFH bin width, sfh_dt) (12-04-12)
	TimeBin = (STEPS*(Halo[Gal[p].HaloNr].SnapNum-1.0))+nstep; //TimeBin = (STEPS*Gal[p].SnapNum)+nstep; //Bin in Yield tables corresponding to current timestep //TEST!: BRUNO: Snapnum would be +1 too low for a 'jumping' galaxy (14-11-13)

	/*//int checka=0;
	//if (Gal[p].DustColdGasClouds_elements[0] >= 0.) {
	//	checka = 1;
		printf("\nStart of model_dust_yields.c:\n");
		printf("Snapnum %i\n", Gal[p].SnapNum);
		for (ee=0; ee<NUM_ELEMENTS; ee++)
			printf("Gal[p].DustColdGasClouds_elements[%i] = %f\n", ee, Gal[p].DustColdGasClouds_elements[ee]);
		for (j=0;j<RNUM;j++)
			printf("Gal[p].DustColdGasCloudsRings_elements[%i][O] = %f\n", j, Gal[p].DustColdGasCloudsRings_elements[j][O_NUM]);
	//}*/

//**************************************************
//Dust grain destruction from supernova shock waves
//**************************************************

#ifdef DUST_DESTRUCTION
#ifdef H2_AND_RINGS
	//if (p == 0 && Gal[p].SnapNum < 16) printf("\nSnap = %i | Step = %i :\n", Gal[p].SnapNum, nstep);
	//printf("model_dust_yields.c\nGal[p].DustColdGasCloudsRings_elements[0][O_NUM] = %e\n", Gal[p].DustColdGasCloudsRings_elements[0][O_NUM]);
	DustColdGasCloudsRings_elements_test1 = Gal[p].DustColdGasCloudsRings_elements[0][O_NUM];
	for (j=0;j<RNUM;j++) {
		//printf("model_dust_yields.c\nGal[p].DustColdGasCloudsRings_elements[%i][O_NUM] = %e | Gal[p].MetalsColdGasRings[%i][0] = %e | Gal[p].MetalsColdGasRings[%i][1] = %e | Gal[p].MetalsColdGasRings[%i][2] = %e\n",
		//		j, Gal[p].DustColdGasCloudsRings_elements[j][O_NUM], j, Gal[p].MetalsColdGasRings[j][0], j, Gal[p].MetalsColdGasRings[j][1], j, Gal[p].MetalsColdGasRings[j][2]);

		total_dust_before=0.0;
		total_dust_after=0.0;
		if ((Gal[p].MetalsColdGasRings[j][0]+Gal[p].MetalsColdGasRings[j][1]+Gal[p].MetalsColdGasRings[j][2]) > 0.0) {
#else //H2_AND_RINGS
	if ((Gal[p].MetalsColdGas[0]+Gal[p].MetalsColdGas[1]+Gal[p].MetalsColdGas[2]) > 0.0) {
#endif //H2_AND_RINGS
		//For dust destruction we follow the prescription of McKee1989.
        //Calculate destruction timescale and destruction fraction: ROB: N.B. dust in ColdGas assumed to only be destroyed by SNe in disc (not bulge or ICL).
		/* ROB: (27-11-21): DiskSNIIRate_current_ts and DiskSNIaRate_current_ts only give the SNII rate in the DISC in the LAST TIMESTEP of that snapshot
		 * Gal[p].DiskSNIIRate, on the other hand, outputs the average SNII rate for the disc from the whole snapshot (similarly to SFR).
		*/
#ifdef H2_AND_RINGS
		R_SN = (DiskSNIIRate_current_ts[j] + DiskSNIaRate_current_ts[j]) / (UnitTime_in_s / SEC_PER_YEAR); //To convert from 1/code_time_units to 1/yr. NOTE: I have put the denominator there in brackets, so the correct conversion is UnitTime_in_s / SEC_PER_YEAR, not UnitTime_in_s * SEC_PER_YEAR. (09-02-22)
		//tdes = (Gal[p].ColdGasRings[j]*(1.0e10/Hubble_h))/(M_CLEARED * F_SN * R_SN);
		frac_Cb = (Gal[p].DustColdGasCloudsRings_elements[j][Cb_NUM]+Gal[p].DustColdGasDiffRings_elements[j][Cb_NUM])/(Gal[p].DustColdGasCloudsRings_elements[j][Cb_NUM]+Gal[p].DustColdGasCloudsRings_elements[j][Si_NUM]+Gal[p].DustColdGasDiffRings_elements[j][Cb_NUM]+Gal[p].DustColdGasDiffRings_elements[j][Si_NUM]); //Fraction of C+Si in dust that is C
		m_cleared = calc_cleared_mass(frac_Cb);
		tdes = (Gal[p].ColdGasRings[j]*(1.0e10/Hubble_h))/(m_cleared * F_SN * R_SN);
		//printf("frac_Cb = %f | m_cleared = %f | R_SN = %e | Mcold = %e | tdes = %e | tdes_old = %e\n", frac_Cb, m_cleared, R_SN, Gal[p].ColdGasRings[j]*(1.0e10/Hubble_h), tdes, (Gal[p].ColdGasRings[j]*(1.0e10/Hubble_h))/(M_CLEARED * F_SN * R_SN));
		//printf("Snapnum %i | Ring %i | (UnitTime_in_s * SEC_PER_YEAR) = %e | UnitTime_in_years = %f\n", Gal[p].SnapNum, j, (UnitTime_in_s / SEC_PER_YEAR), UnitTime_in_years);
		Gal[p].t_des[j] = tdes; //Storing the destruction timescale for each ring as a galaxy property, for output.
#else //H2_AND_RINGS
		R_SN = (DiskSNIIRate_current_ts + DiskSNIaRate_current_ts) / (UnitTime_in_s / SEC_PER_YEAR); //To convert from 1/code_time_units to 1/yr
		//tdes = (Gal[p].ColdGas*(1.0e10/Hubble_h))/(M_CLEARED * F_SN * R_SN);
		frac_Cb = (Gal[p].DustColdGasClouds_elements[Cb_NUM]+Gal[p].DustColdGasDiff_elements[j][Cb_NUM])/(Gal[p].DustColdGasClouds_elements[Cb_NUM]+Gal[p].DustColdGasClouds_elements[Si_NUM]+Gal[p].DustColdGasDiff_elements[Cb_NUM]+Gal[p].DustColdGasDiff_elements[Si_NUM]); //Fraction of C+Si in dust that is C
		m_cleared = calc_cleared_mass(frac_Cb);
		tdes = (Gal[p].ColdGas*(1.0e10/Hubble_h))/(m_cleared * F_SN * R_SN);
		Gal[p].t_des = tdes; //Note, this will be the destruction timescale of the last timestep of each snapshot, rather than the average across the whole snapshot. (09-02-22)
#endif //H2_AND_RINGS

		if (tdes != tdes)
			survive_frac = 1.0; //This accounts for cases where there is no dust to destroy, which causes tdes = -nan
			//printf("check1\n");}
		else
			survive_frac = exp(-dt*UnitTime_in_years/tdes);

		//We assume that the SNR will destroy equal amounts of dust in cold clouds and 
		//the diffuse medium, but all those will end up as diffuse gas, I guess.  
		//Then some will be reaccreted onto cold clouds.
	    //Simplest approximation is just to destroy the same fraction in each.
#ifdef FULL_DUST_RATES
		for (ee=0; ee<NUM_ELEMENTS; ee++) {
#ifdef H2_AND_RINGS
			total_dust_before += Gal[p].DustColdGasDiffRings_elements[j][ee] + Gal[p].DustColdGasCloudsRings_elements[j][ee];
#else //H2_AND_RINGS
			total_dust_before += Gal[p].DustColdGasDiff_elements[ee] + Gal[p].DustColdGasClouds_elements[ee];
#endif //H2_AND_RINGS
		}
#endif //FULL_DUST_RATES
		
		//Update dust elements due to SN dust destruction:
		for (ee=0; ee<NUM_ELEMENTS; ee++) {
#ifdef H2_AND_RINGS
			Gal[p].DustColdGasDiff_elements[ee] += Gal[p].DustColdGasDiffRings_elements[j][ee]*(survive_frac-1.); //This will subtract the original DustColdGasDiffRings_elements[j][ee] and add the new "survival" fraction. (08-02-22)
			Gal[p].DustColdGasClouds_elements[ee] += Gal[p].DustColdGasCloudsRings_elements[j][ee]*(survive_frac-1.);
			Gal[p].DustColdGasDiffRings_elements[j][ee] = Gal[p].DustColdGasDiffRings_elements[j][ee]*survive_frac;
	    	Gal[p].DustColdGasCloudsRings_elements[j][ee] = Gal[p].DustColdGasCloudsRings_elements[j][ee]*survive_frac;
#else //H2_AND_RINGS
	    	Gal[p].DustColdGasDiff_elements[ee] = Gal[p].DustColdGasDiff_elements[ee]*survive_frac;
	    	Gal[p].DustColdGasClouds_elements[ee] = Gal[p].DustColdGasClouds_elements[ee]*survive_frac;
#endif //H2_AND_RINGS
	    }
		//printf("Gal[p].DustColdGasCloudsRings_elements[%i][O_NUM] = %e | survive_frac = %e\n", j, Gal[p].DustColdGasCloudsRings_elements[j][O_NUM], survive_frac);
		//if (Gal[p].DustColdGasCloudsRings_elements[j][O_NUM] < 0.0) exit(0);
		//printf("Snapnum %i | Ring %i | ColdGasRings = %e | tot_ColdGasRings = %e | ColdGas = %e | R_SN = %e | tdes = %e | survive_frac = %f | Gal[p].t_des = %e\n",
		//		Gal[p].SnapNum, j, Gal[p].ColdGasRings[j], tot_ColdGasRings, Gal[p].ColdGas, R_SN, tdes, survive_frac, Gal[p].t_des);

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
	    for (ee=0; ee<NUM_ELEMENTS; ee++) {
#ifdef H2_AND_RINGS
	    	total_dust_after += Gal[p].DustColdGasDiffRings_elements[j][ee] + Gal[p].DustColdGasCloudsRings_elements[j][ee];
#else //H2_AND_RINGS
	    	total_dust_after += Gal[p].DustColdGasDiff_elements[ee] + Gal[p].DustColdGasClouds_elements[ee];
#endif //H2_AND_RINGS
	    }
#endif //FULL_DUST_RATES
	    mass_checks(p,"Dust from destruction: model_dust_yields.c",__LINE__);

#ifdef FULL_DUST_RATES
		Gal[p].DustColdGasRates[4] += (total_dust_before - total_dust_after)/(deltaT * UnitTime_in_years);
		//if (p == 0 && Gal[p].SnapNum < 16) printf("Ring = %i | total_dust_before = %e | total_dust_after = %e | Gal[p].DustColdGasRates[4] = %e\n", j, total_dust_before, total_dust_after, Gal[p].DustColdGasRates[4]);
#endif //FULL_DUST_RATES
		}
#ifdef H2_AND_RINGS
			//printf("Gal[p].DustColdGasCloudsRings_elements[%i][O_NUM] = %e\n", j, Gal[p].DustColdGasCloudsRings_elements[j][O_NUM]);
			DustColdGasCloudsRings_elements_test2 = Gal[p].DustColdGasCloudsRings_elements[0][O_NUM];
		}
#endif //H2_AND_RINGS

#ifdef DUST_HOTGAS //Needs to be outside the rings loop above.
		//*********************************************************************
		//Destruction of dust in the HotGas from thermal sputtering (15-03-22):
		//*********************************************************************
		//Thermal sputtering model from Triani+20:
		//NOTE: Need to check Hubble_h factors for all of these constants! (See h_params.h)
		HotGasDust = 0.0; //[Msun]
		for (ee=0; ee<NUM_ELEMENTS; ee++)
			HotGasDust += Gal[p].DustHotGas_elements[ee];

		if (Gal[p].HotGas > 0.0 && HotGasDust > 0.0) {
			/*
			//TEST: Destroying all dust in the HotGas, to return same result as not modelling it at all (but with modified DustColdGasRates and DustHotGasRates).
			frac_destroyed_HotGas = 1.0;
			for (ee=0; ee<NUM_ELEMENTS; ee++) {
#ifdef FULL_DUST_RATES
				Gal[p].DustHotGasRates += (Gal[p].DustHotGas_elements[ee] * frac_destroyed_HotGas)/(deltaT * UnitTime_in_years);
#endif //FULL_DUST_RATES
				Gal[p].DustHotGas_elements[ee] -= Gal[p].DustHotGas_elements[ee] * frac_destroyed_HotGas;
			}*/

			Tvir = 0.5 * MUMH * pow2(UnitVelocity_in_cm_per_s*Gal[p].Vvir) / BOLTZMANN; //[Kelvin] //ROB: Should prob be ~10^4 to ~10^6 K (see EAGLE analysis in Tumlinson+18)??
			HotGasDensity = (Gal[p].HotGas*UnitMass_in_g) / ((4./3.) * PI * pow3(Gal[p].Rvir*UnitLength_in_cm)); //[g/cm^3] //ROB: Should prob be ~10^-26 to ~10^-29 g/cm^3 (see EAGLE analysis in Tumlinson+18)??
			//Delta_grainRadius = (sputConst * (HotGasDensity/PROTONMASS) * (1./(pow(sputT0/Tvir,sputOmega)+1.))) * 1.e5 * (dt * UnitTime_in_s); //[micron] //Eqn.22 Triani+20 (except here, Delta_grainRadius is positive) //ERROR --> Multiple of 1.0e5 here is to convert cm to microns. (1 CM IS ACTUALLY 1e4 MICRONS) //HotGasDensity in g/cm^3, VirialTemp in K
			Delta_grainRadius = (sputConst * (HotGasDensity/PROTONMASS) * (1./(pow(sputT0/Tvir,sputOmega)+1.))) * 1.e4 * (dt * UnitTime_in_s); //[micron] //Eqn.22 Triani+20 (except here, Delta_grainRadius is positive) //Multiple of 1.0e4 here is to convert cm to microns. //HotGasDensity in g/cm^3, VirialTemp in K
			//printf("\nTvir = %e | T0 = %e | HotGasDensity = %e | Delta_grainRadius = %f\n", Tvir, sputT0, HotGasDensity, Delta_grainRadius); //if (p == 0 && Gal[p].SnapNum < 16)
			if (Delta_grainRadius < 0.0) {
				printf("PROBLEM: Delta_grainRadius = %f | Setting Delta_grainRadius = 0.0...\n", Delta_grainRadius);
				Delta_grainRadius = 0.0;
			}
			if (Delta_grainRadius > CharacGrainRadius)
				Delta_grainRadius = CharacGrainRadius;

			//ROB: (16-03-22): Should have a galaxy struct variable which stores current dust grain radius, for use here like this:
			//Gal[p].DustHotGas_grainRadius -= Delta_grainRadius;
			//tau_sput = sputCharacTime * (Gal[p].DustHotGas_grainRadius/CharacGrainRadius) * (CharacDensity/HotGasDensity) * ((sputT0/Tvir)^(sputOmega)+1.); //[Gyr] //Eqn.23 Triani+20

			tau_sput = sputCharacTime * ((CharacGrainRadius-Delta_grainRadius)/CharacGrainRadius) * (CharacDensity/HotGasDensity) * (pow(sputT0/Tvir,sputOmega)+1.); //[Gyr] //Eqn.23 Triani+20
			Gal[p].t_sput = tau_sput; //[Gyr]
			if (tau_sput == 0.0) //This is the case where all dust has been sputtered away (i.e. Delta_grainRadius >= Gal[p].DustHotGas_grainRadius)
				Delta_Msput = HotGasDust;
			else {
				Delta_Msput = (HotGasDust/(tau_sput/3.)) * ((dt * UnitTime_in_Megayears) / 1.e3); //[Msun] //Eqn.24 Triani+20 (except here, Delta_Msput is positive). Divide by 1.e3 to convert timestep width from Myr to Gyr.
				if (Delta_Msput < 0.0) {
					printf("PROBLEM: Delta_Msput = %e | Setting Delta_Msput = 0.0...\n", Delta_Msput);
					Delta_Msput = 0.0;
				}
				if (Delta_Msput > HotGasDust)
					Delta_Msput = HotGasDust;
			}
			//printf("\nTvir = %e [K] | T0 = %e [K] | HotGasDensity = %e [g/cm^3] | Delta_grainRadius = %f [micron] |\ntau_sput = %e [Gyr] | HotGasDust = %e [Msun] | Delta_Msput = %e [Msun]\n", Tvir, sputT0, HotGasDensity, Delta_grainRadius, tau_sput, HotGasDust, Delta_Msput);

			for (ee=0; ee<NUM_ELEMENTS; ee++) {
				if (Delta_Msput == HotGasDust)
					Gal[p].DustHotGas_elements[ee] = 0.0; //The case where all dust is sputtered away (doing it this way prevents very small precision errors)
				else
					Gal[p].DustHotGas_elements[ee] -= Delta_Msput * (Gal[p].DustHotGas_elements[ee]/HotGasDust); //This will subtract the mass of HotGas dust destroyed by sputtering evenly from each element.
			}

#ifdef FULL_DUST_RATES
				Gal[p].DustHotGasRates += Delta_Msput/(deltaT * UnitTime_in_years); //[1/yr]
#endif //FULL_DUST_RATES
			//printf("Tvir = %e | HotGasDensity = %e | Delta_grainRadius = %f | tau_sput = %e | HotGasDust_elements = [%e,%e,%e...%e,%e,%e] | Delta_Msput = %e\n",
			//		Tvir, HotGasDensity, Delta_grainRadius, tau_sput, Gal[p].DustHotGas_elements[2], Gal[p].DustHotGas_elements[4], Gal[p].DustHotGas_elements[5], Gal[p].DustHotGas_elements[7], Gal[p].DustHotGas_elements[9], Gal[p].DustHotGas_elements[10], Delta_Msput);
		}
#endif //DUST_HOTGAS
#endif //DUST_DESTRUCTION

	/*//if (checka == 1.) {
		printf("***** After dust destruction:\n");
		printf("Snapnum %i\n", Gal[p].SnapNum);
		for (ee=0; ee<NUM_ELEMENTS; ee++)
			printf("Gal[p].DustColdGasClouds_elements[%i] = %f\n", ee, Gal[p].DustColdGasClouds_elements[ee]);
		for (j=0;j<RNUM;j++)
			printf("Gal[p].DustColdGasCloudsRings_elements[%i][O] = %f\n", j, Gal[p].DustColdGasCloudsRings_elements[j][O_NUM]);
	//}*/

  //*****
  //(28-01-22): partition_gas_and_dust_elements() was moved here from just above grain growth. Didn't make any difference to the outputs, but is more consistent when checking the diff/cloud mass fractions for AGB/SN dust production.
  //(02-02-22): Note that partition_gas_and_dust_elements() now has it's own internal loop over rings, so make sure it is only called outside of a loop over rings in this function (i.e. here at the top).
  //(09-03-22): partition_gas_and_dust_elements() has now been moved to the end of update_h2fraction(), so that the cloud and diffuse gas components are always in sync with the current H2 fraction.

  //partition_gas_and_dust_elements(p);
  mass_checks(p,"After partition_gas_and_dust_elements: model_dust_yields.c",__LINE__);

	/*//if (checka == 1.) {
		printf("***** After partition_gas_and_dust_elements():\n");
		printf("Snapnum %i\n", Gal[p].SnapNum);
		for (ee=0; ee<NUM_ELEMENTS; ee++)
			printf("Gal[p].DustColdGasClouds_elements[%i] = %f\n", ee, Gal[p].DustColdGasClouds_elements[ee]);
		for (j=0;j<RNUM;j++)
			printf("Gal[p].DustColdGasCloudsRings_elements[%i][O] = %f\n", j, Gal[p].DustColdGasCloudsRings_elements[j][O_NUM]);
	//}*/

  for (i=0;i<=Gal[p].sfh_ibin;i++) { //LOOP OVER SFH BINS

	  //ROB: (02-02-22): For incorporating rings, I've added a big loop over rings starting here (as is done in model_yields.c).
#ifdef H2_AND_RINGS
	  DustColdGasCloudsRings_elements_test2a = Gal[p].DustColdGasCloudsRings_elements[0][O_NUM];
	  for(j=0;j<RNUM;j++) {
		  DustColdGasCloudsRings_elements_test2b = Gal[p].DustColdGasCloudsRings_elements[0][O_NUM];

#endif

//*****************************************
//DUST ENRICHMENT FROM AGB DISC STARS INTO COLD GAS:
//*****************************************

#ifdef DUST_AGB
#ifdef H2_AND_RINGS
	  if ((Gal[p].sfh_DiskMassRings[j][i] > 0.0) && (Gal[p].MetalsColdGasRings[j][2] > 0.0)) {
		  DiskSFR = Gal[p].sfh_DiskMassRings[j][i]/Gal[p].sfh_dt[i];
#else
	  if ((Gal[p].sfh_DiskMass[i] > 0.0) && (Gal[p].MetalsColdGas[2] > 0.0)) {
		  DiskSFR = Gal[p].sfh_DiskMass[i]/Gal[p].sfh_dt[i];
#endif
		  //pre-calculations to speed up the code
		  DiskSFR_physical_units = DiskSFR * (1.0e10/Hubble_h); //Note: This is NOT in physical units (i.e. NOT in Msun/yr, but in Msun/[code_time_units]). But this is ok, as code_time_units cancel out when multiplying by timestep_width to get 'step_width_times_DiskSFR_physical_units' on the line below ('DiskSFR_physical_units' is never used itself).
		  step_width_times_DiskSFR_physical_units = timestep_width * DiskSFR_physical_units; //ROB: This is the same as DiskSFRxStep_Phys in moedl_yields.c. (17-01-22)

		/* ROB: The Zi and Zi_disp variables are not needed in model_dust_yields.c, as they calculated in model_yields.c and stored in the global variables Zi_saved and Zi_disp_saved instead. (04-01-22)
		 * 	// Note: This approach is ok, but only if update_dust_mass() is called after update_yields_and_return_mass() in main.c, and both are looking at the same galaxy.
		 * 	// Alternative approach would be to run find_initial_metallicity_dust() here, to re-calculate Zi and Zi_disp again, independently of model_yields.c.
		 */
    	
		  //interpolates yields from lookup tables we produced in dust_yield_integrals.c
		  for (k=0;k<AGB_DUST_TYPE_NUM;k++) {
#ifdef H2_AND_RINGS
	    	NormAGBDustYieldRate_actual[k] = NormAGBDustYieldRate[TimeBin][i][Zi_saved[j][i]][k] + ((NormAGBDustYieldRate[TimeBin][i][Zi_saved[j][i]+1][k] - NormAGBDustYieldRate[TimeBin][i][Zi_saved[j][i]][k])*Zi_disp_saved[j][i]);
#else
	    	NormAGBDustYieldRate_actual[k] = NormAGBDustYieldRate[TimeBin][i][Zi_saved[i]][k] + ((NormAGBDustYieldRate[TimeBin][i][Zi_saved[i]+1][k] - NormAGBDustYieldRate[TimeBin][i][Zi_saved[i]][k])*Zi_disp_saved[i]);
#endif
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

		//DustColdGasCloudsRings_elements_test2c = Gal[p].DustColdGasCloudsRings_elements[0][O_NUM];

		//Element Conversion -----------------------------------------------------------------------------------
		//Conversion of dust species (i.e. Ferrosilite) into Actual elements to store
		//in correct arrays (i.e. Forsterite -> Mg/Si/O)
		//All the following conversions are done by mass fraction
		//Corrections added at places to avoid dust masses going beyond the mass in gas-phase elements left to form dust from
		//(given that the newly-formed elements have already been added to the gas in model_yields.c).
		for (ee=0; ee<NUM_ELEMENTS; ee++) {
#ifdef H2_AND_RINGS
			H2frac = Gal[p].H2fractionRings[j];
#else
			H2frac = Gal[p].H2fraction;
#endif //H2_AND_RINGS
#ifndef MAINELEMENTS
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
#endif //MAINELEMENTS
			//else if (ee == O_NUM) { //O
			if (ee == O_NUM) { //O
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

#ifdef H2_AND_RINGS
			//Check how much gas there is actually available to form dust:
			ColdGasDiff_avail = Gal[p].ColdGasDiffRings_elements[j][ee] - Gal[p].DustColdGasDiffRings_elements[j][ee]; //Total diffuse gas available to form dust
			ColdGasClouds_avail = Gal[p].ColdGasCloudsRings_elements[j][ee] - Gal[p].DustColdGasCloudsRings_elements[j][ee]; //Total cloud gas available to form dust
			AGBAllElementsDiff_avail = (1.0-fwind_AGB) * DiskAGBAllElements_ts[j][ee] * (1.0-Gal[p].H2fractionRings[j]); //Newly-ejected element mass into diffuse gas available to form dust
			AGBAllElementsClouds_avail = (1.0-fwind_AGB) * DiskAGBAllElements_ts[j][ee] * Gal[p].H2fractionRings[j]; //Newly-ejected element mass into clouds available to form dust
			//DustColdGasCloudsRings_elements_test3 = Gal[p].DustColdGasCloudsRings_elements[0][O_NUM];
			/*if (New_dust_clouds < 0.0 || ColdGasClouds_avail < 0.0 || AGBAllElementsClouds_avail < 0.0) {
				printf("\nGal[p].DustColdGasCloudsRings_elements[%i][%i] = %e | Gal[p].ColdGasCloudsRings_elements[%i][%i] = %e | New_dust_clouds = %e | ColdGasClouds_avail = %e\n",
						j, ee, Gal[p].DustColdGasCloudsRings_elements[j][ee], j, ee, Gal[p].ColdGasCloudsRings_elements[j][ee], New_dust_clouds, ColdGasClouds_avail);
				printf("AGBAllElementsClouds_avail = %e | fwind_AGB = %e | DiskAGBAllElements_ts[%i][%i] = %e | Gal[p].H2fractionRings[%i] = %e\n",
						AGBAllElementsClouds_avail, fwind_AGB, j, ee, DiskAGBAllElements_ts[j][ee], j, Gal[p].H2fractionRings[j]);
			}*/
			//Add newly-formed dust to ColdGas:
			Gal[p].DustColdGasDiffRings_elements[j][ee]  += min(New_dust_diff, min(ColdGasDiff_avail, AGBAllElementsDiff_avail));
			Gal[p].DustColdGasCloudsRings_elements[j][ee]  += min(New_dust_clouds, min(ColdGasClouds_avail, AGBAllElementsClouds_avail));
			//Gal[p].DustColdGasDiffRings_elements[j][ee]  += max(0.0, min(New_dust_diff, min(ColdGasDiff_avail, AGBAllElementsDiff_avail)));
			//Gal[p].DustColdGasCloudsRings_elements[j][ee]  += max(0.0, min(New_dust_clouds, min(ColdGasClouds_avail, AGBAllElementsClouds_avail)));

			//Adding the same small increments to the global quantities within this loop over rings should work here: (02-02-22):
			Gal[p].DustColdGasDiff_elements[ee]  += min(New_dust_diff, min(ColdGasDiff_avail, AGBAllElementsDiff_avail));
			Gal[p].DustColdGasClouds_elements[ee]  += min(New_dust_clouds, min(ColdGasClouds_avail, AGBAllElementsClouds_avail));
			//Gal[p].DustColdGasDiff_elements[ee]  += max(0.0, min(New_dust_diff, min(ColdGasDiff_avail, AGBAllElementsDiff_avail)));
			//Gal[p].DustColdGasClouds_elements[ee]  += max(0.0, min(New_dust_clouds, min(ColdGasClouds_avail, AGBAllElementsClouds_avail)));
			/*if (New_dust_clouds < 0.0 || ColdGasClouds_avail < 0.0 || AGBAllElementsClouds_avail < 0.0) {
				printf("Gal[p].DustColdGasCloudsRings_elements[%i][%i] = %e | Gal[p].ColdGasCloudsRings_elements[%i][%i] = %e\n",
						j, ee, Gal[p].DustColdGasCloudsRings_elements[j][ee], j, ee, Gal[p].ColdGasCloudsRings_elements[j][ee]);
			}*/
#else
			//Check how much gas there is actually available to form dust:
			ColdGasDiff_avail = Gal[p].ColdGasDiff_elements[ee] - Gal[p].DustColdGasDiff_elements[ee]; //Total diffuse gas available to form dust
			ColdGasClouds_avail = Gal[p].ColdGasClouds_elements[ee] - Gal[p].DustColdGasClouds_elements[ee]; //Total cloud gas available to form dust
			AGBAllElementsDiff_avail = (1.0-fwind_AGB) * DiskAGBAllElements_ts[ee] * (1.0-Gal[p].H2fraction); //Newly-ejected element mass into diffuse gas available to form dust
			AGBAllElementsClouds_avail = (1.0-fwind_AGB) * DiskAGBAllElements_ts[ee] * Gal[p].H2fraction; //Newly-ejected element mass into clouds available to form dust
			//Add newly-formed dust to ColdGas:
			Gal[p].DustColdGasDiff_elements[ee]  += min(New_dust_diff, min(ColdGasDiff_avail, AGBAllElementsDiff_avail));
			Gal[p].DustColdGasClouds_elements[ee]  += min(New_dust_clouds, min(ColdGasClouds_avail, AGBAllElementsClouds_avail));
#endif //H2_AND_RINGS
		}
	} //if ( (Gal[p].sfh_DiskMass[i] > 0.0) && (Gal[p].MetalsColdGas[2] > 0.0) )
    
    mass_checks(p,"Dust from AGBs: model_dust_yields.c",__LINE__);
#endif //DUST_AGB


//*****************************************
//DUST ENRICHMENT FROM SNII FROM DISK STARS INTO COLD GAS:
//*****************************************

#ifdef DUST_SNII
#ifdef H2_AND_RINGS
    //DustColdGasCloudsRings_elements_test4 = Gal[p].DustColdGasCloudsRings_elements[0][O_NUM];
    if ((Gal[p].sfh_DiskMassRings[j][i] > 0.0) && (Gal[p].MetalsColdGasRings[j][0] > 0.0)) {
#ifdef FULL_DUST_RATES
    	//This is estimating the dust in various compounds (say silicates) using the amount of the particular element (say silicon) produced in a process
#ifndef MAINELEMENTS
    	Gal[p].DustColdGasRates[1] += (SNII_prevstep_Cold_Si[j][i] * eta_SNII_Sil * A_Sil_dust/A_Si)/(deltaT * UnitTime_in_years);
		Gal[p].DustColdGasRates[1] += (SNII_prevstep_Cold_Si[j][i] * eta_SNII_SiC * A_SiC_dust/A_Si)/(deltaT * UnitTime_in_years);
		Gal[p].DustColdGasRates[1] += (SNII_prevstep_Cold_Cb[j][i] * eta_SNII_Cb  * A_Cb_dust/A_Cb) /(deltaT * UnitTime_in_years);
#endif //MAINELEMENTS
		Gal[p].DustColdGasRates[1] += (SNII_prevstep_Cold_Fe[j][i] * eta_SNII_Fe  * A_Fe_dust/A_Fe )/(deltaT * UnitTime_in_years);
#endif //FULL_DUST_RATES
#else //H2_AND_RINGS
    if ((Gal[p].sfh_DiskMass[i] > 0.0) && (Gal[p].MetalsColdGas[0] > 0.0)) {
#ifdef FULL_DUST_RATES
#ifndef MAINELEMENTS
    	Gal[p].DustColdGasRates[1] += (SNII_prevstep_Cold_Si[i] * eta_SNII_Sil * A_Sil_dust/A_Si)/(deltaT * UnitTime_in_years);
		Gal[p].DustColdGasRates[1] += (SNII_prevstep_Cold_Si[i] * eta_SNII_SiC * A_SiC_dust/A_Si)/(deltaT * UnitTime_in_years);
		Gal[p].DustColdGasRates[1] += (SNII_prevstep_Cold_Cb[i] * eta_SNII_Cb  * A_Cb_dust/A_Cb) /(deltaT * UnitTime_in_years);
#endif //MAINELEMENTS
		Gal[p].DustColdGasRates[1] += (SNII_prevstep_Cold_Fe[i] * eta_SNII_Fe  * A_Fe_dust/A_Fe )/(deltaT * UnitTime_in_years);
#endif //FULL_DUST_RATES
#endif //H2_AND_RINGS

		//Create dust (based on the prescription of Zhukovska2008)---------------------------
		//SNII_prevstep_x is calculated in model_yields.c
		//It is the amount of a specific metal (i.e. Si) returned to gas phase from SFH_BIN i. //produced in the last timestep

		//ROB: Here, the mass of newly-formed dust is determined by the mass of the "key element" that has just been returned by SNe.
		//So, e.g. if the key element makes up 25% of the total dust mass, then the mass of dust formed is 1/0.25 = 4 times greater than the mass of the key element formed * the conversion efficiency (eta).
		//This is why the A_dust/A_key_element terms are included here.
#ifndef MAINELEMENTS
#ifdef H2_AND_RINGS
		double Dust_Silicates = SNII_prevstep_Cold_Si[j][i] * eta_SNII_Sil * A_Sil_dust/A_Si;
		double Dust_Iron      = SNII_prevstep_Cold_Fe[j][i] * eta_SNII_Fe  * A_Fe_dust/A_Fe;
		double Dust_SiC	      = SNII_prevstep_Cold_Si[j][i] * eta_SNII_SiC * A_SiC_dust/A_Si;
		double Dust_Carbon    = SNII_prevstep_Cold_Cb[j][i] * eta_SNII_Cb  * A_Cb_dust/A_Cb;
#else
		double Dust_Silicates = SNII_prevstep_Cold_Si[i] * eta_SNII_Sil * A_Sil_dust/A_Si;
		double Dust_Iron      = SNII_prevstep_Cold_Fe[i] * eta_SNII_Fe  * A_Fe_dust/A_Fe;
		double Dust_SiC	      = SNII_prevstep_Cold_Si[i] * eta_SNII_SiC * A_SiC_dust/A_Si;
		double Dust_Carbon    = SNII_prevstep_Cold_Cb[i] * eta_SNII_Cb  * A_Cb_dust/A_Cb;
#endif //H2_AND_RINGS
#else //MAINELEMENTS
#ifdef H2_AND_RINGS
		double Dust_Silicates = 0.0;
		double Dust_Iron      = SNII_prevstep_Cold_Fe[j][i] * eta_SNII_Fe  * A_Fe_dust/A_Fe;
		double Dust_SiC	      = 0.0;
		double Dust_Carbon    = 0.0;
#else
		double Dust_Silicates = 0.0;
		double Dust_Iron      = SNII_prevstep_Cold_Fe[i] * eta_SNII_Fe  * A_Fe_dust/A_Fe;
		double Dust_SiC	      = 0.0;
		double Dust_Carbon    = 0.0;
#endif //H2_AND_RINGS
#endif //MAINELEMENTS
            
		//Element conversion -----------------------------------------------------------------
		//Conversion of dust species (i.e. Silicates) into Actual elements to store
		//in correct arrays (i.e. Silicates -> Mg/Si/Fe/O)
		//All the following conversions are done by mass fraction
		//SNII Silicates -------------------
		//ROB: (04-01-22) The mass fractions used here don't make sense to me.
		//Assuming that the silicate mass is just olivines + pyroxenes (and f_ol = 0.32), I get A_Sil_dust = 121.62, close to the quoted value used here (in h_params.h).
		//But, e.g. A_Si / A_Sil_dust = 0.2310, whereas 0.210432 is used here...
		for (ee=0; ee<NUM_ELEMENTS; ee++) {
#ifdef H2_AND_RINGS
			H2frac = Gal[p].H2fractionRings[j];
#else
			H2frac = Gal[p].H2fraction;
#endif //H2_AND_RINGS
#ifndef MAINELEMENTS
				if (ee == Cb_NUM) { //Cb
					New_dust_diff = ((Dust_SiC * SILICONCARBIDE_Cb_FRAC) + (Dust_Carbon * 1.0)) * (1.0 - H2frac);
					New_dust_clouds = ((Dust_SiC * SILICONCARBIDE_Cb_FRAC) + (Dust_Carbon * 1.0)) * H2frac;
				}
				else if (ee == Si_NUM) { //Si
					New_dust_diff = ((Dust_Silicates * SILICATES_Si_FRAC) + (Dust_SiC * SILICONCARBIDE_Si_FRAC)) * (1.0 - H2frac);
					New_dust_clouds = ((Dust_Silicates * SILICATES_Si_FRAC) + (Dust_SiC * SILICONCARBIDE_Si_FRAC)) * H2frac;
				}
#endif //MAINELEMENTS
				//else if (ee == O_NUM) { //O
				if (ee == O_NUM) { //O
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
#ifdef H2_AND_RINGS
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
#else //H2_AND_RINGS
				//Check how much gas there is actually available to form dust:
				ColdGasDiff_avail = Gal[p].ColdGasDiff_elements[ee] - Gal[p].DustColdGasDiff_elements[ee]; //Total diffuse gas available to form dust
				ColdGasClouds_avail = Gal[p].ColdGasClouds_elements[ee] - Gal[p].DustColdGasClouds_elements[ee]; //Total cloud gas available to form dust
				SNIIAllElementsDiff_avail = (1.0-fwind_SNII) * DiskSNIIAllElements_ts[ee] * (1.0-Gal[p].H2fraction); //Newly-ejected element mass into diffuse gas available to form dust
				SNIIAllElementsClouds_avail = (1.0-fwind_SNII) * DiskSNIIAllElements_ts[ee] * Gal[p].H2fraction; //Newly-ejected element mass into clouds available to form dust
				//Add newly-formed dust to ColdGas:
				Gal[p].DustColdGasDiff_elements[ee]  += min(New_dust_diff, min(ColdGasDiff_avail, SNIIAllElementsDiff_avail));
				Gal[p].DustColdGasClouds_elements[ee]  += min(New_dust_clouds, min(ColdGasClouds_avail, SNIIAllElementsClouds_avail));
#endif //H2_AND_RINGS
			} //for (ee=0; ee<NUM_ELEMENTS; ee++)

			mass_checks(p,"Dust from SNe-II: model_dust_yields.c",__LINE__);
	}//if ((Gal[p].sfh_DiskMass[i] > 0.0) && (Gal[p].MetalsColdGas[0] > 0.0))
#endif //DUST_SNII
	
	
//*****************************************
//DUST ENRICHMENT FROM SNIA FROM DISK STARS INTO COLD GAS:
//*****************************************
#ifdef DUST_SNIA
#ifdef H2_AND_RINGS
    if ((Gal[p].sfh_DiskMassRings[j][i] > 0.0) && (Gal[p].MetalsColdGasRings[j][1] > 0.0)) {
		double Dust_Iron = SNIa_prevstep_Cold_Fe[j][i] * eta_SNIa_Fe  * A_Fe_dust/A_Fe;
        #ifdef FULL_DUST_RATES
		    Gal[p].DustColdGasRates[2]  += (SNIa_prevstep_Cold_Fe[j][i] * eta_SNIa_Fe  * A_Fe_dust/A_Fe)/(deltaT * UnitTime_in_years);
        #endif
#else
    if ((Gal[p].sfh_DiskMass[i] > 0.0) && (Gal[p].MetalsColdGas[1] > 0.0)) {
		double Dust_Iron = SNIa_prevstep_Cold_Fe[i] * eta_SNIa_Fe  * A_Fe_dust/A_Fe;
        #ifdef FULL_DUST_RATES		
		    Gal[p].DustColdGasRates[2]  += (SNIa_prevstep_Cold_Fe[i] * eta_SNIa_Fe  * A_Fe_dust/A_Fe)/(deltaT * UnitTime_in_years);
        #endif	
#endif //H2_AND_RINGS
		
#ifdef H2_AND_RINGS
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
#else //H2_AND_RINGS
		Gal[p].DustColdGasDiff_elements[Fe_NUM] += min(Dust_Iron * 1.0 * (1.0-Gal[p].H2fraction),
												   min(Gal[p].ColdGasDiff_elements[Fe_NUM] - Gal[p].DustColdGasDiff_elements[Fe_NUM],
													   (1.0-fwind_SNIa) * DiskSNIaAllElements_ts[Fe_NUM] * (1.0-Gal[p].H2fraction)));
		Gal[p].DustColdGasClouds_elements[Fe_NUM] += min(Dust_Iron * 1.0 * Gal[p].H2fraction,
													 min(Gal[p].ColdGasClouds_elements[Fe_NUM] - Gal[p].DustColdGasClouds_elements[Fe_NUM],
														 (1.0-fwind_SNIa) * DiskSNIaAllElements_ts[Fe_NUM] * Gal[p].H2fraction));
#endif //H2_AND_RINGS

		mass_checks(p,"Dust from SNe-Ia: model_dust_yields.c",__LINE__);
		
	}
#endif //DUST_SNIA

#ifdef H2_AND_RINGS
	  } //for(j=0;j<RNUM;j++)
#endif

} //loop over SFH bins

	/*//if (checka == 1) {
		printf("***** After AGB, SN-II, and SN-Ia dust production:\n");
		printf("Snapnum %i\n", Gal[p].SnapNum);
		for (int itee=0; itee<NUM_ELEMENTS; itee++)
			printf("Gal[p].DustColdGasClouds_elements[%i] = %f\n", itee, Gal[p].DustColdGasClouds_elements[itee]);
		for (int itjj=0;itjj<RNUM;itjj++)
			printf("Gal[p].DustColdGasCloudsRings_elements[%i][O] = %f\n", itjj, Gal[p].DustColdGasCloudsRings_elements[itjj][O_NUM]);
	//}*/

//*****************************************
//DUST GRAIN GROWTH INSIDE MOLECULAR CLOUDS
//*****************************************
//partition_gas_and_dust_elements(p);

#ifdef DUST_GROWTH
#ifdef H2_AND_RINGS
for(j=0;j<RNUM;j++) {
   	if (((Gal[p].MetalsColdGasRings[j][0]+Gal[p].MetalsColdGasRings[j][1]+Gal[p].MetalsColdGasRings[j][2])>0.0) && (Gal[p].ColdGasRings[j] > 0.0)) {
#else
if (((Gal[p].MetalsColdGas[0]+Gal[p].MetalsColdGas[1]+Gal[p].MetalsColdGas[2])>0.0) && (Gal[p].ColdGas > 0.0)) {
#endif //H2_AND_RINGS

        mass_checks(p,"Dust from growth: model_dust_yields.c",__LINE__);
        
        //The number of CO molecules that can be produced from available Carbon and Oxygen
        //in the clouds. By default, assuming only 30% C is locked up as CO in clouds.
        //ROB (10-01-22): This sets the number of CO molecules in the gas phase of the ISM to the minimum of either:
        //		(a) number of carbon atoms in gas phase
        //		(b) Cmax_CO (= 0.3): Fraction of total carbon atoms in gas and dust (this is the max amount of carbon allowed in CO, set in input file)
        //		(c) number of oxygen atoms in gas phase
#ifndef MAINELEMENTS
#ifdef H2_AND_RINGS
        num_CO = min((Gal[p].ColdGasCloudsRings_elements[j][Cb_NUM]-Gal[p].DustColdGasCloudsRings_elements[j][Cb_NUM])/A_Cb,
                     min(Gal[p].ColdGasCloudsRings_elements[j][Cb_NUM]/A_Cb*Cmax_CO, (Gal[p].ColdGasCloudsRings_elements[j][O_NUM]-Gal[p].DustColdGasCloudsRings_elements[j][O_NUM])/A_O));

        //Mass of Carbon and Oxygen available in the gas+dust of clouds for grain growth:
        //double Cb_clouds = Gal[p].ColdGasClouds_elements[Cb_NUM] - num_CO*A_Cb; //ROB: Note, Cb_clouds isn't used anywhere... (17-01-22)
        O_clouds = Gal[p].ColdGasCloudsRings_elements[j][O_NUM] - num_CO*A_O;
#else
        num_CO = min((Gal[p].ColdGasClouds_elements[Cb_NUM]-Gal[p].DustColdGasClouds_elements[Cb_NUM])/A_Cb,
					min(Gal[p].ColdGasClouds_elements[Cb_NUM]/A_Cb*Cmax_CO, (Gal[p].ColdGasClouds_elements[O_NUM]-Gal[p].DustColdGasClouds_elements[O_NUM])/A_O));
        O_clouds = Gal[p].ColdGasClouds_elements[O_NUM] - num_CO*A_O;
#endif //H2_AND_RINGS
#else //MAINELEMENTS
        //NOTE: No carbon is tracked when MAINELEMENTS is on, so no CO is calculable:
        num_CO = 0.0;
#ifdef H2_AND_RINGS
        O_clouds = Gal[p].ColdGasCloudsRings_elements[j][O_NUM];
#else
        O_clouds = Gal[p].ColdGasClouds_elements[O_NUM];
#endif //H2_AND_RINGS
#endif //MAINELEMENTS
        
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
         * I have also re-written these "number of atoms per typical olivine/pyroxene silicate molecule" variables for O, Mg, Si, and Fe calculations accordingly.
         */
#ifndef MAINELEMENTS
        float N_O = OLIVINE_NUMFRAC*4. + (1.-OLIVINE_NUMFRAC)*3.; //number of O atoms per typical silicate molecule
		float N_Mg = OLIVINE_NUMFRAC*2.*OLIVINE_Mg_NUMFRAC + (1.-OLIVINE_NUMFRAC)*PYROXENE_Mg_NUMFRAC;
		float N_Si = OLIVINE_NUMFRAC*1. + (1.-OLIVINE_NUMFRAC)*1.;
		float N_Fe = OLIVINE_NUMFRAC*2.*(1-OLIVINE_Mg_NUMFRAC) + (1.-OLIVINE_NUMFRAC)*(1.-PYROXENE_Mg_NUMFRAC);

		//ROB: This is a nested min() function to find number of olivine/pyroxene silicate molecules formed, given the mass available of the least-abundant constituent element in the gas phase of the molecular clouds:
#ifdef H2_AND_RINGS
		float num_Silicates = max(0., min((O_clouds-Gal[p].DustColdGasCloudsRings_elements[j][O_NUM])/(N_O*A_O),
        								 min((Gal[p].ColdGasCloudsRings_elements[j][Mg_NUM]-Gal[p].DustColdGasCloudsRings_elements[j][Mg_NUM])/(N_Mg*A_Mg),
        									min((Gal[p].ColdGasCloudsRings_elements[j][Si_NUM]-Gal[p].DustColdGasCloudsRings_elements[j][Si_NUM])/(N_Si*A_Si),
        									   (Gal[p].ColdGasCloudsRings_elements[j][Fe_NUM]-Gal[p].DustColdGasCloudsRings_elements[j][Fe_NUM])/(N_Fe*A_Fe)))));
#else
		float num_Silicates = max(0., min((O_clouds-Gal[p].DustColdGasClouds_elements[O_NUM])/(N_O*A_O),
										 min((Gal[p].ColdGasClouds_elements[Mg_NUM]-Gal[p].DustColdGasClouds_elements[Mg_NUM])/(N_Mg*A_Mg),
											min((Gal[p].ColdGasClouds_elements[Si_NUM]-Gal[p].DustColdGasClouds_elements[Si_NUM])/(N_Si*A_Si),
											   (Gal[p].ColdGasClouds_elements[Fe_NUM]-Gal[p].DustColdGasClouds_elements[Fe_NUM])/(N_Fe*A_Fe)))));
#endif //H2_AND_RINGS
#else //MAINELEMENTS
		//For simplicity, we assume there are no silicates at all (as Si is not tracked) when MAINELEMENTS is on: (ROB 07-02-22):
		float num_Silicates = 0.0;
#endif //MAINELEMENTS
        
        //*********************************************
        //Iron Oxide:
        //*********************************************
        /* ROB (10-01-22): I didn't like the way the number of iron oxide molecules was calculated here.
         * There is a hidden normalisation constant of N_Fe_IronOxide floating around, because N_Oxide is actually the number ratio of N_O_IronOxide/N_Fe_IronOxide.
         * Therefore, I have re-written this bit to be exactly the same as for silicates above, with variables like N_O_IronOxide and num_iron_oxide properly
         * representing the unnormalised "number" of molecules.
         */
        //float N_O_IronOxide = (HEMATITE_NUMFRAC*4. + (1.-HEMATITE_NUMFRAC)*3.);
        //float N_Fe_IronOxide = (HEMATITE_NUMFRAC*3. + (1.-HEMATITE_NUMFRAC)*2.);
        float N_O_IronOxide = (HEMATITE_NUMFRAC*3. + (1.-HEMATITE_NUMFRAC)*4.); //First term is from hematite (Fe2 O3), second term is from magnetite [(Fe^2+) (Fe^3+)2 O4].
        float N_Fe_IronOxide = (HEMATITE_NUMFRAC*2. + (1.-HEMATITE_NUMFRAC)*3.); //First term is from hematite (Fe2 O3), second term is from magnetite [(Fe^2+) (Fe^3+)2 O4].

        //*********************************************
        //Calculate the f_max values for C and O (i.e. the maximum condensation fractions):
#ifdef H2_AND_RINGS
        float num_iron_oxide = max(0., min((O_clouds-num_Silicates*N_O*A_O-Gal[p].DustColdGasCloudsRings_elements[j][O_NUM])/(N_O_IronOxide*A_O),
										  (Gal[p].ColdGasCloudsRings_elements[j][Fe_NUM]-num_Silicates*N_Fe*A_Fe-Gal[p].DustColdGasCloudsRings_elements[j][Fe_NUM])/(N_Fe_IronOxide*A_Fe)));
        float f_O_max =  (num_Silicates*N_O*A_O + num_iron_oxide*N_O_IronOxide*A_O + Gal[p].DustColdGasCloudsRings_elements[j][O_NUM])/Gal[p].ColdGasCloudsRings_elements[j][O_NUM]; //Mass fraction of all oxygen in clouds that is in dust (and not CO, by virtue of the fact that O_clouds is used when calculating num_Silicates and num_IronOxide).
#else
        float num_iron_oxide = max(0., min((O_clouds-num_Silicates*N_O*A_O-Gal[p].DustColdGasClouds_elements[O_NUM])/(N_O_IronOxide*A_O),
										  (Gal[p].ColdGasClouds_elements[Fe_NUM]-num_Silicates*N_Fe*A_Fe-Gal[p].DustColdGasClouds_elements[Fe_NUM])/(N_Fe_IronOxide*A_Fe)));
        float f_O_max =  (num_Silicates*N_O*A_O + num_iron_oxide*N_O_IronOxide*A_O + Gal[p].DustColdGasClouds_elements[O_NUM])/Gal[p].ColdGasClouds_elements[O_NUM];
#endif //H2_AND_RINGS
#ifndef MAINELEMENTS
#ifdef H2_AND_RINGS
        float f_Cb_max =  1.0 - num_CO*A_Cb/Gal[p].ColdGasCloudsRings_elements[j][Cb_NUM]; //Mass fraction of all carbon in clouds not in CO (i.e. max available for dust formation)
        //ROB: Why isn't Gal[p].DustColdGasCloudsRings_elements[j][Cb_NUM] also considered here? (07-03-22)
#else
        float f_Cb_max =  1.0 - num_CO*A_Cb/Gal[p].ColdGasClouds_elements[Cb_NUM]; //Mass fraction of all carbon in clouds not in CO
#endif //H2_AND_RINGS
#else //MAINELEMENTS
        float f_Cb_max =  1.0;
#endif //MAINELEMENTS

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
#ifdef H2_AND_RINGS
#ifndef MAINELEMENTS
        Gal[p].f_cmax[j][Cb_NUM] = max(0., min(1., f_Cb_max));
		if (isnan(Gal[p].f_cmax[j][Cb_NUM])) {Gal[p].f_cmax[j][Cb_NUM] = 1.;} //f_Cb_max could legitimately be -nan, if the mass of carbon in ring j is 0.0.
		if (f_Cb_max < 0.0 || f_Cb_max > 1.0)
			printf("***** WARNING: model_dust_yields.c: Predicted f_Cb_max has unrealistic value: f_Cb_max = %f. Therefore, it's set to: %f\n", f_Cb_max, Gal[p].f_cmax[j][Cb_NUM]);
		/*if (f_Cb_max < 0.0 || f_Cb_max > 1.0 || isnan(f_Cb_max)) {
		    printf("***** WARNING: f_Cb_max has unrealistic value: f_Cb_max = %f. Therefore, Gal[p].f_cmax[%i][Cb_NUM] set to %f\n", f_Cb_max, j, Gal[p].f_cmax[j][Cb_NUM]);
		    printf("Gal[p].MetalsColdGasRings = %e | num_CO = %e | Gal[p].ColdGasCloudsRings_elements = %e\n", Gal[p].MetalsColdGasRings[j][0]+Gal[p].MetalsColdGasRings[j][1]+Gal[p].MetalsColdGasRings[j][2], num_CO, Gal[p].ColdGasCloudsRings_elements[j][Cb_NUM]);
		}*/
#endif
		Gal[p].f_cmax[j][O_NUM] = max(0., min(1., f_O_max));
		if (isnan(Gal[p].f_cmax[j][O_NUM])) {Gal[p].f_cmax[j][O_NUM] = 1.;} //f_O_max could legitimately be -nan, if the mass of oxygen in ring j is 0.0.
		if (f_O_max < 0.0 || f_O_max > 1.0) {
			printf("***** WARNING: model_dust_yields.c: Predicted f_O_max has unrealistic value: f_O_max = %f. Therefore, it's set to: %f\n", f_O_max, Gal[p].f_cmax[j][O_NUM]);
			/*printf("Ring: %i | %e | %e | %e | %e | %e | %e | %e\n",
					j, Gal[p].DustColdGasCloudsRings_elements[j][H_NUM], Gal[p].DustColdGasCloudsRings_elements[j][He_NUM], Gal[p].DustColdGasCloudsRings_elements[j][N_NUM],
					Gal[p].DustColdGasCloudsRings_elements[j][Cb_NUM], Gal[p].DustColdGasCloudsRings_elements[j][O_NUM], Gal[p].DustColdGasCloudsRings_elements[j][Mg_NUM],
					Gal[p].DustColdGasCloudsRings_elements[j][Fe_NUM]);
			printf("DustColdGasCloudsRings_elements[0][O_NUM] (at start of model_dust_yields.c) = %e\n", DustColdGasCloudsRings_elements_test1);
			printf("DustColdGasCloudsRings_elements[0][O_NUM] (after SN destruction) = %e\n", DustColdGasCloudsRings_elements_test2);
			printf("DustColdGasCloudsRings_elements[0][O_NUM] (2a) = %e\n", DustColdGasCloudsRings_elements_test2a);
			printf("DustColdGasCloudsRings_elements[0][O_NUM] (2b) = %e\n", DustColdGasCloudsRings_elements_test2b);
			printf("DustColdGasCloudsRings_elements[0][O_NUM] (2c) = %e\n", DustColdGasCloudsRings_elements_test2c);
			printf("DustColdGasCloudsRings_elements[0][O_NUM] (before AGB production) = %e\n", DustColdGasCloudsRings_elements_test3);
			printf("DustColdGasCloudsRings_elements[0][O_NUM] (before SN-II production) = %e\n", DustColdGasCloudsRings_elements_test4);
			printf("DustColdGasCloudsRings_elements[0][O_NUM] (now) = %e\n", Gal[p].DustColdGasCloudsRings_elements[0][O_NUM]);
			exit(0);*/
					//num_Silicates, N_O, A_O, num_iron_oxide, N_O_IronOxide, Gal[p].DustColdGasCloudsRings_elements[j][O_NUM], Gal[p].ColdGasCloudsRings_elements[j][O_NUM]);
		}
#else
#ifndef MAINELEMENTS
        Gal[p].f_cmax[Cb_NUM] = max(0., min(1., f_Cb_max));
		if (isnan(Gal[p].f_cmax[Cb_NUM])) {Gal[p].f_cmax[Cb_NUM] = 1.;}
#endif
		Gal[p].f_cmax[O_NUM] = max(0., min(1., f_O_max));
		if (isnan(Gal[p].f_cmax[O_NUM])) {Gal[p].f_cmax[O_NUM] = 1.;}
#endif //H2_AND_RINGS
        

        //*********************************************
        //Calculate total mass in clouds and dust mass in clouds for updating Gal[p].t_acc, and dust mass in diffuse gas for total gas growth rates:
#ifdef H2_AND_RINGS
		for (ee=0; ee<NUM_ELEMENTS; ee++) {
			Clouds_tot[j] += Gal[p].ColdGasCloudsRings_elements[j][ee];
			DustClouds_init[j] += Gal[p].DustColdGasCloudsRings_elements[j][ee];
			DustDiff_init[j] += Gal[p].DustColdGasDiffRings_elements[j][ee];
		}

		if (DustClouds_init[j] > 0.0) {
			//Gal[p].t_acc += Dust_tAcc0*(Clouds_tot/DustClouds_init);
			tacc = Dust_tAcc0*(Clouds_tot[j]/DustClouds_init[j]);
		}
		else tacc = 1e15;

#else //H2_AND_RINGS
		for (ee=0; ee<NUM_ELEMENTS; ee++) {
			Clouds_tot += Gal[p].ColdGasClouds_elements[ee];
			DustClouds_init += Gal[p].DustColdGasClouds_elements[ee];
			DustDiff_init += Gal[p].DustColdGasDiff_elements[ee];
		}

		if (DustClouds_init == 0) {
            tacc = 1e15; //ROB: This should be it's default value anyway... (07-02-22)
        }
        else {
            tacc = Dust_tAcc0*(Clouds_tot/DustClouds_init);
        }
		Gal[p].t_acc = tacc;
#endif //H2_AND_RINGS


		//*****
		//Update dust element fractions:
		float dt_bins = (dt*UnitTime_in_years); //Width of one timestep in years
#ifdef H2_AND_RINGS
        for(ee=0;ee<NUM_ELEMENTS;ee++) {
			Gal[p].f_i[j][ee] = min(1., Gal[p].DustColdGasDiffRings_elements[j][ee]/Gal[p].ColdGasDiffRings_elements[j][ee]); //initial fraction of elements locked-up in dust in the diffuse medium
			Gal[p].f_c[j][ee] = min(1., Gal[p].DustColdGasCloudsRings_elements[j][ee]/Gal[p].ColdGasCloudsRings_elements[j][ee]);
			drop_fnan(j, ee, p);
			update_fractions(dt_bins, tacc, Dust_tExch, Gal[p].H2fractionRings[j], p, ee, j); //updates the fraction of elements locked-up in clouds and diffuse gas, using t_acc, t_exch, f_cmax, initial f_c and f_i, and the H2 fraction.
			drop_fnan(j, ee, p);
			//Updating the amount of dust in clouds after this grain growth:
			Gal[p].DustColdGasCloudsRings_elements[j][ee] = Gal[p].f_c[j][ee]*Gal[p].ColdGasCloudsRings_elements[j][ee];
			//Gal[p].DustColdGasClouds_elements[ee] += Gal[p].f_c[j][ee]*Gal[p].ColdGasCloudsRings_elements[j][ee];
		}
#else //H2_AND_RINGS
        for(ee=0;ee<NUM_ELEMENTS;ee++) {
			Gal[p].f_i[ee] = min(1., Gal[p].DustColdGasDiff_elements[ee]/Gal[p].ColdGasDiff_elements[ee]);
			Gal[p].f_c[ee] = min(1., Gal[p].DustColdGasClouds_elements[ee]/Gal[p].ColdGasClouds_elements[ee]);
			drop_fnan(ee, p);
			update_fractions(dt_bins, tacc, Dust_tExch, Gal[p].H2fraction, p, ee);
			drop_fnan(ee, p);
			//Updating the amount of dust after this mini-step:
			Gal[p].DustColdGasClouds_elements[ee] = Gal[p].f_c[ee]*Gal[p].ColdGasClouds_elements[ee];
		}
#endif //H2_AND_RINGS
            
        //*****
		//Recalculate properties again, after updating dust fractions:
#ifndef MAINELEMENTS
#ifdef H2_AND_RINGS
        num_CO = min((Gal[p].ColdGasCloudsRings_elements[j][Cb_NUM]-Gal[p].DustColdGasCloudsRings_elements[j][Cb_NUM])/A_Cb,
                     min(Gal[p].ColdGasCloudsRings_elements[j][Cb_NUM]/A_Cb*Cmax_CO, (Gal[p].ColdGasCloudsRings_elements[j][O_NUM]-Gal[p].DustColdGasCloudsRings_elements[j][O_NUM])/A_O));
        O_clouds = Gal[p].ColdGasCloudsRings_elements[j][O_NUM] - num_CO*A_O;
        num_Silicates = max(0., min((O_clouds-Gal[p].DustColdGasCloudsRings_elements[j][O_NUM])/(N_O*A_O),
									 min((Gal[p].ColdGasCloudsRings_elements[j][Mg_NUM]-Gal[p].DustColdGasCloudsRings_elements[j][Mg_NUM])/(N_Mg*A_Mg),
										min((Gal[p].ColdGasCloudsRings_elements[j][Si_NUM]-Gal[p].DustColdGasCloudsRings_elements[j][Si_NUM])/(N_Si*A_Si),
										   (Gal[p].ColdGasCloudsRings_elements[j][Fe_NUM]-Gal[p].DustColdGasCloudsRings_elements[j][Fe_NUM])/(N_Fe*A_Fe)))));
#else //H2_AND_RINGS
        num_CO = min((Gal[p].ColdGasClouds_elements[Cb_NUM]-Gal[p].DustColdGasClouds_elements[Cb_NUM])/A_Cb,
					min(Gal[p].ColdGasClouds_elements[Cb_NUM]/A_Cb*Cmax_CO, (Gal[p].ColdGasClouds_elements[O_NUM]-Gal[p].DustColdGasClouds_elements[O_NUM])/A_O));
        O_clouds = Gal[p].ColdGasClouds_elements[O_NUM] - num_CO*A_O;
        num_Silicates = max(0., min((O_clouds-Gal[p].DustColdGasClouds_elements[O_NUM])/(N_O*A_O),
									 min((Gal[p].ColdGasClouds_elements[Mg_NUM]-Gal[p].DustColdGasClouds_elements[Mg_NUM])/(N_Mg*A_Mg),
										min((Gal[p].ColdGasClouds_elements[Si_NUM]-Gal[p].DustColdGasClouds_elements[Si_NUM])/(N_Si*A_Si),
										   (Gal[p].ColdGasClouds_elements[Fe_NUM]-Gal[p].DustColdGasClouds_elements[Fe_NUM])/(N_Fe*A_Fe)))));
#endif //H2_AND_RINGS
#else //MAINELEMENTS
        num_CO = 0.0;
#ifdef H2_AND_RINGS
        O_clouds = Gal[p].ColdGasCloudsRings_elements[j][O_NUM];
#else //H2_AND_RINGS
        O_clouds = Gal[p].ColdGasClouds_elements[O_NUM];
#endif //H2_AND_RINGS
        num_Silicates = 0.0;
#endif //MAINELEMENTS

#ifdef H2_AND_RINGS
        num_iron_oxide = max(0., min((O_clouds-num_Silicates*N_O*A_O-Gal[p].DustColdGasCloudsRings_elements[j][O_NUM])/(N_O_IronOxide*A_O),
        							 (Gal[p].ColdGasCloudsRings_elements[j][Fe_NUM]-num_Silicates*N_Fe*A_Fe-Gal[p].DustColdGasCloudsRings_elements[j][Fe_NUM])/(N_Fe_IronOxide*A_Fe)));
        f_O_max = (num_Silicates*N_O*A_O + num_iron_oxide*N_O_IronOxide*A_O + Gal[p].DustColdGasCloudsRings_elements[j][O_NUM])/Gal[p].ColdGasCloudsRings_elements[j][O_NUM];
#else
        num_iron_oxide = max(0., min((O_clouds-num_Silicates*N_O*A_O-Gal[p].DustColdGasClouds_elements[O_NUM])/(N_O_IronOxide*A_O),
        							 (Gal[p].ColdGasClouds_elements[Fe_NUM]-num_Silicates*N_Fe*A_Fe-Gal[p].DustColdGasClouds_elements[Fe_NUM])/(N_Fe_IronOxide*A_Fe)));
        f_O_max = (num_Silicates*N_O*A_O + num_iron_oxide*N_O_IronOxide*A_O + Gal[p].DustColdGasClouds_elements[O_NUM])/Gal[p].ColdGasClouds_elements[O_NUM];
#endif //H2_AND_RINGS

#ifndef MAINELEMENTS
#ifdef H2_AND_RINGS
        f_Cb_max =  1.0 - num_CO*A_Cb/Gal[p].ColdGasCloudsRings_elements[j][Cb_NUM];
        Gal[p].f_cmax[j][Cb_NUM] = max(0., min(1., f_Cb_max));
        if (isnan(Gal[p].f_cmax[j][Cb_NUM])) {Gal[p].f_cmax[j][Cb_NUM] = 1.;}
#else //H2_AND_RINGS
        f_Cb_max =  1.0 - num_CO*A_Cb/Gal[p].ColdGasClouds_elements[Cb_NUM];
        Gal[p].f_cmax[Cb_NUM] = max(0., min(1., f_Cb_max));
        if (isnan(Gal[p].f_cmax[Cb_NUM])) {Gal[p].f_cmax[Cb_NUM] = 1.;}
#endif //H2_AND_RINGS
#else //MAINELEMENTS
        f_Cb_max =  1.0;
#endif //MAINELEMENTS

#ifdef H2_AND_RINGS
		Gal[p].f_cmax[j][O_NUM] = max(0., min(1., f_O_max));
		if (isnan(Gal[p].f_cmax[j][O_NUM])) {Gal[p].f_cmax[j][O_NUM] = 1.;}
#else //H2_AND_RINGS
		Gal[p].f_cmax[O_NUM] = max(0., min(1., f_O_max));
		if (isnan(Gal[p].f_cmax[O_NUM])) {Gal[p].f_cmax[O_NUM] = 1.;}
#endif //H2_AND_RINGS

#ifdef H2_AND_RINGS
		for (ee=0; ee<NUM_ELEMENTS; ee++)
			DustClouds_tot[j] += Gal[p].DustColdGasCloudsRings_elements[j][ee];

		tacc = Dust_tAcc0*(Clouds_tot[j]/DustClouds_tot[j]);
#else //H2_AND_RINGS
		for (ee=0; ee<NUM_ELEMENTS; ee++)
			DustClouds_tot += Gal[p].DustColdGasClouds_elements[ee];

		t_acc = Dust_tAcc0*(Clouds_tot/DustClouds_tot);
		Gal[p].t_acc = tacc;
#endif //H2_AND_RINGS

#ifdef H2_AND_RINGS
		for(ee=0;ee<NUM_ELEMENTS;ee++) {
			Gal[p].DustColdGasDiffRings_elements[j][ee] = Gal[p].f_i[j][ee] * Gal[p].ColdGasDiffRings_elements[j][ee];
			Gal[p].DustColdGasCloudsRings_elements[j][ee] = Gal[p].f_c[j][ee] * Gal[p].ColdGasCloudsRings_elements[j][ee]; //ROB: Not sure this line is needed, as f_c and ColdGasCloudsRings_elements haven't change since the last update above. (07-03-22)
			//Gal[p].DustColdGasDiff_elements[ee] += Gal[p].f_i[j][ee] * Gal[p].ColdGasDiffRings_elements[j][ee];
			//Gal[p].DustColdGasClouds_elements[ee] += Gal[p].f_c[j][ee] * Gal[p].ColdGasCloudsRings_elements[j][ee];
			//DustDiff_new[j] += Gal[p].ColdGasCloudsRings_elements[j][ee];
			DustDiff_new[j] += Gal[p].DustColdGasDiffRings_elements[j][ee];
			DustClouds_new[j] += Gal[p].DustColdGasCloudsRings_elements[j][ee];
		}
		DustDiff_Growth = DustDiff_new[j] - DustDiff_init[j];
		DustClouds_Growth = DustClouds_new[j] - DustClouds_init[j];
#else //H2_AND_RINGS
        for(ee=0;ee<NUM_ELEMENTS;ee++) {
			Gal[p].DustColdGasDiff_elements[ee] = Gal[p].f_i[ee] * Gal[p].ColdGasDiff_elements[ee];
			Gal[p].DustColdGasClouds_elements[ee] = Gal[p].f_c[ee] * Gal[p].ColdGasClouds_elements[ee];
			//DustDiff_new += Gal[p].ColdGasClouds_elements[ee];
			DustDiff_new += Gal[p].DustColdGasDiff_elements[ee];
			DustClouds_new += Gal[p].DustColdGasClouds_elements[ee];
		}
        DustDiff_Growth = DustDiff_new - DustDiff_init;
        DustClouds_Growth = DustClouds_new - DustClouds_init;
#endif //H2_AND_RINGS

#ifdef FULL_DUST_RATES
		Gal[p].DustColdGasRates[3] += (DustDiff_Growth + DustClouds_Growth)/(deltaT * UnitTime_in_years);
		//if (Gal[p].SnapNum > 50) printf("SnapNum = %i | DustDiff_Growth = %e | DustClouds_Growth = %e | SnapshotWidth = %e [Myr] | increment = %e | DustColdGasRates[3] = %e\n", Gal[p].SnapNum, DustDiff_Growth, DustClouds_Growth, (deltaT * UnitTime_in_years)/1.e6, (DustDiff_Growth + DustClouds_Growth)/(deltaT * UnitTime_in_years), Gal[p].DustColdGasRates[3]);
#endif
        
		mass_checks(p,"End of dust from growth: model_dust_yields.c",__LINE__);
        
} //if (((Gal[p].MetalsColdGas[0]+Gal[p].MetalsColdGas[1]+Gal[p].MetalsColdGas[2])>0.0) && (Gal[p].ColdGas > 0.0))
#ifdef H2_AND_RINGS
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

/*
//Calculate final global t_des from all rings:
double tot_SNRate=0., tot_ColdGasRings=0.;
double frac_Cb_final, m_cleared_final;
for (int jjj=0;jjj<RNUM;jjj++) {
	tot_SNRate += DiskSNIIRate_current_ts[jjj] + DiskSNIaRate_current_ts[jjj];
	tot_ColdGasRings += Gal[p].ColdGasRings[jjj];
}
//Gal[p].t_des = (tot_ColdGasRings*(1.0e10/Hubble_h))/(M_CLEARED * F_SN * (tot_SNRate/(UnitTime_in_s / SEC_PER_YEAR))); //Note, this will be outputted as the (global) destruction timescale of the last timestep of each snapshot, rather than the average across the whole snapshot. (09-02-22)
frac_Cb_final = (Gal[p].DustColdGasCloudsRings_elements[j][Cb_NUM]+Gal[p].DustColdGasDiffRings_elements[j][Cb_NUM])/(Gal[p].DustColdGasCloudsRings_elements[j][Cb_NUM]+Gal[p].DustColdGasCloudsRings_elements[j][Si_NUM]+Gal[p].DustColdGasDiffRings_elements[j][Cb_NUM]+Gal[p].DustColdGasDiffRings_elements[j][Si_NUM]); //Note, this will be outputted as the (global) destruction timescale of the last timestep of each snapshot, rather than the average across the whole snapshot. (09-02-22)
m_cleared_final = calc_cleared_mass(frac_Cb_final);
Gal[p].t_des = (tot_ColdGasRings*(1.0e10/Hubble_h))/(m_cleared_final * F_SN * (tot_SNRate/(UnitTime_in_s / SEC_PER_YEAR))); //Note, this will be outputted as the (global) destruction timescale of the last timestep of each snapshot, rather than the average across the whole snapshot. (09-02-22)
*/

//Calculate final global t_acc from all rings:
double Clouds_tot_allRings=0., DustClouds_tot_allRings=0.;
for(j=0;j<RNUM;j++) {
	Clouds_tot_allRings += Clouds_tot[j];
	DustClouds_tot_allRings += DustClouds_tot[j];
}
Gal[p].t_acc = Dust_tAcc0*(Clouds_tot_allRings/DustClouds_tot_allRings);
#endif //H2_AND_RINGS
#endif //DUST_GROWTH

/*//if (checka == 1.) {
		printf("***** End of model_dust_yields.c:\n");
		printf("Snapnum %i\n", Gal[p].SnapNum);
		for (ee=0; ee<NUM_ELEMENTS; ee++)
			printf("Gal[p].DustColdGasClouds_elements[%i] = %f\n", ee, Gal[p].DustColdGasClouds_elements[ee]);
		for (j=0;j<RNUM;j++)
			printf("Gal[p].DustColdGasCloudsRings_elements[%i][O] = %f\n", j, Gal[p].DustColdGasCloudsRings_elements[j][O_NUM]);
	//}*/
    
}//update_dust_mass()

