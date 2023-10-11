/*
 * dust_integrals.c
 *
 * A dust version of Rob's yields_integrals/model_yields
 *
 * Pre-calculates the normalised ejecta rates at every timestep, assuming 1 Msun populations.
 * Multiply by SFR from SFH bins (and interpolate between default metallicities) to obtain
 * true dust ejecta rates (done in model_dust.c).
 *
 *
 * Created: Oct2016
 * Author: ScottClay
 *
 * Nov 2017 - Cleaned up code and added comments
 *
 * Edited from: Nov2021
 * Author: robyates
 *
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "allvars.h"
#include "proto.h"

#ifdef DETAILED_DUST

void init_integrated_dust_yields()
{
	int ii, jj, kk, ll;

	for(ii=0;ii<STEPS*(LastDarkMatterSnapShot+1);ii++) {
		for(jj=0;jj<SFH_NBIN;jj++) {
			for(kk=0;kk<ARRAY_Z_NUM;kk++) {
				for(ll=0;ll<AGB_DUST_TYPE_NUM;ll++) {
					NormAGBDustYieldRate[ii][jj][kk][ll]=0.0;									
}}}}}

void integrate_dust_yields()
{
	double previoustime, newtime, deltaT;
	int snap, step,i,mb;
	double timet;

	int Mi_lower, Mi_upper,Zi_Dust;
	int Mi_lower_AGB, Mi_upper_AGB;
	int Mi_lower_Dust, Mi_upper_Dust;
	int width_in_timesteps, mbmax; //Number of current timesteps that fit in any given SFH bin, and the number of mini bins considered for any given SFH bin (max. = 30, for memory considerations)
	double dt, t_lower, t_upper;	
	
	for(snap=0;snap<(LastDarkMatterSnapShot+1);snap++) //LOOP OVER SNAPSHOTS
	{
	    previoustime = NumToTime(snap); //Time to z=0 from start of current snapshot [in code units]
	    newtime = NumToTime(snap+1); //Time to z=0 from end of current snapshot [in code units]
	    deltaT = previoustime - newtime; //snapshot width

	    for(step=0;step<STEPS;step++) //LOOP OVER TIMESTEPS
	    {
	    	dt = deltaT/STEPS;  //Time-width of a timestep in current snapshot [in code units]
	    	timet = previoustime - (step + 0.5) * dt; //Time from middle of the current timestep to z=0 [in code units]
	        for (i=0;i<=SFH_ibin[snap][step];i++) //LOOP OVER SFH BINS
	        {
	   	  	//New method: sub-dividing SFH bins into SFH_NMINIBIN number of 'mini bins': Can later choose inside code which mini bin the characteristic SF time is in:
	        	width_in_timesteps = SFH_dt[snap][step][i]/dt; //Width of SFH bin in number of current timesteps [in code units] //NB: Typecasting a float to an integer here (width_in_timesteps is and integer).
	        	if(width_in_timesteps < 1) width_in_timesteps = 1;
	        	mbmax = width_in_timesteps;
	        	for (mb=1;mb<=mbmax;mb++) //LOOP OVER MINI BINS (New method)
	        	{
	        		//From lower/upper edge of mini-bin to middle of current timestep:
	        		t_lower = (SFH_t[snap][step][i] + (SFH_dt[snap][step][i]) - (mb*(SFH_dt[snap][step][i]/mbmax)) - timet) * UnitTime_in_years/Hubble_h; //IN YEARS //Time from low-z (lower) edge of SFH mini-bin j to middle of current timestep
	        		t_upper = (SFH_t[snap][step][i] + (SFH_dt[snap][step][i]) - (mb*(SFH_dt[snap][step][i]/mbmax)) + (SFH_dt[snap][step][i]/mbmax) - timet) * UnitTime_in_years/Hubble_h; //IN YEARS //Time from high-z (upper) edge of SFH mini-bin j to middle of current timestep

	        	  int Zi;
	        	  for (Zi=0;Zi<ARRAY_Z_NUM;Zi++) //LOOP OVER POSSIBLE INITIAL METALLICITIES
	        	  {
	        	  	
	        	  	//Mi_upper goes to get the bin number, so an integer, from the lifetime array
	        		//Mi_lower = bin number for lifetime array for lowest mass star possible to die
	        		//at this time
	        		Mi_lower = find_initial_mass_dust(t_upper, Zi); //Mass bin (lifetime arrays) corresp. to lowest mass of star to 'die' in current timestep, from SFH bin i.
	        		Mi_upper = find_initial_mass_dust(t_lower, Zi); //Mass bin (lifetime arrays) corresp. to lowest mass of star to 'die' in current timestep, from SFH bin i.
					
					//*****************************************
					//AGB Dust Yields
					//*****************************************
					
					//finds the metallicity in dust table that corresponds to metallicity bin in the lifetimes table 
					Zi_Dust = find_initial_metallicity_comp_dust(Zi, i, 4); 
					
					//Mi_lower_AGB = integer mass bin from lifetime arrays
					//max/min functions below shift Mi_lower_AGB to correspond to lifetime
					//bin closest to limits if we have a limit problem
					Mi_lower_AGB = max_Mi_lower_dust(Mi_lower,4); //Mass bin (lifetime arrays) corresp. to lowest mass of star to 'die' in current timestep, from SFH bin i.
					Mi_upper_AGB = min_Mi_upper_dust(Mi_upper,4); //Mass bin (lifetime arrays) corresp. to lowest mass of star to 'die' in current timestep, from SFH bin i.

					Mi_lower_Dust = find_agb_mass_bin_dust(lifetimeMasses[Mi_lower_AGB]);
					Mi_upper_Dust = find_agb_mass_bin_dust(lifetimeMasses[Mi_upper_AGB]);

					//Actual integration
					int j;
					for (j=Mi_lower_Dust;j<Mi_upper_Dust;j++)
					{       
						int k;
						for (k=0;k<AGB_DUST_TYPE_NUM;k++) { 
							NormAGBDustYieldRate[(STEPS*snap)+step][i][Zi][k] += (AGBDustMasses[j+1]-AGBDustMasses[j]) * ((AGBDustCreated[Zi_Dust][j][k] + AGBDustCreated[Zi_Dust][j+1][k])/2.0);
							//if (NormAGBDustYieldRate[(STEPS*snap)+step][i][Zi][k] < 0.0) printf("step = %i | SFH bin = %i, Z = %i | k = %i | NormAGBDustYieldRate = %e\n",(STEPS*snap)+step,i,Zi,k,NormAGBDustYieldRate[(STEPS*snap)+step][i][Zi][k]);
							}
					}	//for (j=Mi_lower_AGB;j<=Mi_upper_AGB;j++)
	        	} //for (Zi=0;Zi<ARRAY_Z_NUM;Zi++)
	          } //for (j=0;j<=width_in_timesteps;j++) //MINI_BINS	
	        } //for (i=0;i<=SFH_ibin_structure[(SFH_NBIN*snap)+step];i++)
	    } //for(step=0;step<STEPS;step++)
	} //for(snap=0;snap<(LastDarkMatterSnapShot+1);snap++)
	printf("Dust Yield integrals calculated.\n");
}

int find_initial_metallicity_comp_dust(int Zi, int sfh_bin, int table_type)
{
	int i, Zi_bin;
	double Z_in;

	Zi_bin = -1; //This is the bin we want. The closest in dust tables to lifetime metallicity
	i = 0;
//#ifdef BINARYC
//	Z_in = bcMetallicities[Zi];
//#else
	Z_in = lifetimeMetallicities[Zi]; //Actual lifetime metallicity 
//#endif //BINARYC

	switch (table_type)
	{
		case 4: //AGB dust metallicity table
			while (Zi_bin == -1)
			{
				if (AGBDustMetallicities[i] < Z_in)
				{
					i++;
					if (i == AGB_DUST_METAL_NUM) Zi_bin = i-1;
				}
				else Zi_bin = i;
			}
			break;
	}
	return Zi_bin;
}

int find_initial_mass_dust(double lifetime, int Zi_bin)
{
    if (lifetime == 0.0) return LIFETIME_MASS_NUM-1; //If the bin 'touches now', then return max mass (ie: star of shortest lifetime) ie: bin for 120Msun
    else if (lifetime > lifetimes[Zi_bin][0]) return 0; //If true lifetime is longer than max lifetime in table (shouldn't be), then return element number 0
    else {
		int Mi_bin;

		Mi_bin = -1;
		int i;
		i = 0;
		while (Mi_bin == -1) {
			if (lifetimes[Zi_bin][i] > lifetime) {
				i++;
				if (i == LIFETIME_MASS_NUM) Mi_bin = i; //If lifetime is shorter than min lifetime from table, then just return max mass (120 Msun)
			}
			else Mi_bin = i;
		}
		return Mi_bin-1; //This returns element number i for lifetimeMasses[Zi][i] array BELOW the true initial mass corresponding to t_lower or t_upper.
    }
}

int max_Mi_lower_dust(int Mi_lower, int channel_type)
{
	switch (channel_type)
		{
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

int min_Mi_upper_dust(int Mi_upper, int channel_type)
{
	switch (channel_type)
		{
			case 4: //AGB mass limits
					if (lifetimeMasses[Mi_upper] < AGB_MAX_MASS) return Mi_upper;
					else
					{
						int i;
						i = LIFETIME_MASS_NUM-1;
						do { i--; }
						while (lifetimeMasses[i] > AGB_MAX_MASS);
						return i;
					}
					break;
			default: printf("Wrong ejection mode chosen in min_Mi_upper: Use either 2 (SNe-II), 3 (SNe-Ia) or 4 (AGB winds)"); exit(1);
		}
}


int find_agb_mass_bin_dust(double masslimit)
{
	if (masslimit == AGB_MAX_MASS) return AGB_DUST_MASS_NUM-1; //If AGB_MAX_MASS == mass limit = 7Msol, this is max for agb dust tables too. So this works. If you change AGB_Max_Mass - this will need changing
	else
	{
	int Mi_bin;

	Mi_bin = -1;
	int i;
	i = 0;
	while (Mi_bin == -1)
	{
		if (AGBDustMasses[i] < masslimit)
		{
			i++;
			if (i == AGB_DUST_MASS_NUM) Mi_bin = i-1; //If mass is greater than max mass for AGB winds (shouldn't be), then just return max mass (5.0 Msun)
		}
		else Mi_bin = i;
	}
	return Mi_bin;
	}
}

/*void find_actual_ejecta_limits_dust(int channel_type, double Mi_lower_actual, double Mi_upper_actual, int Mi_lower, int Mi_upper, int Zi,
double* Yields_lower_actual, double* Yields_upper_actual)
{
	switch (channel_type)
	{
	case 4: //AGB
		if (Mi_lower == 0)
		{
		    int k;
		    for (k=0;k<AGB_DUST_TYPE_NUM;k++)
		    {
		    	Yields_lower_actual[k] = AGBDustCreated[Zi][0][k];
		    }
		}
		else
		{
		    int k;
		    for (k=0;k<AGB_DUST_TYPE_NUM;k++)
		    {
				Yields_lower_actual[k] = AGBDustCreated[Zi][Mi_lower][k] + ((AGBDustCreated[Zi][Mi_lower+1][k]-AGBDustCreated[Zi][Mi_lower][k]) * ((Mi_lower_actual-AGBDustMasses[Mi_lower])/(AGBDustMasses[Mi_lower+1]-AGBDustMasses[Mi_lower])));
		    }
		}

		if (Mi_upper == AGB_MASS_NUM-1)
		//if (Mi_upper == AGB_DUST_MASS_NUM-1)
		{
		    int k;
		    for (k=0;k<AGB_DUST_TYPE_NUM;k++)
		    {
		    	Yields_upper_actual[k] = AGBDustCreated[Zi][AGB_DUST_MASS_NUM-1][k];
		    }
		}
		else
		{
		    int k;
		    for (k=0;k<AGB_DUST_TYPE_NUM;k++)
		    {
		    	Yields_upper_actual[k] = AGBDustCreated[Zi][Mi_upper][k] + ((AGBDustCreated[Zi][Mi_upper+1][k]-AGBDustCreated[Zi][Mi_upper][k]) * ((Mi_upper_actual-AGBDustMasses[Mi_upper])/(AGBDustMasses[Mi_upper+1]-AGBDustMasses[Mi_upper])));
		     }
		}
		break;
	}
}*/

#endif //DETAILED_EMRICHMENT

