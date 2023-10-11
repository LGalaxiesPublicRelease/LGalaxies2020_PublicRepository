/*
 * model_yields.c
 *
 *  Created on: 18.11.2011
 *      Author: robyates
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "allvars.h"
#include "proto.h"

//INDIVIDUAL ELEMENTS
//All: [H][He][Cb][N][O][Ne][Mg][Si][S][Ca][Fe] or //Only [H][He][O][Mg][Fe]


void update_yields_and_return_mass(int p, int centralgal, double dt, int nstep)
{
	int Zi, igal, ii, mm;
	double timestep_width; //Width of current timestep in CODE UNITS
	int TimeBin; //Bin in Yield arrays corresponding to current timestep
	double Zi_disp;
	double sfh_time;
	//double time_to_ts; //Time from high-z (upper) edge of SFH bin to middle of current timestep (used for massive SNII to hot) [in Myrs]
	//double tcut; //Maximum lifetime of stars that have their ejected put straight into the HotGas [in Myrs]
#ifdef METALRICHWIND
#ifdef GASDENSITYFWIND
	double ColdGasSurfaceDensity;
#endif
#endif
	/*
	//ROB: Now made global variables in h_variables.h:
	double fwind_SNII; //Required for SN-II metal-rich wind implementation
	double fwind_SNIa; //Required for SN-Ia metal-rich wind implementation
	double fwind_AGB; //Required for AGB metal-rich wind implementation*/
	double DiskSFRxStep, DiskSFRxStep_Phys, BulgeSFRxStep, BulgeSFRxStep_Phys, ICMSFRxStep, ICMSFRxStep_Phys;
	double DiskMetallicity, BulgeMetallicity, ICMMetallicity;
	//double TotalMassReturnedToColdDiskGas, TotalMassReturnedToHotGas, TotalMetalsReturnedToHotGas; //ROB: Now declared in h_varibale.h, so they can be called from main.c (08-11-21).
	double SNIIEjectaMass, SNIaEjectaMass, AGBEjectaMass;
	double SNIIRate, SNIaRate, AGBRate;
	double SNIIUnProcessedMetals, SNIaUnProcessedMetals, AGBUnProcessedMetals, SNIIAllMetals, SNIaAllMetals, AGBAllMetals;
#ifdef INDIVIDUAL_ELEMENTS
	int ee;
	double SNIIAllElements[NUM_ELEMENTS], SNIaAllElements[NUM_ELEMENTS], AGBAllElements[NUM_ELEMENTS];
	double SNIIUnProcessedElements[NUM_ELEMENTS], SNIaUnProcessedElements[NUM_ELEMENTS], AGBUnProcessedElements[NUM_ELEMENTS];
	double DiskMetallicityElement_Phys[NUM_ELEMENTS];
	double BulgeMetallicityElement_Phys[NUM_ELEMENTS];
	double ICMMetallicityElement_Phys[NUM_ELEMENTS];
#ifdef DETAILED_DUST
#ifdef H2_AND_RINGS
	double coldgas_init[RNUM][NUM_ELEMENTS];
#else
	double coldgas_init[NUM_ELEMENTS];
#endif //H2_AND_RINGS
	double coldgas_add;
#endif //DETAILED_DUST
#endif //INDIVIDUAL_ELEMENTS
	int n; //Iterator used for loop over NOUT when updating MassWeightedAge
	double AgeCorrectionDisk[NOUT];
	double AgeCorrectionBulge[NOUT];
#ifdef H2_AND_RINGS
	int jj;
	double fractionRings[RNUM];
	double fractionRingsBulge[RNUM];
	//double TotalMassReturnedToColdDiskGasr[RNUM]; //ROB: Now declared in h_varibale.h, so they can be called from main.c (08-11-21).
#endif

	//Set required properties to zero:
	TotalMassReturnedToHotGas=0.0;
	TotalMetalsReturnedToHotGas=0.0;
	TotalMassReturnedToColdDiskGas=0.0;
#ifdef DETAILED_DUST
#ifdef H2_AND_RINGS
    for(jj=0;jj<RNUM;jj++) {
    	TotalMassReturnedToColdDiskGasr[jj]=0.0;
    	DiskSNIIRate_current_ts[jj]=0.0; //Instantaneous SNII rate in the disc (stored for use in model_dust_yields.c)
		DiskSNIaRate_current_ts[jj]=0.0; //Instantaneous SNIa rate in the disc (stored for use in model_dust_yields.c)
		DiskAGBRate_current_ts[jj]=0.0; //Instantaneous AGB rate in the disc (stored for use in model_dust_yields.c)
		for(ee=0;ee<NUM_ELEMENTS;ee++) {
			DiskSNIIAllElements_ts[jj][ee]=0.0; //Total element masses ejected by all SNe-II in this timestep (stored for use in model_dust_yields.c)
			DiskSNIaAllElements_ts[jj][ee]=0.0; //Total element masses ejected by all SNe-Ia in this timestep (stored for use in model_dust_yields.c)
			DiskAGBAllElements_ts[jj][ee]=0.0; //Total element masses ejected by all AGBs in this timestep (stored for use in model_dust_yields.c)
		}
    }
#else //H2_AND_RINGS
    DiskSNIIRate_current_ts=0.0; //Instantaneous SNII rate in the disc (stored for use in model_dust_yields.c)
    DiskSNIaRate_current_ts=0.0; //Instantaneous SNIa rate in the disc (stored for use in model_dust_yields.c)
    DiskAGBRate_current_ts=0.0; //Instantaneous AGB rate in the disc (stored for use in model_dust_yields.c)
    for(ee=0;ee<NUM_ELEMENTS;ee++) {
    	DiskSNIIAllElements_ts[ee]=0.0; //Total element masses ejected by all SNe-II in this timestep (stored for use in model_dust_yields.c)
    	DiskSNIaAllElements_ts[ee]=0.0; //Total element masses ejected by all SNe-Ia in this timestep (stored for use in model_dust_yields.c)
    	DiskAGBAllElements_ts[ee]=0.0; //Total element masses ejected by all AGBs in this timestep (stored for use in model_dust_yields.c)
    }
#endif //H2_AND_RINGS
#else //DETAILED_DUST
#ifdef H2_AND_RINGS
    for(jj=0;jj<RNUM;jj++)
    	TotalMassReturnedToColdDiskGasr[jj]=0.0;
#endif //H2_AND_RINGS
#endif //DETAILED_DUST

    for(n=0;n<NOUT;n++) {
    	AgeCorrectionDisk[n] = 0.0;
    	AgeCorrectionBulge[n] = 0.0;
    }

	timestep_width = dt; //Width of current timestep in CODE UNITS (units cancel out when dividing by SFH bin width, sfh_dt) (12-04-12)
	TimeBin = (STEPS*(Halo[Gal[p].HaloNr].SnapNum-1.0))+nstep; //TimeBin = (STEPS*Gal[p].SnapNum)+nstep; //Bin in Yield tables corresponding to current timestep //TEST!: BRUNO: Snapnum would be +1 too low for a 'jumping' galaxy (14-11-13)
	//timet = NumToTime((Halo[Gal[p].HaloNr].SnapNum-1.0)) - (nstep + 0.5) * dt; //Time from middle of the current timestep to z=0 (used here for MassWeightAge corrections)
	//NB: NumToTime(Gal[p].SnapNum) is the time to z=0 from start of current snapshot
	//    nstep is the number of the current timestep (0-19)
	//    dt is the width of one timestep within current snapshot


	//for stars dying that enrich the Hot gas directly
    if(Gal[p].Type==2)
      igal=Gal[p].CentralGal;
    else
      igal=p;
      
#ifdef DETAILED_DUST
    for (ee=0;ee<NUM_ELEMENTS;ee++) {
#ifdef H2_AND_RINGS
    	for(jj=0;jj<RNUM;jj++)
    		coldgas_init[jj][ee] = Gal[p].ColdGasRings_elements[jj][ee];
#else
    	coldgas_init[ee] = Gal[p].ColdGas_elements[ee];
#endif //H2_AND_RINGS
    }
#endif  

    //**********
    for (ii=0;ii<=Gal[p].sfh_ibin;ii++) //LOOP OVER SFH BINS
    {
    	sfh_time=Gal[p].sfh_t[ii]+(0.5*Gal[p].sfh_dt[ii]);
    	//time_to_ts = ((sfh_time+(0.5*Gal[p].sfh_dt[i])) - timet)*(UnitTime_in_years/Hubble_h)/1.0e6; //Time from high-z (upper) edge of SFH bin to middle of current timestep [in Myrs]
    	//tcut = 2.0*((Gal[p].Rvir/Gal[p].Vvir)/0.0001); //Maximum lifetime of stars that have their ejected put straight into the HotGas [in Myrs]

#ifdef DETAILED_DUST
// These variables store the amount of metals produced by SNe in a timestep for use in model_dustyields.c
//ROB: Looks more like they store the mass of an element ejected into the ColdGas in the current timestep by disc SNe with progenitors born in each of the SFH bins. (29-11-21)
#ifdef H2_AND_RINGS
    	for(jj=0;jj<RNUM;jj++) {
			SNII_prevstep_Cold_Si[jj][ii] = 0.0;
			SNII_prevstep_Cold_Fe[jj][ii] = 0.0;
			SNII_prevstep_Cold_Cb[jj][ii] = 0.0;
			SNIa_prevstep_Cold_Fe[jj][ii] = 0.0;
    	}
#else
    	SNII_prevstep_Cold_Si[ii] = 0.0;
		SNII_prevstep_Cold_Fe[ii] = 0.0;
		SNII_prevstep_Cold_Cb[ii] = 0.0;
		SNIa_prevstep_Cold_Fe[ii] = 0.0;
#endif //H2_AND_RINGS
#endif //DETAILED_DUST

    	mass_checks(p,"model_yields.c",__LINE__);


    //******************************************************
    //ENRICHMENT FROM DISK STARS (INTO COLD GAS & HOT GAS):
    //******************************************************
#ifdef H2_AND_RINGS
    for(jj=0;jj<RNUM;jj++) 	
    	if (Gal[p].DiskMassRings[jj] > 0.0)//ROB: Discs can be destroyed (i.e. converted in to bulges). So only calculate enrichment from stars born in the disc if there is still a disc
    		if (Gal[p].sfh_DiskMassRings[jj][ii] > 0.0)    	
    		{		
#else
    if (Gal[p].DiskMass > 0.0)//ROB: Discs can be destroyed (i.e. converted in to bulges). So only calculate enrichment from stars born in the disc if there is still a disc at the current timestep. This if statement has scope until the end of the following if statement.
    	if (Gal[p].sfh_DiskMass[ii] > 0.0)
    	{
#endif
    

#ifdef METALRICHWIND
#ifdef GASDENSITYFWIND
#ifndef H2_AND_RINGS
    	ColdGasSurfaceDensity = (Gal[p].ColdGas*(1.0e10)*Hubble_h)/(3.14159265*Gal[p].ColdGasRadius*Gal[p].ColdGasRadius*(1.0e6*1.0e6)); //in Msun/pc^2
    	fwind_SNII = min(1.0, 1.0/(ColdGasSurfaceDensity/NORMGASDENSITY)); //Fraction of SN-II ejecta put directly into HotGas. When ColdGasSurfaceDensity <= NORMGASDENSITY (i.e. 10.0 Msun/pc), then fwind_SNII = 1.0.
    	//fwind_SNIa = min(1.0, 1.0/(ColdGasSurfaceDensity/NORMGASDENSITY)); //Fraction of SN-Ia ejecta put directly into HotGas. When ColdGasSurfaceDensity <= NORMGASDENSITY (i.e. 10.0 Msun/pc), then fwind_SNIa = 1.0.
    	//fwind_AGB = min(1.0, 1.0/(ColdGasSurfaceDensity/NORMGASDENSITY)); //Fraction of AGB ejecta put directly into HotGas. When ColdGasSurfaceDensity <= NORMGASDENSITY (i.e. 10.0 Msun/pc), then fwind_AGB = 1.0.
    	fwind_SNIa = FracZSNIatoHot;
    	fwind_AGB = FracZAGBtoHot;
#else
    	if(jj==0)
    		//ColdGasSurfaceDensity = (Gal[p].ColdGasRings[jj]*(1.0e10))/(4.0*3.14159265*RingRadius[jj]*RingRadius[jj]*(1.0e6*1.0e6)); //CHECK h FACTORS!! in Msun/pc^2
    		ColdGasSurfaceDensity = (Gal[p].ColdGasRings[jj]*(1.0e10)*Hubble_h)/(3.14159265*RingRadius[jj]*RingRadius[jj]*(1.0e6*1.0e6)); //in Msun/pc^2
    	else
    		ColdGasSurfaceDensity = (Gal[p].ColdGasRings[jj]*(1.0e10)*Hubble_h)/(3.14159265*((RingRadius[jj]*RingRadius[jj])-(RingRadius[jj-1]*RingRadius[jj-1]))*(1.0e6*1.0e6)); //in Msun/pc^2

    	fwind_SNII = min(1.0, NORMGASDENSITY / log10(ColdGasSurfaceDensity)); //Fraction of SN-II ejecta put directly into HotGas.
    	fwind_SNIa = min(1.0, NORMGASDENSITY / log10(ColdGasSurfaceDensity)); //Fraction of SN-Ia ejecta put directly into HotGas.
    	fwind_AGB = 0.0; //Fraction of AGB ejecta put directly into HotGas.

#endif	//H2_AND_RINGS
//#endif	//GASDENSITYFWIND
//#ifndef GASDENSITYFWIND
#else //GASDENSITYWIND
#ifdef MVIRDEPFWIND
    	fwind_SNII = min(1.0, (FracZSNIItoHot*10.5) / log10(1.e10*Gal[p].Mvir)); //min(1.0, (FracZSNIItoHot*10.) * (Hubble_h/Gal[p].Mvir)); //Sets direct ejection into HotGas by SNe-II to be dependent on gravitating mass (i.e. Mvir). Coefficient chose to be 9.0 to roughly match a value of fwind_SNII = 0.9 for low-mass galaxies [log(M*/Msun)~9.0 --> log(Mvir/Msun)~11.0]. 02-03-20
    	fwind_SNIa = min(1.0, (FracZSNIatoHot*10.5) / log10(1.e10*Gal[p].Mvir));
    	fwind_AGB = min(1.0, (FracZAGBtoHot*10.5) / log10(1.e10*Gal[p].Mvir));
#else
    	fwind_SNII = FracZSNIItoHot; //FracZtoHot;
    	fwind_SNIa = FracZSNIatoHot;
    	fwind_AGB = FracZAGBtoHot;
        //fwind_SNII = min(1.0, 9. * (Hubble_h/Gal[p].Mvir));
        //printf("fwind_SNII = %f, %f, %.2e\n", fwind_SNII, Hubble_h, (Gal[p].Mvir*1.e10)/Hubble_h);
#endif
#endif //GASDENSITYFWIND

#else //METALRICHWIND
     	fwind_SNII = 0.0; //For all stellar ejecta (from disk) to ColdGas
     	fwind_SNIa = 0.0;
     	fwind_AGB = 0.0;
#endif //METALRICHWIND

    	//pre-calculations to speed up the code
	    //timestep_width and dt units cancel out
#ifdef H2_AND_RINGS
    	DiskSFRxStep = timestep_width * Gal[p].sfh_DiskMassRings[jj][ii]/Gal[p].sfh_dt[ii];
    	DiskSFRxStep_Phys = DiskSFRxStep * (1.0e10/Hubble_h); //ROB: I think this is actually just the *mass* of stars formed in the disc in SFH bin ii [in Msun].

    	DiskMetallicity = 0.;
    	for(mm=0;mm<NUM_METAL_CHANNELS;mm++)
    		DiskMetallicity += Gal[p].sfh_MetalsDiskMassRings[jj][ii][mm];
    	DiskMetallicity /= Gal[p].sfh_DiskMassRings[jj][ii];

    	//printf("DiskMetallicity = %f\n", DiskMetallicity);
    	//***** FOR TESTING FIXED-METALLICITY YIELDS: ********
    	//DiskMetallicity = 0.0134; //Solar metallicity (Asplund+09)
    	//****************************************************

#ifdef INDIVIDUAL_ELEMENTS
    	for (ee=0;ee<NUM_ELEMENTS;ee++)
    		DiskMetallicityElement_Phys[ee] = Gal[p].sfh_DiskMass_elementsRings[jj][ii][ee] / (Gal[p].sfh_DiskMassRings[jj][ii]*1.0e10/Hubble_h);
#endif

#else //H2_AND_RINGS
    	DiskSFRxStep = timestep_width * Gal[p].sfh_DiskMass[ii]/Gal[p].sfh_dt[ii];
    	DiskSFRxStep_Phys = DiskSFRxStep * (1.0e10/Hubble_h);

    	DiskMetallicity = 0.;
    	for(mm=0;mm<NUM_METAL_CHANNELS;mm++)
    		DiskMetallicity += Gal[p].sfh_MetalsDiskMass[ii][mm];
    	DiskMetallicity /= Gal[p].sfh_DiskMass[ii];


#ifdef INDIVIDUAL_ELEMENTS
    	for (ee=0;ee<NUM_ELEMENTS;ee++)
    		DiskMetallicityElement_Phys[ee] = Gal[p].sfh_DiskMass_elements[ii][ee] / (Gal[p].sfh_DiskMass[ii]*1.0e10/Hubble_h);
#endif

#endif //H2_AND_RINGS

#ifdef BINARYC
    	Zi = find_initial_metallicity(DiskMetallicity, 1, 5);
#else //BINARYC
    	Zi = find_initial_metallicity(DiskMetallicity, 1, 1);
#endif //BINARYC
    	//Interpolate the disk metallicity on the lifetimeMetallicities tables:
    	Zi_disp = (DiskMetallicity - lifetimeMetallicities[Zi])/(lifetimeMetallicities[Zi+1] - lifetimeMetallicities[Zi]);
    	if (Zi_disp < 0.0) Zi_disp = 0.0; //Don't want to extrapolate yields down below lifetimeMetallicities[0]=0.0004. Instead, assume constant yield below this metallicity.

#ifdef DETAILED_DUST //For use in update_dust_mass() later on
#ifdef H2_AND_RINGS
    	Zi_saved[jj][ii] = Zi;
		Zi_disp_saved[jj][ii] = Zi_disp;
#else
		Zi_saved[ii] = Zi;
		Zi_disp_saved[ii] = Zi_disp;
#endif //H2_AND_RINGS
#endif

#ifndef INDIVIDUAL_ELEMENTS
    	compute_actual_eject_rates(TimeBin, ii, Zi, Zi_disp, Gal[p].sfh_ibin, DiskSFRxStep, DiskSFRxStep_Phys, DiskMetallicity, timestep_width,
				   &SNIIEjectaMass, &SNIIAllMetals, &SNIIUnProcessedMetals,
				   &SNIaEjectaMass, &SNIaAllMetals, &SNIaUnProcessedMetals,
				   &AGBEjectaMass, &AGBAllMetals, &AGBUnProcessedMetals,
				   &SNIIRate, &SNIaRate, &AGBRate);
#else
    	compute_actual_eject_rates(TimeBin, ii, Zi, Zi_disp, Gal[p].sfh_ibin, DiskSFRxStep, DiskSFRxStep_Phys, DiskMetallicity, DiskMetallicityElement_Phys, timestep_width,
				   &SNIIEjectaMass, &SNIIAllMetals, &SNIIUnProcessedMetals,
				   &SNIaEjectaMass, &SNIaAllMetals, &SNIaUnProcessedMetals,
				   &AGBEjectaMass, &AGBAllMetals, &AGBUnProcessedMetals,
				   &SNIIRate, &SNIaRate, &AGBRate,
				   SNIIAllElements, SNIIUnProcessedElements, SNIaAllElements, SNIaUnProcessedElements, AGBAllElements, AGBUnProcessedElements);
#endif

    	//Store for use in dust_yields.c:
#ifdef DETAILED_DUST
#ifdef H2_AND_RINGS
    	DiskSNIIRate_current_ts[jj] += SNIIRate;
    	DiskSNIaRate_current_ts[jj] += SNIaRate;
    	DiskAGBRate_current_ts[jj] += AGBRate; //Not currently used, but calculated anyway. (29-11-21)
    	for(ee=0;ee<NUM_ELEMENTS;ee++) {
    		DiskSNIIAllElements_ts[jj][ee] += SNIIAllElements[ee];
    		DiskSNIaAllElements_ts[jj][ee] += SNIaAllElements[ee];
    		DiskAGBAllElements_ts[jj][ee] += AGBAllElements[ee];
    	}
#else
    	DiskSNIIRate_current_ts += SNIIRate;
    	DiskSNIaRate_current_ts += SNIaRate;
    	DiskAGBRate_current_ts += AGBRate; //Not currently used, but calculated anyway. (29-11-21)
    	for(ee=0;ee<NUM_ELEMENTS;ee++) {
    		DiskSNIIAllElements_ts[ee] += SNIIAllElements[ee];
    		DiskSNIaAllElements_ts[ee] += SNIaAllElements[ee];
    		DiskAGBAllElements_ts[ee] += AGBAllElements[ee];
    	}
#endif //H2_AND_RINGS
#endif //DETAILED_DUST

//If mass return or unprocessed metals larger than what is currently in the stellar population, re-adjust
    	double SNIImetals, SNIametals, AGBmetals;
#ifdef H2_AND_RINGS
		SNIImetals = Gal[p].MetalsDiskMassRings[jj][0];
		SNIametals = Gal[p].MetalsDiskMassRings[jj][1];
		AGBmetals = Gal[p].MetalsDiskMassRings[jj][2];
#else
		SNIImetals = Gal[p].MetalsDiskMass[0];
		SNIametals = Gal[p].MetalsDiskMass[1];
		AGBmetals = Gal[p].MetalsDiskMass[2];
#endif
    	if(SNIIUnProcessedMetals>SNIImetals)
    		SNIIUnProcessedMetals=SNIImetals;
    	if(SNIaUnProcessedMetals>SNIametals)
    		SNIaUnProcessedMetals=SNIametals;
    	if(AGBUnProcessedMetals>AGBmetals)
    		AGBUnProcessedMetals=AGBmetals;




    	mass_checks(p,"model_yields.c",__LINE__);


    	//************************
    	// UPDATE GAS COMPONENTS:
    	//************************

    	//Hot Gas:
     	//Gal[igal].HotGas += (fwind_SNII * SNIIEjectaMass) + (fwind_SNIa * SNIaEjectaMass);
     	Gal[igal].HotGas += (fwind_SNII * SNIIEjectaMass) + (fwind_SNIa * SNIaEjectaMass) + (fwind_AGB * AGBEjectaMass);
     	Gal[igal].MetalsHotGas[0] += fwind_SNII * SNIIAllMetals;
     	Gal[igal].MetalsHotGas[1] += fwind_SNIa * SNIaAllMetals;
     	Gal[igal].MetalsHotGas[2] += fwind_AGB * AGBAllMetals;
     	//Gal[igal].MetalsReheatingRate += ((fwind_SNII * SNIIAllMetals)+(fwind_SNIa * SNIaAllMetals)+(fwind_AGB * AGBAllMetals)) / (dt*STEPS); //ROB: Adding the rate of direct metal transfer from SNe/AGBs to the HotGas in MetalsReheatingRate here. (18-08-20)
     	mass_checks(p,"model_yields.c",__LINE__);
     	//Cold Gas:
     	//Gal[p].ColdGas += (1.0-fwind_SNII)*SNIIEjectaMass + (1.0-fwind_SNIa)*SNIaEjectaMass + AGBEjectaMass;
     	Gal[p].ColdGas += (1.0-fwind_SNII)*SNIIEjectaMass + (1.0-fwind_SNIa)*SNIaEjectaMass + (1.0-fwind_AGB)*AGBEjectaMass;
     	Gal[p].MetalsColdGas[0] += (1.0-fwind_SNII) * SNIIAllMetals;
     	Gal[p].MetalsColdGas[1] += (1.0-fwind_SNIa) * SNIaAllMetals;
     	//Gal[p].MetalsColdGas[2] += AGBAllMetals;
     	Gal[p].MetalsColdGas[2] += (1.0-fwind_AGB) * AGBAllMetals;

#ifdef H2_AND_RINGS
     	/*Gal[p].ColdGasRings[jj] += ((1.0-fwind_SNII) * SNIIEjectaMass + (1.0-fwind_SNIa) * SNIaEjectaMass + AGBEjectaMass);
     	Gal[p].MetalsColdGasRings[jj][0] += (1.0-fwind_SNII) * SNIIAllMetals;
     	Gal[p].MetalsColdGasRings[jj][1] += (1.0-fwind_SNIa) * SNIaAllMetals;
     	Gal[p].MetalsColdGasRings[jj][2] +=  AGBAllMetals;*/
     	Gal[p].ColdGasRings[jj] += ((1.0-fwind_SNII) * SNIIEjectaMass + (1.0-fwind_SNIa) * SNIaEjectaMass + (1.0-fwind_AGB) * AGBEjectaMass);
     	Gal[p].MetalsColdGasRings[jj][0] += (1.0-fwind_SNII) * SNIIAllMetals;
     	Gal[p].MetalsColdGasRings[jj][1] += (1.0-fwind_SNIa) * SNIaAllMetals;
     	Gal[p].MetalsColdGasRings[jj][2] += (1.0-fwind_AGB) * AGBAllMetals;
#endif
    	//Total:
     	//TotalMassReturnedToHotGas += (fwind_SNII * SNIIEjectaMass) + (fwind_SNIa * SNIaEjectaMass);
     	//TotalMassReturnedToColdDiskGas += (1.0-fwind_SNII) * SNIIEjectaMass + (1.0-fwind_SNIa) * SNIaEjectaMass + AGBEjectaMass; //Only use energy from SNe that eject into ColdGas to reheat
     	TotalMassReturnedToHotGas += (fwind_SNII * SNIIEjectaMass) + (fwind_SNIa * SNIaEjectaMass) + (fwind_AGB * AGBEjectaMass);
     	TotalMetalsReturnedToHotGas += (fwind_SNII * SNIIAllMetals) + (fwind_SNIa * SNIaAllMetals) + (fwind_AGB * AGBAllMetals);
     	TotalMassReturnedToColdDiskGas += (1.0-fwind_SNII) * SNIIEjectaMass + (1.0-fwind_SNIa) * SNIaEjectaMass + (1.0-fwind_AGB) * AGBEjectaMass; //Only use energy from SNe that eject into ColdGas to reheat
#ifdef H2_AND_RINGS
     	//TotalMassReturnedToColdDiskGasr[jj] += ((1.0-fwind_SNII) * SNIIEjectaMass + (1.0-fwind_SNIa) * SNIaEjectaMass + AGBEjectaMass);
     	TotalMassReturnedToColdDiskGasr[jj] += ((1.0-fwind_SNII) * SNIIEjectaMass + (1.0-fwind_SNIa) * SNIaEjectaMass + (1.0-fwind_AGB) * AGBEjectaMass);
#endif //H2_AND_RINGS



#ifdef INDIVIDUAL_ELEMENTS
     	for(ee=0;ee<NUM_ELEMENTS;ee++)
     	  {
     	    //Gal[igal].HotGas_elements[ee] += fwind_SNII * SNIIAllElements[ee] + fwind_SNIa * SNIaAllElements[ee];
     	    //Gal[p].ColdGas_elements[ee] += (1.0-fwind_SNII) * SNIIAllElements[ee] + (1.0-fwind_SNIa) * SNIaAllElements[ee] + AGBAllElements[ee];
     		Gal[igal].HotGas_elements[ee] += fwind_SNII * SNIIAllElements[ee] + fwind_SNIa * SNIaAllElements[ee] + fwind_AGB * AGBAllElements[ee];
     		Gal[p].ColdGas_elements[ee] += (1.0-fwind_SNII) * SNIIAllElements[ee] + (1.0-fwind_SNIa) * SNIaAllElements[ee] + (1.0-fwind_AGB) * AGBAllElements[ee];
#ifdef H2_AND_RINGS
     	    Gal[p].ColdGasRings_elements[jj][ee] += (1.0-fwind_SNII) * SNIIAllElements[ee] + (1.0-fwind_SNIa) * SNIaAllElements[ee] + (1.0-fwind_AGB) * AGBAllElements[ee];
#endif // H2_AND_RINGS
     	  }
#ifdef DETAILED_DUST
#ifdef H2_AND_RINGS
#ifndef MAINELEMENTS
     		SNII_prevstep_Cold_Cb[jj][ii] += (1.0-fwind_SNII) * SNIIAllElements[Cb_NUM];
     		SNII_prevstep_Cold_Si[jj][ii] += (1.0-fwind_SNII) * SNIIAllElements[Si_NUM];
#endif //MAINELEMENTS
     		SNII_prevstep_Cold_Fe[jj][ii] += (1.0-fwind_SNII) * SNIIAllElements[Fe_NUM];
     		SNIa_prevstep_Cold_Fe[jj][ii] += (1.0-fwind_SNIa) * SNIaAllElements[Fe_NUM];
#else //H2_AND_RINGS
#ifndef MAINELEMENTS
     		SNII_prevstep_Cold_Cb[ii] += (1.0-fwind_SNII) * SNIIAllElements[Cb_NUM];
     		SNII_prevstep_Cold_Si[ii] += (1.0-fwind_SNII) * SNIIAllElements[Si_NUM];
#endif //MAINELEMENTS
     		SNII_prevstep_Cold_Fe[ii] += (1.0-fwind_SNII) * SNIIAllElements[Fe_NUM];
     		SNIa_prevstep_Cold_Fe[ii] += (1.0-fwind_SNIa) * SNIaAllElements[Fe_NUM];
#endif //H2_AND_RINGS
#endif //DETAILED_DUST
#endif //INDIVIDUAL_ELEMENTS

     	mass_checks(p,"model_yields.c",__LINE__);


     	//*****************************
     	// UPDATE DISK MASS COMPONENTS:
     	//*****************************
     	/*ROB (13-02-13): All the mass/metals/elements in the stars that die in this timestep are lost from the stellar component.
     	 * i.e. All the mass/metals/elements in the stars at birth are removed...
     	 *...Some goes to the gas (+ newly synthesised component), the rest goes into the 'stellar remnants' which are not tracked and do not contribute to the stellar component's mass/metals/elements budget.*/
     	Gal[p].DiskMass -= SNIIEjectaMass + SNIaEjectaMass + AGBEjectaMass;
     	Gal[p].MetalsDiskMass[0] -= SNIIUnProcessedMetals;
     	Gal[p].MetalsDiskMass[1] -= SNIaUnProcessedMetals;
     	Gal[p].MetalsDiskMass[2] -= AGBUnProcessedMetals;
     	Gal[p].DiskSNIIRate += SNIIRate/STEPS; //in 1/(time code units). This records the average rate in each snapshot [i.e. SUM(SNIIRate per timestep)/STEPS], so is set to 0.0 at the start of each snapshot in model_misc.c or main.c (I forget which!). (22-05-20)

#ifdef H2_AND_RINGS
     	Gal[p].DiskMassRings[jj]-= (SNIIEjectaMass + SNIaEjectaMass + AGBEjectaMass);
     	Gal[p].MetalsDiskMassRings[jj][0]-= SNIIUnProcessedMetals;
     	Gal[p].MetalsDiskMassRings[jj][1]-= SNIaUnProcessedMetals;
     	Gal[p].MetalsDiskMassRings[jj][2]-= AGBUnProcessedMetals;
#endif
	    mass_checks(p,"model_yields.c",__LINE__);
#ifdef INDIVIDUAL_ELEMENTS

     	for(ee=0;ee<NUM_ELEMENTS;ee++)
     	  {
     	    Gal[p].DiskMass_elements[ee] -= SNIIUnProcessedElements[ee]+SNIaUnProcessedElements[ee]+AGBUnProcessedElements[ee];
#ifdef H2_AND_RINGS
     	    Gal[p].DiskMassRings_elements[jj][ee] -= (SNIIUnProcessedElements[ee]+SNIaUnProcessedElements[ee]+AGBUnProcessedElements[ee]);
#endif
     	  }
#endif //INDIVIDUAL_ELEMENTS

     	//Update ages:
     	for(n=0;n<NOUT;n++)
     	  {
    	    AgeCorrectionDisk[n] += (sfh_time-NumToTime(ListOutputSnaps[n]))*(SNIIEjectaMass + SNIaEjectaMass + AGBEjectaMass);
    	    if (AgeCorrectionDisk[n] < 0.0) AgeCorrectionDisk[n] = 0.0;
    	  }

    	mass_checks(p,"model_yields.c",__LINE__);
    } //if (Gal[p].sfh_DiskMass[i] > 0.0) -> all disk properties updated













    //******************************************************
    //ENRICHMENT FROM BULGE STARS (INTO COLD GAS & HOT GAS):
    //******************************************************
#ifdef H2_AND_RINGS
    	for(jj=0;jj<RNUM;jj++)
    		if (Gal[p].sfh_BulgeMassRings[jj][ii] > 0.0)
    		{
    			BulgeSFRxStep = timestep_width * Gal[p].sfh_BulgeMassRings[jj][ii]/Gal[p].sfh_dt[ii];
    			BulgeSFRxStep_Phys = BulgeSFRxStep * (1.0e10/Hubble_h);
    			//printf("sfr=%0.5e\n",DiskSFRxStep_Phys);
    			BulgeMetallicity = 0.;
    			for(mm=0;mm<NUM_METAL_CHANNELS;mm++)
    				BulgeMetallicity += Gal[p].sfh_MetalsBulgeMassRings[jj][ii][mm];
    			BulgeMetallicity /= Gal[p].sfh_BulgeMassRings[jj][ii];

    	    	//***** FOR TESTING FIXED-METALLICITY YIELDS: ********
    			//BulgeMetallicity = 0.0134; //Solar metallicity (Asplund+09)
    	    	//****************************************************

#ifdef INDIVIDUAL_ELEMENTS
    			for (ee=0;ee<NUM_ELEMENTS;ee++)
    				BulgeMetallicityElement_Phys[ee] = Gal[p].sfh_BulgeMass_elementsRings[jj][ii][ee] / (Gal[p].sfh_BulgeMassRings[jj][ii]*1.0e10/Hubble_h);
#endif

#else //H2_AND_RINGS

    	if (Gal[p].sfh_BulgeMass[ii] > 0.0)
    	  {
    		BulgeSFRxStep = timestep_width * Gal[p].sfh_BulgeMass[ii]/Gal[p].sfh_dt[ii];
    		BulgeSFRxStep_Phys = BulgeSFRxStep * (1.0e10/Hubble_h);

    		BulgeMetallicity = 0.;
    		for(mm=0;mm<NUM_METAL_CHANNELS;mm++)
    			BulgeMetallicity += Gal[p].sfh_MetalsBulgeMass[ii][mm];
    		BulgeMetallicity /= Gal[p].sfh_BulgeMass[ii];


#ifdef INDIVIDUAL_ELEMENTS
    		for (ee=0;ee<NUM_ELEMENTS;ee++)
    			BulgeMetallicityElement_Phys[ee] = Gal[p].sfh_BulgeMass_elements[ii][ee] / (Gal[p].sfh_BulgeMass[ii]*1.0e10/Hubble_h);
#endif

#endif //H2_AND_RINGS

#ifdef BINARYC
    	Zi = find_initial_metallicity(BulgeMetallicity, 1, 5);
#else //BINARYC
    	Zi = find_initial_metallicity(BulgeMetallicity, 1, 2);
#endif //BINARYC
    	//Interpolate the bulge luminosity on the lifetimeMetallicities tables:
    	Zi_disp = (BulgeMetallicity - lifetimeMetallicities[Zi])/(lifetimeMetallicities[Zi+1] - lifetimeMetallicities[Zi]);
    	if (Zi_disp < 0.0) Zi_disp = 0.0; //Don't want to extrapolate yields down below lifetimeMetallicities[0]=0.0004. Instead, assume constant yield below this metallicity.


#ifndef INDIVIDUAL_ELEMENTS
    	compute_actual_eject_rates(TimeBin, ii, Zi, Zi_disp, Gal[p].sfh_ibin, BulgeSFRxStep, BulgeSFRxStep_Phys, BulgeMetallicity, timestep_width,
				   &SNIIEjectaMass, &SNIIAllMetals, &SNIIUnProcessedMetals,
				   &SNIaEjectaMass, &SNIaAllMetals, &SNIaUnProcessedMetals,
				   &AGBEjectaMass, &AGBAllMetals, &AGBUnProcessedMetals,
				   &SNIIRate, &SNIaRate, &AGBRate);

#else
    	compute_actual_eject_rates(TimeBin, ii, Zi, Zi_disp, Gal[p].sfh_ibin, BulgeSFRxStep, BulgeSFRxStep_Phys, BulgeMetallicity, BulgeMetallicityElement_Phys, timestep_width,
				   &SNIIEjectaMass, &SNIIAllMetals, &SNIIUnProcessedMetals,
				   &SNIaEjectaMass, &SNIaAllMetals, &SNIaUnProcessedMetals,
				   &AGBEjectaMass, &AGBAllMetals, &AGBUnProcessedMetals,
				   &SNIIRate, &SNIaRate, &AGBRate,
				   SNIIAllElements, SNIIUnProcessedElements, SNIaAllElements, SNIaUnProcessedElements, AGBAllElements, AGBUnProcessedElements);
#endif

    	//If mass return or unprocessed metals larger than what is currently in the stellar population, re-adjust
    	double SNIImetals, SNIametals, AGBmetals;

#ifdef H2_AND_RINGS
    	SNIImetals = Gal[p].MetalsBulgeMassRings[jj][0];
    	SNIametals = Gal[p].MetalsBulgeMassRings[jj][1];
    	AGBmetals = Gal[p].MetalsBulgeMassRings[jj][2];
#else
    	SNIImetals = Gal[p].MetalsBulgeMass[0];
    	SNIametals = Gal[p].MetalsBulgeMass[1];
    	AGBmetals = Gal[p].MetalsBulgeMass[2];
#endif
    	if(SNIIUnProcessedMetals>SNIImetals)
    		SNIIUnProcessedMetals=SNIImetals;
    	if(SNIaUnProcessedMetals>SNIametals)
    		SNIaUnProcessedMetals=SNIametals;
    	if(AGBUnProcessedMetals>AGBmetals)
    		AGBUnProcessedMetals=AGBmetals;



    	//*****************************
    	// UPDATE HOT GAS COMPONENTS:
    	//*****************************
#ifndef BULGE_TO_COLD
    	Gal[igal].HotGas += SNIIEjectaMass + SNIaEjectaMass + AGBEjectaMass;
    	//If there was no hotgas left in the galaxy it will probably be stripped next step.
    	//Give it a fake HotRadius for now to avoid crash at mass checks
    	if(Gal[igal].HotGas>0. && Gal[igal].HotRadius==0)
    	  Gal[igal].HotRadius=1.e-10;
    	TotalMassReturnedToHotGas += SNIIEjectaMass + SNIaEjectaMass + AGBEjectaMass;
    	//TotalMetalsReturnedToHotGas += SNIIAllMetals + SNIaAllMetals + AGBAllMetals;
    	Gal[igal].MetalsHotGas[0] += SNIIAllMetals;
    	Gal[igal].MetalsHotGas[1] += SNIaAllMetals;
    	Gal[igal].MetalsHotGas[2] += AGBAllMetals;

#ifdef INDIVIDUAL_ELEMENTS
    	for(ee=0;ee<NUM_ELEMENTS;ee++)
    	  Gal[igal].HotGas_elements[ee] += SNIIAllElements[ee]+SNIaAllElements[ee]+AGBAllElements[ee];
#endif //INDIVIDUAL_ELEMENTS



#else //BULGE_TO_COLD
    	Gal[p].ColdGas += SNIIEjectaMass + SNIaEjectaMass + AGBEjectaMass;
    	Gal[p].MetalsColdGas[0] += SNIIAllMetals;
    	Gal[p].MetalsColdGas[1] += SNIaAllMetals;
    	Gal[p].MetalsColdGas[2] += AGBAllMetals;

#ifdef H2_AND_RINGS
    	Gal[p].ColdGasRings[jj] += (SNIIEjectaMass + SNIaEjectaMass + AGBEjectaMass);
    	Gal[p].MetalsColdGasRings[jj][0] += SNIIAllMetals;
    	Gal[p].MetalsColdGasRings[jj][1] += SNIaAllMetals;
    	Gal[p].MetalsColdGasRings[jj][2] += AGBAllMetals;
#endif

     	TotalMassReturnedToColdDiskGas += SNIIEjectaMass + SNIaEjectaMass + AGBEjectaMass;
#ifdef H2_AND_RINGS
     	TotalMassReturnedToColdDiskGasr[jj]+= (SNIIEjectaMass + SNIaEjectaMass + AGBEjectaMass);
#endif

#ifdef INDIVIDUAL_ELEMENTS
     	for(ee=0;ee<NUM_ELEMENTS;ee++)
     	  {
     	    Gal[p].ColdGas_elements[ee] += SNIIAllElements[ee]+SNIaAllElements[ee]+AGBAllElements[ee];
#ifdef H2_AND_RINGS
     	    Gal[p].ColdGasRings_elements[jj][ee]  += (SNIIAllElements[ee] + SNIaAllElements[ee] + AGBAllElements[ee]);
#endif//H2_AND_RINGS
     	  }
#endif //INDIVIDUAL_ELEMENTS

#endif //BULGE_TO_COLD







    	//*****************************
    	// UPDATE BULGE MASS COMPONENTS:
    	//*****************************
    	Gal[p].BulgeMass -= SNIIEjectaMass + SNIaEjectaMass + AGBEjectaMass;
    	Gal[p].MetalsBulgeMass[0] -= SNIIUnProcessedMetals;
    	Gal[p].MetalsBulgeMass[1] -= SNIaUnProcessedMetals;
    	Gal[p].MetalsBulgeMass[2] -= AGBUnProcessedMetals;
    	Gal[p].BulgeSNIIRate += SNIIRate/STEPS; //in 1/(time code units).

#ifdef H2_AND_RINGS
    	Gal[p].BulgeMassRings[jj]-= (SNIIEjectaMass + SNIaEjectaMass + AGBEjectaMass);
    	Gal[p].MetalsBulgeMassRings[jj][0]-= SNIIUnProcessedMetals;
    	Gal[p].MetalsBulgeMassRings[jj][1]-= SNIaUnProcessedMetals;
    	Gal[p].MetalsBulgeMassRings[jj][2]-= AGBUnProcessedMetals;
#endif

#ifdef INDIVIDUAL_ELEMENTS
    	for(ee=0;ee<NUM_ELEMENTS;ee++)
    	  {
    	    Gal[p].BulgeMass_elements[ee] -= SNIIUnProcessedElements[ee]+SNIaUnProcessedElements[ee]+AGBUnProcessedElements[ee];
#ifdef H2_AND_RINGS
    	    Gal[p].BulgeMassRings_elements[jj][ee] -= (SNIIUnProcessedElements[ee]+SNIaUnProcessedElements[ee]+AGBUnProcessedElements[ee]);
#endif
    	  }
#endif //INDIVIDUAL_ELEMENTS

    	//Update ages:
        for(n=0;n<NOUT;n++)
        {
        	AgeCorrectionBulge[n] += (sfh_time-NumToTime(ListOutputSnaps[n]))*(SNIIEjectaMass + SNIaEjectaMass + AGBEjectaMass);
        	if (AgeCorrectionBulge[n] < 0.0) AgeCorrectionBulge[n] = 0.0;
        }

        mass_checks(p,"model_yields.c",__LINE__);
    } //if (Gal[p].sfh_BulgeMass[i] > 0.0) //all BULGE properties updated









    //*****************************************
    //ENRICHMENT FROM ICL STARS INTO HOT GAS:
    //*****************************************

    if (Gal[p].sfh_ICM[ii] > 0.0)
      {
    	//pre-calculations to speed up the code
    	//Note: This is NOT really an SFR, as no stars are formed in the ICM. Rather, this is a star-transfer rate from satellite disruption to the stellar halo.
    	ICMSFRxStep = timestep_width * Gal[p].sfh_ICM[ii]/Gal[p].sfh_dt[ii];
    	ICMSFRxStep_Phys = ICMSFRxStep * (1.0e10/Hubble_h) ;
    	ICMMetallicity=0.;
    	for(mm=0;mm<NUM_METAL_CHANNELS;mm++)
    	  ICMMetallicity += Gal[p].sfh_MetalsICM[ii][mm];
    	ICMMetallicity /= Gal[p].sfh_ICM[ii];

    	//***** FOR TESTING FIXED-METALLICITY YIELDS: ********
    	//ICMMetallicity = 0.0134; //Solar metallicity (Asplund+09)
    	//****************************************************

#ifdef INDIVIDUAL_ELEMENTS
    	for (ee=0;ee<NUM_ELEMENTS;ee++)
    	  ICMMetallicityElement_Phys[ee] = Gal[p].sfh_ICM_elements[ii][ee] / (Gal[p].sfh_ICM[ii]*1.0e10/Hubble_h);
#endif

#ifdef BINARYC
    	Zi = find_initial_metallicity(ICMMetallicity, 1, 5);
#else //BINARYC
    	Zi = find_initial_metallicity(ICMMetallicity, 1, 3);
#endif //BINARYC
    	//Interpolate the ICM metallicity on the lifetimeMetallicities tables:
    	Zi_disp = (ICMMetallicity - lifetimeMetallicities[Zi])/(lifetimeMetallicities[Zi+1] - lifetimeMetallicities[Zi]);
    	if (Zi_disp < 0.0) Zi_disp = 0.0; //Don't want to extrapolate yields down below lifetimeMetallicities[0]=0.0004. Instead, assume constant yield below this metallicity.


#ifndef INDIVIDUAL_ELEMENTS
    	compute_actual_eject_rates(TimeBin, ii, Zi, Zi_disp, Gal[p].sfh_ibin, ICMSFRxStep, ICMSFRxStep_Phys, ICMMetallicity, timestep_width,
				   &SNIIEjectaMass, &SNIIAllMetals, &SNIIUnProcessedMetals,
				   &SNIaEjectaMass, &SNIaAllMetals, &SNIaUnProcessedMetals,
				   &AGBEjectaMass, &AGBAllMetals, &AGBUnProcessedMetals,
				   &SNIIRate, &SNIaRate, &AGBRate);

#else
    	compute_actual_eject_rates(TimeBin, ii, Zi, Zi_disp, Gal[p].sfh_ibin, ICMSFRxStep, ICMSFRxStep_Phys, ICMMetallicity, ICMMetallicityElement_Phys, timestep_width,
				   &SNIIEjectaMass, &SNIIAllMetals, &SNIIUnProcessedMetals,
				   &SNIaEjectaMass, &SNIaAllMetals, &SNIaUnProcessedMetals,
				   &AGBEjectaMass, &AGBAllMetals, &AGBUnProcessedMetals,
				   &SNIIRate, &SNIaRate, &AGBRate,
				   SNIIAllElements, SNIIUnProcessedElements, SNIaAllElements, SNIaUnProcessedElements, AGBAllElements, AGBUnProcessedElements);
#endif

    	//If mass return or unprocessed metals larger than what is currently in the stellar population, re-adjust
    	double SNIImetals, SNIametals, AGBmetals;

    	SNIImetals = Gal[p].MetalsICM[0];
    	SNIametals = Gal[p].MetalsICM[1];
    	AGBmetals = Gal[p].MetalsICM[2];

    	if(SNIIUnProcessedMetals>SNIImetals)
    		SNIIUnProcessedMetals=SNIImetals;
    	if(SNIaUnProcessedMetals>SNIametals)
    		SNIaUnProcessedMetals=SNIametals;
    	if(AGBUnProcessedMetals>AGBmetals)
    		AGBUnProcessedMetals=AGBmetals;

    	//*****************************
    	// UPDATE HOT GAS COMPONENTS:
    	//*****************************
    	Gal[igal].HotGas += SNIIEjectaMass + SNIaEjectaMass + AGBEjectaMass;
    	//If there was no hotgas left in the galaxy it will probably be stripped next step.
    	//Give it a fake HotRadius for now to avoid crash at mass checks
    	if(Gal[igal].HotGas>0. && Gal[igal].HotRadius==0)
    	  Gal[igal].HotRadius=1.e-10;

    	TotalMassReturnedToHotGas += SNIIEjectaMass + SNIaEjectaMass + AGBEjectaMass;
    	//TotalMetalsReturnedToHotGas += SNIIAllMetals + SNIaAllMetals + AGBAllMetals;
    	Gal[igal].MetalsHotGas[0] += SNIIAllMetals;
    	Gal[igal].MetalsHotGas[1] += SNIaAllMetals;
    	Gal[igal].MetalsHotGas[2] += AGBAllMetals;

#ifdef INDIVIDUAL_ELEMENTS
    	for(ee=0;ee<NUM_ELEMENTS;ee++)
    	  Gal[igal].HotGas_elements[ee] += SNIIAllElements[ee] + SNIaAllElements[ee] + AGBAllElements[ee];
#endif //INDIVIDUAL_ELEMENTS


    	//*****************************
    	// UPDATE ICM COMPONENTS:
    	//*****************************
    	Gal[p].ICM -= SNIIEjectaMass + SNIaEjectaMass + AGBEjectaMass;
    	Gal[p].MetalsICM[0] -= SNIIUnProcessedMetals;
    	Gal[p].MetalsICM[1] -= SNIaUnProcessedMetals;
    	Gal[p].MetalsICM[2] -= AGBUnProcessedMetals;
    	Gal[p].ICMSNIIRate += SNIIRate/STEPS; //in 1/(time code units).

#ifdef INDIVIDUAL_ELEMENTS
    	for(ee=0;ee<NUM_ELEMENTS;ee++)
    	  Gal[p].ICM_elements[ee] -= SNIIUnProcessedElements[ee] + SNIaUnProcessedElements[ee] + AGBUnProcessedElements[ee];
#endif //INDIVIDUAL_ELEMENTS

    	//Update ages:
        //for(n=0;n<NOUT;n++)
        //{
        //	AgeCorrectionICM[n] += (sfh_time-NumToTime(ListOutputSnaps[n]))*(ICMSFRxStep * NormMassEjecRateAllTypes);
        //}

    	mass_checks(p,"model_yields.c",__LINE__);
      } //if (Gal[p].sfh_ICM[i] > 0.0) //all ICM properties updated



    } //for (i=0;i<=Gal[p].sfh_ibin;i++) //MAIN LOOP OVER SFH BINS








    //Update Mass-weighted ages:
    for(n=0;n<NOUT;n++)
      {
    	Gal[p].MassWeightAge[n] -= (AgeCorrectionDisk[n]+AgeCorrectionBulge[n]);
      }

#ifdef DETAILED_DUST
    //Apportion the newly-released metals into the cloud and diffuse phases, in proportion to the H2 fraction:
    //partition_gas_and_dust_elements(p);
#ifdef H2_AND_RINGS
    for(ee=0;ee<NUM_ELEMENTS;ee++) {
    	for(jj=0;jj<RNUM;jj++) {
			coldgas_add = Gal[p].ColdGasRings_elements[jj][ee] - coldgas_init[jj][ee];
			Gal[p].ColdGasCloudsRings_elements[jj][ee] += Gal[p].H2fractionRings[jj] * coldgas_add;
			Gal[p].ColdGasDiffRings_elements[jj][ee] += (1. - Gal[p].H2fractionRings[jj]) * coldgas_add;
			Gal[p].ColdGasClouds_elements[ee] += Gal[p].H2fractionRings[jj] * coldgas_add;
			Gal[p].ColdGasDiff_elements[ee] += (1. - Gal[p].H2fractionRings[jj]) * coldgas_add;
    	}
    }
#else
    for(ee=0;ee<NUM_ELEMENTS;ee++) {
    	coldgas_add = Gal[p].ColdGas_elements[ee] - coldgas_init[ee];
    	Gal[p].ColdGasClouds_elements[ee] += Gal[p].H2fraction * coldgas_add;
    	Gal[p].ColdGasDiff_elements[ee] += (1. - Gal[p].H2fraction) * coldgas_add;
    }
#endif //H2_AND_RINGS
#endif //DETAILED_DUST

#ifndef DETAILED_DUST
	// IF DETAILED_DUST is switched ON
	// SN-Feedback must be called AFTER BOTH model_yields AND model_dustyields
	// and now appears inside main.c if this is the case!
#ifdef FEEDBACK_COUPLED_WITH_MASS_RETURN

    if(TotalMassReturnedToColdDiskGas>0.) {
    	if(TotalMassReturnedToColdDiskGas>Gal[p].ColdGas)
    		TotalMassReturnedToColdDiskGas=Gal[p].ColdGas;
    	//ROB: Shouldn't we do the same check for TotalMassReturnedToColdDiskGasr here too? (29-11-21)
#ifndef H2_AND_RINGS
    	SN_feedback(p, centralgal, TotalMassReturnedToColdDiskGas, "ColdGas", dt);
#else
    	SN_feedback(p, centralgal, TotalMassReturnedToColdDiskGas, TotalMassReturnedToColdDiskGasr, "ColdGas", dt);
#endif
    	Gal[p].MassReturnRateToColdGas += TotalMassReturnedToColdDiskGas / (dt * STEPS); //ROB: Store mass return rate by SNe and stellar winds to the gas phases.
    }

    //this mass will only result in ejection, no reheating
    if(TotalMassReturnedToHotGas>0.) {
    	if(TotalMassReturnedToHotGas>Gal[p].HotGas)
    		TotalMassReturnedToHotGas=Gal[p].HotGas;
#ifndef H2_AND_RINGS
    	SN_feedback(p, centralgal, TotalMassReturnedToHotGas, "HotGas", dt);
#else
    	double HotGasRings[RNUM];
		for(jj=0;jj<RNUM;jj++)
			HotGasRings[jj]=0.;
			//HotGasRings[jj]=TotalMassReturnedToHotGas*fractionRingsBulge[jj];
		SN_feedback(p, centralgal, TotalMassReturnedToHotGas, HotGasRings, "HotGas", dt);
#endif
		Gal[p].MassReturnRateToHotGas += TotalMassReturnedToHotGas / (dt * STEPS); //ROB: Store mass return rate by SNe and stellar winds to the gas phases (from disk stars, bulge stars, and ICL stars).
      }

    if(TotalMetalsReturnedToHotGas>0.) {
    	//NOTE: Only metals returned to the HotGas from disc stars & SNe are considered in TotalMetalsReturnedToHotGas (i.e. not from the bulge or ICL):
    	Gal[p].MetalsReturnRateToHotGas += TotalMetalsReturnedToHotGas / (dt * STEPS); //ROB: Store mass return rate by SNe and stellar winds to the gas phases (from disk stars, bulge stars, and ICL stars).
      }

#endif //FEEDBACK_COUPLED_WITH_MASS_RETURN
#endif //DETAILED_DUST

      mass_checks(p,"model_yields.c",__LINE__);

      /*printf("End of model_yields.c:\n");
      	for (ee=0; ee<NUM_ELEMENTS; ee++)
      		printf("Gal[p].DustColdGasClouds_elements[%i] = %f\n", ee, Gal[p].DustColdGasClouds_elements[ee]);*/

      /*double diskspinpar=sqrt(Gal[p].DiskSpin[0] * Gal[p].DiskSpin[0] +
     			     Gal[p].DiskSpin[1] * Gal[p].DiskSpin[1] +
     			     Gal[p].DiskSpin[2] * Gal[p].DiskSpin[2] );
      int ii;
      if (Gal[p].ColdGas+TotalMassReturnedToColdDiskGas > 1.e-8)
        for(ii=0;ii<3;ii++)
          Gal[p].ColdGasSpin[ii]=((Gal[p].ColdGasSpin[ii])*(Gal[p].ColdGas)+1./sqrt(3)*diskspinpar*TotalMassReturnedToColdDiskGas)/(Gal[p].ColdGas+TotalMassReturnedToColdDiskGas);*/

      if (Gal[p].ColdGas+TotalMassReturnedToColdDiskGas > 1.e-8)
	for (ii = 0; ii < 3; ii++)
	  Gal[p].ColdGasSpin[ii]=((Gal[p].ColdGasSpin[ii])*(Gal[p].ColdGas)+TotalMassReturnedToColdDiskGas*Gal[p].DiskSpin[ii])/(Gal[p].ColdGas+TotalMassReturnedToColdDiskGas);

      if (DiskRadiusModel == 0)
	{
	  Gal[p].ColdGasRadius = get_gas_disk_radius(p);
	  Gal[p].DiskRadius = get_stellar_disk_radius(p);
	}
}


#ifndef INDIVIDUAL_ELEMENTS
void compute_actual_eject_rates(int TimeBin, int ii, int Zi, double Zi_disp, int sfh_ibin, double SFRxStep, double SFRxStep_Phys, double Metallicity, double timestep_width,
				 double *SNIIEjectaMass, double *SNIIAllMetals, double *SNIIUnProcessedMetals,
				 double *SNIaEjectaMass, double *SNIaAllMetals, double *SNIaUnProcessedMetals,
				 double *AGBEjectaMass, double *AGBAllMetals, double *AGBUnProcessedMetals,
				 double *SNIIRate, double *SNIaRate, double *AGBRate)
#else
void compute_actual_eject_rates(int TimeBin, int ii, int Zi, double Zi_disp, int sfh_ibin, double SFRxStep, double SFRxStep_Phys, double Metallicity, double *MetallicityElement_Phys, double timestep_width,
				 double *SNIIEjectaMass, double *SNIIAllMetals, double *SNIIUnProcessedMetals,
				 double *SNIaEjectaMass, double *SNIaAllMetals, double *SNIaUnProcessedMetals,
				 double *AGBEjectaMass, double *AGBAllMetals, double *AGBUnProcessedMetals,
				 double *SNIIRate, double *SNIaRate, double *AGBRate,
				 double *SNIIAllElements, double *SNIIUnProcessedElements,
				 double *SNIaAllElements, double *SNIaUnProcessedElements,
				 double *AGBAllElements, double *AGBUnProcessedElements)
#endif
{
  double NormSNIIMassEjecRate_actual, NormSNIaMassEjecRate_actual, NormAGBMassEjecRate_actual;
  double NormSNIIMetalEjecRate_actual, NormSNIaMetalEjecRate_actual, NormAGBMetalEjecRate_actual;
  double NormSNIINum_actual, NormSNIaNum_actual, NormAGBNum_actual;
#ifdef INDIVIDUAL_ELEMENTS
  int ee;
  double NormSNIIYieldRate_actual[NUM_ELEMENTS], NormSNIaYieldRate_actual[NUM_ELEMENTS], NormAGBYieldRate_actual[NUM_ELEMENTS];
#endif

//pre-calculations to speed up the code
  NormSNIIMassEjecRate_actual = NormSNIIMassEjecRate[TimeBin][ii][Zi] + ((NormSNIIMassEjecRate[TimeBin][ii][Zi+1] - NormSNIIMassEjecRate[TimeBin][ii][Zi])*Zi_disp);
  NormSNIaMassEjecRate_actual = NormSNIaMassEjecRate[TimeBin][ii][Zi] + ((NormSNIaMassEjecRate[TimeBin][ii][Zi+1] - NormSNIaMassEjecRate[TimeBin][ii][Zi])*Zi_disp);
  NormAGBMassEjecRate_actual = NormAGBMassEjecRate[TimeBin][ii][Zi] + ((NormAGBMassEjecRate[TimeBin][ii][Zi+1] - NormAGBMassEjecRate[TimeBin][ii][Zi])*Zi_disp);
  NormSNIIMetalEjecRate_actual = NormSNIIMetalEjecRate[TimeBin][ii][Zi] + ((NormSNIIMetalEjecRate[TimeBin][ii][Zi+1] - NormSNIIMetalEjecRate[TimeBin][ii][Zi])*Zi_disp);
  NormSNIaMetalEjecRate_actual = NormSNIaMetalEjecRate[TimeBin][ii][Zi] + ((NormSNIaMetalEjecRate[TimeBin][ii][Zi+1] - NormSNIaMetalEjecRate[TimeBin][ii][Zi])*Zi_disp);
  NormAGBMetalEjecRate_actual = NormAGBMetalEjecRate[TimeBin][ii][Zi] + ((NormAGBMetalEjecRate[TimeBin][ii][Zi+1] - NormAGBMetalEjecRate[TimeBin][ii][Zi])*Zi_disp);
  NormSNIINum_actual = NormSNIINum[TimeBin][ii][Zi] + ((NormSNIINum[TimeBin][ii][Zi+1] - NormSNIINum[TimeBin][ii][Zi])*Zi_disp);
  NormSNIaNum_actual = NormSNIaNum[TimeBin][ii][Zi] + ((NormSNIaNum[TimeBin][ii][Zi+1] - NormSNIaNum[TimeBin][ii][Zi])*Zi_disp);
  NormAGBNum_actual = NormAGBNum[TimeBin][ii][Zi] + ((NormAGBNum[TimeBin][ii][Zi+1] - NormAGBNum[TimeBin][ii][Zi])*Zi_disp);

#ifdef INDIVIDUAL_ELEMENTS //Work out the actual yield of element k, by interpolating between the yields in the look-up table created by yield_integrals.c.
  for (ee=0;ee<NUM_ELEMENTS;ee++)
    {
      NormSNIIYieldRate_actual[ee] = NormSNIIYieldRate[TimeBin][ii][Zi][ee] + ((NormSNIIYieldRate[TimeBin][ii][Zi+1][ee] - NormSNIIYieldRate[TimeBin][ii][Zi][ee])*Zi_disp);
      NormSNIaYieldRate_actual[ee] = NormSNIaYieldRate[TimeBin][ii][Zi][ee] + ((NormSNIaYieldRate[TimeBin][ii][Zi+1][ee] - NormSNIaYieldRate[TimeBin][ii][Zi][ee])*Zi_disp);
      NormAGBYieldRate_actual[ee] = NormAGBYieldRate[TimeBin][ii][Zi][ee] + ((NormAGBYieldRate[TimeBin][ii][Zi+1][ee] - NormAGBYieldRate[TimeBin][ii][Zi][ee])*Zi_disp);
    }
#endif

#ifdef FAST_TESTING_MODE
    	/* if(NormAGBMassEjecRate_actual>0.1)
    	          NormAGBMassEjecRate_actual=0.01;

    	        if(NormSNIIMassEjecRate_actual+NormSNIaMassEjecRate_actual+NormAGBMassEjecRate_actual>1.)
    	          {
    	            NormSNIIMassEjecRate_actual=0.43;
    	            NormSNIaMassEjecRate_actual=0.3;
    	            NormAGBMassEjecRate_actual=0.1;
    	          }*/
#endif

#ifdef INSTANTANEOUS_RECYCLE //to recover results from instantaneous recycling approximation
  reset_ejection_rates(ii, sfh_ibin,
		       &NormSNIIMassEjecRate_actual, &NormSNIIMetalEjecRate_actual,
		       &NormSNIaMassEjecRate_actual, &NormAGBMassEjecRate_actual,
		       &NormSNIaMetalEjecRate_actual, &NormAGBMetalEjecRate_actual);
#endif //INSTANTANEOUS_RECYCLE


  //SNII
  *SNIIEjectaMass = SFRxStep * NormSNIIMassEjecRate_actual;
#ifdef INSTANTANEOUS_RECYCLE
  *SNIIAllMetals = SFRxStep * (NormSNIIMetalEjecRate_actual + Metallicity * NormSNIIMassEjecRate_actual);
#else
#ifdef PORTINARI
  *SNIIAllMetals = SFRxStep * (NormSNIIMetalEjecRate_actual + Metallicity * NormSNIIMassEjecRate_actual);
#endif
#ifdef CHIEFFI
  //ROB: No unsynth component required for SN-II ejecta, when using the Chieffi & Limongi (2007) yield tables/
  *SNIIAllMetals = SFRxStep * NormSNIIMetalEjecRate_actual;
#endif
#endif //INSTANTANEOUS_RECYCLE
  *SNIIUnProcessedMetals = SFRxStep * Metallicity * NormSNIIMassEjecRate_actual;
  ////*SNIIRate = (SFRxStep/timestep_width) * NormSNIINum_actual; //This is the SNII rate in [1/(10^10/h)(time code units)] in this timestep. So should divide by ((UnitTime_in_s / SEC_PER_YEAR) * (1.e10/Hubble_h)) in save.c to get "per year" before outputting).
  *SNIIRate = (SFRxStep_Phys/timestep_width) * NormSNIINum_actual; //This is the SNII rate in [1/(time code units)] in this timestep. So should divide by (UnitTime_in_s / SEC_PER_YEAR) in save.c to get "per year" before outputting).

  //SNIa
  *SNIaEjectaMass = SFRxStep * NormSNIaMassEjecRate_actual;
  //SNIa yields are written in a different way to SNIIa and AGB so the following line is correct
  *SNIaAllMetals = SFRxStep * NormSNIaMetalEjecRate_actual;
  *SNIaUnProcessedMetals = SFRxStep * Metallicity * NormSNIaMassEjecRate_actual;
  *SNIaRate = (SFRxStep_Phys/timestep_width) * NormSNIaNum_actual; //This is the SNIa rate in [1/(time code units)] in this timestep. So should divide by (UnitTime_in_s / SEC_PER_YEAR) in save.c to get "per year" before outputting).

  //AGB
  *AGBEjectaMass = SFRxStep * NormAGBMassEjecRate_actual;
  *AGBAllMetals = SFRxStep * (NormAGBMetalEjecRate_actual + (Metallicity * NormAGBMassEjecRate_actual));
  *AGBUnProcessedMetals = SFRxStep * Metallicity * NormAGBMassEjecRate_actual;
  *AGBRate = (SFRxStep_Phys/timestep_width) * NormAGBNum_actual; //This is the AGB rate in [1/(time code units)] in this timestep. So should divide by (UnitTime_in_s / SEC_PER_YEAR) in save.c to get "per year" before outputting).

#ifdef INDIVIDUAL_ELEMENTS
  for(ee=0;ee<NUM_ELEMENTS;ee++)
    {
#ifdef PORTINARI
      SNIIAllElements[ee] = SFRxStep_Phys * (NormSNIIYieldRate_actual[ee] + MetallicityElement_Phys[ee] * NormSNIIMassEjecRate_actual);
#endif
#ifdef CHIEFFI
      SNIIAllElements[ee] = SFRxStep_Phys * NormSNIIYieldRate_actual[ee];
#endif
      SNIIUnProcessedElements[ee] = SFRxStep_Phys * MetallicityElement_Phys[ee] * NormSNIIMassEjecRate_actual;

      SNIaAllElements[ee] = SFRxStep_Phys * NormSNIaYieldRate_actual[ee];
      SNIaUnProcessedElements[ee] = SFRxStep_Phys * MetallicityElement_Phys[ee] * NormSNIaMassEjecRate_actual;

      AGBAllElements[ee] = SFRxStep_Phys * (NormAGBYieldRate_actual[ee] + MetallicityElement_Phys[ee] * NormAGBMassEjecRate_actual);
      AGBUnProcessedElements[ee] = SFRxStep_Phys * MetallicityElement_Phys[ee] * NormAGBMassEjecRate_actual;
    }
#endif

}








int find_initial_metallicity(double metallicity, int table_type, int component_type)
{
	int i, Zi_bin;
	Zi_bin = -1;
	i = 0;

	switch (table_type)
	{
	case 1: //Lifetime metallicity table
	  while (Zi_bin == -1)
		{
		  if (lifetimeMetallicities[i] < metallicity)
		    {
			  i++;
			  if (i == LIFETIME_Z_NUM) Zi_bin = i; //If galaxy's Z is higher than max Z from table, then just take max Z from table
		    }
		  else Zi_bin = i;
		}
	  break;
#ifndef BINARYC
	case 2: //SN-II metallicity table
	  while (Zi_bin == -1)
	    {
		  if (SNIIMetallicities[i] < metallicity)
		  {
			  i++;
			  if (i == 5) Zi_bin = i; //If galaxy's Z is higher than max Z from table, then just take max Z from table
		  }
		  else Zi_bin = i;
	    }
	  break;
	  //case 3 //SNIa yields are NOT metallicity dependent
	case 4: //AGB metallicity table
	  while (Zi_bin == -1)
		{
		  if (AGBMetallicities[i] < metallicity)
		    {
			  i++;
			  if (i == 3) Zi_bin = i; //If galaxy's Z is higher than max Z from table, then just take max Z from table
		    }
		  else Zi_bin = i;
		}
	  break;
#else //BINARYC
	case 5: //Binary_c metallicity table
	  while (Zi_bin == -1)
		{
		  if (bcMetallicities[i] < metallicity)
		    {
			  i++;
			  if (i == BC_Z_NUM) Zi_bin = i; //If galaxy's Z is higher than max Z from table, then just take max Z from table
		    }
		  else Zi_bin = i;
		}
	  break;
#endif //BINARYC
	}

	if (Zi_bin == 0 ) return Zi_bin;
	else return Zi_bin-1;

}


#ifdef INSTANTANEOUS_RECYCLE //to recover results from instantaneous recycling approximation
 void reset_ejection_rates(int ii, int sfh_ibin,
		 double *NormSNIIMassEjecRate_actual, double *NormSNIIMetalEjecRate_actual,
		 double *NormSNIaMassEjecRate_actual, double *NormAGBMassEjecRate_actual,
		 double *NormSNIaMetalEjecRate_actual, double *NormAGBMetalEjecRate_actual)
 {
    	if(ii==sfh_ibin)
    	{
    	    *NormSNIIMassEjecRate_actual = RecycleFraction;
    	    *NormSNIIMetalEjecRate_actual = Yield;
    	}
    	else
    	{
    		*NormSNIIMassEjecRate_actual = 0.0;
    		*NormSNIIMetalEjecRate_actual = 0.0;
    	}
    	*NormSNIaMassEjecRate_actual = 0.0;
    	*NormAGBMassEjecRate_actual =  0.0;
    	*NormSNIaMetalEjecRate_actual = 0.0;
    	*NormAGBMetalEjecRate_actual =  0.0;
 }
#endif //INSTANTANEOUS_RECYCLE

