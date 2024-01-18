/*
 * model_yields.c
 *
 *  Created on: 18.11.2011
 *      Author: robyates
 *
 * Updates:
 * 10-05-22: Modified to account for binary_c yields, which need to be integrated over in time for AGBs, SNe-II, and SNe-Ia.
 * 11-10-23: Cleaned-up for uploading to GitHub.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "allvars.h"
#include "proto.h"


void update_yields_and_return_mass(int p, int centralgal, double dt, int nstep)
{
	int Zi, igal, ii, mm;
	double timestep_width; //Width of current timestep in CODE UNITS
	int TimeBin; //Bin in Yield arrays corresponding to current timestep
	double Zi_disp;
	double sfh_time;
#ifdef METALRICHWIND
#ifdef GASDENSITYFWIND
	double ColdGasSurfaceDensity;
#endif
#endif
	double DiskSFRxStep, DiskSFRxStep_Phys, BulgeSFRxStep, BulgeSFRxStep_Phys, ICMSFRxStep, ICMSFRxStep_Phys;
	double DiskMetallicity, BulgeMetallicity, ICMMetallicity;
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
			DiskSNIIAllElements_ts[jj][ee]=0.0; //Total element masses ejected by disc SNe-II in this timestep (stored for use in model_dust_yields.c)
			DiskSNIaAllElements_ts[jj][ee]=0.0; //Total element masses ejected by disc SNe-Ia in this timestep (stored for use in model_dust_yields.c)
			DiskAGBAllElements_ts[jj][ee]=0.0; //Total element masses ejected by disc AGBs in this timestep (stored for use in model_dust_yields.c)
			BulgeSNIIAllElements_ts[jj][ee]=0.0; //Total element masses ejected by bulge SNe-II in this timestep (stored for use in model_dust_yields.c)
			BulgeSNIaAllElements_ts[jj][ee]=0.0; //Total element masses ejected by bulge SNe-Ia in this timestep (stored for use in model_dust_yields.c)
			BulgeAGBAllElements_ts[jj][ee]=0.0; //Total element masses ejected by bulge AGBs in this timestep (stored for use in model_dust_yields.c)
			ICMSNIIAllElements_ts[jj][ee]=0.0; //Total element masses ejected by halo SNe-II in this timestep (stored for use in model_dust_yields.c)
			ICMSNIaAllElements_ts[jj][ee]=0.0; //Total element masses ejected by halo SNe-Ia in this timestep (stored for use in model_dust_yields.c)
			ICMAGBAllElements_ts[jj][ee]=0.0; //Total element masses ejected by halo AGBs in this timestep (stored for use in model_dust_yields.c)
		}
    }
#else //H2_AND_RINGS
    DiskSNIIRate_current_ts=0.0; //Instantaneous SNII rate in the disc (stored for use in model_dust_yields.c)
    DiskSNIaRate_current_ts=0.0; //Instantaneous SNIa rate in the disc (stored for use in model_dust_yields.c)
    DiskAGBRate_current_ts=0.0; //Instantaneous AGB rate in the disc (stored for use in model_dust_yields.c)
    for(ee=0;ee<NUM_ELEMENTS;ee++) {
    	DiskSNIIAllElements_ts[ee]=0.0; //Total element masses ejected by disc SNe-II in this timestep (stored for use in model_dust_yields.c)
    	DiskSNIaAllElements_ts[ee]=0.0; //Total element masses ejected by disc SNe-Ia in this timestep (stored for use in model_dust_yields.c)
    	DiskAGBAllElements_ts[ee]=0.0; //Total element masses ejected by disc AGBs in this timestep (stored for use in model_dust_yields.c)
    	BulgeSNIIAllElements_ts[ee]=0.0; //Total element masses ejected by bulge SNe-II in this timestep (stored for use in model_dust_yields.c)
    	BulgeSNIaAllElements_ts[ee]=0.0; //Total element masses ejected by bulge SNe-Ia in this timestep (stored for use in model_dust_yields.c)
    	BulgeAGBAllElements_ts[ee]=0.0; //Total element masses ejected by bulge AGBs in this timestep (stored for use in model_dust_yields.c)
    	ICMSNIIAllElements_ts[ee]=0.0; //Total element masses ejected by halo SNe-II in this timestep (stored for use in model_dust_yields.c)
		ICMSNIaAllElements_ts[ee]=0.0; //Total element masses ejected by halo SNe-Ia in this timestep (stored for use in model_dust_yields.c)
		ICMAGBAllElements_ts[ee]=0.0; //Total element masses ejected by halo AGBs in this timestep (stored for use in model_dust_yields.c)
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

    for (ii=0;ii<=Gal[p].sfh_ibin;ii++) //LOOP OVER SFH BINS
    {
    	sfh_time=Gal[p].sfh_t[ii]+(0.5*Gal[p].sfh_dt[ii]);

#ifdef DETAILED_DUST
// These variables store the mass of an element ejected into the ColdGas in the current timestep by disc SNe with progenitors born in each of the SFH bins. (29-11-21)
#ifdef H2_AND_RINGS
    	for(jj=0;jj<RNUM;jj++) {
			SNII_prevstep_Cold_Si[jj][ii] = 0.0;
			SNII_prevstep_Cold_Fe[jj][ii] = 0.0;
			SNII_prevstep_Cold_Cb[jj][ii] = 0.0;
			SNIa_prevstep_Cold_Fe[jj][ii] = 0.0;
			SNII_prevstep_Hot_bulge_Si[jj][ii] = 0.0;
			SNII_prevstep_Hot_bulge_Fe[jj][ii] = 0.0;
			SNII_prevstep_Hot_bulge_Cb[jj][ii] = 0.0;
			SNIa_prevstep_Hot_bulge_Fe[jj][ii] = 0.0;
			SNII_prevstep_Hot_ICM_Si[jj][ii] = 0.0;
			SNII_prevstep_Hot_ICM_Fe[jj][ii] = 0.0;
			SNII_prevstep_Hot_ICM_Cb[jj][ii] = 0.0;
			SNIa_prevstep_Hot_ICM_Fe[jj][ii] = 0.0;
    	}
#else
    	SNII_prevstep_Cold_Si[ii] = 0.0;
		SNII_prevstep_Cold_Fe[ii] = 0.0;
		SNII_prevstep_Cold_Cb[ii] = 0.0;
		SNIa_prevstep_Cold_Fe[ii] = 0.0;
		SNII_prevstep_Hot_bulge_Si[ii] = 0.0;
		SNII_prevstep_Hot_bulge_Fe[ii] = 0.0;
		SNII_prevstep_Hot_bulge_Cb[ii] = 0.0;
		SNIa_prevstep_Hot_bulge_Fe[ii] = 0.0;
		SNII_prevstep_Hot_ICM_Si[ii] = 0.0;
		SNII_prevstep_Hot_ICM_Fe[ii] = 0.0;
		SNII_prevstep_Hot_ICM_Cb[ii] = 0.0;
		SNIa_prevstep_Hot_ICM_Fe[ii] = 0.0;
#endif //H2_AND_RINGS
#endif //DETAILED_DUST

    	mass_checks(p,"model_yields.c",__LINE__);


    //******************************************************
    //ENRICHMENT FROM DISC STARS (INTO COLD GAS & HOT GAS):
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
    	fwind_SNII = min(1.0, 1.0/(ColdGasSurfaceDensity/NORMGASDENSITY)); //Fraction of SN-II ejecta put directly into HotGas. When ColdGasSurfaceDensity <= NORMGASDENSITY, then fwind_SNII = 1.0.
    	//fwind_SNIa = min(1.0, 1.0/(ColdGasSurfaceDensity/NORMGASDENSITY)); //Fraction of SN-Ia ejecta put directly into HotGas. When ColdGasSurfaceDensity <= NORMGASDENSITY, then fwind_SNIa = 1.0.
    	//fwind_AGB = min(1.0, 1.0/(ColdGasSurfaceDensity/NORMGASDENSITY)); //Fraction of AGB ejecta put directly into HotGas. When ColdGasSurfaceDensity <= NORMGASDENSITY, then fwind_AGB = 1.0.
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
#else //GASDENSITYWIND
#ifdef MVIRDEPFWIND
    	fwind_SNII = min(1.0, (FracZSNIItoHot*10.5) / log10(1.e10*Gal[p].Mvir)); //Sets direct ejection into HotGas by SNe-II to be dependent on gravitating mass (i.e. Mvir). 02-03-20
    	fwind_SNIa = min(1.0, (FracZSNIatoHot*10.5) / log10(1.e10*Gal[p].Mvir));
    	fwind_AGB = min(1.0, (FracZAGBtoHot*10.5) / log10(1.e10*Gal[p].Mvir));
#else
    	fwind_SNII = FracZSNIItoHot;
    	fwind_SNIa = FracZSNIatoHot;
    	fwind_AGB = FracZAGBtoHot;
#endif //MVIRDEPFWIND
#endif //GASDENSITYFWIND
#else //METALRICHWIND
     	fwind_SNII = 0.0; //For all stellar ejecta (from disk) to ColdGas
     	fwind_SNIa = 0.0;
     	fwind_AGB = 0.0;
#endif //METALRICHWIND

    	//pre-calculations to speed up the code
	    //timestep_width and dt units cancel out
#ifdef H2_AND_RINGS
    	DiskSFRxStep = timestep_width * Gal[p].sfh_DiskMassRings[jj][ii]/Gal[p].sfh_dt[ii]; //[in code units].
    	DiskSFRxStep_Phys = DiskSFRxStep * (1.0e10/Hubble_h); //[in Msun].

    	DiskMetallicity = 0.;
    	for(mm=0;mm<NUM_METAL_CHANNELS;mm++)
    		DiskMetallicity += Gal[p].sfh_MetalsDiskMassRings[jj][ii][mm];
    	DiskMetallicity /= Gal[p].sfh_DiskMassRings[jj][ii];

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
    	Zi = find_initial_metallicity(DiskMetallicity, 5); //Find the metallicity bin BELOW the true metallicity in the SFH bin
    	Zi_disp = (DiskMetallicity - bcMetallicities[Zi])/(bcMetallicities[Zi+1] - bcMetallicities[Zi]); //Fractional position of the SFH metallicity between the bcMetallicities bins.
#else //BINARYC
    	Zi = find_initial_metallicity(DiskMetallicity, 1); //Find the metallicity bin BELOW the true metallicity in the SFH bin
    	Zi_disp = (DiskMetallicity - lifetimeMetallicities[Zi])/(lifetimeMetallicities[Zi+1] - lifetimeMetallicities[Zi]); //Fractional position of the SFH metallicity between the bcMetallicities bins.
#endif //BINARYC
    	if (Zi_disp < 0.0) Zi_disp = 0.0; //Don't want to extrapolate yields down below minimum metallicity considered (e.g. lifetimeMetallicities[0]=0.0004 or bcMetallicities[0]=0.0001). Instead, assume constant yield below this metallicity.

#ifdef DETAILED_DUST //For use in update_dust_mass() later on
#ifdef H2_AND_RINGS
    	Zi_disk_saved[jj][ii] = Zi;
		Zi_disk_disp_saved[jj][ii] = Zi_disp;
#else
		Zi_disk_saved[ii] = Zi;
		Zi_disk_disp_saved[ii] = Zi_disp;
#endif //H2_AND_RINGS
#endif

#ifndef INDIVIDUAL_ELEMENTS
    	compute_actual_eject_rates(TimeBin, ii, Zi, Zi_disp, Gal[p].sfh_ibin, DiskSFRxStep, DiskSFRxStep_Phys, DiskMetallicity, timestep_width,
								   &SNIIEjectaMass, &SNIIAllMetals, &SNIIUnProcessedMetals,
								   &SNIaEjectaMass, &SNIaAllMetals, &SNIaUnProcessedMetals,
								   &AGBEjectaMass, &AGBAllMetals, &AGBUnProcessedMetals,
								   &SNIIRate, &SNIaRate, &AGBRate);
#ifdef H2_AND_RINGS
    	rescale_ejected_material(p, jj, "DiskMass", SNIIEjectaMass, SNIaEjectaMass, AGBEjectaMass,
								 &SNIIUnProcessedMetals, &SNIaUnProcessedMetals, &AGBUnProcessedMetals,
								 &SNIIAllMetals, &SNIaAllMetals, &AGBAllMetals);
#else //H2_AND_RINGS
		rescale_ejected_material(p, "DiskMass", SNIIEjectaMass, SNIaEjectaMass, AGBEjectaMass,
								 &SNIIUnProcessedMetals, &SNIaUnProcessedMetals, &AGBUnProcessedMetals,
								 &SNIIAllMetals, &SNIaAllMetals, &AGBAllMetals);
#endif //H2_AND_RINGS
#else //INDIVIDUAL_ELEMENTS
    	compute_actual_eject_rates(TimeBin, ii, Zi, Zi_disp, Gal[p].sfh_ibin, DiskSFRxStep, DiskSFRxStep_Phys, DiskMetallicity, DiskMetallicityElement_Phys, timestep_width,
								   &SNIIEjectaMass, &SNIIAllMetals, &SNIIUnProcessedMetals,
								   &SNIaEjectaMass, &SNIaAllMetals, &SNIaUnProcessedMetals,
								   &AGBEjectaMass, &AGBAllMetals, &AGBUnProcessedMetals,
								   &SNIIRate, &SNIaRate, &AGBRate,
								   SNIIAllElements, SNIIUnProcessedElements, SNIaAllElements, SNIaUnProcessedElements, AGBAllElements, AGBUnProcessedElements);
#ifdef H2_AND_RINGS
    	//ROB: (08-07-22): New way to scale unprocessed components:
    	rescale_ejected_material(p, jj, "DiskMass", SNIIEjectaMass, SNIaEjectaMass, AGBEjectaMass,
								 &SNIIUnProcessedMetals, &SNIaUnProcessedMetals, &AGBUnProcessedMetals,
								 &SNIIAllMetals, &SNIaAllMetals, &AGBAllMetals,
								 SNIIAllElements, SNIaAllElements, AGBAllElements,
								 SNIIUnProcessedElements, SNIaUnProcessedElements, AGBUnProcessedElements);
#else //H2_AND_RINGS
		rescale_ejected_material(p, "DiskMass", SNIIEjectaMass, SNIaEjectaMass, AGBEjectaMass,
								 &SNIIUnProcessedMetals, &SNIaUnProcessedMetals, &AGBUnProcessedMetals,
								 &SNIIAllMetals, &SNIaAllMetals, &AGBAllMetals,
								 SNIIAllElements, SNIaAllElements, AGBAllElements,
								 SNIIUnProcessedElements, SNIaUnProcessedElements, AGBUnProcessedElements);
#endif //H2_AND_RINGS
#endif //INDIVIDUAL_ELEMENTS

    	//Store for use in dust_yields.c:
#ifdef DETAILED_DUST
#ifdef H2_AND_RINGS
    	DiskSNIIRate_current_ts[jj] += SNIIRate;
    	DiskSNIaRate_current_ts[jj] += SNIaRate;
    	DiskAGBRate_current_ts[jj] += AGBRate; //Not currently used, but calculated anyway. (29-11-21)
#ifdef INDIVIDUAL_ELEMENTS
    	for(ee=0;ee<NUM_ELEMENTS;ee++) {
    		DiskSNIIAllElements_ts[jj][ee] += SNIIAllElements[ee]; //TOTAL element masses ejected by disc SNe-II in this timestep (from all SFH bins). The fwind factors do not need to be included here (they are included later in model_dust_yields.c, where necessary).
    		DiskSNIaAllElements_ts[jj][ee] += SNIaAllElements[ee];
    		DiskAGBAllElements_ts[jj][ee] += AGBAllElements[ee];
    	}
#endif //INDIVIDUAL_ELEMENTS
#else //H2_AND_RINGS
    	DiskSNIIRate_current_ts += SNIIRate;
    	DiskSNIaRate_current_ts += SNIaRate;
    	DiskAGBRate_current_ts += AGBRate; //Not currently used, but calculated anyway. (29-11-21)
#ifdef INDIVIDUAL_ELEMENTS
    	for(ee=0;ee<NUM_ELEMENTS;ee++) {
    		DiskSNIIAllElements_ts[ee] += SNIIAllElements[ee];
    		DiskSNIaAllElements_ts[ee] += SNIaAllElements[ee];
    		DiskAGBAllElements_ts[ee] += AGBAllElements[ee];
    	}
#endif //INDIVIDUAL_ELEMENTS
#endif //H2_AND_RINGS
#endif //DETAILED_DUST

    	mass_checks(p,"model_yields.c",__LINE__);

    	//************************
    	// UPDATE GAS COMPONENTS:
    	//************************

    	//Hot Gas:
     	Gal[igal].HotGas += (fwind_SNII * SNIIEjectaMass) + (fwind_SNIa * SNIaEjectaMass) + (fwind_AGB * AGBEjectaMass);
     	Gal[igal].MetalsHotGas[0] += fwind_SNII * SNIIAllMetals;
     	Gal[igal].MetalsHotGas[1] += fwind_SNIa * SNIaAllMetals;
     	Gal[igal].MetalsHotGas[2] += fwind_AGB * AGBAllMetals;
     	mass_checks(p,"model_yields.c",__LINE__);
     	//Cold Gas:
     	Gal[p].ColdGas += (1.0-fwind_SNII)*SNIIEjectaMass + (1.0-fwind_SNIa)*SNIaEjectaMass + (1.0-fwind_AGB)*AGBEjectaMass;
     	Gal[p].MetalsColdGas[0] += (1.0-fwind_SNII) * SNIIAllMetals;
     	Gal[p].MetalsColdGas[1] += (1.0-fwind_SNIa) * SNIaAllMetals;
     	Gal[p].MetalsColdGas[2] += (1.0-fwind_AGB) * AGBAllMetals;

#ifdef H2_AND_RINGS
     	Gal[p].ColdGasRings[jj] += ((1.0-fwind_SNII) * SNIIEjectaMass + (1.0-fwind_SNIa) * SNIaEjectaMass + (1.0-fwind_AGB) * AGBEjectaMass);
     	Gal[p].MetalsColdGasRings[jj][0] += (1.0-fwind_SNII) * SNIIAllMetals;
     	Gal[p].MetalsColdGasRings[jj][1] += (1.0-fwind_SNIa) * SNIaAllMetals;
     	Gal[p].MetalsColdGasRings[jj][2] += (1.0-fwind_AGB) * AGBAllMetals;
#endif //H2_AND_RINGS

    	//Total:
     	TotalMassReturnedToHotGas += (fwind_SNII * SNIIEjectaMass) + (fwind_SNIa * SNIaEjectaMass) + (fwind_AGB * AGBEjectaMass);
     	TotalMetalsReturnedToHotGas += (fwind_SNII * SNIIAllMetals) + (fwind_SNIa * SNIaAllMetals) + (fwind_AGB * AGBAllMetals);
     	TotalMassReturnedToColdDiskGas += (1.0-fwind_SNII) * SNIIEjectaMass + (1.0-fwind_SNIa) * SNIaEjectaMass + (1.0-fwind_AGB) * AGBEjectaMass; //Only use energy from SNe that eject into ColdGas to reheat
#ifdef H2_AND_RINGS
     	TotalMassReturnedToColdDiskGasr[jj] += ((1.0-fwind_SNII) * SNIIEjectaMass + (1.0-fwind_SNIa) * SNIaEjectaMass + (1.0-fwind_AGB) * AGBEjectaMass);
#endif //H2_AND_RINGS

#ifdef INDIVIDUAL_ELEMENTS
     	for(ee=0;ee<NUM_ELEMENTS;ee++)
     	  {
     	    Gal[igal].HotGas_elements[ee] += fwind_SNII * SNIIAllElements[ee] + fwind_SNIa * SNIaAllElements[ee] + fwind_AGB * AGBAllElements[ee];
     		Gal[p].ColdGas_elements[ee] += (1.0-fwind_SNII) * SNIIAllElements[ee] + (1.0-fwind_SNIa) * SNIaAllElements[ee] + (1.0-fwind_AGB) * AGBAllElements[ee];
#ifdef H2_AND_RINGS
     	    Gal[p].ColdGasRings_elements[jj][ee] += (1.0-fwind_SNII) * SNIIAllElements[ee] + (1.0-fwind_SNIa) * SNIaAllElements[ee] + (1.0-fwind_AGB) * AGBAllElements[ee];
#endif // H2_AND_RINGS
     	  }

#ifdef DETAILED_DUST
#ifdef H2_AND_RINGS
#ifndef MAINELEMENTS
     		SNII_prevstep_Cold_Cb[jj][ii] += SNIIAllElements[Cb_NUM]; //Element mass from Disc SNe-II that is ejected in TOTAL from each SFH bin (ii). Used to calculate SN-II dust production mass and rate from Disc SNe-II in model_dust_yields.c.
     		SNII_prevstep_Cold_Si[jj][ii] += SNIIAllElements[Si_NUM];
#endif //MAINELEMENTS
     		SNII_prevstep_Cold_Fe[jj][ii] += SNIIAllElements[Fe_NUM];
     		SNIa_prevstep_Cold_Fe[jj][ii] += SNIaAllElements[Fe_NUM];
#else //H2_AND_RINGS
#ifndef MAINELEMENTS
     		SNII_prevstep_Cold_Cb[ii] += SNIIAllElements[Cb_NUM];
     		SNII_prevstep_Cold_Si[ii] += SNIIAllElements[Si_NUM];
#endif //MAINELEMENTS
     		SNII_prevstep_Cold_Fe[ii] += SNIIAllElements[Fe_NUM];
     		SNIa_prevstep_Cold_Fe[ii] += SNIaAllElements[Fe_NUM];
#endif //H2_AND_RINGS
#endif //DETAILED_DUST
#endif //INDIVIDUAL_ELEMENTS

     	mass_checks(p,"model_yields.c",__LINE__);


     	//*****************************
     	// UPDATE DISK MASS COMPONENTS:
     	//*****************************
     	Gal[p].DiskMass -= SNIIEjectaMass + SNIaEjectaMass + AGBEjectaMass;
     	Gal[p].MetalsDiskMass[0] -= SNIIUnProcessedMetals;
     	Gal[p].MetalsDiskMass[1] -= SNIaUnProcessedMetals;
     	Gal[p].MetalsDiskMass[2] -= AGBUnProcessedMetals;
     	Gal[p].DiskSNIIRate += SNIIRate/STEPS; //in 1/(time code units). This records the average rate in each snapshot [i.e. SUM(SNIIRate per timestep)/STEPS], so is set to 0.0 at the start of each snapshot in model_misc.c or main.c (I forget which!). (22-05-20)
     	Gal[p].DiskSNIaRate += SNIaRate/STEPS; //in 1/(time code units).

#ifdef H2_AND_RINGS
     	Gal[p].DiskMassRings[jj] -= SNIIEjectaMass + SNIaEjectaMass + AGBEjectaMass;
     	Gal[p].MetalsDiskMassRings[jj][0] -= SNIIUnProcessedMetals;
     	Gal[p].MetalsDiskMassRings[jj][1] -= SNIaUnProcessedMetals;
     	Gal[p].MetalsDiskMassRings[jj][2] -= AGBUnProcessedMetals;
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
    			BulgeMetallicity = 0.;
    			for(mm=0;mm<NUM_METAL_CHANNELS;mm++)
    				BulgeMetallicity += Gal[p].sfh_MetalsBulgeMassRings[jj][ii][mm];
    			BulgeMetallicity /= Gal[p].sfh_BulgeMassRings[jj][ii];

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
    	Zi = find_initial_metallicity(BulgeMetallicity, 5); //Find the metallicity bin BELOW the true metallicity in the SFH bin
    	Zi_disp = (BulgeMetallicity - bcMetallicities[Zi])/(bcMetallicities[Zi+1] - bcMetallicities[Zi]); //Fractional position of the SFH metallicity between the bcMetallicities bins.
#else //BINARYC
    	Zi = find_initial_metallicity(BulgeMetallicity, 1); //Find the metallicity bin BELOW the true metallicity in the SFH bin
    	Zi_disp = (BulgeMetallicity - lifetimeMetallicities[Zi])/(lifetimeMetallicities[Zi+1] - lifetimeMetallicities[Zi]); //Fractional position of the SFH metallicity between the bcMetallicities bins.
#endif //BINARYC
    	if (Zi_disp < 0.0) Zi_disp = 0.0; //Don't want to extrapolate yields down below lifetimeMetallicities[0]=0.0004. Instead, assume constant yield below this metallicity.

#ifdef DETAILED_DUST //For use in update_dust_mass() later on
#ifdef H2_AND_RINGS
    	Zi_bulge_saved[jj][ii] = Zi;
		Zi_bulge_disp_saved[jj][ii] = Zi_disp;
#else
		Zi_bulge_saved[ii] = Zi;
		Zi_bulge_disp_saved[ii] = Zi_disp;
#endif //H2_AND_RINGS
#endif

#ifndef INDIVIDUAL_ELEMENTS
    	compute_actual_eject_rates(TimeBin, ii, Zi, Zi_disp, Gal[p].sfh_ibin, BulgeSFRxStep, BulgeSFRxStep_Phys, BulgeMetallicity, timestep_width,
								   &SNIIEjectaMass, &SNIIAllMetals, &SNIIUnProcessedMetals,
								   &SNIaEjectaMass, &SNIaAllMetals, &SNIaUnProcessedMetals,
								   &AGBEjectaMass, &AGBAllMetals, &AGBUnProcessedMetals,
								   &SNIIRate, &SNIaRate, &AGBRate);
#ifdef H2_AND_RINGS
    	rescale_ejected_material(p, jj, "BulgeMass", SNIIEjectaMass, SNIaEjectaMass, AGBEjectaMass,
								 &SNIIUnProcessedMetals, &SNIaUnProcessedMetals, &AGBUnProcessedMetals,
								 &SNIIAllMetals, &SNIaAllMetals, &AGBAllMetals);
#else //H2_AND_RINGS
		rescale_ejected_material(p, "BulgeMass", SNIIEjectaMass, SNIaEjectaMass, AGBEjectaMass,
								 &SNIIUnProcessedMetals, &SNIaUnProcessedMetals, &AGBUnProcessedMetals,
								 &SNIIAllMetals, &SNIaAllMetals, &AGBAllMetals);
#endif //H2_AND_RINGS
#else //INDIVIDUAL_ELEMENTS
    	compute_actual_eject_rates(TimeBin, ii, Zi, Zi_disp, Gal[p].sfh_ibin, BulgeSFRxStep, BulgeSFRxStep_Phys, BulgeMetallicity, BulgeMetallicityElement_Phys, timestep_width,
								   &SNIIEjectaMass, &SNIIAllMetals, &SNIIUnProcessedMetals,
								   &SNIaEjectaMass, &SNIaAllMetals, &SNIaUnProcessedMetals,
								   &AGBEjectaMass, &AGBAllMetals, &AGBUnProcessedMetals,
								   &SNIIRate, &SNIaRate, &AGBRate,
								   SNIIAllElements, SNIIUnProcessedElements, SNIaAllElements, SNIaUnProcessedElements, AGBAllElements, AGBUnProcessedElements);
#ifdef H2_AND_RINGS
    	rescale_ejected_material(p, jj, "BulgeMass", SNIIEjectaMass, SNIaEjectaMass, AGBEjectaMass,
								 &SNIIUnProcessedMetals, &SNIaUnProcessedMetals, &AGBUnProcessedMetals,
								 &SNIIAllMetals, &SNIaAllMetals, &AGBAllMetals,
								 SNIIAllElements, SNIaAllElements, AGBAllElements,
								 SNIIUnProcessedElements, SNIaUnProcessedElements, AGBUnProcessedElements);
#else //H2_AND_RINGS
		rescale_ejected_material(p, "BulgeMass", SNIIEjectaMass, SNIaEjectaMass, AGBEjectaMass,
								 &SNIIUnProcessedMetals, &SNIaUnProcessedMetals, &AGBUnProcessedMetals,
								 &SNIIAllMetals, &SNIaAllMetals, &AGBAllMetals,
								 SNIIAllElements, SNIaAllElements, AGBAllElements,
								 SNIIUnProcessedElements, SNIaUnProcessedElements, AGBUnProcessedElements);
#endif //H2_AND_RINGS
#endif //INDIVIDUAL_ELEMENTS

//Store for use in dust_yields.c:
#ifdef DETAILED_DUST
#ifdef DUST_HOTGAS
#ifdef H2_AND_RINGS
		//We have not included dust destruction from SN shocks in the bulge yet, so these variables are not yet needed: (01-03-23):
		//BulgeSNIIRate_current_ts[jj] += SNIIRate;
		//BulgeSNIaRate_current_ts[jj] += SNIaRate;
		//BulgeAGBRate_current_ts[jj] += AGBRate; //Not currently used, but calculated anyway. (29-11-21)
#ifdef INDIVIDUAL_ELEMENTS
		for(ee=0;ee<NUM_ELEMENTS;ee++) {
			BulgeSNIIAllElements_ts[jj][ee] += SNIIAllElements[ee];
			BulgeSNIaAllElements_ts[jj][ee] += SNIaAllElements[ee];
			BulgeAGBAllElements_ts[jj][ee] += AGBAllElements[ee];
		}
#endif //INDIVIDUAL_ELEMENTS
#else //H2_AND_RINGS
		//BulgeSNIIRate_current_ts += SNIIRate;
		//BulgeSNIaRate_current_ts += SNIaRate;
		//BulgeAGBRate_current_ts += AGBRate; //Not currently used, but calculated anyway. (29-11-21)
#ifdef INDIVIDUAL_ELEMENTS
		for(ee=0;ee<NUM_ELEMENTS;ee++) {
			BulgeSNIIAllElements_ts[ee] += SNIIAllElements[ee];
			BulgeSNIaAllElements_ts[ee] += SNIaAllElements[ee];
			BulgeAGBAllElements_ts[ee] += AGBAllElements[ee];
		}
#endif //INDIVIDUAL_ELEMENTS
#endif //H2_AND_RINGS
#endif //DUST_HOTGAS
#endif //DETAILED_DUST

		mass_checks(p,"model_yields.c",__LINE__);


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
    	//TotalMetalsReturnedToHotGas += SNIIAllMetals + SNIaAllMetals + AGBAllMetals; //NOTE: Only metals returned to the HotGas from disc stars & SNe are considered in TotalMetalsReturnedToHotGas (i.e. not from the bulge or ICL):
    	Gal[igal].MetalsHotGas[0] += SNIIAllMetals;
    	Gal[igal].MetalsHotGas[1] += SNIaAllMetals;
    	Gal[igal].MetalsHotGas[2] += AGBAllMetals;

#ifdef INDIVIDUAL_ELEMENTS

    	for(ee=0;ee<NUM_ELEMENTS;ee++) {
    	  Gal[igal].HotGas_elements[ee] += SNIIAllElements[ee]+SNIaAllElements[ee]+AGBAllElements[ee];
    	}

#ifdef DETAILED_DUST
#ifdef H2_AND_RINGS
#ifndef MAINELEMENTS
     		SNII_prevstep_Hot_bulge_Cb[jj][ii] += SNIIAllElements[Cb_NUM];
     		SNII_prevstep_Hot_bulge_Si[jj][ii] += SNIIAllElements[Si_NUM];
#endif //MAINELEMENTS
     		SNII_prevstep_Hot_bulge_Fe[jj][ii] += SNIIAllElements[Fe_NUM];
     		SNIa_prevstep_Hot_bulge_Fe[jj][ii] += SNIaAllElements[Fe_NUM];
#else //H2_AND_RINGS
#ifndef MAINELEMENTS
     		SNII_prevstep_Hot_bulge_Cb[ii] += SNIIAllElements[Cb_NUM];
     		SNII_prevstep_Hot_bulge_Si[ii] += SNIIAllElements[Si_NUM];
#endif //MAINELEMENTS
     		SNII_prevstep_Hot_bulge_Fe[ii] += SNIIAllElements[Fe_NUM];
     		SNIa_prevstep_Hot_bulge_Fe[ii] += SNIaAllElements[Fe_NUM];
#endif //H2_AND_RINGS
#endif //DETAILED_DUST

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
    	Gal[p].BulgeSNIaRate += SNIaRate/STEPS; //in 1/(time code units).

#ifdef H2_AND_RINGS
    	Gal[p].BulgeMassRings[jj] -= SNIIEjectaMass + SNIaEjectaMass + AGBEjectaMass;
    	Gal[p].MetalsBulgeMassRings[jj][0] -= SNIIUnProcessedMetals;
    	Gal[p].MetalsBulgeMassRings[jj][1] -= SNIaUnProcessedMetals;
    	Gal[p].MetalsBulgeMassRings[jj][2] -= AGBUnProcessedMetals;
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

#ifdef INDIVIDUAL_ELEMENTS
    	for (ee=0;ee<NUM_ELEMENTS;ee++)
    	  ICMMetallicityElement_Phys[ee] = Gal[p].sfh_ICM_elements[ii][ee] / (Gal[p].sfh_ICM[ii]*1.0e10/Hubble_h);
#endif

#ifdef BINARYC
    	Zi = find_initial_metallicity(ICMMetallicity, 5); //Find the metallicity bin BELOW the true metallicity in the SFH bin
    	Zi_disp = (ICMMetallicity - bcMetallicities[Zi])/(bcMetallicities[Zi+1] - bcMetallicities[Zi]); //Fractional position of the SFH metallicity between the bcMetallicities bins.
#else //BINARYC
    	Zi = find_initial_metallicity(ICMMetallicity, 1); //Find the metallicity bin BELOW the true metallicity in the SFH bin
    	Zi_disp = (ICMMetallicity - lifetimeMetallicities[Zi])/(lifetimeMetallicities[Zi+1] - lifetimeMetallicities[Zi]); //Fractional position of the SFH metallicity between the bcMetallicities bins.
#endif //BINARYC
    	if (Zi_disp < 0.0) Zi_disp = 0.0; //Don't want to extrapolate yields down below lifetimeMetallicities[0]=0.0004. Instead, assume constant yield below this metallicity.

#ifdef DETAILED_DUST //For use in update_dust_mass() later on
#ifdef H2_AND_RINGS
    	Zi_ICM_saved[jj][ii] = Zi;
		Zi_ICM_disp_saved[jj][ii] = Zi_disp;
#else
		Zi_ICM_saved[ii] = Zi;
		Zi_ICM_disp_saved[ii] = Zi_disp;
#endif //H2_AND_RINGS
#endif

#ifndef INDIVIDUAL_ELEMENTS
    	compute_actual_eject_rates(TimeBin, ii, Zi, Zi_disp, Gal[p].sfh_ibin, ICMSFRxStep, ICMSFRxStep_Phys, ICMMetallicity, timestep_width,
				   &SNIIEjectaMass, &SNIIAllMetals, &SNIIUnProcessedMetals,
				   &SNIaEjectaMass, &SNIaAllMetals, &SNIaUnProcessedMetals,
				   &AGBEjectaMass, &AGBAllMetals, &AGBUnProcessedMetals,
				   &SNIIRate, &SNIaRate, &AGBRate);

#ifdef H2_AND_RINGS
    	rescale_ejected_material(p, jj, "ICM", SNIIEjectaMass, SNIaEjectaMass, AGBEjectaMass,
    									 &SNIIUnProcessedMetals, &SNIaUnProcessedMetals, &AGBUnProcessedMetals,
    									 &SNIIAllMetals, &SNIaAllMetals, &AGBAllMetals);
#else //H2_AND_RINGS
    	rescale_ejected_material(p, "ICM", SNIIEjectaMass, SNIaEjectaMass, AGBEjectaMass,
										 &SNIIUnProcessedMetals, &SNIaUnProcessedMetals, &AGBUnProcessedMetals,
										 &SNIIAllMetals, &SNIaAllMetals, &AGBAllMetals);
#endif //H2_AND_RINGS

#else //INDIVIDUAL_ELEMENTS
    	compute_actual_eject_rates(TimeBin, ii, Zi, Zi_disp, Gal[p].sfh_ibin, ICMSFRxStep, ICMSFRxStep_Phys, ICMMetallicity, ICMMetallicityElement_Phys, timestep_width,
				   &SNIIEjectaMass, &SNIIAllMetals, &SNIIUnProcessedMetals,
				   &SNIaEjectaMass, &SNIaAllMetals, &SNIaUnProcessedMetals,
				   &AGBEjectaMass, &AGBAllMetals, &AGBUnProcessedMetals,
				   &SNIIRate, &SNIaRate, &AGBRate,
				   SNIIAllElements, SNIIUnProcessedElements, SNIaAllElements, SNIaUnProcessedElements, AGBAllElements, AGBUnProcessedElements);

#ifdef H2_AND_RINGS
    	rescale_ejected_material(p, jj, "ICM", SNIIEjectaMass, SNIaEjectaMass, AGBEjectaMass,
								 &SNIIUnProcessedMetals, &SNIaUnProcessedMetals, &AGBUnProcessedMetals,
								 &SNIIAllMetals, &SNIaAllMetals, &AGBAllMetals,
								 SNIIAllElements, SNIaAllElements, AGBAllElements,
								 SNIIUnProcessedElements, SNIaUnProcessedElements, AGBUnProcessedElements);
#else //H2_AND_RINGS
		rescale_ejected_material(p, "ICM", SNIIEjectaMass, SNIaEjectaMass, AGBEjectaMass,
								 &SNIIUnProcessedMetals, &SNIaUnProcessedMetals, &AGBUnProcessedMetals,
								 &SNIIAllMetals, &SNIaAllMetals, &AGBAllMetals,
								 SNIIAllElements, SNIaAllElements, AGBAllElements,
								 SNIIUnProcessedElements, SNIaUnProcessedElements, AGBUnProcessedElements);
#endif //H2_AND_RINGS
#endif //INDIVIDUAL_ELEMENTS

//Store for use in dust_yields.c:
#ifdef DETAILED_DUST
#ifdef DUST_HOTGAS
#ifdef H2_AND_RINGS
		//We have not included dust destruction from SN shocks in the stellar halo yet, so these variables are not yet needed: (01-03-23):
		//ICMSNIIRate_current_ts[jj] += SNIIRate;
		//ICMSNIaRate_current_ts[jj] += SNIaRate;
		//ICMAGBRate_current_ts[jj] += AGBRate; //Not currently used, but calculated anyway. (29-11-21)
#ifdef INDIVIDUAL_ELEMENTS
		for(ee=0;ee<NUM_ELEMENTS;ee++) {
			ICMSNIIAllElements_ts[jj][ee] += SNIIAllElements[ee];
			ICMSNIaAllElements_ts[jj][ee] += SNIaAllElements[ee];
			ICMAGBAllElements_ts[jj][ee] += AGBAllElements[ee];
		}
#endif //INDIVIDUAL_ELEMENTS
#else //H2_AND_RINGS
		//ICMSNIIRate_current_ts += SNIIRate;
		//ICMSNIaRate_current_ts += SNIaRate;
		//ICMAGBRate_current_ts += AGBRate; //Not currently used, but calculated anyway. (29-11-21)
#ifdef INDIVIDUAL_ELEMENTS
		for(ee=0;ee<NUM_ELEMENTS;ee++) {
			ICMSNIIAllElements_ts[ee] += SNIIAllElements[ee];
			ICMSNIaAllElements_ts[ee] += SNIaAllElements[ee];
			ICMAGBAllElements_ts[ee] += AGBAllElements[ee];
		}
#endif //INDIVIDUAL_ELEMENTS
#endif //H2_AND_RINGS
#endif //DUST_HOTGAS
#endif //DETAILED_DUST

		mass_checks(p,"model_yields.c",__LINE__);


    	//*****************************
    	// UPDATE HOT GAS COMPONENTS:
    	//*****************************
    	Gal[igal].HotGas += SNIIEjectaMass + SNIaEjectaMass + AGBEjectaMass;
    	//If there was no hotgas left in the galaxy it will probably be stripped next step.
    	//Give it a fake HotRadius for now to avoid crash at mass checks
    	if(Gal[igal].HotGas>0. && Gal[igal].HotRadius==0)
    	  Gal[igal].HotRadius=1.e-10;

    	TotalMassReturnedToHotGas += SNIIEjectaMass + SNIaEjectaMass + AGBEjectaMass;
    	//TotalMetalsReturnedToHotGas += SNIIAllMetals + SNIaAllMetals + AGBAllMetals; //NOTE: Only metals returned to the HotGas from disc stars & SNe are considered in TotalMetalsReturnedToHotGas (i.e. not from the bulge or ICL):
    	Gal[igal].MetalsHotGas[0] += SNIIAllMetals;
    	Gal[igal].MetalsHotGas[1] += SNIaAllMetals;
    	Gal[igal].MetalsHotGas[2] += AGBAllMetals;

#ifdef INDIVIDUAL_ELEMENTS
    	for(ee=0;ee<NUM_ELEMENTS;ee++)
    	  Gal[igal].HotGas_elements[ee] += SNIIAllElements[ee] + SNIaAllElements[ee] + AGBAllElements[ee];

#ifdef DETAILED_DUST
#ifdef H2_AND_RINGS
#ifndef MAINELEMENTS
     		SNII_prevstep_Hot_ICM_Cb[jj][ii] += SNIIAllElements[Cb_NUM];
     		SNII_prevstep_Hot_ICM_Si[jj][ii] += SNIIAllElements[Si_NUM];
#endif //MAINELEMENTS
     		SNII_prevstep_Hot_ICM_Fe[jj][ii] += SNIIAllElements[Fe_NUM];
     		SNIa_prevstep_Hot_ICM_Fe[jj][ii] += SNIaAllElements[Fe_NUM];
#else //H2_AND_RINGS
#ifndef MAINELEMENTS
     		SNII_prevstep_Hot_ICM_Cb[ii] += SNIIAllElements[Cb_NUM];
     		SNII_prevstep_Hot_ICM_Si[ii] += SNIIAllElements[Si_NUM];
#endif //MAINELEMENTS
     		SNII_prevstep_Hot_ICM_Fe[ii] += SNIIAllElements[Fe_NUM];
     		SNIa_prevstep_Hot_ICM_Fe[ii] += SNIaAllElements[Fe_NUM];
#endif //H2_AND_RINGS
#endif //DETAILED_DUST
#endif //INDIVIDUAL_ELEMENTS


    	//*****************************
    	// UPDATE ICM COMPONENTS:
    	//*****************************
    	Gal[p].ICM -= SNIIEjectaMass + SNIaEjectaMass + AGBEjectaMass;
    	Gal[p].MetalsICM[0] -= SNIIUnProcessedMetals;
    	Gal[p].MetalsICM[1] -= SNIaUnProcessedMetals;
    	Gal[p].MetalsICM[2] -= AGBUnProcessedMetals;
    	Gal[p].ICMSNIIRate += SNIIRate/STEPS; //in 1/(time code units).
    	Gal[p].ICMSNIaRate += SNIaRate/STEPS; //in 1/(time code units).

#ifdef INDIVIDUAL_ELEMENTS
    	for(ee=0;ee<NUM_ELEMENTS;ee++)
    	  Gal[p].ICM_elements[ee] -= SNIIUnProcessedElements[ee] + SNIaUnProcessedElements[ee] + AGBUnProcessedElements[ee];
#endif //INDIVIDUAL_ELEMENTS

    	//Update ages: //NOTE: Galaxy mass-weighted ages do not take account of halo star ages currently, only those in the disc and bulge. (11-10-23)
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
	// SN feedback must be called AFTER BOTH model_yields AND model_dust_yields
	// and now appears inside main.c if this is the case.
#ifdef FEEDBACK_COUPLED_WITH_MASS_RETURN

    if(TotalMassReturnedToColdDiskGas>0.) {
    	if(TotalMassReturnedToColdDiskGas>Gal[p].ColdGas)
    		TotalMassReturnedToColdDiskGas=Gal[p].ColdGas;
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
		SN_feedback(p, centralgal, TotalMassReturnedToHotGas, HotGasRings, "HotGas", dt);
#endif
		Gal[p].MassReturnRateToHotGas += TotalMassReturnedToHotGas / (dt * STEPS); //ROB: Store mass return rate by SNe and stellar winds to the gas phases (from disk stars, bulge stars, and ICL stars).
      }

    if(TotalMetalsReturnedToHotGas>0.) {
    	//NOTE: Only metals returned to the HotGas from disc stars & SNe are considered in TotalMetalsReturnedToHotGas (i.e. not from the bulge or ICL):
    	Gal[p].MetalsReturnRateToHotGas += TotalMetalsReturnedToHotGas / (dt * STEPS); //ROB: Store mass return rate by disc SNe and stellar winds to the gas phases.
      }

#endif //FEEDBACK_COUPLED_WITH_MASS_RETURN
#endif //DETAILED_DUST

      mass_checks(p,"model_yields.c",__LINE__);

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
//This function interpolates between the normalised ejecta masses from the discreet metallicity bins generated in yields_integrals.c,
//to find the true normalised ejecta masses at the actual metallicity of the SFH bin.
//Analogous to calc_ejecta_limits() in yields_integrals.c, but interpolating across metallicity bins, rather than mass or time bins. (10-05-22)
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
  *SNIIEjectaMass = SFRxStep * NormSNIIMassEjecRate_actual; //N.B. the masses in the P98 yield tables are *total* ejected masses (newly synth + unproc), unlike the metals and elements which are only newly synth. (04-07-22)
#ifdef INSTANTANEOUS_RECYCLE
  *SNIIAllMetals = SFRxStep * (NormSNIIMetalEjecRate_actual + Metallicity * NormSNIIMassEjecRate_actual);
#else //INSTANTANEOUS_RECYCLE
#if defined(CHIEFFI) || defined(BINARYC)
  //ROB: No unsynth component required for SN-II ejecta, when using the Chieffi & Limongi (2007) or binary_c yield tables
  *SNIIAllMetals = SFRxStep * NormSNIIMetalEjecRate_actual;
#elif defined(PORTINARI)
  *SNIIAllMetals = SFRxStep * (NormSNIIMetalEjecRate_actual + Metallicity * NormSNIIMassEjecRate_actual);
#endif
#endif //INSTANTANEOUS_RECYCLE
  *SNIIUnProcessedMetals = SFRxStep * Metallicity * NormSNIIMassEjecRate_actual;
  *SNIIRate = (SFRxStep_Phys/timestep_width) * NormSNIINum_actual; //This is the SNII rate in [1/(time code units)] in this timestep. So should divide by (UnitTime_in_s / SEC_PER_YEAR) in save.c to get "per year" before outputting).

  //SNIa
  *SNIaEjectaMass = SFRxStep * NormSNIaMassEjecRate_actual;
  //SNIa yields are written in a different way to SNIIa and AGB, so the following line is correct:
  *SNIaAllMetals = SFRxStep * NormSNIaMetalEjecRate_actual;
  *SNIaUnProcessedMetals = SFRxStep * Metallicity * NormSNIaMassEjecRate_actual;
  *SNIaRate = (SFRxStep_Phys/timestep_width) * NormSNIaNum_actual; //This is the SNIa rate in [1/(time code units)] in this timestep. So should divide by (UnitTime_in_s / SEC_PER_YEAR) in save.c to get "per year" before outputting).

  //AGB
  *AGBEjectaMass = SFRxStep * NormAGBMassEjecRate_actual;
#ifdef BINARYC
  //ROB: No unsynth component required for AGB ejecta, when using the binary_c yield tables
  *AGBAllMetals = SFRxStep * NormAGBMetalEjecRate_actual;
#else //BINARYC
  *AGBAllMetals = SFRxStep * (NormAGBMetalEjecRate_actual + (Metallicity * NormAGBMassEjecRate_actual));
#endif //BINARYC
  *AGBUnProcessedMetals = SFRxStep * Metallicity * NormAGBMassEjecRate_actual;
  *AGBRate = (SFRxStep_Phys/timestep_width) * NormAGBNum_actual; //This is the AGB rate in [1/(time code units)] in this timestep. So should divide by (UnitTime_in_s / SEC_PER_YEAR) in save.c to get "per year" before outputting).

#ifdef INDIVIDUAL_ELEMENTS
  for(ee=0;ee<NUM_ELEMENTS;ee++) {
#if defined(CHIEFFI) || defined(BINARYC)
      SNIIAllElements[ee] = SFRxStep_Phys * NormSNIIYieldRate_actual[ee];
#elif defined(PORTINARI)
      SNIIAllElements[ee] = SFRxStep_Phys * (NormSNIIYieldRate_actual[ee] + MetallicityElement_Phys[ee] * NormSNIIMassEjecRate_actual);
#endif
      SNIIUnProcessedElements[ee] = SFRxStep_Phys * MetallicityElement_Phys[ee] * NormSNIIMassEjecRate_actual;
      SNIaAllElements[ee] = SFRxStep_Phys * NormSNIaYieldRate_actual[ee];
      SNIaUnProcessedElements[ee] = SFRxStep_Phys * MetallicityElement_Phys[ee] * NormSNIaMassEjecRate_actual;

#ifdef BINARYC
      AGBAllElements[ee] = SFRxStep_Phys * NormAGBYieldRate_actual[ee];
#else //BINARYC
      AGBAllElements[ee] = SFRxStep_Phys * (NormAGBYieldRate_actual[ee] + MetallicityElement_Phys[ee] * NormAGBMassEjecRate_actual);

#endif //BINARYC
      AGBUnProcessedElements[ee] = SFRxStep_Phys * MetallicityElement_Phys[ee] * NormAGBMassEjecRate_actual;
  }
#endif //INDIVIDUAL_ELEMENTS
}








int find_initial_metallicity(double metallicity, int table_type) //, int component_type
{
	int i, Zi_bin;
	Zi_bin = -1;
	i = 0;

	switch (table_type)
	{
#ifndef BINARYC
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




#ifdef INDIVIDUAL_ELEMENTS
#ifdef H2_AND_RINGS
void rescale_ejected_material(int p, int jj, char component[], double SNIIEjectaMass, double SNIaEjectaMass, double AGBEjectaMass,
							  double* SNIIUnProcessedMetals, double* SNIaUnProcessedMetals, double* AGBUnProcessedMetals,
							  double* SNIIAllMetals, double* SNIaAllMetals, double* AGBAllMetals,
							  double* SNIIAllElements, double* SNIaAllElements, double* AGBAllElements,
							  double* SNIIUnProcessedElements, double* SNIaUnProcessedElements, double* AGBUnProcessedElements) {
#else //H2_AND_RINGS
void rescale_ejected_material(int p, char component[], double SNIIEjectaMass, double SNIaEjectaMass, double AGBEjectaMass,
							  double* SNIIUnProcessedMetals, double* SNIaUnProcessedMetals, double* AGBUnProcessedMetals,
							  double* SNIIAllMetals, double* SNIaAllMetals, double* AGBAllMetals,
							  double* SNIIAllElements, double* SNIaAllElements, double* AGBAllElements,
							  double* SNIIUnProcessedElements, double* SNIaUnProcessedElements, double* AGBUnProcessedElements) {
#endif //H2_AND_RINGS
#else //INDIVIDUAL_ELEMENTS
#ifdef H2_AND_RINGS
void rescale_ejected_material(int p, int jj, char component[],
							  double* SNIIUnProcessedMetals, double* SNIaUnProcessedMetals, double* AGBUnProcessedMetals,
							  double* SNIIAllMetals, double* SNIaAllMetals, double* AGBAllMetals) {
#else //H2_AND_RINGS
void rescale_ejected_material(int p, char component[],
							  double* SNIIUnProcessedMetals, double* SNIaUnProcessedMetals, double* AGBUnProcessedMetals,
							  double* SNIIAllMetals, double* SNIaAllMetals, double* AGBAllMetals) {
#endif //H2_AND_RINGS
#endif //INDIVIDUAL_ELEMENTS
/*
 * ROB: (07-07-22): New way to scale unprocessed components:
 * First, correct element ejected masses so that total elements ejected matches total mass ejected
 * Second, correct metal ejected & unprocessed masses so that total metals ejected matches total metal elements ejected
 * Third, redistribute unprocessed masses from AGB and SNe-Ia to SNe-II, if necessary
*/
	//#####
	//0) Assign metal sub-components (SNe-II, SNe-Ia, or AGBs) for the correct mass component (DiskMass, BulgeMass, or ICM):
	double SNIImetals, SNIametals, AGBmetals;
#ifdef H2_AND_RINGS
	if (strcmp(component,"DiskMass") == 0) {
		SNIImetals = Gal[p].MetalsDiskMassRings[jj][0];
		SNIametals = Gal[p].MetalsDiskMassRings[jj][1];
		AGBmetals = Gal[p].MetalsDiskMassRings[jj][2];
	}
	else if (strcmp(component,"BulgeMass") == 0) {
		SNIImetals = Gal[p].MetalsBulgeMassRings[jj][0];
		SNIametals = Gal[p].MetalsBulgeMassRings[jj][1];
		AGBmetals = Gal[p].MetalsBulgeMassRings[jj][2];
	}
	else if (strcmp(component,"ICM") == 0) {
		SNIImetals = Gal[p].MetalsICM[0];
		SNIametals = Gal[p].MetalsICM[1];
		AGBmetals = Gal[p].MetalsICM[2];
	}
	else  {
		printf("\nERROR: rescale_ejected_material(): Component %c is invalid. Choose from DiskMass, BulgeMass, or ICM.\n", component);
		exit(1);
	}
#else //H2_AND_RINGS
	if (strcmp(component,"DiskMass") == 0) {
		SNIImetals = Gal[p].MetalsDiskMass[0];
		SNIametals = Gal[p].MetalsDiskMass[1];
		AGBmetals = Gal[p].MetalsDiskMass[2];
	}
	else if (strcmp(component,"BulgeMass") == 0) {
		SNIImetals = Gal[p].MetalsBulgeMass[0];
		SNIametals = Gal[p].MetalsBulgeMass[1];
		AGBmetals = Gal[p].MetalsBulgeMass[2];
	}
	else if (strcmp(component,"ICM") == 0) {
		SNIImetals = Gal[p].MetalsICM[0];
		SNIametals = Gal[p].MetalsICM[1];
		AGBmetals = Gal[p].MetalsICM[2];
	}
	else  {
		printf("\nERROR: rescale_ejected_material(): Component %c is invalid. Choose from DiskMass, BulgeMass, or ICM.\n", component);
		exit(1);
	}
#endif //H2_AND_RINGS

#ifdef INDIVIDUAL_ELEMENTS
	int ee;
	double TotSNIIElementsEjected = 0.0, TotSNIaElementsEjected = 0.0, TotAGBElementsEjected = 0.0;
	double SNIIElementCorrectionFactor = 1.0, SNIaElementCorrectionFactor = 1.0, AGBElementCorrectionFactor = 1.0;
	double TotMetalSNIIElementsEjected = 0.0, TotMetalSNIaElementsEjected = 0.0, TotMetalAGBElementsEjected = 0.0;
	double SNIIMetalCorrectionFactor = 1.0, SNIaMetalCorrectionFactor = 1.0, AGBMetalCorrectionFactor = 1.0;
	double TotMetalSNIIUnProcessedElements = 0.0, TotMetalSNIaUnProcessedElements = 0.0, TotMetalAGBUnProcessedElements = 0.0;
	double SNIIMetalExcessFactor, SNIaMetalExcessFactor1, SNIaMetalExcessFactor2, AGBMetalExcessFactor1, AGBMetalExcessFactor2;
	double SNIaElementExcessFactor[NUM_ELEMENTS], AGBElementExcessFactor[NUM_ELEMENTS];
	for(ee=0;ee<NUM_ELEMENTS;ee++) {
		SNIaElementExcessFactor[ee] = 0.0;
		AGBElementExcessFactor[ee] = 0.0;
	}

	//#####
	//1) Correct element unprocessed masses so that total elements ejected matches total mass ejected:
	for(ee=0;ee<NUM_ELEMENTS;ee++) {
		TotSNIIElementsEjected += SNIIAllElements[ee];
		TotSNIaElementsEjected += SNIaAllElements[ee];
		TotAGBElementsEjected += AGBAllElements[ee];
	}
	if (TotSNIIElementsEjected > 0.0)
		SNIIElementCorrectionFactor = (SNIIEjectaMass*1.e10/Hubble_h)/TotSNIIElementsEjected;
	if (TotSNIaElementsEjected > 0.0)
		SNIaElementCorrectionFactor = (SNIaEjectaMass*1.e10/Hubble_h)/TotSNIaElementsEjected;
	if (TotAGBElementsEjected > 0.0)
		AGBElementCorrectionFactor = (AGBEjectaMass*1.e10/Hubble_h)/TotAGBElementsEjected;
	for(ee=0;ee<NUM_ELEMENTS;ee++) {
		SNIIAllElements[ee] *= SNIIElementCorrectionFactor;
		SNIaAllElements[ee] *= SNIaElementCorrectionFactor;
		AGBAllElements[ee] *= AGBElementCorrectionFactor;
	}

	//#####
	//2) Correct metal unprocessed masses so that total metals ejected matches total metal elements ejected:
	for(ee=2;ee<NUM_ELEMENTS;ee++) { //Metal elements only
		TotMetalSNIIElementsEjected += SNIIAllElements[ee];
		TotMetalSNIaElementsEjected += SNIaAllElements[ee];
		TotMetalAGBElementsEjected += AGBAllElements[ee];
		TotMetalSNIIUnProcessedElements += SNIIUnProcessedElements[ee];
		TotMetalSNIaUnProcessedElements += SNIaUnProcessedElements[ee];
		TotMetalAGBUnProcessedElements += AGBUnProcessedElements[ee];
	}

	if (*SNIIAllMetals > 0.0)
		SNIIMetalCorrectionFactor = TotMetalSNIIElementsEjected/(*SNIIAllMetals*1.e10/Hubble_h);
	if (*SNIaAllMetals > 0.0)
		SNIaMetalCorrectionFactor = TotMetalSNIaElementsEjected/(*SNIaAllMetals*1.e10/Hubble_h);
	if (*AGBAllMetals > 0.0)
		AGBMetalCorrectionFactor = TotMetalAGBElementsEjected/(*AGBAllMetals*1.e10/Hubble_h);
	*SNIIUnProcessedMetals = TotMetalSNIIUnProcessedElements*(Hubble_h/1.e10);
	*SNIIAllMetals *= SNIIMetalCorrectionFactor;
	*SNIaUnProcessedMetals = TotMetalSNIaUnProcessedElements*(Hubble_h/1.e10);
	*SNIaAllMetals *= SNIaMetalCorrectionFactor;
	*AGBUnProcessedMetals = TotMetalAGBUnProcessedElements*(Hubble_h/1.e10);
	*AGBAllMetals *= AGBMetalCorrectionFactor;

	//#####
	//3) Check that removed (i.e. unprocessed) material doesn't exceed the amount of material currently in the stellar components:
	if (AGBmetals < *AGBUnProcessedMetals && *AGBUnProcessedMetals > 0.0) {
		AGBMetalExcessFactor1 = AGBmetals / *AGBUnProcessedMetals;
		AGBMetalExcessFactor2 = *AGBUnProcessedMetals - AGBmetals; //in code units
		*AGBUnProcessedMetals *= AGBMetalExcessFactor1;
		for(ee=2;ee<NUM_ELEMENTS;ee++) {
			AGBElementExcessFactor[ee] = AGBUnProcessedElements[ee] - (AGBUnProcessedElements[ee]*AGBMetalExcessFactor1);
			AGBUnProcessedElements[ee] *= AGBMetalExcessFactor1;
		}
		*SNIIUnProcessedMetals += AGBMetalExcessFactor2;
		for(ee=2;ee<NUM_ELEMENTS;ee++)
			SNIIUnProcessedElements[ee] += AGBElementExcessFactor[ee];
	}
	if (SNIametals < *SNIaUnProcessedMetals && *SNIaUnProcessedMetals > 0.0) {
		SNIaMetalExcessFactor1 = SNIametals / *SNIaUnProcessedMetals;
		SNIaMetalExcessFactor2 = *SNIaUnProcessedMetals - SNIametals; //in code units
		*SNIaUnProcessedMetals *= SNIaMetalExcessFactor1;
		for(ee=2;ee<NUM_ELEMENTS;ee++) {
			SNIaElementExcessFactor[ee] = SNIaUnProcessedElements[ee] - (SNIaUnProcessedElements[ee]*SNIaMetalExcessFactor1);
			SNIaUnProcessedElements[ee] *= SNIaMetalExcessFactor1;
		}
		*SNIIUnProcessedMetals += SNIaMetalExcessFactor2;
		for(ee=2;ee<NUM_ELEMENTS;ee++)
			SNIIUnProcessedElements[ee] += SNIaElementExcessFactor[ee];
	}
	if (SNIImetals < *SNIIUnProcessedMetals && *SNIIUnProcessedMetals > 0.0) {
		SNIIMetalExcessFactor = SNIImetals / *SNIIUnProcessedMetals;
		*SNIIUnProcessedMetals *= SNIIMetalExcessFactor;
		for(ee=2;ee<NUM_ELEMENTS;ee++)
			SNIIUnProcessedElements[ee] *= SNIIMetalExcessFactor;
	}
#else //INDIVIDUAL_ELEMENTS
	//OLD METHOD: If mass return or unprocessed metals larger than what is currently in the stellar population, re-adjust:
	if(SNIIUnProcessedMetals>SNIImetals)
		SNIIUnProcessedMetals=SNIImetals;

	if(SNIaUnProcessedMetals>SNIametals)
		SNIaUnProcessedMetals=SNIametals;

	if(AGBUnProcessedMetals>AGBmetals)
		AGBUnProcessedMetals=AGBmetals;
#endif //INDIVIDUAL_ELEMENTS
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

