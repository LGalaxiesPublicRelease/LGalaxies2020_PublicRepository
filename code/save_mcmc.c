#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "allvars.h"
#include "proto.h"

#include "mcmc_vars.h"
#include "mcmc_proto.h"


#ifdef MCMC
/** @brief Writes galaxies into a structure to be used by the MCMC */

void save_galaxy_for_mcmc(int gal_index)
{
  //printf("ID=%lld snap=%d sampleID=%d\n", HaloIDs[HaloGal[gal_index].HaloNr].FirstHaloInFOFgroup,
  //		HaloGal[gal_index].SnapNum, SampleIDs[treenr]);
  int snap, fof, ii;
  float low_mass_limit=9.0, high_mass_limit=12.0, log10_StellarMass;
  double log10_Hubble_h, totmetals;

#ifdef MR_PLUS_MRII
  if(Switch_MR_MRII==1)
    {
      low_mass_limit=9.5;
      high_mass_limit=14.0;
    }
  else
    if(Switch_MR_MRII==2)
      {
	low_mass_limit=6.0;
	high_mass_limit=9.5;
      }
#else
  low_mass_limit=7.27;
  high_mass_limit=13.0;
#ifdef MRII
  low_mass_limit=6.0;
  high_mass_limit=11.27;
#endif
#endif


  log10_Hubble_h=log10(Hubble_h);

  //if(HaloGal[gal_index].SnapNum==53)
  //	printf("id=%lld\n",HaloIDs[HaloGal[gal_index].HaloNr].FirstHaloInFOFgroup);
  for(snap=0;snap<NOUT;snap++)
    {
      if((HaloGal[gal_index].DiskMass+HaloGal[gal_index].BulgeMass)>0.)
	{
#ifndef HALOMODEL
	  log10_StellarMass=log10(1E10 * (HaloGal[gal_index].DiskMass+HaloGal[gal_index].BulgeMass)*Hubble_h);
#else
	  log10_StellarMass=log10(1E10 * (HaloGal[gal_index].DiskMass+HaloGal[gal_index].BulgeMass)/Hubble_h);
#endif
	}
      else
	log10_StellarMass=-99.;

      //THE ERROR IS NOW INCLUDED IN mcmc_likelihood.c
      //StellarMass+=gassdev(&MCMCseed)*0.08*(1+MCMCConstraintsZZ[snap]);

      for(fof=0;fof<NFofsInSample[snap]; fof++)
	if( log10_StellarMass > low_mass_limit &&  log10_StellarMass < high_mass_limit &&
	    HaloIDs[HaloGal[gal_index].HaloNr].FirstHaloInFOFgroup == MCMC_FOF[fof].FoFID[snap])
	  {
			  //MCMC_GAL[TotMCMCGals[snap]].StellarMass[snap] = log10(1E10 * ((HaloGal[gal_index].DiskMass+HaloGal[gal_index].BulgeMass)/Hubble_h));
		  	//if((double)((int)((MCMCConstraintsZZ[snap]*10)+0.5)/10.)>0.9)
		  	//MCMC_GAL[TotMCMCGals[snap]].StellarMass[snap] += ran3(&MCMCseed)*0.04*(1+MCMCConstraintsZZ[snap]);
	    MCMC_GAL[TotMCMCGals[snap]].Type[snap] = HaloGal[gal_index].Type;

	    MCMC_GAL[TotMCMCGals[snap]].log10_StellarMass[snap] = log10_StellarMass;
#ifndef HALOMODEL
	    MCMC_GAL[TotMCMCGals[snap]].StellarMass[snap]=1E10 * (HaloGal[gal_index].DiskMass+HaloGal[gal_index].BulgeMass)*Hubble_h;
#else
	    MCMC_GAL[TotMCMCGals[snap]].StellarMass[snap]=1E10 * (HaloGal[gal_index].DiskMass+HaloGal[gal_index].BulgeMass)/Hubble_h;
#endif

	    MCMC_GAL[TotMCMCGals[snap]].ColdGas[snap] = 1E10 * HaloGal[gal_index].ColdGas*Hubble_h;
#ifdef H2_AND_RINGS
	    MCMC_GAL[TotMCMCGals[snap]].HI[snap] = 1E10 * (HaloGal[gal_index].ColdGas*Hubble_h*(1.-HaloGal[gal_index].H2fraction));
#else
	    MCMC_GAL[TotMCMCGals[snap]].HI[snap] = 1E10 * HaloGal[gal_index].ColdGas*Hubble_h*0.54;
#endif


	    MCMC_GAL[TotMCMCGals[snap]].DiskMass[snap] = 1E10 * HaloGal[gal_index].DiskMass * Hubble_h;
	    MCMC_GAL[TotMCMCGals[snap]].BulgeMass[snap] = 1E10 * HaloGal[gal_index].BulgeMass * Hubble_h;
	    MCMC_GAL[TotMCMCGals[snap]].BlackHoleMass[snap] = 1E10 * HaloGal[gal_index].BlackHoleMass; //black hole in units of h^-1
	    totmetals=0.;
	    for(ii=0;ii<NUM_METAL_CHANNELS;ii++)
	      totmetals += HaloGal[gal_index].MetalsDiskMass[ii]+HaloGal[gal_index].MetalsBulgeMass[ii];
	    MCMC_GAL[TotMCMCGals[snap]].MetalsStellarMass[snap] = 1E10 * (totmetals)*Hubble_h; //MassinMetals in units of h^-2
	    totmetals=0.;
	    for(ii=0;ii<NUM_METAL_CHANNELS;ii++)
	      totmetals += HaloGal[gal_index].MetalsColdGas[ii];
	    MCMC_GAL[TotMCMCGals[snap]].MetalsColdGas[snap] =1E10 * totmetals*Hubble_h; //MassinMetals in units of h^-2
#ifdef INDIVIDUAL_ELEMENTS
	    double N_H=0., N_O=0.;
#ifdef H2_AND_RINGS
            //ELEMENTS and RINGS
	    int ii;
	    double SDSS_aperture=30.;
	    for (ii=0;ii<RNUM;ii++)
	      {
		if(RingRadius[ii]/Hubble_h*1000.<SDSS_aperture && HaloGal[gal_index].H2fractionRings[ii]>0.)
		  {
		    N_H+=HaloGal[gal_index].ColdGasRings_elements[ii][0]/1.;
#ifndef MAINELEMENTS
		    N_O+=HaloGal[gal_index].ColdGasRings_elements[ii][4]/16.;
#else
		    N_O+=HaloGal[gal_index].ColdGasRings_elements[ii][2]/16.;
#endif
		  }
	      }

#else
           //ELEMENTS and NO RINGS
	    N_H=HaloGal[gal_index].ColdGas_elements[0]/1.;
#ifndef MAINELEMENTS
	    N_O=HaloGal[gal_index].ColdGas_elements[4]/16.;
#else
	    N_O=HaloGal[gal_index].ColdGas_elements[2]/16.;
#endif

#endif
            //calculate metallicity
	    if((N_O>0.) && (N_H>0.))
	      MCMC_GAL[TotMCMCGals[snap]].GasMetallicity[snap]=1.e12*(N_O/N_H);
	    else
	      MCMC_GAL[TotMCMCGals[snap]].GasMetallicity[snap]=0.;


#else // INDIVIDUAL_ELEMENTS


	    double Metals=0., Mass=0.;
#ifdef H2_AND_RINGS
//NO ELEMENTS and RINGS
	    int jj;
	    double SDSS_aperture=30.;
	    for (jj=0;jj<RNUM;jj++)
		if(RingRadius[jj]/Hubble_h*1000.<SDSS_aperture && HaloGal[gal_index].H2fractionRings[jj]>0.)
		  {
		    for(ii=0;ii<NUM_METAL_CHANNELS;ii++)
		      Metals += HaloGal[gal_index].MetalsColdGasRings[jj][ii];
		    Mass += HaloGal[gal_index].ColdGasRings[jj];
		    //printf("metallicity=%0.5f\n",Metals/Mass/0.0134*pow(10.,8.69));
		  }
#else
//NO ELEMENTS and No RINGS
	    for(ii=0;ii<NUM_METAL_CHANNELS;ii++)
	      Metals += HaloGal[gal_index].MetalsColdGas[ii];
	    Mass += HaloGal[gal_index].ColdGas;
#endif
	    //calculate metallicity
	    MCMC_GAL[TotMCMCGals[snap]].GasMetallicity[snap]=Metals/Mass/0.0114*pow(10.,8.69);
#endif //INDIVIDUAL_ELEMENTS

	    //in units of Solar Masses yr^-1 h-2
	    MCMC_GAL[TotMCMCGals[snap]].Sfr[snap] = HaloGal[gal_index].Sfr * UnitMass_in_g/UnitTime_in_s*SEC_PER_YEAR/SOLAR_MASS*Hubble_h*Hubble_h;

	    MCMC_GAL[TotMCMCGals[snap]].StellarHalfMassRadius[snap] = HaloGal[gal_index].StellarHalfMassRadius*1000.; //Kpc/h


#ifdef COMPUTE_SPECPHOT_PROPERTIES
#ifdef POST_PROCESS_MAGS
    		struct GALAXY_OUTPUT galaxy_output;

    		//#ifndef HALOMODEL
    		//in case of postprocess magnitudes they are only calculates here, inside prepare
    		prepare_galaxy_for_output(snap, &HaloGal[gal_index], &galaxy_output);


    		MCMC_GAL[TotMCMCGals[snap]].MagU[snap] = galaxy_output.MagDust[0]-5.*log10_Hubble_h;
    		////MCMC_GAL[TotMCMCGals[snap]].MagB[snap] = galaxy_output.MagDust[1]-5.*log10_Hubble_h;
    		MCMC_GAL[TotMCMCGals[snap]].MagV[snap] = galaxy_output.MagDust[1]-5.*log10_Hubble_h;
    		MCMC_GAL[TotMCMCGals[snap]].MagJ[snap] = galaxy_output.MagDust[2]-5.*log10_Hubble_h;
    		////MCMC_GAL[TotMCMCGals[snap]].MagK[snap] = galaxy_output.MagDust[4]-5.*log10_Hubble_h;
    		MCMC_GAL[TotMCMCGals[snap]].Magu[snap] = galaxy_output.MagDust[3]-5.*log10_Hubble_h;
    		////MCMC_GAL[TotMCMCGals[snap]].Magg[snap] = galaxy_output.MagDust[6]-5.*log10_Hubble_h;
    		MCMC_GAL[TotMCMCGals[snap]].Magr[snap] = galaxy_output.MagDust[4]-5.*log10_Hubble_h;
    		////MCMC_GAL[TotMCMCGals[snap]].Magi[snap] = galaxy_output.MagDust[8]-5.*log10_Hubble_h;
    		////MCMC_GAL[TotMCMCGals[snap]].Magz[snap] = galaxy_output.MagDust[9]-5.*log10_Hubble_h;
    		//#endif
    		/*MCMC_GAL[TotMCMCGals[snap]].MagV[snap] = galaxy_output.MagDust[0]-5.*log10_Hubble_h;
    		MCMC_GAL[TotMCMCGals[snap]].MagB[snap] = galaxy_output.MagDust[1]-5.*log10_Hubble_h;
    		MCMC_GAL[TotMCMCGals[snap]].MagK[snap] = galaxy_output.MagDust[2]-5.*log10_Hubble_h;*/
#else

    		MCMC_GAL[TotMCMCGals[snap]].MagU[snap] = lum_to_mag(HaloGal[gal_index].LumDust[0][snap])-5.*log10_Hubble_h;
    		MCMC_GAL[TotMCMCGals[snap]].MagV[snap] = lum_to_mag(HaloGal[gal_index].LumDust[1][snap])-5.*log10_Hubble_h;
    		MCMC_GAL[TotMCMCGals[snap]].MagJ[snap] = lum_to_mag(HaloGal[gal_index].LumDust[2][snap])-5.*log10_Hubble_h;
    		MCMC_GAL[TotMCMCGals[snap]].Magu[snap] = lum_to_mag(HaloGal[gal_index].LumDust[3][snap])-5.*log10_Hubble_h;
    		MCMC_GAL[TotMCMCGals[snap]].Magr[snap] = lum_to_mag(HaloGal[gal_index].LumDust[4][snap])-5.*log10_Hubble_h;

    		//MCMC_GAL[TotMCMCGals[snap]].MagV[snap] = lum_to_mag(HaloGal[gal_index].LumDust[0][snap])-5.*log10_Hubble_h;
    		//MCMC_GAL[TotMCMCGals[snap]].MagB[snap] = lum_to_mag(HaloGal[gal_index].LumDust[1][snap])-5.*log10_Hubble_h;
    		//MCMC_GAL[TotMCMCGals[snap]].MagK[snap] = lum_to_mag(HaloGal[gal_index].LumDust[2][snap])-5.*log10_Hubble_h;

#endif //POST_PROCESS_MAGS
#endif //COMPUTE_SPECPHOT_PROPERTIES

#ifdef HALOMODEL
    		if(snap==0)
    		{
    			MCMC_GAL[TotMCMCGals[snap]].fofid[snap] = fof;
    			//MCMC_GAL[TotMCMCGals[snap]].M_Crit200[snap] = log10(Halo[HaloGal[gal_index].HaloNr].Len*PartMass*1.e10);
    			MCMC_GAL[TotMCMCGals[snap]].M_Crit200[snap] = log10(Halo[HaloGal[gal_index].HaloNr].M_Crit200*1.e10);
    			MCMC_GAL[TotMCMCGals[snap]].M_Mean200[snap] = log10(Halo[HaloGal[gal_index].HaloNr].M_Mean200*1.e10);
#ifdef MCRIT
    			MCMC_GAL[TotMCMCGals[snap]].M_Mean200[snap] = log10(Halo[HaloGal[gal_index].HaloNr].M_Crit200*1.e10);
#endif
    			MCMC_GAL[TotMCMCGals[snap]].x[snap] = HaloGal[gal_index].Pos[0];
    			MCMC_GAL[TotMCMCGals[snap]].y[snap] = HaloGal[gal_index].Pos[1];
    			MCMC_GAL[TotMCMCGals[snap]].z[snap] = HaloGal[gal_index].Pos[2];
    			MCMC_GAL[TotMCMCGals[snap]].Type[snap] = HaloGal[gal_index].Type;
    			MCMC_GAL[TotMCMCGals[snap]].ngal[snap] = 0;
    		}
#endif

    		MCMC_GAL[TotMCMCGals[snap]].Weight[snap] = MCMC_FOF[fof].Weight[snap];

    		//#ifdef MR_PLUS_MRII
    		////

    		//NOW GET PROPERTIES FOR FOF GROUPS - done ofr the particular fof that current galaxy resides in

    		++MCMC_FOF[fof].NGalsInFoF[snap];

#ifdef HALOMODEL
    		if(HaloGal[gal_index].Type==0)
    		{
    			MCMC_FOF[fof].IndexOfCentralGal[snap]=TotMCMCGals[snap];
    			//MCMC_FOF[fof].M_Crit200[snap] = log10(Halo[HaloGal[gal_index].HaloNr].Len*PartMass*1.e10);
    			MCMC_FOF[fof].M_Crit200[snap] = log10(Halo[HaloGal[gal_index].HaloNr].M_Crit200*1.e10);
    			MCMC_FOF[fof].M_Mean200[snap] = log10(Halo[HaloGal[gal_index].HaloNr].M_Mean200*1.e10);
#ifdef MCRIT
    			MCMC_FOF[fof].M_Mean200[snap] = log10(Halo[HaloGal[gal_index].HaloNr].M_Crit200*1.e10);
#endif
    		}
#endif

    		++TotMCMCGals[snap];

    		if(TotMCMCGals[snap] > MCMCAllocFactor)
    			terminate("Maximum number of galaxies in MCMC structure reached. Increase MCMCSmartFactor\n");
    	  }
      }

    return;
}
#endif


