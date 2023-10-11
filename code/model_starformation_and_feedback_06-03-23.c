#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>

#include "allvars.h"
#include "proto.h"

/** @file recipe_starformation_and_feedback.c
 *  @brief recipe_starformation_and_feedback.c computes the amount of stars
 *         formed from the cold gas, the amount of gas reheated from cold to hot
 *         and the amount of gas ejected from hot to external.
 *
 * The routine is divided in two parts, star formation and SN feedback, with a
 * number of different implementations controlled by input parameters.
 *
 *
 *  0 -\f$M_{\rm{crit}}=3.8\times 10^9
 *     \left(\frac{V_{\rm{max}}}{200\,\rm{km s}^{-1}}\right)
 *     \left(\frac{r_{\rm{disk}}}{10\,\rm{kpc}}\right)M_{\odot}\f$
 *     (Eq. 16 Guo2010) (StarFormationModel = 0), \n
 *        - same as 0 but using \f$V_{\rm{max}}\f$ or \f$V_{\rm{max,infall}}\f$
 *          instead of \f$V_{\rm{vir}}\f$ and allowing SF in satellites. *

 *
 * There are 2 options for the <B>SN Feedback Recipe</B>:
 *
 * 0 - \f$\epsilon_{\rm{disk}}=\epsilon
 *      \biggl[0.5+\left(\frac{V_{\rm{max}}}{70km/s}\right)^{-\beta_1}\biggr]\f$,
 *     \f$\epsilon_{\rm{halo}}=\eta
 *      \biggl[0.5+\left(\frac{V_{\rm{max}}}{70km/s}\right)^{-\beta_2}\biggr]\f$
 *     (Eqs. 19 & 21 Guo2010)(FeedbackReheatingModel = 0)
 *
 *
 * Also, Guo2010 alowed for type 1 satellite to have gas cycles and receive
 * gas from their own satellites when these are outside Rvir of the type 0.
 * */


/** @brief Main recipe, calculates the fraction of cold gas turned into stars due
  *        to star formation; the fraction of mass instantaneously recycled and
  *        returned to the cold gas; the fraction of gas reheated from cold to hot,
  *        ejected from hot to external and returned from ejected to hot due to
  *        SN feedback.   */
void starformation(int p, int centralgal, double time, double dt, int nstep)
{
  /* Variables: reff-Rdisk, tdyn=Rdisk/Vmax, strdot=Mstar_dot, stars=strdot*dt*/
  double tdyn, strdot=0., stars, cold_crit;
  int ii;
#ifdef H2_AND_RINGS
  double strdotr[RNUM], starsRings[RNUM];
  double sfe, cold_crit_rate, SigmaGas, SigmaGasratio;
  //double sfe, vmax, cold_crit_rate, SigmaGas, SigmaGasratio, sigmah50=0.0;
  int j;
#endif
#ifdef COMPUTE_SPECPHOT_PROPERTIES
#ifndef POST_PROCESS_MAGS
  double metallicitySF;
#endif
#endif
  double Vmax, gas_radius;

  if(Gal[p].Type == 0)
    Vmax=Gal[p].Vmax;
  else
    Vmax=Gal[p].InfallVmax;

  gas_radius = Gal[p].ColdGasRadius;
  //the h factor in Gal[p].ColdGasRadius goes into cold_crit
  tdyn = gas_radius / Vmax;
  cold_crit = SfrColdCrit * Vmax/200. * gas_radius*100.;

  //standard star formation law (Croton2006, Delucia2007, Guo2010, Henriques2015)
  if(StarFormationModel == 0)
    {
      if(Gal[p].ColdGas > cold_crit)
	strdot = SfrEfficiency * (Gal[p].ColdGas - cold_crit) / tdyn;
      else
	strdot = 0.0;
    }

#ifdef H2_AND_RINGS
	for (int ee=0;ee<NUM_ELEMENTS;ee++)
		if (Gal[p].EjectedMass_elements[ee] < 0.0)
			printf("IN starformation() 1: Gal[%i].EjectedMass_elements[%i] = %e\n", p, ee, Gal[p].EjectedMass_elements[ee]);

  update_h2fraction(p);

  sfe=SfrEfficiency*UnitTime_in_years/Hubble_h; //convert from yr-1 into code units of time // the unit of sfe here is (km/s)/(Mpc/h)=1.022e-3 h^{+1} /Gyr
  //sfe=SfrEfficiency*UnitTime_in_years/Hubble_h/pow((1+ZZ[Gal[p].SnapNum]),1.0); //convert from yr-1 into code units of time

  if(SFRtdyn==1)
    sfe= (sfe/tdyn)/UnitTime_in_years*Hubble_h; // for star formation rate proportional to 1/t_dyn
    //sfe= sfe/1.8/tdyn; // for star formation rate proportional to 1/t_dyn

  for(j=0;j<RNUM;j++)
    {
      if(StarFormationModel == 0)
    	  strdotr[j]=strdot*Gal[p].ColdGasRings[j]/Gal[p].ColdGas;
      else if(StarFormationModel == 2)
      {
    	  if(Gal[p].Type == 0)
    		  cold_crit_rate = SfrColdCrit * Gal[p].Vmax/200. * Gal[p].ColdGasRadius/Gal[p].ColdGas*100.;
    	  else
    		  cold_crit_rate = SfrColdCrit * Gal[p].InfallVmax/200. * Gal[p].ColdGasRadius/Gal[p].ColdGas*100.;

    	  if(cold_crit_rate < 1 && cold_crit_rate>=0)
    		  strdotr[j] = sfe * Gal[p].ColdGasRings[j] * (1 - cold_crit_rate);
    	  else strdotr[j] = 0.0;
      }

      else if(StarFormationModel == 3) /*The star formation law in Krumholz et al. 2009*/
      {
    	  double SigmaGas0=85.0, SF_Law_pow=0.33;

    	  if(j==0)
    		  SigmaGas = Gal[p].ColdGasRings[j] / (M_PI* RingRadius[j]*RingRadius[j])/WARM_PHASE_FACTOR*Clumpingfactor;
    	  else
    		  SigmaGas = Gal[p].ColdGasRings[j] / (M_PI*(RingRadius[j]*RingRadius[j]-RingRadius[j-1]*RingRadius[j-1]))/WARM_PHASE_FACTOR*Clumpingfactor;

    	  /* convert from 10^10 M_sun/h / (Mpc/h)^2 to (M_sun/pc^2) */
    	  SigmaGas=SigmaGas*0.01*Hubble_h;
    	  SigmaGasratio=pow(SigmaGas/SigmaGas0,SF_Law_pow);

    	  if(SigmaGasratio<1.0 && SigmaGasratio>0.0)
    		  strdotr[j] = sfe/SigmaGasratio * Gal[p].ColdGasRings[j]*Gal[p].H2fractionRings[j]/WARM_PHASE_FACTOR ;
    	  //Only cold H2 component is proportional to star formation rate.
    	  else {
    		  if(SigmaGasratio>=1.0)
    			  strdotr[j] = sfe*SigmaGasratio * Gal[p].ColdGasRings[j]*Gal[p].H2fractionRings[j]/WARM_PHASE_FACTOR ;
    		  else strdotr[j]=0.0;
    	  }
      }

      else if(StarFormationModel == 4)	/*The star formation law in Fu et al. 2010*/
      {
    	  //double sigma_H2, N_sf=1.0, sigma2_crit=70, area;

    	  if(Gal[p].H2fractionRings[j]>=0.0) {
    		  strdotr[j] = sfe * Gal[p].ColdGasRings[j]*Gal[p].H2fractionRings[j] / WARM_PHASE_FACTOR;
    		  /*if(j >= 6 && j <= 8) {
	    	  	  strdotr[j] = 0.0; //ROB: +++++BARS TEST!: What happens if no SF is allowed to occur between 940pc and 7.53kpc, due to bars? (05-05-20)
	      	  }*/
    	  }
    	  else strdotr[j]=0.0;
      }
      else strdotr[j] = 0.0;
    }

  for (j=0,strdot=0;j<RNUM;j++)
    strdot+=strdotr[j];

#endif // H2_AND_RINGS



  /* Note that Units of dynamical time are Mpc/Km/s - no conversion on dt needed
   * be mentioned 3.06e19 to 3.15e19 */

  if(strdot < 0.0)
    strdot =0.;
    //terminate("***error stars<0.0***\n");

  stars = strdot * dt;

#ifdef H2_AND_RINGS
  for(j=0;j<RNUM;j++)
    {
      if(strdotr[j] < 0.0)
	strdotr[j] =0.;
      starsRings[j] = strdotr[j] * dt;
      if(starsRings[j] <0.0) starsRings[j] =0.0;
    }
#endif

  //otherwise this check is done inside update_stars_due_to_reheat for stars+reheatedmass!>coldgas
#ifdef FEEDBACK_COUPLED_WITH_MASS_RETURN
  if(stars > Gal[p].ColdGas)
  	stars = Gal[p].ColdGas;
#ifdef H2_AND_RINGS
  for(j=0;j<RNUM;j++)
    if(starsRings[j]>Gal[p].ColdGasRings[j])
      starsRings[j]=Gal[p].ColdGasRings[j];
#endif
#endif

  mass_checks(p,"model_starformation_and_feedback.c",__LINE__);
  mass_checks(centralgal,"model_starformation_and_feedback.c",__LINE__);

  /* Store the value of the metallicity of the cold phase when SF occurs */

#ifdef COMPUTE_SPECPHOT_PROPERTIES
#ifndef POST_PROCESS_MAGS
  metallicitySF=0.;
  if (Gal[p].ColdGas > 0.)
    {
      for(ii=0;ii<NUM_METAL_CHANNELS;ii++)
	metallicitySF+= Gal[p].MetalsColdGas[ii];
      metallicitySF/=Gal[p].ColdGas;

    }
#endif
#endif

//if FEEDBACK_COUPLED_WITH_MASS_RETURN feedback happens only when stars die,
//there is no need to balance it with SF
#ifndef FEEDBACK_COUPLED_WITH_MASS_RETURN
  if (stars > 0.)
#ifndef H2_AND_RINGS
    update_stars_due_to_reheat(p, centralgal, &stars);
#else
    update_stars_due_to_reheat(p, centralgal, &stars, starsRings);
#endif
#endif //FEEDBACK_COUPLED_WITH_MASS_RETURN


/* if update_stars_due_to_reheat is commented out, uncomment this lines:
   if(stars > Gal[p].ColdGas)
    stars = Gal[p].ColdGas;
 #ifdef H2_AND_RINGS
   for(j=0;j<RNUM;j++)
     if(starsRings[j]>Gal[p].ColdGasRings[j])
       starsRings[j]=Gal[p].ColdGasRings[j];
 #endif*/


  mass_checks(p,"model_starformation_and_feedback.c",__LINE__);
  mass_checks(centralgal,"model_starformation_and_feedback.c",__LINE__);

  /*  update the star formation rate */
   /*Sfr=stars/(dt*steps)=strdot*dt/(dt*steps)=strdot/steps -> average over the STEPS*/
   Gal[p].Sfr += stars / (dt * STEPS);
   Gal[p].SfrInst = stars / dt; //*****ROB*****//
 #ifdef H2_AND_RINGS
   for(j=0;j<RNUM;j++)
   {
     Gal[p].SfrRings[j] += starsRings[j] / (dt * STEPS);
     Gal[p].SfrInstRings[j] = starsRings[j] / dt; //*****ROB*****//
   }
 #endif

	for (int ee=0;ee<NUM_ELEMENTS;ee++)
		if (Gal[p].EjectedMass_elements[ee] < 0.0)
			printf("IN starformation() 2: Gal[%i].EjectedMass_elements[%i] = %e\n", p, ee, Gal[p].EjectedMass_elements[ee]);

  // update_from_star_formation can only be called
  // after SN_feedback recipe since stars need to be re_set once the reheated mass is known
  // (star formation and feedback share the same fraction of cold gas)
  if (stars > 0.)
 #ifndef H2_AND_RINGS
    update_from_star_formation(p, stars, "insitu", nstep); // false indicates not a burst
 #else
    update_from_star_formation(p, stars, starsRings, "insitu", nstep); // false indicates not a burst
 #endif
	for (int ee=0;ee<NUM_ELEMENTS;ee++)
		if (Gal[p].EjectedMass_elements[ee] < 0.0)
			printf("IN starformation() 3: Gal[%i].EjectedMass_elements[%i] = %e\n", p, ee, Gal[p].EjectedMass_elements[ee]);

#ifdef COMPUTE_SPECPHOT_PROPERTIES
#ifndef POST_PROCESS_MAGS
  //  Update the luminosities due to the stars formed
  if (stars > 0.0)
    add_to_luminosities(p, stars, time, dt, metallicitySF);
#endif
#endif

  mass_checks(p,"model_starformation_and_feedback.c",__LINE__);
  mass_checks(centralgal,"model_starformation_and_feedback.c",__LINE__);

  update_massweightage(p, stars, time);

#ifndef FEEDBACK_COUPLED_WITH_MASS_RETURN
  /* ifdef FEEDBACK_COUPLED_WITH_MASS_RETURN feedback is only called when stars die,
   * inside DETAILED_METALS_AND_MASS_RETURN */
  if (stars > 0.)
#ifndef H2_AND_RINGS
    SN_feedback(p, centralgal, stars, "ColdGas", dt);
#else
    SN_feedback(p, centralgal, stars, starsRings, "ColdGas", dt);
#endif
#endif

  mass_checks(p,"model_starformation_and_feedback.c",__LINE__);
  mass_checks(centralgal,"model_starformation_and_feedback.c",__LINE__);


  if(DiskInstabilityModel==0)
    {
   //if(Gal[p].ColdGas > 0.0)
    //     check_disk_instability_gas(p,dt);
    if(Gal[p].DiskMass > 0.0)
         check_disk_instability(p,dt);
    }

  if (DiskRadiusModel == 0)
    Gal[p].DiskRadius=get_stellar_disk_radius(p);

  mass_checks(p,"model_starformation_and_feedback.c",__LINE__);
  mass_checks(centralgal,"model_starformation_and_feedback.c",__LINE__);

}




#ifndef H2_AND_RINGS
void update_stars_due_to_reheat(int p, int centralgal, double *stars)
#else
void update_stars_due_to_reheat(int p, int centralgal, double *stars, double starsRings[])
#endif
{
  double reheated_mass, frac, Radius_low=0., totmetals;
  int ii;

#ifndef H2_AND_RINGS
  totmetals=0.;
  for(ii=0;ii<NUM_METAL_CHANNELS;ii++)
    totmetals+= Gal[p].MetalsColdGas[ii];
  reheated_mass=compute_SN_reheat(p, centralgal, *stars, Gal[p].ColdGas, totmetals, Radius_low, Gal[p].ColdGasRadius);
  if((*stars + reheated_mass) > Gal[p].ColdGas)
    {
      frac = Gal[p].ColdGas / (*stars + reheated_mass);
      *stars *= frac;
    }
#else
  int jj;
  for(jj=0;jj<RNUM;jj++)
    {
      if(jj>0)
    	Radius_low=RingRadius[jj-1];

      totmetals=0.;
      for(ii=0;ii<NUM_METAL_CHANNELS;ii++)
	totmetals+= Gal[p].MetalsColdGasRings[jj][ii];

      reheated_mass=compute_SN_reheat(p, centralgal, starsRings[jj], Gal[p].ColdGasRings[jj], totmetals, Radius_low, RingRadius[jj]);
      if((starsRings[jj] + reheated_mass) > Gal[p].ColdGasRings[jj])
	{
	  frac = Gal[p].ColdGasRings[jj] / (starsRings[jj] + reheated_mass);
  	  starsRings[jj] *= frac;
  	}
       *stars+=starsRings[jj];
    }
#endif
}




/** @brief Updates the different components due to star formation: mass
  *        and metals in stars and cold gas and stellar spin. */
//void update_from_star_formation(int p, double time, double stars, double metallicity)
#ifndef H2_AND_RINGS
void update_from_star_formation(int p, double stars, char type_of_event[], int nstep)
#else
void update_from_star_formation(int p, double stars, double starsRings[], char type_of_event[], int nstep)
#endif
{
  int ii;
  double stars_to_add=0., NonRecycledFraction=0.;
#ifdef H2_AND_RINGS
  int jj;
  double stars_to_addr[RNUM], fractionRings[RNUM];
#else
  double fraction;
#endif
#ifdef DETAILED_DUST
  double fractionClouds, fractionDiff=0., cloud_mass, diff_mass; //ROB: These were set as floats before. (02-02-22)
#ifdef H2_AND_RINGS
  double cloud_massRings[RNUM], diff_massRings[RNUM];
  double fractionCloudsRings[RNUM], fractionDiffRings[RNUM];
#endif
  double previoustime, newtime, deltaT;
  previoustime = NumToTime(Gal[p].SnapNum);
  newtime = NumToTime(Gal[p].SnapNum+1);
  deltaT = previoustime - newtime;
#endif

  if(Gal[p].ColdGas <= 0. || stars <= 0.) {
    printf("Gal[p].ColdGas <= 0. || stars <= 0., Coldgas=%0.5e stars=%0.5e, in function update_from_star_formation, model_starformation_and_feedback.c line:%d\n",Gal[p].ColdGas, stars,__LINE__);
    exit(0);
  }

  /* If DETAILED_METALS_AND_MASS_RETURN, no longer an assumed instantaneous
   * recycled fraction. Mass is returned over time via SNe and AGB winds.
   * Update the Stellar Spin when forming stars */
#ifndef DETAILED_METALS_AND_MASS_RETURN
  NonRecycledFraction=(1 - RecycleFraction);
#else
  NonRecycledFraction=1.;
#endif

#ifndef H2_AND_RINGS
  stars_to_add=NonRecycledFraction*stars;
#else
  for(jj=0;jj<RNUM;jj++)
    {
      stars_to_addr[jj]=NonRecycledFraction * starsRings[jj];
      stars_to_add+=stars_to_addr[jj];
    }
#endif //H2_AND_RINGS

  if (Gal[p].DiskMass+stars_to_add > 1.e-8)
    for (ii = 0; ii < 3; ii++)
      Gal[p].DiskSpin[ii]=((Gal[p].DiskSpin[ii])*(Gal[p].DiskMass)+stars_to_add*Gal[p].ColdGasSpin[ii])/(Gal[p].DiskMass+stars_to_add);

  mass_checks(p,"model_starformation_and_feedback.c",__LINE__);

#ifdef DETAILED_DUST
  /* Here, we calculate how much of the gas and dust in the clouds and diffuse sub-components are transferred to the stellar disc
   * during star formation. All the cloud mass is used-up first. If the total mass of newly-formed stars exceeds the cloud
   * sub-component mass, diffuse gas is also used-up.
   */
  int ee;
  //TOTAL:
  cloud_mass = 0.;
  diff_mass = 0.;
  for(ee=0;ee<NUM_ELEMENTS;ee++) {
	  cloud_mass += Gal[p].ColdGasClouds_elements[ee]/(1E10/Hubble_h); //in code units (i.e. 1e10/Hubble_h)
	  diff_mass += Gal[p].ColdGasDiff_elements[ee]/(1E10/Hubble_h); //in code units (i.e. 1e10/Hubble_h)
  }
  fractionClouds = stars_to_add/cloud_mass;
  if (fractionClouds > 1.) {
      fractionClouds = 1.;
      fractionDiff =  (stars_to_add - cloud_mass)/diff_mass;
      if (fractionDiff > 1.)
    	  fractionDiff = 1.;
  }
  else fractionDiff = 0.;
#ifdef H2_AND_RINGS
  //RINGS:
  for (jj=0;jj<RNUM;jj++) {
	  cloud_massRings[jj] = 0.0;
	  diff_massRings[jj] = 0.0;
	  fractionCloudsRings[jj] = 0.0;
	  fractionDiffRings[jj] = 0.0;
  }
  for(ee=0;ee<NUM_ELEMENTS;ee++) {
  	  for (jj=0;jj<RNUM;jj++) {
  		if(Gal[p].ColdGasRings_elements[jj][ee]>0.) {
  			cloud_massRings[jj] += Gal[p].ColdGasCloudsRings_elements[jj][ee]/(1E10/Hubble_h); //in code units (i.e. 1e10/Hubble_h)
  			diff_massRings[jj] += Gal[p].ColdGasDiffRings_elements[jj][ee]/(1E10/Hubble_h); //in code units (i.e. 1e10/Hubble_h)
  		}
  	  }
  }
  for (jj=0;jj<RNUM;jj++) {
	  if (cloud_massRings[jj] > 0.0) {
		  fractionCloudsRings[jj] = stars_to_addr[jj]/cloud_massRings[jj];
		  if (fractionCloudsRings[jj] > 1.0) {
			  fractionCloudsRings[jj] = 1.0;
			  fractionDiffRings[jj] =  (stars_to_addr[jj] - cloud_massRings[jj])/diff_massRings[jj];
			  if (fractionDiffRings[jj] > 1.0)
				  fractionDiffRings[jj] = 1.0;
		  }
		  else fractionDiffRings[jj] = 0.0;
	  }
	  /*else { //ROB: Now pre-set to zero in the rings loop above. (08-03-22)
		  fractionCloudsRings[jj] = 0.0;
		  fractionDiffRings[jj] = 0.0;
	  }*/
  }
#endif //H2_AND_RINGS
#endif //DETAILED_DUST

  /*double tot_ele_clouds=0.0;
  double tot_ele_diff=0.0;
  double met_ele_clouds=0.0;
  double met_ele_diff=0.0;
  for (int eee=0;eee<NUM_ELEMENTS;eee++) {
	  tot_ele_clouds += Gal[p].ColdGasClouds_elements[eee];
	  tot_ele_diff += Gal[p].ColdGasDiff_elements[eee];
	  if (eee > 1) {
		  met_ele_clouds += Gal[p].ColdGasClouds_elements[eee];
		  met_ele_diff += Gal[p].ColdGasDiff_elements[eee];
	  }
  }
  if (p == 0 && Gal[p].SnapNum < 16) printf("1b) Before transfer: Snap = %i | Step = %i | Z_cloud = %e | Z_diff = %e | fracClouds = %f | fracDiff = %f | \n  stars_to_addr = [%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f] | \n  fracCloudsRings = [%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f] | \n  fracDiffRings = [%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f]\n",
		  	  	  	  Gal[p].SnapNum, nstep, met_ele_clouds/tot_ele_clouds, met_ele_diff/tot_ele_diff, fractionClouds, fractionDiff,
					  stars_to_addr[0], stars_to_addr[1], stars_to_addr[2], stars_to_addr[3], stars_to_addr[4], stars_to_addr[5], stars_to_addr[6], stars_to_addr[7], stars_to_addr[8], stars_to_addr[9], stars_to_addr[10], stars_to_addr[11],
					  fractionCloudsRings[0], fractionCloudsRings[1], fractionCloudsRings[2], fractionCloudsRings[3], fractionCloudsRings[4], fractionCloudsRings[5], fractionCloudsRings[6], fractionCloudsRings[7], fractionCloudsRings[8], fractionCloudsRings[9], fractionCloudsRings[10], fractionCloudsRings[11],
					  fractionDiffRings[0], fractionDiffRings[1], fractionDiffRings[2], fractionDiffRings[3], fractionDiffRings[4], fractionDiffRings[5], fractionDiffRings[6], fractionDiffRings[7], fractionDiffRings[8], fractionDiffRings[9], fractionDiffRings[10], fractionDiffRings[11]);
  */

  /*for (int ee=0;ee<NUM_ELEMENTS;ee++)
  			for (int jj=0;jj<RNUM;jj++)
  				if (Gal[p].DiskMassRings_elements[jj][ee] < 0.0) printf("BEFORE transfer_material_with_rings(): Gal[%lld].DiskMassRings_elements[%i][%i] = %e\n", p, jj, ee, Gal[p].DiskMassRings_elements[jj][ee]);*/

#ifdef H2_AND_RINGS
  /*  Update Gas and Metals from star formation */
  for(jj=0;jj<RNUM;jj++)
   if(Gal[p].ColdGasRings[jj]>0.)
     fractionRings[jj]=stars_to_addr[jj]/Gal[p].ColdGasRings[jj];
   else
     fractionRings[jj]=0.;
#ifdef DETAILED_DUST
  transfer_material_with_rings(p,"DiskMass",p,"ColdGas",fractionRings,fractionCloudsRings,fractionDiffRings,"model_starformation_and_feedback.c", __LINE__);
  for(ee=0;ee<NUM_ELEMENTS;ee++)
	  Gal[p].DustColdGasRates[4] += (Gal[p].DustColdGasDiff_elements[ee]*fractionDiff + Gal[p].DustColdGasClouds_elements[ee]*fractionClouds)/(deltaT * UnitTime_in_years);
#else //DETAILED_DUST
  transfer_material_with_rings(p,"DiskMass",p,"ColdGas",fractionRings,"model_starformation_and_feedback.c", __LINE__);
#endif //DETAILED_DUST
#else //H2_AND_RINGS
  fraction=stars_to_add/Gal[p].ColdGas;
#ifdef DETAILED_DUST
  transfer_material(p,"DiskMass",p,"ColdGas",fraction, fractionClouds, fractionDiff,"model_starformation_and_feedback.c", __LINE__);
  //Update dust destruction rate:
  //Gal[p].DustColdGasRates[4] += (elements_total(elements_add(elements_init(),Gal[p].DustColdGasDiff_elements,fraction_diffuse)) + elements_total(elements_add(elements_init(),Gal[p].DustColdGasClouds_elements,fraction_clouds)))/(deltaT * UnitTime_in_years);
  for(ee=0;ee<NUM_ELEMENTS;ee++)
	  Gal[p].DustColdGasRates[4] += (Gal[p].DustColdGasDiff_elements[ee]*fractionDiff + Gal[p].DustColdGasClouds_elements[ee]*fractionClouds)/(deltaT * UnitTime_in_years);
#else
  transfer_material(p,"DiskMass",p,"ColdGas",fraction, "model_starformation_and_feedback.c", __LINE__);
#endif
#endif

  /*for (int ee=0;ee<NUM_ELEMENTS;ee++)
    			for (int jj=0;jj<RNUM;jj++)
    				if (Gal[p].DiskMassRings_elements[jj][ee] < 0.0) printf("AFTER transfer_material_with_rings(): Gal[%lld].DiskMassRings_elements[%i][%i] = %e\n", p, jj, ee, Gal[p].DiskMassRings_elements[jj][ee]);*/

  mass_checks(p,"model_starformation_and_feedback.c",__LINE__);

  /*tot_ele_clouds=0.0;
  tot_ele_diff=0.0;
  met_ele_clouds=0.0;
  met_ele_diff=0.0;
  for (int eee=0;eee<NUM_ELEMENTS;eee++) {
    tot_ele_clouds += Gal[p].ColdGasClouds_elements[eee];
    tot_ele_diff += Gal[p].ColdGasDiff_elements[eee];
    if (eee > 1) {
	    met_ele_clouds += Gal[p].ColdGasClouds_elements[eee];
	    met_ele_diff += Gal[p].ColdGasDiff_elements[eee];
    }
  }
  if (p == 0 && Gal[p].SnapNum < 16) printf("1c) After transfer: Snap = %i | Step = %i | Z_cloud = %e | Z_diff = %e\n",
		  	  	  	  Gal[p].SnapNum, nstep, met_ele_clouds/tot_ele_clouds, met_ele_diff/tot_ele_diff);
  */

/*#ifdef DETAILED_DUST
  transfer_dust_from_starformation(p, fractionDiff, fractionClouds);
  //transfer_dust_from_starformation_with_rings(p, fractionDiffRings, fractionCloudsRings);
#endif*/


#ifdef TRACK_MASSGROWTH_CHANNELS
  //for this calculation we want just the long lived mass and
  //take the instantaneous recycling aproximation even
  //for the detailed chemical enrichment because it is not
  //possible to know which component to eject mass from afterwards
  double long_lived_mass;
  long_lived_mass=stars_to_add;
#ifdef DETAILED_METALS_AND_MASS_RETURN
  long_lived_mass*=(1 - RecycleFraction);
#endif

  if(strcmp(type_of_event,"insitu")==0)
      {
        Gal[p].MassFromInSitu+=long_lived_mass;
#ifdef STAR_FORMATION_HISTORY
#ifdef TRACK_SFH_MASSGROWTH_CHANNELS
        Gal[p].sfh_MassFromInSitu[Gal[p].sfh_ibin]+=long_lived_mass;
#endif
#endif
      }

  if(strcmp(type_of_event,"merger")==0)
    {
      Gal[p].MassFromBursts+=long_lived_mass;
#ifdef STAR_FORMATION_HISTORY
#ifdef TRACK_SFH_MASSGROWTH_CHANNELS
      Gal[p].sfh_MassFromBursts[Gal[p].sfh_ibin]+=long_lived_mass;
#endif
#endif
    }
#endif //TRACK_MASSGROWTH_CHANNELS

#ifdef TRACK_BURST
  if(strcmp(type_of_event,"merger")==0)
    {
      Gal[p].BurstMass+=stars_to_add;
#ifdef STAR_FORMATION_HISTORY
      Gal[p].sfh_BurstMass[Gal[p].sfh_ibin]+=stars_to_add;
#endif
    }
#endif


  if (FeedbackReheatingModel == 0 || FeedbackReheatingModel == 1)
    {
 /* stars (instead of star_to_add) used because the Yield is defined as a
  * fraction of all stars formed, not just long lived */
#ifdef DETAILED_METALS_AND_MASS_RETURN
#ifdef METALS_SELF
      Gal[p].MetalsHotGasSelf.str.type2 += Yield * FracZSNIItoHot * stars;
#endif
#else //IF NOT DETAILED_METALS_AND_MASS_RETURN
      //This part is not used if OPT+=DELAYED_ENRICHMENT_AND MASS_RETURN as yield
      //and recycling fraction are not fixed:
#ifndef H2_AND_RINGS
      Gal[p].MetalsColdGas[0] += Yield * (1.-FracZSNIItoHot) * stars;
#else
      for(jj=0;jj<RNUM;jj++)
	{
	  Gal[p].MetalsColdGasRings[jj][0] += Yield* (1.-FracZSNIItoHot) * starsRings[jj];
	  Gal[p].MetalsColdGas[0] += Yield* (1.-FracZSNIItoHot) * starsRings[jj];
	}
#endif
      Gal[Gal[p].CentralGal].MetalsHotGas[0] += Yield * FracZSNIItoHot * stars;
      Gal[Gal[p].CentralGal].HotGas += Yield * FracZSNIItoHot * stars;
#ifdef METALS_SELF
      Gal[p].MetalsHotGasSelf[0] += Yield * FracZSNIItoHot * stars;
#endif
#endif //DETAILED_METALS_AND_MASS_RETURN
    }

  mass_checks(p,"model_starformation_and_feedback.c",__LINE__);

  if (DiskRadiusModel == 0)
    Gal[p].DiskRadius = get_stellar_disk_radius(p);


}






/* there are two modes for supernova feedback corresponding to when the mass returning
 * by dying stars is returned to the cold gas - reheat and ejection; and when the mass
 * is returned to the hot gas - only ejection.*/
#ifndef H2_AND_RINGS
void SN_feedback(int p, int centralgal, double stars, char feedback_location[], double dt)
#else
void SN_feedback(int p, int centralgal, double stars, double starsRings[], char feedback_location[], double dt)
#endif
{
  double EjectVmax, EjectVvir, SN_Energy, Reheat_Energy, ReScaled_EnergySNcode;
  double reheated_mass=0., ejected_mass=0., totmetals;
  double Radius_low;
  int ii;
  /* SN FEEDBACK RECIPES */
#ifdef H2_AND_RINGS
  int jj, mm;
  double reheated_massr[RNUM], totmetalsRings;
#endif

	for (int ee=0;ee<NUM_ELEMENTS;ee++)
		if (Gal[p].EjectedMass_elements[ee] < 0.0)
			printf("IN SN_feedback() 1: Gal[%i].EjectedMass_elements[%i] = %e\n", p, ee, Gal[p].EjectedMass_elements[ee]);

  Radius_low=0.;
  //REHEAT
#ifndef H2_AND_RINGS
  //when FEEDBACK_COUPLED_WITH_MASS_RETURN some mass goes into HOTGAS and does not produce reheating
  if (strcmp(feedback_location,"HotGas")==0)
    reheated_mass=0;
  else
    {
      totmetals=0.;
      for(ii=0;ii<NUM_METAL_CHANNELS;ii++)
     	totmetals+= Gal[p].MetalsColdGas[ii];
      reheated_mass=compute_SN_reheat(p, centralgal, stars, Gal[p].ColdGas, totmetals, Radius_low, Gal[p].ColdGasRadius);
    }
#else
  reheated_mass=0.0;
  //stars=0;
  for(jj=0;jj<RNUM;jj++)
    {
      if(jj>0)
	Radius_low=RingRadius[jj-1];

      //when FEEDBACK_COUPLED_WITH_MASS_RETURN some mass goes into HOTGAS and does not produce reheating
      if (strcmp(feedback_location,"HotGas")==0)
	  reheated_massr[jj]=0.;
      else
	{
	  totmetals=0.;
	  for(ii=0;ii<NUM_METAL_CHANNELS;ii++)
	    totmetals+= Gal[p].MetalsColdGasRings[jj][ii];
	  reheated_massr[jj]=compute_SN_reheat(p, centralgal, starsRings[jj], Gal[p].ColdGasRings[jj], totmetals, Radius_low, RingRadius[jj]);
	}

      reheated_mass+=reheated_massr[jj];
      //stars+=starsRings[jj];
    }

  //reheated_mass > Gal[p].ColdGas might happen due to precision
  if(reheated_mass > Gal[p].ColdGas)
     reheated_mass = Gal[p].ColdGas;

#endif

	for (int ee=0;ee<NUM_ELEMENTS;ee++)
		if (Gal[p].EjectedMass_elements[ee] < 0.0)
			printf("IN SN_feedback() 2: Gal[%i].EjectedMass_elements[%i] = %e\n", p, ee, Gal[p].EjectedMass_elements[ee]);

  // Determine ejection (for FeedbackEjectionModel 0 we have the dependence on Vmax) Guo2010 - eq 22
  //EJECT
  if(FeedbackEagleScaling == 1)
     {
       //ReScaled_EnergySNcode=EnergySNcode*EAGLE2015_rescale_of_EnergySN (ColdGas, MetalsColdGas, Radius_low, Radius_high);
       //ReScaled_EnergySNcode=EnergySNcode/pow((1 + ZZ[Halo[Gal[p].HaloNr].SnapNum]),4.0);
       double SigmaGas;
       SigmaGas = Gal[p].ColdGas / (M_PI*(Gal[p].ColdGasRadius*Gal[p].ColdGasRadius-Radius_low*Gal[p].ColdGasRadius));
       // convert from 10^10 M_sun/h / (Mpc/h)^2 to (M_sun/pc^2)
       SigmaGas=SigmaGas*0.01*Hubble_h;

       ReScaled_EnergySNcode=EnergySNcode*SigmaGas;

       if(ReScaled_EnergySNcode>EnergySNcode)
          ReScaled_EnergySNcode=EnergySNcode;
     }
   else
     ReScaled_EnergySNcode=EnergySNcode;


  if (Gal[Gal[p].CentralGal].Type == 0)
    {
      EjectVmax=Gal[centralgal].Vmax;
      EjectVvir=Gal[centralgal].Vvir;// main halo Vvir
    }
  else
    {
      EjectVmax=Gal[Gal[p].CentralGal].InfallVmax;
      EjectVvir=Gal[Gal[p].CentralGal].Vvir; //central subhalo Vvir
    }

  if(FeedbackEjectionModel == 0)
    {
      ejected_mass = (FeedbackEjectionEfficiency* (EtaSNcode * ReScaled_EnergySNcode) * stars *
	  min(1./FeedbackEjectionEfficiency, .5+1/pow(EjectVmax/EjectPreVelocity,EjectSlope)) -
	  reheated_mass*EjectVvir*EjectVvir) /(EjectVvir*EjectVvir);
      //ejected_mass = EtaSNcode * ReScaled_EnergySNcode;
    }
  else if(FeedbackEjectionModel == 1)//the ejected material is assumed to have V_SN
    {
	  //ROB: Could replace this with (EnergySNII * SNIIRate)+(EnergySNIa * SNIaRate)+(EnergyAGB * AGBRate)? (15-05-22)
	  SN_Energy = .5 * stars * (EtaSNcode * ReScaled_EnergySNcode);
      Reheat_Energy = .5 * reheated_mass * EjectVvir * EjectVvir;

      ejected_mass = (SN_Energy - Reheat_Energy)/(0.5 * FeedbackEjectionEfficiency*(EtaSNcode * ReScaled_EnergySNcode));

      //if VSN^2<Vvir^2 nothing is ejected
      if(FeedbackEjectionEfficiency*(EtaSNcode * ReScaled_EnergySNcode)<EjectVvir*EjectVvir)
	  ejected_mass =0.0;
    }

  // Finished calculating mass exchanges, so just check that none are negative
  if (reheated_mass < 0.0) reheated_mass = 0.0;
  if (ejected_mass < 0.0) ejected_mass = 0.0;

  Gal[p].ReheatingRate += reheated_mass / (dt*STEPS); //ROB: storing ReheatingRate for outputting (31-03-20)
  Gal[p].EjectionRate += ejected_mass / (dt*STEPS); //ROB: storing EjectionRate for outputting (31-03-20)

	for (int ee=0;ee<NUM_ELEMENTS;ee++)
		if (Gal[p].EjectedMass_elements[ee] < 0.0)
			printf("IN SN_feedback() 3: Gal[%i].EjectedMass_elements[%i] = %e\n", p, ee, Gal[p].EjectedMass_elements[ee]);

  //ROB: Calculate the metal mass reheated:
#ifdef H2_AND_RINGS
  for(jj=0;jj<RNUM;jj++) {
      if(Gal[p].ColdGasRings[jj]>0.) {
    	  totmetalsRings=0.;
    	  for(mm=0;mm<NUM_METAL_CHANNELS;mm++) {
    		 if(Gal[p].MetalsColdGasRings[jj][mm]>0.)
    			 totmetalsRings += Gal[p].MetalsColdGasRings[jj][mm];
    	  }
    	  Gal[p].MetalsReheatingRate += (reheated_massr[jj] * (totmetalsRings/Gal[p].ColdGasRings[jj])) / (dt*STEPS); //reheated_massr[jj]/(Gal[p].ColdGasRings[jj]);
      }
  }
#endif

  /* Update For Feedback */
  /* update cold, hot, ejected gas fractions and respective metallicities
   * there are a number of changes introduced by Guo2010 concerning where
   * the gas ends up */

  //ejected_mass = 0.01*Gal[centralgal].HotGas;
  if (reheated_mass + ejected_mass > 0.)
    {
#ifndef H2_AND_RINGS
      update_from_feedback(p, centralgal, reheated_mass, ejected_mass);
#else
      update_from_feedback(p, centralgal, reheated_mass, ejected_mass,  reheated_massr);
#endif
    }

	for (int ee=0;ee<NUM_ELEMENTS;ee++)
		if (Gal[p].EjectedMass_elements[ee] < 0.0)
			printf("IN SN_feedback() 4: Gal[%i].EjectedMass_elements[%i] = %e\n", p, ee, Gal[p].EjectedMass_elements[ee]);

  mass_checks(p,"model_starformation_and_feedback.c",__LINE__);

}







/*#ifdef H2_AND_RINGS
double compute_SN_reheat(int p, int centralgal, double stars, double ColdGas, int jj)
#else
double compute_SN_reheat(int p, int centralgal, double stars, double ColdGas)
#endif*/
double compute_SN_reheat(int p, int centralgal, double stars, double ColdGas, double MetalsColdGas, double Radius_low, double Radius_high)
{
  double reheated_mass=0.;
  double ReScaled_EnergySNcode;
  double MergeCentralVvir=0.;

  /* In Guo2010 type 1s can eject, reincorporate gas and get gas from their
   * own satellites (is not sent to the type 0 galaxy as in Delucia2007),
   * for gas flow computations:
   * If satellite is inside Rvir of main halo, Vvir of main halo used
   * If it is outside, the Vvir of its central subhalo is used. */

  if(FeedbackEagleScaling == 1)
    {
      //ReScaled_EnergySNcode=EnergySNcode*EAGLE2015_rescale_of_EnergySN (ColdGas, MetalsColdGas, Radius_low, Radius_high);
      //ReScaled_EnergySNcode=EnergySNcode/pow((1 + ZZ[Halo[Gal[p].HaloNr].SnapNum]),4.0);

      double SigmaGas;
#ifdef H2_AND_RINGS
      SigmaGas = ColdGas/(M_PI *(Radius_high*Radius_high-Radius_low*Radius_low))/WARM_PHASE_FACTOR*Clumpingfactor;
#else
      SigmaGas = ColdGas / (M_PI*(Radius_high*Radius_high-Radius_low*Radius_low));
#endif
      // convert from 10^10 M_sun/h / (Mpc/h)^2 to (M_sun/pc^2)
      SigmaGas=SigmaGas*0.01*Hubble_h;

      ReScaled_EnergySNcode=EnergySNcode*SigmaGas;

      if(ReScaled_EnergySNcode>EnergySNcode)
         ReScaled_EnergySNcode=EnergySNcode;
    }
  else
    ReScaled_EnergySNcode=EnergySNcode;

  //REHEAT
  if(ColdGas>0.)
    {
  // Feedback depends on the circular velocity of the host halo
  // Guo2010 - eq 18 & 19
      if(FeedbackReheatingModel == 0)
	{
	  if (Gal[Gal[p].CentralGal].Type == 0)
	    reheated_mass = FeedbackReheatingEpsilon * stars * (.5+1./pow(Gal[Gal[p].CentralGal].Vmax/ReheatPreVelocity,ReheatSlope));
	  else
	    reheated_mass = FeedbackReheatingEpsilon * stars * (.5+1./pow(Gal[Gal[p].CentralGal].InfallVmax/ReheatPreVelocity,ReheatSlope));

	  if(FeedbackReheatingDeansityScaling==1)
	    {
	      double SigmaGas;
#ifdef H2_AND_RINGS
	      SigmaGas = ColdGas/(M_PI *(Radius_high*Radius_high-Radius_low*Radius_low))/WARM_PHASE_FACTOR*Clumpingfactor;
#else
	      SigmaGas = ColdGas / (M_PI*(Radius_high*Radius_high-Radius_low*Radius_low));
#endif
	      // convert from 10^10 M_sun/h / (Mpc/h)^2 to (M_sun/pc^2)
	      SigmaGas=SigmaGas*0.01*Hubble_h;

	      if(SigmaGas>0.)
		reheated_mass*=pow(SigmaGas*0.05,2);
	     // reheated_mass/=pow((1 + ZZ[Halo[Gal[p].HaloNr].SnapNum]),5.0); - much lower reheating, much higher metals
	     // reheated_mass*=pow((1 + ZZ[Halo[Gal[p].HaloNr].SnapNum]),2.0); // much higher reheating, much lower stellar metals, gas metals high-z unaffected
	    }

	  if (reheated_mass * Gal[Gal[p].CentralGal].Vvir * Gal[Gal[p].CentralGal].Vvir > stars * (EtaSNcode * ReScaled_EnergySNcode))
	    reheated_mass = stars * (EtaSNcode * ReScaled_EnergySNcode) / (Gal[Gal[p].CentralGal].Vvir * Gal[Gal[p].CentralGal].Vvir);
	}
    }
  else
    reheated_mass=0.;

  if(reheated_mass > ColdGas)
    reheated_mass = ColdGas;

 /*
      double rd, ringtot;
    int jj;
     rd=Gal[p].ColdGasRadius;
    ringtot=1.-(1+RingRadius[RNUM-1]/rd)/exp(RingRadius[RNUM-1]/rd);
    reheated_massr[0]=(1-(1+RingRadius[0]/rd)/exp(RingRadius[0]/rd))/ringtot*reheated_mass;
    //aux_mass=reheated_massr[0];
    for(jj=1; jj<RNUM; jj++)
      {
      reheated_massr[jj]= ((1+RingRadius[jj-1]/rd)/exp(RingRadius[jj-1]/rd)-(1+RingRadius[jj]/rd)/exp(RingRadius[jj]/rd))/ringtot*reheated_mass;
      //aux_mass+=reheated_massr[jj];
      }*/

  return reheated_mass;



}

double EAGLE2015_rescale_of_EnergySN (double ColdGas, double MetalsColdGas, double Radius_low, double Radius_high)
{
double fth_min=0.3, fth_max=3.0, fth;
double metallicity_Z, z_solar=0.0127;
double nH_birth, nH_grams=1.6737236e-24, nH_0=0.67;
double n=2*log(10);

metallicity_Z=MetalsColdGas/ColdGas;
nH_birth=ColdGas*UnitMass_in_g/nH_grams/(4./3.*M_PI*UnitLength_in_cm*(Radius_high*Radius_high-Radius_low*Radius_low));

fth=fth_min+(fth_max-fth_min)/( 1 + pow(metallicity_Z/(0.1*z_solar),n) * pow(nH_birth/nH_0,-1.*n) );

return  fth;
}

/** @brief Updates cold, hot and external gas components due to SN
 *         reheating and ejection. */
#ifndef H2_AND_RINGS
void update_from_feedback(int p, int centralgal, double reheated_mass, double ejected_mass)
#else
void update_from_feedback(int p, int centralgal, double reheated_mass, double ejected_mass, double reheated_massr[])
#endif
{
  double dis=0.;
  double MassRemain=0.;
  double fraction;
  float fracSNFB;
  int merger_centre=0;
#ifdef H2_AND_RINGS
  double fractionRings[RNUM], tmpfractionRings[RNUM], MassRemainRings[RNUM];
  int jj;
#endif
#ifdef DETAILED_DUST
  int ee;
  double previoustime, newtime, deltaT;
  previoustime = NumToTime(Gal[p].SnapNum);
  newtime = NumToTime(Gal[p].SnapNum+1);
  deltaT = previoustime - newtime;
#endif

  if(Gal[p].ColdGas > 0.) {
    /* REHEAT if galaxy is a type 1 or a type 2 orbiting a type 1 with hot gas
     * being stripped, some of the reheated and ejected masses goes to the type
     * 0 and some stays in the type 1 */
#ifdef H2_AND_RINGS
    for(jj=0;jj<RNUM;jj++)
      if(Gal[p].ColdGasRings[jj]>0.)
    	  fractionRings[jj]=reheated_massr[jj]/(Gal[p].ColdGasRings[jj]);
      else
    	  fractionRings[jj]=0.;
#endif

    if(Gal[p].Type ==0)
      {
    	fracSNFB = ((float)reheated_mass)/((float)Gal[p].ColdGas);
#ifdef H2_AND_RINGS
#ifdef DETAILED_DUST
    	transfer_material_with_rings(p,"HotGas",p,"ColdGas",fractionRings,fractionRings,fractionRings,"model_starformation_and_feedback.c", __LINE__);
    	//Update dust destruction rate:
    	for(ee=0;ee<NUM_ELEMENTS;ee++) {
#ifdef DUST_HOTGAS
    		Gal[p].DustColdGasRates[4] += 0.0; //No dust is destroyed just by being moved to the HotGas when DUST_HOTGAS is on.
#else //DUST_HOTGAS
    		Gal[p].DustColdGasRates[4] += (Gal[p].DustColdGasDiff_elements[ee]*fracSNFB + Gal[p].DustColdGasClouds_elements[ee]*fracSNFB)/(deltaT * UnitTime_in_years);
#endif //DUST_HOTGAS
    	}
#else //DETAILED_DUST
		  for (int ee=0;ee<NUM_ELEMENTS;ee++)
				if (Gal[p].HotGas_elements[ee] < 0.0)
					printf("IN update_from_feedback() 1: Gal[%i].HotGas_elements[%i] = %e | Gal[%i].EjectedMass_elements[%i] = %e\n",
							p, ee, Gal[p].HotGas_elements[ee], p, ee, Gal[p].EjectedMass_elements[ee]);

    	transfer_material_with_rings(p,"HotGas",p,"ColdGas",fractionRings,"model_starformation_and_feedback.c", __LINE__);

		  for (int ee=0;ee<NUM_ELEMENTS;ee++)
				if (Gal[p].HotGas_elements[ee] < 0.0)
					printf("IN update_from_feedback() 2: Gal[%i].HotGas_elements[%i] = %e\n", p, ee, Gal[p].HotGas_elements[ee]);
	//transfer_material_with_rings(p,"ReheatedGas",p,"ColdGas",fractionRings,"model_starformation_and_feedback.c", __LINE__);
#endif //DETAILED_DUST

#else //H2_AND_RINGS
	//fracSNFB = ((float)reheated_mass)/((float)Gal[p].ColdGas);
#ifdef DETAILED_DUST
	//transfer_dust_to_hot(p,fraction_sn_feedback);
	transfer_material(p,"HotGas",p,"ColdGas",fracSNFB,fracSNFB,fracSNFB,"model_starformation_and_feedback.c", __LINE__);
	//Update dust destruction rate:
	for(ee=0;ee<NUM_ELEMENTS;ee++) {
#ifdef DUST_HOTGAS
    		Gal[p].DustColdGasRates[4] += 0.0; //No dust is destroyed just by being moved to the HotGas when DUST_HOTGAS is on.
#else //DUST_HOTGAS
		Gal[p].DustColdGasRates[4] += (Gal[p].DustColdGasDiff_elements[ee]*fracSNFB + Gal[p].DustColdGasClouds_elements[ee]*fracSNFB)/(deltaT * UnitTime_in_years);
#endif //DUST_HOTGAS
	}
#else //DETAILED_DUST
	transfer_material(p,"HotGas",p,"ColdGas",fracSNFB,"model_starformation_and_feedback.c", __LINE__);
	//transfer_material(p,"ReheatedGas",p,"ColdGas",((float)reheated_mass)/((float)Gal[p].ColdGas),"model_starformation_and_feedback.c", __LINE__);
#endif //DETAILED_DUST
#endif //H2_AND_RINGS
      }

    //For satellite galaxies compute how much gas stays in the galaxy and how much goes to central companion
    else if(Gal[p].Type <3)
      {
	if(Gal[p].Type ==1)
	  merger_centre=centralgal;
	else if(Gal[p].Type ==2)
	  merger_centre=Gal[p].CentralGal;

	if(HotGasOnType2Galaxies==0)
	  //if no hot gas in type 2's, share gas between  0 and 1
	  dis=separation_gal(centralgal,Gal[p].CentralGal)/(1+ZZ[Halo[Gal[centralgal].HaloNr].SnapNum]);
	else if(HotGasOnType2Galaxies==1)
	  //if hot gas in type 2's, share gas between itself and mergercentre
	  dis=separation_gal(merger_centre,p)/(1+ZZ[Halo[Gal[centralgal].HaloNr].SnapNum]);

	//compute share of reheated mass
	if ((dis<Gal[centralgal].Rvir && Gal[Gal[p].CentralGal].Type == 1 &&  HotGasOnType2Galaxies==0) ||
	    (dis<Gal[merger_centre].Rvir && HotGasOnType2Galaxies==1))
	  {
	    //mass that remains on type1 (the rest goes to type 0)
	    //for reheat - MassRemain, for eject - ejected_mass
	    MassRemain=reheated_mass*Gal[p].HotRadius/Gal[p].Rvir;
	    ejected_mass = ejected_mass*Gal[p].HotRadius/Gal[p].Rvir;
	    if (MassRemain > reheated_mass)
	      MassRemain = reheated_mass;
#ifdef H2_AND_RINGS
	    for(jj=0;jj<RNUM;jj++)
	      {
		MassRemainRings[jj]=reheated_massr[jj]*Gal[p].HotRadius/Gal[p].Rvir;
		if (MassRemainRings[jj] > reheated_massr[jj])
		  MassRemainRings[jj] = reheated_massr[jj];
	      }
#endif
	  }
	else
	  {
	    MassRemain=reheated_mass;
#ifdef H2_AND_RINGS
	    for(jj=0;jj<RNUM;jj++)
	      MassRemainRings[jj] = reheated_massr[jj];
#endif
	  }

      //needed due to precision issues, since we first remove MassRemain and
      //then (reheated_mass-MassRemain) from the satellite into the type 0 and
      //type 1 the fraction might not add up on the second call since
      //Gal[p].ColdGas is a float and reheated_mass & MassRemain are doubles
      if((MassRemain + reheated_mass)>Gal[p].ColdGas)
	MassRemain=Gal[p].ColdGas-reheated_mass;

#ifdef H2_AND_RINGS
      for(jj=0;jj<RNUM;jj++)
	if((MassRemainRings[jj] + reheated_massr[jj])>Gal[p].ColdGasRings[jj])
	  MassRemainRings[jj] = Gal[p].ColdGasRings[jj]-reheated_massr[jj];
#endif

      mass_checks(p,"model_starformation_and_feedback.c",__LINE__);
      //transfer MassRemain

      if(reheated_mass>0.)
	{
      fracSNFB = MassRemain/Gal[p].ColdGas;
#ifdef H2_AND_RINGS
	  //for(jj=0;jj<RNUM;jj++)
	  //  tmpfractionRings[jj]=fractionRings[jj]*(MassRemain/reheated_mass);
	  for(jj=0;jj<RNUM;jj++)
	    if(Gal[p].ColdGasRings[jj]>0.)
	      tmpfractionRings[jj]=MassRemainRings[jj]/(Gal[p].ColdGasRings[jj]);
	    else
	      tmpfractionRings[jj]=0.;

#ifdef DETAILED_DUST
	  if(HotGasOnType2Galaxies==0)
	    {
	    //transfer to itself if type 1, merger centre if type 2
	    if(Gal[p].CentralGal==p)
	      transfer_material_with_rings(Gal[p].CentralGal,"HotGas",p,"ColdGas",tmpfractionRings,tmpfractionRings,tmpfractionRings,"model_starformation_and_feedback.c", __LINE__);
	    else
	      transfer_material_with_rings(Gal[p].CentralGal,"HotGas",p,"ColdGas",tmpfractionRings,tmpfractionRings,tmpfractionRings,"model_starformation_and_feedback.c", __LINE__);
	    }
	  else if(HotGasOnType2Galaxies==1) //tranfer to itself
	    transfer_material_with_rings(p,"HotGas",p,"ColdGas",tmpfractionRings,tmpfractionRings,tmpfractionRings,"model_starformation_and_feedback.c", __LINE__);

	  //Update dust destruction rate:
	  for(ee=0;ee<NUM_ELEMENTS;ee++) {
#ifdef DUST_HOTGAS
		  Gal[p].DustColdGasRates[4] += 0.0; //No dust is destroyed just by moving to the HotGas when DUST_HOTGAS is on.
#else //DUST_HOTGAS
	  	  Gal[p].DustColdGasRates[4] += (Gal[p].DustColdGasDiff_elements[ee]*fracSNFB + Gal[p].DustColdGasClouds_elements[ee]*fracSNFB)/(deltaT * UnitTime_in_years);
#endif //DUST_HOTGAS
	  }
#else //DETAILED_DUST
	  if(HotGasOnType2Galaxies==0)
	  	    {
	  	    //tranfer to itself if type 1, merger centre if type 2
	  	    if(Gal[p].CentralGal==p)
	  	      transfer_material_with_rings(Gal[p].CentralGal,"HotGas",p,"ColdGas", tmpfractionRings,"model_starformation_and_feedback.c", __LINE__);
	  	      //transfer_material_with_rings(Gal[p].CentralGal,"ReheatedGas",p,"ColdGas", tmpfractionRings,"model_starformation_and_feedback.c", __LINE__);
	  	    else
	  	      transfer_material_with_rings(Gal[p].CentralGal,"HotGas",p,"ColdGas", tmpfractionRings,"model_starformation_and_feedback.c", __LINE__);
	  	    }
	  else if(HotGasOnType2Galaxies==1) //transfer to itself
	  	    transfer_material_with_rings(p,"HotGas",p,"ColdGas",tmpfractionRings,"model_starformation_and_feedback.c", __LINE__);
	  	    //transfer_material_with_rings(p,"ReheatedGas",p,"ColdGas",tmpfractionRings,"model_starformation_and_feedback.c", __LINE__);
#endif //DETAILED_DUST
#else //H2_AND_RINGS
	  //fracSNFB = MassRemain/Gal[p].ColdGas;
#ifdef DETAILED_DUST	  
	  if(HotGasOnType2Galaxies==0) {
	  	    //tranfer to itself if type 1, merger centre if type 2. ROB: Why bother with this if/else if the identical mass transfer is executed? (26-11-21)
	  	    if(Gal[p].CentralGal==p)
	  			transfer_material(Gal[p].CentralGal,"HotGas",p,"ColdGas",fracSNFB,fracSNFB,fracSNFB,"model_starformation_and_feedback.c", __LINE__);
	  	     	// transfer_material(Gal[p].CentralGal,"ReheatedGas",p,"ColdGas", MassRemain/Gal[p].ColdGas,"model_starformation_and_feedback.c", __LINE__);
	  	    else
	  			transfer_material(Gal[p].CentralGal,"HotGas",p,"ColdGas",fracSNFB,fracSNFB,fracSNFB,"model_starformation_and_feedback.c", __LINE__);
	  }
	  else if(HotGasOnType2Galaxies==1)  //transfer to itself
	  	    transfer_material(p,"HotGas",p,"ColdGas",fracSNFB,fracSNFB,fracSNFB,"model_starformation_and_feedback.c", __LINE__);
	  //Update dust destruction rate:
	  for(ee=0;ee<NUM_ELEMENTS;ee++) {
#ifdef DUST_HOTGAS
		  Gal[p].DustColdGasRates[4] += 0.0; //No dust is destroyed just by moving to the HotGas when DUST_HOTGAS is on.
#else //DUST_HOTGAS
		  Gal[p].DustColdGasRates[4] += (Gal[p].DustColdGasDiff_elements[ee]*fracSNFB + Gal[p].DustColdGasClouds_elements[ee]*fracSNFB)/(deltaT * UnitTime_in_years);
#endif //DUST_HOTGAS
	  }
#else //DETAILED_DUST
	  if(HotGasOnType2Galaxies==0) {
	    //tranfer to itself if type 1, merger centre if type 2
	      if(Gal[p].CentralGal==p)
			transfer_material(Gal[p].CentralGal,"HotGas",p,"ColdGas", fracSNFB,"model_starformation_and_feedback.c", __LINE__);
	     	// transfer_material(Gal[p].CentralGal,"ReheatedGas",p,"ColdGas", fracSNFB,"model_starformation_and_feedback.c", __LINE__);
	      else
			transfer_material(Gal[p].CentralGal,"HotGas",p,"ColdGas", fracSNFB,"model_starformation_and_feedback.c", __LINE__);
	    }
	  else if(HotGasOnType2Galaxies==1)  //transfer to itself
	    transfer_material(p,"HotGas",p,"ColdGas", fracSNFB,"model_starformation_and_feedback.c", __LINE__);
	    //transfer_material(p,"ReheatedGas",p,"ColdGas", fracSNFB,"model_starformation_and_feedback.c", __LINE__);
#endif //DETAILED_DUST
#endif
	}

      mass_checks(p,"model_starformation_and_feedback.c",__LINE__);

      //transfer reheated_mass-MassRemain from galaxy to the type 0
      if (reheated_mass > MassRemain)
	if(Gal[p].ColdGas > 0.)
	  {
	    //if the reheat to itself left cold gas below limit do not reheat
	    //to central

		  ///with rings ColdGas is a double and using (float) might cause
		  //(float)(reheated_mass-MassRemain)/Gal[p].ColdGas to be >1
		  fraction=(float)(reheated_mass-MassRemain)/Gal[p].ColdGas;
		  if (fraction>1.+PRECISION_LIMIT) {
			  printf("(reheated_mass-MassRemain)/Gal[p].ColdGas=%f\n",fraction);
			  terminate("reheated too much gas");
		  }
		  fraction = min(fraction,1.);
#ifdef H2_AND_RINGS
	    //cannot use tmpfractionRings defined from fractionRings since
	    // Gal[p].ColdGasRings has changed from MassRemain above
	    for(jj=0;jj<RNUM;jj++)
	      if(Gal[p].ColdGasRings[jj]>0.)
		fractionRings[jj]=(reheated_massr[jj]-MassRemainRings[jj])/Gal[p].ColdGasRings[jj];
		//fractionRings[jj]=(reheated_massr[jj]/Gal[p].ColdGasRings[jj])*	((reheated_mass-MassRemain)/reheated_mass);
	      else
		fractionRings[jj]=0.;

#ifdef DETAILED_DUST
	  if(HotGasOnType2Galaxies==0)	 //transfer to type 0
		  transfer_material_with_rings(centralgal,"HotGas",p,"ColdGas",fractionRings,fractionRings,fractionRings,"model_starformation_and_feedback.c", __LINE__);
	  else if(HotGasOnType2Galaxies==1)//transfer to merger centre
	      transfer_material_with_rings(merger_centre,"HotGas",p,"ColdGas",fractionRings,fractionRings,fractionRings,"model_starformation_and_feedback.c", __LINE__);
	  //Update dust destruction rate:
	  for(ee=0;ee<NUM_ELEMENTS;ee++) {
#ifdef DUST_HOTGAS
		  Gal[p].DustColdGasRates[4] += 0.0; //No dust is destroyed just by moving to the HotGas when DUST_HOTGAS is on.
#else //DUST_HOTGAS
		  Gal[p].DustColdGasRates[4] += (Gal[p].DustColdGasDiff_elements[ee]*fraction + Gal[p].DustColdGasClouds_elements[ee]*fraction)/(deltaT * UnitTime_in_years);
#endif //DUST_HOTGAS
	  }
#else //DETAILED_DUST

	  if(HotGasOnType2Galaxies==0)	 //transfer to type 0
		  transfer_material_with_rings(centralgal,"HotGas",p,"ColdGas",fractionRings,"model_starformation_and_feedback.c", __LINE__);
	  else if(HotGasOnType2Galaxies==1)//transfer to merger centre
	      transfer_material_with_rings(merger_centre,"HotGas",p,"ColdGas",fractionRings,"model_starformation_and_feedback.c", __LINE__);
#endif //DETAILED_DUST
#else //H2_AND_RINGS
#ifdef DETAILED_DUST
	  if(HotGasOnType2Galaxies==0) //transfer to type 0
	    transfer_material(centralgal,"HotGas",p,"ColdGas",fraction,fraction,fraction,"model_starformation_and_feedback.c", __LINE__);
	  else if(HotGasOnType2Galaxies==1) //transfer to merger centre
	    transfer_material(merger_centre,"HotGas",p,"ColdGas",fraction,fraction,fraction,"model_starformation_and_feedback.c", __LINE__);
	  //Update dust destruction rate:
	  for(ee=0;ee<NUM_ELEMENTS;ee++) {
#ifdef DUST_HOTGAS
		  Gal[p].DustColdGasRates[4] += 0.0; //No dust is destroyed just by moving to the HotGas when DUST_HOTGAS is on.
#else //DUST_HOTGAS
		  Gal[p].DustColdGasRates[4] += (Gal[p].DustColdGasDiff_elements[ee]*fraction + Gal[p].DustColdGasClouds_elements[ee]*fraction)/(deltaT * UnitTime_in_years);
#endif //DUST_HOTGAS
	  }
#else
	  if(HotGasOnType2Galaxies==0) //transfer to type 0
	  	    transfer_material(centralgal,"HotGas",p,"ColdGas",fraction,"model_starformation_and_feedback.c", __LINE__);
	  	  else if(HotGasOnType2Galaxies==1) //transfer to merger centre
	  	    transfer_material(merger_centre,"HotGas",p,"ColdGas",fraction,"model_starformation_and_feedback.c", __LINE__);
#endif
#endif //H2_AND_RINGS
	  }

      //} //Gal[p].Type !=2
      mass_checks(p,"model_starformation_and_feedback.c",__LINE__);
    }//types

  }//if(Gal[p].ColdGas > 0.)

  mass_checks(p,"model_starformation_and_feedback.c",__LINE__);

   //**********
   //DO EJECTION OF GAS
  if ( (Gal[Gal[p].CentralGal].HotGas > 0. && HotGasOnType2Galaxies==0) ||
       (Gal[p].HotGas > 0. && HotGasOnType2Galaxies==1) ) {
    
    if(HotGasOnType2Galaxies==0) {
      if (ejected_mass > Gal[Gal[p].CentralGal].HotGas)
	//either eject own gas or merger_centre gas for ttype 2's
	ejected_mass = Gal[Gal[p].CentralGal].HotGas;
      fraction=ejected_mass/Gal[Gal[p].CentralGal].HotGas;
      
    }
    else if(HotGasOnType2Galaxies==1) {
      if (ejected_mass > Gal[p].HotGas && HotGasOnType2Galaxies==1)
	ejected_mass = Gal[p].HotGas;  //always eject own gas
      fraction=ejected_mass/Gal[p].HotGas;
    }

    //**********
    // If type 1, or type 2 orbiting type 1 near type 0
    if (Gal[Gal[p].CentralGal].Type == 1) {
		if (FateOfSatellitesGas == 0) {
#ifdef DETAILED_DUST
			if(HotGasOnType2Galaxies==0)
				transfer_material(Gal[p].CentralGal,"EjectedMass",Gal[p].CentralGal,"HotGas",fraction,0.0,fraction,"model_starformation_and_feedback.c", __LINE__); //Diff_fraction is set to fraction so that some dust in the HotGas is also destroyed, when DUST_HOTGAS is on.
			else if(HotGasOnType2Galaxies==1)
				transfer_material(Gal[p].CentralGal,"EjectedMass",p,"HotGas", fraction,0.0,fraction,"model_starformation_and_feedback.c", __LINE__); //Diff_fraction is set to fraction so that some dust in the HotGas is also destroyed, when DUST_HOTGAS is on.
#ifdef FULL_DUST_RATES
#ifdef DUST_HOTGAS
//Change for EjectaDust:
#ifndef DUST_EJECTEDMASS
			for(ee=0;ee<NUM_ELEMENTS;ee++)
				Gal[p].DustHotGasRates[3] += (Gal[p].DustHotGas_elements[ee]*fraction)/(deltaT * UnitTime_in_years); //If DUST_EJECTEDMASS is OFF, all dust transferred to the Ejecta component is assumed to be fully and instantaneously destroyed.
#endif //DUST_EJECTEDMASS
#endif //DUST_HOTGAS
#endif //FULL_DUST_RATES
#else //DETAILED_DUST
			  for (int ee=0;ee<NUM_ELEMENTS;ee++)
					if (Gal[Gal[p].CentralGal].HotGas_elements[ee] < 0.0)
						printf("IN update_from_feedback() 3: Gal[%i].HotGas_elements[%i] = %e\n", Gal[p].CentralGal, ee, Gal[Gal[p].CentralGal].HotGas_elements[ee]);


			if(HotGasOnType2Galaxies==0)
				transfer_material(Gal[p].CentralGal,"EjectedMass",Gal[p].CentralGal,"HotGas",fraction,"model_starformation_and_feedback.c", __LINE__);
			else if(HotGasOnType2Galaxies==1)
				transfer_material(Gal[p].CentralGal,"EjectedMass",p,"HotGas", fraction,"model_starformation_and_feedback.c", __LINE__);

			  for (int ee=0;ee<NUM_ELEMENTS;ee++)
					if (Gal[Gal[p].CentralGal].HotGas_elements[ee] < 0.0)
						printf("IN update_from_feedback() 4: Gal[%i].HotGas_elements[%i] = %e\n", Gal[p].CentralGal, ee, Gal[Gal[p].CentralGal].HotGas_elements[ee]);
#endif //DETAILED_DUST
		}
		else if (FateOfSatellitesGas == 1) {
#ifdef DETAILED_DUST
			if (dis < Gal[centralgal].Rvir)
				transfer_material(centralgal,"HotGas",Gal[p].CentralGal,"HotGas",fraction,0.0,fraction,"model_starformation_and_feedback.c", __LINE__); //Diff_fraction is set to fraction so that some dust in the HotGas of Gal[p].CentralGal is transferred to HotGas of centralgal, when DUST_HOTGAS is on.
			else {
				transfer_material(Gal[p].CentralGal,"EjectedMass",Gal[p].CentralGal,"HotGas",fraction,0.0,fraction,"model_starformation_and_feedback.c", __LINE__); //Diff_fraction is set to fraction so that some dust in the HotGas is also destroyed, when DUST_HOTGAS is on and DUST_EJECTEDMASS is off, or transferred to EjectedMass when both DUST_HOTGAS and DUST_EJECTEDMASS are on.
//Change for EjectaDust:
#ifdef FULL_DUST_RATES
#ifdef DUST_HOTGAS
#ifndef DUST_EJECTEDMASS
				for(ee=0;ee<NUM_ELEMENTS;ee++)
					Gal[p].DustHotGasRates[3] += (Gal[p].DustHotGas_elements[ee]*fraction)/(deltaT * UnitTime_in_years); //If DUST_EJECTEDMASS is OFF, all dust transferred to the Ejecta component is assumed to be fully and instantaneously destroyed.
#endif //DUST_EJECTEDMASS
#endif //DUST_HOTGAS
#endif //FULL_DUST_RATES
			}
#else
			if (dis < Gal[centralgal].Rvir)
				transfer_material(centralgal,"HotGas",Gal[p].CentralGal,"HotGas",fraction,"model_starformation_and_feedback.c", __LINE__);
			else
				transfer_material(Gal[p].CentralGal,"EjectedMass",Gal[p].CentralGal,"HotGas",fraction,"model_starformation_and_feedback.c", __LINE__);
#endif
		}
    }
    //**********
    // If galaxy type 0 or type 2 merging into type 0
    else {
#ifdef DETAILED_DUST
/*
      if(HotGasOnType2Galaxies==0)
    	  transfer_material(centralgal,"EjectedMass",Gal[p].CentralGal,"HotGas",fraction,0.0,0.0,"model_starformation_and_feedback.c", __LINE__);
      else if(HotGasOnType2Galaxies==1)
    	  transfer_material(centralgal,"EjectedMass",p,"HotGas",fraction,0.0,0.0,"model_starformation_and_feedback.c", __LINE__);
*/
//Change for EjectaDust:
      if(HotGasOnType2Galaxies==0)
          	  transfer_material(centralgal,"EjectedMass",Gal[p].CentralGal,"HotGas",fraction,0.0,fraction,"model_starformation_and_feedback.c", __LINE__);
      else if(HotGasOnType2Galaxies==1)
          	  transfer_material(centralgal,"EjectedMass",p,"HotGas",fraction,0.0,fraction,"model_starformation_and_feedback.c", __LINE__);

//Change for EjectaDust:
#ifdef FULL_DUST_RATES
#ifdef DUST_HOTGAS
#ifndef DUST_EJECTEDMASS
			for(ee=0;ee<NUM_ELEMENTS;ee++)
				Gal[p].DustHotGasRates[3] += (Gal[p].DustHotGas_elements[ee]*fraction)/(deltaT * UnitTime_in_years); //If DUST_EJECTEDMASS is OFF, all dust transferred to the Ejecta component is assumed to be fully and instantaneously destroyed.
#endif //DUST_EJECTEDMASS
#endif //DUST_HOTGAS
#endif //FULL_DUST_RATES

#else //DETAILED_DUST
      if(HotGasOnType2Galaxies==0) {
          for (int ee=0;ee<NUM_ELEMENTS;ee++)
        		if (Gal[Gal[p].CentralGal].HotGas_elements[ee] < 0.0) {
        			printf("IN update_from_feedback() 5: Gal[%i].HotGas_elements[%i] = %e\n", Gal[p].CentralGal, ee, Gal[Gal[p].CentralGal].HotGas_elements[ee]);
        		}

    	  transfer_material(centralgal,"EjectedMass",Gal[p].CentralGal,"HotGas",fraction,"model_starformation_and_feedback.c", __LINE__);

          for (int ee=0;ee<NUM_ELEMENTS;ee++)
        		if (Gal[centralgal].EjectedMass_elements[ee] < 0.0) {
        			printf("IN update_from_feedback() 6: Gal[%i].EjectedMass_elements[%i] = %e\n", centralgal, ee, Gal[centralgal].EjectedMass_elements[ee]);
        			//printf("Gal[%i].HotGas_elements[%i] = %e | fraction = %e\n", Gal[p].CentralGal, ee, Gal[Gal[p].CentralGal].HotGas_elements[ee], fraction);
        		}
      }
      else if(HotGasOnType2Galaxies==1)
          transfer_material(centralgal,"EjectedMass",p,"HotGas",fraction,"model_starformation_and_feedback.c", __LINE__);
#endif //DETAILED_DUST
    }

  }//(Gal[Gal[p].CentralGal].HotGas > 0.)

#ifdef H2_AND_RINGS
  update_h2fraction(p);
#endif

  mass_checks(p,"model_starformation_and_feedback.c",__LINE__);

}


//Age in Mpc/Km/s/h - code units
void update_massweightage(int p, double stars, double time)
{
  int outputbin;
  double age;

  for(outputbin = 0; outputbin < NOUT; outputbin++)
    {
      age = time - NumToTime(ListOutputSnaps[outputbin]);
#ifdef DETAILED_METALS_AND_MASS_RETURN
      Gal[p].MassWeightAge[outputbin] += age * stars;
#else
      Gal[p].MassWeightAge[outputbin] += age * stars * (1. - RecycleFraction);
#endif
    }
}


/** @brief Checks for disk stability using the
 *         Mo, Mao & White (1998) criteria as in Irodotou2018 */

void check_disk_instability_gas(int p, double dt)
{

  double Mcrit, fraction, unstable_mass, mass, BH_unstable_mass;
#ifdef H2_AND_RINGS
  double radius, vmax;
  double dmass, fractionRings[RNUM];
  int jj;
#endif

  mass = Gal[p].ColdGas;

#ifndef H2_AND_RINGS
  /* check stellar disk -> eq 34 Guo2010*/
  if (Gal[p].Type != 0)
    Mcrit = Gal[p].InfallVmax * Gal[p].InfallVmax * Gal[p].ColdGasRadius / G;
  else
    Mcrit = Gal[p].Vmax * Gal[p].Vmax * Gal[p].ColdGasRadius / G;
#else
   if (Gal[p].Type != 0)
     vmax=Gal[p].InfallVmax;
   else
     vmax=Gal[p].Vmax;

   Mcrit = vmax * vmax *  get_gas_disk_radius(p) / G;
#endif

   mass_checks(p,"model_starformation_and_feedback.c",__LINE__);

   unstable_mass = mass - Mcrit;

   if(BHGrowthInDiskInstabilityModel == 1)
       if(unstable_mass > 0.0)
	 {
	   //mass to be transferred to the black hole
	   BH_unstable_mass = unstable_mass * BlackHoleGrowthRate  / (1.0 + pow2((BlackHoleCutoffVelocity / Gal[p].Vvir)));

	   if(BH_unstable_mass>unstable_mass)
	     BH_unstable_mass = unstable_mass;
	   unstable_mass -= BH_unstable_mass;

	   fraction = BH_unstable_mass/mass * BlackHoleGrowthRate  / (1.0 + pow2((BlackHoleCutoffVelocity / Gal[p].Vvir)));

#ifdef H2_AND_RINGS
	   for(jj=0;jj<RNUM;jj++)
	     fractionRings[jj]=0.;

	   dmass = BH_unstable_mass;
	   for(jj=0;jj<RNUM;jj++)
     	     {
	       //mass is transfered first from the inner rings
	       //until the necessary mass is achieved
	       if(dmass>Gal[p].ColdGasRings[jj])
		 {
		   dmass-=Gal[p].ColdGasRings[jj];
		   fractionRings[jj]=1.;
		 }
	       else break;
     	     }

	   	    //ROB: Note: Because a "break" was used above, jj is stuck on the first ring to have ColdGasRings < dmass. Therefore, we can calculate the fractionRings for this ring as follows, to finish-off othe remaining dmass:
	   	    //check needed in case there is a ring with 0 mass in the middle
            if(Gal[p].ColdGasRings[jj]>0.)
              fractionRings[jj]=dmass/Gal[p].ColdGasRings[jj];
            else
              fractionRings[jj]=0.;

#ifdef DETAILED_DUST
            transfer_material_with_rings(p,"BlackHoleMass",p,"ColdGas",fractionRings,fractionRings,fractionRings,"model_starformation_and_feedback.c", __LINE__);
#else
            transfer_material_with_rings(p,"BlackHoleMass",p,"ColdGas",fractionRings,"model_starformation_and_feedback.c", __LINE__);
#endif //DETAILED_DUST
            mass_checks(p,"model_starformation_and_feedback.c",__LINE__);
#else //H2_AND_RINGS
#ifdef DETAILED_DUST
            transfer_material(p,"BlackHoleMass",p,"ColdGas",fraction,fraction,fraction,"model_starformation_and_feedback.c", __LINE__);
#else
            transfer_material(p,"BlackHoleMass",p,"ColdGas",fraction,"model_starformation_and_feedback.c", __LINE__);
#endif
 #endif
            Gal[p].QuasarAccretionRate += BH_unstable_mass/mass*Gal[p].ColdGas / (dt*STEPS);

	 }// if(unstable_mass > 0.0)


   /* add excess stars to the bulge */
   if(unstable_mass > 0.0)
     {
#ifdef H2_AND_RINGS
       for(jj=0;jj<RNUM;jj++)
	 fractionRings[jj]=0.;

       dmass=unstable_mass;

       for(jj=0;jj<RNUM;jj++)
	   {
	     //mass is transfered first from the inner rings
	     //until the necessary mass is achieved
	     if(dmass>Gal[p].ColdGasRings[jj])
	       {
		 dmass-=Gal[p].ColdGasRings[jj];
		 fractionRings[jj]=1.;
	       }
	     else break;
	   }

       //ROB: Note: Because a "break" was used above, jj is stuck on the first ring to have ColdGasRings < dmass. Therefore, we can calculate the fractionRings for this ring as follows, to finish-off the remaining dmass:
       //check needed in case there is a ring with 0 mass in the middle
       if(Gal[p].ColdGasRings[jj]>0.)
	 fractionRings[jj]=dmass/Gal[p].ColdGasRings[jj];
       else
	 fractionRings[jj]=0.;

#ifdef DETAILED_DUST
       transfer_material_with_rings(p,"DiskMass",p,"ColdGas",fractionRings,fractionRings,fractionRings,"model_starformation_and_feedback.c", __LINE__);
#else
       transfer_material_with_rings(p,"DiskMass",p,"ColdGas",fractionRings,"model_starformation_and_feedback.c", __LINE__);
#endif //DETAILED_DUST
       mass_checks(p,"model_starformation_and_feedback.c",__LINE__);
#else //H2_AND_RINGS
#ifdef DETAILED_DUST
       transfer_material(p,"DiskMass",p,"ColdGas",unstable_mass/mass,unstable_mass/mass,unstable_mass/mass, "model_starformation_and_feedback.c", __LINE__);
#else
       transfer_material(p,"DiskMass",p,"ColdGas",unstable_mass/mass, "model_starformation_and_feedback.c", __LINE__);
#endif
#endif

     }// if(unstable_mass > 0.0)
   mass_checks(p,"model_starformation_and_feedback.c",__LINE__);
} //end check_disk_instability_gas


/** @brief Checks for disk stability using the
 *         Mo, Mao & White (1998) criteria */

void check_disk_instability(int p, double dt)
{

  double Mcrit, fraction, stars, diskmass;
#ifdef H2_AND_RINGS
  double rstar, vmax;
  int j;
#endif
/** @brief Calculates the stability of the stellar disk as discussed
 *         in Mo, Mao & White (1998). For unstable stars, the required
 *         amount is transfered to the bulge to make the disk stable again.
 *         Mass, metals and luminosities updated. After Guo2010 the bulge
 *         size is followed and needs to be updated.
 *         Eq 34 & 35 in Guo2010 are used. */

  diskmass = Gal[p].DiskMass;

#ifndef H2_AND_RINGS
  /* check stellar disk -> eq 34 Guo2010*/
  if (Gal[p].Type != 0)
    Mcrit = Gal[p].InfallVmax * Gal[p].InfallVmax * Gal[p].DiskRadius / G;
  else
    Mcrit = Gal[p].Vmax * Gal[p].Vmax * Gal[p].DiskRadius / G;
#else

   if (Gal[p].Type != 0)
     vmax=Gal[p].InfallVmax;
   else
     vmax=Gal[p].Vmax;

   Mcrit = vmax * vmax * get_stellar_disk_radius(p) / G;
#endif

   mass_checks(p,"model_starformation_and_feedback.c",__LINE__);

   stars = diskmass - Mcrit;
   fraction = stars / diskmass;

   /* add excess stars to the bulge */
   if(stars > 0.0)
     {
       /* to calculate the bulge size */
       update_bulgesize_from_disk_instability(p,stars);


       //The bulge will be formed in the same place as the disk was, so the disk rings
       //are transferred directly into bulge rings
#ifdef H2_AND_RINGS
       double dmass, fractionRings[RNUM];
#ifdef DETAILED_DUST
       double fractionCloudsRings[RNUM], fractionDiffRings[RNUM];
#endif
       for(j=0;j<RNUM;j++) {
    	   fractionRings[j]=0.;
#ifdef DETAILED_DUST
    	   fractionCloudsRings[j]=0.;
    	   fractionDiffRings[j]=0.;
#endif
       }

       dmass=stars;
       j=0; //avoid non-definded j if dmass<1e-6
       // if(dmass>1.0e-6)
       for(j=0;j<RNUM;j++) {
    	   //mass is transfered first from the inner rings
    	   //until the necessary mass is achieved
    	   if(dmass>Gal[p].DiskMassRings[j]) {
    		   dmass-=Gal[p].DiskMassRings[j];
    		   fractionRings[j]=1.;
    	   }
    	   else break;
	   }

       //check needed in case there is a ring with 0 mass in the middle
       if(Gal[p].DiskMassRings[j]>0.)
    	   fractionRings[j]=dmass/Gal[p].DiskMassRings[j];
       else
    	   fractionRings[j]=0.;

#ifdef DETAILED_DUST
       transfer_material_with_rings(p,"BulgeMass",p,"DiskMass",fractionRings,fractionCloudsRings,fractionDiffRings,"model_starformation_and_feedback.c", __LINE__);
#else
       transfer_material_with_rings(p,"BulgeMass",p,"DiskMass",fractionRings,"model_starformation_and_feedback.c", __LINE__);
#endif //DETAILED_DUST
       mass_checks(p,"model_starformation_and_feedback.c",__LINE__);
#else
#ifdef DETAILED_DUST
       transfer_material(p,"BulgeMass",p,"DiskMass",fraction,0.0,0.0,"model_starformation_and_feedback.c", __LINE__);
#else
       transfer_material(p,"BulgeMass",p,"DiskMass",fraction, "model_starformation_and_feedback.c", __LINE__);
#endif
#endif

#ifdef BULGESIZE_DEBUG
       mass_checks(p,"model_starformation_and_feedback.c",__LINE__);
       if ((Gal[p].BulgeMass > TINY_MASS && Gal[p].BulgeSize < TINY_LENGTH)||
	   (Gal[p].BulgeMass < TINY_MASS && Gal[p].BulgeSize > TINY_LENGTH)) {
	   printf("BulgeMass=%g BulgeSize=%g\n",Gal[p].BulgeMass,Gal[p].BulgeSize);
	   terminate("bulgesize wrong in disk instablility\n");
       }
#endif

/*
       //burst of star formation from the instability, same as in mergers, with diskmass transferred in instability = mass of satellite
       // and total disk mass = mass of central
       double frac, mass_ratio = fraction;
       frac = collisional_starburst_recipe(mass_ratio, p, p, time, dt*STEPS);
       bulgesize_from_merger(mass_ratio, p, p, Gal[p].BulgeMass+Gal[p].DiskMass, Gal[p].BulgeMass, Gal[p].ColdGas,
			     Gal[p].DiskMass*fraction, 0, Gal[p].ColdGas*fraction, frac,
			     get_gas_disk_radius(p)/3., get_stellar_disk_radius(p)/3., get_gas_disk_radius(p)/3., get_stellar_disk_radius(p)/3.);

       if(mass_ratio > ThreshMajorMerger)
          make_bulge_from_burst(p);*/


     }// if(stars > 0.0)
   mass_checks(p,"model_starformation_and_feedback.c",__LINE__);
} //end check_disk_instability


/** @brief Introduced in Guo2010 to track the change in size of bulges
 *         after their growth due to disk instabilities. */

void update_bulgesize_from_disk_instability(int p, double stars)
{      
  double bulgesize, diskmass, fint, massfrac;
  int  j;


/** @brief Updates bulge from disk instability -> stars represents the mass
  *        transfered to the bulge, which occupies a size in the bulge equal
  *        to that occupied in the disk. */


   // alpha_inter=2.0/C=0.5 (alpha larger than in mergers since
   // the existing and newly formed bulges are concentric)
   fint=4.0;

   // disk instabilities are assumed to NOT remove angular momentum.
   // Since the mass of the disk changes, the spin is changed by the same
   // amount to keep angular momentum constant


  diskmass=Gal[p].DiskMass;
  massfrac=stars/diskmass;
  for (j = 0; j <3 ; j++)
    {
      if(massfrac==1) //everything transferred to the bulge
	Gal[p].DiskSpin[j]=0;
      else
	Gal[p].DiskSpin[j]=Gal[p].DiskSpin[j]/(1-massfrac);
    }
  if (DiskRadiusModel == 0)
    Gal[p].DiskRadius = get_stellar_disk_radius(p);



#ifndef H2_AND_RINGS

#ifdef BULGESIZE_DEBUG
  double orisize;
  orisize=Gal[p].BulgeSize;
#endif




  //GET BULGE SIZE - Eq. 35 in Guo2010
  /* size of newly formed bulge, which consists of the stellar mass transfered
   * from the disk. This is calculated using bulge_from_disk which receives
   * Delta_M/DiskMass and returns Rb/Rd. From eq 35 and since
   * DiskMass=2PISigma(Rd)^2 we see that
   * Delta_M/DiskMass=1-(1+Rb/Rd)*exp(-Rb/Rd), so function bulge_from_disk
   * avoids calculating the slow "ln" function */
  bulgesize=bulge_from_disk(stars/diskmass)*Gal[p].DiskRadius/3.;
  if(Gal[p].BulgeMass < TINY_MASS) {
      /* if previous Bulge Mass = 0
       * -> bulge size is given directly from newly formed bulge */
      Gal[p].BulgeSize=bulgesize;
  } else {
      /* combine the old with newly formed bulge and calculate the
       * bulge size assuming energy conservation as for mergers but
       * using alpha=2. - eq 33 */
      Gal[p].BulgeSize=(Gal[p].BulgeMass+stars)*(Gal[p].BulgeMass+stars) /
	  (Gal[p].BulgeMass*Gal[p].BulgeMass/Gal[p].BulgeSize + stars*stars/bulgesize + fint*Gal[p].BulgeMass*stars/(Gal[p].BulgeSize+bulgesize));
  }

  /* Added by PAT to see if the cause of low bulgesizes could be propagation of
   * tightly-bound bulges.  These originate in unfeasibly small disks. */
#define BULGESIZE_MIN 1e-4
  Gal[p].BulgeSize=max(Gal[p].BulgeSize,BULGESIZE_MIN);

#ifdef BULGESIZE_DEBUG
  if((Gal[p].BulgeMass + stars > TINY_MASS && Gal[p].BulgeSize < TINY_LENGTH)
     || (Gal[p].BulgeMass + stars < TINY_MASS && Gal[p].BulgeSize > TINY_LENGTH)) {
      printf("Original DiskMass=%e, DiskSize=%e\nOriginal BulgeMass=%e, BulgeSize=%e\nTransferred stars=%e, bulgesize=%e\nFinal BulgeMass=%e, BulgeSize=%e\n",
	     Gal[p].DiskMass, Gal[p].DiskRadius, Gal[p].BulgeMass, orisize, stars, bulgesize, Gal[p].BulgeMass+stars, Gal[p].BulgeSize);
      terminate("bulgesize or mass wrong in disk instablility");
    }
#endif
  
#else //H2_AND_RINGS

  /*size of new formed bulge, which consist of the stellar mass trasfered from the disk*/
  /*combine the old bulge with the new materials and caculate the bulge size assuming energy conservation */
  diskmass=stars;
  //j=0;
 // if(diskmass>1.0e-6)
  //  {
      for(j=0;j<RNUM;j++)
  	{
	  //mass is transfered first from the inner rings first
	  //until the necessary mass is achieved
	  if(diskmass>Gal[p].DiskMassRings[j])
	      diskmass-=Gal[p].DiskMassRings[j];
	  else break;
  	}
      if(j==RNUM)
	bulgesize=RingRadius[RNUM-1];
      else
	{
	  if(j==0)
	    bulgesize=diskmass/Gal[p].DiskMassRings[j]*RingRadius[j];
	  else
	    bulgesize=diskmass/Gal[p].DiskMassRings[j]*RingRadius[j]+(1-diskmass/Gal[p].DiskMassRings[j])*RingRadius[j-1];
	}
  //  }
  //else bulgesize=0.5*RingRadius[0];

  if(Gal[p].BulgeMass <1.e-9)
    Gal[p].BulgeSize=bulgesize;
  else
    Gal[p].BulgeSize=(Gal[p].BulgeMass+stars)*(Gal[p].BulgeMass+stars)/
    (Gal[p].BulgeMass*Gal[p].BulgeMass/Gal[p].BulgeSize+stars*stars/bulgesize+fint*Gal[p].BulgeMass*stars/(Gal[p].BulgeSize+bulgesize));

    /*if ((Gal[p].BulgeMass+stars > 1.e-8 && Gal[p].BulgeSize == 0.0)||(Gal[p].BulgeMass+stars == 0 && Gal[p].BulgeSize >1.e-8))
      {
        printf("bulgesize wrong in disk instablility. Diskmass %f, bulgemass %f, bulgesize %f, coldgas %f, masstransfer %f transsize %f\n",
	       Gal[p].DiskMass, Gal[p].BulgeMass, Gal[p].BulgeSize, Gal[p].ColdGas, stars, bulgesize);
        exit(0);
      }*/

#endif //H2_AND_RINGS

}


/** @brief Calculates the size of the disk that contains the
 *         mass transfered to the bulge. */
double bulge_from_disk(double frac)
{
  double x1,x2,x0,value;
/** @brief Calculates the size of the disk that contains the
 *         mass transfered to the bulge. The bulge is assumed
 *         to form with the same size. avoid doing "ln" from eq 35*/
  x1=0.0;
  x2=1.;
  while ((func_size(x2,frac) * func_size(x1,frac))>0) {
    x1=x2;
    x2=x2*2;
  }
  x0=x1+(x2-x1)/2.;
  value=func_size(x0,frac);
  if (value < 0) 
    value = -value;

  while(value>0.00001) {
    if(func_size(x0,frac)*func_size(x2,frac)>0)
      x2=x0;
    else
      x1=x0;
    x0=x1+(x2-x1)/2.;
    value=func_size(x0,frac);
    if (value < 0) 
      value = -value;
  }
    
  return x0;
}


double func_size(double x, double a)
{
  return  exp(-x)*(1+x)-(1-a);
}  

