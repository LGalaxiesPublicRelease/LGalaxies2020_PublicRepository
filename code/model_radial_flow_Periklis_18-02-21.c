#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "allvars.h"
#include "proto.h"

void gas_inflow(int p, double time)
{
  if(RadialFlowModel==0)
    {
      double r_in, r_out, gas_old[RNUM], newarea, r1, r2, frac, alpha, inflowfrac, rgas, Velocity, SurfaceDensity;
      int j, index, ii;
      double gasmetal_old[RNUM][NUM_METAL_CHANNELS];
#ifdef DETAILED_METALS_AND_MASS_RETURN
#ifdef INDIVIDUAL_ELEMENTS
      int kk;
      double gasElements_old[RNUM][NUM_ELEMENTS];
#endif
#endif

      if(Gal[p].ColdGas<1.0e-6)
	rgas=Gal[p].ColdGasRadius;
      else
	{
	  //rgas=0.5*RingRadius[0]*Gal[p].ColdGasRings[0];
	  //for(j=1;j<RNUM;j++) rgas+=(0.5*(RingRadius[j-1]+RingRadius[j])*Gal[p].ColdGasRings[j]);
	  //rgas=rgas/Gal[p].ColdGas/2.0;      //2.0=mean radius/scale length for exponential disk
	  rgas=Gal[p].ColdGasRadius;
	}
      //rgas=Gal[p].ColdGasRadius/3.0;
	 
      //inflowfrac=1.00;
      //alpha=1-1.34e3*time;

      
      for (j=0; j<RNUM; j++)
	{
	  inflowfrac = 1.00;
	  gas_old[j]=Gal[p].ColdGasRings[j]*inflowfrac;
	  for(ii=0;ii<NUM_METAL_CHANNELS;ii++)
	    gasmetal_old[j][ii]=(Gal[p].MetalsColdGasRings[j][ii] * inflowfrac);
#ifdef INDIVIDUAL_ELEMENTS
	  for (kk=0;kk<NUM_ELEMENTS;kk++)
	    gasElements_old[j][kk]=Gal[p].ColdGasRings_elements[j][kk]*inflowfrac;
#endif

	  Gal[p].ColdGasRings[j]*=(1-inflowfrac);
	  for(ii=0;ii<NUM_METAL_CHANNELS;ii++)
	    Gal[p].MetalsColdGasRings[j][ii]=(Gal[p].MetalsColdGasRings[j][ii] * (1.-inflowfrac));
#ifdef INDIVIDUAL_ELEMENTS
	  for (kk=0;kk<NUM_ELEMENTS;kk++)
	    Gal[p].ColdGasRings_elements[j][kk]=Gal[p].ColdGasRings_elements[j][kk]*(1.-inflowfrac);
#endif
	}

      //Msun/pc^2
      /*SurfaceDensity = Gal[p].ColdGas*1e10 / (Gal[p].ColdGasRadius*Gal[p].ColdGasRadius*1e12);
	if(Gal[p].DiskMass>0.)
	{
	//Velocity = pow(1.+SurfaceDensity,0.05);
	Velocity = GasInflowVel*pow(SurfaceDensity,0.1);
	//printf("vel=%e vel_before=%e\n",Velocity,GasInflowVel*pow(SurfaceDensity,0.1));
	}
	else
	Velocity = 0.;*/
      Velocity = GasInflowVel ;

      //printf("%e %e %e\n",Velocity,Gal[p].ColdGasRadius, Gal[p].ColdGas);
      for (j=0; j<RNUM; j++)
	{
	  alpha=1-Velocity*time/Hubble_h; //time unit: (Mpc/h)/(km/s)
	  r_out=RingRadius[j]*alpha;
	  if(j==0) r_in=0.0;
	  else r_in=RingRadius[j-1]*alpha;
	  if(r_in<1.0e-8) r_in=1.0e-8;

	  //constant inflow velocity
	  //vgas=3.0; // km/s
	  //r_out=RingRadius[j]-vgas*time; //time unit: (Mpc/h)/(km/s); radius unit: Mpc/h
	  //if(j==0) r_in=0.0;
	  //else r_in=RingRadius[j-1]-vgas*time;
	  //constant inflow velocity
	  if(r_out<=RingRadius[0])
	    {
	      index=0;
	      Gal[p].ColdGasRings[0]+=gas_old[j];
	      for(ii=0;ii<NUM_METAL_CHANNELS;ii++)
		Gal[p].MetalsColdGasRings[0][ii] += gasmetal_old[j][ii];
#ifdef INDIVIDUAL_ELEMENTS
	      for (kk=0;kk<NUM_ELEMENTS;kk++)
		Gal[p].ColdGasRings_elements[0][kk]+=gasElements_old[j][kk];
#endif
	    }
	  else
	    {
	      newarea=r_out*r_out-r_in*r_in;
	      //index=floor(log(1000*r_out)/log(1.2)+3); //The inverse of RingRadius[i]= pow(1.2,i-3)/1000
	      for(index=0;RingRadius[index]<=r_out&&index<RNUM-1;index++);	//ring number of the r_out
	      index--;
	      r1=RingRadius[index]; r2=r_out;
	      while(r1>=r_in)
		{
		  frac=(r2*r2-r1*r1)/newarea;
		  Gal[p].ColdGasRings[index+1]+=frac*gas_old[j];
		  for(ii=0;ii<NUM_METAL_CHANNELS;ii++)
		    Gal[p].MetalsColdGasRings[index+1][ii] += (gasmetal_old[j][ii]*frac);
#ifdef INDIVIDUAL_ELEMENTS
		  for (kk=0;kk<NUM_ELEMENTS;kk++)
		    Gal[p].ColdGasRings_elements[index+1][kk]+=(gasElements_old[j][kk]*frac);
#endif
		  index--;
		  r2=r1;
		  if(index>-1) r1=RingRadius[index];
		  else r1=0.0;
		}
	      r1=r_in;
	      frac=(r2*r2-r1*r1)/newarea;
	      Gal[p].ColdGasRings[index+1]+=frac*gas_old[j];
	      for(ii=0;ii<NUM_METAL_CHANNELS;ii++)
		Gal[p].MetalsColdGasRings[index+1][ii] += (gasmetal_old[j][ii] * frac);
#ifdef INDIVIDUAL_ELEMENTS
	      for (kk=0;kk<NUM_ELEMENTS;kk++)
		Gal[p].ColdGasRings_elements[index+1][kk]+=(gasElements_old[j][kk]*frac);
#endif
	    }
	}
      if (DiskRadiusModel == 0)
	Gal[p].ColdGasRadius=get_gas_disk_radius(p);



    }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  else{
    double r_in, r_out, gas_old[RNUM], rgas,inflowfrac;//newarea, r1, r2, frac, alpha, inflowfrac, rgas, Velocity, SurfaceDensity;
    double delta_mus[RNUM], deltas[RNUM], widths[RNUM], accretedMass[RNUM], accretedGasFraction[RNUM], r_ring_norm;
    double spread_rs[RNUM][7];
    double dt;
    double mtot,rsub[RNUM], mtot_metals[NUM_METAL_CHANNELS], mtot_elements[NUM_ELEMENTS];
    double r2,r3;
    double Mcool,rd,ringtot;
    int j, k, l ,index, ii,jj;
    double gasmetal_old[RNUM][NUM_METAL_CHANNELS];
#ifdef DETAILED_METALS_AND_MASS_RETURN
#ifdef INDIVIDUAL_ELEMENTS
    int kk;
    double gasElements_old[RNUM][NUM_ELEMENTS];
#endif
#endif


    const double fracs[]={0.02,0.14,0.34,0.34,0.14,0.02};

    dt=time*1e6/1.023; //time in Myr or rather Myr/h, if time in (Mpc/h)/(km/s)=(Mpc/h)/(1.023*pc/Myr)
    

    if(Gal[p].ColdGas==0)
      {
	rgas=0.;
	//printf("1,rgas=%.3f\n",rgas*1e6);
      }
    else if(Gal[p].ColdGas<1.0e-6)
      {
	rgas=RingRadius[0]/2.;
	//printf("2,rgas=%.3f\n",rgas*1e6);
      }
    else if(Gal[p].ColdGas>=1.0e-6)
      {
	if(Gal[p].ColdGasRings[0]>0.)
	  rgas=0.5*RingRadius[0]*Gal[p].ColdGasRings[0];
	for(jj=1;jj<RNUM;jj++)
	  if(Gal[p].ColdGasRings[jj]>0.)
	    rgas+=(0.5*(RingRadius[jj-1]+RingRadius[jj])*Gal[p].ColdGasRings[jj]);

	rgas=3.0*rgas/Gal[p].ColdGas/2.0;      //2.0=mean radius/scale length for exponential disk
	
	//printf("3,rgas=%.3f,m=%.5f,m0=%.5f\n",rgas*1e3,Gal[p].ColdGas,Gal[p].ColdGasRings[0],Gal[p].ColdGasRadius);
      }

    
    for (j=0; j<RNUM; j++)
      {
	inflowfrac=1.0;
	gas_old[j]=Gal[p].ColdGasRings[j]; //saving the gas masses at n in the new array
	for(ii=0;ii<NUM_METAL_CHANNELS;ii++)
	  gasmetal_old[j][ii]=Gal[p].MetalsColdGasRings[j][ii];
#ifdef INDIVIDUAL_ELEMENTS
	for (kk=0;kk<NUM_ELEMENTS;kk++)
	  gasElements_old[j][kk]=Gal[p].ColdGasRings_elements[j][kk];
#endif

	//Gal[p].ColdGasRings[j]*=(1-inflowfrac); //setting the variables to 0 to use at n+1?

	for(ii=0;ii<NUM_METAL_CHANNELS;ii++)
	  Gal[p].MetalsColdGasRings[j][ii]=(Gal[p].MetalsColdGasRings[j][ii] * (1.-inflowfrac));
#ifdef INDIVIDUAL_ELEMENTS
	for (kk=0;kk<NUM_ELEMENTS;kk++)
	  Gal[p].ColdGasRings_elements[j][kk]=Gal[p].ColdGasRings_elements[j][kk]*(1.-inflowfrac);
#endif
      }



    Mcool=Gal[p].CoolingGas;//total mass of cooled gas to the disc
    //printf("%.3f \n",Mcool*1e10);
    //if (Mcool>Gal[p].HotGas)
    // Mcool = Gal[p].HotGas;

    //need to distribute the mass in the rings, here I can change the profile
    //  add the fraction 1/STEPS of the total cooling gas to the cold disk 
    if(Mcool > 0.0)
      {

	rd=get_initial_disk_radius(Gal[p].HaloNr, p)/3.;
	
	ringtot=1.-(1+RingRadius[RNUM-1]/rd)/exp(RingRadius[RNUM-1]/rd); //total area up until the last ring radius

	accretedMass[0]=((1-(1+RingRadius[0]/rd)/exp(RingRadius[0]/rd))/ringtot)*Mcool;
	for(jj=1; jj<RNUM; jj++)
	  accretedMass[jj]= (((1+RingRadius[jj-1]/rd)/exp(RingRadius[jj-1]/rd)-(1+RingRadius[jj]/rd)/exp(RingRadius[jj]/rd))/ringtot)*Mcool;


	
	//ringtot=RingRadius[RNUM-1];
	//accretedMass[0]=(pow(RingRadius[0],3)/pow(ringtot,3))*Mcool;
	//for(jj=1; jj<RNUM; jj++)
	//  accretedMass[jj]=((pow(RingRadius[jj],3)-pow(RingRadius[jj-1],3))/pow(ringtot,3))*Mcool;
	
	
      }
    else
      {
	for(jj=0; jj<RNUM; jj++)
	  accretedMass[jj]=0.0;
      }
  

    for(j=0; j<RNUM; j++)
      {
	//r_ring_norm should be the centres of the rings
	// my recipes are for kpc, km/s, kpc^3/Gyr
	// in the model time is in (Mpc/h)/(km/s) and radii in kpc

	if(rgas>0. && j==0)
	  r_ring_norm=RingRadius[0]/2./(rgas);
	else if(rgas>0. && j>0)
	  r_ring_norm=(RingRadius[j]+RingRadius[j-1])/2./(rgas);
	else
	  r_ring_norm=1;

	//if(p==0)
	// printf("j=%d,rgas=%.3f,rring=%.3f,rnorm=%.3f \n",j,rgas*1e3,RingRadius[j]*1e3,r_ring_norm);


	if(gas_old[j]>0. && accretedMass[j]>=0.){
	  accretedGasFraction[j]=accretedMass[j]/gas_old[j]*60./dt; //unitless
	  //printf("1,%.3f \n",accretedGasFraction[j]);
	}
	else{
	  accretedGasFraction[j]=0.;
	  //printf("2,%.3f,\n",accretedGasFraction[j]);
	}

	if(gas_old[j]>0. && accretedGasFraction[j]<1. && r_ring_norm<1.)
	  {
	    deltas[j] = 35.8*pow(r_ring_norm,1.1) + 21.2*accretedGasFraction[j] - 2.5; //in kpc^3/Gyr
	    widths[j] = pow(deltas[j]*(dt*1e-3),1./3.); //in kpc, I hope
	    //printf("1,%.3f,%.3f,%.3f,%.3f \n",widths[j],deltas[j],accretedGasFraction[j],r_ring_norm);
	  }
	else
	  {
	    deltas[j]= 0;
	    widths[j]= 0;
	  }

	if(r_ring_norm<0.75 && gas_old[j]>0.)// && accretedGasFraction[j]<1.)// && accretedGasFraction[j]<1.)
	  delta_mus[j]=  -1.7 - 6.8*accretedGasFraction[j]; //in km/s
	else if (r_ring_norm>0.75 && r_ring_norm<1. && gas_old[j]>0.)// && accretedGasFraction[j]<1.) 
	  delta_mus[j]= -15.9*r_ring_norm - 6.8*accretedGasFraction[j] + 10.2; //in km/s
	else
	  delta_mus[j]= 0;
	
	//if(r_ring_norm <1. && deltas[j]>1 && delta_mus[j]>1. )
	printf("t=%.3f,d=%.3f,dm=%.3f,rnorm=%.3f,R=%.3f,facc=%.3f\n",dt,deltas[j],delta_mus[j],r_ring_norm,RingRadius[j]*1e3,accretedGasFraction[j]);


	for(k=0; k<7; k++)
	  {
	    //units with radii in kpc
	    spread_rs[j][k] = r_ring_norm*Gal[p].ColdGasRadius*1e3 + delta_mus[j]*dt*1e-3 + (k-3)*widths[j]; 

	  }
      }


    for(j=0; j<RNUM; j++)
      {
	Gal[p].ColdGasRings[j]+= 0;


	if(j==0)
	  {
	    r_in=0;
	    r_out=RingRadius[0];
	  }
	else
	  {
	    r_in=RingRadius[j-1];
	    r_out=RingRadius[j];      
	  }


	for(l=0; l<RNUM; l++)
	  {
	    mtot=0;
	    for(ii=0;ii<NUM_METAL_CHANNELS;ii++)
	      mtot_metals[ii]=0;
#ifdef INDIVIDUAL_ELEMENTS
	    for (kk=0;kk<NUM_ELEMENTS;kk++)
	      mtot_elements[kk]=0;
#endif

	    for(k=0; k<6; k++)
	    
	      {
		if((spread_rs[l][k] > r_out) || (spread_rs[l][k+1] < r_in))
		  {
		    mtot+=0;
		    for(ii=0;ii<NUM_METAL_CHANNELS;ii++)
		      mtot_metals[ii]+=0;
#ifdef INDIVIDUAL_ELEMENTS
		    for(kk=0;kk<NUM_ELEMENTS;kk++)
		      mtot_elements[kk]+=0;
#endif
		  }
		else
		  {
		    r2 = max(r_in,spread_rs[l][k]);
		    r3 = min(r_out,spread_rs[l][k+1]);
		  
		    mtot+= gas_old[l]*fracs[k]* fabs(r2*r2-r3*r3) /fabs(spread_rs[l][k+1]*spread_rs[l][k+1]-spread_rs[l][k]*spread_rs[l][k]);

		    for(ii=0;ii<NUM_METAL_CHANNELS;ii++)
		      mtot_metals[ii]+=gasmetal_old[l][ii]*fracs[k]* fabs(r2*r2-r3*r3) /fabs(spread_rs[l][k+1]*spread_rs[l][k+1]-spread_rs[l][k]*spread_rs[l][k]);
#ifdef INDIVIDUAL_ELEMENTS
		    for (kk=0;kk<NUM_ELEMENTS;kk++)
		      mtot_elements[kk]+= gasElements_old[l][kk]*fracs[k]* fabs(r2*r2-r3*r3) /fabs(spread_rs[l][k+1]*spread_rs[l][k+1]-spread_rs[l][k]*spread_rs[l][k]);
#endif		  
		  
		  }
	      }
 
	  }

	Gal[p].ColdGasRings[j]+=mtot;

	for(ii=0;ii<NUM_METAL_CHANNELS;ii++)
	  Gal[p].MetalsColdGasRings[j][ii]+= mtot_metals[ii];
#ifdef INDIVIDUAL_ELEMENTS
	for (kk=0;kk<NUM_ELEMENTS;kk++)
	  Gal[p].ColdGasRings_elements[j][kk]+=mtot_elements[kk];
#endif            
      } 
    if (DiskRadiusModel == 0)
      Gal[p].ColdGasRadius=get_gas_disk_radius(p);

    
  }
  

}




