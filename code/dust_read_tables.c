/*
 * recipe_dust_tables.c
 *
 *  Created on: Oct2016
 *      Author: scottclay
 *
 *  Edited from: Nov2021
 *  	Author: robyates
 *
 * Similar to Rob's yields_read_tables
 * Reads in the mass and metallicity dependent dust yield tables for use in the code
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "allvars.h"
#include "proto.h"

#ifdef DETAILED_DUST

void read_dust_yield_tables(void)
{
	//------------------------------------------
	//READ AGB MASS LIST:
	//------------------------------------------
	FILE *fd1;
	char buf1[100];
	int i1;
	float m1;
	static char *name1 = "Zhukovska_AGByields_Mass2.txt";

	sprintf(buf1, "./DustTables/%s", name1);

	if(!(fd1 = fopen(buf1, "r")))
        {
          printf("file `%s' not found.\n", buf1);
          exit(0);
        }

	for(i1=0; i1<AGB_DUST_MASS_NUM; i1++)
        {
	  fscanf(fd1, "%f", &m1);
	  AGBDustMasses[i1] = m1;
        }
	fclose(fd1);
	
	//------------------------------------------
	//READ AGB METALLICITY LIST:
	//------------------------------------------
	FILE *fd2;
	char buf2[100];
	int i2;
	float m2;
	static char *name2 = "Zhukovska_AGByields_Metallicity3.txt";

	sprintf(buf2, "./DustTables/%s", name2);

	if(!(fd2 = fopen(buf2, "r")))
        {
          printf("file `%s' not found.\n", buf2);
          exit(0);
        }

	for(i2=0; i2<AGB_DUST_METAL_NUM; i2++)
        {
	  fscanf(fd2, "%f", &m2);
	  AGBDustMetallicities[i2] = m2;
        }
	fclose(fd2);
	
	
	//------------------------------------------
	//READ AGB dust yields - ALL AGB star types C/M/S LIST:
	//------------------------------------------
	FILE *fd3;
	char buf3[100];
	int i3,j3,k3;
	float m3;
	static char *name3 = "Zhukovska_AGByields3.txt";

	sprintf(buf3, "./DustTables/%s", name3);

	if(!(fd3 = fopen(buf3, "r")))
        {
          printf("file `%s' not found.\n", buf3);
          exit(0);
        }

	for(i3=0; i3<AGB_DUST_METAL_NUM; i3++)
        {
        for (j3=0; j3<AGB_DUST_MASS_NUM; j3++)
        	{
        		for(k3=0; k3<AGB_DUST_TYPE_NUM; k3++)
        			{
        			fscanf(fd3, "%f", &m3);
        			AGBDustCreated[i3][j3][k3] = m3 * Chabrier_IMF_dust(AGBDustMasses[j3]);
	        	}
			}
		}
	fclose(fd3);
	printf("AGB dust yield tables read.\n");
}

double Chabrier_IMF_dust(double M)
{
	double e,phi;

	//Coefficient values, normalising the mass-weighted IMF as a function of M for, giving a normalised mass fraction.
	//FOR x = 2.3:
	double A, B;
	//For an IMF normalised over 0.1 --> 120.0 Msun:
	if (IMF_MAX_MASS == 120.0)
	{
		A = 0.842984;
		B = 0.235480;
	}
	//For an IMF normalised over 0.1 --> 100.0 Msun:
	else if (IMF_MAX_MASS == 100.0)
	{
		A = 0.852023;
		B = 0.238004;
	}
	//For an IMF normalised over 0.1 --> 80.0 Msun:
	else if (IMF_MAX_MASS == 80.0)
	{
		A = 0.864027;
		B = 0.241357;
	}
	//For an IMF normalised over 0.1 --> 70.0 Msun:
		else if (IMF_MAX_MASS == 70.0)
		{
			A = 0.871761;
			B = 0.243518;
		}
	//For an IMF normalised over 0.1 --> 60.0 Msun:
	else if (IMF_MAX_MASS == 60.0)
	{
		A = 0.881259;
		B = 0.246171;
	}
	//For an IMF normalised over 0.1 --> 50.0 Msun:
	else if (IMF_MAX_MASS == 50.0)
	{
		A = 0.893355;
		B = 0.249550;
	}
	//For an IMF normalised over 0.1 --> 40.0 Msun:
	else if (IMF_MAX_MASS == 40.0)
	{
		A = 0.909581;
		B = 0.254083;
	}
	//For an IMF normalised over 0.1 --> 30.0 Msun:
	else if (IMF_MAX_MASS == 30.0)
	{
		A = 0.933161;
		B = 0.260669;
	}
	//For an IMF normalised over 0.1 --> 25.0 Msun:
	else if (IMF_MAX_MASS == 25.0)
	{
		A = 0.950066;
		B = 0.265152;
	}
	else {printf("Chabrier_IMF(): Normalization constants for IMF not known. Check upper mass limit (IMF_MAX_MASS)."); exit(1);}

	//FOR x = 2.0:
	//For an IMF normalised over 0.1 --> 120.0 Msun:
	/*const double A = 0.551390;
	const double B = 0.154026;
	A = 0.551390;
	B = 0.154026;*/

	//const double x = 2.3; //Normal Chabrier IMF x = 2.3. Top-heavy IMF e.g. x = 2.0
	const double mc = 0.079;
	const double sigma = 0.69;

	if(M >= 1.0)
	{
		//phi = B*M*pow(M,-x); //ROB: For consistency, I have replaced the slope value defined just above as "x" with the global value defined in h_variables.h as "IMF_SLOPE". (11-02-19)
		phi = B*M*pow(M,-IMF_SLOPE);
	}
	else
	{
		e = (1./M)*exp(-pow(log10(M)-log10(mc),2.)/(2.*pow(sigma,2.)));
		phi = A*M*e;
	}

	return phi/M;
}

#endif
