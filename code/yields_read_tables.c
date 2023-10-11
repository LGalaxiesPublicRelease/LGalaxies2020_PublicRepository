#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "allvars.h"
#include "proto.h"

void read_yield_tables(void)
{
#ifdef DETAILED_DUST
	// Check that the requisite elements are loadable for DETAILED_DUST to run:
#if !defined(Cb_NUM) || !defined(O_NUM) || !defined(Mg_NUM) || !defined(Si_NUM) || !defined(Fe_NUM)
	printf("***** ERROR: read_yield_tables.c: To run with DETAILED_DUST on, the following elements must be included in the yield tables: C, O, Mg, Si, Fe.");
	printf("***** Check h_metals.h to see which elements are currently included.");
	exit(0);
#endif
#endif //DETAILED_DUST

#ifdef BINARYC
	static char *bc_yields_filename = "ensemble_output_";
	static char *bc_yields_directory = "./YieldTables/binary_c_yields/default/binaryStars/";
	static char *bc_SNIa_yields_directory = "/";

	//------------------------------------------
	//READ TIMESTEPS LIST:
	//------------------------------------------
	FILE *fd1;
	char buf1[100];
	int i1;
	double m1;
	static char *name1 = "logt_timesteps.txt";

	sprintf(buf1, "%s%s%s", bc_yields_directory, bc_yields_filename, name1);

	if(!(fd1 = fopen(buf1, "r"))) {
		printf("file `%s' not found.\n", buf1);
		exit(0);
	}

	for(i1=0; i1<BC_TIME_NUM; i1++) {
		fscanf(fd1, "%lf", &m1); //Store at double precision
		bcTimes[i1] = pow(10,m1)*1.e6; //[yr] Convert log times to linear times, and Myr to yr.
	}
	fclose(fd1);

	//------------------------------------------
	//READ METALLICITY LIST:
	//------------------------------------------
	FILE *fd2;
	char buf2[100];
	int i2;
	float m2;
	static char *name2 = "metallicities.txt";

	sprintf(buf2, "%s%s%s", bc_yields_directory, bc_yields_filename, name2);

	if(!(fd2 = fopen(buf2, "r"))) {
		printf("file `%s' not found.\n", buf2);
		exit(0);
	}

	for(i2=0; i2<BC_Z_NUM; i2++) {
		fscanf(fd2, "%f", &m2); //Store at float precision, to match the metallicity list from the lifetime tables below.
		bcMetallicities[i2] = m2;
	}
	fclose(fd2);

	//------------------------------------------
	//READ AGB EJECTED MASS LISTS:
	//------------------------------------------
	FILE *fd3;
	char buf3[100];
	int i3,j3;
	double m3;
	static char *name3[] = {
	  	"AGB_Z0.0001_EjectedMasses.txt",
	  	"AGB_Z0.001_EjectedMasses.txt",
		"AGB_Z0.004_EjectedMasses.txt",
		"AGB_Z0.008_EjectedMasses.txt",
		"AGB_Z0.01_EjectedMasses.txt",
	  	"AGB_Z0.03_EjectedMasses.txt"
		};

	for(j3=0; j3<BC_Z_NUM; j3++) {
		sprintf(buf3, "%s%s%s", bc_yields_directory, bc_yields_filename, name3[j3]);

		if(!(fd3 = fopen(buf3, "r"))) {
			printf("file `%s' not found.\n", buf3);
			exit(0);
		}

		for(i3=0; i3<BC_TIME_NUM; i3++) {
			fscanf(fd3, "%lf", &m3); //Store at double precision
			bcAGBEjectedMasses[j3][i3] = m3 / (bcTimes[i3] * log(10.)); //Convert from log(t/Myr) bins to linear t/yr bins. (For details, see below)
			/*
			 * Note:
			 * m3 = d(M/Msun)/dlog(t/Myr) * 1/Msun  //Default units
			 * m3 / ((bcTimes[i3] / 1.e6) * log(10.)) = d(M/Msun)/d(t/Myr) * 1/Msun  //N.B. bcTimes[] is in yr
			 * m3 / (bcTimes[i3] * log(10.)) = d(M/Msun)/d(t/yr) * 1/Msun  //Units required for yields_integrals.c
			*/
		}
		fclose(fd3);
	}

	//------------------------------------------
	//READ SNII EJECTED MASS LISTS:
	//------------------------------------------
	FILE *fd4;
	char buf4[100];
	int i4,j4;
	double m4;
	static char *name4[] = {
		"SNII_Z0.0001_EjectedMasses.txt",
		"SNII_Z0.001_EjectedMasses.txt",
		"SNII_Z0.004_EjectedMasses.txt",
		"SNII_Z0.008_EjectedMasses.txt",
		"SNII_Z0.01_EjectedMasses.txt",
		"SNII_Z0.03_EjectedMasses.txt"
		};

	for(j4=0; j4<BC_Z_NUM; j4++) {
		sprintf(buf4, "%s%s%s", bc_yields_directory, bc_yields_filename, name4[j4]);

		if(!(fd4 = fopen(buf4, "r"))) {
			printf("file `%s' not found.\n", buf4);
			exit(0);
		}

		for(i4=0; i4<BC_TIME_NUM; i4++) {
			fscanf(fd4, "%lf", &m4); //Store at double precision
			bcSNIIEjectedMasses[j4][i4] = m4 / (bcTimes[i4] * log(10.)); //Convert from log(t/Myr) bins to linear t/yr bins. (For details, see below)
			/*
			 * Note:
			 * m3 = d(M/Msun)/dlog(t/Myr) * 1/Msun  //Default units
			 * m3 / ((bcTimes[i3] / 1.e6) * log(10.)) = d(M/Msun)/d(t/Myr) * 1/Msun  //N.B. bcTimes[] is in yr
			 * m3 / (bcTimes[i3] * log(10.)) = d(M/Msun)/d(t/yr) * 1/Msun  //Units required for yields_integrals.c
			*/
		}
		fclose(fd4);
	}

	//------------------------------------------
	//READ SNIA EJECTED MASS LISTS:
	//------------------------------------------
	FILE *fd5;
	char buf5[100];
	int i5,j5;
	double m5;
	static char *name5[] = {
		"SNIa_Z0.0001_EjectedMasses.txt",
		"SNIa_Z0.001_EjectedMasses.txt",
		"SNIa_Z0.004_EjectedMasses.txt",
		"SNIa_Z0.008_EjectedMasses.txt",
		"SNIa_Z0.01_EjectedMasses.txt",
		"SNIa_Z0.03_EjectedMasses.txt"
		};

	for(j5=0; j5<BC_Z_NUM; j5++) {
		sprintf(buf5, "%s%s%s%s", bc_yields_directory, bc_SNIa_yields_directory, bc_yields_filename, name5[j5]);

		if(!(fd5 = fopen(buf5, "r"))) {
			printf("file `%s' not found.\n", buf5);
			exit(0);
		}

		for(i5=0; i5<BC_TIME_NUM; i5++) {
			fscanf(fd5, "%lf", &m5); //Store at double precision
			bcSNIaEjectedMasses[j5][i5] = (BC_SNIA_SCALE_FACTOR * m5) / (bcTimes[i5] * log(10.)); //Convert from log(t/Myr) bins to linear t/yr bins. (For details, see below)
			/*
			 * Note:
			 * m3 = d(M/Msun)/dlog(t/Myr) * 1/Msun  //Default units
			 * m3 / ((bcTimes[i3] / 1.e6) * log(10.)) = d(M/Msun)/d(t/Myr) * 1/Msun  //N.B. bcTimes[] is in yr
			 * m3 / (bcTimes[i3] * log(10.)) = d(M/Msun)/d(t/yr) * 1/Msun  //Units required for yields_integrals.c
			*/
		}
		fclose(fd5);
	}

	//------------------------------------------
	//READ AGB TOTAL METALS LISTS:
	//------------------------------------------
	FILE *fd6;
	char buf6[100];
	int i6,j6;
	double m6;
	static char *name6[] = {
	  	"AGB_Z0.0001_TotalMetals.txt",
	  	"AGB_Z0.001_TotalMetals.txt",
		"AGB_Z0.004_TotalMetals.txt",
		"AGB_Z0.008_TotalMetals.txt",
		"AGB_Z0.01_TotalMetals.txt",
	  	"AGB_Z0.03_TotalMetals.txt"
		};

	for(j6=0; j6<BC_Z_NUM; j6++) {
		sprintf(buf6, "%s%s%s", bc_yields_directory, bc_yields_filename, name6[j6]);

		if(!(fd6 = fopen(buf6, "r"))) {
			printf("file `%s' not found.\n", buf6);
			exit(0);
		}

		for(i6=0; i6<BC_TIME_NUM; i6++) {
			fscanf(fd6, "%lf", &m6); //Store at double precision
			bcAGBTotalMetals[j6][i6] = m6 / (bcTimes[i6] * log(10.)); //Convert from log(t/Myr) bins to linear t/yr bins.
		}
		fclose(fd6);
	}

	//------------------------------------------
	//READ SNII TOTAL METALS LISTS:
	//------------------------------------------
	FILE *fd7;
	char buf7[100];
	int i7,j7;
	double m7;
	static char *name7[] = {
		"SNII_Z0.0001_TotalMetals.txt",
		"SNII_Z0.001_TotalMetals.txt",
		"SNII_Z0.004_TotalMetals.txt",
		"SNII_Z0.008_TotalMetals.txt",
		"SNII_Z0.01_TotalMetals.txt",
		"SNII_Z0.03_TotalMetals.txt"
		};

	for(j7=0; j7<BC_Z_NUM; j7++) {
		sprintf(buf7, "%s%s%s", bc_yields_directory, bc_yields_filename, name7[j7]);

		if(!(fd7 = fopen(buf7, "r"))) {
			printf("file `%s' not found.\n", buf7);
			exit(0);
		}

		for(i7=0; i7<BC_TIME_NUM; i7++) {
			fscanf(fd7, "%lf", &m7); //Store at double precision
			bcSNIITotalMetals[j7][i7] = m7 / (bcTimes[i7] * log(10.)); //Convert from log(t/Myr) bins to linear t/yr bins.
		}
		fclose(fd7);
	}

	//------------------------------------------
	//READ SNIA TOTAL METALS LISTS:
	//------------------------------------------
	FILE *fd8;
	char buf8[100];
	int i8,j8;
	double m8;
	static char *name8[] = {
		"SNIa_Z0.0001_TotalMetals.txt",
		"SNIa_Z0.001_TotalMetals.txt",
		"SNIa_Z0.004_TotalMetals.txt",
		"SNIa_Z0.008_TotalMetals.txt",
		"SNIa_Z0.01_TotalMetals.txt",
		"SNIa_Z0.03_TotalMetals.txt"
		};

	for(j8=0; j8<BC_Z_NUM; j8++) {
		sprintf(buf8, "%s%s%s%s", bc_yields_directory, bc_SNIa_yields_directory, bc_yields_filename, name8[j8]);

		if(!(fd8 = fopen(buf8, "r"))) {
			printf("file `%s' not found.\n", buf8);
			exit(0);
		}

		for(i8=0; i8<BC_TIME_NUM; i8++) {
			fscanf(fd8, "%lf", &m8); //Store at double precision
			bcSNIaTotalMetals[j8][i8] = (BC_SNIA_SCALE_FACTOR * m8) / (bcTimes[i8] * log(10.)); //Convert from log(t/Myr) bins to linear t/yr bins.
		}
		fclose(fd8);
	}

	//------------------------------------------
	//READ AGB YIELDS LISTS:
	//------------------------------------------
	FILE *fd9;
	char buf9[100];
	int i9,j9,k9;
	double m9;
	static char *name9[] = {
	  	"AGB_Z0.0001_Yields.txt",
	  	"AGB_Z0.001_Yields.txt",
		"AGB_Z0.004_Yields.txt",
		"AGB_Z0.008_Yields.txt",
		"AGB_Z0.01_Yields.txt",
	  	"AGB_Z0.03_Yields.txt"
		};

	for(j9=0; j9<BC_Z_NUM; j9++) {
		sprintf(buf9, "%s%s%s", bc_yields_directory, bc_yields_filename, name9[j9]);

		if(!(fd9 = fopen(buf9, "r"))) {
			printf("file `%s' not found.\n", buf9);
			exit(0);
		}

		for(k9=0; k9<NUM_ELEMENTS; k9++) {
			for(i9=0; i9<BC_TIME_NUM; i9++) {
				fscanf(fd9, "%lf", &m9); //Store at double precision
				bcAGBYields[j9][k9][i9] = m9 / (bcTimes[i9] * log(10.)); //Convert from log(t/Myr) bins to linear t/yr bins.
				if (bcAGBYields[j9][k9][i9] < 0.0) {
					printf("***** WARNING!: [yields_read_tables.c] bcAGBYields[%i][%i][%i] is negative (= %e | raw yield = %e). Forcing it to 0.0. *****\n",
				            j9, k9, i9, bcAGBYields[j9][k9][i9], m9);
					bcAGBYields[j9][k9][i9] = 0.0;
				}

			}
		}
		fclose(fd9);
	}

	//------------------------------------------
	//READ SNII YIELDS LISTS:
	//------------------------------------------
	FILE *fd10;
	char buf10[100];
	int i10,j10,k10;
	double m10;
	static char *name10[] = {
		"SNII_Z0.0001_Yields.txt",
		"SNII_Z0.001_Yields.txt",
		"SNII_Z0.004_Yields.txt",
		"SNII_Z0.008_Yields.txt",
		"SNII_Z0.01_Yields.txt",
		"SNII_Z0.03_Yields.txt"
		};

	for(j10=0; j10<BC_Z_NUM; j10++) {
		sprintf(buf10, "%s%s%s", bc_yields_directory, bc_yields_filename, name10[j10]);

		if(!(fd10 = fopen(buf10, "r"))) {
			printf("file `%s' not found.\n", buf10);
			exit(0);
		}

		for(k10=0; k10<NUM_ELEMENTS; k10++) {
			for(i10=0; i10<BC_TIME_NUM; i10++) {
				fscanf(fd10, "%lf", &m10); //Store at double precision
				bcSNIIYields[j10][k10][i10] = m10 / (bcTimes[i10] * log(10.)); //Convert from log(t/Myr) bins to linear t/yr bins.
				if (bcSNIIYields[j10][k10][i10] < 0.0) {
					printf("***** WARNING!: [yields_read_tables.c] bcSNIIYields[%i][%i][%i] is negative (= %e | raw yield = %e). *****\n",
							j10, k10, i10, bcSNIIYields[j10][k10][i10], m10);
				}
			}
		}
		fclose(fd10);
	}

	//------------------------------------------
	//READ SNIa YIELDS LISTS:
	//------------------------------------------
	FILE *fd11;
	char buf11[100];
	int i11,j11,k11;
	double m11;
	static char *name11[] = {
		"SNIa_Z0.0001_Yields.txt",
		"SNIa_Z0.001_Yields.txt",
		"SNIa_Z0.004_Yields.txt",
		"SNIa_Z0.008_Yields.txt",
		"SNIa_Z0.01_Yields.txt",
		"SNIa_Z0.03_Yields.txt"
		};

	for(j11=0; j11<BC_Z_NUM; j11++) {
		sprintf(buf11, "%s%s%s%s", bc_yields_directory, bc_SNIa_yields_directory, bc_yields_filename, name11[j11]);

		if(!(fd11 = fopen(buf11, "r"))) {
			printf("file `%s' not found.\n", buf11);
			exit(0);
		}

		for(k11=0; k11<NUM_ELEMENTS; k11++) {
			for(i11=0; i11<BC_TIME_NUM; i11++) {
				fscanf(fd11, "%lf", &m11); //Store at double precision
				bcSNIaYields[j11][k11][i11] = (BC_SNIA_SCALE_FACTOR * m11) / (bcTimes[i11] * log(10.)); //Convert from log(t/Myr) bins to linear t/yr bins.
				if (bcSNIaYields[j11][k11][i11] < 0.0) {
					printf("***** WARNING!: [yields_read_tables.c] bcSNIaYields[%i][%i][%i] is negative (= %e | raw yield = %e). *****\n",
							j11, k11, i11, bcSNIaYields[j11][k11][i11], m11);
				}
			}
		}
		fclose(fd11);
	}

	//------------------------------------------
	//READ AGB RATES LISTS:
	//------------------------------------------
	FILE *fd12;
	char buf12[100];
	int i12,j12;
	double m12;
	static char *name12[] = {
	  	"AGB_Z0.0001_Rates.txt",
	  	"AGB_Z0.001_Rates.txt",
		"AGB_Z0.004_Rates.txt",
		"AGB_Z0.008_Rates.txt",
		"AGB_Z0.01_Rates.txt",
	  	"AGB_Z0.03_Rates.txt"
		};

	for(j12=0; j12<BC_Z_NUM; j12++) {
		sprintf(buf12, "%s%s%s", bc_yields_directory, bc_yields_filename, name12[j12]);

		if(!(fd12 = fopen(buf12, "r"))) {
			printf("file `%s' not found.\n", buf12);
			exit(0);
		}

		for(i12=0; i12<BC_TIME_NUM; i12++) {
			fscanf(fd12, "%lf", &m12); //Store at double precision
			bcAGBRates[j12][i12] = m12 / (bcTimes[i12] * log(10.)); //[in dNum/dt 1/Msun] I.e. converted from log(t/Myr) bins to linear t/yr bins.
		}
		fclose(fd12);
	}

	//------------------------------------------
	//READ SNII RATES LISTS:
	//------------------------------------------
	FILE *fd13;
	char buf13[100];
	int i13,j13;
	double m13;
	static char *name13[] = {
		"SNII_Z0.0001_Rates.txt",
		"SNII_Z0.001_Rates.txt",
		"SNII_Z0.004_Rates.txt",
		"SNII_Z0.008_Rates.txt",
		"SNII_Z0.01_Rates.txt",
		"SNII_Z0.03_Rates.txt"
		};

	for(j13=0; j13<BC_Z_NUM; j13++) {
		sprintf(buf13, "%s%s%s", bc_yields_directory, bc_yields_filename, name13[j13]);

		if(!(fd13 = fopen(buf13, "r"))) {
			printf("file `%s' not found.\n", buf13);
			exit(0);
		}

		for(i13=0; i13<BC_TIME_NUM; i13++) {
			fscanf(fd13, "%lf", &m13); //Store at double precision
			bcSNIIRates[j13][i13] = m13 / (bcTimes[i13] * log(10.)); //[in dNum/dt 1/Msun] I.e. converted from log(t/Myr) bins to linear t/yr bins.
		}
		fclose(fd13);
	}

	//------------------------------------------
	//READ SNIA RATES LISTS:
	//------------------------------------------
	FILE *fd14;
	char buf14[100];
	int i14,j14;
	double m14;
	static char *name14[] = {
		"SNIa_Z0.0001_Rates.txt",
		"SNIa_Z0.001_Rates.txt",
		"SNIa_Z0.004_Rates.txt",
		"SNIa_Z0.008_Rates.txt",
		"SNIa_Z0.01_Rates.txt",
		"SNIa_Z0.03_Rates.txt"
		};

	for(j14=0; j14<BC_Z_NUM; j14++) {
		sprintf(buf14, "%s%s%s%s", bc_yields_directory, bc_SNIa_yields_directory, bc_yields_filename, name14[j14]);

		if(!(fd14 = fopen(buf14, "r"))) {
			printf("file `%s' not found.\n", buf14);
			exit(0);
		}

		for(i14=0; i14<BC_TIME_NUM; i14++) {
			fscanf(fd14, "%lf", &m14); //Store at double precision
			bcSNIaRates[j14][i14] = (BC_SNIA_SCALE_FACTOR * m14) / (bcTimes[i14] * log(10.)); //[d/dlog(t/Myr) 1/Msun --> d/d(t/yr) 1/Msun] I.e. converted from log(t/Myr) bins to linear t/yr bins.
		}
		fclose(fd14);
	}

	printf("Binary_c tables read.\n");
	//exit(0);

#ifdef PARALLEL
	if ( ThisTask == 0 )
#endif
#endif //BINARYC

//#ifdef DETAILED_DUST
#if !defined(BINARYC) || defined(DETAILED_DUST)
	//N.B. 26-08-22: If DETAILED_DUST is on, an estimate of single AGB star lifetimes is needed when calculating AGB dust ejecta masses in dust_integrals.c
	//Therefore, we need to load the single-star "lifetimes" tables from Portinari+98 for this, even when BINARYC is on.
	//These lifetime tables are not used to calculate chemical element yields when BINARYC is on, only when it is off.

	//------------------------------------------
	//READ LIFETIME MASS LIST:
	//------------------------------------------
	FILE *fd1L;
	char buf1L[100];
	int i1L;
	float m1L;
	static char *name1L = "stripped_interp_LifetimeMasses.txt";

	sprintf(buf1L, "./YieldTables/%s", name1L);

	if(!(fd1L = fopen(buf1L, "r")))
        {
          printf("file `%s' not found.\n", buf1L);
          exit(0);
        }

	for(i1L=0; i1L<LIFETIME_MASS_NUM; i1L++)
        {
	  fscanf(fd1L, "%f", &m1L);
	  lifetimeMasses[i1L] = m1L;
        }
	fclose(fd1L);

	//------------------------------------------
	//READ LIFETIME METALLICITY LIST:
	//------------------------------------------
	FILE *fd2L;
	char buf2L[100];
	int i2L;
	float m2L;
	static char *name2L = "stripped_LifetimeMetallicities.txt";

	sprintf(buf2L, "./YieldTables/%s", name2L);

	if(!(fd2L = fopen(buf2L, "r")))
        {
          printf("file `%s' not found.\n", buf2L);
          exit(0);
        }

	for(i2L=0; i2L<LIFETIME_Z_NUM; i2L++)
        {
	  fscanf(fd2L, "%f", &m2L);
	  lifetimeMetallicities[i2L] = m2L;
        }
	fclose(fd2L);

	//------------------------------------------
	//READ LIFETIME TABLE:
	//------------------------------------------
	FILE *fd3L;
	char buf3L[100];
	int i3L,j3L;
	float m3L; //,mass;
	//float get_mass(float time, float Z0);
	static char *name3L = "stripped_interp_Lifetimes.txt";
	sprintf(buf3L, "./YieldTables/%s", name3L);

	if(!(fd3L = fopen(buf3L, "r")))
        {
          printf("file `%s' not found.\n", buf3L);
          exit(0);
        }

	for(i3L=0; i3L<LIFETIME_Z_NUM; i3L++)
        {
          for(j3L=0; j3L<LIFETIME_MASS_NUM; j3L++)
          {
        	  fscanf(fd3L, "%f", &m3L);
        	  lifetimes[i3L][j3L] = m3L;
          }
        }
	fclose(fd3L);
#ifdef PARALLEL
	if ( ThisTask == 0 )
#endif
	  printf("Lifetime tables read.\n");

#endif //!defined(BINARYC) || defined(DETAILED_DUST)

//#else //BINARYC
#ifndef BINARYC

	//------------------------------------------
	//READ AGB MASS LIST:
	//------------------------------------------
	FILE *fd4;
	char buf4[100];
	int i4;
	float m4;
	static char *name4 = "stripped_interp_AGBMasses.txt";

	sprintf(buf4, "./YieldTables/%s", name4);

	if(!(fd4 = fopen(buf4, "r")))
        {
          printf("file `%s' not found.\n", buf4);
          exit(0);
        }

	for(i4=0; i4<AGB_MASS_NUM; i4++)
        {
	  fscanf(fd4, "%f", &m4);
	  AGBMasses[i4] = m4;
	    }
	fclose(fd4);

	//------------------------------------------
	//READ AGB METALLICITY LIST:
	//------------------------------------------
	FILE *fd5;
	char buf5[100];
	int i5;
	float m5;
	static char *name5 = "stripped_AGBMetallicities.txt";

	sprintf(buf5, "./YieldTables/%s", name5);

	if(!(fd5 = fopen(buf5, "r")))
        {
          printf("file `%s' not found.\n", buf5);
          exit(0);
        }

	for(i5=0; i5<AGB_Z_NUM; i5++)
        {
	  fscanf(fd5, "%f", &m5);
	  AGBMetallicities[i5] = m5;
        }
	fclose(fd5);

	//------------------------------------------
	//READ AGB EJECTED MASS LISTS:
	//------------------------------------------
	FILE *fd6;
	char buf6[100];
	int i6,j6;
	float m6;
	static char *name6[] = {
  	"convolved_stripped_interp_AGB_Z004_EjectedMasses.txt",
  	"convolved_stripped_interp_AGB_Z008_EjectedMasses.txt",
  	"convolved_stripped_interp_AGB_Z019_EjectedMasses.txt"
	};

	for(i6 = 0; i6 < AGB_Z_NUM; i6++)
    {
		sprintf(buf6, "./YieldTables/%s", name6[i6]);

		if(!(fd6 = fopen(buf6, "r")))
        {
            printf("file `%s' not found.\n", buf6);
            exit(0);
        }

        for(j6=0; j6<AGB_MASS_NUM; j6++)
        {
        	  fscanf(fd6, "%f", &m6);
        	  AGBEjectedMasses[i6][j6] = m6 * Chabrier_IMF(AGBMasses[j6]);
        }
        fclose(fd6);
	 }

	//------------------------------------------
	//READ AGB TOTAL METALS LISTS: This is the sum of all the yields NOT including H and He:
	//------------------------------------------
	FILE *fd7;
	char buf7[100];
	int i7,j7;
	float m7;
	static char *name7[] = {
  	"convolved_stripped_interp_AGB_Z004_TotalMetals.txt",
  	"convolved_stripped_interp_AGB_Z008_TotalMetals.txt",
  	"convolved_stripped_interp_AGB_Z019_TotalMetals.txt"
	};

	for(i7 = 0; i7 < AGB_Z_NUM; i7++)
    	{
	  sprintf(buf7, "./YieldTables/%s", name7[i7]);

	  if(!(fd7 = fopen(buf7, "r")))
          {
            printf("file `%s' not found.\n", buf7);
            exit(0);
          }

            for(j7=0; j7<AGB_MASS_NUM; j7++)
	    {
	      fscanf(fd7, "%f", &m7);
	      AGBTotalMetals[i7][j7] = m7 * Chabrier_IMF(AGBMasses[j7]);
	    }
	  fclose(fd7);
	}

	//------------------------------------------
	//READ AGB YIELD TABLES:
	//------------------------------------------
	FILE *fd8;
	char buf8[100];
	int k8,i8,j8;
	float m8;
	static char *name8[] = {
  	"convolved_stripped_interp_AGB_Z004_Yields.txt",
  	"convolved_stripped_interp_AGB_Z008_Yields.txt",
  	"convolved_stripped_interp_AGB_Z019_Yields.txt"
	};

	for(k8 = 0; k8 < AGB_Z_NUM; k8++)
    {
	  sprintf(buf8, "./YieldTables/%s", name8[k8]);

	  if(!(fd8 = fopen(buf8, "r")))
          {
            printf("file `%s' not found.\n", buf8);
            exit(0);
          }

	  for(i8=0; i8<NUM_ELEMENTS; i8++) //Number of element species (inc H and He) = 11
          {
          for(j8=0; j8<AGB_MASS_NUM; j8++)
            {
            	fscanf(fd8, "%f", &m8);
            	AGBYields[k8][i8][j8] = m8 * Chabrier_IMF(AGBMasses[j8]);
            }
          }
    }

#ifdef PARALLEL
	if ( ThisTask == 0 )
#endif
	  printf("AGB yield tables read.\n");

	//------------------------------------------
	//READ SN-II MASS LIST:
	//------------------------------------------
	FILE *fd9;
	char buf9[100];
	int i9;
	float m9;
#ifdef PORTINARI
	static char *name9 = "stripped_interp_SNIIMasses.txt";
#endif
#ifdef CHIEFFI
	static char *name9 = "stripped_interp_CL_SNIIMasses.txt";
#endif

	sprintf(buf9, "./YieldTables/%s", name9);

	if(!(fd9 = fopen(buf9, "r")))
        {
          printf("file `%s' not found.\n", buf9);
          exit(0);
        }

	for(i9=0; i9<SNII_MASS_NUM; i9++)
        {
	  fscanf(fd9, "%f", &m9);
	  SNIIMasses[i9] = m9;
        }
	fclose(fd9);

	//------------------------------------------
	//READ SN-II METALLICITY LIST:
	//------------------------------------------
	FILE *fd10;
	char buf10[100];
	int i10;
	float m10;
#ifdef PORTINARI
	static char *name10 = "stripped_SNIIMetallicities.txt";
#endif
#ifdef CHIEFFI
	static char *name10 = "stripped_CL_SNIIMetallicities.txt";
#endif

	sprintf(buf10, "./YieldTables/%s", name10);

	if(!(fd10 = fopen(buf10, "r")))
        {
          printf("file `%s' not found.\n", buf10);
          exit(0);
        }

	for(i10=0; i10<SNII_Z_NUM; i10++)
        {
			fscanf(fd10, "%f", &m10);
			SNIIMetallicities[i10] = m10;
        }
	fclose(fd10);

	//------------------------------------------
	//READ SN-II EJECTED MASS LISTS:
	//------------------------------------------
	FILE *fd11;
	char buf11[100];
	int i11,j11;
	float m11;
#ifdef PORTINARI
	static char *name11[] = {
  	"convolved_stripped_interp_SNII_Z0004_EjectedMasses.txt",
	"convolved_stripped_interp_SNII_Z004_EjectedMasses.txt",
  	"convolved_stripped_interp_SNII_Z008_EjectedMasses.txt",
  	"convolved_stripped_interp_SNII_Z02_EjectedMasses.txt",
	"convolved_stripped_interp_SNII_Z05_EjectedMasses.txt"
	};
#endif
#ifdef CHIEFFI
	static char *name11[] = {
  	"convolved_stripped_interp_CL_SNII_Z0_EjectedMasses.txt",
	"convolved_stripped_interp_CL_SNII_Z000001_EjectedMasses.txt",
  	"convolved_stripped_interp_CL_SNII_Z0001_EjectedMasses.txt",
  	"convolved_stripped_interp_CL_SNII_Z001_EjectedMasses.txt",
	"convolved_stripped_interp_CL_SNII_Z006_EjectedMasses.txt",
  	"convolved_stripped_interp_CL_SNII_Z02_EjectedMasses.txt"
	};
#endif

	for(i11 = 0; i11 < SNII_Z_NUM; i11++)
    	{
	  sprintf(buf11, "./YieldTables/%s", name11[i11]);

	  if(!(fd11 = fopen(buf11, "r")))
          {
            printf("file `%s' not found.\n", buf11);
            exit(0);
          }
          for(j11=0; j11<SNII_MASS_NUM; j11++)
	  {
	    fscanf(fd11, "%f", &m11);
	    SNIIEjectedMasses[i11][j11] = m11 * Chabrier_IMF(SNIIMasses[j11]);
	  }
	  fclose(fd11);
	}

	//------------------------------------------
	//READ SN-II TOTAL METALS LISTS:
	//------------------------------------------
	FILE *fd12;
	char buf12[100];
	int i12,j12;
	float m12;
#ifdef PORTINARI
	static char *name12[] = {
  	"convolved_stripped_interp_SNII_Z0004_TotalMetals.txt",
	"convolved_stripped_interp_SNII_Z004_TotalMetals.txt",
  	"convolved_stripped_interp_SNII_Z008_TotalMetals.txt",
  	"convolved_stripped_interp_SNII_Z02_TotalMetals.txt",
	"convolved_stripped_interp_SNII_Z05_TotalMetals.txt"
	};
#endif
#ifdef CHIEFFI
	static char *name12[] = {
  	"convolved_stripped_interp_CL_SNII_Z0_TotalMetals.txt",
	"convolved_stripped_interp_CL_SNII_Z000001_TotalMetals.txt",
  	"convolved_stripped_interp_CL_SNII_Z0001_TotalMetals.txt",
  	"convolved_stripped_interp_CL_SNII_Z001_TotalMetals.txt",
	"convolved_stripped_interp_CL_SNII_Z006_TotalMetals.txt",
  	"convolved_stripped_interp_CL_SNII_Z02_TotalMetals.txt"
	};
#endif

	for(i12 = 0; i12 < SNII_Z_NUM; i12++)
    {
	  sprintf(buf12, "./YieldTables/%s", name12[i12]);

	  if(!(fd12 = fopen(buf12, "r")))
      {
            printf("file `%s' not found.\n", buf12);
            exit(0);
      }

      for(j12=0; j12<SNII_MASS_NUM; j12++)
	  {
	    fscanf(fd12, "%f", &m12);
	    SNIITotalMetals[i12][j12] = m12 * Chabrier_IMF(SNIIMasses[j12]);
	  }
	  fclose(fd12);
	}

	//------------------------------------------
	//READ SN-II YIELD TABLES:
	//------------------------------------------
	FILE *fd13;
	char buf13[100];
	int k13,i13,j13;
	float m13;
#ifdef PORTINARI
	static char *name13[] = {
  	"convolved_stripped_interp_SNII_Z0004_Yields.txt",
	"convolved_stripped_interp_SNII_Z004_Yields.txt",
  	"convolved_stripped_interp_SNII_Z008_Yields.txt",
  	"convolved_stripped_interp_SNII_Z02_Yields.txt",
	"convolved_stripped_interp_SNII_Z05_Yields.txt"
	};
#endif
#ifdef CHIEFFI
	static char *name13[] = {
  	"convolved_stripped_interp_CL_SNII_Z0_Yields.txt",
	"convolved_stripped_interp_CL_SNII_Z000001_Yields.txt",
  	"convolved_stripped_interp_CL_SNII_Z0001_Yields.txt",
  	"convolved_stripped_interp_CL_SNII_Z001_Yields.txt",
	"convolved_stripped_interp_CL_SNII_Z006_Yields.txt",
  	"convolved_stripped_interp_CL_SNII_Z02_Yields.txt"
	};
#endif

	for(k13 = 0; k13 < SNII_Z_NUM; k13++)
    	{
	  sprintf(buf13, "./YieldTables/%s", name13[k13]);

	  if(!(fd13 = fopen(buf13, "r")))
          {
            printf("file `%s' not found.\n", buf13);
            exit(0);
          }

	  for(i13=0; i13<NUM_ELEMENTS; i13++) //Number of element species (inc H and He) = 11
          {
            for(j13=0; j13<SNII_MASS_NUM; j13++) //Number of initial masses = 85 (for Portinari yields) (11 in initial yield tables)
	    {
	      fscanf(fd13, "%f", &m13);
	      SNIIYields[k13][i13][j13] = m13 * Chabrier_IMF(SNIIMasses[j13]);
	    }
	  }
	}

#ifdef PARALLEL
	if ( ThisTask == 0 )
#endif
	  printf("SN-II yield tables read.\n");

#ifndef DTD
	//------------------------------------------
	//READ SN-Ia MASS LIST:
	//------------------------------------------
	FILE *fd17;
	char buf17[100];
	int i17;
	float m17;
	static char *name17 = "stripped_interp_SNIaMasses.txt";

	sprintf(buf17, "./YieldTables/%s", name17);

	if(!(fd17 = fopen(buf17, "r")))
        {
          printf("file `%s' not found.\n", buf17);
          exit(0);
        }

	for(i17=0; i17<SNIA_MASS_NUM; i17++)
        {
	  fscanf(fd17, "%f", &m17);
	  SNIaMasses[i17] = m17;
        }
	fclose(fd17);

	//------------------------------------------
	//READ SN-Ia TOTAL METALS:
	//------------------------------------------
	FILE *fd15;
	char buf15[100];
	int i15;
	float m15;
	static char *name15 = "convolved_stripped_interp_SNIa_TotalMetals.txt";

	sprintf(buf15, "./YieldTables/%s", name15);

	 if(!(fd15 = fopen(buf15, "r")))
     {
        printf("file `%s' not found.\n", buf15);
        exit(0);
     }

	 for(i15 = 0; i15 < SNIA_MASS_NUM; i15++)
	 {
    	fscanf(fd15, "%f", &m15);
    	SNIaTotalMetals[i15] = m15 * Chabrier_IMF(SNIaMasses[i15]);
	 }

		//------------------------------------------
		//READ SN-Ia EJECTED MASSES:
		//------------------------------------------
		FILE *fd16;
		char buf16[100];
		int i16;
		float m16;
		static char *name16 = "convolved_stripped_interp_SNIa_EjectedMasses.txt";

		sprintf(buf16, "./YieldTables/%s", name16);

		 if(!(fd16 = fopen(buf16, "r")))
	     {
	        printf("file `%s' not found.\n", buf16);
	        exit(0);
	     }

		 for(i16 = 0; i16 < SNIA_MASS_NUM; i16++)
		 {
	    	fscanf(fd16, "%f", &m16);
	    	SNIaEjectedMasses[i16] = m16 * Chabrier_IMF(SNIaMasses[i16]);
		 }

	//------------------------------------------
	//READ SN-Ia YIELD TABLES:
	//------------------------------------------
		FILE *fd14;
		char buf14[100];
		int i14,j14;
		float m14;
		static char *name14 = "convolved_stripped_interp_SNIa_Yields.txt";

		sprintf(buf14, "./YieldTables/%s", name14);

		 if(!(fd14 = fopen(buf14, "r")))
	     {
	        printf("file `%s' not found.\n", buf14);
	        exit(0);
	     }

		 for(i14=0; i14<42; i14++) //Number of element species (inc H and He) = 42
	     {
	         for(j14=0; j14<SNIA_MASS_NUM; j14++) //Number of initial masses = 48 (30 in initial table)
		     {
		      fscanf(fd14, "%f", &m14);
		      SNIaYields[i14][j14] = m14 * Chabrier_IMF(SNIaMasses[j14]);
		    }
		  }

#ifdef PARALLEL
		 if ( ThisTask == 0 )
#endif
		   printf("SN-Ia yield tables read.\n");

#else
		//------------------------------------------
		//READ SN-Ia YIELD TABLES:
		//------------------------------------------
			FILE *fd14;
			char buf14[100];
			int i14,j14;
			float m14;
			static char *name14 = "SNIaYields.txt";

			sprintf(buf14, "./YieldTables/%s", name14);

			 if(!(fd14 = fopen(buf14, "r")))
		     {
		        printf("file `%s' not found.\n", buf14);
		        exit(0);
		     }

			 for(i14=0; i14<42; i14++) //Number of element species (inc H and He) = 42
		     {
			      fscanf(fd14, "%f", &m14);
			      SNIaYields[i14] = m14;
			  }
#ifdef PARALLEL
			 if ( ThisTask == 0 )
#endif
			   printf("SN-Ia yield tables read.\n");

#endif
#endif //BINARYC


#ifdef WRITE_YIELD_DATA
//------------------------------------------
//READ SOLAR ABUNDANCES LIST:
//------------------------------------------
FILE *fd15;
char buf15[100];
int i15;
float m15;
static char *name15 = "solar_abundances_Asplund09.txt";

sprintf(buf15, "./YieldTables/BurstYields/%s", name15);

if(!(fd15 = fopen(buf15, "r")))
	{
	  printf("file `%s' not found.\n", buf15);
	  exit(0);
	}

for(i15=0; i15<11; i15++)
	{
  fscanf(fd15, "%f", &m15);
  solarAbundMassRatios[i15] = m15;
	}
fclose(fd15);
printf("Solar abundance ratios read.\n");
#endif //WRITE_YIELD_DATA

}


#ifndef BINARYC
double Chabrier_IMF(double M)
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

	const double mc = 0.079;
	const double sigma = 0.69;
	if(M >= 1.0)
	{
		phi = B*M*pow(M,-IMF_SLOPE);
	}
	else
	{
		e = (1./M)*exp(-pow(log10(M)-log10(mc),2.)/(2.*pow(sigma,2.)));
		phi = A*M*e;
	}

	return phi/M;
}
#endif //BINARYC

