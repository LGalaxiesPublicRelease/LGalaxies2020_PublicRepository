// Variables but some parameters hidden in here too.
// May be better to split by function

// I have no idea why the following variable is not declared extern,
// but compilation breaks if you do so.
char *inputFile; // Name of the input file
extern double MinGalOutputMass;

extern int FirstFile;		/* first and last file for processing */
extern int LastFile;

extern int Ntrees;		/* number of trees in current file */
extern double AllocValue_MaxHaloGal;
extern double AllocValue_MaxGal;
extern double AllocValue_MaxGalTree;

extern int MaxGal;		/* Maximum number of galaxies allowed for Gal[] array */
extern int NHaloGal, MaxHaloGal;
extern int NGalTree, MaxGalTree;
extern int *HaloGalHeap;
extern int IndexStored;

extern int LastSnapShotNr;

extern int LastDarkMatterSnapShot;
#ifdef MR_PLUS_MRII //OPTION for MCMC
extern int LastDarkMatterSnapShot_MR;
extern int LastDarkMatterSnapShot_MRII;
#endif


extern char SpecPhotDir[512];
extern char PhotPrefix[50];
extern char SpecPhotIMF[50];
extern char McFile[512];
extern char FileWithFilterNames[512];
extern char CoolFunctionsDir[512];
extern char CosmologyTablesDir[512];
extern char OutputDir[512];
/* in case a second parameter is given as argument to the code, this will be taken as a
 * temporary outputdir to allow fast I/O. OutputDir will be replaced by this directory
 * and in the end everything will be moved to the FinalOutputDir (original OutputDir
 * given in input.par )*/
extern char FinalOutputDir[512];
extern char FileNameGalaxies[512];
extern char SimulationDir[512];
extern char FileWithOutputRedshifts[512];

extern char FileWithZList[512];
//variable used to scale to a different cosmology
extern char FileWithZList_OriginalCosm[512];
#ifdef MR_PLUS_MRII  //OPTION for MCMC
extern char FileWithZList_MR[512];
extern char FileWithZList_OriginalCosm_MR[512];
extern char FileWithZList_MRII[512];
extern char FileWithZList_OriginalCosm_MRII[512];
#endif

extern double ScalePos;
extern double ScaleMass;

#ifdef SPECIFYFILENR
extern char   FileNrDir[512];
extern int    ListInputFilrNr[111];
#endif

extern int TotHalos;
extern int TotGalaxies[NOUT];
extern int *TreeNgals[NOUT];

extern int *FirstHaloInSnap;

extern int *TreeNHalos;
extern int *TreeFirstHalo;

extern void *TreeAuxData;


extern double MaxMemSize;

extern size_t AllocatedBytes;
extern size_t HighMarkBytes;
extern size_t FreeBytes;

extern int ThisTask, NTask;

#ifdef GALAXYTREE
extern int GalCount;
extern int TotGalCount;
#ifdef NORMALIZEDDB
extern int TotGalSFHBinCount;
#endif
#endif

/* Cosmological parameters */
extern double BaryonFrac;
extern double Sigma8;
extern double Omega;
extern double OmegaLambda;
extern double Hubble_h;
extern double Omega_OriginalCosm;
extern double OmegaLambda_OriginalCosm;
extern double Hubble_h_OriginalCosm;
//SIMULATION RELATED
extern double PartMass;
extern double BoxSize;
extern double PartMass_OriginalCosm;
extern double BoxSize_OriginalCosm;
#ifdef MR_PLUS_MRII  //OPTION for MCMC
extern double PartMass_MR;
extern double BoxSize_MR;
extern double PartMass_OriginalCosm_MR;
extern double BoxSize_OriginalCosm_MR;
extern double PartMass_MRII;
extern double BoxSize_MRII;
extern double PartMass_OriginalCosm_MRII;
extern double BoxSize_OriginalCosm_MRII;
#endif


/* flags */
extern int StarFormationModel;
#ifdef H2_AND_RINGS
extern int H2FractionRecipe;
extern int SFRtdyn;
#endif
extern int FeedbackReheatingModel;
extern int FeedbackReheatingDeansityScaling;
extern int FeedbackEjectionModel;
extern int FeedbackEagleScaling;
extern int FateOfSatellitesGas;
extern int ReIncorporationModel;
#ifdef EXCESS_MASS
extern int InfallModel;
#endif
extern int ReionizationModel;
extern int BlackHoleGrowth;
extern int AGNRadioModeModel;
extern int DiskRadiusModel;
extern int DiskInstabilityModel;
extern int BHGrowthInDiskInstabilityModel;
extern int HotGasStripingModel;
extern int HotGasOnType2Galaxies;
extern int StarBurstModel;
extern int BulgeFormationInMinorMergersOn;
extern int MetallicityOption;

/* parameters */
extern double Reionization_z0;
extern double Reionization_zr;
extern double Yield;
extern double RecycleFraction;
extern double ThreshMajorMerger;
extern double MergerTimeMultiplier;
extern double RamPressureStrip_CutOffMass;
extern double RamPressureRadiusThreshold;
extern double SfrEfficiency;
extern double SfrColdCrit;
extern double SfrBurstEfficiency;
extern double SfrBurstSlope;
extern double AgnEfficiency;
extern double BlackHoleGrowthRate;
extern double BlackHoleDisruptGrowthRate;
extern double BlackHoleSeedMass;
extern double BlackHoleAccretionRate;
extern double BlackHoleCutoffVelocity;
extern double FeedbackReheatingEpsilon;
extern double ReheatPreVelocity;
extern double ReheatSlope;
extern double FeedbackEjectionEfficiency;
extern double EjectPreVelocity;
extern double EjectSlope;
extern double ReIncorporationFactor;
extern double ReincZpower;
extern double ReincVelocitypower;
extern double FracZSNIItoHot;
extern double FracZSNIatoHot;
extern double FracZAGBtoHot;
extern double A_FACTOR;
#ifdef H2_AND_RINGS
extern double Clumpingfactor;
extern double Warmphasefactor;
extern double GasInflowVel;
#endif
#ifdef FEEDBACK_COUPLED_WITH_MASS_RETURN
extern double EnergySNcode, EnergySN;
extern double EnergySNIIcode, EnergySNII;
extern double EnergySNIacode, EnergySNIa;
extern double EnergyAGBcode, EnergyAGB;
#else
extern double EnergySNcode, EnergySN;
#endif
extern double EtaSNcode, EtaSN;

#ifdef DETAILED_DUST
extern double Dust_tExch;
extern double Dust_tAcc0;
extern double Cmax_CO;
extern double FracDustSNIItoHot;
extern double FracDustSNIatoHot;
extern double FracDustAGBtoHot;
#endif

#ifdef H2_AND_RINGS
extern double RingRadius[RNUM];
extern double RingArea[RNUM];
extern double InverseRingArea[RNUM];
#endif

extern double
	UnitLength_in_cm,
	UnitTime_in_s,
	UnitVelocity_in_cm_per_s,
	UnitMass_in_g,
	RhoCrit,
	UnitPressure_in_cgs,
	UnitDensity_in_cgs,
	UnitCoolingRate_in_cgs,
	UnitEnergy_in_cgs,
	UnitTime_in_Megayears, //Using time as stored in the code, this gives Myr/h
	UnitTime_in_years,
	G,
	Hubble,
	a0, ar;

//extern int NUM_ELEMENTS;

extern int ListOutputSnaps[NOUT];
extern float ListOutputRedshifts[NOUT];

extern double ZZ[MAXSNAPS];
extern double AA[MAXSNAPS];
//variable used to scale to a different cosmology
extern double AA_OriginalCosm[MAXSNAPS];

extern double Age[MAXSNAPS];

extern int    Zlistlen;

extern gsl_rng *random_generator;


extern int    NumMergers;


/*  tabulated stuff */
#ifdef STAR_FORMATION_HISTORY
/* SFH_ is the reference structure for storing the star formation histories in
 * logarithmic bins. It is computed in init.c generating a binning structure for
 * each snapshot/time step. In the code galaxy structures are adjusted with respect
 * to this structure at each step. */
extern double SFH_t[MAXSNAPS][STEPS][SFH_NBIN]; //Time to present at the lower edge of the bin (code units)
extern double SFH_dt[MAXSNAPS][STEPS][SFH_NBIN]; //Time width of the bin (code units)
extern int SFH_Nbins[MAXSNAPS][STEPS][SFH_NBIN]; //Number of bins merged in each bin (only useful for the merging algorithm)
extern int SFH_ibin[MAXSNAPS][STEPS]; //Last active bin
#ifdef DETAILED_METALS_AND_MASS_RETURN
extern double tau_t[STEPS*MAXSNAPS]; //Time-to-z=0 of every timestep in the code. (Used for SNe rates in yield_integrals.c)
extern double tau_dt[STEPS*MAXSNAPS];//Width of every timestep in the code. (Used for SNe rates in yield_integrals.c)
#endif
#endif //STAR_FORMATION_HISTORY

#ifdef DETAILED_METALS_AND_MASS_RETURN
//#define NUM_CHANNELS 3 //Not needed. Can use NUM_METAL_CHANNELS (defined in h_metals.h) instead. Number of channels through which elements are assumed to be ejected: SNe-II, SNe-Ia, AGB.
#ifdef BINARYC
#define BC_TIME_NUM 52 //Number of time bins in the default binary_c ensemble outputs
#define BC_Z_NUM 6 //For old binary_c ensembles: BC_Z_NUM 6
#define BC_LOGT_BINWIDTH 0.1 //Width of binary_c log(t) timebins in dex
#define BC_SNIA_SCALE_FACTOR 1. //10. //1. //Default: 1.0 //An arbitrary scale factor to multiply the binary_c SNIa masses, metals, yields, and rates by, in order to match observed z~0 SN-Ia rates from Graur+18.
#endif //BINARYC

#if !defined(BINARYC) || defined(DETAILED_DUST)
#define LIFETIME_MASS_NUM 150
#define LIFETIME_Z_NUM 6
#endif //!defined(BINARYC) || defined(DETAILED_DUST)

#ifndef BINARYC
//Number of interpolated points within the mass ranges for the four types of yield table:
#define LIFETIME_MASS_NUM 150
#define LIFETIME_Z_NUM 6
#define AGB_MASS_NUM 59 //55 //ROB: 59, when going from 0.85 to 7 Msun
#define AGB_Z_NUM 3
#ifdef PORTINARI
#define SNII_MASS_NUM 85  //ROB: 85, from 6 <= M[Msun] <= 120. Change SNII_MIN_MASS and SNII_MAX_MASS for shorter ranges. (NB: 85 is the number of interpolated steps for the max range of 6-120 Msun found in the SNII yield table files in /YieldTables. This full max range is always read in, but SNII_MIN_MASS and SNII_MAX_MASS limit the range actually used in yields_integrals.c.)
#define SNII_Z_NUM 5
#endif //PORTINARI
#ifdef CHIEFFI
#define SNII_MASS_NUM 81 //ROB: 56 if 7 <= M[Msun] <= 50. 81 if 7 <= M[Msun] <= 120. (NB: You can set SNII_MASS_NUM 81, and SNII_MAX_MASS 50. But DON'T put SNII_MASS_NUM > 81 ever!)
#define SNII_Z_NUM 6
#endif //CHIEFFI
#define SNIA_MASS_NUM 83 //48 //Number increased after extending range to cover M2 masses (07-02-12)
#endif //BINARYC

//Set the depth of the metallicity dimension (ARRAY_Z_NUM) used for all loops/arrays generated in yields_integrals.c to either the binary_c metallicity grid or original lifetime metallicity grid:
#ifdef BINARYC
#define ARRAY_Z_NUM BC_Z_NUM
#else //BINARYC
#define ARRAY_Z_NUM LIFETIME_Z_NUM
#endif //BINARYC

/*
//Mass ranges for the different modes of ejection:
#if !defined(BINARYC) || defined(DETAILED_DUST) //These are needed for calculating AGB dust ejecta masses when DETAILED_DUST is on, even when using BINARYC chemical element yields.
#define IMF_MIN_MASS 0.1 //NOTE: This is not used anywhere explicitly in the code yet. Instances of the minimum star mass *may* be implicit/hard-coded as 0.1 Msun in other parts of L-Galaxies. (19-05-20)
#define IMF_MAX_MASS 120.0 //40.0 //Maximum star mass assumed to exist when normalising the Chabrier IMF in yields_read_tables.c and obtaining KALPHA and F316 in yields_integrals.c. This new parameter now allows max. star mass to be different (i.e. higher) than SNII_MAX_MASS. (18-05-20)
#define AGB_MIN_MASS 0.85
#define AGB_MAX_MASS 7.0 //6.0
#endif //!defined(BINARYC) || defined(DETAILED_DUST)
*/

#ifdef BINARYC //These are needed for calculating AGB dust ejecta masses when DETAILED_DUST is on, even when using BINARYC chemical element yields.
#define IMF_MIN_MASS 0.1 //NOTE: This is not used anywhere explicitly in the code yet. Instances of the minimum star mass *may* be implicit/hard-coded as 0.1 Msun in other parts of L-Galaxies. (19-05-20)
#define IMF_MAX_MASS 80.0 //120.0 //Maximum star mass assumed to exist when normalising the Chabrier IMF to get AGB dust yields in dust_read_tables.c
#else //BINARYC
#define IMF_MIN_MASS 0.1 //NOTE: This is not used anywhere explicitly in the code yet. Instances of the minimum star mass *may* be implicit/hard-coded as 0.1 Msun in other parts of L-Galaxies. (19-05-20)
#define IMF_MAX_MASS 120.0 //120.0 //40.0 //Maximum star mass assumed to exist when normalising the Chabrier IMF in yields_read_tables.c and dust_read_tables.c, and when obtaining KALPHA and F316 in yields_integrals.c.
#endif //BINARYC

#if !defined(BINARYC) || defined(DETAILED_DUST) //These are needed for calculating AGB dust ejecta masses when DETAILED_DUST is on, even when using BINARYC chemical element yields.
#define AGB_MIN_MASS 0.85
#define AGB_MAX_MASS 7.0 //6.0
#endif
#ifndef BINARYC
#define SNIA_MIN_MASS 3.0
#define SNIA_MAX_MASS 16.0
#define SNIA_MIN_TIME 35.0*1.0e6
#define SNIA_MAX_TIME 21.0*1.0e9
#ifdef PORTINARI
#define SNII_MIN_MASS 7.0 //6.0
#define SNII_MAX_MASS 120.0 //120.0 //25.0 //NOTE: This is no longer the same as IMF_MAX_MASS (see above). (18-05-20)
#endif //PORTINARI
#ifdef CHIEFFI
#define SNII_MIN_MASS 7.0
#define SNII_MAX_MASS 120.0 //50.0 //NOTE: This is no longer the same as IMF_MAX_MASS (see above). (18-05-20)
#endif //CHIEFFI
#endif //BINARYC

int ELETOBIGCOUNTA;
int FRACCOUNTA;

//Arrays that yield tables are written to:
#ifdef BINARYC
//char bcElements[11][5];
double bcTimes[BC_TIME_NUM];
double bcMetallicities[BC_Z_NUM];
double bcAGBEjectedMasses[BC_Z_NUM][BC_TIME_NUM];
double bcSNIIEjectedMasses[BC_Z_NUM][BC_TIME_NUM];
double bcSNIaEjectedMasses[BC_Z_NUM][BC_TIME_NUM];
double bcAGBTotalMetals[BC_Z_NUM][BC_TIME_NUM];
double bcSNIITotalMetals[BC_Z_NUM][BC_TIME_NUM];
double bcSNIaTotalMetals[BC_Z_NUM][BC_TIME_NUM];
double bcAGBYields[BC_Z_NUM][NUM_ELEMENTS][BC_TIME_NUM];
double bcSNIIYields[BC_Z_NUM][NUM_ELEMENTS][BC_TIME_NUM];
double bcSNIaYields[BC_Z_NUM][NUM_ELEMENTS][BC_TIME_NUM];
double bcAGBRates[BC_Z_NUM][BC_TIME_NUM];
double bcSNIIRates[BC_Z_NUM][BC_TIME_NUM];
double bcSNIaRates[BC_Z_NUM][BC_TIME_NUM];
#endif //BINARYC

//#ifdef DETAILED_DUST
#if !defined(BINARYC) || defined(DETAILED_DUST)
double lifetimeMasses[LIFETIME_MASS_NUM];
double lifetimeMetallicities[LIFETIME_Z_NUM];
double lifetimes[LIFETIME_Z_NUM][LIFETIME_MASS_NUM];
#endif //!defined(BINARYC) || defined(DETAILED_DUST)

#ifndef BINARYC
double AGBMasses[AGB_MASS_NUM]; //Initial star masses [Msun]
double AGBMetallicities[AGB_Z_NUM]; //Initial star metallicities [Msun]
double AGBEjectedMasses[AGB_Z_NUM][AGB_MASS_NUM]; //Total mass ejected [Msun]
double AGBTotalMetals[AGB_Z_NUM][AGB_MASS_NUM]; //Total metal YIELD ejected [Msun]
double AGBYields[AGB_Z_NUM][NUM_ELEMENTS][AGB_MASS_NUM]; //YIELD ejected, for each element [Msun]
double SNIIMasses[SNII_MASS_NUM];
double SNIIMetallicities[SNII_Z_NUM];
double SNIIEjectedMasses[SNII_Z_NUM][SNII_MASS_NUM];
double SNIITotalMetals[SNII_Z_NUM][SNII_MASS_NUM];
double SNIIYields[SNII_Z_NUM][NUM_ELEMENTS][SNII_MASS_NUM];
#ifndef DTD
double SNIaMasses[SNIA_MASS_NUM];
double SNIaEjectedMasses[SNIA_MASS_NUM];
double SNIaTotalMetals[SNIA_MASS_NUM];
double SNIaYields[42][SNIA_MASS_NUM];
#else //DTD
double SNIaYields[42];
#endif //DTD
#endif //BINARYC

//Integrated yields arrays:
double NormSNIIMassEjecRate[STEPS*MAXSNAPS][SFH_NBIN][ARRAY_Z_NUM];
double NormSNIIMetalEjecRate[STEPS*MAXSNAPS][SFH_NBIN][ARRAY_Z_NUM];
double NormAGBMassEjecRate[STEPS*MAXSNAPS][SFH_NBIN][ARRAY_Z_NUM];
double NormAGBMetalEjecRate[STEPS*MAXSNAPS][SFH_NBIN][ARRAY_Z_NUM];
double NormSNIaMassEjecRate[STEPS*MAXSNAPS][SFH_NBIN][ARRAY_Z_NUM];
double NormSNIaMetalEjecRate[STEPS*MAXSNAPS][SFH_NBIN][ARRAY_Z_NUM];
#ifdef INDIVIDUAL_ELEMENTS
double NormSNIIYieldRate[STEPS*MAXSNAPS][SFH_NBIN][ARRAY_Z_NUM][NUM_ELEMENTS];
double NormAGBYieldRate[STEPS*MAXSNAPS][SFH_NBIN][ARRAY_Z_NUM][NUM_ELEMENTS];
double NormSNIaYieldRate[STEPS*MAXSNAPS][SFH_NBIN][ARRAY_Z_NUM][NUM_ELEMENTS];
#endif //INDIVIDUAL_ELEMENTS

#ifdef WRITE_YIELD_DATA
double solarAbundMassRatios[11]; //Array to store solar abundance ratios from Asplund+09 (see yields_read_tables.c)
//Yield arrays for a 1Msun burst in the first SFH minibin (for writing to a text file from yields_integrals.c):
double NormSNIIYieldRate_burst[STEPS*MAXSNAPS][ARRAY_Z_NUM][NUM_ELEMENTS];
double NormSNIaYieldRate_burst[STEPS*MAXSNAPS][ARRAY_Z_NUM][NUM_ELEMENTS];
double NormAGBYieldRate_burst[STEPS*MAXSNAPS][ARRAY_Z_NUM][NUM_ELEMENTS];
//Arrays for total mass ejected per Msun formed from every SFH bin for each timestep [Msun/Msun]:
double TotNormSNIIEjMass[STEPS*MAXSNAPS][ARRAY_Z_NUM];
double TotNormSNIaEjMass[STEPS*MAXSNAPS][ARRAY_Z_NUM];
double TotNormAGBEjMass[STEPS*MAXSNAPS][ARRAY_Z_NUM];
//Arrays for total metals mass ejected per Msun formed from every SFH bin for each timestep [Msun/Msun]:
double TotNormSNIIMetMass[STEPS*MAXSNAPS][ARRAY_Z_NUM];
double TotNormSNIaMetMass[STEPS*MAXSNAPS][ARRAY_Z_NUM];
double TotNormAGBMetMass[STEPS*MAXSNAPS][ARRAY_Z_NUM];
//Arrays for total SN number per Msun formed from every SFH bin for each timestep [Msun/Msun]:
double TotNormSNIINum[STEPS*MAXSNAPS][ARRAY_Z_NUM];
double TotNormSNIaNum[STEPS*MAXSNAPS][ARRAY_Z_NUM];
double TotNormAGBNum[STEPS*MAXSNAPS][ARRAY_Z_NUM];
#endif //WRITE_YIELD_DATA

//Arrays used to plot SNe rates from SFH bins (yield_integrals.c):
double TheSFH[SFH_NBIN];
double SNII_Rate[STEPS*MAXSNAPS][ARRAY_Z_NUM]; //Rate of SNe-II exploding in this timestep [units: 1/yr]
double SNIa_Rate[STEPS*MAXSNAPS][ARRAY_Z_NUM];
double AGB_Rate[STEPS*MAXSNAPS][ARRAY_Z_NUM];
double Mi_lower_lastTS[SFH_NBIN][ARRAY_Z_NUM][2]; //Used in yields_integrals.c to calculate the number of SNe-II from a 1Msun SP born in the first timestep. ONLY WORKS FOR SFHBIN[0]!. 2-element array (for every SFH bin at every metallicity) which stores the lower mass threshold for the 1st minibin of the previous timestep. This is then used as the upper mass threshold of the last minibin in the current timestep, to prevent any mass range over/under-lapping during integration due to dt and dt_SFH changing relative to each other over time. (20-05-20)

//Temp SNIa arrays for testing revamp of yields_integrals.c (08-04-22):
double SNIa_Rate2[STEPS*MAXSNAPS][ARRAY_Z_NUM];
double NormSNIaNum2[STEPS*MAXSNAPS][SFH_NBIN][ARRAY_Z_NUM];
double NormSNIaMassEjecRate2[STEPS*MAXSNAPS][SFH_NBIN][ARRAY_Z_NUM];
double NormSNIaMassEjecRateOld[STEPS*MAXSNAPS][SFH_NBIN][ARRAY_Z_NUM];
double NormSNIaMetalEjecRate2[STEPS*MAXSNAPS][SFH_NBIN][ARRAY_Z_NUM];
double NormSNIaYieldRate2[STEPS*MAXSNAPS][SFH_NBIN][ARRAY_Z_NUM][NUM_ELEMENTS];

//Arrays used to plot SNe rates from SFH-timesteps (calc_SNe_rates.c):
double TheSFH2[STEPS*MAXSNAPS];
double SNIIRate2[STEPS*MAXSNAPS][ARRAY_Z_NUM];
double SNIaRate2[STEPS*MAXSNAPS][ARRAY_Z_NUM];
double AGBRate2[STEPS*MAXSNAPS][ARRAY_Z_NUM];

//Arrays used to calculate SN rates in-code (i.e. in model_yields.c):
double NormSNIINum[STEPS*MAXSNAPS][SFH_NBIN][ARRAY_Z_NUM]; //Number of SNe-II exploding (per Msun formed) in this timestep per Msun [units: # / Msun]. (It's normalised this way because we have normalised the Chabrier IMF to 1Msun already).
double NormSNIaNum[STEPS*MAXSNAPS][SFH_NBIN][ARRAY_Z_NUM]; //Number of SNe-Ia exploding (per Msun formed) in this timestep per Msun [units: # / Msun]. (It's normalised this way because we have normalised the Chabrier IMF to 1Msun already).
double NormAGBNum[STEPS*MAXSNAPS][SFH_NBIN][ARRAY_Z_NUM]; //Number of AGB stars dying (per Msun formed) in this timestep per Msun [units: # / Msun]. (It's normalised this way because we have normalised the Chabrier IMF to 1Msun already).

//Arrays used to store disc SN rates and newly-ejected element masses (for use in model_dust_yields.c):
#ifdef DETAILED_DUST
#ifdef H2_AND_RINGS
double DiskSNIIRate_current_ts[RNUM]; //SNII rate in the disc in this current timestep [units: 1/(time_code_units)]
double DiskSNIaRate_current_ts[RNUM]; //SNIa rate in the disc in this current timestep [units: 1/(time_code_units)]
double DiskAGBRate_current_ts[RNUM]; //AGB rate in the disc in this current timestep [units: 1/(time_code_units)]
#ifdef INDIVIDUAL_ELEMENTS
double DiskSNIIAllElements_ts[RNUM][NUM_ELEMENTS]; //Total element masses ejected by disc SNe-II in this timestep [units: Msun]
double DiskSNIaAllElements_ts[RNUM][NUM_ELEMENTS]; //Total element masses ejected by disc SNe-Ia in this timestep [units: Msun]
double DiskAGBAllElements_ts[RNUM][NUM_ELEMENTS]; //Total element masses ejected by disc AGBs in this timestep [units: Msun]
#ifdef DUST_HOTGAS
double BulgeSNIIAllElements_ts[RNUM][NUM_ELEMENTS]; //Total element masses ejected by bulge SNe-II in this timestep [units: Msun]
double BulgeSNIaAllElements_ts[RNUM][NUM_ELEMENTS]; //Total element masses ejected by bulge SNe-Ia in this timestep [units: Msun]
double BulgeAGBAllElements_ts[RNUM][NUM_ELEMENTS]; //Total element masses ejected by bulge AGBs in this timestep [units: Msun]
double ICMSNIIAllElements_ts[RNUM][NUM_ELEMENTS]; //Total element masses ejected by halo SNe-II in this timestep [units: Msun]
double ICMSNIaAllElements_ts[RNUM][NUM_ELEMENTS]; //Total element masses ejected by halo SNe-Ia in this timestep [units: Msun]
double ICMAGBAllElements_ts[RNUM][NUM_ELEMENTS]; //Total element masses ejected by halo AGBs in this timestep [units: Msun]
#endif //DUST_HOTGAS
#endif //INDIVIDUAL_ELEMENTS
#else //H2_AND_RINGS
double DiskSNIIRate_current_ts; //SNII rate in the disc in this current timestep [units: 1/(time_code_units)]
double DiskSNIaRate_current_ts; //SNIa rate in the disc in this current timestep [units: 1/(time_code_units)]
double DiskAGBRate_current_ts; //AGB rate in the disc in this current timestep [units: 1/(time_code_units)]
#ifdef INDIVIDUAL_ELEMENTS
double DiskSNIIAllElements_ts[NUM_ELEMENTS]; //Total element masses ejected by disc SNe-II in this timestep [units: Msun]
double DiskSNIaAllElements_ts[NUM_ELEMENTS]; //Total element masses ejected by disc SNe-Ia in this timestep [units: Msun]
double DiskAGBAllElements_ts[NUM_ELEMENTS]; //Total element masses ejected by disc AGBs in this timestep [units: Msun]
#ifdef DUST_HOTGAS
double BulgeSNIIAllElements_ts[NUM_ELEMENTS]; //Total element masses ejected by bulge SNe-II in this timestep [units: Msun]
double BulgeSNIaAllElements_ts[NUM_ELEMENTS]; //Total element masses ejected by bulge SNe-Ia in this timestep [units: Msun]
double BulgeAGBAllElements_ts[NUM_ELEMENTS]; //Total element masses ejected by bulge AGBs in this timestep [units: Msun]
double ICMSNIIAllElements_ts[NUM_ELEMENTS]; //Total element masses ejected by halo SNe-II in this timestep [units: Msun]
double ICMSNIaAllElements_ts[NUM_ELEMENTS]; //Total element masses ejected by halo SNe-Ia in this timestep [units: Msun]
double ICMAGBAllElements_ts[NUM_ELEMENTS]; //Total element masses ejected by halo AGBs in this timestep [units: Msun]
#endif //DUST_HOTGAS
#endif //INDIVIDUAL_ELEMENTS
#endif //H2_AND_RINGS
#endif //DETAILED_DUST

//Params now made global (for dust model) so that they can be called in main.c for SN feedback while also being updated in model_yields.c (08-11-21):
double TotalMassReturnedToColdDiskGas;
double TotalMassReturnedToHotGas;
double TotalMetalsReturnedToHotGas;
double TotalMassReturnedToColdDiskGasr[RNUM];

//Actual wind efficiencies used in model_yields.c (required to be global for use in model_dust_yields.c):
double fwind_SNII; //Required for SN-II metal-rich wind implementation
double fwind_SNIa; //Required for SN-Ia metal-rich wind implementation
double fwind_AGB; //Required for AGB metal-rich wind implementation

//IMF parameters (for chemical enrichment):
#define IMF_SLOPE 2.3 //2.6 //2.0 //2.15 //High-mass slope of Chabrier IMF. (2.3 = normal Chabrier IMF. <2.3 = top-heavy, >2.3 = bottom-heavy.)
//For consistency with yields_integrals.c, IMF_SLOPE must equal either: 2.0, 2.15, 2.3, or 2.6.

//SNIa parameters:
//#define A_FACTOR 0.04 //NOW DEFINED IN THE INPUT.PAR FILES (11-11-21) //0.028 //0.035 //Fraction of mass from all objects between SNIa_MIN_MASS and SNIA_MAX_MASS that comes from SN-Ia. //0.028 preferred in Yates+13.
//#define FRAC2HOT 0.9 //Fraction of material released by disk stars that goes straight into the HotGas. Rest goes in ColdGas.
#ifdef DTD
//#define KALPHA 1.4765 //1.59203 //Now set in yield_integrals.c
//#define	F316 0.0384 //Integral of the IMF (by number) from 3.0 - 16.0 Msun //Now set in yield_integrals.c
#define SNIAEJECMASS 1.2300971 //Total mass (and total metals) ejected by a SNIa explosion in Msun //Value form original yield table (42 elements): 1.3740855. //Value when only considering 11 elements: 1.2300971
#ifdef BIMODALDTD
	//#define DTD_NORM 0.903206 //(26Myrs - 21Gyrs)
	//#define DTD_NORM 0.896668 //For P98 Z=0.02 lifetimes (35Myrs - 17Gyrs)
	#define DTD_NORM 0.900348 //(35Myrs - 21Gyrs)
#endif //BIMODALDTD
#ifdef CUSTOMDTD
	//#define DTD_NORM 0.524836 //(26Myrs - 21Gyrs)
	//#define DTD_NORM 0.606746 //For P98 Z=0.02 lifetimes (35Myrs - 17Gyrs, for ~32% in prompt component)
	#define DTD_NORM 0.610431 //(35Myrs - 21Gyrs, for ~32% in prompt component)
#endif //CUSTOMDTD
#ifdef GAUSSIANDTD
	#define DTD_NORM = 1.0 //For P98 Z=0.02 lifetimes (35Myrs - 17Gyrs) and (35Myrs - 21Gyrs)
	#define TAUCHARAC 1.0 //Characteristic delay time for SNe-Ia (i.e. peak of Gaussian distribution) in Gyrs //default: 2.0
	#define SIGMA_TD 0.2*TAUCHARAC //0.2 for narrow-DTD, 0.5 for wide_DTD
#endif //GAUSSIANDTD
#ifdef POWERLAWDTD
	//#define DTD_NORM 7.21863 // (26Myrs - 21Gyrs)
	#define DTD_NORM 6.72574 //6.72544 // (35Myrs - 21Gyrs)
	//#define DTD_NORM 6.56087 //For P98 Z=0.02 lifetimes (35Myrs - 17Gyrs)
	//#define DTD_NORM 6.22432 // (35Myrs - 11Gyrs)
	//#define DTD_NORM 6.02197 // (50Myrs - 17Gyrs)
	//#define DTD_NORM 6.35503 // (40Myrs - 17Gyrs)
	#define DTD_SLOPE -1.12 //Slope of power law, according to Maoz et al. (2012)
#endif //POWERLAWDTD
#ifdef RUITERDTD
	//#define DTD_NORM 1.09545 //(26Myrs - 21Gyrs)
	#define DTD_NORM 1.09545 //(35Myrs - 21Gyrs)
	//#define DTD_NORM 1.08422 //For P98 Z=0.02 lifetimes (35Myrs - 17Gyrs)
	#define TAUCHARAC 0.5 //Peak of Gaussian (prompt) component [in Gyrs]
	#define SIGMA_TD 0.2*TAUCHARAC //Width of Gaussian (prompt) component
	#define DTD_SLOPE -2.0 //Slope of power law (delayed) component (see Ruiter et al. 2012)
#endif //RUITERDTD
#endif //DTD

#endif //DETAILED_METALS_AND_MASS_RETURN

#ifdef DETAILED_DUST
	// Numbers for reading in dust tables in dustyields_read_tables.c
	#define AGB_DUST_MASS_NUM 28
	#define AGB_DUST_METAL_NUM 4
	#define AGB_DUST_TYPE_NUM 11

	float AGBDustMasses[AGB_DUST_MASS_NUM]; //Initial star masses [Msun]
	float AGBDustMetallicities[AGB_DUST_METAL_NUM]; //Initial star metallicities [Msun]
	float AGBDustCreated[AGB_DUST_METAL_NUM][AGB_DUST_MASS_NUM][AGB_DUST_TYPE_NUM]; //Total mass ejected [Msun]

	//float NormAGBDustYieldRate[(STEPS*MAXSNAPS)][SFH_NBIN][LIFETIME_Z_NUM][AGB_DUST_TYPE_NUM];
	float NormAGBDustYieldRate[(STEPS*MAXSNAPS)][SFH_NBIN][ARRAY_Z_NUM][AGB_DUST_TYPE_NUM];

	//Variables required for partition_gas_and_dust_elements():
	double ColdGasDiff_elements_old[NUM_ELEMENTS]; //Stores the element masses contained in the "diffuse gas" sub-component of ColdGas from the previous timestep. For use in shuffle_ISM().
	double ColdGasClouds_elements_old[NUM_ELEMENTS]; //Same as above, but for the "molecular clouds" sub-component.
	double ColdGasDiff_elements_change[NUM_ELEMENTS]; //Stores the difference in element masses between the current and previous timestep for the "diffuse gas" sub-component of ColdGas. For use in shuffle_ISM().
	double ColdGasClouds_elements_change[NUM_ELEMENTS]; //Same as above, but for the "molecular clouds" sub-component.
#ifdef H2_AND_RINGS
	double ColdGasDiffRings_elements_old[RNUM][NUM_ELEMENTS];
	double ColdGasCloudsRings_elements_old[RNUM][NUM_ELEMENTS];
	double ColdGasDiffRings_elements_change[RNUM][NUM_ELEMENTS];
	double ColdGasCloudsRings_elements_change[RNUM][NUM_ELEMENTS];
	// For storing metallicity in model_yields.c
	int Zi_disk_saved[RNUM][SFH_NBIN];
	double Zi_disk_disp_saved[RNUM][SFH_NBIN];
	int Zi_bulge_saved[RNUM][SFH_NBIN];
	double Zi_bulge_disp_saved[RNUM][SFH_NBIN];
	int Zi_ICM_saved[RNUM][SFH_NBIN];
	double Zi_ICM_disp_saved[RNUM][SFH_NBIN];
	// Variables to hold amount of metals created in prev time step in model_yields.c
	float SNII_prevstep_Cold_Cb[RNUM][SFH_NBIN];
	float SNII_prevstep_Cold_Si[RNUM][SFH_NBIN];
	float SNII_prevstep_Cold_Fe[RNUM][SFH_NBIN];
	float SNIa_prevstep_Cold_Fe[RNUM][SFH_NBIN];
	float SNII_prevstep_Hot_bulge_Cb[RNUM][SFH_NBIN];
	float SNII_prevstep_Hot_bulge_Si[RNUM][SFH_NBIN];
	float SNII_prevstep_Hot_bulge_Fe[RNUM][SFH_NBIN];
	float SNIa_prevstep_Hot_bulge_Fe[RNUM][SFH_NBIN];
	float SNII_prevstep_Hot_ICM_Cb[RNUM][SFH_NBIN];
	float SNII_prevstep_Hot_ICM_Si[RNUM][SFH_NBIN];
	float SNII_prevstep_Hot_ICM_Fe[RNUM][SFH_NBIN];
	float SNIa_prevstep_Hot_ICM_Fe[RNUM][SFH_NBIN];
#else
	int Zi_disk_saved[SFH_NBIN];
	double Zi_disk_disp_saved[SFH_NBIN];
	int Zi_bulge_saved[SFH_NBIN];
	double Zi_bulge_disp_saved[SFH_NBIN];
	int Zi_ICM_saved[SFH_NBIN];
	double Zi_ICM_disp_saved[SFH_NBIN];
	float SNII_prevstep_Cold_Cb[SFH_NBIN];
	float SNII_prevstep_Cold_Si[SFH_NBIN];
	float SNII_prevstep_Cold_Fe[SFH_NBIN];
	float SNIa_prevstep_Cold_Fe[SFH_NBIN];
	float SNII_prevstep_Hot_bulge_Cb[SFH_NBIN];
	float SNII_prevstep_Hot_bulge_Si[SFH_NBIN];
	float SNII_prevstep_Hot_bulge_Fe[SFH_NBIN];
	float SNIa_prevstep_Hot_bulge_Fe[SFH_NBIN];
	float SNII_prevstep_Hot_ICM_Cb[SFH_NBIN];
	float SNII_prevstep_Hot_ICM_Si[SFH_NBIN];
	float SNII_prevstep_Hot_ICM_Fe[SFH_NBIN];
	float SNIa_prevstep_Hot_ICM_Fe[SFH_NBIN];
#endif //H2_AND_RINGS
#endif //DETAILED_DUST

#ifdef COMPUTE_SPECPHOT_PROPERTIES
// SSP PHOT_TABLES - magnitues of starburst population as a function of age

#ifdef M05
#define SSP_NAGES 220		// Age grid of the SSP tables
#define SSP_NMETALLICITES 4			// Number of Metalicities used
#ifdef SPEC_PHOTABLES_ON_THE_FLY
#define SSP_NLambda 1221
#endif
#endif

#ifdef BC03
#define SSP_NAGES 221		// Age grid of the SSP tables
#define SSP_NMETALLICITES 6			// Number of Metalicities used
#ifdef SPEC_PHOTABLES_ON_THE_FLY
#define SSP_NLambda 1221
#endif
#endif

#ifdef CB07
#define SSP_NAGES 221		// Age grid of the SSP tables
#define SSP_NMETALLICITES 6			// Number of Metalicities used
#ifdef SPEC_PHOTABLES_ON_THE_FLY
#define SSP_NLambda 1238
#endif
#endif

//table containing the Metallicity grid of the SSP tables (converted to log10)
extern double SSP_logMetalTab[SSP_NMETALLICITES];
//table containing the Age grid of the SSP tables (originally in years, converted to log10(internal time units 1e12 Yrs/h))
extern double SSP_logAgeTab[SSP_NAGES];
//table containing redshift (different from the one in the code when scaling to future times)
extern double RedshiftTab[MAXSNAPS];
extern double LumTables[NMAG][SSP_NMETALLICITES][MAXSNAPS][SSP_NAGES];
extern double FilterLambda[NMAG+1];//wavelength of each filter + 1 for V-band

#ifdef SPEC_PHOTABLES_ON_THE_FLY
#define MAX_NLambdaFilter 1000
extern int NLambdaFilter[NMAG];
//VEGA
#define NLambdaVega 3303
#endif

//DUST EXTINCTION
#define ExpTauBCBulge 0.5	// constant extinction for young stars in bulges.
#define MUWIDTH  0.2
#define MUCENTER 0.3
extern long mu_seed;

#endif //COMPUTE_SPECPHOT_PROPERTIES


extern size_t HighMark;

#ifdef UPDATETYPETWO
extern int NtotHalos, TotIds, Nids, TotSnaps, OffsetIDs;
extern int *CountIDs_halo, *OffsetIDs_halo, *CountIDs_snaptree, *OffsetIDs_snaptree;
extern long long *IdList;
extern float *PosList, *VelList;
#endif


extern int Hashbits;
extern int NumWrittenInParallel;
extern double ScaleFactor;	// factor by which to multiply a position to get its ph index (after floring)


#ifdef USE_MEMORY_TO_MINIMIZE_IO
extern char *ptr_auxdata, *ptr_treedata, *ptr_dbids, *ptr_galaxydata, *ptr_galsnapdata[NOUT];
extern size_t offset_auxdata, offset_treedata, offset_dbids;
extern size_t offset_galaxydata, maxstorage_galaxydata, filled_galaxydata;
extern size_t offset_galsnapdata[NOUT], maxstorage_galsnapdata[NOUT], filled_galsnapdata[NOUT];
#endif


extern float Reion_z[46],Reion_Mc[46];

extern FILE *tree_file;
extern FILE *treeaux_file;
extern FILE *treedbids_file;
extern FILE *FdGalTree;
extern FILE *FdGalTreeSFH;
extern FILE *FdGalDumps[NOUT];

/*H2 fraction table*/
#ifdef H2_AND_RINGS
//#define RHO_LEN 101
//#define Z_LEN 13
#define RHO_LEN 420
#define Z_LEN 6
extern double H2Fraction[LENZ][LENSIGMAH];
extern double H2Fraction_Zgrid[LENZ];
extern double H2Fraction_SigmaHgrid[LENSIGMAH];
#endif

#ifdef HDF5_OUTPUT
int b[NOUT];
struct GALAXY_OUTPUT galaxy_output_hdf5[NOUT][NRECORDS_APP];
#endif //HDF5_OUTPUT

#ifdef DEBUG
extern FILE *FdGalDebug;
#ifdef H2_AND_RINGS
#define NDebugProps 23
extern char DebugProperties[NDebugProps][100];
#else
#define NDebugProps 11
extern char DebugProperties[NDebugProps][100];
#endif
#endif
