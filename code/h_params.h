//To understand the units in the code read through set_units in init.c!!!
#define  GRAVITY     6.672e-8
#define  SOLAR_MASS  1.989e33
#define  SOLAR_LUM   3.826e33
#define  RAD_CONST   7.565e-15
#define  AVOGADRO    6.0222e23
#define  BOLTZMANN   1.3806e-16
#define  GAS_CONST   8.31425e7
#define  SPEEDOFLIGHT 2.9979e10 // ROB: Changed from "C" to "SPEEDOFLIGHT" (09-04-21)
#define  PLANCK      6.6262e-27
#define  CM_PER_MPC  3.085678e24 // TODO this is read in from input.par
#define  PROTONMASS  1.6726e-24 //[grams]
#define  HUBBLE      3.2407789e-18   /* in h/sec */
#define  MUMH        0.59341*PROTONMASS // Could be a variable but we never change it!
#define	 PI			 3.14159 //ROB: Added on 16-03-22: My mind is blown that this wasn't already here?! How have we not needed pi before? It was only defined in mcmc_halomodel.c

//To understand the units in the code read through set_units in init.c!!!
#define  SEC_PER_MEGAYEAR   3.155e13
#define  SEC_PER_YEAR       3.155e7

#define MIN_ALLOC_NUMBER       1000
#define ALLOC_INCREASE_FACTOR  1.1
#define ALLOC_DECREASE_FACTOR  0.7
#define PRECISION_LIMIT 1.e-7
#define TINY_MASS 1.e-8 //ROB: Is this supposed to be in code units (i.e. 1.e10/Hubble_h)?
#define TINY_LENGTH 1.e-6

#ifdef GALAXYTREE
#undef  NOUT
#define NOUT MAXSNAPS
#endif

//WATCH OUT! In the case of MCMC running both MR and MRII the larger value is used to "allocate" all the arrays
//inside the code its LastDarkMatterSnapShot+1 that defines the extent of the loops
//(in MCMC MR_plus_MRII mode this are not always identical)
#ifdef MRII
#define  MAXSNAPS  68     /* Number of snapshots in the dark matter simulation */
#else

#ifdef PHOENIX
#define  MAXSNAPS  72
#else

#ifdef NIFTY
#define  MAXSNAPS  108
#else

#ifdef CATERPILLAR
#define  MAXSNAPS  320
#else

#define  MAXSNAPS  64  //NORMAL MILLENNIUM

#endif //CATERPILLAR
#endif //NIFTY
#endif //PHOENIX
#endif //else MRII

#define  MAXGALFAC 2.3 /*1.5/2.3 - maximum fraction of satellite without a halo (for memory allocation)*/


#ifdef FAST_TESTING_MODE
#define  STEPS 10
#else
#define  STEPS 20 //60
#endif

#ifdef H2_AND_RINGS
#define RNUM 12          /* radially divide one disk into RNUM */
//#define RNUM 30          /* radially divide one disk into RNUM */
#define WARM_PHASE_FACTOR 1.3 // to use when deriving HI and H2 from total cold gas because 1/3 of it is ionized)
#define LENSIGMAH 40
#define LENZ 6
//#define LENSIGMAH 101
//#define LENZ 13
#endif

#ifdef METALRICHWIND
#ifdef GASDENSITYFWIND
#define NORMGASDENSITY 1.5 //90. //100. //10. //20. //Msun/pc^2 //ISM gas surface density to normalise to when calculating density-dependent direct ejection into HotGas (18-05-18)
#endif
#endif

#ifdef INDIVIDUAL_ELEMENTS
#define PRISTINE_H_FRACTION 0.75 //Fraction of infalling gas assumed to be hydrogen
#define PRISTINE_HE_FRACTION 0.25 //Fraction of infalling gas assumed to be helium
#endif

#ifdef DETAILED_DUST
//Dust rate array dimensionality:
#ifdef FULL_DUST_RATES
#define NUM_COLDGAS_DUST_RATES 5 // [AGB, SNII, SNIa, GrainGrowth, SNShockDestruction] // Needed for HDF5 table creation
#define NUM_HOTGAS_DUST_RATES 4 // [AGB, SNII, SNIa, SputteringDestruction] // Dust production in the bulge and halo, plus destruction in the HotGas via sputtering
#endif //FULL_DUST_RATES

//Element atomic weights:
#define A_H 1.008
#define A_He 4.003
#define A_Cb 12.01 //0
#define A_N 14.007 //1   N considered to be volatile
#define A_O 15.999 //2
#define A_Ne 20.1797 //3  Ne does not form dust, non-reactive
#define A_Mg 24.305 //4
#define A_Si 28.085 //5
#define A_S 32.06 //6   S considered to be volatile
#define A_Ca 40.078 //7
#define A_Fe 55.845 //8

//Taken from Zhukovska2008:
#define A_Sil_dust 121.62 //121.41
#define A_Fe_dust  55.845 //55.85 //Note: Original value had a different precision to A_Fe above.
#define A_SiC_dust 40.10 //SiC dust does not form in the ISM
#define A_Cb_dust  12.01

//Dust condensation efficiency parameters (from Zhukvska+08):
//From Vijayan+19: "These efficiency parameters are defined considering the effects of the reverse shock
//and are therefore smaller than they would be for initial dust condensation."
#define eta_SNII_Sil 0.00035
#define eta_SNII_Fe  0.001
#define eta_SNII_SiC 0.0003
#define eta_SNII_Cb  0.15
#define eta_SNIa_Fe  0.005

//Defining the mass fractions of each dust type that is made up of each element (for use in model_dust_yields.c):
//Note: ROB: Given that the atomic weights above are only to 3 d.c.s, shouldn't we truncate these dust-type fractions to 3 d.c.s too?
//Forsterite (Mg2SiO4):
//#define FORSTERITE_ELEMENT_NUM 3
//#define FORSTERITE_ELEMENTS [4,6,7] //O, Mg, Si
#define FORSTERITE_O_FRAC 0.454874 //i.e. (4*A_O) / ((2*A_Mg)+A_Si+(4*A_O)) = 0.454869
#define FORSTERITE_Mg_FRAC 0.345504
#define FORSTERITE_Si_FRAC 0.199622
//Fayalrite (Fe2SiO4):
#define FAYALITE_O_FRAC 0.314063
#define FAYALITE_Si_FRAC 0.137827
#define FAYALITE_Fe_FRAC 0.548110
//Enstatite (MgSi03):
#define ENSTATITE_O_FRAC 0.478124
#define ENSTATITE_Mg_FRAC 0.243050
#define ENSTATITE_Si_FRAC 0.279768
//Ferrosilite (Fe2Si206):
#define FERROSILITE_O_FRAC 0.363819
#define FERROSILITE_Si_FRAC 0.212884
#define FERROSILITE_Fe_FRAC 0.423297
//Quartz (SiO4):
#define QUARTZ_O_FRAC 0.694998
#define QUARTZ_Si_FRAC 0.305002
//Silicon Carbide (SiC):
#define SILICONCARBIDE_Cb_FRAC 0.299547
#define SILICONCARBIDE_Si_FRAC 0.700453
//Olivine ([Mg,Fe]2SiO4):
//#define OLIVINE_Mg_FRAC
//Silicates:
#define SILICATES_O_FRAC 0.439462 //0.419567
#define SILICATES_Mg_FRAC 0.205116 //0.091053
#define SILICATES_Si_FRAC 0.237635 //0.210432
#define SILICATES_Fe_FRAC 0.117787 //0.278948

#define OLIVINE_Mg_NUMFRAC 0.8 //The average fraction of Mg+Fe assumed to be Mg in Olivine (i.e. parameter x in Zhukovska+08)
#define PYROXENE_Mg_NUMFRAC 0.8 //The average fraction of Mg+Fe assumed to be Mg in Pyroxene (i.e. parameter x in Zhukovska+08)
#define OLIVINE_NUMFRAC 0.32 //The average fraction of Olivine+Pyroxene assumed to be Olivine (i.e. parameter f_ol in Zhukovska+08)
#define HEMATITE_NUMFRAC 0.5 //The average fraction of Iron oxide assumed to be Hematite (rather than Magnetite)

#ifdef DUST_DESTRUCTION
#define M_CLEARED 1200.0 //Msol //Default assumed value for mass of gas that is cleared of dust by an average SNe (see Vijayan+19)
#define M_CLEARED_Cb 1180.0 //[Msun] //From Hu+19 table 2 (for n_H = 0.3 cm^-3)
#define M_CLEARED_Si 1660.0 //[Msun] //From Hu+19 table 2 (for n_H = 0.3 cm^-3)
#define F_SN 0.36 //Dimensionless
//#ifdef DUST_HOTGAS
#if defined(DUST_HOTGAS) || defined(DUST_EJECTEDMASS)
#define sputConst 3.2e-18 //[cm^4 s^-1]
#define m_proton 1.6726219e-24 //Mass of a proton [grams]
#define sputT0 2.e6 //Critical temp above which the sputtering rate flattens [Kelvin]
#define sputOmega 2.5 //Controls the sputtering rate in the low-temperature regime
#define sputCharacTime 0.17 //Characteristic sputtering time [Gyr]
#define CharacGrainRadius 0.1 //Assumed initial dust grain radius [micron]
#define CharacDensity 1.e-27 //Characteristic gas density for sputtering timescale [g/cm^3]
#endif //defined(DUST_HOTGAS) || defined(DUST_EJECTEDMASS)
#endif //DUST_DESTRUCTION
#endif //DETAILED_DUST

#define  ALLOCPARAMETER 50.  /* new definition !!! THIS HAS TO BE 50 !!! DONT EVER EVER EVER CHANGE !!! */

#ifdef STAR_FORMATION_HISTORY
#define SFH_NMERGE 3 //54 //15 // SFH_NMERGE=Nmax+1 (Nmax used in Shamshiri2014)
//#define SFH_NMERGE 10  //  SFH_NMERGE=Nmax+1 (Nmax used in Shamshiri2014)
#ifdef NIFTY
#define SFH_NBIN 22 //  NIFTY - 108 snapshots
#else
#ifdef CATERPILLAR
#define SFH_NBIN 24 //  CATERPILLAR - 320 snapshots
#else
//#define SFH_NBIN 1280 //STEPS*MAXSNAPS
#define SFH_NBIN 20 //370 //110
#endif //NIFTY
#endif //CATERPILLAR
#endif //STAR_FORMATION_HISTORY

