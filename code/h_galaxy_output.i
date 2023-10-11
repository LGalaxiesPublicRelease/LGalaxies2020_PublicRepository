# 0 "./code/h_galaxy_output.h"
# 1 "/cygdrive/c/Users/ry22aas/robyates/Astro/L-Galaxies/L-Galaxies2020_plusBinaries/L-Galaxies2020_plusBinaries_version//"
# 0 "<built-in>"
# 0 "<command-line>"
# 1 "./code/h_galaxy_output.h"
/*
 * Galaxy structure for output.
 *
 * NOTE: due to the way that the HDF5 builder routines work, variables should be commented out
 * with slash-star ... star-slash and not //.  Otherwise they are included in the properties list.
 *
 */
# 41 "./code/h_galaxy_output.h"
#pragma pack(1)
struct GALAXY_OUTPUT {
# 72 "./code/h_galaxy_output.h"
//    int   NumDisr; //***** ROB: Geoff's disruption counter (24-03-20) *****
    int Type; // None //Galaxy type: 0 for central galaxies of a main halo, 1 for central galaxies in sub-halos, 2 for satellites without halo.

    int HaloIndex; // None // ?Unique ID of MPA halo containing this galaxy





    float HaloM_Crit200; // 1e10 Msun/h // M200 cf critical last time this halo was a type 0
# 91 "./code/h_galaxy_output.h"
    int SnapNum; // None //The snapshot number where this galaxy was identified.
    float LookBackTimeToSnap; // yr // The time from a given snapshot to z=0
    float CentralMvir; // 10^10/h Msun // virial mass of background (FOF) halo containing this galaxy
    float CentralRvir; // Mpc/h // Proper[?] R200 cf critical of background (FOF) halo containing this galaxy
    float DistanceToCentralGal[3]; // Mpc/h // Proper[?] components of the distance between this galaxy and the galaxy at the centre of the FoF group.
    float Pos[3]; // 1/h Mpc // Comoving galaxy/subhalo position
    float Vel[3]; // km/s // Galaxy/subhalo peculiar velocity
    int Len; // None // Number of particles in the associated subhalo  
    /* properties of subhalo at the last time this galaxy was a central galaxy */
    float Mvir; // 10^10/h Msun // M200 cf critical of the halo last time galaxy was type 0
    float Rvir; // Mpc/h // R200 cf critical of the halo last time galaxy was type 0
    float Vvir; // km/s // Virial velocity of the halo last time galaxy was type 0
    float Vmax; // km/s //Maximum rotational velocity of the subhalo, or the last value for type 2's galaxies.
    float ColdGasSpin[3]; // Mpc/h km/s // The specific angular momentum of the cold gas disk
    float DiskSpin[3]; // Mpc/h km/s // The specific angular momentum of the stellar disk
# 131 "./code/h_galaxy_output.h"
    /* baryonic reservoirs */
    float ColdGas; // 10^10/h Msun // Mass in cold gas.

    float ColdGasRings[RNUM]; // 1e10 Msun/h //  Mass of clod gas in each annulur ring.
    float H2fraction; // None //  Fraction of ColdGas in the form of H_2
    float H2fractionRings[RNUM]; // None //  H2 fraction within each annular ring.

    float StellarMass; // 10^10/h Msun // Total mass in stars in the disk and the bulge combined
    float DiskMass; // 10^10/h Msun // Mass of stars in the disk
    float BulgeMass; // 10^10/h Msun // Mass of stars in the bulge

    float DiskMassRings[RNUM]; // 1e10 Msun/h // Mass of stars within each annular ring
    float BulgeMassRings[RNUM]; // 1e10 Msun/h // Mass of stars within each annular ring

    float HotGas; // 10^10/h Msun // Mass in hot gas
    /*float ReheatedGas; // 10^10/h Msun // Mass in reheated gas*/
    float EjectedMass; // 10^10/h Msun // Mass in ejected gas



    float BlackHoleMass; // 10^10/h Msun // Mass of central black hole
    /* float BlackHoleGas; // 10^10/h Msun // Mass in BH accretion disk */
    /* ICL magnitude and mass*/
    float ICM; //10^10/h Msun //Total mass in metals in intra-cluster stars, for type 0,1
# 165 "./code/h_galaxy_output.h"
    float MetalsColdGas[NUM_METAL_CHANNELS]; // 10^10/h Msun // Mass in metals in cold gas.

    float MetalsColdGasRings[RNUM][NUM_METAL_CHANNELS]; // 10^10/h Msun // Mass in metals in cold gas in each annular ring

    float MetalsStellarMass[NUM_METAL_CHANNELS]; // 10^10/h Msun // Mass in metals in the disk
    float MetalsDiskMass[NUM_METAL_CHANNELS]; // 10^10/h Msun // Mass in metals in the disk
    float MetalsBulgeMass[NUM_METAL_CHANNELS]; // 10^10/h Msun // Mass in metals in the bulge

    float MetalsDiskMassRings[RNUM][NUM_METAL_CHANNELS]; // 10^10/h Msun // Mass in metals in stars in each annular ring
    float MetalsBulgeMassRings[RNUM][NUM_METAL_CHANNELS]; // 10^10/h Msun // Mass in metals in stars in each annular ring

    float MetalsHotGas[NUM_METAL_CHANNELS]; // 10^10/h Msun // Mass in metals in the hot gas
    /* float MetalsReheatedGas[NUM_METAL_CHANNELS]; // 10^10/h Msun // Mass in metals in the Reheated gas */
    float MetalsEjectedMass[NUM_METAL_CHANNELS]; // 10^10/h Msun // Mass in metals in the ejected gas



    float MetalsICM[NUM_METAL_CHANNELS]; // 10^10/h Msun // Mass in metals in intra-cluster stars, for type 0,1



    float DiskSNIIRate; //yr^-1 //Rate of SN-II explosions from the disc (averaged over a whole snapshot)
    float BulgeSNIIRate; //yr^-1 //Rate of SN-II explosions from the bulge (averaged over a whole snapshot)
    float ICMSNIIRate; //yr^-1 //Rate of SN-II explosions from the halo stars (averaged over a whole snapshot)
    float DiskSNIaRate; //yr^-1 //Rate of SN-Ia explosions from the disc (averaged over a whole snapshot)
    float BulgeSNIaRate; //yr^-1 //Rate of SN-Ia explosions from the bulge (averaged over a whole snapshot)
    float ICMSNIaRate; //yr^-1 //Rate of SN-Ia explosions from the halo stars (averaged over a whole snapshot)
# 216 "./code/h_galaxy_output.h"
    float Sfr; // Msun/yr // Star formation rate
    float SfrInst; // Msun/yr // Instantaneous star formation rate (i.e. the sfr in the final timestep of each snapshot) considering secular SF only (i.e. not including starbursts) //*****ROB*****//

    float SfrRings[RNUM]; // Msun/yr // Star formation rate within each annular ring
    float SfrInstRings[RNUM]; // Msun/yr // Instantaneous star formation rate within each annular ring (i.e. the sfr in the final timestep of each snapshot) considering secular SF only (i.e. not including starbursts) //*****ROB*****//

    float SfrBulge; // Msun/yr // Star formation rate in bulge.



    float BulgeSize; // Mpc/h // Half mass radius of bulge
    float DiskRadius; // Mpc/h // Size of the stellar disk, 3x the scale length.
    float ColdGasRadius; // Mpc/h // Size of the gas disk, 3x the scale length.
    float StellarHalfMassRadius; // Mpc/h // stellar Half mass radius
# 242 "./code/h_galaxy_output.h"
       /* magnitudes in various bands */
# 280 "./code/h_galaxy_output.h"
    float MassWeightAge; //10^9yr //The age of this galaxy weighted by mass of its components.
# 318 "./code/h_galaxy_output.h"
    //All: [H][He][Cb][N][O][Ne][Mg][Si][S][Ca][Fe] or //Only [H][He][O][Mg][Fe]
# 330 "./code/h_galaxy_output.h"
    float DiskMass_elements[NUM_ELEMENTS]; // Msun // Mass of elements locked up in stars in disk.
    float BulgeMass_elements[NUM_ELEMENTS]; // Msun // Mass of elements locked up in stars in bulge.

    float DiskMassRings_elements[RNUM][NUM_ELEMENTS]; // Msun // Mass of elements locked up in stars in each annular ring.
    float BulgeMassRings_elements[RNUM][NUM_ELEMENTS]; // Msun // Mass of elements locked up in stars in each annular ring.

    float ColdGas_elements[NUM_ELEMENTS]; // Msun // Mass of elements locked up in stars in cold gas.

    float ColdGasRings_elements[RNUM][NUM_ELEMENTS]; // Msun // Mass of elements locked up in cold gas in each annular ring.
# 350 "./code/h_galaxy_output.h"
    float HotGas_elements[NUM_ELEMENTS]; // Msun // Mass of elements locked up in hot gas.
    /* float ReheatedGas_elements[NUM_ELEMENTS]; // Msun // Mass of elements locked up in reheated gas. */
    float ICM_elements[NUM_ELEMENTS]; // Msun // Mass of elements locked up in stars in the ICM
    float EjectedMass_elements[NUM_ELEMENTS]; // Msun // Mass of elements locked up in ejected gas.







  float t_des[RNUM];

  float t_sput_HotGas;


  float t_sput_EjectedMass;

  //struct DustRates DustColdGasRates; // ? // Rates of creation and destruction of dust
  float DustColdGasRates[NUM_COLDGAS_DUST_RATES]; // ? // Rates of creation and destruction of dust in the ColdGas

  float DustHotGasRates[NUM_HOTGAS_DUST_RATES]; // ? // Rates of creation and destruction of dust in the Hotgas


  float DustEjectedMassRates; // Rate of destruction of dust in the EjectedMass



  float t_acc[RNUM];





  float f_cmax[RNUM][NUM_ELEMENTS]; // unitless // max fraction of molecular media //ROB: I think this is actually the maximum mass fraction of element ee in clouds that is in dust.




  float ColdGasDiff_elements[NUM_ELEMENTS]; // Msun // Mass of elements in the diffused phase (in ColdGas)
  float ColdGasClouds_elements[NUM_ELEMENTS]; // Msun // Mass of elements in the cloud phase (in ColdGas)
  float DustColdGasDiff_elements[NUM_ELEMENTS]; // Msun // Mass of elements in the diffuse phase locked up in dust (in ColdGas)
  float DustColdGasClouds_elements[NUM_ELEMENTS]; // Msun // Mass of elements in the cloud phase locked up in dust (in ColdGas)

  float DustHotGas_elements[NUM_ELEMENTS]; // Msun // Mass of elements locked up in dust (in HotGas)


  float DustEjectedMass_elements[NUM_ELEMENTS]; // Msun // Mass of elements locked up in dust (in EjectedMass)


    float ColdGasDiffRings_elements[RNUM][NUM_ELEMENTS]; // Msun // Mass of elements in diffuse gas in each annular ring
    float ColdGasCloudsRings_elements[RNUM][NUM_ELEMENTS]; // Msun // Mass of elements in clouds in each annular ring
    float DustColdGasDiffRings_elements[RNUM][NUM_ELEMENTS]; // Msun // Mass of elements in diffuse gas locked up in dust in each annular ring
    float DustColdGasCloudsRings_elements[RNUM][NUM_ELEMENTS]; // Msun // Mass of elements in diffuse gas locked up in dust in each annular ring


};

// next only of interest to DB output, which generally requires complete tree

struct SFH_BIN {
    long long GalID; // None // ID of the galaxy
    short snapnum; // None // snapnum of the galaxy, repeated here for faster lookups of times etc
    short sfh_ibin; // None //Index of highest bin currently in use
    /* float sfh_time; // yr // Lookback time at the middle of bin. */
    /* float sfh_dt; // yr // time width of bin. */
    float sfh_DiskMass; // 1e10 Msun/h // SFH of disk
    float sfh_BulgeMass; // 1e10 Msun/h // SFH of bulge

    float sfh_DiskMassRings[RNUM]; // 10^10 Msun/h // Star formation history in the disk RINGS.
    float sfh_BulgeMassRings[RNUM]; // 10^10 Msun/h // Star formation history in the bulge RINGS.

    float sfh_ICM; // 1e10 Msun/h // SFH of ICM
# 431 "./code/h_galaxy_output.h"
    float sfh_MetalsDiskMass[NUM_METAL_CHANNELS]; // 1e10 Msun/h // Metals locked up in stars in disk.
    float sfh_MetalsBulgeMass[NUM_METAL_CHANNELS]; // 1e10 Msun/h // Metals locked up in stars in bulge.
    float sfh_MetalsICM[NUM_METAL_CHANNELS]; // 1e10 Msun/h // Metals locked up in stars in ICM.
};

struct SFH_Time {
    int snapnum; // None // snapnum
    int bin; // None // index of current bin
    /* proposal: in output write the start of the bin and its end, rather than centre and dt */
    double lookbacktime; // yr // lookback time to centre of current bin
    double dt; // yr // width of the current bin
    int nbins; // None // number of highest resolution bins used to create current bin
};

#pragma pack()
