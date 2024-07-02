CREATE TABLE GALAXIES (
 type INTEGER NOT NULL 
,  haloIndex INTEGER NOT NULL 
,  haloM_Crit200 REAL NOT NULL 
,  snapNum INTEGER NOT NULL 
,  lookBackTimeToSnap REAL NOT NULL 
,  centralMvir REAL NOT NULL 
,  centralRvir REAL NOT NULL 
,  distanceToCentralGal_1 REAL NOT NULL 
,  distanceToCentralGal_2 REAL NOT NULL 
,  distanceToCentralGal_3 REAL NOT NULL 
,  pos_1 REAL NOT NULL 
,  pos_2 REAL NOT NULL 
,  pos_3 REAL NOT NULL 
,  vel_1 REAL NOT NULL 
,  vel_2 REAL NOT NULL 
,  vel_3 REAL NOT NULL 
,  len INTEGER NOT NULL 
,  mvir REAL NOT NULL 
,  rvir REAL NOT NULL 
,  vvir REAL NOT NULL 
,  vmax REAL NOT NULL 
,  coldGasSpin_1 REAL NOT NULL 
,  coldGasSpin_2 REAL NOT NULL 
,  coldGasSpin_3 REAL NOT NULL 
,  diskSpin_1 REAL NOT NULL 
,  diskSpin_2 REAL NOT NULL 
,  diskSpin_3 REAL NOT NULL 
,  coldGas REAL NOT NULL 
,  h2fraction REAL NOT NULL 
,  stellarMass REAL NOT NULL 
,  diskMass REAL NOT NULL 
,  bulgeMass REAL NOT NULL 
,  hotGas REAL NOT NULL 
,  ejectedMass REAL NOT NULL 
,  blackHoleMass REAL NOT NULL 
,  iCM REAL NOT NULL 
,  diskSNIIRate REAL NOT NULL 
,  bulgeSNIIRate REAL NOT NULL 
,  iCMSNIIRate REAL NOT NULL 
,  diskSNIaRate REAL NOT NULL 
,  bulgeSNIaRate REAL NOT NULL 
,  iCMSNIaRate REAL NOT NULL 
,  sfr REAL NOT NULL 
,  sfrBulge REAL NOT NULL 
,  bulgeSize REAL NOT NULL 
,  diskRadius REAL NOT NULL 
,  coldGasRadius REAL NOT NULL 
,  stellarHalfMassRadius REAL NOT NULL 
,  massWeightAge REAL NOT NULL 
,  t_sput_HotGas REAL NOT NULL 
,  t_sput_EjectedMass REAL NOT NULL 
,  dustEjectedMassRates REAL NOT NULL 
)
CREATE TABLE SFH_Times (
 snapnum INTEGER NOT NULL 
,  none // NOT NULL 
,  none // NOT NULL 
,  yr // NOT NULL 
,  yr // NOT NULL 
 -- size = 20
)
