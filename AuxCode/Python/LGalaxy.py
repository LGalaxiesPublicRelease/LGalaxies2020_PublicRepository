# numpy dtype for LGAL_GAL_STRUCT
import numpy
struct_dtype = numpy.dtype([
('Type',numpy.int32,1),
('HaloIndex',numpy.int32,1),
('HaloM_Mean200',numpy.float32,1),
('HaloM_Crit200',numpy.float32,1),
('HaloM_TopHat',numpy.float32,1),
('HaloPos',numpy.float32,3),
('HaloVel',numpy.float32,3),
('HaloVelDisp',numpy.float32,1),
('HaloVmax',numpy.float32,1),
('HaloSpin',numpy.float32,3),
('SnapNum',numpy.int32,1),
('LookBackTimeToSnap',numpy.float32,1),
('CentralMvir',numpy.float32,1),
('CentralRvir',numpy.float32,1),
('DistanceToCentralGal',numpy.float32,3),
('Pos',numpy.float32,3),
('Vel',numpy.float32,3),
('Len',numpy.int32,1),
('Mvir',numpy.float32,1),
('Rvir',numpy.float32,1),
('Vvir',numpy.float32,1),
('Vmax',numpy.float32,1),
('ColdGasSpin',numpy.float32,3),
('DiskSpin',numpy.float32,3),
('InfallVmax',numpy.float32,1),
('InfallVmaxPeak',numpy.float32,1),
('InfallSnap',numpy.int32,1),
('InfallHotGas',numpy.float32,1),
('HotRadius',numpy.float32,1),
('OriMergTime',numpy.float32,1),
('MergTime',numpy.float32,1),
('ColdGas',numpy.float32,1),
('ColdGasRings',numpy.float32,RNUM),
('H2fraction',numpy.float32,1),
('H2fractionRings',numpy.float32,RNUM),
('StellarMass',numpy.float32,1),
('DiskMass',numpy.float32,1),
('BulgeMass',numpy.float32,1),
('DiskMassRings',numpy.float32,RNUM),
('BulgeMassRings',numpy.float32,RNUM),
('HotGas',numpy.float32,1),
('EjectedMass',numpy.float32,1),
('BlackHoleMass',numpy.float32,1),
('ICM',numpy.float32,1),
('MassFromInSitu',numpy.float32,1),
('MassFromMergers',numpy.float32,1),
('MassFromBursts',numpy.float32,1),
('MetalsColdGas',numpy.float32,NUM_METAL_CHANNELS),
('MetalsColdGasRings',numpy.float32,RNUM),
('MetalsStellarMass',numpy.float32,NUM_METAL_CHANNELS),
('MetalsDiskMass',numpy.float32,NUM_METAL_CHANNELS),
('MetalsBulgeMass',numpy.float32,NUM_METAL_CHANNELS),
('MetalsDiskMassRings',numpy.float32,RNUM),
('MetalsBulgeMassRings',numpy.float32,RNUM),
('MetalsHotGas',numpy.float32,NUM_METAL_CHANNELS),
('MetalsEjectedMass',numpy.float32,NUM_METAL_CHANNELS),
('MetalsICM',numpy.float32,NUM_METAL_CHANNELS),
('DiskSNIIRate',numpy.float32,1),
('BulgeSNIIRate',numpy.float32,1),
('ICMSNIIRate',numpy.float32,1),
('MassReturnRateToColdGas',numpy.float32,1),
('MassReturnRateToHotGas',numpy.float32,1),
('MetalsReturnRateToHotGas',numpy.float32,1),
('PrimordialAccretionRate',numpy.float32,1),
('CoolingRadius',numpy.float32,1),
('ReheatingRate',numpy.float32,1),
('MetalsReheatingRate',numpy.float32,1),
('EjectionRate',numpy.float32,1),
('CoolingRate',numpy.float32,1),
('CoolingRate_beforeAGN',numpy.float32,1),
('QuasarAccretionRate',numpy.float32,1),
('RadioAccretionRate',numpy.float32,1),
('Sfr',numpy.float32,1),
('SfrInst',numpy.float32,1),
('SfrRings',numpy.float32,RNUM),
('SfrInstRings',numpy.float32,RNUM),
('SfrBulge',numpy.float32,1),
('XrayLum',numpy.float32,1),
('BulgeSize',numpy.float32,1),
('DiskRadius',numpy.float32,1),
('ColdGasRadius',numpy.float32,1),
('StellarHalfMassRadius',numpy.float32,1),
('StellarHalfLightRadius',numpy.float32,1),
('CosInclination',numpy.float32,1),
('DisruptOn',numpy.int32,1),
('MergeOn',numpy.int32,1),
('MagDust',numpy.float32,5),
('Mag',numpy.float32,5),
('MagBulge',numpy.float32,5),
('MassWeightAge',numpy.float32,1),
('rbandWeightAge',numpy.float32,1),
('sfh_ibin',numpy.int32,1),
('sfh_numbins',numpy.int32,1),
('sfh_DiskMass',numpy.float32,SFH_NBIN),
('sfh_BulgeMass',numpy.float32,SFH_NBIN),
('sfh_DiskMassRings',numpy.float32,RNUM),
('sfh_BulgeMassRings',numpy.float32,RNUM),
('sfh_ICM',numpy.float32,SFH_NBIN),
('sfh_MetalsDiskMass',numpy.float32,SFH_NBIN),
('sfh_MetalsBulgeMass',numpy.float32,SFH_NBIN),
('sfh_MetalsICM',numpy.float32,SFH_NBIN),
('sfh_DiskMass_elements',numpy.float32,SFH_NBIN),
('sfh_BulgeMass_elements',numpy.float32,SFH_NBIN),
('sfh_ICM_elements',numpy.float32,SFH_NBIN),
('DiskMass_elements',numpy.float32,NUM_ELEMENTS),
('BulgeMass_elements',numpy.float32,NUM_ELEMENTS),
('DiskMassRings_elements',numpy.float32,RNUM),
('BulgeMassRings_elements',numpy.float32,RNUM),
('ColdGas_elements',numpy.float32,NUM_ELEMENTS),
('ColdGasRings_elements',numpy.float32,RNUM),
('HotGas_elements',numpy.float32,NUM_ELEMENTS),
('ICM_elements',numpy.float32,NUM_ELEMENTS),
('EjectedMass_elements',numpy.float32,NUM_ELEMENTS),
('ending','i4',0)
])
properties_used = {}
for el in struct_dtype.names:
	properties_used[el] = False
