# python structure for modified L-Galaxies 2020 model (Yates+20b)
import numpy
LGalaxiesStruct = numpy.dtype([
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
('ColdGasRings',numpy.float32,12),
('H2fraction',numpy.float32,1),
('H2fractionRings',numpy.float32,12),
('StellarMass',numpy.float32,1),
('DiskMass',numpy.float32,1),
('BulgeMass',numpy.float32,1),
('DiskMassRings',numpy.float32,12),
('BulgeMassRings',numpy.float32,12),
('HotGas',numpy.float32,1),
('EjectedMass',numpy.float32,1),
('BlackHoleMass',numpy.float32,1),
('ICM',numpy.float32,1),
('MassFromInSitu',numpy.float32,1),
('MassFromMergers',numpy.float32,1),
('MassFromBursts',numpy.float32,1),
('MetalsColdGas',numpy.float32,3),
('MetalsColdGasRings',numpy.float32,[12,3]),
('MetalsStellarMass',numpy.float32,3),
('MetalsDiskMass',numpy.float32,3),
('MetalsBulgeMass',numpy.float32,3),
('MetalsDiskMassRings',numpy.float32,[12,3]),
('MetalsBulgeMassRings',numpy.float32,[12,3]),
('MetalsHotGas',numpy.float32,3),
('MetalsEjectedMass',numpy.float32,3),
('MetalsICM',numpy.float32,3),
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
('SfrRings',numpy.float32,12),
('SfrInstRings',numpy.float32,12),
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
('sfh_DiskMass',numpy.float32,20),
('sfh_BulgeMass',numpy.float32,20),
('sfh_DiskMassRings',numpy.float32,[12,20]),
('sfh_BulgeMassRings',numpy.float32,[12,20]),
('sfh_ICM',numpy.float32,20),
('sfh_MetalsDiskMass',numpy.float32,[20,3]),
('sfh_MetalsBulgeMass',numpy.float32,[20,3]),
('sfh_MetalsDiskMassRings',numpy.float32,12,20,3),
('sfh_MetalsBulgeMassRings',numpy.float32,12,20,3),
('sfh_MetalsICM',numpy.float32,[20,3]),
('sfh_DiskMass_elements',numpy.float32,[20,11]),
('sfh_BulgeMass_elements',numpy.float32,[20,11]),
('sfh_DiskMass_elementsRings',numpy.float32,12,20,11),
('sfh_BulgeMass_elementsRings',numpy.float32,12,20,11),
('sfh_ICM_elements',numpy.float32,[20,11]),
('DiskMass_elements',numpy.float32,11),
('BulgeMass_elements',numpy.float32,11),
('DiskMassRings_elements',numpy.float32,[12,11]),
('BulgeMassRings_elements',numpy.float32,[12,11]),
('ColdGas_elements',numpy.float32,11),
('ColdGasRings_elements',numpy.float32,[12,11]),
('HotGas_elements',numpy.float32,11),
('ICM_elements',numpy.float32,11),
('EjectedMass_elements',numpy.float32,11),
('ending','i4',0)
])
properties_used = {}
for el in LGalaxiesStruct.names:
	properties_used[el] = True
