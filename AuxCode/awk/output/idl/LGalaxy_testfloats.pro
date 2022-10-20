;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
FUNCTION LGalaxy_testfloats, LGs, nstart 
; test whether floats are NaN or too small for SQLServer
; assumes the existence of a function testFloat accepting an array of floats
 badranges = []
 bad = 0
 sel = testFloat(LGs.HaloM_Mean200)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'HaloM_Mean200 --- ', nstart+sel
     print, 'HaloM_Mean200 --- ', LGs[sel].HaloM_Mean200
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.HaloM_Crit200)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'HaloM_Crit200 --- ', nstart+sel
     print, 'HaloM_Crit200 --- ', LGs[sel].HaloM_Crit200
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.HaloM_TopHat)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'HaloM_TopHat --- ', nstart+sel
     print, 'HaloM_TopHat --- ', LGs[sel].HaloM_TopHat
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.HaloPos(0))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'HaloPos[0] --- ', nstart+sel
     print, 'HaloPos --- ', LGs[sel].HaloPos(0)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.HaloPos(1))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'HaloPos[1] --- ', nstart+sel
     print, 'HaloPos --- ', LGs[sel].HaloPos(1)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.HaloPos(2))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'HaloPos[2] --- ', nstart+sel
     print, 'HaloPos --- ', LGs[sel].HaloPos(2)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.HaloVel(0))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'HaloVel[0] --- ', nstart+sel
     print, 'HaloVel --- ', LGs[sel].HaloVel(0)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.HaloVel(1))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'HaloVel[1] --- ', nstart+sel
     print, 'HaloVel --- ', LGs[sel].HaloVel(1)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.HaloVel(2))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'HaloVel[2] --- ', nstart+sel
     print, 'HaloVel --- ', LGs[sel].HaloVel(2)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.HaloVelDisp)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'HaloVelDisp --- ', nstart+sel
     print, 'HaloVelDisp --- ', LGs[sel].HaloVelDisp
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.HaloVmax)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'HaloVmax --- ', nstart+sel
     print, 'HaloVmax --- ', LGs[sel].HaloVmax
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.HaloSpin(0))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'HaloSpin[0] --- ', nstart+sel
     print, 'HaloSpin --- ', LGs[sel].HaloSpin(0)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.HaloSpin(1))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'HaloSpin[1] --- ', nstart+sel
     print, 'HaloSpin --- ', LGs[sel].HaloSpin(1)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.HaloSpin(2))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'HaloSpin[2] --- ', nstart+sel
     print, 'HaloSpin --- ', LGs[sel].HaloSpin(2)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.LookBackTimeToSnap)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'LookBackTimeToSnap --- ', nstart+sel
     print, 'LookBackTimeToSnap --- ', LGs[sel].LookBackTimeToSnap
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.CentralMvir)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'CentralMvir --- ', nstart+sel
     print, 'CentralMvir --- ', LGs[sel].CentralMvir
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.CentralRvir)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'CentralRvir --- ', nstart+sel
     print, 'CentralRvir --- ', LGs[sel].CentralRvir
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.DistanceToCentralGal(0))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'DistanceToCentralGal[0] --- ', nstart+sel
     print, 'DistanceToCentralGal --- ', LGs[sel].DistanceToCentralGal(0)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.DistanceToCentralGal(1))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'DistanceToCentralGal[1] --- ', nstart+sel
     print, 'DistanceToCentralGal --- ', LGs[sel].DistanceToCentralGal(1)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.DistanceToCentralGal(2))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'DistanceToCentralGal[2] --- ', nstart+sel
     print, 'DistanceToCentralGal --- ', LGs[sel].DistanceToCentralGal(2)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Pos(0))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Pos[0] --- ', nstart+sel
     print, 'Pos --- ', LGs[sel].Pos(0)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Pos(1))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Pos[1] --- ', nstart+sel
     print, 'Pos --- ', LGs[sel].Pos(1)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Pos(2))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Pos[2] --- ', nstart+sel
     print, 'Pos --- ', LGs[sel].Pos(2)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Vel(0))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Vel[0] --- ', nstart+sel
     print, 'Vel --- ', LGs[sel].Vel(0)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Vel(1))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Vel[1] --- ', nstart+sel
     print, 'Vel --- ', LGs[sel].Vel(1)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Vel(2))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Vel[2] --- ', nstart+sel
     print, 'Vel --- ', LGs[sel].Vel(2)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Mvir)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Mvir --- ', nstart+sel
     print, 'Mvir --- ', LGs[sel].Mvir
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Rvir)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Rvir --- ', nstart+sel
     print, 'Rvir --- ', LGs[sel].Rvir
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Vvir)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Vvir --- ', nstart+sel
     print, 'Vvir --- ', LGs[sel].Vvir
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Vmax)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Vmax --- ', nstart+sel
     print, 'Vmax --- ', LGs[sel].Vmax
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ColdGasSpin(0))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ColdGasSpin[0] --- ', nstart+sel
     print, 'ColdGasSpin --- ', LGs[sel].ColdGasSpin(0)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ColdGasSpin(1))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ColdGasSpin[1] --- ', nstart+sel
     print, 'ColdGasSpin --- ', LGs[sel].ColdGasSpin(1)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ColdGasSpin(2))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ColdGasSpin[2] --- ', nstart+sel
     print, 'ColdGasSpin --- ', LGs[sel].ColdGasSpin(2)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.DiskSpin(0))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'DiskSpin[0] --- ', nstart+sel
     print, 'DiskSpin --- ', LGs[sel].DiskSpin(0)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.DiskSpin(1))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'DiskSpin[1] --- ', nstart+sel
     print, 'DiskSpin --- ', LGs[sel].DiskSpin(1)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.DiskSpin(2))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'DiskSpin[2] --- ', nstart+sel
     print, 'DiskSpin --- ', LGs[sel].DiskSpin(2)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.InfallVmax)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'InfallVmax --- ', nstart+sel
     print, 'InfallVmax --- ', LGs[sel].InfallVmax
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.InfallVmaxPeak)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'InfallVmaxPeak --- ', nstart+sel
     print, 'InfallVmaxPeak --- ', LGs[sel].InfallVmaxPeak
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.InfallHotGas)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'InfallHotGas --- ', nstart+sel
     print, 'InfallHotGas --- ', LGs[sel].InfallHotGas
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.HotRadius)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'HotRadius --- ', nstart+sel
     print, 'HotRadius --- ', LGs[sel].HotRadius
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.OriMergTime)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'OriMergTime --- ', nstart+sel
     print, 'OriMergTime --- ', LGs[sel].OriMergTime
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MergTime)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MergTime --- ', nstart+sel
     print, 'MergTime --- ', LGs[sel].MergTime
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ColdGas)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ColdGas --- ', nstart+sel
     print, 'ColdGas --- ', LGs[sel].ColdGas
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.H2fraction)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'H2fraction --- ', nstart+sel
     print, 'H2fraction --- ', LGs[sel].H2fraction
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.StellarMass)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'StellarMass --- ', nstart+sel
     print, 'StellarMass --- ', LGs[sel].StellarMass
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.DiskMass)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'DiskMass --- ', nstart+sel
     print, 'DiskMass --- ', LGs[sel].DiskMass
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.BulgeMass)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'BulgeMass --- ', nstart+sel
     print, 'BulgeMass --- ', LGs[sel].BulgeMass
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.HotGas)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'HotGas --- ', nstart+sel
     print, 'HotGas --- ', LGs[sel].HotGas
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.EjectedMass)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'EjectedMass --- ', nstart+sel
     print, 'EjectedMass --- ', LGs[sel].EjectedMass
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.BlackHoleMass)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'BlackHoleMass --- ', nstart+sel
     print, 'BlackHoleMass --- ', LGs[sel].BlackHoleMass
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ICM)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ICM --- ', nstart+sel
     print, 'ICM --- ', LGs[sel].ICM
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MassFromInSitu)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MassFromInSitu --- ', nstart+sel
     print, 'MassFromInSitu --- ', LGs[sel].MassFromInSitu
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MassFromMergers)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MassFromMergers --- ', nstart+sel
     print, 'MassFromMergers --- ', LGs[sel].MassFromMergers
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MassFromBursts)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MassFromBursts --- ', nstart+sel
     print, 'MassFromBursts --- ', LGs[sel].MassFromBursts
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.DiskSNIIRate)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'DiskSNIIRate --- ', nstart+sel
     print, 'DiskSNIIRate --- ', LGs[sel].DiskSNIIRate
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.BulgeSNIIRate)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'BulgeSNIIRate --- ', nstart+sel
     print, 'BulgeSNIIRate --- ', LGs[sel].BulgeSNIIRate
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ICMSNIIRate)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ICMSNIIRate --- ', nstart+sel
     print, 'ICMSNIIRate --- ', LGs[sel].ICMSNIIRate
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MassReturnRateToColdGas)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MassReturnRateToColdGas --- ', nstart+sel
     print, 'MassReturnRateToColdGas --- ', LGs[sel].MassReturnRateToColdGas
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MassReturnRateToHotGas)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MassReturnRateToHotGas --- ', nstart+sel
     print, 'MassReturnRateToHotGas --- ', LGs[sel].MassReturnRateToHotGas
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MetalsReturnRateToHotGas)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MetalsReturnRateToHotGas --- ', nstart+sel
     print, 'MetalsReturnRateToHotGas --- ', LGs[sel].MetalsReturnRateToHotGas
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.PrimordialAccretionRate)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'PrimordialAccretionRate --- ', nstart+sel
     print, 'PrimordialAccretionRate --- ', LGs[sel].PrimordialAccretionRate
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.CoolingRadius)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'CoolingRadius --- ', nstart+sel
     print, 'CoolingRadius --- ', LGs[sel].CoolingRadius
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ReheatingRate)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ReheatingRate --- ', nstart+sel
     print, 'ReheatingRate --- ', LGs[sel].ReheatingRate
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MetalsReheatingRate)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MetalsReheatingRate --- ', nstart+sel
     print, 'MetalsReheatingRate --- ', LGs[sel].MetalsReheatingRate
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.EjectionRate)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'EjectionRate --- ', nstart+sel
     print, 'EjectionRate --- ', LGs[sel].EjectionRate
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.CoolingRate)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'CoolingRate --- ', nstart+sel
     print, 'CoolingRate --- ', LGs[sel].CoolingRate
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.CoolingRate_beforeAGN)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'CoolingRate_beforeAGN --- ', nstart+sel
     print, 'CoolingRate_beforeAGN --- ', LGs[sel].CoolingRate_beforeAGN
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.QuasarAccretionRate)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'QuasarAccretionRate --- ', nstart+sel
     print, 'QuasarAccretionRate --- ', LGs[sel].QuasarAccretionRate
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.RadioAccretionRate)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'RadioAccretionRate --- ', nstart+sel
     print, 'RadioAccretionRate --- ', LGs[sel].RadioAccretionRate
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Sfr)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Sfr --- ', nstart+sel
     print, 'Sfr --- ', LGs[sel].Sfr
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.SfrBulge)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'SfrBulge --- ', nstart+sel
     print, 'SfrBulge --- ', LGs[sel].SfrBulge
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.XrayLum)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'XrayLum --- ', nstart+sel
     print, 'XrayLum --- ', LGs[sel].XrayLum
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.BulgeSize)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'BulgeSize --- ', nstart+sel
     print, 'BulgeSize --- ', LGs[sel].BulgeSize
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.DiskRadius)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'DiskRadius --- ', nstart+sel
     print, 'DiskRadius --- ', LGs[sel].DiskRadius
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ColdGasRadius)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ColdGasRadius --- ', nstart+sel
     print, 'ColdGasRadius --- ', LGs[sel].ColdGasRadius
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.StellarHalfMassRadius)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'StellarHalfMassRadius --- ', nstart+sel
     print, 'StellarHalfMassRadius --- ', LGs[sel].StellarHalfMassRadius
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.StellarHalfLightRadius)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'StellarHalfLightRadius --- ', nstart+sel
     print, 'StellarHalfLightRadius --- ', LGs[sel].StellarHalfLightRadius
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.CosInclination)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'CosInclination --- ', nstart+sel
     print, 'CosInclination --- ', LGs[sel].CosInclination
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagDust(0))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagDust[0] --- ', nstart+sel
     print, 'MagDust --- ', LGs[sel].MagDust(0)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagDust(1))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagDust[1] --- ', nstart+sel
     print, 'MagDust --- ', LGs[sel].MagDust(1)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagDust(2))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagDust[2] --- ', nstart+sel
     print, 'MagDust --- ', LGs[sel].MagDust(2)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagDust(3))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagDust[3] --- ', nstart+sel
     print, 'MagDust --- ', LGs[sel].MagDust(3)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagDust(4))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagDust[4] --- ', nstart+sel
     print, 'MagDust --- ', LGs[sel].MagDust(4)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Mag(0))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Mag[0] --- ', nstart+sel
     print, 'Mag --- ', LGs[sel].Mag(0)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Mag(1))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Mag[1] --- ', nstart+sel
     print, 'Mag --- ', LGs[sel].Mag(1)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Mag(2))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Mag[2] --- ', nstart+sel
     print, 'Mag --- ', LGs[sel].Mag(2)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Mag(3))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Mag[3] --- ', nstart+sel
     print, 'Mag --- ', LGs[sel].Mag(3)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Mag(4))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Mag[4] --- ', nstart+sel
     print, 'Mag --- ', LGs[sel].Mag(4)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagBulge(0))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagBulge[0] --- ', nstart+sel
     print, 'MagBulge --- ', LGs[sel].MagBulge(0)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagBulge(1))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagBulge[1] --- ', nstart+sel
     print, 'MagBulge --- ', LGs[sel].MagBulge(1)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagBulge(2))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagBulge[2] --- ', nstart+sel
     print, 'MagBulge --- ', LGs[sel].MagBulge(2)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagBulge(3))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagBulge[3] --- ', nstart+sel
     print, 'MagBulge --- ', LGs[sel].MagBulge(3)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagBulge(4))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagBulge[4] --- ', nstart+sel
     print, 'MagBulge --- ', LGs[sel].MagBulge(4)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MassWeightAge)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MassWeightAge --- ', nstart+sel
     print, 'MassWeightAge --- ', LGs[sel].MassWeightAge
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.rbandWeightAge)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'rbandWeightAge --- ', nstart+sel
     print, 'rbandWeightAge --- ', LGs[sel].rbandWeightAge
     badranges=[badranges,sel]
 endif
if(bad) then begin 
     print, 'badranges found: ',badranges
endif
return, badranges
end
