;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
PRO LGalaxy_zerofloats, LGs 
; test whether floats are NaN or too small for SQLServer
; if so, set offending values to 0
; assumes the existence of a function testFloat accepting an array of floats
 sel = testFloat(LGs.HaloM_Mean200)
 if(sel(0) gt -1) then begin
     LGs[sel].HaloM_Mean200 = 0
 endif
 sel = testFloat(LGs.HaloM_Crit200)
 if(sel(0) gt -1) then begin
     LGs[sel].HaloM_Crit200 = 0
 endif
 sel = testFloat(LGs.HaloM_TopHat)
 if(sel(0) gt -1) then begin
     LGs[sel].HaloM_TopHat = 0
 endif
 sel = testFloat(LGs.HaloPos(0))
 if(sel(0) gt -1) then begin
     LGs[sel].HaloPos(0) = 0
 endif
 sel = testFloat(LGs.HaloPos(1))
 if(sel(0) gt -1) then begin
     LGs[sel].HaloPos(1) = 0
 endif
 sel = testFloat(LGs.HaloPos(2))
 if(sel(0) gt -1) then begin
     LGs[sel].HaloPos(2) = 0
 endif
 sel = testFloat(LGs.HaloVel(0))
 if(sel(0) gt -1) then begin
     LGs[sel].HaloVel(0) = 0
 endif
 sel = testFloat(LGs.HaloVel(1))
 if(sel(0) gt -1) then begin
     LGs[sel].HaloVel(1) = 0
 endif
 sel = testFloat(LGs.HaloVel(2))
 if(sel(0) gt -1) then begin
     LGs[sel].HaloVel(2) = 0
 endif
 sel = testFloat(LGs.HaloVelDisp)
 if(sel(0) gt -1) then begin
     LGs[sel].HaloVelDisp = 0
 endif
 sel = testFloat(LGs.HaloVmax)
 if(sel(0) gt -1) then begin
     LGs[sel].HaloVmax = 0
 endif
 sel = testFloat(LGs.HaloSpin(0))
 if(sel(0) gt -1) then begin
     LGs[sel].HaloSpin(0) = 0
 endif
 sel = testFloat(LGs.HaloSpin(1))
 if(sel(0) gt -1) then begin
     LGs[sel].HaloSpin(1) = 0
 endif
 sel = testFloat(LGs.HaloSpin(2))
 if(sel(0) gt -1) then begin
     LGs[sel].HaloSpin(2) = 0
 endif
 sel = testFloat(LGs.LookBackTimeToSnap)
 if(sel(0) gt -1) then begin
     LGs[sel].LookBackTimeToSnap = 0
 endif
 sel = testFloat(LGs.CentralMvir)
 if(sel(0) gt -1) then begin
     LGs[sel].CentralMvir = 0
 endif
 sel = testFloat(LGs.CentralRvir)
 if(sel(0) gt -1) then begin
     LGs[sel].CentralRvir = 0
 endif
 sel = testFloat(LGs.DistanceToCentralGal(0))
 if(sel(0) gt -1) then begin
     LGs[sel].DistanceToCentralGal(0) = 0
 endif
 sel = testFloat(LGs.DistanceToCentralGal(1))
 if(sel(0) gt -1) then begin
     LGs[sel].DistanceToCentralGal(1) = 0
 endif
 sel = testFloat(LGs.DistanceToCentralGal(2))
 if(sel(0) gt -1) then begin
     LGs[sel].DistanceToCentralGal(2) = 0
 endif
 sel = testFloat(LGs.Pos(0))
 if(sel(0) gt -1) then begin
     LGs[sel].Pos(0) = 0
 endif
 sel = testFloat(LGs.Pos(1))
 if(sel(0) gt -1) then begin
     LGs[sel].Pos(1) = 0
 endif
 sel = testFloat(LGs.Pos(2))
 if(sel(0) gt -1) then begin
     LGs[sel].Pos(2) = 0
 endif
 sel = testFloat(LGs.Vel(0))
 if(sel(0) gt -1) then begin
     LGs[sel].Vel(0) = 0
 endif
 sel = testFloat(LGs.Vel(1))
 if(sel(0) gt -1) then begin
     LGs[sel].Vel(1) = 0
 endif
 sel = testFloat(LGs.Vel(2))
 if(sel(0) gt -1) then begin
     LGs[sel].Vel(2) = 0
 endif
 sel = testFloat(LGs.Mvir)
 if(sel(0) gt -1) then begin
     LGs[sel].Mvir = 0
 endif
 sel = testFloat(LGs.Rvir)
 if(sel(0) gt -1) then begin
     LGs[sel].Rvir = 0
 endif
 sel = testFloat(LGs.Vvir)
 if(sel(0) gt -1) then begin
     LGs[sel].Vvir = 0
 endif
 sel = testFloat(LGs.Vmax)
 if(sel(0) gt -1) then begin
     LGs[sel].Vmax = 0
 endif
 sel = testFloat(LGs.ColdGasSpin(0))
 if(sel(0) gt -1) then begin
     LGs[sel].ColdGasSpin(0) = 0
 endif
 sel = testFloat(LGs.ColdGasSpin(1))
 if(sel(0) gt -1) then begin
     LGs[sel].ColdGasSpin(1) = 0
 endif
 sel = testFloat(LGs.ColdGasSpin(2))
 if(sel(0) gt -1) then begin
     LGs[sel].ColdGasSpin(2) = 0
 endif
 sel = testFloat(LGs.DiskSpin(0))
 if(sel(0) gt -1) then begin
     LGs[sel].DiskSpin(0) = 0
 endif
 sel = testFloat(LGs.DiskSpin(1))
 if(sel(0) gt -1) then begin
     LGs[sel].DiskSpin(1) = 0
 endif
 sel = testFloat(LGs.DiskSpin(2))
 if(sel(0) gt -1) then begin
     LGs[sel].DiskSpin(2) = 0
 endif
 sel = testFloat(LGs.InfallVmax)
 if(sel(0) gt -1) then begin
     LGs[sel].InfallVmax = 0
 endif
 sel = testFloat(LGs.InfallVmaxPeak)
 if(sel(0) gt -1) then begin
     LGs[sel].InfallVmaxPeak = 0
 endif
 sel = testFloat(LGs.InfallHotGas)
 if(sel(0) gt -1) then begin
     LGs[sel].InfallHotGas = 0
 endif
 sel = testFloat(LGs.HotRadius)
 if(sel(0) gt -1) then begin
     LGs[sel].HotRadius = 0
 endif
 sel = testFloat(LGs.OriMergTime)
 if(sel(0) gt -1) then begin
     LGs[sel].OriMergTime = 0
 endif
 sel = testFloat(LGs.MergTime)
 if(sel(0) gt -1) then begin
     LGs[sel].MergTime = 0
 endif
 sel = testFloat(LGs.ColdGas)
 if(sel(0) gt -1) then begin
     LGs[sel].ColdGas = 0
 endif
 sel = testFloat(LGs.H2fraction)
 if(sel(0) gt -1) then begin
     LGs[sel].H2fraction = 0
 endif
 sel = testFloat(LGs.StellarMass)
 if(sel(0) gt -1) then begin
     LGs[sel].StellarMass = 0
 endif
 sel = testFloat(LGs.DiskMass)
 if(sel(0) gt -1) then begin
     LGs[sel].DiskMass = 0
 endif
 sel = testFloat(LGs.BulgeMass)
 if(sel(0) gt -1) then begin
     LGs[sel].BulgeMass = 0
 endif
 sel = testFloat(LGs.HotGas)
 if(sel(0) gt -1) then begin
     LGs[sel].HotGas = 0
 endif
 sel = testFloat(LGs.EjectedMass)
 if(sel(0) gt -1) then begin
     LGs[sel].EjectedMass = 0
 endif
 sel = testFloat(LGs.BlackHoleMass)
 if(sel(0) gt -1) then begin
     LGs[sel].BlackHoleMass = 0
 endif
 sel = testFloat(LGs.ICM)
 if(sel(0) gt -1) then begin
     LGs[sel].ICM = 0
 endif
 sel = testFloat(LGs.MassFromInSitu)
 if(sel(0) gt -1) then begin
     LGs[sel].MassFromInSitu = 0
 endif
 sel = testFloat(LGs.MassFromMergers)
 if(sel(0) gt -1) then begin
     LGs[sel].MassFromMergers = 0
 endif
 sel = testFloat(LGs.MassFromBursts)
 if(sel(0) gt -1) then begin
     LGs[sel].MassFromBursts = 0
 endif
 sel = testFloat(LGs.DiskSNIIRate)
 if(sel(0) gt -1) then begin
     LGs[sel].DiskSNIIRate = 0
 endif
 sel = testFloat(LGs.BulgeSNIIRate)
 if(sel(0) gt -1) then begin
     LGs[sel].BulgeSNIIRate = 0
 endif
 sel = testFloat(LGs.ICMSNIIRate)
 if(sel(0) gt -1) then begin
     LGs[sel].ICMSNIIRate = 0
 endif
 sel = testFloat(LGs.MassReturnRateToColdGas)
 if(sel(0) gt -1) then begin
     LGs[sel].MassReturnRateToColdGas = 0
 endif
 sel = testFloat(LGs.MassReturnRateToHotGas)
 if(sel(0) gt -1) then begin
     LGs[sel].MassReturnRateToHotGas = 0
 endif
 sel = testFloat(LGs.MetalsReturnRateToHotGas)
 if(sel(0) gt -1) then begin
     LGs[sel].MetalsReturnRateToHotGas = 0
 endif
 sel = testFloat(LGs.PrimordialAccretionRate)
 if(sel(0) gt -1) then begin
     LGs[sel].PrimordialAccretionRate = 0
 endif
 sel = testFloat(LGs.CoolingRadius)
 if(sel(0) gt -1) then begin
     LGs[sel].CoolingRadius = 0
 endif
 sel = testFloat(LGs.ReheatingRate)
 if(sel(0) gt -1) then begin
     LGs[sel].ReheatingRate = 0
 endif
 sel = testFloat(LGs.MetalsReheatingRate)
 if(sel(0) gt -1) then begin
     LGs[sel].MetalsReheatingRate = 0
 endif
 sel = testFloat(LGs.EjectionRate)
 if(sel(0) gt -1) then begin
     LGs[sel].EjectionRate = 0
 endif
 sel = testFloat(LGs.CoolingRate)
 if(sel(0) gt -1) then begin
     LGs[sel].CoolingRate = 0
 endif
 sel = testFloat(LGs.CoolingRate_beforeAGN)
 if(sel(0) gt -1) then begin
     LGs[sel].CoolingRate_beforeAGN = 0
 endif
 sel = testFloat(LGs.QuasarAccretionRate)
 if(sel(0) gt -1) then begin
     LGs[sel].QuasarAccretionRate = 0
 endif
 sel = testFloat(LGs.RadioAccretionRate)
 if(sel(0) gt -1) then begin
     LGs[sel].RadioAccretionRate = 0
 endif
 sel = testFloat(LGs.Sfr)
 if(sel(0) gt -1) then begin
     LGs[sel].Sfr = 0
 endif
 sel = testFloat(LGs.SfrBulge)
 if(sel(0) gt -1) then begin
     LGs[sel].SfrBulge = 0
 endif
 sel = testFloat(LGs.XrayLum)
 if(sel(0) gt -1) then begin
     LGs[sel].XrayLum = 0
 endif
 sel = testFloat(LGs.BulgeSize)
 if(sel(0) gt -1) then begin
     LGs[sel].BulgeSize = 0
 endif
 sel = testFloat(LGs.DiskRadius)
 if(sel(0) gt -1) then begin
     LGs[sel].DiskRadius = 0
 endif
 sel = testFloat(LGs.ColdGasRadius)
 if(sel(0) gt -1) then begin
     LGs[sel].ColdGasRadius = 0
 endif
 sel = testFloat(LGs.StellarHalfMassRadius)
 if(sel(0) gt -1) then begin
     LGs[sel].StellarHalfMassRadius = 0
 endif
 sel = testFloat(LGs.StellarHalfLightRadius)
 if(sel(0) gt -1) then begin
     LGs[sel].StellarHalfLightRadius = 0
 endif
 sel = testFloat(LGs.CosInclination)
 if(sel(0) gt -1) then begin
     LGs[sel].CosInclination = 0
 endif
 sel = testFloat(LGs.MagDust(0))
 if(sel(0) gt -1) then begin
     LGs[sel].MagDust(0) = 0
 endif
 sel = testFloat(LGs.MagDust(1))
 if(sel(0) gt -1) then begin
     LGs[sel].MagDust(1) = 0
 endif
 sel = testFloat(LGs.MagDust(2))
 if(sel(0) gt -1) then begin
     LGs[sel].MagDust(2) = 0
 endif
 sel = testFloat(LGs.MagDust(3))
 if(sel(0) gt -1) then begin
     LGs[sel].MagDust(3) = 0
 endif
 sel = testFloat(LGs.MagDust(4))
 if(sel(0) gt -1) then begin
     LGs[sel].MagDust(4) = 0
 endif
 sel = testFloat(LGs.Mag(0))
 if(sel(0) gt -1) then begin
     LGs[sel].Mag(0) = 0
 endif
 sel = testFloat(LGs.Mag(1))
 if(sel(0) gt -1) then begin
     LGs[sel].Mag(1) = 0
 endif
 sel = testFloat(LGs.Mag(2))
 if(sel(0) gt -1) then begin
     LGs[sel].Mag(2) = 0
 endif
 sel = testFloat(LGs.Mag(3))
 if(sel(0) gt -1) then begin
     LGs[sel].Mag(3) = 0
 endif
 sel = testFloat(LGs.Mag(4))
 if(sel(0) gt -1) then begin
     LGs[sel].Mag(4) = 0
 endif
 sel = testFloat(LGs.MagBulge(0))
 if(sel(0) gt -1) then begin
     LGs[sel].MagBulge(0) = 0
 endif
 sel = testFloat(LGs.MagBulge(1))
 if(sel(0) gt -1) then begin
     LGs[sel].MagBulge(1) = 0
 endif
 sel = testFloat(LGs.MagBulge(2))
 if(sel(0) gt -1) then begin
     LGs[sel].MagBulge(2) = 0
 endif
 sel = testFloat(LGs.MagBulge(3))
 if(sel(0) gt -1) then begin
     LGs[sel].MagBulge(3) = 0
 endif
 sel = testFloat(LGs.MagBulge(4))
 if(sel(0) gt -1) then begin
     LGs[sel].MagBulge(4) = 0
 endif
 sel = testFloat(LGs.MassWeightAge)
 if(sel(0) gt -1) then begin
     LGs[sel].MassWeightAge = 0
 endif
 sel = testFloat(LGs.rbandWeightAge)
 if(sel(0) gt -1) then begin
     LGs[sel].rbandWeightAge = 0
 endif
end
