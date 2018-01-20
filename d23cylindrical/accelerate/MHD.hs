{-# LANGUAGE Rank2Types #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE TypeSynonymInstances #-}
{-# LANGUAGE FlexibleInstances #-}
module MHD where 

import qualified Prelude as P

import Data.Array.Accelerate as Acc 
import Data.Array.Accelerate.Linear
import Advect
import Control.Lens 
import Limits

type PressureField sh = Array sh Double 
type MagneticField sh = Array sh (V3 Double)
type MomentumField sh = Array sh (V3 Double) 
type DensityField sh = Array sh Double
type EnergyField sh = Array sh Double
type ForceField sh = Array sh (V3 Double) 
type VelocityField sh = Array sh (V3 Double)
type SpeedField sh = Array sh Double 

magneticenergy :: Shape sh => Acc (MagneticField sh) -> Acc (EnergyField sh)
magneticenergy mag = Acc.map (\b -> (quadrance b) / 2 ) mag

kineticenergy :: Shape sh => Acc (MomentumField sh) -> Acc (DensityField sh) -> Acc (EnergyField sh)
kineticenergy mom den = Acc.zipWith (\m d -> (quadrance m) / 2.0 / d) mom den
 
velocity :: Shape sh => Acc (MomentumField sh) -> Acc (DensityField sh) -> Acc (VelocityField sh) 
velocity mom den = Acc.zipWith (\m d -> m ^/ d) mom den 

pressure :: Shape sh => Exp Double -> Acc (EnergyField sh) -> Acc (EnergyField sh) -> Acc (EnergyField sh) -> Acc (PressureField sh) 
pressure gamma total kinetic magnetic = 
    Acc.zipWith3 (\ t k m -> (gamma - 1)*(t - k - m)) total kinetic magnetic

totalPressure :: Shape sh => Acc (PressureField sh) -> Acc (EnergyField sh) -> Acc (PressureField sh) 
totalPressure p em = Acc.zipWith (+) p em 

csound2 :: Shape sh => Exp Double -> Acc (PressureField sh) -> Acc (DensityField sh) -> Acc (SpeedField sh) 
csound2 gamma pressure density = Acc.zipWith (\p d -> gamma * p / d) pressure density

monoatomic_gamma :: Exp Double 
monoatomic_gamma = constant (5/3)

cmax :: Shape sh =>Exp (V3 Double) ->Acc (DensityField sh) -> Acc (MomentumField sh) -> Acc (EnergyField sh)  -> Acc (MagneticField sh) -> Acc (SpeedField sh)
cmax dir den mom ener mag = Acc.zipWith (+) currvel sound 
                        where
                            gamma = monoatomic_gamma
                            magd = Acc.map (\b -> dot b dir) mag
                            vd = Acc.map (\v -> dot v dir) (velocity mom den) 
                            currvel = Acc.zipWith (\v d -> abs (v/d)) vd den
                            c2 = csound2 gamma (pressure gamma ener (kineticenergy mom den) (magneticenergy mag)) den
                            alfen2 = Acc.zipWith (\ b d -> (quadrance b) / d) mag den
                            cfast2 = Acc.zipWith (+) c2 alfen2
                            sound = Acc.zipWith4 (\c cf b d ->  sqrt (0.5 * cf + sqrt(cf - 4*c*(b /d)))) cfast2 c2 magd den

type MHD sh = (DensityField sh, MomentumField sh, EnergyField sh, MagneticField sh)

mhdflux :: Shape sh => (Exp (V3 Double)) -> Acc (MHD sh) -> Acc (MHD sh)
mhdflux dir mhd = lift (fden,fmom,fener,fmag)
                            where
                                (den,mom,ener,mag) = unlift mhd
                                gamma = monoatomic_gamma --plasma has no molecules
                                emag = magneticenergy mag 
                                ekin = kineticenergy mom den
                                ptot = totalPressure (pressure gamma ener ekin emag) emag
                                vel = velocity mom den
                                fden = Acc.zipWith (\v d -> (dot dir v) * d) vel den 
                                fmom = Acc.zipWith4 (\v m b p -> (dot dir v)*^m ^-^ (dot dir b)*^b ^+^ p*^dir) vel mom mag ptot
                                fener =  Acc.zipWith5 (\v e b d em -> (dot v dir) * (e + d + em) - (dot b dir) * (dot b v)) vel ener mag den emag
                                fmag = Acc.zipWith (\b v -> (dot dir v)*^b ^-^ (dot b dir)*^v) mag vel


projectMHDArray :: Shape sh => (Exp Double -> Exp Double -> Exp Double-> Exp Double) -> Acc (MHD sh) -> Acc (MHD sh) -> Acc (MHD sh) -> Acc (MHD sh, MHD sh) 
projectMHDArray limiter p c n = lift (u,d)
                        where
                            (pden,pmom,pener,pmag) = unlift p
                            (cden,cmom,cener,cmag) = unlift c
                            (nden,nmom,nener,nmag) = unlift n
                            projden = projectScalarArray limiter pden cden nden
                            projmom = projectVectorArray limiter pmom cmom nmom
                            projener = projectScalarArray limiter pener cener nener
                            projmag = projectVectorArray limiter pmag cmag nmag
                            u = (afst projden,afst projmom,afst projener,afst projmag)
                            d = (asnd projden,asnd projmom,asnd projener,asnd projmag)
