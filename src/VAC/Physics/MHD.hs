{-# LANGUAGE Rank2Types #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE TypeSynonymInstances #-}
{-# LANGUAGE FlexibleInstances #-}
module MHD where 

import qualified Prelude as P

import Types
import Data.Array.Accelerate as Acc 
import Data.Array.Accelerate.Linear
import Control.Lens 
import Limits

type Pressure = Precision 
type Magnetic = V3 Precision
type Momentum = V3 Precision
type Density = Precision
type Energy = Precision
type Force = V3 Precision
type Velocity = V3 Precision
type Speed = Precision 

magneticenergy :: Exp Magnetic -> Exp Energy
magneticenergy mag = (quadrance mag) / 2

kineticenergy :: Exp Momentum -> Exp Density -> Exp Energy
kineticenergy mom den = quadrance mom / 2.0 / den
 
velocity ::  Exp Momentum -> Exp Density -> Exp Velocity
velocity mom den = mom ^/ den  

momentum :: Exp Velocity -> Exp Density -> Exp Momentum
momentum v d = v^*d

pressure :: Exp Precision -> Exp Energy -> Exp Energy -> Exp Energy -> Exp Pressure
pressure gamma total kinetic magnetic = (gamma - 1)*(total - kinetic - magnetic)

thermalenergy :: Exp Precision -> Exp Pressure -> Exp Energy
thermalenergy gamma pressure = pressure / (gamma - 1)

totalPressure :: Exp Pressure -> Exp Energy -> Exp Pressure 
totalPressure p em = p + em 

csound2 ::  Exp Precision -> Exp Pressure -> Exp Density -> Exp Speed 
csound2 gamma pressure density = gamma * pressure / density

monoatomic_gamma :: Exp Precision 
monoatomic_gamma = constant (5/3)

cmax :: Exp (V3 Precision) -> Exp Density -> Exp Momentum -> Exp Energy  -> Exp Magnetic -> Exp Speed
cmax dir den mom ener mag = currvel + sound 
                        where
                            gamma = monoatomic_gamma
                            magd = dot dir $ mag 
                            vd = dot dir $ (velocity mom den) 
                            currvel = abs (vd/den)
                            c2 = csound2 gamma (pressure gamma ener (kineticenergy mom den) (magneticenergy mag)) den
                            alfen2 = (quadrance mag) / den
                            cfast2 = c2 + alfen2
                            sound = sqrt (0.5 * c2 + sqrt(c2 - 4*cfast2*(magd /den)))

type MHD = (Density, Momentum, Energy, Magnetic)

mhdzero :: MHD 
mhdzero = (0.0, V3 0.0 0.0 0.0,  0.0, V3 0.0 0.0 0.0)

mhdspeed :: Exp (V3 Precision) -> Exp MHD -> Exp Speed 
mhdspeed dir mhd = cmax dir den mom ener mag where 
        (den,mom,ener,mag) = unlift mhd 

mhdscale :: Exp Precision -> Exp MHD -> Exp MHD 
mhdscale amt mhd = lift (amt*den,amt*^mom,amt*ener,amt*^mag)
                        where (den,mom,ener,mag) = unlift mhd

mhdmerge :: Exp MHD -> Exp MHD -> Exp MHD 
mhdmerge a b = lift (aden+bden, amom ^+^ bmom, aener+bener, amag^+^bmag) 
                        where
                            aden :: Exp Density --not sure why these types are needed
                            aener :: Exp Energy
                            (aden,amom, aener,amag) = unlift a 
                            (bden,bmom, bener,bmag) = unlift b


mhdflux :: Exp (V3 Precision) -> Exp MHD -> Exp MHD
mhdflux dir mhd = lift (fden,fmom,fener,fmag)
                            where
                                (den,mom,ener,mag) = unlift mhd
                                gamma = monoatomic_gamma --plasma has no molecules
                                emag = magneticenergy mag 
                                ekin = kineticenergy mom den
                                ptot = totalPressure (pressure gamma ener ekin emag) emag
                                vel = velocity mom den
                                vd = dot dir vel
                                bd = dot dir mag 
                                fden = vd * den  
                                fmom = vd*^mom ^-^ bd*^mag ^+^ ptot*^dir
                                fener =  vd * (ener + den + emag) - bd * (dot mag vel)
                                fmag = (vd*^mag) ^-^ (bd*^vel)


--perform forwards and backwards projection with conservative variables
conservativeMHDProject :: (Exp Precision -> Exp Precision -> Exp Precision-> Exp Precision) -> Exp MHD -> Exp MHD -> Exp MHD -> Exp (MHD,MHD)
conservativeMHDProject limiter p c n = lift (u,d)
                        where
                            (pden,pmom,pener,pmag) = unlift p
                            (cden,cmom,cener,cmag) = unlift c
                            (nden,nmom,nener,nmag) = unlift n
                            projden = projectScalar limiter pden cden nden
                            projmom = projectVector limiter pmom cmom nmom
                            projener = projectScalar limiter pener cener nener
                            projmag = projectVector limiter pmag cmag nmag
                            u = (fst projden,fst projmom,fst projener,fst projmag)
                            d = (snd projden,snd projmom,snd projener,snd projmag)


--perform forwards and backwards projections with primative (pressure velocity) variables
primativeMHDProject :: (Exp Precision -> Exp Precision -> Exp Precision-> Exp Precision) -> Exp MHD -> Exp MHD -> Exp MHD -> Exp (MHD,MHD)
primativeMHDProject limiter p c n = lift (u,d)
                        where
                            gamma = monoatomic_gamma
                            (pden,pmom,pener,pmag) = unlift p
                            (cden,cmom,cener,cmag) = unlift c
                            (nden,nmom,nener,nmag) = unlift n
                            projden = projectScalar limiter pden cden nden
                            projmag = projectVector limiter pmag cmag nmag
                            projvel = projectVector limiter (velocity pmom pden) (velocity cmom cden) (velocity nmom nden)
                            uvel = fst projvel
                            dvel = snd projvel 
                            umom = momentum uvel (fst projden)
                            dmom = momentum dvel (snd projden) 
                            ppress = pressure gamma pener (kineticenergy pmom pden) (magneticenergy pmag)
                            cpress = pressure gamma cener (kineticenergy cmom cden) (magneticenergy cmag)
                            npress = pressure gamma nener (kineticenergy nmom nden) (magneticenergy nmag)
                            projpres  = projectScalar limiter ppress cpress npress 
                            upress = fst projpres 
                            dpress = snd projpres
                            uthermener = thermalenergy gamma (fst projpres)
                            dthermener = thermalenergy gamma (snd projpres)
                            utotener = uthermener + (kineticenergy umom (fst projden)) + magneticenergy (fst projmag)
                            dtotener = dthermener + (kineticenergy dmom (snd projden)) + magneticenergy (snd projmag)
                            u = (fst projden,umom,utotener,fst projmag)
                            d = (snd projden,dmom,dtotener,snd projmag)