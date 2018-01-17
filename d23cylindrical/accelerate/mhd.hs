{-# LANGUAGE Rank2Types #-}
module MHD where 

import Linear
import Control.Lens 

type Momentum a = V3 a
type Density a = a 
type Velocity a = V3 a

type Energy a = a
type Pressure a = a 

type Fluid a = (Momentum a, Density a) 
type PrimFluid a = (Velocity a, Density a)

velocity :: Fractional a => Fluid a -> Velocity a 
velocity (momentum, density) = momentum ^/ density

density :: Fractional a => Fluid a -> Density a 
density (_,density) = density

primFluid :: Fractional a => Fluid a -> PrimFluid a 
primFluid fluid = (velocity fluid, density) where (_,density) = fluid 

consFluid :: Fractional a => PrimFluid a -> Fluid a 
consFluid (velocity, density) = (velocity ^* density, density) 

--quadrance is their funny name for length squared
ekin:: Fractional a => Fluid a -> Energy a 
ekin (momentum, density) = lensq / 2.0 / density
                            where lensq = quadrance momentum

type Magnetic a = V3 a

type MHD1 a = (Fluid a, Energy a, Magnetic a)
type PrimMHD1 a = (PrimFluid a, Pressure a, Magnetic a)

emag:: Fractional a => Magnetic a -> Energy a 
emag mag = lensq / 2.0 where lensq = quadrance mag


pressurePlasma1 :: Fractional a => a -> MHD1 a -> Pressure a 
pressurePlasma1 gamma (fluid, energy, mag) = (gamma - 1.0) * energy - ekin fluid - emag mag 

pressureTotal1 :: Fractional a => a -> MHD1 a -> Pressure a  
pressureTotal1 gamma mhdFluid = pressurePlasma1 gamma mhdFluid + emag mag 
                                    where (_,_,mag) = mhdFluid 

primMHD1 :: Fractional a => a -> MHD1 a -> PrimMHD1 a
primMHD1 gamma mhdFluid = (primfluid, pressure, magnetic)
                            where
                                (fluid,_,magnetic) = mhdFluid
                                primfluid = primFluid fluid
                                pressure = pressurePlasma1 gamma mhdFluid

consMHD1 :: Fractional a => a -> PrimMHD1 a -> MHD1 a
consMHD1 gamma primMHDFluid = (fluid, energy, magnetic) 
                                where 
                                    (primFluid, pressure, magnetic) = primMHDFluid
                                    fluid = consFluid primFluid
                                    energy = pressure / gamma + ekin fluid + emag magnetic


squareSound:: Fractional a => a -> MHD1 a -> a 
squareSound gamma mhd = gamma * (pressurePlasma1 gamma mhd) / density fluid 
                            where (fluid, _, _) = mhd 

cmax :: Floating a => a -> (Lens' (V3 a) a) -> MHD1 a -> a
cmax gamma l mhdFluid = currvel + sound 
                        where
                            (fluid,energy,magnetic) = mhdFluid
                            bd = magnetic ^. l --magnetic field along direction of travel
                            vd = (velocity fluid) ^. l --velocity along direction of travel
                            dens = density fluid
                            currvel = abs (vd / dens)
                            csound2 = squareSound gamma mhdFluid
                            cfast2 = csound2 + (quadrance magnetic) / dens
                            sound = sqrt (0.5 * cfast2 + sqrt(cfast2 - 4*csound2*(bd / dens)))

-- perp sets the direction of the vector along the lens to zero
-- this means that the remainder of the vector is perpendicular to 
-- the lens axis
perp :: Num a => (Lens' (f a) a) -> f a -> f a
perp l v = v & l .~ 0


flux :: Fractional a => a -> (Lens' (V3 a) a) -> MHD1 a -> MHD1 a 
flux gamma dir mhd = ((momflux,denflux),enerflux,magflux)
                            where 
                                (fluid, ener, mag) = mhd
                                (mom,den) = fluid
                                v = mom ^/ den -- velocity vector
                                vi = v^.dir  -- velocity along direction of derivative
                                bi = mag ^.dir  -- magnetic component along direction of derivative
                                n = zero & dir .~ 1.0 -- unit vector along direction of derivative
                                pt = pressureTotal1 gamma mhd
                                denflux = vi*den
                                momflux = vi * den *^ v ^-^ bi *^ mag + den *^ n 
                                enerflux = vi*ener - bi  * (dot mag v) + vi*den
                                magflux = vi *^ mag - bi *^ v


--TODO: state projection and flux functions

