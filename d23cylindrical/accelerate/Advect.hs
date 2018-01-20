{-# LANGUAGE Rank2Types #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE FunctionalDependencies #-}

module Advect where 

import Data.Array.Accelerate as Acc
import Data.Array.Accelerate.Linear

class Advect a b | a -> b where
    projection :: (Exp Double -> Exp Double -> Exp Double-> Exp Double) -> Acc a -> Acc a -> Acc a -> Acc (a,a)
    flux :: Exp (V3 Double) -> Acc a -> Acc b 

-- the TVD method looks at the flux and state on either side of an interface
-- and then determines the downstream flux
-- for hancock we just use the upstream and downstream flux as the upstream and 
-- downstream flux respectively
hancock :: forall a b. Exp (V3 Double) -> (a,a) -> (b,b) -> (a,a)
hancock _ (uf,df) _ = (uf,df) 


--
--ktscheme :: (Ord a,Floating a,Advect b) => (b->a) -> b -> b -> b
--ktscheme sc l r = scale (fromRational 0.5) $ fr + fl - scale s (r-l)
--                    where 
--                        s = max (sc l) (sc r)
--                        fr = flux r
--                        fl = flux l
--
--tvdmusclf :: Exp (V3 Double) -> (a,a) -> (b,b) -> (a,a) 
--tvdmusclf dir (uf, df) (us,ds) = (f,f)
--                    where 
--                        f = uf - df + scale (speed (average us ds)) (us - ds)
--

--these will have to be mapped over
--fluxsum :: (Floating a, Advect b) => AreaPair a -> FluxPair b -> State b 
--fluxsum (au,ad) (fu,fd) = (scale au fu) - (scale ad fd) 

--compute the state derivative given a cell geometry and the fluxes
--accumulate :: (Floating a,Advect b) => Geometry a -> VoxelFlux b -> State b
--accumulate geom flux = scale (recip volume) totalflux
--        where 
--            (faceareas,volume) = geom
--            dflux = liftA2 fluxsum faceareas flux 
--            totalflux = sum dflux

--two step performs a basic two step time integration
--it takes a geometry over which 
--twostep :: (Floating a,Advect b) => (State b -> VoxelFlux b) -> (State b -> VoxelFlux b) -> Geometry a -> Time a -> State b -> State b
--twostep predict main geom dt start = end  
--        where 
--        du_start = accumulate geom (predict start)
--        mid = start + scale (dt/2) start
--        du_mid = accumulate geom (main mid) 
--        end = start + scale dt du_mid

