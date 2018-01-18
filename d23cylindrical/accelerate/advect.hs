{-# LANGUAGE Rank2Types #-}
module Advect where 

import Control.Applicative
import Control.Lens 
import Linear

type AreaPair a = (a,a)
type Volume a = a
type Geometry a = (V3 (AreaPair a),Volume a)
type State a = a
type ProjectedState a = (a,a)  
type FluxPair a = (a,a)
type VoxelFlux a = V3 (FluxPair a)
type Time a = a 

class Num a => Advect a where 
    --given a cell state, compute the boundry flux (need to include direction)
    flux :: a -> a
    zeroflux :: a 
    --given a cell state, compute the spectral radius (max eignevalue) of df(u)/du
    --typically maximum speed along direction 
    speed :: Floating b => (Lens' (V3 b) b) -> a -> b
    primative :: a -> a 
    conserve :: a -> a
    scale :: Floating b => b -> a -> a 

extend :: (Advect a) => ((a,a,a)->a) -> (a,a,a) -> (a,a)
extend limiter (p,c,n) = (conserve lp,conserve rp) 
                        where
                            slope = limiter (primative p,primative c,primative n)
                            step = scale (fromRational 0.5) slope
                            lp = c - step 
                            rp = c + step

hancock :: (Floating b, Advect a) => (Lens' (V3 b) b) -> (a,a) -> (a,a) -> (a,a) 
hancock _ (_,upwindr) (dwnwindl,_) =  (flux upwindr,flux dwnwindl) 

ktscheme :: (Ord a,Floating a,Advect b) => (b->a) -> b -> b -> b
ktscheme sc l r = scale (fromRational 0.5) $ fr + fl - scale s (r-l)
                    where 
                        s = max (sc l) (sc r)
                        fr = flux r
                        fl = flux l

tvdmusclf :: (Ord b,Floating b, Advect a) => (Lens' (V3 b) b) -> (a,a) -> (a,a) -> (a,a) 
tvdmusclf dir (ul, ur) (dl,dr) = (uf,df)
    where 
          sc = speed dir 
          uf = ktscheme sc ul ur 
          df = ktscheme sc dl dr

--these will have to be mapped over
fluxsum :: (Floating a, Advect b) => AreaPair a -> FluxPair b -> State b 
fluxsum (au,ad) (fu,fd) = (scale au fu) - (scale ad fd) 

--compute the state derivative given a cell geometry and the fluxes
accumulate :: (Floating a,Advect b) => Geometry a -> VoxelFlux b -> State b
accumulate geom flux = scale (recip volume) totalflux
        where 
            (faceareas,volume) = geom
            dflux = liftA2 fluxsum faceareas flux 
            totalflux = sum dflux

--two step performs a basic two step time integration
--it takes a geometry over which 
twostep :: (Floating a,Advect b) => (State b -> VoxelFlux b) -> (State b -> VoxelFlux b) -> Geometry a -> Time a -> State b -> State b
twostep predict main geom dt start = end  
        where 
        du_start = accumulate geom (predict start)
        mid = start + scale (dt/2) start
        du_mid = accumulate geom (main mid) 
        end = start + scale dt du_mid