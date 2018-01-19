{-# LANGUAGE Rank2Types #-}
{-# LANGUAGE ScopedTypeVariables #-}
module Advect where 

import Control.Applicative
import Control.Lens 
--import Linear

import Data.Array.Accelerate as Acc
import Data.Array.Accelerate.Linear

type AreaPair a = (a,a)
type Volume a = a
type Geometry a = (V3 (AreaPair a),Volume a)
type State a = a
type ProjectedState a = (a,a)  
type FluxPair a = (a,a)
type VoxelFlux a = V3 (FluxPair a)
type Time a = a 

data Step a b c = Step {
    projstate ::  (a c,a c,a c) -> (a c,a c), 
    toflux::  Lens' (V3 c) c -> a c -> b c,
    zeroflux:: b c, 
    shock ::  Lens' (V3 c) c -> (a c,a c) -> b c,
    merge ::  (b c,b c,b c) -> (b c,b c,b c) -> (b c,b c),
    fluxfix ::  V3 (b c, b c) -> V3 (b c, b c),
    derivative :: V3 (b c) -> a c}


advect:: forall a b c. Lens' (V3 c) c -> Step a b c -> ( a c,a c, a c, a c, a c) -> (b c, b c) 
advect lens prog (pp,p,c,n,nn) = res
        where
            proj :: (a c, a c, a c) -> (a c, a c)
            proj =  (projstate prog)
            (_,pds) = proj (pp,p,c) --upwind reconstructed states 
            (cus,cds) = proj (p,c,n)  --current reconstructed states
            (nus,_) = proj (c,n,nn) --downwind reconstructed states 
            getshock = shock prog
            usf = getshock lens (pds,cus) --shock fluxes associated with boundry states upstream
            dsf = getshock lens (cds,nus) --shock fluxes associated with boundry states downstream
            flux = toflux prog $ lens 
            res = lift $ (merge prog) (flux pds, usf,flux cus) (flux cds,dsf,flux nus) 


advect1 :: Step a b c ->  Boundary (Array DIM1 (a c)) -> Acc (Array DIM1 (a c)) -> Acc (Array DIM1 (VoxelFlux (b c)))
advect1 prog boundry input = Acc.map filter fluxbundle
            where 
                xflux :: Acc (Array DIM1 (FluxPair (b c)))
                xflux = stencil (lift . (advect _x prog) . unlift) boundry input
                zf :: b c
                zf = (zeroflux prog)
                toflux :: Exp (FluxPair (b c)) -> Exp (V3 (FluxPair (b c)))
                toflux xf = lift $ V3 xf zf zf
                fluxbundle :: Acc (Array DIM1 (VoxelFlux (b c)))
                fluxbundle = Acc.map toflux xflux
                filter:: Exp (VoxelFlux (b c)) -> Exp (VoxelFlux (b c))
                filter = (fluxfix prog) 



--class Advect where 
--    speed :: Floating b => (Lens' (V3 b) b) -> a -> b
--
--class Num a => Scale a where 
--    scale :: Floating b => b -> a -> a 
--
--class Primative a where
--    primative :: State b => a -> b 
--    conserve :: State b => b -> a
--
--extend :: (Primative a,Scale a) => ((a,a,a)->a) -> (a,a,a) -> (a,a)
--extend limiter (p,c,n) = (conserve lp,conserve rp) 
--                        where
--                            slope = limiter (primative p,primative c,primative n)
--                            step = scale (fromRational 0.5) slope
--                            lp = c - step 
--                            rp = c + step
--
--
--hancockmerge :: (a,a,a) -> (a,a,a) -> (a,a) 
--hancockmerge (ul,um,ur) (dl,dm,dr) = (ur,dl) 
--
--hancockshock :: Lens' (V3 c) c -> a c -> b c 
--hancockshock _ _ = zeroflux
--
--hancock :: (Floating b, Advect a) => (Lens' (V3 b) b) -> (a,a) -> (a,a) -> (a,a) 
--hancock _ (_,upwindr) (dwnwindl,_) =  (flux upwindr,flux dwnwindl) 
--
--ktscheme :: (Ord a,Floating a,Advect b) => (b->a) -> b -> b -> b
--ktscheme sc l r = scale (fromRational 0.5) $ fr + fl - scale s (r-l)
--                    where 
--                        s = max (sc l) (sc r)
--                        fr = flux r
--                        fl = flux l
--
--tvdmusclf :: (Ord b,Floating b, Advect a) => (Lens' (V3 b) b) -> (a,a) -> (a,a) -> (a,a) 
--tvdmusclf dir (ul, ur) (dl,dr) = (uf,df)
--    where 
--          sc = speed dir 
--          uf = ktscheme sc ul ur 
--          df = ktscheme sc dl dr

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

