{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE ConstraintKinds #-}
module Advect where 

import qualified Prelude as P 

import Control.Lens
import Control.Applicative
import Data.Array.Accelerate as Acc
import Data.Array.Accelerate.Linear

--class (Elt diff, Elt flux) => Accumulate flux diff | diff -> flux where 
--    accumulate :: Exp Double -> Exp flux -> Exp diff 
--
--class Elt diff => Merge diff  where 
--    merge :: Exp diff -> Exp diff -> Exp diff  

type Projector a = Exp a -> Exp a -> Exp a -> Exp (a,a)

proj :: Elt a => Projector a -> Stencil5 a -> Exp ((a,a),(a,a))
proj prj (pp,p,c,n,nn) = lift $ (us,ds) 
    where 
    pd =  snd $ prj pp p c
    cen =   prj p c n 
    nu =  fst $ prj c n nn
    us = (pd,fst cen)
    ds = (snd cen,nu)

stencil1D :: Elt a => Projector a -> Stencil5 a -> Exp (V1 ((a,a),(a,a)))
stencil1D prj vals = lift . V1 $ proj prj vals

stencil2D :: Elt a => Projector a ->Stencil5x5 a -> Exp (V2 ((a,a),(a,a)))
stencil2D prj (ppy,py,c,ny,nny) = lift $ V2 (d1^._x) dn 
            where 
                d1 = unlift $ stencil1D prj c 
                dn = proj prj (ppy^._3,py^._3,c^._3,ny^._3,nny^._3) 

stencil3D :: Elt a => Projector a -> Stencil5x5x5 a -> Exp (V3 ((a,a),(a,a)))
stencil3D prj (ppz,pz,c,nz,nnz) = lift $ V3 (d2^._x) (d2^._y) dn 
            where
                d2 = unlift $ stencil2D prj c
                dn = proj prj (ppz^._3^._3,pz^._3^._3,c^._3^._3,nz^._3^._3,nnz^._3^._3)

-- A fluxer takes a quantity dA which is a vector with a magnitude equal to the cell boundry area
-- and the direction normal to the surface in the direction of the downstream cell 
-- it also takes the quanties of the left and right state at the cell boundry.
-- The result of this computation is the left and right fluxes which will be the 
-- same for a conservative method.
type Fluxer vec state flux = Exp (vec Double) -> Exp (state,state) -> Exp (flux,flux) 

type Patch order dir  = order (dir Double, dir Double) 
type Geom order dir  = (Patch order dir,Double)


type Blox f a = (Elt (f a), Unlift Exp (f (Exp a)), Plain (f (Exp a)) ~ f a)

createFlux :: forall v v2 state flux.
       ( Applicative v, Elt state, Elt flux
       , Blox v ((state,state), (state,state))
       , Blox v (v2 Double, v2 Double)
       , Blox v (flux,flux)
       , Blox v2 Double
       ) => Fluxer v2 state flux -> Exp (Patch v v2) -> Exp (v ((state,state),(state,state))) -> Exp (v (flux,flux))
createFlux fluxer cellgeom projected_state = lift fluxes where
                    norms :: v (Exp (v2 Double, v2 Double))
                    norms = unlift cellgeom
                    leftnorm :: v (Exp (v2 Double))
                    leftnorm = P.fmap fst norms 
                    rightnorm :: v (Exp (v2 Double))
                    rightnorm = P.fmap snd norms
                    states :: v (Exp ((state,state),(state,state)))
                    states = unlift projected_state
                    leftstates :: v (Exp (state,state))
                    leftstates = P.fmap fst states
                    rightstates :: v (Exp (state,state))
                    rightstates = P.fmap snd states 
                    leftfluxes :: v (Exp (flux,flux))
                    leftfluxes = liftA2 fluxer leftnorm leftstates
                    rightfluxes :: v (Exp (flux,flux))
                    rightfluxes = liftA2 (fluxer) rightnorm rightstates
                    uf :: v (Exp flux) 
                    uf = P.fmap fst leftfluxes
                    df :: v (Exp flux) 
                    df = P.fmap snd rightfluxes
                    fluxes :: v (Exp (flux,flux))
                    fluxes = liftA2 (\x y -> lift (x,y)) uf df


----core of the advection flux calculation algorithm
--advect :: (Interp state, Flux state fl) => Exp (V3 Double, V3 Double) -> (Exp state, Exp state, Exp state, Exp state,Exp state) -> Exp (fl,fl) 
--advect areas state = lift (uflux,dflux) where
--                    (pp,p,c,n,nn) = state
--                    (_,dp) = step pp p c
--                    (uc,dc) =  step p c n
--                    (un,_) =  step c n nn
--                    unorm :: Exp (V3 Double)
--                    dnorm :: Exp (V3 Double) 
--                    (unorm,dnorm) =  unlift $ areas
--                    ustate =  lift (dp,uc)
--                    dstate =  lift (dc,un) 
--                    (_,uflux) =  unlift $ flux unorm ustate
--                    (dflux,_) =   unlift $ flux dnorm dstate

--advect1D :: (Interp state, Flux state flux) => Acc (Array DIM1 (V3 (V3 Double,V3 Double))) -> Boundary (Array DIM1 state) -> Acc (Array DIM1 state) -> Acc (Array DIM1 (V1 (flux,flux))) 
--advect1D geom bound input = fl where
--                            area1 :: Acc (Array DIM1 (V3 Double, V3 Double))
--                            area1 = Acc.map (\v -> v^._x) geom
--                            sten = stencil stencil1D bound 
--                            st = sten input
--                            fl1 = advect area1 st
--                            fl = Acc.map (lift.V1) fl1 



