{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE ConstraintKinds #-}
{-# LANGUAGE FlexibleContexts #-}
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
stencil2D prj (ppx,px,c,nx,nnx) = lift $ V2 dn (d1^._x) 
            where 
                d1 = unlift $ stencil1D prj c 
                dn = proj prj (ppx^._3,px^._3,c^._3,nx^._3,nnx^._3) 

stencil3D :: Elt a => Projector a -> Stencil5x5x5 a -> Exp (V3 ((a,a),(a,a)))
stencil3D prj (ppx,px,c,nx,nnx) = lift $ V3 dn (d2^._x) (d2^._y) 
            where
                d2 = unlift $ stencil2D prj c
                dn = proj prj (ppx^._3^._3,px^._3^._3,c^._3^._3,nx^._3^._3,nnx^._3^._3)

-- A fluxer takes a quantity dA which is a vector with a magnitude equal to the cell boundry area
-- and the direction normal to the surface in the direction of the downstream cell 
-- it also takes the quanties of the left and right state at the cell boundry.
-- The result of this computation is the left and right fluxes which will be the 
-- same for a conservative method.
type Fluxer vec state flux = Exp (vec Double) -> Exp (state,state) -> Exp (flux,flux) 

type Patch order  = order (order Double, order Double) 
type Geom order  = (Patch order,Double)


createFlux :: forall v state flux. 
                (P.Monad v,Elt state, Elt flux,Elt (v Double),
                Box v ((state,state),(state,state)),
                Box v (v Double, v Double),
                Box v (flux,flux),Box v Double) =>
                Fluxer v state flux -> Exp (Patch v) -> Exp (v ((state,state),(state,state))) -> Exp (v (flux,flux))
createFlux fluxer cellgeom projected_state = lift fluxes where
                    norms :: v (Exp (v Double, v Double))
                    norms = unlift cellgeom
                    states :: v (Exp ((state,state),(state,state)))
                    states = unlift projected_state
                    fluxes :: v (Exp (flux,flux))
                    fluxes = do 
                        n <- norms 
                        s <- states
                        let leftn = fst n
                        let rightn = snd n
                        let lefts = fst s
                        let rights = snd s
                        let leftf = fluxer leftn lefts
                        let rightf = fluxer rightn rights
                        let uf = snd leftf
                        let df = fst rightf
                        P.return $ lift (uf,df)



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



