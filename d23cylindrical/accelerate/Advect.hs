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
import Types 


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



cellcomp :: forall v state flux diff. 
            (P.Monad v,Elt state, Elt flux,Elt diff,Elt (v Double),
            Box v ((state,state),(state,state)),
            Box v (v Double, v Double),
            Elt (Patch v), 
            Box v (flux,flux),Box v Double) => Fluxer v state flux -> Differ v flux diff -> Exp (Geom v) -> Exp (v ((state,state),(state,state))) -> Exp diff
cellcomp fluxer differ geom states = derivative
                        where 
                            patch :: Exp (Patch v) 
                            patch = fst geom 
                            vol :: Exp Double
                            vol = snd geom
                            fluxes :: Exp (v (flux,flux))
                            fluxes = createFlux fluxer patch states 
                            derivative = differ vol fluxes 








