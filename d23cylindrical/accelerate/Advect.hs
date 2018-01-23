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

-- the TVD method looks at the flux and state on either side of an interface
-- and then determines the downstream flux
-- for hancock we just use the upstream and downstream flux as the upstream and 
-- downstream flux respectively
hancock :: (Elt flux, Elt state) => FluxFunc v state flux -> Fluxer v state flux 
hancock f dir states = lift $ (f dir $ fst states, f dir $ snd states)


average :: (Elt flux, Elt state) =>  Merger flux -> Scaler flux -> FluxFunc v state flux -> Fluxer v state flux
average m s f dir states = lift (avg,avg) where
                                left = f dir $ fst states 
                                right = f dir $ snd states
                                avg = s 0.5 $ m left right 


diff :: forall flux diff v. 
    (P.Foldable v,P.Monad v, Elt flux,
    Elt diff, Box v (flux,flux),Box v flux,
    Box v diff)=>
    Merger diff -> Flow flux diff -> Differ v flux diff 
diff m f v fluxes = P.foldl1 m net where 
        dims = unlift fluxes :: v (Exp (flux,flux))
        net :: v (Exp diff) 
        net = do 
            d <- dims
            let inflow = f (recip v) $ fst d 
            let outflow = f (negate $ recip v) $ snd d
            P.return $ m inflow outflow


accum :: Merger state -> Scaler state -> Accumulator state state
accum m s t d i =  m i ds where ds = s t d


proj :: Elt a => Projector a -> Stencil5 a -> Exp ((a,a),(a,a))
proj prj (pp,p,c,n,nn) = lift $ (us,ds) 
    where 
    pd =  snd $ prj pp p c
    cen =   prj p c n 
    nu =  fst $ prj c n nn
    us = (pd,fst cen)
    ds = (snd cen,nu)

stencil1D :: Elt state => Projector state -> Stencil5 state -> Exp (V1 ((state,state),(state,state)))
stencil1D prj vals = lift . V1 $ proj prj vals

stencil2D :: Elt state => Projector state ->Stencil5x5 state -> Exp (V2 ((state,state),(state,state)))
stencil2D prj (ppx,px,c,nx,nnx) = lift $ V2 dn (d1^._x) 
            where 
                d1 = unlift $ stencil1D prj c 
                dn = proj prj (ppx^._3,px^._3,c^._3,nx^._3,nnx^._3) 

stencil3D :: Elt state => Projector state -> Stencil5x5x5 state -> Exp (V3 ((state,state),(state,state)))
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
            Box v (flux,flux),Box v Double) => Fluxer v state flux -> Differ v flux diff -> Exp (Cell v) -> Exp (v ((state,state),(state,state))) -> Exp diff
cellcomp fluxer differ geom states = derivative
                        where 
                            patch :: Exp (Patch v) 
                            patch = fst geom 
                            vol :: Exp Double
                            vol = snd geom
                            fluxes :: Exp (v (flux,flux))
                            fluxes = createFlux fluxer patch states 
                            derivative = differ vol fluxes 


advection3D :: forall  state flux diff. 
    (Elt state, Elt flux,Elt diff)=> Projector state -> Fluxer V3 state flux -> Differ V3 flux diff -> Acc (Array DIM3 (Cell V3)) -> Acc (Array DIM3 state) -> Acc (Array DIM3 diff)
advection3D proj flux diff geom input = derivative 
                            where 
                                projected_states = stencil (stencil3D proj) wrap input :: Acc (Array DIM3 (V3 ((state,state),(state,state))))
                                cellcomputer g s = cellcomp flux diff g s :: Exp diff
                                derivative = Acc.zipWith cellcomputer geom projected_states :: Acc (Array DIM3 diff) 




