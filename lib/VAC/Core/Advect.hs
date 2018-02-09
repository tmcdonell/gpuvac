{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE ConstraintKinds #-}
{-# LANGUAGE FlexibleContexts #-}
module VAC.Core.Advect where 

import qualified Prelude as P 

import Control.Lens
import Control.Applicative
import Data.Array.Accelerate as Acc
import Data.Array.Accelerate.Linear
import VAC.Core.Types as Types

-- the TVD method looks at the flux and state on either side of an interface
-- and then determines the downstream flux
-- for hancock we just use the upstream and downstream flux as the upstream and 
-- downstream flux respectively
hancock :: (Elt flux, Elt state) => FluxFunc state flux -> Fluxer state flux 
hancock f dir states = lift $ (f dir $ fst states, f dir $ snd states)


average :: (Elt flux, Elt state) =>  Merger flux -> Scaler flux -> FluxFunc state flux -> Fluxer state flux
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


onAxis :: Elt a => Projector a -> Stencil5 a -> Exp ((a,a),(a,a))
onAxis prj (pp,p,c,n,nn) = lift $ (us,ds) 
    where 
    pd =  snd $ prj pp p c
    cen =   prj p c n 
    nu =  fst $ prj c n nn
    us = (pd,fst cen)
    ds = (snd cen,nu)

stencilX :: Elt a => Projector a -> Stencil3x3x5 a -> Exp ((a,a),(a,a))
stencilX prj ( (_,(_,pp,_),_)
              ,(_,(_,p,_),_)
              ,(_,(_,c,_),_)
              ,(_,(_,n,_),_)
              ,(_,(_,nn,_),_) ) = onAxis prj (pp,p,c,n,nn)

stencilY :: Elt a => Projector a -> Stencil3x5x3 a -> Exp ((a,a),(a,a))
stencilY prj (_,((_,pp,_),(_,p,_),(_,c,_),(_,n,_),(_,nn,_)),_) = onAxis prj (pp,p,c,n,nn) 

stencilZ :: Elt a => Projector a -> Stencil5x3x3 a -> Exp ((a,a),(a,a))
stencilZ prj (_,(_,(pp,p,c,n,nn),_),_) = onAxis prj (pp,p,c,n,nn)

createFlux :: forall state flux.(Elt state,Elt flux) => Fluxer state flux -> Exp (V3 Precision,V3 Precision) -> Exp ((state,state),(state,state))  -> Exp (flux,flux)
createFlux fluxer geom projected_state = lift fluxes where 
        (ustates,dstates) = unlift projected_state :: (Exp (state,state),Exp (state,state))
        (unorm,dnorm) = unlift geom :: (Exp (V3 Precision),Exp (V3 Precision))
        uflux = snd $ fluxer unorm ustates :: Exp flux
        dflux = fst $ fluxer dnorm dstates :: Exp flux
        fluxes = (uflux,dflux)

advect:: forall state stencil flux sh. (Stencil sh state stencil,Elt state,Elt flux) => (stencil -> Exp ((state,state),(state,state))) -> Fluxer state flux -> Acc (Array sh PatchPair) -> Acc (Array sh state) -> Acc (Array sh (flux,flux))
advect sten fluxer geom input = fluxes where 
        states = stencil sten wrap input :: Acc (Array sh ((state,state),(state,state)))
        fluxes = Acc.zipWith (createFlux fluxer) geom states :: Acc (Array sh (flux,flux))

advection3D :: forall state flux diff. 
    (Elt state, Elt flux,Elt diff)=> Projector state -> Fluxer state flux -> Differ V3 flux diff -> Acc (Array DIM3 (Cell V3)) -> Acc (Array DIM3 state) -> Acc (Array DIM3 diff)
advection3D proj fluxer diff geom input = derivative 
                            where
                                (faces,volume) = unzip geom :: (Acc (Array DIM3 (Faces V3)),Acc (Array DIM3 Precision))
                                xfluxes = advect (stencilX proj) fluxer (map (\v -> v^._x) faces) input :: Acc (Array DIM3 (flux,flux))
                                yfluxes = advect (stencilY proj) fluxer (map (\v -> v^._y) faces) input :: Acc (Array DIM3 (flux,flux))
                                zfluxes = advect (stencilZ proj) fluxer (map (\v -> v^._z) faces) input :: Acc (Array DIM3 (flux,flux))
                                tobox a b c= lift $ V3 a b c :: Exp (V3 (flux,flux))
                                fluxes = Acc.zipWith3 tobox xfluxes yfluxes zfluxes :: Acc (Array DIM3 (V3 (flux,flux)))
                                derivative = Acc.zipWith diff volume fluxes :: Acc (Array DIM3 diff)




