{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE ConstraintKinds #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE RankNTypes #-}
module VAC.Core.Advect where 

import qualified Prelude as P 

import Control.Lens
import Control.Applicative
import Data.Array.Accelerate as Acc
import Data.Array.Accelerate.Linear
import VAC.Core.Slicing as Slicing
import VAC.Core.Types as Types

-- the TVD method looks at the flux and state on either side of an interface
-- and then determines the downstream flux
-- for hancock we just use the upstream and downstream flux as the upstream and 
-- downstream flux respectively
hancock :: (Elt flux, Elt state) => FluxFunc state flux -> Fluxer state flux 
hancock f dir stateL stateR = lift $ (f dir stateL, f dir stateR)


average :: (Elt flux, Elt state) =>  Merger flux -> Scaler flux -> FluxFunc state flux -> Fluxer state flux
average m s f dir stateL stateR = lift (avg,avg) where
                                left = f dir stateL
                                right = f dir stateR
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


onDim :: forall state flux sh. (Shape sh,Elt state,Elt flux) => Lens' (Exp sh) (Exp Int) 
    -> Projector state -> Fluxer state flux -> Acc (Array sh (V3 Precision)) 
    -> Acc (Array sh state)  -> (Acc (Array sh flux),Acc (Array sh flux))
onDim dim proj fluxr geom states  = (left dim fluxR,right dim fluxL)
    where
        --previous, current, and next states along the chosen dimension 
        (prev,curr,next) = unlift $ sten3 dim states :: (Acc (Array sh state), Acc (Array sh state), Acc (Array sh state))
        --up upwind and downwind projected states at the current index
        (upwind, dwnwind) = unzip $ Acc.zipWith3 proj prev curr next :: (Acc (Array sh state),Acc (Array sh state))
        --trim two of either side of the geometry to account for the slicing
        (fluxL,fluxR) = unzip $ Acc.zipWith3 fluxr (trim dim 2 geom) (left dim dwnwind) (right dim upwind) :: (Acc (Array sh flux), Acc (Array sh flux))
        
        


advection3D :: forall state flux diff. 
    (Elt state, Elt flux,Elt diff)=> Projector state -> Fluxer state flux 
    -> Differ V3 flux diff -> Acc (Array DIM3 (Cell V3)) -> Acc (Array DIM3 state) 
    -> Acc (Array DIM3 diff)
advection3D proj fluxer diff geom input = derivative 
                            where
                                (faces,volume) = unzip geom :: (Acc (Array DIM3 (Faces V3)),Acc (Array DIM3 Precision))
                                xfluxes = onDim _1 proj fluxer xfaces input 
                                yfluxes = onDim _2 proj fluxer yfaces input 
                                zfluxes = onDim _3 proj fluxer zfaces input 
                                
                                tobox a b c= lift $ V3 a b c :: Exp (V3 (flux,flux))
                                fluxes = Acc.zipWith3 tobox xfluxes yfluxes zfluxes :: Acc (Array DIM3 (V3 (flux,flux)))
                                derivative = Acc.zipWith diff volume fluxes :: Acc (Array DIM3 diff)




