--Copyright 2018 GeneralFusion Inc.
--
--Licensed under the Apache License, Version 2.0 (the "License");
--you may not use this file except in compliance with the License.
--You may obtain a copy of the License at
--
--    http://www.apache.org/licenses/LICENSE-2.0
--
--Unless required by applicable law or agreed to in writing, software
--distributed under the License is distributed on an "AS IS" BASIS,
--WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
--See the License for the specific language governing permissions and
--limitations under the License.
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE ConstraintKinds #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE RankNTypes #-}
module VAC.Core.Advect where 

import qualified Prelude as P 


import Data.Array.Accelerate as Acc
import Data.Array.Accelerate.Control.Lens
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
    -> Acc (Array sh state)  -> Acc (Array sh (flux,flux))
onDim dim proj fluxr geom states  = zip (left dim fluxR) (right dim fluxL)
    where
        --previous, current, and next states along the chosen dimension 
        (prev,curr,next) = unlift $ sten3 dim states :: (Acc (Array sh state), Acc (Array sh state), Acc (Array sh state))
        --up upwind and downwind projected states at the current index
        (upwind, dwnwind) = unzip $ Acc.zipWith3 proj prev curr next :: (Acc (Array sh state),Acc (Array sh state))
        --trim two of either side of the geometry to account for the slicing
        (fluxL,fluxR) = unzip $ Acc.zipWith3 fluxr (trim dim 2 geom) (left dim dwnwind) (right dim upwind) :: (Acc (Array sh flux), Acc (Array sh flux))
        
        
modshape :: Shape sh => sh -> sh -> sh
modshape Z Z = Z
modshape (e:.len) (b:.el) = (modshape e b):.(mod (el-2) (len-4))

advection3D :: forall state flux diff. 
    (Elt state, Elt flux,Elt diff)=> Projector state -> Fluxer state flux 
    -> Differ V3 flux diff -> Acc Geometry3D -> Acc (Array DIM3 state) 
    -> Acc (Array DIM3 diff)
advection3D proj fluxer diff geom input = backpermute (shape input) (modshape (shape input)) derivative 
                            where
                                (volume,faces1,faces2,faces3) = unlift geom :: (Acc (Array DIM3 Precision),Acc (Array DIM3 (V3 Precision)),Acc (Array DIM3 (V3 Precision)),Acc (Array DIM3 (V3 Precision)))
                                fluxfunc :: Lens' (Exp DIM3) (Exp Int) -> Acc (Array DIM3 (V3 Precision)) -> Acc (Array DIM3 (flux,flux))
                                fluxfunc lens faces = onDim lens proj fluxer faces input
                                fluxes1 = fluxfunc _1 faces1 :: Acc (Array DIM3 (flux,flux))
                                fluxes2 = fluxfunc _2 faces2 :: Acc (Array DIM3 (flux,flux))
                                fluxes3 = fluxfunc _3 faces3 :: Acc (Array DIM3 (flux,flux))
                                
                                tobox a b c= lift $ V3 a b c :: Exp (V3 (flux,flux))
                                fluxes = Acc.zipWith3 tobox fluxes1 fluxes2 fluxes3 :: Acc (Array DIM3 (V3 (flux,flux)))
                                derivative = Acc.zipWith diff volume fluxes :: Acc (Array DIM3 diff)
                               



