{-# LANGUAGE Rank2Types #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE FunctionalDependencies #-}
module Integration where

import qualified Prelude as P 
import Data.Array.Accelerate as Acc

type Accumulator diff state = Exp Double -> Exp diff -> Exp state -> Exp state

twostep :: (Shape sh, Elt diff, Elt state) => Accumulator diff state -> (Acc (Array sh state) -> Acc (Array sh diff)) -> (Acc (Array sh state) -> Acc (Array sh diff)) -> Exp Double -> Acc (Array sh state) -> Acc (Array sh state) 
twostep integrate predict advance time start = final 
            where
                step t d u = Acc.zipWith (\di st -> integrate t di st) d u 
                du_start = predict start
                mid = step (time/2) du_start start  
                du_mid = advance mid
                final = step time du_mid start 