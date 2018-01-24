{-# LANGUAGE Rank2Types #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE FunctionalDependencies #-}
module Integration where

import qualified Prelude as P 
import Data.Array.Accelerate as Acc

import Types 

accum :: Merger state -> Scaler state -> Accumulator state state
accum m s t d i =  m i ds where ds = s t d

twostep :: (Shape sh, Elt diff, Elt state) =>  (Acc (Array sh state) -> Acc (Array sh diff)) -> (Acc (Array sh state) -> Acc (Array sh diff)) -> Accumulator diff state -> (Acc (Array sh state) -> Acc (Array sh state)) -> Simulator sh state
twostep predict advance integrate boundry time start = final 
            where
                step t d u = boundry $ Acc.zipWith (\di st -> integrate t di st) d u 
                du_start = predict start
                mid = step (time/2) du_start start 
                du_mid = advance mid
                final = step time du_mid start 