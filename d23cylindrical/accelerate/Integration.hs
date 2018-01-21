{-# LANGUAGE Rank2Types #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE FunctionalDependencies #-}
module Integration where

import qualified Prelude as P 
import Data.Array.Accelerate as Acc


class (Elt diff, Elt state) => Integrate diff state | state -> diff where 
    integrate :: Exp Double -> Exp diff -> Exp state -> Exp state 

twostep :: (Shape sh, Integrate diff state) => (Acc (Array sh state) -> Acc (Array sh diff)) -> (Acc (Array sh state) -> Acc (Array sh diff)) -> Exp Double -> Acc (Array sh state) -> Acc (Array sh state) 
twostep predict advance time start = final 
            where
                step t d u = Acc.zipWith (\di st -> integrate t di st) d u 
                du_start = predict start
                mid = step (time/2) du_start start  
                du_mid = advance mid
                final = step time du_mid start 