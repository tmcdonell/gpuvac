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
{-# LANGUAGE Rank2Types #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE FunctionalDependencies #-}
module VAC.Core.Integration where

import qualified Prelude as P 
import Data.Array.Accelerate as Acc

import VAC.Core.Types as Types

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