{-# LANGUAGE BangPatterns    #-}
{-# LANGUAGE TemplateHaskell #-}
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

import Sim
import VAC.Physics.MHD                                    as MHD
import VAC.Core.Types

import Data.Array.Accelerate                              as A
import Data.Array.Accelerate.Linear
import Data.Array.Accelerate.LLVM.PTX                     as PTX
-- import Data.Array.Accelerate.LLVM.Native                  as CPU

import Prelude                                            as P


startcond :: Array DIM3 MHD
startcond = $(PTX.runQ initial)
-- startcond = CPU.runN initial

cell :: Array DIM3 (Cell V3)
cell = $(PTX.runQ voxels)
-- cell = CPU.runN voxels

step :: Array DIM3 (Cell V3) -> Array DIM3 MHD -> Array DIM3 MHD
step = $(PTX.runQ testsimulator)
-- step = CPU.runN testsimulator

values :: Array DIM3 MHD
values = loop 0 startcond
  where
    lIMIT = 100

    loop :: Int -> Array DIM3 MHD -> Array DIM3 MHD
    loop !i !p
      | i P.< lIMIT = loop (i+1) (step cell p)
      | otherwise   = p

main :: IO ()
main = putStrLn . show $ indexArray values (Z:.0:.0:.0)

