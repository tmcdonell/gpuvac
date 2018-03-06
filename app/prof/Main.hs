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
{-# LANGUAGE TemplateHaskell #-}
import Sim
import VAC.Physics.MHD as MHD
import VAC.Core.Types 

import Data.List as L
import Data.Array.Accelerate as Acc
import Data.Array.Accelerate.LLVM.PTX as PTX 

startcond :: Array DIM3 MHD
startcond = $(PTX.runQ initial)

cell :: Geometry3D 
cell = PTX.runN geometry

step :: Geometry3D -> Array DIM3 MHD -> Array DIM3 MHD
step = $(PTX.runQ testsimulator)

values :: Array DIM3 MHD
values = L.foldl' (\p _ -> step cell p) startcond [1..100]

main :: IO ()
main = putStrLn . show $ indexArray values (Z:.0:.0:.0)
