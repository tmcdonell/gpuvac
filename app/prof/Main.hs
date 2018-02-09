{-# LANGUAGE TemplateHaskell #-}
import Sim
import VAC.Physics.MHD as MHD
import VAC.Core.Types 

import Data.List as L
import Data.Array.Accelerate as Acc
import Data.Array.Accelerate.Linear
import Data.Array.Accelerate.LLVM.PTX as PTX 

startcond :: Array DIM3 MHD
startcond = $(PTX.runQ initial)

cell :: Array DIM3 (Cell V3)
cell = $(PTX.runQ voxels)

step :: Array DIM3 (Cell V3) -> Array DIM3 MHD -> Array DIM3 MHD
step = $(PTX.runQ testsimulator)

values :: Array DIM3 MHD
values = L.foldl' (\p _ -> step cell p) startcond [1..100]

main :: IO ()
main = putStrLn . show $ indexArray values (Z:.0:.0:.0)
