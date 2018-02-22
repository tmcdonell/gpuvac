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
values = L.foldl' (\p _ -> step cell p) startcond [1..1]

main :: IO ()
main = putStrLn . show $ indexArray values (Z:.0:.0:.0)
