{-# LANGUAGE TemplateHaskell #-}
import Sim
import VAC.Physics.MHD as MHD

import Data.List as L
import Data.Array.Accelerate as Acc
import Data.Array.Accelerate.Linear 
import Data.Array.Accelerate.LLVM.PTX as PTX 

startcond = PTX.runN initial

step = PTX.run1 testsimulator

values = L.foldl' (\p _ -> step p) startcond [1..1000]

main = putStrLn . show $ indexArray values (Z:.0:.0:.0)
