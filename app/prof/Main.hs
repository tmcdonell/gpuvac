{-# LANGUAGE TemplateHaskell #-}
import Sim
import VAC.Physics.MHD as MHD

import Data.List as L
import Data.Array.Accelerate as Acc
import Data.Array.Accelerate.LLVM.PTX as PTX 

startcond :: Array DIM3 MHD
startcond = PTX.runN initial

step :: Array DIM3 MHD -> Array DIM3 MHD
step = PTX.run1 testsimulator

values :: Array DIM3 MHD
values = L.foldl' (\p _ -> step p) startcond [1..100]

main :: IO ()
main = putStrLn . show $ indexArray values (Z:.0:.0:.0)
