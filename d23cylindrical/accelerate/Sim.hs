module Sim where 

import qualified Prelude as P 

import Types 
import Advect
import Limits
import MHD
import Integration
import Geometry 

import Data.Array.Accelerate as Acc
import Data.Array.Accelerate.Linear

projector :: Projector MHD
projector = primativeMHDProject dwlimiter2

differ :: Differ V3 MHD MHD
differ = diff mhdmerge mhdscale

predictorFlux :: Fluxer V3 MHD MHD 
predictorFlux = hancock mhdflux

predictor:: Acc (Array DIM3 (Cell V3)) -> Acc (Array DIM3 MHD) -> Acc (Array DIM3 MHD)
predictor = advection3D projector predictorFlux differ 

simulatorFlux :: Fluxer V3 MHD MHD 
simulatorFlux = average mhdmerge mhdscale mhdflux

simulator:: Acc (Array DIM3 (Cell V3)) -> Acc (Array DIM3 MHD) -> Acc (Array DIM3 MHD)
simulator = advection3D projector simulatorFlux differ 

accumulator :: Accumulator MHD MHD 
accumulator = accum mhdmerge mhdscale

constructSimulator :: (Acc (Array DIM3 MHD) -> Acc (Array DIM3 MHD)) -> Acc (Array DIM3 (Cell V3))-> Simulator DIM3 MHD 
constructSimulator constructboundary cells = twostep (predictor $ cells) (simulator $ cells) accumulator constructboundary 

geometry :: Acc (Array DIM3 (V3 Double,Cell V3))
geometry = generateGeometry (cylindrical3D (0.1,1) (0,1) (0,2*pi)) (constant (Z:.50:.350:.50))

(locations,cells) = unzip geometry :: (Acc (Array DIM3 (V3 Double)),  Acc (Array DIM3 (Cell V3)))

testsimulator ::  Acc (Array DIM3 MHD) -> Acc (Array DIM3 MHD)
testsimulator input = constructSimulator P.id cells (constant 0) input

initial :: Acc (Array DIM3 MHD) 
initial = Acc.map (\_ -> constant mhdzero) locations
