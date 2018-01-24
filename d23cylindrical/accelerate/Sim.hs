module Sim where 

import Types 
import Advect
import Limits
import MHD
import Integration

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

constructSimulator :: Acc (Array DIM3 (Cell V3)) -> (Acc (Array DIM3 MHD) -> Acc (Array DIM3 MHD)) -> Simulator DIM3 MHD 
constructSimulator cells constructboundary = twostep (predictor $ cells) (simulator $ cells) accumulator constructboundary 
