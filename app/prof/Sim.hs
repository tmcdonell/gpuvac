module Sim where 

import qualified Prelude as P 

import VAC.Core.Types as Types
import VAC.Core.Advect as Advect
import VAC.Core.Limits as Limits
import VAC.Physics.MHD as MHD
import VAC.Core.Integration as Integration
import VAC.Core.Geometry as Geometry

import Data.Array.Accelerate as Acc
import Data.Array.Accelerate.Linear

dimensions :: Exp DIM3 
dimensions = constant (Z:.50:.350:.100)

projector :: Projector MHD
projector = primativeMHDProject dwlimiter2

differ :: Differ V3 MHD MHD
differ = diff mhdmerge mhdscale

predictorFlux :: Fluxer MHD MHD 
predictorFlux = hancock mhdflux

predictor:: Acc Geometry3D -> Acc (Array DIM3 MHD) -> Acc (Array DIM3 MHD)
predictor = advection3D projector predictorFlux differ

simulatorFlux :: Fluxer MHD MHD 
simulatorFlux = average mhdmerge mhdscale mhdflux

simulator:: Acc Geometry3D -> Acc (Array DIM3 MHD) -> Acc (Array DIM3 MHD)
simulator = advection3D projector simulatorFlux differ

accumulator :: Accumulator MHD MHD 
accumulator = accum mhdmerge mhdscale

constructSimulator :: (Acc (Array DIM3 MHD) -> Acc (Array DIM3 MHD)) -> Acc Geometry3D -> Simulator DIM3 MHD 
constructSimulator constructboundary geom = twostep (predictor $ geom) (simulator $ geom) accumulator constructboundary 

geometry :: Acc Geometry3D
geometry = cylindrical3D (0.1,1) (0,1) dimensions 

testsimulator ::  Acc Geometry3D -> Acc (Array DIM3 MHD) -> Acc (Array DIM3 MHD)
testsimulator cells input = constructSimulator P.id cells (constant 1e-9) input

initial :: Acc (Array DIM3 MHD) 
initial = Acc.generate dimensions (\_ -> constant mhdzero)
