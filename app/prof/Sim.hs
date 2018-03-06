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
testsimulator cells input = constructSimulator P.id cells (constant 0.0) input

initial :: Acc (Array DIM3 MHD) 
initial = Acc.generate dimensions (\_ -> constant mhdzero)
