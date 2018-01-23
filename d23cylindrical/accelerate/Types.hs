module Types where 

import Data.Array.Accelerate as Acc
import Data.Array.Accelerate.Linear

type Limiter = Exp Double -> Exp Double -> Exp Double -> Exp Double 

-- type Patch is a vector of vector pairs, the first vector encodes which coordinate 
-- frame we are travelling along, while the inner tuple contains the vectors encoding
-- the upstream and downstream surface normals
type Patch order  = order (order Double, order Double) 

-- A Cell contains a Patch which has all the appropriate area details as well as a 
-- volume term. 
type Cell order  = (Patch order,Double)

-- A projector allows us to project a cell state to the cell boundaries. It takes
-- the previous cell state, the current cell and the next cell along one of the 
-- coordinate dimensions. It uses this information to compute a local slope and projection
-- this function is where any manner of TVD or slope limiting is injected
type Projector a = Exp a -> Exp a -> Exp a -> Exp (a,a)

-- A fluxer takes a quantity dA which is a vector with a magnitude equal to the cell boundry area
-- and the direction normal to the surface in the direction of the downstream cell 
-- it also takes the quanties of the upstream and downstream state at the cell boundry.
-- The result of this computation is the upstream and downstream fluxes at the boundry 
-- which should be the same for a conservative method. The advection method, such as 
-- hancock or kt would be encoded within this function
type Fluxer vec state flux = Exp (vec Double) -> Exp (state,state) -> Exp (flux,flux) 

-- A Differ is responsible for taking the cell volume, as well as a vector of 
-- flux pairs and reducing that to a time derivative. 
type Differ vec flux diff = Exp Double -> Exp (vec (flux,flux)) -> Exp diff 

-- An Accumulator is responsible for consuming a time derivative and a time delta
-- it takes the change in time, the derivative in time and the initial state
-- the result should be essentially final = intial + dt * du_dt. However, this allows 
-- for time derivative to be of a different shape and type to your state variable.
type Accumulator diff state = Exp Double -> Exp diff -> Exp state -> Exp state