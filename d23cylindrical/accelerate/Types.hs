module Types where 

import Data.Array.Accelerate as Acc
import Data.Array.Accelerate.Linear

type Projector a = Exp a -> Exp a -> Exp a -> Exp (a,a)

-- A fluxer takes a quantity dA which is a vector with a magnitude equal to the cell boundry area
-- and the direction normal to the surface in the direction of the downstream cell 
-- it also takes the quanties of the left and right state at the cell boundry.
-- The result of this computation is the left and right fluxes which will be the 
-- same for a conservative method.
type Fluxer vec state flux = Exp (vec Double) -> Exp (state,state) -> Exp (flux,flux) 

type Patch order  = order (order Double, order Double) 
type Geom order  = (Patch order,Double)

type Differ vec flux diff = Exp Double -> Exp (vec (flux,flux)) -> Exp diff 

type Accumulator diff state = Exp Double -> Exp diff -> Exp state -> Exp state