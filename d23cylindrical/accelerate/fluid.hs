module Fluid where 

import Linear

class Num a => Fluid a where 
    pressure :: Floating b => a -> b
    sound :: Floating b => a -> V3 b
    velocity :: Floating b => a -> V3 b 

data Hydro a = Hydro a (a,V3 a,a)

instance Floating a => Fluid (Hydrodynamic a) where
    pressure :: Floating a => (Hydro a) -> a
    pressure (Hydro gamma (density,momentum, energy)) = (gamma - 1)*(energy - quadrance momentum / density / 2.0)
    sound :: Floating a => (Hydro a) -> V3 a
    sound h = let  (Hydro gamma (density,_,_) = h  in gamma * (pressure h) / density
    velocity :: Floating a => (Hydro a) -> V3 a 
    velocity (Hydro _ (density, momentum, _)) = momentum ^/ density 

