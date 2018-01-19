module Physics where

import Data.Array.Accelerate as Acc 
import Data.Array.Accelerate.Linear

type PressureField sh = Array sh Double 
type MagneticField sh = Array sh (V3 Double)
type MomentumField sh = Array sh (V3 Double) 
type DensityField sh = Array sh Double
type EnergyField sh = Array sh Double
type ForceField sh = Array sh (V3 Double) 

class Energy a where 
    energy :: Exp a -> Exp Double

instance Energy (Exp Magnet) where
    energy (Magnet v) = (quadrance v) / 2.0 

magneticenergy :: Shape sh => Acc (MagneticField sh) -> Acc (EnergyField sh)
magneticenergy mag = Acc.map (\b -> (quadrance b) / 2 ) mag

kineticenergy :: Shape sh => Acc (MomentumField sh) -> Acc (DensityField sh) -> Acc (EnergyField sh)
kineticenergy mom den = Acc.zipWith (\m d -> (quadrance m) / 2.0 / d) mom den 

