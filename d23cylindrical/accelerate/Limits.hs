module Limits where

import Control.Applicative
import qualified Prelude as P
import Data.Array.Accelerate as Acc 
import Data.Array.Accelerate.Linear

--slope limited derivative
dwlimiter2 :: Exp Double -> Exp Double -> Exp Double -> Exp Double
dwlimiter2 p c n = 
    2.0 * s * maximum 
        where 
            nc = n-c
            np = n-p
            cp = c-p
            s = signum nc
            maximum = max 0 minimum 
            minimum = min absolute $ min previous average
            absolute = abs nc
            previous = s * cp 
            average = 0.25 * s * np 

minmod :: Exp Double -> Exp Double
minmod r = max  0 $ min r 1

ratio :: Exp Double -> Exp Double -> Exp Double -> Exp Double
ratio p c n =  (c - p) / (n - c)

slope :: (Exp Double -> Exp Double) -> Exp Double -> Exp Double -> Exp Double -> Exp Double
slope limiter p c n = (limiter $ ratio p c n) * (n-c) 

vectorlimit::(Exp Double -> Exp Double-> Exp Double ->Exp Double) -> Exp (V3 Double) -> Exp (V3 Double) -> Exp (V3 Double) -> Exp (V3 Double)
vectorlimit limiter p c n = lift tuple
                            where
                                up :: V3 (Exp Double) 
                                up = unlift p 
                                uc :: V3 (Exp Double) 
                                uc = unlift c 
                                un :: V3 (Exp Double) 
                                un = unlift n 
                                tuple :: V3 (Exp Double)
                                tuple = liftA3 limiter up uc un

--projectScalar :: Shape sh => (Exp Double,Exp Double,Exp Double)->Exp Double) -> Acc (Array sh Double, Array sh Double, Array sh Double) -> Acc (Array sh Double, Array sh Double)
--projectScalar limiter pcn = lift (u,d) 
--                                where 
--                                    (p,c,n) = unlift pcn 
--                                    Acc.zipWith3 (,,) 


