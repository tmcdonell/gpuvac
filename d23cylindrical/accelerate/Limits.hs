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

scalarProjector :: (Exp Double -> Exp Double -> Exp Double -> Exp Double) -> Exp Double -> Exp Double -> Exp Double -> Exp (Double, Double)
scalarProjector limiter p c n = lift (u,d) 
                        where 
                            s = limiter p c n
                            u = c - 0.5*s 
                            d = c + 0.5*s 

vectorProjector :: (Exp Double -> Exp Double -> Exp Double -> Exp Double) -> Exp (V3 Double) -> Exp (V3 Double) -> Exp (V3 Double) -> Exp (V3 Double,V3 Double)
vectorProjector limiter p c n = lift (P.fmap fst tuple, P.fmap snd tuple) 
                                        where
                                            up :: V3 (Exp Double) 
                                            up = unlift p 
                                            uc :: V3 (Exp Double) 
                                            uc = unlift c 
                                            un :: V3 (Exp Double) 
                                            un = unlift n 
                                            tuple :: V3 (Exp (Double, Double))
                                            tuple = liftA3 (scalarProjector limiter) up uc un

arrayProjection ::(Elt a, Shape sh) => (Exp a -> Exp a -> Exp a -> Exp (a,a)) -> Acc (Array sh a) -> Acc (Array sh a) -> Acc (Array sh a) -> Acc (Array sh a, Array sh a)
arrayProjection projector p c n = lift $ unzip proj where 
                        proj = Acc.zipWith3 projector p c n 


projectScalarArray :: Shape sh => (Exp Double -> Exp Double -> Exp Double-> Exp Double) -> Acc (Array sh Double) -> Acc (Array sh Double) -> Acc (Array sh Double) -> Acc (Array sh Double, Array sh Double)
projectScalarArray limiter p c n = arrayProjection (scalarProjector limiter) p c n 

projectVectorArray :: Shape sh => (Exp Double -> Exp Double -> Exp Double-> Exp Double) -> Acc (Array sh (V3 Double)) -> Acc (Array sh (V3 Double)) -> Acc (Array sh (V3 Double)) -> Acc (Array sh (V3 Double), Array sh (V3 Double))
projectVectorArray limiter p c n = arrayProjection (vectorProjector limiter) p c n 


--
--ktscheme :: (Ord a,Floating a,Advect b) => (b->a) -> b -> b -> b
--ktscheme sc l r = scale (fromRational 0.5) $ fr + fl - scale s (r-l)
--                    where 
--                        s = max (sc l) (sc r)
--                        fr = flux r
--                        fl = flux l
--
--tvdmusclf :: Exp (V3 Double) -> (a,a) -> (b,b) -> (a,a) 
--tvdmusclf dir (uf, df) (us,ds) = (f,f)
--                    where 
--                        f = uf - df + scale (speed (average us ds)) (us - ds)
--