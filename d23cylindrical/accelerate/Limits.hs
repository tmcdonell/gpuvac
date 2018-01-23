{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE ConstraintKinds #-}
{-# LANGUAGE FlexibleContexts #-}
module Limits where

import Control.Applicative
import qualified Prelude as P
import Data.Array.Accelerate as Acc 
import Data.Array.Accelerate.Linear

import Types

--slope limited derivative
dwlimiter2 :: Limiter
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

vectorlimit::Limiter -> Exp (V3 Double) -> Exp (V3 Double) -> Exp (V3 Double) -> Exp (V3 Double)
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

projectScalar :: Limiter -> Exp Double -> Exp Double -> Exp Double -> Exp (Double, Double)
projectScalar limiter p c n = lift (u,d) 
                        where 
                            s = limiter p c n
                            u = c - 0.5*s 
                            d = c + 0.5*s 

projectVector :: Limiter -> Exp (V3 Double) -> Exp (V3 Double) -> Exp (V3 Double) -> Exp (V3 Double,V3 Double)
projectVector limiter p c n = lift (P.fmap fst tuple, P.fmap snd tuple) 
                                        where
                                            up :: V3 (Exp Double) 
                                            up = unlift p 
                                            uc :: V3 (Exp Double) 
                                            uc = unlift c 
                                            un :: V3 (Exp Double) 
                                            un = unlift n 
                                            tuple :: V3 (Exp (Double, Double))
                                            tuple = liftA3 (projectScalar limiter) up uc un



