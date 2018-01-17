{-# LANGUAGE Rank2Types #-}
module Advect where 

import Control.Lens 
import Linear

class Floating a => Advect a where 
    --given a cell state, compute the boundry flux (need to include direction)
    flux :: a -> a
    --given a cell state, compute the spectral radius (max eignevalue) of df(u)/du
    --typically maximum speed along direction 
    speed :: Floating b => (Lens' (V3 a) b) -> a -> b
    primative :: a -> a 
    conserve :: a -> a
    scale :: Floating b => b -> a -> a 

extend :: (Advect a) => ((a,a,a)->a) -> (a,a,a) -> (a,a)
extend limiter (p,c,n) = (conserve lp,conserve rp) 
                        where
                            slope = limiter (primative p,primative c,primative n)
                            lp = c - (fromRational 0.5) * slope
                            rp = c + (fromRational 0.5) * slope


hancock :: (Advect a,Floating b) => ((a,a,a) -> a) -> (b,b,b) -> (a,a,a) -> a
hancock limiter (lA,volume,rA) state = scale (recip volume) $ (scale lA (flux l)) - (scale rA (flux r))
        where
            (_,c,_) = state  
            (l,r) = extend limiter state
            