module Limits where

--slope limited derivative
dwlimiter2 (p,c,n) = 
        2.0 * s * maximum 
            where 
                s = signum (n-c)
                maximum = max 0 minimum 
                minimum = min absolute $ min previous average
                absolute = abs (n-c)
                previous = s * (c-p) 
                average = 0.25 * s * (n-p) 

minmod :: (Ord a,Num a) => a -> a
minmod r = max  0 $ min r 1

ratio :: Fractional a => (a,a,a) -> a
ratio (p,c,n) =  (c - p) / (n - c)

slope :: (Fractional a, Ord a, Num a) => (a -> a) -> (a,a,a) -> a
slope limiter stencil = (limiter $ ratio stencil) * (n-c) 
                            where (_,c,n) = stencil