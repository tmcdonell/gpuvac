{-# LANGUAGE Rank2Types #-}
module Advect where 

import Control.Lens 
import Linear

class Num a => Advect a where 
    --given a cell state, compute the boundry flux (need to include direction)
    flux :: a -> a
    --given a cell state, compute the spectral radius (max eignevalue) of df(u)/du
    --typically maximum speed along direction 
    speed :: Floating b => (Lens' (V3 b) b) -> a -> b
    primative :: a -> a 
    conserve :: a -> a
    scale :: Floating b => b -> a -> a 

extend :: (Advect a) => ((a,a,a)->a) -> (a,a,a) -> (a,a)
extend limiter (p,c,n) = (conserve lp,conserve rp) 
                        where
                            slope = limiter (primative p,primative c,primative n)
                            step = scale (fromRational 0.5) slope
                            lp = c - step 
                            rp = c + step


hancock :: (Floating b, Advect a) => (Lens' (V3 b) b) -> (a,a) -> (a,a) -> (a,a) 
hancock _ (_,upwindr) (dwnwindl,_) =  (flux upwindr,flux dwnwindl) 

ktscheme :: (Ord a,Floating a,Advect b) => (b->a) -> b -> b -> b
ktscheme sc l r = scale (fromRational 0.5) $ fr + fl - scale s (r-l)
                    where 
                        s = max (sc l) (sc r)
                        fr = flux r
                        fl = flux l

tvdmusclf :: (Ord b,Floating b, Advect a) => (Lens' (V3 b) b) -> (a,a) -> (a,a) -> (a,a) 
tvdmusclf dir (ul, ur) (dl,dr) = (uf,df)
    where 
          sc = speed dir 
          uf = ktscheme sc ul ur 
          df = ktscheme sc dl dr