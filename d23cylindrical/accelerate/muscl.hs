{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE ScopedTypeVariables #-}
module MUSCL where

import qualified Data.Array.Accelerate.LLVM.PTX     as GPU
import Data.Array.Accelerate as A
import qualified Prelude as P

minmod :: (Ord a, P.Num (Exp a)) => Exp a -> Exp a
minmod r = max  0 $ min r 1

ratio :: P.Fractional a => (a,a,a) -> a
ratio (p,c,n) =  (c - p) / (n - c)

correctflux :: P.Fractional a => (t -> a) -> (t -> t -> a) -> t -> t -> a
correctflux toflx transport left right = 0.5 * (toflx right + toflx left - transport left right)

totransport :: P.Num a => (a -> a -> a) -> a -> a -> a
totransport cmpt left right = (cmpt left right) * (right - left)

step :: P.Fractional t => t -> (t -> t) -> (t -> t) -> (t -> t -> t) -> (t, t, t, t, t) -> t
step dt limiter toflx cmpt (pp,p,c,n,nn) = c + dt * (upwind_flux - downwind_flux) 
    where 
        limit_prev = limiter $ ratio (pp,p,c)
        limit_curr = limiter $ ratio (p,c,n)
        limit_next = limiter $ ratio (c,n,nn)
        downwind_left = c + 0.5 * limit_curr * (c - p) 
        downwind_right = n - 0.5 * limit_next * (n - c)
        upwind_left = p + 0.5 * limit_prev * (c - p) 
        upwind_right = c - 0.5 * limit_curr * (n - c)
        downwind_flux = correctflux toflx (totransport cmpt) downwind_left downwind_right
        upwind_flux = correctflux toflx (totransport cmpt) upwind_left upwind_right

simple1d :: (P.Fractional (Exp a), Ord a) => (Exp a, Exp a, Exp a, Exp a, Exp a) -> Exp a
simple1d sten = step 0.1 minmod (\x->x) (\x y -> 0.5*x + 0.5*y) sten

advance:: (Stencil t a (Exp b, Exp b, Exp b, Exp b, Exp b),P.Fractional (Exp b), Ord b) => Acc (Array t a) -> Acc (Array t b)
advance input = stencil simple1d clamp input

input :: Array DIM1 Double
input = fromList (Z:.100) [0..]

res :: Array DIM1 Double
res = GPU.run1 advance input
