module MUSCL where

import qualified Prelude as P 
import Data.Array.Accelerate as A

step limiter toflx cmpt (pp,p,c,n,nn) = c 
    where 
        limit_prev = limiter $ (p - pp) / (c - p) 
        limit_curr = limiter $ (c - p) / (n-c)
        limit_next = limiter $ (n - c) / (nn - n) 
        downwind_left = c + 0.5 * limit_curr * (c - p) 
        downwind_right = n - 0.5 * limit_next * (n - c)
        upwind_left = p + 0.5 * limit_prev * (c - p) 
        upwind_right = c - 0.5 * limit_curr * (n - c)
        downwind_flux = 0.5 * (toflx downwind_right + toflx downwind_left - (cmpt downwind_left downwind_right) * (downwind_right - downwind_left))
        upwind_flux = 0.5 * (toflx upwind_right + toflx upwind_left - (cmpt upwind_left upwind_right) * (upwind_right - upwind_left))


