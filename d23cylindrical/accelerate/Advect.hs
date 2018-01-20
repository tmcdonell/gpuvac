{-# LANGUAGE Rank2Types #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE FunctionalDependencies #-}

module Advect where 

import qualified Prelude as P 

import Control.Lens
import Data.Array.Accelerate as Acc
import Data.Array.Accelerate.Linear

class Interp state where 
    interp :: Exp state -> Exp state -> Exp state -> Exp (state,state) 

class Flux state flux | state -> flux where 
    flux :: Exp (V3 Double) -> Exp (state,state) -> Exp (flux,flux) 

class Accumulate flux diff | diff -> flux where 
    accumulate :: Exp Double -> Exp flux -> Exp diff 

class Merge diff  where 
    merge :: Exp diff -> Exp diff -> Exp diff  

class Integrate diff state | state -> diff where 
    integrate :: Exp Double -> Exp diff -> Exp state -> Exp state 

type Geom = (V3 (V3 Double, V3 Double),Double) 

--stencil1D (pp,p,c,n,nn) = lift (pp,p,c,n,nn) 

interpStep :: (Interp state,Shape sh,Elt state) => Acc (Array sh state) ->Acc (Array sh state) ->Acc (Array sh state) -> (Acc (Array sh state), Acc (Array sh state))
interpStep p c n = unzip $ Acc.zipWith3 interp p c n 


--core of the advection flux calculation algorithm
advect :: (Shape sh, Interp state, Flux state fl,Elt state,Elt fl) => Acc (Array sh (V3 Double, V3 Double)) -> Acc (Array sh (state,state,state,state,state)) -> Acc (Array sh (fl,fl)) 
advect areas state = zip uflux dflux where 
                    (pp,p,c,n,nn) = unzip5 state
                    (_,dp) = interpStep pp p c
                    (uc,dc) = interpStep p c n
                    (un,_) = interpStep c n nn
                    (unorm,dnorm) = unzip areas
                    ustate = zip dp uc
                    dstate = zip dc un
                    (_,uflux) = unzip $ Acc.zipWith flux unorm ustate
                    (dflux,_) = unzip $ Acc.zipWith flux dnorm dstate


-- the TVD method looks at the flux and state on either side of an interface
-- and then determines the downstream flux
-- for hancock we just use the upstream and downstream flux as the upstream and 
-- downstream flux respectively
hancock :: forall a b. Exp (V3 Double) -> (a,a) -> (b,b) -> (a,a)
hancock _ (uf,df) _ = (uf,df) 



--compute the state derivative given a cell geometry and the fluxes
--accumulate :: (Floating a,Advect b) => Geometry a -> VoxelFlux b -> State b
--accumulate geom flux = scale (recip volume) totalflux
--        where 
--            (faceareas,volume) = geom
--            dflux = liftA2 fluxsum faceareas flux 
--            totalflux = sum dflux

--two step performs a basic two step time integration
--it takes a geometry over which 
--twostep :: (Floating a,Advect b) => (State b -> VoxelFlux b) -> (State b -> VoxelFlux b) -> Geometry a -> Time a -> State b -> State b
--twostep predict main geom dt start = end  
--        where 
--        du_start = accumulate geom (predict start)
--        mid = start + scale (dt/2) start
--        du_mid = accumulate geom (main mid) 
--        end = start + scale dt du_mid

