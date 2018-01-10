
import qualified Prelude as P
import Data.Array.Accelerate          as Acc

import qualified Data.Array.Accelerate.LLVM.PTX     as GPU
import qualified Data.Array.Accelerate.LLVM.Native     as CPU

type Snapshot = Array DIM1 Double
type State = (Snapshot,Snapshot)

range :: Int -> Snapshot
range length = fromList (Z:.length) [0..]

sten :: Stencil3 Double -> Exp Double
sten (a,b,c) = P.sum $ P.zipWith (*) [1.0,-2.0,1.0] [a,b,c] 

boundry :: Boundary Snapshot
boundry = function (\_ -> constant 0.0)

applysten :: Acc Snapshot -> Acc Snapshot
applysten = stencil sten boundry


advance:: Acc State -> Acc State
advance state = 
            lift (next,current)
        where 
            next ::  Acc Snapshot
            next = zipWith3 (\a b c -> a+b+c) diff (map ((*) 2.0) current) (map ( (*) $ negate 1.0) prev)
            diff ::  Acc Snapshot
            diff = applysten current
            current ::  Acc Snapshot
            current = afst state
            prev ::  Acc Snapshot
            prev = asnd state


programRun :: State -> State 
programRun = GPU.run1 advance 

programAdvance :: State -> Int -> State 
programAdvance initial count = P.foldr (\_ state -> programRun state) initial [0..count-1]

start :: State
start = (range 1000, range 1000)

--res = GPU.run (derivative 10000000)
