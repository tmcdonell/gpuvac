{-# LANGUAGE ScopedTypeVariables #-}

import qualified Data.List as L 
import qualified Prelude as P

import qualified Foreign.CUDA.Driver                    as CUDA
import Data.Array.Accelerate          as A

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

for2 :: Arrays a => (Acc a -> Acc a) -> Acc a -> Acc a
for2 f i = (f.f) i 

for4 :: Arrays a => (Acc a -> Acc a) -> Acc a -> Acc a
for4 f i = (for2 (for2 f)) i 

for8 :: Arrays a => (Acc a -> Acc a) -> Acc a -> Acc a
for8 f i = (for4 (for2 f)) i 

for10 :: Arrays a => (Acc a -> Acc a) -> Acc a -> Acc a
for10 f i = (for8 f .for2 f) i 

for100 :: Arrays a => (Acc a -> Acc a) -> Acc a -> Acc a
for100 f i = (for10 (for10 f)) i 

for10000 :: Arrays a => (Acc a -> Acc a) -> Acc a -> Acc a
for10000 f i = (for100 (for100 f)) i 

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
programRun = GPU.run1 $ for10000 advance 

programAdvance :: State -> Int -> State 
programAdvance initial count = L.foldl' (\state _ -> programRun state) initial [0..count-1]

-- need to do multiple iterations on GPU
start :: State
start = (range 1000, range 1000)

main :: P.IO ()
main = do
  CUDA.initialise []
  dev <- CUDA.device 0
  ctx <- CUDA.create dev [CUDA.SchedYield]  -- http://tmcdonell-bot.github.io/accelerate-travis-buildbot/cuda-0.9.0.0/Foreign-CUDA-Driver-Context-Base.html#t:ContextFlag
  ptx <- GPU.createTargetFromContext ctx
  let r = GPU.run1With ptx (for10000 advance)
  let output = P.show $ r start
  P.putStrLn output  

