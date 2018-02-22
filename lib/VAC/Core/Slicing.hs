{-# LANGUAGE RankNTypes #-}
module VAC.Core.Slicing where

import Data.Array.Accelerate as A
import Data.Array.Accelerate.Control.Lens

slitOn :: (Shape sh, Elt e)
    => Lens' (Exp sh) (Exp Int)
    -> Exp Int
    -> Exp Int
    -> Acc (Array sh e)
    -> Acc (Array sh e)
slitOn dim m n acc =
  let
      m'  = the (unit m)
      n'  = the (unit n)
      sh' = over dim (\i -> n' `A.min` ((i-m') `A.max` 0)) (shape acc)
  in
  backpermute sh' (& dim +~ m') acc

trim :: (Shape sh, Elt e) 
    => Lens' (Exp sh) (Exp Int)
    -> Exp Int
    -> Acc (Array sh e) 
    -> Acc (Array sh e) 
trim dim n acc = 
    let 
        m = (shape acc) ^.dim
    in 
    slitOn dim n (m-n) acc

dim3Array :: Acc (Array DIM3 Float)
dim3Array = enumFromN (constant (Z:.10:.10:.10)) 0

test :: Acc (Array DIM3 Float)
test = trim _3 (constant 2) dim3Array

sten3 :: (Shape sh, Elt e)
    => Lens' (Exp sh) (Exp Int)
    -> Acc (Array sh e) 
    -> Acc (Array sh e, Array sh e, Array sh e) 
sten3 dim acc = 
    let 
        n = (shape acc)^.dim
    in
    lift (slitOn dim 0 (n-2) acc, slitOn dim 1 (n-1) acc, slitOn dim 2 n acc)

left :: (Shape sh, Elt e)
    => Lens' (Exp sh) (Exp Int)
    -> Acc (Array sh e) 
    -> Acc (Array sh e)
left dim acc = 
    let 
        n = (shape acc)^.dim
    in
    slitOn dim 0 (n-1) acc


right :: (Shape sh, Elt e)
    => Lens' (Exp sh) (Exp Int)
    -> Acc (Array sh e) 
    -> Acc (Array sh e)
right dim acc = 
    let 
        n = (shape acc)^.dim
    in
    slitOn dim 1 n acc

sten2 :: (Shape sh, Elt e)
    => Lens' (Exp sh) (Exp Int)
    -> Acc (Array sh e) 
    -> Acc (Array sh e, Array sh e) 
sten2 dim acc = 
    let 
        m = (shape acc)^.dim
        n= m-1
    in
    lift (slitOn dim 0 (n-1) acc, slitOn dim 1 n acc)