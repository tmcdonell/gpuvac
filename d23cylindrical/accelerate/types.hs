{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE ScopedTypeVariables #-}
module Types where

import Data.Array.Accelerate
import Data.Array.Accelerate.Linear

type R = Double 

type Momentum = V3 R
type Density = R 

type Fluid = (Momentum, Density) 

type Energy = R 
type Magnetic = V3 R

type SingleFluidMHD = (Fluid,Energy,Magnetic)
