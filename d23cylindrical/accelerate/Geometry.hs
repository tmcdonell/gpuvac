module Geometry where 

import qualified Prelude as P 

import Data.Array.Accelerate as Acc
import Data.Array.Accelerate.Linear
import Types 


cartesian3D :: (Real,Real) -> (Real,Real) -> (Real,Real) -> Exp DIM3 -> Exp DIM3 -> Exp (V3 Real,Cell V3)
cartesian3D (xmin,xmax) (ymin,ymax) (zmin,zmax) dimensions index = lift (position, geom)
                                    where
                                        (xpos,ypos,zpos) = unlift $ unindex3 index :: (Exp Int,Exp Int,Exp Int) 
                                        (xsize,ysize,zsize) = unlift $ unindex3 dimensions :: (Exp Int,Exp Int,Exp Int) 
                                        dx = constant $ xmax - xmin :: Exp Real 
                                        dy = constant $ ymax - ymin :: Exp Real 
                                        dz = constant $ zmax - zmin :: Exp Real 
                                        dV = dx*dy*dz :: Exp Real 
                                        dAx = lift $ V3 (dy*dz) 0 0 :: Exp (V3 Real) 
                                        dAy = lift $ V3 0 (dx*dz) 0 :: Exp (V3 Real) 
                                        dAz = lift $ V3 0 0 (dx*dy) :: Exp (V3 Real) 
                                        patch = lift $ V3 (dAx,dAx) (dAy,dAy) (dAz,dAz) :: Exp (V3 (V3 Real, V3 Real))
                                        geom = lift $ (patch,dV) :: Exp (Cell V3) 
                                        x = dx * (fromIntegral xpos / fromIntegral xsize) + dx/2.0 + constant xmin :: Exp Real 
                                        y = dy * (fromIntegral ypos / fromIntegral ysize) + dy/2.0 + constant ymin :: Exp Real
                                        z = dz * (fromIntegral zpos / fromIntegral zsize) + dz/2.0 + constant zmin:: Exp Real 
                                        position = lift $ V3 x y z :: Exp (V3 Real) 


cylindrical3D :: (Real,Real) -> (Real,Real) -> (Real,Real) -> Exp DIM3 -> Exp DIM3 -> Exp (V3 Real,Cell V3)
cylindrical3D (rmin,rmax) (zmin,zmax) (tmin,tmax) dimensions index = lift (position,geom)  
                                        where 
                                            (rpos,zpos,tpos) = unlift $ unindex3 index :: (Exp Int,Exp Int,Exp Int) 
                                            (rsize,zsize,tsize) = unlift $ unindex3 dimensions :: (Exp Int,Exp Int,Exp Int) 
                                            dr = constant $ (rmax - rmin) :: Exp Real 
                                            dz = constant $ (zmax - zmin) :: Exp Real 
                                            dt = constant $ (tmax - tmin) :: Exp Real 
                                            rs = dr * (fromIntegral rpos / fromIntegral rsize) + constant rmin :: Exp Real 
                                            r = rs + dr / 2.0 :: Exp Real 
                                            re = rs + dr :: Exp Real 
                                            z = dz * (fromIntegral zpos / fromIntegral zsize) + dz/2.0 + constant zmin:: Exp Real
                                            ts = dt * (fromIntegral tpos / fromIntegral tsize) + constant tmin:: Exp Real 
                                            t = ts + dt/2.0 :: Exp Real 
                                            te = ts + dt :: Exp Real 
                                            dAz' = dr * r * dt :: Exp Real
                                            dV = dAz'*dz :: Exp Real
                                            dAr' rp = V3 ((cos t) * (dz*dt*rp)) ((sin t) * (dz*dt*rp)) 0.0 :: V3 (Exp Real)
                                            dAru = lift $ dAr' rs  :: Exp (V3 Real)
                                            dArd = lift $ dAr' re :: Exp (V3 Real)
                                            dAt' tp = lift $ V3 (-1*(sin tp)*dr*dz) ((cos tp)*dr*dz) 0.0 :: Exp (V3 Real)
                                            dAtu = dAt' ts :: Exp (V3 Real)
                                            dAtd = dAt' te :: Exp (V3 Real)
                                            dAz = lift $ V3 0 0 (negate dAz') :: Exp (V3 Real) 
                                            patch = lift $ V3 (dAru,dArd) (dAz,dAz) (dAtu,dAtd) :: Exp (V3 (V3 Real, V3 Real))
                                            geom = lift $ (patch,dV) :: Exp (Cell V3) 
                                            position = lift $ V3 r z t :: Exp (V3 Real) 


generateGeometry :: (Elt a, Shape sh) => (Exp sh -> Exp sh -> Exp a) -> Exp sh -> Acc (Array sh a) 
generateGeometry gener size = generate size (gener size) 