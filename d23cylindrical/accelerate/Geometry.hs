module Geometry where 

import qualified Prelude as P 

import Data.Array.Accelerate as Acc
import Data.Array.Accelerate.Linear

type Patch order  = order (order Double, order Double) 
type Geom order  = (Patch order,Double)

cartesian3D :: (Double,Double) -> (Double,Double) -> (Double,Double) -> Exp DIM3 -> Exp DIM3 -> Exp (V3 Double,Geom V3)
cartesian3D (xmin,xmax) (ymin,ymax) (zmin,zmax) dimensions index = lift (position, geom)
                                    where
                                        (xpos,ypos,zpos) = unlift $ unindex3 index :: (Exp Int,Exp Int,Exp Int) 
                                        (xsize,ysize,zsize) = unlift $ unindex3 dimensions :: (Exp Int,Exp Int,Exp Int) 
                                        dx = constant $ xmax - xmin :: Exp Double 
                                        dy = constant $ ymax - ymin :: Exp Double 
                                        dz = constant $ zmax - zmin :: Exp Double 
                                        dV = dx*dy*dz :: Exp Double 
                                        dAx = lift $ V3 (dy*dz) 0 0 :: Exp (V3 Double) 
                                        dAy = lift $ V3 0 (dx*dz) 0 :: Exp (V3 Double) 
                                        dAz = lift $ V3 0 0 (dx*dy) :: Exp (V3 Double) 
                                        patch = lift $ V3 (dAx,dAx) (dAy,dAy) (dAz,dAz) :: Exp (V3 (V3 Double, V3 Double))
                                        geom = lift $ (patch,dV) :: Exp (Geom V3) 
                                        x = dx * (fromIntegral xpos / fromIntegral xsize) + dx/2.0 + constant xmin :: Exp Double 
                                        y = dy * (fromIntegral ypos / fromIntegral ysize) + dy/2.0 + constant ymin :: Exp Double
                                        z = dz * (fromIntegral zpos / fromIntegral zsize) + dz/2.0 + constant zmin:: Exp Double 
                                        position = lift $ V3 x y z :: Exp (V3 Double) 


cylindrical3D :: (Double,Double) -> (Double,Double) -> (Double,Double) -> Exp DIM3 -> Exp DIM3 -> Exp (V3 Double,Geom V3)
cylindrical3D (rmin,rmax) (zmin,zmax) (tmin,tmax) dimensions index = lift (position,geom)  
                                        where 
                                            (rpos,zpos,tpos) = unlift $ unindex3 index :: (Exp Int,Exp Int,Exp Int) 
                                            (rsize,zsize,tsize) = unlift $ unindex3 dimensions :: (Exp Int,Exp Int,Exp Int) 
                                            dr = constant $ (rmax - rmin) :: Exp Double 
                                            dz = constant $ (zmax - zmin) :: Exp Double 
                                            dt = constant $ (tmax - tmin) :: Exp Double 
                                            rs = dr * (fromIntegral rpos / fromIntegral rsize) + constant rmin :: Exp Double 
                                            r = rs + dr / 2.0 :: Exp Double 
                                            re = rs + dr :: Exp Double 
                                            z = dz * (fromIntegral zpos / fromIntegral zsize) + dz/2.0 + constant zmin:: Exp Double
                                            ts = dt * (fromIntegral tpos / fromIntegral tsize) + constant tmin:: Exp Double 
                                            t = ts + dt/2.0 :: Exp Double 
                                            te = ts + dt :: Exp Double 
                                            dAz' = dr * r * dt :: Exp Double
                                            dV = dAz'*dz :: Exp Double
                                            dAr' rp = V3 ((cos t) * (dz*dt*rp)) ((sin t) * (dz*dt*rp)) 0.0 :: V3 (Exp Double)
                                            dAru = lift $ dAr' rs  :: Exp (V3 Double)
                                            dArd = lift $ dAr' re :: Exp (V3 Double)
                                            dAt' tp = lift $ V3 (-1*(sin tp)*dr*dz) ((cos tp)*dr*dz) 0.0 :: Exp (V3 Double)
                                            dAtu = dAt' ts :: Exp (V3 Double)
                                            dAtd = dAt' te :: Exp (V3 Double)
                                            dAz = lift $ V3 0 0 (negate dAz') :: Exp (V3 Double) 
                                            patch = lift $ V3 (dAru,dArd) (dAz,dAz) (dAtu,dAtd) :: Exp (V3 (V3 Double, V3 Double))
                                            geom = lift $ (patch,dV) :: Exp (Geom V3) 
                                            position = lift $ V3 r z t :: Exp (V3 Double) 


generateGeometry :: (Elt a, Shape sh) => (Exp sh -> Exp sh -> Exp a) -> Exp sh -> Acc (Array sh a) 
generateGeometry gener size = generate size (gener size) 