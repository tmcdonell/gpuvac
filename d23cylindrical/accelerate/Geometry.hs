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
                                        x = dx * (fromIntegral xpos / fromIntegral xsize) + dx/2.0 :: Exp Double 
                                        y = dy * (fromIntegral ypos / fromIntegral ysize) + dy/2.0 :: Exp Double
                                        z = dz * (fromIntegral zpos / fromIntegral zsize) + dz/2.0 :: Exp Double 
                                        position = lift $ V3 x y z :: Exp (V3 Double) 


cylindrical3D :: (Double,Double) -> (Double,Double) -> (Double,Double) -> Exp sh -> Exp (V3 Double,Geom V3)
cylindrical3D (rmin,rmax) (zmin,zmax) (tmin,tmax) dimensions index = lift (position,geom)  
                                        where 
                                            (rpos,zpos,tpos) = unlift $ unindex3 index :: (Exp Int,Exp Int,Exp Int) 
                                            (rsize,zsize,tsize) = unlift $ unindex3 dimensions :: (Exp Int,Exp Int,Exp Int) 
                                            dr = constant $ rmax - rmin :: Exp Double 
                                            dz = constant $ zmax - zmin :: Exp Double 
                                            dt = constant $ tmax - tmin :: Exp Double 
                                            rmin = dr * (fromIntegral rpos / fromIntegral rsize) :: Exp Double 
                                            r = r + dr / 2.0 :: Exp Double 
                                            rmax = r + dr :: Exp Double 
                                            z = dz * (fromIntegral zpos / fromIntegral zsize) + dy/2.0 :: Exp Double
                                            tmin = dt * (fromIntegral tpos / fromIntegral tsize) :: Exp Double 
                                            t = tmin + dt/2.0 :: Exp Double 
                                            tmax = tmin + dt :: Exp Double 
                                            dAz' = dt/ 2.0 * (rmax^2.0 - rmin^2.0) :: Exp Double
                                            dV = dAz'*dz :: Exp Double
                                            dAr' r = dz *dt * r :: Exp Double -> Exp Double 
                                            dAt' = dr * dz :: Exp Double
                                            dAru = lift $ V3 ((cos t) * (dAr' rmin)) ((sin t) * (dAr' rmin)) 0.0 :: Exp (V3 Double)
                                            dArd = lift $ V3 ((cos t) * (dAr' rmax)) ((sin t) * (dAr' rmax)) 0.0 :: Exp (V3 Double)
                                            dAtu = lift $ V3 (-1*(sin tmin)*dAt') ((cos tmin)*dAt') 0.0 :: Exp (V3 Double)
                                            dAtd = lift $ V3 (-1*(sin tmax)*dAt') ((cos tmax)*dAt') 0.0 :: Exp (V3 Double)
                                            dAz = lift $ V3 0 0 (negate dAz') :: Exp (V3 Double) 
                                            patch = lift $ V3 (dAru,dArd) (dAz,dAz) (dAtu,dAtd) :: Exp (V3 (V3 Double, V3 Double))
                                            geom = lift $ (patch,dV) :: Exp (Geom V3) 
                                            position = lift $ V3 r z t :: Exp (V3 Double) 


generateGeometry :: (Elt a, Shape sh) => (Exp sh -> Exp sh -> Exp a) -> Exp sh -> Acc (Array sh a) 
generateGeometry gener size = generate size (gener size) 