module VAC.Core.Geometry where 

import qualified Prelude as P 

import Data.Array.Accelerate as Acc
import Data.Array.Accelerate.Linear
import VAC.Core.Types as Types


cartesian3D :: (Precision,Precision) -> (Precision,Precision) -> (Precision,Precision) -> Exp DIM3 -> Exp DIM3 -> Exp (V3 Precision,Cell V3)
cartesian3D (xmin,xmax) (ymin,ymax) (zmin,zmax) dimensions index = lift (position, geom)
                                    where
                                        (xpos,ypos,zpos) = unlift $ unindex3 index :: (Exp Int,Exp Int,Exp Int) 
                                        (xsize,ysize,zsize) = unlift $ unindex3 dimensions :: (Exp Int,Exp Int,Exp Int) 
                                        dx = constant $ xmax - xmin :: Exp Precision 
                                        dy = constant $ ymax - ymin :: Exp Precision 
                                        dz = constant $ zmax - zmin :: Exp Precision 
                                        dV = dx*dy*dz :: Exp Precision 
                                        dAx = lift $ V3 (dy*dz) 0 0 :: Exp (V3 Precision) 
                                        dAy = lift $ V3 0 (dx*dz) 0 :: Exp (V3 Precision) 
                                        dAz = lift $ V3 0 0 (dx*dy) :: Exp (V3 Precision) 
                                        patch = lift $ V3 (dAx,dAx) (dAy,dAy) (dAz,dAz) :: Exp (V3 (V3 Precision, V3 Precision))
                                        geom = lift $ (patch,dV) :: Exp (Cell V3) 
                                        x = dx * (fromIntegral xpos / fromIntegral xsize) + dx/2.0 + constant xmin :: Exp Precision 
                                        y = dy * (fromIntegral ypos / fromIntegral ysize) + dy/2.0 + constant ymin :: Exp Precision
                                        z = dz * (fromIntegral zpos / fromIntegral zsize) + dz/2.0 + constant zmin:: Exp Precision 
                                        position = lift $ V3 x y z :: Exp (V3 Precision) 


cylindrical3D :: (Precision,Precision) -> (Precision,Precision) -> (Precision,Precision) -> Exp DIM3 -> Exp DIM3 -> Exp (V3 Precision,Cell V3)
cylindrical3D (rmin,rmax) (zmin,zmax) (tmin,tmax) dimensions index = lift (position,geom)  
                                        where 
                                            (rpos,zpos,tpos) = unlift $ unindex3 index :: (Exp Int,Exp Int,Exp Int) 
                                            (rsize,zsize,tsize) = unlift $ unindex3 dimensions :: (Exp Int,Exp Int,Exp Int) 
                                            dr = constant $ (rmax - rmin) :: Exp Precision 
                                            dz = constant $ (zmax - zmin) :: Exp Precision 
                                            dt = constant $ (tmax - tmin) :: Exp Precision 
                                            rs = dr * (fromIntegral rpos / fromIntegral rsize) + constant rmin :: Exp Precision 
                                            r = rs + dr / 2.0 :: Exp Precision 
                                            re = rs + dr :: Exp Precision 
                                            z = dz * (fromIntegral zpos / fromIntegral zsize) + dz/2.0 + constant zmin:: Exp Precision
                                            ts = dt * (fromIntegral tpos / fromIntegral tsize) + constant tmin:: Exp Precision 
                                            t = ts + dt/2.0 :: Exp Precision 
                                            te = ts + dt :: Exp Precision 
                                            dAz' = dr * r * dt :: Exp Precision
                                            dV = dAz'*dz :: Exp Precision
                                            dAr' rp = V3 ((cos t) * (dz*dt*rp)) ((sin t) * (dz*dt*rp)) 0.0 :: V3 (Exp Precision)
                                            dAru = lift $ dAr' rs  :: Exp (V3 Precision)
                                            dArd = lift $ dAr' re :: Exp (V3 Precision)
                                            dAt' tp = lift $ V3 (-1*(sin tp)*dr*dz) ((cos tp)*dr*dz) 0.0 :: Exp (V3 Precision)
                                            dAtu = dAt' ts :: Exp (V3 Precision)
                                            dAtd = dAt' te :: Exp (V3 Precision)
                                            dAz = lift $ V3 0 0 (negate dAz') :: Exp (V3 Precision) 
                                            patch = lift $ V3 (dAru,dArd) (dAz,dAz) (dAtu,dAtd) :: Exp (V3 (V3 Precision, V3 Precision))
                                            geom = lift $ (patch,dV) :: Exp (Cell V3) 
                                            position = lift $ V3 r z t :: Exp (V3 Precision) 


generateGeometry :: (Elt a, Shape sh) => (Exp sh -> Exp sh -> Exp a) -> Exp sh -> Acc (Array sh a) 
generateGeometry gener size = generate size (gener size) 