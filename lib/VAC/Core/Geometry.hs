--Copyright 2018 GeneralFusion Inc.
--
--Licensed under the Apache License, Version 2.0 (the "License");
--you may not use this file except in compliance with the License.
--You may obtain a copy of the License at
--
--    http://www.apache.org/licenses/LICENSE-2.0
--
--Unless required by applicable law or agreed to in writing, software
--distributed under the License is distributed on an "AS IS" BASIS,
--WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
--See the License for the specific language governing permissions and
--limitations under the License.
module VAC.Core.Geometry where 

import qualified Prelude as P 

import Data.Array.Accelerate as Acc
import Data.Array.Accelerate.Linear hiding (fromQuaternion)
import Linear.Matrix (fromQuaternion)

import VAC.Core.Types as Types

delta :: (Precision, Precision) -> Exp Int -> Exp Precision
delta (rmin,rmax) n = (constant $ rmax - rmin) / fromIntegral n 
                      

position :: (Precision,Precision) -> Exp Int -> Exp Int -> Exp Precision 
position (rmin,rmax) n i = dx*(fromIntegral i) / (fromIntegral n) + dx/2.0 + constant rmin
                                where 
                                    dx = delta (rmin,rmax) n 

cartesianVoxel :: Exp Precision -> Exp Precision -> Exp Precision -> Exp (Cell V3) 
cartesianVoxel dx dy dz = lift $ (patch,dV)
                            where
                                patch = lift $ V3 (dAx,dAx) (dAy,dAy) (dAz,dAz) :: Exp (Faces V3)
                                dAx = lift $ V3 (dy*dz) 0 0 :: Exp (V3 Precision)
                                dAy = lift $ V3 0 (dx*dz) 0 :: Exp (V3 Precision)
                                dAz = lift $ V3 0 0 (dx*dy) :: Exp (V3 Precision)
                                dV = dx*dy*dz :: Exp Precision


cartesian3D :: (Precision,Precision) -> (Precision,Precision) -> (Precision,Precision) -> Exp DIM3 -> Exp DIM3 -> Exp (V3 Precision,Cell V3)
cartesian3D xrange yrange zrange dimensions index = lift (location, geom)
                                    where
                                        (xpos,ypos,zpos) = unlift $ unindex3 index :: (Exp Int,Exp Int,Exp Int) 
                                        (xsize,ysize,zsize) = unlift $ unindex3 dimensions :: (Exp Int,Exp Int,Exp Int) 
                                        dx = delta xrange xsize :: Exp Precision
                                        dy = delta yrange ysize :: Exp Precision
                                        dz = delta zrange zsize :: Exp Precision
                                        geom = cartesianVoxel dx dy dz :: Exp (Cell V3)
                                        x = position xrange xsize xpos :: Exp Precision 
                                        y = position yrange ysize ypos :: Exp Precision
                                        z = position zrange zsize zpos :: Exp Precision
                                        location = lift $ V3 x y z :: Exp (V3 Precision)

--creates a complete cylindrical voxel centered on theta = 0
cylindricalVoxel :: Exp Precision -> Exp Precision -> Exp Precision -> Exp Precision -> Exp (Cell V3)
cylindricalVoxel dr dz dt r = lift (patch,dV) where
                                    dAr' rc = V3 (rc*dt*dz) 0.0 0.0 
                                    rp = r - dr / 2.0 
                                    rn = r + dr / 2.0
                                    dAt' t = V3 (-1*(sin t)*dr*dz) ((cos tp)*dr*dz) 0.0
                                    tp = -1.0*dt/2.0
                                    tn = dt/2.0
                                    dAz = V3 0.0 0.0 (negate $ dr*dt*r)
                                    patch = V3 (dAr' rp, dAr' rn) (dAz,dAz) (dAt' tp,dAt' tn) 
                                    dV = dr*dt*r*dz

rotateAboutZ :: Exp Precision -> Quaternion (Exp Precision) 
rotateAboutZ amount = Quaternion (cos halfAngle) $ unlift imag where
            axis = constant $ V3 0.0 0.0 1.0 :: Exp (V3 Precision) 
            halfAngle = amount / 2.0 :: Exp Precision
            imag = axis ^* sin halfAngle :: Exp (V3 Precision) 

rotatePatchPair :: Exp (M33 Precision) -> Exp (V3 Precision,V3 Precision) -> Exp (V3 Precision, V3 Precision)
rotatePatchPair rot patches = lift (rot !* p1, rot !* p2) where (p1,p2) = unlift patches

rotateCell :: Quaternion (Exp Precision) -> Exp (Cell V3) -> Exp (Cell V3) 
rotateCell rotquat cell = lift $ (newpatch,volume)  where 
        (patches, volume) = unlift cell :: (Exp (Faces V3), Exp Precision)
        rotmat = lift $ fromQuaternion rotquat :: Exp (M33 Precision)
        patchpairs = unlift patches :: V3 (Exp PatchPair)
        mvfunc pair = rotatePatchPair rotmat pair :: Exp PatchPair 
        newpairs = P.fmap mvfunc patchpairs :: V3 (Exp PatchPair)
        newpatch = lift newpairs :: Exp (Faces V3) 

sweepCell ::Exp Precision -> Exp (Cell V3) -> Exp (Cell V3)
sweepCell amount input = rotateCell (rotateAboutZ amount) input

cylindrical3D :: (Precision,Precision) -> (Precision,Precision) -> Exp DIM3 -> Exp DIM3 -> Exp (V3 Precision,Cell V3)
cylindrical3D rrange zrange dimensions index = lift (pos,newvox)  
                                        where
                                            trange = (-pi,pi)
                                            (rpos,zpos,tpos) = unlift $ unindex3 index :: (Exp Int,Exp Int,Exp Int) 
                                            (rsize,zsize,tsize) = unlift $ unindex3 dimensions :: (Exp Int,Exp Int,Exp Int) 
                                            dr = delta rrange rsize :: Exp Precision 
                                            dz = delta zrange zsize :: Exp Precision 
                                            dt = delta trange tsize :: Exp Precision 
                                            r = position rrange rsize rpos :: Exp Precision
                                            z = position zrange zsize zpos :: Exp Precision 
                                            t = position trange tsize tpos - dt / 2.0 :: Exp Precision -- we want 0 aligned with first cell
                                            vox = cylindricalVoxel dr dz dt r :: Exp (Cell V3) 
                                            newvox = sweepCell t vox :: Exp (Cell V3) 
                                            pos = lift $ V3 r z t :: Exp (V3 Precision)


generateGeometry :: (Elt a, Shape sh) => (Exp sh -> Exp sh -> Exp a) -> Exp sh -> Acc (Array sh a) 
generateGeometry gener size = generate size $ gener size
