module VAC.Core.Geometry where 

import qualified Prelude as P 
import Data.Array.Accelerate.Control.Lens
import Data.Array.Accelerate as Acc
import Data.Array.Accelerate.Linear
import Linear.Matrix (fromQuaternion)

import VAC.Core.Types as Types

delta :: (Precision, Precision) -> Exp Int -> Exp Precision
delta (rmin,rmax) n = (constant $ rmax - rmin) / fromIntegral n 

cartesian1D :: (Precision,Precision) -> Exp DIM1 -> Acc Geometry1D
cartesian1D xrange dim = 
        let 
            lenx = dim ^. _1
            dx = delta xrange lenx
            padlen = dim & _1 .~ (lenx+1)
        in 
            lift (fill dim dx,fill padlen 1)



cartesian2D :: (Precision, Precision) -> (Precision,Precision) -> Exp DIM2 -> Acc Geometry2D
cartesian2D xrange yrange dim = 
        let 
            lenx = dim ^. _1
            leny = dim ^. _2
            dx = delta xrange lenx
            dy = delta yrange leny
            padx = dim & _1 .~ (lenx+1)
            pady = dim & _2 .~ (leny+1)
            nx = (zero & _x .~ dy) :: Exp (V2 Precision)
            ny = (zero & _y .~ dx) :: Exp (V2 Precision)
        in
            lift (fill dim (dx*dy),fill padx nx, fill pady ny)


cartesian3D :: (Precision, Precision) -> (Precision,Precision) -> (Precision,Precision) -> Exp DIM3 -> Acc Geometry3D
cartesian3D xrange yrange zrange dim = 
        let 
            lenx = dim ^. _1
            leny = dim ^. _2
            lenz = dim ^. _3
            dx = delta xrange lenx
            dy = delta yrange leny
            dz = delta zrange lenz
            padx = dim & _1 .~ (lenx+1)
            pady = dim & _2 .~ (leny+1)
            padz = dim & _3 .~ (lenz+1)
            nx =  zero & _x .~ (dy*dz) :: Exp (V3 Precision)
            ny =  zero & _y .~ (dx*dz) :: Exp (V3 Precision)
            nz =  zero & _z .~ (dx*dy) :: Exp (V3 Precision)
        in
            lift (fill dim (dx*dy*dz),fill padx nx, fill pady ny,fill padz nz)


cylindrical1D :: (Precision,Precision) -> Exp DIM1 -> Acc Geometry1D
cylindrical1D rrange dim = 
    let 
        lenr = dim ^. _1
        dr = delta rrange lenr
        padr = dim & _1 .~ (lenr+1)
        dz = constant 1.0 :: Exp Precision
        dt = constant $ 2*pi :: Exp Precision
        rmin = constant $ rrange ^. _1 :: Exp Precision
        getR :: Exp DIM1 -> Exp Precision
        getR sh = rmin + dr*(fromIntegral idx) where idx = (sh ^. _1) :: Exp Int
        voxelR = map (\v -> v + dr/2.0) $ generate dim getR
        voxelV = fill dim (dr*dz*dt)
        nR = map (\v -> v * dt *dz) $ generate padr getR
    in 
        lift (zipWith (*) voxelR voxelV, nR) 


cylindrical2D :: (Precision,Precision) -> (Precision,Precision) -> Exp DIM2 -> Acc Geometry2D
cylindrical2D rrange zrange dim = 
    let 
        lenr = dim ^. _1
        lenz = dim ^. _2
        dr = delta rrange lenr
        dz = delta zrange lenz
        padr = dim & _1 .~ (lenr+1)
        padz = dim & _2 .~ (lenz+1)
        dt = constant $ 2*pi :: Exp Precision
        rmin = constant $ rrange ^. _1 :: Exp Precision
        getR :: Exp DIM2 -> Exp Precision
        getR sh = rmin + dr*(fromIntegral idx) where idx = (sh ^. _1) :: Exp Int
        voxelV = map (\v -> (v + dr/2.0)*dr*dz*dt) $ generate dim getR
        nR = map (\v -> zero & _x .~ (v * dt *dz)) $ generate padr getR
        nZ = map (\v -> zero & _y .~ (v * dt * dr)) $ generate padz getR
    in 
        lift (voxelV, nR,nZ) 

rotateAboutZ :: Exp Precision -> Quaternion (Exp Precision) 
rotateAboutZ amount = Quaternion (cos halfAngle) $ unlift imag where
            axis = constant $ V3 0.0 0.0 1.0 :: Exp (V3 Precision) 
            halfAngle = amount / 2.0 :: Exp Precision
            imag = axis ^* sin halfAngle :: Exp (V3 Precision) 

cylindrical3D :: (Precision,Precision) -> (Precision,Precision) -> Exp DIM3 -> Acc Geometry3D
cylindrical3D rrange zrange dim = 
    let 
        lenr = dim ^. _1
        lenz = dim ^. _2
        lent = dim ^. _3
        dr = delta rrange lenr
        dz = delta zrange lenz
        dt = delta (0,2*pi) lent
        padr = dim & _1 .~ (lenr+1)
        padz = dim & _2 .~ (lenz+1)
        padt = dim & _3 .~ (lent+1)

        rmin = constant $ rrange ^. _1 :: Exp Precision
        getR :: Exp DIM3 -> Exp Precision
        getR sh = rmin + dr*(fromIntegral idx) where idx = (sh ^. _1) :: Exp Int

        getT :: Exp DIM3 -> Exp Precision
        getT sh = dt * (fromIntegral (sh ^. _3))

        voxelV = map (\v -> (v + dr/2.0)*dr*dz*dt) $ generate dim getR
        nR = map (\v -> zero & _x .~ (v * dt *dz)) $ generate padr getR
        nZ = map (\v -> zero & _y .~ (v * dt * dr)) $ generate padz getR
        nT = fill padt (zero & _z .~ (dr*dz))

        rotT = map (\v -> lift $ fromQuaternion $ rotateAboutZ v) $ generate padt getT
        rotR = map (\v -> lift $ fromQuaternion $ rotateAboutZ (v+dt/2)) $ generate padr getT
        nT' = zipWith (!*) rotT nT
        nR' = zipWith (!*) rotR nR
    in 
        lift (voxelV, nR',nZ,nT')

