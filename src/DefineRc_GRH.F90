#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

module DefineRc_GRH

  implicit none
  contains

    subroutine ICPertFLRW_GRH_DefineRc ( x, y, z, &
                                     Rc, dxRc, dyRc, dzRc, &
                                     dxdxRc, dxdyRc, dxdzRc, dydyRc, dydzRc, dzdzRc, &
                                     dxdxdxRc, dxdxdyRc, dxdxdzRc, dxdydyRc, dxdydzRc, dxdzdzRc, &
                                     dydydyRc, dydydzRc, dydzdzRc, dzdzdzRc, &
                                     dxdxdxdxRc, dxdxdxdyRc, dxdxdxdzRc, dxdxdydyRc, dxdxdydzRc, &
                                     dxdxdzdzRc, dxdydydyRc, dxdydydzRc, dxdydzdzRc, dxdzdzdzRc, &
                                     dydydydyRc, dydydydzRc, dydydzdzRc, dydzdzdzRc, dzdzdzdzRc )
      DECLARE_CCTK_PARAMETERS
      DECLARE_CCTK_FUNCTIONS
      ! Input
      CCTK_REAL :: x, y, z

      ! Output
      CCTK_REAL :: Rc, dxRc, dyRc, dzRc
      CCTK_REAL :: dxdxRc, dxdyRc, dxdzRc, dydyRc, dydzRc, dzdzRc
      CCTK_REAL :: dxdxdxRc, dxdxdyRc, dxdxdzRc, dxdydyRc, dxdydzRc, dxdzdzRc
      CCTK_REAL :: dydydyRc, dydydzRc, dydzdzRc, dzdzdzRc
      CCTK_REAL :: dxdxdxdxRc, dxdxdxdyRc, dxdxdxdzRc, dxdxdydyRc, dxdxdydzRc
      CCTK_REAL :: dxdxdzdzRc, dxdydydyRc, dxdydydzRc, dxdydzdzRc, dxdzdzdzRc
      CCTK_REAL :: dydydydyRc, dydydydzRc, dydydzdzRc, dydzdzdzRc, dzdzdzdzRc
      
      ! exp
      CCTK_REAL :: tsteepx, tsteepy, tsteepz
      CCTK_REAL :: expfacx, expfacy, expfacz
      CCTK_REAL :: xts, yts, zts
      CCTK_REAL :: dxxts, dyyts, dzzts
      CCTK_REAL :: dxdxxts, dydyyts, dzdzzts
      CCTK_REAL :: dxdxdxxts, dydydyyts, dzdzdzzts
      CCTK_REAL :: dxdxdxdxxts, dydydydyyts, dzdzdzdzzts
      CCTK_REAL :: expx, expy, expz
      CCTK_REAL :: dxexpx, dyexpy, dzexpz
      CCTK_REAL :: dxdxexpx, dydyexpy, dzdzexpz
      CCTK_REAL :: dxdxdxexpx, dydydyexpy, dzdzdzexpz
      CCTK_REAL :: dxdxdxdxexpx, dydydydyexpy, dzdzdzdzexpz

      ! Local variables
      integer, parameter :: dp = 8
      logical :: SinRc, CosRc, ExpRc
      ! sin
      integer :: m
      real(dp), parameter :: pi = 4._dp*atan(1._dp)
      CCTK_REAL :: twopi, kx, ky, kz, kxy
      CCTK_REAL :: sinx, siny, sinz, sinxy
      CCTK_REAL :: cosx, cosy, cosz, cosxy
      twopi = 2._dp * pi

      SinRc = CCTK_EQUALS (ICPertFLRW_GRH_Rcprofile, "sin")
      CosRc = CCTK_EQUALS (ICPertFLRW_GRH_Rcprofile, "cos")
      ExpRc = CCTK_EQUALS (ICPertFLRW_GRH_Rcprofile, "exp")

      ! Initialise them to zero
      Rc = 0._dp
      ! Single derivatives
      dxRc = 0._dp
      dyRc = 0._dp
      dzRc = 0._dp
      ! Double derivatives
      dxdxRc = 0._dp
      dxdyRc = 0._dp
      dxdzRc = 0._dp
      dydyRc = 0._dp
      dydzRc = 0._dp
      dzdzRc = 0._dp
      ! Triple derivatives
      dxdxdxRc = 0._dp
      dxdxdyRc = 0._dp
      dxdxdzRc = 0._dp
      dxdydyRc = 0._dp
      dxdydzRc = 0._dp
      dxdzdzRc = 0._dp
      dydydyRc = 0._dp
      dydydzRc = 0._dp
      dydzdzRc = 0._dp
      dzdzdzRc = 0._dp
      ! Quadruple derivatives
      dxdxdxdxRc = 0._dp
      dxdxdxdyRc = 0._dp
      dxdxdxdzRc = 0._dp
      dxdxdydyRc = 0._dp
      dxdxdydzRc = 0._dp
      dxdxdzdzRc = 0._dp
      dxdydydyRc = 0._dp
      dxdydydzRc = 0._dp
      dxdydzdzRc = 0._dp
      dxdzdzdzRc = 0._dp
      dydydydyRc = 0._dp
      dydydydzRc = 0._dp
      dydydzdzRc = 0._dp
      dydzdzdzRc = 0._dp
      dzdzdzdzRc = 0._dp

      if (SinRc) then
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !           Sinusoidal
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do m = 1, 20
          kx = twopi / ICPertFLRW_GRH_lambda_x(m)
          ky = twopi / ICPertFLRW_GRH_lambda_y(m)
          kz = twopi / ICPertFLRW_GRH_lambda_z(m)

          sinx = sin(x * kx + ICPertFLRW_GRH_phi_x(m))
          siny = sin(y * ky + ICPertFLRW_GRH_phi_y(m))
          sinz = sin(z * kz + ICPertFLRW_GRH_phi_z(m))
          cosx = cos(x * kx + ICPertFLRW_GRH_phi_x(m))
          cosy = cos(y * ky + ICPertFLRW_GRH_phi_y(m))
          cosz = cos(z * kz + ICPertFLRW_GRH_phi_z(m))
  
          Rc = Rc + ICPertFLRW_GRH_Amp_x(m) * sinx &
                  + ICPertFLRW_GRH_Amp_y(m) * siny &
                  + ICPertFLRW_GRH_Amp_z(m) * sinz

          ! Single derivatives
          dxRc = dxRc + ICPertFLRW_GRH_Amp_x(m) * kx * cosx
          dyRc = dyRc + ICPertFLRW_GRH_Amp_y(m) * ky * cosy
          dzRc = dzRc + ICPertFLRW_GRH_Amp_z(m) * kz * cosz

          ! Double derivatives
          dxdxRc = dxdxRc - ICPertFLRW_GRH_Amp_x(m) * (kx**2._dp) * sinx
          dydyRc = dydyRc - ICPertFLRW_GRH_Amp_y(m) * (ky**2._dp) * siny
          dzdzRc = dzdzRc - ICPertFLRW_GRH_Amp_z(m) * (kz**2._dp) * sinz

          ! Triple derivatives
          dxdxdxRc = dxdxdxRc - ICPertFLRW_GRH_Amp_x(m) * (kx**3._dp) * cosx
          dydydyRc = dydydyRc - ICPertFLRW_GRH_Amp_y(m) * (ky**3._dp) * cosy
          dzdzdzRc = dzdzdzRc - ICPertFLRW_GRH_Amp_z(m) * (kz**3._dp) * cosz
 
          ! Quadruple derivatives
          dxdxdxdxRc = dxdxdxdxRc + ICPertFLRW_GRH_Amp_x(m) * (kx**4._dp) * sinx
          dydydydyRc = dydydydyRc + ICPertFLRW_GRH_Amp_y(m) * (ky**4._dp) * siny
          dzdzdzdzRc = dzdzdzdzRc + ICPertFLRW_GRH_Amp_z(m) * (kz**4._dp) * sinz
        enddo
      else if (CosRc) then
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !           Cosine
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do m = 1, 20
          kx = twopi / ICPertFLRW_GRH_lambda_x(m)
          ky = twopi / ICPertFLRW_GRH_lambda_y(m)
          kz = twopi / ICPertFLRW_GRH_lambda_z(m)
          kxy = twopi / ICPertFLRW_GRH_lambda_xy(m)

          cosx = cos(x * kx)
          cosy = cos(y * ky)
          cosz = cos(z * kz)
          sinx = sin(x * kx)
          siny = sin(y * ky)
          sinz = sin(z * kz)
          cosxy = cos((x - y) * kxy)
          sinxy = sin((x - y) * kxy)

          Rc = Rc + ICPertFLRW_GRH_Amp_x(m) * cosx &
                  + ICPertFLRW_GRH_Amp_y(m) * cosy &
                  + ICPertFLRW_GRH_Amp_z(m) * cosz &
                  + ICPertFLRW_GRH_Amp_xy(m) * cosxy

          ! Single derivatives
          dxRc = dxRc - ICPertFLRW_GRH_Amp_x(m) * kx * sinx &
                      - ICPertFLRW_GRH_Amp_xy(m) * kxy * sinxy
          dyRc = dyRc - ICPertFLRW_GRH_Amp_y(m) * ky * siny &
                      + ICPertFLRW_GRH_Amp_xy(m) * kxy * sinxy
          dzRc = dzRc - ICPertFLRW_GRH_Amp_z(m) * kz * sinz

          ! Double derivatives
          dxdxRc = dxdxRc - ICPertFLRW_GRH_Amp_x(m) * (kx**2._dp) * cosx &
                          - ICPertFLRW_GRH_Amp_xy(m) * (kxy**2._dp) * cosxy
          dxdyRc = dxdyRc + ICPertFLRW_GRH_Amp_xy(m) * (kxy**2._dp) * cosxy
          dydyRc = dydyRc - ICPertFLRW_GRH_Amp_y(m) * (ky**2._dp) * cosy &
                          - ICPertFLRW_GRH_Amp_xy(m) * (kxy**2._dp) * cosxy
          dzdzRc = dzdzRc - ICPertFLRW_GRH_Amp_z(m) * (kz**2._dp) * cosz

          ! Triple derivatives
          dxdxdxRc = dxdxdxRc + ICPertFLRW_GRH_Amp_x(m) * (kx**3._dp) * sinx &
                              + ICPertFLRW_GRH_Amp_xy(m) * (kxy**3._dp) * sinxy
          dxdxdyRc = dxdxdyRc - ICPertFLRW_GRH_Amp_xy(m) * (kxy**3._dp) * sinxy
          dxdydyRc = dxdydyRc + ICPertFLRW_GRH_Amp_xy(m) * (kxy**3._dp) * sinxy
          dydydyRc = dydydyRc + ICPertFLRW_GRH_Amp_y(m) * (ky**3._dp) * siny &
                              - ICPertFLRW_GRH_Amp_xy(m) * (kxy**3._dp) * sinxy
          dzdzdzRc = dzdzdzRc + ICPertFLRW_GRH_Amp_z(m) * (kz**3._dp) * sinz

          ! Quadruple derivatives
          dxdxdxdxRc = dxdxdxdxRc + ICPertFLRW_GRH_Amp_x(m) * (kx**4._dp) * cosx &
                                  + ICPertFLRW_GRH_Amp_xy(m) * (kxy**4._dp) * cosxy
          dxdxdxdyRc = dxdxdxdyRc - ICPertFLRW_GRH_Amp_xy(m) * (kxy**4._dp) * cosxy
          dxdxdydyRc = dxdxdydyRc + ICPertFLRW_GRH_Amp_xy(m) * (kxy**4._dp) * cosxy
          dxdydydyRc = dxdydydyRc - ICPertFLRW_GRH_Amp_xy(m) * (kxy**4._dp) * cosxy
          dydydydyRc = dydydydyRc + ICPertFLRW_GRH_Amp_y(m) * (ky**4._dp) * cosy &
                                  + ICPertFLRW_GRH_Amp_xy(m) * (kxy**4._dp) * cosxy
          dzdzdzdzRc = dzdzdzdzRc + ICPertFLRW_GRH_Amp_z(m) * (kz**4._dp) * cosz
        enddo

      else if (ExpRc) then
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !           Exponential
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          tsteepx = 2._dp * ICPertFLRW_GRH_steepness_x
          tsteepy = 2._dp * ICPertFLRW_GRH_steepness_y
          tsteepz = 2._dp * ICPertFLRW_GRH_steepness_z

          expfacx = (- 1._dp / 2._dp) * (ICPertFLRW_GRH_variance_x) ** (- tsteepx)
          expfacy = (- 1._dp / 2._dp) * (ICPertFLRW_GRH_variance_y) ** (- tsteepy)
          expfacz = (- 1._dp / 2._dp) * (ICPertFLRW_GRH_variance_z) ** (- tsteepz)

          xts = x ** tsteepx
          yts = y ** tsteepy
          zts = z ** tsteepz

          expx = exp( expfacx * xts)
          expy = exp( expfacy * yts)
          expz = exp( expfacz * zts)

          Rc = ICPertFLRW_GRH_exp_amplitude * expx * expy * expz

          ! x^2s derivarives
          dxxts = 0._dp
          dxdxxts = 0._dp
          dxdxdxxts = 0._dp
          dxdxdxdxxts = 0._dp
          if (tsteepz >= 2._dp) then
              dxxts = tsteepx * (x ** (tsteepx - 1._dp))
              dxdxxts = tsteepx * (tsteepx - 1._dp) * (x ** (tsteepx - 2._dp))
              if (tsteepz >= 4._dp) then
                  dxdxdxxts = (tsteepx * (tsteepx - 1._dp) * (tsteepx - 2._dp)    &
                               * (x ** (tsteepx - 3._dp)))
                  dxdxdxdxxts = (tsteepx * (tsteepx - 1._dp) * (tsteepx - 2._dp)  &
                                 * (tsteepx - 3._dp) * (x ** (tsteepx - 4._dp)))
              endif
          endif

          ! y^2s derivarives
          dyyts = 0._dp
          dydyyts = 0._dp
          dydydyyts = 0._dp
          dydydydyyts = 0._dp
          if (tsteepz >= 2._dp) then
              dyyts = tsteepy * (y ** (tsteepy - 1._dp))
              dydyyts = tsteepy * (tsteepy - 1._dp) * (y ** (tsteepy - 2._dp))
              if (tsteepz >= 4._dp) then
                  dydydyyts = (tsteepy * (tsteepy - 1._dp) * (tsteepy - 2._dp)    &
                               * (y ** (tsteepy - 3._dp)))
                  dydydydyyts = (tsteepy * (tsteepy - 1._dp) * (tsteepy - 2._dp)  &
                                 * (tsteepy - 3._dp) * (y ** (tsteepy - 4._dp)))
              endif
          endif

          ! z^2s derivarives
          dzzts = 0._dp
          dzdzzts = 0._dp
          dzdzdzzts = 0._dp
          dzdzdzdzzts = 0._dp
          if (tsteepz >= 2._dp) then
              dzzts = tsteepz * (z ** (tsteepz - 1._dp))
              dzdzzts = tsteepz * (tsteepz - 1._dp) * (z ** (tsteepz - 2._dp))
              if (tsteepz >= 4._dp) then
                  dzdzdzzts = (tsteepz * (tsteepz - 1._dp) * (tsteepz - 2._dp)    &
                               * (z ** (tsteepz - 3._dp)))
                  dzdzdzdzzts = (tsteepz * (tsteepz - 1._dp) * (tsteepz - 2._dp)  &
                                 * (tsteepz - 3._dp) * (z ** (tsteepz - 4._dp)))
              endif
          endif
          
          ! Single derivatives
          dxexpx = expfacx * dxxts * expx
          dyexpy = expfacy * dyyts * expy
          dzexpz = expfacz * dzzts * expz
          dxRc = ICPertFLRW_GRH_exp_amplitude * dxexpx * expy * expz
          dyRc = ICPertFLRW_GRH_exp_amplitude * expx * dyexpy * expz
          dzRc = ICPertFLRW_GRH_exp_amplitude * expx * expy * dzexpz

          ! Double derivatives
          dxdxexpx = expfacx * (dxdxxts * expx + dxxts * dxexpx)
          dydyexpy = expfacy * (dydyyts * expy + dyyts * dyexpy)
          dzdzexpz = expfacz * (dzdzzts * expz + dzzts * dzexpz)
          dxdxRc = ICPertFLRW_GRH_exp_amplitude * dxdxexpx * expy * expz
          dxdyRc = ICPertFLRW_GRH_exp_amplitude * dxexpx * dyexpy * expz
          dxdzRc = ICPertFLRW_GRH_exp_amplitude * dxexpx * expy * dzexpz
          dydyRc = ICPertFLRW_GRH_exp_amplitude * expx * dydyexpy * expz
          dydzRc = ICPertFLRW_GRH_exp_amplitude * expx * dyexpy * dzexpz
          dzdzRc = ICPertFLRW_GRH_exp_amplitude * expx * expy * dzdzexpz

          ! Triple derivatives
          dxdxdxexpx = expfacx * (dxdxdxxts * expx           &
                                  + 2._dp * dxdxxts * dxexpx &
                                  + dxxts * dxdxexpx)
          dydydyexpy = expfacy * (dydydyyts * expy           &
                                  + 2._dp * dydyyts * dyexpy &
                                  + dyyts * dydyexpy)
          dzdzdzexpz = expfacz * (dzdzdzzts * expz           &
                                  + 2._dp * dzdzzts * dzexpz &
                                  + dzzts * dzdzexpz)
          dxdxdxRc = ICPertFLRW_GRH_exp_amplitude * dxdxdxexpx * expy * expz
          dxdxdyRc = ICPertFLRW_GRH_exp_amplitude * dxdxexpx * dyexpy * expz
          dxdxdzRc = ICPertFLRW_GRH_exp_amplitude * dxdxexpx * expy * dzexpz
          dxdydyRc = ICPertFLRW_GRH_exp_amplitude * dxexpx * dydyexpy * expz
          dxdydzRc = ICPertFLRW_GRH_exp_amplitude * dxexpx * dyexpy * dzexpz
          dxdzdzRc = ICPertFLRW_GRH_exp_amplitude * dxexpx * expy * dzdzexpz
          dydydyRc = ICPertFLRW_GRH_exp_amplitude * expx * dydydyexpy * expz
          dydydzRc = ICPertFLRW_GRH_exp_amplitude * expx * dydyexpy * dzexpz
          dydzdzRc = ICPertFLRW_GRH_exp_amplitude * expx * dyexpy * dzdzexpz
          dzdzdzRc = ICPertFLRW_GRH_exp_amplitude * expx * expy * dzdzdzexpz

          ! Quadruple derivatives
          dxdxdxdxexpx = expfacx * (dxdxdxdxxts * expx + 3._dp * dxdxdxxts * dxexpx    &
                                    + 3._dp * dxdxxts * dxdxexpx + dxxts * dxdxdxexpx)
          dydydydyexpy = expfacy * (dydydydyyts * expy + 3._dp * dydydyyts * dyexpy    &
                                    + 3._dp * dydyyts * dydyexpy + dyyts * dydydyexpy)
          dzdzdzdzexpz = expfacz * (dzdzdzdzzts * expz + 3._dp * dzdzdzzts * dzexpz    &
                                    + 3._dp * dzdzzts * dzdzexpz + dzzts * dzdzdzexpz)
          dxdxdxdxRc = ICPertFLRW_GRH_exp_amplitude * dxdxdxdxexpx * expy * expz
          dxdxdxdyRc = ICPertFLRW_GRH_exp_amplitude * dxdxdxexpx * dyexpy * expz
          dxdxdxdzRc = ICPertFLRW_GRH_exp_amplitude * dxdxdxexpx * expy * dzexpz
          dxdxdydyRc = ICPertFLRW_GRH_exp_amplitude * dxdxexpx * dydyexpy * expz
          dxdxdydzRc = ICPertFLRW_GRH_exp_amplitude * dxdxexpx * dyexpy * dzexpz
          dxdxdzdzRc = ICPertFLRW_GRH_exp_amplitude * dxdxexpx * expy * dzdzexpz
          dxdydydyRc = ICPertFLRW_GRH_exp_amplitude * dxexpx * dydydyexpy * expz
          dxdydydzRc = ICPertFLRW_GRH_exp_amplitude * dxexpx * dydyexpy * dzexpz
          dxdydzdzRc = ICPertFLRW_GRH_exp_amplitude * dxexpx * dyexpy * dzdzexpz
          dxdzdzdzRc = ICPertFLRW_GRH_exp_amplitude * dxexpx * expy * dzdzdzexpz
          dydydydyRc = ICPertFLRW_GRH_exp_amplitude * expx * dydydydyexpy * expz
          dydydydzRc = ICPertFLRW_GRH_exp_amplitude * expx * dydydyexpy * dzexpz
          dydydzdzRc = ICPertFLRW_GRH_exp_amplitude * expx * dydyexpy * dzdzexpz
          dydzdzdzRc = ICPertFLRW_GRH_exp_amplitude * expx * dyexpy * dzdzdzexpz
          dzdzdzdzRc = ICPertFLRW_GRH_exp_amplitude * expx * expy * dzdzdzdzexpz

      endif
  end subroutine ICPertFLRW_GRH_DefineRc
end module DefineRc_GRH
