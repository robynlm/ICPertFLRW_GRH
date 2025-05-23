! Examples found in:
! EinsteinInitialData/Hydro_RNSID/src/rnsid_rfr.c  ==> Hydro_rnsid_init
! EinsteinInitialData/NoExcision/src/overwriteBSSN.F90  ==> NoExcision_OverwriteBSSN
! EinsteinInitialData/TOVSolver/src/tov.c  ==> TOV_C_Exact
! GRHydro_InitData/src/GRHydro_TOV.F90

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine ICPertFLRW_GRH_TimelevelCopy (CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  integer :: i,j,k
  integer :: timelevels
  
  ! ==================== Metric
  call CCTK_ActiveTimeLevelsVN (timelevels, cctkGH, "ADMBase::gxx")
  if (timelevels > 1) then
    do k = 1, cctk_lsh(3)
      do j = 1, cctk_lsh(2)
        do i = 1, cctk_lsh(1)
          gxx_p(i,j,k) = gxx(i,j,k)
          gxy_p(i,j,k) = gxy(i,j,k)
          gxz_p(i,j,k) = gxz(i,j,k)
          gyy_p(i,j,k) = gyy(i,j,k)
          gyz_p(i,j,k) = gyz(i,j,k)
          gzz_p(i,j,k) = gzz(i,j,k)
        enddo
      enddo
    enddo
  end if
  if (timelevels > 2) then
    do k = 1, cctk_lsh(3)
      do j = 1, cctk_lsh(2)
        do i = 1, cctk_lsh(1)
          gxx_p_p(i,j,k) = gxx(i,j,k)
          gxy_p_p(i,j,k) = gxy(i,j,k)
          gxz_p_p(i,j,k) = gxz(i,j,k)
          gyy_p_p(i,j,k) = gyy(i,j,k)
          gyz_p_p(i,j,k) = gyz(i,j,k)
          gzz_p_p(i,j,k) = gzz(i,j,k)
        enddo
      enddo
    enddo
  end if

  ! ==================== Curvature
  call CCTK_ActiveTimeLevelsVN (timelevels, cctkGH, "ADMBase::kxx")
  if (timelevels > 1) then
    do k = 1, cctk_lsh(3)
      do j = 1, cctk_lsh(2)
        do i = 1, cctk_lsh(1)
          kxx_p(i,j,k) = kxx(i,j,k)
          kxy_p(i,j,k) = kxy(i,j,k)
          kxz_p(i,j,k) = kxz(i,j,k)
          kyy_p(i,j,k) = kyy(i,j,k)
          kyz_p(i,j,k) = kyz(i,j,k)
          kzz_p(i,j,k) = kzz(i,j,k)
        enddo
      enddo
    enddo
  end if
  if (timelevels > 2) then
    do k = 1, cctk_lsh(3)
      do j = 1, cctk_lsh(2)
        do i = 1, cctk_lsh(1)
          kxx_p_p(i,j,k) = kxx(i,j,k)
          kxy_p_p(i,j,k) = kxy(i,j,k)
          kxz_p_p(i,j,k) = kxz(i,j,k)
          kyy_p_p(i,j,k) = kyy(i,j,k)
          kyz_p_p(i,j,k) = kyz(i,j,k)
          kzz_p_p(i,j,k) = kzz(i,j,k)
        enddo
      enddo
    enddo
  end if

  ! ==================== Lapse
  call CCTK_ActiveTimeLevelsVN (timelevels, cctkGH, "ADMBase::alp")
  if (timelevels > 1) then
    do k = 1, cctk_lsh(3)
      do j = 1, cctk_lsh(2)
        do i = 1, cctk_lsh(1)
          alp_p(i,j,k) = alp(i,j,k)
          dtalp_p(i,j,k) = dtalp(i,j,k)
        enddo
      enddo
    enddo
  end if
  if (timelevels > 2) then
    do k = 1, cctk_lsh(3)
      do j = 1, cctk_lsh(2)
       do i = 1, cctk_lsh(1)
         alp_p_p(i,j,k) = alp(i,j,k)
         dtalp_p_p(i,j,k) = dtalp(i,j,k)
        enddo
      enddo
    enddo
  end if

  ! ==================== Shift
  call CCTK_ActiveTimeLevelsVN (timelevels, cctkGH, "ADMBase::betax")
  if (timelevels > 1) then
    do k = 1, cctk_lsh(3)
      do j = 1, cctk_lsh(2)
        do i = 1, cctk_lsh(1)
          betax_p(i,j,k) = betax(i,j,k)
          betay_p(i,j,k) = betay(i,j,k)
          betaz_p(i,j,k) = betaz(i,j,k)
          dtbetax_p(i,j,k) = dtbetax(i,j,k)
          dtbetay_p(i,j,k) = dtbetay(i,j,k)
          dtbetaz_p(i,j,k) = dtbetaz(i,j,k)
        enddo
      enddo
    enddo
  end if
  if (timelevels > 2) then
    do k = 1, cctk_lsh(3)
      do j = 1, cctk_lsh(2)
        do i = 1, cctk_lsh(1)
          betax_p_p(i,j,k) = betax(i,j,k)
          betay_p_p(i,j,k) = betay(i,j,k)
          betaz_p_p(i,j,k) = betaz(i,j,k)
          dtbetax_p_p(i,j,k) = dtbetax(i,j,k)
          dtbetay_p_p(i,j,k) = dtbetay(i,j,k)
          dtbetaz_p_p(i,j,k) = dtbetaz(i,j,k)
        enddo
      enddo
    enddo
  end if
      
  ! ==================== Fluid
  call CCTK_ActiveTimeLevelsVN (timelevels, cctkGH, "HydroBase::rho")
  if (timelevels > 1) then
    do k = 1, cctk_lsh(3)
      do j = 1, cctk_lsh(2)
        do i = 1, cctk_lsh(1)
          rho_p(i,j,k) = rho(i,j,k)
          eps_p(i,j,k) = eps(i,j,k)
          press_p(i,j,k) = press(i,j,k)
          vel_p(i,j,k,:) = vel(i,j,k,:)
          w_lorentz_p(i,j,k) = w_lorentz(i,j,k)
        enddo
      enddo
    enddo
  end if
  if (timelevels > 2) then
    do k = 1, cctk_lsh(3)
      do j = 1, cctk_lsh(2)
        do i = 1, cctk_lsh(1)
          rho_p_p(i,j,k) = rho(i,j,k)
          eps_p_p(i,j,k) = eps(i,j,k)
          press_p_p(i,j,k) = press(i,j,k)
          vel_p_p(i,j,k,:) = vel(i,j,k,:)
          w_lorentz_p(i,j,k) = w_lorentz(i,j,k)
        enddo
      enddo
    enddo
  end if

end subroutine ICPertFLRW_GRH_TimelevelCopy
