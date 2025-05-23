#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine ICPertFLRW_GRH_DeltaT (CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  integer, parameter :: dp = 8
  real(dp), parameter :: pi = 4._dp*atan(1._dp)
  CCTK_REAL :: H0, a0, tau0_EdS, eoswfac, cf_omega
  CCTK_REAL :: tau, aa, aa0, redshift
  logical :: want_pflrw_spacetime, want_proper_time

  character :: msg*1000

  want_pflrw_spacetime = CCTK_EQUALS (initial_data, "ICPertFLRW_GRH")
  want_proper_time = CCTK_EQUALS (ICPertFLRW_GRH_time, "proper")

  if (ICPertFLRW_GRH_need_to_fix_dtfac .EQ. 1) then
    ! This condition is so this calculation is
    ! not repeated when using mesh refinement
    if (want_proper_time) then
      if (want_pflrw_spacetime) then

        ! Define the pressure to energy density ratio
        cf_omega = 0._dp
        eoswfac = 3._dp * (1._dp + cf_omega)

        ! Cosmological measurements today
        H0 = ICPertFLRW_GRH_h * ICPertFLRW_GRH_c / 2997.9_dp ! Units are Mpc
        tau0_EdS = 2._dp / ( eoswfac * H0 )

        ! Define scale factor in terms of the reference redshift
        if (ICPertFLRW_GRH_z_comoving_ref .LT. 0._dp) then
          aa = 1._dp
        else
          a0 = 1._dp + ICPertFLRW_GRH_z_comoving_ref
        endif

        tau = cctk_time
        aa0 = ( tau / tau0_EdS )**(2._dp/eoswfac)
        redshift = (1._dp / aa0) - 1._dp

        if (ICPertFLRW_GRH_z_comoving_ref .LT. 0._dp) then
          a0 = aa * (1._dp + redshift)
        else
          aa = a0 / (1._dp + redshift)
        endif

        ! This overwrites the Time::dtfac parameter
        call CCTK_INFO("Define time spacing: dtfac = ICPertFLRW_GRH_time_fac * scale factor")
        dtfac = ICPertFLRW_GRH_time_fac * aa
        write (msg, '("The initial scale factor is ",g16.6)') aa
        call CCTK_INFO (msg)
        write (msg, '("dtfac is ",g16.6)') dtfac
        call CCTK_INFO (msg)

        ICPertFLRW_GRH_need_to_fix_dtfac = 0

      endif  !if you want pflrw spacetime
    endif  !if you want proper time
  endif !if ICPertFLRW_GRH_need_to_fix_dtfac
end subroutine ICPertFLRW_GRH_DeltaT




