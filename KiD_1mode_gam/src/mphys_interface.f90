! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module to interface with choice of microphysics schemes
!

module mphys_interface

  Use mphys_thompson07, only: mphys_thompson07_interface
  Use mphys_thompson09, only: mphys_thompson09_interface
  Use mphys_morr_two_moment, only: mphys_morrison_interface
#if UM_MICRO ==1
  Use mphys_um7_3, only: mphys_um7_3_interface
#endif
#if SHIPWAY_MICRO == 1
  Use mphys_4A, only: mphys_4A_interface
#endif
  Use mphys_amp, only:mphys_amp_interface
  Use mphys_tau_bin, only:mphys_tau_bin_interface
  Use switches

contains

  subroutine mphys_column(scheme_id)
    
    integer, intent(in) :: scheme_id

    select case (scheme_id)
   case(imphys_thompson09) ! Greg Thompson's mphys scheme
       call mphys_thompson09_interface
    case(imphys_thompson07) ! Greg Thompson's mphys scheme
       call mphys_thompson07_interface
    case(imphys_morr_two_moment) ! Hugh Morrisons's mphys scheme
       call mphys_morrison_interface
#if UM_MICRO == 1

    case(imphys_um7_3)      ! mphys scheme from um version 7.3
       call mphys_um7_3_interface
!print*, "test1"
#endif

    case(imphys_tau_bin)    ! TAU bin mphys scheme
       call mphys_tau_bin_interface
!print*, "test2"
#if SHIPWAY_MICRO == 1
    case(imphys_4A)         ! Shipway 4A scheme
       call mphys_4A_interface
!print*, "test3"
#endif
    case(imphys_amp)
       call mphys_amp_interface
!print*, "test4"
    end select

  end subroutine mphys_column

end module mphys_interface

