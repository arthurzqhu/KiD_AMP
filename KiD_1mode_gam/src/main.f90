! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Driver for 1D Kinematic Driver (KiD) model
!
! Author: Ben Shipway
!
! For version details see header_data.f90
!

Module main

  Use typeKind
  Use parameters, only : dt, dg_dt, nx, nz, num_h_moments, h_shape
  Use namelists, only : read_namelist, l_ppe, n_ppe
  Use runtime, only : time, time_step, n_times
  Use switches
  Use set_profiles,  only : read_profiles
  Use interpolation, only : interpolate_input, interpolate_forcing
  Use diagnostics,   only : save_diagnostics_1d, save_diagnostics_2d, &
       write_diagnostics, query_dgstep
  Use derived_fields, only : calc_derived_fields
  Use advection_interface, only : advect_column
  Use mphys_interface, only : mphys_column
  Use stepfields, only : step_column
  Use divergence, only : diverge_column
  Use micro_prm, only: check_bintype,set_constants,nkr,npm,l_sctab_mod,diag_dt3,diag_dt4
  Use column_variables
  Use ppe_fun, only: load_latinhc
  Implicit none

contains

   subroutine main_loop
      real(8), allocatable :: temp_field(:,:)
      integer :: itime      ! loop counter for time
      real(8) :: t1, t2
      character(len=100) :: sigfig_fmt

      ! Start by reading in namelists
      allocate(temp_field(nz,0:nx+1))
      if (l_namelists) call read_namelist

      call check_bintype
      call set_constants

      ! the initial condition will be sampled from LH regardless of s_sample_dist
      if (l_ppe) call load_latinhc(n_ppe)

      ! Set up the initial fields and forcing
      if (l_input_file)then
         call read_profiles(input_file)
      else
         call read_profiles(icase)
      end if

      call interpolate_input(ifiletype)

      if (icase .ne. 501) then
         call interpolate_forcing
      endif

      call calc_derived_fields

      ! Do we want to do diagnostics on this timestep?
      call query_dgstep

      if ( nx == 1 ) then
         call save_diagnostics_1d
      else
         call save_diagnostics_2d
      endif

      call CPU_TIME(t1)
      diag_dt3 = 0.
      diag_dt4 = 0.

      do itime=1,n_times
         time=time+dt
         time_step=time_step+1

         ! Do we want to do diagnostics on this timestep?
         call query_dgstep

         if (icase .ne. 501) then
            call interpolate_forcing
         endif

         call calc_derived_fields

         if (l_advect)then
            call advect_column(scheme_id=0)
         end if

         if (l_diverge)then
            call diverge_column
         end if

         if (l_mphys)then
            call mphys_column(scheme_id=imphys)
         end if

         call step_column

         if ( nx == 1 ) then
            call save_diagnostics_1d
         else
            call save_diagnostics_2d
         endif
      end do
      call CPU_TIME(t2)
      print*, 'main loop time', t2-t1
      print*, 'invert_moments time', diag_dt3, ' s'
      ! print*, 'ml_moment2state time', diag_dt4, ' s'

      if (l_write_dgs) call write_diagnostics

      ! ! write lutable files for single category AMP:
      ! if (l_sctab_mod) call write_sctab

   end subroutine main_loop

End Module main
