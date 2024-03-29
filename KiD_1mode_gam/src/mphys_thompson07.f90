! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module containing interface to greg thompson's microphysics (2007) code
!
module mphys_thompson07

  Use parameters, only : num_h_moments, num_h_bins, nspecies, nz, dt &
       , h_names, mom_units, max_char_len, nx
  Use column_variables
  Use physconst, only : p0, r_on_cp, pi

  Use module_mp_thompson07
  Use diagnostics, only: save_dg, i_dgtime

  Implicit None

  !Logical switches 
  logical :: micro_unset=.True.
  integer:: ih, imom
  character(max_char_len) :: name, units

contains

  Subroutine mphys_thompson07_interface

    real :: pptrain, pptsnow, pptgraul, pptice
    real :: pptrain_2d(nx), pptsnow_2d(nx), pptgraul_2d(nx) &
         , pptice_2d(nx)    
    real :: t1d(nz), p1d(nz), dz1d(nz),qv1d(nz),qc1d(nz) &
         ,qr1d(nz), qi1d(nz), ni1d(nz), qs1d(nz), qg1d(nz)

    integer :: kts, kte, i, j, k

    kts=1
    kte=nz
    j=1

    do i=1,nx
       pptrain = 0.
       pptsnow = 0.
       pptgraul = 0.
       pptice = 0.
       do k=1,nz
          t1d(k) = (theta(k,i) + (dtheta_adv(k,i)+dtheta_div(k,i))*dt )*exner(k,i)
          p1d(k) = p0*exner(k,i)**(1./r_on_cp)
          dz1d(k) = dz(k)
          qv1d(k) = qv(k,i) + (dqv_adv(k,i)+dqv_div(k,i))*dt
          
          qc1d(k) = hydrometeors(k,i,1)%moments(1,1) & 
               + (dhydrometeors_adv(k,i,1)%moments(1,1) &
               + dhydrometeors_div(k,i,1)%moments(1,1))*dt 
          
          qr1d(k) = hydrometeors(k,i,2)%moments(1,1) & 
               + (dhydrometeors_adv(k,i,2)%moments(1,1) &
               + dhydrometeors_div(k,i,2)%moments(1,1))*dt 
          
!          qi1d(k) = hydrometeors(k,i,3)%moments(1,1) & 
!               + (dhydrometeors_adv(k,i,3)%moments(1,1) &
!               + dhydrometeors_div(k,i,3)%moments(1,1))*dt 
!          
!          ni1d(k) = hydrometeors(k,i,3)%moments(1,2) & 
!               + (dhydrometeors_adv(k,i,3)%moments(1,2) &
!               + dhydrometeors_div(k,i,3)%moments(1,2))*dt 
!          
!          qs1d(k) = hydrometeors(k,i,4)%moments(1,1) & 
!               + (dhydrometeors_adv(k,i,4)%moments(1,1) &
!               + dhydrometeors_div(k,i,4)%moments(1,1))*dt 
!          
!          qg1d(k) = hydrometeors(k,i,5)%moments(1,1) & 
!               + (dhydrometeors_adv(k,i,5)%moments(1,1) &
!               + dhydrometeors_div(k,i,5)%moments(1,1))*dt 
          
       end do

       ! Initialise microphysics 
       if (micro_unset)then
          call thompson07_init
          micro_unset=.False.
       end if
       
       
       call mp_thompson(qv1d, qc1d, qi1d, qr1d, qs1d, qg1d, ni1d, &
            t1d, p1d, dz1d,                                       &
            pptrain, pptsnow, pptgraul, pptice,                   &
            kts, kte, dt, i, j)
       

       if (nx == 1) then
          ! Save some diagnostics
          imom=1
          !rain ppt
          ih=2
          name='surface_ppt_for_'//trim(h_names(ih))
          units=trim(mom_units(imom))//' m'
          call save_dg(pptrain, name, i_dgtime,  units, dim='time')
          !ice ppt
          ih=3
          name='surface_ppt_for_'//trim(h_names(ih))
          units=trim(mom_units(imom))//' m'
          call save_dg(pptice, name, i_dgtime,  units, dim='time')
          !snow ppt
          ih=4
          name='surface_ppt_for_'//trim(h_names(ih))
          units=trim(mom_units(imom))//' m'
          call save_dg(pptsnow, name, i_dgtime,  units, dim='time')
          !graupel ppt
          ih=5
          name='surface_ppt_for_'//trim(h_names(ih))
          units=trim(mom_units(imom))//' m'
          call save_dg(pptgraul, name, i_dgtime,  units, dim='time')
          !total ppt
          name='total_surface_ppt'
          units=trim(mom_units(imom))//' m'
          call save_dg((pptice+pptrain+pptsnow+pptgraul), name, i_dgtime, &
               units, dim='time') 
       else
          pptrain_2d(i) = pptrain
          pptice_2d(i) = pptice
          pptsnow_2d(i) = pptsnow
          pptgraul_2d(i) = pptgraul
       endif

       ! back out tendencies
       do k=1,nz
          dtheta_mphys(k,i)=(t1d(k)/exner(k,i)-theta(k,i))/dt & 
               - ( dtheta_adv(k,i)+dtheta_div(k,i))
          
          dqv_mphys(k,i)=(qv1d(k)-qv(k,i))/dt & 
               - ( dqv_adv(k,i)+dqv_div(k,i))
          
          dhydrometeors_mphys(k,i,1)%moments(1,1)= &
               (qc1d(k)-hydrometeors(k,i,1)%moments(1,1))/dt & 
               - (dhydrometeors_adv(k,i,1)%moments(1,1)  &
               + dhydrometeors_div(k,i,1)%moments(1,1))
          
          dhydrometeors_mphys(k,i,2)%moments(1,1)= &
               (qr1d(k)-hydrometeors(k,i,2)%moments(1,1))/dt & 
               - (dhydrometeors_adv(k,i,2)%moments(1,1)  &
               + dhydrometeors_div(k,i,2)%moments(1,1))

!          dhydrometeors_mphys(k,i,3)%moments(1,1)= &
!               (qi1d(k)-hydrometeors(k,i,3)%moments(1,1))/dt & 
!               - (dhydrometeors_adv(k,i,3)%moments(1,1)  &
!               + dhydrometeors_div(k,i,3)%moments(1,1))
!
!          dhydrometeors_mphys(k,i,3)%moments(1,2)= &
!               (ni1d(k)-hydrometeors(k,i,3)%moments(1,2))/dt & 
!               - (dhydrometeors_adv(k,i,3)%moments(1,2)  &
!               + dhydrometeors_div(k,i,3)%moments(1,2))
!          
!          dhydrometeors_mphys(k,i,4)%moments(1,1)= &
!               (qs1d(k)-hydrometeors(k,i,4)%moments(1,1))/dt & 
!               - (dhydrometeors_adv(k,i,4)%moments(1,1)  &
!               + dhydrometeors_div(k,i,4)%moments(1,1))
!          
!          dhydrometeors_mphys(k,i,5)%moments(1,1)= &
!               (qg1d(k)-hydrometeors(k,i,5)%moments(1,1))/dt & 
!               - (dhydrometeors_adv(k,i,5)%moments(1,1)  &
!               + dhydrometeors_div(k,i,5)%moments(1,1))
       end do

    end do

    if (nx > 1) then
       ! Save some mean scalar diagnostics
       imom=1
       !rain ppt
       ih=2
       name='surface_ppt_for_'//trim(h_names(ih))
       units=trim(mom_units(imom))//' m'
       call save_dg(pptrain_2d/nx, name, i_dgtime,  units, dim='time')
       !ice ppt
       ih=3
       name='surface_ppt_for_'//trim(h_names(ih))
       units=trim(mom_units(imom))//' m'
       call save_dg(pptice_2d/nx, name, i_dgtime,  units, dim='time')
       !snow ppt
       ih=4
       name='surface_ppt_for_'//trim(h_names(ih))
       units=trim(mom_units(imom))//' m'
       call save_dg(pptsnow_2d/nx, name, i_dgtime,  units, dim='time')
       !graupel ppt
       ih=5
       name='surface_ppt_for_'//trim(h_names(ih))
       units=trim(mom_units(imom))//' m'
       call save_dg(pptgraul_2d/nx, name, i_dgtime,  units, dim='time')
       !total ppt
       name='total_surface_ppt'
       units=trim(mom_units(imom))//' m'
       call save_dg((pptice_2d+pptrain_2d+pptsnow_2d+pptgraul_2d)/nx, &
            name, i_dgtime, units, dim='time') 
       
       ! save the all horizontal columns
       imom=1
       !rain ppt
       ih=2
       name='surface_ppt_for_'//trim(h_names(ih))
       units=trim(mom_units(imom))//' m'
       call save_dg(pptrain_2d, name, i_dgtime,  units, dim='time')
       !ice ppt
       ih=3
       name='surface_ppt_for_'//trim(h_names(ih))
       units=trim(mom_units(imom))//' m'
       call save_dg(pptice_2d, name, i_dgtime,  units, dim='time')
       !snow ppt
       ih=4
       name='surface_ppt_for_'//trim(h_names(ih))
       units=trim(mom_units(imom))//' m'
       call save_dg(pptsnow_2d, name, i_dgtime,  units, dim='time')
       !graupel ppt
       ih=5
       name='surface_ppt_for_'//trim(h_names(ih))
       units=trim(mom_units(imom))//' m'
       call save_dg(pptgraul_2d, name, i_dgtime,  units, dim='time')
       !total ppt
       name='total_surface_ppt'
       units=trim(mom_units(imom))//' m'
       call save_dg((pptice_2d+pptrain_2d+pptsnow_2d+pptgraul_2d), &
            name, i_dgtime, units, dim='time')        

    endif

  end Subroutine mphys_thompson07_interface

end module mphys_thompson07
