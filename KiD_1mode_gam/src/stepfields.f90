! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Increment the fields
!
!
module stepfields

  Use runtime
  Use typeKind

  Implicit None

  integer:: k, j

 contains

   subroutine step_column
     Use parameters, only : dt, nspecies, &
          num_h_moments, num_h_bins, num_aero_moments, num_aero_bins, nz, &
          nx
     Use physconst, only : Ru => R, p0, r_on_cp
     Use switches
     Use column_variables
     Use namelists, only : ampORbin, bintype
     Use micro_prm, only : pmomsc, diams, nkr, D_min, D_max, npmc, tempvar_debug3


     !local variables
     integer :: ih, imom, ibin
     real(wp) :: pos_tmp, srf_tmp(0:nx+1) ! temporary storage
     logical :: sp_test ! test for rounding error
     character(100) :: fmt ! format string
     real :: field(nz,0:nx+1)
     real(wp) :: moment_tentative, M_prev, M_curr, DIAM_CHAR(num_h_moments(1)), &
        mom3, mom0, momx, momy, DIAM_CHAR0, DIAM_CHARX, DIAM_CHARY
     integer :: mom_number, mom_previous_number, imomx, imomy
     logical :: l_print

     !------------------------------
     ! Manipulate mask if necessary
     !------------------------------
     if (.not. l_periodic_bound .and. nx>1)then
       field_mask(:,0)=0.
       field_mask(:,1)=0.
       field_mask(:,nx)=0.
       field_mask(:,nx+1)=0.
     end if
     

     !-----
     !theta
     !-----
     if (.not. l_fix_theta)then
        srf_tmp(:)=theta(1,:)
        ! Advective part
        if (.not. l_noadv_theta) then
           theta(:,:)=theta(:,:)+dtheta_adv(:,:)*field_mask(:,:)*dt
        endif
        ! mphys part
        if (.not. l_nomphys_theta) then
           theta(:,:)=theta(:,:)+dtheta_mphys(:,:)*field_mask(:,:)*dt
        endif
        ! divergence part
        theta(:,:)=theta(:,:)+dtheta_div(:,:)*field_mask(:,:)*dt
        ! forced part
        theta(:,:)=theta(:,:)+Tforce(:,:)*exner(:,:)*field_mask(:,:)*dt
        ! surface value
        if (isurface_fixed==isurface_fixed)theta(1,:)=srf_tmp(:)
     end if
     !-----
     !qv
     !-----
     if (.not. l_fix_qv)then
        srf_tmp(:)=qv(1,:)
        ! Advective part
        if (l_advect) then
        if (.not. l_noadv_qv)then
           if (maxval(dqv_adv(:,:))/=0. .or. minval(dqv_adv(:,:))/=0.) then
           end if
           if (l_posadv_qv)then
              qv(:,:)=qv(:,:)+max(0.,dqv_adv(:,:))*field_mask(:,:)*dt
           else
              qv(:,:)=qv(:,:)+dqv_adv(:,:)*field_mask(:,:)*dt
           end if
        end if
        endif
        ! mphys part
        if (.not. l_nomphys_qv) then
           qv(:,:)=qv(:,:)+dqv_mphys(:,:)*field_mask(:,:)*dt
        endif
        if (l_diverge) then
        ! divergence part
        qv(:,:)=qv(:,:)+dqv_div(:,:)*field_mask(:,:)*dt
        ! forced part
        qv(:,:)=qv(:,:)+qforce(:,:)*field_mask(:,:)*dt
        endif
        ! surface value
        if (isurface_fixed==isurface_fixed)qv(1,:)=srf_tmp(:)

     end if
     !------------
     !aerosol
     !------------
     ! Advective part
     if (.not. l_fix_aerosols)then
        do ih=1,naerosol
           do imom=1,num_aero_moments(ih)
              do ibin=1,num_aero_bins(ih)
                 do j=1,nx
                 do k=1,nz
                    aerosol(k,j,ih)%moments(ibin,imom)=           &
                         aerosol(k,j,ih)%moments(ibin,imom) + (   &
                         daerosol_adv(k,j,ih)%moments(ibin,imom)  &
                         )*field_mask(k,j)*dt
                 end do
                 end do
              end do
           end do
        end do
     end if
     ! mphys part
     do ih=1,naerosol
        do imom=1,num_aero_moments(ih)
           do ibin=1,num_aero_bins(ih)
              do j=1,nx
              do k=1,nz
                 aerosol(k,j,ih)%moments(ibin,imom)=              &
                      aerosol(k,j,ih)%moments(ibin,imom) + (      &
                      + daerosol_mphys(k,j,ih)%moments(ibin,imom) &
                      )*field_mask(k,j)*dt
              end do
              end do
           end do
        end do
     end do
     ! divergence part
     if (.not. l_fix_aerosols)then
        do ih=1,naerosol
           do imom=1,num_aero_moments(ih)
              do ibin=1,num_aero_bins(ih)
                 do j=1,nx
                 do k=1,nz
                    aerosol(k,j,ih)%moments(ibin,imom)=            &
                         aerosol(k,j,ih)%moments(ibin,imom) + (    &
                         + daerosol_div(k,j,ih)%moments(ibin,imom) &
                         )*field_mask(k,j)*dt
                 end do
                 end do
              end do
           end do
        end do
     end if
     !------------
     !hydrometeors
     !------------
     ! Advective part

     l_print = .false.
     if (.not. l_noadv_hydrometeors)then
        do j=1,nx
        do k=1,nz
           do ih=1,nspecies
              do ibin=1,num_h_bins(ih)
                 do imom=1,num_h_moments(ih)
                    hydrometeors(k,j,ih)%moments(ibin,imom) = &
                         hydrometeors(k,j,ih)%moments(ibin,imom) + (   &
                         dhydrometeors_adv(k,j,ih)%moments(ibin,imom)  &
                         )*field_mask(k,j)*dt
                 end do
                 end do
              end do
           end do
        end do

                 ! ! constraining the moment ratios 0 -> <largest> -> <intermediate>

                 ! if (ampORbin .eq. 'amp' .and. npmc == 4 &
                 !    .and. all(hydrometeors(k,j,ih)%moments(ibin,1:4)>0.)) then
                 !    mom3 = hydrometeors(k,j,ih)%moments(ibin,1)
                 !    mom0 = hydrometeors(k,j,ih)%moments(ibin,2)
                 !    momx = hydrometeors(k,j,ih)%moments(ibin,3)
                 !    momy = hydrometeors(k,j,ih)%moments(ibin,4)
                 !    imomx = pmomsc(3)
                 !    imomy = pmomsc(4)


                 !    ! review this part
                 !    if (pmomsc(3) > 3) then
                 !       DIAM_CHAR0 = (mom3/mom0)**(1./3.)
                 !       DIAM_CHARX = (momx/mom3)**(1./(imomx-3.))
                 !       DIAM_CHARY = (momy/momx)**(1./(imomy-imomx))
                 !    elseif (pmomsc(3) < 3 .and. pmomsc(4) > 3) then
                 !       DIAM_CHAR0 = (momx/mom0)**(1./imomx)
                 !       DIAM_CHARX = (mom3/momx)**(1./(3.-imomx))
                 !       DIAM_CHARY = (momy/mom3)**(1./(imomy-3.))
                 !    elseif (pmomsc(4) < 3) then
                 !       DIAM_CHAR0 = (momx/mom0)**(1./imomx)
                 !       DIAM_CHARX = (momy/momx)**(1./(imomy-imomx))
                 !       DIAM_CHARY = (mom3/momy)**(1./(3.-imomy))
                 !    endif

                 !    if (DIAM_CHAR0 .ne. DIAM_CHAR0) then
                 !       print*, 'nan in DIAM_CHAR0 during adv'
                 !       print*, k,ih,ibin
                 !       print*, 'moms', hydrometeors(k,j,ih)%moments(ibin,:)
                 !       stop
                 !    endif

                 !    if (DIAM_CHARX .ne. DIAM_CHARX) then
                 !       print*, 'nan in DIAM_CHARX during adv'
                 !       print*, k,ih,ibin
                 !       print*, 'moms', hydrometeors(k,j,ih)%moments(ibin,:)
                 !       stop
                 !    endif

                 !    if (DIAM_CHARY .ne. DIAM_CHARY) then
                 !       print*, 'nan in DIAM_CHARY during adv'
                 !       print*, k,ih,ibin
                 !       print*, 'moms', hydrometeors(k,j,ih)%moments(ibin,:)
                 !       print*, 'from adv interface:', tempvar_debug3
                 !       stop
                 !    endif


                 !    if (DIAM_CHAR0 .ne. DIAM_CHAR0) then
                 !       print*, 'mom3, mom0', mom3, mom0, momx, momy
                 !       print*, 'hyd0, dhyd0', k, hydrometeors(k,j,ih)%moments(ibin,2), &
                 !          dhydrometeors_adv(k,j,ih)%moments(ibin,2)
                 !       print*, 'from adv intf:', tempvar_debug3
                 !       stop
                 !    endif

                 !    if (DIAM_CHAR0 < D_min) then
                 !       DIAM_CHAR0 = D_min
                 !       print*, 'DIAM_CHAR0 < D_min', DIAM_CHAR0, D_min
                 !    endif
                 !    if (DIAM_CHAR0 > D_max*0.997) then
                 !       DIAM_CHAR0 = D_max*0.997
                 !       print*, 'DIAM_CHAR0 > D_max', DIAM_CHAR0, D_max
                 !    endif
                 !    if (DIAM_CHARX < DIAM_CHAR0*1.001) then
                 !       DIAM_CHARX = DIAM_CHAR0*1.001
                 !       print*, 'DIAM_CHARX < DIAM_CHAR0', DIAM_CHARX, DIAM_CHAR0
                 !    endif
                 !    if (DIAM_CHARX > D_max*0.998) then
                 !       DIAM_CHARX = D_max*0.998
                 !       print*, 'DIAM_CHARX > D_max', DIAM_CHARX, D_max
                 !    endif
                 !    if (DIAM_CHARY < DIAM_CHARX*1.001) then
                 !       DIAM_CHARY = DIAM_CHARX*1.001
                 !       print*, 'DIAM_CHARY < DIAM_CHARX', DIAM_CHARY, DIAM_CHARX
                 !    endif
                 !    if (DIAM_CHARY > D_max) then
                 !       DIAM_CHARY = D_max*0.999
                 !       print*, 'DIAM_CHARY > D_max', DIAM_CHARY, D_max
                 !    endif
                 !    ! if (DIAM_CHAR0 .eq. DIAM_CHAR0) then
                 !    !    print*, DIAM_CHAR0, DIAM_CHARX, DIAM_CHARY
                 !    !    stop
                 !    ! endif

                 !    if (pmomsc(3) > 3) then
                 !       mom0 = mom3/(DIAM_CHAR0**3.)
                 !       momx = DIAM_CHARX**(imomx-3.)*mom3
                 !       momy = DIAM_CHARY**(imomy-imomx)*momx
                 !    elseif (pmomsc(3) < 3 .and. pmomsc(4) > 3) then
                 !       momx = DIAM_CHARX**(imomx-3.)*mom3
                 !       mom0 = momx/(DIAM_CHAR0**imomx)
                 !       momy = DIAM_CHARY**(imomy-3.)*mom3
                 !    elseif (pmomsc(4) < 3) then
                 !       momy = DIAM_CHARY**(imomy-3.)*mom3
                 !       momx = DIAM_CHARX**(imomx-imomy)*momy
                 !       mom0 = momx/(DIAM_CHAR0**imomx)
                 !    endif
                 !    hydrometeors(k,j,ih)%moments(ibin,2) = mom0
                 !    hydrometeors(k,j,ih)%moments(ibin,3) = momx
                 !    hydrometeors(k,j,ih)%moments(ibin,4) = momy
                 ! endif


                 ! if (any(hydrometeors(k,j,ih)%moments(ibin,:)<0)) then
                 !    print*, k,j,ih
                 !    print*, hydrometeors(k,j,ih)%moments(ibin,:)
                 !    stop 'adv moment<0'
                 ! endif

     end if
     ! if (l_print) stop

     ! mphys part
     do ih=1,nspecies
        do imom=1,num_h_moments(ih)
           do ibin=1,num_h_bins(ih)
              do j=1,nx
                 do k=1,nz
                    hydrometeors(k,j,ih)%moments(ibin,imom)=              &
                       hydrometeors(k,j,ih)%moments(ibin,imom) + (      &
                       dhydrometeors_mphys(k,j,ih)%moments(ibin,imom)   &
                       )*field_mask(k,j)*dt
                 end do
              end do
           end do
        end do
     end do

     ! divergence part
     if (.not. l_nodiv_hydrometeors)then
        do ih=1,nspecies
           do imom=1,num_h_moments(ih)
              do ibin=1,num_h_bins(ih)
                 do j=1,nx
                 do k=1,nz
                    hydrometeors(k,j,ih)%moments(ibin,imom)=            &
                         hydrometeors(k,j,ih)%moments(ibin,imom) + (    &
                         dhydrometeors_div(k,j,ih)%moments(ibin,imom)   &
                         )*field_mask(k,j)*dt
                 end do
                 end do
              end do
           end do
        end do
     end if
     !      ! forced part - not yet implemented
     !         do ih=1,nspecies
     !            do imom=1,num_h_moments(ih)
     !               do ibin=1,num_h_bins(ih)
     !                  hydrometeors(k,ih)%moments(ibin,imom)=              &
     !                       hydrometeors(k,ih)%moments(ibin,imom) + (      &
     !                       + dhydrometeors_force(k,ih)%moments(ibin,imom) &
     !                       )*field_mask(k,j)*dt
     !               end do
     !            end do
     !         end do

     ! Positivity check - to prevent rounding errors giving negative
     ! values - !!! NB this violates conservation !!!
     if (l_force_positive)then
     do ih=1,nspecies
        do imom=1,num_h_moments(ih)
           do ibin=1,num_h_bins(ih)
              do j=1,nx
              do k=1,nz
                 pos_tmp=max(0., hydrometeors(k,j,ih)%moments(ibin&
                      &,imom))
                 if (hydrometeors(k,j,ih)%moments(ibin,imom)<0.) then
                    sp_test=1.e3*SPACING(MAX( &
                         dt*dhydrometeors_mphys(k,j,ih)%moments(ibin,imom) &
                         ,dt*dhydrometeors_div(k,j,ih)%moments(ibin,imom) &
                         ,dt*dhydrometeors_adv(k,j,ih)%moments(ibin,imom) &
                         ,hydrometeors(k,j,ih)%moments(ibin,imom) &
                         - dt*(                                   &
                         dhydrometeors_mphys(k,j,ih)%moments(ibin,imom) &
                         +dhydrometeors_div(k,j,ih)%moments(ibin,imom) &
                         +dhydrometeors_adv(k,j,ih)%moments(ibin,imom)))) &
                         < ABS(hydrometeors(k,j,ih)%moments(ibin,imom))
                    if (sp_test)then
                       fmt='( A, /, A, T33, E11.3, /, &
                            &  A, T33, E11.3, /, &
                            &  A, T18, I2,  /, &
                            &  A, T18, I2,  /, &
                            &  A, T18, I2 )'
                       write(*,fmt)  'Warning some &
                            &negative numbers have been generated which do&
                            & not look like rounding error.', 'Estimat&
                            &ed rounding error: ', 10.*SPACING(MAX( &
                            dt*dhydrometeors_mphys(k,j,ih)%moments(ibin,imom) &
                            ,dt*dhydrometeors_div(k,j,ih)%moments(ibin,imom) &
                            ,dt*dhydrometeors_adv(k,j,ih)%moments(ibin,imom) &
                            ,hydrometeors(k,j,ih)%moments(ibin,imom) &
                            - dt*(                                   &
                            dhydrometeors_mphys(k,j,ih)%moments(ibin,imom) &
                            +dhydrometeors_div(k,j,ih)%moments(ibin,imom) &
                            +dhydrometeors_adv(k,j,ih)%moments(ibin&
                            &,imom)))),   &
                            'Negative value of hydrometeor: ', &
                            hydrometeors(k,j,ih)%moments(ibin,imom), &
                            'species: ',ih, &
                            'moment: ',imom, &
                            'bin: ',ibin
                       write(*, '(''Timestep: '', I4, '' z level: '', &
                            &I4)') time_step, k
                       write(*,*) ''
                       write(*,*)   &
                             dt*dhydrometeors_mphys(k,j,ih)%moments(ibin,imom) &
                             ,dt*dhydrometeors_div(k,j,ih)%moments(ibin,imom) &
                             ,dt*dhydrometeors_adv(k,j,ih)%moments(ibin,imom) &
                             ,hydrometeors(k,j,ih)%moments(ibin,imom) &
                             - dt*(                                   &
                             dhydrometeors_mphys(k,j,ih)%moments(ibin,imom) &
                             +dhydrometeors_div(k,j,ih)%moments(ibin,imom) &
                             +dhydrometeors_adv(k,j,ih)%moments(ibin&
                             &,imom))
                    end if
                 end if
                 hydrometeors(k,j,ih)%moments(ibin,imom)=pos_tmp
              end do
              end do
           end do
        end do
     end do
     end if
     !----------------
     ! update pressure
     !----------------
     if (l_pupdate) then
        do j=1,nx
           exner(:,j)=(rho(:)*Ru*theta(:,j)/p0)**(r_on_cp/(1.-r_on_cp))
        end do
     endif

   end subroutine step_column

 end module stepfields
