! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module containing interface to greg thompson's microphysics (2007) code
!
module mphys_amp

  Use parameters, only : num_h_moments, num_h_bins, h_shape, nspecies, nz, dt &
       , h_names, mom_units, max_char_len, nx
  Use column_variables
  Use physconst, only : p0, r_on_cp, pi, rhow

  Use module_hujisbm
  Use micro_prm
  Use diagnostics, only: save_dg, i_dgtime, my_save_dg_bin_dp
  Use switches, only: l_advect,l_diverge, l_noadv_theta, l_noadv_qv
  Use namelists, only: bintype, ampORbin, l_noadv_hydrometeors
  use switches, only: zctrl
  Implicit None

  !Logical switches
  logical :: micro_unset=.True.
  integer:: ih, imom
  character(max_char_len) :: name, units

contains

  Subroutine mphys_amp_interface
    use parameters, only: flag_count,max_nbins
    integer :: i, j, k, imom, rain_alt
    real, dimension(nz,nx) :: t2d, p2d, qv2d
    real(8), dimension(nz,nx,num_h_moments(1)) :: Mpc2d
    real(8), dimension(nz,nx,num_h_moments(1)) :: Mpr2d
    real(8),dimension(nz,nx,10) :: mc,mr
    real(8), save, dimension(nz,nx,2) :: guessc2d,guessr2d
    real(8), dimension(nz,nx,max_nbins) :: aer2d,dropsm2d,dropsn2d,dropsinitm2d,dropsinitn2d
    real, dimension(nz) :: field
    real(8),dimension(nz,flag_count) :: fieldflag
    real(8),dimension(nz,nx) :: dm
    real(8), dimension(nz) :: fielddp
    real(8), dimension(nz) :: dm_c, dm_r, dm_w
    real(8), dimension(nz,nx) :: fielddp2d, dm_c2d, dm_r2d, dm_w2d
    real(8), dimension(nz,max_nbins) :: fieldbin
    real(8), dimension(nz,nx,max_nbins) :: fieldbin2d
    real(8), dimension(nz,nx,2,flag_count) :: flag
    real(8), dimension(num_h_moments(2)) :: mr_s
    real(8) :: inf=huge(fielddp(1))
    real(8), dimension(nz,nx) :: reff, m3w
    real(8), dimension(nx) :: opdep,albd
    real(8), dimension(nz) :: lwc
    real(8), dimension(nz) :: nmask ! mask for non-NaN values
    real(8) :: m2w, m0w
    character(1) :: Mnum
    real(8) :: dnr_s,rm_s !source dnr and rain mass

    do i=1,nx
       do k=1,nz
          p2d(k,i) = p0*exner(k,i)**(1./r_on_cp)

          t2d(k,i) = theta(k,i)*exner(k,i)
          qv2d(k,i) = qv(k,i)

          !if (.not. l_noadv_qv) qv2d(k,i) = qv(k,i) + dqv_adv(k,i)*dt             
          !if (.not. l_noadv_theta) t2d(k,i) = (theta(k,i) + dtheta_adv(k,i)*dt )*exner(k,i)
          !if (l_diverge) then
          !   qv2d(k,i) = qv(k,i) + dqv_div(k,i)*dt             
          !   t2d(k,i) = (theta(k,i) + dtheta_div(k,i)*dt )*exner(k,i)
          !endif

          if (ampORbin .eq. 'amp') then
             do imom=1,num_h_moments(1)
               Mpc2d(k,i,imom) = hydrometeors(k,i,1)%moments(1,imom)
               !if (.not. l_noadv_hydrometeors) Mpc2d(k,i,imom)=Mpc2d(k,i,imom) &
               !     + dhydrometeors_adv(k,i,1)%moments(1,imom)*dt
               !if (l_diverge) Mpc2d(k,i,imom)=Mpc2d(k,i,imom) &
               !     + dhydrometeors_div(k,i,1)%moments(1,imom)*dt
             enddo
             if (any(Mpc2d(k,i,:)==0.)) Mpc2d(k,i,:)=0.
             do imom=1,num_h_moments(2)
               Mpr2d(k,i,imom) = hydrometeors(k,i,2)%moments(1,imom)
               !if (.not. l_noadv_hydrometeors) Mpr2d(k,i,imom)=Mpr2d(k,i,imom) &
               !     + dhydrometeors_adv(k,i,2)%moments(1,imom)*dt
               !if (l_diverge) Mpr2d(k,i,imom)=Mpr2d(k,i,imom) &
               !     + dhydrometeors_div(k,i,2)%moments(1,imom)*dt 
             enddo
             if (any(Mpr2d(k,i,:)==0.)) Mpr2d(k,i,:)=0.
          else !bin
             Mpc2d=0.;Mpr2d=0.
             do j=1,nkr
                dropsm2d(k,i,j)=hydrometeors(k,i,1)%moments(j,1)/col
                !if (.not. l_noadv_hydrometeors) dropsm2d(k,i,j)=dropsm2d(k,i,j) &
                !   + dhydrometeors_adv(k,i,1)%moments(j,1)*dt/col
                !if (l_diverge) dropsm2d(k,i,j)=dropsm2d(k,i,j) &
                !   + dhydrometeors_div(k,i,1)%moments(j,1)*dt/col
                if (bintype .eq. 'tau')  then
                   dropsn2d(k,i,j)=hydrometeors(k,i,1)%moments(j,2)/col
                   !if (.not. l_noadv_hydrometeors) dropsn2d(k,i,j)=dropsn2d(k,i,j) &
                   !   + dhydrometeors_adv(k,i,1)%moments(j,2)*dt/col
                   !if (l_diverge) dropsn2d(k,i,j)=dropsn2d(k,i,j) &
                   !   + dhydrometeors_div(k,i,1)%moments(j,2)*dt/col
                end if
             enddo
          endif
       enddo
   enddo

   ! Initialise microphysics
   if (micro_unset)then
      if (ampORbin .eq. 'amp') then
         guessc2d(:,:,1) = h_shape(1) !shape parameter
         guessr2d(:,:,1) = h_shape(2)
         guessc2d(:,:,2) = 0.001         !characteristic diameter dn
         guessr2d(:,:,2) = 0.001
         call amp_init(aer2d,Mpc2d,Mpr2d,guessc2d,guessr2d)
      elseif (ampORbin .eq. 'bin') then
         if (bintype .eq. 'sbm') then
            call sbm_init(aer2d,dropsm2d)
         elseif (bintype .eq. 'tau') then
            call tau_init(aer2d,dropsm2d,dropsn2d)
         endif
      endif
      micro_unset=.False.
   endif

   aer2d = 0.

   ! set rain source if there is one
   if (rain_source(1)>0.) then 
       ! set the zctrl(1) to be where the rain source is, which happens to be nz-1
       ! cant be nz because tau does not sediment moisture from the topmost layer
       rain_alt=int(zctrl(1)/((zctrl(1))/nz))-1
       rm_s=rain_source(1)**3*pi/6*rhow*rain_source(2)
       dnr_s=(rm_s*6./3.14159/rhow/rain_source(2)*gamma(h_shape(2))/gamma(h_shape(2)+3))**(1./3.)

       do j=1,nx
           if (bintype .eq. 'sbm') then
               CALL init_dist_sbm(dble(0.),h_shape(1),dble(0.),rm_s,h_shape(2),dnr_s,&
                    diams,dropsm2d(rain_alt,j,:))
           elseif (bintype .eq. 'tau') then
               CALL init_dist_tau(dble(0.),h_shape(1),dble(0.),rm_s,h_shape(2),dnr_s,&
                    dropsm2d(rain_alt,j,:),dropsn2d(rain_alt,j,:))
           endif
       enddo

       if (ampORbin .eq. 'amp') then
           guessr2d(rain_alt,:,2)=dnr_s
       
           pmomsr=(/3,imomr1,imomr2/)
           do j=1,nx
               do i=1,num_h_moments(2)
                   if (bintype .eq. 'sbm') then
                       mr_s(i)=sum(dropsm2d(rain_alt,j,split_bins+1:nkr)/xl(split_bins+1:nkr)&
                             *diams(split_bins+1:nkr)**pmomsr(i))*col*rhow
                   elseif (bintype .eq. 'tau') then
                       mr_s(i)=sum(dropsm2d(rain_alt,j,split_bins+1:nkr)/binmass(split_bins+1:nkr)&
                             *diams(split_bins+1:nkr)**pmomsc(i))*col
           
                   end if
               enddo
               Mpr2d(rain_alt,j,1:num_h_moments(2))=mr_s(1:num_h_moments(2))
           enddo 
       endif
   endif

   if (ampORbin .eq. 'amp') then
      dropsm2d=0.
      dropsn2d=0.
      dropsinitm2d=0.
      dropsinitn2d=0.

      call mp_amp(Mpc2d,Mpr2d,guessc2d,guessr2d, &
           p2d,t2d,qv2d,aer2d,dropsm2d,dropsn2d,mc,&
           mr,flag,dropsinitm2d,dropsinitn2d)

      !dropsm2d=dropsinitm2d
      !dropsn2d=dropsinitn2d
      !print*, Mpc2d(46,1,:)
      !print*, Mpr2d(46,1,:)
   elseif (ampORbin .eq. 'bin') then
      if (bintype .eq. 'sbm') then
         call mp_sbm(dropsm2d,p2d,t2d,qv2d,aer2d,mc,mr)
         !print*, dropsm2d(46,:,:)
      elseif (bintype .eq. 'tau') then
         call mp_tau(dropsm2d,dropsn2d,t2d,qv2d,mc,mr)
      endif
   endif

   !if (ampORbin .eq. 'bin') then
   !   if (bintype .eq. 'sbm') num_h_moments=(/1,1/)
   !   if (bintype .eq. 'tau') num_h_moments=(/2,1/)
   !endif

  ! back out tendencies
dqv_mphys = 0.
dtheta_mphys = 0.

do imom=1,num_h_moments(1)
  do i=1,num_h_bins(1)
   dhydrometeors_mphys(:,:,1)%moments(i,imom)=0.
  enddo
enddo

do imom=1,num_h_moments(2)
  do i=1,num_h_bins(2)
   dhydrometeors_mphys(:,:,2)%moments(i,imom)=0.
  enddo
enddo

do i=1,nx
   do k=1,nz
      dtheta_mphys(k,i)=(t2d(k,i)/exner(k,i)-theta(k,i))/dt
      dqv_mphys(k,i)=(qv2d(k,i)-qv(k,i))/dt

      !if (.not. l_noadv_qv) dqv_mphys(k,i)=dqv_mphys(k,i)-dqv_adv(k,i)      
      !if (.not. l_noadv_theta) dtheta_mphys(k,i)=dtheta_mphys(k,i)-dtheta_adv(k,i)
      !if (l_diverge) then
      !   dtheta_mphys(k,i)=dtheta_mphys(k,i)-dtheta_div(k,i)
      !   dqv_mphys(k,i)=dqv_mphys(k,i)-dqv_div(k,i)
      !endif

      if (ampORbin .eq. 'amp') then ! when bulk
         do imom=1,num_h_moments(1)
            dhydrometeors_mphys(k,i,1)%moments(1,imom)= &
                 (Mpc2d(k,i,imom)-hydrometeors(k,i,1)%moments(1,imom))/dt
           !if (.not. l_noadv_hydrometeors) dhydrometeors_mphys(k,i,1)%moments(1,imom)= &
           !              dhydrometeors_mphys(k,i,1)%moments(1,imom) &
           !              -dhydrometeors_adv(k,i,1)%moments(1,imom)
           !if (l_diverge) dhydrometeors_mphys(k,i,1)%moments(1,imom)= &
           !               dhydrometeors_mphys(k,i,1)%moments(1,imom) &
           !               -dhydrometeors_div(k,i,1)%moments(1,imom)
         enddo
         do imom=1,num_h_moments(2)
            dhydrometeors_mphys(k,i,2)%moments(1,imom)= &
                 (Mpr2d(k,i,imom)-hydrometeors(k,i,2)%moments(1,imom))/dt
           !if (.not. l_noadv_hydrometeors) dhydrometeors_mphys(k,i,2)%moments(1,imom)= &
           !              dhydrometeors_mphys(k,i,2)%moments(1,imom) &
           !              -dhydrometeors_adv(k,i,2)%moments(1,imom)
           !if (l_diverge) dhydrometeors_mphys(k,i,2)%moments(1,imom)= &
           !               dhydrometeors_mphys(k,i,2)%moments(1,imom) &
           !               -dhydrometeors_div(k,i,2)%moments(1,imom)
         enddo

      else ! when bin
         do j=1,nkr
            dhydrometeors_mphys(k,i,1)%moments(j,1)= &
                   (dropsm2d(k,i,j)*col-hydrometeors(k,i,1)%moments(j,1))/dt
            !if (.not. l_noadv_hydrometeors) dhydrometeors_mphys(k,i,1)%moments(j,1)= &
            !              dhydrometeors_mphys(k,i,1)%moments(j,1) &
            !              -dhydrometeors_adv(k,i,1)%moments(j,1)
            !if (l_diverge) dhydrometeors_mphys(k,i,1)%moments(j,1)= &
            !               dhydrometeors_mphys(k,i,1)%moments(j,1) &
            !               -dhydrometeors_div(k,i,1)%moments(j,1)
            if (bintype .eq. 'tau') then
               dhydrometeors_mphys(k,i,1)%moments(j,2)= &
                     (dropsn2d(k,i,j)*col-hydrometeors(k,i,1)%moments(j,2))/dt
               !if (.not. l_noadv_hydrometeors) dhydrometeors_mphys(k,i,1)%moments(j,2)= &
               !              dhydrometeors_mphys(k,i,1)%moments(j,2) &
               !              -dhydrometeors_adv(k,i,1)%moments(j,2)
               !if (l_diverge) dhydrometeors_mphys(k,i,1)%moments(j,2)= &
               !               dhydrometeors_mphys(k,i,1)%moments(j,2) &
               !               -dhydrometeors_div(k,i,1)%moments(j,2)
            endif
         enddo

      endif

   enddo
enddo

! Save some diagnostics
!fitting flag
if (ampORbin .eq. 'amp') then
   if (nx==1) then
   field(:)=flag(:,nx,1,1)
   call save_dg(field,'oflagc', i_dgtime,units='unitless', dim='z')
   field(:)=flag(:,nx,2,1)
   call save_dg(field,'oflagr', i_dgtime,units='unitless', dim='z')

   field(:)=flag(:,nx,1,2)
   call save_dg(field,'flagoobc', i_dgtime,units='unitless', dim='z')
   field(:)=flag(:,nx,2,2)
   call save_dg(field,'flagoobr', i_dgtime,units='unitless', dim='z')

!   field(:)=flag(:,nx,1,3)
!   call save_dg(field,'fitting_flag_cloud_x_intol', i_dgtime,units='unitless', dim='z')
!   field(:)=flag(:,nx,2,3)
!   call save_dg(field,'fitting_flag_rain_x_intol', i_dgtime,units='unitless', dim='z')

!   field(:)=flag(:,nx,1,4)
!   call save_dg(field,'fitting_flag_cloud_y_intol', i_dgtime,units='unitless', dim='z')
!   field(:)=flag(:,nx,2,4)
!   call save_dg(field,'fitting_flag_rain_y_intol', i_dgtime,units='unitless', dim='z')
   else
   endif

endif

   ! if (imomc1.ne.3) then
   !       field(:)=flag(:,nx,1)
   !       call save_dg(field,'fitting_flag_cloud', i_dgtime,units='unitless', dim='z')
   !       field(:)=flag(:,nx,2)
   !       call save_dg(field,'fitting_flag_rain', i_dgtime,units='unitless', dim='z')
   ! endif

!diagnosed moments
do i=1,10
  write(Mnum,'(I1)') i-1

  name='diagM'//Mnum//'_cloud'
  units='m^'//Mnum
  if (nx==1) then
     fielddp(:)=mc(:,1,i)
     call save_dg(fielddp,name,i_dgtime,units,dim='z')
  else
     fielddp2d(:,:)=mc(:,:,i)
     call save_dg(fielddp2d,name,i_dgtime,units,dim='z,x')
  endif

  name='diagM'//Mnum//'_rain'
  if (nx==1) then
     fielddp(:)=mr(:,1,i)
     call save_dg(fielddp,name,i_dgtime,units,dim='z')
  else
     fielddp2d(:,:)=mr(:,:,i)
     call save_dg(fielddp2d,name,i_dgtime,units,dim='z,x')
  endif
enddo

! diagnose mass mean diameter and effective radius
! mc and mr here might not be additive because it's updated after diagnosing from
! the DSD after mphys -ahu

do j=1,nx
   do k=1,nz
      m3w(k,j)=mc(k,j,4)+mr(k,j,4)
      m2w=mc(k,j,3)+mr(k,j,3)
      m0w=mc(k,j,1)+mr(k,j,1)
      reff(k,j)=m3w(k,j)/m2w*0.5
      if (nx==1) then
         dm_w(k)=(m3w(k,j)/m0w)**(1./3.)
         dm_c(k)=(mc(k,j,4)/mc(k,j,1))**(1./3.)
         dm_r(k)=(mr(k,j,4)/mr(k,j,1))**(1./3.)
         if ( (dm_w(k) .ne. dm_w(k)) .or. (dm_w(k)>inf) ) dm_w(k)=0.
         if ( (dm_c(k) .ne. dm_c(k)) .or. (dm_c(k)>inf) ) dm_c(k)=0.
         if ( (dm_r(k) .ne. dm_r(k)) .or. (dm_r(k)>inf) ) dm_r(k)=0.
      else
         dm_w2d(k,j)=(m3w(k,j)/m0w)**(1./3.)
         dm_c2d(k,j)=(mc(k,j,4)/mc(k,j,1))**(1./3.)
         dm_r2d(k,j)=(mr(k,j,4)/mr(k,j,1))**(1./3.)
         if ( (dm_w2d(k,j) .ne. dm_w2d(k,j)) .or. (dm_w2d(k,j)>inf) ) dm_w2d(k,j)=0.
         if ( (dm_c2d(k,j) .ne. dm_c2d(k,j)) .or. (dm_c2d(k,j)>inf) ) dm_c2d(k,j)=0.
         if ( (dm_r2d(k,j) .ne. dm_r2d(k,j)) .or. (dm_r2d(k,j)>inf) ) dm_r2d(k,j)=0.
      endif
      if ( (reff(k,j) .ne. reff(k,j)) .or. (reff(k,j)>inf) ) reff(k,j)=0.
   enddo
enddo

if (nx==1) then
   call save_dg(dm_w,'Dm',i_dgtime,'micron',dim='z')
   call save_dg(dm_c,'Dm',i_dgtime,'micron',dim='z')
   call save_dg(dm_r,'Dm',i_dgtime,'micron',dim='z')
   call save_dg(reff(:,1)*1.e6,'reff',i_dgtime,'micron',dim='z')
else
   call save_dg(dm_w2d,'Dm',i_dgtime,'micron',dim='z,x')
   call save_dg(dm_c2d,'Dm',i_dgtime,'micron',dim='z,x')
   call save_dg(dm_r2d,'Dm',i_dgtime,'micron',dim='z,x')
   call save_dg(reff*1.e6,'reff',i_dgtime,'micron',dim='z,x')
endif

!diagnose optical depth and albedo
units='unitless'
opdep(:)=0.0
do j=1,nx
   do k=1,nz
      lwc(k)=m3w(k,j)*pi/6*rhow*rho(k)
      if (reff(k,j).ne.0. .and. reff(k,j)<inf) then
         opdep(j)=opdep(j)+3*lwc(k)/(2*rhow*reff(k,j))*dz(k)
      endif
   enddo
   albd(j)=opdep(j)/(opdep(j)+13.3)
enddo

if (nx==1) then
   call save_dg(opdep(1), 'opt_dep', i_dgtime,  units,dim='time')
   call save_dg(albd(1), 'albedo', i_dgtime,  units,dim='time')
else
   call save_dg(opdep(:), 'opt_dep', i_dgtime,  units,dim='x')
   call save_dg(albd(:), 'albedo', i_dgtime,  units,dim='x')
   call save_dg(sum(opdep(1:nx))/nx, 'mean_opt_dep', i_dgtime,  units,dim='time')
   call save_dg(sum(albd(1:nx))/nx, 'mean_albedo', i_dgtime,  units,dim='time')
endif

!bin distributions
if (ampORbin .eq. 'bin') then
    fieldbin(:,:)=dropsm2d(:,nx,:)
    name='mass_dist'
    units='kg/kg/ln(r)'
    call save_dg('bin',fieldbin,name,i_dgtime,units)

    if (bintype .eq. 'tau') then
        fieldbin(:,:)=dropsn2d(:,nx,:)
        name='num_dist'
        units='1/kg/ln(r)'
        call save_dg('bin',fieldbin,name,i_dgtime,units)
    endif

elseif (ampORbin .eq. 'amp') then
    fieldbin(:,:)=dropsinitm2d(:,nx,:)
    name='mass_dist_init'
    units='kg/kg/ln(r)'
    call save_dg('bin',fieldbin,name,i_dgtime,units)

    if (bintype .eq. 'tau') then
        fieldbin(:,:)=dropsinitn2d(:,nx,:)
        name='num_dist_init'
        units='1/kg/ln(r)'
        call save_dg('bin',fieldbin,name,i_dgtime,units)
    end if

endif

!parameters
do j=1,nx
   if (nx==1) then
      call save_dg(guessc2d(:,nx,1),'nu_c', i_dgtime,units='unitless', dim='z')
      call save_dg(guessr2d(:,nx,1),'nu_r', i_dgtime,units='unitless', dim='z')
      call save_dg(guessc2d(:,nx,2),'Dn_c', i_dgtime,units='m', dim='z')
      call save_dg(guessr2d(:,nx,2),'Dn_r', i_dgtime,units='m', dim='z')
   else
      call save_dg(guessc2d(:,:,1),'nu_c', i_dgtime,units='unitless', dim='z,x')
      call save_dg(guessr2d(:,:,1),'nu_r', i_dgtime,units='unitless', dim='z,x')
      call save_dg(guessc2d(:,:,2),'Dn_c', i_dgtime,units='m', dim='z,x')
      call save_dg(guessr2d(:,:,2),'Dn_r', i_dgtime,units='m', dim='z,x')
   endif
enddo

end Subroutine mphys_amp_interface

end module mphys_amp
