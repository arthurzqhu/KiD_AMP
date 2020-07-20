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
  Use physconst, only : p0, r_on_cp, pi

  Use module_hujisbm
  Use micro_prm
  Use diagnostics, only: save_dg, i_dgtime, my_save_dg_bin_dp
  Use switches, only: l_advect,l_diverge
  Use namelists, only: bintype, ampORbin

  Implicit None

  !Logical switches
  logical :: micro_unset=.True.
  integer:: ih, imom
  character(max_char_len) :: name, units

contains

  Subroutine mphys_amp_interface
    use parameters, only: flag_count,max_nbins
    integer :: i, j, k, imom
    real, dimension(nz,nx) :: t2d, p2d, qv2d
    real(8), dimension(nz,nx,num_h_moments(1)) :: Mpc2d
    real(8), dimension(nz,nx,num_h_moments(2)) :: Mpr2d
    real(8),dimension(nz,nx,10) :: mc,mr
    real(8), save, dimension(nz,nx,2) :: guessc2d,guessr2d
    real, dimension(nz,nx,max_nbins) :: aer2d,dropsm2d,dropsn2d,dropsinit2d
    real, dimension(nz) :: field
    real(8),dimension(nz,flag_count) :: fieldflag
    real(8), dimension(nz) :: fielddp
    real(8), dimension(nz,max_nbins) :: fielddp2d
    real(8), dimension(nz,nx,2,flag_count) :: flag
    character(1) :: Mnum

    do i=1,nx
       do k=1,nz
          !if (l_advect .and. l_diverge) then
          !  t2d(k,i) = (theta(k,i) + (dtheta_adv(k,i)+dtheta_div(k,i))*dt )*exner(k,i)
          !  qv2d(k,i) = qv(k,i) + (dqv_adv(k,i)+dqv_div(k,i))*dt
          !elseif (l_advect) then
          !  t2d(k,i) = (theta(k,i) + dtheta_adv(k,i)*dt )*exner(k,i)
          !  qv2d(k,i) = qv(k,i) + dqv_adv(k,i)*dt
          !elseif (l_diverge) then
          !  t2d(k,i) = (theta(k,i) + dtheta_div(k,i)*dt )*exner(k,i)
          !  qv2d(k,i) = qv(k,i) + dqv_div(k,i)*dt
          !else
            t2d(k,i) = theta(k,i)*exner(k,i)
            qv2d(k,i) = qv(k,i)
          !endif
          p2d(k,i) = p0*exner(k,i)**(1./r_on_cp)

          if (ampORbin .eq. 'amp') then
             do imom=1,num_h_moments(1)
                Mpc2d(k,i,imom) = hydrometeors(k,i,1)%moments(1,imom)
         !       Mpc2d(k,i,imom) = hydrometeors(k,i,1)%moments(1,imom) &
         !            + (dhydrometeors_adv(k,i,1)%moments(1,imom) &
         !            + dhydrometeors_div(k,i,1)%moments(1,imom))*dt
             enddo
             if (any(Mpc2d(k,i,:)==0.)) Mpc2d(k,i,:)=0.
             do imom=1,num_h_moments(2)
                Mpr2d(k,i,imom) = hydrometeors(k,i,2)%moments(1,imom)
         !       Mpr2d(k,i,imom) = hydrometeors(k,i,2)%moments(1,imom) &
         !            + (dhydrometeors_adv(k,i,2)%moments(1,imom) &
         !            + dhydrometeors_div(k,i,2)%moments(1,imom))*dt
             enddo
             if (any(Mpr2d(k,i,:)==0.)) Mpr2d(k,i,:)=0.
          else !bin
             Mpc2d=0.;Mpr2d=0.
             do j=1,nkr
                dropsm2d(k,i,j)=hydrometeors(k,i,1)%moments(j,1)/col
         !       dropsm2d(k,i,j)=hydrometeors(k,i,1)%moments(j,1) &
         !            + (dhydrometeors_adv(k,i,1)%moments(j,1) &
         !            + dhydrometeors_div(k,i,1)%moments(j,1))*dt
                 if (bintype .eq. 'tau')  then
                     dropsn2d(k,i,j)=hydrometeors(k,i,1)%moments(j,2)/col
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
   if (ampORbin .eq. 'amp') then
      dropsm2d=0.
      dropsinit2d=0.
!print*,'s',Mpc2d(18,1,:)
      call mp_amp(Mpc2d,Mpr2d,guessc2d,guessr2d, &
           p2d,t2d,qv2d,aer2d,dropsm2d,mc,mr,flag,dropsinit2d)
!print*,'e',Mpc2d(18,1,:)
   elseif (ampORbin .eq. 'bin') then
      if (bintype .eq. 'sbm') then
         call mp_sbm(dropsm2d,p2d,t2d,qv2d,aer2d,mc,mr)
      elseif (bintype .eq. 'tau') then
         call mp_tau(dropsm2d,dropsn2d,t2d,qv2d,mc,mr)
      endif
   endif

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
      !if (l_advect) dtheta_mphys(k,i)=dtheta_mphys(k,i)-dtheta_adv(k,i)
      !if (l_diverge) dtheta_mphys(k,i)=dtheta_mphys(k,i)-dtheta_div(k,i)

      dqv_mphys(k,i)=(qv2d(k,i)-qv(k,i))/dt
      !if (l_advect) dqv_mphys(k,i)=dqv_mphys(k,i)-dqv_adv(k,i)
      !if (l_diverge) dqv_mphys(k,i)=dqv_mphys(k,i)-dqv_div(k,i)

      if (ampORbin .eq. 'amp') then ! when bulk
         do imom=1,num_h_moments(1)
            dhydrometeors_mphys(k,i,1)%moments(1,imom)= &
                 (Mpc2d(k,i,imom)-hydrometeors(k,i,1)%moments(1,imom))/dt
         !   if (l_advect) dhydrometeors_mphys(k,i,1)%moments(1,imom)= &
         !                 dhydrometeors_mphys(k,i,1)%moments(1,imom) &
         !                 -dhydrometeors_adv(k,i,1)%moments(1,imom)
         !   if (l_diverge) dhydrometeors_mphys(k,i,1)%moments(1,imom)= &
         !                  dhydrometeors_mphys(k,i,1)%moments(1,imom) &
         !                  -dhydrometeors_div(k,i,1)%moments(1,imom)
         enddo
         do imom=1,num_h_moments(2)
            dhydrometeors_mphys(k,i,2)%moments(1,imom)= &
                 (Mpr2d(k,i,imom)-hydrometeors(k,i,2)%moments(1,imom))/dt
         !   if (l_advect) dhydrometeors_mphys(k,i,2)%moments(1,imom)= &
         !                 dhydrometeors_mphys(k,i,2)%moments(1,imom) &
         !                 -dhydrometeors_adv(k,i,2)%moments(1,imom)
         !   if (l_diverge) dhydrometeors_mphys(k,i,2)%moments(1,imom)= &
         !                  dhydrometeors_mphys(k,i,2)%moments(1,imom) &
         !                  -dhydrometeors_div(k,i,2)%moments(1,imom)
         enddo

      else ! when bin
         do j=1,nkr
            dhydrometeors_mphys(k,i,1)%moments(j,1)= &
                 (dropsm2d(k,i,j)*col-hydrometeors(k,i,1)%moments(j,1))/dt
         !   if (l_advect) dhydrometeors_mphys(k,i,1)%moments(j,1)= &
         !                 dhydrometeors_mphys(k,i,1)%moments(j,1) &
         !                 -dhydrometeors_adv(k,i,1)%moments(j,1)
         !   if (l_diverge) dhydrometeors_mphys(k,i,1)%moments(j,1)= &
         !                  dhydrometeors_mphys(k,i,1)%moments(j,1) &
         !                  -dhydrometeors_div(k,i,1)%moments(j,1)
            if (bintype .eq. 'tau') then
                dhydrometeors_mphys(k,i,1)%moments(j,2)= &
                     (dropsn2d(k,i,j)*col-hydrometeors(k,i,1)%moments(j,2))/dt
            endif
         enddo

      endif

   enddo
enddo

! Save some diagnostics
!fitting flag
!if (imomc1.ne.3) then
!   field(:)=flag(:,nx,1,1)
!   call save_dg(field,'old_fitting_flag_cloud', i_dgtime,units='unitless', dim='z')
!   field(:)=flag(:,nx,2,1)
!   call save_dg(field,'old_fitting_flag_rain', i_dgtime,units='unitless', dim='z')
!
!   field(:)=flag(:,nx,1,2)
!   call save_dg(field,'fitting_flag_cloud_oob', i_dgtime,units='unitless', dim='z')
!   field(:)=flag(:,nx,2,2)
!   call save_dg(field,'fitting_flag_rain_oob', i_dgtime,units='unitless', dim='z')
!
!   field(:)=flag(:,nx,1,3)
!   call save_dg(field,'fitting_flag_cloud_x_intol', i_dgtime,units='unitless', dim='z')
!   field(:)=flag(:,nx,2,3)
!   call save_dg(field,'fitting_flag_rain_x_intol', i_dgtime,units='unitless', dim='z')
!
!   field(:)=flag(:,nx,1,4)
!   call save_dg(field,'fitting_flag_cloud_y_intol', i_dgtime,units='unitless', dim='z')
!   field(:)=flag(:,nx,2,4)
!   call save_dg(field,'fitting_flag_rain_y_intol', i_dgtime,units='unitless', dim='z')
!endif

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
     fielddp(:)=mc(:,nx,i)
     call save_dg(fielddp,name,i_dgtime,units,dim='z')

     name='diagM'//Mnum//'_rain'
     fielddp(:)=mr(:,nx,i)
     call save_dg(fielddp,name,i_dgtime,units,dim='z')
   enddo

!bin distributions
if (ampORbin .eq. 'bin') then
    fielddp2d(:,:)=dropsm2d(:,nx,:)
    name='mass_dist'
    units='kg/kg/ln(r)'
    call save_dg('bin',fielddp2d,name,i_dgtime,units)

    if (bintype .eq. 'tau') then
        fielddp2d(:,:)=dropsn2d(:,nx,:)
        name='num_dist'
        units='1/kg/ln(r)'
        call save_dg('bin',fielddp2d,name,i_dgtime,units)
    endif

elseif (ampORbin .eq. 'amp') then
    fielddp2d(:,:)=dropsinit2d(:,nx,:)
    name='drops_init'
    units='kg/kg/ln(r)'
    call save_dg('bin',fielddp2d,name,i_dgtime,units)
endif

!parameters
  field(:)=guessc2d(:,nx,1)
  call save_dg(field,'shape_parameter_cloud', i_dgtime,units='unitless', dim='z')
  field(:)=guessr2d(:,nx,1)
  call save_dg(field,'shape_parameter_rain', i_dgtime,units='unitless', dim='z')
  field(:)=guessc2d(:,nx,2)
  call save_dg(field,'characteristic_diameter_cloud', i_dgtime,units='m', dim='z')
  field(:)=guessr2d(:,nx,2)
  call save_dg(field,'characteristic_diameter_rain', i_dgtime,units='m', dim='z')


end Subroutine mphys_amp_interface

end module mphys_amp
