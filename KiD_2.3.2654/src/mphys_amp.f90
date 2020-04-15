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

  Implicit None

  !Logical switches
  logical :: micro_unset=.True.
  integer:: ih, imom
  character(max_char_len) :: name, units

contains

  Subroutine mphys_amp_interface
    use parameters, only: flag_count
    integer :: i, j, k, imom
    real, dimension(nz,nx) :: t2d, p2d, qv2d
    real(8), dimension(nz,nx,num_h_moments(1)) :: Mpc2d
    real(8), dimension(nz,nx,num_h_moments(2)) :: Mpr2d
    real(8),dimension(nz,nx,10) :: mc,mr
    real(8), save, dimension(nz,nx,2) :: guessc2d,guessr2d
    real, dimension(nz,nx,nkr) :: aer2d,drops2d,dropsinit2d
    real, dimension(nz) :: field
    real(8),dimension(nz,flag_count) :: fieldflag
    real(8), dimension(nz) :: fielddp
    real(8), dimension(nz,nkr) :: fielddp2d
    real(8), dimension(nz,nx,2,flag_count) :: flag
    real(8), dimension(nz,nx,2) :: oMxM3,oMyM3,nMxM3,nMyM3
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

          if (imomc1.ne.3) then
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
          else
             Mpc2d=0.;Mpr2d=0.
             do j=1,nkr
                drops2d(k,i,j)=hydrometeors(k,i,1)%moments(j,1)
         !       drops2d(k,i,j)=hydrometeors(k,i,1)%moments(j,1) &
         !            + (dhydrometeors_adv(k,i,1)%moments(j,1) &
         !            + dhydrometeors_div(k,i,1)%moments(j,1))*dt
             enddo
          endif
       enddo
   enddo

   ! Initialise microphysics
   if (micro_unset)then
      if (imomc1.ne.3) then
         guessc2d(:,:,1) = h_shape(1) !shape parameter
         guessr2d(:,:,1) = h_shape(2)
         guessc2d(:,:,2) = 0.001         !characteristic diameter dn
         guessr2d(:,:,2) = 0.001
         call amp_init(aer2d,Mpc2d,Mpr2d,guessc2d,guessr2d)
      else
         call sbm_init(aer2d,drops2d)
      endif
      micro_unset=.False.
   endif

   aer2d = 0.
   if (imomc1 .ne. 3) then
      drops2d=0.
      dropsinit2d=0.
!print*,'s',Mpc2d(18,1,:)
      call mp_amp(Mpc2d,Mpr2d,guessc2d,guessr2d, &
           p2d,t2d,qv2d,aer2d,drops2d,mc,mr,flag,dropsinit2d,oMxM3,oMyM3,nMxM3,&
           nMyM3)
!print*,'e',Mpc2d(18,1,:)
   else
      call mp_sbm(drops2d,p2d,t2d,qv2d,aer2d,mc,mr)
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

         if (imomc1.ne.3) then
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
         else
            !do j=1,nkr
               dhydrometeors_mphys(k,i,1)%moments(j,1)= &
                    (drops2d(k,i,j)-hydrometeors(k,i,1)%moments(j,1))/dt
            !   if (l_advect) dhydrometeors_mphys(k,i,1)%moments(j,1)= &
            !                 dhydrometeors_mphys(k,i,1)%moments(j,1) &
            !                 -dhydrometeors_adv(k,i,1)%moments(j,1)
            !   if (l_diverge) dhydrometeors_mphys(k,i,1)%moments(j,1)= &
            !                  dhydrometeors_mphys(k,i,1)%moments(j,1) &
            !                  -dhydrometeors_div(k,i,1)%moments(j,1)
            !enddo

         endif
      enddo
   enddo

! Save some diagnostics
!fitting flag
  if (imomc1.ne.3) then
     field(:)=flag(:,nx,1,1)
     call save_dg(field,'old_fitting_flag_cloud', i_dgtime,units='unitless', dim='z')
     field(:)=flag(:,nx,2,1)
     call save_dg(field,'old_fitting_flag_rain', i_dgtime,units='unitless', dim='z')

     field(:)=flag(:,nx,1,2)
     call save_dg(field,'fitting_flag_cloud_oob', i_dgtime,units='unitless', dim='z')
     field(:)=flag(:,nx,2,2)
     call save_dg(field,'fitting_flag_rain_oob', i_dgtime,units='unitless', dim='z')

     field(:)=flag(:,nx,1,3)
     call save_dg(field,'fitting_flag_cloud_x_intol', i_dgtime,units='unitless', dim='z')
     field(:)=flag(:,nx,2,3)
     call save_dg(field,'fitting_flag_rain_x_intol', i_dgtime,units='unitless', dim='z')

     field(:)=flag(:,nx,1,4)
     call save_dg(field,'fitting_flag_cloud_y_intol', i_dgtime,units='unitless', dim='z')
     field(:)=flag(:,nx,2,4)
     call save_dg(field,'fitting_flag_rain_y_intol', i_dgtime,units='unitless', dim='z')
  endif

  ! Saving ratios
  field(:)=oMxM3(:,nx,1)
  call save_dg(field,'oMxM3_c',i_dgtime,units='n/a', dim='z')
  field(:)=oMxM3(:,nx,2)
  call save_dg(field,'oMxM3_r',i_dgtime,units='n/a', dim='z')
  field(:)=oMyM3(:,nx,1)
  call save_dg(field,'oMyM3_c',i_dgtime,units='n/a', dim='z')
  field(:)=oMyM3(:,nx,2)
  call save_dg(field,'oMyM3_r',i_dgtime,units='n/a', dim='z')
  field(:)=nMxM3(:,nx,1)
  call save_dg(field,'nMxM3_c',i_dgtime,units='n/a', dim='z')
  field(:)=nMxM3(:,nx,2)
  call save_dg(field,'nMxM3_r',i_dgtime,units='n/a', dim='z')
  field(:)=nMyM3(:,nx,1)
  call save_dg(field,'nMyM3_c',i_dgtime,units='n/a', dim='z')
  field(:)=nMyM3(:,nx,2)
  call save_dg(field,'nMyM3_r',i_dgtime,units='n/a', dim='z')

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
  fielddp2d(:,:)=drops2d(:,nx,:)
  name='drops'
  units='kg/kg/ln(r)'
  call save_dg('bin',fielddp2d,name,i_dgtime,units)

  fielddp2d(:,:)=dropsinit2d(:,nx,:)
  name='drops_init'
  units='kg/kg/ln(r)'
  call save_dg('bin',fielddp2d,name,i_dgtime,units)

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
