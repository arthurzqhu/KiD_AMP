
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module to interface with choice of advection schemes
!

module advection_interface

  Use runtime, only : time
  Use parameters, only : num_h_moments, num_h_bins, nz, nspecies &
            ,num_aero_moments,num_aero_bins, aero_mom_init, dt, nx &
            ,max_char_len
  Use column_variables
  Use ultf2d_mod, only : ultf2d_interface 
  Use switches


  Implicit none

    real(wp), allocatable :: field_adv(:,:)
    real(wp), allocatable :: field(:,:)

 contains

  subroutine advect_column(scheme_id) 
    integer, intent(in), optional :: scheme_id
    integer :: cscheme_id
!    real(wp), allocatable :: field(:,:)
    integer :: ih, imom, ibin, k, j
    real(8) :: norm_factor, scale_factor=1

    if (present(scheme_id))then
       cscheme_id=scheme_id
    else
       cscheme_id=0
    end if
    ! Advect each model variable in column...

    allocate(field(nz,0:nx+1))
    allocate(field_adv(nz,0:nx+1))

    if (.not. l_fix_theta)then
  
       field(:,:)=theta(:,:)
    ! Temperature perturbation
       call generic_advection(  &
            &  field            &
            & ,field_adv        &
            & ,cscheme_id       &
            & ,wth_surf         & 
            )

       dtheta_adv(:,:)=field_adv(:,:)
    else
       dtheta_adv(:,:)=0.0
    endif

    ! vapour
    if (.not. l_fix_qv)then
       
       field(:,:)=qv(:,:)

       norm_factor=sum(field)
       if (norm_factor .ne. 0.) field=scale_factor*field/norm_factor

       call generic_advection(      &
            &  field                &
            & ,field_adv            &
            & ,cscheme_id           &
            & ,wqv_surf             &
            )

       if (norm_factor .ne. 0.) then 
          field=field*norm_factor/scale_factor
          field_adv=field_adv*norm_factor/scale_factor
       endif

       do j=1,nx
          do k=1,nz
             if (field(k,j)+dt*field_adv(k,j) < 0.0) &
                  field_adv(k,j)=-0.99999*(field(k,j)/dt)
          enddo
       end do
       dqv_adv(:,:)=field_adv(:,:) 
    else
       dqv_adv(:,:)=0.0
    endif

    field(:,:)=ss(:,:)

    norm_factor=sum(field)
    if (norm_factor .ne. 0.) field=scale_factor*field/norm_factor

    call generic_advection(      &
         &  field                &
         & ,field_adv            &
         & ,scheme_id)

    if (norm_factor .ne. 0.) then 
       field=field*norm_factor/scale_factor
       field_adv=field_adv*norm_factor/scale_factor
    endif
   
    dss_adv(:,:)=field_adv(:,:)

    ! Aerosol
    do ih=1,naerosol
       do ibin=1,num_aero_bins(ih)
          do imom=1,num_aero_moments(ih)
             do j=1,nx
                do k=1,nz
                   field(k,j)=aerosol(k,j,ih)%moments(ibin,imom)
                end do
             enddo

             norm_factor=sum(field)
             if (norm_factor .ne. 0.) field=scale_factor*field/norm_factor

             call generic_advection(            &
                  &  field                      &
                  & ,field_adv                  &
                  & ,scheme_id)

             if (norm_factor .ne. 0.) then 
                field=field*norm_factor/scale_factor
                field_adv=field_adv*norm_factor/scale_factor
             endif

             do j = 1,nx
                do k=1,nz
                   if (field(k,j)+dt*field_adv(k,j) < 0.0) &
                        field_adv(k,j)=-0.99999*(field(k,j)/dt)
                   daerosol_adv(k,j,ih)%moments(ibin,imom)=field_adv(k,j)
                 enddo
             end do             
          end do
       end do
    end do    

    ! Hydrometeors
    do ih=1,nspecies
       do ibin=1,num_h_bins(ih)
          do imom=1,num_h_moments(ih)

             do j=1,nx
                do k=1,nz
                   field(k,j)=hydrometeors(k,j,ih)%moments(ibin,imom)
                end do
             enddo

             norm_factor=sum(field)

             if (norm_factor /= norm_factor) then
               print*, 'ih, ibin, imom', ih, ibin, imom
               print*, 'nan in hydrometeors, check dhydrometeors_* at the previous timestep of', time
               stop
             endif

             if (norm_factor .ne. 0.) field=scale_factor*field/norm_factor

             call generic_advection(            &
                &  field                      &
                & ,field_adv                  &
                & ,scheme_id)

             ! ! print*, 'fieldadv', field_adv(1,1), norm_factor, scale_factor
             ! if (norm_factor .ne. norm_factor .or. norm_factor > HUGE(norm_factor)) then
             !    print*, 'field', field
             !    stop
             ! endif

             if (norm_factor .ne. 0.) then 
                field=field*norm_factor/scale_factor
                field_adv=field_adv*norm_factor/scale_factor
             endif

             do j=1,nx
                do k=1,nz
                   ! if (field(k,j)<0) then
                   !    print*, k, ih, ibin, imom, field(k,j)
                   !    print*, 'hyd', hydrometeors(k,j,ih)%moments(ibin,imom)
                   !    print*, 'dmphys', dhydrometeors_mphys(k,j,ih)%moments(ibin,imom)*dt
                   !    stop 'moment<0'
                   ! endif
                   if (field(k,j)+dt*field_adv(k,j) < 0.0) &
                      field_adv(k,j)=-.999999*(field(k,j)/dt)
                   dhydrometeors_adv(k,j,ih)%moments(ibin,imom)=field_adv(k,j)

                enddo
             end do
             ! print*, 'adv', time, hydrometeors(1,1,1)%moments(1,1), hydrometeors(1,1,2)%moments(1,1)

          end do 
       end do
    end do

    deallocate(field_adv)
    deallocate(field)

  end subroutine advect_column

  subroutine generic_advection(field, field_adv, scheme_id, surface_flux)


    real(wp) :: field(:,:)
    integer, intent(in) :: scheme_id
    real(wp), intent(out) :: field_adv(:,:)
    real(wp), intent(in), optional :: surface_flux(:)
    ! local variables
    real(wp) :: csurface_flux(0:nx+1)

    if (present(surface_flux))then
       csurface_flux=surface_flux
    else
       csurface_flux=0.
    end if

    !declarations for 2D ultimate advection (scheme_id=0)
    select case (scheme_id)
    case(0) ! 2D Ultimate advection scheme from LEM

       call ultf2d_interface(field, v_half, w_half, z, z_half, rho, rho_half, &
            field_adv, csurface_flux)

    case(1) ! simple first order SL advection
       call sl_advection(field, field_adv)
    end select

  end subroutine generic_advection

  subroutine sl_advection(field, field_adv)
    !
    ! A simple semi-lagrangian scheme
    !

    real(wp), intent(in) :: field(:,:)
    real(wp), intent(out) :: field_adv(:,:)

    ! local variables
    integer :: k
    real(wp) :: z_depart

    field_adv=field
    Print*, ' SL advection not yet implemented'
    STOP

    ! Calculate departure point
    
  end subroutine sl_advection
          

end module advection_interface
