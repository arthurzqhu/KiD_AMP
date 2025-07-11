! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Save any diagnostics here...
!
!

Module diagnostics

  Use typeKind
  Use parameters
  Use column_variables
  Use runtime
  Use physconst, only: r_on_cp, p0, pi
  Use common_physics, only : qsaturation, qisaturation
  Use header_data
  Use switches, only: l_diverge, wctrl

  Use namelists, only: KiD_outdir, KiD_outfile, fileNameOut,ampORbin,nmom_diag, imomc1, imomc2, &
    l_save_tend, l_save_adv, l_save_div, l_save_mphys, rain_source
  use micro_prm, only: QtoM3, pmomsc

  Implicit none

  type, public :: dgID
     character(max_char_len) :: name
     integer            :: ID=unset_integer
  end type dgID

  type, public :: diagBinData   ! Holds bin data (e.g. sizes)
     character(max_char_len) :: name, units, longname
     character(max_char_len) :: dim
     real(wp), pointer :: data(:)
  end type diagBinData

  type, public :: diagScalarTS   ! Holds scalar timeseries data
     character(max_char_len) :: name, units, longname
     character(max_char_len) :: dim
     real(wp), pointer :: data(:)
  end type diagScalarTS

  type, public :: diag1DTS   ! Holds 1D timeseries data
     character(max_char_len) :: name, units, longname
     character(max_char_len) :: dim
     real(wp), pointer :: data(:,:)
  end type diag1DTS

  type, public :: diag1D_nxTS   ! Holds 1D timeseries data in nx (only used for 2D)
     character(max_char_len) :: name, units, longname
     character(max_char_len) :: dim
     real(wp), pointer :: data(:,:)
  end type diag1D_nxTS

  type, public :: diag2DTS   ! Holds 2D timeseries data
     character(max_char_len) :: name, units, longname
     character(max_char_len) :: dim
     real(wp), pointer :: data(:,:,:)
  end type diag2DTS

  type, public :: diag_bin2DTS   ! Holds 2D bin-height timeseries data
     character(max_char_len) :: name, units, longname
     character(max_char_len) :: dim
     real(wp), pointer :: data(:,:,:)
  end type diag_bin2DTS

  type, public :: diag_proc2DTS   ! Holds 2D procrate-height timeseries data
      character(max_char_len) :: name, units, longname
      character(max_char_len) :: dim
      real(wp), pointer :: data(:,:,:)
  end type diag_proc2DTS

  type, public :: diag_flag2DTS   ! Holds 2D flag-height timeseries data
      character(max_char_len) :: name, units, longname
      character(max_char_len) :: dim
      real(wp), pointer :: data(:,:,:)
  end type diag_flag2DTS

  type, public :: dgIDarray
     type(dgID) :: dgIDs(max_dgs)
     integer    :: nIds=0
     integer    :: lastAccessed=0
  end type dgIDarray

  type(diag1DTS), target :: instant_column(max_dgs) ! Instantaneous fields
  type(diag1D_nxTS), target :: instant_nx(max_dgs) ! Instaneous nx field for 2D sim
  type(diag2DTS), target :: instant_2D(max_dgs) ! Instantaneous fields
  type(diagScalarTS), target :: scalars(max_dgs) ! Instantaneous scalars
  type(diagBinData), target :: binData(max_dgs) ! bin data
  type(diag_bin2DTS), target :: instant_bindgs(max_dgs) ! Instantaneous bin diags
  type(diag_proc2DTS), target :: instant_procdgs(max_dgs) ! Instantaneous proc rates diags
  type(diag_flag2DTS), target :: instant_flagdgs(max_dgs) ! Instantaneous flag diags

  type(dgIDarray), save, target :: ID_instant_column
  type(dgIDarray), save, target :: ID_instant_nx
  type(dgIDarray), save, target :: ID_instant_2D
  type(dgIDarray), save, target :: ID_scalars
  type(dgIDarray), save, target :: ID_binData
  type(dgIDarray), save, target :: ID_instant_bindgs
  type(dgIDarray), save, target :: ID_instant_procdgs
  type(dgIDarray), save, target :: ID_instant_flagdgs

  integer :: maxn_dgtimes=0

  interface allocate_dgs
     module procedure allocate_dgs_1DTS, allocate_dgs_2DTS, &
          allocate_dgs_ScalarTS, allocate_dgs_bindata, &
          allocate_dgs_bindgs, allocate_dgs_flagdgs
  end interface

  interface save_dg

     module procedure save_dg_1D_sp, save_dg_1D_point_sp,            &
        save_dg_1D_range_sp, save_dg_scalar_sp, save_dg_1D_range_dp, &
        save_dg_1D_dp, save_dg_1D_point_dp, save_dg_scalar_dp,       &
        save_dg_2D_point_sp, save_dg_2D_point_dp,                    &
        save_dg_2D_sp, save_dg_2D_dp,                                &
        save_dg_bin_sp, save_dg_bin_dp!, save_dg_flag_sp,             &
        ! save_dg_flag_dp, save_dg_bin_point_sp, save_dg_bin_point_dp

  end interface

  interface save_dg_2d

     module procedure save_dg_2D_point_sp, save_dg_2D_point_dp,      &
        save_dg_2D_sp, save_dg_2D_dp

  end interface

  interface save_binData

     module procedure save_binData_sp, save_binData_dp

  end interface

  interface print_dg
     module procedure print_dg_1D_point
  end interface

  character(*), parameter :: fmta= "(A5, i2.2, A1, i2.2, A1, i4.4, A7"
  character(*), parameter :: fmtb= ", i2.2, A1, i2.2, A1, i2.2, A5, A4)"
  character(*), parameter :: fmt=fmta//fmtb

  integer :: k, j

  integer :: k_here, i_here, j_here

contains

  subroutine query_dgstep

    if(mod(time,dg_dt)<dt.and.time>dgstart)then
       l_dgstep=.true.
       i_dgtime=i_dgtime+1
    else
       l_dgstep=.false.
    end if

  end subroutine query_dgstep


  subroutine save_diagnostics_1d
    use global_fun, only: momk
    real(wp), allocatable :: field(:)
    real(wp), allocatable :: field_flag(:,:)
    real(wp), allocatable :: field_nx(:)
    real(wp), allocatable :: field_bin_c(:)
    real(wp), allocatable :: field_bin_r(:)
    real(wp), allocatable :: field_bin(:,:)
    real(wp), allocatable :: field_moms(:,:)
    real(wp), allocatable :: m3_arr(:), m0_arr(:), mx_arr(:)
    character(max_char_len) :: name, units, dims
    integer :: k, ih, imom, ibin, ift

    if (.not. l_dgstep) return

    !
    ! Instantaneous (for 1-D run) or horizontally averaged (for 2d run)
    ! column data
    !
    allocate(field(nz))
    allocate(field_flag(nz,flag_count))
    allocate(field_nx(0:nx+1))
    allocate(field_bin_c(nz))
    allocate(field_bin_r(nz))
    allocate(field_bin(nz,max_nbins))
    allocate(field_moms(nz,nmom_diag))

    ! set dimensions for the diagnostic column or grid output, i.e.
    ! 1-D set-up is dimensioned with column - 'z'

    dims = 'z'

    ! pressure
    field(:)=pmb(:,nx)
    call save_dg(field, 'pressure', i_dgtime,  units='mb',dim=dims)
    !vertical velocity
    field(:)=w(:,nx)
    call save_dg(field, 'w', i_dgtime, units='m/s',dim=dims)
    !vapour
    field(:)=qv(:,nx)
    call save_dg(field, 'vapour', i_dgtime,  units='kg/kg',dim=dims)
    !theta
    field(:)=theta(:,nx)
    call save_dg(field,'theta', i_dgtime,  units='K',dim=dims)
    !temperature(K)
    field(:)=TdegK(:,nx)
    call save_dg(field, 'temperature', i_dgtime, units='K',dim=dims)

    ! ! aerosol
    ! do ih=1,naerosol
    !    do imom=1,num_aero_moments(ih)
    !       do k=1,nz
    !          field(k)=sum(aerosol(k,nx,ih)%moments(:,imom))
    !       end do
    !       name=trim(aero_names(ih))//'_'//trim(aero_mom_names(imom))
    !       units=trim(aero_mom_units(imom))
    !       call save_dg(field, name, i_dgtime, units,dim=dims)
    !    end do
    ! end do

    !hydrometeors

    do ih=1,nspecies
       do imom=1,num_h_moments(ih)
          name=trim(h_names(ih))//'_'//trim(mom_names(imom))
          units=trim(mom_units(imom))
          do k=1,nz
             if (ampORbin .eq. 'bin') then
                if (ih==1) then ! cloud
                  field(k)=sum(hydrometeors(k,nx,1)%moments(1:split_bins,imom))
                elseif (ih==2) then ! rain
                  field(k)=sum(hydrometeors(k,nx,1)%moments(split_bins+1:max_nbins,imom))
                endif
                if (imom==1) field(k)=field(k)*QtoM3
             else
                field(k)=sum(hydrometeors(k,nx,ih)%moments(:,imom))
             endif
          end do
          call save_dg(field, name, i_dgtime,  units,dim=dims)

       end do

    end do

    ! output extra moments

    if (ampORbin .eq. 'bin') then
      do ih = 1,nspecies
        do imom=3,4
          name=trim(h_names(ih))//'_'//trim(mom_names(imom))
          units='m^k/m2'
          do k=1,nz
            if (ih==1) then
              m3_arr = hydrometeors(k,nx,1)%moments(1:split_bins,1)
              m0_arr = hydrometeors(k,nx,1)%moments(1:split_bins,2)
            elseif(ih==2) then
              m3_arr = hydrometeors(k,nx,1)%moments(split_bins+1:max_nbins,1)
              m0_arr = hydrometeors(k,nx,1)%moments(split_bins+1:max_nbins,2)
            endif
            if (ih==1) field(k)=momk(m3_arr, m0_arr, dble(pmomsc(imom)), 1)
            if (ih==2) field(k)=momk(m3_arr, m0_arr, dble(pmomsc(imom)), split_bins+1)
          enddo

          call save_dg(field, name, i_dgtime,  units,dim=dims)

        enddo
      enddo
    endif


    !4D flags
    !do ift=1,4 !iflagtype
    ! call save_dg('flag',field_flag,'fitting_flag',i_dgtime,units='unitless',dim=dims)
    !end do

    !========================================================
    ! Tendency terms
    !========================================================
    if (l_save_tend) then

    if (l_save_adv) then
    !
    ! Advection
    !
    !-----------
    ! theta
    field(:)=dtheta_adv(:,nx)
    call save_dg(field(:), 'dtheta_adv', i_dgtime,  units='K/s',dim=dims)
    ! qv
    field(:)=dqv_adv(:,nx)
    call save_dg(field, 'dqv_adv', i_dgtime,  units='kg/kg/s',dim=dims)

    ! ! aerosols
    ! do ih=1,naerosol
    !    do imom=1,num_aero_moments(ih)
    !       do k=1,nz
    !          field(k)=sum(daerosol_adv(k,nx,ih)%moments(:,imom))
    !       end do
    !       name=trim(aero_names(ih))//'_'//trim(aero_mom_names(imom))//'_adv'
    !       units=trim(aero_mom_units(imom)//'/s')
    !       call save_dg(field , name, i_dgtime, units,dim=dims)
    !    end do
    ! end do

    !hydrometeors
    do ih=1,nspecies
       do imom=1,num_h_moments(ih)
          name=trim(h_names(ih))//'_'//trim(mom_names(imom))//'_adv'
          units=trim(mom_units(imom))//'/s'
          do k=1,nz
             if (ampORbin .eq. 'bin')then
                if (ih==1) then
                    field(k)=sum(dhydrometeors_adv(k,nx,1)%moments(1:split_bins,imom))
                elseif (ih==2) then
                    field(k)=sum(dhydrometeors_adv(k,nx,1)%moments(split_bins+1:max_nbins,imom))
                endif
                if (imom==1) field(k)=field(k)*QtoM3
             else
                 field(k)=sum(dhydrometeors_adv(k,nx,ih)%moments(:,imom))
             endif
          end do
          call save_dg(field, name, i_dgtime,  units,dim=dims)
       end do
    end do
    endif

    if (l_save_div) then
    !
    !-----------
    !
    ! Divergence
    !
    !-----------
    if (l_diverge)then
       ! theta
       field(:)=dtheta_div(:,nx)
       call save_dg(field, 'dtheta_div', i_dgtime,  units='K/s',dim=dims)
       ! qv tendency
       field(:)=dqv_div(:,nx)
       call save_dg(field, 'dqv_div', i_dgtime,  units='kg/kg/s',dim=dims)

       !!aerosol
       !do ih=1,naerosol
       !   do imom=1,num_aero_moments(ih)
       !      do k=1,nz
       !         field(k)=sum(daerosol_div(k,nx,ih)%moments(:,imom))
       !      end do
       !      name=trim(aero_names(ih))//'_'//trim(aero_mom_names(imom))//'_div'
       !      units=trim(aero_mom_units(imom)//'/s')
       !      call save_dg(field, name, i_dgtime, units,dim=dims)
       !   end do
       !end do

       !hydrometeors
       do ih=1,nspecies
          do imom=1,num_h_moments(ih)
             name=trim(h_names(ih))//'_'//trim(mom_names(imom))//'_div'
             units=trim(mom_units(imom)//'/s')


             do k=1,nz
                if (ampORbin .eq. 'bin')then
                   if (ih==1) then
                       field(k)=sum(dhydrometeors_div(k,nx,1)%moments(1:split_bins,imom))
                   elseif (ih==2) then
                       field(k)=sum(dhydrometeors_div(k,nx,1)%moments(split_bins+1:max_nbins,imom))
                   else
                        ! need in case there's ice -ahu
                   endif
                   if (imom==1) field(k)=field(k)*QtoM3
                else
                   field(k)=sum(dhydrometeors_div(k,nx,ih)%moments(:,imom))
                endif
             end do
             call save_dg(field, name, i_dgtime,  units,dim=dims)
          end do
       end do
    end if

    endif
    
    if (l_save_mphys) then
    !-------------
    !
    ! Microphysics
    !
    !-------------
    !
    ! theta
    field(:)=dtheta_mphys(:,nx)
    call save_dg(field, 'dtheta_mphys', i_dgtime,  units='K/s',dim=dims)
    ! qv tendency
    field(:)=dqv_mphys(:,nx)
    call save_dg(field, 'dqv_mphys', i_dgtime,  units='kg/kg/s',dim=dims)

    !!aerosol
    !do ih=1,naerosol
    !   do imom=1,num_aero_moments(ih)
    !      do k=1,nz
    !         field(k)=sum(daerosol_mphys(k,nx,ih)%moments(:,imom))
    !      end do
    !      name=trim(aero_names(ih))//'_'//trim(aero_mom_names(imom))//'_mphys'
    !      units=trim(aero_mom_units(imom)//'/s')
    !      call save_dg(field, name, i_dgtime, units,dim=dims)
    !   end do
    !end do

    !hydrometeors
    do ih=1,nspecies
       do imom=1,num_h_moments(ih)
          name=trim(h_names(ih))//'_'//trim(mom_names(imom))//'_mphys'
          units=trim(mom_units(imom)//'/s')

          do k=1,nz
              if (ampORbin .eq. 'bin')then
                 if (ih==1) then
                     field(k)=sum(dhydrometeors_mphys(k,nx,1)%moments(1:split_bins,imom))
                 elseif (ih==2) then
                     field(k)=sum(dhydrometeors_mphys(k,nx,1)%moments(split_bins+1:max_nbins,imom))
                 else
                     ! need in case there's ice -ahu
                 endif
                 if (imom==1) field(k)=field(k)*QtoM3
              else
                 field(k)=sum(dhydrometeors_mphys(k,nx,ih)%moments(:,imom))
              endif
          end do
          call save_dg(field, name, i_dgtime,  units,dim=dims)
       end do
    end do
    endif
    !----------------
    !
    ! Imposed forcing
    !
    !----------------

    !! theta
    !field(:)=Tforce(:,nx)*exner(:,nx)
    !call save_dg(field, 'dtheta_force', i_dgtime,  units='K/s',dim=dims)
    !! qv
    !field(:)=qforce(:,nx)
    !call save_dg(field, 'dqv_force', i_dgtime,  units='kg/kg/s',dim=dims)
    !!hydrometeors
    !do ih=1,nspecies
    !   do imom=1,num_h_moments(ih)
    !      name=trim(h_names(ih))//'_'//trim(mom_names(imom))//'_force'
    !      units=trim(mom_units(imom))//'/s'
    !      do k=1,nz
    !          if (ampORbin .eq. 'bin')then
    !             if (ih==1) then
    !                 field(k)=sum(dhydrometeors_force(k,nx,1)%moments(1:split_bins,imom))*QtoM3
    !             elseif (ih==2) then
    !                 field(k)=sum(dhydrometeors_force(k,nx,1)%moments(split_bins+1:max_nbins,imom))*QtoM3
    !             else
    !                 ! need in case there's ice -ahu
    !             endif
    !          else
    !             field(k)=sum(dhydrometeors_force(k,nx,ih)%moments(:,imom))
    !          endif
    !      end do
    !      call save_dg(field, name, i_dgtime,  units,dim=dims)
    !   end do
    !end do
    endif ! l_save_tend

    !--------------
    !
    ! Miscellaneous
    !
    !--------------
    !Other diags
    !RH liquid
     do k=1,nz
        if (z(k) < 20000.)then
           field(k)= 100.*qv(k,nx)/ &
                qsaturation(TdegK(k,nx),pmb(k,nx))
        else
           field(k)=unset_real
        end if
     end do
     call save_dg(field(:)*field_mask(:,nx), 'RH', i_dgtime, '%', dims)

    !RH ice
     do k=1,nz
        if (z(k) < 20000.)then
           field(k)= 100.*qv(k,nx)/ &
                qisaturation(TdegK(k,nx),pmb(k,nx))
        else
           field(k)=unset_real
        end if
     end do
     call save_dg(field(:)*field_mask(:,nx), 'RH_ice', i_dgtime, '%', dims)



    !-------------------------
    !
    ! Instantaneous scalars
    !
    !-------------------------
    call save_dg(time, 'time', i_dgtime,  units='s',dim='time')

    !-------------------------
    !
    ! Column integrated values
    !
    !-------------------------

    !hydrometeors
    do ih=1,nspecies
       do imom=1,num_h_moments(ih)
          name=trim(h_names(ih))//'_'//trim(mom_names(imom))//'_path'
          units='m^k/m2'

          do k=1,nz
              if (ampORbin .eq. 'bin')then
                 if (ih==1) then
                     field(k)=sum(rho(k)*dz(k)*hydrometeors(k,nx,1)%moments(1:split_bins,imom))
                 elseif (ih==2) then
                     field(k)=sum(rho(k)*dz(k)*hydrometeors(k,nx,1)%moments(split_bins+1:max_nbins,imom))
                 endif
                 if (imom==1) field(k)=field(k)*QtoM3
              else
                  field(k)=sum(rho(k)*dz(k)*hydrometeors(k,nx,ih)%moments(:,imom))
              endif
          end do

          if (ih==2 .and. imom==1) then
            if (sum(field) /= sum(field)) then
              print*, field
              stop
            endif
          endif

          ! ! if (ih==2) print*, i_dgtime, sum(field)
          call save_dg(sum(field), name, i_dgtime,  units,dim='time')
          ! print*, 'path', ih, imom, sum(field)

         ! print*, 'ih, imom, path', ih, imom, field(25)
         ! if (sum(field)>0 .and. imom==4) then
         !   ! print*, hydrometeors(25,1,1)%moments(:,:)
         !   stop
         ! endif


       end do

    end do

    ! output extra moments


    if (ampORbin .eq. 'bin') then
      do ih = 1,nspecies
        field(:) = 0.
        do imom=3,4
          name=trim(h_names(ih))//'_'//trim(mom_names(imom))//'_path'
          units='m^k/m2'
          do k=1,nz
            if (ih==1) then
              m3_arr = hydrometeors(k,nx,1)%moments(1:split_bins,1)
              m0_arr = hydrometeors(k,nx,1)%moments(1:split_bins,2)
            elseif(ih==2) then
              m3_arr = hydrometeors(k,nx,1)%moments(split_bins+1:max_nbins,1)
              m0_arr = hydrometeors(k,nx,1)%moments(split_bins+1:max_nbins,2)
            endif
            if (ih==1) field(k)=rho(k)*dz(k)*momk(m3_arr, m0_arr, dble(pmomsc(imom)), 1)
            if (ih==2) field(k)=rho(k)*dz(k)*momk(m3_arr, m0_arr, dble(pmomsc(imom)), split_bins+1)

          enddo
          call save_dg(sum(field), name, i_dgtime,  units,dim='time')

          ! print*, 'path', ih, imom, sum(field)

         ! print*, 'ih, imom, layer', ih, imom, field(25)
         ! if (sum(field)>0 .and. imom==4) then
         !   stop
         ! endif

        enddo

      enddo

    endif

    deallocate(field_bin)
    deallocate(field_moms)
    deallocate(field)
    deallocate(field_nx)
    deallocate(field_bin_c)
    deallocate(field_bin_r)

    n_dgtimes=i_dgtime

    ! Write out progress indicator in 2% chunks.
    if ( mod( int(100*i_dgtime/float(maxn_dgtimes)), 20) < &
         mod( int(100*(i_dgtime-1)/float(maxn_dgtimes)), 20) ) &
         write (unit=6,fmt='(T3, i3, a)') int(100*i_dgtime&
         &/float(maxn_dgtimes)), '% completed...'

  end subroutine save_diagnostics_1d

  subroutine save_diagnostics_2d

    ! arrays for the 2-D diagnostics
    real(wp), allocatable :: field_nx(:)
    real(wp), allocatable :: field_2D(:,:)
    real(wp), allocatable :: field_bin_c_2d(:,:)
    real(wp), allocatable :: field_bin_r_2d(:,:)
    ! arrays to calculate the scalar diags
    real(wp), allocatable :: field(:)
    real(wp), allocatable :: field_bin_c(:)
    real(wp), allocatable :: field_bin_r(:)

    character(max_char_len) :: name, units, dims
    integer :: k, ih, imom, ibin

    if (.not. l_dgstep) return

    !
    ! Instantaneous grid diags for a 2-D run
    !
    allocate(field_nx(0:nx+1))
    allocate(field_2D(nz,1:nx))
    allocate(field_bin_c_2d(nz,1:nx))
    allocate(field_bin_r_2d(nz,1:nx))
    !
    ! Instantaneous horizontally averged diags for a scalar diags
    !
    allocate(field(nz))
    allocate(field_bin_c(nz))
    allocate(field_bin_r(nz))

    ! set dimensions for the diagnostic column or grid output, i.e.
    ! 2-D set-up is dimensioned with grid - 'z,x'

    dims = 'z,x'

    ! pressure
    field_2d(:,1:nx) = pmb(:,1:nx)
    call save_dg(field_2d, 'pressure', i_dgtime,  units='mb',dim=dims)
    !temperature(K)
    field_2d(:,1:nx) = TdegK(:,1:nx)
    call save_dg(field_2d, 'temperature', i_dgtime, units='K', dim=dims)
    !theta(K)
    field_2d(:,1:nx) = theta(:,1:nx)
    call save_dg(field_2d, 'theta', i_dgtime, units='K', dim=dims)
    !vapour
    field_2d(:,1:nx) = qv(:,1:nx)
    call save_dg(field_2d, 'vapor', i_dgtime, units='kg/kg', dim=dims)
    !vertical velocity
    field_2d(:,1:nx) = w(:,1:nx)
    call save_dg(field_2d, 'w', i_dgtime, units='m/s', dim=dims)
    !horizontal velocity
    field_2d(:,1:nx) = v(:,1:nx)
    call save_dg(field_2d, 'v_2d', i_dgtime, units='m/s', dim=dims)

    !!aerosol
    !do ih=1,naerosol
    !   do imom=1,num_aero_moments(ih)
    !      do k=1,nz
    !         do j=1,nx
    !            field_2D(k,j)=sum(aerosol(k,j,ih)%moments(:,imom))
    !         end do
    !      end do
    !      name=trim(aero_names(ih))//'_'//trim(aero_mom_names(imom))
    !      units=trim(aero_mom_units(imom))
    !      call save_dg(field_2D, name, i_dgtime, units,dim='z,x')
    !   end do
    !end do

    !hydrometeors

    do ih=1,nspecies
       do imom=1,num_h_moments(ih)
          name=trim(h_names(ih))//'_'//trim(mom_names(imom))
          units=trim(mom_units(imom))
          do k=1,nz
             do j=1,nx
                if (ampORbin .eq. 'bin') then
                   if (ih==1) then
                      field_2D(k,j)=sum(hydrometeors(k,j,1)%moments(1:split_bins,imom))
                   elseif (ih==2) then
                      field_2D(k,j)=sum(hydrometeors(k,j,1)%moments(split_bins+1:max_nbins,imom))
                   endif
                   if (imom==1) field_2d(k,j)=field_2d(k,j)*QtoM3
                else
                   field_2D(k,j)=sum(hydrometeors(k,j,ih)%moments(:,imom))
                endif
             end do
          end do
          call save_dg(field_2D(1:nz,1:nx), name, i_dgtime,  units,dim=dims)
       enddo
    enddo

    !========================================================
    ! Tendency terms
    !========================================================
    !
    ! Advection
    !
    !-----------
    ! theta
    field_2d(:,1:nx)=dtheta_adv(:,1:nx)
    call save_dg(field_2d, 'dtheta_adv', i_dgtime,  units='K/s',dim=dims)
    ! qv
    field_2d(:,1:nx)=dqv_adv(:,1:nx)
    call save_dg(field_2d, 'dqv_adv', i_dgtime,  units='kg/kg/s',dim=dims)

    ! ! aerosols
    ! do ih=1,naerosol
    !    do imom=1,num_aero_moments(ih)
    !       do k=1,nz
    !          do j=1,nx
    !             field_2d(k,j)=sum(daerosol_adv(k,j,ih)%moments(:,imom))
    !          end do
    !       end do
    !       name=trim(aero_names(ih))//'_'//trim(aero_mom_names(imom))//'_adv'
    !       units=trim(aero_mom_units(imom))//'/s'
    !       call save_dg(field_2D(1:nz,1:nx), name, i_dgtime,  units,dim='z,x')
    !    end do
    ! end do

    !hydrometeors
    do ih=1,nspecies
       do imom=1,num_h_moments(ih)
          name=trim(h_names(ih))//'_'//trim(mom_names(imom))//'_adv'
          units=trim(mom_units(imom))//'/s'
          do k=1,nz
             do j=1,nx
                if (ampORbin .eq. 'bin') then
                   if (ih==1) then
                      field_2D(k,j)=sum(dhydrometeors_adv(k,j,1)%moments(1:split_bins,imom))
                   elseif (ih==2) then
                      field_2D(k,j)=sum(dhydrometeors_adv(k,j,1)%moments(split_bins+1:max_nbins,imom))
                   endif
                   if (imom==1) field_2d(k,j)=field_2d(k,j)*QtoM3
                else
                   field_2D(k,j)=sum(dhydrometeors_adv(k,j,ih)%moments(:,imom))
                endif
             end do
          end do
          call save_dg(field_2d, name, i_dgtime,  units,dim=dims)
       end do
    end do
    !
    !-----------
    !
    ! Divergence
    !
    !-----------
    if (l_diverge)then
       ! theta
       field_2d(:,1:nx)=dtheta_div(:,1:nx)
       call save_dg(field_2d, 'dtheta_div', i_dgtime,  units='K/s',dim=dims)
       ! qv tendency
       field_2d(:,1:nx)=dqv_div(:,1:nx)
       call save_dg(field_2d, 'dqv_div', i_dgtime,  units='kg/kg/s',dim=dims)

       !!aerosol
       !do ih=1,naerosol
       !   do imom=1,num_aero_moments(ih)
       !      do k=1,nz
       !         do j=1,nx
       !            field_2D(k,j)=sum(daerosol_div(k,j,ih)%moments(:,imom))
       !         end do
       !      end do
       !      name=trim(aero_names(ih))//'_'//trim(aero_mom_names(imom))//'_div'
       !      units=trim(aero_mom_units(imom))//'/s'
       !      call save_dg(field_2D(1:nz,1:nx), name, i_dgtime,  units,dim='z,x')
       !   end do
       !end do

       !hydrometeors
       do ih=1,nspecies
          do imom=1,num_h_moments(ih)
             name=trim(h_names(ih))//'_'//trim(mom_names(imom))//'_div'
             units=trim(mom_units(imom)//'/s')
             do k=1,nz
                do j=1,nx
                   if (ampORbin .eq. 'bin') then
                      if (ih==1) then
                         field_2D(k,j)=sum(dhydrometeors_div(k,j,1)%moments(1:split_bins,imom))
                      elseif (ih==2) then
                         field_2D(k,j)=sum(dhydrometeors_div(k,j,1)%moments(split_bins+1:max_nbins,imom))
                      endif
                      if (imom==1) field_2d(k,j)=field_2d(k,j)*QtoM3
                   else
                      field_2D(k,j)=sum(dhydrometeors_div(k,j,ih)%moments(:,imom))
                   endif
                enddo
             end do
                call save_dg(field_2d, name, i_dgtime,  units,dim=dims)
          end do
       end do
    end if
    !-------------
    !
    ! Microphysics
    !
    !-------------
    !
    ! theta
    field_2d(:,1:nx)=dtheta_mphys(:,1:nx)
    call save_dg(field_2d, 'dtheta_mphys', i_dgtime,  units='K/s',dim='z')
    ! qv tendency
    field_2d(:,1:nx)=dqv_mphys(:,1:nx)
    call save_dg(field_2d, 'dqv_mphys', i_dgtime,  units='kg/kg/s',dim='z')

    !!aerosol
    !do ih=1,naerosol
    !   do imom=1,num_aero_moments(ih)
    !      do k=1,nz
    !         do j=1,nx
    !            field_2D(k,j)=sum(daerosol_mphys(k,j,ih)%moments(:,imom))
    !         end do
    !      end do
    !      name=trim(aero_names(ih))//'_'//trim(aero_mom_names(imom))//'_mphys'
    !      units=trim(aero_mom_units(imom)//'/s')
    !      call save_dg(field_2D(1:nz,1:nx), name, i_dgtime,  units,dim='z,x')
    !   end do
    !end do

    !hydrometeors
    do ih=1,nspecies
       do imom=1,num_h_moments(ih)
          name=trim(h_names(ih))//'_'//trim(mom_names(imom))//'_mphys'
          units=trim(mom_units(imom)//'/s')
          do k=1,nz
             do j=1,nx
                if (ampORbin .eq. 'bin') then
                   if (ih==1) then
                      field_2D(k,j)=sum(dhydrometeors_mphys(k,j,1)%moments(1:split_bins,imom))
                   elseif (ih==2) then
                      field_2D(k,j)=sum(dhydrometeors_mphys(k,j,1)%moments(split_bins+1:max_nbins,imom))
                   endif
                   if (imom==1) field_2d(k,j)=field_2d(k,j)*QtoM3
                else
                   field_2D(k,j)=sum(dhydrometeors_mphys(k,j,ih)%moments(:,imom))
                endif
             end do
           end do
           call save_dg(field_2d, name, i_dgtime,  units,dim=dims)
      end do
    end do
    !----------------
    !
    ! Imposed forcing
    !
    !----------------
    ! theta
    field_2D(:,1:nx)=Tforce(:,1:nx)*exner(:,1:nx)
    call save_dg(field_2d, 'dtheta_force', i_dgtime,  units='K/s',dim=dims)
    ! qv
    field_2D(:,1:nx)=qforce(:,1:nx)*exner(:,1:nx)
    call save_dg(field_2d, 'dqv_force', i_dgtime,  units='kg/kg/s',dim=dims)
    !hydrometeors
    do ih=1,nspecies
       do imom=1,num_h_moments(ih)
          name=trim(h_names(ih))//'_'//trim(mom_names(imom))//'_force'
          units=trim(mom_units(imom))//'/s'
          do k=1,nz
             do j=1,nx
                if (ampORbin .eq. 'bin') then
                   if (ih==1) then
                      field_2D(k,j)=sum(dhydrometeors_force(k,j,1)%moments(1:split_bins,imom))
                   elseif (ih==2) then
                      field_2D(k,j)=sum(dhydrometeors_force(k,j,1)%moments(split_bins+1:max_nbins,imom))
                   endif
                   if (imom==1) field_2d(k,j)=field_2d(k,j)*QtoM3
                else
                   field_2D(k,j)=sum(dhydrometeors_force(k,j,ih)%moments(:,imom))
                endif
             end do
           end do
           call save_dg(field_2d, name, i_dgtime,  units,dim=dims)
       end do
     end do
    !--------------
    !
    ! Miscellaneous
    !
    !--------------
    !Other diags
    !RH liquid
     do k=1,nz
        do j=1,nx
           if (z(k) < 20000.)then
              field_2D(k,j)= 100.*qv(k,j)/ &
                   qsaturation(TdegK(k,j),pmb(k,j))
           else
              field_2D(k,j)=unset_real
           end if
        end do
     enddo
     call save_dg(field_2D(1:nz,1:nx)*field_mask(1:nz,1:nx), 'RH', i_dgtime, units, dim='z,x')

    !RH ice
     do k=1,nz
        do j=1,nx
           if (z(k) < 20000.)then
              field_2D(k,j)= 100.*qv(k,j)/ &
                   qisaturation(TdegK(k,j),pmb(k,j))
           else
              field_2D(k,j)=unset_real
           end if
        end do
     enddo
     call save_dg(field_2D(1:nz,1:nx)*field_mask(1:nz,1:nx), 'RH_ice', i_dgtime, units, dim='z,x')

    !-------------------------
    !
    ! Instantaneous scalars
    !
    !-------------------------
    call save_dg(time, 'time', i_dgtime,  units='s',dim='time')

    !-------------------------
    !
    ! Horizontally averaged
    ! column integrated values
    !
    !-------------------------
    !hydrometeors
    do ih=1,nspecies
       do imom=1,num_h_moments(ih)
          name='mean_'//trim(h_names(ih))//'_'//trim(mom_names(imom))//'_path'
          units=trim(mom_units(imom))//' kg/m2'
          do j=1,nx
            do k=1,nz
              if (ampORbin .eq. 'bin')then
                 if (ih==1) then
                    field_2D(k,j)=sum(rho(k)*dz(k)*hydrometeors(k,j,1)%moments(1:split_bins,imom))
                 elseif (ih==2) then
                    field_2D(k,j)=sum(rho(k)*dz(k)*hydrometeors(k,j,1)%moments(split_bins+1:max_nbins,imom))
                 endif
                 if (imom==1) field_2d(k,j)=field_2d(k,j)*QtoM3
              else
                 field_2D(k,j)=sum(rho(k)*dz(k)*hydrometeors(k,j,ih)%moments(:,imom))
              endif
            end do
            field_nx(j) = sum(field_2D(:,j))
          end do
          ! horizontally averaged column integrated value
          call save_dg(sum(field_nx(1:nx))/nx, name, i_dgtime,  units,dim='time')
          ! column integrated values for each column
          name=trim(h_names(ih))//'_'//trim(mom_names(imom))//'_path'
          call save_dg(field_nx(1:nx), name, i_dgtime,  units,dim='x')
       end do
    end do

    deallocate(field_nx)
    deallocate(field_2D)
    deallocate(field_bin_c_2d)
    deallocate(field_bin_r_2d)
    deallocate(field)
    deallocate(field_bin_c)
    deallocate(field_bin_r)

    n_dgtimes=i_dgtime

    ! Write out progress indicator in 2% chunks.
    if ( mod( int(100*i_dgtime/float(maxn_dgtimes)), 2) < &
         mod( int(100*(i_dgtime-1)/float(maxn_dgtimes)), 2) ) &
         write (unit=6,fmt='(T3, i3, a)') int(100*i_dgtime&
         &/float(maxn_dgtimes)), '% completed...'

  end subroutine save_diagnostics_2d

  subroutine allocate_dgs_1DTS(dgStore)

    type(diag1DTS), intent(inout) :: dgStore
    integer :: n_offset, n_dgtimes

    n_offset = dgstart/dt-1
    n_dgtimes = n_times !- n_offset

    if (nx == 1) then

       maxn_dgtimes=max(maxn_dgtimes, int(n_dgtimes*dt/dg_dt)+1)
       allocate(dgStore%data(nz, maxn_dgtimes))
       dgStore%data=unset_real

    else

       maxn_dgtimes=max(maxn_dgtimes, int(n_dgtimes*dt/dg_dt)+1)
       allocate(dgStore%data(1:nx, maxn_dgtimes))
       dgStore%data=unset_real

    endif

  end subroutine allocate_dgs_1DTS

  subroutine allocate_dgs_2DTS(dgStore)

    type(diag2DTS), intent(inout) :: dgStore
    integer :: n_offset, n_dgtimes
    n_offset = dgstart/dt-1
    n_dgtimes = n_times !- n_offset

    maxn_dgtimes=max(maxn_dgtimes, int(n_dgtimes*dt/dg_dt)+1)
    allocate(dgStore%data(nz, 1:nx, maxn_dgtimes))
    dgStore%data=unset_real

  end subroutine allocate_dgs_2DTS

  subroutine allocate_dgs_bindgs(dgStore)

    type(diag_bin2dTS), intent(inout) :: dgStore
    integer :: n_offset, n_dgtimes

    n_offset = dgstart/dt-1
    n_dgtimes = n_times !- n_offset

    maxn_dgtimes=max(maxn_dgtimes, int(n_dgtimes*dt/dg_dt)+1)
    allocate(dgStore%data(nz, max_nbins, maxn_dgtimes))
    dgStore%data=unset_real

  end subroutine allocate_dgs_bindgs

  subroutine allocate_procdgs(dgStore)

    type(diag_proc2dTS), intent(inout) :: dgStore
    integer :: n_offset, n_dgtimes

    n_offset = dgstart/dt-1
    n_dgtimes = n_times !- n_offset

    maxn_dgtimes=max(maxn_dgtimes, int(n_dgtimes*dt/dg_dt)+1)
    allocate(dgStore%data(nz, nmom_diag, maxn_dgtimes))
    dgStore%data=unset_real

  end subroutine allocate_procdgs

  subroutine allocate_dgs_flagdgs(dgStore)

    type(diag_flag2DTS), intent(inout) :: dgStore
    integer :: n_offset, n_dgtimes

    n_offset = dgstart/dt-1
    n_dgtimes = n_times !- n_offset

    maxn_dgtimes=max(maxn_dgtimes, int(n_dgtimes*dt/dg_dt)+1)
    allocate(dgStore%data(nz, flag_count, maxn_dgtimes))
    dgStore%data=unset_real

  end subroutine allocate_dgs_flagdgs

  subroutine allocate_dgs_ScalarTS(dgStore)

    type(diagScalarTS), intent(inout) :: dgStore
    integer :: n_offset, n_dgtimes

    n_offset = dgstart/dt-1
    n_dgtimes = n_times !- n_offset
    maxn_dgtimes=max(maxn_dgtimes, int(n_dgtimes*dt/dg_dt)+1)
    allocate(dgStore%data(maxn_dgtimes))
    dgStore%data=unset_real
!print*,'EOA'
  end subroutine allocate_dgs_ScalarTS

  subroutine allocate_dgs_bindata(dgStore)

    type(diagBinData), intent(inout) :: dgStore

    allocate(dgStore%data(max_nbins))
    dgStore%data=unset_real

  end subroutine allocate_dgs_bindata

  subroutine getUniqueId(name, dg_index, ivar)

    character(*), intent(in) :: name
    type(dgIDarray), intent(inout) :: dg_index
    integer, intent(out) :: ivar

    integer :: i
    logical ::l_new ! logical to indicate a new variable

    ! First check to see if name already defined
    ! Start checking last accessed index and then
    ! loop around all other values (linear search)
    l_new=.true.
    do i=0,dg_index%nids-1
       ivar=mod( dg_index%lastAccessed + i - 1, dg_index%nids)+1
       if (trim(name)==trim(dg_index%dgIds(ivar)%name))then
          dg_index%lastAccessed=ivar
          l_new=.false.
          exit
       end if
    end do

    if (l_new)then
       ! If we haven't found it then assign the next free slot
       dg_index%nids=dg_index%nids+1
       ivar=dg_index%nids
       dg_index%dgIds(ivar)%name=name
       dg_index%lastAccessed=ivar
    end if

  end subroutine getUniqueId

  subroutine save_dg_1d_sp(field, name, itime, units, dim, longname)

    real(sp), intent(in) :: field(:)
    character(*), intent(in) :: name
    integer, intent(in) :: itime
    character(*), intent(in), optional :: units, dim, longname

    !local variables for column diagnostics
    type(diag1DTS), pointer :: dg(:)
    type(dgIDarray), pointer  :: dg_index

    character(max_char_len):: cunits, cdim, clongname
    integer :: ivar

    if (.not. l_dgstep) return

    ! We're assuming diagnostics are instant for now
    ! could put in an optional argument later to do
    ! averaged, accumulated, etc. later.
    dg=>instant_column
    dg_index=>ID_instant_column

    if (present(units))then
       cunits=units
    else
       cunits='Not set'
    end if

    if (present(dim))then
       cdim=dim
    else
       cdim='z'
    end if

    if (present(longname))then
       clongname=longname
    else
       clongname=name
    end if

    call getUniqueId(name, dg_index, ivar)

    if (.not. associated(dg(ivar)%data)) then
       call allocate_dgs(dg(ivar))
       dg(ivar)%name=name
       dg(ivar)%units=trim(cunits)
       dg(ivar)%dim=trim(cdim)
       dg(ivar)%longname=trim(clongname)
    end if

    dg(ivar)%data(:,itime)=field(:)

  end subroutine save_dg_1d_sp

  subroutine save_dg_1d_range_sp(ks, ke, field, name, itime, units, &
     dim, longname)

    integer, intent(in) :: ks, ke  ! start and end of range
    real(sp), intent(in) :: field(:)
    character(*), intent(in) :: name
    integer, intent(in) :: itime
    character(*), intent(in), optional :: units, dim, longname

    !local variables
    type(diag1DTS), pointer :: dg(:)
    type(dgIDarray), pointer  :: dg_index
    character(max_char_len):: cunits, cdim, clongname
    integer :: ivar

    if (.not. l_dgstep) return

    ! We're assuming diagnostics are instant for now
    ! could put in an optional argument later to do
    ! averaged, accumulated, etc. later.
    dg=>instant_column
    dg_index=>ID_instant_column


    if (present(units))then
       cunits=units
    else
       cunits='Not set'
    end if

    if (present(dim))then
       cdim=dim
    else
       cdim='z'
    end if

    if (present(longname))then
       clongname=longname
    else
       clongname=name
    end if

    call getUniqueId(name, dg_index, ivar)

    if (.not. associated(dg(ivar)%data)) then
       call allocate_dgs(dg(ivar))
       dg(ivar)%name=name
       dg(ivar)%units=trim(cunits)
       dg(ivar)%dim=trim(cdim)
       dg(ivar)%longname=trim(clongname)
    end if
    dg(ivar)%data(ks:ke,itime)=field(ks:ke)

  end subroutine save_dg_1d_range_sp

  subroutine save_dg_1d_point_sp(k, value, name, itime, units, &
       dim, longname)

    integer, intent(in) :: k ! index of value to input
    real(sp), intent(in) :: value
    character(*), intent(in) :: name
    integer, intent(in) :: itime
    character(*), intent(in), optional :: units, dim, longname

    !local variables
    type(diag1DTS), pointer :: dg(:)
    type(dgIDarray), pointer  :: dg_index
    character(max_char_len):: cunits, cdim, clongname
    integer :: ivar

    if (.not. l_dgstep) return

    ! We're assuming diagnostics are instant for now
    ! could put in an optional argument later to do
    ! averaged, accumulated, etc. later.
    dg=>instant_column
    dg_index=>ID_instant_column

    if (present(units))then
       cunits=units
    else
       cunits='Not set'
    end if

    if (present(dim))then
       cdim=dim
    else
       cdim='z'
    end if

    if (present(longname))then
       clongname=longname
    else
       clongname=name
    end if

    call getUniqueId(name, dg_index, ivar)

    if (.not. associated(dg(ivar)%data)) then
       call allocate_dgs(dg(ivar))
       dg(ivar)%name=name
       dg(ivar)%units=trim(cunits)
       dg(ivar)%dim=trim(cdim)
       dg(ivar)%longname=trim(clongname)
    end if
    dg(ivar)%data(k,itime)=value

  end subroutine save_dg_1d_point_sp

  !subroutine save_dg_bin_point_sp(bin, k, field, name, itime, units, &
  !   dim, longname)

  !  integer, intent(in) :: k ! index of value to input
  !  real(sp), intent(in) :: field(:)
  !  character(*), intent(in) :: name
  !  character(*), intent(in) :: bin
  !  integer, intent(in) :: itime
  !  character(*), intent(in), optional :: units, dim, longname

  !  !local variables
  !  type(diag_bin2DTS), pointer :: dg(:)
  !  type(dgIDarray), pointer  :: dg_index
  !  character(max_char_len):: cunits, cdim, clongname
  !  integer :: ivar

  !  if (bin=='bin') then
  !    if (.not. l_dgstep) return

  !    ! We're assuming diagnostics are instant for now
  !    ! could put in an optional argument later to do
  !    ! averaged, accumulated, etc. later.
  !    dg=>instant_bindgs
  !    dg_index=>ID_instant_bindgs

  !    if (present(units))then
  !      cunits=units
  !    else
  !      cunits='Not set'
  !    end if

  !    if (present(dim))then
  !      cdim=dim
  !    else
  !      cdim='z'
  !    end if

  !    if (present(longname))then
  !      clongname=longname
  !    else
  !      clongname=name
  !    end if

  !    call getUniqueId(name, dg_index, ivar)

  !    if (.not. associated(dg(ivar)%data)) then
  !      call allocate_dgs(dg(ivar))
  !      dg(ivar)%name=name
  !      dg(ivar)%units=trim(cunits)
  !      dg(ivar)%dim=trim(cdim)
  !      dg(ivar)%longname=trim(clongname)
  !    end if

  !    dg(ivar)%data(k,:,itime)=field(:)
  !  end if

  !end subroutine save_dg_bin_point_sp

  !subroutine save_dg_bin_point_dp(bin, k, field, name, itime, units, &
  !   dim, longname)

  !  integer, intent(in) :: k ! index of value to input
  !  real(dp), intent(in) :: field(:)
  !  character(*), intent(in) :: name
  !  character(*), intent(in) :: bin
  !  integer, intent(in) :: itime
  !  character(*), intent(in), optional :: units, dim, longname

  !  !local variables
  !  type(diag_bin2DTS), pointer :: dg(:)
  !  type(dgIDarray), pointer  :: dg_index
  !  character(max_char_len):: cunits, cdim, clongname
  !  integer :: ivar

  !  if (bin=='bin') then
  !    if (.not. l_dgstep) return

  !    ! We're assuming diagnostics are instant for now
  !    ! could put in an optional argument later to do
  !    ! averaged, accumulated, etc. later.
  !    dg=>instant_bindgs
  !    dg_index=>ID_instant_bindgs

  !    if (present(units))then
  !      cunits=units
  !    else
  !      cunits='Not set'
  !    end if

  !    if (present(dim))then
  !      cdim=dim
  !    else
  !      cdim='z'
  !    end if

  !    if (present(longname))then
  !      clongname=longname
  !    else
  !      clongname=name
  !    end if

  !    call getUniqueId(name, dg_index, ivar)

  !    if (.not. associated(dg(ivar)%data)) then
  !      call allocate_dgs(dg(ivar))
  !      dg(ivar)%name=name
  !      dg(ivar)%units=trim(cunits)
  !      dg(ivar)%dim=trim(cdim)
  !      dg(ivar)%longname=trim(clongname)
  !    end if

  !    dg(ivar)%data(k,:,itime)=field(:)
  !  end if

  !end subroutine save_dg_bin_point_dp

  subroutine save_dg_2d_sp(field, name, itime, units, &
       dim, longname)

    real(sp), intent(in) :: field(:,:)
    character(*), intent(in) :: name
    integer, intent(in) :: itime
    character(*), intent(in), optional :: units, dim, longname

    !local variables
    type(diag2DTS), pointer :: dg(:)
    type(dgIDarray), pointer  :: dg_index
    character(max_char_len):: cunits, cdim, clongname
    integer :: ivar

    if (.not. l_dgstep) return

    ! We're assuming diagnostics are instant for now
    ! could put in an optional argument later to do
    ! averaged, accumulated, etc. later.
    dg=>instant_2D
    dg_index=>ID_instant_2D

    if (present(units))then
       cunits=units
    else
       cunits='Not set'
    end if

    if (present(dim))then
       cdim=dim
    else
       cdim='z'
    end if

    if (present(longname))then
       clongname=longname
    else
       clongname=name
    end if

    call getUniqueId(name, dg_index, ivar)

    if (.not. associated(dg(ivar)%data)) then
       call allocate_dgs(dg(ivar))
       dg(ivar)%name=name
       dg(ivar)%units=trim(cunits)
       dg(ivar)%dim=trim(cdim)
       dg(ivar)%longname=trim(clongname)
    end if

    dg(ivar)%data(:,:,itime)=field(:,:)

  end subroutine save_dg_2d_sp

  subroutine save_dg_2d_point_sp(k, l, value, name, itime, units, &
       dim, longname)

    integer, intent(in) :: k ,l! indices of value to input
    real(sp), intent(in) :: value
    character(*), intent(in) :: name
    integer, intent(in) :: itime
    character(*), intent(in), optional :: units, dim, longname

    !local variables
    type(diag2DTS), pointer :: dg(:)
    type(dgIDarray), pointer  :: dg_index
    character(max_char_len):: cunits, cdim, clongname
    integer :: ivar

    if (.not. l_dgstep) return

    ! We're assuming diagnostics are instant for now
    ! could put in an optional argument later to do
    ! averaged, accumulated, etc. later.
    dg=>instant_2D
    dg_index=>ID_instant_2D

    if (present(units))then
       cunits=units
    else
       cunits='Not set'
    end if

    if (present(dim))then
       cdim=dim
    else
       cdim='z'
    end if

    if (present(longname))then
       clongname=longname
    else
       clongname=name
    end if

    call getUniqueId(name, dg_index, ivar)

    if (.not. associated(dg(ivar)%data)) then
       call allocate_dgs(dg(ivar))
       dg(ivar)%name=name
       dg(ivar)%units=trim(cunits)
       dg(ivar)%dim=trim(cdim)
       dg(ivar)%longname=trim(clongname)
    end if
    dg(ivar)%data(k,l,itime)=value

  end subroutine save_dg_2d_point_sp

  subroutine save_dg_scalar_sp(value, name, itime, units, dim, longname)

    real(sp), intent(in) :: value
    character(*), intent(in) :: name
    integer, intent(in) :: itime
    character(*), intent(in), optional :: units, dim, longname

    !local variables
    type(diagScalarTS), pointer:: dg(:)
    type(dgIDarray), pointer :: dg_index
    character(max_char_len):: cunits, cdim, clongname
    integer :: ivar

    if (.not. l_dgstep) return

    dg=>scalars
    dg_index=>ID_scalars

    if (present(units))then
       cunits=units
    else
       cunits='Not set'
    end if

    if (present(dim))then
       cdim=dim
    else
       cdim='time'
    end if

    if (present(longname))then
       clongname=longname
    else
       clongname=name
    end if

    call getUniqueId(name, dg_index, ivar)

    if (.not. associated(dg(ivar)%data)) then
       call allocate_dgs(dg(ivar))
       dg(ivar)%name=name
       dg(ivar)%units=trim(cunits)
       dg(ivar)%dim=trim(cdim)
       dg(ivar)%longname=trim(clongname)
    end if
    dg(ivar)%data(itime)=value

  end subroutine save_dg_scalar_sp

  subroutine save_dg_1d_dp(field, name, itime, units, dim, longname)

    real(dp), intent(in) :: field(:)
    character(*), intent(in) :: name
    integer, intent(in) :: itime
    character(*), intent(in), optional :: units, dim, longname

    !local variables
    type(diag1DTS), pointer :: dg(:)
    type(dgIDarray), pointer  :: dg_index
    character(max_char_len):: cunits, cdim, clongname
    integer :: ivar


    if (.not. l_dgstep) return

    ! We're assuming diagnostics are instant for now
    ! could put in an optional argument later to do
    ! averaged, accumulated, etc. later.
    dg=>instant_column
    dg_index=>ID_instant_column

    if (present(units))then
       cunits=units
    else
       cunits='Not set'
    end if

    if (present(dim))then
       cdim=dim
    else
       cdim='z'
    end if

    if (present(longname))then
       clongname=longname
    else
       clongname=name
    end if

    call getUniqueId(name, dg_index, ivar)

    if (.not. associated(dg(ivar)%data)) then
       call allocate_dgs(dg(ivar))
       dg(ivar)%name=name
       dg(ivar)%units=trim(cunits)
       dg(ivar)%dim=trim(cdim)
       dg(ivar)%longname=trim(clongname)
    end if
    dg(ivar)%data(:,itime)=field(:)

  end subroutine save_dg_1d_dp

  subroutine save_dg_1d_range_dp(ks, ke, field, name, itime, units, &
       dim, longname)

    integer, intent(in) :: ks, ke  ! start and end of range
    real(dp), intent(in) :: field(:)
    character(*), intent(in) :: name
    integer, intent(in) :: itime
    character(*), intent(in), optional :: units, dim, longname

    !local variables
    type(diag1DTS), pointer :: dg(:)
    type(dgIDarray), pointer  :: dg_index
    character(max_char_len):: cunits, cdim, clongname
    integer :: ivar

    if (.not. l_dgstep) return

    ! We're assuming diagnostics are instant for now
    ! could put in an optional argument later to do
    ! averaged, accumulated, etc. later.
    dg=>instant_column
    dg_index=>ID_instant_column


    if (present(units))then
       cunits=units
    else
       cunits='Not set'
    end if

    if (present(dim))then
       cdim=dim
    else
       cdim='z'
    end if

    if (present(longname))then
       clongname=longname
    else
       clongname=name
    end if

    call getUniqueId(name, dg_index, ivar)

    if (.not. associated(dg(ivar)%data)) then
       call allocate_dgs(dg(ivar))
       dg(ivar)%name=name
       dg(ivar)%units=trim(cunits)
       dg(ivar)%dim=trim(cdim)
       dg(ivar)%longname=trim(clongname)
    end if
    dg(ivar)%data(ks:ke,itime)=field(ks:ke)

  end subroutine save_dg_1d_range_dp

  subroutine save_dg_1d_point_dp(k, value, name, itime, units, &
       dim, longname)

    integer, intent(in) :: k ! index of value to input
    real(dp), intent(in) :: value
    character(*), intent(in) :: name
    integer, intent(in) :: itime
    character(*), intent(in), optional :: units, dim, longname

    !local variables
    type(diag1DTS), pointer :: dg(:)
    type(dgIDarray), pointer  :: dg_index
    character(max_char_len):: cunits, cdim, clongname
    integer :: ivar

    if (.not. l_dgstep) return

    ! We're assuming diagnostics are instant for now
    ! could put in an optional argument later to do
    ! averaged, accumulated, etc. later.
    dg=>instant_column
    dg_index=>ID_instant_column

    if (present(units))then
       cunits=units
    else
       cunits='Not set'
    end if

    if (present(dim))then
       cdim=dim
    else
       cdim='z'
    end if

    if (present(longname))then
       clongname=longname
    else
       clongname=name
    end if

    call getUniqueId(name, dg_index, ivar)

    if (.not. associated(dg(ivar)%data)) then
       call allocate_dgs(dg(ivar))
       dg(ivar)%name=name
       dg(ivar)%units=trim(cunits)
       dg(ivar)%dim=trim(cdim)
       dg(ivar)%longname=trim(clongname)
    end if
    dg(ivar)%data(k,itime)=value

  end subroutine save_dg_1d_point_dp

  subroutine save_dg_2d_dp(field, name, itime, units, &
       dim, longname)

    real(dp), intent(in) :: field(:,:)
    character(*), intent(in) :: name
    integer, intent(in) :: itime
    character(*), intent(in), optional :: units, dim, longname

    !local variables
    type(diag2DTS), pointer :: dg(:)
    type(dgIDarray), pointer  :: dg_index
    character(max_char_len):: cunits, cdim, clongname
    integer :: ivar

    if (.not. l_dgstep) return

    ! We're assuming diagnostics are instant for now
    ! could put in an optional argument later to do
    ! averaged, accumulated, etc. later.
    dg=>instant_2D
    dg_index=>ID_instant_2D

    if (present(units))then
       cunits=units
    else
       cunits='Not set'
    end if

    if (present(dim))then
       cdim=dim
    else
       cdim='z'
    end if

    if (present(longname))then
       clongname=longname
    else
       clongname=name
    end if

    call getUniqueId(name, dg_index, ivar)

    if (.not. associated(dg(ivar)%data)) then
       call allocate_dgs(dg(ivar))
       dg(ivar)%name=name
       dg(ivar)%units=trim(cunits)
       dg(ivar)%dim=trim(cdim)
       dg(ivar)%longname=trim(clongname)
    end if
    dg(ivar)%data(:,:,itime)=field(:,:)

  end subroutine save_dg_2d_dp

  subroutine save_dg_2d_point_dp(k, l, value, name, itime, units, &
       dim, longname)

    integer, intent(in) :: k ,l! indices of value to input
    real(dp), intent(in) :: value
    character(*), intent(in) :: name
    integer, intent(in) :: itime
    character(*), intent(in), optional :: units, dim, longname

    !local variables
    type(diag2DTS), pointer :: dg(:)
    type(dgIDarray), pointer  :: dg_index
    character(max_char_len):: cunits, cdim, clongname
    integer :: ivar

    if (.not. l_dgstep) return

    ! We're assuming diagnostics are instant for now
    ! could put in an optional argument later to do
    ! averaged, accumulated, etc. later.
    dg=>instant_2D
    dg_index=>ID_instant_2D

    if (present(units))then
       cunits=units
    else
       cunits='Not set'
    end if

    if (present(dim))then
       cdim=dim
    else
       cdim='z'
    end if

    if (present(longname))then
       clongname=longname
    else
       clongname=name
    end if

    call getUniqueId(name, dg_index, ivar)

    if (.not. associated(dg(ivar)%data)) then
       call allocate_dgs(dg(ivar))
       dg(ivar)%name=name
       dg(ivar)%units=trim(cunits)
       dg(ivar)%dim=trim(cdim)
       dg(ivar)%longname=trim(clongname)
    end if
    dg(ivar)%data(k,l,itime)=value

  end subroutine save_dg_2d_point_dp


  subroutine save_dg_scalar_dp(value, name, itime, units, dim, longname)

    real(dp), intent(in) :: value
    character(*), intent(in) :: name
    integer, intent(in) :: itime
    character(*), intent(in), optional :: units, dim, longname

    !local variables
    type(diagScalarTS), pointer:: dg(:)
    type(dgIDarray), pointer :: dg_index
    character(max_char_len):: cunits, cdim, clongname
    integer :: ivar

    if (.not. l_dgstep) return

    dg=>scalars
    dg_index=>ID_scalars

    if (present(units))then
       cunits=units
    else
       cunits='Not set'
    end if

    if (present(dim))then
       cdim=dim
    else
       cdim='time'
    end if

    if (present(longname))then
       clongname=longname
    else
       clongname=name
    end if

    call getUniqueId(name, dg_index, ivar)

    if (.not. associated(dg(ivar)%data)) then
       call allocate_dgs(dg(ivar))
       dg(ivar)%name=name
       dg(ivar)%units=trim(cunits)
       dg(ivar)%dim=trim(cdim)
       dg(ivar)%longname=trim(clongname)
    end if
    dg(ivar)%data(itime)=value

  end subroutine save_dg_scalar_dp

  subroutine save_dg_bin_sp(bin, field, name, itime, units, &
     dim, longname)

    real(sp), intent(in) :: field(:,:)
    character(*), intent(in) :: name
    character(*), intent(in) :: bin
    integer, intent(in) :: itime
    character(*), intent(in), optional :: units, dim, longname

    !local variables
    type(diag_bin2DTS), pointer :: dg(:)
    type(dgIDarray), pointer  :: dg_index
    character(max_char_len):: cunits, cdim, clongname
    integer :: ivar

    if (bin=='bin') then
      if (.not. l_dgstep) return

      ! We're assuming diagnostics are instant for now
      ! could put in an optional argument later to do
      ! averaged, accumulated, etc. later.
      dg=>instant_bindgs
      dg_index=>ID_instant_bindgs

      if (present(units))then
        cunits=units
      else
        cunits='Not set'
      end if

      if (present(dim))then
        cdim=dim
      else
        cdim='z'
      end if

      if (present(longname))then
        clongname=longname
      else
        clongname=name
      end if

      call getUniqueId(name, dg_index, ivar)

      if (.not. associated(dg(ivar)%data)) then
        call allocate_dgs(dg(ivar))
        dg(ivar)%name=name
        dg(ivar)%units=trim(cunits)
        dg(ivar)%dim=trim(cdim)
        dg(ivar)%longname=trim(clongname)
      end if

      dg(ivar)%data(:,:,itime)=field(:,:)
    end if

  end subroutine save_dg_bin_sp

  subroutine save_dg_bin_dp(bin, field, name, itime, units, &
     dim, longname)

    real(dp), intent(in) :: field(:,:)
    character(*), intent(in) :: name
    character(*), intent(in) :: bin
    integer, intent(in) :: itime
    character(*), intent(in), optional :: units, dim, longname

    !local variables
    type(diag_bin2DTS), pointer :: dg(:)
    type(dgIDarray), pointer  :: dg_index
    character(max_char_len):: cunits, cdim, clongname
    integer :: ivar

   if (bin=='bin') then
      if (.not. l_dgstep) return
      ! print*, 'dp called bin'
      ! We're assuming diagnostics are instant for now
      ! could put in an optional argument later to do
      ! averaged, accumulated, etc. later.
      dg=>instant_bindgs
      dg_index=>ID_instant_bindgs

      if (present(units))then
        cunits=units
      else
        cunits='Not set'
      end if

      if (present(dim))then
        cdim=dim
      else
        cdim='z'
      end if

      if (present(longname))then
        clongname=longname
      else
        clongname=name
      end if

      call getUniqueId(name, dg_index, ivar)

      if (.not. associated(dg(ivar)%data)) then
        call allocate_dgs(dg(ivar))
        dg(ivar)%name=name
        dg(ivar)%units=trim(cunits)
        dg(ivar)%dim=trim(cdim)
        dg(ivar)%longname=trim(clongname)
      end if

      dg(ivar)%data(:,:,itime)=field(:,:)
  end if

  end subroutine save_dg_bin_dp

  subroutine save_proc_dp(field, name, itime, units, &
     dim, longname)

    real(dp), intent(in) :: field(:,:)
    character(*), intent(in) :: name
    integer, intent(in) :: itime
    character(*), intent(in), optional :: units, dim, longname

    !local variables
    type(diag_proc2DTS), pointer :: dg(:)
    type(dgIDarray), pointer  :: dg_index
    character(max_char_len):: cunits, cdim, clongname
    integer :: ivar

    if (.not. l_dgstep) return
    ! print*, 'dp called nmom'
    ! We're assuming diagnostics are instant for now
    ! could put in an optional argument later to do
    ! averaged, accumulated, etc. later.
    dg=>instant_procdgs
    dg_index=>ID_instant_procdgs

    if (present(units))then
      cunits=units
    else
      cunits='Not set'
    end if

    if (present(dim))then
      cdim=dim
    else
      cdim='z'
    end if

    if (present(longname))then
      clongname=longname
    else
      clongname=name
    end if

    call getUniqueId(name, dg_index, ivar)

    if (.not. associated(dg(ivar)%data)) then
      call allocate_procdgs(dg(ivar))
      dg(ivar)%name=name
      dg(ivar)%units=trim(cunits)
      dg(ivar)%dim=trim(cdim)
      dg(ivar)%longname=trim(clongname)
    end if

    dg(ivar)%data(:,:,itime)=field(:,:)

  end subroutine save_proc_dp

  !subroutine save_dg_flag_sp(field, flag, name, itime, units, &
  !   dim, longname)

  !  real(sp), intent(in) :: field(:,:)
  !  character(*), intent(in) :: name
  !  character(*), intent(in) :: flag
  !  integer, intent(in) :: itime
  !  character(*), intent(in), optional :: units, dim, longname

  !  !local variables
  !  type(diag_flag2DTS), pointer :: dg(:)
  !  type(dgIDarray), pointer  :: dg_index
  !  character(max_char_len):: cunits, cdim, clongname
  !  integer :: ivar

  !  if (flag=='flag') then
  !    if (.not. l_dgstep) return

  !    dg=>instant_flagdgs
  !    dg_index=>ID_instant_flagdgs

  !    if (present(units))then
  !      cunits=units
  !    else
  !      cunits='Not set'
  !    end if

  !    if (present(dim))then
  !      cdim=dim
  !    else
  !      cdim='z'
  !    end if

  !    if (present(longname))then
  !      clongname=longname
  !    else
  !      clongname=name
  !    end if

  !    call getUniqueId(name, dg_index, ivar)

  !    if (.not. associated(dg(ivar)%data)) then
  !      call allocate_dgs(dg(ivar))
  !      dg(ivar)%name=name
  !      dg(ivar)%units=trim(cunits)
  !      dg(ivar)%dim=trim(cdim)
  !      dg(ivar)%longname=trim(clongname)
  !    end if

  !    dg(ivar)%data(:,:,itime)=field(:,:)
  !  end if

  !end subroutine save_dg_flag_sp

  !subroutine save_dg_flag_dp(field, flag, name, itime, units, &
  !   dim, longname)

  !  real(dp), intent(in) :: field(:,:)
  !  character(*), intent(in) :: name
  !  character(*), intent(in) :: flag
  !  integer, intent(in) :: itime
  !  character(*), intent(in), optional :: units, dim, longname

  !  !local variables
  !  type(diag_flag2DTS), pointer :: dg(:)
  !  type(dgIDarray), pointer  :: dg_index
  !  character(max_char_len):: cunits, cdim, clongname
  !  integer :: ivar

  !  if (flag=='flag') then
  !    if (.not. l_dgstep) return

  !    dg=>instant_flagdgs
  !    dg_index=>ID_instant_flagdgs

  !    if (present(units))then
  !      cunits=units
  !    else
  !      cunits='Not set'
  !    end if

  !    if (present(dim))then
  !      cdim=dim
  !    else
  !      cdim='z'
  !    end if

  !    if (present(longname))then
  !      clongname=longname
  !    else
  !      clongname=name
  !    end if

  !    call getUniqueId(name, dg_index, ivar)

  !    if (.not. associated(dg(ivar)%data)) then
  !      call allocate_dgs(dg(ivar))
  !      dg(ivar)%name=name
  !      dg(ivar)%units=trim(cunits)
  !      dg(ivar)%dim=trim(cdim)
  !      dg(ivar)%longname=trim(clongname)
  !    end if

  !    dg(ivar)%data(:,:,itime)=field(:,:)
  !  end if

  !end subroutine save_dg_flag_dp

  subroutine save_binData_sp(field, name, units, &
       longname)
    ! Save the bin data (e.g. sizes or masses)
    real(sp), intent(in) :: field(:)
    character(*), intent(in) :: name
    character(*), intent(in), optional :: units, longname

    !local variables
    type(diagBinData), pointer :: dg(:)
    type(dgIDarray), pointer  :: dg_index
    character(max_char_len):: cunits, cdim, clongname
    integer :: ivar

    dg=>binData
    dg_index=>ID_binData

    if (present(units))then
       cunits=units
    else
       cunits='Not set'
    end if

    cdim='nbins'

    if (present(longname))then
       clongname=longname
    else
       clongname=name
    end if

    call getUniqueId(name, dg_index, ivar)

    if (.not. associated(dg(ivar)%data)) then
       call allocate_dgs(dg(ivar))
       dg(ivar)%name=name
       dg(ivar)%units=trim(cunits)
       dg(ivar)%dim=trim(cdim)
       dg(ivar)%longname=trim(clongname)
    end if
    dg(ivar)%data(:)=field(:)

  end subroutine save_binData_sp

  subroutine save_binData_dp(field, name, units, &
       longname)
    ! Save the bin data (e.g. sizes or masses)
    real(dp), intent(in) :: field(:)
    character(*), intent(in) :: name
    character(*), intent(in), optional :: units, longname

    !local variables
    type(diagBinData), pointer :: dg(:)
    type(dgIDarray), pointer  :: dg_index
    character(max_char_len):: cunits, cdim, clongname
    integer :: ivar

    dg=>binData
    dg_index=>ID_binData

    if (present(units))then
       cunits=units
    else
       cunits='Not set'
    end if

    cdim='nbins'

    if (present(longname))then
       clongname=longname
    else
       clongname=name
    end if

    call getUniqueId(name, dg_index, ivar)

    if (.not. associated(dg(ivar)%data)) then
       call allocate_dgs(dg(ivar))
       dg(ivar)%name=name
       dg(ivar)%units=trim(cunits)
       dg(ivar)%dim=trim(cdim)
       dg(ivar)%longname=trim(clongname)
    end if
    dg(ivar)%data(:)=field(:)

  end subroutine save_binData_dp

  subroutine my_save_dg_bin_dp(bin, field, name, itime, units, &
     dim, longname)

    real(dp), intent(in) :: field(:,:)
    character(*), intent(in) :: name
    character(*), intent(in) :: bin
    integer, intent(in) :: itime
    character(*), intent(in), optional :: units, dim, longname

    !local variables
    type(diag_bin2DTS), pointer :: dg(:)
    type(dgIDarray), pointer  :: dg_index
    character(max_char_len):: cunits, cdim, clongname
    integer :: ivar

    if (.not. l_dgstep) return

    ! We're assuming diagnostics are instant for now
    ! could put in an optional argument later to do
    ! averaged, accumulated, etc. later.
    dg=>instant_bindgs
    dg_index=>ID_instant_bindgs

    if (present(units))then
       cunits=units
    else
       cunits='Not set'
    end if

    if (present(dim))then
       cdim=dim
    else
       cdim='z'
    end if

    if (present(longname))then
       clongname=longname
    else
       clongname=name
    end if

    call getUniqueId(name, dg_index, ivar)

    if (.not. associated(dg(ivar)%data)) then
       call allocate_dgs(dg(ivar))
       dg(ivar)%name=name
       dg(ivar)%units=trim(cunits)
       dg(ivar)%dim=trim(cdim)
       dg(ivar)%longname=trim(clongname)
    end if

    dg(ivar)%data(:,:,itime)=field(:,:)

  end subroutine my_save_dg_bin_dp

  subroutine write_diagnostics

    Use netcdf
    Use switches, only : icase
    Use namelists

    character(max_char_len) :: outfile, nmlfile, outdir, name

    ! Local variables
    integer :: ivar, itmp=1, offset=3, n_dgs, k, iq, tim, imom
    character(4) :: char4

    ! netcdf variables
    integer :: status, ncid, zid, xid, timeid, binsid, flagsid, nmomid
    integer :: tzdim(2), txdim(2), tzldim(3), tzxdim(3), tzfdim(3), dim, &
               dim2d(2), dim3d(3), tzmdim(3)
    integer, allocatable :: varid(:)
    real(8), allocatable :: field3d(:,:,:)

    write(6,*) 'Integration complete.'

    if (fileNameOut=='')then
      write(char4,'(I4.4)')icase
      outdir='output/'
      if (trim(KiD_outdir)/='')outdir=KiD_outdir
      if (trim(KiD_outfile)/='') then
        outfile = trim(outdir)//trim(KiD_outfile)
      else
        outfile=trim(outdir)//trim(modelName)//'_m-'//trim(mphys_id)// &
          '_b-'//trim(ampORbin)//'+'//trim(bintype)// '_u-'// &
          trim(username)//'_c-'//char4//'_v-0001.nc'
      endif
      offset=10
      call sanitize(outfile)
      status=nf90_eexist
      do while(status==nf90_eexist)
        ! set to "overwrite" mode for convenience -ahu
        status=nf90_create(outfile, nf90_clobber, ncid)
        ! status=nf90_create(outfile, nf90_noclobber, ncid)
        if (status==nf90_eexist)then
          write(char4,'(I4.4)')itmp
        else
          exit
        end if
        itmp=itmp+1
        outfile=trim(outfile(1:len_trim(outfile)-offset))//'_v-'//char4//'.nc'
      end do
    else
      outfile=trim(fileNameOut)
      status=nf90_create(outfile, nf90_clobber, ncid)
    end if

    write(6,*) 'Writing data to: ', trim(outfile)


    !------------------------
    ! Write global attributes
    !------------------------

    call date_and_time(zone=dt_zone,values=dt_values)

    write ( dateString, fmt )  'Date ',dt_values(3), '/', dt_values(2), '/',&
         & dt_values(1), '; time ', dt_values(5), ':', dt_values(6),&
         & ':', dt_values(7), dt_zone,' UTC'

    status=nf90_put_att(ncid, nf90_global, 'title',              &
         ' Microphysics dataset produced by '//               &
         trim(modelName)//' model version '//trim(version)// &
         '('//trim(revision)//').')

    status=nf90_put_att(ncid, nf90_global, 'creation date',      &
         dateString)

    status=nf90_put_att(ncid, nf90_global, 'User', username)
    status=nf90_put_att(ncid, nf90_global, 'Institution', institution)
    status=nf90_put_att(ncid, nf90_global, 'Microphysics ID', mphys_id)
    status=nf90_put_att(ncid, nf90_global, 'Advection ID', advection_id)
    status=nf90_put_att(ncid, nf90_global, 'references', references)
    status=nf90_put_att(ncid, nf90_global, 'comments', comments)
    status=nf90_put_att(ncid, nf90_global, 'Na', aero_N_init/1e6)
    status=nf90_put_att(ncid, nf90_global, 'w', wctrl(1))

    if (Dm_init>0.) then
      status=nf90_put_att(ncid, nf90_global, 'dm', Dm_init)
    endif
    if (Nd_init>0.) then
      status=nf90_put_att(ncid, nf90_global, 'Nd', Nd_init)
    endif
    if (rain_source(1)>0.) then
      status=nf90_put_att(ncid, nf90_global, 'dm', rain_source(1))
    endif

    status=nf90_def_dim(ncid, 'time', int(n_dgtimes, kind=incdfp), timeid)
    call check_ncstatus(status)
    status=nf90_def_dim(ncid, 'z', int(nz, kind=incdfp), zid)
    call check_ncstatus(status)
    status=nf90_def_dim(ncid, 'nmom', int(nmom_diag, kind=incdfp), nmomid)
    call check_ncstatus(status)
    status=nf90_def_dim(ncid, 'bins', int(max_nbins, kind=incdfp), binsid)
    call check_ncstatus(status)
    status=nf90_def_dim(ncid, 'x', int(nx, kind=incdfp), xid)
    call check_ncstatus(status)
    status=nf90_def_dim(ncid, 'flags', int(flag_count, kind=incdfp), flagsid)
    call check_ncstatus(status)

    tzdim=(/ timeid, zid /)
    txdim=(/ timeid, xid /)
    tzldim=(/ timeid, binsid, zid /)
    tzmdim=(/ timeid, nmomid, zid /)
    tzxdim=(/ timeid, xid, zid /)
    tzfdim=(/timeid, flagsid, zid/)

    ! Do the dim variables
    n_dgs=2
    allocate(varid(n_dgs))
    ivar=1
    status=nf90_def_var(ncid, 'z', NF90_DOUBLE, zid, varid(ivar))
    call check_ncstatus(status)
    status=nf90_put_att(ncid, varid(ivar),'units', 'm')
    call check_ncstatus(status)
    ivar=2
    status=nf90_def_var(ncid, 'x', NF90_DOUBLE, xid, varid(ivar))
    call check_ncstatus(status)
    status=nf90_put_att(ncid, varid(ivar),'units', 'm')
    call check_ncstatus(status)

    status=nf90_enddef(ncid)

    ivar=1
    status=nf90_put_var(ncid, varid(ivar), z)
    call check_ncstatus(status)

    ivar=2
    if (nx==1) x(:)=0.
    status=nf90_put_var(ncid, varid(ivar), x(1:nx))
    call check_ncstatus(status)

    deallocate(varid)

    ! Do the bin data variables
    n_dgs=ID_binData%nids
    allocate(varid(n_dgs))
    status=nf90_redef(ncid)
    dim=binsid

    do ivar=1,n_dgs
       name=binData(ivar)%name
       call sanitize(name)
       status=nf90_def_var(ncid, name &
            , NF90_DOUBLE, dim, varid(ivar))
       call check_ncstatus(status, binData(ivar)%name)
       status=nf90_put_att(ncid, varid(ivar),        &
            'units', binData(ivar)%units)
       status=nf90_put_att(ncid, varid(ivar),        &
            'dim', binData(ivar)%dim)
       status=nf90_put_att(ncid, varid(ivar),        &
            'missing_value', unset_real)
    end do

    status=nf90_enddef(ncid)

    do ivar=1,n_dgs
       status=nf90_put_var(ncid, varid(ivar), &
            binData(ivar)%data(1:max_nbins))
       call check_ncstatus(status)
    end do


    deallocate(varid)

    ! Do scalars (including time)
    n_dgs=ID_Scalars%nids
    allocate(varid(n_dgs))
    status=nf90_redef(ncid)
    dim=timeid

    do ivar=1,n_dgs
       name=scalars(ivar)%name
       call sanitize(name)
       status=nf90_def_var(ncid, name &
            , NF90_DOUBLE, dim, varid(ivar))
       call check_ncstatus(status, scalars(ivar)%name)
       status=nf90_put_att(ncid, varid(ivar),        &
            'units', scalars(ivar)%units)
       status=nf90_put_att(ncid, varid(ivar),        &
            'dim', scalars(ivar)%dim)
       status=nf90_put_att(ncid, varid(ivar),        &
            'missing_value', unset_real)
    end do

    status=nf90_enddef(ncid)

    do ivar=1,n_dgs
       status=nf90_put_var(ncid, varid(ivar), &
            scalars(ivar)%data(1:n_dgtimes))
       call check_ncstatus(status)
    end do


    deallocate(varid)

    if (nx == 1) then

       ! Do the instantaneous diags (1D)
       n_dgs=ID_instant_column%nids
       allocate(varid(n_dgs))
       status=nf90_redef(ncid)
       dim2d=tzdim

       do ivar=1,n_dgs
       ! print *, 'ivar', ivar

          status=nf90_def_var(ncid, instant_column(ivar)%name &
               , NF90_DOUBLE, dim2d, varid(ivar))
            ! print*, ncid, instant_column(ivar)%name
            ! print*, NF90_DOUBLE, dim2d, varid(ivar)
          call check_ncstatus(status)
          status=nf90_put_att(ncid, varid(ivar),        &
               'units', instant_column(ivar)%units)
          status=nf90_put_att(ncid, varid(ivar),        &
               'dim', instant_column(ivar)%dim)
          status=nf90_put_att(ncid, varid(ivar),        &
               'missing_value', unset_real)
       end do

       status=nf90_enddef(ncid)

       do ivar=1,n_dgs
          status=nf90_put_var(ncid, varid(ivar), &
               transpose(instant_column(ivar)%data(1:nz,1:n_dgtimes)))
          call check_ncstatus(status)
       end do

       deallocate(varid)

       ! Do the instantaneous bin diags (2D)
       n_dgs=ID_instant_bindgs%nids
       allocate(varid(n_dgs))
       status=nf90_redef(ncid)
       dim3d=tzldim

       do ivar=1,n_dgs
          status=nf90_def_var(ncid, instant_bindgs(ivar)%name &
               , NF90_DOUBLE, dim3d, varid(ivar))
          call check_ncstatus(status)
          status=nf90_put_att(ncid, varid(ivar),        &
               'units', instant_bindgs(ivar)%units)
          status=nf90_put_att(ncid, varid(ivar),        &
               'dim', instant_bindgs(ivar)%dim)
          status=nf90_put_att(ncid, varid(ivar),        &
               'missing_value', unset_real)
       end do

       status=nf90_enddef(ncid)

       allocate(field3d(n_dgtimes, max_nbins, nz))

       do ivar=1,n_dgs
          do tim=1,n_dgtimes
             do iq=1, max_nbins
                do k = 1,nz
                   field3d(tim,iq,k)= instant_bindgs(ivar)%data(k,iq,tim)
                enddo
             enddo
          enddo
          status=nf90_put_var(ncid, varid(ivar), field3d)

          call check_ncstatus(status)

       enddo

       deallocate(field3d)

       deallocate(varid)


       ! Do the instantaneous proc diags (2D)
       n_dgs=ID_instant_procdgs%nids
       allocate(varid(n_dgs))
       status=nf90_redef(ncid)
       dim3d=tzmdim

       do ivar=1,n_dgs
          status=nf90_def_var(ncid, instant_procdgs(ivar)%name &
               , NF90_DOUBLE, dim3d, varid(ivar))
          call check_ncstatus(status)
          status=nf90_put_att(ncid, varid(ivar),        &
               'units', instant_procdgs(ivar)%units)
          status=nf90_put_att(ncid, varid(ivar),        &
               'dim', instant_procdgs(ivar)%dim)
          status=nf90_put_att(ncid, varid(ivar),        &
               'missing_value', unset_real)
       end do

       status=nf90_enddef(ncid)

       allocate(field3d(n_dgtimes, nmom_diag, nz))

       do ivar=1,n_dgs
          do tim=1,n_dgtimes
             do imom=1, nmom_diag
                do k = 1,nz
                   field3d(tim,imom,k)= instant_procdgs(ivar)%data(k,imom,tim)
                enddo
             enddo
          enddo
          status=nf90_put_var(ncid, varid(ivar), field3d)

          call check_ncstatus(status)

       enddo

       deallocate(field3d)

       deallocate(varid)


       ! Do the instantaneous flag diags (2D)
       n_dgs=ID_instant_flagdgs%nids
       allocate(varid(n_dgs))
       status=nf90_redef(ncid)
       dim3d=tzfdim

       do ivar=1,n_dgs
          status=nf90_def_var(ncid, instant_flagdgs(ivar)%name &
               , NF90_DOUBLE, dim3d, varid(ivar))
          call check_ncstatus(status)
          status=nf90_put_att(ncid, varid(ivar),        &
               'units', instant_flagdgs(ivar)%units)
          status=nf90_put_att(ncid, varid(ivar),        &
               'dim', instant_flagdgs(ivar)%dim)
          status=nf90_put_att(ncid, varid(ivar),        &
               'missing_value', unset_real)
       end do

       status=nf90_enddef(ncid)

       allocate(field3d(n_dgtimes, flag_count, nz))

       do ivar=1,n_dgs
          do tim=1,n_dgtimes
             do iq=1, flag_count
                do k = 1,nz
                   field3d(tim,iq,k)= instant_flagdgs(ivar)%data(k,iq,tim)
                enddo
             enddo
          enddo
          status=nf90_put_var(ncid, varid(ivar), field3d)

          call check_ncstatus(status)

       enddo

       deallocate(field3d)

       deallocate(varid)

    else

       ! Do the instantaneous diags (2D)
       n_dgs=ID_instant_2D%nids
       allocate(varid(n_dgs))
       status=nf90_redef(ncid)
       dim3d=tzxdim

       do ivar=1,n_dgs
          status=nf90_def_var(ncid, instant_2D(ivar)%name &
               , NF90_DOUBLE, dim3d, varid(ivar))
          call check_ncstatus(status)
          status=nf90_put_att(ncid, varid(ivar),        &
               'units', instant_2D(ivar)%units)
          status=nf90_put_att(ncid, varid(ivar),        &
               'dim', instant_2D(ivar)%dim)
          status=nf90_put_att(ncid, varid(ivar),        &
               'missing_value', unset_real)
       end do

       status=nf90_enddef(ncid)

       allocate(field3d(n_dgtimes, nx, nz))

       do ivar=1,n_dgs
          do tim=1,n_dgtimes
             do j = 1, nx
                do k = 1,nz
                   field3d(tim,j,k)= field_mask(k,j)*instant_2D(ivar)%data(k,j,tim) &
                        +(1.-field_mask(k,j))*unset_real
                enddo
             enddo
          enddo
          status=nf90_put_var(ncid, varid(ivar), field3d)

          call check_ncstatus(status)

       enddo

       deallocate(field3d)

       deallocate(varid)

       ! Do the instantaneous diags (nx vs time)
       n_dgs=ID_instant_column%nids
       allocate(varid(n_dgs))
       status=nf90_redef(ncid)
       dim2d=txdim

       do ivar=1,n_dgs
          status=nf90_def_var(ncid, instant_column(ivar)%name &
               , NF90_DOUBLE, dim2d, varid(ivar))
          call check_ncstatus(status)
          status=nf90_put_att(ncid, varid(ivar),        &
               'units', instant_column(ivar)%units)
          status=nf90_put_att(ncid, varid(ivar),        &
               'dim', instant_column(ivar)%dim)
          status=nf90_put_att(ncid, varid(ivar),        &
               'missing_value', unset_real)
       end do

       status=nf90_enddef(ncid)

       do ivar=1,n_dgs
          status=nf90_put_var(ncid, varid(ivar), &
               transpose(instant_column(ivar)%data(1:nx,1:n_dgtimes)))
          call check_ncstatus(status)
       end do

       deallocate(varid)

    endif

    status=nf90_close(ncid)
    call check_ncstatus(status)


    ! We now write out the full namelist data to file
    nmlfile=outfile(1:len_trim(outfile)-3)//'.nml'
    open(99, file=nmlfile)
    write(99, '(A)') 'Namelist data corresponding to data file: '//&
         trim(outfile)
    write(99, '(A)') dateString
    write(99, mphys)
    write(99, control)
    write(99, case)
    write(99, switch)
    write(99, addcontrol)
    close(99)


  end subroutine write_diagnostics

  subroutine check_ncstatus(status, name)
    ! Wrapper for checking netcdf command status
    ! Could do more with it if wanted
    Use netcdf
    integer, intent(in) :: status
    character(*), intent(in), optional :: name

    if (status /= nf90_noerr)then
       write(6,*) 'Netcdf error: ', nf90_strerror(status)
       if (present(name)) &
            write(6,*) 'Error occured with variable: ', trim(name)
       stop
    end if

  end subroutine check_ncstatus

  subroutine print_dg_1D_point(name, itime, k, dg_obj, dg_index)

    character(*), intent(in) :: name
    integer, intent(in) :: itime, k
    type(diag1DTS), intent(in) :: dg_obj(:)
    type(dgIDarray), intent(inout) :: dg_index

    !local
    integer :: ivar

    call getUniqueId(name, dg_index, ivar)

    if (associated(dg_obj(ivar)%data))then
       print*, name, time, k, dg_obj(ivar)%data(k,itime)
    end if

  end subroutine print_dg_1D_point

  subroutine sanitize(name)
    ! function to replace character that netcdf doesn't
    ! like with alternatives.
    !
    character(*), intent(inout) :: name
    integer :: i

    do i = 1, len_trim(name)
       if( name(i:i) == '[' ) name(i:i) = '_'
       if( name(i:i) == ']' ) name(i:i) = '_'
       if( name(i:i) == '=' ) name(i:i) = '-'
       if( name(i:i) == ' ' ) name(i:i) = '_'
    end do

  end subroutine sanitize


end Module diagnostics
