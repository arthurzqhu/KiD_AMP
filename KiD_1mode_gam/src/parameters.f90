! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Some parameters and arrays needed by different schemes
! Note that some of these are no longer parameters, but can be
! passed in through a namelist (see namelists.f90)
!
!

#ifndef DEF_NX
#define DEF_NX 1
#endif

#ifndef DEF_NZ
#define DEF_NZ 120
#endif

module parameters


  Use typeKind
  Implicit None

  ! control parameters
  real               :: dt=1.        ! Time step
  integer, parameter :: nz=DEF_NZ    ! Number of height levels
  integer, parameter :: nx=DEF_NX    ! Number of horizontal gridpoints

  ! microphysics specific parameters
  integer, parameter :: max_nmoments=6
  integer, parameter :: max_nbins=34
  integer, parameter :: flag_count=4
  integer, parameter :: nspecies=2
#if SHIPWAY_MICRO == 1
  integer, parameter :: naerosol=9 ! number of aerosol species
#else
  integer, parameter :: naerosol=1 ! number of aerosol species
#endif

  ! ppe and latin hypercube
  double precision, allocatable :: lsample(:,:)

  ! diagnostic parameters
  integer, parameter ::     &
       max_dgs       = 300  & ! Maximum number of dg variables in
                              ! a diagnostic array
      ,max_char_len  = 200    ! Maximum variable name length

  real(wp)           :: dg_dt  = 60.
                        ! time interval between diagnostic calcs
                        ! should be an integer multiple of dt

  integer            :: diaglevel = 5 ! switch for level of diagnostic output

  real(wp)           :: dgstart = 0.0 ! start time for diagnostics

 ! other parameters
  real(wp),parameter :: unset_real=-999
  integer,parameter  :: unset_integer=-999

  ! Microphysics parameters
  integer :: num_aero_moments(naerosol) = 1
  integer :: num_aero_bins(naerosol)= 1
  real(wp) :: aero_mom_init(max_nmoments)
  real(wp) :: aero_N_init(naerosol)=0.0
  real(wp) :: aero_rd_init(naerosol)=0.0
  real(wp) :: aero_sig_init(naerosol)=0.0
  real(wp) :: Dm_init=0.

  integer :: num_h_moments(nspecies)= (/  &
        1 & ! cloud
       ,1 & ! rain
       /)
      ! ,2 & ! ice
      ! ,1 & ! snow
      ! ,1 & ! graupel
      ! /)
   real :: h_shape(nspecies) = (/1., 1./)

  ! If using bulk microphysics, simply set these to 1.
  integer :: num_h_bins(nspecies)= 1

  real(wp) :: mom_init(max_nmoments)
  character(10) :: h_names(nspecies)= &
       (/ 'cloud     '   &
       ,  'rain      '   &
       /)
       !,  'ice       '   &
       !,  'snow      '   &
       !,  'graupel   '  /)

  character(10) :: mom_names(max_nmoments)
  character(10) :: mom_units(max_nmoments)= &
       (/ 'kg/kg     '   &
       ,  'm^k/kg    '   &
       ,  'm^k/kg    '   &
       ,  'm^k/kg    '   &
       ,  'm^k/kg    '   &
       ,  'm^k/kg    '  /)

  character(10) :: aero_names(naerosol)='aero'! &
     !(/ 'aitken    '   &
     !,  'accum     '   &
     !,  'coarse    '   &
!#if SHIPWAY_MICRO == 1
!     ,  'active    '   &
!     ,  'active_r  '   &
!     ,  'dust      '   &
!     ,  'dust_ice  '   &
!     ,  'dust_cloud'   &
!     ,  'sol_ice   '   &
!#endif
!     /)

  character(10) :: aero_mom_names(max_nmoments)
  character(10) :: aero_mom_units(max_nmoments)
   ! if using bin microphysics set variable below to
   ! determine the last bin of cloud
   ! if not using bin, this value is ignored
   integer :: split_bins

   type pymc_filedirs_type
     character(len=200) :: nevp_dir
     character(len=200) :: condevp_dir
     character(len=200) :: coal_dir
     character(len=200) :: sed_dir
   end type



end module parameters
