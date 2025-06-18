! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module to deal with namelists for input
!
module namelists

  Use parameters
  Use header_data, only : mphys_id
  Use switches
  Use switches_bin
  ! Use micro_prm, only:imomc1,imomc2,imomr1,imomr2,donucleation, &
  !                     docondensation,docollisions,dosedimentation, &
  !                     cloud_init,rain_init

#if SHIPWAY_MICRO == 1
  ! Temporary for adding in 4a switches
  Use mphys_switches, only: iopt_act, option, aerosol_option        &
     , l_aaut, l_aacc, l_aevp, l_ased, iopt_tidy, l_warm            &
     , l_inuc, iopt_rcrit, iopt_inuc                                &
     , l_evaporation, l_rain, l_sed, l_boussinesq, diag_mu_option   &
     , l_sed_3mdiff, l_cons, l_abelshipway, l_condensation          &
     , l_active_inarg2000, l_oneway, l_newoptions
  Use mphys_parameters, only: p1, p2, p3, sp1, sp2, sp3             &
     , max_step_length, max_sed_length
#endif

  implicit none
  integer:: imomc1,imomc2,imomr1,imomr2,n_cat,nmom_diag
  real(8) :: ss_init
  real(8), dimension(4):: cloud_init,rain_init, rain_source ! rain_source=Dm,N
  real(8), dimension(2) :: cloud_layer_DN
  real(8) :: ovc_factor=0.0 ! overcorrection factor
  real(8) :: rhctrl ! actually saturation ratio not relative humidity
  real(8), allocatable :: moments_diag(:)
  logical :: docollisions, docondensation, donucleation, dosedimentation, &
             dobreakup, l_truncated, l_init_test, log_predictNc, l_use_nn, &
             l_boss_partition_liq, l_ppe, l_getrates, l_boss_save_dsd, &
             l_ppe_nevp, l_ppe_condevp, l_ppe_coal, l_ppe_sed, doevap
  character(200) :: param_val_fpath = "../../CloudBOSS/boss_slc_param_values.csv"
  character(200) :: param_infl_sigma_fpath = "../../CloudBOSS/boss_slc_param_sigma.csv"
  character(200) :: param_val_fpath_2cat = "../../CloudBOSS/boss_2cat_param_values.csv"
  character(20)  :: s_sample_dist
  character(200) :: posterior_path, lhs_path
  integer :: n_init

! integer switch for type of BOSS autoconversion (q)
!1 = 1-term auto wo rain dependence
!2 = 2-term auto wo rain dependence
!3 = 2-term auto w rain dependence
integer, public :: iautoq=1

! integer switch for type of BOSS Nc fall speed 
! 1 = default P3 method 
! 2 = Stokes velocity of mean mass radius cloud (mc)
integer, public :: ivTnc = 1

! integer switch for type of BOSS qc fall speed
! 1 = default P3 method 
! 2 = power law based on mc
integer, public :: ivTqc = 1

! integer switch for type of BOSS Nr fall speed 
! 1 = default P3 method 
! 2 = power law based on mr 
integer, public :: ivTnr = 1

! integer switch for type of BOSS qr fall speed
! 1 = default P3 method 
! 2 = power law based on mr
integer, public :: ivTqr = 1

!!!!!!!!!!!!!!!!!!!!!!!
! P3 limiters in nml 
! defaults set to current P3 values 
! kl add 101923

! minimum number mean cloud drop diameter [m]
real,public :: dNc_min = 1e-6
! maximum number mean cloud drop diameter [m]
real,public :: dNc_max = 40e-6
! minimum number mean raindrop diameter [m]
real,public :: dNr_min = 10e-6
! maximum number mean raindrop diameter [m]
real,public :: dNr_max = 2e-3

! cloud number fall speed parameters 
! for ivTnc = 1
! no parameters used 
! for ivTnc = 2
! vTnc = min(vT_stokes(mc),vTncmax)
real, public :: vTncmax = 0.5 ! [m/s]

! cloud mass fall speed parameters 
! for ivTqc = 1
! no parameters used 
! for ivTqc = 2
! vTqc = min(10^log_a_vTqc * mr^b_vTqc,vTqcmax)
real, public :: vTqcmax = 0.5 ! [m/s]

! rain number fall speed parameters 
! for ivTnr = 1
! no parameters used 
! for ivTnr = 2
! vTnr = min(10^log_a_vTnr * mr^b_vTnr,vTnrmax)
real, public :: vTnrmax = 10. ! [m/s]

! rain mass fall speed parameters 
! for ivTqr = 1
! no parameters used 
! for ivTqr = 2
! vTqrmax = min(10^log_a_vTqr * mr^b_vTqr,vTqrmax)
real, public :: vTqrmax = 10. ! [m/s]
integer, public :: idraw = 1, rand_seed

! ppe variables
integer, public :: irealz, n_ppe
real, public :: deflation_factor = 1., Na_min, Na_max, w_min, w_max, Dm_min, Dm_max, &
  Nd_min, Nd_max
logical :: l_save_tend, l_save_adv, l_save_div, l_save_mphys
character(200) :: kidpath
character(200) :: nevp_dir, condevp_dir, coal_dir, sed_dir

  namelist/mphys/num_h_moments, num_h_bins, h_shape, mom_init, &
       h_names, mom_names, mom_units,num_aero_moments,num_aero_bins, &
       aero_mom_init, aero_N_init, aero_sig_init, aero_rd_init, aero_names, &
       imomc1,imomc2,imomr1,imomr2,doevap,donucleation,docondensation,docollisions, &
       dosedimentation,dobreakup,cloud_init,rain_init,ss_init, rain_source, &
       param_val_fpath, param_infl_sigma_fpath, log_predictNc, param_val_fpath_2cat,iautoq,ivTnc,ivTqc, &
       ivTnr,ivTqr,dNc_min,dNc_max,dNr_min,dNr_max,vTncmax,vTqcmax,vTnrmax,vTqrmax,&
       n_cat,idraw,rand_seed,nmom_diag,cloud_layer_DN

  namelist/control/dt, dg_dt, mphys_scheme, mphys_var &
       , wctrl, zctrl, tctrl, pctrl_z, pctrl_v, pctrl_T, ipctrl &
       , xctrl, lhf_ctrl, shf_ctrl, diaglevel, dgstart, rhctrl
  namelist/case/input_file, l_input_file, ifiletype, icase

  namelist/ppe/l_ppe, s_sample_dist, posterior_path, irealz, &
               deflation_factor, Na_min, Na_max, w_min, w_max, &
               n_ppe, lhs_path, l_ppe_nevp, l_ppe_condevp, l_ppe_coal, l_ppe_sed, &
               Dm_min, Dm_max, n_init, nevp_dir, condevp_dir, coal_dir, sed_dir, &
               Nd_min, Nd_max

  namelist/switch/l_mphys, l_advect, l_diverge, l_pupdate &
       , l_fix_qv, l_nomphys_qv, l_noadv_qv, l_posadv_qv &
       , l_fix_theta, l_nomphys_theta, l_noadv_theta  &
       , l_noadv_hydrometeors, l_nodiv_hydrometeors, l_sediment &
       , isurface, l_noadv_aerosols, l_nodiv_aerosols, l_fix_aerosols &
       , l_sed_ult, l_diverge_advection, l_periodic_bound  &
       , l_force_positive, l_noevaporation, l_nocondensation &
       , l_truncated, l_init_test, l_use_nn, l_boss_partition_liq &
       , l_getrates, l_boss_save_dsd, l_save_tend, l_save_adv, l_save_div, l_save_mphys

  logical :: iiwarm=.false.
  character(200) :: KiD_outdir=''
  character(200) :: KiD_outfile=''
  character(3) :: bintype='' ! underlying 'sbm' or 'tau'
  character(3) :: ampORbin='' ! run 'bin' as a standalone or with 'amp' on top of it
  character(1) :: initprof='' ! 'c' as a constant column of water, 'i' linearly increasing ...
  ! towards the cloud top, where cloud and rain mass = cloud_init(1) and rain_init(1) 
  logical :: mp_proc_dg, extralayer, l_hist_run

  namelist/addcontrol/iiwarm, KiD_outdir, KiD_outfile, ovc_factor, &
          mp_proc_dg, bintype, ampORbin, l_coll_coal, initprof, extralayer, &
          l_hist_run, moments_diag, kidpath &
#if SHIPWAY_MICRO == 1
     ! Shipway 4A ...
     , option, l_evap, l_sed_3mdiff &
     , max_step_length, max_sed_length, diag_mu, l_rsource, l_raut &
     , l_evaporation, l_rain, l_sed, l_boussinesq, diag_mu_option   &
     , p1, p2, p3 &
     , sp1, sp2, sp3 &
     , l_abelshipway, l_cons &
     , l_cond, l_condensation, iopt_act &
     , aerosol_option, l_aaut, l_aacc, l_aevp, l_ased &
     , iopt_tidy, l_warm, l_inuc, iopt_rcrit   &
     , l_active_inarg2000, iopt_inuc, l_cu_cold, l_oneway, l_newoptions &
#endif
     ! Thompson 09...
     , l_reuse_thompson_lookup

  ! Namelist input...

  character(200) :: fileName=''
  character(200) :: fileNameIn=''
  character(200) :: fileNameOut=''
  character(200) :: namelistIn='namelists/input.nml'
  character(200) :: fileIn=''
  character(200) :: fileOut=''
  logical :: fexist

  namelist/namelistToUse/fileIn, fileOut

contains

  subroutine read_namelist
    !
    ! Read the namelists from file
    !

#if COMMANDLINE == 1
    ! This bit is F2003 - If your compiler doesnt support
    ! it you need to comment the line out you can then specify
    ! which namelist to use throught namelists/input.nml file
    ! or else use command line processing that works with your
    ! compiler (nearly all will do it but not in a portable way).
    write(*,*) 'Querying command line'
    CALL GET_COMMAND_ARGUMENT(1,fileNameIn)
    CALL GET_COMMAND_ARGUMENT(2,fileNameOut)
#endif

    if (trim(fileNameIn)=='')then  ! Not input at the command line
                                   ! so use input.nml
#ifdef DEF_CASE
      write(namelistIn, '(A,A,A)') 'namelists/','DEF_CASE','_input.nml'
#endif

      write(*,*) 'Unable to determine input file from command line, so querying ', trim(namelistIn), ' instead...'
      inquire(file=namelistIn, exist=fexist)
      if (fexist) then
        open(2, file=namelistIn)
        read(2, namelistToUse)
        close(2)
        write(*, namelistToUse)
      end if
      fileNameIn  = fileIn
      if (trim(fileOut)/='')fileNameOut = fileOut
    end if

    if (trim(fileNameIn)/='')fileName=fileNameIn

    write(6,*) 'Using namelist: ', trim(fileName)

    open(1, file=fileName)
!    rewind(1)
    read(1,mphys)
!    rewind(1)
    read(1,case)
!    rewind(1)
    read(1,control)
!    rewind(1)
    read(1,ppe)
!    rewind(1)
    read(1,switch)
!    rewind(1)
allocate(moments_diag(nmom_diag))
    read(1,addcontrol)
    close(1)

    select case(mphys_scheme)
    case('lem2.4')
       imphys=imphys_lem2_4
       mphys_id='LEM2.4'
    case('tau_bin')
       imphys=imphys_tau_bin
       mphys_id='TAU_bin'
    case('thompson') ! NB same as thompson09
       imphys=imphys_thompson09
       mphys_id='thompson09'
    case('thompson09')
       imphys=imphys_thompson09
       mphys_id='thompson09'
    case('thompson06')
       imphys=imphys_thompson06
       mphys_id='thompson06'
    case('thompson07')
       imphys=imphys_thompson07
       mphys_id='thompson07'
    case('morr_two_moment')
       imphys=imphys_morr_two_moment
       mphys_id='morr_two_moment'
    case('um7_3')
       imphys=imphys_um7_3
       mphys_id='um7_3'
    case('wsm6')
       imphys=imphys_wsm6
       mphys_id='wsm6'
    case('wdm6')
       imphys=imphys_wdm6
       mphys_id='wdm6'
    case('4A')
       imphys=imphys_4A
       mphys_id='4A'
    case('amp')
       imphys=imphys_amp
       mphys_id='amp'
    case('boss')
       imphys=imphys_boss
       mphys_id='boss'
     case default
       print*, 'Mphys scheme not recognized: ', mphys_scheme
       print*, 'Did you mean:'
       print*, '   lem2.4?'
       print*, '   um7_3?'
       print*, '   tau_bin?'
       print*, '   thompson09?'
       print*, '   thompson07?'
       print*, '   morr_two_moment?'
       print*, '   wsm6?'
       print*, '   4A?'
       print*, '(NB not all available in release version)'
       stop
    end select

    if (trim(input_file)=='')l_input_file=.False.

    if (.not. l_input_file)ifiletype=itest_case

  end subroutine read_namelist

end module namelists
