module mphys_boss

use parameters, only : num_h_moments, num_h_bins, h_shape, nspecies, nz, dt &
  , h_names, mom_units, max_char_len, nx, dg_dt, dgstart
use column_variables
use common_physics, only : qsaturation
use physconst, only : p0, r_on_cp, pi, rhow
use namelists, only: log_predictNc, l_truncated, l_init_test, n_cat, l_use_nn &
  , l_boss_partition_liq, l_boss_save_dsd

use module_hujisbm
use micro_prm
use diagnostics, only: save_dg, i_dgtime, my_save_dg_bin_dp
use switches, only: l_advect,l_diverge, l_noadv_theta, l_noadv_qv
use switches, only: zctrl
use runtime, only: time
use module_mp_boss
! use global_fun

real, dimension(nz,nx) :: th_new,th_old,qv_new,qv_old,qitot,qirim,nitot,birim,ssat,uzpl,pres,&
  dzq,diag_ze,diag_effc,diag_effi,diag_vmi,diag_di,diag_rhoi,SCF_out
real, dimension(nx) :: prt_liq,prt_sol
real :: clbfact_dep, clbfact_sub,scpf_pfrac,scpf_resfact
character(len=16), parameter :: model = 'KiD'
logical :: typeDiags_ON,debug_on,scpf_on,micro_unset=.true.,l_diag
real, dimension(nz,max_nbins) :: fieldbin
real(8), dimension(nz) :: fielddp
real(8), save, dimension(nz,nx,2) :: guessc,guessr


contains
  
subroutine mphys_boss_interface(npmcr)
use mphys_tau_bin_declare, only: xkgmean, dgmean
use module_bin_init, only: bin_init, data

integer :: j,k,imom,it
integer :: npmcr(nspecies), p3stat
real, dimension(nz,nx,npmcr(1)) :: qcs
real, dimension(nz,nx,npmcr(2)) :: qrs
real, dimension(nz) :: Dm_c, Dm_r, Dm_w
real, dimension(nx) :: prt_drzl,prt_rain,prt_crys,prt_snow,prt_grpl,prt_pell,prt_hail,prt_sndp
real, dimension(nz,nx,max_nbins) :: ffcd_mass_befm, ffcd_num_befm, ffcd_mass_aftm, ffcd_num_aftm
real(8),dimension(nz,nx,1:10) :: mc,mr
character(max_char_len) :: name, units
character(1) :: Mnum
logical :: prep_stop = .false.

npm = maxval(npmcr)

scpf_pfrac = 0.
scpf_resfact = 0.
scpf_on = .false.
clbfact_dep = 1.
clbfact_sub = 1.
debug_on = .false.
l_diag = (mod(time,dg_dt)<dt) .and. (time>=dgstart)

it = time/dt
qcs(:,:,:) = 0.
qrs(:,:,:) = 0.

do imom = 1,npm
  do j=1,nx
    do k=1,nz
        qcs(k,j,imom) = hydrometeors(k,j,1)%moments(1,imom)
        if (n_cat==2 .or. (l_boss_partition_liq .and. mod(time-dt,dg_dt)<dt.and.time>dgstart+dt)) then
          qrs(k,j,imom) = hydrometeors(k,j,2)%moments(1,imom)
          ! if (k==80 .and. imom==1) print*, time, l_diag, qrs(k,j,imom),qcs(k,j,imom)
          ! if (k==80 .and. imom==1) stop
        ! else
        endif
        ! if (k==80 .and. imom==1) print*, sngl(time), l_diag, qrs(k,j,imom),qcs(k,j,imom)
    enddo
  enddo
enddo
    ! print*, 'before', time, hydrometeors(1,1,1)%moments(1,1), hydrometeors(1,1,2)%moments(1,1)

if (n_cat==1 .and. l_boss_partition_liq .and. mod(time-dt,dg_dt)<dt.and.time>dgstart+dt) then
  qcs = qcs+qrs
  qrs = 0.
endif

! qcs(:,:,1) now mass
qcs(:,:,1) = qcs(:,:,1)*M3toq
if (n_cat==2) qrs(:,:,1) = qrs(:,:,1)*M3toq
th_old(:,1) = theta(:,1)
qv_old(:,1) = qv(:,1)
th_new(:,1) = theta(:,1)
qv_new(:,1) = qv(:,1)
qitot = 0.
qirim = 0.
nitot = 0.
birim = 0.
! note: code for prediction of ssat not currently avaiable, set 2D array to 0
ssat = 0.
uzpl(:,1) = w(:,1)
pres(:,1) = p0*exner(:,1)**(1./r_on_cp)
dzq(:,1) = dz
typeDiags_ON = .false.

if (micro_unset .and. l_boss_partition_liq) then

  if (n_cat==2) then
    pmomsc(1:3)=(/3,imomc1,imomc2/)
    pmomsr(1:3)=(/3,imomr1,imomr2/)
  else
    pmomsc(1:6)=(/3, 0, imomc1, imomc2, imomr1, imomr2/)
    pmomsr(1:6)=(/3, 0, imomc1, imomc2, imomr1, imomr2/)
  endif

  guessc(:,:,1) = h_shape(1) !shape parameter
  guessr(:,:,1) = h_shape(2)
  guessc(:,:,2) = dnc_def        !characteristic diameter dn
  guessr(:,:,2) = dnr_def

  call set_micro
  call bin_init !initialises the cloud bin categories
  call data     !reads in and sets the coll-coal kernal

  diams=dgmean
  do imom = 1,10
    do ib = 1,nkr
      pdiams(imom, ib) = diams(ib)**(dble(imom)-1.)
    enddo
  enddo


endif


if (n_cat==1) then ! single cat
  if (micro_unset) then
    call boss_slc_init(lookup_file_dir='.', nCat=1, trplMomI=.false., model='KiD', &
      stat=p3stat, abort_on_err=.false., dowr=.true.,qcs=qcs,its=1,ite=nx,kts=1,kte=nz,&
      npm=npm)
    micro_unset = .false.
    ! l_use_nn = .true.
    ! if (l_use_nn) &
    !   call moment2state_net % load("/Users/arthurhu/Downloads/moments_to_state_fixed_mu_nn.txt")
  endif

  ! print*, 'bef', qcs(80,1,1)
  ! if (mod(time-dt,dg_dt)<dt.and.time>dgstart+dt) print*, 'hyd', hydrometeors(80,1,1)%moments(1,1), hydrometeors(80,1,2)%moments(1,1)
  call boss_slc_main(qcs,npm,th_old,th_new,qv_old,qv_new,dt,qitot,qirim,&
    nitot,birim,ssat,uzpl,pres,dzq,it,prt_liq,prt_sol,1,nx,1, &
    nz,1,diag_ze,diag_effc,diag_effi,diag_vmi,diag_di,diag_rhoi,&
    log_predictNc,typeDiags_ON,trim(model),clbfact_dep,clbfact_sub,debug_on,scpf_on,&
    scpf_pfrac,scpf_resfact,SCF_out)
  ! print*, 'aft', qcs(80,1,1)
  ! print*, qcs(:,:,2)

else ! two cat
  if (micro_unset) then
    call boss_2cat_init(lookup_file_dir='.', nCat=1, trplMomI=.false., &
      model='KiD',stat=p3stat,abort_on_err=.false., dowr=.true.,iparam=4, &
      qcs=qcs,qrs=qrs,its=1,ite=nx,kts=1,kte=nz,npm=npm)
    micro_unset = .false.
  endif

  ! print*, 'bef', qrs(80,1,1:2)
  call boss_2cat_main(qcs(:,:,1),qcs(:,:,2),qrs(:,:,1),qrs(:,:,2),th_old,th_new, &
    qv_old,qv_new,dt,qitot,qirim,nitot,birim,ssat,uzpl,pres,dzq,it,prt_liq,prt_sol, &
    1,nx,1,nz,1,diag_ze,diag_effc,diag_effi,diag_vmi,diag_di,diag_rhoi,4,log_predictNc, &
    typeDiags_ON,trim(model),clbfact_dep,clbfact_sub,debug_on,scpf_on,scpf_pfrac,&
    scpf_resfact,SCF_out)
  ! print*, 'aft', qrs(80,1,1:2)

endif

! qcs(:,:,1) now 3rd moment
qcs(:,:,1) = qcs(:,:,1)*QtoM3
if (n_cat==2) qrs(:,:,1) = qrs(:,:,1)*QtoM3

! back out tendencies
dqv_mphys(:,1) = (qv_new(:,1) - qv_old(:,1))/dt
dtheta_mphys(:,1) = (th_new(:,1) - th_old(:,1))/dt

do imom = 1, num_h_moments(1)
  do j=1,nx
    do k=1,nz
      if (qcs(k,j,imom) .ne. qcs(k,j,imom)) qcs(k,j,imom) = 0.
      dhydrometeors_mphys(k,j,1)%moments(1,imom) = &
        (qcs(k,j,imom) - hydrometeors(k,j,1)%moments(1,imom))/dt
      ! if (dhydrometeors_mphys(k,j,1)%moments(1,imom) .ne. dhydrometeors_mphys(k,j,1)%moments(1,imom)) then
      !   print*, 'nan found at k, imom', k, imom
      !   print*, 'qcs', qcs(k-5:k+5,j,imom)
      !   stop
      ! endif
    enddo
  enddo
end do
! print*, qcs(15,1,:)

do imom = 1, num_h_moments(2)
  do j=1,nx
    do k=1,nz
      if (n_cat==2) then
        dhydrometeors_mphys(k,j,2)%moments(1,imom) = &
          (qrs(k,j,imom) - hydrometeors(k,j,2)%moments(1,imom))/dt
      else
        dhydrometeors_mphys(k,j,2)%moments(1,imom) = 0.
      endif
    enddo
  enddo
enddo

! save surface rain rate
if (l_diag) then
  call save_dg(sum(prt_liq(1:nx))/nx,'mean_surface_ppt',i_dgtime,'kg m-2 s-1',dim='time')
endif

if (l_boss_partition_liq .and. l_diag) then

  call invert_mom_boss(qcs, qrs, ffcd_mass_aftm, ffcd_num_aftm, npmcr)

  do j=1,nx
    do k=1,nz
      call calcmoms(dble(ffcd_mass_aftm(k,j,:)), dble(ffcd_num_aftm(k,j,:)), 10, mc(k,j,:), mr(k,j,:))
    enddo
  enddo
  ! print*, 'qcs', time, qcs(80,1,1)
  ! print*, 'mcs', ime, sngl(mc(80,1,4)),sngl(mr(80,1,4))

  Dm_c = (mc(:,1,4)/mc(:,1,1))**(1/3.)
  Dm_r = (mr(:,1,4)/mr(:,1,1))**(1/3.)
  Dm_w = ((mc(:,1,4)+mr(:,1,4))/(mc(:,1,1)+mr(:,1,1)))**(1/3.)
  call save_dg(Dm_c, 'Dm_c',i_dgtime,'m',dim='z')
  call save_dg(Dm_r, 'Dm_r',i_dgtime,'m',dim='z')
  call save_dg(Dm_w, 'Dm_w',i_dgtime,'m',dim='z')

  do imom=1,10
    write(Mnum,'(I1)') imom-1
    name='diagM'//Mnum//'_cloud'
    units='m^'//Mnum
    fielddp(:)=mc(:,1,imom)
    call save_dg(fielddp,name,i_dgtime,units,dim='z')

    name='diagM'//Mnum//'_rain'
    fielddp(:)=mr(:,1,imom)
    call save_dg(fielddp,name,i_dgtime,units,dim='z')
  enddo

  if (l_boss_save_dsd) then
    fieldbin(:,:)=ffcd_mass_aftm(:,nx,:)
    name='mass_dist_init'
    units='kg/kg/ln(r)'
    call save_dg('bin',fieldbin,name,i_dgtime,units)

    fieldbin(:,:)=ffcd_num_aftm(:,nx,:)
    name='num_dist_init'
    units='1/kg/ln(r)'
    call save_dg('bin',fieldbin,name,i_dgtime,units)
  endif

  ! if it's diagnosed then update moments with mc, mr instead

  if (n_cat==1) then
    do imom = 1, num_h_moments(1)
      ip = pmomsc(imom)+1
      do j=1,nx
        do k=1,nz
          dhydrometeors_mphys(k,j,1)%moments(1,imom) = &
            (mc(k,j,ip) - hydrometeors(k,j,1)%moments(1,imom))/dt
        enddo
      enddo
    end do
    ! print*, 'dmc', mc(80,1,4), hydrometeors(80,1,1)%moments(1,1)

    do imom = 1, num_h_moments(2)
      ip = pmomsr(imom)+1
      do j=1,nx
        do k=1,nz
          ! to ensure that the moment set is not broken by the DSD diagnosis
          qrs(k,j,imom) = qcs(k,j,imom) - mc(k,j,ip)
          dhydrometeors_mphys(k,j,2)%moments(1,imom) = &
            (qrs(k,j,imom) - hydrometeors(k,j,2)%moments(1,imom))/dt
        enddo
      enddo
    enddo
  endif

endif

! also diagnose cloud vs rain

end subroutine mphys_boss_interface

subroutine invert_mom_boss(qcs, qrs, ffcd_mass, ffcd_num, npmcr)

use module_mp_amp, only: invert_moments, invert_moments_nn

integer, intent(in) :: npmcr(nspecies)
real, dimension(nz,nx,npmcr(1)), intent(in) :: qcs
real, dimension(nz,nx,npmcr(2)), intent(in) :: qrs
real, dimension(nz,nx,max_nbins), intent(out) :: ffcd_mass, ffcd_num
integer :: j,k,imom,it
double precision :: mom_pred(num_h_moments(1), n_cat)
double precision :: gam_param(2,2), sqerr, flag(nz,nx,n_cat)

! NOTE: make sure qcs(:,:,1) is M3 not q

ffcd_mass(:,:,:) = 0.
ffcd_num(:,:,:) = 0.

do j=1,nx
  do k=1,nz
    mom_pred(1:npm,1) = qcs(k,j,1:npm)
    if (n_cat==2) mom_pred(1:npm,2) = qrs(k,j,1:npm)

    ! mom_pred(1,1) = mom_pred(1,1)*QtoM3
    ! if (n_cat==2) mom_pred(1,2) = mom_pred(1,2)*QtoM3

    gam_param(1,1) = guessc(k,j,2)
    gam_param(2,1) = guessc(k,j,1)
    gam_param(1,2) = guessr(k,j,2)
    gam_param(2,2) = guessr(k,j,1)
    if (sum(mom_pred(1,:)) <= total_m3_th) then
      ffcd_mass(k,j,:) = 0.
      ffcd_num(k,j,:) = 0.
    else
      call invert_moments(mom_pred, gam_param, ffcd_mass(k,j,:), &
        ffcd_num(k,j,:),flag(k,j,:),sqerr)
    endif
    guessc(k,j,2) = gam_param(1,1)
    guessc(k,j,1) = gam_param(2,1)
    guessr(k,j,2) = gam_param(1,2)
    guessr(k,j,1) = gam_param(2,2)
  enddo
enddo

end subroutine invert_mom_boss

end module mphys_boss
