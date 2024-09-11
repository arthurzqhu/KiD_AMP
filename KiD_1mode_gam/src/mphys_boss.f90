module mphys_boss

use parameters, only : num_h_moments, num_h_bins, h_shape, nspecies, nz, dt &
  , h_names, mom_units, max_char_len, nx
use column_variables
use common_physics, only : qsaturation
use physconst, only : p0, r_on_cp, pi, rhow
use namelists, only: log_predictNc

use module_hujisbm
use micro_prm
use diagnostics, only: save_dg, i_dgtime, my_save_dg_bin_dp
use switches, only: l_advect,l_diverge, l_noadv_theta, l_noadv_qv
use switches, only: zctrl
use runtime, only: time
use module_mp_boss
use global_fun

integer :: j,k,imom,it
real, dimension(nz,nx) :: th_new,th_old,qv_new,qv_old,qitot,qirim,nitot,birim,ssat,uzpl,pres,&
  dzq,diag_ze,diag_effc,diag_effi,diag_vmi,diag_di,diag_rhoi,SCF_out
real, dimension(nx) :: prt_liq,prt_sol
real :: clbfact_dep, clbfact_sub,scpf_pfrac,scpf_resfact
character(len=16), parameter :: model = 'KiD'
logical :: typeDiags_ON,debug_on,scpf_on,micro_unset=.true.


contains
  
subroutine mphys_boss_interface(npmcr)

integer :: npmcr(nspecies), p3stat
real, dimension(nz,nx,npmcr(1)) :: qcs
real, dimension(nz,nx,npmcr(2)) :: qrs
real, dimension(nz,nx) :: Dm_c, Dm_r, Dm_w
real, dimension(nx) :: prt_drzl,prt_rain,prt_crys,prt_snow,prt_grpl,prt_pell,prt_hail,prt_sndp
logical :: prep_stop = .false.

npm = maxval(npmcr)

scpf_pfrac = 0.
scpf_resfact = 0.
scpf_on = .false.
clbfact_dep = 1.
clbfact_sub = 1.
debug_on = .false.

it = time/dt
qcs(:,:,:) = 0.
qrs(:,:,:) = 0.

do imom = 1,npm
  do j=1,nx
    do k=1,nz
      qcs(k,j,imom) = hydrometeors(k,j,1)%moments(1,imom)
    enddo
  enddo
enddo

if (npm<4) then 
  do imom = 1,npm
    do j=1,nx
      do k=1,nz
        qrs(k,j,imom) = hydrometeors(k,j,2)%moments(1,imom)
      enddo
    enddo
  enddo
endif

! qcs(:,:,1) now mass
qcs(:,:,1) = qcs(:,:,1)*pio6rw
if (npm<4) qrs(:,:,1) = qrs(:,:,1)*pio6rw
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

if (npm>=4) then ! single cat
  if (micro_unset) then
    call boss_slc_init(lookup_file_dir='.', nCat=1, trplMomI=.false., model='KiD', &
      stat=p3stat, abort_on_err=.false., dowr=.true.,qcs=qcs,its=1,ite=nx,kts=1,kte=nz,&
      npm=npm)
    micro_unset = .false.
  endif

  call boss_slc_main(qcs,npm,th_old,th_new,qv_old,qv_new,dt,qitot,qirim,&
    nitot,birim,ssat,uzpl,pres,dzq,it,prt_liq,prt_sol,1,nx,1, &
    nz,1,diag_ze,diag_effc,diag_effi,diag_vmi,diag_di,diag_rhoi,&
    log_predictNc,typeDiags_ON,trim(model),clbfact_dep,clbfact_sub,debug_on,scpf_on,&
    scpf_pfrac,scpf_resfact,SCF_out)

else ! two cat
  if (micro_unset) then
    call boss_2cat_init(lookup_file_dir='.', nCat=1, trplMomI=.false., &
      model='KiD',stat=p3stat,abort_on_err=.false., dowr=.true.,iparam=4, &
      qcs=qcs,qrs=qrs,its=1,ite=nx,kts=1,kte=nz,npm=npm)
    micro_unset = .false.
  endif

  call boss_2cat_main(qcs(:,:,1),qcs(:,:,2),qrs(:,:,1),qrs(:,:,2),th_old,th_new, &
    qv_old,qv_new,dt,qitot,qirim,nitot,birim,ssat,uzpl,pres,dzq,it,prt_liq,prt_sol, &
    1,nx,1,nz,1,diag_ze,diag_effc,diag_effi,diag_vmi,diag_di,diag_rhoi,4,log_predictNc, &
    typeDiags_ON,trim(model),clbfact_dep,clbfact_sub,debug_on,scpf_on,scpf_pfrac,&
    scpf_resfact,SCF_out)

endif

! qcs(:,:,1) now 3rd moment
qcs(:,:,1) = qcs(:,:,1)*ipio6rw
if (npm < 4) qrs(:,:,1) = qrs(:,:,1)*ipio6rw

! back out tendencies
dqv_mphys(:,1) = (qv_new(:,1) - qv_old(:,1))/dt
dtheta_mphys(:,1) = (th_new(:,1) - th_old(:,1))/dt

do imom = 1, num_h_moments(1)
  do j=1,nx
    do k=1,nz
      dhydrometeors_mphys(k,j,1)%moments(1,imom) = &
        (qcs(k,j,imom) - hydrometeors(k,j,1)%moments(1,imom))/dt
    enddo
  enddo
end do

! if (dhydrometeors_mphys(20,1,1)%moments(1,1)>0) then
!   print*, 'qc in mphys_boss', qcs(20,1,1)*pio6rw
!   print*, 'dqc in mphys_boss', dhydrometeors_mphys(20,1,1)%moments(1,1)*pio6rw
! endif

if (npm<4) then
  do j=1,nx
    do k=1,nz
      do imom = 1, num_h_moments(1)
        dhydrometeors_mphys(k,j,2)%moments(1,imom) = &
          (qrs(k,j,imom) - hydrometeors(k,j,2)%moments(1,imom))/dt
      enddo
    enddo
  enddo
endif

! save surface rain rate
call save_dg(sum(prt_liq(1:nx))/nx*1e3,'mean_surface_ppt',i_dgtime,'kg m-2 s-1',dim='time')

Dm_c = qcs(:,:,1)/qcs(:,:,2)
Dm_r = qrs(:,:,1)/qrs(:,:,2)
Dm_w = (qcs(:,:,1)+qrs(:,:,1))/(qcs(:,:,2)+qrs(:,:,2))
call save_dg(Dm_c, 'Dm_c',i_dgtime,'m',dim='time')
call save_dg(Dm_r, 'Dm_r',i_dgtime,'m',dim='time')
call save_dg(Dm_w, 'Dm_w',i_dgtime,'m',dim='time')

! invert moments just bc we can

end subroutine mphys_boss_interface

end module mphys_boss
