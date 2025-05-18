module ppe_fun

use netcdf

contains

! ----------------------------------------------------------
subroutine read_netcdf(variable, filename, varname, group_name)

  implicit none

  ! Variable declarations
  integer :: ncid         ! NetCDF file ID
  integer :: varid        ! Variable ID
  integer :: grpid        ! Group ID
  integer :: parentid     ! Parent ID - could be ncid or grpid depending on whether there's group
  integer :: status       ! Return value for error handling
  integer :: dimids(2)    ! Dimension IDs for the variable
  integer :: var_dims(2)  ! Dimensions of the variable
  integer :: ndims, nvars, natts, unlimdimid
  integer :: nx, ny, xtype, i
  integer, allocatable, dimension(:) :: varids
  real(8), dimension(:,:), allocatable :: variable   ! Array to hold the variable data

  character(len=200) :: filename
  character(len=*) :: varname
  character(len=100), optional :: group_name
  logical :: l_group

  if (allocated(variable)) deallocate(variable)
  l_group = present(group_name)

  ! Open the NetCDF file (read-only mode)
  status = nf90_open(filename, nf90_nowrite, ncid)
  if (status /= nf90_noerr) then
    print *, 'Error opening NetCDF file in read_netcdf'
    stop
  end if

  if (l_group) then
    ! Get the group ID
    status = nf90_inq_grp_ncid(ncid, group_name, grpid)
    if (status /= NF90_NOERR) stop "Group "//group_name//" not found"
    parentid = grpid
  else
    parentid = ncid
  endif


  status = nf90_inq_varid(parentid, trim(varname), varid)
  if (status /= nf90_noerr) then
    print *, 'Error getting variable ID: ', varname, nf90_strerror(status)
    status = nf90_close(ncid)
    stop
  end if

  status = NF90_INQUIRE_VARIABLE(parentid, varid, xtype=xtype, ndims=ndims, &
      dimids=dimids, natts=natts)

  ! Get dimension lengths (assuming a 2D variable here)
  status = nf90_inquire_dimension(parentid, dimids(1), len=nx)
  if (status /= nf90_noerr) then
    print *, 'Error getting dimension length for dim 1:', nf90_strerror(status)
    status = nf90_close(ncid)
    stop
  endif

  status = nf90_inquire_dimension(parentid, dimids(2), len=ny)
  if (status /= nf90_noerr) then
    print *, 'Error getting dimension length for dim 2:', nf90_strerror(status)
    status = nf90_close(ncid)
    stop
  endif

  ! Allocate the array to hold the data
  if (.not. allocated(variable)) allocate(variable(nx, ny))

  ! Read the variable data
  status = nf90_get_var(parentid, varid, variable)
  if (status /= nf90_noerr) then
    print *, 'Error reading variable data: ', status
    status = nf90_close(ncid)
    stop
  end if

  ! Close the NetCDF file
  status = nf90_close(ncid)
  if (status /= nf90_noerr) then
    print *, 'Error closing NetCDF file', status
    stop
  end if

  ! Output the data (or process as needed)
  ! print *, 'Data read successfully:', variable

end subroutine read_netcdf

! ----------------------------------------------------------

subroutine load_latinhc(n_ppe)
use parameters, only: lsample, aero_N_init, Dm_init
use switches, only: wctrl
use namelists, only: irealz, Na_max, Na_min, w_max, w_min, lhs_path, &
                     l_ppe_nevp, l_ppe_condevp, l_ppe_coal, l_ppe_sed, &
                     Dm_min, Dm_max, n_init
use micro_prm, only: n_param_nevp, n_param_condevp, n_param_coal, n_param_sed, npp
use module_lhc, only: write_lhc_nc
  
integer, intent(in) :: n_ppe
character(len=200) :: filename
character(len=100) :: varname
character(len=6) :: n_ppe_str, npp_str
logical :: file_exist

npp = 0
if (l_ppe_nevp) npp = npp + n_param_nevp
if (l_ppe_condevp) npp = npp + n_param_condevp
if (l_ppe_coal) npp = npp + n_param_coal
if (l_ppe_sed) npp = npp + n_param_sed

write(npp_str,'(I0)') npp+n_init
write(n_ppe_str,'(I0)') n_ppe

filename = trim(lhs_path)//'/lhs_out_p'// trim(npp_str)//'_n'//trim(n_ppe_str)//'.nc'
inquire(file = filename, exist = file_exist)
if (.not. file_exist) call write_lhc_nc(npp+n_init, n_ppe)
print*, 'loading LHS:', filename
varname = 'lhs_sample'
call read_netcdf(lsample, filename, varname)

if (Dm_min>0.) then
  Dm_init = lsample(1, irealz) * (Dm_max - Dm_min) + Dm_min
elseif (Na_min>0. .and. w_min>0.) then
  ! draw a aero_N_init between Na_min and Na_max
  aero_N_init(1) = lsample(1, irealz) * (Na_max - Na_min) + Na_min
  ! draw a wctrl(1) between w_min and w_max
  wctrl(1) = lsample(2, irealz) * (w_max - w_min) + w_min
endif

end subroutine load_latinhc
! ----------------------------------------------------------
subroutine get_perturbed_params(params_save, pvalue_mean, pvalue_isd)

use parameters, only: lsample
use micro_prm, only: n_param_nevp, n_param_condevp, n_param_coal, n_param_sed, npp
use namelists, only: irealz, l_ppe_nevp, l_ppe_condevp, l_ppe_coal, l_ppe_sed, deflation_factor, &
  n_init
integer :: ilsample, inevp, icondevp, icoal, ised, iparam_indproc, iparam_allproc, n_param
real, allocatable, dimension(:) :: pvalue_mean, pvalue_isd, params_save
double precision :: nudge_diff

! icondevp = 1; icoal = 1; ised = 1; 
ilsample = 1; inevp = 0; icondevp = 0; icoal = 0; ised = 0
n_param = n_param_nevp + n_param_condevp + n_param_coal + n_param_sed

if (l_ppe_nevp) then
  do iparam_indproc = 1, n_param_nevp
    ilsample = iparam_indproc + inevp
    iparam_allproc = iparam_indproc
    nudge_diff = (lsample(ilsample+n_init,irealz)-.5)*2*pvalue_isd(iparam_allproc)*deflation_factor
    params_save(iparam_allproc) = pvalue_mean(iparam_allproc) + nudge_diff
  enddo
  ! push the starting index of other processes if nevp params are perturbed
  icondevp = ilsample; icoal = ilsample; ised = ilsample
endif

if (l_ppe_condevp) then
  do iparam_indproc = 1, n_param_condevp
    ilsample = iparam_indproc + icondevp
    iparam_allproc = iparam_indproc + n_param_nevp
    nudge_diff = (lsample(ilsample+n_init,irealz)-.5)*2*pvalue_isd(iparam_allproc)*deflation_factor
    params_save(iparam_allproc) = pvalue_mean(iparam_allproc) + nudge_diff
  enddo
  icoal = ilsample; ised = ilsample
endif

if (l_ppe_coal) then
  pvalue_isd(n_param_nevp+n_param_condevp+1) = 10 ! a0coal
  pvalue_isd(n_param_nevp+n_param_condevp+2) = 30 ! mtrans
  pvalue_isd(n_param_nevp+n_param_condevp+3:n_param_nevp+n_param_condevp+n_param_coal)=2 ! the exponents

  do iparam_indproc = 1, n_param_coal
    ilsample = iparam_indproc + icoal
    iparam_allproc = iparam_indproc + n_param_nevp + n_param_condevp
    nudge_diff = (lsample(ilsample+n_init,irealz)-.5)*2*pvalue_isd(iparam_allproc)*deflation_factor
    params_save(iparam_allproc) = pvalue_mean(iparam_allproc) + nudge_diff
  enddo
  ised = ilsample
endif

if (l_ppe_sed) then
  pvalue_isd(n_param_nevp+n_param_condevp+n_param_coal+3:n_param) = 5 ! the exponents
  pvalue_isd(n_param_nevp+n_param_condevp+n_param_coal+1) = 50 ! afall
  pvalue_isd(n_param_nevp+n_param_condevp+n_param_coal+2) = 30 ! mlim
  pvalue_isd(n_param-1) = 30 ! mlim2

  do iparam_indproc = 1, n_param_sed
    ilsample = iparam_indproc + ised
    iparam_allproc = iparam_indproc + n_param_nevp + n_param_condevp + n_param_coal
    nudge_diff = (lsample(ilsample+n_init,irealz)-.5)*2*pvalue_isd(iparam_allproc)*deflation_factor
    params_save(iparam_allproc) = pvalue_mean(iparam_allproc) + nudge_diff
  enddo
endif

  
end subroutine get_perturbed_params
! ----------------------------------------------------------
subroutine get_perturbed_params_custom(params_save, posterior_binmeans, posterior_cdf, n_bins)
  
use parameters, only: lsample
use micro_prm, only: n_param_nevp, n_param_condevp, n_param_coal, n_param_sed, npp
use namelists, only: irealz, l_ppe_nevp, l_ppe_condevp, l_ppe_coal, l_ppe_sed, n_init
integer :: ilsample, inevp, icondevp, icoal, ised, ibin, iparam_allproc, iparam_indproc
double precision, allocatable :: posterior_binmeans(:,:), posterior_cdf(:,:)
real, allocatable :: params_save(:)
double precision :: nudge_diff

ilsample = 1; inevp = 0; icondevp = 0; icoal = 0; ised = 0

! pick a number between 0 and 1 from lhs and find it in the cdf
if (l_ppe_nevp) then
  do iparam_indproc = 1, n_param_nevp
    ilsample = iparam_indproc + inevp
    iparam_allproc = iparam_indproc
    do ibin = 1, n_bins
      if (lsample(ilsample+n_init, irealz) <= posterior_cdf(ibin, ilsample)) then
        params_save(iparam_allproc) = posterior_binmeans(ibin, ilsample)
        exit
      endif
    enddo
  enddo
  icondevp = ilsample; icoal = ilsample; ised = ilsample
endif

if (l_ppe_condevp) then
  do iparam_indproc = 1, n_param_condevp
    ilsample = iparam_indproc + icondevp
    iparam_allproc = iparam_indproc + n_param_nevp
    do ibin = 1, n_bins
      if (lsample(ilsample+n_init, irealz) <= posterior_cdf(ibin, ilsample)) then
        params_save(iparam_allproc) = posterior_binmeans(ibin, ilsample)
        exit
      endif
    enddo
  enddo
  icoal = ilsample; ised = ilsample
endif

if (l_ppe_coal) then
  do iparam_indproc = 1, n_param_coal
    ilsample = iparam_indproc + icoal
    iparam_allproc = iparam_indproc + n_param_nevp + n_param_condevp
    do ibin = 1, n_bins
      if (lsample(ilsample+n_init, irealz) <= posterior_cdf(ibin, ilsample)) then
        params_save(iparam_allproc) = posterior_binmeans(ibin, ilsample)
        exit
      endif
    enddo
  enddo
  ised = ilsample
endif

if (l_ppe_sed) then
  print*, 'n_param_sed', n_param_sed
  do iparam_indproc = 1, n_param_sed
    ilsample = iparam_indproc + ised
    iparam_allproc = iparam_indproc + n_param_nevp + n_param_condevp + n_param_coal
    do ibin = 1, n_bins
      if (lsample(ilsample+n_init, irealz) <= posterior_cdf(ibin, ilsample)) then
        params_save(iparam_allproc) = posterior_binmeans(ibin, ilsample)
        exit
      endif
    enddo
  enddo
endif

end subroutine get_perturbed_params_custom

! ----------------------------------------------------------

subroutine get_posterior(posterior_binmeans, posterior_cdf, n_bins, npp)
use namelists, only: custom_dens_path, custom_bins_path
use csv_module

double precision, allocatable, dimension(:,:) :: posterior_dens, posterior_binedges, posterior_binmeans, posterior_cdf
integer :: ilsample, ibin, npp, iparam
logical :: stat_ok
type(csv_file) :: dens_csv, bins_csv
  
allocate(posterior_dens(n_bins, npp))
allocate(posterior_binedges(n_bins+1, npp))
allocate(posterior_binmeans(n_bins, npp))
allocate(posterior_cdf(n_bins, npp))

! read the file
call dens_csv%read(trim(custom_dens_path),header_row=1,status_ok=stat_ok)
call bins_csv%read(trim(custom_bins_path),header_row=1,status_ok=stat_ok)

! get MCMC posterior distributions
do iparam = 1, npp
  do ibin = 1, n_bins+1
    call bins_csv%get(ibin,iparam,posterior_binedges(ibin, iparam),stat_ok)
    if (ibin > n_bins) cycle
    call dens_csv%get(ibin,iparam,posterior_dens(ibin, iparam),stat_ok)
  enddo
enddo
posterior_binmeans = (posterior_binedges(2:n_bins+1, :) + posterior_binedges(1:n_bins, :))/2

! build the cdf
posterior_cdf(1, :) = posterior_dens(1, :)
do ibin = 2, n_bins
  posterior_cdf(ibin, :) = posterior_cdf(ibin-1, :) + posterior_dens(ibin, :)
enddo

! in case CDF doesn't sum up to 1
do iparam = 1, npp
  posterior_cdf(:, iparam) = posterior_cdf(:, iparam)/posterior_cdf(n_bins, iparam)
enddo

end subroutine get_posterior
! ----------------------------------------------------------
subroutine get_perturbed_params_pymc(params_save, pymc_filedirs)

use parameters, only: pymc_filedirs_type, lsample
use micro_prm, only: n_param_nevp, n_param_condevp, n_param_coal, n_param_sed
use namelists, only: n_init, l_ppe_nevp, l_ppe_condevp, l_ppe_coal, l_ppe_sed, irealz

real, allocatable, dimension(:) :: params_save
real :: params_nevp(n_param_nevp), params_condevp(n_param_condevp), params_coal(n_param_coal), &
  params_sed(n_param_sed)
type(pymc_filedirs_type) :: pymc_filedirs
real(8) :: large_val = 1.e3
real(8), allocatable, dimension(:) :: max_std
integer :: n_param, inevp, icondevp, icoal, ised

inevp = n_init; icondevp = n_init; icoal = n_init; ised = n_init;
n_param = n_param_nevp + n_param_condevp + n_param_coal + n_param_sed


if (l_ppe_nevp) then
  allocate(max_std(n_param_nevp))
  max_std(:) = large_val
  call load_pymc(params_nevp, max_std, lsample(inevp+1:inevp+n_param_nevp, irealz), n_param_nevp, pymc_filedirs%nevp_dir)
  deallocate(max_std)
  params_save(1:n_param_nevp) = params_nevp
  icondevp = n_param_nevp; icoal = n_param_nevp; ised = n_param_nevp
endif

if (l_ppe_condevp) then
  allocate(max_std(n_param_condevp))
  max_std(:) = large_val
  call load_pymc(params_condevp, max_std, lsample(icondevp+1:icondevp+n_param_condevp, irealz), n_param_condevp, pymc_filedirs%condevp_dir)
  deallocate(max_std)
  params_save(n_param_nevp+1:n_param_nevp+n_param_condevp) = params_condevp
  icoal = icondevp + n_param_condevp; ised = icondevp + n_param_condevp
endif


if (l_ppe_coal) then
  allocate(max_std(n_param_coal))
  max_std = [10,30,2,2,2,2,2,2,2,2,2,2]
  call load_pymc(params_coal, max_std, lsample(icoal+1:icoal+n_param_coal, irealz), n_param_coal, pymc_filedirs%coal_dir)
  deallocate(max_std)
  params_save(n_param_nevp+n_param_condevp+1:n_param_nevp+n_param_condevp+n_param_coal) = &
    params_coal
  ised = icoal + n_param_coal
endif

if (l_ppe_sed) then
  allocate(max_std(n_param_sed))
  max_std = [50,30,2,2,2,2,2,2,2,2,30,2]
  call load_pymc(params_sed, max_std, lsample(ised+1:ised+n_param_sed, irealz), n_param_sed, pymc_filedirs%sed_dir)
  deallocate(max_std)
  params_save(n_param_nevp+n_param_condevp+n_param_coal+1:n_param) = &
    params_sed
endif


end subroutine get_perturbed_params_pymc
! ----------------------------------------------------------
subroutine load_pymc(params_loaded, max_std, trunc_lhs, n_param, pymc_filedir)

use netcdf

real, dimension(n_param) :: params_loaded
real(8), allocatable, dimension(:,:,:) :: params_proc_sampled
real(8), allocatable, dimension(:,:) :: params_read
character(len=200) :: pymc_filedir
character(len=100) :: group_name
character(len=10) :: varname
integer :: ncid, varid, cid, did, pid
integer :: nchain, ndraw, n_param
real(8), allocatable :: posterior(:,:,:), a0coal(:,:)
real(8) :: max_std(n_param), trunc_lhs(n_param)
character(len=20), allocatable, dimension(:) :: vnames

! put all the offline posterior samples to the `posterior` variable
group_name = 'posterior'
vnames = get_pnames_from_netcdf(pymc_filedir, group_name)

if (size(vnames) /= n_param) then
  print*, size(vnames), n_param
  print*, vnames
  STOP 'inconsistent n_param'
endif

! allocate(params_proc_sampled(30000,8,12))

do i = 1, n_param
  call read_netcdf(params_read, pymc_filedir, trim(vnames(i)), group_name)
  if (.not. allocated(params_proc_sampled)) then
    ndraw = size(params_read,1)
    nchain = size(params_read,2)
    allocate(params_proc_sampled(ndraw, nchain, n_param))
  endif
  params_proc_sampled(:,:,i) = params_read
enddo

call draw_mvnormal(params_proc_sampled, nchain, ndraw, n_param, max_std, trunc_lhs, params_loaded)

end subroutine load_pymc
! ----------------------------------------------------------
subroutine draw_mvnormal(posterior, nchain, ndraw, n_param, max_std, trunc_lhs, theta)
  use namelists, only: irealz

  implicit none
  integer, intent(in)  :: nchain, ndraw, n_param
  real(8),    intent(in)  :: posterior(nchain, ndraw, n_param)
  real,    intent(out) :: theta(n_param)

  real(8) :: mean(n_param), cov(n_param,n_param), L(n_param,n_param)
  real(8) :: z(n_param, 1), perturbed_val(n_param, 1)
  integer :: total, i, j, k, info
  real(8) :: u1, u2, inflStd(n_param), s(n_param), max_std(n_param), &
    inflCov(n_param, n_param), cappedCov(n_param,n_param), trunc_lhs(n_param)

  !— 1) collapse chains, compute mean —
  total = nchain*ndraw
  mean = 0.0
  do i=1,nchain; do j=1,ndraw
     mean = mean + posterior(i,j,:)
  end do; end do
  mean = mean / total

  !— 2) compute covariance —
  cov = 0.0
  do i=1,nchain; do j=1,ndraw
    cov = cov + matmul(reshape(posterior(i,j,:)-mean,(/n_param,1/)), &
                       reshape(posterior(i,j,:)-mean,(/1,n_param/)))
  end do; end do
  cov = cov / (total-1)
  cappedCov = cov

  !- 3) inflate covariance
  inflCov = cov*ndraw
  do i = 1, n_param
    if (inflCov(i,i) > 0.) then
      inflStd(i) = sqrt(inflCov(i,i))
      s(i)       = min(1., max_std(i) / inflStd(i))
    else
      inflStd(i) = 0.
      s(i)       = 0.
    end if
  end do
  do i = 1, n_param
    do j = 1, n_param
      cappedCov(i,j) = s(i) * s(j) * inflCov(i,j)
      ! if (i/=j) cappedCov(i,j) = 0.
    end do
  end do

  !— 4) Cholesky factorization cov = L·Lᵀ (lower triangular) —
  L = cappedCov
  call dpotrf('L', n_param, L, n_param, info)
  if (info /= 0) then
     print *, "Cholesky failed, info=", info
     stop
  end if

  !— 5) transform: theta = mean + L·z —
  z(:, 1) = (trunc_lhs-.5)*2!*pvalue_isd(iparam_allproc)*deflation_factor
  perturbed_val = matmul(L, z)/5
  ! print*, perturbed_val
  ! print*, 'cov', cov
  ! print*, 'cappedCov', cappedCov
  ! stop
  theta = sngl(mean + perturbed_val(:, 1))
  
end subroutine draw_mvnormal
! ----------------------------------------------------------
function get_pnames_from_netcdf(filename, group_name) result(vnames)

  use global_fun, only: contains_any
  use netcdf

  integer :: grpid        ! Group ID
  integer :: ndims, nvars, natts, unlimdimid, status
  integer :: n_param, iv
  integer, allocatable, dimension(:) :: varids
  logical, allocatable, dimension(:) :: l_pnames
  character(len=200) :: filename
  character(len=20)  :: group_name
  character(len=20), allocatable, dimension(:) :: vnames, vn_skip, vnames_all


  ! Open the NetCDF file (read-only mode)
  status = nf90_open(filename, nf90_nowrite, ncid)
  if (status /= nf90_noerr) then
    print *, 'Error opening NetCDF file in get_pnames_from_netcdf: ', trim(nf90_strerror(status))
    stop
  end if

  ! Get the group ID
  status = nf90_inq_grp_ncid(ncid, group_name, grpid)
  if (status /= NF90_NOERR) stop "Group "//group_name//" not found"

  allocate(vn_skip(5))
  vn_skip(1) = "chain"
  vn_skip(2) = "draw"
  vn_skip(3) = "sum_sq"
  vn_skip(4) = "log_pppd"
  vn_skip(5) = "neg"

  status = nf90_inquire(grpid, ndims, nvars, natts, unlimdimid)
  allocate(varids(nvars))
  allocate(vnames_all(nvars))
  allocate(l_pnames(nvars))
  status = nf90_inq_varids(grpid, nvars, varids)

  !--- loop and print each name ---
  do i = 1, nvars
    status = NF90_INQUIRE_VARIABLE(grpid, varids(i), name=vnames_all(i))
    l_pnames(i) = .not. contains_any(vnames_all(i), vn_skip)
    if (status /= NF90_NOERR) then
      print *, ' error reading var ', status, varids(i), trim(nf90_strerror(status))
      stop
    end if
  end do

  n_param = count(l_pnames)
  allocate(vnames(n_param))
  iv = 1
  do i = 1, nvars
    if (l_pnames(i)) then
      vnames(iv) = vnames_all(i)
      iv = iv + 1
    endif
  enddo

end function get_pnames_from_netcdf
! ----------------------------------------------------------

end module ppe_fun
