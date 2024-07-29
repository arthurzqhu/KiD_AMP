module module_mp_amp
  use micro_prm, only: pmomsc, pmomsr, dnc_def, dnr_def, npm, n_cat, nkr, &
                       diams, pio6rw, ipio6rw, pdiams, total_m3_th, max_nbins, col, &
                       nuterm31, nuterm32, nutermx1, nutermx2, nutermw1, nutermw2, tORf, &
                       D_min, D_max, tempvar_debug, rtemp_debug
  use module_hujisbm, only: xl
  use mphys_tau_bin_declare, only: xkgmean
  use namelists, only: bintype, l_truncated, l_init_test
  use parameters, only: split_bins, h_shape
  use physconst, only: pi
  use global_fun, only: gammq
  implicit none
  
  integer :: icat, scl, ecl, lcl_catbound1(2), lcl_catbound2(2), imomw, imomx, imomy, imomz, &
    imom_consv
  double precision :: M3p(2), M0p, Mwp, Mxp(2), Myp, Mzp, relaxw, &
    relaxx, relaxy, nu_def, m1frac, dummy, nuterm_consv1, nuterm_consv2, relax0
  integer, parameter :: n_mom_diag=10
  double precision :: dntab(50,50,2) !, debug_arr(nkr,5)
  logical :: l_improved
  type, private :: guess_param
    double precision :: & 
      ! guesses of dn, nu for two modes l1 and l2, not the most concise implementation
      ! but probably the most readable
      l1prm(2), &
      l2prm(2), & 
      ! error information
      error(4), sqerr, &
      ! mass fraction of liquid category 1
      l1_mfrac
  end type

  ! constants
  double precision, parameter :: dn_min=1e-9, dn_max=1e-3, nu_max=100, nu_min=-1

contains

subroutine  invert_moments(mom_pred, gam_param, amp_distm, amp_distn, flags, sqerr)

double precision, dimension(npm, n_cat), intent(inout) :: mom_pred
double precision, dimension(2,2), intent(inout) :: gam_param
real, dimension(nkr), intent(out) :: amp_distm, amp_distn
double precision, intent(out) :: flags(n_cat)
double precision, dimension(nkr) :: amp_distm_dble, amp_distn_dble
double precision :: minmxm3, maxmxm3, mw, mx, m3, m0
type(guess_param) :: guess_best, guess_prev, guess_try1, guess_try2, &
                     guess_try3, guess_start, guess_old, guess_int
double precision :: dn(2), rand1(4), rand2(4), Mconsv_M3, tol, &
                    err_threshold, sqerr_threshold, sqerr, t3sqerr, t3err(4)
integer :: exit_signal, ntry, itry, ibguess,i,i1,i2
double precision :: dn1t(201), dn2t(251), ers(201,251), gstep, ers1(201,251), ers2(201,251), &
  gstep1, gstep2

err_threshold = 100
sqerr_threshold = 100

amp_distm_dble(:) = 0.
amp_distn_dble(:) = 0.
exit_signal = 0
l_improved = .false.
ibguess = 0

! get the index boundaries for the two categories
lcl_catbound1 = (/1,split_bins+1/)
lcl_catbound2 = (/split_bins,nkr/)

if ( n_cat == 2 ) then
  imomx = pmomsc(2)
  if ( npm == 3 ) imomy = pmomsc(3)
elseif ( n_cat == 1 ) then
  imomw = pmomsc(3)
  imomx = pmomsc(4)
  if ( npm >= 5 ) then
    imomy = pmomsc(5)
  endif
  if ( npm >= 6 ) then
    imomz = pmomsc(6)
  endif
endif

! previous guesses

if ( any(gam_param(1,:) .le. 1e-20) ) then
  !initiate
  gam_param(1,1) = dnc_def
  gam_param(1,2) = dnr_def
  gam_param(2,1) = h_shape(1)
  gam_param(2,2) = h_shape(2)
endif

if ( n_cat==1 ) then

  ! initialize variables
  relax0 = 1.
  relaxw = 1.
  if (npm >= 5) then
    relaxx = 1.
  endif
  if (npm >= 6) then
    relaxy = 1.
  endif

  tol = 1e-4
  exit_signal = 0
  
  M3p(1) = mom_pred(1,1)
  M0p = mom_pred(2,1)
  Mwp = mom_pred(3,1)
  Mxp(1) = mom_pred(4,1)

  if (npm >= 5) then
    Myp = mom_pred(5,1)
  endif
  if (npm >= 6) then
    Mzp = mom_pred(6,1)
  endif

  if (M3p(1) < total_m3_th) then
    amp_distm(:) = 0.
    amp_distn(:) = 0.
    return
  endif

  dn(1) = gam_param(1,1)
  dn(2) = gam_param(1,2)

  guess_prev%l1prm(1) = dn(1)
  guess_prev%l2prm(1) = dn(2)
  guess_prev%l1prm(2) = gam_param(2,1) ! nu1
  guess_prev%l2prm(2) = gam_param(2,2) ! nu2
  ! print*, guess_prev
  ! stop

  ! a moment ratio needed for moment conservation
  if ( npm == 4 ) then
    Mconsv_M3 = Mxp(1)/M3p(1)
    imom_consv = imomx
    nuterm_consv1 = nutermx1
    nuterm_consv2 = nutermx2
  elseif ( npm == 5 ) then
    Mconsv_M3 = Myp/M3p(1)
    imom_consv = imomy
  elseif ( npm == 6 ) then
    Mconsv_M3 = Mzp/M3p(1)
    imom_consv = imomz
  endif

  ! calculate the error from the initial guesses
  CALL calcerr_1cat(guess_prev, Mconsv_M3)
  guess_best = guess_prev
  guess_int = guess_best

  if (guess_prev%l1prm(1).ne.0 .or. guess_prev%l2prm(1).ne.0) then
    ! put initial guesses into find_params
    guess_try1 = guess_prev
    if ( npm == 4 ) then
      call find_params(fcn_4m, guess_try1, tol, exit_signal)
    elseif ( npm == 5 ) then
      call find_params(fcn_5m, guess_try1, tol, exit_signal)
    elseif ( npm == 6 ) then
      call find_params(fcn_6m, guess_try1, tol, exit_signal)
    endif
    call update_param_guess(guess_best, guess_try1, mark_improve=.false.)
  endif

  ! guess_start is fixed so that the gradient descent algorithm doesn't
  ! drift us away from the global minimum
  guess_start = guess_best
  
  ! if guess_try1 is not good enough then start random guess
  if ( exit_signal .ne. 1 ) then
    ! snap back to default if outside a reasonable range

    if ( guess_best%l1prm(1)>=dn_max .or. guess_best%l1prm(1)<=dn_min ) then
      ! print*, 'dnc', guess_best%l1prm(1)
      guess_best%l1prm(1) = dnc_def
    endif
    if ( guess_best%l2prm(1)>=dn_max .or. guess_best%l2prm(1)<=dn_min ) then
      ! print*, 'dnr', guess_best%l2prm(1)
      guess_best%l2prm(1) = dnr_def
    endif
    if ( npm >= 5 ) then
      if ( guess_best%l1prm(2)>=nu_max .or. guess_best%l1prm(2)<=nu_min ) guess_best%l1prm(2) = h_shape(1)
    endif
    if ( npm >= 6 ) then
      if ( guess_best%l2prm(2)>=nu_max .or. guess_best%l2prm(2)<=nu_min ) guess_best%l2prm(2) = h_shape(2)
    endif


    ntry = 0
    call random_init(.true., .true.)
    do itry = 1,10000
       relax0 = 1.
       relaxw = 1.
       if ( npm >= 5 ) relaxx = 1.
       if ( npm >= 6 ) relaxy = 1.
      ! print*, 'itry', itry
      call RANDOM_NUMBER(rand1); call RANDOM_NUMBER(rand2)
      guess_try2%l1prm(1) = guess_start%l1prm(1)* &
        exp(sqrt(-2*0.25*log(rand1(1)))*cos(2.*3.14159*rand2(1)))
      guess_try2%l2prm(1) = guess_start%l2prm(1)* &
        exp(sqrt(-2*0.25*log(rand1(2)))*cos(2.*3.14159*rand2(2)))

      if ( npm == 4 ) then
        guess_try2%l1prm(2) = h_shape(1)
        guess_try2%l2prm(2) = h_shape(2)
      endif
      if ( npm >= 5 ) then
        guess_try2%l1prm(2) = guess_start%l1prm(2)+ &
          sqrt(-2*9*log(rand1(3)))*cos(2.*3.14159*rand2(3))
      endif
      if ( npm >= 6 ) then
        guess_try2%l2prm(2) = guess_start%l2prm(2)+ &
          sqrt(-2*9*log(rand1(4)))*cos(2.*3.14159*rand2(4))
      endif

      call calcerr_1cat(guess_try2, Mconsv_M3) ! make sure the sqerr is not nan

      call update_param_guess(guess_best, guess_try2, mark_improve=.true.)

      if (guess_try2%sqerr < 1) then
        ! start the gradient descent if the guess is good enough
        ntry = ntry + 1

        if ( npm == 4 ) then
          call find_params(fcn_4m, guess_try2, tol, exit_signal)
        elseif ( npm == 5 ) then
          call find_params(fcn_5m, guess_try2, tol, exit_signal)
        elseif ( npm == 6 ) then
          call find_params(fcn_6m, guess_try2, tol, exit_signal)
        endif

        call update_param_guess(guess_best, guess_try2, mark_improve=.true.)
        if ( exit_signal .eq. 1 ) goto 10 ! good enough -> stop searching
      endif

      !this may not be necessary, we'll see
      if (abs(guess_try2%error(2))<0.25 .and. guess_try2%sqerr<sqerr_threshold) then
        guess_int = guess_try2
        ! ibguess = 1
      endif

      if ( ntry >= 5 .or. itry >=5000) then
        !If we've tried a lot and failed, then let's try relaxing the error
        !on momw but not M0
        guess_try3 = guess_int
        if ( npm == 4 ) then
          call find_params(fcn_4m, guess_try3, tol, exit_signal)
        elseif ( npm == 5 ) then
          call find_params(fcn_5m, guess_try3, tol, exit_signal)
        elseif ( npm == 6 ) then
          call find_params(fcn_6m, guess_try3, tol, exit_signal)
        endif

        ! if ( abs(guess_try3%error(2))>tol ) relaxw = max(tol/abs(guess_try3%error(2)), tol)
        if ( abs(guess_try3%error(1))>tol ) relax0 = max(tol/abs(guess_try3%error(1)), tol)
        if ( npm >= 5 .and. abs(guess_try3%error(3))>tol ) &
          relaxx = max(tol/abs(guess_try3%error(3)), tol)
        if ( npm >= 6 .and. abs(guess_try3%error(4))>tol ) &
          relaxy = max(tol/abs(guess_try3%error(4)), tol)

        if ( npm == 4 ) then
          call find_params(fcn_4m, guess_try3, tol, exit_signal)
        elseif ( npm == 5 ) then
          call find_params(fcn_5m, guess_try3, tol, exit_signal)
        elseif ( npm == 6 ) then
          call find_params(fcn_6m, guess_try3, tol, exit_signal)
        endif

        t3err(:) = 0.
        t3err(1:2) = guess_try3%error(1:2)
        ! t3err(2) = guess_try3%error(2)/relaxw
        t3err(1) = guess_try3%error(1)/relax0
        if ( npm >= 5 ) then
          t3err(3) = guess_try3%error(3)/relaxx
        endif
        if ( npm >= 6 ) then
          t3err(4) = guess_try3%error(4)/relaxy
        endif
        t3sqerr = sqrt(sum(t3err(1:(npm-2))**2))
        ! call update_param_guess(guess_best, guess_try3, mark_improve=.true.)

        if (t3sqerr <= guess_best%sqerr .or. (abs(guess_best%error(2))>.01 .and. &
           abs(t3err(2))<.01 )) then
           guess_best = guess_try3
           l_improved = .true.
        endif

        if ( exit_signal .eq. 1 .or. exit_signal .eq. 4 .or. abs(guess_best%error(2))<=.01 &
          .or. ntry >= 10 ) goto 10
       ! ibguess = 0
       ! relaxw = 1.
       if ( npm >= 5 ) relaxx = 1.
       if ( npm >= 6 ) relaxy = 1.

      endif
    enddo
   endif

  10 continue

  if (l_improved) guess_start = guess_best

  call calcdist(guess_start, amp_distm_dble)

  if ( guess_start%l1prm(1)>=dn_max .or. guess_start%l1prm(1)<=dn_min ) then
    guess_start%l1prm(1) = dnc_def
  endif
  if ( guess_start%l2prm(1)>=dn_max .or. guess_start%l2prm(1)<=dn_min ) then
    guess_start%l2prm(1) = dnr_def
  endif
  if ( npm >= 5 ) then
    if ( guess_start%l1prm(2)>=nu_max .or. guess_start%l1prm(2)<=nu_min ) guess_start%l1prm(2) = h_shape(1)
  endif
  if ( npm >= 6 ) then
    if ( guess_start%l2prm(2)>=nu_max .or. guess_start%l2prm(2)<=nu_min ) guess_start%l2prm(2) = h_shape(2)
  endif

  gam_param(1:2,1) = guess_start%l1prm(1:2)
  gam_param(1:2,2) = guess_start%l2prm(1:2)
  sqerr = guess_start%sqerr

  if (guess_start%sqerr>tol) then
    flags = 1
  else
    flags = 0
  endif

! --------------------- two category AMP -------------------
elseif ( n_cat==2 ) then
  ! for 2M (3M will be implemented later if necessary)
  ! initialize and zeros variables
  tol = 1.d-8
  guess_best%l1prm = (/dble(dnc_def), dble(h_shape(1))/)
  guess_best%l2prm = (/dble(dnr_def), dble(h_shape(2))/)

  do icat=1,2
    M3p(icat) = mom_pred(1, icat)
    Mxp(icat) = mom_pred(2, icat)

    if (icat == 1) nu_def = h_shape(1)
    if (icat == 2) nu_def = h_shape(2)

    if (M3p(icat) < total_m3_th) cycle
    if (Mxp(icat) <= 0.) cycle

    scl = lcl_catbound1(icat)
    ecl = lcl_catbound2(icat)
    flags(icat) = 0

    Mconsv_M3 = Mxp(icat)/M3p(icat)

    ! force Mconsv_M3 into the proper range if it's not already
    minmxm3 = min(diams(scl)**(imomx-3.), diams(ecl)**(imomx-3.))
    maxmxm3 = max(diams(scl)**(imomx-3.), diams(ecl)**(imomx-3.))
    if ( Mconsv_M3 <= minmxm3 ) then
      Mconsv_M3 = minmxm3*1.0001
      mom_pred(2, icat) = Mconsv_M3 * mom_pred(1, icat)
      flags(icat) = 1
    elseif ( Mconsv_M3 >= maxmxm3 ) then
      Mconsv_M3 = maxmxm3*0.9999
      mom_pred(2, icat) = Mconsv_M3 * mom_pred(1, icat)
      flags(icat) = 1
    endif

    dn(icat) = (Mconsv_M3*gamma(nu_def+3.)/gamma(nu_def+imomx))**(1./(imomx-3.))

    if (dn(icat) .ne. dn(icat) .or. dn(icat) > huge(dn(icat))) then
      print*, 'Mxp(icat), M3p(icat)', Mxp(icat), M3p(icat)
      print*, 'Mconsv_M3', Mconsv_M3
      print*, 'nu_def', nu_def
      print*, 'imomx', imomx
    endif

    if (icat == 1) guess_best%l1prm(1) = dn(1)
    if (icat == 2) guess_best%l2prm(1) = dn(2)

    call find_params(fcn_2m, guess_best, tol, exit_signal)

    if ( abs(guess_best%error(1)) > tol ) flags(icat) = 1
  enddo

  call calcdist(guess_best, amp_distm_dble)

  if ( guess_best%l1prm(1)>=dn_max .or. guess_best%l1prm(1)<=dn_min ) then
    guess_best%l1prm(1) = dnc_def
  endif

  if ( guess_best%l2prm(1)>=dn_max .or. guess_best%l2prm(1)<=dn_min ) then
    guess_best%l2prm(1) = dnr_def
  endif

  gam_param(1:2,1) = guess_best%l1prm(1:2)
  gam_param(1:2,2) = guess_best%l2prm(1:2)
  sqerr = guess_best%sqerr

endif

amp_distm = real(amp_distm_dble)
if (bintype .eq. 'tau') amp_distn = amp_distm/xkgmean

! print*, 'itry', itry
! print*, 'gam_param', gam_param(1,:)
! stop
! print*, 'error', guess_start%error
! print*, 'relaxw', relaxw

! print*, 'sum(amp_distm)',sum(amp_distm)*col*ipio6rw
! print*, 'amp_distm', amp_distm
! print*, 'guess_best%l1prm(1:2)', guess_best%l1prm(1:2)
! print*, 'guess_best%l2prm(1:2)', guess_best%l2prm(1:2)

end subroutine invert_moments

! --------------- gamma parameters -> discrete distribution -----------------
subroutine calcdist(guess_best, mass_dist)

type(guess_param), intent(in) :: guess_best
integer :: ntry, lcl
double precision, dimension(nkr), intent(out) :: mass_dist
double precision :: md_part(nkr,2), dn(2), nu(2)
double precision :: m3, Mconsv_M3, exptermc, exptermr, n01, n02, ml1, ml2, m0, mw, mx, n0, &
  exptermrc(nkr)

dn(1) = guess_best%l1prm(1)!*1e-6
dn(2) = guess_best%l2prm(1)!*1e-6
nu(1) = guess_best%l1prm(2)
nu(2) = guess_best%l2prm(2)

mass_dist(1:nkr) = 0.

if (n_cat == 1) then

  if ( npm == 4 ) then
    Mconsv_M3 = Mxp(1)/M3p(1)
  elseif ( npm == 5 ) then
    Mconsv_M3 = Myp/M3p(1)
  elseif ( npm == 6 ) then
    Mconsv_M3 = Mzp/M3p(1)
  endif

  ! print*, 'dn1, dn2', dn
  ! print*, 'm1frac', guess_best%l1_mfrac

  if (l_truncated) then
    call incgamma_1cat(dn(1), nu(1), dn(2), nu(2), Mconsv_M3, mass_dist)
    mass_dist = mass_dist*M3p(1)
  else

    ml1 = M3p(1)*guess_best%l1_mfrac
    ml2 = M3p(1) - ml1
    n01 = ml1/gamma(nu(1)+3)*log(2.)/3*pio6rw
    n02 = ml2/gamma(nu(2)+3)*log(2.)/3*pio6rw


    do lcl=1,nkr
      if (ml1>0. .and. dn(1)>0.) then
        exptermc=exp(-1.*diams(lcl)/dn(1))
        mass_dist(lcl) = n01*exptermc*(diams(lcl)/dn(1))**(nu(1)+3)
      endif
      if (ml2>0. .and. dn(2)>0.) then
        exptermr=exp(-1.*diams(lcl)/dn(2))
        mass_dist(lcl) = mass_dist(lcl) + n02*exptermr*(diams(lcl)/dn(2))**(nu(2)+3)
      endif
    enddo
    mass_dist = mass_dist/col
  endif

  return

elseif (n_cat == 2) then
  ! should work for both 2M and 3M
  md_part(1:nkr,1:2) = 0.
  do icat = 1,2

    if (M3p(icat) <= total_m3_th) goto 112 ! skip loop. can also do cycle...?

    if (l_truncated) then

      scl = lcl_catbound1(icat)
      ecl = lcl_catbound2(icat)
      call incgamma_2cat(nu(icat), dn(icat), md_part(:, icat))
      call calc_mom(m3, md_part(:,icat), 3)
      if (m3==0. .or. dn(icat)==0.) then
        print*, 'scl, ecl', scl, ecl
        print*, 'md_part', md_part(:,icat)
        print*, 'cannot find rxfinal and mass = ',M3p(icat)*pio6rw, Mxp(icat)
        print*, 'm3, nu, dn, icat', m3, nu, dn, icat!, dnbounds(:,icat)
        stop
      endif
      md_part(:, icat) = md_part(:, icat)*M3p(icat)/m3
    else

      n0 = M3p(icat)/gamma(nu(1)+3)*log(2.)/3*pio6rw
      do lcl=1,nkr
        if (M3p(icat)>0. .and. dn(1)>0.) then
          exptermrc(lcl)=exp(-1.*diams(lcl)/dn(icat))
          md_part(lcl,icat) = n0*exptermrc(lcl)*(diams(lcl)/dn(icat))**(nu(icat)+3)/col
        endif
      enddo

    endif


112 continue
  enddo

  mass_dist(1:nkr) = md_part(1:nkr, 1) + md_part(1:nkr, 2)

  return

endif
  
end subroutine calcdist

! --------------- incgamma_2cat -----------------
subroutine incgamma_2cat(nu, dn, md)
! construct a normalized partial gamma distribution (md) given nu and dn

double precision, intent(in) :: nu, dn
double precision, intent(out) :: md(nkr)
double precision :: expterm, term1(max_nbins), term2
integer :: lcl

md(:) = 0

if (scl >= ecl) then
   print*, 'scl, ecl', scl, ecl
   stop 'scl ecl not initiated correctly'
endif

  do lcl=scl,ecl !Loop over bins
     term1(lcl)=-diams(lcl)/dn+(nu+3.)*log(diams(lcl)/dn)
  enddo
  term2 = -maxval(term1(scl:ecl))

  do lcl=scl,ecl
     md(lcl)=exp(term1(lcl)+term2)
  enddo

end subroutine incgamma_2cat

! --------------- incgamma_1cat -----------------
subroutine incgamma_1cat(dn1, nu1, dn2, nu2, Mconsv_M3, md)
double precision, intent(in) :: Mconsv_M3
double precision, intent(out) :: md(nkr)
double precision, intent(inout) :: dn1, nu1, dn2, nu2
double precision :: md1(nkr), md2(nkr), m31, m32, mcs_1, mcs_2, gterm1a(nkr), gterm1b &
                  , gterm2a(nkr), gterm2b
double precision :: infinity, dummy, frac
integer :: lcl

infinity = HUGE(dummy)

! make sure dn1 and dn2 are in the right order
if (dn1 > dn2) then
  dummy = dn1
  dn1 = dn2
  dn2 = dummy
endif

! get the base function shapes
do lcl = 1,nkr
  gterm1a(lcl)=-diams(lcl)/dn1+(nu1+3)*log(diams(lcl)/dn1)
  gterm2a(lcl)=-diams(lcl)/dn2+(nu2+3)*log(diams(lcl)/dn2)
enddo

gterm1b = maxval(gterm1a(1:nkr))*(-1.)
gterm2b = maxval(gterm2a(1:nkr))*(-1.)

do lcl=1,nkr
  md1(lcl)=exp(gterm1a(lcl)+gterm1b)
  md2(lcl)=exp(gterm2a(lcl)+gterm2b)
enddo

! debug_arr(:,1) = md1
! debug_arr(:,2) = md2

call calc_mom(m31, md1, 3)
call calc_mom(m32, md2, 3)

if (m31.eq.0) then
  mcs_1=0.
else
  call calc_mom(mcs_1, md1, imom_consv)
  mcs_1 = mcs_1/m31
endif

if (m32.eq.0) then
  mcs_2=0.
else
  call calc_mom(mcs_2, md2, imom_consv)
  mcs_2 = mcs_2/m32
endif

if (dn1 <= 0. .or. m31>infinity .or. mcs_1>infinity) then
  m1frac=0.
  mcs_1=0;m31=1
  md = md2/m32
  ! dn1 = dnc_def
elseif (dn2 <= 0. .or. m32>infinity .or. mcs_2>infinity) then
  m1frac=1.
  mcs_2=0;m32=1
  md = md1/m31
  ! dn2 = dnr_def
else
  m1frac=(Mconsv_M3-mcs_2)/(mcs_1-mcs_2)
  if (m1frac<0) then
    m1frac=0.
    ! dn1 = dnc_def
  endif
  if (m1frac>1) then
    m1frac=1.
    ! dn2 = dnr_def
  endif
  md = md1/m31*m1frac+ md2/m32*(1.-m1frac)
endif

return
end subroutine incgamma_1cat

! --------------- get_mom_2m -------------------------
subroutine get_mom_2m(nu,dn,m3,mx)
double precision, intent(in) :: nu, dn
double precision, intent(out) :: m3, mx
double precision :: nuterm3, nutermx

if (icat == 1) then
  nuterm3 = nuterm31
  nutermx = nutermw1 ! this is not a typo
else
  nuterm3 = nuterm32
  nutermx = nutermw2 ! this is not a typo
endif

m3 = dn**3*nuterm3
mx = dn**imomx*nutermx
end subroutine get_mom_2m

! --------------- get_mom_from_param -----------------
subroutine get_mom_from_param(dn1, dn2, nu1, nu2, Mconsv_M3, m3, m0, mw)
double precision, intent(in) :: Mconsv_M3
double precision, intent(out) :: m3, m0, mw!, mx, my, mz
double precision, intent(inout) :: dn1, dn2
double precision :: md1(nkr), md2(nkr), m31, m32, mcs_1, mcs_2
double precision :: infinity, mw1frac, mw2frac, mcsm31, mcsm32, mcsmw1, mcsmw2, M3Mw, m3mw1, m3mw2
double precision :: m3f1, m3f2, mwf1, mwf2, nu1, nu2, n01, n02, mw1, mw2, m01, m02
double precision :: gincf_01, gincf_02, gincf_31, gincf_32, gincf_w1, gincf_w2, gincf_csv1, gincf_csv2
integer :: lcl

if (npm > 4) then
   nuterm31 = gamma(nu1+3)/gamma(nu1)
   nuterm32 = gamma(nu2+3)/gamma(nu2)
   nutermw1 = gamma(nu1+imomw)/gamma(nu1)
   nutermw2 = gamma(nu2+imomw)/gamma(nu2)
   nuterm_consv1 = gamma(nu1+imom_consv)/gamma(nu1)
   nuterm_consv2 = gamma(nu2+imom_consv)/gamma(nu2)
endif

! make sure dn1 and dn2 are in the right order
if (dn1 > dn2) then
  dummy = dn1
  dn1 = dn2
  dn2 = dummy
endif

gincf_01   = 1. !gammq(nu1, D_min/dn1)-gammq(nu1, D_max/dn1)
gincf_02   = 1. !gammq(nu2, D_min/dn2)-gammq(nu2, D_max/dn2)
gincf_31   = 1. !gammq(nu1+3, D_min/dn1)-gammq(nu1+3, D_max/dn1)
gincf_32   = 1. !gammq(nu2+3, D_min/dn2)-gammq(nu2+3, D_max/dn2)
gincf_w1   = 1. !gammq(nu1+imomw, D_min/dn1)-gammq(nu1+imomw, D_max/dn1)
gincf_w2   = 1. !gammq(nu2+imomw, D_min/dn2)-gammq(nu2+imomw, D_max/dn2)
gincf_csv1 = 1. !gammq(nu1+imom_consv, D_min/dn1)-gammq(nu1+imom_consv, D_max/dn1)
gincf_csv2 = 1. !gammq(nu2+imom_consv, D_min/dn2)-gammq(nu2+imom_consv, D_max/dn2)

infinity = HUGE(dummy)
M3Mw = M3p(1)/Mwp

101 continue

! shortcut for 4M. for 5-6M nu terms need to be recalculated for each iteration
! nutermk1 = gamma(nu1+k), where k = moment order
m31 = dn1**3*nuterm31*gincf_31
m32 = dn2**3*nuterm32*gincf_32
! mw1 = dn1**imomw*nutermw1*gincf_w1
! mw2 = dn2**imomw*nutermw2*gincf_w2
mcs_1 = dn1**imom_consv*nuterm_consv1*gincf_csv1
mcs_2 = dn2**imom_consv*nuterm_consv2*gincf_csv2

mcsm31 = mcs_1/m31
mcsm32 = mcs_2/m32
if (dn1 <= 0. .or. m31>infinity .or. mcs_1>infinity) then 
   m1frac = 0.
   ! dn1 = dnc_def
   m3 = m32
   ! goto 102
   ! go to 101
elseif (dn2 <= 0. .or. m32>infinity .or. mcs_2>infinity) then
   m1frac = 1.
   ! dn2 = dnr_def
   m3 = m31
   ! goto 102
   ! go to 101
else
   m1frac=(Mconsv_M3-mcsm32)/(mcsm31-mcsm32)
   ! print*, 'm1frac', m1frac
   ! m32frac = (mcsm31-Mconsv_M3)/(mcsm31-mcsm32)
   if (m1frac < 0.) then
      m1frac = 0.
      m3 = m32
      ! dn1 = dnc_def
   elseif (m1frac > 1.) then
      m1frac = 1.
      m3 = m31
      ! dn2 = dnr_def
   ! elseif (m1frac < 0.01 .and. m1frac > 0.) then
   !    m3 = m32/(1-m1frac)
      ! print*, 'm3 from m32', m32/(1-m1frac)
      ! print*, 'm3 from m31', m31/m1frac
   else
      m3f1 = m31/m1frac
      m3f2 = m32/(1-m1frac)
      m3 = min(m3f1, m3f2)
      ! m3 = m1frac*m31 + (1-m1frac)*m32
      ! print*, m3, m3f1, m3f2
   endif
endif

m3f1 = m3*m1frac
m3f2 = m3*(1-m1frac)
n01 = m3f1/(dn1**(nu1+3)*gamma(nu1+3)*gincf_31)
n02 = m3f2/(dn2**(nu2+3)*gamma(nu2+3)*gincf_32)

m01 = n01*(dn1**nu1*gamma(nu1))*gincf_01
mw1 = n01*(dn1**(nu1+imomw)*gamma(nu1+imomw))*gincf_w1
m02 = n02*(dn2**nu2*gamma(nu2))*gincf_02
mw2 = n02*(dn2**(nu2+imomw)*gamma(nu2+imomw))*gincf_w2
m0 = m01+m02
mw = mw1+mw2


return
end subroutine get_mom_from_param

! --------------- calc_mom -------------------
subroutine calc_mom(mom_val, mass_dist, mom_order)
! calculate moment value (mom_val) from mass_dist given an order of moment (mom_order)

integer, intent(in) :: mom_order
double precision, intent(in) :: mass_dist(nkr)
double precision, intent(out) :: mom_val

if (bintype .eq. 'tau') then
   mom_val = sum(mass_dist(1:nkr)/xkgmean*pdiams(mom_order+1,1:nkr))*col
elseif (bintype .eq. 'sbm') then
   mom_val = sum(mass_dist(1:nkr)/xl*pdiams(mom_order+1,1:nkr))*col*1000
endif

end subroutine calc_mom

! --------------- calc_pmoms -------------------
subroutine calc_pmoms(momc, momr, distm, distn)
use, intrinsic :: ieee_arithmetic

real, dimension(nkr), intent(in) :: distm, distn
real, dimension(n_mom_diag), intent(out) :: momc, momr
integer :: ibin, imom
real, dimension(nkr) :: diag_m, diag_D

diag_m = distm/distn
diag_D = (diag_m*ipio6rw)**(1./3.)

do ibin = 1,nkr
  if ( (diag_D(ibin) .ne. diag_D(ibin)) .or. (.not. ieee_is_finite(diag_D(ibin))) ) &
    diag_D(ibin) = diams(ibin)
enddo

do imom=0,n_mom_diag-1
  momc(imom+1) = sum( distn(1:split_bins-1)*diag_D(1:split_bins-1)**(real(imom)) )
  momr(imom+1) = sum( distn(split_bins:nkr)*diag_D(split_bins:nkr)**(real(imom)) )
enddo

end subroutine calc_pmoms
! --------------- fcn_2m -----------------
subroutine fcn_2m(n, dn, err_vec, iflag)
integer :: n, iflag
double precision :: dn(n), err_vec(n)
double precision :: nu, md_part(1:nkr), mx, m3
 
if (dn(1)<0. .or. dn(1)>1e4) then
  err_vec=100
else
  if (icat == 1) nu = h_shape(1)
  if (icat == 2) nu = h_shape(2)
  md_part(:) = 0.

  if (l_truncated) then
    call incgamma_2cat(nu, dn(1), md_part)
    call calc_mom(m3, md_part, 3)
    call calc_mom(mx, md_part, imomx)
  else
    call get_mom_2m(nu,dn(1),m3,mx)
  endif

  if (m3>0.) then
    err_vec(1) = log10( (Mxp(icat)/M3p(icat))*(m3/mx) )
  else
    ! print*, 'nonpositive m3', m3
    ! print*, 'Mxp(icat), M3p(icat)', Mxp(icat), M3p(icat)
    ! print*,  'scl, ecl', scl, ecl
    ! print*, 'dn,nu', dn,nu
    ! print*, 'imomx', imomx
    ! print*, 'md_part', md_part
    ! print*, 'icat', icat
    err_vec(1)=100.
    ! stop
  endif
endif

return
end subroutine fcn_2m

! --------------- fcn_4m -----------------
subroutine fcn_4m(n, dn, err_vec, iflag)
integer :: n, iflag
double precision :: dn(n), err_vec(n)
double precision :: nu1, dn1, nu2, dn2, Mconsv_M3, md(nkr), mw, m3, m0

nu1 = h_shape(1)
dn1 = dn(1)
nu2 = h_shape(2)
dn2 = dn(2)
Mconsv_M3 = Mxp(1)/M3p(1)
md(:) = 0.

if (l_truncated) then
  call incgamma_1cat(dn1, nu1, dn2, nu2, Mconsv_M3, md)
  call calc_mom(m3, md, 3)
  call calc_mom(m0, md, 0)
  call calc_mom(mw, md, imomw)
else
  call get_mom_from_param(dn1, dn2, nu1, nu2, Mconsv_M3, m3, m0, mw)
  ! ahu 2024-06-30
  dn = (/dn1,dn2/)
endif

err_vec(1) = log10( (M0p/M3p(1))/(m0/m3) )*relax0
err_vec(2) = log10( (Mwp/M3p(1))/(mw/m3) )

! print*, 'err_vec', err_vec
! print*, 'relaxw', relaxw

if (any(md<0)) then
  err_vec = err_vec*1000
endif

end subroutine fcn_4m

! --------------- fcn_5m -----------------
subroutine fcn_5m(n, prms, err_vec, iflag)
integer :: n, iflag
double precision :: prms(n), err_vec(n)
double precision :: nu1, dn1, nu2, dn2, Mconsv_M3, md(nkr), mw, m3, m0, mx, my

dn1 = prms(1)
dn2 = prms(2)
nu1 = prms(3)
Mconsv_M3 = Myp/M3p(1)
md(:) = 0.
call incgamma_1cat(dn1, nu1, dn2, nu2, Mconsv_M3, md)

call calc_mom(m3, md, 3)
call calc_mom(m0, md, 0)
call calc_mom(mw, md, imomw)
call calc_mom(my, md, imomy)

err_vec(1) = log10( (M0p/M3p(1))/(m0/m3) )*relax0
err_vec(2) = log10( (Mwp/M3p(1))/(mw/m3) )
err_vec(3) = log10( (Myp/M3p(1))/(mx/m3) )*relaxx

if (any(md<0)) then
  err_vec = err_vec*1000
endif

end subroutine fcn_5m

! --------------- fcn_6m -----------------
subroutine fcn_6m(n, prms, err_vec, iflag)
integer :: n, iflag
double precision :: prms(n), err_vec(n)
double precision :: nu1, dn1, nu2, dn2, Mconsv_M3, md(nkr), mw, m3, m0, mz, my

dn1 = prms(1)
dn2 = prms(2)
nu1 = prms(3)
nu2 = prms(4)
Mconsv_M3 = Mzp/M3p(1)
md(:) = 0.
call incgamma_1cat(dn1, nu1, dn2, nu2, Mconsv_M3, md)

call calc_mom(m3, md, 3)
call calc_mom(m0, md, 0)
call calc_mom(mw, md, imomw)
call calc_mom(my, md, imomy)
call calc_mom(mz, md, imomz)

err_vec(1) = log10( (M0p/M3p(1))/(m0/m3) )*relax0
err_vec(2) = log10( (Mwp/M3p(1))/(mw/m3) )
err_vec(3) = log10( (Myp/M3p(1))/(my/m3) )*relaxx
err_vec(4) = log10( (Mzp/M3p(1))/(mz/m3) )*relaxy

if (any(md<0)) then
  err_vec = err_vec*1000
endif

end subroutine fcn_6m

! --------------- calcerr_1cat ------------------
subroutine calcerr_1cat(guess_eval, Mconsv_M3)
type(guess_param), intent(inout) :: guess_eval
double precision, intent(in) :: Mconsv_M3
double precision :: dn1, nu1, dn2, nu2, md(nkr)
double precision :: m3, m0, mw, mx, my, mz, ratio

md(:) = 0.
guess_eval%error(:) = 0.
guess_eval%sqerr = 0.
guess_eval%l1_mfrac = 0.
dn1 = guess_eval%l1prm(1)
nu1 = guess_eval%l1prm(2)
dn2 = guess_eval%l2prm(1)
nu2 = guess_eval%l2prm(2)

if (l_truncated) then
  call incgamma_1cat(dn1, nu1, dn2, nu2, Mconsv_M3, md)
  call calc_mom(m3, md, 3)
  call calc_mom(m0, md, 0)
  call calc_mom(mw, md, imomw)
else
  call get_mom_from_param(dn1, dn2, nu1, nu2, Mconsv_M3, m3, m0, mw)
endif

! if (l_init_test) then
!   print*, 'm3,m0,mw', m3, m0, mw
!   print*, 'm3/m0 vs. M3p/M0p',m3/m0,M3p(1)/M0p
!   print*, 'm3/mw vs. M3p/Mwp',m3/mw,M3p(1)/Mwp
! endif

ratio = (M0p/M3p(1)) * (m3/m0)
if (ratio>0) then
  guess_eval%error(1) = log10(ratio)
else
  guess_eval%error(1) = log10(-ratio)*1000
endif

! print*, 'm0 ratio', ratio

ratio = (Mwp/M3p(1)) * (m3/mw)
if (ratio>0) then
  guess_eval%error(2) = log10(ratio)
else
  guess_eval%error(2) = log10(-ratio)*1000
endif

! print*, 'mw ratio', ratio

guess_eval%error(3) = 0.
guess_eval%error(4) = 0.

if ( npm >= 5 ) then
  call calc_mom(mx, md, imomx)
  ratio = (Mxp(1)/M3p(1)) * (m3/mx)
  if (ratio>0) then
    guess_eval%error(3) = log10(ratio)
  else
    guess_eval%error(3) = log10(-ratio)*1000
  endif
endif

if ( npm >= 6 ) then
  call calc_mom(my, md, imomy)
  ratio = (Myp/M3p(1)) * (m3/my)
  if (ratio>0) then
    guess_eval%error(4) = log10(ratio)
  else
    guess_eval%error(4) = log10(-ratio)*1000
  endif
endif

if (any(md<0)) guess_eval%error=guess_eval%error*1000
guess_eval%sqerr = sqrt(sum(guess_eval%error**2))
guess_eval%l1_mfrac = m1frac

end subroutine calcerr_1cat

! --------------- update_param_guess -----------------
subroutine update_param_guess(guess_best, guess_candidate, mark_improve)
! replace guess_best with guess_candidate if its sqerr (squared error)
! is smaller
type(guess_param), intent(inout) :: guess_best
type(guess_param), intent(in) :: guess_candidate
logical, optional :: mark_improve

if ( .not. present(mark_improve) ) mark_improve = .false.
if ( guess_candidate%sqerr < guess_best%sqerr ) then
  guess_best = guess_candidate
  if ( mark_improve ) l_improved = .true.
endif

end subroutine update_param_guess
! --------------- minpack interface -----------------
subroutine find_params(fcn, guess_best, tol, exit_signal)
external fcn
type(guess_param), intent(inout) :: guess_best
double precision, intent(in) :: tol
integer, intent(out) :: exit_signal
integer :: n
double precision :: guess(4), err_vec(4) ! maximum possible for 6M U-AMP
integer, parameter :: lwa=50 ! should be large enough, doesn't have to be precise
double precision, dimension(lwa) :: wa

guess_best%error(:) = 0.
guess_best%sqerr = 0.

if (npm == 2) n=1
if (npm == 3 .or. npm == 4) n=2
if (npm == 5) n=3
if (npm == 6) n=4

if ( npm < 4 ) then
  ! only get dn if 2m, get both dn and nu if 3m
  if (icat == 1) guess(1:n) = guess_best%l1prm(1:n)
  if (icat == 2) guess(1:n) = guess_best%l2prm(1:n)
elseif ( npm >= 4 ) then
  guess(1) = guess_best%l1prm(1)
  guess(2) = guess_best%l2prm(1)
endif

if ( npm >= 5 ) then
  guess(3) = guess_best%l1prm(2)
endif
if ( npm >= 6 ) then
  guess(4) = guess_best%l2prm(2)
endif

! print*, 'prior guesses', guess(1:n)
call hybrd1(fcn, n, guess(1:n), err_vec(1:n), tol, exit_signal, wa, lwa)
! print*, 'post guesses', guess(1:n)

if ( npm < 4 ) then
  ! only get dn if 2m, get both dn and nu if 3m
  if (icat == 1) guess_best%l1prm(1:n) = guess(1:n)
  if (icat == 2) guess_best%l2prm(1:n) = guess(1:n)
elseif ( npm >= 4 ) then
  guess_best%l1prm(1) = guess(1)
  guess_best%l2prm(1) = guess(2)
endif

if ( npm >= 5 ) then
  guess_best%l1prm(2) = guess(3)
endif
if ( npm >= 6 ) then
  guess_best%l2prm(2) = guess(4)
endif

! print*, 'posterior guesses', guess

guess_best%error(1:n) = err_vec(1:n)
guess_best%sqerr = sqrt(sum(err_vec(1:n)**2))
guess_best%l1_mfrac = m1frac
! print*, 'err_vec(1:n)', err_vec(1:n)
! print*, 'sqerr', sqrt(sum(err_vec(1:n)**2))

end subroutine find_params
! --------------- end of minpack interface -----------------

end module
