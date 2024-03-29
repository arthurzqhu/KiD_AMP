! -------- 2 category AMP -----------
! fcn_2p: {{{
subroutine fcn_2p(n,x,fvec,iflag)
use micro_prm, only:ihyd,nubounds,dnbounds,nkr,diams,M3p,Mxp,Myp,momx,momy,skr,ekr,relax,nkr
use parameters, only: max_nbins
implicit none
integer n,iflag
double precision x(n),fvec(n)

double precision rx, nu, dn, m3, mx, my
double precision, dimension(max_nbins):: md

!The purpose of this routine is to calculate the number/mass distribution
!of drops given parameters of the gamma distribution
!f(D)=n0*(D/dn)^(nu+3-1)*exp(-D/dn)
!The +3 is added because this is really a mass distribution, but we want
!the shape parameter for a number distribution.
!n0=(total mass)/gamma(nu+3)
!This function is used in conjunction with the equation solving routines.
!fvec characterizes how close the given parameters are to a solution.

!Some bounds on the parameters
if (x(2)<=0. .or. x(2)>dnbounds(2,ihyd) &
   .or. x(1)<nubounds(1,ihyd)*2. .or. x(1)>nubounds(2,ihyd)*2.) then
!If the parameters exceed the bounds, assign a high error
  fvec(1)=100.
  fvec(2)=100.
else
  rx=1.e-10 !This is an arbitrary value of rx
  nu=x(1) !Shape parameter
  dn=x(2) !Scale diameter multiplied by 1e6 to make its magnitude
                !approximately the same as nu

  md=0. !Initialize the mass distribution
  call incgamma_norm(rx,nu,dn,skr,ekr,md)
!  call incjohnsonsb(nu,dn*1.e6,skr,ekr,md)

  !Calculate the moments - 0th, 3rd, and xth
  call calcmom(m3,md,3)
  call calcmom(mx,md,momx)
  call calcmom(my,md,momy)

  !Calculate the errors in the moments
  if (m3>0.) then
    fvec(1)=log10((Mxp/M3p)*(m3/mx))
    fvec(2)=log10((Myp/M3p)*(m3/my))*relax**2.
  else
    fvec=100.
  endif
endif
return
End subroutine fcn_2p

! }}}
! fcn_1p: {{{
subroutine fcn_1p(n,x,fvec,iflag)
use micro_prm, only:diams,Mxp,M3p,skr,ekr,nug,momx,nkr
use parameters, only: max_nbins
implicit none
integer n,iflag
double precision x(n),fvec(n)

double precision rx, nu, dn, m3, mx
double precision, dimension(max_nbins):: md

!The purpose of this routine is to calculate the number/mass distribution
!of drops given parameters of the gamma distribution
!f(D)=n0*(D/dn)^(nu+3-1)*exp(-D/dn)
!The +3 is added because this is really a mass distribution, but we want
!the shape parameter for a number distribution.
!n0=(total mass)/gamma(nu+3)
!This function is used in conjunction with the equation solving routines.
!fvec characterizes how close the given parameters are to a solution.

!Some bounds on the parameters
if (x(1)<=0. .or. x(1)>1.e10) then
  !If the parameters exceed the bounds, assign a high error
  fvec(1)=100.
else
  rx=1.e-10 !This is an arbitrary value of rx
  nu=nug !nug is the global value of nu
  dn=x(1) !x(1) is the scale diameter multiplied by 1e6 to make its magnitude
                !approximately the same as nu

  md=0. !Initialize the distribution
  call incgamma_norm(rx,nu,dn,skr,ekr,md)

  !Calculate the moments - xth and 3rd
  call calcmom(m3,md,3)!sum(md(1:nkr)/xl(1:nkr)*diams(1:nkr)**3)!*col*1000.
  call calcmom(mx,md,momx)!sum(md(1:nkr)/xl(1:nkr)*diams(1:nkr)**momx)!*col*1000.

  !Calculate the errors in the moments
  if (m3>0.) then
    fvec(1)=log10((Mxp/M3p)*(m3/mx))
  else
    fvec(1)=100.
  endif
endif
return
End subroutine fcn_1p

! }}} 
! calcerr: {{{
subroutine calcerr(rx,nu,dn,error)
use micro_prm, only:diams,M3p,Mxp,Myp,col,nkr,skr,ekr,relax,momx,momy
use parameters, only: max_nbins
implicit none
double precision :: ratio,m3,mx,my,md(max_nbins),error(2),rx,nu,dn

  md = 0.
  call incgamma_norm(rx,nu,dn,skr,ekr,md)
!  call incjohnsonsb(nu,dn*1.e6,skr,ekr,md)

  !Calculate the moments - 3rd, xth, and yth
  call calcmom(m3,md,3) 
  call calcmom(mx,md,momx) 
  call calcmom(my,md,momy)

  !Calculate the errors in the moments
  ratio = (Mxp/mx)*(m3/M3p)
  error(1)=log10(ratio)

  ratio = (Myp/my)*(m3/M3p)
  error(2)=log10(ratio)*relax**2.

End subroutine calcerr


! }}}
! calcmom: {{{
! can now be completely replaced by dist2momv() but am too lazy to do that - ahu
subroutine calcmom(mom,md,momnum)

use micro_prm, only: diams,col,nkr,binmass
use namelists, only: bintype
use parameters, only: max_nbins
use module_hujisbm, only:xl
implicit none

double precision :: mom,md(max_nbins)
integer :: momnum

!x1000 for sbm to convert xl from g to kg
!tau binmass already in kg
if (bintype .eq. 'sbm') then
    mom=sum(md(1:nkr)/xl(1:nkr)*diams(1:nkr)**momnum)*col*1000.
elseif (bintype .eq. 'tau') then
    mom=sum(md(1:nkr)/binmass(1:nkr)*diams(1:nkr)**momnum)*col
endif
end subroutine calcmom

! }}}
! calc_nummom: {{{
! can now be completely replaced by dist2momv() but am too lazy to do that - ahu
subroutine calc_nummom(mom,md,momnum)

use micro_prm, only: diams,col,nkr,binmass
use namelists, only: bintype
use parameters, only: max_nbins
use module_hujisbm, only:xl
implicit none

double precision :: mom,md(max_nbins)
integer :: momnum

!x1000 for sbm to convert xl from g to kg
!tau binmass already in kg
if (bintype .eq. 'sbm') then
    mom=sum(md(1:nkr)*diams(1:nkr)**momnum)*col
elseif (bintype .eq. 'tau') then
    mom=sum(md(1:nkr)*diams(1:nkr)**momnum)*col
endif
end subroutine calc_nummom

!}}}
! incgamma: {{{
subroutine incgamma(rx,nu,dn,ia,iz,md)
use micro_prm, only:diams,nkr
use parameters, only: max_nbins 
use namelists, only: bintype
use mphys_tau_bin_declare, only:dgmean
implicit none

double precision :: rx,nu,dn,n0,expterm,md(max_nbins)
integer :: ia,iz,kr

  n0=rx
  do kr=ia,iz !Loop over bins
    if (bintype .eq. 'sbm') then
        expterm=exp(-1.*diams(kr)/dn)
        md(kr)=n0*expterm*(diams(kr)/dn)**(nu+3.)
    elseif (bintype .eq. 'tau') then
        !expterm=exp(-1.*diams(kr)/dn)
        !md(kr)=n0*expterm*(diams(kr)/dn)**(nu+3.)
        md(kr)=exp(log(rx)-diams(kr)/dn+(nu+3.)*log(diams(kr)/dn))
    endif
    !Check for NaNs
!    if (md(kr).ne.md(kr) .or. md(kr)*0.0 .ne.0.0 &
!       .or. md(kr)/md(kr) .ne. 1.0) md(kr)=0.
  enddo

return
End subroutine incgamma

! }}}
! incgamma_norm: {{{
subroutine incgamma_norm(rx,nu,dn,ia,iz,md)
use micro_prm, only:diams,nkr,ihyd
use parameters, only: max_nbins
use namelists, only: bintype
use mphys_tau_bin_declare, only:dgmean
implicit none

double precision :: rx,nu,dn,n0,expterm,md(max_nbins),term1(max_nbins),term2
integer :: ia,iz,kr

!Adele - substantially different than before. Making use of exp and log 
!properties to ensure that the value of md(ia)=1

  n0=rx
  do kr=ia,iz !Loop over bins
     term1(kr)=-diams(kr)/dn+(nu+3.)*log(diams(kr)/dn)
  enddo
  term2 = maxval(term1(ia:iz))*-1.

  do kr=ia,iz
     md(kr)=exp(term1(kr)+term2)
  enddo

return
End subroutine incgamma_norm

! }}}
! incjohnsonsb: {{{
subroutine incjohnsonsb(s1,s2,ia,iz,md)
use micro_prm, only:diams,nkr
use parameters, only: max_nbins
implicit none

double precision :: s1,s2,minx,ranx,z,expterm,md(max_nbins)
integer :: ia,iz,kr

minx=0.
ranx=diams(nkr)*2. ! might need to change it based on sbm or tau -ahu

  do kr=ia,iz !Loop over bins
    z=(diams(kr)-minx)/ranx
    expterm=exp(-0.5*(s1+s2*log(z/(1.-z)))**2)
    md(kr)=diams(kr)/(z*(1-z))*expterm
    !Check for NaNs
    if (md(kr).ne.md(kr) .or. md(kr)*0.0 .ne.0.0 &
       .or. md(kr)/md(kr) .ne. 1.0) md(kr)=0.
  enddo

return
End subroutine incjohnsonsb

! }}}
! calcdist: {{{
subroutine calcdist(x,md)
use micro_prm, only:nkr,diams,M3p,col,skr,ekr,rxfinal,ihyd,dnbounds
use parameters, only: max_nbins
use diagnostics, only: i_dgtime
implicit none
double precision x(2)
integer n

double precision rx, nu, dn, m3
double precision, dimension(max_nbins):: md


!Use this function to force correct 3rd moment
!Very similar to fcn_2p
  nu=x(1)
  dn=x(2)
  rx = 1
  n=0
 do
  n=n+1
  md=0.

  call incgamma_norm(rx,nu,dn,skr,ekr,md)
  call calcmom(m3,md,3)

  if ((m3==0. .or. dn==0.) .and. M3p<1.e-20) then
    rxfinal=0.
    md = 0.
    exit
  elseif (m3==0. .or. dn==0.) then
    print*,'cannot find rxfinal and mass = ',M3p*3.1415/6.*1000
    print*,m3,rx,nu,dn,ihyd,dnbounds(:,ihyd)
    stop
  endif
  if (M3p/m3 > 1.e10) then !Might need similar condition for M3p/m3< 1.e-10
    rx=rx*1000.
  else
    rxfinal = M3p/m3*rx
    md = md*rxfinal/rx

    exit
  endif
  if (n>1000) then
    print*,'calcdist stuck',m3,M3p
    stop
  endif
 enddo

return
End subroutine calcdist
! }}}
! searchparamsG: {{{
Subroutine searchparamsG(guess,ihyd,md,flag)

use micro_prm, only:nkr,npm
use parameters, only:flag_count, max_nbins

implicit none
real(8) :: guess(2)
real(8) :: md(max_nbins)
integer :: ihyd
real(8),dimension(flag_count) :: flag

if (npm==3) then !if rain is 2M, then cloud is also 2M
   CALL searchparams3M(guess,ihyd,md,flag)
else
   CALL searchparams2M(guess,md,flag(1))
endif

End subroutine searchparamsG

! }}}
! searchparams3M: {{{
Subroutine searchparams3M(guess,ihyd,md,flag)

use micro_prm, only:relax,Mp,M3p,Mxp,Myp,momx,nkr, &
               nutab,dntab,minmaxmx,mintab,maxtab,ntab, &
               cloud_mr_th, rain_mr_th,minmaxmy
use parameters, only:flag_count, max_nbins
use namelists, only:ovc_factor
use diagnostics, only:i_dgtime
use, intrinsic :: ieee_arithmetic, only: IEEE_Value, IEEE_QUIET_NAN
use, intrinsic :: iso_fortran_env, only: real32

implicit none
real(8) :: guess(2),oguess(2),vals(2),ovals(2),tol
real(8) :: MxM3,MyM3,irm,iry,wgtm,wgty1,wgty2
real(8) :: minmy1,maxmy1,minmy2,maxmy2,sqval,osqval,minsqval,md(max_nbins),min12,max12
real(8) :: minMx3,maxMx3,minMy3,maxMy3,oMxM3,oMyM3,rdummy,rdummy2
integer :: i,info,n,ihyd,im1,im2,iy1a,iy2a,ix1b,ix2b,dummy,dummy2!,flag
real(8), dimension(flag_count) :: flag
integer, parameter :: lwa=max_nbins
real(8),dimension(lwa) :: wa
external :: fcn_2p
real(real32) :: nan

nan = IEEE_VALUE(nan, IEEE_QUIET_NAN)

flag(:)=0.
!This subroutine drives the search for PDF parameters that satisfy the
!moments.

M3p=Mp(1);Mxp=Mp(2);Myp=Mp(3)

relax=1. !Initial value of relax. Do not change.
tol=1.d-6 !Tolerance of minimization routines
n=2 !Number of equations to solve. Do not change.

info=0 !Message sent back from minimization routines.

oguess=guess !oguess is the current best guess
ovals=0.
CALL calcerr(100.d0,oguess(1),oguess(2),ovals)
osqval = sqrt(sum(ovals**2))

!First try - previous values
if (guess(2).ne.0) then
  CALL hybrd1(fcn_2p,n,guess,vals,tol,info,wa,lwa)
  sqval = sqrt(sum(vals**2))
  if (sqval .lt. osqval) then
     ovals = vals
     osqval = sqval
     oguess = guess
  endif
endif

!If first try was wildly off, or no 1st guess exists - use look up table
if (guess(2).eq.0 .or. sqval>tol) then
  !Ratio of xth moment to 3rd moment.
  MxM3 = Mxp/M3p
  !Ratio of yth moment to 3rd moment.
  MyM3 = Myp/M3p

  oMxM3 = MxM3
  oMyM3 = MyM3

  !Need to see if these ratios are in the solution space
  !If not, adjust Mx and/or My until they are in the solution space

  !Lower and upper limits of Mx/M3 can be calculated,
  !but I threw those limits out as being too hard to deal with
  !Actual practical limits are stored in the minmaxmx table
  minMx3 = minmaxmx(1,momx+1,ihyd)
  maxMx3 = minmaxmx(2,momx+1,ihyd)

  ! test skipping this part...

  !Check now to see if MxM3 is out of allowable range and adjust Mx
  if (MxM3 < minMx3 .or. MxM3 > maxMx3) then
     minMy3 = minmaxmy(1,ihyd)
     maxMy3 = minmaxmy(2,ihyd)
     if (MyM3>=minMy3 .and. MyM3<=maxMy3) then !Mx/M3 is not in range, My/M3 is in range
        dummy = minloc(abs(mintab(:,ihyd)-MyM3),DIM=1,MASK=(mintab(:,ihyd)-MyM3)<0)
        rdummy = minMx3*(maxMx3/minMx3)**((real(dummy)-1.)/(real(ntab)-1.)) 
        dummy2 = minloc(abs(maxtab(:,ihyd)-MyM3),DIM=1,MASK=(maxtab(:,ihyd)-MyM3)>0)
        rdummy2 = minMx3*(maxMx3/minMx3)**((real(dummy2)-1.)/(real(ntab)-1.)) 
        minMx3 = min(rdummy,rdummy2)
        maxMx3 = max(rdummy,rdummy2)
     endif
     if (MxM3< minMx3) then
         MxM3 = min(minMx3*(1.+ovc_factor),maxMx3)
     elseif (MxM3 > maxMx3) then
         MxM3 = max(maxMx3*(1.-ovc_factor),minMx3)
     endif
     Mxp = MxM3 * M3p
     Mp(2)=Mxp
     !Reset to global min and max values of Mx/M3
     minMx3 = minmaxmx(1,momx+1,ihyd)
     maxMx3 = minmaxmx(2,momx+1,ihyd)
  endif

  !Find where we are on a log10 scale bewteen those two points
  !in terms of the number of points in our look up tables
  irm = log10(MxM3/minMx3)/log10(maxMx3/minMx3)*(real(ntab)-1.)+1.
  irm = max(1.,min(real(ntab)-1.,irm))
  im1 = int(floor(irm))
  wgtm = irm-floor(irm)
  im2=im1+1

  !Use im1 and im2 to find min and max value of My/M3
  !ihyd1=cloud, ihyd2=rain

  minmy1 = mintab(im1,ihyd)
  maxmy1 = maxtab(im1,ihyd)
  minmy2 = mintab(im2,ihyd)
  maxmy2 = maxtab(im2,ihyd)

  !Check now to see if MxM3 is out of allowable range and adjust Mx
  min12=minmy2*wgtm+minmy1*(1-wgtm)
  max12=maxmy2*wgtm+maxmy1*(1-wgtm)

  if (MyM3 < min12) then
    MyM3 = min(min12*(1.+ovc_factor),max12)
    Myp = MyM3 * M3p
  elseif (MyM3 > max12) then
    MyM3 = max(max12*(1.-ovc_factor),min12)
    Myp = MyM3 * M3p
  endif

  Mp(3)=Myp

  !Again, find where we are on the log scale between the upper and lower
  !limits. Do this twice, once for each of the nearest min and max pairs
  iry=log10(MyM3/minmy1)/log10(maxmy1/minmy1)*(real(ntab)-1.)+1.
  iry = max(1.,min(real(ntab),iry))
  iy1a = floor(iry)
  wgty1=iry-iy1a
  iy2a=min(ntab,iy1a+1)

  iry=log10(MyM3/minmy2)/log10(maxmy2/minmy2)*(real(ntab)-1.)+1.
  iry = max(1.,min(real(ntab),iry))
  ix1b = floor(iry)
  wgty2=iry-ix1b
  ix2b=min(ntab,ix1b+1)

  flag(2) = sqrt( log10(MxM3/oMxM3)**2 + log10(MyM3/oMyM3)**2 )
  !Find best-guess parameters in the look up tables
  guess(1) = (1.-wgtm)*((1.-wgty1)*nutab(iy1a,im1,ihyd)+wgty1*nutab(iy2a,im1,ihyd)) + &
             wgtm*((1.-wgty2)*nutab(ix1b,im2,ihyd)+wgty2*nutab(ix2b,im2,ihyd))
  guess(2) = (1.-wgtm)*((1.-wgty1)*dntab(iy1a,im1,ihyd)+wgty1*dntab(iy2a,im1,ihyd)) + &
             wgtm*((1.-wgty2)*dntab(ix1b,im2,ihyd)+wgty2*dntab(ix2b,im2,ihyd))

  !Use the best-guess. See what happens
  CALL hybrd1(fcn_2p,n,guess,vals,tol,info,wa,lwa)

  sqval = sqrt(sum(vals**2))
  if (sqval .lt. osqval) then
     ovals = vals
     osqval = sqval
     oguess = guess
  endif
endif

!If not close enough, jiggle the values and relax tolerance on Mx
!Start with current best guess
guess = oguess
vals = ovals
sqval = osqval
if (sqval>tol) then
  i=0
  minsqval=sqrt(sum(ovals**2))
  do
    if(abs(vals(1))>tol) guess = guess*(1.01)**(-1*(i+1))
    CALL hybrd1(fcn_2p,n,guess,vals,tol,info,wa,lwa)
    sqval = sqrt(sum(vals**2))

    !Keep track of best guess. Don't worry about error on My
    !for determining best guess
    if (sqval<abs(minsqval)) then
      minsqval=sqval
      ovals = vals
      osqval=sqval
      oguess=guess
    endif

    !If fitting was successful, or Mx/M3 error low enough
    !or we need to give up, or if increase in error such
    !that it's not even worth our time to keep trying
    if (info.eq.1 .or. abs(vals(1)).le.tol .or. i>4 &
         .or. abs(vals(1))>100.*abs(ovals(1))) then
       !print*,'fail'
       exit
     endif

    relax = relax*0.1
    i=i+1
  enddo

  guess=oguess !We're done. Set guess to best guess

endif

!Set flag to 1 if fitting didn't work as well as we wished
if (abs(vals(1))>tol .or. abs(vals(2))>tol) flag(1) = 1
if (abs(vals(1))>tol) flag(3)=1
if (abs(vals(2))>tol) flag(4)=1
!if (abs(vals(1))>tol .or. abs(vals(2))>tol) flag(1)=1

!Set flag to -1 or nan if cloud or rain mass didn't reach the threshold
if ((ihyd==1 .and. M3p<cloud_mr_th) .or. (ihyd==2 .and. M3p<rain_mr_th)) then
   flag(1) = -1
   flag(2:flag_count) = nan
end if

if (guess(1) .ne. guess(1) .or. guess(2) .ne. guess(2)) then
    print*,'ihyd=',ihyd
    print*, 'guess=',guess
endif
!Force third moment to have no error and calculate final distribution

CALL calcdist(guess,md)

return

End subroutine searchparams3M

! }}}
! searchparams2M: {{{
Subroutine searchparams2M(guess,md,flag)

use micro_prm, only:Mp,M3p,Mxp,ntab,skr,ekr,nkr,momx,nug,diams,total_m3_th
use parameters, only: max_nbins
use global_fun
implicit none
real(8) :: guess(2),oguess(2),vals(1),tol,guessin(1)
real(8) :: MxM3,minmxm3,maxmxm3
real(8) :: md(max_nbins)
integer :: info,n!,flag
integer, parameter :: lwa=max_nbins
real(8),dimension(lwa) :: wa
real(8) :: flag, M3t, M0t
external :: fcn_1p

flag = 0
!This subroutine drives the search for PDF parameters that satisfy the
!moments.
M3p=Mp(1)
Mxp=Mp(2)

if (M3p <= total_m3_th) then
   md = 0.
   return
endif


tol=1.d-8 !Tolerance of minimization routines
n=1 !Number of equations to solve. Do not change.
info=0 !Message sent back from minimization routines.

!Find min and max values of MxM3
minmxm3=min(diams(skr)**(momx-3.),diams(ekr)**(momx-3.))
maxmxm3=max(diams(skr)**(momx-3.),diams(ekr)**(momx-3.))

!Generic 2nd predicted moment
MxM3 = Mxp/M3p
if (MxM3<=minmxm3) then
   MxM3=minmxm3*1.0001
   Mxp=MxM3*M3p
   flag=1
elseif (MxM3>=maxmxm3) then
   MxM3=maxmxm3*0.9999
   Mxp=MxM3*M3p
   flag=1
endif

oguess=guess
nug=guess(1)!not actually a guess. A parameter.
guess(2) = (MxM3*gamma(nug+3.)/gamma(nug+momx))**(1./(momx-3.))
vals = 2.
guessin=guess(2)
CALL hybrd1(fcn_1p,n,guessin,vals,tol,info,wa,lwa)
! print*, 'after hybrd1', guessin

if (vals(1)==1.) then
  guess=oguess
else
  guess(2)=guessin(1)
endif

!Set flag to 1 if fitting didn't work as well as we wished
if (abs(vals(1))>tol) flag = 1

!For 3M fitting, we jiggle the values and relax tolerance on Mx
!No sense in doing that here. The fitting is much easier, and should
!always give the same answer.
CALL calcdist(guess,md)

! call calcmom(M3t, md, 3)
! call calcmom(M0t, md, 0)
! print*, 'guess meand from md', guess(2), get_meandiam(M3t, M0t)

return

End subroutine searchparams2M
! }}}

! -------- single category AMP -----------
! dn below is different from above (a factor of 1e6)

! fcn_1p_sc: {{{
subroutine fcn_1p_sc(n,x,fvec,iflag)
use micro_prm, only:diams,M0temp,M3temp,skr,ekr,nug,nkr
use parameters, only: max_nbins
implicit none
integer n,iflag
double precision x(n),fvec(n)

double precision rx, nu, dn, m3, m0
double precision, dimension(max_nbins):: md

!The purpose of this routine is to calculate the number/mass distribution
!of drops given parameters of the gamma distribution
!f(D)=n0*(D/dn)^(nu+3-1)*exp(-D/dn)
!The +3 is added because this is really a mass distribution, but we want
!the shape parameter for a number distribution.
!n0=(total mass)/gamma(nu+3)
!This function is used in conjunction with the equation solving routines.
!fvec characterizes how close the given parameters are to a solution.

!Some bounds on the parameters
if (x(1)<=0. .or. x(1)>1.e10) then
  !If the parameters exceed the bounds, assign a high error
  fvec(1)=100.
else
  rx=1.e-10 !This is an arbitrary value of rx
  nu=nug !nug is the global value of nu
  dn=x(1) !x(1) is the scale diameter multiplied by 1e6 to make its magnitude
                !approximately the same as nu

  md=0. !Initialize the distribution
  skr = 1
  ekr = nkr
  call incgamma_norm(rx,nu,dn,skr,ekr,md)

  !Calculate the moments - xth and 3rd
  call calcmom(m3,md,3)!sum(md(1:nkr)/xl(1:nkr)*diams(1:nkr)**3)!*col*1000.
  call calcmom(m0,md,0)!sum(md(1:nkr)/xl(1:nkr)*diams(1:nkr)**momx)!*col*1000.

  !Calculate the errors in the moments
  if (m3>0.) then
    fvec(1)=log10((M0temp/M3temp)*(m3/m0))
  else
    fvec(1)=100.
  endif
endif
return
End subroutine fcn_1p_sc

! }}} 

! guessparams2M :{{{
subroutine guessparams2M(Ma, Mb, moma, momb, guess_post, tol, vals)
! ma = m3, mb = mx

use micro_prm, only: skr, ekr, diams, nug
use parameters, only: max_nbins
implicit none
real(8), intent(in) :: Ma, Mb, tol
integer, intent(in) :: moma, momb
real(8), intent(out) :: guess_post(2)
real(8), intent(out) :: vals
real(8) :: MbMa, guess_int, minmbma, maxmbma, guessin
integer, parameter :: lwa=max_nbins
real(8),dimension(lwa) :: wa
integer :: info, n
external :: fcn_1p_sc

n = 1
info = 0
MbMa = Mb/Ma
guessin = (MbMa*gamma(nug+dble(moma))/&
   gamma(nug+dble(momb)))**(1./dble(momb-moma))
vals = 2.
CALL hybrd1(fcn_1p_sc,n,guessin,vals,tol,info,wa,lwa)

guess_post = (/nug, guessin/)

end subroutine guessparams2M

! }}}
! is_unimodal: {{{
subroutine is_unimodal(Ma, Mb, moma, momb, guess_prior, Mz_pred, momz, l_unimodal, &
                       guess_post, print_flag_opt)

use micro_prm, only: tol_unimod_frac
use sort_mod, only: are_close
implicit none

real(8), intent(in) :: Ma, Mb, guess_prior(2), Mz_pred
integer, intent(in) :: moma, momb, momz
real(8), intent(out) :: guess_post(2)
logical, intent(out) :: l_unimodal
real(8) :: nu, dn, vals, Mz_est
logical, optional :: print_flag_opt
logical :: print_flag

if (.not. present(print_flag_opt)) then
   print_flag = .false.
else
   print_flag = print_flag_opt
endif

call guessparams2M(Ma, Mb, moma, momb, guess_prior, guess_post, 1.d-8, vals)

nu = guess_post(1)
dn = guess_post(2)*1e-6
Mz_est = Ma * dn**(momz-3) * gamma(nu+momz)/gamma(nu+moma)
if (print_flag) print*, 'est:', Mz_est, ', pred:', Mz_pred

l_unimodal = are_close(Mz_pred, Mz_est, 'ratio', tol_unimod_frac)

end subroutine is_unimodal

! }}}

! get_relative_moment_vals: {{{

Subroutine  get_relative_moment_vals(nu1,dn1,nu2,dn2,m3,m0,mw,mx,my)
use parameters, only: num_h_moments, max_nbins
use micro_prm, only: nuterm31, nuterm32, nutermz1, nutermz2, nutermw1, nutermw2, &
   M3p,Mwp,Mzp,momz,momw,dnc_def,dnr_def, nkr, diams
real(8), intent(inout) :: nu1,dn1,nu2,dn2
real(8) :: mz1m31, mz2m32, m31mw1, m32mw2, m3frac, mwfrac, &
   MzM3,M3Mw,m31,m32,mw1,mw2,mz1,mz2, infinity, m32frac, mw2frac, &
   gterm1b, gterm2b, m3f1, m3f2, mwf1, mwf2
real(8), dimension(max_nbins) :: md1(max_nbins), md2(max_nbins), &
   gterm1a(max_nbins), gterm2a(max_nbins)
real(8), intent(out) :: m3,m0,mw
real(8), intent(out), optional :: mx, my

infinity = HUGE(.01)

! print*, 'grmv dn1, dn2', dn1, dn2

MzM3 = Mzp/M3p
M3Mw = M3p/Mwp

101 continue

m31 = dn1**3*nuterm31
m32 = dn2**3*nuterm32
mw1 = dn1**momw*nutermw1
mw2 = dn2**momw*nutermw2
mz1 = dn1**momz*nutermz1
mz2 = dn2**momz*nutermz2

! do kr=1,nkr
!    gterm1a(kr)=-diams(kr)/dn1+(nu1+3)*log(diams(kr)/dn1)
!    gterm2a(kr)=-diams(kr)/dn2+(nu2+3)*log(diams(kr)/dn2)
! enddo

! gterm1b = maxval(gterm1a(1:nkr))*(-1.)
! gterm2b = maxval(gterm2a(1:nkr))*(-1.)

! do kr=1,nkr
!    md1(kr)=exp(gterm1a(kr)+gterm1b)
!    md2(kr)=exp(gterm2a(kr)+gterm2b)
! enddo

! call calcmom(m31, md1, 3)
! call calcmom(m32, md2, 3)
! call calcmom(mw1, md1, momw)
! call calcmom(mw2, md2, momw)
! call calcmom(mz1, md1, momz)
! call calcmom(mz2, md2, momz)

! print*, 'truncated m31 m32', m31, m32
! print*, 'full m31, m32', dn1*nuterm31, dn2*nuterm32
! stop

m0 = 1.

mz1m31 = mz1/m31
mz2m32 = mz2/m32
if (dn1 <= 0. .or. m31>infinity .or. mz1>infinity) then 
   m3frac = 0.
   dn1 = dnc_def
   ! m3 = m32
   go to 101
elseif (dn2 <= 0. .or. m32>infinity .or. mz2>infinity) then
   m3frac = 1.
   dn2 = dnr_def
   ! m3 = m31
   go to 101
else
   m3frac=(MzM3-mz2m32)/(mz1m31-mz2m32)
   ! m32frac = (mz1m31-MzM3)/(mz1m31-mz2m32)
   if (m3frac < 0.) then
      m3frac = 0.
      m3 = m32
   elseif (m3frac > 1.) then
      m3frac = 1.
      m3 = m31
   ! elseif (m3frac < 0.01 .and. m3frac > 0.) then
   !    m3 = m32/(1-m3frac)
      ! print*, 'm3 from m32', m32/(1-m3frac)
      ! print*, 'm3 from m31', m31/m3frac
   else
      m3f1 = m31/m3frac
      m3f2 = m32/(1-m3frac)
      m3 = min(m3f1, m3f2)
      ! print*, m3, m3f1, m3f2
   endif
endif

m31mw1 = m31/mw1
m32mw2 = m32/mw2
if (dn1<=0. .or. m31>infinity .or. mw1>infinity) then
   mwfrac = 0.
   dn1 = dnc_def
   ! mw = mw2
   go to 101
elseif (dn2<=0. .or. m32>infinity .or. mw2>infinity) then
   mwfrac = 1.
   dn2 = dnr_def
   ! mw = mw1
   go to 101
else
   mwfrac = (M3Mw-m32mw2)/(m31mw1-m32mw2)
   ! mw2frac = (m31mw1-M3Mw)/(m31mw1-m32mw2)
   if (mwfrac <= 0.) then
      mwfrac = 0.
      mw = mw2
   elseif (mwfrac >= 1.) then
      mwfrac = 1.
      mw = mw1
   ! elseif (mwfrac < 0.001 .and. mwfrac > 0.) then
   !    mw = mw2/(1-mwfrac)
      ! print*, 'mw from mw2', mw2/(1-mwfrac)
      ! print*, 'mw from mw1', mw1/mwfrac
      ! stop
   else
      ! mw = mw1/mwfrac
      mwf1 = mw1/mwfrac
      mwf2 = mw2/(1-mwfrac)
      mw = min(mwf1, mwf2)
      ! print*, mw, mwf1, mwf2
      ! stop
   endif
endif

! print*, 'grmv m31, m32, m3frac', m31, m32, m3frac
! print*, 'grmv mwfrac', mw1, mw2, mwfrac
! print*, 'm32frac', m32frac
! print*, 'mw2frac', mw2frac

if (num_h_moments(1) == 6) then
   stop '6M U-AMP not yet implemented.'
endif

End Subroutine get_relative_moment_vals 
! }}}

! fcn_6p: {{{
subroutine fcn_6p(n,x,fvec,iflag)
use micro_prm, only:ihyd,nubounds,dnbounds,nkr,diams, &
                    M0p,M3p,Mxp,Myp,Mwp,Mzp,momw,momx,momy,momz, &
                    col,skr,ekr,rxfinal,relaxw,relaxx,relaxy
use parameters, only: max_nbins
implicit none
integer n,iflag
double precision x(n),fvec(n)

double precision rx, nu1, dn1, nu2, dn2, frac, m3, m0, mw, mx, my, mz, m3frac
double precision, dimension(max_nbins):: md
integer kr

!The purpose of this routine is to calculate the number/mass distribution
!of drops given parameters of the gamma distribution
!f(D)=n0*(D/dn)^(nu+3-1)*exp(-D/dn)
!The +3 is added because this is really a mass distribution, but we want
!the shape parameter for a number distribution.
!n0=(total mass)/gamma(nu+3)
!This function is used in conjunction with the equation solving routines.
!fvec characterizes how close the given parameters are to a solution.

  rx=1.e-10 !This is an arbitrary value of rx
  nu1=x(1) !Shape parameter
  dn1=x(2) !Scale diameter multiplied by 1e6 to make its magnitude 
                !approximately the same as nu
  nu2=x(3) !Shape parameter
  dn2=x(4) !Scale diameter multiplied by 1e6 to make its magnitude 
                !approximately the same as nu
  frac = Mzp/M3p

  md=0. !Initialize the mass distribution
  call incgamma2(rx,nu1,dn1,nu2,dn2,frac,md,m3frac)

  !Calculate the moments - 0th, 3rd, wth, xth, yth, and zth
  call calcmom(m3, md, 3)
  call calcmom(m0, md, 0)
  call calcmom(mw, md, momw)
  call calcmom(mx, md, momx)
  call calcmom(my, md, momy)
 
  fvec(1)=((Mwp/M3p)*(m3/mw)-1.)*relaxw
  fvec(2)=((Mxp/M3p)*(m3/mx)-1.)*relaxx
  fvec(3)=((Myp/M3p)*(m3/my)-1.)*relaxy
  fvec(4)=((M0p/M3p)*(m3/m0)-1.)

  if (any(md<0)) then
     fvec=fvec*1000
  endif
return
End subroutine fcn_6p
! }}}

! fcn_4p: {{{
subroutine fcn_4p(n,x,fvec,iflag)
use micro_prm, only:ihyd,nubounds,dnbounds,nkr,diams, &
    M3p,M0p,Mwp,Mzp,momw,momz,col,skr,ekr,rxfinal,relaxw,nug1,nug2, &
    nuterm31, nuterm32, nutermz1, nutermz2, nutermw1, nutermw2
use parameters, only: max_nbins
implicit none
integer n,iflag
double precision x(n),fvec(n)

double precision rx, nu1, dn1, nu2, dn2, MzM3, m3, m0, mw, mz
double precision m31_ana, m32_ana, m3_ana, m01_ana, m02_ana, m0_ana,&
                 mw1_ana, mw2_ana, mw_ana, mwfrac, m0_frac, m3frac, &
                 M3Mw, M31Mw1, M32Mw2, M3M0, M31M01, M32M02, &
                 m31mw1_ana, m32mw2_ana
double precision, dimension(max_nbins):: md
integer kr

!The purpose of this routine is to calculate the number/mass distribution
!of drops given parameters of the gamma distribution
!f(D)=n0*(D/dn)^(nu+3-1)*exp(-D/dn)
!The +3 is added because this is really a mass distribution, but we want
!the shape parameter for a number distribution.
!n0=(total mass)/gamma(nu+3)
!This function is used in conjunction with the equation solving routines.
!fvec characterizes how close the given parameters are to a solution.

  rx=1.e-10 !This is an arbitrary value of rx
  nu1=nug1 !Shape parameter
  dn1=x(1) !Scale diameter multiplied by 1e6 to make its magnitude 
                !approximately the same as nu
  nu2=nug2 !Shape parameter
  dn2=x(2) !Scale diameter multiplied by 1e6 to make its magnitude 
                !approximately the same as nu
  MzM3 = Mzp/M3p

  md=0. !Initialize the mass distribution
  ! print*, 'fcn_4p'
  ! call get_relative_moment_vals(nu1, dn1, nu2, dn2, m3, m0, mw)

  ! call get_relative_moment_vals(nu1, dn1, nu2, dn2, m3_ana, m0_ana, mw_ana)
  call incgamma2(rx,nu1,dn1,nu2,dn2,MzM3,md,m3frac)
  !Calculate the moments - 0th, 3rd, wth, xth, yth, and zth
  call calcmom(m3, md, 3)
  call calcmom(m0, md, 0)
  call calcmom(mw, md, momw)

  fvec(1)=log10((M0p/M3p)*(m3/m0))
  fvec(2)=log10((Mwp/M3p)*(m3/mw))*relaxw

  if (any(md<0)) then
     fvec=fvec*1000
  endif
return
End subroutine fcn_4p

! }}}

! calcerr2: {{{
subroutine calcerr2(rx,nu1,dn1,nu2,dn2,frac,error)
use micro_prm, only:nkr,diams,M0p,M3p,Mwp,Mxp,Myp,Mzp,col,relax,momw,momx,momy,momz
use parameters, only: max_nbins, num_h_moments

implicit none
double precision :: ratio,m3,m0,mw,mx,my,mz,md(max_nbins)
double precision :: error(4),rx,nu1,dn1,nu2,dn2,frac,m3frac
integer :: kr

  md = 0.

  ! call get_relative_moment_vals(nu1, dn1, nu2, dn2, m3, m0, mw)

  ! print*, 'calcerr2'
  call incgamma2(rx,nu1,dn1,nu2,dn2,frac,md,m3frac)
  !Calculate the moments - 3rd, xth, and yth
  call calcmom(m3, md, 3)
  call calcmom(m0, md, 0)
  call calcmom(mw, md, momw)
  call calcmom(mx, md, momx)
  call calcmom(my, md, momy)
  call calcmom(mz, md, momz)


  !Calculate the errors in the moments

  if (num_h_moments(1) == 4) then 
     ratio = (M0p/M3p)*(m3/m0)
     error(1) = log10(ratio)
     if (ratio<0) error(1) = log10(-ratio)*1000

     ratio = (Mwp/M3p)*(m3/mw)
     error(2) = log10(ratio)
     if (ratio<0) error(2) = log10(-ratio)*1000

     error(3)=0.
     error(4)=0.
  else
     ratio = (Mwp/mw)*(m3/M3p)
     error(1)=(ratio)-1.

     ratio = (Mxp/mx)*(m3/M3p)
     error(2)=(ratio)-1.

     ratio = (Myp/my)*(m3/M3p)
     error(3)=(ratio)-1.

     ratio = (M0p/m0)*(m3/M3p)
     error(4)=(ratio)-1.
   endif

!print*,'calc2',m3/m0,M3p/M0p
  if (any(md<0)) error=error*1000

End subroutine calcerr2

! }}}
! incgamma2: {{{
subroutine incgamma2(rx,nu1,dn1,nu2,dn2,MzM3,md,m3frac)
use micro_prm, only:diams,nkr,col,momz,tempvar_debug,l_printflag,tempvar_debug2, &
   Mp, dnc_def, dnr_def, nuterm31, nuterm32, nutermz1, nutermz2, momw
use parameters, only: max_nbins
implicit none

double precision :: rx,nu1,dn1,nu2,dn2,m3frac,n0,expterm1,expterm2
double precision :: gterm1a(max_nbins),gterm1b,gterm2a(max_nbins),gterm2b
double precision :: md(max_nbins), md1(max_nbins), md2(max_nbins), m31,m32,mz1,mz2,MzM3
double precision :: infinity
double precision :: m31_ana, m32_ana, mz1m31_ana, mz2m32_ana, t1, t2, t3, t4, &
   m3frac_ana, mw1, mw2, dummy
integer :: kr, dn_change_counter

dn_change_counter = 0
infinity = HUGE(rx)

! 100 continue

n0=rx

if (dn1 > dn2) then
   dummy = dn1
   dn1 = dn2
   dn2 = dummy
endif

do kr=1,nkr
   gterm1a(kr)=-diams(kr)/dn1+(nu1+3)*log(diams(kr)/dn1)
   gterm2a(kr)=-diams(kr)/dn2+(nu2+3)*log(diams(kr)/dn2)
enddo

gterm1b = maxval(gterm1a(1:nkr))*(-1.)
gterm2b = maxval(gterm2a(1:nkr))*(-1.)

do kr=1,nkr
   md1(kr)=exp(gterm1a(kr)+gterm1b)
   md2(kr)=exp(gterm2a(kr)+gterm2b)
enddo
call calcmom(m31, md1, 3)
call calcmom(m32, md2, 3)

if (m31.eq.0) then
   mz1=0.
else
   call calcmom(mz1, md1, momz)
   mz1 = mz1/m31
   ! mz1m31_ana = dn1**momz*nutermz1/m31_ana
endif

if (m32.eq.0) then
   mz2=0.
else
   call calcmom(mz2, md2, momz)
   mz2 = mz2/m32
   ! mz2m32_ana = dn2**momz*nutermz2/m32_ana
endif

! print*, 'incg m31, m32', m31, m32
! print*, 'incg m31_ana, m32_ana', m31, m32
if (dn1 <= 0. .or. m31>infinity .or. mz1>infinity) then 
   m3frac=0.
   mz1=0;m31=1
   md = md2/m32
   dn1 = dnc_def
elseif (dn2 <= 0. .or. m32>infinity .or. mz2>infinity) then
   m3frac=1.
   mz2=0;m32=1
   md = md1/m31
   dn2 = dnr_def
else
   m3frac=(MzM3-mz2)/(mz1-mz2)
   ! m3frac=(MzM3-mz2m32_ana)/(mz1m31_ana-mz2m32_ana)
   ! if (l_printflag) print*, m3frac

   if (m3frac < 0) then
      ! give up nudging the dn has been nudged up by a lot
      m3frac = 0.
      dn1 = dnc_def
      ! if (dn_change_counter > 0) then
      ! else
      !    dn_change_counter = dn_change_counter + 1
      !    dn2 = dn2*1.01
      !    go to 100
      ! endif
   elseif (m3frac > 1) then
      ! print*, m3frac, dn1, dn2
      ! print*, Mp
      m3frac = 1.
      dn2 = dnr_def
      ! might need to do a similar treatment to one of the dns, 
      ! but i'll do it when the problem occurs
      ! dn2 = 0.
   endif

   md = md1/m31*m3frac + md2/m32*(1.-m3frac)
endif


! if (m3frac>0. .and. m3frac<1.) then
!    print*, 'new/old', (t3-t2)/(t2-t1)
!    print*, 'm3frac', m3frac
!    print*, 'm3frac_ana', m3frac_ana
! endif

! tempvar_debug(1:9) = (/dn1, dn2, m31, m32, mz1, mz2, MzM3, (MzM3-mz2)/(mz1-mz2), m3frac/)
! tempvar_debug2 = md1/m31*(MzM3-mz2)/(mz1-mz2)+ md2/m32*(1.-(MzM3-mz2)/(mz1-mz2))


! tempvar_debug = (/dn1, dn2, m31, m32, m3frac/)


return
End subroutine incgamma2

! }}}
! calcdist2: {{{
subroutine calcdist2(n,x,md)
use micro_prm, only:nkr,diams,M3p,Mzp,col,skr,ekr,rxfinal,fracfinal,momx,momy
use parameters, only: max_nbins
implicit none
integer n
double precision x(n+2) ! feel like it should have a +2 here...

double precision rx, nu1, dn1, nu2, dn2, frac,m3, m3frac
double precision, dimension(max_nbins):: md

!Use this function to force correct 3rd moment 
!Very similar to fcn_6p
rx=100.
nu1=x(1)
dn1=x(2)
nu2=x(3)
dn2=x(4)
frac=Mzp/M3p

md=0.
call incgamma2(rx,nu1,dn1,nu2,dn2,frac,md,m3frac)

rxfinal = M3p
md = md*M3p

fracfinal=frac

return
End subroutine calcdist2

! }}}
! ! searchParamsByFrac: {{{
! Subroutine searchParamsByFrac(frac_m3, frac_m0, md, npm, flag)

! use micro_prm, only:nkr,momx,momy
! use parameters, only:max_nbins

! implicit none
! real(8) :: frac_m3, frac_m0
! real(8) :: md(max_nbins)
! integer :: npm,flag

! if (npm==4) then !if rain is 2M, then cloud is also 2M
!    CALL searchParamsBF4M(frac_m3,frac_m0,md,flag)
! else
!    CALL searchParamsBF6M(frac_m3,frac_m0,md,flag)
! endif

! End subroutine searchParamsByFrac

! ! }}}
! searchparamsG2: {{{
Subroutine searchparamsG2(guess1,guess2,md,npm,flag)

use micro_prm, only:nkr,momx,momy
use parameters, only:max_nbins

implicit none
real(8) :: guess1(2),guess2(2),bigguess(4)
real(8) :: ld,ud,md(max_nbins)
integer :: ihyd,npm
real(8) :: flag

bigguess(1:2)=guess1
bigguess(3:4)=guess2
if (npm==4) then !if rain is 2M, then cloud is also 2M
   CALL searchparams4M(bigguess,md,flag)
else
   CALL searchparams6M(bigguess,md,flag)
endif

guess1=bigguess(1:2)
guess2=bigguess(3:4)

End subroutine searchparamsG2

! }}}
! searchparams6M: {{{
Subroutine searchparams6M(guess,md,flag)

use micro_prm, only:relaxw,relaxx,relaxy,Mp,M3p,M0p,Mwp,Mxp,Myp,Mzp,momx,nkr, &
               nutab,dntab,minmaxmx,mintab,maxtab,ntab,dnc_def,dnr_def
use parameters, only: max_nbins

implicit none
real(8) :: guess(4),oguess(4),bguess(4),vals(4),ovals(4),bvals(4),tol
real(8) :: M0M3,MwM3,MxM3,MyM3,MzM3,irm,iry,wgtm,wgty1,wgty2
real(8) :: minmy1,maxmy1,minmy2,maxmy2,sqval,osqval,bsqval,minsqval,md(max_nbins),min12,max12
real(8) :: minMx3,maxMx3,temp
real(8) :: rand1(4),rand2(4)
integer :: i,info,n,ihyd,im1,im2,iy1a,iy2a,ix1b,ix2b,flag,ntry,ibguess,iimp
integer, parameter :: lwa=max_nbins
real(8),dimension(lwa) :: wa
external :: fcn_6p

flag = 0
!This subroutine drives the search for PDF parameters that satisfy the 
!moments. 
M3p=Mp(1);M0p=Mp(2);Mwp=Mp(3);Mxp=Mp(4);Myp=Mp(5);Mzp=Mp(6)
MzM3=Mzp/M3p

relaxw=1. !Initial value of relax. Do not change.
relaxx=1.
relaxy=1.

tol=1.d-4 !Tolerance of minimization routines
n=4 !Number of equations to solve. Do not change.

info=0 !Message sent back from minimization routines.

oguess=guess !oguess is the current best guess
ovals=0.
bvals=100
bsqval=100. !back-up values
bguess=guess
ibguess=0
iimp=0

CALL calcerr2(100.d0,oguess(1),oguess(2),oguess(3), &
                     oguess(4),MzM3,ovals)
osqval = sqrt(sum(ovals**2))

vals=100.
!First try - previous values
if (guess(2).ne.0) then
  CALL hybrd1(fcn_6p,n,guess,vals,tol,info,wa,lwa)
  sqval = sqrt(sum(vals**2))
  if (sqval .lt. osqval) then! .and. guess(2)>0. .and. guess(4)>0.) then
     ovals = vals
     osqval = sqval
     oguess = guess
  endif
endif

guess=oguess

!If first try was wildly off, or no 1st guess exists - use random first guesses
if (guess(2).eq.0 .or. info .ne. 1) then

   if (oguess(1).ge.100) oguess(1)=4.
   if (oguess(3).ge.100) oguess(3)=4.
   if (oguess(2).ge.10000 .or. oguess(2) .le.1e-3) oguess(2)=dnc_def
   if (oguess(4).ge.10000 .or. oguess(4) .le.1e-3) oguess(4)=dnr_def

   ntry = 0
   i=1
   do i=1,10000
      CALL random_number(rand1)
      CALL random_number(rand2)

      guess(1)=oguess(1)+sqrt(-2*9*log(rand1(1)))*cos(2.*3.14159*rand2(1)) !st. dev. of 5
      guess(3)=oguess(3)+sqrt(-2*9*log(rand1(3)))*cos(2.*3.14159*rand2(3)) !st. dev. of 5
      guess(2)=oguess(2)*exp(sqrt(-2*0.25*log(rand1(2)))*cos(2.*3.14159*rand2(2)))
      guess(4)=oguess(4)*exp(sqrt(-2*0.25*log(rand1(4)))*cos(2.*3.14159*rand2(4)))

      CALL calcerr2(100.d0,guess(1),guess(2),guess(3), &
                      guess(4),MzM3,vals)
  
      sqval = sqrt(sum(vals**2.))

      if (sqval .lt. osqval)then! .and. guess(2)>0. .and.guess(4)>0.) then
         ovals = vals
         osqval = sqval
         oguess = guess
         iimp = 1
      endif

      if (sqval/sqval.ne.1) then
         print*,guess,oguess
         stop 'sqval NaN'
      endif

      if (sqval < 1) then !if it looks like a good first guess, try and see
         !Use the best-guess. See what happens
         ntry = ntry + 1

         CALL hybrd1(fcn_6p,n,guess,vals,tol,info,wa,lwa)
         sqval = sqrt(sum(vals**2))

         if (sqval .lt. osqval)then ! .and. guess(2)>0. .and.guess(4)>0.) then
            ovals = vals
            osqval = sqval
            oguess = guess
            iimp = 1
         endif
         !if (info.eq.1) print*,'info 1'
         if (info.eq.1) exit
      endif

      if (abs(vals(4))<0.25 .and. sqval<bsqval) then
            bguess = guess
            bsqval = sqval
            bvals = vals
            ibguess = 1
      endif

      if ((ntry.ge.5.or.i>=5000) .and. ibguess == 1) then
         !If we've tried a lot and failed, then let's try relaxing error 
         !on momw,momx,momy but not M0.
         guess = bguess
         CALL hybrd1(fcn_6p,n,guess,bvals,tol,info,wa,lwa)
         if (abs(bvals(1))>tol) relaxw=max(tol/abs(bvals(1)),1.e-4)
         if (abs(bvals(2))>tol) relaxx=max(tol/abs(bvals(2)),1.e-4)
         if (abs(bvals(3))>tol) relaxy=max(tol/abs(bvals(3)),1.e-4)
         CALL hybrd1(fcn_6p,n,guess,vals,tol,info,wa,lwa)
         vals(1)=vals(1)/relaxw
         vals(2)=vals(2)/relaxx
         vals(3)=vals(3)/relaxy

         sqval = sqrt(sum(vals**2))

         if (sqval .le. osqval .or. (abs(ovals(4))>.01 .and. abs(vals(4))<.01))then 
! .and. guess(2)>0. .and.guess(4)>0.) then
           ovals = vals
           osqval = sqval
           oguess = guess
           iimp = 1
         endif

!        if (info.eq.1) print*,'bg info 1'
         if (info.eq.1 .and. info.eq.4) exit

!        if(abs(ovals(4))<.01 .and. ntry>=10) print*,'bg ovals4 .01',ovals(4),ntry
!        if(abs(ovals(4))<.01) print*,'bg ovals4 .01',ovals(4),ntry,ovals
         if(abs(ovals(4))<.01 .and. ntry>=10) exit
         ibguess=0
         relaxw=1.
         relaxx=1.
         relaxy=1.
      endif
!print*,'in loop',i,osqval,oguess(1),iimp
   enddo
!else
! print*,'first guess worked'
endif

if (iimp==1) then
!iimp==1 if guess improved with the random guessing prodecure
 guess = oguess
 sqval = osqval
 vals = ovals
endif
!print*,'numiter',i,iimp,vals(4)

!Set flag to 1 if fitting didn't work as well as we wished
if (sqval>tol) flag = 1

!Force third moment to have no error and calculate final distribution
CALL calcdist2(n,guess,md)
!print*,'guess',guess
! print*, 'tries', i
return

End subroutine searchparams6M

! }}}
! searchparams4M: {{{
Subroutine searchparams4M(guess,md,flag)

use micro_prm
use parameters, only: max_nbins
use global_fun

implicit none
real(8) :: guess(4),oguess(4),bguess(4),vals(4),ovals(4),bvals(4),tol,tguess(2),&
   lut_guess(2), pguess(4), mguess(4), mvals(4), mval1_mat(1001,10001), mval2_mat(1001,10001)
real(8) :: M0M3,MwM3,MxM3,MyM3,MzM3,irm,iry,wgtm,wgty1,wgty2
real(8) :: minmy1,maxmy1,minmy2,maxmy2,sqval,osqval,bsqval,minsqval,md(max_nbins),min12,max12,md1(max_nbins),md2(max_nbins)
real(8) :: minMx3,maxMx3,temp
real(8) :: rand1(4),rand2(4)
integer :: i,j,info,n,im1,im2,iy1a,iy2a,ix1b,ix2b,ntry,ibguess,iimp,idx_dns, itrial
integer, parameter :: lwa=max_nbins
real(8),dimension(lwa) :: wa
external :: fcn_4p
real(8), dimension(:), allocatable :: recorded_dns, recorded_moms
logical :: l_sctab_mod_curr, l_recorded, l_getlut
real(8):: error_4m, lut_skip
real(8) :: guess2m(2), M3t, M0t, dummy, Mwt, Mzt
real(8), dimension(max_nbins) :: dist1, dist2, distfull
logical :: l_unimodal, momw_close, momz_close
real(8) :: m31, m32, mz1, mz2, frac, vals2m, m3, tol2m
real(8) :: fguess(2), m3frac, m0frac
real(8), allocatable :: Mw_err(:,:), Mz_err(:,:)
real(8) :: flag, t1, t2, t3

call CPU_TIME(t1)
if (.not. allocated(tempvar_debug)) allocate(tempvar_debug(13))
flag = 0.
i = 0
l_failed1 = .false.
l_sctab_mod_curr = .false.
l_getlut = .false.
! l_massnoguess = .false.
tol2m = 1.d-8

!This subroutine drives the search for PDF parameters that satisfy the 
!moments. 
M3p=Mp(1);M0p=Mp(2);Mwp=Mp(3);Mzp=Mp(4)
Mxp=0.;Myp=0.
MzM3=Mzp/M3p
nug1 = guess(1)
nug2 = guess(3)

if (M3p <= total_m3_th) then
   md = 0.
   return
endif

relaxw=1. !Initial value of relax. Do not change.
tol=1.d-4 !Tolerance of minimization routines
n=2 !Number of equations to solve. Do not change.
info=0 !Message sent back from minimization routines.

oguess=guess !oguess is the current best guess
pguess=guess !pguess always store the previous guess and should be a constant
ovals=0.
vals=0
bvals=100
bsqval=100.
bguess=guess
ibguess=0
iimp=0.
error_4m=0.
ntry = 0


! ahu: map the error space around the previous guess and the correct Dn
! map a factor of 0.95 - 1.05 for dn1 and dn2 with precision of 0.001 and 0.0005 
! (for easier postprocess)

! mguess = guess ! m = map
! do i = 1,1001
!    print*, i
!    mguess(2) = guess(2)*(0.95+(i-1)/10000.)
!    do j = 1,10001
!       mguess(4) = guess(4)*(0.75+(j-1)/20000.)
!       call calcerr2(100.d0, mguess(1), mguess(2), mguess(3), mguess(4), MzM3, mvals)
!       mval1_mat(i,j) = mvals(1)
!       mval2_mat(i,j) = mvals(2)
!    enddo
! enddo
! open(55, file = "./output/M0_gerr.txt")
! sigfig_fmt = '(1001e'//itoa(7+8)//'.'//itoa(7)//')'
! write(55, trim(sigfig_fmt)) mval1_mat
! close(55)
! open(55, file = "./output/Mw_gerr.txt")
! sigfig_fmt = '(1001e'//itoa(7+8)//'.'//itoa(7)//')'
! write(55, trim(sigfig_fmt)) mval2_mat
! close(55)


CALL calcerr2(100.d0,oguess(1),oguess(2),oguess(3), &
                     oguess(4),MzM3,ovals)
! if (l_printflag) print*, 'initial oguess', oguess((/2,4/))
! if (l_printflag) print*, 'initial dn1 frac', tempvar_debug(8)
! if (l_printflag) print*, 'initial m31, m32, mz1, mz2, MzM3', tempvar_debug(3:7)
! ! if (l_printflag) print*, 'md with uncortd frac', tempvar_debug2
! if (l_printflag) print*, 'initial ovals', ovals(1:2)

osqval = sqrt(sum(ovals(1:2)**2))
! if (l_printflag) print*, i, tguess, ovals(1:2)

vals=100.

!First try - previous values
if (guess(2).ne.0 .or. guess(4).ne.0) then
   tguess(1) = guess(2)
   tguess(2) = guess(4)
   CALL hybrd1(fcn_4p,n,tguess,vals(1:2),tol,info,wa,lwa)

   ! ! reorder dn if they are not ordered correctly
   ! if (tguess(2) < tguess(1)) then
   !    dummy = tguess(1)
   !    tguess(1) = tguess(2)
   !    tguess(2) = dummy
   ! endif

   sqval = sqrt(sum(vals(1:2)**2))
   ! print*, 'sqval,osqval', sqval, osqval
   if (sqval .lt. osqval) then
     ovals = vals
     osqval = sqval
     oguess = (/guess(1),tguess(1),guess(3),tguess(2)/)
     bguess = oguess
   endif
endif
if (l_printflag) print*, 'vals after 1st hybrd', vals(1:2)
if (l_printflag) print*, 'guess after 1st hybrd', tguess

if (osqval .ne. osqval) osqval = 1.e10 ! an arbitrary large number

guess = oguess

! print*, 'guesses', tguess((/1,2/))

call CPU_TIME(t2)
! If first two tries were wildly off, or no 1st guess exists - try backtracking
! the previous moment ratios and start from there. 
! if that doesn't work, then random guesses
! do itrial = 1,2
if (info .ne. 1) then
   ! if (itrial == 1) oguess = guess
   ! if (itrial == 2) oguess = pguess
   
   if (l_getlut) then
      oguess(2) = lut_guess(1)
      oguess(4) = lut_guess(2)
   else
      if (l_printflag) print*, 'oguess after 1st hybrd', oguess((/2,4/))
      if (oguess(2).ge.DnRangeMax .or. oguess(2) .le. DnRangeMin) oguess(2)=dnc_def
      if (oguess(4).ge.DnRangeMax .or. oguess(4) .le. DnRangeMin) oguess(4)=dnr_def
   endif

   do i=1,10000
      CALL random_number(rand1)
      CALL random_number(rand2)

      tguess(1)=guess(2)*exp(sqrt(-2*0.25*log(rand1(2)))*cos(2.*3.14159*rand2(2)))
      tguess(2)=guess(4)*exp(sqrt(-2*0.25*log(rand1(4)))*cos(2.*3.14159*rand2(4)))

      CALL calcerr2(100.d0,guess(1),tguess(1),guess(3), &
                      tguess(2),MzM3,vals)
      vals(3:4)=0.
      sqval = sqrt(sum(vals(1:2)**2.))

      if (sqval .lt. osqval) then! .and. guess(2)>0. .and.guess(4)>0.) then
         ovals = vals
         osqval = sqval
         oguess = (/guess(1),tguess(1),guess(3),tguess(2)/)
         iimp = 1
      endif

      if (sqval/sqval.ne.1) then
         print*, 'tguess', tguess
         print*, 'pguess', pguess
         print*, 'Mp', Mp
         print*, 'vals', vals
         print*, 'sqval NaN'
         md = 0.
         stop
         ! return
      endif

      if (sqval < 1) then
         !Use the best-guess. See what happens
         ntry = ntry + 1
         CALL hybrd1(fcn_4p,n,tguess,vals(1:2),tol,info,wa,lwa)
         sqval = sqrt(sum(vals(1:2)**2))
         if (sqval .lt. osqval) then! .and. guess(2)>0. .and.guess(4)>0.) then
            ovals = vals
            osqval = sqval
            oguess = (/guess(1),tguess(1),guess(3),tguess(2)/)
            iimp = 1
         endif
         if (info.eq.1) then
            exit
         endif
      endif

      if (abs(vals(2))<0.25 .and. sqval<bsqval) then
         bguess = (/guess(1),tguess(1),guess(3),tguess(2)/)
         bsqval = sqval
         bvals = vals
         ibguess = 1
      endif

      if (ntry.ge.5.or.i>=5000) then
        !If we've tried a lot and failed, then let's try relaxing the error
        !on momw but not M0
        tguess(1)=bguess(2)
        tguess(2)=bguess(4)
        
        CALL hybrd1(fcn_4p,n,tguess,bvals(1:2),tol,info,wa,lwa)

        if (abs(bvals(2))>tol) relaxw=max(tol/abs(bvals(2)),1.e-4)

        CALL hybrd1(fcn_4p,n,tguess,vals(1:2),tol,info,wa,lwa)
        vals(2)=vals(2)/relaxw

        sqval = sqrt(sum(vals(1:2)**2.))

        if (sqval .le. osqval .or. (abs(ovals(2))>.01 .and. abs(vals(1))<.01)) then
           ! if (all(tguess>0)) then
              ovals = vals
              osqval = sqval
              oguess = (/guess(1),tguess(1),guess(3),tguess(2)/)
              iimp = 1
           ! endif
        endif

        if (info.eq.1 .or. info .eq. 4) exit

        if(abs(ovals(2))<.01 .or. ntry>=10) exit
        ibguess=0
        relaxw=1.
      endif

      
   enddo 

else
endif
! enddo

call CPU_TIME(t3)

diag_dt1 = diag_dt1 + t2-t1
diag_dt2 = diag_dt2 + t3-t2
itries = itries + i


! print*, '2nd/1st', (t3-t2)/(t2-t1)
! print*, 'total tries', i


if (iimp==1) then
!iimp==1 if guess is improved with the random guessing procedure
   guess = oguess
   sqval = osqval
   vals = ovals
endif

! if (i>=1) then
!    ! print*, 'tries', i
!    ! print*, 'guess(2)/pguess(2)', guess(2)/pguess(2)
!    ! print*, 'guess(4)/pguess(4)', guess(4)/pguess(4)
!    PF_change_arr = (/dble(i),guess(2)/pguess(2),guess(4)/pguess(4),sqval,dble(info)/)
! endif

! if (any(guess<0)) then
!    print*, 'guess<0'
!    print*, 'stop'
!    print*, 'iimp', iimp
!    stop
! endif


tguess = guess((/2,4/))

! make sure guesses are in the right order

! if (guess(2)>guess(4)) then
!    tguess = guess(1:2)
!    guess(1:2) = guess(3:4)
!    guess(3:4) = tguess
! endif

if (guess(2) == 0. .and. guess(4) == 0.) then
   ! if (tempvar_debug(1) .ne. 0. .or. tempvar_debug(2) .ne. 0) then
   !    l_massnoguess = .true.
   !    ! print*, 'Dns', tempvar_debug(1:2)
   !    ! print*, 'tguesses', tguess
   !    ! print*, 'frac', tempvar_debug(8)
   !    ! print*, 'Mp', Mp(1:4)
   !    ! print*, 'meand', get_meandiam(Mp(1), Mp(2))
   !    ! print*, 'sqval', sqval
   !    ! print*, 'osqval', osqval
   !    ! print*, 'improved?', iimp
   !    ! stop
   ! endif
   md=0.
   return
endif


! l_printflag = .true.

!Force third moment to have no error and calculate final distribution
CALL calcdist2(n,guess,md)

if (guess(2) > DnRangeMax .or. guess(2) < DnRangeMin) guess(2) = dnc_def
if (guess(4) > DnRangeMax .or. guess(4) < DnRangeMin) guess(4) = dnr_def

if (vals(2)>tol) flag=1
if (vals(1)>tol) flag=2
if (vals(1)>tol .and. vals(2)>tol) flag=3

! if (l_printflag) print*, 'final error', vals(1:2)
! ! print*, 'frac', tempvar_debug(9)
! call calcmom(m3t, md, 3)
! call calcmom(m0t, md, 0)
! call calcmom(mwt, md, momw)
! call calcmom(mzt, md, momz)
! tempvar_debug(10:13) = (/m3t, m0t, mwt, mzt/)
! if (l_printflag) print*, 'info', info
! if (l_printflag) print*, 'tries', i
! if (l_printflag) print*, 'ntries', ntry
! if (l_printflag) print*, 'dn1 fraction', tempvar_debug(9)
! if (l_printflag) print*, 'guesses', guess((/2,4/))
! if (l_printflag) print*, 'moments predicted', Mp(1:4)
! if (l_printflag) print*, 'moments after pf ', (/m3t, m0t, mwt, mzt/)


!Set flag to 1 if fitting didn't work as well as we wished

! print*, 'dn1 fraction', tempvar_debug(9)
! print*, 'guesses', guess((/2,4/))
! print*, 'total_m3_th', total_m3_th

! if (sqval>tol) then
!    ! print*, 'Mp predicted', Mp(1:4)
!    ! print*, 'Mp after pf ', (/m3t,m0t,mwt,mzt/)
!    flag = 1.
! endif

! if (l_printflag) then
!    ! print*, 'predicted moms', Mp(1:4)
!    print*, 'moms predicted', Mp(1:4)
!    print*, 'moms after pf ', (/m3t, m0t, mwt, mzt/)
!    print*, 'meand predicted', get_meandiam(Mp(1), Mp(2))
!    print*, 'meand after pf ', get_meandiam(m3t, m0t)
!   nug=guess(1)
!   guess2m(1) = nug
!   call guessparams2M(Mp(1), Mp(4), 3, 9, guess2m, tol2m, vals2m)
!   call calcdist(guess2m, md)
!   call calcmom(m3t, md, 3)
!   call calcmom(m0t, md, 0)
!   call calcmom(mwt, md, momw)
!   call calcmom(mzt, md, momz)
!   ! print*, 'md assuming 1 mode', md
!   print*, 'moms assuming 1 mode', (/m3t, m0t, mwt, mzt/)
!   print*, 'meand assuming 1 mode', get_meandiam(m3t, m0t)
!   ! stop
! endif

! if (info .eq. 4) then
!    momw_close = are_close(mwt, Mp(3), 'ratio', tol_rat_opt=1d0)
!    momz_close = are_close(mzt, Mp(4), 'ratio', tol_rat_opt=1d0)
!    ! print*, 'momw, momz similar', momw_close, momz_close
! endif

return

End subroutine searchparams4M
! }}}
