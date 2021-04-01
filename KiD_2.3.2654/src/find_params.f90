subroutine fcn_2p(n,x,fvec,iflag)
use micro_prm, only:ihyd,nubounds,dnbounds,nkr,diams,M3p,Mxp,Myp,momx,momy,skr,ekr,relax,nkr
use module_hujisbm, only:xl
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
!-----------------------------------------------
subroutine fcn_1p(n,x,fvec,iflag)
use micro_prm, only:diams,Mxp,M3p,skr,ekr,nug,momx,nkr
use module_hujisbm, only:xl
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
!------------------------------------------------------------------------
subroutine calcerr(rx,nu,dn,error)
use micro_prm, only:diams,M3p,Mxp,Myp,col,nkr,skr,ekr,relax,momx,momy
use module_hujisbm, only:xl
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
!------------------------------------------------------------------------
subroutine calcmom(mom,md,momnum)

use micro_prm, only: diams,col,nkr,binmass
use mphys_tau_bin_declare, only:dgmean,xk,xkgmean
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
!------------------------------------------------------------------------
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
!------------------------------------------------------------------------
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
!------------------------------------------------------------------------
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
!----------------------------------------------------------------------
subroutine calcdist(x,md)
use micro_prm, only:nkr,diams,M3p,col,skr,ekr,rxfinal,ihyd,dnbounds
use module_hujisbm, only:xl
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
  rx = 1.e-10
  n=0
 do
  n=n+1
  md=0.

  call incgamma_norm(rx,nu,dn,skr,ekr,md)
  !call incjohnsonsb(nu,dn*1.e6,skr,ekr,md)
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
!--------------------------------------------------------
Subroutine searchparamsG(guess,ihyd,md,flag)

use micro_prm, only:nkr,npm
use parameters, only:flag_count, max_nbins

implicit none
real(8) :: guess(2)
real(8) :: md(max_nbins)
integer :: ihyd
real(8),optional,dimension(flag_count) :: flag

if (npm==3) then !if rain is 2M, then cloud is also 2M
   CALL searchparams3M(guess,ihyd,md,flag)
else
   CALL searchparams2M(guess,md,flag)
endif

End subroutine searchparamsG
!--------------------------------------------------------
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
  guess(2) = (1.-wgtm)*((1.-wgty1)*dntab(iy1a,im1,ihyd)+wgty1*dntab(ix2b,im1,ihyd)) + &
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
    if (info.eq.1 .or. i>4) then
      flag(1)=vals(1)
      !if (info .ne.1)flag(1)=1
      !if (vals(1)==100.)flag(1)=2
      exit
    endif

    relax = relax*0.1
    i=i+1
  enddo

  guess=oguess !We're done. Set guess to best guess

endif

!Set flag to 1 if fitting didn't work as well as we wished
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
!--------------------------------------------------------
Subroutine searchparams2M(guess,md,flag)

use micro_prm, only:Mp,M3p,Mxp,ntab,skr,ekr,nkr,momx,nug,diams
use parameters, only: max_nbins
implicit none
real(8) :: guess(2),oguess(2),vals(1),tol,guessin(1)
real(8) :: MxM3,minmxm3,maxmxm3
real(8) :: md(max_nbins)
integer :: info,n,flag
integer, parameter :: lwa=max_nbins
real(8),dimension(lwa) :: wa
external :: fcn_1p

flag = 0
!This subroutine drives the search for PDF parameters that satisfy the
!moments.
M3p=Mp(1)
Mxp=Mp(2)

tol=1.d-8 !Tolerance of minimization routines
n=1 !Number of equations to solve. Do not change.
info=0 !Message sent back from minimization routines.

!Find min and max values of MxM3
minmxm3=min(diams(skr)**(momx-3),diams(ekr)**(momx-3))
maxmxm3=max(diams(skr)**(momx-3),diams(ekr)**(momx-3))

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
guess(2) = (MxM3*gamma(nug+3.)/gamma(nug+momx))**(1./(momx-3.))*1.e6
vals = 2.
guessin=guess(2)
CALL hybrd1(fcn_1p,n,guessin,vals,tol,info,wa,lwa)

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

return

End subroutine searchparams2M
