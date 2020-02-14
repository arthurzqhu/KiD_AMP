
subroutine fcn_2p(n,x,fvec,iflag)
use micro_prm, only:ihyd,nubounds,dnbounds,nkr,diams,M3p,Mxp,Myp,momx,momy,skr,ekr,relax
use module_hujisbm, only:xl
implicit none
integer n,iflag
double precision x(n),fvec(n)

double precision rx, nu, dn, m3, mx, my
double precision, dimension(nkr):: md

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
   .or. x(1)<nubounds(1,ihyd) .or. x(1)>nubounds(2,ihyd)) then
!If the parameters exceed the bounds, assign a high error
  fvec(1)=100.
  fvec(2)=100.
else
  rx=1.e-10 !This is an arbitrary value of rx
  nu=x(1) !Shape parameter
  dn=x(2)*1.e-6 !Scale diameter multiplied by 1e6 to make its magnitude
                !approximately the same as nu

  md=0. !Initialize the mass distribution
  call incgamma_norm(rx,nu,dn,skr,ekr,md)
!  call incjohnsonsb(nu,dn*1.e6,skr,ekr,md)

  !Calculate the moments - 0th, 3rd, and xth
  m3=sum(md/xl*diams**3)! don't need this part since will cancel in ratio *col*1000.
  mx=sum(md/xl*diams**momx)!*col*1000.
  my=sum(md/xl*diams**momy)!*col*1000.

  !Calculate the errors in the moments
  if (m3>0.) then
    fvec(1)=log10((Mxp/M3p)*(m3/mx))*relax**2.
    fvec(2)=log10((Myp/M3p)*(m3/my))
  else
    fvec=100.
  endif
endif
return
End subroutine fcn_2p
!-----------------------------------------------
subroutine fcn_1p(n,x,fvec,iflag)
use micro_prm, only:nkr,diams,Mxp,M3p,skr,ekr,nug,momx
use module_hujisbm, only:xl
implicit none
integer n,iflag
double precision x(n),fvec(n)

double precision rx, nu, dn, m3, mx
double precision, dimension(nkr):: md

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
  dn=x(1)*1.e-6 !x(1) is the scale diameter multiplied by 1e6 to make its magnitude
                !approximately the same as nu
  md=0. !Initialize the distribution
  call incgamma_norm(rx,nu,dn,skr,ekr,md)

  !Calculate the moments - xth and 3rd
  m3=sum(md/xl*diams**3)!*col*1000.
  mx=sum(md/xl*diams**momx)!*col*1000.

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

implicit none
double precision :: ratio,m3,mx,my,md(nkr),error(2),rx,nu,dn

  md = 0.
  call incgamma_norm(rx,nu,dn,skr,ekr,md)
!  call incjohnsonsb(nu,dn*1.e6,skr,ekr,md)

  !Calculate the moments - 3rd, xth, and yth
  m3=sum(md/xl*diams**3.)*col*1000.
  mx=sum(md/xl*diams**momx)*col*1000.
  my=sum(md/xl*diams**momy)*col*1000.

  !Calculate the errors in the moments
  ratio = (Mxp/mx)*(m3/M3p)
  !if (ratio>1.) ratio = 1./ratio
  !error(1)=(1.-ratio)*relax**2
  error(1)=log10(ratio)*relax**2.

  ratio = (Myp/my)*(m3/M3p)
  error(2)=log10(ratio)

End subroutine calcerr
!------------------------------------------------------------------------
subroutine incgamma(rx,nu,dn,ia,iz,md)
use micro_prm, only:diams,nkr
implicit none

double precision :: rx,nu,dn,n0,expterm,md(nkr)
integer :: ia,iz,kr

  n0=rx/gamma(nu+3)
  do kr=ia,iz !Loop over bins
    expterm=exp(-1.*diams(kr)/dn)
    md(kr)=n0*expterm*(diams(kr)/dn)**(nu+3.)
    !Check for NaNs
!    if (md(kr).ne.md(kr) .or. md(kr)*0.0 .ne.0.0 &
!       .or. md(kr)/md(kr) .ne. 1.0) md(kr)=0.
  enddo

return
End subroutine incgamma
!------------------------------------------------------------------------
subroutine incgamma_norm(rx,nu,dn,ia,iz,md)
use micro_prm, only:diams,nkr
implicit none

double precision :: rx,nu,dn,n0,expterm,md(nkr)
integer :: ia,iz,kr

  n0=rx
  do kr=ia,iz !Loop over bins
    expterm=exp(-1.*diams(kr)/dn)
    md(kr)=n0*expterm*(diams(kr)/dn)**(nu+3.)
    !Check for NaNs
!    if (md(kr).ne.md(kr) .or. md(kr)*0.0 .ne.0.0 &
!       .or. md(kr)/md(kr) .ne. 1.0) md(kr)=0.
  enddo

return
End subroutine incgamma_norm
!------------------------------------------------------------------------
subroutine incjohnsonsb(s1,s2,ia,iz,md)
use micro_prm, only:diams,nkr
implicit none

double precision :: s1,s2,minx,ranx,z,expterm,md(nkr)
integer :: ia,iz,kr

minx=0.
ranx=diams(nkr)*2.

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
implicit none
double precision x(2)
integer n

double precision rx, nu, dn, m3
double precision, dimension(nkr):: md

!Use this function to force correct 3rd moment
!Very similar to fcn_2p
  rx=1.e-10
  nu=x(1)
  dn=x(2)*1.e-6
  n=0
 do
  n=n+1
  md=0.
  call incgamma_norm(rx,nu,dn,skr,ekr,md)
  !call incjohnsonsb(nu,dn*1.e6,skr,ekr,md)

  m3=sum(md/xl*diams**3)*col*1000.

  if ((m3==0. .or. dn==0.) .and. M3p<1.e-20) then
    rxfinal=0.
    md = 0.
    exit
  elseif (m3==0. .or. dn==0.) then
    print*,'cannot find rxfinal and mass = ',M3p*3.1415/6.*1000
    print*,m3,nu,dn,ihyd,dnbounds(:,ihyd)
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

implicit none
real(8) :: guess(2)
real(8) :: md(nkr)
integer :: ihyd,flag

if (npm==3) then !if rain is 2M, then cloud is also 2M
   CALL searchparams3M(guess,ihyd,md,flag)
else
   CALL searchparams2M(guess,md,flag)
endif

End subroutine searchparamsG
!--------------------------------------------------------
Subroutine searchparams3M(guess,ihyd,md,flag)

use micro_prm, only:relax,Mp,M3p,Mxp,Myp,momx,nkr, &
               nutab,dntab,minmaxmx,mintab,maxtab,ntab

implicit none
real(8) :: guess(2),oguess(2),vals(2),ovals(2),tol
real(8) :: MxM3,MyM3,irm,iry,wgtm,wgty1,wgty2
real(8) :: minmy1,maxmy1,minmy2,maxmy2,sqval,osqval,minsqval,md(nkr),min12,max12
real(8) :: minMx3,maxMx3
integer :: i,info,n,ihyd,im1,im2,iy1a,iy2a,ix1b,ix2b,flag
integer, parameter :: lwa=33
real(8),dimension(lwa) :: wa
external :: fcn_2p

flag = 0
!This subroutine drives the search for PDF parameters that satisfy the
!moments.
M3p=Mp(1);Mxp=Mp(2);Myp=Mp(3)

relax=1. !Initial value of relax. Do not change.
tol=1.d-8 !Tolerance of minimization routines
n=2 !Number of equations to solve. Do not change.

info=0 !Message sent back from minimization routines.

oguess=guess !oguess is the current best guess
ovals=0.
!print*,oguess
CALL calcerr(100.d0,oguess(1),oguess(2)*1.e-6,ovals)
osqval = sqrt(sum(ovals**2))

!First try - previous values
if (guess(2).ne.0) then
  CALL hybrd1(fcn_2p,n,guess,vals,tol,info,wa,lwa)
  sqval = sqrt(sum(vals**2))
  if (sqval .lt. osqval) then
  !if (abs(vals(1)) .lt. abs(ovals(1))) then
     ovals = vals
     osqval = sqval
     oguess = guess
  endif
endif

!If first try was wildly off, or no 1st guess exists - use look up table
if (guess(2).eq.0 .or. abs(vals(1))>1.0e-4) then
  !Ratio of xth moment to 3rd moment.
  MxM3 = Mxp/M3p
  !Ratio of yth moment to 3rd moment.
  MyM3 = Myp/M3p

  !Need to see if these ratios are in the solution space
  !If not, adjust Mx and/or My until they are in the solution space

  !Lower and upper limits of Mx/M3 can be calculated,
  !but I threw those limits out as being too hard to deal with
  !Actual limits are stored in the minmaxmx table
  minMx3 = minmaxmx(1,momx+1,ihyd)
  maxMx3 = minmaxmx(2,momx+1,ihyd)

  !Check now to see if MxM3 is out of allowable range and adjust Mx
  if (MxM3 < minMx3) then
    flag=-abs(minMx3/MxM3)
    MxM3 = minMx3*(1.+0.0)
    Mxp = MxM3 * M3p
  elseif (MxM3 > maxMx3) then
    flag=-abs(MxM3/maxMx3)
    MxM3 = maxMx3*(1.-0.0)
    Mxp = MxM3 * M3p
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
  min12=min(minmy1,minmy2)
  max12=max(maxmy1,maxmy2)
  if (MyM3 < min12) then
    flag=-abs(min12/myM3)
    MyM3 = min12*(1.+0.0)
    Myp = MyM3 * M3p
  elseif (MyM3 > max12) then
    flag=-abs(myM3/max12)
    MyM3 = max12*(1.-0.0)
    Myp = MyM3 * M3p
  endif

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

  !Find best-guess parameters in the look up tables
  guess(1) = (1.-wgtm)*((1.-wgty1)*nutab(iy1a,im1,ihyd)+wgty1*nutab(iy2a,im1,ihyd)) + &
             wgtm*((1.-wgty2)*nutab(ix1b,im2,ihyd)+wgty2*nutab(ix2b,im2,ihyd))
  guess(2) = (1.-wgtm)*((1.-wgty1)*dntab(iy1a,im1,ihyd)+wgty1*dntab(ix2b,im1,ihyd)) + &
             wgtm*((1.-wgty2)*dntab(ix1b,im2,ihyd)+wgty2*dntab(ix2b,im2,ihyd))

  !Use the best-guess. See what happens
  CALL hybrd1(fcn_2p,n,guess,vals,tol,info,wa,lwa)
  sqval = sqrt(sum(vals**2))
!print*,sqval,osqval
  if (sqval .lt. osqval) then
     ovals = vals
     osqval = sqval
     oguess = guess
  endif
  !print*,'c',guess, vals
endif

!If not close enough, jiggle the values and relax tolerance on Mx
!Start with current best guess
guess = oguess
vals = ovals
sqval = osqval
i=0
!print*,'b',vals,guess
if (abs(vals(1))>tol.and.abs(vals(1))<1000.) then
  i=0
  minsqval=ovals(1)
  do
    if(abs(vals(1))>tol) guess = guess*(1.01)**(-1*(i+1))
    CALL hybrd1(fcn_2p,n,guess,vals,tol,info,wa,lwa)
    sqval = sqrt(sum(vals**2))

    !Keep track of best guess. Don't worry about error on My
    !for determining best guess
    if (abs(vals(1))<abs(minsqval)) then
      minsqval=vals(1)
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
if (abs(vals(1))>tol .and. flag>=0) flag=3
if (abs(vals(2))>tol .and. flag>=0) flag=2
if (abs(vals(1))>tol .and. abs(vals(2))>tol .and. flag>=0) flag=4

!Force third moment to have no error and calculate final distribution
!print*,'a',guess
CALL calcdist(guess,md)
!print*,ihyd,guess
return

End subroutine searchparams3M
!--------------------------------------------------------
Subroutine searchparams2M(guess,md,flag)

use micro_prm, only:Mp,M3p,Mxp,ntab,skr,ekr,nkr,momx,nug,diams

implicit none
real(8) :: guess(2),oguess(2),vals,tol,guessin(1)
real(8) :: MxM3,minmxm3,maxmxm3
real(8) :: md(nkr)
integer :: info,n,flag
integer, parameter :: lwa=33
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

if (vals==1.) then
  guess=oguess
else
  guess(2)=guessin(1)
endif

!Set flag to 1 if fitting didn't work as well as we wished
if (abs(vals)>tol) flag = 1

!For 3M fitting, we jiggle the values and relax tolerance on Mx
!No sense in doing that here. The fitting is much easier, and should
!always give the same answer.
CALL calcdist(guess,md)
return

End subroutine searchparams2M
