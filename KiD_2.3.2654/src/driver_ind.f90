subroutine amp_init(aer2d,Mpc2d,Mpr2d,guessc2d,guessr2d)

use switches, only: zctrl
use parameters, only:nx,nz,num_h_moments,h_shape,max_nbins
use column_variables, only: z
use module_hujisbm
use micro_prm

implicit none
character(4)::momstr
character(len=100)::lutfolder
real(8), dimension(nz,nx,num_h_moments(1)) :: Mpc2d
real(8), dimension(nz,nx,num_h_moments(2)) :: Mpr2d
real(8), dimension(nz,nx,2) :: guessc2d,guessr2d
real(8),dimension(max_nbins) :: ffcd_mass, ffcd_num
real(8),dimension(num_h_moments(1)) :: mc
real(8),dimension(num_h_moments(2)) :: mr
real :: dnc,dnr
real, dimension(nz,nx,max_nbins) :: aer2d
integer :: i,k

!Set some parameters
ndtcoll=1
ipris=0; ihail=0; igraup=0;iceprocs=0;imbudget=1

pmomsc=(/3,imomc1,imomc2/)
pmomsr=(/3,imomr1,imomr2/)
!if (imomc1==imomc2) then
!  npm=2 !Run 2M scheme, Num. Pred. Moments = 2
!  print*,2
!else
!  npm=3 !Run 3M scheme, Num. Pred. Moments = 3
!  print*,3
!endif
npm=num_h_moments(1)

!Open output files, or decide that we don't need to run this combination of moments
if (imomc1>=7 .and. imomc1 .ne. imomc2) then
 print*,"not a valid moment combination",imomc1,imomc2
 stop
else !Run AMP
  !Open and read lookup tables
  if (npm==3) then
    lutfolder='./src/input_data/lutables/'

    write(momstr,'(A,I1,A,I1)') 'M',imomc1,'M',imomc2
    open(17,file=trim(lutfolder)//'cloud_nu_'//momstr//'.txt')
    read(17,*) nutab(:,:,1); close(17)
    open(17,file=trim(lutfolder)//'cloud_dn_'//momstr//'.txt')
    read(17,*) dntab(:,:,1); close(17)
    open(17,file=trim(lutfolder)//'cloud_minmaxmx.txt')
    read(17,*) minmaxmx(:,:,1); close(17)
    open(17,file=trim(lutfolder)//'cloud_minmy_'//momstr//'.txt')
    read(17,*) mintab(:,1); close(17)
    open(17,file=trim(lutfolder)//'cloud_maxmy_'//momstr//'.txt')
    read(17,*) maxtab(:,1); close(17)
    !Find min and max of parameters from the tables
    dnbounds(1,1)=minval(dntab(:,:,1)); dnbounds(2,1)=maxval(dntab(:,:,1))
    nubounds(1,1)=minval(nutab(:,:,1)); nubounds(2,1)=maxval(nutab(:,:,1))

    write(momstr,'(A,I1,A,I1)') 'M',imomr1,'M',imomr2
    open(17,file=trim(lutfolder)//'rain_nu_'//momstr//'.txt')
    read(17,*) nutab(:,:,2); close(17)
    open(17,file=trim(lutfolder)//'rain_dn_'//momstr//'.txt')
    read(17,*) dntab(:,:,2); close(17)
    open(17,file=trim(lutfolder)//'rain_minmaxmx.txt')
    read(17,*) minmaxmx(:,:,2); close(17)
    open(17,file=trim(lutfolder)//'rain_minmy_'//momstr//'.txt')
    read(17,*) mintab(:,2); close(17)
    open(17,file=trim(lutfolder)//'rain_maxmy_'//momstr//'.txt')
    read(17,*) maxtab(:,2); close(17)
    dnbounds(1,2)=minval(dntab(:,:,2)); dnbounds(2,2)=maxval(dntab(:,:,2))
    nubounds(1,2)=minval(nutab(:,:,2)); nubounds(2,2)=maxval(nutab(:,:,2))
  endif
endif

!call microphysics initialization routines
if (bintype .eq. 'sbm') then
    CALL micro_init_sbm()
    CALL micro_init_sbm2()
    CALL kernalsdt()
elseif (bintype .eq. 'tau') then
    CALL micro_init_tau()
endif

!Set up initial distribution and set moment values and parameter guesses
if (npm==3) then
  dnc=dntab(50,50,1); dnr=dntab(50,50,2)
else
  dnc=0.;dnr=0.
endif

if(cloud_init(1)>0.)dnc = (cloud_init(1)*6./3.14159/1000./cloud_init(2)*gamma(h_shape(1))/gamma(h_shape(1)+3))**(1./3.)
if(rain_init(1)>0.)dnr = (rain_init(1)*6./3.14159/1000./rain_init(2)*gamma(h_shape(2))/gamma(h_shape(2)+3))**(1./3.)
if(cloud_init(1)>0. .or. rain_init(1)>0.) then
    if (bintype .eq. 'sbm') then
        CALL init_dist_sbm(cloud_init(1),h_shape(1),dnc,rain_init(1),&
                           h_shape(2),dnr,diams,ffcd_mass)
    elseif (bintype .eq. 'tau') then
        CALL init_dist_tau(cloud_init(1),h_shape(1),dnc,rain_init(1),%
                           h_shape(2),dnr,diams,ffcd_mass,ffcd_num)
    endif
else
   ffcd_mass(:)=0.
   if (bintype .eq. 'tau') ffcd_num(:)=0.
endif

guessc2d(:,:,2)=dnc
guessr2d(:,:,2)=dnr

do i=1,num_h_moments(1)
    mc(i)=sum(ffcd_mass(1:krdrop)/xl(1:krdrop)*diams(1:krdrop)**pmomsc(i))*col*1000.
enddo
do i=1,num_h_moments(2)
    mr(i)=sum(ffcd_mass(krdrop+1:nkr)/xl(krdrop+1:nkr)*diams(krdrop+1:nkr)**pmomsr(i))*col*1000.
enddo

do i=1,nx
   do k=1,nz
      if (z(k)>=zctrl(2) .and. z(k)<=zctrl(3)) then
         Mpc2d(k,i,1:num_h_moments(1))=mc(1:num_h_moments(1))
         Mpr2d(k,i,1:num_h_moments(2))=mr(1:num_h_moments(2))
      endif
   enddo
enddo

end subroutine amp_init
!------------------------------------------------------------
subroutine sbm_init(aer2d,dropsm2d)

use switches, only: zctrl
use parameters, only:nx,nz,num_h_moments,h_shape,max_nbins
use column_variables, only: z
use module_hujisbm
use micro_prm

implicit none
real(8),dimension(max_nbins) :: ffcd
real :: dnc,dnr
real, dimension(nz,nx,max_nbins) :: aer2d,dropsm2d
integer :: i,k

!Set some parameters
ndtcoll=1
ipris=0; ihail=0; igraup=0;iceprocs=0;imbudget=1

!call microphysics initialization routines
CALL micro_init_sbm()
CALL micro_init_sbm2()
CALL kernalsdt()

!Set up initial distribution and set moment values and parameter guesses
dnc=0.;dnr=0.
if(cloud_init(1)>0.)dnc = (cloud_init(1)*6./3.14159/1000. &
                          /cloud_init(2)*gamma(h_shape(1)) &
                          /gamma(h_shape(1)+3))**(1./3.)
if(rain_init(1)>0.)dnr = (rain_init(1)*6./3.14159/1000. &
                         /rain_init(2)*gamma(h_shape(2)) &
                         /gamma(h_shape(2)+3))**(1./3.)

CALL init_dist_sbm(cloud_init(1),h_shape(1),dnc,rain_init(1),h_shape(2),dnr,diams,ffcd)

do i=1,nx
   do k=1,nz
      if (z(k)>=zctrl(2) .and. z(k)<=zctrl(3)) then
         dropsm2d(k,i,:)=ffcd
      endif
   enddo
enddo

end subroutine sbm_init

!------------------------------------------------------------
subroutine tau_init(aer2d,dropsm2d,dropsn2d)

use switches, only: zctrl
use parameters, only:nx,nz,num_h_moments,h_shape,max_nbins
use column_variables, only: z
use micro_prm

implicit none
!real(8),dimension(NQP) :: tcd ! tau composite distribution -ahu
!NQP=AERO_BIN+Ln2+LK+LK+1=71
real :: dnc,dnr
real(8),dimension(max_nbins) :: ffcd_mass,ffcd_num
integer :: i,k
real, dimension(nz,nx,max_nbins) :: aer2d,dropsm2d,dropsn2d

CALL micro_init_tau()
!Set up initial distribution and set moment values and parameter guesses
dnc=0.;dnr=0.
if(cloud_init(1)>0.)dnc = (cloud_init(1)*6./3.14159/1000. &
                          /cloud_init(2)*gamma(h_shape(1)) &
                          /gamma(h_shape(1)+3))**(1./3.)
if(rain_init(1)>0.)dnr = (rain_init(1)*6./3.14159/1000. &
                         /rain_init(2)*gamma(h_shape(2)) &
                         /gamma(h_shape(2)+3))**(1./3.)

CALL init_dist_tau(cloud_init(1),h_shape(1),dnc,rain_init(1),h_shape(2),&
                   dnr,diams,ffcd_mass,ffcd_num)

do i=1,nx
   do k=1,nz
      if (z(k)>=zctrl(2) .and. z(k)<=zctrl(3)) then
         dropsm2d(k,i,:)=ffcd_mass
         dropsn2d(k,i,:)=ffcd_num
      endif
   enddo
enddo


end subroutine tau_init
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine mp_amp(Mpc,Mpr,guessc,guessr,press,tempk,qv,fncn,ffcd_mass,&
    ffcd_num,mc,mr,flag,ffcd_massinit,ffcd_numinit)

use module_hujisbm
use micro_prm
use parameters, only: nx,nz,num_h_moments,flag_count,max_nbins
use mphys_tau_bin_declare, only: XKmean
use, intrinsic :: ieee_arithmetic, only: IEEE_Value, IEEE_QUIET_NAN
use, intrinsic :: iso_fortran_env, only: real32

implicit none
integer:: i,j,k,ip
real(8),dimension(nz,nx,2,flag_count)::flag

real, dimension(nz,nx)::tempk,press,qv
real, dimension(nz,nx,max_nbins)::ffcd_mass,ffcd_num,fncn,ffcd_massinit,&
                                  ffcd_numinit
real(8),dimension(nz,nx,num_h_moments(1)) :: Mpc
real(8),dimension(nz,nx,num_h_moments(2)) :: Mpr
real(8),dimension(nz,nx,2) :: guessc,guessr
real(8),dimension(nz,nx,10) :: mc,mr,mc0,mr0
real(8),dimension(max_nbins) :: ffcloud_mass,ffrain_mass,ffcloud_num,ffrain_num
real(8) :: dummy,realpmom,newdiam
real(real32) :: nan

nan = IEEE_VALUE(nan, IEEE_QUIET_NAN)

do k=1,nz
 do j=1,nx
   flag(k,j,:,:) = 0
   !Find gamma PDF parameters
   !searchparamsG returns the distribution bins
   !-----------CLOUD---------------------------------
   ihyd=1
   Mp(1:num_h_moments(1))=Mpc(k,j,:) !Mp is the global variable
   skr=1; ekr=krdrop
   momx=pmomsc(2) !momx is a global variable
   momy=pmomsc(3) !pmomsc(1)=3 always
!if(k==18)print*,k,Mp
   if (Mp(1)>0.) then
      CALL searchparamsG(guessc(k,j,:),ihyd,ffcloud_mass,flag(k,j,1,:))
      if (bintype .eq. 'tau') then
          ! assuming a single droplet size for all drops in a bin for now
          ! CHECK UNITS! -ahu!!!
          ffcloud_num=ffcloud_mass/XKmean
      endif
   else
      ffcloud_mass=0.
      ffcloud_num=0.
      flag(k,j,1,1)=-1
      flag(k,j,1,2:flag_count)=nan
   endif
   !----------RAIN--------------------------------
   ihyd=2
   Mp(1:num_h_moments(2))=Mpr(k,j,:)
   skr=krdrop+1; ekr=nkr
   momx = pmomsr(2)
   momy=pmomsr(3) !pmomsr(1)=3 always
   if (Mp(1)>0.) then
      CALL searchparamsG(guessr(k,j,:),ihyd,ffrain_mass,flag(k,j,2,:))
      if (bintype .eq. 'tau') then
          ffrain_num=ffrain_mass/XKmean
      endif
   else
      ffrain_mass=0.
      ffrain_num=0.
      flag(k,j,2,1)=-1
      flag(k,j,2,2:flag_count)=nan
   endif
   !----------SUM the Cloud and Rain Distributions-----------
   ffcd_mass(k,j,:) = ffcloud_mass+ffrain_mass
   ffcd_num(k,j,:) = ffcloud_num+ffrain_num
   ffcd_massinit(k,j,:)=ffcd_mass(k,j,:)
   ffcd_numinit(k,j,:)=ffcd_num(k,j,:)

if(ffcd_mass(k,j,1).ne.ffcd_mass(k,j,1))then
    print*,'NaNs in ffcd_mass'
    print*,'NaN',k,j,ffcloud_mass(1),ffrain_mass(1)
    print*,Mpc(k,j,:),Mpr(k,j,:)
    stop
endif

if(ffcd_num(k,j,1).ne.ffcd_num(k,j,1))then
    print*,'NaNs in ffcd_num'
    print*,'NaN',k,j,ffcloud_num(1),ffrain_num(1)
    stop
endif

   !Calculate moments - most of the time they're the same as what we used to find parameters
   !But in the case that parameters couldn't be found, we want to know what the actual moment
   !values are of our distributions
   call calcmom(ffcd_mass,mc0(k,j,:),mr0(k,j,:))
   ! do i=1,10
   !   mc0(k,j,i)=sum(ffcloud_mass(1:krdrop)/xl(1:krdrop)*diams(1:krdrop)**(i-1))*col*1000.
   !   mr0(k,j,i)=sum(ffrain_mass(krdrop+1:nkr)/xl(krdrop+1:nkr)*diams(krdrop+1:nkr)**(i-1))*col*1000.
   ! enddo
 enddo
enddo
!print*, Mpc(25,1,1),mc0(25,1,4),flag(25,1,1)
!print*,Mpc(25,1,2),mc0(25,1,momx+1)
!------CALL MICROPHYSICS--------------------

if (bintype .eq. 'sbm') then
    call micro_proc_sbm(press,tempk,qv,fncn,ffcd_mass)
elseif (bintype .eq. 'tau') then
    call micro_proc_tau(tempk,qv,ffcd_mass,ffcd_num)
endif

!---------CALC MOMENTS-----------------------
do k=1,nz
 do j=1,nx
     call calcmom(ffcd_mass,mc(k,j,:),mr(k,j,:))
   ! do i=1,10
   !    mc(k,j,i)=sum(ffcd_mass(k,j,1:krdrop)/xl(1:krdrop)*diams(1:krdrop)**(i-1))*col*1000.
   !    mr(k,j,i)=sum(ffcd_mass(k,j,krdrop+1:nkr)/xl(krdrop+1:nkr)*diams(krdrop+1:nkr)**(i-1))*col*1000.
   ! enddo

   !---------UPDATE MOMENTS---------------------------------
   !Mp=predicted value before microphysics. m0=value after finding parameters
   !If m0 .ne. Mp, then finding parameters failed
   !If Mp<m0, it's possible to predict too large a change to M
   !To avoid, this, scale the change by Mp/m0.
   !mc=diagnosed value after microphysics. If there are no issues
   !then m becomes the new Mp
   do i=1,npm  !Loop over the moments
      ip=pmomsc(i)+1 !pmomsc contains the predicted moments. +1 to get correct index in mc and mc0
      if (mc0(k,j,ip) > 0.) then
         dummy=(mc(k,j,ip)-mc0(k,j,ip))*Mpc(k,j,i)/mc0(k,j,ip)+Mpc(k,j,i)
         if(dummy<0.) then
           Mpc(k,j,i) = 0.
         else
           Mpc(k,j,i)=dummy
         endif
         if (i>=2) then !check new mean size
           realpmom = real(pmomsc(i))
           newdiam = (Mpc(k,j,1)/Mpc(k,j,i))**(1./(3.-realpmom))
           if (newdiam < diams(1)) then
              Mpc(k,j,i) = Mpc(k,j,1)/(1.01*diams(1))**(3.-realpmom)
           elseif (newdiam > diams(krdrop)) then
              Mpc(k,j,i) = Mpc(k,j,1)/(0.99*diams(krdrop))**(3.-realpmom)
           endif
         endif
         !Could also think about checking diameter in terms of my/mx for 3M ...
         !Would reduce acceptable space, but would not be the same as checking
         !for the actual solution space. Feels like there ought to be a way to
         !easily check for the actual solution space without resorting to look
         !up tables. Could also just use the look up table solution space here.
      else
         Mpc(k,j,i)=mc(k,j,ip)
      endif
      mc(k,j,ip)=Mpc(k,j,i)

      ip=pmomsr(i)+1
      if (mr0(k,j,ip)>0.) then
         dummy=(mr(k,j,ip)-mr0(k,j,ip))*Mpr(k,j,i)/mr0(k,j,ip)+Mpr(k,j,i)
         if(dummy<0.) then
           Mpr(k,j,i)=0.
         else
           Mpr(k,j,i)=dummy
         endif
         if (i>=2) then !check new mean size
           realpmom = real(pmomsr(i))
           newdiam = (Mpr(k,j,1)/Mpr(k,j,i))**(1./(3.-realpmom))
           if (newdiam < diams(krdrop+1)) then
              Mpr(k,j,i) = Mpr(k,j,1)/(1.01*diams(krdrop+1))**(3.-realpmom)
           elseif (newdiam > diams(nkr)) then
              Mpr(k,j,i) = Mpr(k,j,1)/(0.99*diams(nkr))**(3.-realpmom)
           endif
         endif
      else
         Mpr(k,j,i)=mr(k,j,ip)
      endif
      mr(k,j,ip)=Mpr(k,j,i)
    enddo

    if (any(mc(k,j,:)==0.)) then
      mc(k,j,:)=0;Mpc(k,j,:)=0
    endif
    if (any(mr(k,j,:)==0.)) then
      mr(k,j,:)=0;Mpr(k,j,:)=0
    endif
 enddo
enddo
end subroutine mp_amp
!---------------------------------------------------------------------
subroutine mp_sbm(ffcd,press,tempk,qv,fncn,mc,mr)

use module_hujisbm
use micro_prm
use parameters, only: nx,nz,num_h_moments,max_nbins

implicit none
integer:: i,j,k,ip

real, dimension(nz,nx)::tempk,press,qv
real, dimension(nz,nx,max_nbins)::ffcd,fncn
real(8),dimension(nz,nx,10) :: mc,mr

!------CALL MICROPHYSICS--------------------
call micro_proc_sbm(press,tempk,qv,fncn,ffcd)

!---------CALC MOMENTS-----------------------
do k=1,nz
 do j=1,nx
   do i=1,10
      mc(k,j,i)=sum(ffcd(k,j,1:krdrop)/xl(1:krdrop)*diams(1:krdrop)**(i-1))*col*1000.
      mr(k,j,i)=sum(ffcd(k,j,krdrop+1:nkr)/xl(krdrop+1:nkr)*diams(krdrop+1:nkr)**(i-1))*col*1000.
   enddo
 enddo
enddo
end subroutine mp_sbm

!---------------------------------------------------------------------
subroutine mp_tau(ffcd_mass2d,ffcd_num2d,tempk,qv,mc,mr)

use module_hujisbm
use micro_prm
use parameters, only: nx,nz,num_h_moments,max_nbins
use namelists, only: l_advect
use mphys_tau_bin_declare, only: lk_cloud,x_bin
implicit none
integer:: i,j,k,ip

real, dimension(nz,nx):: tempk,qv
real, dimension(nz,nx,max_nbins)::ffcd_mass2d,ffcd_num2d
real(8),dimension(nz,nx,10) :: mc,mr ! moments

!------CALL MICROPHYSICS--------------------
call micro_proc_tau(tempk,qv,ffcd_mass2d,ffcd_num2d)

!---------CALC MOMENTS-----------------------
do k=1,nz
 do j=1,nx
   do i=1,10
      mc(k,j,i)=sum(ffcd_mass2d(k,j,1:lk_cloud)/x_bin(2:lk_cloud+1)*diams(1:lk_cloud)**(i-1))*col*1.e6
      mr(k,j,i)=sum(ffcd_mass2d(k,j,lk_cloud+1:nkr)/x_bin(lk_cloud+2:nkr+1)*diams(lk_cloud+1:nkr)**(i-1))*col*1.e6
   enddo
 enddo
enddo
! there might be a better way to calculate moments, but will leave it like that for now. -ahu!!!

end subroutine mp_tau

subroutine calcmom(ffcd_mass,mc,mr)

use micro_prm
use mphys_tau_bin_declare, only: lk_cloud,x_bin
Use namelists, only: bintype
implicit none
integer :: i,j,k,ip
real, dimension(max_nbins)::ffcd_mass
real(8),dimension(10) :: mc,mr ! moments

do i=1,10
    if (bintype .eq. 'sbm') then
        mc=sum(ffcd_mass(1:krdrop)/xl(1:krdrop)*diams(1:krdrop)**(i-1))*col*1000.
        mr=sum(ffcd_mass(krdrop+1:nkr)/xl(krdrop+1:nkr)*diams(krdrop+1:nkr)**(i-1))*col*1000.
    elseif (bintype .eq. 'tau') then
        mc=sum(ffcd_mass(1:lk_cloud)/x_bin(2:lk_cloud+1)*diams(1:lk_cloud)**(i-1))*col*1.e6
        mr=sum(ffcd_mass(lk_cloud+1:nkr)/x_bin(lk_cloud+2:nkr+1)*diams(lk_cloud+1:nkr)**(i-1))*col*1.e6
    endif
enddo

end subroutine calcmom
