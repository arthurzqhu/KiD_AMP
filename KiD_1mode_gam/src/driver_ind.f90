subroutine amp_init(aer2d,Mpc2d,Mpr2d,guessc2d,guessr2d)

use switches, only: zctrl
use parameters, only:nx,nz,num_h_moments,h_shape,max_nbins,split_bins
use column_variables, only: z
use module_hujisbm
use micro_prm
use mphys_tau_bin_declare, only: lk_cloud,xk,xkgmean,dgmean
use namelists, only: bintype, l_truncated, l_use_nn
use global_fun, only: get_meandiam
use physconst, only: pi
use mphys_tau_bin_declare, only: diam
! use module_mp_boss

implicit none
character(4)::momstr
character(2)::rdrop_binstr
character(len=100) :: lutfolder
real(8), dimension(nz,nx,num_h_moments(1)) :: Mpc2d
real(8), dimension(nz,nx,num_h_moments(2)) :: Mpr2d
real(8), dimension(nz,nx,2) :: guessc2d,guessr2d
real(8),dimension(max_nbins) :: ffcd_mass, ffcd_num
real(8),dimension(nz,nx,max_nbins) :: dropsm2d, dropsn2d
real(8),dimension(num_h_moments(1)) :: mc
real(8),dimension(num_h_moments(2)) :: mr
real(8) :: z_cbi,z_cti,d_cloudi
real(8) :: dnc,dnr
real(8), dimension(nz,nx,max_nbins) :: aer2d
real(8), dimension(max_nbins) :: diag_m, diag_D !diagnosed mass and diam of each bin
integer :: i,k,imom,ip,ib
double precision :: infinity, ratio_fix
integer, parameter :: max_rows = 120  ! Maximum number of rows in the CSV file
integer, parameter :: max_cols = 34   ! Maximum number of columns in the CSV file
integer :: row, colm, iostat
character(1) :: delimiter = ','       ! The delimiter used in the CSV file
character(100) :: filename, filenamen ! The name of the CSV file
integer :: num_rows = 0               ! The number of rows in the CSV file
integer :: num_cols = 0               ! The number of columns in the CSV file

infinity = HUGE(dnr)

if (l_truncated) then
   tORf = 'truncated'
else
   tORf = 'fullgam'
endif

!Set some parameters
ndtcoll=1
ipris=0; ihail=0; igraup=0;iceprocs=0;imbudget=1
z_cbi=zctrl(2)
z_cti=zctrl(3)
d_cloudi=z_cti-z_cbi

npmc=num_h_moments(1)
npmr=num_h_moments(2)
npm = maxval(num_h_moments)
dnc = guessc2d(1,1,2)
dnr = guessr2d(1,1,2)

if (npmc < 4) then
   pmomsc(1:3)=(/3,imomc1,imomc2/)
   pmomsr(1:3)=(/3,imomr1,imomr2/)
else
   pmomsc(1:6)=(/3, 0, imomc1, imomc2, imomr1, imomr2/)
   pmomsr(1:6)=(/3, 0, imomc1, imomc2, imomr1, imomr2/)
endif

! 3m initialization: {{{
!Open output files, or decide that we don't need to run this combination of moments
if (npmc < 4 .and. imomc1>=7 .and. imomc1 .ne. imomc2) then
 print*,"not a valid moment combination",imomc1,imomc2
 stop
else !Run AMP
  write(rdrop_binstr,'(I0)') split_bins
  !Open and read lookup tables
  lutfolder='./src/input_data/'//trim(bintype)//'_lutables_'//rdrop_binstr//'/'
endif

! if (npmc >= 4) then
!   if (bintype .eq. 'tau') then
!      lutfolder='./src/input_data/'//trim(bintype)//'_lutables_34/'
!   elseif (bintype .eq. 'sbm') then
!      lutfolder='./src/input_data/'//trim(bintype)//'_lutables_33/'
!   endif
! endif

if (npmc==3) then
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
   dnc=dntab(50,50,1)
endif

if (npmr==3) then
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
   dnr=dntab(50,50,2)
endif

if (npmr>=3 .or. npmc>=3) then
   dntab = dntab * 1.e-6
   !Find min and max of parameters from the tables
   dnbounds(1,1)=minval(dntab(:,:,1)); dnbounds(2,1)=maxval(dntab(:,:,1))
   nubounds(1,1)=minval(nutab(:,:,1)); nubounds(2,1)=maxval(nutab(:,:,1))
   minmaxmy(1,1)=minval(mintab(:,1)); minmaxmy(2,1)=maxval(maxtab(:,1))
   dnbounds(1,2)=minval(dntab(:,:,2)); dnbounds(2,2)=maxval(dntab(:,:,2))
   nubounds(1,2)=minval(nutab(:,:,2)); nubounds(2,2)=maxval(nutab(:,:,2))
   minmaxmy(1,2)=minval(mintab(:,2)); minmaxmy(2,2)=maxval(maxtab(:,2))
endif

! }}}

!call microphysics initialization routines
if (bintype .eq. 'sbm') then
    CALL micro_init_sbm()
    CALL micro_init_sbm2()
    CALL kernalsdt()
elseif (bintype .eq. 'tau') then
    CALL micro_init_tau()
endif

if (bintype .eq. 'tau') then
  D_min = diam(1)/100
  D_max = diam(nkr+1)/100
  ! D_min = diams(1)
  ! D_max = diams(nkr)
  ! print*, D_min, D_max
  ! stop
endif

do imom = 1,10
  do ib = 1,nkr
    pdiams(imom, ib) = diams(ib)**(dble(imom)-1.)
  enddo
enddo

if (initprof .eq. 'c') then
   if(cloud_init(1)>0.) then
      dnc = (cloud_init(1)*QtoM3/cloud_init(2)* & 
             gamma(h_shape(1))/gamma(h_shape(1)+3))**(1./3.)
   endif
   if(rain_init(1)>0.) then
      dnr = (rain_init(1)*QtoM3/rain_init(2)* &
             gamma(h_shape(2))/gamma(h_shape(2)+3))**(1./3.)
   endif
   if(cloud_init(1)>0. .or. rain_init(1)>0.) then
       if (bintype .eq. 'sbm') then
           CALL init_dist_sbm(cloud_init(1),h_shape(1),dnc,rain_init(1),&
                              h_shape(2),dnr,diams,ffcd_mass)
       elseif (bintype .eq. 'tau') then
           CALL init_dist_tau(cloud_init(1),h_shape(1),dnc,rain_init(1),&
                              h_shape(2),dnr,ffcd_mass,ffcd_num)
       endif
   else
      ffcd_mass(:)=0.
      if (bintype .eq. 'tau') ffcd_num(:)=0.
   endif
   guessc2d(:,:,2)=dnc
   guessr2d(:,:,2)=dnr
   ratio_fix = (cloud_init(1)+rain_init(1))/sum(ffcd_mass*col)
   ffcd_mass = ffcd_mass*ratio_fix
   ratio_fix = (cloud_init(2)+rain_init(2))/sum(ffcd_num*col)
   ffcd_num = ffcd_num*ratio_fix
   ! stop
   do i=1,num_h_moments(1)
     if (bintype .eq. 'sbm') then
       mc(i)=sum(ffcd_mass(1:split_bins)/xl(1:split_bins)*diams(1:split_bins)**pmomsc(i))*col*1000.
     elseif (bintype .eq. 'tau') then
     select case (pmomsc(i))
     case (0)
       mc(i)=sum(ffcd_num(1:split_bins))*col
     case default
       mc(i)=sum(ffcd_mass(1:split_bins)/binmass(1:split_bins)*diams(1:split_bins)**pmomsc(i))*col
     end select 
     endif
   enddo
   do i=1,num_h_moments(2)
     if (bintype .eq. 'sbm') then
       mr(i)=sum(ffcd_mass(split_bins+1:nkr)/xl(split_bins+1:nkr)*diams(split_bins+1:nkr)**pmomsr(i))*col*1000.
     elseif (bintype .eq. 'tau') then
     select case (pmomsr(i))
     case (0)
       mr(i)=sum(ffcd_num(split_bins+1:nkr))*col
     case default
       mr(i)=sum(ffcd_mass(split_bins+1:nkr)/binmass(split_bins+1:nkr)*diams(split_bins+1:nkr)**pmomsc(i))*col
     end select 
     end if
   enddo
endif

! print*, mc+mr
! print*, dnc, dnr
! print*, cloud_init(1)*QtoM3/(dnc**(h_shape(1)+3)*gamnu1_3)*dnc**(h_shape(1)+pmomsc(3))*gamnu1_w
! ! print*, cloud_init(2)
! stop

do i=1,nx
   do k=1,nz
      if (z(k)>=zctrl(2) .and. z(k)<=zctrl(3)) then
         if (initprof .eq. 'i') then
            if (cloud_init(1)>0.) then
               dnc = (cloud_init(1)*(z(k)-z_cbi)/d_cloudi*6./3.14159/1000./&
                  cloud_init(2)*gamma(h_shape(1))/gamma(h_shape(1)+3))**(1./3.)
               if (gamma(h_shape(1)+3) > infinity) dnc = 0.
            endif
            if (rain_init(1)>0.) then 
               dnr = (rain_init(1)*(z(k)-z_cbi)/d_cloudi*6./3.14159/1000./&
                  rain_init(2)*gamma(h_shape(2))/gamma(h_shape(2)+3))**(1./3.)
               if (gamma(h_shape(2)+3) > infinity) dnr = 0.
            endif

            if (cloud_init(1)>0. .or. rain_init(1)>0.) then
                if (bintype .eq. 'sbm') then
                    CALL init_dist_sbm(cloud_init(1)*(z(k)-z_cbi)/d_cloudi,&
                       h_shape(1),dnc,rain_init(1)*(z(k)-z_cbi)/d_cloudi,&
                       h_shape(2),dnr,diams,ffcd_mass)
                elseif (bintype .eq. 'tau') then
                    CALL init_dist_tau(cloud_init(1)*(z(k)-z_cbi)/d_cloudi,&
                       h_shape(1),dnc,rain_init(1)*(z(k)-z_cbi)/d_cloudi,&
                       h_shape(2),dnr,ffcd_mass,ffcd_num)
                endif
            else
               ffcd_mass(:)=0.
               if (bintype .eq. 'tau') ffcd_num(:)=0.
            endif

            do imom=1,num_h_moments(1)
              if (bintype .eq. 'sbm') then
                mc(imom)=sum(ffcd_mass(1:split_bins)/xl(1:split_bins)*&
                   diams(1:split_bins)**pmomsc(imom))*col*1000.
              elseif (bintype .eq. 'tau') then
                mc(imom)=sum(ffcd_mass(1:split_bins)/binmass(1:split_bins)*&
                   diams(1:split_bins)**pmomsc(imom))*col
              endif
            enddo

            do imom=1,num_h_moments(2)
              if (bintype .eq. 'sbm') then
                mr(imom)=sum(ffcd_mass(split_bins+1:nkr)/xl(split_bins+1:nkr)*&
                   diams(split_bins+1:nkr)**pmomsr(imom))*col*1000.
              elseif (bintype .eq. 'tau') then
                mr(imom)=sum(ffcd_mass(split_bins+1:nkr)/binmass(split_bins+1:nkr)*&
                   diams(split_bins+1:nkr)**pmomsr(imom))*col
              end if
            enddo

! little tricky here so this implementation is not definitive - different moments increase at 
! different rates. what's done here is that the increase rate is proportional to the order of
! the moment, given that 3rd moment (mass) increases linearly and 0th moment (number) is 
! constant. -ahu
!            Mpc2d(k,i,1)=mc(1)*(z(k)-z_cbi)/d_cloudi
!            Mpr2d(k,i,1)=mr(1)*(z(k)-z_cbi)/d_cloudi
!
!            do imom=2,npmc
!               ip=pmomsc(imom)+1
!               Mpc2d(k,i,imom)=mc(imom)*((z(k)-z_cbi)/d_cloudi)**(pmomsc(imom)/3.)
!            enddo
!
!            do imom=2,npmr
!               ip=pmomsr(imom)+1
!               Mpr2d(k,i,imom)=mr(imom)*((z(k)-z_cbi)/d_cloudi)**(pmomsr(imom)/3.)
!            enddo
!
!            print*, k,Mpc2d(k,i,1)
         endif

         ! update Mpc2d
         if (npmc < 4) then
            Mpc2d(k,i,1:num_h_moments(1)) = mc(1:num_h_moments(1))
            Mpr2d(k,i,1:num_h_moments(2)) = mr(1:num_h_moments(2))
         else
            Mpc2d(k,i,1:num_h_moments(1)) = mc(1:num_h_moments(1)) + mr(1:num_h_moments(2))
         endif
      endif
   enddo
enddo

if (extralayer) then
   dnr = (2d-4*6./3.14159/1000./1.e5* & 
          gamma(h_shape(2))/gamma(h_shape(2)+3))**(1./3.)
   if (bintype .eq. 'sbm') then
       CALL init_dist_sbm(0d0,h_shape(1),dnc,2d-4,&
                          h_shape(2),dnr,diams,ffcd_mass)
   elseif (bintype .eq. 'tau') then
       CALL init_dist_tau(0d0,h_shape(1),dnc,2d-4,&
                          h_shape(2),dnr,ffcd_mass,ffcd_num)
   endif
   do i=1,num_h_moments(1)
     if (bintype .eq. 'sbm') then
       mc(i)=sum(ffcd_mass(1:split_bins)/xl(1:split_bins)*diams(1:split_bins)**pmomsc(i))*col*1000.
     elseif (bintype .eq. 'tau') then
       mc(i)=sum(ffcd_mass(1:split_bins)/binmass(1:split_bins)*diams(1:split_bins)**pmomsc(i))*col
     endif
   enddo
   do i=1,num_h_moments(2)
     if (bintype .eq. 'sbm') then
       mr(i)=sum(ffcd_mass(split_bins+1:nkr)/xl(split_bins+1:nkr)*diams(split_bins+1:nkr)**pmomsr(i))*col*1000.
     elseif (bintype .eq. 'tau') then
       mr(i)=sum(ffcd_mass(split_bins+1:nkr)/binmass(split_bins+1:nkr)*diams(split_bins+1:nkr)**pmomsc(i))*col
     end if
   enddo
   do k = 10,37
      guessr2d(k,:,2)=dnr
      if (npmc < 4) then
         Mpc2d(k,1,1:num_h_moments(1)) = mc(1:num_h_moments(1))
         Mpr2d(k,1,1:num_h_moments(2)) = mr(1:num_h_moments(2))
      else
         Mpc2d(k,1,1:num_h_moments(1)) = mc(1:num_h_moments(1))+mr(1:num_h_moments(2))
      endif
   enddo
endif

! print*, mc+mr
! stop



if (l_hist_run) then
   ! Open the CSV file
   filename='./hist_run_prof/sbm_a200w8_400s.txt'
   ! if (bintype .eq. 'tau') filenamen='./hist_run_prof/'//trim(bintype)//'n.txt'
   open(10, file=filename, status='old', action='read', iostat=iostat)
   if (iostat /= 0) then
      print *, "Error opening file:", filename
      stop
   endif

   ! if (bintype .eq. 'tau') then
   !    open(11, file=filenamen, status='old', action='read', iostat=iostat)
   !    if (iostat /= 0) then
   !       print *, "Error opening file:", filenamen
   !       stop
   !    endif
   ! endif

   ! Read the data from the CSV file
   do while (.true.)
      read(10, *, iostat=iostat) dropsm2d(num_rows+1,1,:)
      ! if (bintype .eq. 'tau') read(11, *, iostat=iostat) dropsn2d(num_rows+1,1,:)
      if (iostat /= 0) exit ! Exit the loop when there is no more data to read
      num_rows = num_rows + 1
      num_cols = size(dropsm2d, 3)
      if (num_rows > max_rows) then
         print *, "Error: Maximum number of rows exceeded"
         stop
      endif
      if (num_cols > max_cols) then
         print *, "Error: Maximum number of columns exceeded"
         stop
      endif
   enddo

   ! Close the CSV file
   close(10)
   ! if (bintype .eq. 'tau') close(11)
   do k=1,nz
      do i=1,nx
         if (bintype .eq. 'tau') dropsn2d(k,i,:) = dropsm2d(k,i,:)/binmass
      enddo
   enddo

   do i=1,nx
      do k=1,nz

         if (bintype .eq. 'tau') then
            diag_m=binmass
            diag_D=(diag_m*QtoM3)**(1./3.)
            do ib=1,nkr
               if ((diag_D(ib) .ne. diag_D(ib)) .or. (diag_D(ib)>infinity)) diag_D(ib)=diams(ib)
            end do
         endif

         do imom=1,num_h_moments(1)
            if (bintype .eq. 'sbm') then
               mc(imom)=sum(dropsm2d(k,i,1:split_bins)/xl(1:split_bins)*&
                  diams(1:split_bins)**pmomsc(imom))*col*1000.
            elseif (bintype .eq. 'tau') then
               mc(imom)=sum(dropsn2d(k,i,1:split_bins)*diag_D(1:split_bins)**pmomsc(imom))*col
            endif
         enddo
         do imom=1,num_h_moments(2)
            if (bintype .eq. 'sbm') then
               mr(imom)=sum(dropsm2d(k,i,split_bins+1:nkr)/xl(split_bins+1:nkr)*&
                  diams(split_bins+1:nkr)**pmomsc(imom))*col*1000.
            elseif (bintype .eq. 'tau') then
               mr(imom)=sum(dropsn2d(k,i,split_bins+1:nkr)*diag_D(split_bins+1:nkr)**pmomsc(imom))*col
            endif
         enddo

         if (npmc < 4) then
            Mpc2d(k,i,1:num_h_moments(1)) = mc(1:num_h_moments(1))
            Mpr2d(k,i,1:num_h_moments(2)) = mr(1:num_h_moments(2))
         else
            Mpc2d(k,i,1:num_h_moments(1)) = mc(1:num_h_moments(1)) + mr(1:num_h_moments(2))
         endif

      enddo
   enddo
endif ! l_hist_run

if (l_use_nn) &
  call moment2state_net % load("/Users/arthurhu/Downloads/moments_to_state_fixed_mu_nn.txt")

! momx = imomc1
! momy = imomc2
! if (npm>=4) call read_boss_slc_param

end subroutine amp_init

! sbm_init: {{{
!------------------------------------------------------------
subroutine sbm_init(aer2d,dropsm2d)

use switches, only: zctrl
use parameters, only:nx,nz,num_h_moments,h_shape,max_nbins
use column_variables, only: z
use module_hujisbm
use micro_prm

implicit none
real(8),dimension(max_nbins) :: ffcd
real(8) :: dnc,dnr
real(8), dimension(nz,nx,max_nbins) :: aer2d,dropsm2d
real(8) :: z_cbi,z_cti,d_cloudi
integer :: i,k
real(8) :: infinity
integer, parameter :: max_rows = 120  ! Maximum number of rows in the CSV file
integer, parameter :: max_cols = 34   ! Maximum number of columns in the CSV file
integer :: row, colm, iostat
character(1) :: delimiter = ','       ! The delimiter used in the CSV file
character(100) :: filename            ! The name of the CSV file
integer :: num_rows = 0               ! The number of rows in the CSV file
integer :: num_cols = 0               ! The number of columns in the CSV file

infinity = huge(dnc)

!Set some parameters
ndtcoll=1
ipris=0; ihail=0; igraup=0;iceprocs=0;imbudget=1
z_cbi=zctrl(2)
z_cti=zctrl(3)
d_cloudi=z_cti-z_cbi

!call microphysics initialization routines
CALL micro_init_sbm()
CALL micro_init_sbm2()
CALL kernalsdt()

!Set up initial distribution and set moment values and parameter guesses
dnc=0.;dnr=0.
if (initprof .eq. 'c') then
   if(cloud_init(1)>0.) dnc = (cloud_init(1)*6./3.14159/1000. &
                              /cloud_init(2)*gamma(h_shape(1)) &
                              /gamma(h_shape(1)+3))**(1./3.)
   if(rain_init(1)>0.) dnr = (rain_init(1)*6./3.14159/1000. &
                             /rain_init(2)*gamma(h_shape(2)) &
                             /gamma(h_shape(2)+3))**(1./3.)
   if (gamma(h_shape(1)+3) > infinity) dnc = 0.
   if (gamma(h_shape(2)+3) > infinity) dnr = 0.
   CALL init_dist_sbm(cloud_init(1),h_shape(1),dnc,rain_init(1),h_shape(2),dnr,diams,ffcd)
endif

if (l_hist_run) then
   ! Open the CSV file
   filename='./hist_run_prof/sbm_a200w8_400s.txt'
   open(10, file=filename, status='old', action='read', iostat=iostat)
   if (iostat /= 0) then
      print *, "Error opening file:", filename
      stop
   endif

   ! Read the data from the CSV file
   do while (.true.)
      read(10, *, iostat=iostat) dropsm2d(num_rows+1,1,:)
      if (iostat /= 0) exit ! Exit the loop when there is no more data to read
      num_rows = num_rows + 1
      num_cols = size(dropsm2d, 3)
      if (num_rows > max_rows) then
         print *, "Error: Maximum number of rows exceeded"
         stop
      endif
      if (num_cols > max_cols) then
         print *, "Error: Maximum number of columns exceeded"
         stop
      endif
   enddo

   ! Close the CSV file
   close(10)
endif

if (extralayer) then
   dnr = (2d-4*6./3.14159/1000./1.e5* & 
          gamma(h_shape(2))/gamma(h_shape(2)+3))**(1./3.)
   CALL init_dist_sbm(0d0,h_shape(1),dnc,2d-4,&
                      h_shape(2),dnr,diams,ffcd)
   do k=10,37
      dropsm2d(k,1,:) = ffcd
   enddo
endif

do i=1,nx
   do k=1,nz
      if (z(k)>=zctrl(2) .and. z(k)<=zctrl(3)) then
         if (initprof .eq. 'c') then
            dropsm2d(k,i,:)=ffcd
         elseif (initprof .eq. 'i') then
            !calculate dnc/r for each layer with water
            if(cloud_init(1)>0.) dnc = (cloud_init(1)*(z(k)-z_cbi)/d_cloudi*6./3.14159/1000. &
                                       /cloud_init(2)*gamma(h_shape(1)) &
                                       /gamma(h_shape(1)+3))**(1./3.)
            if(rain_init(1)>0.) dnr = (rain_init(1)*(z(k)-z_cbi)/d_cloudi*6./3.14159/1000. &
                                      /rain_init(2)*gamma(h_shape(2)) &
                                      /gamma(h_shape(2)+3))**(1./3.)
            if (gamma(h_shape(1)+3) > infinity) dnc = 0.
            if (gamma(h_shape(2)+3) > infinity) dnr = 0.
            CALL init_dist_sbm(cloud_init(1)*(z(k)-z_cbi)/d_cloudi,h_shape(1),&
               dnc,rain_init(1)*(z(k)-z_cbi)/d_cloudi,h_shape(2),dnr,diams,ffcd)

            dropsm2d(k,i,:)=ffcd
            !dropsm2d(k,i,:)=ffcd*(z(k)-z_cbi)/d_cloudi
         endif
      endif
   enddo
enddo

end subroutine sbm_init
! }}}

! tau_init: {{{
!------------------------------------------------------------
subroutine tau_init(aer2d,dropsm2d,dropsn2d)

use switches, only: zctrl
use parameters, only:nx,nz,num_h_moments,h_shape,max_nbins
use column_variables, only: z
use micro_prm
use mphys_tau_bin_declare, only: xk, x_bin, xkgmean
! use module_mp_boss
use namelists, only: moments_diag, nmom_diag

implicit none
!real(8),dimension(NQP) :: tcd ! tau composite distribution -ahu
!NQP=AERO_BIN+Ln2+LK+LK+1=71
real(8) :: dnc,dnr
real(8),dimension(max_nbins) :: ffcd_mass,ffcd_num
integer :: i,k
real(8), dimension(nz,nx,max_nbins) :: aer2d,dropsm2d,dropsn2d
real(8) :: z_cbi,z_cti,d_cloudi
real(8) :: infinity
integer, parameter :: max_rows = 120  ! Maximum number of rows in the CSV file
integer, parameter :: max_cols = 34   ! Maximum number of columns in the CSV file
integer :: row, colm, iostat
character(1) :: delimiter = ','       ! The delimiter used in the CSV file
character(100) :: filename,filenamen  ! The name of the CSV file
integer :: num_rows = 0               ! The number of rows in the CSV file
integer :: num_cols = 0               ! The number of columns in the CSV file
integer :: ib


infinity = huge(dnc)

z_cbi=zctrl(2)
z_cti=zctrl(3)
d_cloudi=z_cti-z_cbi

CALL micro_init_tau()
!Set up initial distribution and set moment values and parameter guesses
dnc=0.;dnr=0.
if (initprof .eq. 'c')  then
   if(cloud_init(1)>0.)dnc = (cloud_init(1)*6./3.14159/1000. &
                             /cloud_init(2)*gamma(h_shape(1)) &
                             /gamma(h_shape(1)+3))**(1./3.)
   if(rain_init(1)>0.)dnr = (rain_init(1)*6./3.14159/1000. &
                            /rain_init(2)*gamma(h_shape(2)) &
                            /gamma(h_shape(2)+3))**(1./3.)
   if (gamma(h_shape(1)+3) > infinity) dnc = 0.
   if (gamma(h_shape(2)+3) > infinity) dnr = 0.
   CALL init_dist_tau(cloud_init(1),h_shape(1),dnc,rain_init(1),h_shape(2),&
                      dnr,ffcd_mass,ffcd_num)
endif

if (l_hist_run) then
   ! Open the CSV file
   filename='./hist_run_prof/sbm_a200w8_400s.txt'
   ! filenamen='./hist_run_prof/'//trim(bintype)//'n.txt'
   open(10, file=filename, status='old', action='read', iostat=iostat)
   if (iostat /= 0) then
      print *, "Error opening file:", filename
      stop
   endif

   ! open(11, file=filenamen, status='old', action='read', iostat=iostat)
   ! if (iostat /= 0) then
   !    print *, "Error opening file:", filenamen
   !    stop
   ! endif

   ! Read the data from the CSV file
   do while (.true.)
      read(10, *, iostat=iostat) dropsm2d(num_rows+1,1,:)
      ! read(11, *, iostat=iostat) dropsn2d(num_rows+1,1,:)
      if (iostat /= 0) exit ! Exit the loop when there is no more data to read
      num_rows = num_rows + 1
      num_cols = size(dropsm2d, 3)
      if (num_rows > max_rows) then
         print *, "Error: Maximum number of rows exceeded"
         stop
      endif
      if (num_cols > max_cols) then
         print *, "Error: Maximum number of columns exceeded"
         stop
      endif
   enddo

   ! Close the CSV file
   close(10)
endif

do k=1,nz
   do i=1,nx
      if (bintype .eq. 'tau') dropsn2d(k,i,:) = dropsm2d(k,i,:)/binmass
   enddo
enddo

do i=1,nx
   do k=1,nz
      if (z(k)>=zctrl(2) .and. z(k)<=zctrl(3)) then
         if (initprof .eq. 'i') then
            if(cloud_init(1)>0.)dnc = (cloud_init(1)*(z(k)-z_cbi)/d_cloudi*6./3.14159/1000. &
                                      /cloud_init(2)*gamma(h_shape(1)) &
                                      /gamma(h_shape(1)+3))**(1./3.)
            if(rain_init(1)>0.)dnr = (rain_init(1)*(z(k)-z_cbi)/d_cloudi*6./3.14159/1000. &
                                     /rain_init(2)*gamma(h_shape(2)) &
                                     /gamma(h_shape(2)+3))**(1./3.)
            if (gamma(h_shape(1)+3) > infinity) dnc = 0.
            if (gamma(h_shape(2)+3) > infinity) dnr = 0.
            CALL init_dist_tau(cloud_init(1)*(z(k)-z_cbi)/d_cloudi,&
               h_shape(1),dnc,rain_init(1)*(z(k)-z_cbi)/d_cloudi,&
               h_shape(2),dnr,ffcd_mass,ffcd_num)
         endif

         dropsm2d(k,i,:)=ffcd_mass
         dropsn2d(k,i,:)=ffcd_num
      endif
   enddo
enddo

if (extralayer) then
   dnr = (2d-4*6./3.14159/1000./1.e5* & 
          gamma(h_shape(2))/gamma(h_shape(2)+3))**(1./3.)
   CALL init_dist_tau(0d0, h_shape(1), dnc, 2d-4, h_shape(2), dnr, ffcd_mass, ffcd_num)

   ! ib = 31
   ! ffcd_mass(:) = 0.
   ! ffcd_num(:) = 0.
   ! ffcd_num(ib) = 1d5/col
   ! ffcd_mass(ib) = xkgmean(ib)*ffcd_num(ib)

   ! do k=50,50
   do k=10,37
      dropsm2d(k,1,:) = ffcd_mass
      dropsn2d(k,1,:) = ffcd_num
   enddo
endif

! momx = imomc1
! momy = imomc2
! call read_boss_slc_param

end subroutine tau_init
! }}}

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine mp_amp(Mpc,Mpr,guessc,guessr,press,tempk,qv,fncnr8,ffcdr8_mass,&
    ffcdr8_num,mc,mr,flag,ffcdr8_massinit,ffcdr8_numinit,ffcdr8_massf,ffcdr8_numf)

use module_hujisbm
use micro_prm
use parameters, only: nx,nz,num_h_moments,flag_count,max_nbins
use mphys_tau_bin_declare, only: xk,lk_cloud,xkgmean,xkk1
use physconst, only: pi
use diagnostics, only: i_dgtime
use, intrinsic :: ieee_arithmetic, only: IEEE_Value, IEEE_QUIET_NAN
use, intrinsic :: iso_fortran_env, only: real32
use global_fun
use column_variables, only: hydrometeors
use module_mp_amp, only: invert_moments, invert_moments_nn
use namelists, only: l_truncated, l_init_test, l_use_nn, n_cat

implicit none
integer:: i,j,k,ip
real(8),dimension(nz,nx,2,flag_count)::flag,flag_dummy

real, dimension(nz,nx)::tempk,press,qv,mc_temp,mr_temp
real(8), dimension(nz,nx,max_nbins)::ffcdr8_mass,ffcdr8_num,fncnr8,ffcdr8_massinit,&
                                  ffcdr8_numinit, ffcdr8_massf, ffcdr8_numf
real, dimension(nz,nx,max_nbins)::ffcd_mass,ffcd_num,fncn
real(8),dimension(nz,nx,num_h_moments(1)) :: Mpc
real(8),dimension(nz,nx,num_h_moments(2)) :: Mpr
real(8),dimension(nz,nx,2) :: guessc,guessr
real(8),dimension(nz,nx,1:10) :: mc,mr,mc0,mr0
real(8),dimension(max_nbins) :: ffcloud_mass,ffrain_mass,ffcloud_num,ffrain_num
real(8) :: tnum_gam, true_num, correct_ratio, true_mass, tmass_gam, tm1_gam, true_m1
real(8) :: dummy,realpmom,newdiam,Dmin,Dmax,Dskip
real(real32) :: nan
real(8) :: dn1_arr(100), dn2_arr(100), powers(100)
real(8) :: mom_pred(npm, 2), gam_param(2,2), sqerr

diag_dt1 = 0.
diag_dt2 = 0.
! diag_dt3 = 0.
itries = 0

nan = IEEE_VALUE(nan, IEEE_QUIET_NAN)

call CPU_TIME(diag_dt1)
do k=1,nz
 do j=1,nx

   flag(k,j,:,:) = 0.
   !Find gamma PDF parameters

   if (bintype .eq. 'sbm') then
      Dmin = diams(1) 
      Dskip = diams(1)*1.0001
      Dmax = diams(nkr)
   endif
   if (bintype .eq. 'tau') then 
      Dmin = xkk1(1)*2
      Dskip = dgmean(1)
      Dmax = dgmean(nkr)
   endif

   if (l_init_test) then
     mpc(k,j,1:4) = (/5.1056581531715761E-006, 362680242.18796563, 1.5132587792741748E-010, 1.7143553902324072E-019/)
     ! Mpc(k,j,1:2) = (/1e-007, 1e8/)
     ! Mpr(k,j,1:2) = (/1.4000821838411837E-007,   70.996367122565701/)
     ! Mpc(k,j,1:4) = (/1.1958954365427025E-005, 12326993.134452719, 1.5975858216584202E-008, 2.3280344892313758E-011/)
     ! Mpc(k,j,1:4) = (/1.1271152298852096E-005, 30032998.830784883, 1.4230874205695079E-008, 1.9600799200941356E-011/)
     ! Mpc(k,j,1:4) = (/5.9659897103853119E-003, 8929380865.9754963, 8.2108958909978291E-006, 1.2169634741639012E-008/)
   endif

   mom_pred(:,1) = Mpc(k,j,:)
   if (npmc < 4) then
      mom_pred(:,2) = Mpr(k,j,:)
   endif

   gam_param(1,1) = guessc(k,j,2)
   gam_param(2,1) = guessc(k,j,1)
   gam_param(1,2) = guessr(k,j,2)
   gam_param(2,2) = guessr(k,j,1)

   if (sum(mom_pred(1,:)) <= total_m3_th) then
      ffcd_mass(k,j,:) = 0.
      ffcd_num(k,j,:) = 0.
      flag(k,j,:,1) = -1
      flag(k,j,:,2:flag_count) = nan
   else
     if (l_init_test) print*, 'mom_pred', mom_pred(:,1)+mom_pred(:,2)

     if (l_use_nn .and. npm>=4) then
       call invert_moments_nn(mom_pred, gam_param, ffcd_mass(k,j,:), ffcd_num(k,j,:))
     else
       ! print*, k, Mpc(k,1,1:2)
       call invert_moments(mom_pred, gam_param, ffcd_mass(k,j,:), ffcd_num(k,j,:), flag(k,j,:,1), sqerr)
     endif

     if (ffcd_mass(k,j,1).ne.ffcd_mass(k,j,1)) then
     print*, k, mom_pred
       print*, 'ffcd nan'
       print*, 'gam_param', gam_param
       ! print*, (ffcd_mass(k,j,:))
       stop
     endif

     if (l_init_test) print*, 'cloud m3', sum(ffcd_mass(k,j,1:split_bins))*col*QtoM3
     if (l_init_test) print*, 'rain m3', sum(ffcd_mass(k,j,split_bins+1:nkr))*col*QtoM3
     if (l_init_test) print*, 'cloud m0', sum(ffcd_num(k,j,1:split_bins))*col
     if (l_init_test) print*, 'rain m0', sum(ffcd_num(k,j,split_bins+1:nkr))*col
     if (l_init_test .and. (.not. l_use_nn)) print*, 'gam_param after', gam_param
     ! if (l_init_test) stop
     ! if (l_init_test) print*, 'sqerr', sqerr

   endif

   Mpc(k,j,:) = mom_pred(:,1)
   if (npmc < 4) then
      Mpr(k,j,:) = mom_pred(:,2)
   endif

   ! print*, 'ffcd', ffcd_mass(k,j,:)
   ! print*, 'sum(ffcd)', sum(ffcd_mass(k,j,:))*col*QtoM3

   ! open(20, file = 'ffcd_'//trim(tORf)//'.txt')
   ! write(20, '(34e15.7)') real(ffcd_mass(k,j,:))
   ! write(20, '(34e15.7)') real(ffcd_num(k,j,:))
   ! close(20)

   guessc(k,j,2) = gam_param(1,1)
   guessc(k,j,1) = gam_param(2,1)
   guessr(k,j,2) = gam_param(1,2)
   guessr(k,j,1) = gam_param(2,2)

   ffcdr8_massinit(k,j,:) = dble(ffcd_mass(k,j,:))
   ffcdr8_numinit(k,j,:) = dble(ffcd_num(k,j,:))

   if (n_cat == 2) then
      call calcmoms(ffcdr8_massinit(k,j,:),ffcdr8_numinit(k,j,:),10,mc0(k,j,1:10),mr0(k,j,1:10))
   elseif (n_cat == 1) then
      call calcmoms_sc(ffcdr8_massinit(k,j,:), ffcdr8_numinit(k,j,:), 10, mc0(k,j,1:10))
   endif

   if (npm == 4) then
     if (mc0(k,j,4)>1 .or. mc0(k,j,4).ne.mc0(k,j,4)) then
       print*, 'mpc(k,j,:)',mpc(k,j,:)
       print*, 'ffcd',ffcd_mass(k,j,:), ffcd_num(k,j,:)
       print*, 'mc0(k,j,pmomsc(1:4)+1)',mc0(k,j,pmomsc(1:4)+1)
       print*, 'gam_param(1,:)', gam_param(1,:)
       stop 'mass too high'
     endif
   endif

   ! if (mpc(k,j,1)>0) then
   !   ! print*, 'pmomsc(1:4)+1',pmomsc(1:4)+1
   !   print*, 'before', k, sngl(mpc(k,j,:))
   !   print*, 'after ', k, sngl(mc0(k,j,pmomsc(1:4)+1))
   !   print*, ''
   !   ! stop
   ! endif

   ! print*, 'sqerr', sqerr
   ! print*, 'gam_param(1,:)',gam_param(1,:)
   ! print*, 'Mpc', Mpc(k,j,:)
   ! print*, 'mc0', mc0(k,j,(/4,1,5,6/))
   ! print*, 'ffcd', ffcd_mass(k,j,:)
   ! if (l_init_test) print*, 'mc0', mc0(k,j,(/4,1,5,6/))
   ! if (l_init_test) print*, 'mr0', mr0(k,j,(/4,1,5,6/))
   ! stop


 enddo
enddo
call CPU_TIME(diag_dt2)
diag_dt3 = diag_dt3 + diag_dt2-diag_dt1

     ! print*, 'invert_moments', diag_dt3
   if (l_init_test) stop
! stop

! if (i_dgtime>=272 .and. i_dgtime<=273) then
!   k=8
!   print*, 'rain', mr0(k,1,4)
!   print*, guessc(k,1,2)
!   print*, guessr(k,1,2)
!   print*, 'Mpc', Mpc(k,1,:)
!   print*, 'Mpr', Mpr(k,1,:)
!   print*, 'ffcd', ffcd_mass(k,1,:)
!   print*, 'amp_distm_dble', tempvar_debug(2:35)
!   print*, 'amp_distm', rtemp_debug(1:34)
!   ! print*, 'mdpart', tempvar_debug(35:68)
!   ! print*, 'mdist', tempvar_debug(69:102)
! endif
! if (i_dgtime == 273) stop

! if (i_dgtime>=184 .and. i_dgtime <= 185) then
!   k=28
!   print*, ''
!   print*, 'Mpc', Mpc(k,1,:)
!   print*, 'cloud mass', sum(ffcd_mass(k,1,1:split_bins))*col
!   print*, 'rain mass', sum(ffcd_mass(k,1,split_bins+1:nkr))*col
!   print*, 'mc0', mc0(k,1,(/4,1,5,6/))
!   print*, 'guess', guessc(k,1,2), guessr(k,1,2)
!   ! print*, 'ffcd', ffcd_mass(k,j,:)
!   print*, ''
!   ! print*, 'mc0+mr0', sum(ffcd_mass(k,1,1:nkr))*col
!   ! if (sum(ffcd_mass(k,1,split_bins+1:nkr))*col>1e-3) then
!   !   stop
!   ! endif
! endif
! if (i_dgtime == 185) stop

!------CALL MICROPHYSICS--------------------

! ffcd_mass=real(ffcdr8_mass)

fncn=real(fncnr8)
! print*, 'ffcd before mphys', ffcd_mass(40,1,:)
if (bintype .eq. 'sbm') then
   call micro_proc_sbm(press,tempk,qv,fncn,ffcd_mass)
elseif (bintype .eq. 'tau') then
   ! ffcd_num=real(ffcdr8_num)
   call micro_proc_tau(tempk,qv,ffcd_mass,ffcd_num)
   ffcdr8_num=dble(ffcd_num)
endif

ffcdr8_mass=dble(ffcd_mass)
fncnr8=dble(fncn)

ffcdr8_massf = ffcdr8_mass
ffcdr8_numf = ffcdr8_num

!---------CALC MOMENTS-----------------------
do k=1,nz
  do j=1,nx
    if (npm < 4) then
      ! 2-cat AMP: {{{
      call calcmoms(ffcdr8_mass(k,j,:),ffcdr8_num(k,j,:),10,mc(k,j,:),mr(k,j,:))
      !---------UPDATE MOMENTS---------------------------------
      !Mp=predicted value before microphysics. m0=value after finding parameters
      !If m0 .ne. Mp, then finding parameters failed
      !If Mp<m0, it's possible to predict too large a change to M
      !To avoid, this, scale the change by Mp/m0.
      !mc=diagnosed value after microphysics. If there are no issues
      !then m becomes the new Mp

      do i=1,npmc  !Loop over the moments
        ip=pmomsc(i)+1 !pmomsc contains the predicted moments. +1 to get correct index in mc and mc0
        if (mc0(k,j,ip) > 0.) then
          dummy=mc(k,j,ip)
          ! dummy=(mc(k,j,ip)-mc0(k,j,ip))*Mpc(k,j,i)/mc0(k,j,ip)+Mpc(k,j,i)
          if(dummy<0.) then
            Mpc(k,j,i) = 0.
          else
            Mpc(k,j,i)=dummy
          endif
        else
          Mpc(k,j,i)=mc(k,j,ip)
        endif
        mc(k,j,ip)=Mpc(k,j,i)
      enddo

      do i=1,npmr
        ip=pmomsr(i)+1
        if (mr0(k,j,ip)>0.) then
          dummy=mr(k,j,ip)
          ! dummy=(mr(k,j,ip)-mr0(k,j,ip))*Mpr(k,j,i)/mr0(k,j,ip)+Mpr(k,j,i)
          if(dummy<0.) then
            Mpr(k,j,i)=0.
          else
            Mpr(k,j,i)=dummy
          endif
          if (i>=2) then !check new mean size
            realpmom = real(pmomsr(i))
            newdiam = (Mpr(k,j,1)/Mpr(k,j,i))**(1./(3.-realpmom))
            if (newdiam < diams(split_bins+1)) then
              Mpr(k,j,i) = Mpr(k,j,1)/(1.01*diams(split_bins+1))**(3.-realpmom)
            elseif (newdiam > diams(nkr)) then
              Mpr(k,j,i) = Mpr(k,j,1)/(0.99*diams(nkr))**(3.-realpmom)
            endif
          endif
        else
          Mpr(k,j,i)=mr(k,j,ip)
        endif
        ! mr(k,j,ip)=Mpr(k,j,i)
      enddo

      ! if (Mpr(k,j,1) .ne. Mpr(k,j,1)) then
      !    print*, 'mr', mr(k,j,:)
      !    print*, 'mr0', mr0(k,j,:)
      !    print*, 'mpr', Mpr(k,j,:)
      !    stop
      ! endif

      if (any(mc(k,j,:)==0.)) then
        mc(k,j,:)=0;Mpc(k,j,:)=0
      endif
      if (any(mr(k,j,:)==0.)) then
        ! mr(k,j,:)=0;Mpr(k,j,:)=0
      endif
      ! }}}
    else
      ! 1-cat AMP: {{{
      call calcmoms_sc(ffcdr8_mass(k,j,:), ffcdr8_num(k,j,:), 10, mc(k,j,:))
      do i=1,npm  !Loop over the moments
        ip=pmomsc(i)+1 !pmomsc contains the predicted moments. +1 to get correct index in mc and mc0
        if (mc0(k,j,ip) > 0.) then
          dummy=mc(k,j,ip)!
          ! dummy = (mc(k,j,ip)-mc0(k,j,ip))*Mpc(k,j,i)/mc0(k,j,ip)+Mpc(k,j,i)
          if(dummy<0.) then
            Mpc(k,j,i) = 0.
          else
            ! if (l_truncated) then
              Mpc(k,j,i)=dummy
            ! else
              ! if (dummy .ne. dummy) then
              !   print*, 'mc', mc(k,j,:)
              !   stop
              ! endif
              ! Mpc(k,j,i)=dummy
              ! Mpc(k,j,i)=mc(k,j,ip)-mc0(k,j,ip) + Mpc(k,j,i)
            ! endif
          endif

          ! take a look at this part
          if (i>=2) then !check new mean size
            realpmom = real(pmomsc(i))
            newdiam = (Mpc(k,j,1)/Mpc(k,j,i))**(1./(3.-realpmom))
            if (newdiam < diams(1)) then
              Mpc(k,j,i) = Mpc(k,j,1)/(1.01*diams(1))**(3.-realpmom)
            elseif (newdiam > diams(nkr)) then
              Mpc(k,j,i) = Mpc(k,j,1)/(0.99*diams(nkr))**(3.-realpmom)
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
      enddo

      ! if (any(mc(k,j,:)==0.)) then
      if (any(mc(k,j,:)==0.)) then
        mc(k,j,:)=0; Mpc(k,j,:)=0
      else ! update mc and mr as diagnosed cloud/rain moments
        call calcmoms(ffcdr8_mass(k,j,:), ffcdr8_num(k,j,:), 10, mc(k,j,:), mr(k,j,:))
      endif


      ! }}}
    endif

  enddo
enddo

end subroutine mp_amp
!---------------------------------------------------------------------
subroutine mp_sbm(ffcdr8,press,tempk,qv,fncnr8,mc,mr)

use module_hujisbm
use micro_prm
use parameters, only: nx,nz,num_h_moments,max_nbins
use global_fun

implicit none
integer:: i,j,k,ip
real, dimension(nz,nx)::tempk,press,qv
real, dimension(nz,nx,max_nbins)::ffcd,fncn
real(8), dimension(nz,nx,max_nbins)::ffcdr8,fncnr8
real(8), dimension(max_nbins)::dummy
real(8),dimension(nz,nx,10) :: mc,mr

!------CALL MICROPHYSICS--------------------
ffcd=real(ffcdr8)
fncn=real(fncnr8)
! print*, 'ffcd before mphys', ffcd(40,1,:)
call micro_proc_sbm(press,tempk,qv,fncn,ffcd)
! print*, 'ffcd after mphys', ffcd(40,1,:)
ffcdr8=dble(ffcd)
fncnr8=dble(fncn)

!---------CALC MOMENTS-----------------------
do k=1,nz
 do j=1,nx
  call calcmoms(ffcdr8(k,j,:),dummy,10,mc(k,j,:),mr(k,j,:))
 enddo
enddo

! print*, 'meand after mphys', get_meandiam(mc(40,1,4)+mr(40,1,4), mc(40,1,1)+mr(40,1,1))
! print*, 'mom45 after mphys', mc(40,1,5)+mr(40,1,5), mc(40,1,6)+mr(40,1,6)

end subroutine mp_sbm

!---------------------------------------------------------------------
subroutine mp_tau(ffcdr8_mass2d,ffcdr8_num2d,tempk,qv,mc,mr)

use module_hujisbm
use micro_prm
use parameters, only: nx,nz,num_h_moments,max_nbins
use namelists, only: l_advect
use mphys_tau_bin_declare, only: lk_cloud,xk,xkgmean
use global_fun
implicit none
integer:: i,j,k,ip

real, dimension(nz,nx):: tempk,qv
real(8), dimension(nz,nx,max_nbins)::ffcdr8_mass2d,ffcdr8_num2d
real, dimension(nz,nx,max_nbins)::ffcd_mass2d,ffcd_num2d
real(8),dimension(nz,nx,10) :: mc,mr ! moments

!------CALL MICROPHYSICS--------------------
ffcd_mass2d=real(ffcdr8_mass2d)
ffcd_num2d=real(ffcdr8_num2d)
! print*, 'ffcd before mphys', ffcd_mass2d(40,1,:)
call micro_proc_tau(tempk,qv,ffcd_mass2d,ffcd_num2d)
! print*, 'ffcd after mphys', ffcd_mass2d(40,1,:)
ffcdr8_mass2d=dble(ffcd_mass2d)
ffcdr8_num2d=dble(ffcd_num2d)
!---------CALC MOMENTS-----------------------
do k=1,nz
   do j=1,nx
      call calcmoms(ffcdr8_mass2d(k,j,:),ffcdr8_num2d(k,j,:),10,mc(k,j,:),mr(k,j,:))
   enddo
enddo

! print*, 'meand after mphys', get_meandiam(mc(40,1,4)+mr(40,1,4), mc(40,1,1)+mr(40,1,1))
! print*, 'mom45 after mphys', mc(40,1,5)+mr(40,1,5), mc(40,1,6)+mr(40,1,6)

! there might be a better way to calculate moments, but will leave it like that for now. -ahu
end subroutine mp_tau

subroutine calcmoms(ffcdr8_mass,ffcdr8_num,momnum,mc,mr)

use module_hujisbm
use micro_prm
Use namelists, only: bintype
Use physconst, only: pi
use parameters, only: split_bins

implicit none
integer :: i,ib,ip,momnum
real(8), dimension(max_nbins)::ffcdr8_mass
real(8), dimension(max_nbins) :: ffcdr8_num
real(8), dimension(nkr) :: diag_m, diag_D !diagnosed mass and diam of each bin
real(8), dimension(momnum) :: mc,mr ! moments
real(8) :: inf=huge(mc(1))

if (bintype .eq. 'tau') then
   diag_m=ffcdr8_mass/ffcdr8_num
   diag_D=(diag_m*QtoM3)**(1./3.)
   do ib=1,nkr
      if ((diag_D(ib) .ne. diag_D(ib)) .or. (diag_D(ib)>inf)) diag_D(ib)=diams(ib)
   end do
endif

do i=1,momnum
   if (bintype .eq. 'sbm') then
      mc(i)=sum(ffcdr8_mass(1:split_bins)/xl(1:split_bins)*diams(1:split_bins)**(i-1.))*col*1000.
      mr(i)=sum(ffcdr8_mass(split_bins+1:nkr)/xl(split_bins+1:nkr)*diams(split_bins+1:nkr)**(i-1.))*col*1000.
   elseif (bintype .eq. 'tau') then
      mc(i)=sum(ffcdr8_num(1:split_bins)*diag_D(1:split_bins)**(i-1))*col
      mr(i)=sum(ffcdr8_num(split_bins+1:nkr)*diag_D(split_bins+1:nkr)**(i-1))*col
   endif
enddo

end subroutine calcmoms


subroutine calcmoms_sc(ffcdr8_mass,ffcdr8_num,momnum,mc)

use module_hujisbm
use micro_prm
Use namelists, only: bintype
Use physconst, only: pi

implicit none
integer :: i,ib,ip,momnum
real(8), dimension(max_nbins)::ffcdr8_mass
real(8), dimension(max_nbins) :: ffcdr8_num
real(8), dimension(nkr) :: diag_m, diag_D !diagnosed mass and diam of each bin
real(8), dimension(momnum) :: mc,mr ! moments
real(8) :: inf=huge(mc(1))

if (bintype .eq. 'tau') then
   diag_m=ffcdr8_mass/ffcdr8_num
   diag_D=(diag_m/(1000.*pi/6))**(1./3.)
   do ib=1,nkr
      if ((diag_D(ib) .ne. diag_D(ib)) .or. (diag_D(ib)>inf)) then
         diag_D(ib)=diams(ib)
      endif
   end do
endif

do i=1,momnum
   if (bintype .eq. 'sbm') then
      mc(i)=sum(ffcdr8_mass(1:nkr)/xl*diams**(i-1.))*col*1000.
   elseif (bintype .eq. 'tau') then
      mc(i)=sum(ffcdr8_num(1:nkr)*diag_D**(i-1.))*col
   endif
enddo

end subroutine calcmoms_sc
