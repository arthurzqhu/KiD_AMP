!A. Igel - 5/2013 - Adapting for CLOUD model
! This module adapts the SBM from wrf, originally from A. Khain.

Subroutine micro_proc_sbm(press,tempk,qv,fncn,ffcd)

use module_hujisbm
use micro_prm
use physconst, only: pi
use parameters, only: nx,nz,dt,aero_N_init,max_nbins,split_bins
Use diagnostics, only: save_dg, i_dgtime, save_binData
use column_variables, only: z_half
use namelists, only: mp_proc_dg
IMPLICIT NONE

REAL::pcgs,rhocgs
REAL, DIMENSION(nz,nx)::tempk,press,qv
REAL, DIMENSION(nz,nx,max_nbins)::ffcd,fncn
! SBM VARIABLES
REAL,DIMENSION (max_nbins) :: FF1IN,FF3IN,FF4IN,FF5IN,&
                        FF1R,FF3R,FF4R,FF5R,FCCN,FIN
REAL,DIMENSION (max_nbins,icemax) :: FF2IN,FF2R
REAL, DIMENSION(nz) :: rhocgs_z,pcgs_z,zcgs_z
REAL, DIMENSION(nz,max_nbins) :: vr1z,ff1z

REAL :: SUP2_OLD,dthalf

! tmp vars to calculate process rate. not necessary but make the code look cleaner -ahu
real, dimension(max_nbins) :: ff1r_prev
REAL, DIMENSION(nz,max_nbins) :: ff1z_prev
real :: proc_tmass, proc_cmass, proc_rmass

DOUBLE PRECISION :: ss_max=0.003d0
double precision del1in_limit
DOUBLE PRECISION DEL1NR,DEL2NR,DEL12R,DEL12RD,ES1N,ES2N,EW1N,EW1PN
DOUBLE PRECISION DELSUP1,DELSUP2,DELDIV1,DELDIV2
DOUBLE PRECISION TT,QQ,TTA,QQA,PP,QQt
DOUBLE PRECISION DIV1,DIV2,DIV3,DIV4,DEL1IN,DEL2IN,DEL1AD,DEL2AD
!For bulk nucleation
REAL FACTZ,CONCCCN,CONCDROP
!Flags
INTEGER ISYM1,ISYM3,ISYM4,ISYM5
INTEGER, DIMENSION (3):: ISYM2
INTEGER DIFFU
!!! For CCN regeneration
real fccn0(max_nbins)
real :: ndrop, subtot  ! for diagnostic CCN
!Functions
real :: sum_pris, sum_snow
INTEGER :: k,i,j
integer itimestep, kr, ikl,myflag

!Zero Out Budgets
CALL zero_budgets()
dthalf=dt/2.
do i=1,nx
    do k=1,nz
         !Adele - check units on everything
         PP = press(k,i)*10.
         pcgs = PP
         TT=tempk(k,i)
         QQ=qv(k,i)
         rhocgs=PP/TT/287./10000.
         !Convert distributions from d(mixing ratio)/dln(r) to d(#/cm3)/d(mass)
         ! CCN

         !Convert distributions from d(mixing ratio)/dln(r) to
         !d(#/cm3)/d(mass)/mass. This is a meaningless quantity, but it is
         !treated this way through the remainder of the bin micro code.
         !So be it. Adele

         DO KR=1,NKR
            fccn(kr)=0.
            !FCCN(KR)=fncn(k,i,kr)*rhocgs/xccn(kr)
         END DO
         ! LIQUID
         DO KR=1,NKR
            FF1R(KR)=ffcd(k,i,kr)*rhocgs/xl(kr)/xl(kr)/3.0
         END DO

         FF2R=0.;FF3R=0.;FF4R=0.;FF5R=0.
         ! ICE
         !IF (ICEPROCS.EQ.1)THEN
         !   DO KR=1,NKR
               ! COLUMNS
         !      IF (IPRIS == 1 .or. IPRIS >= 4) THEN
         !         FF2R(KR,1)=ffic(kr)*rhocgs/xi(kr,1)/xi(kr,1)/3.0
         !      ENDIF
               ! PLATES
         !      IF (IPRIS == 2 .or. IPRIS >= 4) THEN
         !         FF2R(KR,2)=ffip(kr)*rhocgs/xi(kr,2)/xi(kr,2)/3.0
         !      ENDIF
               ! DENDRITES
         !      IF (IPRIS >= 3) THEN
         !         FF2R(KR,3)=ffid(kr)*rhocgs/xi(kr,3)/xi(kr,3)/3.0
         !      ENDIF
               ! SNOW/RAMS AGGREGATES
         !      FF3R(KR)=ffsn(kr)*rhocgs/xs(kr)/xs(kr)/3.0
               ! GRAUPEL
         !      IF (IGRAUP > 0) THEN
         !         FF4R(KR)=ffgl(kr)*rhocgs/xg(kr)/xg(kr)/3.0
         !      ENDIF
               ! HAIL
         !      IF (IHAIL > 0) THEN
         !         FF5R(KR)=ffhl(kr)*rhocgs/xh(kr)/xh(kr)/3.0
         !      ENDIF
         !   END DO
         !ENDIF
!-----------------------------------------------
       IF (docondensation) THEN
!Condensation and Nucleation Preparation
            !old values, 1 liquid, 2 ice
            ES1N=AA1_MY*DEXP(-BB1_MY/TT)
            ES2N=AA2_MY*DEXP(-BB2_MY/TT)
            EW1N=QQ*PP/(0.622+0.378*QQ)
            DIV1=EW1N/ES1N        !Saturation ratio wrt liquid
            DEL1IN=EW1N/ES1N-1.   !supersat ratio wrt liquid
            DIV2=EW1N/ES2N        !Saturation ratio wrt ice
            DEL2IN=EW1N/ES2N-1.   !supersat ratio wrt liquid
            SUP2_OLD=DEL2IN
!---------------------------------------------------
            !Latent heating budget prep
           ! CALL LHV_BUDGET(FF1R,FF2R,FF3R,FF4R,FF5R, &
           !                 EXNER,RHOCGS,'beg')
            FF1IN=FF1R
            FF2IN=0.
           ! FF2IN=FF2R
           ! FF3IN=FF3R
           ! FF4IN=FF4R
           ! FF5IN=FF5R

            IF (donucleation) THEN
              IF (DEL1IN.GT.0.) THEN
                 !aero_N_init=#/kg, concccn=#/g
                 concccn = aero_N_init(1)*1.e-3*(100.*del1in)**bccn
                 concdrop = sum(ff1in*xl)*3.*col/rhocgs
                 IF (concccn>concdrop) THEN
                    ff1in(1)=ff1in(1)+(concccn-concdrop)*rhocgs/(3.0*col*xl(1))

                    if (mp_proc_dg) then
                       proc_tmass = (ff1in(1)-ff1r(1))/(rhocgs/xl(kr)/xl(kr)/3.0)*col
                       if (nx==1) then
                          call save_dg(k,proc_tmass/dt,'dm_nuc',i_dgtime,units='kg/kg/s',dim='z')
                       else
                          call save_dg(k,i,proc_tmass/dt,'dm_nuc',i_dgtime,units='kg/kg/s',dim='z')
                       endif


                    endif
                 ENDIF
              ENDIF

           !  IF (DEL1IN.GT.0.OR.(ICEPROCS.EQ.1.and.DEL2IN.GT.0))THEN !If supersat wrt liq or ice
           !  !mo Limit SS for the droplet nucleation
           !  del1in_limit = min(DEL1IN,ss_max)
           !  !NUCLEATION
           !  CALL JERNUCL01(FF1IN,FF2IN,FCCN &
           !               ,TT,rhocgs &
           !               ,DEL1IN_limit,DEL2IN &
           !               ,SUP2_OLD &
           !               ,dtlt,fin &
           !               ,nuccldrt,nuccldct &
           !               ,nucicert,nucicect &
           !               ,inucifnrt,inucifnct)
           ! ENDIF
           ENDIF
!-------------------------------------------------------------------
!Condensation/Deposition/Evaporation/Sublimation
           FF1R=FF1IN
           !FF2R=FF2IN
           !Flags to determine if deposition/sublimation is necessary
           ISYM1=0;ISYM2=0;ISYM3=0;ISYM4=0;ISYM5=0
           if (any(FF1R.gt.1.e-6))ISYM1=1
               !if (iceprocs.eq.1) then
               !   if(any(FF2R.gt.1.e-6))ISYM2=1
               !   if(any(FF3R.gt.1.e-6))ISYM3=1
               !   if(any(FF4R.gt.1.e-6))ISYM4=1
               !   if(any(FF5R.gt.1.e-6))ISYM5=1
               !endif

               !CALL VAP_BUDGET(FF1R,FF2R,FF3R,FF4R,FF5R,RHOCGS,'beg')

               !If liquid
               IF (ISYM1.EQ.1.AND.((TT-273.15).GT.-0.187.OR. &
                  (sum(ISYM2).EQ.0.AND.ISYM3.EQ.0.AND.ISYM4.EQ.0.AND.ISYM5.EQ.0)))THEN

                  if (mp_proc_dg) ff1r_prev=ff1r

                  CALL ONECOND1(TT,QQ,PP,rhocgs &
                          ,VR1,pcgs &
                          ,DEL1IN,DEL2IN,DIV1,DIV2 &
                          ,FF1R,FF1IN,XL,RLEC,RO1BL &
                          ,AA1_MY,BB1_MY,AA2_MY,BB2_MY &
                          ,C1_MEY,C2_MEY &
                          ,COL,DTCOND,ICEMAX,NKR,ISYM1 &
                          ,ISYM2,ISYM3,ISYM4,ISYM5,1,1,1,0.,1.,1)

                  ! diag condensation rate
                  if (mp_proc_dg) then
                     proc_cmass=0.
                     proc_rmass=0.
                     proc_tmass=0.

                     do kr=1,split_bins
                        proc_cmass = proc_cmass + (ff1r(kr)-ff1r_prev(kr))/(rhocgs/xl(kr)/xl(kr)/3.0)*col
                        proc_tmass = proc_tmass + (ff1r(kr)-ff1r_prev(kr))/(rhocgs/xl(kr)/xl(kr)/3.0)*col
                     enddo

                     do kr=split_bins,nkr
                        proc_rmass = proc_rmass + (ff1r(kr)-ff1r_prev(kr))/(rhocgs/xl(kr)/xl(kr)/3.0)*col
                        proc_tmass = proc_tmass + (ff1r(kr)-ff1r_prev(kr))/(rhocgs/xl(kr)/xl(kr)/3.0)*col
                     enddo

                     if (nx==1) then
                        call save_dg(k,proc_cmass/dt,'dm_cloud_ce',i_dgtime,units='kg/kg/s',dim='z')
                        call save_dg(k,proc_rmass/dt,'dm_rain_ce',i_dgtime,units='kg/kg/s',dim='z')
                        call save_dg(k,proc_tmass/dt,'dm_ce',i_dgtime,units='kg/kg/s',dim='z')
                     else
                        call save_dg(k,i,proc_cmass/dt,'dm_cloud_ce',i_dgtime,units='kg/kg/s',dim='z,x')
                        call save_dg(k,i,proc_rmass/dt,'dm_rain_ce',i_dgtime,units='kg/kg/s',dim='z,x')
                        call save_dg(k,i,proc_tmass/dt,'dm_ce',i_dgtime,units='kg/kg/s',dim='z,x')
                     endif

                  endif
                  !If ice
                  !ELSE IF(ISYM1.EQ.0.AND.(TT-273.15).LE.-0.187.AND. &
                  !   (any(ISYM2.EQ.1).OR.ISYM3.EQ.1.OR.ISYM4.EQ.1.OR.ISYM5.EQ.1))THEN
                  !   print*,'not supported yet'
                  !   stop
                    ! CALL ONECOND2(TT,QQ,PP,rhocgs &
                    !      ,VR2,VR3,VR4,VR5,pcgs &
                    !      ,DEL1IN,DEL2IN,DIV1,DIV2 &
                    !      ,FF2R,FF2IN,XI,RIEC,RO2BL &
                    !      ,FF3R,FF3IN,XS,RSEC,RO3BL &
                    !      ,FF4R,FF4IN,XG,RGEC,RO4BL &
                    !      ,FF5R,FF5IN,XH,RHEC,RO5BL &
                    !      ,AA1_MY,BB1_MY,AA2_MY,BB2_MY &
                    !      ,COL,DTCOND,ICEMAX,NKR &
                    !      ,ISYM2,ISYM3,ISYM4,ISYM5)
                  !If mixed-phase
                  !ELSE IF(ISYM1.EQ.1.AND.(TT-273.15).LE.-0.187.AND. &
                  !   (any(ISYM2.EQ.1).OR.ISYM3.EQ.1.OR.ISYM4.EQ.1.OR.ISYM5.EQ.1))THEN
                  !   print*,'not supported yet'
                  !   stop
                  !   CALL ONECOND3(TT,QQ,PP,rhocgs &
                  !        ,VR1,VR2,VR3,VR4,VR5,pcgs &
                  !        ,DEL1IN,DEL2IN,DIV1,DIV2 &
                  !        ,FF1R,FF1IN,XL,RLEC,RO1BL &
                  !        ,FF2R,FF2IN,XI,RIEC,RO2BL &
                  !        ,FF3R,FF3IN,XS,RSEC,RO3BL &
                  !        ,FF4R,FF4IN,XG,RGEC,RO4BL &
                  !        ,FF5R,FF5IN,XH,RHEC,RO5BL &
                  !        ,AA1_MY,BB1_MY,AA2_MY,BB2_MY &
                  !        ,COL,DTCOND,ICEMAX,NKR &
                  !        ,ISYM1,ISYM2,ISYM3,ISYM4,ISYM5)
                  !END IF

                  !CALL VAP_BUDGET(FF1R,FF2R,FF3R,FF4R,FF5R,RHOCGS,'end')
           ENDIF
       endif !if condensation
  !-----------------------------------------------------------------------------------------
  !------------------COLLISION-COALESCENCE--------------------------------------------------
       IF (docollisions) then
            IF ((iceprocs.eq.0 .and. any(FF1R>0.)) &
                 .or. any(FF1R>0.).or.any(FF2R>0.).or.any(FF3R>0.))  THEN
               !Graupel and hail do not self-collect, so don't need to call collection if they
               !are the only species present
               cld2raint=0.

               if (mp_proc_dg) ff1r_prev=ff1r
               
               CALL COAL_BOTT_NEW(FF1R,FF2R,FF3R, &
                    FF4R,FF5R,TT,QQ,PP,rhocgs,dthalf,TCRIT,TTCOAL, &
                    cld2raint,rimecldt, &
                    aggregatet,rimecldsnowt, &
                    rimecldaggrt,rimecldgraut, &
                    rimecldhailt,rain2icet, &
                    rain2snt,rain2agt, &
                    rain2grt,rain2hat)

               if (mp_proc_dg) then

                  proc_cmass=0.
                  proc_rmass=0.
!                  proc_tmass=0.

                  do kr=1,split_bins
                     proc_cmass = proc_cmass + (ff1r(kr)-ff1r_prev(kr))/(rhocgs/xl(kr)/xl(kr)/3.0)*col
!                     proc_tmass = proc_tmass + (ff1r(kr)-ff1r_prev(kr))/(rhocgs/xl(kr)/xl(kr)/3.0)*col
                  enddo

                  do kr=split_bins,nkr
                     proc_rmass = proc_rmass + (ff1r(kr)-ff1r_prev(kr))/(rhocgs/xl(kr)/xl(kr)/3.0)*col
!                     proc_tmass = proc_tmass + (ff1r(kr)-ff1r_prev(kr))/(rhocgs/xl(kr)/xl(kr)/3.0)*col
                  enddo

                  if (nx==1) then
                     call save_dg(k,proc_cmass/dt,'dm_cloud_coll',i_dgtime,units='kg/kg/s',dim='z')
                     call save_dg(k,proc_rmass/dt,'dm_rain_coll',i_dgtime,units='kg/kg/s',dim='z')
                     !call save_dg(k,proc_tmass/dt,'dm_coll',i_dgtime,units='kg/kg/s',dim='z')
                  else
                     call save_dg(k,i,proc_cmass/dt,'dm_cloud_coll',i_dgtime,units='kg/kg/s',dim='z,x')
                     call save_dg(k,i,proc_rmass/dt,'dm_rain_coll',i_dgtime,units='kg/kg/s',dim='z,x')
                     !call save_dg(k,i,proc_tmass/dt,'dm_coll',i_dgtime,units='kg/kg/s',dim='z,x')
                  endif

               endif

            END IF
       END IF
  !-----------------------------------------------------------------------------------------
  !update temperature and water vapor
       tempk(k,i)=TT
       qv(k,i)=QQ
  !Prepare for sedimentation
       rhocgs_z(k)=rhocgs
       pcgs_z(k)=pcgs
       zcgs_z(k)=z_half(k)*100. !convert to cgs
       vr1z(k,:)=vr1(:)*dsqrt(1.d6/pcgs)
       ff1z(k,:)=ff1r(:)/rhocgs
        ! if (iceprocs.eq.1) then
        !    vr2z=
        !    vr3z=
        !    vr4z=
        !    vr5z=
        ! endif
   enddo !end loop over k

   ! +-----------------------------+
   ! Hydrometeor Sedimentation
   ! +-----------------------------+
   ! ... Drops ...
   if(dosedimentation) then
      if (mp_proc_dg) ff1z_prev=ff1z

      call falfluxhucm_z(ff1z,vr1z,rhocgs_z,pcgs_z,zcgs_z,dt,1,nz,nkr)

      if (mp_proc_dg) then
         proc_tmass=0.

         ! probably not the most efficient way to implement -ahu
         do k=1,nz
            do kr=1,nkr
               proc_tmass = proc_tmass + (ff1z(k,kr)-ff1z_prev(k,kr))/(1./xl(kr)/xl(kr)/3.0)*col
            enddo

            if (nx==1) then
               call save_dg(k,proc_tmass/dt,'dm_sed',i_dgtime,units='kg/kg/s',dim='z')
            else
               call save_dg(k,i,proc_tmass/dt,'dm_sed',i_dgtime,units='kg/kg/s',dim='z,x')
            endif
         enddo

      endif


   endif
   ! if(iceprocs == 1)then
   !    call FALFLUXHUCM_Z(ffx_z,VRX,RHOCGS_z,PCGS_z,ZCGS_z,DT,kts,kte,nkr)
   !    call FALFLUXHUCM_Z(ffx_z,VRX,RHOCGS_z,PCGS_z,ZCGS_z,DT,kts,kte,nkr)
   ! end if ! if (iceprocs == 1)
!------------------------------------------------------
   do k=1,nz
     ! Convert back to f(diam)
     DO KR=1,NKR
        fncn(k,i,KR)=0.
        !fncn(k,i,KR)=FCCN(KR)/rhocgs*xccn(kr)
        ffcd(k,i,KR)=ff1z(k,kr)*xl(kr)*xl(kr)*3.0
     END DO
     !IF (ICEPROCS.EQ.1)THEN
     !   DO KR=1,NKR
     !     if (ipris == 1 .or. ipris >= 4) &
     !        ffic(kr)=FF2R(KR,1)/rhocgs*xi(kr,1)*xi(kr,1)*3.0
     !     if (ipris == 2 .or. ipris >= 4) &
     !        ffip(kr)=FF2R(KR,2)/rhocgs*xi(kr,2)*xi(kr,2)*3.0
     !     if (ipris >= 3) &
     !        ffid(kr)=FF2R(KR,3)/rhocgs*xi(kr,3)*xi(kr,3)*3.0
     !     ffsn(kr)=FF3R(KR)/rhocgs*xs(kr)*xs(kr)*3.0
     !     if (igraup > 0) &
     !        ffgl(kr)=FF4R(KR)/rhocgs*xg(kr)*xg(kr)*3.0
     !     if (ihail > 0) &
     !        ffhl(kr)=FF5R(KR)/rhocgs*xh(kr)*xh(kr)*3.0
     !   END DO
     !END IF
   enddo
enddo !end loop over i
RETURN
END SUBROUTINE micro_proc_sbm

!---------------------------------------------------------------------------------
Subroutine micro_init_sbm

use micro_prm
use module_hujisbm
use column_variables, only:rho

IMPLICIT NONE
INTEGER,parameter :: hujisbm_unit1 = 22
integer :: kr
REAL PI, rhocgs, x0ccn
data pi/3.141592654/
! dtime - timestep of integration (calculated in main program) :
! ax - coefficient used for masses calculation
! ima(i,j) - k-category number, c(i,j) - courant number

dlnr=dlog(2.d0)/(3.d0*scal)

!--- Read in various lookup tables

900   FORMAT(6E13.5)
! MASSES :
OPEN(UNIT=hujisbm_unit1,FILE="./src/input_data/sbm_input/masses.asc", &
     FORM="FORMATTED")
READ(hujisbm_unit1,900) XL(otn),XI(otn,oti),XS(otn),XG(otn),XH(otn)
CLOSE(hujisbm_unit1)
! BULKRADIUS

OPEN(UNIT=hujisbm_unit1,FILE="./src/input_data/sbm_input/bulkradii.asc_s_0_03_0_9", &
     FORM="FORMATTED")
READ(hujisbm_unit1,*) RADXXO(otn,oth)
CLOSE(hujisbm_unit1)


DO KR=1,NKR
   DIAMS(KR)=RADXXO(KR,1)*2.*0.01
ENDDO

!Choose aerosol type
!!!!!!AEROMEDRAD currently uninitialized
!IF (aerotype == 1) THEN
!ammonium sulfate
!   ROCCN0 = 1.769
!   RO_SOLUTE = 1.769
!   MWAERO = 132.0
!   IONS = 3
!   RP_G = aeromedrad*100.
!   SIG_G = 1.8
!ELSEIF (aerotype == 2) THEN
!NaCl
!   ROCCN0 = 2.165
!   RO_SOLUTE = 2.165
!   MWAERO = 58.4
!   IONS = 2
!   RP_G = aeromedrad*100.
!   SIG_G = 1.8
!ELSE
!   PRINT*, 'INVALID AERO_CHEM TYPE. MUST EQUAL 1 OR 2', aerotype
!   STOP
!ENDIF

!Find CCN density, mass, and diameter for each bin
!X0CCN =XL(ICCN)/(2.**(NKR-1)) !Minimum ccn mass
!DO KR=1,NKR
!   ROCCN(KR)=ROCCN0
!   XCCN(KR)=X0CCN*2.**(KR-1)
!   RCCN(KR)=(3.*XCCN(KR)/4./3.141593/ROCCN(KR))**(1./3.)
!ENDDO

! Functional form
!DO KR=1,NKR
!  fccnr0(kr) = 1./(sqrt(2.*pi)*log(sig_g))    &
!               *exp(-(log(rccn(kr)/rp_g))**2./2.0/(log(sig_g))**2.)
!ENDDO
!fccnr0 = fccnr0 / (sum(fccnr0) * col)

!Calculations are done in CGS units
!do i=1,nx
! do k=1,nz
!   rhocgs = rho(k)/1000.
!   DO KR=1,NKR
!     fncn(k,i,kr)=FCCNR0(KR)/rhocgs*xccn(kr)*naero
!   ENDDO
! enddo
!enddo

return
END SUBROUTINE micro_init_sbm
!---------------------------------------------------------------------------------
Subroutine micro_init_sbm2()

use micro_prm
use module_hujisbm
use parameters, only:dt

IMPLICIT NONE
INTEGER KR
INTEGER,parameter :: hujisbm_unit1 = 22
REAL PI,NDTCOLL_REAL
data PI/3.141592654/

!--- Read in various lookup tables
! CAPACITIES :
! Capacities are used for the condensation rates
OPEN(UNIT=hujisbm_unit1,FILE="./src/input_data/sbm_input/capacity.asc",  &
    FORM="FORMATTED")
900   FORMAT(6E13.5)
READ(hujisbm_unit1,900) RLEC(otn),RIEC(otn,oti),RSEC(otn),RGEC(otn),RHEC(otn)
CLOSE(hujisbm_unit1)

! MASSES :
OPEN(UNIT=hujisbm_unit1,FILE="./src/input_data/sbm_input/masses.asc", &
     FORM="FORMATTED")
READ(hujisbm_unit1,900) XL(otn),XI(otn,oti),XS(otn),XG(otn),XH(otn)
CLOSE(hujisbm_unit1)


! TERMINAL VELOCITY :
OPEN(UNIT=hujisbm_unit1,FILE="./src/input_data/sbm_input/termvels.asc", &
     FORM="FORMATTED")
READ(hujisbm_unit1,*) VR1(otn),VR2(otn,oti),VR3(otn),VR4(otn),VR5(otn)
CLOSE(hujisbm_unit1)

! KERNELS DEPENDING ON PRESSURE :
OPEN(UNIT=hujisbm_unit1,FILE="./src/input_data/sbm_input/kernels_z.asc",  &
     FORM="FORMATTED")
READ(hujisbm_unit1,900) YWLL_1000MB(otn,otn),YWLL_750MB(otn,otn),YWLL_500MB(otn,otn)
CLOSE(hujisbm_unit1)

! KERNELS NOT DEPENDING ON PRESSURE :
OPEN(UNIT=hujisbm_unit1,FILE="./src/input_data/sbm_input/kernels.asc_s_0_03_0_9",  &
     FORM="FORMATTED")
READ(hujisbm_unit1,900) &
   YWLL(otn,otn),YWLI(otn,otn,oti),YWLS(otn,otn),YWLG(otn,otn),YWLH(otn,otn), &
   YWIL(otn,otn,oti),YWII(otn,otn,oti,oti),YWIS(otn,otn,oti),YWIG(otn,otn,oti),YWIH(otn,otn,oti), &
   YWSL(otn,otn),YWSI(otn,otn,oti),YWSS(otn,otn),YWSG(otn,otn),YWSH(otn,otn), &
   YWGL(otn,otn),YWGI(otn,otn,oti),YWGS(otn,otn),YWGG(otn,otn),YWGH(otn,otn), &
   YWHL(otn,otn),YWHI(otn,otn,oti),YWHS(otn,otn),YWHG(otn,otn),YWHH(otn,otn)
close (hujisbm_unit1)

! BULKDENSITY :
OPEN(UNIT=hujisbm_unit1,FILE="./src/input_data/sbm_input/bulkdens.asc_s_0_03_0_9", &
     FORM="FORMATTED")
READ(hujisbm_unit1,900) RO1BL(otn),RO2BL(otn,oti),RO3BL(otn),RO4BL(otn),RO5BL(otn)
CLOSE(hujisbm_unit1)

! BULKRADIUS
OPEN(UNIT=hujisbm_unit1,FILE="./src/input_data/sbm_input/bulkradii.asc_s_0_03_0_9", &
     FORM="FORMATTED")
READ(hujisbm_unit1,*) RADXXO(otn,oth)
CLOSE(hujisbm_unit1)

do kr=1,nkr
   xl_mg(kr)=xl(kr)*1.e3
   xs_mg(kr)=xs(kr)*1.e3
   xg_mg(kr)=xg(kr)*1.e3
   xh_mg(kr)=xh(kr)*1.e3
   xi1_mg(kr)=xi(kr,1)*1.e3
   xi2_mg(kr)=xi(kr,2)*1.e3
   xi3_mg(kr)=xi(kr,3)*1.e3
enddo

!Initialize some other parameters
DTCOND=dt/REAL(NCOND)

call courant_bott

CALL BREAKINIT
!  collision was called every NDTCOLL*dt
ndtcoll_real = real(ndtcoll)
call kernals(ndtcoll_real)

return
END SUBROUTINE micro_init_sbm2
!--------------------------------------------------
FUNCTION sum_pris(ff2r,rhocgs)

use micro_prm, only:nkr, icemax, col3, krpris, ipris
use parameters, only: max_nbins
use module_hujisbm, only:xi,sum_mass

implicit none
real :: sum_pris,rhocgs
real, dimension(max_nbins,icemax)::ff2r

IF (IPRIS >=4) THEN
    sum_pris = sum_mass(ff2r(:,1),XI(:,1),rhocgs,1,KRPRIS(1)-1) &
              + sum_mass(ff2r(:,2),XI(:,2),rhocgs,1,KRPRIS(2)-1) &
              + sum_mass(ff2r(:,3),XI(:,3),rhocgs,1,KRPRIS(3)-1)

ELSEIF (IPRIS==1 .or. IPRIS==2 .or. IPRIS==3) THEN
  sum_pris = sum_mass(ff2r(:,IPRIS),XI(:,IPRIS),rhocgs,1,KRPRIS(IPRIS)-1)
ELSE
  sum_pris = 0.
ENDIF

END FUNCTION sum_pris
!--------------------------------------------------
FUNCTION sum_snow(ff2r,rhocgs)

use micro_prm, only:nkr, icemax, col3, krpris, ipris
use module_hujisbm, only:xi,sum_mass
use parameters, only: max_nbins
implicit none
real:: sum_snow,rhocgs
real, dimension(max_nbins,icemax)::ff2r

IF (IPRIS >=4) THEN
    sum_snow = sum_mass(ff2r(:,1),XI(:,1),rhocgs,KRPRIS(1),NKR) &
             + sum_mass(ff2r(:,2),XI(:,2),rhocgs,KRPRIS(2),NKR) &
             + sum_mass(ff2r(:,3),XI(:,3),rhocgs,KRPRIS(3),NKR)

ELSEIF (IPRIS==1 .or. IPRIS==2 .or. IPRIS==3) THEN
  sum_snow = sum_mass(ff2r(:,IPRIS),XI(:,IPRIS),rhocgs,KRPRIS(IPRIS),NKR)
ELSE
  sum_snow = 0.
ENDIF

END FUNCTION sum_snow
!----------------------------------------------------
Subroutine zero_budgets()

use micro_prm

implicit none

!-------Zero out micro budgets-----------------------
if (imbudget >= 1) then
   latheatvap = 0.
   if (iceprocs == 1) latheatfrz = 0.
endif

   if (imbudget >= 1) then
      vapliqt = 0.
      latheatvapt = 0.
      if (iceprocs == 1) then
         latheatfrzt = 0.
         vapicet = 0.
      endif
   endif
   if (imbudget >= 2) then
      vapcldt = 0.
      vapraint = 0.
      if (iceprocs == 1) then
         if (ipris > 0) vapprist = 0.
         if (ipris > 0) vapsnowt = 0.
         vapaggrt = 0.
         if (igraup > 0) vapgraut = 0.
         if (ihail > 0) vaphailt = 0.
      endif
   endif
END SUBROUTINE zero_budgets

!-------------------------------------------------------------------
Subroutine vap_budget(ff1r,ff2r,ff3r,ff4r,ff5r,dens,str)

use micro_prm
use module_hujisbm, only:xl,xs,xg,xh,sum_mass
use parameters, only: max_nbins
implicit none

real, dimension(max_nbins) :: ff1r,ff3r,ff4r,ff5r
real, dimension(max_nbins,icemax) :: ff2r
real :: dens,plusminus,sum_pris,sum_snow
character(len=3) :: str

!plusminus = -1 for budget prep, = +1 for budget finish
if (str .eq. 'beg') then
   plusminus = -1.
elseif (str .eq. 'end') then
   plusminus = 1.
else
   print*, 'Invalid string for budgets'
   stop
endif

if (imbudget >= 1) then
   vapliqt = vapliqt + &
      plusminus * sum_mass(ff1r,xl(:),dens,1,nkr)

   if(iceprocs.eq.1) &
      vapicet = vapicet + plusminus * &
          sum_mass(ff2r(:,1)+ff2r(:,2)+ff2r(:,3)+ff3r+ff4r+ff5r,xs,dens,1,nkr)
endif
if (imbudget >= 2) then
    vapcldt = vapcldt +  &
                 plusminus * sum_mass(ff1r,xl,dens,1,krdrop-1)
    vapraint = vapraint + &
                 plusminus * sum_mass(ff1r,xl,dens,krdrop-1,nkr)
    if(iceprocs.eq.1) then
       if (ipris > 0) &
       vapprist = vapprist + &
                      plusminus * sum_pris(ff2r,dens)
       if (ipris > 0) &
       vapsnowt = vapsnowt + &
                      plusminus * sum_snow(ff2r,dens)
       vapaggrt = vapaggrt + &
                      plusminus * sum_mass(ff3r,xs,dens,1,nkr)
       if (igraup > 0) &
       vapgraut = vapgraut + &
                      plusminus * sum_mass(ff4r,xg,dens,1,nkr)
       if (ihail > 0) &
       vaphailt = vaphailt + &
                      plusminus * sum_mass(ff5r,xh,dens,1,nkr)
    endif
endif

return
END SUBROUTINE vap_budget
!-------------------------------------------------------------------
Subroutine lhv_budget(ff1r,ff2r,ff3r,ff4r,ff5r, &
                      pitot,dens,str)

use micro_prm
use module_hujisbm, only:xl,xs,sum_mass
use parameters, only: max_nbins
implicit none

real, dimension(max_nbins) :: ff1r,ff3r,ff4r,ff5r
real, dimension(max_nbins,icemax) :: ff2r
real :: pitot,dens,fac,plusminus
character(len=3) :: str

!plusminus = -1 for budget prep, = +1 for budget finish
if (str .eq. 'beg') then
   plusminus = -1.
elseif (str .eq. 'end') then
   plusminus = 1.
else
   print*, 'Invalid string for budgets'
   stop
endif

if (lhrtheta) then
   fac = 1./pitot
else
   fac = cpi
endif

if (imbudget >= 1) then
   latheatvap = latheatvap + &
      plusminus * alvl * fac * sum_mass(ff1r,xl(:),dens,1,nkr)

   if(iceprocs.eq.1) &
      latheatvap = latheatvap + &
         plusminus * alvi * fac * &
         sum_mass(ff2r(:,1)+ff2r(:,2)+ff2r(:,3)+ff3r+ff4r+ff5r,xs,dens,1,nkr)
endif

if (plusminus .eq. 1) then
   latheatvapt = latheatvapt + &
                                       latheatvap
endif
return
END SUBROUTINE lhv_budget
!-------------------------------------------------------------------
Subroutine lhf_budget(ff1r,pitot,dens,str)

use micro_prm
use module_hujisbm, only:xl,sum_mass
use parameters, only: max_nbins
implicit none

real, dimension(max_nbins) :: ff1r
real :: pitot,dens,fac,plusminus
character(len=3) :: str

!plusminus = -1 for budget prep, = +1 for budget finish
if (str .eq. 'beg') then
   plusminus = -1.
elseif (str .eq. 'end') then
   plusminus = 1.
else
   print*, 'Invalid string for budgets'
   stop
endif

if (lhrtheta) then
   fac = 1./pitot
else
   fac = cpi
endif

if (imbudget >= 1) then
   latheatfrz = latheatfrz + &
      plusminus * (-1.) * alli * fac * sum_mass(ff1r,xl(:),dens,1,nkr)
endif

if (plusminus .eq. 1) then
   latheatfrzt = latheatfrzt + &
                                       latheatfrz
endif
return
END SUBROUTINE lhf_budget
!-------------------------------------------------------------------
Subroutine init_dist_sbm(rxc,gnuc,dnc,rxr,gnur,dnr,diams,ffcd)
use micro_prm, only: nkr
use parameters, only: max_nbins
implicit none

integer :: kr
real :: rxc,gnuc,dnc,rxr,gnur,dnr
real(8):: n0c,exptermc,n0r,exptermr
real(8), dimension(max_nbins) :: ffcd,diams

!print*, 'inside init dist',rx,gnu,dn
!Setting up a mass distribution, not a number distribution
!So increase gnu by 3

ffcd=0.
n0c=rxc/gamma(gnuc+3)
n0r=rxr/gamma(gnur+3)

do kr=1,nkr
  if (rxc>0.) then
    exptermc=exp(-1.*diams(kr)/dnc)
    ffcd(kr) = n0c*exptermc*(diams(kr)/dnc)**(gnuc+3)
  endif
  if (rxr>0.) then
     exptermr=exp(-1.*diams(kr)/dnr)
     ffcd(kr) = ffcd(kr) + n0r*exptermr*(diams(kr)/dnr)**(gnur+3)
  endif
!If ffcd(kr) is NaN, set to zero
!    if (ffcd(kr).ne.ffcd(kr) .or. ffcd(kr)*0.0.ne.0.0 &
!       .or. ffcd(kr)/ffcd(kr).ne. 1.0) ffcd(kr)=0.
enddo

!print*, 'ffcd',ffcd

return
!Diagnose shape parameter
!  if (mom<0) then
!    !Shape parameter is given by mx
!    gnu=mx
!  elseif (mom==1) then
!    xx=cx**(-2./3.)*mx/m3**(1./3.)
!    if (xx>=1.) xx=0.99
!    x3=xx**3.
!    gnu=(-3.*x3-xx**1.5*sqrt(x3+8.))/(2.*(x3-1))
!  elseif (mom==2) then
!    xx=cx**(-1./3.)*mx/m3**(2./3.)
!    if (xx>=1.) xx=0.99
!    x3=xx**3.
!    gnu=(-4.*x3-sqrt(8.*x3+1.)+1)/(2.*(x3-1))
!  elseif (mom==4) then
!    xx=cx**(1./3.)*mx/m3**(4./3.)
!    if (xx<=1.) xx=1.01
!    x3=xx**3.
!    xa=(45*x3**2.+27*x3+(0.,1.)*sqrt(3.)*xx**3.)**(1./3.)
!    gnu=real(-1.*(x3-3.)/(x3-1.)+xa/(3.**(2./3.)*(x3-1.))+(3.*x3**2.+33.*x3)/(3.**(4./3.)*(x3-1.)*xa))
!  elseif (mom==6) then
!    xx=cx*mx/m3**(2.)
!    if (xx<=1.) xx=1.01
!    xa=(729.*xx**2.+sqrt(4.*(-3.*xx**2.-75.*xx-3.)**3.+(729.*xx**2.+729.*xx)**2.)+729.*xx)**(1./3.)
!    gnu=-1.*(xx-4.)/(xx-1.)+xa/(3.*2.**(1./3.)*(xx-1.))+2.**(1./3.)*(3.*xx**2.+75*xx+3)/(3.*xa*(xx-1.))
!  endif
END SUBROUTINE init_dist_sbm

!---------------------------------------------------------------------------------
subroutine micro_proc_tau(tempk,qv,ffcd_mass2d,ffcd_num2d)
use parameters, only: nz, nx, dt, max_nbins
use common_physics, only: qsaturation
Use diagnostics, only: save_dg, i_dgtime, save_binData
Use switches, only: l_sediment, mphys_var, l_fix_aerosols
Use switches_bin
Use namelists, only: dosedimentation, docollisions, docondensation, &
                     donucleation, dobreakup, l_coll_coal, l_break, mp_proc_dg
use column_variables, only: dtheta_adv, dtheta_div, dqv_adv, dqv_div, dss_adv, &
                            dss_div, daerosol_adv, daerosol_div, &
                            dhydrometeors_adv, dhydrometeors_div, exner, &
                            dtheta_mphys,dqv_mphys, daerosol_mphys, &
                            dhydrometeors_mphys, dss_mphys, ss, aerosol, &
                            w_half, dz_half, rho, dz, theta, hydrometeors
use mphys_tau_bin_declare, only: JMINP, JMAXP, LK, ICDKG_BIN, ICDNC_BIN, KKP,&
                                 NQP, IAERO_BIN, ICDKG_BIN, ICDNC_BIN, iqv, &
                                 iqss, ln2, nqp, rprefrcp, prefrcp, prefn, dzn, &
                                 rhon, rdz_on_rhon, tref, IMICROBIN, IMICROBIN, &
                                 CCNNEWTOT, CCNNEWAVG, CCNORIGTOT, eps, &
                                 AMKORIG, ANKORIG, QSATPW, JJP, CCNOLD, CCN, &
                                 AMKOLD, ANKOLD, DS, DUS, CDNCEVAP, am1_diag, &
                                 an1_diag, AM0, AN0, AMN, IINHOM_mix, DUS1, &
                                 tevap_bb, tevap_squires, t_mix, t_mix_approx, &
                                 da_no, da_no_rat, r_int, r_bar, ssat, evap_a, &
                                 DIEMC, DIENC, DIEMD, DIEND, ANK, AMK, AN1OLD, &
                                 AM1OLD, AN0, AM0, AN1, AM1, AMKCC, ANKCC,AN2, &
                                 DG1, SG1, dcrit, dqn_act, XK, IRAINBIN, &
                                 IMICROBIN, rmass_cw, QL_SED, QLN_SED, dD, &
                                 xkk1, xkmean, IRAINP, lk_cloud, CCNORIGAVG, &
                                 dqn_reg, l_dodgs, dth_dt, dq_dt
use namelists, only: aero_N_init,l_advect,l_diverge,ampORbin
use micro_prm, only: col, qindices, q_lem, th_lem, sq_lem,sth_lem, w_lem, &
                     ffcdprev_mass, ffcdprev_num
use physconst, only : p0, this_r_on_cp=>r_on_cp, pi
use module_mp_tau_bin, only: GET_FORCING, COND_new, EVAP_new, REBIN, SXY, &
                             SCONC, BREAK, MICROcheck, REFFCALC, BIN_SEDIMENT, &
                             XACT

implicit none

!generic 1d and 2d KiD arrays for diagnostics
real, dimension(KKP) :: field ! field for 1-D KiD diagnostics
real, dimension(KKP,JMINP:JMAXP):: field_2d ! field for 2-D KiD diagnostics

!variables for fixed radiative flux
!CALL RAD_MOD
 REAL ::                                                          &
&  RADFAC                                                          &
            ! factor to change from flux to heating rate
& ,QL_1D                                                           &
            ! q_l
& ,SWF_NET                                                         &
            ! Net SW flux (W/m2)
& ,SWF_DN                                                          &
            ! Downward SW flux (W/m2)
& ,SWF_UP    ! Upward SW flux (W/m2)
REAL, DIMENSION(jjp) ::                                           &
&  Z_BL_TOP  ! BL top used in radiation (m)
REAL, DIMENSION(jjp,kkp) :: fnt_lw, sth_lw
!     variables for fixed cooling rate
REAL :: cool_kday,sth_lwmax,qttol
INTEGER,DIMENSION(JMINP:JMAXP) ::  K_BL_TOP
! Subprogram arguments
! IN
INTEGER :: I
REAL, DIMENSION(JMINP:JMAXP,KKP) ::                       &
&    TH       ! potential temperature perturbation
REAL, DIMENSION(JMINP:JMAXP,KKP,NQP) ::                   &
&    Q        ! moisture fields

! INOUT
REAL, DIMENSION(JMINP:JMAXP,KKP) ::                    &
&    STH      ! potential temperature tendency
REAL, DIMENSION(JMINP:JMAXP,KKP,NQP) ::                &
&    SQ        ! moisture fields' tendency
!
integer ijj,ikk
!
!END of RADIATION DECLARATIONS
!
! Local variables
!

REAL,DIMENSION(JMINP:JMAXP,KKP) :: DQLDT,DNQLDT,QLOLD,QLNEW           &
&                             , NQLOLD,NQLNEW, RAINOLD,RAINNOOLD   &
&                             ,RAINNEW,RAINNONEW,DRAINDT,DRAINNODT &
&                             , QSATMIX                            &
                                       !saturation mixing ratio
&                             , RH, RH2                            &
                                       !relative humidity
&                             , DTHDT                              &
                                     ! microphysical tendency
                                     ! in potential temperature
&                             , DQVDT                              &
                                     ! microphysical tendency
                                     ! in water vapour
&                             , TBASE                              &
                                     !value of T around which
                                     !saturation functions base
&                             ,QVOLD,QTOT
!
REAL,DIMENSION(JMINP:JMAXP,KKP,LN2) :: CCNfrac_loc
REAL :: THOLD                                                     &
              !potential temp pertubation
& ,      THNEW                                                     &
              !potential temp pertubation after bin microphysics
& ,      PMB                                                       &
              !reference pressure in millibars
& ,      vapour                                                    &
               !total water vapour for the whole grid
& ,      water                                                     &
              ! total water for the whole grid
& , u0, T_OLD, totloss, dccndt,TH2T,qln_tot
REAL :: TOTCCNNUC, TOTEVAP,totevap2,TOTCCNREG,fein_ccnfrac
REAL, DIMENSION(JMINP:JMAXP, KKP) :: cell_ccnreg
REAL, DIMENSION(JMINP:JMAXP, KKP) :: t_0, qst_0
REAL :: cloudrop
REAL, DIMENSION(JMINP:JMAXP, KKP) :: ds_0_temp ! supersat @ beginning of the timestep
REAL :: ds_0      ! dynamic term for supersaturation total trans
REAL :: QmassFLD,QnumFLD
real :: t
!
!local integers
!
INTEGER, DIMENSION(JMINP:JMAXP) :: KQLINV
INTEGER :: KQLMAX
character(2) :: str2

!vars in CLOUDBIN subroutine
REAL,PARAMETER :: AL=597  !latent heat of evaporation and cond
REAL,PARAMETER :: CPBIN=0.24 !specific heat of water vapour
REAL, DIMENSION(JMINP:JMAXP,KKP) :: CCNTOT
REAL :: CN1
REAL DM !the condensated mass calced in EVAP and COND routines
real dm_cloud_evap ! the evaporated mass
real dm_rain_evap ! the evaporated mass
REAL DDDD!supersat percentage
REAL DAMKDT!tendency in bin resolved mass cloud due to microphysics
REAL DANKDT!tendency in bin resolved number cloud due to microhysics
REAL QVNEW!vapour mixing ratio after the nucleation, before code
          !proceeds to evap or cond, after evap or cond, this var
          !QV + change due to nucleation and cond or evap. It is
          !used to calc DQVDT
REAL,DIMENSION(JMINP:JMAXP,KKP) :: TAU                                &
                                   !microphys forcing
&                ,VSW                                              &
&                ,EN                                               &
                      !estimated ETA for end of timestep
&                ,EA   !average ETA for the timestep
REAL :: QST_nuc,ds_force,temp_nuc,mpos
      !added 28/04/08 for new inhom mixing method
REAL,DIMENSION(JMINP:JMAXP,KKP) :: tau_dum
REAL ::  delm

REAL,DIMENSION(JMINP:JMAXP,KKP) :: NCC,MCC
real, dimension(kkp) :: auto_con_mass, d_rmass, d_rnum, d_cmass, d_cnum
real :: rmass_tot_orig, rmass_tot_new
real :: proc_cmass, proc_rmass, proc_tmass, proc_cnum, proc_rnum, proc_tnum
real :: dtcalc
INTEGER :: LT, IT, inhom_evap
INTEGER :: ccn_pos, cloud_pos, count, loop_count

integer :: j, k, l, iq, ih, imom
real :: rdt
real, dimension(nz,nx) :: tempk, qv
real, dimension(nz,nx,max_nbins) :: ffcd_mass2d, ffcd_num2d
integer, parameter :: offset=1 ! 1 = no microphysics on the bottom level



!press & tempk currently not used

rdt = 1./dt
i=1

! prepare for the tau mphys
rprefrcp(1+offset:kkp)=exner(1:kkp-offset,nx) ! I think this is upside-down in LEM
! AH - 04/03/10, line below leads to divide by 0
!      corrected by setting array to 2:kkp. Problem highlighted
!      by Theotonio Pauliquevis
! prefrcp(:)=1./rprefrcp(:)
prefrcp(1+offset:kkp)=1./rprefrcp(1+offset:kkp)

prefn(1+offset:kkp)=p0*exner(1:kkp-offset,nx)**(1./this_r_on_cp)

dzn(1+offset:kkp)=dz_half(1:kkp-offset)

rhon(1+offset:kkp)=rho(1:kkp-offset)
rdz_on_rhon(1+offset:kkp)=1./(dz(1:kkp-offset)*rhon(1+offset:kkp))
! Reference temperature (this is fixed in lem, but
! shouldn't make a difference for microphysics if we
! just set it to be the current profile (i.e. th'=0)
tref(1+offset:kkp)=theta(1:kkp-offset,nx)*exner(1:kkp-offset,nx)




do j = jminp,jmaxp
    q_lem(j,1+offset:kkp,iqv) = qv(1:kkp-offset,j)
    q_lem(j,1+offset:kkp,iqss) = ss(1:kkp-offset,j)
    do iq=1,ln2
        ih=qindices(IAERO_BIN(iq))%ispecies
        imom=qindices(IAERO_BIN(iq))%imoment
        do k=1,nz-offset
            q_lem (j, k+offset, IAERO_BIN(iq)) = aerosol(k,j,ih)%moments(iq,imom)
        end do
    enddo
    do iq=1,lk
        ! mass bins
        ih=qindices(ICDKG_BIN(iq))%ispecies
        imom=qindices(ICDKG_BIN(iq))%imoment
        do k=1,nz-offset
            q_lem (j, k+offset, ICDKG_BIN(iq)) = ffcd_mass2d(k,j,iq)*col!*(pi/6*1000.)
        end do
        ! number bins
        ih=qindices(ICDNC_BIN(iq))%ispecies
        imom=qindices(ICDNC_BIN(iq))%imoment
        do k=1,nz-1
            q_lem(j,k+offset,icdnc_bin(iq)) = ffcd_num2d(k,j,iq)*col
        end do
    enddo
    th_lem (j, :) = 0.0
    w_lem(j,1+offset:kkp)=w_half(1:kkp-offset,j)
enddo

! -------------- set the tendency due to adv and div as input -ahu ------------
!if (l_advect .or. l_diverge) then
  do j=jminp,jmaxp
     sth_lem(j,1+offset:kkp)=dtheta_adv(1:kkp-offset,j)+dtheta_div(1:kkp-offset,j)
     sq_lem(j,1+offset:kkp,iqv)=dqv_adv(1:kkp-offset,j)+dqv_div(1:kkp-offset,j)

     sq_lem(j,1+offset:kkp,iqss)=dss_adv(1:kkp-offset,j)+dss_div(1:kkp-offset,j)

     do iq=1,ln2
        ih=qindices(iaero_bin(iq))%ispecies
        imom=qindices(iaero_bin(iq))%imoment
        do k=1,nz-offset
           sq_lem(j,k+offset,iaero_bin(iq))=(daerosol_adv(k,j,ih)%moments(iq,imom) &
                + daerosol_div(k,j,ih)%moments(iq,imom))
        end do
     enddo
     do iq=1,lk
        ih=qindices(icdkg_bin(iq))%ispecies
        imom=qindices(icdkg_bin(iq))%imoment
        do k=1,nz-offset
            if (ampORbin .eq. 'bin') then
                sq_lem(j,k+offset,icdkg_bin(iq))=(dhydrometeors_adv(k,j,ih)%moments(iq,imom) &
                    + dhydrometeors_div(k,j,ih)%moments(iq,imom))!*(pi/6*1000.)
            elseif (ampORbin .eq. 'amp') then
                sq_lem(j,k+offset,icdkg_bin(iq))=0.
                ! (ffcd_mass2d(k,j,iq) - ffcdprev_mass(k,j,iq))*col/dt
            endif
        end do
        ih=qindices(icdnc_bin(iq))%ispecies
        imom=qindices(icdnc_bin(iq))%imoment
        do k=1,nz-offset
            if (ampORbin .eq. 'bin') then
                sq_lem(j,k+offset,icdnc_bin(iq))=(dhydrometeors_adv(k,j,ih)%moments(iq,imom) &
                    + dhydrometeors_div(k,j,ih)%moments(iq,imom))
            elseif (ampORbin .eq. 'amp') then
                sq_lem(j,k+offset,icdnc_bin(iq))=0.
                ! (ffcd_num2d(k,j,iq) - ffcdprev_num(k,j,iq))*col/dt
            endif
        end do
     end do
   end do

   DO k = 1+offset, nz
      DO j = jminp,jmaxp
         DO IQ = 1, LK
            CALL ADVECTcheck(j,k,iq,DT,Q_lem(j,k,ICDKG_BIN(iq)),              &
  &                 Q_lem(j,k,ICDNC_BIN(iq)),SQ_lem(j,k,ICDKG_BIN(iq)),  &
  &                 SQ_lem(j,k,ICDNC_BIN(iq)))
         ENDDO
      ENDDO
   ENDDO

!end if

!print*, q_lem(j,1+offset:kkp,iqss)

! ------------------------------ mphys starts ----------------------------------

TH = th_lem
Q = q_lem
STH = sth_lem
SQ = sq_lem

IF(IMICROBIN == 1.AND.IRAINP == 0) THEN

!Calculate the domain total bin resolved CCN number for regen
!when Ln2 >  3

TOTCCNNUC = 0.0
DO IQ = 1, Ln2
   CCNNEWTOT(IQ) = 0.0
   CCNNEWAVG(IQ) = 0.0
   DO K = 2, KKP
      DO J = JMINP, JMAXP
         CCNNEWTOT(IQ) =  CCNNEWTOT(IQ) +(Q(J,K,IAERO_BIN(IQ))+    &
              (SQ(J,K,IAERO_BIN(IQ))*DT))
      ENDDO
   ENDDO
   TOTCCNNUC = TOTCCNNUC + (CCNORIGTOT(IQ) - CCNNEWTOT(IQ))
ENDDO

ENDIF

eps = epsilon(1.0d0)

do K = 2, KKP
do J = JMINP,JMAXP
   DO IQ=1,LK
!AMKORIG and ANKORIG for use in positivity check of mass and number
!change due to bin microphysics
      AMKORIG(J,K,IQ)=Q(J,K,ICDKG_BIN(IQ)) + &
           (SQ(J,K,ICDKG_BIN(IQ))*DT)
      ANKORIG(J,K,IQ)=Q(J,K,ICDNC_BIN(IQ)) + &
           (SQ(J,K,ICDNC_BIN(IQ))*DT)

      if (AMKORIG(J,K,IQ) < eps .or.         &
           ANKORIG(J,K,IQ) < eps .or.        &
           AMKORIG(J,K,IQ) > ANKORIG(J,K,IQ)) then ! a check for large mass? -ahu
         AMKORIG(J,K,IQ) = 0.0
         ANKORIG(J,K,IQ) = 0.0
      endif

   ENDDO

ENDDO
ENDDO

if (dosedimentation) then
   CALL BIN_SEDIMENT(I,DT,AMKORIG,ANKORIG,Q,SQ,RDT)

   if (mp_proc_dg) then
      do k=2,kkp
         do j=jminp,jmaxp
            proc_tmass=0.0
            proc_tnum=0.0
            do l=1,lk
                proc_tmass=proc_tmass+ql_sed(j,k,l)
                proc_tnum=proc_tnum+qln_sed(j,k,l)
            enddo
   
            call save_dgproc(proc_tmass/dt,proc_tnum/dt, &
                   'dm_sed', 'dn_sed', k, j)
         enddo 
      enddo
   endif

endif                    ! sedimentation calculation


totevap = 0.0
totevap2 = 0.0
totccnreg = 0.0
DS_0 = 0.0
QLOLD(:,:) = 0.0
NQLOLD(:,:) = 0.0
QLNEW(:,:) = 0.0
NQLNEW(:,:) = 0.0

DO K=2,KKP
    PMB = 0.01*PREFN(K)
    DO J=JMINP,JMAXP
    ! 1. set base values
    ! a) calc ss at beginning of DT for Get_forcing
    T_0(j,k) = TREF(K) + TH(j,k)*RPREFRCP(K)
    QST_0(j,k) = QSATURATION(t_0(j,k),PMB)
    DS_0 = q(j,k,IQSS) + (SQ(J,K,IQSS)*DT)

    ! b) Set base values after dynamics for use in microphys
    THOLD=TH(J,K) + (STH(J,K)*DT)
    QVOLD(J,K)=Q(J,K,IQV) + (SQ(J,K,IQV)*DT)

    DO IQ = 1, LK
       QLOLD(J,K)=QLOLD(J,K)+(Q(J,K,ICDKG_BIN(IQ))               &
        +(SQ(J,K,ICDKG_BIN(IQ))*DT))
       NQLOLD(J,K)=NQLOLD(J,K)+(Q(J,K,ICDNC_BIN(IQ))               &
        +(SQ(J,K,ICDNC_BIN(IQ))*DT))
    ENDDO

    TBASE(J,K)=TREF(K) + THOLD*RPREFRCP(K)

    !
    ! 2. calculate saturation functions
    !
    QSATPW(J,K)  = QSATURATION(TBASE(J,K),PMB)

    ! 3. call cloudbin to calculate the bin microphysics (i.e.
    !   nucleation, condensation, evaporation, collection, breakup)

    ! CALL CLOUDBIN(I,J,K,Q,SQ,AMKORIG,ANKORIG,QSATPW,RH,TBASE,TREF,    &
    !           DQVDT(J,K),DT,RDT,PMB,QVOLD(J,K),QLOLD(J,K),    &
    !           totevap,totccnreg,DS_0)

    ! Set values of bins before microphysics, these values do not include the
    ! effect of dynamics from this timestep. I do not think this is needed
    ! as all I want is a tendency due ti microphysics. The tendency for aerosol
    ! is calculated in SUBNUC, as this is where nucleation is calculated

    if (jjp == 1) then
       call save_dg(k,rhon(k),'density', i_dgtime, units='kg/m3',dim='z')
    endif

    QVNEW = QVOLD(J,K)
    DO L=1,LN2
      CCNOLD(J,K,L)=Q(J,K,IAERO_BIN(L)) + (SQ(J,K,IAERO_BIN(L))*DT)
      CCN(J,K,L) = CCNOLD(J,K,L) !this line is required so that
      !CCN does not equal zero if DS < 0, i.e. for calc of aerosol
                              !tendency
    ENDDO

    DO L=1,LK
      AMKOLD(J,K,L)=AMKORIG(J,K,L)
      ANKOLD(J,K,L)=ANKORIG(J,K,L)
    END DO



    !****************************************************************
    !          DS() IS THE SUPERSATURATED QUANTITY
    !****************************************************************
    DS(J,K)=QVNEW-QSATPW(J,K)
    DUS(J,K)=DS(J,K)/QSATPW(J,K)
    QST_nuc = QSATPW(J,K) !used in activation calculation
    TEMP_NUC = TBASE(J,K)
    DS_FORCE = 0.0

    IF (DT.GT.4.0) THEN
     LT = CEILING(DT/1.0)
     DTcalc = DT/REAL(LT)
    ELSE
     LT=1
     DTcalc=DT
    ENDIF

    CDNCEVAP(j,k) = 0.0
    am1_diag(j,k) = 0.0
    an1_diag(j,k) = 0.0

    !AH 0410 - moving all unit conversion outside
    !        the calls for get_force, cond and evap_new
    !        to minimise precision errors, also don't change
    !        amkold or ankold
    DO L=1,LK
     IF(AMKOLD(J,K,L) > eps .AND.       &
          ANKOLD(J,K,L) > eps )THEN
        am0(l)=(amkold(j,k,l)*Rhon(k))/1.E3
        an0(l)=(ankold(j,k,l)*Rhon(k))/1.E6
        AMN(L)=am0(L)/an0(L)
     ELSE
        AMN(L)=0.
        am0(l)=0.
        an0(l)=0.
     ENDIF
     am1_diag(j,k) = am1_diag(j,k)+am0(l)
     an1_diag(j,k) = an1_diag(j,k)+an0(l)
    ENDDO

    DO it = 1,lt

       IF(IINHOM_mix == 0) then
          VSW(J,K) = MIN(1.,-1.E+10*MAX(DS(j,k),DS_0))
          CALL GET_FORCING(J,K,QSATPW,TBASE,DS_0,DS,                      &
               AM0,AN0,DTcalc,TAU,EN,EA,VSW(J,K))

          inhom_evap = 0

       ENDIF

       DM=0.0
       dm_cloud_evap = 0.0
       dm_rain_evap = 0.0
       NCC(J,K)=0.0
       MCC(J,K)=0.0
       DUS1(J,K)=0.0
       tevap_bb(j,k) = 0.0
       tevap_squires(j,k) = 0.0
       t_mix(j,k) = 0.0
       t_mix_approx(J,K) = 0.0
       da_no(J,K) = 0.0
       da_no_rat(J,K) = 0.0
       r_int(j,k) = 0.0
       r_bar(j,k) = 0.0
       ssat(j,k) = 0.0
       evap_a(j,K) = 0.0

       ssat(J,K) = EA(J,K)/QSATPW(J,K)

       DO L=1,LK
          DIEMC(L)=0.0
          DIENC(L)=0.0
          DIEMD(L)=0.0
          DIEND(L)=0.0
          ANK(J,k,L) = 0.0
          AMK(j,k,l) = 0.0
       ENDDO
       IF (IINHOM_MIX == 0) then
          DS_force = tau(j,k)
          tau_dum(j,k) = tau(j,k)
       endif

! 3.1 cond/evap

       IF (DS_force > eps .and. docondensation) THEN
!*****************************************************************
!        CONDENSATION
!     COND RECEIVES MKOLD,NKOLD RETURNS MK,NK
!****************************************************************
          AN1OLD(J,K) = 0.0
          AM1OLD(j,k) = 0.0
          DO L=1,LK
             AN1OLD(J,K) =AN1OLD(J,K)+AN0(l)
             AM1OLD(j,K) = AM1OLD(j,K)+am0(l)
          ENDDO


!***************************************************************
          CALL COND_new(J,K,DM,TBASE,QSATPW,RH,Q,SQ,DT,RDT,PMB,QVNEW,      &
               TAU_dum,it,LT)
!***************************************************************

!!!!!!!!!! AMK and ANK now have weird units !!!!!!!!!! -ahu

          CALL REBIN(J,K)

          AN1(J,K) = 0.0
          AM1(j,k) = 0.0
          DO L=1,LK
             AN1(J,K) = AN1(J,K) + (ANK(J,K,L))
             AM1(J,K) = AM1(j,k) + (AMK(j,k,l))
          ENDDO

          IF(ABS((AN1(J,K))-(AN1OLD(J,K))) >  1.) THEN
             print *, 'cond conservation prob after rebin'
             print *, 'AN1OLD', AN1OLD(j,k),k,j
             print *, 'AN1', AN1(j,k),k,j
          ENDIF

       ELSEIF(DS_force < -eps .and. docondensation) then !DS_force < 0.0
!****************************************************************
    !                     EVAPORATION
    !     EVAP RECEIVES MKD,NKD RETURNS MK,NK
    !***********************************************************
          AN1OLD(J,K) = 0.0

          DO L=1,LK
             AN1OLD(J,K) =AN1OLD(J,K)+AN0(l)
          ENDDO

          IF (AN1OLD(J,K) > eps) then

             CALL EVAP_new(I,J,K,DM,TBASE,QSATPW,RH,Q,SQ,DT,RDT,PMB,      &
                  QVNEW,TAU_dum,EA,it,LT,inhom_evap,delm)

             CALL REBIN(J,K)
             !
             AN1(J,K) = 0.0
             DO L=1,LK
                AN1(J,K) = AN1(J,K) + (ANK(J,K,L))
             ENDDO

          ELSE
             DO L=1,LK
                AMK(J,K,L)=AM0(L)
                ANK(J,K,L)=AN0(L)
             ENDDO
          ENDIF

        ELSE

           DO L=1,LK
              AMK(J,K,L)=am0(L)
              ANK(J,K,L)=an0(L)
           ENDDO
       ENDIF

         !update fields
         an1old(j,k) = 0.0
         an1(j,k) = 0.0
         DO L=1,LK
            ! n0 and m0 in cjs units for next lt iteration if
            ! needed
            AN0(L)=ank(j,k,L)
            AM0(L)=amk(j,k,L)
            IF(AM0(L) >  eps .AND.                       &
                 AN0(L) > eps)THEN
               AMN(L)=AM0(L)/AN0(L)
            ELSE
               AMN(L)=0.
            ENDIF
    !
    ! AH 0410 - convert mass and number back to kg/kg and #/kg
    !           respectively for source term calc and thermodynamic
    !           updates
    !

            amk(j,k,l)=amk(j,k,l)*1.e3/rhon(k)
            ank(j,k,l)=ank(j,k,l)*1.e6/rhon(k)
!!!!!!!!!! AMK and ANK now have the right units !!!!!!!!!! -ahu

            if(amk(j,k,l) < eps.or.                 &
                 ank(j,k,l) < eps) then
               amk(j,k,l) = 0.0
               ank(j,k,l) = 0.0
            endif
            an1old(j,k) = an1old(j,k) + ankorig(j,k,l)
            an1(j,k) = an1(j,k) + ank(j,k,l)

         ENDDO

    ! AH 0410 - work out total mass and number change due to cond
    !           and evap
         dm = 0.0
         do l = 1,lk
            dm = dm + (amk(j,k,l) - amkorig(j,k,l))
         enddo

    ! save mphys tendency due to condensation for each bin - takes up a lot of space -ahu
      !if (mp_proc_dg) then
      !  call save_dg('bin',k,amk(j,k,:)-amkorig(j,k,:),'dhyd_cond_mass',i_dgtime,units='kg/kg/s',dim='z')
      !  call save_dg('bin',k,ank(j,k,:)-ankorig(j,k,:),'dhyd_cond_num',i_dgtime,units='#/kg/s',dim='z')
      !endif

    ! AH 0410 - Update thermodynamic field once evap and cond
    !           finished
         if ( it .ne. lt )  DS_0 = EN(j,k)

         TBASE(J,K)=TBASE(J,K)+AL*DM/CPBIN
         QVNEW = QVNEW - DM

if (tbase(j,k)<0) then
    print*, 'negative temperature', j,k,tbase(j,k)
    print*, 'dm',dm
    stop
endif

         QSATPW(J,K)=QSATURATION(TBASE(J,K),PMB)
         DS(J,K)=QVNEW-QSATPW(J,K)
    ENDDO               !end of iteration over LT

    !subtract total number after evaporation, from total number before evap
    !this value tells the number of droplets that have evaporated completely
    !
      IF (AN1(J,K) >  AN1OLD(J,K)) THEN
         CDNCEVAP(J,K) = 0.0
      ELSE
         CDNCEVAP(J,K) = (AN1OLD(J,K) - (AN1(J,K)))
      ENDIF
      
      ! this part is probably redundant, will delete later if so. -ahu
      if (mp_proc_dg .and. docondensation) then
         if (jjp > 1) then
            !    diagnostics for total cond/evap rate liquid water (change in mass)
            call save_dg(k,j,dm/dt,'dm_ce', i_dgtime, &
                 units='kg/kg/s',dim='z,x')
            !    cond/evap rate of cloud
            do l = 1,lk_cloud
               dm_cloud_evap = dm_cloud_evap + (amk(j,k,l) - amkorig(j,k,l))
            enddo
            call save_dg(k,j,dm_cloud_evap/dt,'dm_cloud_ce', i_dgtime, &
                 units='kg/kg/s',dim='z,x')
            !     cond/evap rate of rain
            do l = lk_cloud+1, lk
               dm_rain_evap = dm_rain_evap + (amk(j,k,l) - amkorig(j,k,l))
            enddo
            call save_dg(k,j,dm_rain_evap/dt,'dm_rain_ce', i_dgtime, &
                 units='kg/kg/s',dim='z,x')
         else
            !    diagnostics for total cond/evap rate liquid water (change in mass)
            call save_dg(k,dm/dt,'dm_ce', i_dgtime, &
                 units='kg/kg/s',dim='z')
            !    cond/evap rate of cloud
            do l = 1,lk_cloud
               dm_cloud_evap = dm_cloud_evap + (amk(j,k,l) - amkorig(j,k,l))
            enddo
            call save_dg(k,dm_cloud_evap/dt,'dm_cloud_ce', i_dgtime, &
                 units='kg/kg/s',dim='z')
            !     cond/evap rate of rain
            do l = lk_cloud+1, lk
               dm_rain_evap = dm_rain_evap + (amk(j,k,l) - amkorig(j,k,l))
            enddo
            call save_dg(k,dm_rain_evap/dt,'dm_rain_ce', i_dgtime, &
                 units='kg/kg/s',dim='z')
         endif
      endif
! 3.2 do activation after cond/evap, using updated supersat for timestep
if (donucleation) then
    EA(J,K) = DS(J,K)

    IF(EA(J,K) >  0.0) THEN
      AN1OLD(J,K) = 0.0
      AM1OLD(J,K) = 0.0
      CCNTOT(J,K) = 0.0
      DO L = 1, LK
        AN1OLD(J,K)  = AN1OLD(J,K) +ANK(J,K,L)
        AM1OLD(J,K)  = AM1OLD(J,K) +AMK(J,K,L)
      ENDDO

      DO L = 1, LN2
        CCNTOT(J,K) =  CCNTOT(J,K) + CCNOLD(J,K,L)
      ENDDO

      AN1(J,K) = 0.0
      AM1(J,K) = 0.0
      amkcc(1) = 0.0
      ankcc(1) = 0.0

      DO L = 1, LK
        AN1(J,K) = AN1(J,K) + ANK(J,K,L)
      ENDDO

      AM1(J,K) = AMK(J,K,1)
      DDDD=(EA(J,K)/QSATPW(J,K))*100

      ccn_pos = 1
      cloud_pos = 1

      CN1 = CCN(j,k,ccn_pos)/rhon(k) ! convert to /kg

      if (.not. l_fix_aerosols) then
         AN2(J,K) = MAX(0.,                                            &
    &       XACT(tbase(j,k),DDDD,DG1,SG1,dcrit(j,k))*                  &
            CN1)

         !adjust CCN to show the removal of CCN by activation
          CCN(j,k,ccn_pos) = (CN1 - AN2(J,K))*rhon(k) ! convert to /m^3

      else
         AN2(J,K) = MAX(0.,                                            &
    &       XACT(tbase(j,k),DDDD,DG1,SG1,dcrit(j,k))*                  &
            CN1-AN1(J,K))

    !     &       CCN(J,K,ccn_pos)-AN1(J,K))

      endif

      dqn_act(J,K) = AN2(J,K)/dt

      ANK(J,K,1) = ANK(J,K,1) + AN2(J,K)
      ! the 0.25 factor attempts to accomodate for the fact
      ! that not all CCN will grow to the size of the 1st bin

      AMK(J,K,1) = AMK(J,K,1) + AN2(J,K)*XK(1)*0.25

      AMK(J,K,1) = MAX(ANK(J,K,1)*XK(1)*1.01, AMK(J,K,1))

      MCC(J,K) = AMK(J,K,1) - AM1(J,K)
      AM1(j,k) = 0.0
      DO L = 1, LK
        AM1(J,K) = AM1(J,K) + AMK(J,K,L)
      ENDDO

    !calculate the new termodynamic fields following activation
      TBASE(J,K)=TBASE(J,K)+AL/CPBIN*MCC(J,K)
      QVNEW = QVNEW - MCC(J,K)
      QSATPW(J,K)=QSATURATION(TBASE(J,K),PMB)
    ENDIF                          ! end of do activation

    if (mp_proc_dg) then
       if (nx==1) then
          call save_dg(k,MCC(j,k),'dm_nuc',i_dgtime,units='kg/kg/s',dim='z')
       else
          call save_dg(k,j,MCC(j,k),'dm_nuc',i_dgtime,units='kg/kg/s',dim='z')
       endif
    endif
endif

    AN1(j,k) = 0.0
    DO L = 1, LK
      AN1(J,K) = AN1(J,K) + ANK(J,K,L)
    ENDDO

    ! AH 0410 - to update supersat for advection. The update of ss is
    !           performed here not in stepfields (code commented out in
    !           stepfields)
    !
    DS(J,K) =  QVNEW - QSATPW(J,K)
    q(j,k,iqss) = DS(j,k)

!*************************************************************
!        UPDATING CHANGE IN MASS CAUSED BY MICROPHYSICS
!             AFTER  CONDENSATION-EVAPORATION
!*************************************************************

        DO L=1,LK
            DIEMC(L)=AMK(J,K,L)-AMKORIG(J,K,L)!AMKORIG is the mass prior to
                !bin micro, therefore DIEMC (old variable name) is
                !mass resulting from bin micro
            DIENC(L)=ANK(J,K,L)-ANKORIG(J,K,L)!ANKORIG is the number prior to
                !bin micro, therefore DIENC (old variable name) is
                !number resulting from bin micro
        ENDDO

        if (IRAINBIN == 1.AND.IMICROBIN == 1) then

! 3.3
!************************************************************
!            COLLECTION + BREAKUP
!************************************************************
!
! AH 0410 - The collision-coalescence code uses the original format
!           of the TAU model (so is an older version than that used
!           in RAMS). AMKORIG and ANKORIG, which are the same as
!           AMKOLD and ANKOLD) are used to initialise calculations in
!           SXY, SCONC and BREAK. Unit conversions are performed
!           in the respective subroutines (this is not neccesary and
!           should be moved out so that AN0 and AM0 are used
!
            rmass_cw(k) = 0.0

            loop_count = 0
            count = 0
            DO L=1,LK
              AMK(J,K,L)=AMKORIG(J,K,L)
              ANK(J,K,L)=ANKORIG(J,K,L)
            ENDDO

            AM1(J,K)=0.0
            AN1(J,K)=0.0
            DO L=1,LK
              AM1(J,K)=AM1(J,K)+AMKORIG(J,K,L)
              AN1(J,K)=AN1(J,K)+ANKORIG(J,K,L)
            ENDDO


            IF (AM1(J,K) >  lk*eps .AND. AN1(J,K) > lk*eps ) THEN

               if (docollisions) then


                  CALL SXY(J,K,DT)
                  CALL SCONC(J,K,DT)

                  if (l_break .or. dobreakup) then
                     CALL BREAK(J,K,DT)
                  endif

               endif

!*************************************************************
!        UPDATING CHANGE IN MASS CAUSED BY MICROPHYSICS
!             AFTER  COLLECTION & BREAKUP
!*************************************************************
               DO L=1,LK
                  DIEMD(L)=AMK(J,K,L)-AMKORIG(J,K,L)
                  DIEND(L)=ANK(J,K,L)-ANKORIG(J,K,L)
               ENDDO

               if (mp_proc_dg .and. docollisions) then

                  d_cmass(k) = 0.0
                  d_cnum(k) = 0.0
                  d_rnum(k) = 0.0
                  d_rmass(k) = 0.0

                  do l = lk_cloud+1, lk
                     d_rmass(k) = d_rmass(k) + DIEMD(L)
                     d_rnum(k) = d_rnum(k) + DIEND(L)
                  enddo

                  do l = 1,lk_cloud
                     d_cmass(k) = d_cmass(k) + DIEMD(L)
                     d_cnum(k) = d_cnum(k) + DIEND(L)
                  enddo

                  if (jjp == 1) then
                     call save_dg(k,d_cmass(k)/dt,'dm_cloud_coll', i_dgtime, &
                          units='kg/kg/s',dim='z')
                     call save_dg(k,d_cnum(k)/dt,'dn_cloud_coll', i_dgtime, &
                          units='#/kg/s',dim='z')
                     call save_dg(k,d_rmass(k)/dt,'dm_rain_coll', i_dgtime, &
                          units='kg/kg/s',dim='z')
                     call save_dg(k,d_rnum(k)/dt,'dn_rain_coll', i_dgtime, &
                          units='#/kg/s',dim='z')
                  else

                     call save_dg(k,d_cmass(k)/dt,'dm_cloud_coll', i_dgtime, &
                          units='kg/kg/s',dim='z,x')
                     call save_dg(k,d_cnum(k)/dt,'dn_cloud_coll', i_dgtime, &
                          units='#/kg/s',dim='z,x')
                     call save_dg(k,d_rmass(k)/dt,'dm_rain_coll', i_dgtime, &
                          units='kg/kg/s',dim='z,x')
                     call save_dg(k,d_rnum(k)/dt,'dn_rain_coll', i_dgtime, &
                          units='#/kg/s',dim='z,x')
                  endif

               endif

            ELSE
               DO L=1,LK
                  DIEMD(L) = 0.0
                  DIEND(L)= 0.0
               ENDDO
           ENDIF
        ENDIF

!***********************************************************
!     UPDATING MASS & CONCENTRATION AFTER MICROPHYSICS
!***********************************************************
!      AM1(J,K)=0.0
!      AN1(J,K)=0.0

        DO L=1,LK
            IF(IRAINBIN == 1.AND.IMICROBIN == 1) THEN
                !add microphysical change resulting from  nuc,cond/evap,
                !collision-coalscence, breakup and sedimentation
                AMK(J,K,L)=DIEMC(L)+AMKORIG(J,K,L)+DIEMD(L)+QL_SED(J,K,L)
                ANK(J,K,L)=DIENC(L)+ANKORIG(J,K,L)+DIEND(L)+QLN_SED(J,K,L)
            ELSE!only add microphysical change resulting from  nuc,cond/evap
                AMK(J,K,L)=DIEMC(L)+AMKORIG(J,K,L)
                ANK(J,K,L)=DIENC(L)+ANKORIG(J,K,L)
            ENDIF
        ENDDO
!
!     Calculate the change in "rain bin" mass and number
!     (i.e. bin greater than  100 microns diameter)
!      if (j == 0) then
!       rmass_tot_orig = 0.0
!       rmass_tot_new = 0.0
!       !  d_rmass(k) = 0.0
!       auto_con_mass(k) = 0.0
!
!       do l=lk_cloud+1,lk
!           rmass_tot_orig = rmass_tot_orig + AMKORIG(J,K,L)
!           rmass_tot_new = rmass_tot_new + AMK(J,K,L)
!       enddo
!
!       ! d_rmass(k) = rmass_tot_new - rmass_tot_orig
!
!       auto_con_mass(k) =  d_rmass(k) - rmass_cw(k)
!
!       if (jjp == 1) then
!          call save_dg(k,auto_con_mass(k)/dt,'drain_auto', i_dgtime, &
!               units='kg/kg/s',dim='z')
!
!          call save_dg(k,rmass_cw(k)/dt,'drain_rain', i_dgtime, &
!               units='kg/kg/s',dim='z')
!       else
!          call save_dg(k,j,auto_con_mass(k)/dt,'drain_auto', i_dgtime, &
!               units='kg/kg/s',dim='z,x')
!
!          call save_dg(k,j,rmass_cw(k)/dt,'drain_rain', i_dgtime, &
!               units='kg/kg/s',dim='z,x')
!
!       endif

!      endif

!**********************************************************
!     Update mass and conc, tendencies due to micro and end
!**********************************************************
!First update change in aerosol (if aerosol is not fixed)

        if (.not. l_fix_aerosols) then
           DO L = 1,LN2
              DCCNDT=(CCN(J,K,L)-CCNOLD(J,K,L))*RDT
              IF ((Q(J,K,IAERO_BIN(L))+(SQ(J,K,IAERO_BIN(L))+DCCNDT)*DT)     &
        &         <  0.0) THEN
                 DCCNDT=-0.9999*CCNOLD(J,K,L)*RDT
              ENDIF
              SQ(J,K,IAERO_BIN(L))=SQ(J,K,IAERO_BIN(L))+DCCNDT
           ENDDO
        endif
        DO L=1,LK
!AH 0410 - Calculate the change in mass and number due to microphysics i.e.
!          nucleation, cond, evap, collection and sedimentation
!
            IF(IRAINBIN == 1.AND.IMICROBIN == 1) THEN
            ! sum the effects of cond/evap+coll/coal+sedi
            ! on mass and number
             DAMKDT=(DIEMC(L)+DIEMD(L)+QL_SED(J,K,L)) * RDT
             DANKDT=(DIENC(L)+DIEND(L)+QLN_SED(J,K,L)) * RDT
            ELSE
             ! only use the term from cond/evap (and nuc)
             DAMKDT = DIEMC(L) * RDT
             DANKDT = DIENC(L) * RDT
            ENDIF
!AH 0410 - Call microcheck to update the source fields for bin resolved mass
!          and number with DAMKDT and DANKDT
!
            CALL MICROcheck(J,K,L,DT,AMKORIG(J,K,L), &
               ANKORIG(j,k,l),SQ(J,K,ICDKG_BIN(L)), &
               SQ(J,K,ICDNC_BIN(L)),DAMKDT,DANKDT,RDT,Q(J,K,ICDKG_BIN(L)) &
               ,Q(J,K,ICDNC_BIN(L)))

        ENDDO

        dD(1)=xkk1(1)*2.0*1.e6
        do l = 2, lk
           dD(l)=(xkk1(l)-xkk1(l-1))*2.0*1.e6
        end do

        call save_binData((xkk1(:)*2.0*1.e6),'bins_D_lower', &
             units='microns', longname='lower bound of droplet diameter')
        call save_binData(xk(:),'bins_mass_lower', &
             units='kg', longname='lower bound of droplet mass')
        call save_binData(xkmean(:),'bins', &
             units='kg', longname='mean droplet mass')
        call save_binData(((xkmean(:)*6./3141.59)**(1./3.)*1.e6), &
             'bins_D', units='microns' &
             , longname='mean droplet diameter')
        call save_binData(dD(:), 'dD', units='microns' &
             , longname='width of bin')

! Finally calc the microphys change in QV
        DQVDT(j,k)=(QVNEW-QVOLD(J,K))*RDT

! 4. calculate the change in theta due to bin microphysics
        THNEW=(TBASE(J,K) - TREF(K))*PREFRCP(K)
        DTHDT(J,K)=(THNEW-THOLD)*RDT
        totevap = totevap + CDNCEVAP(J,K)

    ENDDO
ENDDO

!
!regeneration for single prognostic CCN variable
if (.not. l_fix_aerosols) then
    DO IQ = 1, LN2
        CCNNEWTOT(IQ) = 0.0
        CCNNEWAVG(IQ) = 0.0
        DO J = JMINP, JMAXP
            DO K = 2, KKP
                CCNNEWTOT(IQ) =  CCNNEWTOT(IQ) + (Q(J,K,IAERO_BIN(IQ))+   &
                &          (SQ(J,K,IAERO_BIN(IQ))*DT))
            ENDDO
        ENDDO
        TOTCCNNUC = TOTCCNNUC + (CCNORIGTOT(IQ) - CCNNEWTOT(IQ))
        CCNNEWAVG(IQ) = CCNNEWTOT(IQ)/(JJP*(KKP))
    ENDDO
    DO IQ = 1, LN2
        DO J = JMINP, JMAXP
            DO K = 2, KKP
                IF(CCNNEWAVG(IQ) <  CCNORIGAVG(IQ)   &
                    & .AND.CDNCEVAP(J,K) >  0.0) THEN
                    DCCNDT = CDNCEVAP(J,K)*RDT
                    SQ(J,K,IAERO_BIN(IQ)) = SQ(J,K,IAERO_BIN(IQ)) + DCCNDT
                    dqn_reg(J,K) = DCCNDT
                    CDNCEVAP(J,K) = 0.0
                ELSE
                    dqn_reg(j,k) = 0.0
                    CDNCEVAP(J,K) = 0.0
              ENDIF
            ENDDO
        ENDDO
    ENDDO

    if (jjp == 1) then
        field(:) = dqn_reg(jjp,:)
        call save_dg(field,'ccn_reg', i_dgtime, &
             units='#/kg/s',dim='z')
    else
        do K = 2,KKP
            do J = JMINP, JMAXP
            field_2d(k, j) = dqn_reg(j,k)
            enddo
        enddo
        call save_dg(field_2d(1:kkp,1:jjp),'ccn_reg', i_dgtime, &
            units='#/kg/s',dim='z,x')
    endif
endif

if (jjp == 1) then
    field(:) = dqn_act(jjp,:)
    call save_dg(field,'ccn_act', i_dgtime, units='#/kg/s',dim='z')
    else
    do  K = 2,KKP
        do J = JMINP, JMAXP
            field_2d(k, j) = dqn_act(j,k)
        enddo
    enddo
    call save_dg(field_2d(1:kkp,1:jjp), 'ccn_act', i_dgtime, units='#/kg/s',&
        dim='z,x')
endif

DO K = 2,KKP
    DO J = JMINP,JMAXP
        STH(J,K) = STH(J,K) + DTHDT(J,K)
        SQ(J,K,IQV) = SQ(J,K,IQV)+DQVDT(J,K)
        DO IQ=1,LK
        !Call microcheck to update the source fields for bin resolved mass
        !and number, and check that small numbers are not causing erroneous
        !values of mass and number that lead to numerical instabilities

            IF((Q(J,K,ICDKG_BIN(IQ))+(SQ(J,K,ICDKG_BIN(IQ))*DT)) < 0.0 &
            .or.(Q(J,K,ICDNC_BIN(IQ))+(SQ(J,K,ICDNC_BIN(IQ))*DT)) < 0.0) THEN
                QLNEW(J,K)=QLNEW(J,K)
                NQLNEW(J,K)=NQLNEW(J,K)
            ELSE
                QLNEW(J,K)=QLNEW(J,K)+(Q(J,K,ICDKG_BIN(IQ))               &
                              +(SQ(J,K,ICDKG_BIN(IQ))*DT))
                NQLNEW(J,K)=NQLNEW(J,K)+(Q(J,K,ICDNC_BIN(IQ))             &
                              +(SQ(J,K,ICDNC_BIN(IQ))*DT))
            ENDIF
        ENDDO
        DQLDT(J,K)=(QLNEW(J,K)-QLOLD(J,K))*RDT
        DNQLDT(J,K)=(NQLNEW(J,K)-NQLOLD(J,K))*RDT
    ENDDO

ENDDO

! 6. calculate effective radius for radiation
CALL REFFCALC(Q,SQ,DT,RDT)

 !
! 7. update diagnostic fields if necessary
!
DO K = 2, KKP
    DO J = JMINP, JMAXP
        IF(l_dodgs)THEN
            dth_dt(j,k) = DTHDT(J,K)
            dq_dt(j,k,iqv) = DQVDT(J,K)
        ENDIF
    ENDDO
ENDDO

!do j=jminp,jmaxp
!    do k=1,kkp
!        do iq=1,lk
!            if (q(j,k,icdkg_bin(iq))>.1) then
!                print*, 'liquid too massive'
                !print*, 'q', j,k,iq,q(j,k,icdkg_bin(iq))
                !print*, 'sq', j,k,iq,sq(j,k,icdkg_bin(iq))
!                stop
!            endif
!        enddo
!    enddo
!enddo

th_lem = TH
q_lem = Q
sth_lem = STH
sq_lem = SQ

! ------------------------------- mphys ends -----------------------------------

! ---------- sq and sth were tendencies due to mphys + adv + div. -------------
! ---------------------- set them as only due to mphys. -----------------------
do j=jminp,jmaxp
!    if (l_advect .or. l_diverge) then
        sth_lem(j,1+offset:kkp)=sth_lem(j,1+offset:kkp)     &
            - (dtheta_adv(1:kkp-offset,j)+dtheta_div(1:kkp-offset,j))
        sq_lem(j,1+offset:kkp,iqv)=sq_lem(j,1+offset:kkp,iqv)   &
            - (dqv_adv(1:kkp-offset,j)+dqv_div(1:kkp-offset,j))

        do iq=1,ln2
            ih=qindices(iaero_bin(iq))%ispecies
            imom=qindices(iaero_bin(iq))%imoment
            do k=1,nz-offset
                sq_lem(j,k+offset,iaero_bin(iq))=sq_lem(j,k+offset,iaero_bin(iq))       &
                    - (daerosol_adv(k,j,ih)%moments(iq,imom)                &
                    + daerosol_div(k,j,ih)%moments(iq,imom))
            end do
        enddo

        do iq=1,lk
            ih=qindices(icdkg_bin(iq))%ispecies
            imom=qindices(icdkg_bin(iq))%imoment
            do k=1,nz-offset
                if (ampORbin .eq. 'bin') then
                    sq_lem(j,k+offset,icdkg_bin(iq))= sq_lem(j,k+offset,icdkg_bin(iq))      &
                        - (dhydrometeors_adv(k,j,ih)%moments(iq,imom)           &
                        + dhydrometeors_div(k,j,ih)%moments(iq,imom))!*(pi/6*1000.)
                elseif (ampORbin .eq. 'amp') then
                    ! sq_lem(j,k+offset,icdkg_bin(iq)) = sq_lem(j,k+offset,icdkg_bin(iq)) &
                    !     - (ffcd_mass2d(k,j,iq) - ffcdprev_mass(k,j,iq))*col/dt
                endif
            end do
            ih=qindices(icdnc_bin(iq))%ispecies
            imom=qindices(icdnc_bin(iq))%imoment
            do k=1,nz-offset
                if (ampORbin .eq. 'bin') then
                    sq_lem(j,k+offset,icdnc_bin(iq))= sq_lem(j,k+offset,icdnc_bin(iq))      &
                        - (dhydrometeors_adv(k,j,ih)%moments(iq,imom)           &
                        + dhydrometeors_div(k,j,ih)%moments(iq,imom))
                elseif (ampORbin .eq. 'amp') then
                    ! sq_lem(j,k+offset,icdnc_bin(iq)) = sq_lem(j,k+offset,icdnc_bin(iq)) &
                    !     - (ffcd_num2d(k,j,iq) - ffcdprev_num(k,j,iq))*col/dt
                endif
            end do
        end do
!    end if
!   For now set no microphysics on the bottom level - this would be
!   better done by having a subterranian level 0 in column variables
    if (offset==1) sq_lem(j,1,1)=0
    tempk(1:kkp-offset,j) = ((tempk(1:kkp-offset,j)/exner(1:kkp-offset,j)) &
                            + sth_lem(j,1+offset:kkp)*dt)*exner(1:kkp-offset,j)
    qv(1:kkp-offset,j) = qv(1:kkp-offset,j) + sq_lem(j,1+offset:kkp,iqv)*dt

!
!   update supersaturation field here (not in step fields)
!
    ss(1:kkp-offset,j) = q_lem(j,1+offset:kkp,iqss)

!   daerosol_mphys is not yet updated in the main amp routine -ahu
    do iq=1,ln2
        ih=qindices(iaero_bin(iq))%ispecies
        imom=qindices(iaero_bin(iq))%imoment
        do k=1,nz-offset
           daerosol_mphys(k,j,ih)%moments(iq,imom) =              &
               sq_lem(j,k+offset,iaero_bin(iq))
        end do
    enddo

    do iq=1,lk
        do k=1,nz-offset
            ffcd_mass2d(k,j,iq) = ffcd_mass2d(k,j,iq) + sq_lem(j,k+offset,icdkg_bin(iq))*dt/col!/(pi/6*1000.)
        end do
        do k=1,nz-offset
            ffcd_num2d(k,j,iq) = ffcd_num2d(k,j,iq) + sq_lem(j,k+offset,icdnc_bin(iq))*dt/col
        end do
    end do
end do

! check big liquid
do k=1,nz
    do j=jminp,jmaxp
        do l=1,lk
            if (q_lem(j,k,icdkg_bin(l))>.1 .or. sq_lem(j,k,icdkg_bin(l))>.1) then
                !print*, 'liquid too massive'
                !print*, 'k=',k, 'l=',l
                !print*, 'q=',q_lem(j,k,icdkg_bin(l))
                !print*, 'sq=',sq_lem(j,k,icdkg_bin(l))
                !stop
            end if
        end do
    end do
end do

! check NaN
do k=1,nz
    do j=jminp,jmaxp
        do iq=1,nqp
            if (q_lem(j,k,iq) .ne. q_lem(j,k,iq)) then
                print*, 'some NaNs here'
                print*,'q', k,iq
                stop
            end if
            if (sq_lem(j,k,iq) .ne. sq_lem(j,k,iq)) then
                print*, 'some NaNs here'
                print*,'sq', k,iq
                stop
            end if
        enddo
    enddo
enddo

do k=1,nz
    !print*, k,sum(q_lem(j,k,icdkg_bin(:)))
enddo

!print*,'sq', sq_lem(1,31,icdkg_bin(:))/col

end subroutine micro_proc_tau

!-------------------------------------------------------------------
subroutine micro_init_tau
use micro_prm
use module_bin_init
use mphys_tau_bin_declare
implicit none

integer :: j,k

! Set up microphysics species
call set_micro
call bin_init !initialises the cloud bin categories
call data     !reads in and sets the coll-coal kernal

diams=dgmean!DIAM(1:max_nbins)*.01
binmass=xkgmean!xk

DO IQ = 1,LN2
 ih=qindices(IAERO_BIN(iq))%ispecies
 imom=qindices(IAERO_BIN(iq))%imoment
 DO k=1,nz-1
   DO j = JMINP , JMAXP
     CCNORIG(j,k+1,IQ) = aerosol(k,j,ih)%moments(iq,imom)
   ENDDO
 ENDDO
ENDDO

DO k = 1, KKP
 DO j = JMINP, JMAXP
   TOTCCNORIG(j,k) = 0.0
   DO IQ = 1, LN2
     TOTCCNORIG(j,k) = TOTCCNORIG(j,k) + CCNORIG(j,k,IQ)
   ENDDO
 ENDDO
ENDDO
DO IQ = 1, Ln2
  CCNORIGTOT(IQ) = 0.0
  CCNORIGAVG(IQ) = 0.0
   DO k = 1, KKP
     DO j = JMINP, JMAXP
       CCNORIGTOT(IQ) = CCNORIGTOT(IQ) + CCNORIG(j,k,IQ)
     ENDDO
   ENDDO
   CCNORIGAVG(IQ) = CCNORIGTOT(IQ)/(JJP*KKP)
ENDDO

!diams = DIAM(1:max_nbins)*.01 !convert to metric and ditch the last dummy element (?) -ahu


end subroutine micro_init_tau

!----------- sub-subroutine for micro_init_tau ------------
subroutine set_micro
Use parameters, only : num_h_moments, num_h_bins, &
     nspecies, mom_names, h_names, mom_units, max_char_len, &
     num_aero_moments,num_aero_bins, aero_mom_init
Use mphys_tau_bin_declare
Use micro_prm
use namelists, only: ampORbin
implicit none

integer :: i

rdt=1./dt
! vapour

iq=1
iqv=iq

if (num_h_bins(1) >= 11 .or. ampORbin .eq. 'amp') then
  IMICROBIN=1
  IRAINBIN=1

  call qcount(iqss, iq)   ! advected supersat.

!      call qcount(iql, iq)   ! total water for diag

  if (num_aero_moments(1) >= 1) then
     do i = 1,ln2
        call qcount(IAERO_BIN(i),iq) ! aerosol bins
     enddo
  end if
  do i = 1,lk
    call qcount(ICDKG_BIN(i),iq) ! cloud mass bins
  enddo
  do i = 1,lk
    call qcount(ICDNC_BIN(i),iq) ! cloud number bins
  enddo

  nqs=aero_bin+ln2+lk+lk+1

  allocate(qindices(nqs))

  ! Set qindices for later use
  if (num_aero_moments(1) >= 1) then
     do i = 1,ln2
        qindices(IAERO_BIN(i))%ispecies=1
        qindices(IAERO_BIN(i))%imoment=1
     enddo
  end if
  if (num_h_bins(1) >= 4)then ! cloud mass
     do i = 1,lk
        qindices(ICDKG_BIN(i))%ispecies=1
        qindices(ICDKG_BIN(i))%imoment=1
        qindices(ICDNC_BIN(i))%ispecies=1
        qindices(ICDNC_BIN(i))%imoment=2
     enddo
  end if

end if ! bin model selected

end subroutine set_micro

! ------------ sub-subroutine for micro_init_tau ------------
subroutine qcount(var, count)
integer, intent(out) :: var
integer, intent(inout) ::  count

count=count+1
var=count
end subroutine qcount

!-------------------------------------------------------------------
Subroutine init_dist_tau(rxc,gnuc,dnc,rxr,gnur,dnr,ffcd_mass,ffcd_num)

use micro_prm, only:nkr, diams, krdrop,binmass
use parameters, only: max_nbins
use mphys_tau_bin_declare, only: DIAM, NQP, xkgmean, dgmean,xk
use physconst, only: pi, rhoW

implicit none

integer :: kr
real:: rxc,gnuc,dnc,rxr,gnur,dnr
real(8):: n0c,exptermc,n0r,exptermr
real(8), dimension(max_nbins) :: ffcd_mass, ffcd_num


!Setting up a mass distribution, not a number distribution
!So increase gnu by 3

ffcd_mass=0.
ffcd_num=0.

n0c=rxc/gamma(gnuc+3)
n0r=rxr/gamma(gnur+3)

do kr=1,nkr
  if (rxc>0.) then
    exptermc=exp(-1.*diams(kr)/dnc)
    ffcd_mass(kr) = n0c*exptermc*(diams(kr)/dnc)**(gnuc+3)
    ffcd_num(kr) = ffcd_mass(kr)/binmass(kr)!(diams(kr)**3*pi/6.*rhoW)
  endif
  if (rxr>0.) then
     exptermr=exp(-1.*diams(kr)/dnr)
     ffcd_mass(kr) = ffcd_mass(kr) + n0r*exptermr*(diams(kr)/dnr)**(gnur+3)
     ffcd_num(kr) = ffcd_num(kr) + ffcd_mass(kr)/binmass(kr)!(diams(kr)**3*pi/6.*rhoW)
!  else
!     ffcd_mass(krdrop:nkr)=0.
  endif

!If ffcd(kr) is NaN, set to zero
!    if (ffcd(kr).ne.ffcd(kr) .or. ffcd(kr)*0.0.ne.0.0 &
!       .or. ffcd(kr)/ffcd(kr).ne. 1.0) ffcd(kr)=0.
enddo
return
!Diagnose shape parameter
!  if (mom<0) then
!    !Shape parameter is given by mx
!    gnu=mx
!  elseif (mom==1) then
!    xx=cx**(-2./3.)*mx/m3**(1./3.)
!    if (xx>=1.) xx=0.99
!    x3=xx**3.
!    gnu=(-3.*x3-xx**1.5*sqrt(x3+8.))/(2.*(x3-1))
!  elseif (mom==2) then
!    xx=cx**(-1./3.)*mx/m3**(2./3.)
!    if (xx>=1.) xx=0.99
!    x3=xx**3.
!    gnu=(-4.*x3-sqrt(8.*x3+1.)+1)/(2.*(x3-1))
!  elseif (mom==4) then
!    xx=cx**(1./3.)*mx/m3**(4./3.)
!    if (xx<=1.) xx=1.01
!    x3=xx**3.
!    xa=(45*x3**2.+27*x3+(0.,1.)*sqrt(3.)*xx**3.)**(1./3.)
!    gnu=real(-1.*(x3-3.)/(x3-1.)+xa/(3.**(2./3.)*(x3-1.))+(3.*x3**2.+33.*x3)/(3.**(4./3.)*(x3-1.)*xa))
!  elseif (mom==6) then
!    xx=cx*mx/m3**(2.)
!    if (xx<=1.) xx=1.01
!    xa=(729.*xx**2.+sqrt(4.*(-3.*xx**2.-75.*xx-3.)**3.+(729.*xx**2.+729.*xx)**2.)+729.*xx)**(1./3.)
!    gnu=-1.*(xx-4.)/(xx-1.)+xa/(3.*2.**(1./3.)*(xx-1.))+2.**(1./3.)*(3.*xx**2.+75*xx+3)/(3.*xa*(xx-1.))
!  endif
END SUBROUTINE init_dist_tau

! ------------ originally part of mphys_tau_bin ------------
SUBROUTINE  ADVECTcheck(j,k,IQ,DT,ZQmass,ZQnum,Sourcemass,&
     Sourcenum)
!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
use mphys_tau_bin_declare, only: X_BIN, rhon
IMPLICIT NONE

!CALL PRAMETR
!CALL RI
!CALL XD
!CALL GRID1
!local variable
REAL :: QmassFLD, ZQmass, Qmass
REAL :: Qnumfield, ZQnum, Qnum
REAL :: Sourcemass
REAL :: Sourcenum
REAL :: AVG, AVGinit
REAL :: DT,RDT
!loop counters
INTEGER j, k, IQ

RDT = 1.0/DT

!First calculate the new field
QmassFLD=(ZQmass+(DT*Sourcemass))
Qnumfield = (ZQnum+(DT*Sourcenum))
!Change units to microphys units (just for consistency

QmassFLD = (QmassFLD*rhon(k))/1.e3
Qnumfield = (Qnumfield*rhon(k))/1.e6
Sourcemass = (Sourcemass*rhon(k))/1.e3
Sourcenum = (Sourcenum*rhon(k))/1.e6

!calculate the average particle size for the bin

IF(Qnumfield  >  0.0)THEN
  AVG =QmassFLD/Qnumfield
  IF(AVG >  (2.*X_BIN(IQ))) THEN
    sourcenum=((QmassFLD/(2.*X_BIN(IQ)-1.e-20))-                        &
&                             ((ZQnum*rhon(k))/1.e6))*RDT
  ENDIF
  IF(AVG <  X_BIN(IQ).AND.AVG >  0.0) THEN
    sourcenum=((QmassFLD/(X_BIN(IQ)+1.e-20))-                           &
&                             ((ZQnum*rhon(k))/1.e6))*RDT
  ENDIF

ENDIF
!
Sourcemass = (Sourcemass*1.e3)/rhon(k)
Sourcenum  = (Sourcenum*1.e6)/rhon(k)

!do positivity check after normalisation as changing SQ by normalisation
!can lead to negatives
IF((ZQmass+(DT*Sourcemass) <  0.0).or.                            &
&                    (ZQnum+(DT*Sourcenum) <  0.0)) THEN
  Sourcemass = -0.99999 * ZQmass/DT
  Sourcenum =  -0.99999 * ZQnum/DT
ENDIF

IF (ABS(Sourcemass) <  1.e-20.OR. ABS(Sourcenum) <  1.e-20) THEN
    Sourcemass = 0.0
    Sourcenum = 0.0
ENDIF

END subroutine ADVECTcheck


subroutine save_dgproc(varmass, varnum, namemass, namenum, k, j)
   use parameters, only : nx
   use namelists, only: bintype
   use diagnostics, only: save_dg, i_dgtime

   implicit none
   real, intent(in) :: varmass
   real, intent(in), optional :: varnum
   character(*), intent(in) :: namemass
   character(*), intent(in), optional :: namenum
   integer, intent(in) :: k, j

   if (nx==1) then
      call save_dg(k,varmass,namemass,i_dgtime,units='kg/kg/s',dim='z')
      if (present(varnum)) call save_dg(k,varnum,namenum,i_dgtime,units='#/kg/s',dim='z')
   else
      call save_dg(k,j,varmass,namemass,i_dgtime,units='kg/kg/s',dim='z,x')
      if (present(varnum)) call save_dg(k,j,varnum,namenum,i_dgtime,units='#/kg/s',dim='z,x')
   endif

end subroutine 
