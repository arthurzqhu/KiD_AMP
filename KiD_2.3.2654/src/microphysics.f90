!A. Igel - 5/2013 - Adapting for CLOUD model
! This module adapts the SBM from wrf, originally from A. Khain.

Subroutine micro_proc_sbm(press,tempk,qv,fncn,ffcd)

use module_hujisbm
use micro_prm
use parameters, only: nx,nz,dt,aero_N_init,max_nbins
use column_variables, only:z_half

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
            DEL1IN=EW1N/ES1N-1.
            DIV2=EW1N/ES2N        !Saturation ratio wrt ice
            DEL2IN=EW1N/ES2N-1.
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
                 concccn = aero_N_init(1)*1.e-6*(100.*del1in)**bccn
                 concdrop = sum(ff1in*xl)*3.*col
                 IF (concccn>concdrop) THEN
                    ff1in(1)=ff1in(1)+(concccn-concdrop)/(3.0*col*xl(1))
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
                  CALL ONECOND1(TT,QQ,PP,rhocgs &
                          ,VR1,pcgs &
                          ,DEL1IN,DEL2IN,DIV1,DIV2 &
                          ,FF1R,FF1IN,XL,RLEC,RO1BL &
                          ,AA1_MY,BB1_MY,AA2_MY,BB2_MY &
                          ,C1_MEY,C2_MEY &
                          ,COL,DTCOND,ICEMAX,NKR,ISYM1 &
                          ,ISYM2,ISYM3,ISYM4,ISYM5,1,1,1,0.,1.,1)

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
               CALL COAL_BOTT_NEW(FF1R,FF2R,FF3R, &
                    FF4R,FF5R,TT,QQ,PP,rhocgs,dthalf,TCRIT,TTCOAL, &
                    cld2raint,rimecldt, &
                    aggregatet,rimecldsnowt, &
                    rimecldaggrt,rimecldgraut, &
                    rimecldhailt,rain2icet, &
                    rain2snt,rain2agt, &
                    rain2grt,rain2hat)
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
   if(dosedimentation) call falfluxhucm_z(ff1z,vr1z,rhocgs_z,pcgs_z,zcgs_z,dt,1,nz,nkr)
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
!print*, XL
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
real:: rxc,gnuc,dnc,rxr,gnur,dnr
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
subroutine micro_proc_tau(thpert,qv,ss,ffcd_mass2d,ffcd_num2d)
use parameters, only: nz, nx, dt, max_nbins
use column_variables, only: dtheta_adv, dtheta_div, dqv_adv, dqv_div, dss_adv, &
                            dss_div, daerosol_adv, daerosol_div, &
                            dhydrometeors_adv, dhydrometeors_div, exner, &
                            dtheta_mphys,dqv_mphys, daerosol_mphys, dhydrometeors_mphys, &
                            dss_mphys
use mphys_tau_bin_declare, only: JMINP, JMAXP, LK, ICDKG_BIN, ICDNC_BIN, KKP,&
                                 NQP, IAERO_BIN, ICDKG_BIN, ICDNC_BIN, iqv, &
                                 iqss, ln2, nqp
use module_mp_tau_bin, only: tau_bin
use namelists, only: aero_N_init,l_advect,l_diverge
use micro_prm, only: col, qindices, q_lem, th_lem, sq_lem,sth_lem, w_lem

integer :: j, k, iq, ih, imom
real :: rdt
real, dimension(nz,nx) :: thpert, qv,ss
real, dimension(nz,nx,max_nbins) :: ffcd_mass2d, ffcd_num2d

!press & tempk currently not used

rdt = 1./dt

! prepare for the tau mphys
do k = 1, nz
  do j = jminp, jmaxp
    q_lem(j,k,iqv) = qv(k,j)
    q_lem(j,k,iqss) = ss(k,j)
    th_lem(j,k) = thpert(k,j)
    do iq = 1,ln2
      q_lem(j,k,iaero_bin(iq)) = aero_N_init(1) ! this (1) might be problematic in the future
                                                ! but will leave it here for now
      !aero_N a constant in TAU standalone, will worry about it later if it's not -ahu
    end do
    do iq=1,lk
      q_lem(j,k,icdkg_bin(iq)) = ffcd_mass2d(k,j,iq)*col
      q_lem(j,k,icdnc_bin(iq)) = ffcd_num2d(k,j,iq)*col
    end do
  end do
end do

if (l_advect .or. l_diverge) then
  do j=jminp,jmaxp
     sth_lem(j,2:kkp)=dtheta_adv(1:kkp-1,j)+dtheta_div(1:kkp-1,j)
     sq_lem(j,2:kkp,iqv)=dqv_adv(1:kkp-1,j)+dqv_div(1:kkp-1,j)

     sq_lem(j,2:kkp,iqss)=dss_adv(1:kkp-1,j)+dss_div(1:kkp-1,j)

     do iq=1,ln2
        ih=qindices(iaero_bin(iq))%ispecies
        imom=qindices(iaero_bin(iq))%imoment
        do k=1,nz-1
           sq_lem(j,k+1,iaero_bin(iq))=(daerosol_adv(k,j,ih)%moments(iq,imom) &
                + daerosol_div(k,j,ih)%moments(iq,imom))
        end do
     enddo
     do iq=1,lk
        ih=qindices(icdkg_bin(iq))%ispecies
        imom=qindices(icdkg_bin(iq))%imoment
        do k=1,nz-1
           sq_lem(j,k+1,icdkg_bin(iq))=dhydrometeors_adv(k,j,ih)%moments(iq,imom) &
                + dhydrometeors_div(k,j,ih)%moments(iq,imom)
        end do
        ih=qindices(icdnc_bin(iq))%ispecies
        imom=qindices(icdnc_bin(iq))%imoment
        do k=1,nz-1
           sq_lem(j,k+1,icdnc_bin(iq))=dhydrometeors_adv(k,j,ih)%moments(iq,imom) &
                + dhydrometeors_div(k,j,ih)%moments(iq,imom)
        end do
     end do
   end do

   DO k = 2, nz
      DO j = jminp,jmaxp
         DO IQ = 1, LK
            CALL ADVECTcheck(j,k,iq,DT,Q_lem(j,k,ICDKG_BIN(iq)),              &
  &                 Q_lem(j,k,ICDNC_BIN(iq)),SQ_lem(j,k,ICDKG_BIN(iq)),  &
  &                 SQ_lem(j,k,ICDNC_BIN(iq)))
         ENDDO
      ENDDO
   ENDDO

end if

call tau_bin(1, th_lem, q_lem, sth_lem, sq_lem, dt, rdt)

!print*, 'q_lem mass1',q_lem(1,30,icdkg_bin(1))
!print*, 'q_lem num1',q_lem(1,30,icdnc_bin(1))
!stop

! output to ffcd after mphys, might not be right -ahu
do j=jminp,jmaxp
    do k=1,nz
        qv(k,j) = q_lem(j,k,iqv)
        ss(k,j) = q_lem(j,k,iqss)
        thpert(k,j) = th_lem(j,k)
        do iq=1,lk
            ffcd_mass2d(k,j,iq) = (q_lem(j,k,icdkg_bin(iq)) + sq_lem(j,k,icdkg_bin(iq))*dt)/col
            ffcd_num2d(k,j,iq) = (q_lem(j,k,icdnc_bin(iq)) + sq_lem(j,k,icdnc_bin(iq))*dt)/col
        enddo
    enddo
enddo

!print*, sq_lem(1,61,1)

do k=1,nz
    do j=jminp,jmaxp
        do iq=1,nqp
            if (q_lem(j,k,iq) .ne. q_lem(j,k,iq)) then
                print*,'q', k,iq
            end if
            if (sq_lem(j,k,iq) .ne. sq_lem(j,k,iq)) then
                print*,'sq', k,iq
            end if
        enddo
    enddo
enddo

if (any(q_lem .ne. q_lem) .or. any(sq_lem .ne. sq_lem)) then
    print*, 'q', q_lem(1,58,:)
    print*, 'sq', sq_lem(1,58,:)
!    print*, ffcd_mass2d(40,1,:)
!    print*, ffcd_num2d(40,1,:)
    print*, '**at least the mphys routine finished**'
stop

end if

!print*, q_lem(1,30,15), sq_lem(1,30,15)
!print*, 'mass',ffcd_mass2d(30,1,1)
!print*, 'num', ffcd_num2d(30,1,1)


do j=jminp,jmaxp
    if (l_advect .or. l_diverge) then
        sth_lem(j,2:kkp)=sth_lem(j,2:kkp)-(dtheta_adv(1:kkp-1,j)+dtheta_div(1:kkp-1,j))
        sq_lem(j,2:kkp,iqv)=sq_lem(j,2:kkp,iqv)-(dqv_adv(1:kkp-1,j)+dqv_div(1:kkp-1,j))

        do iq=1,ln2
            ih=qindices(iaero_bin(iq))%ispecies
            imom=qindices(iaero_bin(iq))%imoment
            do k=1,nz-1
                sq_lem(j,k+1,iaero_bin(iq))=sq_lem(j,k+1,iaero_bin(iq))       &
                    - (daerosol_adv(k,j,ih)%moments(iq,imom)                &
                    + daerosol_div(k,j,ih)%moments(iq,imom))
            end do
        enddo

        do iq=1,lk
            ih=qindices(icdkg_bin(iq))%ispecies
            imom=qindices(icdkg_bin(iq))%imoment
            do k=1,nz-1
                sq_lem(j,k+1,icdkg_bin(iq))= sq_lem(j,k+1,icdkg_bin(iq))      &
                    - (dhydrometeors_adv(k,j,ih)%moments(iq,imom)           &
                    + dhydrometeors_div(k,j,ih)%moments(iq,imom))
            end do
            ih=qindices(icdnc_bin(iq))%ispecies
            imom=qindices(icdnc_bin(iq))%imoment
            do k=1,nz-1
                sq_lem(j,k+1,icdnc_bin(iq))= sq_lem(j,k+1,icdnc_bin(iq))      &
                    - (dhydrometeors_adv(k,j,ih)%moments(iq,imom)           &
                    + dhydrometeors_div(k,j,ih)%moments(iq,imom))
            end do
        end do
    end if
!   For now set no microphysics on the bottom level - this would be
!   better done by having a subterranian level 0 in column variables
    sq_lem(j,1,1)=0

    dtheta_mphys(1:kkp-1,j)=sth_lem(j,2:kkp)

    dqv_mphys(1:kkp-1,j)=sq_lem(j,2:kkp,iqv)

!
!   update supersaturation field here (not in step fields)
!
    ss(1:kkp-1,j) = q_lem(j,2:kkp,iqss)

    do iq=1,ln2
        ih=qindices(iaero_bin(iq))%ispecies
        imom=qindices(iaero_bin(iq))%imoment
        do k=1,nz-1
           daerosol_mphys(k,j,ih)%moments(iq,imom) =              &
               sq_lem(j,k+1,iaero_bin(iq))
        end do
    enddo

    do iq=1,lk
        ih=qindices(icdkg_bin(iq))%ispecies
        imom=qindices(icdkg_bin(iq))%imoment
        do k=1,nz-1
            dhydrometeors_mphys(k,j,ih)%moments(iq,imom) =         &
                   sq_lem(j,k+1,icdkg_bin(iq))
        end do
        ih=qindices(icdnc_bin(iq))%ispecies
        imom=qindices(icdnc_bin(iq))%imoment
        do k=1,nz-1
            dhydrometeors_mphys(k,j,ih)%moments(iq,imom) =        &
                sq_lem(j,k+1,icdnc_bin(iq))
        end do
    end do
end do

end subroutine micro_proc_tau

!-------------------------------------------------------------------
subroutine micro_init_tau
use micro_prm
use module_bin_init
use mphys_tau_bin_declare
use module_bin_init
implicit none

integer :: j,k
rprefrcp(2:kkp)=exner(1:kkp-1,nx) ! I think this is upside-down in LEM
! AH - 04/03/10, line below leads to divide by 0
!      corrected by setting array to 2:kkp. Problem highlighted
!      by Theotonio Pauliquevis
! prefrcp(:)=1./rprefrcp(:)
prefrcp(2:kkp)=1./rprefrcp(2:kkp)

prefn(2:kkp)=p0*exner(1:kkp-1,nx)**(1./this_r_on_cp)

dzn(2:kkp)=dz_half(1:kkp-1)

rhon(2:kkp)=rho(1:kkp-1)
rdz_on_rhon(2:kkp)=1./(dz(1:kkp-1)*rhon(2:kkp))
! Reference temperature (this is fixed in lem, but
! shouldn't make a difference for microphysics if we
! just set it to be the current profile (i.e. th'=0)
tref(2:kkp)=theta(1:kkp-1,nx)*exner(1:kkp-1,nx)


! Set up microphysics species
call set_micro

do j=jminp,jmaxp
   q_lem (j, 2:kkp, iqv) = qv(1:kkp-1, j)
   q_lem (j, 2:kkp, iqss) = ss(1:kkp-1, j)

   do iq=1,ln2
     ih=qindices(IAERO_BIN(iq))%ispecies
     imom=qindices(IAERO_BIN(iq))%imoment
     do k=1,nz-1
         q_lem (j, k+1, IAERO_BIN(iq)) = aerosol(k,j,ih)%moments(iq,imom)
     end do
   enddo

   do iq=1,lk
     ! mass bins
     ih=qindices(ICDKG_BIN(iq))%ispecies
     imom=qindices(ICDKG_BIN(iq))%imoment
     do k=1,nz-1
        q_lem (j, k+1, ICDKG_BIN(iq)) = hydrometeors(k,j,ih)%moments(iq,imom)
     end do
     ! number bins
     ih=qindices(ICDNC_BIN(iq))%ispecies
     imom=qindices(ICDNC_BIN(iq))%imoment
     do k=1,nz-1
        q_lem (j, k+1, ICDNC_BIN(iq)) = hydrometeors(k,j,ih)%moments(iq,imom)
     end do
   enddo
   th_lem (j, :) = 0.0
   w_lem(j,2:kkp)=w_half(1:kkp-1,j)
end do

call bin_init !initialises the cloud bin categories
call data     !reads in and sets the coll-coal kernal

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

end subroutine micro_init_tau

! sub-subroutine for micro_init_tau ------------
subroutine set_micro
Use parameters, only : num_h_moments, num_h_bins, &
     nspecies, mom_names, h_names, mom_units, max_char_len, &
     num_aero_moments,num_aero_bins, aero_mom_init
Use mphys_tau_bin_declare
Use micro_prm
implicit none

integer :: i

rdt=1./dt

! vapour

iq=1
iqv=iq

if (num_h_bins(1) >= 11) then
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

! sub-subroutine for micro_init_tau ------------
subroutine qcount(var, count)
integer, intent(out) :: var
integer, intent(inout) ::  count

count=count+1
var=count
end subroutine qcount

!-------------------------------------------------------------------
Subroutine init_dist_tau(rxc,gnuc,dnc,rxr,gnur,dnr,diams,ffcd_mass,ffcd_num)

use micro_prm, only:nkr
use parameters, only: max_nbins
use mphys_tau_bin_declare, only: DIAM, NQP
use physconst, only: pi, rhoW

implicit none

integer :: kr
real:: rxc,gnuc,dnc,rxr,gnur,dnr
real(8):: n0c,exptermc,n0r,exptermr
real(8), dimension(max_nbins) :: diams, ffcd_mass, ffcd_num

diams = DIAM(1:max_nbins)*.01 !convert to metric and ditch the last dummy element (?) -ahu
!print*, 'inside init dist',rx,gnu,dn
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
    ffcd_num(kr) = ffcd_mass(kr)/(diams(kr)**3*pi/6.*rhoW)
  endif
  if (rxr>0.) then
     exptermr=exp(-1.*diams(kr)/dnr)
     ffcd_mass(kr) = ffcd_mass(kr) + n0r*exptermr*(diams(kr)/dnr)**(gnur+3)
     ffcd_num(kr) = ffcd_num(kr) + ffcd_mass(kr)/(diams(kr)**3*pi/6.*rhoW)
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
