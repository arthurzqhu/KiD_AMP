module module_hujisbm

use micro_prm

! YWLL_1000MB(nkr,nkr) - input array of kernels for pressure 1000mb
! YWLL_750MB(nkr,nkr) - input array of kernels for pressure 750mb
! YWLL_500MB(nkr,nkr) - input array of kernels for pressure 500mb
REAL, SAVE :: &
! CRYSTALS 
   YWLI(NKR,NKR,ICEMAX) &
! MIXTURES
   ,YWIL(NKR,NKR,ICEMAX),YWII(NKR,NKR,ICEMAX,ICEMAX) &
   ,YWIS(NKR,NKR,ICEMAX),YWIG(NKR,NKR,ICEMAX) &
   ,YWIH(NKR,NKR,ICEMAX),YWSI(NKR,NKR,ICEMAX) &
   ,YWGI(NKR,NKR,ICEMAX),YWHI(NKR,NKR,ICEMAX)
!
REAL, DIMENSION(NKR,NKR),SAVE :: &
    YWLL_1000MB,YWLL_750MB,YWLL_500MB,YWLL,YWLS,YWLG,YWLH &
! SNOW :
    ,YWSL,YWSS,YWSG,YWSH &
! GRAUPELS :
    ,YWGL,YWGS,YWGG,YWGH &
! HAIL :
    ,YWHL,YWHS,YWHG,YWHH
REAL,SAVE :: &
       XI(NKR,ICEMAX) &
      ,RADXX(NKR,NHYDR-1),MASSXX(NKR,NHYDR-1),DENXX(NKR,NHYDR-1) &
      ,RADXXO(NKR,NHYDRO),MASSXXO(NKR,NHYDRO),DENXXO(NKR,NHYDRO) &
      ,RIEC(NKR,ICEMAX),COEFIN(NKR),SLIC(NKR,6),TLIC(NKR,2) &
      ,RO2BL(NKR,ICEMAX)
REAL, SAVE :: VR1(NKR),VR2(NKR,ICEMAX),VR3(NKR) &
      ,VR4(NKR),VR5(NKR),VRX(NKR)
REAL,DIMENSION(NKR),SAVE ::  &
       XL,RLEC,XX,XCCN,XS,RSEC &
      ,XG,RGEC,XH,RHEC,RO1BL,RO3BL,RO4BL,RO5BL &
      ,ROCCN,RCCN,DROPRADII
REAL, SAVE ::  ICEN(NKR)
REAL, SAVE ::  FCCNR0(NKR), FCCNR_mp(NKR)

REAL :: C2,C3,C4
double precision,save ::  cwll(nkr,nkr)
double precision,save::  &
      xl_mg(0:nkr),xs_mg(0:nkr),xg_mg(0:nkr),xh_mg(0:nkr) &
     ,xi1_mg(0:nkr),xi2_mg(0:nkr),xi3_mg(0:nkr) &
     ,chucm(nkr,nkr),ima(nkr,nkr) &
     ,cwll_1000mb(nkr,nkr),cwll_750mb(nkr,nkr),cwll_500mb(nkr,nkr) &
     ,cwli_1(nkr,nkr),cwli_2(nkr,nkr),cwli_3(nkr,nkr) &
     ,cwls(nkr,nkr),cwlg(nkr,nkr),cwlh(nkr,nkr) &

     ,cwil_1(nkr,nkr),cwil_2(nkr,nkr),cwil_3(nkr,nkr) &

     ,cwii_1_1(nkr,nkr),cwii_1_2(nkr,nkr),cwii_1_3(nkr,nkr) &
     ,cwii_2_1(nkr,nkr),cwii_2_2(nkr,nkr),cwii_2_3(nkr,nkr) &
     ,cwii_3_1(nkr,nkr),cwii_3_2(nkr,nkr),cwii_3_3(nkr,nkr) &

     ,cwis_1(nkr,nkr),cwis_2(nkr,nkr),cwis_3(nkr,nkr) &
     ,cwig_1(nkr,nkr),cwig_2(nkr,nkr),cwig_3(nkr,nkr) &
     ,cwih_1(nkr,nkr),cwih_2(nkr,nkr),cwih_3(nkr,nkr) &

     ,cwsl(nkr,nkr) &
     ,cwsi_1(nkr,nkr),cwsi_2(nkr,nkr),cwsi_3(nkr,nkr)&
     ,cwss(nkr,nkr),cwsg(nkr,nkr),cwsh(nkr,nkr) &
     ,cwgl(nkr,nkr)&
     ,cwgi_1(nkr,nkr),cwgi_2(nkr,nkr),cwgi_3(nkr,nkr)&
     ,cwgs(nkr,nkr),cwgg(nkr,nkr),cwgh(nkr,nkr) &

     ,cwhl(nkr,nkr) &
     ,cwhi_1(nkr,nkr),cwhi_2(nkr,nkr),cwhi_3(nkr,nkr) &
     ,cwhs(nkr,nkr),cwhg(nkr,nkr),cwhh(nkr,nkr) &
     ,dlnr &
     ,CTURBLL(KRMAX_LL,KRMAX_LL)&
     ,CTURB_LL(K0_LL,K0_LL)&
     ,CTURBGL(KRMAXG_GL,KRMAXL_GL)&
     ,CTURB_GL(K0G_GL,K0L_GL)

DOUBLE PRECISION, save :: &
     BRKWEIGHT(JBREAK),PKIJ(JBREAK,JBREAK,JBREAK), &
     QKJ(JBREAK,JBREAK),ECOALMASSM(NKR,NKR)

contains
 ! +-------------------------------------------------------------+
   SUBROUTINE FALFLUXHUCM_Z(chem_new,VR1,RHOCGS,PCGS,ZCGS,DT, &
   						               kts,kte,nkr)

     IMPLICIT NONE

 	   integer,intent(in) :: kts,kte,nkr
 	   real(kind=r4size),intent(inout) :: chem_new(:,:)
 	   real(kind=r4size),intent(in) :: rhocgs(:),pcgs(:),zcgs(:),VR1(:,:),DT

 	  ! ... Locals
 	  integer :: I,J,K,KR
    real(kind=r4size) :: TFALL,DTFALL,VFALL(KTE),DWFLUX(KTE)
    integer :: IFALL,N,NSUB

 ! FALLING FLUXES FOR EACH KIND OF CLOUD PARTICLES: C.G.S. UNIT
 ! ADAPTED FROM GSFC CODE FOR HUCM
 !  The flux at k=1 is assumed to be the ground so FLUX(1) is the
 ! flux into the ground. DWFLUX(1) is at the lowest half level where
 ! Q(1) etc are defined. The formula for FLUX(1) uses Q(1) etc which
 ! is actually half a grid level above it. This is what is meant by
 ! an upstream method. Upstream in this case is above because the
 ! velocity is downwards.
 ! USE UPSTREAM METHOD (VFALL IS POSITIVE)
       
       DO KR=1,NKR
        IFALL=0
        DO k = kts,kte
           IF(chem_new(K,KR).GE.1.E-20)IFALL=1
        END DO
        IF (IFALL.EQ.1)THEN
         TFALL=1.E10
         DO K=kts,kte
          ! [KS] VFALL(K) = VR1(K,KR)*SQRT(1.E6/PCGS(K))
           VFALL(K) = VR1(K,KR) ! ... [KS] : The pressure effect is taken into account at the beggining of the calculations
           TFALL=AMIN1(TFALL,ZCGS(K)/(VFALL(K)+1.E-20))
         END DO
         IF(TFALL.GE.1.E10)STOP
	 NSUB=(INT(2.0*DT/TFALL)+1)
         DTFALL=DT/NSUB

         DO N=1,NSUB
           DO K=KTS,KTE-1
             DWFLUX(K)=-(RHOCGS(K)*VFALL(K)*chem_new(k,kr)- &
             RHOCGS(K+1)* &
             VFALL(K+1)*chem_new(K+1,KR))/(RHOCGS(K)*(ZCGS(K+1)- &
             ZCGS(K)))
           END DO
 ! NO Z ABOVE TOP, SO USE THE SAME DELTAZ
           DWFLUX(KTE)=-(RHOCGS(KTE)*VFALL(KTE)* &
      &                 chem_new(kte,kr))/(RHOCGS(KTE)*(ZCGS(KTE)-ZCGS(KTE-1)))
           DO K=kts,kte
            chem_new(k,kr)=chem_new(k,kr)+DWFLUX(K)*DTFALL
           END DO
         END DO
        END IF
       END DO

       RETURN
       END SUBROUTINE FALFLUXHUCM_Z
 ! +----------------------------------+
!--------------------------------------------------
FUNCTION sum_mass(speciesbin,xbin,dens,krs,kre)

use micro_prm, only:nkr,col

implicit none

integer :: krs,kre
real, dimension(nkr) :: speciesbin,xbin
real :: sum_mass,dens

sum_mass = col3/dens*sum(speciesbin(krs:kre)*xbin(krs:kre)*xbin(krs:kre))

END FUNCTION
!-------------------------------------------------
SUBROUTINE JERNUCL01(PSI1,PSI2,FCCNR,DTT,ROR,DSUP1,DSUP2 &
 ,SUP2_OLD,DT,finr &
 ,nuccldrt,nuccldct,nucicert,nucicect,inucifnrt,inucifnct)   

IMPLICIT NONE 

REAL SUP2_OLD, &
     FCCNR(NKR),FINR(NKR)
DOUBLE PRECISION DTT,DSUP1,DSUP2
REAL TT,ROR, &
     TPC,SUP1,SUP2,DEL1N,DEL2N,AL1,AL2, &
     TEMP1,TEMP2,TEMP3,A1,B1,A2,B2,DT
real nucicert,nucicect,inucifnrt,inucifnct,nuccldrt,nuccldct
REAL PSI1(NKR),PSI2(NKR,ICEMAX)
REAL alwc
DATA A1,B1,A2,B2/-0.639,0.1296,-2.8,0.262/
DATA TEMP1,TEMP2,TEMP3/-5.,-2.,-20./
DATA AL1/2500./,AL2/2834./
SUP1=DSUP1
SUP2=DSUP2

TT=DTT

DEL1N=100.*SUP1
TPC=TT-273.15

IF(DEL1N.GT.0.AND.TPC.GT.-30.) THEN
   CALL WATER_NUCL (PSI1,FCCNR,XL,TT,SUP1  &
       ,RCCN,DROPRADII,NKR &
       ,nuccldrt,nuccldct,ror)
ENDIF

IF (ICEPROCS.EQ.1)THEN
   DEL2N=100.*SUP2
   IF(TPC.LT.0..AND.TPC.GE.-35..AND.DEL2N.GT.0.) THEN
         CALL ICE_NUCL (PSI2,XI,TT,ROR,SUP2,SUP2_OLD,SUP1,DT,finr &
                        ,nucicect,nucicert,inucifnct,inucifnrt)
   ENDIF
ENDIF

RETURN
END SUBROUTINE JERNUCL01
!======================================================================   
SUBROUTINE WATER_NUCL (PSI1,FCCNR,X1,TT,SUP1 &
,RCCN,DROPRADII,NKR,nuccldrt,nuccldct,dens)

USE MICRO_PRM, ONLY:IMBUDGET
IMPLICIT NONE
INTEGER NDROPMAX,KR,NKR
REAL PSI1(NKR),FCCNR(NKR),X1(NKR)
REAL DROPCONCN(NKR)
REAL RCCN(NKR),DROPRADII(NKR)
REAL TT,SUP1,DX
real nuccldrt,nuccldct,dens

CALL NUCLEATION (SUP1,TT,FCCNR,DROPCONCN  &
     ,NDROPMAX,RCCN,DROPRADII,NKR)

! NEW WATER SIZE DISTRIBUTION FUNCTION (BEGIN)
  DO KR=1,NDROPMAX
     DX=3.*COL*X1(KR)
     PSI1(KR)=PSI1(KR)+DROPCONCN(KR)/DX

     IF (imbudget >= 1) THEN
        nuccldrt = nuccldrt + dropconcn(kr)*x1(kr)/dens
        nuccldct = nuccldct + dropconcn(kr)/dens*1000.
     ENDIF
  ENDDO

RETURN
END SUBROUTINE WATER_NUCL
!======================================================================
!     ICE_NUCL Modifcation History
!     (April 2007 J. Comstock)
!     modified to include classical theory heterogeneous nucleation
!     via the condensation/immersion freezing mode
!     Added ICEFLAG=0 Use Meyers param, ICEFLAG=1 Use Classical Theory

SUBROUTINE ICE_NUCL (PSI2,X2,TT,ROR,SUP2,SUP2_OLD &
                     ,SUP1,DT,finr &
                     ,nucicect,nucicert,inucifnct,inucifnrt)
  IMPLICIT NONE
  INTEGER ITYPE,KR,ICE,NRGI,K1,ki
  REAL DEL2N,SUP1,SUP2,C1,C2,TPC,TT,ROR
  REAL DX
  REAL HELEK1,HELEK2,TPCC,DEL2NN,FF1BN
  REAL FACT,DSUP2N,SUP2_OLD,DELTACD,DELTAF,ADDF
  REAL X2(NKR,ICEMAX),PSI2(NKR,ICEMAX)
  REAL DT,finr(NKR)
  REAL QHET(NKR),NUMHET(NKR)

  real rcn(nkr)
  real inic(nkr)
  real nucicert,nucicect,inucifnrt,inucifnct

  REAL A1,B1,A2,B2
  DATA A1,B1,A2,B2/-0.639,0.1296,-2.8,0.262/
  REAL TEMP1,TEMP2,TEMP3
  DATA TEMP1,TEMP2,TEMP3/-5.,-2.,-20./

  C1=C1_MEY
  C2=C2_MEY
! TYPE OF ICE WITH NUCLEATION (BEGIN)

  TPC=TT-273.15
  ITYPE=0

  IF((TPC.GT.-4.0).OR.(TPC.LE.-8.1.AND.TPC.GT.-12.7).OR.&
     (TPC.LE.-17.8.AND.TPC.GT.-22.4)) THEN
    ITYPE=2
  ELSE
    IF((TPC.LE.-4.0.AND.TPC.GT.-8.1).OR.(TPC.LE.-22.4)) THEN
      ITYPE=1
    ELSE
      ITYPE=3
    ENDIF
  ENDIF

! NEW CRYSTAL SIZE DISTRIBUTION FUNCTION                      (BEGIN)

  ICE=ITYPE

  IF (ICEFLAG .EQ. 0) THEN  !USE MEYERS ICE NUCLEATION SCHEME
     NRGI=2
     IF(TPC.LT.TEMP1) THEN
        DEL2N=100.*SUP2
        DEL2NN=DEL2N
        IF(DEL2N.GT.50.0) DEL2NN=50.
           HELEK1=C1*EXP(A1+B1*DEL2NN)
      ELSE
           HELEK1=0.
      ENDIF

      IF(TPC.LT.TEMP2) THEN
         TPCC=TPC
         IF(TPCC.LT.TEMP3) TPCC=TEMP3
            HELEK2=C2*EXP(A2-B2*TPCC)
      ELSE
         HELEK2=0.
      ENDIF

      FF1BN=HELEK1+HELEK2

      FACT=1.
      DSUP2N=(SUP2-SUP2_OLD)*100.

      SUP2_OLD=SUP2

      IF(DSUP2N.GT.50.) DSUP2N=50.

      DELTACD=FF1BN*B1*DSUP2N

      IF(DELTACD.GE.FF1BN) DELTACD=FF1BN
      IF(DELTACD.GT.0.) THEN
         DELTAF=DELTACD*FACT
         DO KR=1,NRGI-1
            DX=X2(KR,ICE)*COL3
            ADDF=DELTAF/DX
            PSI2(KR,ICE)=PSI2(KR,ICE)+ADDF
            IF (imbudget >=1) THEN
               nucicert = nucicert + DELTAF/ROR
               nucicect = nucicect + DELTAF/X2(KR,ICE)/ROR*1000.
            ENDIF
            IF (imbudget >= 2) THEN
               inucifnrt = inucifnrt + DELTAF/ROR
               inucifnct = inucifnct + DELTAF/X2(KR,ICE)/ROR*1000.
            ENDIF
         ENDDO
      ENDIF
   ELSE   !USE CLASSICAL THEORY NUCLEATION SCHEME

!===================CCN-version ice nucleation (added by J. Fan)=======================

      DO KR=1,NKR
         rcn(kr) = rccn(kr)*1.e-2              ! transform to m from cm
      ENDDO
          
      DO KR=1,NKR
         inic(kr) = finr(kr)*col*1.e6  ! m-3
      ENDDO

       CALL heteronuc_ccn(TT,SUP1,DT,NKR,inic,rcn,QHET,NUMHET)
       ! Distribute nucleated aerosols to appropriate ice size bin (QHET:kg; NUMHET: m-3)

       ! Update IN
       DO KR=1, NKR
          finr(kr) = inic(kr)*1.e-6/col
       ENDDO

       if (sum(NUMHET(:)) > 0.) then
          do ki=1,NKR
             do kr=1,NKR
                if (QHET (ki)*1.e3 .LE. X2(1,ice) )  k1=1
                if (QHET (ki)*1.e3 .GE. X2(nkr,ice) )  k1=nkr
                if ((kr.GT.1) .and. (QHET (ki)*1.e3 .LE. X2(kr,ice)).and.  &
                   (QHET (ki)*1.e3 .GT. X2(kr-1,ice)))  k1=kr
             enddo
             DX=3.*X2(k1,ICE)*COL
             PSI2(K1,ICE) = PSI2(K1,ICE) + NUMHET(ki)*1.e-6/DX
             IF (imbudget >=1) THEN
                nucicert = nucicert + numhet(ki)/DX*COL3*X2(k1,ICE)*X2(k1,ICE)/ROR*1.e-6
                nucicect = nucicect + numhet(ki)/DX*COL3*X2(k1,ICE)/ROR*1.e-3
             ENDIF
             IF (imbudget >=2) THEN
                inucifnrt = inucifnrt + numhet(ki)/DX*COL3*X2(k1,ICE)*X2(k1,ICE)/ROR*1.e-6
                inucifnct = inucifnct + numhet(ki)/DX*COL3*X2(k1,ICE)/ROR*1.e-3
             ENDIF
              
          enddo
       endif
   ENDIF

RETURN
END SUBROUTINE ICE_NUCL

!=====================================================================
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Computes heterogeneous freezing rate, crystal number and mass assuming
! immersion freezing from ccn distribution (assuming ccn in aqueous soln)

! Assumes nucleation is immersion/condensation freezing from solution drops+ccn

subroutine heteronuc_ccn(temp,Sw,dt,INBINMAX,icenuclei, &
           draeros,qhet,numhet)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     This is main subroutine for heterogeneous nucleation
!     described in Morrison et al. 2005 JAS.

!     modified to interface with HUCM Spectral Bin Model (18 Apr 2007 JMC)

!     Germ formation rates due to condensation freezing are described in
!     Khvorostyanov and Curry 2005 JAS.
!     This formulation passes the aerosol (IN) size distribution into subroutine
!     and tracks IN loss due to nucleation scavenging.

!     modified to pass the drop size distribution rather than the IN dry aerosol
!     distribution and removed aerosol growth. (Oct 2007 JMC)

implicit none

!     input:
integer INBINMAX    ! number of aerosol (IN) bins
real temp,        &  ! temperature (K)
    Sw  ,        &  ! water supersaturation ratio (dimensionless)
    dt  ,        &  ! time step (s)
    icenuclei(INBINMAX),& ! ice nuclei (from ccn distribution) (m-3)
    draeros(INBINMAX)    ! dry aerosol radius (m)

!     output:
!     note: numhet and qhet are the total number aerosol nucleated during time step,
!     not nucleation rate!!!
real numhet(INBINMAX)     ! number of aerosol nucleated for aerosol bin (m-3)
real qhet(INBINMAX)       ! mass of single ice crystal nucleated per aerosol bin (kg)

!     internal:
real raeros,          &    ! radius of insoluble substrate (m)
    wraeros,         &    ! wet aerosol radius (m)
    ajlsfrz,         &    ! het frz rate for particle
    probhet,         &    ! probability of het frz for particle
    Tc,              &    ! temperature (C)
    hrel                 ! saturation ratio (unitless)
integer kr                ! counter for aerosol size bin

!     other parameter
real rhow,             &   !density of water
     rhoi                  !density of ice
real pi

Tc = temp-273.15          !! K to C
hrel = Sw + 1.0           ! saturation ratio

!     define constants
pi = 4.*atan(1.)

!      qvaero = 0.85   !aerosol solubility in terms of volume fraction
!                      this quantity is assumed, may want to upgrade formulation
!     wettability
!      wetcoef = 0.9
!     relative area of active sites
!      alf = 0.5e-5
!     misfit strain parameter
!      epsil = 2.5e-2            ! e=2.5% Turnbull vonegut, 1952
!      epsil = 0.1e-2            !! e~0; Turnbull vonegut, 1952
!      epsil=0.01  !e=1%
!     density of water (kg/m3)
rhow = 1000.0
!     density of ice (kg/m3)
rhoi = 900.0
!     aerosol parameter describing composition (Khvorostyanov and Curry 1999)
!      betafr = 0.5

!-----------------------------------------------------------------------
!     numerical integration of size bins

!     initialize nucleation number and mass

do kr = 1,INBINMAX
   numhet(kr) = 0.
   qhet(kr) = 0.
end do

!     main loop around aerosol size bins, starting with largest bin

do kr=INBINMAX,1,-1

!     determine radius of insoluble portion of aerosol raeros
!     size of insoluble substrate determines germ formation rate
!     See Khvorostyanov and Curry 2004

   raeros = draeros(kr) * (1.-qvaero)**(1./3.)  !m

!   print*,'raeros=',raeros
!     call to get germ formation rate

   if (raeros.gt.0.) then
      call jhetfrz2(Tc,hrel,raeros,ajlsfrz)
   else
      ajlsfrz = 0.
   end if

!     calculate probability of aerosol particle freezing within dt

   probhet = 1.-exp(-ajlsfrz*dt)
!   if (probhet > 1.e-3) print *,'probhet',probhet,ajlsfrz,raeros

!     print *,probhet,ajlsfrz,raeros

!     if there is no freezing associated with largest size bin, then
!     exit loop, there will be no nucleation
!     use probability of 1.e-10 as cutoff

   if (kr.eq.INBINMAX) then
      if (probhet.lt.1.e-10) return
   end if

!     nucleate if probability is significant

  if (probhet.ge.1.e-10) then

!     number of ice nucleated for each drop size bin (1/m3)

      numhet(kr) = probhet * AMAX1(icenuclei(kr),0.0)
!    print*,'numhet',kr,numhet(kr),ajlsfrz,probhet,icenuclei(kr)

!----------------------------------------------------------------------------------
!     note, use wetted aerosol for radius to calculate mass
!     of nucleated crystals for size bin (kg m-3)
!     this is done by multiplying nucleation rate (m-3) by
!     4/3*pi*rhoi*r^3
!     rhoi = 900 kg m-3 is assumed bulk density of ice

      call req_haze(temp,Sw,draeros(kr),wraeros)
!      print*, 'wet rad', kr, draeros(kr), wraeros

      qhet(kr) = 4./3.*pi*rhoi*(wraeros)**3.0 !mass of a single particle in (kg)

      icenuclei(kr) = icenuclei(kr) - numhet(kr)

   end if                 !! significant freezing probability
end do                    !! end loop around size bins
return
end subroutine heteronuc_ccn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine jhetfrz2(ttc,relh,raeros,ajlshet)

!     this subroutine calculates germ formation rates on insoluble subsrate for
!     given T, S, etc.
!     Khvorostyanov and Curry, 2005 JAS
!     All units are mks (JMC) 18 Apr 2007

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

implicit none

!     input data

real ttc,     &      ! temperature (C!!!!)
     relh,    &      ! saturation ratio over water
     sig_iw0, &      !surface tension ice/water interface
     raeros          !radius of the insolouble substrate (m)

!     output data

real  ajlshet       !nucleation rate (s-1)
real  sw_th    !threshold saturation

!     internal variables
!!!!!!! Modified by Yi Wang !!!!!!!
!!!! Using c1s (c,one,s) is same as cls (c,l,s)
!!!! Change c1s to cleans

REAL cturn,rhoimean,akbol,rhoi,akplank,cleans,almmeanc  &
          ,a_lm0,a_lm1,a_lm2,ttze,a_lmtn(4),rlog,cor1,cor2   &
          ,almmean,ttk,tt0,humf,expind,sw_thresh,fgerms,Rv
REAL factivh,xxx,fisize,fisizeSS,xmfi,ccc1,ccc2,ccc3, &
     ffshape2,fshape,fgerms0,dfgerms,coefnucl,rgerms

!     Turnb-Vonneg strain const, [NT/m^-2], 24 June 2000, PK97, p.343
CTurn=1.7e10
cleans=1.e19     ! [m^-2], conc. molec. per 1 m^2 of surface, PK97, p.342
              ! PK97, p. 206, was 10^18!!!  !m^-2, number of molec. per 1
              ! [m^2], contact. ice germ
RHOIMEAN=900.0      !kg/m3, mean ice density
RHOI=900.0          !kg/m3, ice density


akbol    =1.380622e-23    ![J/K], Boltzmann constant
akplank  =6.626176e-34    ![J*s], Planck constant
Rv       =4.615e2         !specific gas constant for water vapor (J/kg/K)

sig_iw0=(28.0+0.25*TTC)*1.e-3 !surface tension, [NT/m] ice-water interface


TTZE=273.15

A_LM0    =70.0            ! Lm=70 cal/g !Prup'95
A_LM1    =0.708
A_LM2    =-2.5e-3         !L''m, Prup.'78, p.89
ALMMEANC =79.7 + 0.693*TTC-2.5E-3*TTC**2 !old PK78, corr. 2 June 1998
A_LMTN(1)=79.7 + A_LM1*TTC+A_LM2*(TTC**2) ! Lm, Prup'78, p.89;
A_LMTN(2)=A_LM0 + A_LM1*TTC+A_LM2*TTC**2 ! Pr95, fig.8, PK97, fig. 3-12

!     below:  Temperature correction  fi(T) from
!     Integr. melt. heat  -temp. corr.

RLOG    =alog(TTZE/(TTC+TTZE))
cor1    =-A_LM1/A_LM0*(TTZE+TTC/RLOG)
cor2    = A_LM2/A_LM0*(TTZE**2.0-TTC/RLOG*(TTC-2.0*TTZE)/2.)

A_LMTN(4)=1.+cor1+cor2
A_LMTN(3)=A_LM0 * A_LMTN(4)

!     effective melting heat

ALMMEANC=A_LMTN(1)        ![kcal/kg] or [cal/g]
ALMMEAN=ALMMEANC*4.18e3   ![J/kg]

TTK=TTC+273.15            !Temperature Kelvin
TT0=273.15                !0 C (Kelvin=273.15)
!     param. G=Rv*T/Lm (divided by corr. Lm(TTC), 1 June 1998
humf=0.11*TTK/ALMMEANC


expind=CTurn*epsil**2.0/(RHOI*ALMMEAN)

!     thershold saturation

Sw_thresh=(TTK/TT0*exp(expind))**(1./humf)
Sw_th=Sw_thresh           !for output

!     if sw less than theshold then exit out, no nucleation
!      print*,' Sw_th=',Sw_th
if (Sw_thresh.GE.relh) then !no nucleation, go out
   AJLSHet=1.e-20
   FGERMS=1.e-14
   rgerms=1.e-08
   go to 98765
endif

!     activation energy

!******Activation energy at T=-30 C is 10 kkal/mole (PK97, Fig. 3-11, p.95);
!***   this is =0.694*10^-12 erg; starting point for this fit at T=-30 C.
!***   It decreases with T and coincides with Jensen94 at T=-90 C

FACTIVH= &           ! Activation energy, erg
          0.694E-12*(1.+0.027*(TTC+30.0)) !Linear fit to Prupp.'95  4 Nov 1997
if (TTC.gt.-30.0) FACTIVH=0.385e-12*  &
          exp(-8.423e-3*TTC+6.384e-4*TTC**2.0+7.891e-6*TTC**3.0) !p.96

if (TTC.le.-40.) FACTIVH=  & !for`low T<-40 C
          0.694E-12*(1.+0.027*(TTC+30.0)*exp(0.010*(TTC+30.0))) !MY CHOICE
FACTIVH = FACTIVH * 1.E-7  !convert from erg to Joules

!     rgerms in m
rgerms=2.0*sig_iw0/(ALMMEAN*RHOIMEAN*(ALOG(273.15/TTK)   &
          +humf*ALOG(relh))-CTurn*epsil**2.0) !corr 24 June 2000 for misfit strain

!     this check is to make sure that radius of germ is greater than 0

if (rgerms.gt.0.) go to 12345 !continue if rgerms>0

!*****************go out with AJLSHet=0, if rgerms<0  **********************

if (rgerms.le.0.) then
   AJLSHet=1.e-20         !no nucleation, go out
   FGERMS=1.e-24
   rgerms=1.e-18
!   print*,'rgerms.le.0.'
   go to 98765      !go out with AJLSHet=0, if rgerms<0
endif

12345 xxx=raeros/rgerms

!     shape factor

fisizeSS=(1.0-2.0*wetcoef*xxx+xxx**2.0)**0.5 !12 Nov 1998, PK97, p. 302
fisize=fisizeSS           !18 July 2000, derevn.
xmfi=(xxx-wetcoef)/fisizeSS

ccc1=((1-wetcoef*xxx)/fisizeSS)**3.0
ccc2=(xxx**3.0)*(2.0-3.0*xmfi+xmfi**3.0)
ccc3=3*wetcoef*(xxx**2)*(xmfi-1)
ffshape2=0.5*(1+ccc1+ccc2+ccc3)
fshape=ffshape2
!print*,'fshape=',fshape
FGERMS0=4./3.*3.1416*sig_iw0*(rgerms**2.0)*fshape !Germ energy with shape
!print*,'fgerms0=',fgerms0,sig_iw0,rgerms

!*****!Fletcher's correction for active site, PK97, p.345, !11 April 1999
DFGERMS=alf*(raeros**2.0)*sig_iw0*(1.-wetcoef) !11 April 1999, PK97, p.345
!print*,'dfgerms=',dfgerms

!     correction for active site, PK97, p.345; 24 June 2000

FGERMS=FGERMS0-DFGERMS
if (FGERMS.le.0.) FGERMS=0. ! no zero values ???

!     output to screen
!     if (FGERMS.le.0.) write (*,*) 'FGERMS LT 0.',' TTK=',TTK,
!     # ' relh=',relh,' FGERMS0=',FGERMS0,' DFGERMS=',DFGERMS,'alf=',alf,
!     # ' raeros=',raeros

!     preexpon. fact. PK97, p.342, heterogen. 12 Nov 1998

coefnucl=(akbol*TTK/akplank)*    & !12 Nov 1998, PK97, p. 342
          (4.0*3.1416*raeros**2.0)*cleans
!print*,'coefnucl=',coefnucl,akbol,ttk,akplank,raeros,cleans
!     calculate freezing rate

AJLSHet=           &     !s^-1, heterogen. salt nucleation rate, PK97, p. 342
          coefnucl*exp(-(FACTIVH+FGERMS)/(akbol*TTK)) !KS.98

if (AJLSHet.lt.1.e-20) AJLSHet=1.e-20

98765 return  !go here if rgerms<0
END subroutine jhetfrz2


!!==============================================================

SUBROUTINE REQ_HAZE(temp,sw,draeros,wraeros)

!     this subroutine calculates germ formation rates on insoluble subsrate for
!     given T, S, etc.
!     Khvorostyanov and Curry, 2005 JAS
!     All units are mks (JMC) 18 Apr 2007

!       purpose
!       input   temp            [K]
!               sw              saturation ratio
!               draeros         [m] dry aerosol radius
!               betafr          aerosol parameter describing composition
!               qvaeros         aerosol solubility in terms of volume fraction (unitless)
!       output  wraeros         [m]
!
!       comment
!               dilute solution approx.
!               works for RHW <=100.0
!       Note
!               output wraeros >= draeros
!       source:
!               Hugh Morrison, Jennifer Comstock
!       reference:
!               Khvorostyanov and Curry (1999) JGR
!               http://mathworld.wolfram.com
!===============================================================

IMPLICIT NONE

!      INCLUDE "aerosol_prop_mks.inc"

real      temp            ! temperature [K}
real      sw              ! saturation ratio (water)
real      draeros         ! dry aerosol radius [m]
real      Tc,T0
real      H
real      sigvl,rhow
real      b, A, BB
real      wr100
real      a0,a1,a2,s1,s2
real      Q,Q3,R,R2,D,D12
real      theta
real      wraeros
real      Mw,Rv
!real      rhos2

!print*,'Req_haze...'

!rhos2 = 1770.0        !density of dry aerosol (kg/m3) (NH4)2SO4 Ammonium Sulate KC1999
!Ms2 = 0.132           !molecular weight of dry aerosol (g/mol) (NH4)2SO4 Ammonium Sulfate KC1999
Mw = 1.8e-2           !molecular weight of water (kg/mol)
Rv = 4.615e2              !specific gas constant for water vapor (J/kg/K)
T0 = 273.15
Tc= temp-T0       ! K to C
H = sw+1.0                ! saturation ratio
rhow=  1000.0     ! density of water (kg/m3)

sigvl= 0.0761 - 1.55e-4 * Tc
b =    2.0 * qvaero*(rhos2/rhow)*(Mw/Ms2) !aerosol parameter describing composition,  Eq (11), KC1999
BB =   2.0 * sigvl/(Rv*temp*rhow) !Kelvin parameter, KC1999
A =    b * (draeros)**(2.0+2.0*betafr)
wr100= (A/BB)**0.5   !wet radius at RH=100%

!print*,'BB kel=',BB
!if (h > 0.0) print*,'check para', qvaero, rhos2, Ms2, draeros, betafr,h

if (H .lt. 1.0) then

!     cubic formula, http://mathworld.wolfram.com

   a0     = A/(H-1.)
   a1     = 0.
   a2     = -BB/(H-1.)
!   print*,a0,a1,a2
!     cubic root solution

   Q= (3.*a1   -(a2*a2))/9.0
   R= (9.*a1*a2-27.*a0-2.*a2*a2*a2)/54.0
   Q3     = Q*Q*Q
   R2     = R*R
   D= Q3+R2
!   print*,'Q&R',Q,R
! print*,'Dvalue', D, Q, R, a0, a2, A, BB

   if (D .gt. 0.0) then
      D12 = D**0.5

      if ( (R + D12).gt.0.0) then
         s1 = (R + D12)**(1./3.)
      else
         s1 = -(abs(R + D12)**(1./3.))
      endif

      if ( (R - D12).gt.0.0) then
         s2 = (R - D12)**(1./3.)
      else
         s2 = -(abs(R - D12))**(1./3.)
      endif
!      print*,'s1&s2=',s1,s2

      wraeros = (s1+s2)-a2/3.
!      print*,'D gt 0:',wraeros
   else if (D .lt. 0.0) then
!     3 real solutions
!     choose the first solution (most realistic)
      theta = acos(R/sqrt(-Q3))

      wraeros =   2 * ((-Q)**0.5) * cos(  theta / 3. ) - a2/3.0
!      print*,'D lt 0:',wraeros
   endif
else
!     for water saturated or greater

   wraeros = wr100
!     make sure aerosol size does not exceed value when RH = 100%

   if (wraeros .gt. wr100) then
      wraeros = wr100
   endif

!     make sure wraeros is not smaller than draeros
   wraeros = AMAX1(wraeros,draeros)
endif

RETURN
END SUBROUTINE REQ_HAZE
!==============================================================================
!     Creates ice nuclei distribution from input ccn distribution
!     currently assume some arbitrary fraction of ccn are viable IN

SUBROUTINE MAKE_IN_DIST(ccn,binmax,icenuc)
IMPLICIT NONE
!     Inputs:
integer binmax     !number of ccn and IN bins
real ccn(binmax)   !ccn distribution
!     Outputs:
real icenuc(binmax)  !IN distribution

integer i

do i=1,binmax
   icenuc(i) = ccn(i)*col*fracin*1.e6     ! m-3
enddo

RETURN

END SUBROUTINE MAKE_IN_DIST

!===================end CCN-version ice nucleation=====================================

SUBROUTINE NUCLEATION (SUP1,TT,FCCNR,DROPCONCN  &
      ,NDROPMAX,RCCN,DROPRADII,NKR)
! DROPCONCN(KR), 1/cm^3 - drop bin concentrations, KR=1,...,NKR

! determination of new size spectra due to drop nucleation

IMPLICIT NONE
INTEGER NDROPMAX,IDROP,INEXT,ISMALL,KR,NCRITI
INTEGER IMIN,IMAX,NKR,I
REAL &
 SUP1,TT,RACTMAX,RCRITI,BKOE, &
 AKOE,DEG01
REAL CCNCONC(NKR)
REAL CCNCONC_BFNUCL,CCNCONC_AFNUCL, DEl_CCNCONC
REAL RCCN(NKR),DROPRADII(NKR),FCCNR(NKR)
REAL RACT(NKR),DROPCONCN(NKR)

  DEG01=1./3.


! calculation initial value of NDROPMAX - maximal number of drop bin
! which is activated

! initial value of NDROPMAX

  NDROPMAX=0

  DO KR=1,NKR
! initialization of bin radii of activated drops
     RACT(KR)=0.
! initialization of aerosol(CCN) bin concentrations
     CCNCONC(KR)=0.
! initialization of drop bin concentrations
     DROPCONCN(KR)=0.
  ENDDO

! CCNCONC_BFNUCL - concentration of aerosol particles before
!nucleation

  CCNCONC_BFNUCL=0.
  DO I=1,NKR
     CCNCONC_BFNUCL=CCNCONC_BFNUCL+FCCNR(I)
  ENDDO

  CCNCONC_BFNUCL=CCNCONC_BFNUCL*COL

  IF(CCNCONC_BFNUCL.EQ.0.) THEN
     RETURN    
  ELSE
     CALL BOUNDARY(IMIN,IMAX,FCCNR,NKR)
     CALL CRITICAL (AKOE,BKOE,TT,RCRITI,SUP1,DEG01)
!     print*, 'rcriti',RCRITI,imax,RCCN(IMAX)
 
     IF(RCRITI.GE.RCCN(IMAX))  RETURN
  END IF

! calculation of CCNCONC(I) - aerosol(CCN) bin concentrations;
!     I=IMIN,...,IMAX
! determination of NCRITI - number bin in which is located RCRITI
  IF (IMIN.EQ.1)THEN
   CALL CCNIMIN(IMIN,IMAX,RCRITI,NCRITI,RCCN,CCNCONC, &
      FCCNR,NKR)
   CALL CCNLOOP(IMIN,IMAX,RCRITI,NCRITI,RCCN,CCNCONC, &
      FCCNR,NKR)
  ELSE
   CALL CCNLOOP(IMIN,IMAX,RCRITI,NCRITI,RCCN,CCNCONC, &
      FCCNR,NKR)
  END IF

! calculation CCNCONC_AFNUCL - ccn concentration after nucleation

 CCNCONC_AFNUCL=0.

 DO I=IMIN,IMAX
    CCNCONC_AFNUCL=CCNCONC_AFNUCL+FCCNR(I)
 ENDDO

 CCNCONC_AFNUCL=CCNCONC_AFNUCL*COL

 DEL_CCNCONC=CCNCONC_BFNUCL-CCNCONC_AFNUCL

 CALL ACTIVATE(IMIN,IMAX,AKOE,BKOE,RCCN,RACTMAX,NKR)
 CALL DROPMAX(DROPRADII,RACTMAX,NDROPMAX,NKR)
! put nucleated droplets into the drop bin according to radius
! change in drop concentration due to activation DROPCONCN(IDROP)
 ISMALL=NCRITI

 INEXT=ISMALL

 DO IDROP=1,NDROPMAX
    DROPCONCN(IDROP)=0.
    DO I=ISMALL,IMAX
       IF(RACT(I).LE.DROPRADII(IDROP)) THEN
          DROPCONCN(IDROP)=DROPCONCN(IDROP)+CCNCONC(I)
          INEXT=I+1
       ENDIF
    ENDDO
    ISMALL=INEXT
 ENDDO

RETURN
END SUBROUTINE NUCLEATION
!########################################################################
SUBROUTINE BOUNDARY(IMIN,IMAX,FCCNR,NKR)
! IMIN - left CCN spectrum boundary
IMPLICIT NONE
INTEGER I,IMIN,IMAX,NKR
REAL FCCNR(NKR)

IMIN=0

DO I=1,NKR
   IF(FCCNR(I).NE.0.) THEN
     IMIN=I
     EXIT
   ENDIF
ENDDO

! IMAX - right CCN spectrum boundary
IMAX=0

DO I=NKR,1,-1
   IF(FCCNR(I).NE.0.) THEN
     IMAX=I
     EXIT
   ENDIF
ENDDO
RETURN
END  SUBROUTINE BOUNDARY
!######################################################################
SUBROUTINE CRITICAL (AKOE,BKOE,TT,RCRITI,SUP1,DEG01)
! AKOE & BKOE - constants in Koehler equation
IMPLICIT NONE
REAL AKOE,BKOE,TT,RCRITI,SUP1,DEG01

  AKOE=3.3E-05/TT
  BKOE=ions*4.3/mwaero
  BKOE=BKOE*(4./3.)*3.141593*RO_SOLUTE                  

! RCRITI, cm - critical radius of "dry" aerosol

  RCRITI=(AKOE/3.)*(4./BKOE/SUP1/SUP1)**DEG01
RETURN
END  SUBROUTINE CRITICAL
!#######################################################################
SUBROUTINE CCNIMIN(IMIN,IMAX,RCRITI,NCRITI,RCCN,CCNCONC, &
      FCCNR,NKR)

IMPLICIT NONE
INTEGER IMIN,II,IMAX,NCRITI,NKR
REAL RCRITI
REAL RCCN(NKR),FCCNR(NKR),CCNCONC(NKR)
REAL RCCN_MIN
REAL DLN1,DLN2
! rccn_min - minimum aerosol(ccn) radius
  RCCN_MIN=RCCN(1)/10000.
! calculation of ccnconc(ii)=fccnr(ii)*col - aerosol(ccn) bin
!                                            concentrations,
!                                            ii=imin,...,imax
! determination of ncriti   - number bin in which is located rcriti
! calculation of ccnconc(ncriti)=fccnr(ncriti)*dln1/(dln1+dln2),
! where,    
! dln1=Ln(rcriti)-Ln(rccn_min)
! dln2=Ln(rccn(1)-Ln(rcriti)

IF(RCRITI.LE.RCCN_MIN) THEN
   NCRITI=1
   DO II=NCRITI+1,IMAX
      CCNCONC(II)=COL*FCCNR(II)     
      FCCNR(II)=0.
   ENDDO
ELSEIF(RCRITI.GT.RCCN_MIN.AND.RCRITI.LT.RCCN(IMIN)) THEN
   NCRITI=1
   DO II=NCRITI+1,IMAX
      CCNCONC(II)=COL*FCCNR(II)
      FCCNR(II)=0.
   ENDDO
   DLN1=ALOG(RCRITI)-ALOG(RCCN_MIN)
   DLN2=ALOG(RCCN(1))-ALOG(RCRITI)
   CCNCONC(NCRITI)=DLN2*FCCNR(NCRITI)
   FCCNR(NCRITI)=FCCNR(NCRITI)*DLN1/(DLN1+DLN2)
ENDIF
RETURN
END SUBROUTINE CCNIMIN
!####################################################################
SUBROUTINE CCNLOOP(IMIN,IMAX,RCRITI,NCRITI,RCCN,CCNCONC, &
                   FCCNR,NKR)
IMPLICIT NONE
INTEGER I,IMIN,IMAX,NKR,II,NCRITI
REAL RCRITI,RCCN(NKR),CCNCONC(NKR),FCCNR(NKR)
REAL DLN1,DLN2

  IF(IMIN.GT.1) THEN
    IF(RCRITI.LE.RCCN(IMIN-1)) THEN
      NCRITI=IMIN
      DO II=NCRITI,IMAX
         CCNCONC(II)=COL*FCCNR(II)
         FCCNR(II)=0.
      ENDDO
      RETURN
    ENDIF
    IF(RCRITI.LT.RCCN(IMIN).AND.RCRITI.GT.RCCN(IMIN-1)) THEN
       NCRITI=IMIN
      
       DO II=NCRITI+1,IMAX
          CCNCONC(II)=COL*FCCNR(II)
          FCCNR(II)=0.
       ENDDO
       DLN1=ALOG(RCRITI)-ALOG(RCCN(IMIN-1))
       DLN2=COL-DLN1
       CCNCONC(NCRITI)=DLN2*FCCNR(NCRITI)
       FCCNR(NCRITI)=FCCNR(NCRITI)*DLN1/COL
       RETURN
    ENDIF
  ENDIF

  DO I=IMIN,IMAX-1
     IF(RCRITI.EQ.RCCN(I)) THEN
       NCRITI=I+1
       DO II=I+1,IMAX
          CCNCONC(II)=COL*FCCNR(II)
          FCCNR(II)=0.
       ENDDO
       RETURN
     ENDIF
     IF(RCRITI.GT.RCCN(I).AND.RCRITI.LT.RCCN(I+1)) THEN
       NCRITI=I+1
       IF(I.NE.IMAX-1) THEN
         DO II=NCRITI+1,IMAX
            CCNCONC(II)=COL*FCCNR(II)
            FCCNR(II)=0.
         ENDDO
       ENDIF
       DLN1=ALOG(RCRITI)-ALOG(RCCN(I))
       DLN2=COL-DLN1
       CCNCONC(NCRITI)=DLN2*FCCNR(NCRITI)
       FCCNR(NCRITI)=FCCNR(NCRITI)*DLN1/COL
       RETURN
     END IF
  ENDDO
RETURN
END  SUBROUTINE CCNLOOP
!##################################################################
SUBROUTINE ACTIVATE(IMIN,IMAX,AKOE,BKOE,RCCN,RACTMAX,NKR)
IMPLICIT NONE

INTEGER IMIN,IMAX,NKR
INTEGER I,I0,I1
REAL RCCN(NKR)
REAL  R03,SUPCRITI,RACT(NKR),XKOE
REAL AKOE,BKOE,AKOE23,RACTMAX

! Spectrum of activated drops                                 (begin) 
  DO I=IMIN,IMAX
  ! critical water supersaturations appropriating CCN radii
     XKOE=(4./27.)*(AKOE**3./BKOE)
     AKOE23=AKOE*2./3.
     R03=RCCN(I)**3.
     SUPCRITI=SQRT(XKOE/R03)

! RACT(I) - radii of activated drops, I=IMIN,...,IMAX

     IF(RCCN(I).LE.(0.3E-5)) RACT(I)=AKOE23/SUPCRITI
     IF(RCCN(I).GT.(0.3E-5)) RACT(I)=5.*RCCN(I)
  ENDDO

  I0=IMIN

   DO I=IMIN,IMAX-1
      IF(RACT(I+1).LT.RACT(I)) THEN
         I0=I+1
         EXIT
      ENDIF
   ENDDO

   I1=I0-1

   IF(I0.EQ.IMIN) GOTO 47

   IF(I0.EQ.IMAX) THEN
      RACT(IMAX)=RACT(IMAX-1)
      GOTO 47
   ENDIF

   IF(RACT(IMAX).LE.RACT(I0-1)) THEN
      DO I=I0,IMAX
         RACT(I)=RACT(I0-1)
      ENDDO
      GOTO 47
   ENDIF

   DO I=I0+1,IMAX
      IF(RACT(I).GE.RACT(I0-1)) THEN
         I1=I
         EXIT
      ENDIF
   ENDDO

! line interpolation RACT(I) for I=I0,...,I1

   DO I=I0,I1
      RACT(I)=RACT(I0-1)+(I-I0+1)*(RACT(I1)-RACT(I0-1)) &
                           /(I1-I0+1)
   ENDDO

 47    CONTINUE

 RACTMAX=0.

 DO I=IMIN,IMAX
    RACTMAX=AMAX1(RACTMAX,RACT(I))
 ENDDO
RETURN

END SUBROUTINE ACTIVATE
!-----------------------------------------------------------
SUBROUTINE DROPMAX(DROPRADII,RACTMAX,NDROPMAX,NKR)
IMPLICIT NONE
INTEGER IDROP,NKR,NDROPMAX
REAL RACTMAX,DROPRADII(NKR)
! calculation of NDROPMAX - maximal number of drop bin which
! is activated

  NDROPMAX=1

  DO IDROP=1,NKR
     IF(RACTMAX.LE.DROPRADII(IDROP)) THEN
       NDROPMAX=IDROP
       RETURN
     ENDIF
  ENDDO
RETURN
END SUBROUTINE DROPMAX

!-------------------------------------------------
   SUBROUTINE ONECOND1 &
 & (TT,QQ,PP,ROR &
 & ,VR1,PSINGLE &
 & ,DEL1N,DEL2N,DIV1,DIV2 &
 & ,FF1,PSI1,R1,RLEC,RO1BL &
 & ,AA1_MY,BB1_MY,AA2_MY,BB2_MY &
 & ,C1_MEY,C2_MEY &
 & ,COL,DTCOND,ICEMAX,NKR,ISYM1 &
   ,ISYM2,ISYM3,ISYM4,ISYM5,Iin,Jin,Kin,W_in,DX_in,Itimestep)

        IMPLICIT NONE


 INTEGER NKR,ICEMAX, ISYM1, ISYM2(ICEMAX),ISYM3,ISYM4,ISYM5, Iin, Jin, Kin, &
 sea_spray_no_temp_change_per_grid, Itimestep
 REAL    COL,VR1(NKR),PSINGLE &
      &       ,AA1_MY,BB1_MY,AA2_MY,BB2_MY &
      &       ,DTCOND, W_in,DX_in

 REAL C1_MEY,C2_MEY
 INTEGER I_ABERGERON,I_BERGERON, &
      & KR,ICE,ITIME,KCOND,NR,NRM, &
      & KLIMIT, &
      & KM,KLIMITL
 REAL AL1,AL2,D,GAM,POD, &
      & RV_MY,CF_MY,D_MYIN,AL1_MY,AL2_MY,ALC,DT0LREF,DTLREF, &
      & A1_MYN, BB1_MYN, A2_MYN, BB2_MYN,DT,DTT,XRAD, &
      & TPC1, TPC2, TPC3, TPC4, TPC5, &
      & EPSDEL, EPSDEL2,DT0L, DT0I,&
      & ROR, &
      & CWHUCM,B6,B8L,B8I, &
      & DEL1,DEL2,DEL1S,DEL2S, &
      & TIMENEW,TIMEREV,SFN11,SFN12, &
      & SFNL,SFNI,B5L,B5I,B7L,B7I,DOPL,DOPI,RW,RI,QW,PW, &
      & PI,QI,DEL1N0,DEL2N0,D1N0,D2N0,DTNEWL,DTNEWL1,D1N,D2N, &
      & DEL_R1,DT0L0,DT0I0, &
      & DTNEWL0, &
      & DTNEWL2
 REAL DT_WATER_COND,DT_WATER_EVAP

 INTEGER K
 ! NEW ALGORITHM OF CONDENSATION (12.01.00)

 REAL  FF1_OLD(NKR),SUPINTW(NKR)
       DOUBLE PRECISION DSUPINTW(NKR),DD1N,DB11_MY,DAL1,DAL2
       DOUBLE PRECISION COL3,RORI,TPN,TPS,QPN,QPS,TOLD,QOLD &
      &                  ,FI1_K,FI2_K,FI3_K,FI4_K,FI5_K &
      &                  ,R1_K,R2_K,R3_K,R4_K,R5_K &
      &                  ,FI1R1,FI2R2,FI3R3,FI4R4,FI5R5 &
      &                  ,RMASSLAA,RMASSLBB,RMASSIAA,RMASSIBB &
      &                  ,ES1N,ES2N,EW1N,ARGEXP &
      &                  ,TT,QQ,PP &
      &                  ,DEL1N,DEL2N,DIV1,DIV2 &
      &                  ,OPER2,OPER3,AR1,AR2

        DOUBLE PRECISION DELMASSL1

 ! DROPLETS

         REAL R1(NKR) &
      &           ,RLEC(NKR),RO1BL(NKR) &
      &           ,FI1(NKR),FF1(NKR),PSI1(NKR) &
      &           ,B11_MY(NKR),B12_MY(NKR)

 ! WORK ARRAYS

 ! NEW ALGORITHM OF MIXED PHASE FOR EVAPORATION


  REAL DTIMEO(NKR),DTIMEL(NKR) &
      &           ,TIMESTEPD(NKR)

 ! NEW ALGORITHM (NO TYPE OF ICE)

  REAL :: FL1(NKR), sfndummy(3), R1N(NKR)
 INTEGER :: IDROP
DOUBLE PRECISION :: R1D(NKR),R1ND(NKR)

OPER2(AR1)=0.622/(0.622+0.378*AR1)/AR1
OPER3(AR1,AR2)=AR1*AR2/(0.622+0.378*AR1)

DATA AL1 /2500./, AL2 /2834./, D /0.211/ &
      &      ,GAM /1.E-4/, POD /10./

DATA RV_MY,CF_MY,D_MYIN,AL1_MY,AL2_MY &
      &      /461.5,0.24E-1,0.211E-4,2.5E6,2.834E6/

DATA A1_MYN, BB1_MYN, A2_MYN, BB2_MYN &
      &      /2.53,5.42,3.41E1,6.13/

  DATA TPC1, TPC2, TPC3, TPC4, TPC5 &
      &      /-4.0,-8.1,-12.7,-17.8,-22.4/


  DATA EPSDEL, EPSDEL2 /0.1E-03,0.1E-03/

  DATA DT0L, DT0I /1.E20,1.E20/

  DOUBLE PRECISION :: DEL1_d , DEL2_d, RW_d , PW_d, RI_d, PI_d, D1N_d, D2N_d, &
       VR1_d(NKR)

 sfndummy = 0.0
 B12_MY = 0.0
 B11_MY = 0.0

  I_ABERGERON=0
  I_BERGERON=0
  COL3=3.0*COL
 ITIME=0
 KCOND=0
 DT_WATER_COND=0.4
 DT_WATER_EVAP=0.4
 ITIME=0
 KCOND=0
 DT0LREF=0.2
 DTLREF=0.4

 NR=NKR
 NRM=NKR-1
 DT=DTCOND
 DTT=DTCOND
 XRAD=0.

  CWHUCM=0.
 XRAD=0.
 B6=CWHUCM*GAM-XRAD
 B8L=1./ROR
 B8I=1./ROR
 RORI=1./ROR

 DO KR=1,NKR
    FF1_OLD(KR)=FF1(KR)
    SUPINTW(KR)=0.0
    DSUPINTW(KR)=0.0
 ENDDO

 TPN=TT
 QPN=QQ
 DO KR=1,NKR
     FI1(KR)=FF1(KR)
 END DO

 ! WARM MP (CONDENSATION OR EVAPORATION) (BEGIN)
 TIMENEW=0.
 ITIME=0

 TOLD = TPN
 QOLD = QPN
 R1D = R1
 R1ND = R1D
 SFNL = 0.0
 SFN11 = 0.0

 56  ITIME = ITIME+1
 TIMEREV = DT-TIMENEW
 TIMEREV = DT-TIMENEW
 DEL1 = DEL1N
 DEL2 = DEL2N
 DEL1S = DEL1N
 DEL2S = DEL2N
 TPS = TPN
 QPS = QPN

 IF(ISYM1 == 1)THEN
  FL1 = 0.0
  VR1_d = VR1
  CALL JERRATE_KS &
     (R1D,TPS,PP,VR1_d,RLEC,RO1BL,B11_MY,1,1,fl1,NKR,ICEMAX)
  sfndummy(1)=SFN11
  CALL JERTIMESC_KS(FI1,R1D,SFNDUMMY,B11_MY,B8L,1,NKR,ICEMAX,COL)
  SFN11 = sfndummy(1)
 ENDIF

 SFN12 = 0.0
 SFNL = SFN11 + SFN12
 SFNI = 0.

 B5L=BB1_MY/TPS/TPS
 B5I=BB2_MY/TPS/TPS
 B7L=B5L*B6
 B7I=B5I*B6
 DOPL=1.+DEL1S
 DOPI=1.+DEL2S
 RW=(OPER2(QPS)+B5L*AL1)*DOPL*SFNL
 RI=(OPER2(QPS)+B5L*AL2)*DOPL*SFNI
 QW=B7L*DOPL
 PW=(OPER2(QPS)+B5I*AL1)*DOPI*SFNL
 PI=(OPER2(QPS)+B5I*AL2)*DOPI*SFNI
 QI=B7I*DOPI

 IF(RW.NE.RW .or. PW.NE.PW)THEN
    print*, 'NaN In ONECOND1'
    print*, 'fatal error in ONECOND1 (RW or PW are NaN), model stop'
 ENDIF

 KCOND=10
 IF(DEL1N >= 0.0D0) KCOND=11

   IF(KCOND == 11) THEN
      DTNEWL = DT
      DTNEWL = DT
      DTNEWL = AMIN1(DTNEWL,TIMEREV)
      TIMENEW = TIMENEW + DTNEWL
      DTT = DTNEWL

      IF (DTT < 0.0) then
       print*,"fatal error in ONECOND1-DEL1N>0:(DTT<0), model stop"
       stop
      ENDIF

      DEL1_d = DEL1
      DEL2_d = DEL2
      RW_d = RW
      PW_d = PW
      RI_d = RI
      PI_d = PI

      CALL JERSUPSAT_KS(DEL1_d,DEL2_d,DEL1N,DEL2N, &
                    RW_d,PW_d,RI_d,PI_d, &
                    DTT,D1N_d,D2N_d,0.0,0.0, &
                    ISYM1,ISYM2,ISYM3,ISYM4,ISYM5)
      DEL1 = DEL1_d
      DEL2 = DEL2_d
      RW = RW_d
      PW = PW_d
      RI = RI_d
      PI = PI_d
      D1N = D1N_d
      D2N = D2N_d

      IF(ISYM1 == 1)THEN
       IDROP = ISYM1
       CALL JERDFUN_KS(R1D, R1ND, B11_MY, FI1, PSI1, fl1, D1N, &
                   ISYM1, 1, 1, TPN, IDROP, FR_LIM, FRH_LIM, ICEMAX, NKR, COL, 1, Iin, Jin ,Kin, Itimestep)
      ENDIF

      IF((DEL1.GT.0.AND.DEL1N.LT.0) .AND.ABS(DEL1N).GT.EPSDEL) THEN
        print*,"fatal error in ONECOND1-1 (DEL1.GT.0.AND.DEL1N.LT.0), model stop"
        stop
      ENDIF

    ! IN CASE : KCOND.EQ.11
    ELSE

      ! EVAPORATION - ONLY WATER
      ! IN CASE : KCOND.NE.11
     DTIMEO = DT
      DTNEWL = DT
      DTNEWL = AMIN1(DTNEWL,TIMEREV)
      TIMENEW = TIMENEW + DTNEWL
      DTT = DTNEWL

      IF (DTT < 0.0) then
        print*,"fatal error in ONECOND1-DEL1N<0:(DTT<0), model stop"
        stop
      ENDIF

      DEL1_d = DEL1
      DEL2_d = DEL2
      RW_d = RW
      PW_d = PW
      RI_d = RI
      PI_d = PI
      CALL JERSUPSAT_KS(DEL1_d,DEL2_d,DEL1N,DEL2N, &
        RW_d,PW_d,RI_d,PI_d, &
        DTT,D1N_d,D2N_d,0.0,0.0, &
        ISYM1,ISYM2,ISYM3,ISYM4,ISYM5)
      DEL1 = DEL1_d
      DEL2 = DEL2_d
      RW = RW_d
      PW = PW_d
      RI = RI_d
      PI = PI_d
      D1N = D1N_d
      D2N = D2N_d

      IF(ISYM1 == 1)THEN
        IDROP = ISYM1
        CALL JERDFUN_KS(R1D, R1ND, B11_MY, &
                  FI1, PSI1, fl1, D1N, &
                  ISYM1, 1, 1, TPN, IDROP, FR_LIM, FRH_LIM, ICEMAX, NKR, COL, 1, Iin, Jin ,Kin, Itimestep)
      ENDIF

      IF((DEL1.LT.0.AND.DEL1N.GT.0) &
        .AND.ABS(DEL1N).GT.EPSDEL) THEN
         print*,"fatal error in ONECOND1-2 (DEL1.LT.0.AND.DEL1N.GT.0), model stop"
         stop
      ENDIF

    ENDIF


 RMASSLBB=0.
 RMASSLAA=0.

 ! ... before JERNEWF (ONLY WATER)
 DO K=1,NKR
  FI1_K = FI1(K)
  R1_K = R1(K)
  FI1R1 = FI1_K*R1_K*R1_K
  RMASSLBB = RMASSLBB+FI1R1
 ENDDO
 RMASSLBB = RMASSLBB*COL3*RORI
 IF(RMASSLBB.LE.0.) RMASSLBB=0.
 ! ... after JERNEWF (ONLY WATER)
 DO K=1,NKR
  FI1_K=PSI1(K)
  R1_K=R1(K)
  FI1R1=FI1_K*R1_K*R1_K
  RMASSLAA=RMASSLAA+FI1R1
 END DO
 RMASSLAA=RMASSLAA*COL3*RORI
 IF(RMASSLAA.LE.0.) RMASSLAA=0.

 DELMASSL1 = RMASSLAA - RMASSLBB
 QPN = QPS - DELMASSL1
 DAL1 = AL1
 TPN = TPS + DAL1*DELMASSL1

 IF(ABS(DAL1*DELMASSL1) > 3.0 )THEN
  print*,"ONECOND1-in(start)"
 print*,"I=",Iin,"J=",Jin,"Kin",Kin,"W",w_in,"DX",dx_in
  print*,"DELMASSL1",DELMASSL1,"DT",DTT
  print*,"DEL1N,DEL2N,DEL1,DEL2,D1N,D2N,RW,PW,RI,PI,DT"
  print*,DEL1N,DEL2N,DEL1,DEL2,D1N,D2N,RW,PW,RI,PI,DTT
  print*,"TPS",TPS,"QPS",QPS
 print*,'FI1 before',FI1,'PSI1 after',PSI1
  print*,"ONECOND1-in(end)"
  print*,"(fatal error in ONECOND1-in (ABS(DAL1*DELMASSL1) > 3.0), model stop"
  stop
 ENDIF

 ! ... SUPERSATURATION (ONLY WATER)
 ARGEXP=-BB1_MY/TPN
 ES1N=AA1_MY*DEXP(ARGEXP)
 ARGEXP=-BB2_MY/TPN
 ES2N=AA2_MY*DEXP(ARGEXP)
 EW1N=OPER3(QPN,PP)
 IF(ES1N == 0.0D0)THEN
          DEL1N=0.5
          DIV1=1.5
 ELSE
          DIV1 = EW1N/ES1N
          DEL1N = EW1N/ES1N-1.
 END IF
 IF(ES2N == 0.0D0)THEN
          DEL2N=0.5
          DIV2=1.5
 ELSE
          DEL2N = EW1N/ES2N-1.
          DIV2 = EW1N/ES2N
 END IF
 IF(ISYM1 == 1) THEN
  DO KR=1,NKR
           SUPINTW(KR)=SUPINTW(KR)+B11_MY(KR)*D1N
           DD1N=D1N
           DB11_MY=B11_MY(KR)
           DSUPINTW(KR)=DSUPINTW(KR)+DB11_MY*DD1N
  ENDDO
 ENDIF

 ! ... REPEATE TIME STEP (ONLY WATER: CONDENSATION OR EVAPORATION)
 IF(TIMENEW.LT.DT) GOTO 56

 57  CONTINUE

 IF(ISYM1 == 1) THEN
    CALL JERDFUN_NEW_KS (R1D,R1ND,SUPINTW, &
      FF1_OLD,PSI1, &
      TPN,IDROP,FR_LIM, NKR, COL,1,Iin,Jin,Kin,Itimestep)
 ENDIF ! in case ISYM1/=0

 RMASSLAA=0.0
 RMASSLBB=0.0

 DO K=1,NKR
  FI1_K=FF1_OLD(K)
  R1_K=R1(K)
  FI1R1=FI1_K*R1_K*R1_K
  RMASSLBB=RMASSLBB+FI1R1
 ENDDO
 RMASSLBB=RMASSLBB*COL3*RORI
 IF(RMASSLBB.LT.0.0) RMASSLBB=0.0

 DO K=1,NKR
  FI1_K=PSI1(K)
  R1_K=R1(K)
  FI1R1=FI1_K*R1_K*R1_K
  RMASSLAA=RMASSLAA+FI1R1
 ENDDO
 RMASSLAA=RMASSLAA*COL3*RORI
 IF(RMASSLAA.LT.0.0) RMASSLAA=0.0
 DELMASSL1 = RMASSLAA-RMASSLBB

 QPN = QOLD - DELMASSL1
 DAL1 = AL1
 TPN = TOLD + DAL1*DELMASSL1

 IF(ABS(DAL1*DELMASSL1) > 5.0 )THEN
  print*,"ONECOND1-out (start)"
  print*,"I=",Iin,"J=",Jin,"Kin",Kin,"W",w_in,"DX",dx_in
  print*,"DEL1N,DEL2N,D1N,D2N,RW,PW,RI,PI,DT"
  print*,DEL1N,DEL2N,D1N,D2N,RW,PW,RI,PI,DTT
  print*,"I=",Iin,"J=",Jin,"Kin",Kin
  print*,"TPS=",TPS,"QPS=",QPS,"delmassl1",delmassl1
  print*,"DAL1=",DAL1
  print*,RMASSLBB,RMASSLAA
  print*,"FI1",FI1
  print*,"PSI1",PSI1
  print*,"ONECOND1-out (end)"
  IF(ABS(DAL1*DELMASSL1) > 5.0 )THEN
   print*,"fatal error in ONECOND1-out (ABS(DAL1*DELMASSL1) > 5.0), model stop"
   stop
  ENDIF
 ENDIF

 ! ... SUPERSATURATION
 ARGEXP=-BB1_MY/TPN
 ES1N=AA1_MY*DEXP(ARGEXP)
 ARGEXP=-BB2_MY/TPN
 ES2N=AA2_MY*DEXP(ARGEXP)
 EW1N=OPER3(QPN,PP)
 IF(ES1N == 0.0D0)THEN
   DEL1N=0.5
   DIV1=1.5
  stop "fatal error in ONECOND1 (ES1N.EQ.0), model stop"
 ELSE
    DIV1=EW1N/ES1N
    DEL1N=EW1N/ES1N-1.
 END IF
 IF(ES2N.EQ.0)THEN
    DEL2N=0.5
    DIV2=1.5
    stop "fatal error in ONECOND1 (ES2N.EQ.0), model stop"
 ELSE
    DEL2N=EW1N/ES2N-1.
    DIV2=EW1N/ES2N
 END IF

 TT=TPN
 QQ=QPN
 DO KR=1,NKR
  FF1(KR)=PSI1(KR)
 ENDDO

 RETURN
 END SUBROUTINE ONECOND1
 ! +----------------------------------------------------------------------------+
 SUBROUTINE ONECOND2 &
       & (TT,QQ,PP,ROR  &
       & ,VR2,VR3,VR4,VR5,PSINGLE &
       & ,DEL1N,DEL2N,DIV1,DIV2 &
       & ,FF2,PSI2,R2,RIEC,RO2BL &
       & ,FF3,PSI3,R3,RSEC,RO3BL &
       & ,FF4,PSI4,R4,RGEC,RO4BL &
       & ,FF5,PSI5,R5,RHEC,RO5BL &
       & ,AA1_MY,BB1_MY,AA2_MY,BB2_MY &
       & ,C1_MEY,C2_MEY &
       & ,COL,DTCOND,ICEMAX,NKR &
       & ,ISYM1,ISYM2,ISYM3,ISYM4,ISYM5, &
        Iin,Jin,Kin,W_in,DX_in,Itimestep)

    IMPLICIT NONE

       INTEGER NKR,ICEMAX,ISYM1, Iin, Jin, Kin, Itimestep
       REAL    COL,VR2(NKR,ICEMAX),VR3(NKR),VR4(NKR) &
      &           ,VR5(NKR),PSINGLE &
      &       ,AA1_MY,BB1_MY,AA2_MY,BB2_MY &
      &       ,DTCOND,W_in,DX_in

       REAL C1_MEY,C2_MEY
       INTEGER I_MIXCOND,I_MIXEVAP,I_ABERGERON,I_BERGERON, &
      & KR,ICE,ITIME,ICM,KCOND,NR,NRM,INUC, &
      & ISYM2(ICEMAX),ISYM3,ISYM4,ISYM5,KP,KLIMIT, &
      & KM,ITER,KLIMITL,KLIMITG,KLIMITH,KLIMITI_1,KLIMITI_2,KLIMITI_3, &
      & NCRITI
       REAL AL1,AL2,D,GAM,POD, &
      & RV_MY,CF_MY,D_MYIN,AL1_MY,AL2_MY,ALC,DT0LREF,DTLREF, &
      & A1_MYN, BB1_MYN, A2_MYN, BB2_MYN,DT,DTT,XRAD, &
      & TPC1, TPC2, TPC3, TPC4, TPC5, &
      & EPSDEL, DT0L, DT0I, &
      & ROR, &
      & DEL1NUC,DEL2NUC, &
      & CWHUCM,B6,B8L,B8I,RMASSGL,RMASSGI, &
      & DEL1,DEL2,DEL1S,DEL2S, &
      & TIMENEW,TIMEREV,SFN11,SFN12, &
      & SFNL,SFNI,B5L,B5I,B7L,B7I,DOPL,DOPI,OPERQ,RW,RI,QW,PW, &
      & PI,QI,D1N0,D2N0,DTNEWL,DTNEWL1,D1N,D2N, &
      & DEL_R1,DT0L0,DT0I0,SFN31,SFN32,SFN52, &
      & SFNII1,SFN21,SFN22,DTNEWI3,DTNEWI4,DTNEWI5,DTNEWI2_1, &
      & DTNEWI2_2,DTNEWI1,DEL_R2,DEL_R4,DEL_R5,SFN41,SFN42, &
      & SNF51,DTNEWI2_3,DTNEWI2,DTNEWI_1,DTNEWI_2, &
      & DTNEWL0,DTNEWG1,DTNEWH1,DTNEWI_3, &
      & DTNEWL2,SFN51,SFNII2,DEL_R3,DTNEWI
        REAL DT_WATER_COND,DT_WATER_EVAP,DT_ICE_COND,DT_ICE_EVAP, &
      &  DT_MIX_COND,DT_MIX_EVAP,DT_MIX_BERGERON,DT_MIX_ANTIBERGERON

        INTEGER K

       DOUBLE PRECISION DD1N,DB11_MY,DAL1,DAL2
       DOUBLE PRECISION COL3,RORI,TPN,TPS,QPN,QPS,TOLD,QOLD &
      &                  ,FI1_K,FI2_K,FI3_K,FI4_K,FI5_K &
      &                  ,R1_K,R2_K,R3_K,R4_K,R5_K &
      &                  ,FI1R1,FI2R2,FI3R3,FI4R4,FI5R5 &
      &                  ,RMASSLAA,RMASSLBB,RMASSIAA,RMASSIBB &
      &                  ,ES1N,ES2N,EW1N,ARGEXP &
      &                  ,TT,QQ,PP &
      &                  ,DEL1N,DEL2N,DIV1,DIV2 &
      &                  ,OPER2,OPER3,AR1,AR2

        DOUBLE PRECISION DELTAQ1,DELMASSI1,DELMASSL1

         CHARACTER*70 CPRINT

 ! CRYSTALS

  REAL R2(NKR,ICEMAX) &
      &           ,RIEC(NKR,ICEMAX) &
      &           ,RO2BL(NKR,ICEMAX) &
      &           ,FI2(NKR,ICEMAX),PSI2(NKR,ICEMAX) &
      &           ,FF2(NKR,ICEMAX) &
      &           ,B21_MY(NKR,ICEMAX),B22_MY(NKR,ICEMAX)

 ! SNOW
         REAL R3(NKR) &
      &           ,RSEC(NKR),RO3BL(NKR) &
      &           ,FI3(NKR),FF3(NKR),PSI3(NKR) &
      &           ,B31_MY(NKR),B32_MY(NKR)

 ! GRAUPELS

         REAL R4(NKR) &
      &           ,RGEC(NKR),RO4BL(NKR) &
      &           ,FI4(NKR),FF4(NKR),PSI4(NKR) &
      &           ,B41_MY(NKR),B42_MY(NKR)

 ! HAIL
         REAL R5(NKR) &
      &           ,RHEC(NKR),RO5BL(NKR) &
      &           ,FI5(NKR),FF5(NKR),PSI5(NKR) &
      &           ,B51_MY(NKR),B52_MY(NKR)

 ! CCN

  REAL DTIMEG(NKR),DTIMEH(NKR)

  REAL DEL2D(ICEMAX),DTIMEO(NKR),DTIMEL(NKR) &

      &           ,DTIMEI_1(NKR),DTIMEI_2(NKR),DTIMEI_3(NKR) &
      &           ,SFNI1(ICEMAX),SFNI2(ICEMAX) &
      &           ,TIMESTEPD(NKR) &
      &           ,FI1REF(NKR),PSI1REF(NKR) &
      &           ,FI2REF(NKR,ICEMAX),PSI2REF(NKR,ICEMAX)&
      &           ,FCCNRREF(NKR)

  REAL :: FL1(NKR), sfndummy(3), FL3(NKR), FL4(NKR), FL5(NKR), &
      R2N(NKR,ICEMAX), R3N(NKR), R4N(NKR), R5N(NKR)
  INTEGER :: IDROP, ISYMICE
  DOUBLE PRECISION :: R2D(NKR,ICEMAX),R3D(NKR), R4D(NKR), R5D(NKR), &
        R2ND(NKR,ICEMAX),R3ND(NKR), R4ND(NKR), R5ND(NKR), &
        VR2_d(NKR,ICEMAX), VR3_d(NKR), VR4_d(NKR), VR5_d(NKR)

  OPER2(AR1)=0.622/(0.622+0.378*AR1)/AR1
  OPER3(AR1,AR2)=AR1*AR2/(0.622+0.378*AR1)

  DATA AL1 /2500./, AL2 /2834./, D /0.211/ &
      &      ,GAM /1.E-4/, POD /10./

  DATA RV_MY,CF_MY,D_MYIN,AL1_MY,AL2_MY &
      &      /461.5,0.24E-1,0.211E-4,2.5E6,2.834E6/

  DATA A1_MYN, BB1_MYN, A2_MYN, BB2_MYN &
      &      /2.53,5.42,3.41E1,6.13/

  DATA TPC1, TPC2, TPC3, TPC4, TPC5 &
      &      /-4.0,-8.1,-12.7,-17.8,-22.4/

  DATA EPSDEL/0.1E-03/

  DATA DT0L, DT0I /1.E20,1.E20/

  DOUBLE PRECISION :: DEL1_d, DEL2_d, RW_d, PW_d, RI_d, PI_d, D1N_d, D2N_d

  B22_MY = 0.0
  B32_MY = 0.0
  B42_MY = 0.0
  B52_MY = 0.0

  B21_MY = 0.0
  B31_MY = 0.0
  B41_MY = 0.0
  B51_MY = 0.0

  SFNDUMMY = 0.0
  R2D = R2
  R3D = R3
  R4D = R4
  R5D = R5
  R2ND = R2D
  R3ND = R3D
  R4ND = R4D
  R5ND = R5D

  SFNI1 = 0.0
  SFN31 = 0.0
  SFN41 = 0.0
  SFN51 = 0.0

  I_MIXCOND=0
  I_MIXEVAP=0
  I_ABERGERON=0
  I_BERGERON=0
  COL3=3.0*COL
  ICM=ICEMAX
  ITIME=0
  KCOND=0
  DT_WATER_COND=0.4
  DT_WATER_EVAP=0.4
  DT_ICE_COND=0.4
  DT_ICE_EVAP=0.4
  DT_MIX_COND=0.4
  DT_MIX_EVAP=0.4
  DT_MIX_BERGERON=0.4
  DT_MIX_ANTIBERGERON=0.4
  ICM=ICEMAX
  ITIME=0
  KCOND=0
  DT0LREF=0.2
  DTLREF=0.4

  NR=NKR
  NRM=NKR-1
  DT=DTCOND
  DTT=DTCOND
  XRAD=0.

  CWHUCM=0.
  XRAD=0.
  B6=CWHUCM*GAM-XRAD
  B8L=1./ROR
  B8I=1./ROR
  RORI=1./ROR

  TPN=TT
  QPN=QQ

    DO ICE=1,ICEMAX
    SFNI1(ICE)=0.
    SFNI2(ICE)=0.
    DEL2D(ICE)=0.
    ENDDO

    TIMENEW = 0.
    ITIME = 0

 ! ONLY ICE (CONDENSATION OR EVAPORATION) :

   46 ITIME = ITIME + 1

    TIMEREV=DT-TIMENEW

    DEL1=DEL1N
    DEL2=DEL2N
    DEL1S=DEL1N
    DEL2S=DEL2N
    DEL2D(1)=DEL2N
    DEL2D(2)=DEL2N
    DEL2D(3)=DEL2N
    TPS=TPN
    QPS=QPN
    DO KR=1,NKR
    FI3(KR)=PSI3(KR)
    FI4(KR)=PSI4(KR)
    FI5(KR)=PSI5(KR)
    DO ICE=1,ICEMAX
    FI2(KR,ICE)=PSI2(KR,ICE)
    ENDDO
    ENDDO

    IF(sum(ISYM2) > 0) THEN
      FL1 = 0.0
      VR2_d = VR2
    ! ... ice crystals
     CALL JERRATE_KS (R2D,TPS,PP,VR2_d,RIEC,RO2BL,B21_MY,3,2,fl1,NKR,ICEMAX)

     CALL JERTIMESC_KS (FI2,R2D,SFNI1,B21_MY,B8I,ICM,NKR,ICEMAX,COL)
    ENDIF
    IF(ISYM3 == 1) THEN
      FL3 = 0.0
      VR3_d = VR3
    ! ... snow
     CALL JERRATE_KS (R3D,TPS,PP,VR3_d,RSEC,RO3BL,B31_MY,1,3,fl3,NKR,ICEMAX)

     sfndummy(1) = SFN31
     CALL JERTIMESC_KS(FI3,R3D,SFNDUMMY,B31_MY,B8I,1,NKR,ICEMAX,COL)
       SFN31 = sfndummy(1)
    ENDIF
    IF(ISYM4 == 1) THEN
      FL4 = 0.0
      VR4_d = VR4
    ! ... graupel
     CALL JERRATE_KS(R4D,TPS,PP,VR4_d,RGEC,RO4BL,B41_MY,1,2,fl4,NKR,ICEMAX)

     sfndummy(1) = SFN41
     CALL JERTIMESC_KS(FI4,R4D,SFNDUMMY,B41_MY,B8I,1,NKR,ICEMAX,COL)
       SFN41 = sfndummy(1)
    ENDIF
    IF(ISYM5 == 1) THEN
      FL5 = 0.0
      VR5_d = VR5
    ! ... hail
     CALL JERRATE_KS(R5D,TPS,PP,VR5_d,RHEC,RO5BL,B51_MY,1,2,fl5,NKR,ICEMAX)

     sfndummy(1) = SFN51
     CALL JERTIMESC_KS(FI5,R5D,SFNDUMMY,B51_MY,B8I,1,NKR,ICEMAX,COL)
       SFN51 = sfndummy(1)
    ENDIF


    SFNII1 = SFNI1(1) + SFNI1(2) + SFNI1(3)
    SFN21 = SFNII1 + SFN31 + SFN41 + SFN51
    SFNL = 0.0
    SFN22 = 0.0
    SFNI = SFN21 + SFN22

    B5L=BB1_MY/TPS/TPS
    B5I=BB2_MY/TPS/TPS
    B7L=B5L*B6
    B7I=B5I*B6
    DOPL=1.+DEL1S
    DOPI=1.+DEL2S
    OPERQ=OPER2(QPS)
    RW=(OPERQ+B5L*AL1)*DOPL*SFNL
    QW=B7L*DOPL
    PW=(OPERQ+B5I*AL1)*DOPI*SFNL
    RI=(OPERQ+B5L*AL2)*DOPL*SFNI
    PI=(OPERQ+B5I*AL2)*DOPI*SFNI
    QI=B7I*DOPI

     KCOND=20
     IF(DEL2N > 0.0) KCOND=21

    IF(RW.NE.RW .or. PW.NE.PW)THEN
      print*, 'NaN In ONECOND2'
      stop "fatal error in ONECOND2 (RW or PW are NaN), model stop"
    ENDIF

 ! ... (ONLY ICE)
    IF(KCOND == 21)  THEN
    ! ... ONLY_ICE: CONDENSATION
       DTNEWL = DT
       DTNEWL = AMIN1(DTNEWL,TIMEREV)
       TIMENEW = TIMENEW + DTNEWL
       DTT = DTNEWL

    IF (DTT < 0.0) stop "fatal error in ONECOND2-DEL2N>0:(DTT<0), model stop"

    DEL1_d = DEL1
    DEL2_d = DEL2
    RW_d = RW
    PW_d = PW
    RI_d = RI
    PI_d = PI
    CALL JERSUPSAT_KS(DEL1_d,DEL2_d,DEL1N,DEL2N, &
              RW_d,PW_d,RI_d,PI_d, &
              DTT,D1N_d,D2N_d,0.0,0.0, &
              ISYM1,ISYM2,ISYM3,ISYM4,ISYM5)
    DEL1 = DEL1_d
    DEL2 = DEL2_d
    RW = RW_d
    PW = PW_d
    RI = RI_d
    PI = PI_d
    D1N = D1N_d
    D2N = D2N_d

    IF(sum(ISYM2) > 0)THEN
     IDROP = 0
     FL1 = 0.0
     IF(ISYM2(1) == 1) THEN
       CALL JERDFUN_KS(R2D(:,1), R2ND(:,1), B21_MY(:,1), &
           FI2(:,1), PSI2(:,1), fl1, D2N, &
           ISYM2(1), ICM, 1, TPN, IDROP, FR_LIM, FRH_LIM, ICEMAX, NKR, COL, 21, Iin, Jin ,Kin, Itimestep)
     ENDIF
     IF(ISYM2(2) == 1) THEN
       CALL JERDFUN_KS(R2D(:,2), R2ND(:,2), B21_MY(:,2), &
           FI2(:,2), PSI2(:,2), fl1, D2N, &
           ISYM2(2), ICM, 2, TPN, IDROP, FR_LIM, FRH_LIM, ICEMAX, NKR, COL, 22, Iin, Jin ,Kin, Itimestep)
     ENDIF
     IF(ISYM2(3) == 1) THEN
       CALL JERDFUN_KS(R2D(:,3), R2ND(:,3), B21_MY(:,3), &
           FI2(:,3), PSI2(:,3), fl1, D2N, &
           ISYM2(3), ICM, 3, TPN, IDROP, FR_LIM, FRH_LIM, ICEMAX, NKR, COL, 23, Iin, Jin ,Kin, Itimestep)

     ! IN CASE : ISYM2.NE.0
     ENDIF
    ENDIF

    IF(ISYM3 == 1) THEN
     IDROP = 0
     FL3 = 0.0
     CALL JERDFUN_KS(R3D, R3ND, B31_MY, &
         FI3, PSI3, fl3, D2N, &
         ISYM3, 1, 3, TPN, IDROP, FR_LIM, FRH_LIM, ICEMAX, NKR, COL, 3, Iin, Jin ,Kin, Itimestep)
    ENDIF


    IF(ISYM4 == 1) THEN
      IDROP = 0
      FL4 = 0.0
      CALL JERDFUN_KS(R4D, R4ND, B41_MY, &
         FI4, PSI4, fl4, D2N, &
         ISYM4, 1, 4, TPN, IDROP, FR_LIM, FRH_LIM, ICEMAX, NKR, COL, 4, Iin, Jin ,Kin, Itimestep)
     ! IN CASE : ISYM4.NE.0
    ENDIF

    IF(ISYM5 == 1) THEN
     IDROP = 0
     FL5 = 0.0
     CALL JERDFUN_KS(R5D, R5ND, B51_MY, &
        FI5, PSI5, fl5, D2N, &
        ISYM5, 1, 5, TPN, IDROP, FR_LIM, FRH_LIM, ICEMAX, NKR, COL, 5, Iin, Jin ,Kin, Itimestep)
    ! IN CASE : ISYM5.NE.0
    ENDIF

    IF((DEL2.GT.0.AND.DEL2N.LT.0) &
            .AND.ABS(DEL2N).GT.EPSDEL) THEN
     stop "fatal error in module_mp_fast_sbm (DEL2.GT.0.AND.DEL2N.LT.0), model stop"
    ENDIF

    ELSE
    ! ... IN CASE KCOND.NE.21
    ! ONLY ICE: EVAPORATION
        DTNEWL = DT
        DTNEWL = AMIN1(DTNEWL,TIMEREV)
        TIMENEW = TIMENEW + DTNEWL
        DTT = DTNEWL

      IF (DTT < 0.0) stop "fatal error in ONECOND2-DEL2N<0:(DTT<0), model stop"

      DEL1_d = DEL1
      DEL2_d = DEL2
      RW_d = RW
      PW_d = PW
      RI_d = RI
      PI_d = PI
      CALL JERSUPSAT_KS(DEL1_d,DEL2_d,DEL1N,DEL2N, &
               RW_d,PW_d,RI_d,PI_d, &
                DTT,D1N_d,D2N_d,0.0,0.0, &
                ISYM1,ISYM2,ISYM3,ISYM4,ISYM5)
       DEL1 = DEL1_d
      DEL2 = DEL2_d
      RW = RW_d
      PW = PW_d
      RI = RI_d
      PI = PI_d
      D1N = D1N_d
      D2N = D2N_d

    IF(sum(ISYM2) > 0) THEN
      IDROP = 0
      FL1 = 0.0
      IF(ISYM2(1)==1)THEN
       CALL JERDFUN_KS(R2D(:,1), R2ND(:,1), B21_MY(:,1), &
            FI2(:,1), PSI2(:,1), fl1, D2N, &
            ISYM2(1), ICM, 1, TPN, IDROP, FR_LIM, FRH_LIM, ICEMAX, NKR, COL, 21, Iin, Jin ,Kin, Itimestep)
        ENDIF
      IF(ISYM2(2)==1)THEN
          CALL JERDFUN_KS(R2D(:,2), R2ND(:,2), B21_MY(:,2), &
            FI2(:,2), PSI2(:,2), fl1, D2N, &
         ISYM2(2), ICM, 2, TPN, IDROP, FR_LIM, FRH_LIM, ICEMAX, NKR, COL, 22, Iin, Jin ,Kin, Itimestep)
         ENDIF
      IF(ISYM2(3)==1)THEN
          CALL JERDFUN_KS(R2D(:,3), R2ND(:,3), B21_MY(:,3), &
         FI2(:,3), PSI2(:,3), fl1, D2N, &
            ISYM2(3), ICM, 3, TPN, IDROP, FR_LIM, FRH_LIM, ICEMAX, NKR, COL, 23, Iin, Jin ,Kin, Itimestep)
     ENDIF
    ENDIF

       IF(ISYM3 == 1) THEN
    ! ... SNOW
     IDROP = 0
     FL3 = 0.0
     CALL JERDFUN_KS(R3D, R3ND, B31_MY, &
         FI3, PSI3, fl3, D2N, &
         ISYM3, 1, 3, TPN, IDROP, FR_LIM, FRH_LIM, ICEMAX, NKR, COL, 3, Iin, Jin ,Kin, Itimestep)
    ! IN CASE : ISYM3.NE.0
       ENDIF

     IF(ISYM4 == 1) THEN
     ! ... GRAUPELS (ONLY_ICE: EVAPORATION)
         ! ... New JERDFUN
         IDROP = 0
         FL4 = 0.0
         CALL JERDFUN_KS(R4D, R4ND, B41_MY, &
                         FI4, PSI4, fl4, D2N, &
                         ISYM4, 1, 4, TPN, IDROP, FR_LIM, FRH_LIM, ICEMAX, NKR, COL, 4, Iin, Jin ,Kin, Itimestep)
     ! IN CASE : ISYM4.NE.0
     ENDIF

       IF(ISYM5 == 1) THEN
         ! ... HAIL (ONLY_ICE: EVAPORATION)
           ! ... New JERDFUN
           IDROP = 0
           FL5 = 0.0
           CALL JERDFUN_KS(R5D, R5ND, B51_MY, &
                           FI5, PSI5, fl5, D2N, &
                           ISYM5, 1, 5, TPN, IDROP, FR_LIM, FRH_LIM, ICEMAX, NKR, COL, 5, Iin, Jin ,Kin, Itimestep)
             ! IN CASE : ISYM5.NE.0
       ENDIF

       IF((DEL2.LT.0.AND.DEL2N.GT.0) &
            .AND.ABS(DEL2N).GT.EPSDEL) THEN
             stop "fatal error in module_mp_fast_sbm (DEL2.LT.0.AND.DEL2N.GT.0), model stop"
       ENDIF

    ! IN CASE : KCOND.NE.21
     ENDIF

 ! MASSES
    RMASSIBB=0.0
    RMASSIAA=0.0

    DO K=1,NKR
    DO ICE = 1,ICEMAX
    FI2_K = FI2(K,ICE)
    R2_K = R2(K,ICE)
    FI2R2 = FI2_K*R2_K*R2_K
    RMASSIBB = RMASSIBB + FI2R2
     ENDDO
    FI3_K=FI3(K)
    FI4_K=FI4(K)
    FI5_K=FI5(K)
    R3_K=R3(K)
    R4_K=R4(K)
    R5_K=R5(K)
    FI3R3=FI3_K*R3_K*R3_K
    FI4R4=FI4_K*R4_K*R4_K
    FI5R5=FI5_K*R5_K*R5_K
    RMASSIBB=RMASSIBB+FI3R3
    RMASSIBB=RMASSIBB+FI4R4
    RMASSIBB=RMASSIBB+FI5R5
    ENDDO
    RMASSIBB=RMASSIBB*COL3*RORI
    IF(RMASSIBB.LT.0.0) RMASSIBB=0.0

    DO K=1,NKR
    DO ICE =1,ICEMAX
    FI2_K=PSI2(K,ICE)
    R2_K=R2(K,ICE)
    FI2R2=FI2_K*R2_K*R2_K
    RMASSIAA=RMASSIAA+FI2R2
    ENDDO
    FI3_K = PSI3(K)
    FI4_K = PSI4(K)
    FI5_K = PSI5(K)
    R3_K=R3(K)
    R4_K=R4(K)
    R5_K=R5(K)
    FI3R3=FI3_K*R3_K*R3_K
    FI4R4=FI4_K*R4_K*R4_K
    FI5R5=FI5_K*R5_K*R5_K
    RMASSIAA=RMASSIAA+FI3R3
    RMASSIAA=RMASSIAA+FI4R4
    RMASSIAA=RMASSIAA+FI5R5
    ENDDO
   RMASSIAA = RMASSIAA*COL3*RORI

   IF(RMASSIAA.LT.0.0) RMASSIAA=0.0

   DELMASSI1 = RMASSIAA-RMASSIBB
   QPN = QPS-DELMASSI1
   DAL2 = AL2
   TPN = TPS+DAL2*DELMASSI1

    IF(ABS(DAL2*DELMASSI1) > 5.0 )THEN
      print*,"ONECOND2-out (start)"
      print*,"I=",Iin,"J=",Jin,"Kin",Kin,"W",w_in,"DX",dx_in
      print*,"DEL1N,DEL2N,D1N,D2N,RW,PW,RI,PI,DT"
      print*,DEL1N,DEL2N,D1N,D2N,RW,PW,RI,PI,DTT
      print*,"TPS=",TPS,"QPS=",QPS,"delmassi1",delmassi1
      print*,"DAL1=",DAL2
      print*,RMASSIBB,RMASSIAA
      print*,"FI2_1",FI2(:,1)
      print*,"FI2_2",FI2(:,2)
      print*,"FI2_3",FI2(:,3)
      print*,"FI3",FI3
      print*,"FI4",FI4
      print*,"FI5",FI5
      print*,"PSI2_1",PSI2(:,1)
      print*,"PSI2_2",PSI2(:,2)
      print*,"PSI2_3",PSI2(:,3)
      print*,"PSI3",PSI3
      print*,"PSI4",PSI4
      print*,"PSI5",PSI5
      print*,"ONECOND2-out (end)"
      IF(ABS(DAL2*DELMASSI1) > 5.0 )THEN
      stop "fatal error in ONECOND2-out (ABS(DAL2*DELMASSI1) > 5.0), model stop"
   ENDIF
    ENDIF

 ! ... SUPERSATURATION
    ARGEXP=-BB1_MY/TPN
    ES1N=AA1_MY*DEXP(ARGEXP)
    ARGEXP=-BB2_MY/TPN
    ES2N=AA2_MY*DEXP(ARGEXP)
    EW1N=OPER3(QPN,PP)
    IF(ES1N == 0.0)THEN
     DEL1N=0.5
     DIV1=1.5
     stop "fatal error in ONECOND2 (ES1N.EQ.0), model stop"
    ELSE
     DIV1=EW1N/ES1N
     DEL1N=EW1N/ES1N-1.
    END IF
    IF(ES2N == 0.0)THEN
     DEL2N=0.5
     DIV2=1.5
     stop "fatal error in ONECOND2 (ES2N.EQ.0), model stop"
    ELSE
     DEL2N=EW1N/ES2N-1.
     DIV2=EW1N/ES2N
    END IF

 !  END OF TIME SPLITTING
 ! (ONLY ICE: CONDENSATION OR EVAPORATION)
   IF(TIMENEW.LT.DT) GOTO 46

   TT=TPN
   QQ=QPN
   DO KR=1,NKR
    DO ICE=1,ICEMAX
     FF2(KR,ICE)=PSI2(KR,ICE)
    ENDDO
    FF3(KR)=PSI3(KR)
    FF4(KR)=PSI4(KR)
    FF5(KR)=PSI5(KR)
   ENDDO

   RETURN
   END SUBROUTINE ONECOND2
 ! +----------------------------------------------------------------------------+
         SUBROUTINE ONECOND3 &
        & (TT,QQ,PP,ROR &
        & ,VR1,VR2,VR3,VR4,VR5,PSINGLE &
        & ,DEL1N,DEL2N,DIV1,DIV2 &
        & ,FF1,PSI1,R1,RLEC,RO1BL &
        & ,FF2,PSI2,R2,RIEC,RO2BL &
        & ,FF3,PSI3,R3,RSEC,RO3BL &
        & ,FF4,PSI4,R4,RGEC,RO4BL &
        & ,FF5,PSI5,R5,RHEC,RO5BL &
        & ,AA1_MY,BB1_MY,AA2_MY,BB2_MY &
        & ,C1_MEY,C2_MEY &
        & ,COL,DTCOND,ICEMAX,NKR &
        & ,ISYM1,ISYM2,ISYM3,ISYM4,ISYM5, &
         Iin,Jin,Kin,W_in,DX_in, Itimestep)

        IMPLICIT NONE
        INTEGER ICEMAX,NKR,KR,ITIME,ICE,KCOND,K &
      &           ,ISYM1,ISYM2(ICEMAX),ISYM3,ISYM4,ISYM5, Kin, Iin, Jin, Itimestep
        INTEGER KLIMITL,KLIMITG,KLIMITH,KLIMITI_1, &
      &  KLIMITI_2,KLIMITI_3
        INTEGER I_MIXCOND,I_MIXEVAP,I_ABERGERON,I_BERGERON
        REAL ROR,VR1(NKR),VR2(NKR,ICEMAX),VR3(NKR),VR4(NKR) &
      &           ,VR5(NKR),PSINGLE &
      &           ,AA1_MY,BB1_MY,AA2_MY,BB2_MY &
      &           ,C1_MEY,C2_MEY &
      &           ,COL,DTCOND,W_in,DX_in

 ! DROPLETS

         REAL R1(NKR)&
      &           ,RLEC(NKR),RO1BL(NKR) &
      &           ,FI1(NKR),FF1(NKR),PSI1(NKR) &
      &           ,B11_MY(NKR),B12_MY(NKR)

 ! CRYSTALS

  REAL R2(NKR,ICEMAX) &
      &           ,RIEC(NKR,ICEMAX) &
      &           ,RO2BL(NKR,ICEMAX) &
      &           ,FI2(NKR,ICEMAX),PSI2(NKR,ICEMAX) &
      &           ,FF2(NKR,ICEMAX) &
      &           ,B21_MY(NKR,ICEMAX),B22_MY(NKR,ICEMAX) &
      &           ,RATE2(NKR,ICEMAX),DEL_R2M(NKR,ICEMAX)

 ! SNOW
         REAL R3(NKR) &
      &           ,RSEC(NKR),RO3BL(NKR) &
      &           ,FI3(NKR),FF3(NKR),PSI3(NKR) &
      &           ,B31_MY(NKR),B32_MY(NKR) &
      &           ,DEL_R3M(NKR)

 ! GRAUPELS

         REAL R4(NKR) &
      &           ,RGEC(NKR),RO4BL(NKR) &
      &           ,FI4(NKR),FF4(NKR),PSI4(NKR) &
      &           ,B41_MY(NKR),B42_MY(NKR) &
      &           ,DEL_R4M(NKR)

 ! HAIL
         REAL R5(NKR) &
      &           ,RHEC(NKR),RO5BL(NKR) &
      &           ,FI5(NKR),FF5(NKR),PSI5(NKR) &
      &           ,B51_MY(NKR),B52_MY(NKR) &
      &           ,DEL_R5M(NKR)

       DOUBLE PRECISION DD1N,DB11_MY,DAL1,DAL2
       DOUBLE PRECISION COL3,RORI,TPN,TPS,QPN,QPS,TOLD,QOLD &
      &                  ,FI1_K,FI2_K,FI3_K,FI4_K,FI5_K &
      &                  ,R1_K,R2_K,R3_K,R4_K,R5_K &
      &                  ,FI1R1,FI2R2,FI3R3,FI4R4,FI5R5 &
      &                  ,RMASSLAA,RMASSLBB,RMASSIAA,RMASSIBB &
      &                  ,ES1N,ES2N,EW1N,ARGEXP &
      &                  ,TT,QQ,PP,DEL1N0,DEL2N0 &
      &                  ,DEL1N,DEL2N,DIV1,DIV2 &
      &                  ,OPER2,OPER3,AR1,AR2

        DOUBLE PRECISION DELTAQ1,DELMASSI1,DELMASSL1

        REAL A1_MYN, BB1_MYN, A2_MYN, BB2_MYN
         DATA A1_MYN, BB1_MYN, A2_MYN, BB2_MYN &
      &      /2.53,5.42,3.41E1,6.13/
        REAL B8L,B8I,SFN11,SFN12,SFNL,SFNI
        REAL B5L,B5I,B7L,B7I,B6,DOPL,DEL1S,DEL2S,DOPI,RW,QW,PW, &
      &  RI,PI,QI,SFNI1(ICEMAX),SFNI2(ICEMAX),AL1,AL2
        REAL D1N,D2N,DT0L, DT0I,D1N0,D2N0
        REAL SFN21,SFN22,SFNII1,SFNII2,SFN31,SFN32,SFN41,SFN42,SFN51, &
      &  SFN52
        REAL DEL1,DEL2
        REAL  TIMEREV,DT,DTT,TIMENEW
        REAL DTIMEG(NKR),DTIMEH(NKR),totccn_before,totccn_after

        REAL DEL2D(ICEMAX),DTIMEO(NKR),DTIMEL(NKR) &
      &           ,DTIMEI_1(NKR),DTIMEI_2(NKR),DTIMEI_3(NKR)
        REAL DT_WATER_COND,DT_WATER_EVAP,DT_ICE_COND,DT_ICE_EVAP, &
      &  DT_MIX_COND,DT_MIX_EVAP,DT_MIX_BERGERON,DT_MIX_ANTIBERGERON
        REAL DTNEWL0,DTNEWL1,DTNEWI1,DTNEWI2_1,DTNEWI2_2,DTNEWI2_3, &
      & DTNEWI2,DTNEWI_1,DTNEWI_2,DTNEWI3,DTNEWI4,DTNEWI5, &
      & DTNEWL,DTNEWL2,DTNEWG1,DTNEWH1
        REAL TIMESTEPD(NKR)

        DATA AL1 /2500./, AL2 /2834./
        REAL EPSDEL,EPSDEL2
        DATA EPSDEL, EPSDEL2 /0.1E-03,0.1E-03/

     REAL :: FL1(NKR), FL2(NKR,ICEMAX), FL3(NKR), FL4(NKR), FL5(NKR), SFNDUMMY(3), &
          R1N(NKR), R2N(NKR,ICEMAX), R3N(NKR), R4N(NKR), R5N(NKR)
     INTEGER :: IDROP, ICM, ISYMICE
     DOUBLE PRECISION :: R1D(NKR),R2D(NKR,ICEMAX),R3D(NKR), R4D(NKR), R5D(NKR), &
           R1ND(NKR),R2ND(NKR,ICEMAX),R3ND(NKR), R4ND(NKR), R5ND(NKR)


     DATA DT0L, DT0I /1.E20,1.E20/

     DOUBLE PRECISION :: DEL1_d, DEL2_d , RW_d, PW_d , RI_d , PI_d , D1N_d, D2N_d, &
     VR1_d(NKR), VR2_d(NKR,ICEMAX), VR3_d(NKR), VR4_d(NKR), VR5_d(NKR), &
     TTinput,QQinput,DEL1Ninput,DEL2Ninput

        OPER2(AR1)=0.622/(0.622+0.378*AR1)/AR1
        OPER3(AR1,AR2)=AR1*AR2/(0.622+0.378*AR1)



 TTinput = TT
 QQinput = QQ
 DEL1Ninput = DEL1N
 DEL2Ninput = DEL2N

 B12_MY = 0.0
 B22_MY = 0.0
 B32_MY = 0.0
 B42_MY = 0.0
 B52_MY = 0.0

 B21_MY = 0.0
 B31_MY = 0.0
 B41_MY = 0.0
 B51_MY = 0.0

 ICM = ICEMAX
 R1D = R1
 R2D = R2
 R3D = R3
 R4D = R4
 R5D = R5
 R1ND = R1D
 R2ND = R2D
 R3ND = R3D
 R4ND = R4D
 R5ND = R5D

 VR1_d = VR1
 VR2_d = VR2
 VR3_d = VR3
 VR4_d = VR4
 VR5_d = VR5

 SFN11 = 0.0
 SFNI1 = 0.0
 SFN31 = 0.0
 SFN41 = 0.0
 SFN51 = 0.0

 DT_WATER_COND=0.4
 DT_WATER_EVAP=0.4
 DT_ICE_COND=0.4
 DT_ICE_EVAP=0.4
 DT_MIX_COND=0.4
 DT_MIX_EVAP=0.4
 DT_MIX_BERGERON=0.4
 DT_MIX_ANTIBERGERON=0.4

 I_MIXCOND=0
 I_MIXEVAP=0
 I_ABERGERON=0
 I_BERGERON=0

 ITIME = 0
 TIMENEW = 0.0
 DT = DTCOND
 DTT = DTCOND

 B6=0.
 B8L=1./ROR
 B8I=1./ROR

 RORI=1.D0/ROR
  COL3=3.D0*COL
 TPN=TT
 QPN=QQ

 16  ITIME = ITIME + 1
 IF((TPN-273.15).GE.-0.187) GO TO 17
 TIMEREV = DT - TIMENEW
 DEL1 = DEL1N
 DEL2 = DEL2N
 DEL1S = DEL1N
 DEL2S = DEL2N

 DEL2D(1) = DEL2N
 DEL2D(2) = DEL2N
 DEL2D(3) = DEL2N
 TPS = TPN
 QPS = QPN
 DO KR = 1,NKR
  FI1(KR) = PSI1(KR)
  FI3(KR) = PSI3(KR)
  FI4(KR) = PSI4(KR)
  FI5(KR) = PSI5(KR)
  DO ICE = 1,ICEMAX
   FI2(KR,ICE) = PSI2(KR,ICE)
  ENDDO
 ENDDO

 IF(ISYM1 == 1)THEN
  FL1 = 0.0
  CALL JERRATE_KS &
         (R1D,TPS,PP,VR1_d,RLEC,RO1BL,B11_MY,1,1,fl1,NKR,ICEMAX)

  sfndummy(1) = SFN11
  CALL JERTIMESC_KS(FI1,R1D,SFNDUMMY,B11_MY,B8L,1,NKR,ICEMAX,COL)
  SFN11 = sfndummy(1)
 ENDIF

 IF(sum(ISYM2) > 0) THEN
   FL1 = 0.0
   ! ... ice crystals
    CALL JERRATE_KS (R2D,TPS,PP,VR2_d,RIEC,RO2BL,B21_MY,3,2,fl1,NKR,ICEMAX)
    CALL JERTIMESC_KS (FI2,R2D,SFNI1,B21_MY,B8I,ICM,NKR,ICEMAX,COL)
 ENDIF
 IF(ISYM3 == 1) THEN
   FL3 = 0.0
   ! ... snow
   CALL JERRATE_KS (R3D,TPS,PP,VR3_d,RSEC,RO3BL,B31_MY,1,3,fl3,NKR,ICEMAX)
   sfndummy(1) = SFN31
   CALL JERTIMESC_KS(FI3,R3D,SFNDUMMY,B31_MY,B8I,1,NKR,ICEMAX,COL)
    SFN31 = sfndummy(1)
 ENDIF
 IF(ISYM4 == 1) THEN
   FL4 = 0.0
   ! ... graupel
   CALL JERRATE_KS(R4D,TPS,PP,VR4_d,RGEC,RO4BL,B41_MY,1,2,fl4,NKR,ICEMAX)
   sfndummy(1) = SFN41
   CALL JERTIMESC_KS(FI4,R4D,SFNDUMMY,B41_MY,B8I,1,NKR,ICEMAX,COL)
   SFN41 = sfndummy(1)
 ENDIF
 IF(ISYM5 == 1) THEN
   FL5 = 0.0
   ! ... hail
   CALL JERRATE_KS(R5D,TPS,PP,VR5_d,RHEC,RO5BL,B51_MY,1,2,fl5,NKR,ICEMAX)
   sfndummy(1) = SFN51
   CALL JERTIMESC_KS(FI5,R5D,SFNDUMMY,B51_MY,B8I,1,NKR,ICEMAX,COL)
   SFN51 = sfndummy(1)
 ENDIF

  SFNII1 = SFNI1(1) + SFNI1(2) + SFNI1(3)
  SFN21 = SFNII1 + SFN31 + SFN41 + SFN51
  SFN12 = 0.0
  SFNL = SFN11 + SFN12
  SFN22 = 0.0
  SFNI = SFN21 + SFN22

  B5L=BB1_MY/TPS/TPS
  B5I=BB2_MY/TPS/TPS
  B7L=B5L*B6
  B7I=B5I*B6
  DOPL=1.+DEL1S
  DOPI=1.+DEL2S
  RW=(OPER2(QPS)+B5L*AL1)*DOPL*SFNL
  QW=B7L*DOPL
  PW=(OPER2(QPS)+B5I*AL1)*DOPI*SFNL
  RI=(OPER2(QPS)+B5L*AL2)*DOPL*SFNI
  PI=(OPER2(QPS)+B5I*AL2)*DOPI*SFNI
  QI=B7I*DOPI

  IF(RW.NE.RW .or. PW.NE.PW)THEN
    print*, 'NaN In ONECOND3'
    stop "fatal error in ONECOND3 (RW or PW are NaN), model stop"
  ENDIF

  ! DEL1 > 0, DEL2 < 0    (ANTIBERGERON MIXED PHASE - KCOND=50)
  ! DEL1 < 0 AND DEL2 < 0 (EVAPORATION MIXED_PHASE - KCOND=30)
  ! DEL1 > 0 AND DEL2 > 0 (CONDENSATION MIXED PHASE - KCOND=31)
  ! DEL1 < 0, DEL2 > 0    (BERGERON MIXED PHASE - KCOND=32)

  KCOND=50
  IF(DEL1N .LT. 0.0 .AND. DEL2N .LT. 0.0) KCOND=30
  IF(DEL1N .GT. 0.0 .AND. DEL2N .GT. 0.0) KCOND=31
  IF(DEL1N .LT. 0.0 .AND. DEL2N .GT. 0.0) KCOND=32

  IF(KCOND == 50) THEN
   DTNEWL = DT
    DTNEWL = AMIN1(DTNEWL,TIMEREV)
    TIMENEW = TIMENEW + DTNEWL
    DTT = DTNEWL

   ! ... Incase the Anti-Bregeron regime we do not call diffusional-growth
   PRINT*, "Anti-Bregeron Regime, No DIFFU"
   PRINT*,  DEL1, DEL2, TT, QQ, Kin
   GO TO 17
    ! IN CASE : KCOND = 50
  ENDIF
  IF(KCOND == 31) THEN
  ! ... DEL1 > 0 AND DEL2 > 0 (CONDENSATION MIXED PHASE - KCOND=31)
  ! ... CONDENSATION MIXED PHASE (BEGIN)
   DTNEWL = DT
    DTNEWL = AMIN1(DTNEWL,TIMEREV)
    TIMENEW = TIMENEW + DTNEWL
    DTT = DTNEWL
  ! CONDENSATION MIXED PHASE (END)
 ! IN CASE : KCOND = 31
  ENDIF
   IF(KCOND == 30) THEN
   ! ... DEL1 < 0 AND DEL2 < 0 (EVAPORATION MIXED_PHASE - KCOND=30)
   ! ... EVAPORATION MIXED PHASE (BEGIN)
   DTNEWL = DT
    DTNEWL = AMIN1(DTNEWL,TIMEREV)
    TIMENEW = TIMENEW + DTNEWL
    DTT = DTNEWL
  ! EVAPORATION MIXED PHASE (END)
  ! IN CASE : KCOND = 30
  ENDIF
  IF(KCOND == 32) THEN
   ! ... IF(DEL1N < 0.0 .AND. DEL2N > 0.0) KCOND=32
   ! ... BERGERON MIXED PHASE (BEGIN)
   DTNEWL = DT
    DTNEWL = AMIN1(DTNEWL,TIMEREV)
    TIMENEW = TIMENEW + DTNEWL
    DTT = DTNEWL
  ! BERGERON MIXED PHASE (END)
  ! IN CASE : KCOND = 32
  ENDIF

   IF (DTT < 0.0) stop "fatal error in ONECOND3:(DTT<0), model stop"

  DEL1_d = DEL1
  DEL2_d = DEL2
  RW_d = RW
  PW_d = PW
  RI_d = RI
  PI_d = PI
  CALL JERSUPSAT_KS(DEL1_d,DEL2_d,DEL1N,DEL2N, &
       RW_d,PW_d,RI_d,PI_d, &
       DTT,D1N_d,D2N_d,0.0,0.0, &
       ISYM1,ISYM2,ISYM3,ISYM4,ISYM5)
  DEL1 = DEL1_d
  DEL2 = DEL2_d
  RW = RW_d
  PW = PW_d
  RI = RI_d
  PI = PI_d
  D1N = D1N_d
  D2N = D2N_d

  IF(ISYM1 == 1) THEN
   ! DROPLETS
   ! DROPLET DISTRIBUTION FUNCTION
   IDROP = ISYM1
   FL1 = 0.0
   CALL JERDFUN_KS(R1D, R1ND, B11_MY, &
       FI1, PSI1, fl1, D1N, &
       ISYM1, 1, 1, TPN, IDROP, FR_LIM, FRH_LIM, ICEMAX, NKR, COL, 1, Iin, Jin ,Kin, Itimestep)
   ! IN CASE ISYM1.NE.0
  ENDIF
  IF(sum(ISYM2) > 0) THEN
   ! CRYSTALS
   IDROP = 0
   FL1 = 0.0
   IF(ISYM2(1)==1)THEN
    CALL JERDFUN_KS(R2D(:,1), R2ND(:,1), B21_MY(:,1), &
           FI2(:,1), PSI2(:,1), fl1, D2N, &
        ISYM2(1), ICM, 1, TPN, IDROP, FR_LIM, FRH_LIM, ICEMAX, NKR, COL, 21, Iin, Jin ,Kin, Itimestep)
   ENDIF
   IF(ISYM2(2)==1)THEN
      CALL JERDFUN_KS(R2D(:,2), R2ND(:,2), B21_MY(:,2), &
           FI2(:,2), PSI2(:,2), fl1, D2N, &
        ISYM2(2), ICM, 2, TPN, IDROP, FR_LIM, FRH_LIM, ICEMAX, NKR, COL, 22, Iin, Jin ,Kin, Itimestep)
   ENDIF
   IF(ISYM2(3)==1)THEN
      CALL JERDFUN_KS(R2D(:,3), R2ND(:,3), B21_MY(:,3), &
           FI2(:,3), PSI2(:,3), fl1, D2N, &
        ISYM2(3), ICM, 3, TPN, IDROP, FR_LIM, FRH_LIM, ICEMAX, NKR, COL, 23, Iin, Jin ,Kin, Itimestep)
   ENDIF
  ENDIF

  IF(ISYM3 == 1) THEN
   ! SNOW
   IDROP = 0
   FL3 = 0.0
   CALL JERDFUN_KS(R3D, R3ND, B31_MY, &
       FI3, PSI3, fl3, D2N, &
       ISYM3, 1, 3, TPN, IDROP, FR_LIM, FRH_LIM, ICEMAX, NKR, COL, 3, Iin, Jin ,Kin, Itimestep)
  ! IN CASE ISYM3.NE.0
  ENDIF

  IF(ISYM4 == 1) THEN
  ! GRAUPELS
   IDROP = 0
   FL4 = 0.0
   CALL JERDFUN_KS(R4D, R4ND, B41_MY, &
       FI4, PSI4, fl4, D2N, &
       ISYM4, 1, 4, TPN, IDROP, FR_LIM, FRH_LIM, ICEMAX, NKR, COL, 4, Iin, Jin ,Kin, Itimestep)

  ! IN CASE ISYM4.NE.0
  ENDIF

    IF(ISYM5 == 1) THEN
     ! HAIL
    IDROP = 0
    FL5 = 0.0
    CALL JERDFUN_KS(R5D, R5ND, B51_MY, &
      FI5, PSI5, fl5, D2N, &
      ISYM5, 1, 5, TPN, IDROP, FR_LIM, FRH_LIM, ICEMAX, NKR, COL, 5, Iin, Jin ,Kin, Itimestep)
  ! IN CASE ISYM5.NE.0
  ENDIF

 RMASSLBB=0.D0
 RMASSIBB=0.D0
 RMASSLAA=0.D0
 RMASSIAA=0.D0

 DO K=1,NKR
  FI1_K=FI1(K)
  R1_K=R1(K)
  FI1R1=FI1_K*R1_K*R1_K
  RMASSLBB=RMASSLBB+FI1R1
  DO ICE =1,ICEMAX
  FI2_K=FI2(K,ICE)
  R2_K=R2(K,ICE)
  FI2R2=FI2_K*R2_K*R2_K
  RMASSIBB=RMASSIBB+FI2R2
  ENDDO
   FI3_K=FI3(K)
   FI4_K=FI4(K)
   FI5_K=FI5(K)
   R3_K=R3(K)
   R4_K=R4(K)
   R5_K=R5(K)
   FI3R3=FI3_K*R3_K*R3_K
   FI4R4=FI4_K*R4_K*R4_K
   FI5R5=FI5_K*R5_K*R5_K
   RMASSIBB=RMASSIBB+FI3R3
   RMASSIBB=RMASSIBB+FI4R4
   RMASSIBB=RMASSIBB+FI5R5
   ENDDO
   RMASSIBB=RMASSIBB*COL3*RORI
   IF(RMASSIBB.LT.0.0) RMASSIBB=0.0
   RMASSLBB=RMASSLBB*COL3*RORI
   IF(RMASSLBB.LT.0.0) RMASSLBB=0.0
   DO K=1,NKR
   FI1_K=PSI1(K)
   R1_K=R1(K)
   FI1R1=FI1_K*R1_K*R1_K
   RMASSLAA=RMASSLAA+FI1R1
   DO ICE =1,ICEMAX
   FI2(K,ICE)=PSI2(K,ICE)
   FI2_K=FI2(K,ICE)
   R2_K=R2(K,ICE)
   FI2R2=FI2_K*R2_K*R2_K
   RMASSIAA=RMASSIAA+FI2R2
   ENDDO
   FI3_K=PSI3(K)
   FI4_K=PSI4(K)
   FI5_K=PSI5(K)
   R3_K=R3(K)
   R4_K=R4(K)
   R5_K=R5(K)
   FI3R3=FI3_K*R3_K*R3_K
   FI4R4=FI4_K*R4_K*R4_K
   FI5R5=FI5_K*R5_K*R5_K
   RMASSIAA=RMASSIAA+FI3R3
   RMASSIAA=RMASSIAA+FI4R4
   RMASSIAA=RMASSIAA+FI5R5
   ENDDO
  RMASSIAA=RMASSIAA*COL3*RORI
  IF(RMASSIAA.LE.0.0) RMASSIAA=0.0
  RMASSLAA=RMASSLAA*COL3*RORI
  IF(RMASSLAA.LT.0.0) RMASSLAA=0.0

  DELMASSL1=RMASSLAA-RMASSLBB
  DELMASSI1=RMASSIAA-RMASSIBB
  DELTAQ1=DELMASSL1+DELMASSI1
  QPN=QPS-DELTAQ1
  DAL1=AL1
  DAL2=AL2
  TPN = TPS + DAL1*DELMASSL1+DAL2*DELMASSI1

  IF(ABS(DAL1*DELMASSL1+DAL2*DELMASSI1) > 5.0 )THEN
   print*,"ONECOND3-input-start"
   print*,"TTinput",TTinput,"QQinput",QQinput,"PP",PP
   print*,'DEL1Ninput',DEL1Ninput,'DEL2Ninput',DEL2Ninput
   print*,"ROR",ROR,'VR1',VR1,'PSINGLE',PSINGLE
   print*,'DIV1',DIV1,'DIV2',DIV2
   print*,'R1',R1,'RLEC',RLEC,'RO1BL',RO1BL
   print*,'const',AA1_MY,BB1_MY,AA2_MY,BB2_MY
   print*,'const',C1_MEY,C2_MEY,COL
   print*,'DTCOND',DTCOND,'ICEMAX',ICEMAX,'NKR',NKR
   print*,'ISYM1',ISYM1,'ISYM2',ISYM2,'ISYM3',ISYM3,'ISYM4',ISYM4,'ISYM5',ISYM5
   print*,Iin,Jin,Kin,W_in,DX_in
   print*,"ONECOND3-input-end"

   print*,"ONECOND3-out (start)"
   print*,"I=",Iin,"J=",Jin,"Kin",Kin,"W",w_in,"DX",dx_in
   print*,"DEL1N,DEL2N,D1N,D2N,RW,PW,RI,PI,DT"
   print*,DEL1N,DEL2N,D1N,D2N,RW,PW,RI,PI,DTT
   print*,"TPS=",TPS,"TPN=",TPN,"QPS=",QPS,"delmassl1",delmassl1,"delmassi1",delmassi1
   print*,"DAL2=",DAL2,"DAL1=",DAL1
   print*,RMASSLAA,RMASSLBB
   print*,RMASSIAA,RMASSIBB
   print*,"FI1",FI1
   print*,"FI3",FI3
   print*,"FI4",FI4
   print*,"FI5",FI5
   print*,"PSI1",PSI1
   print*,"R1D",R1D,"R1ND",R1ND
   print*,"PSI3",PSI3
   print*,"R3D",R3D,"R3ND",R3ND
   print*,"PSI4",PSI4
   print*,"R4D",R4D,"R4ND",R4ND
   print*,"PSI5",PSI5
   print*,"R5D",R5D,"R5ND",R5ND
   print*,"ONECOND3-out (end)"
   IF(ABS(DAL1*DELMASSL1+DAL2*DELMASSI1) > 5.0 )THEN
    stop "fatal error in ONECOND3-out (ABS(DAL1*DELMASSL1+DAL2*DELMASSI1) > 5.0), model stop"
   ENDIF
  ENDIF

 ! SUPERSATURATION
  ARGEXP=-BB1_MY/TPN
  ES1N=AA1_MY*DEXP(ARGEXP)
  ARGEXP=-BB2_MY/TPN
  ES2N=AA2_MY*DEXP(ARGEXP)
  EW1N=OPER3(QPN,PP)
  IF(ES1N == 0.0)THEN
   DEL1N=0.5
   DIV1=1.5
   print*,'es1n onecond3 = 0'
   stop "fatal error in ONECOND3 (ES1N.EQ.0), model stop"
  ELSE
   DIV1=EW1N/ES1N
   DEL1N=EW1N/ES1N-1.
  END IF
  IF(ES2N == 0.0)THEN
   DEL2N=0.5
   DIV2=1.5
   print*,'es2n onecond3 = 0'
   stop "fatal error in ONECOND3 (ES2N.EQ.0), model stop"
  ELSE
   DEL2N=EW1N/ES2N-1.
   DIV2=EW1N/ES2N
  END IF
  ! END OF TIME SPLITTING

  IF(TIMENEW < DT) GOTO 16
  17 CONTINUE

  TT=TPN
  QQ=QPN
  DO KR=1,NKR
     FF1(KR)=PSI1(KR)
     DO ICE=1,ICEMAX
        FF2(KR,ICE)=PSI2(KR,ICE)
     ENDDO
     FF3(KR)=PSI3(KR)
     FF4(KR)=PSI4(KR)
     FF5(KR)=PSI5(KR)
  ENDDO

   RETURN
   END SUBROUTINE ONECOND3
!===============================================================
! +---------------------------------------------+
     double precision FUNCTION POLYSVP (T,TYPE)
! ..................................
!  COMPUTE SATURATION VAPOR PRESSURE

!  POLYSVP RETURNED IN UNITS OF PA.
!  T IS INPUT IN UNITS OF K.
!  TYPE REFERS TO SATURATION WITH RESPECT TO LIQUID (0) OR ICE (1)

! REPLACE GOFF-GRATCH WITH FASTER FORMULATION FROM FLATAU ET AL. 1992, TABLE 4 (RIGHT-HAND COLUMN)
! +------------------------------------------------------------------------------------------------+

      IMPLICIT NONE

      REAL DUM
      REAL T
      INTEGER TYPE
! ice
      real a0i,a1i,a2i,a3i,a4i,a5i,a6i,a7i,a8i
      data a0i,a1i,a2i,a3i,a4i,a5i,a6i,a7i,a8i /&
 6.11147274, 0.503160820, 0.188439774e-1, &
        0.420895665e-3, 0.615021634e-5,0.602588177e-7, &
        0.385852041e-9, 0.146898966e-11, 0.252751365e-14/

! liquid
      real a0,a1,a2,a3,a4,a5,a6,a7,a8

! V1.7
      data a0,a1,a2,a3,a4,a5,a6,a7,a8 /&
 6.11239921, 0.443987641, 0.142986287e-1, &
        0.264847430e-3, 0.302950461e-5, 0.206739458e-7, &
        0.640689451e-10,-0.952447341e-13,-0.976195544e-15/
      real dt

! ICE

      IF (TYPE == 1) THEN
         POLYSVP = (10.**(-9.09718*(273.16/T-1.)-3.56654*                &
             LOG10(273.16/T)+0.876793*(1.-T/273.16)+      &
             LOG10(6.1071)))*100.0*10.0
      END IF

! LIQUID

      IF (TYPE == 0) THEN
        POLYSVP = (10.**(-7.90298*(373.16/T-1.)+                        &
              5.02808*LOG10(373.16/T)-         &
              1.3816E-7*(10**(11.344*(1.-T/373.16))-1.)+    &
              8.1328E-3*(10**(-3.49149*(373.16/T-1.))-1.)+    &
              LOG10(1013.246)))*100.0*10.0
         END IF

      RETURN
      END FUNCTION POLYSVP
! +----------------------------------------------------------+
      SUBROUTINE JERRATE_KS (xlS, &
        TP,PP, &
        Vxl,RIEC,RO1BL, &
        B11_MY, &
        ID,IN,fl1,NKR,ICEMAX)

IMPLICIT NONE
! ... Interface
INTEGER,INTENT(IN) :: ID, IN, NKR, ICEMAX
  REAL(KIND=r4size),INTENT(IN) :: RO1BL(NKR,ID),RIEC(NKR,ID),FL1(NKR)
  REAL(KIND=r4size),INTENT(INOUT) :: B11_MY(NKR,ID)
  REAL(KIND=r8size),INTENT(IN) :: PP, TP, xlS(NKR,ID),Vxl(NKR,ID)
! ... Interface
! ... Locals
INTEGER :: KR, nskin(nkr), ICE
REAL(KIND=r4size) :: VENTPLM(NKR), FD1(NKR,ICEMAX),FK1(NKR,ICEMAX), xl_MY1(NKR,ICEMAX), &
              AL1_MY(2),ESAT1(2), TPreal
REAL(KIND=r8size) :: PZERO, TZERO, CONST, D_MY, COEFF_VISCOUS, SHMIDT_NUMBER,     &
        A, B, RVT, SHMIDT_NUMBER03, XLS_KR_ICE, RO1BL_KR_ICE, VXL_KR_ICE, REINOLDS_NUMBER, &
           RESHM, VENTPL, CONSTL, DETL

REAL(KIND=r4size) :: deg01,deg03

! A1L_MY - CONSTANTS FOR "MAXWELL": MKS
REAL(KIND=r8size),parameter:: RV_MY=461.5D4, CF_MY=2.4D3, D_MYIN=0.211D0

! CGS :

! RV_MY, CM*CM/SEC/SEC/KELVIN - INDIVIDUAL GAS CONSTANT
!                               FOR WATER VAPOUR
!RV_MY=461.5D4

! D_MYIN, CM*CM/SEC - COEFFICIENT OF DIFFUSION OF WATER VAPOUR

!D_MYIN=0.211D0

! PZERO, DYNES/CM/CM - REFERENCE PRESSURE

PZERO=1.013D6

! TZERO, KELVIN - REFERENCE TEMPERATURE

TZERO=273.15D0

do kr=1,nkr
if (in==2 .and. fl1(kr)==0.0 .or. in==6 .or. in==3 .and. tp<273.15) then
   nskin(kr) = 2
else !in==1 or in==6 or lef/=0
   nskin(kr) = 1
endif
enddo

! CONSTANTS FOR CLAUSIUS-CLAPEYRON EQUATION :

! A1_MY(1),G/SEC/SEC/CM

! A1_MY(1)=2.53D12

! A1_MY(2),G/SEC/SEC/CM

! A1_MY(2)=3.41D13

! BB1_MY(1), KELVIN

! BB1_MY(1)=5.42D3

! BB1_MY(2), KELVIN

! BB1_MY(2)=6.13D3

! AL1_MY(1), CM*CM/SEC/SEC - LATENT HEAT OF VAPORIZATION

AL1_MY(1)=2.5D10

! AL1_MY(2), CM*CM/SEC/SEC - LATENT HEAT OF SUBLIMATION

AL1_MY(2)=2.834D10

! CF_MY, G*CM/SEC/SEC/SEC/KELVIN - COEFFICIENT OF
!                                  THERMAL CONDUCTIVITY OF AIR
!CF_MY=2.4D3

  DEG01=1.0/3.0
  DEG03=1.0/3.0

CONST=12.566372D0

! coefficient of diffusion

D_MY=D_MYIN*(PZERO/PP)*(TP/TZERO)**1.94D0

! coefficient of viscousity

! COEFF_VISCOUS=0.13 cm*cm/sec

        COEFF_VISCOUS=0.13D0

! Shmidt number

        SHMIDT_NUMBER=COEFF_VISCOUS/D_MY

! Constants used for calculation of Reinolds number

        A=2.0D0*(3.0D0/4.0D0/3.141593D0)**DEG01
        B=A/COEFF_VISCOUS

RVT=RV_MY*TP
  ! ESAT1(IN)=A1_MY(IN)*DEXP(-BB1_MY(IN)/TP)
  ! if (IN==1) then
  !            ESAT1(IN)=EW(TP)
  ! ELSE
  !            ESAT1(IN)=EI(TP)
  ! endif

! ... (KS) - update the saturation vapor pressure
!ESAT1(1)=EW(TP)
    !ESAT1(2)=EI(TP)
TPreal = TP
ESAT1(1) = POLYSVP(TPreal,0)
ESAT1(2) = POLYSVP(TPreal,1)

   DO KR=1,NKR
      VENTPLM(KR)=0.0D0
    ENDDO

SHMIDT_NUMBER03=SHMIDT_NUMBER**DEG03

   DO ICE=1,ID
      DO KR=1,NKR

          xlS_KR_ICE=xlS(KR,ICE)
          RO1BL_KR_ICE=RO1BL(KR,ICE)
          Vxl_KR_ICE=Vxl(KR,ICE)
! Reynolds numbers
          REINOLDS_NUMBER= &
              B*Vxl_KR_ICE*(xlS_KR_ICE/RO1BL_KR_ICE)**DEG03
          RESHM=DSQRT(REINOLDS_NUMBER)*SHMIDT_NUMBER03

          IF(REINOLDS_NUMBER<2.5D0) THEN
            VENTPL=1.0D0+0.108D0*RESHM*RESHM
            VENTPLM(KR)=VENTPL
          ELSE
            VENTPL=0.78D0+0.308D0*RESHM
            VENTPLM(KR)=VENTPL
          ENDIF

        ENDDO
! cycle by KR

! VENTPL_MAX is given in MICRO.PRM include file

     DO KR=1,NKR

        VENTPL=VENTPLM(KR)

        IF(VENTPL>VENTPL_MAX) THEN
          VENTPL=VENTPL_MAX
          VENTPLM(KR)=VENTPL
        ENDIF

        CONSTL=CONST*RIEC(KR,ICE)

        FD1(KR,ICE)=RVT/D_MY/ESAT1(nskin(kr))
        FK1(KR,ICE)=(AL1_MY(nskin(kr))/RVT-1.0D0)*AL1_MY(nskin(kr))/CF_MY/TP

        xl_MY1(KR,ICE)=VENTPL*CONSTL
        ! growth rate
        DETL=FK1(KR,ICE)+FD1(KR,ICE)
        B11_MY(KR,ICE)=xl_MY1(KR,ICE)/DETL

       ENDDO
! cycle by KR

      ENDDO
! cycle by ICE

RETURN
END SUBROUTINE JERRATE_KS

! SUBROUTINE JERRATE
! ................................................................................
SUBROUTINE JERTIMESC_KS (FI1,X1,SFN11, &
       B11_MY,CF,ID,NKR,ICEMAX,COL)

IMPLICIT NONE

! ... Interface
INTEGER,INTENT(IN) :: ID,NKR,ICEMAX
 REAL(KIND=r4size),INTENT(in) :: B11_MY(NKR,ID), FI1(NKR,ID), COL, CF
 REAL(KIND=r8size),INTENT(in) :: X1(NKR,ID)
 REAL(KIND=r4size),INTENT(out) :: SFN11(ID)
! ... Interface

! ... Locals
INTEGER :: ICE, KR
 REAL(KIND=r4size) :: SFN11S, FK, DELM, FUN, B11
! ... Locals

 DO ICE=1,ID
     SFN11S=0.0D0
     SFN11(ICE)=CF*SFN11S
   DO KR=1,NKR
! value of size distribution functions
       FK=FI1(KR,ICE)
! delta-m
       DELM=X1(KR,ICE)*3.0D0*COL
! integral's expression
        FUN=FK*DELM
! values of integrals
        B11=B11_MY(KR,ICE)
     SFN11S=SFN11S+FUN*B11
  ENDDO
! cycle by kr
! correction
  SFN11(ICE)=CF*SFN11S
  ENDDO

! cycle by ice

RETURN
END SUBROUTINE JERTIMESC_KS
! +--------------------------------------------------------+
SUBROUTINE JERSUPSAT_KS (DEL1,DEL2,DEL1N,DEL2N, &
                     RW,PW,RI,PI, &
           DT,DEL1INT,DEL2INT,DYN1,DYN2, &
           ISYM1,ISYM2,ISYM3,ISYM4,ISYM5)

IMPLICIT NONE
! ... Interface
INTEGER,INTENT(INOUT) ::ISYM1, ISYM2(:), ISYM3, ISYM4, ISYM5
  REAL(KIND=r4size),INTENT(IN) ::   DT, DYN1, DYN2
  REAL(KIND=r8size),INTENT(IN) :: DEL1, DEL2
  REAL(KIND=r8size),INTENT(INOUT) :: DEL1N,DEL2N,DEL1INT,DEL2INT,RW, PW, RI, PI
! ... Interface
! ... Locals
   INTEGER :: I, ISYMICE
    REAL(KIND=r8size) :: X, EXPM1, DETER, EXPR, EXPP, A, ALFA, BETA, GAMA, G31, G32, G2, EXPB, EXPG, &
  C11, C21, C12, C22, A1DEL1N, A2DEL1N, A3DEL1N, A4DEL1N, A1DEL1INT, A2DEL1INT, &
A3DEL1INT, A4DEL1INT, A1DEL2N, A2DEL2N, A3DEL2N , A4DEL2N, A1DEL2INT, A2DEL2INT, &
A3DEL2INT, A4DEL2INT, A5DEL2INT
! ... Locals

EXPM1(x)=x+x*x/2.0D0+x*x*x/6.0D0+x*x*x*x/24.0D0+ &
                 x*x*x*x*x/120.0D0

ISYMICE = sum(ISYM2) + ISYM3 + ISYM4 + ISYM5
IF(AMAX1(RW,PW,RI,PI)<=RW_PW_RI_PI_MIN) THEN

    RW = 0.0
    PW = 0.0
    RI = 0.0
    PI = 0.0
    ISYM1 = 0
    ISYMICE = 0

ELSE

    IF(DMAX1(RW,PW)>RW_PW_MIN) THEN

   ! ... (KS) - A zero can pass through, assign a minimum value
   IF(RW < RW_PW_MIN*RW_PW_MIN) RW = 1.0D-20
   IF(PW < RW_PW_MIN*RW_PW_MIN) PW = 1.0D-20
   ! ... (KS) ...................................................

        IF(DMAX1(PI/PW,RI/RW)<=RATIO_ICEW_MIN) THEN
   ! only water
         RI = 0.0
         PI = 0.0
         ISYMICE = 0
      ENDIF

        IF(DMIN1(PI/PW,RI/RW)>1.0D0/RATIO_ICEW_MIN) THEN
   ! only ice
         RW = 0.0
         PW = 0.0
         ISYM1 = 0
       ENDIF

    ELSE
   ! only ice
   RW = 0.0
     PW = 0.0
     ISYM1 = 0

   ENDIF
  ENDIF

 IF(ISYMICE == 0)THEN
  ISYM2 = 0
  ISYM3 = 0
  ISYM4 = 0
  ISYM5 = 0
 ENDIF

    DETER=RW*PI-PW*RI


    IF(RW==0.0 .AND. RI==0.0) THEN

          DEL1N=DEL1+DYN1*DT
          DEL2N=DEL2+DYN2*DT
          DEL1INT=DEL1*DT+DYN1*DT*DT/2.0D0
          DEL2INT=DEL2*DT+DYN2*DT*DT/2.0D0

          GOTO 100

    ENDIF

! solution of equation for supersaturation with
! different DETER values

    IF(RI==0.0) THEN
! only water                                                     (start)

      EXPR=EXP(-RW*DT)
      IF(ABS(RW*DT)>1.0E-6) THEN
        DEL1N=DEL1*EXPR+(DYN1/RW)*(1.0D0-EXPR)
        DEL2N=PW*DEL1*EXPR/RW-PW*DYN1*DT/RW- &
              PW*DYN1*EXPR/(RW*RW)+DYN2*DT+ &
              DEL2-PW*DEL1/RW+PW*DYN1/(RW*RW)
        DEL1INT=-DEL1*EXPR/RW+DYN1*DT/RW+ &
                 DYN1*EXPR/(RW*RW)+DEL1/RW-DYN1/(RW*RW)
        DEL2INT=PW*DEL1*EXPR/(-RW*RW)-PW*DYN1*DT*DT/(2.0D0*RW)+ &
                PW*DYN1*EXPR/(RW*RW*RW)+DYN2*DT*DT/2.0D0+ &
                DEL2*DT-PW*DEL1*DT/RW+PW*DYN1*DT/(RW*RW)+ &
                PW*DEL1/(RW*RW)-PW*DYN1/(RW*RW*RW)
        GOTO 100
! in case DABS(RW*DT)>1.0D-6
      ELSE

! in case DABS(RW*DT)<=1.0D-6

          EXPR=EXPM1(-RW*DT)
          DEL1N=DEL1+DEL1*EXPR+(DYN1/RW)*(0.0D0-EXPR)
          DEL2N=PW*DEL1*EXPR/RW-PW*DYN1*DT/RW- &
                   PW*DYN1*EXPR/(RW*RW)+DYN2*DT+DEL2
          DEL1INT=-DEL1*EXPR/RW+DYN1*DT/RW+DYN1*EXPR/(RW*RW)
          DEL2INT=PW*DEL1*EXPR/(-RW*RW)-PW*DYN1*DT*DT/(2.0D0*RW)+ &
                     PW*DYN1*EXPR/(RW*RW*RW)+DYN2*DT*DT/2.0D0+ &
                     DEL2*DT-PW*DEL1*DT/RW+PW*DYN1*DT/(RW*RW)
          GOTO 100

        ENDIF
! only water                                                    (end)

! in case RI==0.0D0
    ENDIF

    IF(RW==0.0) THEN
! only ice                                                    (start)

      EXPP=EXP(-PI*DT)

      IF(ABS(PI*DT)>1.0E-6) THEN

        DEL2N = DEL2*EXPP+(DYN2/PI)*(1.0D0-EXPP)
        DEL2INT = -DEL2*EXPP/PI+DYN2*DT/PI+ &
                   DYN2*EXPP/(PI*PI)+DEL2/PI-DYN2/(PI*PI)
        DEL1N = +RI*DEL2*EXPP/PI-RI*DYN2*DT/PI- &
                  RI*DYN2*EXPP/(PI*PI)+DYN1*DT+ &
                  DEL1-RI*DEL2/PI+RI*DYN2/(PI*PI)
        DEL1INT = -RI*DEL2*EXPP/(PI*PI)-RI*DYN2*DT*DT/(2.0D0*PI)+ &
                    RI*DYN2*EXPP/(PI*PI*PI)+DYN1*DT*DT/2.0D0+ &
                    DEL1*DT-RI*DEL2*DT/PI+RI*DYN2*DT/(PI*PI)+ &
                    RI*DEL2/(PI*PI)-RI*DYN2/(PI*PI*PI)
        GOTO 100
! in case DABS(PI*DT)>1.0D-6
      ELSE

! in case DABS(PI*DT)<=1.0D-6

          EXPP=EXPM1(-PI*DT)
          DEL2N=DEL2+DEL2*EXPP-EXPP*DYN2/PI
          DEL2INT=-DEL2*EXPP/PI+DYN2*DT/PI+DYN2*EXPP/(PI*PI)
          DEL1N=+RI*DEL2*EXPP/PI-RI*DYN2*DT/PI- &
                    RI*DYN2*EXPP/(PI*PI)+DYN1*DT+DEL1
          DEL1INT=-RI*DEL2*EXPP/(PI*PI)-RI*DYN2*DT*DT/(2.0D0*PI)+ &
                      RI*DYN2*EXPP/(PI*PI*PI)+DYN1*DT*DT/2.0D0+ &
                      DEL1*DT-RI*DEL2*DT/PI+RI*DYN2*DT/(PI*PI)
          GOTO 100

      ENDIF
! only ice                                                      (end)

! in case RW==0.0D0
    ENDIF

    IF(RW/=0.0 .AND. RI/=0.0) THEN

      A=(RW-PI)*(RW-PI)+4.0E0*PW*RI

      IF(A < 0.0) THEN
           PRINT*,   'IN SUBROUTINE JERSUPSAT: A < 0'
            PRINT*,   'DETER'
            PRINT 201, DETER
            PRINT*,   'RW,PW,RI,PI'
            PRINT 204, RW,PW,RI,PI
            PRINT*,   'DT,DYN1,DYN2'
            PRINT 203, DT,DYN1,DYN2
            PRINT*,   'DEL1,DEL2'
            PRINT 202, DEL1,DEL2
           PRINT*,   'STOP 1905:A < 0'
           stop "fatal error: STOP 1905:A < 0, model stop"
       ENDIF
! water and ice                                               (start)
       ALFA=DSQRT((RW-PI)*(RW-PI)+4.0D0*PW*RI)

! 5/8/04 Nir, Beta is negative to the simple solution so
! it will decay

        BETA=0.5D0*(ALFA+RW+PI)
        GAMA=0.5D0*(ALFA-RW-PI)
        G31=PI*DYN1-RI*DYN2
        G32=-PW*DYN1+RW*DYN2
        G2=RW*PI-RI*PW
        IF (G2 == 0.0D0) G2 = 1.0004d-11*1.0003d-11-1.0002d-11*1.0001e-11 ! ... (KS) - 24th,May,2016
        EXPB=DEXP(-BETA*DT)
        EXPG=DEXP(GAMA*DT)

        IF(DABS(GAMA*DT)>1.0E-6) THEN
          C11=(BETA*DEL1-RW*DEL1-RI*DEL2-BETA*G31/G2+DYN1)/ALFA
          C21=(GAMA*DEL1+RW*DEL1+RI*DEL2-GAMA*G31/G2-DYN1)/ALFA
          C12=(BETA*DEL2-PW*DEL1-PI*DEL2-BETA*G32/G2+DYN2)/ALFA
          C22=(GAMA*DEL2+PW*DEL1+PI*DEL2-GAMA*G32/G2-DYN2)/ALFA
          DEL1N=C11*EXPG+C21*EXPB+G31/G2
          DEL1INT=C11*EXPG/GAMA-C21*EXPB/BETA+(C21/BETA-C11/GAMA) &
                  +G31*DT/G2
          DEL2N=C12*EXPG+C22*EXPB+G32/G2
          DEL2INT=C12*EXPG/GAMA-C22*EXPB/BETA+(C22/BETA-C12/GAMA) &
                  +G32*DT/G2
            GOTO 100
! in case DABS(GAMA*DT)>1.0D-6
        ELSE
! in case DABS(GAMA*DT)<=1.0D-6
            IF(ABS(RI/RW)>1.0E-12) THEN
              IF(ABS(RW/RI)>1.0E-12) THEN
                ALFA=DSQRT((RW-PI)*(RW-PI)+4.0D0*PW*RI)
                BETA=0.5D0*(ALFA+RW+PI)
                GAMA=0.5D0*(ALFA-RW-PI)
              IF (GAMA == 0.0D0) GAMA=0.5D0*(2.002d-10-2.001d-10) ! ... (KS) - 24th,May,2016
                EXPG=EXPM1(GAMA*DT)
                EXPB=DEXP(-BETA*DT)

! beta/alfa could be very close to 1 that why I transform it
! remember alfa-beta=gama

                C11=(BETA*DEL1-RW*DEL1-RI*DEL2+DYN1)/ALFA
                C21=(GAMA*DEL1+RW*DEL1+RI*DEL2-GAMA*G31/G2-DYN1)/ALFA
                C12=(BETA*DEL2-PW*DEL1-PI*DEL2+DYN2)/ALFA
                C22=(GAMA*DEL2+PW*DEL1+PI*DEL2-GAMA*G32/G2-DYN2)/ALFA

                A1DEL1N=C11
                A2DEL1N=C11*EXPG
                A3DEL1N=C21*EXPB
                A4DEL1N=G31/G2*(GAMA/ALFA+(GAMA/ALFA-1.0D0)*EXPG)

                DEL1N=A1DEL1N+A2DEL1N+A3DEL1N+A4DEL1N

                A1DEL1INT=C11*EXPG/GAMA
                A2DEL1INT=-C21*EXPB/BETA
                A3DEL1INT=C21/BETA
                A4DEL1INT=G31/G2*DT*(GAMA/ALFA)

                DEL1INT=A1DEL1INT+A2DEL1INT+A3DEL1INT+A4DEL1INT

                A1DEL2N=C12
                A2DEL2N=C12*EXPG
                A3DEL2N=C22*EXPB
                A4DEL2N=G32/G2*(GAMA/ALFA+ &
                       (GAMA/ALFA-1.0D0)* &
                       (GAMA*DT+GAMA*GAMA*DT*DT/2.0D0))

                DEL2N=A1DEL2N+A2DEL2N+A3DEL2N+A4DEL2N

                A1DEL2INT=C12*EXPG/GAMA
                A2DEL2INT=-C22*EXPB/BETA
                A3DEL2INT=C22/BETA
                A4DEL2INT=G32/G2*DT*(GAMA/ALFA)
                A5DEL2INT=G32/G2*(GAMA/ALFA-1.0D0)* &
                                 (GAMA*DT*DT/2.0D0)

                DEL2INT=A1DEL2INT+A2DEL2INT+A3DEL2INT+A4DEL2INT+ &
                        A5DEL2INT

! in case DABS(RW/RI)>1D-12
              ELSE

! in case DABS(RW/RI)<=1D-12

                X=-2.0D0*RW*PI+RW*RW+4.0D0*PW*RI

                ALFA=PI*(1+(X/PI)/2.0D0-(X/PI)*(X/PI)/8.0D0)
                BETA=PI+(X/PI)/4.0D0-(X/PI)*(X/PI)/16.0D0+RW/2.0D0
                GAMA=(X/PI)/4.0D0-(X/PI)*(X/PI)/16.0D0-RW/2.0D0

                EXPG=EXPM1(GAMA*DT)
                EXPB=DEXP(-BETA*DT)

              C11=(BETA*DEL1-RW*DEL1-RI*DEL2+DYN1)/ALFA
              C21=(GAMA*DEL1+RW*DEL1+RI*DEL2-GAMA*G31/G2-DYN1)/ALFA
              C12=(BETA*DEL2-PW*DEL1-PI*DEL2+DYN2)/ALFA
              C22=(GAMA*DEL2+PW*DEL1+PI*DEL2-GAMA*G32/G2-DYN2)/ALFA

                DEL1N=C11+C11*EXPG+C21*EXPB+ &
                         G31/G2*(GAMA/ALFA+(GAMA/ALFA-1)*EXPG)
                DEL1INT=C11*EXPG/GAMA-C21*EXPB/BETA+(C21/BETA)+ &
                           G31/G2*DT*(GAMA/ALFA)
                DEL2N=C12+C12*EXPG+C22*EXPB+G32/G2*(GAMA/ALFA+ &
                        (GAMA/ALFA-1.0D0)* &
                        (GAMA*DT+GAMA*GAMA*DT*DT/2.0D0))
               DEL2INT=C12*EXPG/GAMA-C22*EXPB/BETA+ &
               (C22/BETA)+G32/G2*DT*(GAMA/ALFA)+ &
                G32/G2*(GAMA/ALFA-1.0D0)*(GAMA*DT*DT/2.0D0)

! in case DABS(RW/RI)<=1D-12
          ENDIF
! alfa/beta 2
! in case DABS(RI/RW)>1D-12

            ELSE

! in case DABS(RI/RW)<=1D-12

              X=-2.0D0*RW*PI+PI*PI+4.0D0*PW*RI

              ALFA=RW*(1.0D0+(X/RW)/2.0D0-(X/RW)*(X/RW)/8.0D0)
              BETA=RW+(X/RW)/4.0D0-(X/RW)*(X/RW)/16.0D0+PI/2.0D0
              GAMA=(X/RW)/4.0D0-(X/RW)*(X/RW)/16.0D0-PI/2.0D0

              EXPG=EXPM1(GAMA*DT)
              EXPB=DEXP(-BETA*DT)

              C11=(BETA*DEL1-RW*DEL1-RI*DEL2+DYN1)/ALFA
             C21=(GAMA*DEL1+RW*DEL1+RI*DEL2-GAMA*G31/G2-DYN1)/ALFA
              C12=(BETA*DEL2-PW*DEL1-PI*DEL2+DYN2)/ALFA
             C22=(GAMA*DEL2+PW*DEL1+PI*DEL2-GAMA*G32/G2-DYN2)/ALFA

              DEL1N=C11+C11*EXPG+C21*EXPB+ &
                    G31/G2*(GAMA/ALFA+(GAMA/ALFA-1.0D0)*EXPG)
              DEL1INT=C11*EXPG/GAMA-C21*EXPB/BETA+(C21/BETA)+ &
                      G31/G2*DT*(GAMA/ALFA)
              DEL2N=C12+C12*EXPG+C22*EXPB+G32/G2* &
                    (GAMA/ALFA+ &
                    (GAMA/ALFA-1.0D0)*(GAMA*DT+GAMA*GAMA*DT*DT/2.0D0))
           DEL2INT=C12*EXPG/GAMA-C22*EXPB/BETA+C22/BETA+ &
              G32/G2*DT*(GAMA/ALFA)+ &
              G32/G2*(GAMA/ALFA-1.0D0)*(GAMA*DT*DT/2.0D0)
! alfa/beta
! in case DABS(RI/RW)<=1D-12
     ENDIF
! in case DABS(GAMA*DT)<=1D-6
   ENDIF

! water and ice                                                 (end)

! in case ISYM1/=0.AND.ISYM2/=0

        ENDIF

 100    CONTINUE

  201 FORMAT(1X,D13.5)
  202 FORMAT(1X,2D13.5)
  203 FORMAT(1X,3D13.5)
  204 FORMAT(1X,4D13.5)

        RETURN
        END SUBROUTINE JERSUPSAT_KS

! SUBROUTINE JERSUPSAT
! ....................................................................
 SUBROUTINE JERDFUN_KS (xi,xiN,B21_MY, &
                    FI2,PSI2,fl2,DEL2N, &
                    ISYM2,IND,ITYPE,TPN,IDROP, &
                    FR_LIM,FRH_LIM,ICEMAX,NKR,COL,Ihydro,Iin,Jin,Kin,Itimestep)

 IMPLICIT NONE
! ... Interface
 INTEGER,INTENT(IN) :: ISYM2, IND, ITYPE, NKR, ICEMAX, Ihydro, Iin, Jin ,Kin, Itimestep
 INTEGER,INTENT(INOUT) :: IDROP
 REAL(kind=R4SIZE),INTENT(IN) :: B21_MY(:), FI2(:), FR_LIM(:), FRH_LIM(:), &
            DEL2N, COL
 REAL(kind=R8SIZE),INTENT(IN) :: TPN, xi(:)
 REAL(kind=R8SIZE),INTENT(INOUT) :: xiN(:)
 REAL(kind=R4SIZE),INTENT(INOUT) :: PSI2(:), FL2(:)
! ... Interface

! ... Locals
 INTEGER :: ITYP, KR, NR, ICE, K, IDSD_Negative
 REAL(kind=R8SIZE) :: FL2_NEW(NKR), FI2R(NKR), PSI2R(NKR), C, DEGREE1, DEGREE2, DEGREE3, D, RATEXI, &
                 B, A, xiR(NKR),xiNR(NKR), FR_LIM_KR
! ... Locals


 C = 2.0D0/3.0D0

 DEGREE1 = 1.0D0/3.0D0
 DEGREE2 = C
 DEGREE3 = 3.0D0/2.0D0

 IF(IND > 1) THEN
   ITYP = ITYPE
 ELSE
   ITYP = 1
 ENDIF

 DO KR=1,NKR
    PSI2R(KR) = FI2(KR)
    FI2R(KR) = FI2(KR)
 ENDDO

 NR=NKR

! new size distribution functions                             (start)

 IF(ISYM2 == 1) THEN
   IF(IND==1 .AND. ITYPE==1) THEN
! drop diffusional growth
     DO KR=1,NKR
        D=xi(KR)**DEGREE1
        RATExi=C*DEL2N*B21_MY(KR)/D
        B=xi(KR)**DEGREE2
        A=B+RATExi
        IF(A<0.0D0) THEN
          xiN(KR)=1.0D-50
        ELSE
          xiN(KR)=A**DEGREE3
        ENDIF
     ENDDO
! in case IND==1.AND.ITYPE==1
   ELSE
! in case IND/=1.OR.ITYPE/=1
          DO KR=1,NKR
             RATExi = DEL2N*B21_MY(KR)
             xiN(KR) = xi(KR) + RATExi
          ENDDO
   ENDIF

! recalculation of size distribution functions                (start)

      DO KR=1,NKR
        xiR(KR) = xi(KR)
        xiNR(KR) = xiN(KR)
       FI2R(KR) = FI2(KR)
      END DO

     IDSD_Negative = 0
   CALL JERNEWF_KS &
        (NR,xiR,FI2R,PSI2R,xiNR,ISIGN_3POINT,TPN,IDROP,NKR,COL,IDSD_Negative,Ihydro,Iin,Jin,Kin,Itimestep)
   IF(IDSD_Negative == 1)THEN
    IF(ISIGN_KO_1 == 1) THEN
     ! ... (KS) - we do not use Kovatch-Ouland as separate method
     ! CALL JERNEWF_KO_KS &
       !     (NR,xiR,FI2R,PSI2R,xiNR,NKR,COL)
    ENDIF
   ENDIF

     DO KR=1,NKR
       !IF(ITYPE==5) THEN
       !    FR_LIM_KR=FRH_LIM(KR)
       !ELSE
       !    FR_LIM_KR=FR_LIM(KR)
      !ENDIF
      IF(PSI2R(KR)<0.0D0) THEN
       PRINT*,    'STOP 1506 : PSI2R(KR)<0.0D0, in JERDFUN_KS'
      stop "fatal error in PSI2R(KR)<0.0D0, in JERDFUN_KS, model stop"
      ENDIF
        PSI2(KR) = PSI2R(KR)
     ENDDO
! cycle by ICE
! recalculation of size distribution functions                  (end)
! in case ISYM2/=0
 ENDIF
! new size distribution functions                               (end)

  201 FORMAT(1X,D13.5)
  304   FORMAT(1X,I2,2X,4D13.5)

 RETURN
 END SUBROUTINE JERDFUN_KS
! +----------------------------------------------------------------------------+
  SUBROUTINE JERNEWF_KS &
           (NRX,RR,FI,PSI,RN,I3POINT,TPN,IDROP,NKR,COL,IDSD_Negative,Ihydro, &
              Iin,Jin,Kin,Itimestep)

        IMPLICIT NONE
! ... Interface
  INTEGER,INTENT(IN) :: NRX, I3POINT, NKR, Ihydro, Iin, Jin, Kin, Itimestep
  INTEGER,INTENT(INOUT) :: IDROP, IDSD_Negative
  real(kind=R8SIZE),INTENT(IN) :: TPN
  real(kind=R4SIZE),INTENT(IN) :: COL
  real(kind=R8SIZE),INTENT(INOUT) :: PSI(:), RN(:), FI(:), RR(:)
! ... Interface

! ... Locals
  INTEGER :: KMAX, KR, I, K , NRXP, ISIGN_DIFFUSIONAL_GROWTH, NRX1,  &
              I3POINT_CONDEVAP, IEvap
  real(kind=R8SIZE) :: RNTMP,RRTMP,RRP,RRM,RNTMP2,RRTMP2,RRP2,RRM2, GN1,GN2, &
               GN3,GN1P,GMAT,GMAT2, &
        CDROP(NRX),DELTA_CDROP(NRX),RRS(NRX+1),PSINEW(NRX+1), &
        PSI_IM,PSI_I,PSI_IP, AOLDCON, ANEWCON, AOLDMASS, ANEWMASS

  INTEGER,PARAMETER :: KRDROP_REMAPING_MIN = 6, KRDROP_REMAPING_MAX = 12
! ... Locals

 IF(TPN .LT. 273.15-5.0D0) IDROP=0

! INITIAL VALUES FOR SOME VARIABLES

  NRXP = NRX + 1
!   NRX1 = 24
!   NRX1 = 35
   NRX1 = NKR

   DO I=1,NRX
! RN(I), g - new masses after condensation or evaporation
     IF(RN(I) < 0.0D0) THEN
       RN(I) = 1.0D-50
        FI(I) = 0.0D0
     ENDIF
  ENDDO

! new change 26.10.09                                         (start)
 DO K=1,NRX
    RRS(K)=RR(K)
 ENDDO
! new change 26.10.09                                           (end)

 I3POINT_CONDEVAP = I3POINT

 IEvap = 0
 IF(RN(1) < RRS(1)) THEN
! evaporation
   I3POINT_CONDEVAP = 0
! new change 26.10.09                                         (start)
   IDROP = 0
! new change 26.10.09                                           (end)
   NRX1 = NRX
   IEvap = 1
 ENDIF

 IF(IDROP == 0) I3POINT_CONDEVAP = 0

! new change 26.10.09                                         (start)

 DO K=1,NRX
    PSI(K)=0.0D0
    CDROP(K)=0.0D0
    DELTA_CDROP(K)=0.0D0
    PSINEW(K)=0.0D0
 ENDDO

 RRS(NRXP)=RRS(NRX)*1024.0D0

 PSINEW(NRXP) = 0.0D0

! new change 26.10.09                                           (end)

 ISIGN_DIFFUSIONAL_GROWTH = 0

 DO K=1,NRX
    IF(RN(K).NE.RR(K)) THEN
    ISIGN_DIFFUSIONAL_GROWTH = 1
    GOTO 2000
    ENDIF
 ENDDO

 2000   CONTINUE

 IF(ISIGN_DIFFUSIONAL_GROWTH == 1) THEN

! Kovetz-Olund method                                         (start)

! new change 26.10.09                                         (start)
   DO K=1,NRX1 ! ... [KS] >> NRX1-1
! new change 26.10.09                                           (end)

   IF(FI(K) > 0.0) THEN
     IF(DABS(RN(K)-RR(K)) < 1.0D-16) THEN
          PSINEW(K) = FI(K)*RR(K)
          CYCLE
       ENDIF

     I = 1
     DO WHILE (.NOT.(RRS(I) <= RN(K) .AND. RRS(I+1) >= RN(K)) &
                 .AND.I.LT.NRX1) ! [KS] >> was NRX1-1
                  I = I + 1
       ENDDO

       IF(RN(K).LT.RRS(1)) THEN
          RNTMP=RN(K)
          RRTMP=0.0D0
          RRP=RRS(1)
          GMAT2=(RNTMP-RRTMP)/(RRP-RRTMP)
          PSINEW(1)=PSINEW(1)+FI(K)*RR(K)*GMAT2
     ELSE

        RNTMP=RN(K)
        RRTMP=RRS(I)
        RRP=RRS(I+1)
        GMAT2=(RNTMP-RRTMP)/(RRP-RRTMP)
        GMAT=(RRP-RNTMP)/(RRP-RRTMP)
        PSINEW(I)=PSINEW(I)+FI(K)*RR(K)*GMAT
        PSINEW(I+1)=PSINEW(I+1)+FI(K)*RR(K)*GMAT2
     ENDIF
! in case FI(K).NE.0.0D0
   ENDIF

 3000    CONTINUE

   ENDDO
! cycle by K

   DO KR=1,NRX1
       PSI(KR)=PSINEW(KR)
   ENDDO

   DO KR=NRX1+1,NRX
       PSI(KR)=FI(KR)
   ENDDO
! Kovetz-Olund method                                           (end)

! calculation both new total drop concentrations(after KO) and new
! total drop masses (after KO)

! 3point method                                               (start)
   IF(I3POINT_CONDEVAP == 1) THEN
     DO K=1,NRX1-1
     IF(FI(K) > 0.0) THEN
        IF(DABS(RN(K)-RR(K)).LT.1.0D-16) THEN
            PSI(K) = FI(K)*RR(K)
            GOTO 3001
          ENDIF

          IF(RRS(2).LT.RN(K)) THEN
             I = 2
             DO WHILE &
                     (.NOT.(RRS(I) <= RN(K) .AND. RRS(I+1) >= RN(K)) &
                     .AND.I.LT.NRX1-1)
                    I = I + 1
           ENDDO
             RNTMP=RN(K)

             RRTMP=RRS(I)
             RRP=RRS(I+1)
             RRM=RRS(I-1)

             RNTMP2=RN(K+1)

             RRTMP2=RRS(I+1)
             RRP2=RRS(I+2)
             RRM2=RRS(I)

             GN1=(RRP-RNTMP)*(RRTMP-RNTMP)/(RRP-RRM)/ &
                  (RRTMP-RRM)

             GN1P=(RRP2-RNTMP2)*(RRTMP2-RNTMP2)/ &
                   (RRP2-RRM2)/(RRTMP2-RRM2)

             GN2=(RRP-RNTMP)*(RNTMP-RRM)/(RRP-RRTMP)/ &
                   (RRTMP-RRM)

             GMAT=(RRP-RNTMP)/(RRP-RRTMP)

             GN3=(RRTMP-RNTMP)*(RRM-RNTMP)/(RRP-RRM)/ &
                                           (RRP-RRTMP)
             GMAT2=(RNTMP-RRTMP)/(RRP-RRTMP)

             PSI_IM = PSI(I-1)+GN1*FI(K)*RR(K)

             PSI_I = PSI(I)+GN1P*FI(K+1)*RR(K+1)+&
                   (GN2-GMAT)*FI(K)*RR(K)

             PSI_IP = PSI(I+1)+(GN3-GMAT2)*FI(K)*RR(K)

             IF(PSI_IM > 0.0D0) THEN

               IF(PSI_IP > 0.0D0) THEN

                 IF(I > 2) THEN
! smoothing criteria
                   IF(PSI_IM > PSI(I-2) .AND. PSI_IM < PSI_I &
                     .AND. PSI(I-2) < PSI(I) .OR. PSI(I-2) >= PSI(I)) THEN

                      PSI(I-1) = PSI_IM

                      PSI(I) = PSI(I) + FI(K)*RR(K)*(GN2-GMAT)

                      PSI(I+1) = PSI_IP
! in case smoothing criteria
                   ENDIF
! in case I.GT.2
                 ENDIF

! in case PSI_IP.GT.0.0D0
        ELSE
               EXIT
             ENDIF
! in case PSI_IM.GT.0.0D0
      ELSE
            EXIT
          ENDIF
! in case I.LT.NRX1-2
!         ENDIF

! in case RRS(2).LT.RN(K)
       ENDIF

! in case FI(K).NE.0.0D0
      ENDIF

 3001 CONTINUE

     ENDDO
        ! cycle by K

      ! in case I3POINT_CONDEVAP.NE.0
   ENDIF
! 3 point method                                                (end)

! PSI(K) - new hydrometeor size distribution function

   DO K=1,NRX1
      PSI(K)=PSI(K)/RR(K)
   ENDDO

   DO K=NRX1+1,NRX
   PSI(K)=FI(K)
   ENDDO

   IF(IDROP == 1) THEN
    DO K=KRDROP_REMAPING_MIN,KRDROP_REMAPING_MAX
     CDROP(K)=3.0D0*COL*PSI(K)*RR(K)
    ENDDO
        ! KMAX - right boundary spectrum of drop sdf
         !(KRDROP_REMAP_MIN =< KMAX =< KRDROP_REMAP_MAX)
    DO K=KRDROP_REMAPING_MAX,KRDROP_REMAPING_MIN,-1
       KMAX=K
       IF(PSI(K).GT.0.0D0) GOTO 2011
    ENDDO

  2011  CONTINUE
 ! Andrei's new change 28.04.10                                (start)
    DO K=KMAX-1,KRDROP_REMAPING_MIN,-1
 ! Andrei's new change 28.04.10                                  (end)
     IF(CDROP(K).GT.0.0D0) THEN
      DELTA_CDROP(K)=CDROP(K+1)/CDROP(K)
       IF(DELTA_CDROP(K).LT.COEFF_REMAPING) THEN
        CDROP(K)=CDROP(K)+CDROP(K+1)
        CDROP(K+1)=0.0D0
       ENDIF
     ENDIF
    ENDDO

    DO K=KRDROP_REMAPING_MIN,KMAX
     PSI(K)=CDROP(K)/(3.0D0*COL*RR(K))
    ENDDO

 ! in case IDROP.NE.0
    ENDIF

! new change 26.10.09                                           (end)

! in case ISIGN_DIFFUSIONAL_GROWTH.NE.0
        ELSE
! in case ISIGN_DIFFUSIONAL_GROWTH.EQ.0
       DO K=1,NRX
          PSI(K)=FI(K)
       ENDDO
       ENDIF

  DO KR=1,NRX
      IF(PSI(KR) < 0.0) THEN ! ... (KS)
     IDSD_Negative = 1
     print*, "IDSD_Negative=",IDSD_Negative,"kr",kr
     PRINT*,    'IN SUBROUTINE JERNEWF'
     PRINT*,  'PSI(KR)<0'
     PRINT*,    'BEFORE EXIT'
     PRINT*,    'ISIGN_DIFFUSIONAL_GROWTH'
     PRINT*,     ISIGN_DIFFUSIONAL_GROWTH
     PRINT*,    'I3POINT_CONDEVAP'
     PRINT*,     I3POINT_CONDEVAP
     PRINT*,    'K,RR(K),RN(K),K=1,NRX'
     PRINT*,    (K,RR(K),RN(K),K=1,NRX)
     PRINT*,    'K,RR(K),RN(K),FI(K),PSI(K),K=1,NRX'
     PRINT 304, (K,RR(K),RN(K),FI(K),PSI(K),K=1,NRX)
     PRINT*,  IDROP,Ihydro,Iin,Jin,Kin,Itimestep
     stop "fatal error in SUBROUTINE JERNEWF PSI(KR)<0, < min, model stop"
   ENDIF
  ENDDO

  304   FORMAT(1X,I2,2X,4D13.5)

        RETURN
        END SUBROUTINE JERNEWF_KS
! +------------------------------------------------------------------+
 SUBROUTINE JERDFUN_NEW_KS &
           (xi,xiN,B21_MY, &
           FI2,PSI2, &
           TPN,IDROP,FR_LIM,NKR,COL,Ihydro,Iin,Jin,Kin,Itimestep)

 IMPLICIT NONE

! ... Interface
 INTEGER,INTENT(INOUT) :: IDROP, NKR
 INTEGER,INTENT(IN) :: Ihydro,Iin,Jin,Kin,Itimestep
 REAL(kind=R4SIZE),intent(IN) :: FI2(:), B21_MY(:), FR_LIM(:), COL
 REAL(kind=R8SIZE), INTENT(IN) :: TPN, xi(:)
 REAL(kind=R4SIZE),INTENT(INOUT) :: PSI2(:)
 REAL(kind=R8SIZE),INTENT(INOUT) :: xiN(:)
! ... Interface

! ... Locals
 INTEGER :: NR, KR, IDSD_Negative
 REAL(kind=R8SIZE) :: C, DEGREE1, DEGREE2, DEGREE3, D, RATEXI, B, A, &
                 xiR(NKR),FI2R(NKR),PSI2R(NKR),xiNR(NKR)
! ... Locals

 C=2.0D0/3.0D0

 DEGREE1=C/2.0D0
 DEGREE2=C
 DEGREE3=3.0D0/2.0D0

 NR=NKR

 xiR = xi
 FI2R = FI2
 PSI2R = PSI2
 xiNR = xiN

! new drop size distribution functions                             (start)

! drop diffusional growth

 DO KR=1,NKR
    D = xiR(KR)**DEGREE1
! Andrei's new change of 3.09.10                              (start)
!    RATExi=C*DEL2N*B21_MY(KR)/D
    RATExi = C*B21_MY(KR)/D
! Andrei's new change of 3.09.10                                (end)
    B = xiR(KR)**DEGREE2
    A = B+RATExi
    IF(A<0.0D0) THEN
      xiNR(KR) = 1.0D-50
    ELSE
      xiNR(KR) = A**DEGREE3
    ENDIF
 ENDDO

! recalculation of size distribution functions                (start)

 IDSD_Negative = 0
 CALL JERNEWF_KS &
   (NR,xiR,FI2R,PSI2R,xiNR,ISIGN_3POINT,TPN,IDROP,NKR,COL,IDSD_Negative,Ihydro,Iin,Jin,Kin,Itimestep)
 IF(IDSD_Negative == 1)THEN
  IF(ISIGN_KO_2 == 1) THEN
   ! ... (KS) - we do not use Kovatch-Ouland as separate method
    ! CALL JERNEWF_KO_KS &
      !      (NR,xiR,FI2R,PSI2R,xiNR,NKR,COL)
  ENDIF
 ENDIF

 PSI2 = PSI2R

! recalculation of drop size distribution functions                  (end)
! new drop size distribution functions                          (end)

  201 FORMAT(1X,D13.5)

 RETURN
 END SUBROUTINE JERDFUN_NEW_KS
! +---------------------------------------------------------+
 SUBROUTINE Relaxation_Time(TPS,QPS,PP,ROR,DEL1S,DEL2S, &
                      R1,VR1,FF1in,RLEC,RO1BL, &
                      R2,VR2,FF2in,RIEC,RO2BL, &
                      R3,VR3,FF3in,RSEC,RO3BL, &
                      R4,VR4,FF4in,RGEC,RO4BL, &
                      R5,VR5,FF5in,RHEC,RO5BL, &
                      NKR,ICEMAX,COL,DTdyn,NCOND,DTCOND)

 implicit none
! ... Interface
 integer,intent(in) :: NKR,ICEMAX
 integer,intent(out) :: NCOND
 real(kind=R4SIZE),intent(in) :: R1(:),FF1in(:),RLEC(:),RO1BL(:), &
        R2(:,:),FF2in(:,:),RIEC(:,:),RO2BL(:,:), &
        R3(NKR),FF3in(:),RSEC(:),RO3BL(:), &
        R4(NKR),FF4in(:),RGEC(:),RO4BL(:), &
        R5(NKR),FF5in(:),RHEC(:),RO5BL(:), &
        ROR,COL,DTdyn,VR1(:),VR2(:,:),VR3(:),VR4(:),VR5(:)
  real(kind=R8SIZE),intent(in) :: TPS,QPS,PP,DEL1S,DEL2S
  real(kind=R4SIZE),intent(out) :: DTCOND
! ... Interface
! ... Local
 integer :: ISYM1, ISYM2(ICEMAX), ISYM3, ISYM4, ISYM5, ISYM_SUM, ICM
  real(kind=R8SIZE),parameter :: AA1_MY = 2.53D12, BB1_MY = 5.42D3, AA2_MY = 3.41D13, &
                                 BB2_MY = 6.13E3, AL1 = 2500.0, AL2 = 2834.0
 real(kind=R8SIZE),parameter :: TAU_Min = 0.1 ! [s]
 real(kind=R8SIZE) :: OPER2, AR1, TAU_RELAX, B5L, B5I, &
                 R1D(NKR), R2D(NKR,ICEMAX), R3D(NKR), R4D(NKR), R5D(NKR), &
                       VR1_d(nkr),VR2_d(nkr,icemax),VR3_d(nkr),VR4_d(nkr),VR5_d(nkr)
 real(kind=R4SIZE) :: B11_MY(NKR), B21_MY(NKR,ICEMAX), B31_MY(NKR), &
                    B41_MY(NKR), B51_MY(NKR), FL1(NKR), FL3(NKR), FL4(NKR), FL5(NKR), &
                       SFNDUMMY(3), SFN11, SFNI1(ICEMAX), SFNII1, SFN21, SFN31, SFN41, SFN51, SFNI, SFNL, B8L, B8I, RI, PW, &
                      DOPL, DOPI, TAU_w, TAU_i, phi, RW, PI
! ... Local

  OPER2(AR1)=0.622/(0.622+0.378*AR1)/AR1
    VR1_d = VR1
    VR2_d = VR2
    VR3_d = VR3
    VR4_d = VR4
    VR5_d = VR5


  ISYM1 = 0
  ISYM2 = 0
  ISYM3 = 0
  ISYM4 = 0
  ISYM5 = 0
  IF(sum(FF1in) > 0.0) ISYM1 = 1
  IF(sum(FF2in(:,1)) > 1.0D-10) ISYM2(1) = 1
  IF(sum(FF2in(:,2)) > 1.0D-10) ISYM2(2) = 1
  IF(sum(FF2in(:,3)) > 1.0D-10) ISYM2(3) = 1
  IF(sum(FF3in) > 1.0D-10) ISYM3 = 1
  IF(sum(FF4in) > 1.0D-10) ISYM4 = 1
  IF(sum(FF5in) > 1.0D-10) ISYM5 = 1

  ISYM_SUM = ISYM1 + sum(ISYM2) + ISYM3 + ISYM4  + ISYM5
  IF(ISYM_SUM == 0)THEN
   TAU_RELAX = DTdyn
   NCOND = nint(DTdyn/TAU_RELAX)
      DTCOND = TAU_RELAX
        RETURN
  ENDIF

  R1D = R1
  R2D = R2
  R3D = R3
  R4D = R4
  R5D = R5
  B8L=1./ROR
    B8I=1./ROR
  ICM = ICEMAX
  SFN11 = 0.0
  SFNI1 = 0.0
  SFN31 = 0.0
  SFN41 = 0.0
  SFN51 = 0.0
  B11_MY = 0.0
  B21_MY = 0.0
  B31_MY = 0.0
  B41_MY = 0.0
  B51_MY = 0.0


    ! ... Drops
    IF(ISYM1 == 1)THEN
     FL1 = 0.0
     CALL JERRATE_KS &
              (R1D,TPS,PP,VR1_d,RLEC,RO1BL,B11_MY,1,1,fl1,NKR,ICEMAX)
     sfndummy(1) = SFN11
     SFN11 = sfndummy(1)
    ENDIF
    ! ... IC
    !IF(sum(ISYM2) > 0) THEN
    ! FL1 = 0.0
    ! ! ... ice crystals
    ! CALL JERRATE_KS (R2D,TPS,PP,VR2_d,RIEC,RO2BL,B21_MY,3,2,fl1,NKR,ICEMAX)
    ! CALL JERTIMESC_KS (FF2in,R2D,SFNI1,B21_MY,B8I,ICM,NKR,ICEMAX,COL)
    !ENDIF
      ! ... Snow
      IF(ISYM3 == 1) THEN
     FL3 = 0.0
     ! ... snow
     CALL JERRATE_KS (R3D,TPS,PP,VR3_d,RSEC,RO3BL,B31_MY,1,3,fl3,NKR,ICEMAX)
     sfndummy(1) = SFN31
     CALL JERTIMESC_KS(FF3in,R3D,SFNDUMMY,B31_MY,B8I,1,NKR,ICEMAX,COL)
       SFN31 = sfndummy(1)
      ENDIF
      ! ... Graupel
     IF(ISYM4 == 1) THEN
     FL4 = 0.0
     ! ... graupel
     CALL JERRATE_KS(R4D,TPS,PP,VR4_d,RGEC,RO4BL,B41_MY,1,2,fl4,NKR,ICEMAX)

     sfndummy(1) = SFN41
     CALL JERTIMESC_KS(FF4in,R4D,SFNDUMMY,B41_MY,B8I,1,NKR,ICEMAX,COL)
       SFN41 = sfndummy(1)
    ENDIF
      ! ... Hail
      IF(ISYM5 == 1) THEN
        FL5 = 0.0
        ! ... hail
        CALL JERRATE_KS(R5D,TPS,PP,VR5_d,RHEC,RO5BL,B51_MY,1,2,fl5,NKR,ICEMAX)

        sfndummy(1) = SFN51
        CALL JERTIMESC_KS(FF5in,R5D,SFNDUMMY,B51_MY,B8I,1,NKR,ICEMAX,COL)
        SFN51 = sfndummy(1)
     ENDIF

    SFNII1 = 0.0
    SFN21 = 0.0
    SFNL = 0.0
    SFNI = 0.0
    RI = 0.0
    PW = 0.0
    SFNII1 = SFNI1(1)+SFNI1(2)+SFNI1(3)
    SFN21 = SFNII1 + SFN31 + SFN41 + SFN51
    SFNL = SFN11  ! Liquid
    SFNI = SFN21  ! Total Ice

    B5L=BB1_MY/TPS/TPS
    B5I=BB2_MY/TPS/TPS
    DOPL=1.+ DEL1S
    DOPI=1.+ DEL2S
    RW=(OPER2(QPS)+B5L*AL1)*DOPL*SFNL
    RI=(OPER2(QPS)+B5L*AL2)*DOPL*SFNI
    PW=(OPER2(QPS)+B5I*AL1)*DOPI*SFNL
    PI=(OPER2(QPS)+B5I*AL2)*DOPI*SFNI

      TAU_w = DTdyn
      TAU_i = DTdyn
      phi = (1.0 + DEL2S)/(1.0 + DEL1S)
      if(PW > 0.0 .or. PI > 0.0) TAU_w = (PW + phi*PI)**(-1.0)
      if(RW > 0.0 .or. RI > 0.0) TAU_i =  phi/(RW + RI*phi)
      TAU_RELAX = DTdyn
    IF(PW > 0.0 .or. RI > 0.0) TAU_RELAX = (PW + RI)**(-1.0)/3.0
    IF(PW > 0.0 .and. RI > 0.0) TAU_RELAX = min(TAU_w,TAU_i)/3.0

      if(TAU_RELAX > DTdyn) TAU_RELAX = DTdyn/3.0
    if(TAU_RELAX < TAU_Min) TAU_RELAX = TAU_Min
      IF(PW <= 0.0 .and. RI <= 0.0) TAU_RELAX = DTdyn

    !if(TAU_RELAX < DTdyn .and. IDebug_Print_DebugModule==1)then
    !  print*,"in Relaxation_Time,TAU_RELAX < DTdyn"
   !   print*,TAU_RELAX
    !endif

    !NCOND = nint(DTdyn/TAU_RELAX)
    NCOND = ceiling(DTdyn/TAU_RELAX)
      DTCOND = TAU_RELAX

 RETURN
 END SUBROUTINE Relaxation_Time
SUBROUTINE BREAKINIT
!     GT    : MASS DISTRIBUTION FUNCTION
!     XT_MG : MASS OF BIN IN MG
!     JMAX  : NUMBER OF BINS
IMPLICIT NONE
INTEGER,PARAMETER :: hujisbm_unit1=22
LOGICAL, PARAMETER :: PRINT_diag=.FALSE.
INTEGER AP,IE,JE,KE
PARAMETER (AP = 1)
INTEGER I,J,K,JDIFF
REAL PI,D0
DOUBLE PRECISION M(0:JBREAK),ALM
INTEGER IP,KP,JP,KQ,JQ
CHARACTER*20 FILENAME_P,FILENAME_Q

FILENAME_P = 'coeff_p.asc'
FILENAME_Q = 'coeff_q.asc'

IE = JBREAK
JE = JBREAK
KE = JBREAK
PI    = 3.1415927
D0    = 0.0101593
M(1)  = PI/6.0 * D0**3

!.....IN CGS
!.....SHIFT BETWEEN COAGULATION AND BREAKUP GRID

JDIFF = JMAX - JBREAK

!.....INITIALIZATION
!........CALCULATING THE BREAKUP GRID
ALM  = 2.d0
M(0)  = M(1)/ALM
DO K=1,KE-1
   M(K+1) = M(K)*ALM
ENDDO
DO K=1,KE
   BRKWEIGHT(K) = 2./(M(K)**2 - M(K-1)**2)
ENDDO

!........OUTPUT
!WRITE (*,*) 'COLL_BREAKUP_INI: COAGULATION AND BREAKUP GRID'
!WRITE (*,'(2A5,5A15)') 'ICOAG','IBREAK','XCOAG','DCOAG', &
!    'XBREAK','DBREAK','MWEIGHT'

!........READ DER BREAKUP COEFFICIENTS FROM INPUT FILE

!WRITE (*,*) 'COLL_BREAKUP: READ THE BREAKUP COEFFS'
!WRITE (*,*) '              FILE PKIJ: ', FILENAME_P
OPEN(UNIT=hujisbm_unit1,FILE="./src/input_data/sbm_input/coeff_p.asc",  &
   FORM="FORMATTED",STATUS="OLD")
DO K=1,KE
   DO I=1,IE
      DO J=1,I
         READ(hujisbm_unit1,'(3I6,1E16.8)') KP,IP,JP,PKIJ(KP,IP,JP)
      ENDDO
   ENDDO
ENDDO
CLOSE(hujisbm_unit1)

OPEN(UNIT=hujisbm_unit1,FILE="./src/input_data/sbm_input/coeff_q.asc",  &
   FORM="FORMATTED",STATUS="OLD")
DO K=1,KE
   DO J=1,JE
      READ(hujisbm_unit1,'(2I6,1E16.8)') KQ,JQ,QKJ(KQ,JQ)
   ENDDO
ENDDO
CLOSE(hujisbm_unit1)

DO I=1,JMAX
   DO J=1,JMAX
      ECOALMASSM(I,J)=1.0D0
   ENDDO
ENDDO

DO I=1,JMAX
   DO J=1,JMAX
      ECOALMASSM(I,J)=ECOALMASS(XL(I),XL(J))
   ENDDO
ENDDO
RETURN
END SUBROUTINE BREAKINIT
!##############################################################################
REAL FUNCTION ECOALMASS(ETA,KSI)
IMPLICIT NONE
REAL PI
PARAMETER (PI = 3.1415927)

REAL ETA,KSI
REAL KPI,RHO
REAL DETA,DKSI

PARAMETER (RHO  = 1.0)

!     REAL ECOALDIAM
!     EXTERNAL ECOALDIAM

KPI = 6./PI

DETA = (KPI*ETA/RHO)**(1./3.)
DKSI = (KPI*KSI/RHO)**(1./3.)

ECOALMASS = ECOALDIAM(DETA,DKSI)

RETURN
END FUNCTION ECOALMASS
!#############################################################################

!------------------------------------------------
!     COALESCENCE EFFICIENCY AS FUNC OF DIAMETERS
!------------------------------------------------

REAL FUNCTION ECOALDIAM(DETA,DKSI)
!     IMPLICIT NONE

REAL DETA,DKSI
REAL DGR,DKL,RGR,RKL,P,Q,E,X,QMIN,QMAX
REAL ZERO,ONE,EPS,PI

PARAMETER (ZERO = 0.0)
PARAMETER (ONE  = 1.0)
PARAMETER (EPS  = 1.0E-30)
PARAMETER (PI   = 3.1415927)

!     REAL   ECOALLOWLIST,ECOALOCHS
!     EXTERNAL ECOALLOWLIST,ECOALOCHS

DGR = MAX(DETA,DKSI)
DKL = MIN(DETA,DKSI)

RGR = 0.5*DGR
RKL = 0.5*DKL

P = (RKL / RGR)
Q = (RKL * RGR)**0.5
Q = 0.5 * (RKL + RGR)

qmin = 250e-4
qmax = 400e-4  
if (q.lt.qmin) then
   e = max(ecoalOchs(Dgr,Dkl),ecoalBeard(Dgr,Dkl)) 
elseif (q.ge.qmin.and.q.lt.qmax) then
   x = (q - qmin) / (qmax - qmin)
   e = sin(pi/2.0*x)**2 * ecoalLowList(Dgr,Dkl) &
          + sin(pi/2.0*(1 - x))**2 * ecoalOchs(Dgr,Dkl)
elseif (q.ge.qmax) then
   e = ecoalLowList(Dgr,Dkl)
else
   e  = 1.0
endif

ECOALDIAM  = MAX(MIN(ONE,E),EPS)

RETURN
END FUNCTION  ECOALDIAM
!############################################################################
!--------------------------------------------------
!     COALESCENCE EFFICIENCY (LOW&LIST)
!--------------------------------------------------
REAL FUNCTION ECOALLOWLIST(DGR,DKL)
IMPLICIT NONE
REAL PI,SIGMA,KA,KB,EPSI
REAL DGR,DKL,RGR,RKL
REAL ST,SC,ET,DSTSC,CKE,W1,W2,DC,ECL
REAL QQ0,QQ1,QQ2

PARAMETER (EPSI=1.E-20)

PI = 3.1415927
SIGMA = 72.8
KA = 0.778
KB = 2.61E-4

RGR = 0.5*DGR
RKL = 0.5*DKL

CALL COLLENERGY(DGR,DKL,CKE,ST,SC,W1,W2,DC)

DSTSC = ST-SC
ET = CKE+DSTSC
IF (ET .LT. 50.0) THEN
   QQ0=1.0+(DKL/DGR)
   QQ1=KA/QQ0**2
   QQ2=KB*SIGMA*(ET**2)/(SC+EPSI)
   ECL=QQ1*EXP(-QQ2)
ELSE
   ECL=0.0
ENDIF

ECOALLOWLIST = ECL

RETURN
END FUNCTION ECOALLOWLIST
!###########################################################################
!--------------------------------------------------
!     COALESCENCE EFFICIENCY (BEARD AND OCHS)
!--------------------------------------------------

REAL FUNCTION ECOALOCHS(D_L,D_S)
IMPLICIT NONE
REAL D_L,D_S
REAL PI,SIGMA,N_W,R_S,R_L,DV,P,G,X,E
REAL EPSF,FPMIN
PARAMETER (EPSF  = 1.E-30)
PARAMETER (FPMIN = 1.E-30)

PI = 3.1415927
SIGMA = 72.8

R_S = 0.5 * D_S
R_L = 0.5 * D_L
P   = R_S / R_L

DV  = ABS(VTBEARD(D_L) - VTBEARD(D_S))
IF (DV.LT.FPMIN) DV = FPMIN
N_W = R_S * DV**2 / SIGMA
G   = 2**(3./2.)/(6.*PI) * P**4 * (1.+ P) / ((1.+P**2)*(1.+P**3))
X   = N_W**(0.5) * G
E   = 0.767 - 10.14 * X

ECOALOCHS = E

RETURN
END FUNCTION ECOALOCHS
!##########################################################################
!-----------------------------------------
!     CALCULATING THE COLLISION ENERGY
!-----------------------------------------

SUBROUTINE COLLENERGY(DGR,DKL,CKE,ST,SC,W1,W2,DC)
IMPLICIT NONE

REAL DGR,DKL,DC
REAL K10,PI,SIGMA,RHO
REAL CKE,W1,W2,ST,SC
REAL DGKA3,DGKB3,DGKA2
REAL V1,V2,DV
REAL EPSF,FPMIN
PARAMETER (EPSF  = 1.E-30)
PARAMETER (FPMIN = 1.E-30)

PI    = 3.1415927
RHO   = 1.0
SIGMA = 72.8

K10=RHO*PI/12.0D0

DGR = MAX(DGR,EPSF)
DKL = MAX(DKL,EPSF)

DGKA2=(DGR**2)+(DKL**2)

DGKA3=(DGR**3)+(DKL**3)

IF (DGR.NE.DKL) THEN
   V1 = VTBEARD(DGR)
   V2 = VTBEARD(DKL)
   DV = (V1-V2)
   IF (DV.LT.FPMIN) DV = FPMIN
   DV = DV**2
   IF (DV.LT.FPMIN) DV = FPMIN
   DGKB3=(DGR**3)*(DKL**3)
   CKE = K10 * DV * DGKB3/DGKA3
ELSE
   CKE = 0.0D0
ENDIF
ST = PI*SIGMA*DGKA2
SC = PI*SIGMA*DGKA3**(2./3.)

W1=CKE/(SC+EPSF)
W2=CKE/(ST+EPSF)

DC=DGKA3**(1./3.)

RETURN
END SUBROUTINE COLLENERGY
!######################################################################
!--------------------------------------------------
!     CALCULATING TERMINAL VELOCITY (BEARD-FORMULA)
!--------------------------------------------------

REAL FUNCTION VTBEARD(DIAM)
IMPLICIT NONE

REAL DIAM,AA
REAL ROP,RU,AMT,PP,RL,TT,ETA,DENS,CD,D,A
REAL ALA,GR,SI,BOND,PART,XX,YY,RE,VT
REAL B00,B11,B22,B33,B44,B55,B0,B1,B2,B3,B4,B5,B6
INTEGER ID

DATA B00,B11,B22,B33,B44,B55,B0,B1,B2,B3,B4,B5,B6/-5.00015, &
     5.23778,-2.04914,.475294,-.0542819,.00238449,-3.18657,.992696, &
     -.153193E-2,-.987059E-3,-.578878E-3,.855176E-4,-.327815E-5/

AA   = DIAM/2.0
ROP  = 1.0
RU   = 8.3144E+7
AMT  = 28.9644
ID   = 10000
PP   = FLOAT(ID)*100.
RL   = RU/AMT
TT   = 283.15
ETA  = (1.718+.0049*(TT-273.15))*1.E-4
DENS = PP/TT/RL
ALA  = 6.6E-6*1.01325E+6/PP*TT/293.15
GR   = 979.69
SI   = 76.1-.155*(TT-273.15)

IF (AA.GT.500.E-4) THEN
   BOND = GR*(ROP-DENS)*AA*AA/SI
   PART = (SI**3*DENS*DENS/(ETA**4*GR*(ROP-DENS)))**(1./6.)
   XX = LOG(16./3.*BOND*PART)
   YY = B00+B11*XX+B22*XX*XX+B33*XX**3+B44*XX**4+B55*XX**5
   RE = PART*EXP(YY)
   VT = ETA*RE/2./DENS/AA
ELSEIF (AA.GT.1.E-3) THEN
   CD = 32.*AA*AA*AA*(ROP-DENS)*DENS*GR/3./ETA/ETA
   XX = LOG(CD)
   RE = EXP(B0+B1*XX+B2*XX*XX+B3*XX**3+B4*XX**4+B5*XX**5+B6*XX**6)
   D  = CD/RE/24.-1.
   VT = ETA*RE/2./DENS/AA
ELSE
   A  = 1.+1.26*ALA/AA
   A  = A*2.*AA*AA*GR*(ROP-DENS)/9./ETA
   CD = 12*ETA/A/AA/DENS
   VT = A
ENDIF

VTBEARD = VT

RETURN
END FUNCTION VTBEARD
!########################################################################
!-------------------------------------------------- 
!     Function f. Coalescence-Efficiency 
!     Eq. (7) of Beard and Ochs (1995)
!--------------------------------------------------
REAL FUNCTION ecoalBeard(D_l,D_s) 
 
IMPLICIT NONE 
REAL D_l,D_s
REAL R_s,R_l
REAL rcoeff
REAL epsf
PARAMETER (epsf  = 1.e-30) 
INTEGER its
COMPLEX acoeff(4),x

R_s = 0.5 * D_s
R_l = 0.5 * D_l

rcoeff = 5.07 - log(R_s*1e4) - log(R_l*1e4/200.0)

acoeff(1) = CMPLX(rcoeff)
acoeff(2) = CMPLX(-5.94)
acoeff(3) = CMPLX(+7.27)
acoeff(4) = CMPLX(-5.29)

x = (0.50,0)

CALL LAGUER(acoeff,3,x,its)

EcoalBeard = REAL(x)

RETURN 
END FUNCTION ecoalBeard 
!######################################################################

SUBROUTINE laguer(a,m,x,its)
INTEGER m,its,MAXIT,MR,MT
REAL EPSS
COMPLEX a(m+1),x
PARAMETER (EPSS=2.e-7,MR=8,MT=10,MAXIT=MT*MR)
INTEGER iter,j
REAL abx,abp,abm,err,frac(MR)
COMPLEX dx,x1,b,d,f,g,h,sq,gp,gm,g2
SAVE frac
DATA frac /.5,.25,.75,.13,.38,.62,.88,1./
do iter=1,MAXIT
   its=iter
   b=a(m+1)
   err=abs(b)
   d=cmplx(0.,0.)
   f=cmplx(0.,0.)
   abx=abs(x)
   do j=m,1,-1
      f=x*f+d
      d=x*d+b
      b=x*b+a(j)
      err=abs(b)+abx*err
   enddo
   err=EPSS*err
   if(abs(b).le.err) then
      return
   else
      g=d/b
      g2=g*g
      h=g2-2.*f/b
      sq=sqrt((m-1)*(m*h-g2))
      gp=g+sq
      gm=g-sq
      abp=abs(gp)
      abm=abs(gm)
      if(abp.lt.abm) gp=gm
      if (max(abp,abm).gt.0.) then
         dx=m/gp
      else
         dx=exp(cmplx(log(1.+abx),float(iter)))
      endif
   endif
   x1=x-dx
   if(x.eq.x1)return
   if (mod(iter,MT).ne.0) then
      x=x1
   else
      x=x-dx*frac(iter/MT)
   endif
enddo
return
END SUBROUTINE laguer
!##################################################################
subroutine courant_bott
implicit none
integer k,kk,j,i
double precision x0
! ima(i,j) - k-category number,
! chucm(i,j)   - courant number :
! logarithmic grid distance(dlnr) :
xl_mg(0)=xl_mg(1)/2
do i=1,nkr
   do j=i,nkr
      x0=xl_mg(i)+xl_mg(j)
      do k=j,nkr
         kk=k
         if(xl_mg(k).ge.x0.and.xl_mg(k-1).lt.x0) then
           chucm(i,j)=dlog(x0/xl_mg(k-1))/(3.d0*dlnr)
           if(chucm(i,j).gt.1.-1.d-08) then
             chucm(i,j)=0.
             kk=kk+1
           endif
           ima(i,j)=min(nkr-1,kk-1)

           goto 2000
         endif
      enddo
 2000       continue
      chucm(j,i)=chucm(i,j)
      ima(j,i)=ima(i,j)
   enddo
enddo
return
end subroutine courant_bott
!#####################################################################

SUBROUTINE KERNALS(DTIME)
! KHAIN30/07/99
IMPLICIT NONE
INTEGER I,J
REAL PI
!******************************************************************
data pi/3.141592654/
! dtime - timestep of integration (calculated in main program) :
!Adele - dtime is now the number of time steps between calls to 
!COAL_BOTT_NEW. Kernals are multiplied by the timestep length
!in KERNALSDT
! dlnr - logarithmic grid distance
! ima(i,j) - k-category number, c(i,j) - courant number 
! cw*(i,j) (in cm**3) - multiply help kernel with constant 
! timestep(dt) and logarithmic grid distance(dlnr) :
REAL DTIME
! logarithmic grid distance(dlnr) :
! dlnr=dlog(2.d0)/(3.d0*scal)
! scal is micro.prm file parameter(scal=1.d0 for x(k+1)=x(k)*2)
! calculation of cw*(i,j) (in cm**3) - multiply help kernel 
! with constant timestep(dtime) and logarithmic grid distance(dlnr) :
!     print*,'dlnr in kernal = ',dlnr,dtime
  DO I=1,NKR
     DO J=1,NKR
        CWLL_1000MB(I,J)=DTIME*DLNR*YWLL_1000MB(I,J)
        CWLL_750MB(I,J)=DTIME*DLNR*YWLL_750MB(I,J)
        CWLL_500MB(I,J)=DTIME*DLNR*YWLL_500MB(I,J)

        CWLL(I,J)=DTIME*DLNR*YWLL(I,J)
        CWLS(I,J)=DTIME*DLNR*YWLS(I,J)
        CWLG(I,J)=DTIME*DLNR*YWLG(I,J)
        CWLH(I,J)=DTIME*DLNR*YWLH(I,J)

        CWSL(I,J)=DTIME*DLNR*YWSL(I,J)
        CWSS(I,J)=DTIME*DLNR*YWSS(I,J)
        CWSG(I,J)=DTIME*DLNR*YWSG(I,J)
        CWSH(I,J)=DTIME*DLNR*YWSH(I,J)

        CWGL(I,J)=DTIME*DLNR*YWGL(I,J)
        IF(RADXXO(I,6).LT.2.0D-2) THEN
          IF(RADXXO(J,1).LT.1.0D-3) THEN
            IF(RADXXO(J,1).GE.7.0D-4) THEN
              CWGL(I,J)=DTIME*DLNR*YWGL(I,J)/1.5D0
            ELSE
              CWGL(I,J)=DTIME*DLNR*YWGL(I,J)/3.0D0
            ENDIF
          ENDIF
        ENDIF
        IF(I.LE.14.AND.J.LE.7) CWGL(I,J)=0.0D0
        CWGS(I,J)=DTIME*DLNR*YWGS(I,J)
        CWGG(I,J)=DTIME*DLNR*YWGG(I,J)
        CWGH(I,J)=DTIME*DLNR*YWGH(I,J)

        CWHL(I,J)=DTIME*DLNR*YWHL(I,J)
        CWHS(I,J)=DTIME*DLNR*YWHS(I,J)
        CWHG(I,J)=DTIME*DLNR*YWHG(I,J)
        CWHH(I,J)=DTIME*DLNR*YWHH(I,J)

        CWLI_1(I,J)=DTIME*DLNR*YWLI(I,J,1)
        CWLI_2(I,J)=DTIME*DLNR*YWLI(I,J,2)
        CWLI_3(I,J)=DTIME*DLNR*YWLI(I,J,3)
        
        CWIL_1(I,J)=DTIME*DLNR*YWIL(I,J,1)
        CWIL_2(I,J)=DTIME*DLNR*YWIL(I,J,2)
        CWIL_3(I,J)=DTIME*DLNR*YWIL(I,J,3)

        CWIS_1(I,J)=DTIME*DLNR*YWIS(I,J,1)
        CWIS_2(I,J)=DTIME*DLNR*YWIS(I,J,2)
        CWIS_3(I,J)=DTIME*DLNR*YWIS(I,J,3)

        CWSI_1(I,J)=DTIME*DLNR*YWSI(I,J,1)
        CWSI_2(I,J)=DTIME*DLNR*YWSI(I,J,2)
        CWSI_3(I,J)=DTIME*DLNR*YWSI(I,J,3)

        CWIG_1(I,J)=DTIME*DLNR*YWIG(I,J,1)
        CWIG_2(I,J)=DTIME*DLNR*YWIG(I,J,2)
        CWIG_3(I,J)=DTIME*DLNR*YWIG(I,J,3)

        CWGI_1(I,J)=DTIME*DLNR*YWGI(I,J,1)
        CWGI_2(I,J)=DTIME*DLNR*YWGI(I,J,2)
        CWGI_3(I,J)=DTIME*DLNR*YWGI(I,J,3)

        CWIH_1(I,J)=DTIME*DLNR*YWIH(I,J,1)
        CWIH_2(I,J)=DTIME*DLNR*YWIH(I,J,2)
        CWIH_3(I,J)=DTIME*DLNR*YWIH(I,J,3)

        CWHI_1(I,J)=DTIME*DLNR*YWHI(I,J,1)
        CWHI_2(I,J)=DTIME*DLNR*YWHI(I,J,2)
        CWHI_3(I,J)=DTIME*DLNR*YWHI(I,J,3)

        CWII_1_1(I,J)=DTIME*DLNR*YWII(I,J,1,1)
        CWII_1_2(I,J)=DTIME*DLNR*YWII(I,J,1,2)
        CWII_1_3(I,J)=DTIME*DLNR*YWII(I,J,1,3)

        CWII_2_1(I,J)=DTIME*DLNR*YWII(I,J,2,1)
        CWII_2_2(I,J)=DTIME*DLNR*YWII(I,J,2,2)
        CWII_2_3(I,J)=DTIME*DLNR*YWII(I,J,2,3)

        CWII_3_1(I,J)=DTIME*DLNR*YWII(I,J,3,1)
        CWII_3_2(I,J)=DTIME*DLNR*YWII(I,J,3,2)
        CWII_3_3(I,J)=DTIME*DLNR*YWII(I,J,3,3)
     ENDDO
  ENDDO
  CALL TURBCOEF
  DO J=1,7
     DO I=15,24-J
        CWGL(I,J)=0.0D0
     ENDDO
  ENDDO
  DO I=1,NKR
     DO J=1,NKR
        CWLG(J,I)=CWGL(I,J)
     ENDDO
  ENDDO
  DO I=KRMING_GL,KRMAXG_GL
     DO J=KRMINL_GL,KRMAXL_GL
       IF (ICETURB.EQ.1)THEN
         CWGL(I,J)=CTURBGL(I,J)*CWGL(I,J)
       ELSE
         CWGL(I,J)=CWGL(I,J)
       END IF
     ENDDO
  ENDDO
  DO I=KRMING_GL,KRMAXG_GL
     DO J=KRMINL_GL,KRMAXL_GL
        CWLG(J,I)=CWGL(I,J)
     ENDDO
  ENDDO

RETURN
END SUBROUTINE KERNALS
!-----------------------------------------------------------
SUBROUTINE KERNALSDT
use parameters, only:dt
implicit none
  CWLL_1000MB=CWLL_1000MB*DT
  CWLL_750MB=CWLL_750MB*DT
  CWLL_500MB=CWLL_500MB*DT

  CWLL=CWLL*DT
  CWLS=CWLS*DT
  CWLG=CWLG*DT
  CWLH=CWLH*DT

  CWSL=CWSL*DT
  CWSS=CWSS*DT
  CWSG=CWSG*DT
  CWSH=CWSH*DT

  CWGL=CWGL*DT
  CWGS=CWGS*DT
  CWGG=CWGG*DT
  CWGH=CWGH*DT

  CWHL=CWHL*DT
  CWHS=CWHS*DT
  CWHG=CWHG*DT
  CWHH=CWHH*DT

  CWLI_1=CWLI_1*DT
  CWLI_2=CWLI_2*DT
  CWLI_3=CWLI_3*DT
  
  CWIL_1=CWIL_1*DT
  CWIL_2=CWIL_2*DT
  CWIL_3=CWIL_3*DT

  CWIS_1=CWIS_1*DT
  CWIS_2=CWIS_2*DT
  CWIS_3=CWIS_3*DT

  CWSI_1=CWSI_1*DT
  CWSI_2=CWSI_2*DT
  CWSI_3=CWSI_3*DT

  CWIG_1=CWIG_1*DT
  CWIG_2=CWIG_2*DT
  CWIG_3=CWIG_3*DT

  CWGI_1=CWGI_1*DT
  CWGI_2=CWGI_2*DT
  CWGI_3=CWGI_3*DT

  CWIH_1=CWIH_1*DT
  CWIH_2=CWIH_2*DT
  CWIH_3=CWIH_3*DT

  CWHI_1=CWHI_1*DT
  CWHI_2=CWHI_2*DT
  CWHI_3=CWHI_3*DT

  CWII_1_1=CWII_1_1*DT
  CWII_1_2=CWII_1_2*DT
  CWII_1_3=CWII_1_3*DT

  CWII_2_1=CWII_2_1*DT
  CWII_2_2=CWII_2_2*DT
  CWII_2_3=CWII_2_3*DT

  CWII_3_1=CWII_3_1*DT
  CWII_3_2=CWII_3_2*DT
  CWII_3_3=CWII_3_3*DT

END SUBROUTINE KERNALSDT
!#############################################################3
SUBROUTINE TURBCOEF
IMPLICIT NONE
INTEGER I,J
DOUBLE PRECISION X_KERN,Y_KERN
DOUBLE PRECISION RL_LL(K0_LL),RL_GL(K0L_GL),RG_GL(K0G_GL)
RL_LL(1)=RADXXO(KRMIN_LL,1)*1.E4
RL_LL(2)=10.0D0
RL_LL(3)=20.0D0
RL_LL(4)=30.0D0
RL_LL(5)=40.0D0
RL_LL(6)=50.0D0
RL_LL(7)=60.0D0
RL_LL(8)=RADXXO(KRMAX_LL,1)*1.E4
DO J=1,K0_LL
   DO I=1,K0_LL
      CTURB_LL(I,J)=1.0D0
   ENDDO
ENDDO 
  CTURB_LL(1,1)=4.50D0
  CTURB_LL(1,2)=4.50D0
  CTURB_LL(1,3)=3.00D0
  CTURB_LL(1,4)=2.25D0
  CTURB_LL(1,5)=1.95D0
  CTURB_LL(1,6)=1.40D0
  CTURB_LL(1,7)=1.40D0
  CTURB_LL(1,8)=1.40D0

  CTURB_LL(2,1)=4.50D0
  CTURB_LL(2,2)=4.50D0
  CTURB_LL(2,3)=3.00D0
  CTURB_LL(2,4)=2.25D0
  CTURB_LL(2,5)=1.95D0
  CTURB_LL(2,6)=1.40D0
  CTURB_LL(2,7)=1.40D0
  CTURB_LL(2,8)=1.40D0

  CTURB_LL(3,1)=3.00D0
  CTURB_LL(3,2)=3.00D0
  CTURB_LL(3,3)=2.70D0
  CTURB_LL(3,4)=2.25D0
  CTURB_LL(3,5)=1.65D0
  CTURB_LL(3,6)=1.40D0
  CTURB_LL(3,7)=1.40D0
  CTURB_LL(3,8)=1.40D0

  CTURB_LL(4,1)=2.25D0
  CTURB_LL(4,2)=2.25D0
  CTURB_LL(4,3)=2.25D0
  CTURB_LL(4,4)=1.95D0
  CTURB_LL(4,5)=1.65D0
  CTURB_LL(4,6)=1.40D0
  CTURB_LL(4,7)=1.40D0
  CTURB_LL(4,8)=1.40D0

  CTURB_LL(5,1)=1.95D0
  CTURB_LL(5,2)=1.95D0
  CTURB_LL(5,3)=1.65D0
  CTURB_LL(5,4)=1.65D0
  CTURB_LL(5,5)=1.65D0
  CTURB_LL(5,6)=1.40D0
  CTURB_LL(5,7)=1.40D0
  CTURB_LL(5,8)=1.40D0

  CTURB_LL(6,1)=1.40D0
  CTURB_LL(6,2)=1.40D0
  CTURB_LL(6,3)=1.40D0
  CTURB_LL(6,4)=1.40D0
  CTURB_LL(6,5)=1.40D0
  CTURB_LL(6,6)=1.40D0
  CTURB_LL(6,7)=1.40D0
  CTURB_LL(6,8)=1.40D0

  CTURB_LL(7,1)=1.40D0
  CTURB_LL(7,2)=1.40D0
  CTURB_LL(7,3)=1.40D0
  CTURB_LL(7,4)=1.40D0
  CTURB_LL(7,5)=1.40D0
  CTURB_LL(7,6)=1.40D0
  CTURB_LL(7,7)=1.40D0
  CTURB_LL(7,8)=1.40D0

  CTURB_LL(8,1)=1.40D0
  CTURB_LL(8,2)=1.40D0
  CTURB_LL(8,3)=1.40D0
  CTURB_LL(8,4)=1.40D0
  CTURB_LL(8,5)=1.40D0
  CTURB_LL(8,6)=1.40D0
  CTURB_LL(8,7)=1.40D0
  CTURB_LL(8,8)=1.40D0
DO J=1,K0_LL
   DO I=1,K0_LL
      CTURB_LL(I,J)=(CTURB_LL(I,J)-1.0D0)/1.5D0+1.0D0
   ENDDO
ENDDO
DO I=KRMIN_LL,KRMAX_LL
   DO J=KRMIN_LL,KRMAX_LL
      CTURBLL(I,J)=1.0D0
   ENDDO
ENDDO
DO I=KRMIN_LL,KRMAX_LL
   X_KERN=RADXXO(I,1)*1.0D4
   IF(X_KERN.LT.RL_LL(1)) X_KERN=RL_LL(1)
   IF(X_KERN.GT.RL_LL(K0_LL)) X_KERN=RL_LL(K0_LL) 
   DO J=KRMIN_LL,KRMAX_LL
      Y_KERN=RADXXO(J,1)*1.0D4
      IF(Y_KERN.LT.RL_LL(1)) Y_KERN=RL_LL(1)
      IF(Y_KERN.GT.RL_LL(K0_LL)) Y_KERN=RL_LL(K0_LL)
      CTURBLL(I,J)=F(X_KERN,Y_KERN,RL_LL,RL_LL,CTURB_LL &
                     ,K0_LL,K0_LL) 
   ENDDO
ENDDO
  RL_GL(1) = RADXXO(1,1)*1.E4 
  RL_GL(2) = 8.0D0
  RL_GL(3) = 10.0D0
  RL_GL(4) = 16.0D0
  RL_GL(5) = 20.0D0
  RL_GL(6) = 30.0D0
  RL_GL(7) = 40.0D0
  RL_GL(8) = 50.0D0
  RL_GL(9) = 60.0D0
  RL_GL(10)= 70.0D0
  RL_GL(11)= 80.0D0
  RL_GL(12)= 90.0D0
  RL_GL(13)=100.0D0
  RL_GL(14)=200.0D0
  RL_GL(15)=300.0D0
  RL_GL(16)=RADXXO(24,1)*1.0D4
! TURBULENCE GRAUPEL BULK RADII IN MKM
  RG_GL(1) = RADXXO(1,6)*1.0D4 
  RG_GL(2) = 30.0D0  
  RG_GL(3) = 60.0D0 
  RG_GL(4) = 100.0D0 
  RG_GL(5) = 200.0D0 
  RG_GL(6) = 300.0D0
  RG_GL(7) = 400.0D0
  RG_GL(8) = 500.0D0
  RG_GL(9) = 600.0D0
  RG_GL(10)= 700.0D0
  RG_GL(11)= 800.0D0
  RG_GL(12)= 900.0D0
  RG_GL(13)=1000.0D0
  RG_GL(14)=2000.0D0
  RG_GL(15)=3000.0D0
  RG_GL(16)=RADXXO(33,6)*1.0D4
DO I=KRMING_GL,KRMAXG_GL
   DO J=KRMINL_GL,KRMAXL_GL
      CTURBGL(I,J)=1.0D0
   ENDDO
ENDDO
DO I=1,K0G_GL
   DO J=1,K0L_GL
      CTURB_GL(I,J)=1.0D0
   ENDDO
ENDDO 
IF(IEPS_400.EQ.1) THEN
    CTURB_GL(1,1)=0.0D0
    CTURB_GL(1,2)=0.0D0
    CTURB_GL(1,3)=1.2D0
    CTURB_GL(1,4)=1.3D0
    CTURB_GL(1,5)=1.4D0
    CTURB_GL(1,6)=1.5D0
    CTURB_GL(1,7)=1.5D0
    CTURB_GL(1,8)=1.5D0
    CTURB_GL(1,9)=1.5D0
    CTURB_GL(1,10)=1.5D0
    CTURB_GL(1,11)=1.5D0
    CTURB_GL(1,12)=1.0D0
    CTURB_GL(1,13)=1.0D0
    CTURB_GL(1,14)=1.0D0
    CTURB_GL(1,15)=1.0D0

    CTURB_GL(2,1)=1.0D0
    CTURB_GL(2,2)=1.4D0
    CTURB_GL(2,3)=1.8D0
    CTURB_GL(2,4)=2.2D0
    CTURB_GL(2,5)=2.6D0
    CTURB_GL(2,6)=3.0D0
    CTURB_GL(2,7)=2.85D0
    CTURB_GL(2,8)=2.7D0
    CTURB_GL(2,9)=2.55D0
    CTURB_GL(2,10)=2.4D0
    CTURB_GL(2,11)=2.25D0
    CTURB_GL(2,12)=1.0D0
    CTURB_GL(2,13)=1.0D0
    CTURB_GL(2,14)=1.0D0

    CTURB_GL(3,1)=7.5D0
    CTURB_GL(3,2)=7.5D0
    CTURB_GL(3,3)=4.5D0
    CTURB_GL(3,4)=4.5D0
    CTURB_GL(3,5)=4.65D0
    CTURB_GL(3,6)=4.65D0
    CTURB_GL(3,7)=4.5D0
    CTURB_GL(3,8)=4.5D0
    CTURB_GL(3,9)=4.0D0
    CTURB_GL(3,10)=3.0D0
    CTURB_GL(3,11)=2.0D0
    CTURB_GL(3,12)=1.5D0
    CTURB_GL(3,13)=1.3D0
    CTURB_GL(3,14)=1.0D0
    
    CTURB_GL(4,1)=5.5D0
    CTURB_GL(4,2)=5.5D0
    CTURB_GL(4,3)=4.5D0
    CTURB_GL(4,4)=4.5D0
    CTURB_GL(4,5)=4.65D0
    CTURB_GL(4,6)=4.65D0
    CTURB_GL(4,7)=4.5D0
    CTURB_GL(4,8)=4.5D0
    CTURB_GL(4,9)=4.0D0
    CTURB_GL(4,10)=3.0D0
    CTURB_GL(4,11)=2.0D0
    CTURB_GL(4,12)=1.5D0
    CTURB_GL(4,13)=1.35D0
    CTURB_GL(4,14)=1.0D0
 
    CTURB_GL(5,1)=4.5D0
    CTURB_GL(5,2)=4.5D0
    CTURB_GL(5,3)=3.3D0
    CTURB_GL(5,4)=3.3D0
    CTURB_GL(5,5)=3.3D0
    CTURB_GL(5,6)=3.4D0
    CTURB_GL(5,7)=3.8D0
    CTURB_GL(5,8)=3.8D0
    CTURB_GL(5,9)=3.8D0
    CTURB_GL(5,10)=3.6D0
    CTURB_GL(5,11)=2.5D0
    CTURB_GL(5,12)=2.0D0
    CTURB_GL(5,13)=1.4D0
    CTURB_GL(5,14)=1.0D0
 
    CTURB_GL(6,1)=4.0D0
    CTURB_GL(6,2)=4.0D0
    CTURB_GL(6,3)=2.8D0
    CTURB_GL(6,4)=2.8D0
    CTURB_GL(6,5)=2.85D0
    CTURB_GL(6,6)=2.9D0
    CTURB_GL(6,7)=3.0D0
    CTURB_GL(6,8)=3.1D0
    CTURB_GL(6,9)=2.9D0
    CTURB_GL(6,10)=2.6D0
    CTURB_GL(6,11)=2.5D0
    CTURB_GL(6,12)=2.0D0
    CTURB_GL(6,13)=1.3D0
    CTURB_GL(6,14)=1.1D0

    CTURB_GL(7,1)=3.5D0
    CTURB_GL(7,2)=3.5D0
    CTURB_GL(7,3)=2.5D0
    CTURB_GL(7,4)=2.5D0
    CTURB_GL(7,5)=2.6D0
    CTURB_GL(7,6)=2.7D0
    CTURB_GL(7,7)=2.8D0
    CTURB_GL(7,8)=2.8D0
    CTURB_GL(7,9)=2.8D0
    CTURB_GL(7,10)=2.6D0
    CTURB_GL(7,11)=2.3D0
    CTURB_GL(7,12)=2.0D0
    CTURB_GL(7,13)=1.3D0
    CTURB_GL(7,14)=1.1D0

    CTURB_GL(8,1)=3.25D0
    CTURB_GL(8,2)=3.25D0
    CTURB_GL(8,3)=2.3D0
    CTURB_GL(8,4)=2.3D0
    CTURB_GL(8,5)=2.35D0
    CTURB_GL(8,6)=2.37D0
    CTURB_GL(8,7)=2.55D0
    CTURB_GL(8,8)=2.55D0
    CTURB_GL(8,9)=2.55D0
    CTURB_GL(8,10)=2.3D0
    CTURB_GL(8,11)=2.1D0
    CTURB_GL(8,12)=1.9D0
    CTURB_GL(8,13)=1.3D0
    CTURB_GL(8,14)=1.1D0

    CTURB_GL(9,1)=3.0D0
    CTURB_GL(9,2)=3.0D0
    CTURB_GL(9,3)=3.1D0
    CTURB_GL(9,4)=2.2D0
    CTURB_GL(9,5)=2.2D0
    CTURB_GL(9,6)=2.2D0
    CTURB_GL(9,7)=2.3D0
    CTURB_GL(9,8)=2.3D0
    CTURB_GL(9,9)=2.5D0
    CTURB_GL(9,10)=2.5D0
    CTURB_GL(9,11)=2.2D0
    CTURB_GL(9,12)=1.8D0
    CTURB_GL(9,13)=1.25D0
    CTURB_GL(9,14)=1.1D0

    CTURB_GL(10,1)=2.75D0
    CTURB_GL(10,2)=2.75D0
    CTURB_GL(10,3)=2.0D0
    CTURB_GL(10,4)=2.0D0
    CTURB_GL(10,5)=2.0D0
    CTURB_GL(10,6)=2.1D0
    CTURB_GL(10,7)=2.2D0
    CTURB_GL(10,8)=2.2D0
    CTURB_GL(10,9)=2.3D0
    CTURB_GL(10,10)=2.3D0
    CTURB_GL(10,11)=2.3D0
    CTURB_GL(10,12)=1.8D0
    CTURB_GL(10,13)=1.2D0
    CTURB_GL(10,14)=1.1D0

    CTURB_GL(11,1)=2.6D0
    CTURB_GL(11,2)=2.6D0
    CTURB_GL(11,3)=1.95D0
    CTURB_GL(11,4)=1.95D0
    CTURB_GL(11,5)=1.95D0
    CTURB_GL(11,6)=2.05D0
    CTURB_GL(11,7)=2.15D0
    CTURB_GL(11,8)=2.15D0
    CTURB_GL(11,9)=2.25D0
    CTURB_GL(11,10)=2.25D0
    CTURB_GL(11,11)=1.9D0
    CTURB_GL(11,12)=1.8D0
    CTURB_GL(11,13)=1.2D0
    CTURB_GL(11,14)=1.1D0

    CTURB_GL(12,1)=2.4D0
    CTURB_GL(12,2)=2.4D0
    CTURB_GL(12,3)=1.85D0
    CTURB_GL(12,4)=1.85D0
    CTURB_GL(12,5)=1.85D0
    CTURB_GL(12,6)=1.75D0
    CTURB_GL(12,7)=1.85D0
    CTURB_GL(12,8)=1.85D0
    CTURB_GL(12,9)=2.1D0
    CTURB_GL(12,10)=2.1D0
    CTURB_GL(12,11)=1.9D0
    CTURB_GL(12,12)=1.8D0 
    CTURB_GL(12,13)=1.3D0
    CTURB_GL(12,14)=1.1D0

    CTURB_GL(13,1)=1.67D0
    CTURB_GL(13,2)=1.67D0
    CTURB_GL(13,3)=1.75D0
    CTURB_GL(13,4)=1.83D0
    CTURB_GL(13,5)=1.87D0
    CTURB_GL(13,6)=2.0D0
    CTURB_GL(13,7)=2.1D0
    CTURB_GL(13,8)=2.12D0
    CTURB_GL(13,9)=2.15D0
    CTURB_GL(13,10)=2.18D0
    CTURB_GL(13,11)=2.19D0
    CTURB_GL(13,12)=1.67D0
    CTURB_GL(13,13)=1.28D0
    CTURB_GL(13,14)=1.0D0

    CTURB_GL(14,1)=1.3D0
    CTURB_GL(14,2)=1.3D0
    CTURB_GL(14,3)=1.35D0
    CTURB_GL(14,4)=1.4D0
    CTURB_GL(14,5)=1.6D0
    CTURB_GL(14,6)=1.7D0
    CTURB_GL(14,7)=1.7D0
    CTURB_GL(14,8)=1.7D0
    CTURB_GL(14,9)=1.7D0
    CTURB_GL(14,10)=1.7D0
    CTURB_GL(14,11)=1.7D0
    CTURB_GL(14,12)=1.4D0
    CTURB_GL(14,13)=1.25D0
    CTURB_GL(14,14)=1.0D0

    CTURB_GL(15,1)=1.17D0
    CTURB_GL(15,2)=1.17D0
    CTURB_GL(15,3)=1.17D0
    CTURB_GL(15,4)=1.25D0
    CTURB_GL(15,5)=1.3D0
    CTURB_GL(15,6)=1.35D0
    CTURB_GL(15,7)=1.4D0
    CTURB_GL(15,8)=1.4D0
    CTURB_GL(15,9)=1.45D0
    CTURB_GL(15,10)=1.47D0
    CTURB_GL(15,11)=1.44D0
    CTURB_GL(15,12)=1.3D0
    CTURB_GL(15,13)=1.12D0
    CTURB_GL(15,14)=1.0D0

    CTURB_GL(16,1)=1.17D0
    CTURB_GL(16,2)=1.17D0
    CTURB_GL(16,3)=1.17D0
    CTURB_GL(16,4)=1.25D0
    CTURB_GL(16,5)=1.3D0
    CTURB_GL(16,6)=1.35D0
    CTURB_GL(16,7)=1.4D0
    CTURB_GL(16,8)=1.45D0
    CTURB_GL(16,9)=1.45D0
    CTURB_GL(16,10)=1.47D0
    CTURB_GL(16,11)=1.44D0
    CTURB_GL(16,12)=1.3D0
    CTURB_GL(16,13)=1.12D0
    CTURB_GL(16,14)=1.0D0
ENDIF
IF(IEPS_800.EQ.1) THEN
    CTURB_GL(1,1) =0.00D0
    CTURB_GL(1,2) =0.00D0
    CTURB_GL(1,3) =1.00D0
    CTURB_GL(1,4) =1.50D0
    CTURB_GL(1,5) =1.40D0
    CTURB_GL(1,6) =1.30D0
    CTURB_GL(1,7) =1.20D0
    CTURB_GL(1,8) =1.10D0
    CTURB_GL(1,9) =1.00D0
    CTURB_GL(1,10)=1.00D0
    CTURB_GL(1,11)=1.00D0
    CTURB_GL(1,12)=1.00D0
    CTURB_GL(1,13)=1.00D0
    CTURB_GL(1,14)=1.00D0
    CTURB_GL(1,15)=1.00D0
    CTURB_GL(1,16)=1.00D0

    CTURB_GL(2,1) =0.00D0
    CTURB_GL(2,2) =0.00D0
    CTURB_GL(2,3) =1.00D0
    CTURB_GL(2,4) =2.00D0
    CTURB_GL(2,5) =1.80D0
    CTURB_GL(2,6) =1.70D0
    CTURB_GL(2,7) =1.60D0
    CTURB_GL(2,8) =1.50D0
    CTURB_GL(2,9) =1.50D0
    CTURB_GL(2,10)=1.50D0
    CTURB_GL(2,11)=1.50D0
    CTURB_GL(2,12)=1.50D0
    CTURB_GL(2,13)=1.50D0
    CTURB_GL(2,14)=1.00D0
    CTURB_GL(2,15)=1.00D0
    CTURB_GL(2,16)=1.00D0

    CTURB_GL(3,1) =0.00D0
    CTURB_GL(3,2) =0.00D0
    CTURB_GL(3,3) =4.00D0
    CTURB_GL(3,4) =7.65D0
    CTURB_GL(3,5) =7.65D0
    CTURB_GL(3,6) =8.00D0
    CTURB_GL(3,7) =8.00D0
    CTURB_GL(3,8) =7.50D0
    CTURB_GL(3,9) =6.50D0
    CTURB_GL(3,10)=6.00D0
    CTURB_GL(3,11)=5.00D0
    CTURB_GL(3,12)=4.50D0
    CTURB_GL(3,13)=4.00D0
    CTURB_GL(3,14)=2.00D0
    CTURB_GL(3,15)=1.30D0
    CTURB_GL(3,16)=1.00D0

    CTURB_GL(4,1) =7.50D0
    CTURB_GL(4,2) =7.50D0
    CTURB_GL(4,3) =7.50D0
    CTURB_GL(4,4) =7.65D0
    CTURB_GL(4,5) =7.65D0
    CTURB_GL(4,6) =8.00D0
    CTURB_GL(4,7) =8.00D0
    CTURB_GL(4,8) =7.50D0
    CTURB_GL(4,9) =6.50D0
    CTURB_GL(4,10)=6.00D0
    CTURB_GL(4,11)=5.00D0
    CTURB_GL(4,12)=4.50D0
    CTURB_GL(4,13)=4.00D0
    CTURB_GL(4,14)=2.00D0
    CTURB_GL(4,15)=1.30D0
    CTURB_GL(4,16)=1.00D0
    
    CTURB_GL(5,1) =5.50D0
    CTURB_GL(5,2) =5.50D0
    CTURB_GL(5,3) =5.50D0
    CTURB_GL(5,4) =5.75D0
    CTURB_GL(5,5) =5.75D0
    CTURB_GL(5,6) =6.00D0
    CTURB_GL(5,7) =6.25D0
    CTURB_GL(5,8) =6.17D0
    CTURB_GL(5,9) =5.75D0
    CTURB_GL(5,10)=5.25D0
    CTURB_GL(5,11)=4.75D0
    CTURB_GL(5,12)=4.25D0
    CTURB_GL(5,13)=4.00D0
    CTURB_GL(5,14)=2.00D0
    CTURB_GL(5,15)=1.35D0
    CTURB_GL(5,16)=1.00D0
 
    CTURB_GL(6,1) =4.50D0
    CTURB_GL(6,2) =4.50D0
    CTURB_GL(6,3) =4.50D0
    CTURB_GL(6,4) =4.75D0
    CTURB_GL(6,5) =4.75D0
    CTURB_GL(6,6) =5.00D0
    CTURB_GL(6,7) =5.25D0
    CTURB_GL(6,8) =5.25D0
    CTURB_GL(6,9) =5.00D0
    CTURB_GL(6,10)=4.75D0
    CTURB_GL(6,11)=4.50D0
    CTURB_GL(6,12)=4.00D0
    CTURB_GL(6,13)=3.75D0
    CTURB_GL(6,14)=2.00D0
    CTURB_GL(6,15)=1.40D0
    CTURB_GL(6,16)=1.00D0
 
    CTURB_GL(7,1) =4.00D0
    CTURB_GL(7,2) =4.00D0
    CTURB_GL(7,3) =4.00D0
    CTURB_GL(7,4) =4.00D0
    CTURB_GL(7,5) =4.00D0
    CTURB_GL(7,6) =4.25D0
    CTURB_GL(7,7) =4.50D0
    CTURB_GL(7,8) =4.67D0
    CTURB_GL(7,9) =4.50D0
    CTURB_GL(7,10)=4.30D0
    CTURB_GL(7,11)=4.10D0
    CTURB_GL(7,12)=3.80D0
    CTURB_GL(7,13)=3.50D0
    CTURB_GL(7,14)=2.00D0
    CTURB_GL(7,15)=1.30D0
    CTURB_GL(7,16)=1.10D0

    CTURB_GL(8,1) =3.50D0
    CTURB_GL(8,2) =3.50D0
    CTURB_GL(8,3) =3.50D0
    CTURB_GL(8,4) =3.65D0
    CTURB_GL(8,5) =3.65D0
    CTURB_GL(8,6) =3.80D0
    CTURB_GL(8,7) =4.1D02
    CTURB_GL(8,8) =4.17D0
    CTURB_GL(8,9) =4.17D0
    CTURB_GL(8,10)=4.00D0
    CTURB_GL(8,11)=3.80D0
    CTURB_GL(8,12)=3.67D0
    CTURB_GL(8,13)=3.40D0
    CTURB_GL(8,14)=2.00D0
    CTURB_GL(8,15)=1.30D0
    CTURB_GL(8,16)=1.10D0

    CTURB_GL(9,1) =3.25D0
    CTURB_GL(9,2) =3.25D0
    CTURB_GL(9,3) =3.25D0
    CTURB_GL(9,4) =3.25D0
    CTURB_GL(9,5) =3.25D0
    CTURB_GL(9,6) =3.50D0
    CTURB_GL(9,7) =3.75D0
    CTURB_GL(9,8) =3.75D0
    CTURB_GL(9,9) =3.75D0
    CTURB_GL(9,10)=3.75D0
    CTURB_GL(9,11)=3.60D0
    CTURB_GL(9,12)=3.40D0
    CTURB_GL(9,13)=3.25D0
    CTURB_GL(9,14)=2.00D0
    CTURB_GL(9,15)=1.30D0
    CTURB_GL(9,16)=1.10D0
    
    CTURB_GL(10,1) =3.00D0
    CTURB_GL(10,2) =3.00D0
    CTURB_GL(10,3) =3.00D0
    CTURB_GL(10,4) =3.10D0
    CTURB_GL(10,5) =3.10D0
    CTURB_GL(10,6) =3.25D0
    CTURB_GL(10,7) =3.40D0
    CTURB_GL(10,8) =3.50D0
    CTURB_GL(10,9) =3.50D0
    CTURB_GL(10,10)=3.50D0
    CTURB_GL(10,11)=3.40D0
    CTURB_GL(10,12)=3.25D0
    CTURB_GL(10,13)=3.15D0
    CTURB_GL(10,14)=1.90D0
    CTURB_GL(10,15)=1.30D0
    CTURB_GL(10,16)=1.10D0

    CTURB_GL(11,1) =2.75D0
    CTURB_GL(11,2) =2.75D0
    CTURB_GL(11,3) =2.75D0
    CTURB_GL(11,4) =2.75D0
    CTURB_GL(11,5) =2.75D0
    CTURB_GL(11,6) =3.00D0
    CTURB_GL(11,7) =3.25D0
    CTURB_GL(11,8) =3.25D0
    CTURB_GL(11,9) =3.25D0
    CTURB_GL(11,10)=3.25D0
    CTURB_GL(11,11)=3.25D0
    CTURB_GL(11,12)=3.15D0
    CTURB_GL(11,13)=3.00D0
    CTURB_GL(11,14)=1.80D0
    CTURB_GL(11,15)=1.30D0
    CTURB_GL(11,16)=1.10D0

    CTURB_GL(12,1) =2.60D0
    CTURB_GL(12,2) =2.60D0
    CTURB_GL(12,3) =2.60D0
    CTURB_GL(12,4) =2.67D0
    CTURB_GL(12,5) =2.67D0
    CTURB_GL(12,6) =2.75D0
    CTURB_GL(12,7) =3.00D0
    CTURB_GL(12,8) =3.17D0
    CTURB_GL(12,9) =3.17D0
    CTURB_GL(12,10)=3.17D0
    CTURB_GL(12,11)=3.10D0
    CTURB_GL(12,12)=2.90D0
    CTURB_GL(12,13)=2.80D0
    CTURB_GL(12,14)=1.87D0
    CTURB_GL(12,15)=1.37D0
    CTURB_GL(12,16)=1.10D0

    CTURB_GL(13,1) =2.40D0
    CTURB_GL(13,2) =2.40D0
    CTURB_GL(13,3) =2.40D0
    CTURB_GL(13,4) =2.50D0
    CTURB_GL(13,5) =2.50D0
    CTURB_GL(13,6) =2.67D0
    CTURB_GL(13,7) =2.83D0
    CTURB_GL(13,8) =2.90D0
    CTURB_GL(13,9) =3.00D0
    CTURB_GL(13,10)=2.90D0
    CTURB_GL(13,11)=2.85D0
    CTURB_GL(13,12)=2.80D0
    CTURB_GL(13,13)=2.75D0
    CTURB_GL(13,14)=1.83D0
    CTURB_GL(13,15)=1.30D0
    CTURB_GL(13,16)=1.10D0

    CTURB_GL(14,1) =1.67D0
    CTURB_GL(14,2) =1.67D0
    CTURB_GL(14,3) =1.67D0
    CTURB_GL(14,4) =1.75D0
    CTURB_GL(14,5) =1.75D0
    CTURB_GL(14,6) =1.83D0
    CTURB_GL(14,7) =1.87D0
    CTURB_GL(14,8) =2.00D0
    CTURB_GL(14,9) =2.10D0
    CTURB_GL(14,10)=2.12D0
    CTURB_GL(14,11)=2.15D0
    CTURB_GL(14,12)=2.18D0
    CTURB_GL(14,13)=2.19D0
    CTURB_GL(14,14)=1.67D0
    CTURB_GL(14,15)=1.28D0
    CTURB_GL(14,16)=1.00D0

    CTURB_GL(15,1) =1.30D0
    CTURB_GL(15,2) =1.30D0
    CTURB_GL(15,3) =1.30D0
    CTURB_GL(15,4) =1.35D0
    CTURB_GL(15,5) =1.35D0
    CTURB_GL(15,6) =1.40D0
    CTURB_GL(15,7) =1.60D0
    CTURB_GL(15,8) =1.70D0
    CTURB_GL(15,9) =1.70D0
    CTURB_GL(15,10)=1.70D0
    CTURB_GL(15,11)=1.70D0
    CTURB_GL(15,12)=1.70D0
    CTURB_GL(15,13)=1.70D0
    CTURB_GL(15,14)=1.40D0
    CTURB_GL(15,15)=1.25D0
    CTURB_GL(15,16)=1.00D0

    CTURB_GL(16,1) =1.17D0
    CTURB_GL(16,2) =1.17D0
    CTURB_GL(16,3) =1.17D0
    CTURB_GL(16,4) =1.17D0
    CTURB_GL(16,5) =1.17D0
    CTURB_GL(16,6) =1.25D0
    CTURB_GL(16,7) =1.30D0
    CTURB_GL(16,8) =1.35D0
    CTURB_GL(16,9) =1.40D0
    CTURB_GL(16,10)=1.45D0
    CTURB_GL(16,11)=1.45D0
    CTURB_GL(16,12)=1.47D0
    CTURB_GL(16,13)=1.44D0
    CTURB_GL(16,14)=1.30D0
    CTURB_GL(16,15)=1.12D0
    CTURB_GL(16,16)=1.00D0
ENDIF
IF(IEPS_800.EQ.1.AND.IEPS_1600.EQ.1) THEN
   DO I=1,K0G_GL
      DO J=1,K0L_GL
         CTURB_GL(I,J)=CTURB_GL(I,J)*1.7D0
      ENDDO
    ENDDO 
ENDIF
DO J=1,K0L_GL
   DO I=1,K0G_GL
      CTURB_GL(I,J)=(CTURB_GL(I,J)-1.0D0)/1.5D0+1.0D0
   ENDDO
ENDDO
DO I=KRMING_GL,KRMAXG_GL
   DO J=KRMINL_GL,KRMAXL_GL
      CTURBGL(I,J)=1.
   ENDDO
ENDDO
DO I=KRMING_GL,KRMAXG_GL 
   X_KERN=RADXXO(I,6)*1.0D4
   IF(X_KERN.LT.RG_GL(1)) X_KERN=RG_GL(1)
   IF(X_KERN.GT.RG_GL(K0G_GL)) X_KERN=RG_GL(K0G_GL) 
   DO J=KRMINL_GL,KRMAXL_GL
      Y_KERN=RADXXO(J,1)*1.0D4
      IF(Y_KERN.LT.RL_GL(1)) Y_KERN=RL_GL(1)
      IF(Y_KERN.GT.RL_GL(K0L_GL)) Y_KERN=RL_GL(K0L_GL)
      CTURBGL(I,J)=F(X_KERN,Y_KERN,RG_GL,RL_GL,CTURB_GL &
                     ,K0G_GL,K0L_GL)
   ENDDO
ENDDO
IF(IEPS_800.EQ.1) THEN
   DO I=KRMING_GL,15
      DO J=KRMINL_GL,13
         IF(CTURBGL(I,J).LT.3.0D0) CTURBGL(I,J)=3.0D0
      ENDDO
   ENDDO
ENDIF
IF(IEPS_1600.EQ.1) THEN
   DO I=KRMING_GL,15
      DO J=KRMINL_GL,13
         IF(CTURBGL(I,J).LT.5.1D0) CTURBGL(I,J)=5.1D0
      ENDDO
   ENDDO
ENDIF
DO I=1,33
   DO J=1,24
      IF(I.LE.14.AND.J.EQ.8) CTURBGL(I,J)=1.0D0
      IF(I.GT.14.AND.J.LE.8) CTURBGL(I,J)=1.2D0
   ENDDO
ENDDO     
RETURN
END SUBROUTINE TURBCOEF
!===================================================================
  real function f(x,y,x0,y0,table,k0,kk0)
! two-dimensional linear interpolation of the collision efficiency
! with help table(k0,kk0)

 implicit none
 integer k0,kk0,k,ir,kk,iq
 double precision x,y,p,q,ec,ek
 double precision x0(k0),y0(kk0),table(k0,kk0)

do k=2,k0
   if(x.le.x0(k).and.x.ge.x0(k-1)) then
      ir=k     
   elseif(x.gt.x0(k0)) then
      ir=k0+1
   elseif(x.lt.x0(1)) then
      ir=1
   endif
enddo
do kk=2,kk0
   if(y.le.y0(kk).and.y.ge.y0(kk-1)) iq=kk
enddo
if(ir.lt.k0+1) then
   if(ir.ge.2) then
      p =(x-x0(ir-1))/(x0(ir)-x0(ir-1))
      q =(y-y0(iq-1))/(y0(iq)-y0(iq-1))
      ec=(1.d0-p)*(1.d0-q)*table(ir-1,iq-1)+ &
             p*(1.d0-q)*table(ir,iq-1)+ &
             q*(1.d0-p)*table(ir-1,iq)+ &
                  p*q*table(ir,iq)    
    else
      q =(y-y0(iq-1))/(y0(iq)-y0(iq-1))
      ec=(1.d0-q)*table(1,iq-1)+q*table(1,iq)    
    endif
  else
    q =(y-y0(iq-1))/(y0(iq)-y0(iq-1))
    ek=(1.d0-q)*table(k0,iq-1)+q*table(k0,iq)
    ec=min(ek,1.d0) 
  endif
  f=ec
  return
  end function f
!======================================================================
SUBROUTINE COAL_BOTT_NEW(FF1R,FF2R,FF3R, &
   FF4R,FF5R,TT,QQ,PP,RHO,dthalf,TCRIT,TTCOAL, &
   cld2raint,rimecldt,aggregatet,rimecldsnowt,rimecldaggrt, &
   rimecldgraut,rimecldhailt,rain2icet,rain2snowt,rain2agt,&
   rain2grt,rain2hat)

use micro_prm, only:ipris

IMPLICIT NONE
INTEGER KR,ICE
INTEGER icol_drop,icol_snow,icol_graupel,icol_hail, &
 icol_column,icol_plate,icol_dendrite,icol_drop_brk
double precision  g1(nkr),g2(nkr,icemax),g3(nkr),g4(nkr),g5(nkr)
double precision gdumb(nkr),xl_dumb(0:nkr)
double precision g2_1(nkr),g2_2(nkr),g2_3(nkr)
real cont_fin_drop,dconc,conc_icempl,deldrop,t_new, &
 cont_fin_ice,conc_old,conc_new,cont_init_ice, &
 cont_init_drop,ALWC
REAL    FF1R(NKR),FF2R(NKR,ICEMAX),FF3R(NKR),FF4R(NKR),FF5R(NKR)
REAL DTHALF
REAL TCRIT,TTCOAL
INTEGER I,J,IT,NDIV
REAL RHO
REAL :: cld2raint,rimecldt,aggregatet,rimecldsnowt,rimecldaggrt, &
        rimecldgraut,rimecldhailt,rain2icet,rain2snowt,rain2agt, &
        rain2grt,rain2hat,totalcloud,totalrain,totalprissnow
DOUBLE PRECISION break_drop_bef,break_drop_aft,dtbreakup
DOUBLE PRECISION break_drop_per,mytotal
DOUBLE PRECISION TT,QQ,PP,prdkrn,prdkrn1
parameter (prdkrn1=1.d0)

icol_drop_brk=0
icol_drop=0
icol_snow=0
icol_graupel=0
icol_hail=0
icol_column=0
icol_plate=0
icol_dendrite=0
t_new=tt
CALL MISC1(PP,cwll_1000mb,cwll_750mb,cwll_500mb,cwll,nkr)
! THIS IS FOR BREAKUP
DO I=1,NKR
   DO J=1,NKR
      CWLL(I,J)=ECOALMASSM(I,J)*CWLL(I,J)
   ENDDO
ENDDO

! THIS IS FOR TURBULENCE
IF (LIQTURB.EQ.1)THEN
   DO I=1,KRMAX_LL
      DO J=1,KRMAX_LL
         CWLL(I,J)=CTURBLL(I,J)*CWLL(I,J)
      END DO
   END DO
END IF
CALL MODKRN(TT,QQ,PP,PRDKRN,TTCOAL)
DO KR=1,NKR
   G1(KR)=FF1R(KR)*3.*XL(KR)*XL(KR)*1.E3
   G2(KR,1)=FF2R(KR,1)*3*xi(KR,1)*XI(KR,1)*1.e3
   G2(KR,2)=FF2R(KR,2)*3.*xi(KR,2)*XI(KR,2)*1.e3
   G2(KR,3)=FF2R(KR,3)*3.*xi(KR,3)*XI(KR,3)*1.e3
   G3(KR)=FF3R(KR)*3.*xs(kr)*xs(kr)*1.e3
   G4(KR)=FF4R(KR)*3.*xg(kr)*xg(kr)*1.e3
   G5(KR)=FF5R(KR)*3.*xh(kr)*xh(kr)*1.e3
   g2_1(kr)=g2(KR,1)
   g2_2(KR)=g2(KR,2)
   g2_3(KR)=g2(KR,3)
   if(kr.gt.(nkr-jbreak).and.g1(kr).gt.1.e-17)icol_drop_brk=1
   if(g1(kr).gt.1.e-10)icol_drop=1
   if(g2_1(kr).gt.1.e-10)icol_column=1
   if(g2_2(kr).gt.1.e-10)icol_plate=1
   if(g2_3(kr).gt.1.e-10)icol_dendrite=1
   if(g3(kr).gt.1.e-10)icol_snow=1
   if(g4(kr).gt.1.e-10)icol_graupel=1
   if(g5(kr).gt.1.e-10)icol_hail=1
ENDDO 
!Adele
mytotal=sum(g1)
! calculation of initial hydromteors content in g/cm**3 :
cont_init_drop=0.
cont_init_ice=0.
do kr=1,nkr
   cont_init_drop=cont_init_drop+g1(kr)
   cont_init_ice=cont_init_ice+g3(kr)+g4(kr)+g5(kr)
   do ice=1,icemax
      cont_init_ice=cont_init_ice+g2(kr,ice)
   enddo
enddo
cont_init_drop=col*cont_init_drop*1.e-3
cont_init_ice=col*cont_init_ice*1.e-3
! calculation of alwc in g/m**3
alwc=cont_init_drop*1.e6
! calculation of initial hydromteors content in g/cm**3 :

!----------LIQUID-LIQUID INTERACTIONS-----------------------
totalcloud = sum(g1(1:krdrop-1))*col*1e-3/rho
if (icol_drop.eq.1)then 
  if (imbudget >= 1) then
     cld2raint = cld2raint - totalcloud 
  endif
  call coll_xxx (G1,CWLL,XL_MG,CHUCM,IMA,NKR)
!print*,'after cc',sum(g1)-mytotal
  if (imbudget >= 1) then
     totalcloud = sum(g1(1:krdrop-1))*col*1e-3/rho
     cld2raint = cld2raint + totalcloud 
  endif
!-----------RAINDROP BREAKUP------------------------------
!Adele - turned breakup off
  if(icol_drop_brk.eq.10)then
     ndiv=1
10 continue
     do it = 1,ndiv
       if (ndiv.gt.10000)stop 'ndiv'
       dtbreakup = dthalf/ndiv
       if (it.eq.1)then
          do kr=1,nkr
             gdumb(kr)= g1(kr)*1.D-3
             xl_dumb(kr)=xl_mg(KR)*1.D-3
          end do
          break_drop_bef=0.d0
          do kr=1,nkr
             break_drop_bef=break_drop_bef+g1(kr)*1.D-3
          enddo
       end if
       call breakup(gdumb,xl_dumb,dtbreakup,brkweight, &
            pkij,qkj,nkr,jbreak)
     end do
     break_drop_aft=0.0d0
     do kr=1,nkr
        break_drop_aft=break_drop_aft+gdumb(kr)
     enddo
     break_drop_per=1.-break_drop_aft/break_drop_bef
     if (break_drop_per.gt.1.e-6)then
         ndiv=ndiv*2
print*,break_drop_per,mytotal-sum(g1),mytotal,sum(g1),break_drop_aft,break_drop_bef
         GO TO 10
     else
        do kr=1,nkr
           g1(kr)=gdumb(kr)*1.D3
        end do
     end if
  end if
!print*,'after breakup',sum(g1)-mytotal
  totalrain = sum(g1(krdrop:nkr))*col*1e-3/rho
!-----------------SNOW-LIQUID COLLISIONS---------------------
 if (iceprocs .eq. 1) then
  if (imbudget >= 1) then
     rimecldt = rimecldt - totalcloud
     rain2icet = rain2icet - totalrain
  endif
  if (icol_snow.eq.1)then 
     if (imbudget >= 2) then
        rimecldaggrt = rimecldaggrt - totalcloud 
        rain2agt = rain2agt - totalrain
     endif
     if(tt.lt.tcrit) then
        call coll_xyz (g1,g3,g4,cwls,xl_mg,xs_mg, &
                       chucm,ima,prdkrn1,nkr,0)
     endif
     if(tt.ge.tcrit) then
        call coll_xyz (g1,g3,g5,cwls,xl_mg,xs_mg, &
                       chucm,ima,prdkrn1,nkr,0)
     endif
     if(alwc.lt.alcr) then
        call coll_xyx (g3,g1,cwsl,xs_mg,xl_mg, &
                       chucm,ima,prdkrn1,nkr,1)
     endif
     if(alwc.ge.alcr) then
        call coll_xyz (g3,g1,g4,cwsl,xs_mg,xl_mg, &
                       chucm,ima,prdkrn1,nkr,1)
     endif
     if (imbudget >= 2) then
        totalcloud = sum(g1(1:krdrop-1))*col*1e-3/rho
        rimecldaggrt = rimecldaggrt + totalcloud 
        totalrain = sum(g1(krdrop:nkr))*col*1e-3/rho
        rain2agt = rain2agt + totalrain
     endif
  end if
!-----------------GRAUPEL-LIQUID COLLISIONS----------------
! water - graupel = graupel (t < tcrit ; xl_mg ge xg_mg)
! graupel - water = graupel (t < tcrit ; xg_mg > xl_mg)
! water - graupel = hail (t ge tcrit ; xl_mg ge xg_mg)
! graupel - water = hail (t ge tcrit ; xg_mg > xl_mg)
  if (icol_graupel.eq.1)then 
     if (imbudget >= 2) then
        rimecldgraut = rimecldgraut - totalcloud 
        rain2grt = rain2grt - totalrain
     endif
     if (alwc.lt.alcr_hail) then
        call coll_xyy (g1,g4,cwlg,xl_mg,xg_mg, &
                       chucm,ima,prdkrn1,nkr,0)
! if cwc .gt. alcr_hail, water+graupel - hail for nr .gt. kp_hail
     else
        call coll_xyyz(g1,g4,g5,cwlg,xl_mg,xg_mg,chucm, &
                       ima,prdkrn1,nkr,0,kp_hail)
     endif
 ! ice-multiplication preparation:
     conc_old=0.
     conc_new=0.
     do kr=kr_icempl,nkr
        conc_old=conc_old+col*g1(kr)/xl_mg(kr)
     enddo
     if (alwc.lt.alcr_hail) then
        call coll_xyx (g4,g1,cwgl,xg_mg,xl_mg, &
                       chucm,ima,prdkrn1,nkr,1)
     else
        call coll_xyxz(g4,g1,g5,cwlg,xg_mg,xl_mg,chucm, &
                       ima,prdkrn1,nkr,1,kp_hail)
     endif
     if (imbudget >= 2) then
         totalcloud = sum(g1(1:krdrop-1))*col*1e-3/rho
         rimecldgraut = rimecldgraut + totalcloud 
         totalrain = sum(g1(krdrop:nkr))*col*1e-3/rho
         rain2grt = rain2grt + totalrain
     endif
!----ICE-MULTIPLICATION : Hallet-Mossop processes (1 per 250 collisions)
     if (icempl.eq.1) then
        if (tt.ge.265.15.and.tt.le.270.15) then
           do kr=kr_icempl,nkr
              conc_new=conc_new+col*g1(kr)/xl_mg(kr)
           enddo
           dconc=conc_old-conc_new
           if(tt.le.268.15) then
              conc_icempl=dconc*4.e-3*(265.15-tt)/(265.15-268.15)
           elseif(tt.gt.268.15) then
              conc_icempl=dconc*4.e-3*(tcrit-tt)/(tcrit-268.15)
           endif
           g2_2(1)=g2_2(1)+conc_icempl*xi2_mg(1)/col
         endif
     endif
  endif
!---------------HAIL-LIQUID COLLISIONS-----------------------------
! water - hail = hail (xl_mg ge xh_mg)    (kxyy=2)
! hail - water = hail (xh_mg > xl_mg)     (kxyx=3)
  if(icol_hail.eq.1) then
    if (imbudget >= 2) then
       rimecldhailt = rimecldhailt - totalcloud 
       rain2hat = rain2hat - totalrain
    endif
    call coll_xyy (g1,g5,cwlh,xl_mg,xh_mg, &
                       chucm,ima,prdkrn1,nkr,0)
    call coll_xyx (g5,g1,cwhl,xh_mg,xl_mg, &
                       chucm,ima,prdkrn1,nkr,1)
    if (imbudget >= 2) then
       totalcloud = sum(g1(1:krdrop-1))*col*1e-3/rho
       rimecldhailt = rimecldhailt + totalcloud 
       totalrain = sum(g1(krdrop:nkr))*col*1e-3/rho
       rain2hat = rain2hat + totalrain
    endif
  endif
!--------------COLUMN-LIQUID COLLISIONS--------------------------
! interactions between water and crystals :
! interactions between water and columns :
! water - columns = graupel (t < tcrit ; xl_mg ge xi_mg)    (kxyz=6)
! water - columns = hail (t ge tcrit ; xl_mg ge xi_mg)      (kxyz=7)
! columns - water = columns/graupel (xi_mg > xl_mg)             (kxyx=4); kxyxz=2)
! now: columns - water = columns (xi_mg > xl_mg)             (kxyx=4); kxyxz=2)
  if (imbudget >= 2 .and. ipris > 0) then
     rimecldsnowt = rimecldsnowt - totalcloud 
     rain2snowt = rain2snowt - totalrain
  endif
  if(icol_column.eq.1) then
     if(tt.lt.tcrit) then
        call coll_xyz (g1,g2_1,g4,cwli_1,xl_mg,xi1_mg, &
                         chucm,ima,prdkrn,nkr,0)
     endif
     if(tt.ge.tcrit) then
        call coll_xyz (g1,g2_1,g5,cwli_1,xl_mg,xi1_mg, &
                        chucm,ima,prdkrn,nkr,0)
     endif
     call coll_xyx (g2_1,g1,cwil_1,xi1_mg,xl_mg, &
                         chucm,ima,prdkrn,nkr,1)
  endif
!------------PLATE-LIQUID COLLISIONS---------------------------------
! water - plates = graupel (t < tcrit ; xl_mg ge xi2_mg)    (kxyz=8)
! water - plates = hail (t ge tcrit ; xl_mg ge xi2_mg)      (kxyz=9)
! plates - water = plates/graupel (xi2_mg > xl_mg)              (kxyx=5; kxyxz=3)
!now: plates - water = plates (xi2_mg > xl_mg)              (kxyx=5; kxyxz=3)
  if(icol_plate.eq.1) then
    if(tt.lt.tcrit) then
      call coll_xyz (g1,g2_2,g4,cwli_2,xl_mg,xi2_mg, &
                         chucm,ima,prdkrn,nkr,0)
    endif
    if(tt.ge.tcrit) then
      call coll_xyz (g1,g2_2,g5,cwli_2,xl_mg,xi2_mg, &
                         chucm,ima,prdkrn,nkr,0)
    endif
    call coll_xyx (g2_2,g1,cwil_2,xi2_mg,xl_mg, &
                         chucm,ima,prdkrn,nkr,1)
  endif
!--------------DENDRITE-LIQUID COLLISIONS------------------------------
! water - dendrites = graupel (t < tcrit ; xl_mg ge xi3_mg) (kxyz=10)
! water - dendrites = hail (t ge tcrit ; xl_mg ge xi3_mg)   (kxyz=11)
! dendrites - water = dendrites/graupel (xi3_mg > xl_mg)         (kxyx=6; kxyxz=4)
!now dendrites - water = dendrites (xi3_mg > xl_mg)         (kxyx=6; kxyxz=4)
  if(icol_dendrite.eq.1) then
    if(tt.lt.tcrit) then
      call coll_xyz (g1,g2_3,g4,cwli_3,xl_mg,xi3_mg, &
                         chucm,ima,prdkrn,nkr,0)
    endif
    if(tt.ge.tcrit) then
      call coll_xyz (g1,g2_3,g5,cwli_3,xl_mg,xi3_mg, &
                         chucm,ima,prdkrn,nkr,0)
    endif
    call coll_xyx (g2_3,g1,cwil_3,xi3_mg,xl_mg, &
                         chucm,ima,prdkrn,nkr,1)
  endif
  if (imbudget >= 1) then
     totalcloud = sum(g1(1:krdrop-1))*col*1e-3/rho
     totalrain = sum(g1(krdrop:nkr))*col*1e-3/rho
     rimecldt = rimecldt + totalcloud 
     rain2icet = rain2icet + totalrain
  endif
  if (imbudget >= 2 .and. ipris > 0) then
     rimecldsnowt = rimecldsnowt + totalcloud 
     rain2snowt = rain2snowt + totalrain
  endif
 endif !endif ice present
endif !endif liquid collisions

!------------CRYSTAL-CRYSTAL INTERACTIONS--------------------------
! if(t.le.TTCOAL) - no interactions between crystals
if(tt.gt.TTCOAL .and. iceprocs == 1) then
  if (imbudget >= 1) then
     totalprissnow = sum(g2)*col*1e-3/rho
     aggregatet = rimecldt - totalcloud 
  endif
!------------COLUMN-ICE COLLISIONS--------------------------------
if(icol_column.eq.1) then
  ! columns - columns = snow
  call coll_xxy (g2_1,g3,cwii_1_1,xi1_mg, &
     chucm,ima,prdkrn,nkr)
  ! columns - plates = snow (xi1_mg ge xi2_mg)                (kxyz=12)
  ! plates - columns = snow (xi2_mg > xi1_mg)                 (kxyz=13)
  if(icol_plate.eq.1) then     
   call coll_xyz (g2_1,g2_2,g3,cwii_1_2,xi1_mg,xi2_mg, &
     chucm,ima,prdkrn,nkr,0)
   call coll_xyz (g2_2,g2_1,g3,cwii_2_1,xi2_mg,xi1_mg, &
     chucm,ima,prdkrn,nkr,1)
  end if
  ! columns - dendrites = snow (xi1_mg ge xi3_mg)             (kxyz=14)
  ! dendrites - columns = snow (xi3_mg > xi1_mg)              (kxyz=15)
  if(icol_dendrite.eq.1) then
     call coll_xyz (g2_1,g2_3,g3,cwii_1_3,xi1_mg,xi3_mg, &
     chucm,ima,prdkrn,nkr,0)
     call coll_xyz (g2_3,g2_1,g3,cwii_3_1,xi3_mg,xi1_mg, &
     chucm,ima,prdkrn,nkr,1)
  end if
  ! columns - snow = snow (xi1_mg ge xs_mg)                   (kxyy=3)
  ! snow - columns = snow (xs_mg > xi1_mg)                    (kxyx=7)
  if(icol_snow.eq.1) then
!       call coll_xyy (g2_1,g3,cwis_1,xi1_mg,xs_mg, &
!                      chucm,ima,prdkrn,nkr,0)
   call coll_xyx (g3,g2_1,cwsi_1,xs_mg,xi1_mg, &
     chucm,ima,prdkrn,nkr,1)
  endif          
endif
!-----------------PLATES-ICE COLLISIONS-------------------------
 ! plates - plates = snow
if(icol_plate.eq.1) then
   call coll_xxy (g2_2,g3,cwii_2_2,xi2_mg, &
     chucm,ima,prdkrn,nkr)
  ! plates - dendrites = snow (xi2_mg ge xi3_mg)              (kxyz=17)
  ! dendrites - plates = snow (xi3_mg > xi2_mg)               (kxyz=18)
  if(icol_dendrite.eq.1) then
   call coll_xyz (g2_2,g2_3,g3,cwii_2_3,xi2_mg,xi3_mg, &
     chucm,ima,prdkrn,nkr,0)
   call coll_xyz (g2_3,g2_2,g3,cwii_3_2,xi3_mg,xi2_mg, &
     chucm,ima,prdkrn,nkr,1)
  end if
 ! plates - snow = snow (xi2_mg ge xs_mg)                    (kxyy=4)
 ! snow - plates = snow (xs_mg > xi2_mg)                     (kxyx=12)
  if(icol_snow.eq.1) then
!       call coll_xyy (g2_2,g3,cwis_2,xi2_mg,xs_mg, &
!                      chucm,ima,prdkrn,nkr,0)
   call coll_xyx (g3,g2_2,cwsi_2,xs_mg,xi2_mg, &
     chucm,ima,prdkrn,nkr,1)
   end if
endif
!------------DENDRITES-ICE COLLISIONS---------------------------
 ! dendrites - dendrites = snow
if(icol_dendrite.eq.1) then
  call coll_xxy (g2_3,g3,cwii_3_3,xi3_mg, &
                 chucm,ima,prdkrn,nkr)
  ! dendrites - snow = snow (xi3_mg ge xs_mg)                 (kxyy=5)
  ! snow - dendrites = snow (xs_mg > xi3_mg)                  (kxyx=17)
  if(icol_snow.eq.1) then
!       call coll_xyy (g2_3,g3,cwis_3,xi3_mg,xs_mg, &
!                    chucm,ima,prdkrn,nkr,0)
   call coll_xyx (g3,g2_3,cwsi_3,xs_mg,xi3_mg, &
     chucm,ima,prdkrn,nkr,1)
  end if
endif
!-------------SNOW-SNOW COLLISIONS-------------------------------
  if(icol_snow.ne.0) then
   ! snow - snow = snow
   call coll_xxx_prd (g3,cwss,xs_mg,chucm,ima,prdkrn,nkr)
   ! snow - graupel = snow (xs_mg > xg_mg)                     (kxyx=22)
   ! graupel - snow = graupel (xg_mg ge xs_mg)                 (kxyx=23)

   !JF4 no snow-graupel interactions
   !         if(icol_graupel.eq.1) then
   !          call coll_xyx (g3,g4,cwsg,xs_mg,xg_mg, &
   !     &                chucm,ima,prdkrn,nkr,1)
   !         endif
   !end JF4
  endif
  if (imbudget >= 1) then
     totalprissnow = sum(g2)*col*1e-3/rho
     aggregatet = rimecldt + totalcloud 
  endif
endif

!-------------END OF COLLISION CALCULATIONS-------------------

! calculation of finish hydrometeors contents in g/cm**3 :
cont_fin_drop=0.
cont_fin_ice=0.
do kr=1,nkr
   g2(kr,1)=g2_1(kr)
   g2(kr,2)=g2_2(kr)
   g2(kr,3)=g2_3(kr)
   cont_fin_drop=cont_fin_drop+g1(kr)
   cont_fin_ice=cont_fin_ice+g3(kr)+g4(kr)+g5(kr)
   do ice=1,icemax
      cont_fin_ice=cont_fin_ice+g2(kr,ice)
   enddo
enddo
cont_fin_drop=col*cont_fin_drop*1.e-3
cont_fin_ice=col*cont_fin_ice*1.e-3
deldrop=cont_init_drop-cont_fin_drop
! deldrop in g/cm**3
! resulted value of temperature (rob in g/cm**3) :
if(t_new.le.273.15) then
  if(deldrop.ge.0.) then
    t_new=t_new+320.*deldrop/rho
  else
! if deldrop < 0
    if(abs(deldrop).gt.cont_init_drop*0.05) then
      print *, '*** deldrop < 0 in coal_bott ***'
      print *,'*** coal_bott ***'
      print*,'cont_fin_drop = ',cont_fin_drop,cont_fin_ice
      print*,'cont_init_drop = ',cont_init_drop,cont_init_ice
      stop
    endif
  endif
endif

! recalculation of density function f1,f2,f3,f4,f5 in 1/(g*cm**3) :  
DO KR=1,NKR
   FF1R(KR)=G1(KR)/(3.*XL(KR)*XL(KR)*1.E3)
   FF2R(KR,1)=G2(KR,1)/(3*xi(KR,1)*XI(KR,1)*1.e3)
   FF2R(KR,2)=G2(KR,2)/(3.*xi(KR,2)*XI(KR,2)*1.e3)
   FF2R(KR,3)=G2(KR,3)/(3.*xi(KR,3)*XI(KR,3)*1.e3)
   FF3R(KR)=G3(KR)/(3.*xs(kr)*xs(kr)*1.e3)
   FF4R(KR)=G4(KR)/(3.*xg(kr)*xg(kr)*1.e3)
   FF5R(KR)=G5(KR)/(3.*xh(kr)*xh(kr)*1.e3)
ENDDO

tt=t_new
RETURN
END SUBROUTINE COAL_BOTT_NEW
!--------------------------------------------------------------
SUBROUTINE MISC1(PP,cwll_1000mb,cwll_750mb,cwll_500mb,cwll,nkr)
IMPLICIT NONE
INTEGER kr1,kr2,NKR
DOUBLE PRECISION PP
REAL P_Z
double precision cwll(nkr,nkr),cwll_1,cwll_2,cwll_3 &
     ,cwll_1000mb(nkr,nkr),cwll_750mb(nkr,nkr),cwll_500mb(nkr,nkr)

P_Z=PP
do kr1=1,nkr
   do kr2=1,nkr
      cwll_1=cwll_1000mb(kr1,kr2)
      cwll_2=cwll_750mb(kr1,kr2)
      cwll_3=cwll_500mb(kr1,kr2)
      if(p_z.ge.p1) cwll(kr1,kr2)=cwll_1
      if(p_z.eq.p2) cwll(kr1,kr2)=cwll_2
      if(p_z.eq.p3) cwll(kr1,kr2)=cwll_3
      if(p_z.lt.p1.and.p_z.gt.p2) &
         cwll(kr1,kr2)=cwll_2+ &
         (cwll_1-cwll_2)*(p_z-p2)/(p1-p2) 
      if(p_z.lt.p2.and.p_z.gt.p3) &
         cwll(kr1,kr2)=cwll_3+ &
         (cwll_2-cwll_3)*(p_z-p3)/(p2-p3)
      if(p_z.lt.p3) cwll(kr1,kr2)=cwll_3
   enddo
enddo
RETURN
END SUBROUTINE  MISC1
!----------------------------------------------------------
subroutine coll_xxx (g,ckxx,x,chucm,ima,nkr)
implicit double precision (a-h,o-z)
dimension g(nkr),ckxx(nkr,nkr),x(0:nkr)
dimension chucm(nkr,nkr)
double precision ima(nkr,nkr)
gmin=1.d-60
! lower and upper integration limit ix0,ix1
do i=1,nkr-1
   ix0=i
   if(g(i).gt.gmin) exit
enddo

if(ix0.eq.nkr-1) return

do i=nkr-1,1,-1
   ix1=i
   if(g(i).gt.gmin) exit
enddo

do i=ix0,ix1
   do j=i,ix1
      k=ima(i,j)
      kp=k+1
      x0=ckxx(i,j)*g(i)*g(j)
      x0=min(x0,g(i)*x(j))
      if(j.ne.k) then
         x0=min(x0,g(j)*x(i))
      endif
      gsi=x0/x(j)
      gsj=x0/x(i)
      gsk=gsi+gsj
      g(i)=g(i)-gsi
      if(g(i).lt.0.d0) g(i)=0.d0
      g(j)=g(j)-gsj
      gk=g(k)+gsk
      if(g(j).lt.0.d0.and.gk.lt.gmin) then
         g(j)=0.d0
         g(k)=g(k)+gsi
      endif
      flux=0.d0

      if(gk.gt.gmin) then
         x1=dlog(g(kp)/gk+1.d-15)
         if (x1.eq.0)then
            flux=0  
         else
            flux=gsk/x1*(dexp(0.5*x1)-dexp(x1*(0.5-chucm(i,j))))
            flux=min(flux,gsk)
         end if

         g(k)=gk-flux
         if(gk .lt. flux) flux=gk
         if(g(k).lt.0.d0) g(k)=0.d0
         g(kp)=g(kp)+flux
      endif
   end do
end do
end subroutine coll_xxx
!-------------------------------------------------------------
subroutine coll_xxx_prd (g,ckxx,x,chucm,ima,prdkrn,nkr)
implicit double precision (a-h,o-z)
dimension g(nkr),ckxx(nkr,nkr),x(0:nkr)
dimension chucm(nkr,nkr)
double precision ima(nkr,nkr)

gmin=1.d-60
! lower and upper integration limit ix0,ix1
do i=1,nkr-1
   ix0=i
   if(g(i).gt.gmin) exit
enddo

if(ix0.eq.nkr-1) return

do i=nkr-1,1,-1
   ix1=i
   if(g(i).gt.gmin) exit
enddo

do i=ix0,ix1
   do j=i,ix1
      k=ima(i,j)
      kp=k+1
      x0=ckxx(i,j)*g(i)*g(j)*prdkrn
      x0=min(x0,g(i)*x(j))
      if(j.ne.k) then
        x0=min(x0,g(j)*x(i))
      endif
      gsi=x0/x(j)
      gsj=x0/x(i)
      gsk=gsi+gsj
      g(i)=g(i)-gsi
      if(g(i).lt.0.d0) g(i)=0.d0
      g(j)=g(j)-gsj
      gk=g(k)+gsk
      if(g(j).lt.0.d0.and.gk.lt.gmin) then
         g(j)=0.d0
         g(k)=g(k)+gsi
      endif
      flux=0.d0

      if(gk.gt.gmin) then
         x1=dlog(g(kp)/gk+1.d-15)
         if (x1.eq.0)then
            flux=0  
         else
            flux=gsk/x1*(dexp(0.5*x1)-dexp(x1*(0.5-chucm(i,j))))
            flux=min(flux,gsk)
         end if
         g(k)=gk-flux
         if(gk .lt. flux) flux=gk

         if(g(k).lt.0.d0) g(k)=0.d0
          g(kp)=g(kp)+flux
      endif
   enddo
enddo
return
end subroutine coll_xxx_prd 
!-------------------------------------------------------------------
subroutine modkrn(TT,QQ,PP,PRDKRN,TTCOAL)
implicit none
real epsf,tc,ttt1,ttt,factor,qs2,qq1,dele,f,factor_t
double precision TT,QQ,PP,satq2,t,p
double precision prdkrn
REAL at,bt,ct,dt,temp,a,b,c,d,tc_min,tc_max
real factor_max,factor_min
REAL TTCOAL
data at,bt,ct,dt/0.88333,0.0931878,0.0034793,4.5185186e-05/

satq2(t,p)=3.80e3*(10**(9.76421-2667.1/t))/p
temp(a,b,c,d,tc)=d*tc*tc*tc+c*tc*tc+b*tc+a

IF (QQ.LE.0)QQ=1.E-12
epsf=.5
tc=tt-273.15
if(tc.le.0) then
  ttt1  =temp(at,bt,ct,dt,tc)
  ttt   =ttt1
  qs2   =satq2(tt,pp)
  qq1   =qq*(0.622+0.378*qs2)/(0.622+0.378*qq)/qs2
  dele  =ttt*qq1
  if(tc.ge.-6.) then
    factor = dele
    if(factor.lt.epsf) factor=epsf
    if(factor.gt.1.) factor=1.
  endif                        
  factor_t=factor
  if(tc.ge.-12.5.and.tc.lt.-6.) factor_t=0.5
  if(tc.ge.-17.0.and.tc.lt.-12.5) factor_t=1.
  if(tc.ge.-20.0.and.tc.lt.-17.) factor_t=0.4
  if(tc.lt.-20.) then
    tc_min=ttcoal-273.15
    tc_max=-20.
    factor_max=0.25
    factor_min=0.
    f=factor_min+(tc-tc_min)*(factor_max-factor_min)/  &
                                    (tc_max-tc_min)
    factor_t=f
  endif
  if (factor_t.lt.0)factor_t=0.01
  prdkrn=factor_t
else
  prdkrn=1.d0
end if
RETURN
END SUBROUTINE modkrn 
!-------------------------------------------------------------
subroutine coll_xxy(gx,gy,ckxx,x,chucm,ima,prdkrn,nkr)
implicit double precision (a-h,o-z)
dimension chucm(nkr,nkr)
double precision ima(nkr,nkr)
dimension  gx(nkr),gy(nkr),ckxx(nkr,nkr),x(0:nkr)

gmin=1.d-60
! lower and upper integration limit ix0,ix1
do i=1,nkr-1
   ix0=i
   if(gx(i).gt.gmin) exit
enddo

if(ix0.eq.nkr-1) return

do i=nkr-1,1,-1
   ix1=i
   if(gx(i).gt.gmin) exit
enddo

! collisions
do i=ix0,ix1
   do j=i,ix1
      k=ima(i,j)
      kp=k+1
      x0=ckxx(i,j)*gx(i)*gx(j)*prdkrn
      x0=min(x0,gx(i)*x(j))
      x0=min(x0,gx(j)*x(i))
      gsi=x0/x(j)
      gsj=x0/x(i)
      gsk=gsi+gsj
      gx(i)=gx(i)-gsi
      if(gx(i).lt.0.d0) gx(i)=0.d0
      gx(j)=gx(j)-gsj
      if(gx(j).lt.0.d0) gx(j)=0.d0
      gk=gy(k)+gsk
      flux=0.d0
      if(gk.gt.gmin) then
         x1=dlog(gy(kp)/gk+1.d-15)
         if (x1.eq.0)then
            flux=0  
         else
            flux=gsk/x1*(dexp(0.5*x1)-dexp(x1*(0.5-chucm(i,j))))
            flux=min(flux,gsk)
         end if
         gy(k)=gk-flux
         if(gk .lt. flux) flux=gk
         if(gy(k).lt.0.d0) gy(k)=0.d0
         gy(kp)=gy(kp)+flux
      endif
   enddo
enddo
return
end subroutine coll_xxy
!====================================================================
subroutine coll_xyy(gx,gy,ckxy,x,y,chucm,ima, &
  prdkrn,nkr,indc)
implicit double precision (a-h,o-z)
dimension gy(nkr),gx(nkr),ckxy(nkr,nkr),x(0:nkr),y(0:nkr)
dimension chucm(nkr,nkr)
double precision ima(nkr,nkr)

gmin=1.d-60
! lower and upper integration limit ix0,ix1
do i=1,nkr-1
   ix0=i
   if(gx(i).gt.gmin) exit
enddo

if(ix0.eq.nkr-1) return

do i=nkr-1,1,-1
   ix1=i
   if(gx(i).gt.gmin) exit
enddo

! lower and upper integration limit iy0,iy1
do i=1,nkr-1
   iy0=i
   if(gy(i).gt.gmin) exit
enddo

if(iy0.eq.nkr-1) return

do i=nkr-1,1,-1
   iy1=i
   if(gy(i).gt.gmin) exit
enddo

! collisions :
do i=iy0,iy1
   jmin=i
   if(jmin.eq.(nkr-1)) return
   if(i.lt.ix0) jmin=ix0-indc
   do j=jmin+indc,ix1         
      k=ima(i,j)
      kp=k+1
      x0=ckxy(j,i)*gy(i)*gx(j)*prdkrn
      x0=min(x0,gy(i)*x(j))
      x0=min(x0,gx(j)*y(i))
      gsi=x0/x(j)
      gsj=x0/y(i)
      gsk=gsi+gsj
      gy(i)=gy(i)-gsi
      if(gy(i).lt.0.d0) gy(i)=0.d0
      gx(j)=gx(j)-gsj
      if(gx(j).lt.0.d0) gx(j)=0.d0
      gk=gy(k)+gsk
      flux=0.d0
      if(gk.gt.gmin) then
        x1=dlog(gy(kp)/gk+1.d-15)
        if (x1.eq.0)then
           flux=0  
        else
           flux=gsk/x1*(dexp(0.5*x1)-dexp(x1*(0.5-chucm(i,j))))
           flux=min(flux,gsk)
        end if
        gy(k)=gk-flux
        if(gk .lt. flux) flux=gk
        if(gy(k).lt.0.d0) gy(k)=0.d0
        gy(kp)=gy(kp)+flux
      endif
   enddo
enddo
return
end subroutine coll_xyy
!=================================================================
subroutine coll_xyx(gx,gy,ckxy,x,y,chucm,ima, &
 prdkrn,nkr,indc)
implicit double precision (a-h,o-z)
dimension gy(nkr),gx(nkr),ckxy(nkr,nkr),x(0:nkr),y(0:nkr)
dimension chucm(nkr,nkr)
double precision ima(nkr,nkr)

gmin=1.d-60
! lower and upper integration limit ix0,ix1
do i=1,nkr-1
   ix0=i
   if(gx(i).gt.gmin) exit
enddo

if(ix0.eq.nkr-1) return

do i=nkr-1,1,-1
   ix1=i
   if(gx(i).gt.gmin) exit
enddo

! lower and upper integration limit iy0,iy1
do i=1,nkr-1
   iy0=i
   if(gy(i).gt.gmin) exit
enddo

if(iy0.eq.nkr-1) return

do i=nkr-1,1,-1
   iy1=i
   if(gy(i).gt.gmin) exit
enddo

! collisions :
do i=iy0,iy1
   jmin=i
   if(jmin.eq.(nkr-1)) return
   if(i.lt.ix0) jmin=ix0-indc
   do j=jmin+indc,ix1
      k=ima(i,j)
      kp=k+1
      x0=ckxy(j,i)*gy(i)*gx(j)*prdkrn
      x0=min(x0,gy(i)*x(j))
      if(j.ne.k) then
         x0=min(x0,gx(j)*y(i))
      endif
      gsi=x0/x(j)
      gsj=x0/y(i)
      gsk=gsi+gsj
      gy(i)=gy(i)-gsi
      if(gy(i).lt.0.d0) gy(i)=0.d0
      gx(j)=gx(j)-gsj
      gk=gx(k)+gsk
      if(gx(j).lt.0.d0.and.gk.lt.gmin) then
         gx(j)=0.d0
         gx(k)=gx(k)+gsi
      endif
      flux=0.d0            
      if(gk.gt.gmin) then
         x1=dlog(gx(kp)/gk+1.d-15)
         if (x1.eq.0)then
             flux=0  
         else
             flux=gsk/x1*(dexp(0.5*x1)-dexp(x1*(0.5-chucm(i,j))))
             flux=min(flux,gsk)
         end if
         gx(k)=gk-flux
         if(gk .lt. flux) flux=gk

         if(gx(k).lt.0.d0) gx(k)=0.d0
         gx(kp)=gx(kp)+flux
      endif
   enddo
enddo
return
end subroutine coll_xyx
!=====================================================================
subroutine coll_xyxz(gx,gy,gz,ckxy,x,y,chucm,ima, &
 prdkrn,nkr,indc, kp_bound)
implicit double precision (a-h,o-z)
dimension gy(nkr),gx(nkr),gz(nkr),ckxy(nkr,nkr),x(0:nkr),y(0:nkr)
dimension chucm(nkr,nkr)
double precision ima(nkr,nkr)
integer kp_bound

gmin=1.d-60
! lower and upper integration limit ix0,ix1
do i=1,nkr-1
   ix0=i
   if(gx(i).gt.gmin) exit
enddo

if(ix0.eq.nkr-1) return

do i=nkr-1,1,-1
   ix1=i
   if(gx(i).gt.gmin) exit
enddo

! lower and upper integration limit iy0,iy1
do i=1,nkr-1
   iy0=i
   if(gy(i).gt.gmin) exit
enddo

if(iy0.eq.nkr-1) return

do i=nkr-1,1,-1
   iy1=i
   if(gy(i).gt.gmin) exit
enddo

! collisions :
do i=iy0,iy1
   jmin=i
   if(jmin.eq.(nkr-1)) return
   if(i.lt.ix0) jmin=ix0-indc
   do j=jmin+indc,ix1
      k=ima(i,j)
      kp=k+1
      x0=ckxy(j,i)*gy(i)*gx(j)*prdkrn
      x0=min(x0,gy(i)*x(j))
      if(j.ne.k) then
         x0=min(x0,gx(j)*y(i))
      endif
      gsi=x0/x(j)
      gsj=x0/y(i)
      gsk=gsi+gsj
      gy(i)=gy(i)-gsi
      if(gy(i).lt.0.d0) gy(i)=0.d0
      gx(j)=gx(j)-gsj
      gk=gx(k)+gsk
      if(gx(j).lt.0.d0.and.gk.lt.gmin) then
        gx(j)=0.d0
        gx(k)=gx(k)+gsi
      endif
      flux=0.d0
      if(kp.lt.kp_bound) gkp=gx(kp)
      if(kp.ge.kp_bound) gkp=gz(kp)
      if(gk.gt.gmin) then
         x1=dlog(gkp/gk+1.d-15)
         if (x1.eq.0)then
            flux=0  
         else
            flux=gsk/x1*(dexp(0.5*x1)-dexp(x1*(0.5-chucm(i,j))))
            flux=min(flux,gsk)
         end if
         gx(k)=gk-flux
         if(gk .lt. flux) flux=gk

         if(gx(k).lt.0.d0) gx(k)=0.d0
         if(kp.lt.kp_bound) gx(kp)=gkp+flux
         if(kp.ge.kp_bound) gz(kp)=gkp+flux

      endif
   enddo
enddo
return
end subroutine coll_xyxz
!=====================================================================
subroutine coll_xyz(gx,gy,gz,ckxy,x,y,chucm,ima, &
                    prdkrn,nkr,indc)
implicit double precision (a-h,o-z)
dimension gx(nkr),gy(nkr),gz(nkr),ckxy(nkr,nkr),x(0:nkr),y(0:nkr)
dimension chucm(nkr,nkr)
double precision ima(nkr,nkr)

gmin=1.d-60
! lower and upper integration limit ix0,ix1
do i=1,nkr-1
   ix0=i
   if(gx(i).gt.gmin) exit
enddo

if(ix0.eq.nkr-1) return

do i=nkr-1,1,-1
   ix1=i
   if(gx(i).gt.gmin) exit
enddo
! lower and upper integration limit iy0,iy1
do i=1,nkr-1
   iy0=i
   if(gy(i).gt.gmin) exit
enddo

if(iy0.eq.nkr-1) return

do i=nkr-1,1,-1
   iy1=i
   if(gy(i).gt.gmin) exit
enddo

! collisions :
do i=iy0,iy1
   jmin=i
   if(jmin.eq.(nkr-1)) return
   if(i.lt.ix0) jmin=ix0-indc
    do j=jmin+indc,ix1         
      k=ima(i,j)
      kp=k+1
      x0=ckxy(j,i)*gy(i)*gx(j)*prdkrn
      x0=min(x0,gy(i)*x(j))
      x0=min(x0,gx(j)*y(i))
      gsi=x0/x(j)
      gsj=x0/y(i)
      gsk=gsi+gsj
      gy(i)=gy(i)-gsi
      if(gy(i).lt.0.d0) gy(i)=0.d0
      gx(j)=gx(j)-gsj
      if(gx(j).lt.0.d0) gx(j)=0.d0
      gk=gz(k)+gsk
      flux=0.d0
      if(gk.gt.gmin) then
         x1=dlog(gz(kp)/gk+1.d-15)
         if (x1.eq.0)then
            flux=0  
         else
            flux=gsk/x1*(dexp(0.5*x1)-dexp(x1*(0.5-chucm(i,j))))
            flux=min(flux,gsk)
         end if
         gz(k)=gk-flux
         if(gk .lt. flux) flux=gk
!
         if(gz(k).lt.0.d0) gz(k)=0.d0
         gz(kp)=gz(kp)+flux
      endif
   enddo
enddo
return
end subroutine coll_xyz

!========================================================================
subroutine coll_xyyz(gx,gy,gz,ckxy,x,y,chucm,ima, &
                     prdkrn,nkr,indc, kp_bound)
implicit double precision (a-h,o-z)
dimension gy(nkr),gx(nkr),gz(nkr),ckxy(nkr,nkr),x(0:nkr),y(0:nkr)
dimension chucm(nkr,nkr)
double precision ima(nkr,nkr)
integer kp_bound

gmin=1.d-60
! lower and upper integration limit ix0,ix1
do i=1,nkr-1
   ix0=i
   if(gx(i).gt.gmin) exit
enddo

if(ix0.eq.nkr-1) return

do i=nkr-1,1,-1
   ix1=i
   if(gx(i).gt.gmin) exit
enddo
! lower and upper integration limit iy0,iy1
do i=1,nkr-1
   iy0=i
   if(gy(i).gt.gmin) exit
enddo

if(iy0.eq.nkr-1) return

do i=nkr-1,1,-1
   iy1=i
   if(gy(i).gt.gmin) exit
enddo

gsi=0.
gsj=0.

! collisions :
do i=iy0,iy1
   jmin=i
   if(jmin.eq.(nkr-1)) return
   if(i.lt.ix0) jmin=ix0-indc
   do j=jmin+indc,ix1
      k=ima(i,j)
      kp=k+1
      x0=ckxy(j,i)*gy(i)*gx(j)*prdkrn
      x0=min(x0,gy(i)*x(j))
!Adele - one line
      x0=min(x0,gx(j)*y(i))
!      if(j.ne.k) then
!        x0=min(x0,gx(j)*y(i))
!      endif
      gsi=x0/x(j)
      gsj=x0/y(i)
      gsk=gsi+gsj
      gy(i)=gy(i)-gsi
      if(gy(i).lt.0.d0) gy(i)=0.d0
      gx(j)=gx(j)-gsj
!Adele - one line
      if(gx(j).lt.0.d0) gx(j)=0.d0
      gk=gy(k)+gsk
      flux=0.d0
      if(kp.lt.kp_bound) gkp=gy(kp)
      if(kp.ge.kp_bound) gkp=gz(kp)
      if(gk.gt.gmin) then
        x1=dlog(gkp/gk+1.d-15)
        if (x1.eq.0)then
           flux=0  
        else
           flux=gsk/x1*(dexp(0.5*x1)-dexp(x1*(0.5-chucm(i,j))))
           flux=min(flux,gsk)
        end if
        gy(k)=gk-flux
        if(gk .lt. flux) flux=gk

        if(gy(k).lt.0.d0) gy(k)=0.d0
        if(kp.lt.kp_bound) gy(kp)=gkp+flux
        if(kp.ge.kp_bound) gz(kp)=gkp+flux
      endif
   enddo
enddo
return
end subroutine coll_xyyz

!===============================================================
SUBROUTINE BREAKUP(GT_MG,XT_MG,DT,BRKWEIGHT, &
                   PKIJ,QKJ,JMAX,JBREAK)
!.....INPUT VARIABLES
!
! GT    : MASS DISTRIBUTION FUNCTION
! XT_MG : MASS OF BIN IN MG
! JMAX  : NUMBER OF BINS
! DT    : TIMESTEP IN S

INTEGER JMAX

!.....LOCAL VARIABLES

INTEGER JBREAK,AP,IA,JA,KA,IE,JE,KE
DOUBLE PRECISION EPS,NEGSUM

PARAMETER (AP = 1)
PARAMETER (IA = 1)
PARAMETER (JA = 1)
PARAMETER (KA = 1)
PARAMETER (EPS = 1.D-20)

INTEGER I,J,K,JDIFF
DOUBLE PRECISION GT_MG(JMAX),XT_MG(0:JMAX),DT
DOUBLE PRECISION BRKWEIGHT(JBREAK),PKIJ(JBREAK,JBREAK,JBREAK), &
                 QKJ(JBREAK,JBREAK)
DOUBLE PRECISION HLP(JMAX)
DOUBLE PRECISION FT(JMAX),FA(JMAX)
DOUBLE PRECISION DG(JMAX),DF(JMAX),DBREAK(JBREAK),GAIN,LOSS
REAL PI
PARAMETER (PI = 3.1415927)
IE = JBREAK
JE = JBREAK
KE = JBREAK
!.....IN CGS

!.....SHIFT BETWEEN COAGULATION AND BREAKUP GRID

JDIFF = JMAX - JBREAK

!.....INITIALIZATION

!.....TRANSFORMATION FROM G(LN X) = X**2 F(X) TO F(X)
DO J=1,JMAX
   FT(J) = GT_MG(J) / XT_MG(J)**2
ENDDO

!.....SHIFT TO BREAKUP GRID

DO K=1,KE
   FA(K) = FT(K+JDIFF)
ENDDO

!.....BREAKUP: BLECK'S FIRST ORDER METHOD
!
!     PKIJ: GAIN COEFFICIENTS
!     QKJ : LOSS COEFFICIENTS
!

DO K=1,KE
   GAIN = 0.0
   DO I=1,IE
      DO J=1,I
         GAIN = GAIN + FA(I)*FA(J)*PKIJ(K,I,J)
      ENDDO
   ENDDO
   LOSS = 0.0
   DO J=1,JE
      LOSS = LOSS + FA(J)*QKJ(K,J)
   ENDDO
   DBREAK(K) = BRKWEIGHT(K) * (GAIN - FA(K)*LOSS)
ENDDO

!.....SHIFT RATE TO COAGULATION GRID

DO J=1,JDIFF
   DF(J) = 0.0
ENDDO
DO J=1,KE
   DF(J+JDIFF) = DBREAK(J)
ENDDO
!.....TRANSFORMATION TO MASS DISTRIBUTION FUNCTION G(LN X)

DO J=1,JMAX
   DG(J) = DF(J) * XT_MG(J)**2
ENDDO

!.....TIME INTEGRATION

DO J=1,JMAX
HLP(J) = 0.0
NEGSUM = 0.0
GT_MG(J) = GT_MG(J) + DG(J) * DT
IF (GT_MG(J).LT.0) THEN
   HLP(J) = MIN(GT_MG(J),HLP(J))
   GT_MG(J) = EPS
ENDIF
ENDDO
RETURN
END SUBROUTINE BREAKUP
  
END MODULE module_hujisbm
