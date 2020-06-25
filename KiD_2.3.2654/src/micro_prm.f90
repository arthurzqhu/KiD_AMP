module micro_prm
use parameters, only: max_nbins
use namelists, only:imomc1,imomc2,imomr1,imomr2,donucleation, &
                    docondensation,docollisions,dosedimentation, &
                    cloud_init,rain_init,bintype
!use parameters, only:nx,nz,num_h_moments
implicit none

!*****Variables originally from RAMS*******!
real, parameter ::                    &
        alvl     = 2.50e6             &
    ,   alvi     = 2.834e6            &
    ,   alli     = 0.334e6            &
    ,   cpi      = 1. / 1004.
!real(8) :: press,tempk,rhw
logical :: lhrtheta = .true.
integer :: imbudget,icloud,idriz,irain,ipris,isnow,iaggr,igraup,ihail
!******Variables Needed for COMPUTING BUDGETS ******************************
!For imbudget>=1
real :: latheatvap,latheatfrz,nuccldrt,cld2raint &
,ice2raint,nucicert,vapliqt,vapicet,melticet,rimecldt,aggregatet &
,rain2icet,latheatvapt,latheatfrzt,nuccldct,nucicect

!For imbudget>=2
real :: inuchomrt,inuccontrt,inucifnrt,inucifnct,inuchazrt   &
,vapcldt,vapraint,vapprist,vapsnowt,vapaggrt,vapgraut,vaphailt     &
,vapdrizt,meltprist,meltsnowt,meltaggrt,meltgraut,melthailt         &
,rimecldsnowt,rimecldaggrt,rimecldgraut,rimecldhailt,rain2prt        &
,rain2snt,rain2agt,rain2grt,rain2hat,aggrselfprist                   &
,aggrselfsnowt,aggrprissnowt

!******Variables for SBM*******!
logical :: iprint
! number of bins
integer :: nkr!=33

! for bulk nucleation
integer, parameter :: BULKNUC=0

! COL=Ln2/3 (for xl(k+1)=xl(k)*2.)
real :: COL = 0.231049060186648 !=log(2.)/3.
real :: COL3 = 0.693147180559945  !=log(2.)

! Aerosols
!  add for diagnostic CCN. If diagCCN=.false., CCN is prognostic
LOGICAL, PARAMETER :: diagCCN=.false.

! Parameters are used for calculation of

! RO_SOLUTE - density of aerosol+h20
! ROCCN0, g/cm^3 - density of aerosol
! Molecular weight of dry aerosol
! Number of ions
real :: RO_SOLUTE,ROCCN0,MWAERO,IONS

! ACCN in 1/cm^3, BCCN - coefficients for calculation of FCCNR(KR)
! FCCNR(KR)=1.5*ACCN*BCCN*S_KR**BCCN,
! where S_KR - water supersaturation, %
real, parameter :: ACCN=1.0000E02, BCCN=0.4620E00

! ICCN - set maximum mass of aerosol
! XL(ICCN) - maximum mass of aerosol, is determined as minimum
! mass drop in bin with number ICCN
integer, parameter :: ICCN=1

!c if ICEPROCS=1 it is ice microphysics
!ICEPROCS is set to one by setting any of the ice flags in RAMSIN to a value above 0
integer :: ICEPROCS,IFASTSBM

! ICEFLAG=0 USE MEYERS ICE NUCLEATION SCHEME; ICEFLAG=1 USE Classical Theory
integer :: ICEFLAG

!c ICEMAX - number of ice crystal types
integer, parameter :: ICEMAX=3

!c flags for turbulance in collection kernals
integer, parameter :: LIQTURB=1
integer, parameter :: ICETURB=1
integer, parameter :: NHYDR=5,NHYDRO=7                &
             ,K0_LL=8,KRMIN_LL=1,KRMAX_LL=19,L0_LL=6  &
             ,IEPS_400=1,IEPS_800=0,IEPS_1600=0       & !For turbulent enhancement to collection kernal
             ,K0L_GL=16,K0G_GL=16                     & !Bins that have turbulent enhancement
             ,KRMINL_GL=1,KRMAXL_GL=24                &
             ,KRMING_GL=1,KRMAXG_GL=33                &
             ,KRDROP=14                               & !First Bin number that is rain in FFCD
             ,JMAX=33,JBREAK = 18                       !For rain breakup routines
integer,dimension(ICEMAX), parameter::    &
         KRPRIS=(/14,16,16/)                              !First Bin number that is RAMS snow for each
                                                        !of the 3 ice crystal types

!Coefficient for deposition-condensation freezing, and contact nuclei fomrulations from Meyers, respectively.
!The code came to me (Adele) with C2_MEY=0, corresponding to no contact nuclei. This should stay zero
!unless the nucleation code is fixed. It currently does not treat the contact nuclei correctly.
real, parameter :: C1_MEY=0.00012,C2_MEY=0.

! For drop freezing (iceform=1) or drop evaporation nuclei (iceform = 2) mechanism
integer, parameter :: iceform = 2

!c KRFREEZ set bin boundary between plates and hail with freezing
!c AFREEZE, BFREEZE, BFREEZEMAX - coefficients in formula of
!c                                homogeneous freezing
integer, PARAMETER :: KRFREEZ=21
REAL, PARAMETER :: BFREEZMAX=0.66E0

!c other parameters and thresholds
real, parameter :: AFREEZMY=0.3333E-04,BFREEZMY=0.6600E00
real, parameter :: cloud_mr_th=1d-8, rain_mr_th=1d-10 !cloud/rain mixing ratio threshold, in units of kg/g

!c Parameters are used in algorithm of diffusional growth
!c NCOND determine timestep (DT/NCOND) with diffusional growth
integer, parameter :: NCOND=1
real :: DTCOND

!c Coeffients for diffusional growth
real, parameter :: A1_MYN=2.53,BB1_MYN=5.42, A2_MYN=3.41E1,BB2_MYN=6.13
real, parameter :: AA1_MY=2.53E12,BB1_MY=5.42E3, AA2_MY=3.41E13,BB2_MY=6.13E3
real, parameter :: DEL_BB=BB2_MY-BB1_MY, DEL_BBN=BB2_MYN-BB1_MYN, DEL_BBR=BB1_MYN/DEL_BBN

!c COAGULATION
!Do COAL_BOTT_NEW every NDTCOLL time steps
integer :: NDTCOLL
real :: DTCOLL

!c Kernels (collisions), depend on heights
!c P1=1000mb=1000000 dynes/cm^2
real, parameter :: p1=1000000.0,p2=750000.0,p3=500000.0

!c Parameters used into subroutine coal_bott (collisions):
!c if temperature less than TTCOAL then no calculations of
!c temperature dependence of collision kernels
real, parameter :: TTCOAL=233.15E00

!c TCRIT - temperature boundary between graupel and hail with
!c coagulation
real, parameter :: TCRIT=270.15E00

!c alcr=1.0 g/m**3 :threshold for graupel formation snow-drop collisions
real, parameter :: ALCR=1.5

!c  threshold for hail formation graupel-drop collisions
integer, parameter :: alcr_hail=3

!c threshold bin number for hail formation from graupel-drop collisions
integer, parameter :: kp_hail=25

!c scal=1.0 is valid under condition of mass doubling :
!c xl(k+1)=xl(k)*2
real, parameter :: SCAL=1.0E00

!c Parameters used for ice multiplication :
!c if icempl=0 no ice multiplication;
!c if icempl=1 - ice multiplication is included
integer, parameter :: ICEMPL=1
integer, parameter :: kr_icempl=9

!c 3-point remapping
integer, parameter :: I3POINT=1

! For new ice nucleation from J. Comstock (ICEFLAG =1)

!       version
!               adapted from aerosol_prop_mks.h; from J.Comstock
!       purpose
!               -define properties for aerosol composition used in
!               heterogeneous nucleation (KC1999)
!               -define aerosol size distribution parameters
!       NOTE MKS units

!       =========================================================================

real, parameter :: qvaero=0.4    !aerosol solubility in terms of volume fraction (Qv)
real, parameter :: betafr=0.5     !aerosol parameter describing composition
real, parameter :: Ms2=0.132      !molecular weight of dry aerosol (g/mol) (NH4)2SO4 Ammonium Sulfate KC1999
real, parameter :: rhos2=1770.0   !density of dry aerosol (kg/m3) (NH4)2SO4 Ammonium Sulate KC1999
real, parameter ::  wetcoef = 0.90 !wettability
real, parameter :: alf = 0.1e-5   ! relative area of active sites
real, parameter :: epsil = 2.5e-2 ! misfit strain parameter
real :: fracin ! IN fraction

!       ==========================================================================
!!! For ccn regeneration
real :: ccnreg
real :: inreg
!       ================================================
! For idealized ccn distributions
real :: sig_g, rp_g !cm

!
!double precision, dimension(nkr) :: fncn,ffcd,ffip,ffid,ffic,ffgl,ffhl,ffsn
!    real, dimension(nz,nx,num_h_moments(1)) :: Mpc2d
!    real, dimension(nz,nx,num_h_moments(2)) :: Mpr2d
!    real, dimension(nz,nx,2) :: guessc2d,guessr2d
!integer:: imomc1,imomc2,imomr1,imomr2
integer, dimension(3):: pmomsc,pmomsr
!real, dimension(2):: cloud_init,rain_init
double precision :: aeromedrad, naero=0., relax, nug
integer :: aerotype=1,npm
real :: dtlt
double precision :: Mp(3),M3p,Mxp,Myp,rxfinal
logical :: parcel!, docollisions, docondensation, donucleation, dosedimentation
integer :: skr,ekr,momx,momy,ihyd
integer, parameter :: ntab=100
double precision, dimension(ntab,ntab,2) :: nutab,dntab
double precision, dimension(2,10,2) :: minmaxmx
double precision, dimension(ntab,2) :: mintab,maxtab
double precision, dimension(2,2) :: nubounds, dnbounds
real(8),dimension(max_nbins) :: diams

integer, parameter :: r4size = 4, r8size = 8
INTEGER,PARAMETER :: ISIGN_KO_1 = 0, ISIGN_KO_2 = 0,  ISIGN_3POINT = 1,  &
                      IDebug_Print_DebugModule = 1
DOUBLE PRECISION,PARAMETER::COEFF_REMAPING = 0.0066667D0
DOUBLE PRECISION,PARAMETER::VENTPL_MAX = 5.0D0

DOUBLE PRECISION,PARAMETER::RW_PW_MIN = 1.0D-10
DOUBLE PRECISION,PARAMETER::RI_PI_MIN = 1.0D-10
DOUBLE PRECISION,PARAMETER::RW_PW_RI_PI_MIN = 1.0D-10
DOUBLE PRECISION,PARAMETER::RATIO_ICEW_MIN = 1.0D-4

REAL (KIND=R4SIZE) :: FR_LIM(max_nbins), FRH_LIM(max_nbins)

!******* for operating SBM when max_nbins is set to 34 *******
integer :: idx &
           ,otn(33) = (/(idx,idx=1,33,1)/) & !an array of consecutive integer from one to nkr
           ! had to hard code the array size to 33 since nkr is not a parameter -ahu
           ,oti(icemax) = (/(idx,idx=1,icemax,1)/) &
           ,oth(nhydro) = (/(idx,idx=1,nhydro,1)/)
contains
    subroutine check_bintype
    implicit none
 
    if (bintype .eq. 'sbm') then
        nkr=33
    elseif (bintype .eq. 'tau') then
        nkr=34
    endif

    end subroutine check_bintype

end module micro_prm
