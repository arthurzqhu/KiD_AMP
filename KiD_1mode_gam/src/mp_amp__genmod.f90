        !COMPILER-GENERATED INTERFACE MODULE: Wed Jun 15 06:12:36 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE MP_AMP__genmod
          INTERFACE 
            SUBROUTINE MP_AMP(MPC,MPR,GUESSC,GUESSR,PRESS,TEMPK,QV,FNCN,&
     &FFCDR8_MASS,FFCDR8_NUM,MC,MR,FLAG,FFCDR8_MASSINIT,FFCDR8_NUMINIT)
              REAL(KIND=8) :: MPC(120,1,NUM_H_MOMENTS((1)))
              REAL(KIND=8) :: MPR(120,1,NUM_H_MOMENTS((2)))
              REAL(KIND=8) :: GUESSC(120,1,2)
              REAL(KIND=8) :: GUESSR(120,1,2)
              REAL(KIND=4) :: PRESS(120,1)
              REAL(KIND=4) :: TEMPK(120,1)
              REAL(KIND=4) :: QV(120,1)
              REAL(KIND=8) :: FNCN(120,1,34)
              REAL(KIND=8) :: FFCDR8_MASS(120,1,34)
              REAL(KIND=8) :: FFCDR8_NUM(120,1,34)
              REAL(KIND=8) :: MC(120,1,1:10)
              REAL(KIND=8) :: MR(120,1,1:10)
              REAL(KIND=8) :: FLAG(120,1,2,4)
              REAL(KIND=8) :: FFCDR8_MASSINIT(120,1,34)
              REAL(KIND=8) :: FFCDR8_NUMINIT(120,1,34)
            END SUBROUTINE MP_AMP
          END INTERFACE 
        END MODULE MP_AMP__genmod
