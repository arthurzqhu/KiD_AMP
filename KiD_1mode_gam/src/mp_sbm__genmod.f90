        !COMPILER-GENERATED INTERFACE MODULE: Wed Jun 15 06:12:36 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE MP_SBM__genmod
          INTERFACE 
            SUBROUTINE MP_SBM(FFCDR8,PRESS,TEMPK,QV,FNCNR8,MC,MR)
              REAL(KIND=8) :: FFCDR8(120,1,34)
              REAL(KIND=4) :: PRESS(120,1)
              REAL(KIND=4) :: TEMPK(120,1)
              REAL(KIND=4) :: QV(120,1)
              REAL(KIND=8) :: FNCNR8(120,1,34)
              REAL(KIND=8) :: MC(120,1,10)
              REAL(KIND=8) :: MR(120,1,10)
            END SUBROUTINE MP_SBM
          END INTERFACE 
        END MODULE MP_SBM__genmod
