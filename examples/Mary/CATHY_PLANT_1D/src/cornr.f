C
C*************************** CORNR **************************************
C
C Correction for the Newton-Raphson method to inizialize the void ratio 
C distribution in subroutine CHVELOP 
C
C***********************************************************************
C
      DOUBLE PRECISION FUNCTION CORNR(INDEi,SWi,D,PSI,PNODIi,
     1                                SNODIi)
C
      IMPLICIT NONE
      INCLUDE 'CATHY.H'
      REAL*8   INDEi,SWi,D,PSI,PNODIi,SNODIi
      REAL*8   CBETA,THETA0
      REAL*8   UNO
      PARAMETER (UNO=1.0D0)
      INCLUDE 'SOILCHAR.H'
C
      IF (SWi .LT. UNO) THEN
         CBETA=CBETA0+CANG*D 
         THETA0=PNODIi/(UNO-PNODIi)
         IF (CBETA .GT. UNO) THEN
            CBETA=UNO
         END IF
         CORNR=(INDEi - ABS(THETA0+UNO)**(UNO-CBETA) *
     1          (INDEi*SWi+UNO)**CBETA + UNO)/
     2         (UNO - ABS(THETA0+UNO)**(UNO-CBETA) * CBETA *
     2          (INDEi*SWi+UNO)**(CBETA-UNO) * SWi)
      ELSE
cm       CORNR=(INDEi - (SNODIi*PSI + PNODIi) * (UNO+INDEi))/
cm   1         (UNO - SNODIi * PSI - PNODIi)
         CORNR=0.0D0
      END IF
C
      RETURN
      END











