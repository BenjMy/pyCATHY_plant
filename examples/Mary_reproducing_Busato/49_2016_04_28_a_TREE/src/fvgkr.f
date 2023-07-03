C
C**************************  FVGKR *************************************
C
C  calculate relative hydraulic conductivity using van Genuchten 
C  characteristic equation
C
C***********************************************************************
C
      DOUBLE PRECISION FUNCTION FVGKR(PSI,SE)
C
      IMPLICIT NONE
      INCLUDE 'CATHY.H'
      REAL*8   OMEGA,V1
      REAL*8   PSI,SE
      INCLUDE 'SOILCHAR.H'
C
      IF (PSI .LT. -1.0D-14) THEN
         OMEGA=ABS(SE)**VGMR
         V1=1.0D0 - (ABS(1.0D0 - OMEGA)**VGM)
         FVGKR=DSQRT(SE)*V1*V1
      ELSE
         FVGKR=1.0D0
      END IF
C
      RETURN
      END
