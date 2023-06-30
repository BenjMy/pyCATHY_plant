C
C**************************  FVGDSE ************************************
C
C  calculate derivative of effective saturation with respect to pressure 
C  head using van Genuchten characteristic equation
C
C***********************************************************************
C
      DOUBLE PRECISION FUNCTION FVGDSE(PSI)
C
      IMPLICIT NONE
      INCLUDE 'CATHY.H'
      REAL*8   BETA,B1,B1R
      REAL*8   PSI
      INCLUDE 'SOILCHAR.H'
C
      IF (PSI .LT. -1.0D-14) THEN
         BETA=(ABS(PSI/VGPSAT))**VGN
         B1=BETA+1.0D0
         B1R=1.0D0/B1
         FVGDSE=VGN1*((ABS(PSI)**VGN1)/VGPSN)*(ABS(B1)**VGNR)*B1R*B1R
      ELSE
         FVGDSE=0.0D0
      END IF
C
      RETURN
      END
