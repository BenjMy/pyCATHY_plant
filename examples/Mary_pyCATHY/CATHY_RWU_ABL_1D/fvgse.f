C
C**************************  FVGSE *************************************
C
C  calculate effective saturation using van Genuchten 
C  characteristic equation
C
C***********************************************************************
C
      DOUBLE PRECISION FUNCTION FVGSE(PSI)
C
      IMPLICIT NONE
      INCLUDE 'CATHY.H'
      REAL*8   BETA
      REAL*8   PSI
      INCLUDE 'SOILCHAR.H'
C
      IF (PSI .LT. -1.0D-14) THEN
         BETA=(ABS(PSI/VGPSAT))**VGN
         FVGSE=ABS(1.0D0/(BETA+1.0D0))**VGM
      ELSE
         FVGSE=1.0D0
      END IF
C
      RETURN
      END
