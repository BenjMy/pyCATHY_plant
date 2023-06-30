C
C**************************  FXVDKR ************************************
C
C  calculate derivative of relative hydraulic conductivity with respect
C  to pressure head using extended van Genuchten characteristic 
C  equation
C
C***********************************************************************
C
      DOUBLE PRECISION FUNCTION FXVDKR(PSI)
C
      IMPLICIT NONE
      INCLUDE 'CATHY.H'
      REAL*8   BETA,B1,V1,V2,V3,V4
      REAL*8   PSI
      INCLUDE 'SOILCHAR.H'
C
      IF (PSI .LT. -1.0D-14) THEN
         BETA=(PSI/VGPSAT)**VGN
         B1=BETA+1.0D0
         V1=(B1**VGM) - (BETA**VGM)
         V2=VGPSAT/PSI
         V3=VGN1*BETA*V2*((1.0D0/B1)**VGM52)/VGPSAT
         V4=V2*((2.5D0/B1)*BETA - 2.0D0) - 0.5D0*(B1**VGMM1)
         FXVDKR=V3*V1*V4
      ELSE
         FXVDKR=0.0D0
      END IF
C
      RETURN
      END
