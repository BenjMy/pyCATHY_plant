C
C**************************  FVGDKR ************************************
C
C  calculate derivative of relative hydraulic conductivity with respect
C  to pressure head using van Genuchten characteristic equation.
C  We use the same calculation as for the extended van Genuchten
C  curves (i.e. the derivative is formulated in terms of PSI rather than
C  in terms of SE and DSEDP) in order to avoid division by zero.
C
C***********************************************************************
C
      DOUBLE PRECISION FUNCTION FVGDKR(PSI)
C
      IMPLICIT NONE
      INCLUDE 'CATHY.H'
      REAL*8   BETA,B1,V1,V2,V3,V4
      REAL*8   PSI
      INCLUDE 'SOILCHAR.H'
C
      IF (PSI .LT. -1.0D-14) THEN
         BETA=(ABS(PSI/VGPSAT))**VGN
         B1=BETA+1.0D0
         V1=ABS(B1)**VGM - ABS(BETA)**VGM
         V2=VGPSAT/PSI
         V3=VGN1*BETA*V2*(ABS(1.0D0/B1)**VGM52)/VGPSAT
         V4=V2*((2.5D0/B1)*BETA - 2.0D0) - 0.5D0*(ABS(B1)**VGMM1)
         FVGDKR=V3*V1*V4
      ELSE
         FVGDKR=0.0D0
      END IF
C
      RETURN
      END
