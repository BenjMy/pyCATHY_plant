C
C**************************  FXVDMC ************************************
C
C  calculate derivative of moisture content with respect to pressure 
C  head using extended van Genuchten characteristic equation.
C  for the extended van Genuchten equation this gives the value of the
C  overall storage coefficient, or general storage term.
C
C***********************************************************************
C
      DOUBLE PRECISION FUNCTION FXVDMC(PSI,SS,POR,INOD)
C
      IMPLICIT NONE
      INCLUDE 'CATHY.H'
      INTEGER  INOD
      REAL*8   BETA,B1,B1R
      REAL*8   PSI,SS,POR
      INCLUDE 'SOILCHAR.H'
C
      IF (PSI .LT. VGPNOT(INOD)) THEN
         BETA=(PSI/VGPSAT)**VGN
         B1=BETA+1.0D0
         B1R=1.0D0/B1
         FXVDMC=VGN1*(POR-VGRMC)*((DABS(PSI)**VGN1)/VGPSN)*
     1               (B1**VGNR)*B1R*B1R
      ELSE
         FXVDMC=SS
      END IF
C
      RETURN
      END
