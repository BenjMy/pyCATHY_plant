C
C**************************  CHPIC0 ************************************
C
C  calculate soil moisture characteristics needed for Picard scheme
C  KSLOPE=0 : analytical differentiation of moisture curves
C
C***********************************************************************
C
      SUBROUTINE CHPIC0(N,PTNEW,SNODI,PNODI,SW,CKRW,ETAI,IVGHU)
C
      IMPLICIT NONE
      INCLUDE 'CATHY.H'
      INTEGER  I
      INTEGER  N,IVGHU
      REAL*8   FVGKR,FVGSE,FVGDSE,FXVKR,FXVMC,FXVDMC
      REAL*8   FHUKR2,FHUKR3,FHUSE,FHUDSE,FBCKR,FBCSE,FBCDSE
      REAL*8   PSI,SE
      REAL*8   PTNEW(*),SNODI(*),PNODI(*),SW(*),CKRW(*),ETAI(*)
      INCLUDE 'SOILCHAR.H'
C
      IF (IVGHU .EQ. 0) THEN
         DO I=1,N
            PSI      = PTNEW(I)
            SE       = FVGSE(PSI)
            SW(I)    = VGPNOT(I)*SE + VGRMC/PNODI(I)
            ETAI(I)  = SW(I)*SNODI(I) + PNODI(I)*VGPNOT(I)*FVGDSE(PSI)
            CKRW(I)  = FVGKR(PSI,SE)
         END DO
      ELSE IF (IVGHU .EQ. 1) THEN
         DO I=1,N
            PSI      = PTNEW(I)
            SW(I)    = FXVMC(PSI,SNODI(I),PNODI(I),I)/PNODI(I)
            ETAI(I)  = FXVDMC(PSI,SNODI(I),PNODI(I),I)
            CKRW(I)  = FXVKR(PSI)
         END DO
      ELSE IF (IVGHU .EQ. 2) THEN
         DO I=1,N
            PSI      = PTNEW(I)
            SE       = FHUSE(PSI)
            SW(I)    = HUSWR1*SE + HUSWR
            ETAI(I)  = SW(I)*SNODI(I) + PNODI(I)*HUSWR1*FHUDSE(PSI)
            CKRW(I)  = FHUKR2(PSI,SE)
         END DO
      ELSE IF (IVGHU .EQ. 3) THEN
         DO I=1,N
            PSI      = PTNEW(I)
            SE       = FHUSE(PSI)
            SW(I)    = HUSWR1*SE + HUSWR
            ETAI(I)  = SW(I)*SNODI(I) + PNODI(I)*HUSWR1*FHUDSE(PSI)
            CKRW(I)  = FHUKR3(PSI,SE)
         END DO
      ELSE
         DO I=1,N
            PSI      = PTNEW(I)
            SE       = FBCSE(PSI)
            SW(I)    = BCPORM(I)*SE + BCRMC/PNODI(I)
            ETAI(I)  = SW(I)*SNODI(I) + PNODI(I)*BCPORM(I)*FBCDSE(PSI)
            CKRW(I)  = FBCKR(PSI)
         END DO
      END IF
C
      RETURN
      END
