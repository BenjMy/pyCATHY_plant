C
C**************************  CHPIC4 ************************************
C
C  calculate soil moisture characteristics needed for Picard scheme
C  KSLOPE=4 : localized tangent slope differentiation of moisture 
C  curves.
C  Note that DSETAN contains tangent slope values at each node only for 
C  the case IVGHU=1; for the other IVGHU cases the tangent slope is 
C  constant for all nodes and we use DSETAN(1).
C
C***********************************************************************
C
      SUBROUTINE CHPIC4(N,PTNEW,SNODI,PNODI,SW,CKRW,ETAI,IVGHU)
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
            PSI       = PTNEW(I)
            SE        = FVGSE(PSI)
            SW(I)     = VGPNOT(I)*SE + VGRMC/PNODI(I)
            CKRW(I)   = FVGKR(PSI,SE)
            IF (PSI .GE. PSEL  .AND.  PSI .LE. PSER) THEN
               ETAI(I)= SW(I)*SNODI(I) + PNODI(I)*VGPNOT(I)*DSETAN(1)
            ELSE
               ETAI(I)= SW(I)*SNODI(I) + PNODI(I)*VGPNOT(I)*FVGDSE(PSI)
            END IF
         END DO
      ELSE IF (IVGHU .EQ. 1) THEN
         DO I=1,N
            PSI       = PTNEW(I)
            SW(I)     = FXVMC(PSI,SNODI(I),PNODI(I),I)/PNODI(I)
            CKRW(I)   = FXVKR(PSI)
            IF (PSI .GE. PSEL  .AND.  PSI .LE. PSER) THEN
               ETAI(I)= DSETAN(I)
            ELSE
               ETAI(I)= FXVDMC(PSI,SNODI(I),PNODI(I),I)
            END IF
         END DO
      ELSE IF (IVGHU .EQ. 2) THEN
         DO I=1,N
            PSI       = PTNEW(I)
            SE        = FHUSE(PSI)
            SW(I)     = HUSWR1*SE + HUSWR
            CKRW(I)   = FHUKR2(PSI,SE)
            IF (PSI .GE. PSEL  .AND.  PSI .LE. PSER) THEN
               ETAI(I)= SW(I)*SNODI(I) + PNODI(I)*HUSWR1*DSETAN(1)
            ELSE
               ETAI(I)= SW(I)*SNODI(I) + PNODI(I)*HUSWR1*FHUDSE(PSI)
            END IF
         END DO
      ELSE IF (IVGHU .EQ. 3) THEN
         DO I=1,N
            PSI       = PTNEW(I)
            SE        = FHUSE(PSI)
            SW(I)     = HUSWR1*SE + HUSWR
            CKRW(I)   = FHUKR3(PSI,SE)
            IF (PSI .GE. PSEL  .AND.  PSI .LE. PSER) THEN
               ETAI(I)= SW(I)*SNODI(I) + PNODI(I)*HUSWR1*DSETAN(1)
            ELSE
               ETAI(I)= SW(I)*SNODI(I) + PNODI(I)*HUSWR1*FHUDSE(PSI)
            END IF
         END DO
      ELSE
         DO I=1,N
            PSI       = PTNEW(I)
            SE        = FBCSE(PSI)
            SW(I)     = BCPORM(I)*SE + BCRMC/PNODI(I)
            CKRW(I)   = FBCKR(PSI)
            IF (PSI .GE. PSEL  .AND.  PSI .LE. PSER) THEN
               ETAI(I)= SW(I)*SNODI(I) + PNODI(I)*BCPORM(I)*DSETAN(1)
            ELSE
               ETAI(I)= SW(I)*SNODI(I) + PNODI(I)*BCPORM(I)*FBCDSE(PSI)
            END IF
         END DO
      END IF
C
      RETURN
      END
