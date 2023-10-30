C
C**************************  CHPIC3 ************************************
C
C  calculate soil moisture characteristics needed for Picard scheme
C  KSLOPE=3 : localized chord slope and analytical differentiation of 
C             moisture curves
C
C***********************************************************************
C
      SUBROUTINE CHPIC3(N,PTNEW,PTOLD,SNODI,PNODI,SW,CKRW,ETAI,IVGHU)
C
      IMPLICIT NONE
      INCLUDE 'CATHY.H'
      INTEGER  I
      INTEGER  N,IVGHU
      REAL*8   FVGKR,FVGSE,FVGDSE,FXVKR,FXVMC,FXVDMC
      REAL*8   FHUKR2,FHUKR3,FHUSE,FHUDSE,FBCKR,FBCSE,FBCDSE
      REAL*8   PSI,POLD,DP,SE,SEOLD,MC,MCOLD
      REAL*8   PTNEW(*),PTOLD(*),SNODI(*),PNODI(*),SW(*),CKRW(*),ETAI(*)
      INCLUDE 'SOILCHAR.H'
C
      IF (IVGHU .EQ. 0) THEN
         DO I=1,N
            PSI          = PTNEW(I)
            SE           = FVGSE(PSI)
            SW(I)        = VGPNOT(I)*SE + VGRMC/PNODI(I)
            CKRW(I)      = FVGKR(PSI,SE)
            IF (PSI .GE. PSEL  .AND.  PSI .LE. PSER) THEN
               POLD      = PTOLD(I)
               DP        = PSI - POLD
               IF (DABS(DP) .LT. TOLKSL) THEN
                  ETAI(I)= SW(I)*SNODI(I) +
     1                     PNODI(I)*VGPNOT(I)*FVGDSE(PSI)
               ELSE
                  SEOLD  = FVGSE(POLD)
                  ETAI(I)= SW(I)*SNODI(I) +
     1                     PNODI(I)*VGPNOT(I)*((SE - SEOLD)/DP)
               END IF
            ELSE
               ETAI(I)   = SW(I)*SNODI(I) +
     1                     PNODI(I)*VGPNOT(I)*FVGDSE(PSI)
            END IF
         END DO
      ELSE IF (IVGHU .EQ. 1) THEN
         DO I=1,N
            PSI          = PTNEW(I)
            SW(I)        = FXVMC(PSI,SNODI(I),PNODI(I),I)/PNODI(I)
            CKRW(I)      = FXVKR(PSI)
            IF (PSI .GE. PSEL  .AND.  PSI .LE. PSER) THEN
               POLD      = PTOLD(I)
               DP        = PSI - POLD
               IF (DABS(DP) .LT. TOLKSL) THEN
                  ETAI(I)= FXVDMC(PSI,SNODI(I),PNODI(I),I)
               ELSE
                  MC     = FXVMC(PSI,SNODI(I),PNODI(I),I)
                  MCOLD  = FXVMC(POLD,SNODI(I),PNODI(I),I)
                  ETAI(I)= (MC-MCOLD)/DP
               END IF
            ELSE
               ETAI(I)   = FXVDMC(PSI,SNODI(I),PNODI(I),I)
            END IF
         END DO
      ELSE IF (IVGHU .EQ. 2) THEN
         DO I=1,N
            PSI          = PTNEW(I)
            SE           = FHUSE(PSI)
            SW(I)        = HUSWR1*SE + HUSWR
            CKRW(I)      = FHUKR2(PSI,SE)
            IF (PSI .GE. PSEL  .AND.  PSI .LE. PSER) THEN
               POLD      = PTOLD(I)
               DP        = PSI - POLD
               IF (DABS(DP) .LT. TOLKSL) THEN
                  ETAI(I)= SW(I)*SNODI(I) +
     1                     PNODI(I)*HUSWR1*FHUDSE(PSI)
               ELSE
                  SEOLD  = FHUSE(POLD)
                  ETAI(I)= SW(I)*SNODI(I) +
     1                     PNODI(I)*HUSWR1*((SE - SEOLD)/DP)
               END IF
            ELSE
               ETAI(I)   = SW(I)*SNODI(I) +
     1                     PNODI(I)*HUSWR1*FHUDSE(PSI)
            END IF
         END DO
      ELSE IF (IVGHU .EQ. 3) THEN
         DO I=1,N
            PSI          = PTNEW(I)
            SE           = FHUSE(PSI)
            SW(I)        = HUSWR1*SE + HUSWR
            CKRW(I)      = FHUKR3(PSI,SE)
            IF (PSI .GE. PSEL  .AND.  PSI .LE. PSER) THEN
               POLD      = PTOLD(I)
               DP        = PSI - POLD
               IF (DABS(DP) .LT. TOLKSL) THEN
                  ETAI(I)= SW(I)*SNODI(I) +
     1                     PNODI(I)*HUSWR1*FHUDSE(PSI)
               ELSE
                  SEOLD  = FHUSE(POLD)
                  ETAI(I)= SW(I)*SNODI(I) +
     1                     PNODI(I)*HUSWR1*((SE - SEOLD)/DP)
               END IF
            ELSE
               ETAI(I)   = SW(I)*SNODI(I) +
     1                     PNODI(I)*HUSWR1*FHUDSE(PSI)
            END IF
         END DO
      ELSE
         DO I=1,N
            PSI          = PTNEW(I)
            SE           = FBCSE(PSI)
            SW(I)        = BCPORM(I)*SE + BCRMC/PNODI(I)
            CKRW(I)      = FBCKR(PSI)
            IF (PSI .GE. PSEL  .AND.  PSI .LE. PSER) THEN
               POLD      = PTOLD(I)
               DP        = PSI - POLD
               IF (DABS(DP) .LT. TOLKSL) THEN
                  ETAI(I)= SW(I)*SNODI(I) +
     1                     PNODI(I)*BCPORM(I)*FBCDSE(PSI)
               ELSE
                  SEOLD  = FBCSE(POLD)
                  ETAI(I)= SW(I)*SNODI(I) +
     1                     PNODI(I)*BCPORM(I)*((SE - SEOLD)/DP)
               END IF
            ELSE
               ETAI(I)   = SW(I)*SNODI(I) +
     1                     PNODI(I)*BCPORM(I)*FBCDSE(PSI)
            END IF
         END DO
      END IF
C
      RETURN
      END
