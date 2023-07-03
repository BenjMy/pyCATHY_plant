C
C**************************  PICUNS ************************************
C
C  calculate characteristic curves for unsaturated zone,
C  Picard case
C
C***********************************************************************
C
      SUBROUTINE PICUNS(NLKP,N,NT,KSLOPE,IVGHU,NTRI,TP,TETRA,
     1                  PTNEW,PTOLD,SNODI,PNODI,SW,CKRW,ETAI,
     2                  SWE,CKRWE,ETAE)
C
      IMPLICIT NONE
      INTEGER  NLKP,N,NT,KSLOPE,IVGHU,NTRI
      INTEGER  TP(*),TETRA(5,*)
      REAL*8   PTNEW(*),PTOLD(*),SNODI(*),PNODI(*)
      REAL*8   SW(*),CKRW(*),ETAI(*),SWE(*),CKRWE(*),ETAE(*)
C
      IF (IVGHU .EQ. -1) THEN
         CALL MOISTABPIC(NLKP,NT,NTRI,TETRA,PTNEW,PNODI,SNODI,
     1                   SWE,CKRWE,ETAE)
         CALL ELTNOD(N,NT,TP,TETRA,ETAE,ETAI)
         CALL ELTNOD(N,NT,TP,TETRA,CKRWE,CKRW)
         CALL ELTNOD(N,NT,TP,TETRA,SWE,SW)
      ELSE
         IF (KSLOPE .EQ. 0) THEN
            CALL CHPIC0(N,PTNEW,SNODI,PNODI,SW,CKRW,ETAI,IVGHU)
         ELSE IF (KSLOPE .EQ. 1) THEN
            CALL CHPIC1(N,PTNEW,PTOLD,SNODI,PNODI,SW,CKRW,ETAI,IVGHU)
         ELSE IF (KSLOPE .EQ. 2) THEN
            CALL CHPIC2(N,PTNEW,PTOLD,SNODI,PNODI,SW,CKRW,ETAI,IVGHU)
         ELSE IF (KSLOPE .EQ. 3) THEN
            CALL CHPIC3(N,PTNEW,PTOLD,SNODI,PNODI,SW,CKRW,ETAI,IVGHU)
         ELSE
            CALL CHPIC4(N,PTNEW,SNODI,PNODI,SW,CKRW,ETAI,IVGHU)
         END IF
         CALL NODELT(NT,TETRA,CKRW,CKRWE)
         CALL NODELT(NT,TETRA,ETAI,ETAE)
      END IF
C
      RETURN
      END
