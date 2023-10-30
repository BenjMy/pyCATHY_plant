C
C**************************  STORMB ************************************
C
C  calculate volume of change in storage
C  between the current time level and the previous time level.
C  DSTORE > 0 for net increase in storage.
C
C***********************************************************************
C
      SUBROUTINE STORMB(PNEW,PTIMEP,N,DSTORE,VOLNOD,ETAIP,ETAI)
C
      IMPLICIT NONE
      INTEGER  K
      INTEGER  N
      REAL*8   DSTORE
      REAL*8   PNEW(*),PTIMEP(*),VOLNOD(*),ETAIP(*),ETAI(*)
C
      DSTORE=0.0D0
      DO K=1,N
         DSTORE=DSTORE + VOLNOD(K)*(ETAI(K) + ETAIP(K))*0.5D0*
     1                             (PNEW(K) - PTIMEP(K))
      END DO
C
      RETURN
      END
