C
C**************************  RHSPIC ************************************
C
C  assemble RHS vector, without contribution of the unsaturated
C  zone gravitational term and without the boundary conditions:
C  Picard scheme
C
C***********************************************************************
C
      SUBROUTINE RHSPIC(N,JA,TOPOL,DELTAT,PTNEW,PNEW,
     1                  PTIMEP,COEF1,COEF2,TNOTI)
C
      IMPLICIT  NONE
      INTEGER   I,J,K,M
      INTEGER   N
      INTEGER   JA(*),TOPOL(*)
      REAL*8    RDT,COEF
      REAL*8    DELTAT
      REAL*8    PTNEW(*),PNEW(*),PTIMEP(*),COEF1(*),COEF2(*),TNOTI(*)
C
      CALL INIT0R(N,TNOTI)
      RDT=1.0D0/DELTAT
      DO K=1,N
         I=TOPOL(K)
         J=TOPOL(K+1)-1
         DO M=I,J
            COEF=COEF2(M)*RDT 
            TNOTI(K)=TNOTI(K) - COEF1(M)*PTNEW(JA(M)) - 
     1                          COEF*(PNEW(JA(M)) - PTIMEP(JA(M)))
            IF (M .NE. I) TNOTI(JA(M))=TNOTI(JA(M)) - 
     1                          COEF1(M)*PTNEW(K) - 
     2                          COEF*(PNEW(K) - PTIMEP(K))
         END DO
      END DO
C
      RETURN
      END
