C
C**************************  SHLPIC ************************************
C
C  restore diagonal elements of LHS system matrix corresponding to
C  Dirichlet nodes. 
C  Also, set solution at Dirichlet nodes to the prescribed 
C  values. This is done since the solved solution at
C  Dirichlet nodes may not be exactly equal to the presribed
C  values, due to inaccuracies and roundoff errors which could 
C  arise from the way we imposed Dirichlet conditions (by multiplying
C  diagonal terms by a 'large' number).
C  Picard scheme (symmetric storage).
C
C***********************************************************************
C
      SUBROUTINE SHLPIC(NP,TOPOL,CONTP,PRESC,LHSP,COEF1,PNEW,POLD,
     1                  NSF,NSFNUM,NSFNOD,SFEX,LHSSF,DUPUIT,Z,
     2                  NNOD,IFATM,LHSATM)
C
      IMPLICIT  NONE
      INCLUDE  'CATHY.H'
      INTEGER   I,J,IND,INOD
      INTEGER   NP,NSF,NNOD,DUPUIT
      INTEGER   TOPOL(*),CONTP(*),NSFNUM(*),NSFNOD(NSFMAX,*),SFEX(*)
      INTEGER   IFATM(*)
      REAL*8    PRESC(*),LHSP(*),COEF1(*),PNEW(*),POLD(*),Z(*)
      REAL*8    LHSSF(NSFMAX,*),LHSATM(*)
C
      DO I=1,NP
         J=CONTP(I)
         PNEW(J)=PRESC(I)
         IND=TOPOL(J)
         COEF1(IND)=LHSP(I)
      END DO
      DO I=1,NSF
         DO J=SFEX(I),NSFNUM(I)
            INOD=NSFNOD(I,J)
            IND=TOPOL(INOD)
          IF (DUPUIT.EQ.0) THEN
            PNEW(INOD)=0.0D0
          ELSE
            PNEW(INOD)=0.0D0+Z(NSFNOD(I,SFEX(I)))-Z(INOD)
          END IF
C           WRITE(66666,*)INOD,PNEW(INOD)
            COEF1(IND)=LHSSF(I,J)
         END DO
      END DO
      DO I=1,NNOD
         IF (IFATM(I) .EQ. 1 .OR. IFATM(I) .EQ. 2) THEN
            IND=TOPOL(I)
            PNEW(I)=POLD(I)
            COEF1(IND)=LHSATM(I)
         END IF
      END DO
C
C     WRITE(66666,*)
      RETURN
      END
