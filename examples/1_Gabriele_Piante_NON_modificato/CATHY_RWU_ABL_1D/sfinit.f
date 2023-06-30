C
C**************************  SFINIT ************************************
C
C  calculate seepage face exit points based on initial pressure heads
C
C***********************************************************************
C
      SUBROUTINE SFINIT(NSF,NSFNUM,NSFNOD,SFEX,SFEXP,SFEXIT,DUPUIT,
     1                  PTIMEP,PNEW,SFFLAG,DELTAT,Z)
C
      IMPLICIT NONE
      INCLUDE 'CATHY.H'
      INTEGER  I,J,K,INOD
      INTEGER  NSF,DUPUIT
      INTEGER  NSFNUM(*),NSFNOD(NSFMAX,*),SFEX(*),SFEXP(*),SFEXIT(*)
      INTEGER  SFFLAG(*)
      REAL*8   DELTAT
      REAL*8   PTIMEP(*),PNEW(*),Z(*)
      INCLUDE 'IOUNITS.H'
C
      DO I=1,NSF
         SFEX(I)=NSFNUM(I) + 1
         SFEXP(I)=SFEX(I)
         SFEXIT(I)=SFEX(I)
         DO J=1,NSFNUM(I)
            INOD=NSFNOD(I,J)
            IF (PTIMEP(INOD) .GE. 0.0D0) THEN
               SFEX(I)=J
               SFEXP(I)=SFEX(I)
               SFEXIT(I)=SFEX(I)
               DO K=J,NSFNUM(I)
                  INOD=NSFNOD(I,K)
                IF (DUPUIT.EQ.0) THEN
                    PTIMEP(INOD)=0.0D0
                    PNEW(INOD)=0.0D0
                ELSE
                    PTIMEP(INOD)=0.0D0+Z(NSFNOD(I,SFEX(I)))-Z(INOD)
                    PNEW(INOD)=0.0D0+Z(NSFNOD(I,SFEX(I)))-Z(INOD)
                END IF
C                   WRITE(666,*) INOD,PNEW(INOD)
               END DO
               GO TO 50
            END IF
         END DO
   50    CONTINUE
         IF (SFEX(I) .EQ. NSFNUM(I)+1) THEN
            SFFLAG(1)=SFFLAG(1) + 1
            WRITE(IOUT10,2100) I,0.0D0,DELTAT,0
         ELSE IF (SFEX(I) .EQ. 1) THEN
            SFFLAG(2)=SFFLAG(2) + 1
            WRITE(IOUT10,2200) I,0.0D0,DELTAT,0
         END IF
      END DO
C
      RETURN
 2100 FORMAT(  ' SFFLAG(1) AT SEEPAGE FACE ',I4,
     1         ' (TIME=',1PE8.2,', DELTAT=',1PE8.2,', ITER=',I4,')',
     2       /,11X,'NO EXIT POINT; SEEPAGE FACE COMPLETELY',
     3         ' UNSATURATED')
 2200 FORMAT(  ' SFFLAG(2) AT SEEPAGE FACE ',I4,
     1         ' (TIME=',1PE8.2,', DELTAT=',1PE8.2,', ITER=',I4,')',
     2       /,11X,'EXIT POINT AT TOP OF SEEPAGE FACE; SEEPAGE FACE',
     3         ' COMPLETELY SATURATED')
      END
