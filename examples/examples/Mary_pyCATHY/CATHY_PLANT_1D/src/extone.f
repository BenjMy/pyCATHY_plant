C
C**************************  EXTONE ************************************
C
C  calculate new position of the exit point along each seepage face,
C  and adjust boundary conditions for the seepage face nodes to reflect
C  changes in the position of the exit point.
C  Check for convergence of seepage face exit points. 
C  Case ISFONE=1 : seepage face exit point updating performed by 
C                  checking only the one node above and one node
C                  below the current exit point position
C
C***********************************************************************
C
      SUBROUTINE EXTONE(NSF,NSFNUM,NSFNOD,SFEX,SFEXIT,SFQ,PNEW,
     1                  SFFLAG,ITER,TIME,DELTAT,KSF,DUPUIT,Z)
C
      IMPLICIT  NONE
      INCLUDE  'CATHY.H'
      INTEGER   I,J,K,KEX,INOD,DUPUIT
      INTEGER   NSF,ITER,KSF
      INTEGER   NSFNUM(*),NSFNOD(NSFMAX,*),SFEX(*),SFEXIT(*),SFFLAG(*)
      REAL*8    TIME,DELTAT
      REAL*8    SFQ(NSFMAX,*),PNEW(*),Z(*)
      INCLUDE  'IOUNITS.H'
C
      DO 100 I=1,NSF
         KEX=SFEX(I)
         IF (KEX .LE. NSFNUM(I)) THEN
            J=KEX
            IF (SFQ(I,J) .GE. 0.0D0) THEN
               SFEX(I)=J + 1
               IF (J+1 .GT. NSFNUM(I)) THEN
                  SFEX(I)=NSFNUM(I) + 1
                  SFFLAG(1)=SFFLAG(1) + 1
                  WRITE(IOUT10,2100) I,TIME,DELTAT,ITER
               END IF
               SFQ(I,J)=0.0D0
               DO K=KEX-1,1,-1
                  INOD=NSFNOD(I,K)
                  IF (PNEW(INOD) .GE. 0.0D0) THEN
                     SFFLAG(3)=SFFLAG(3) + 1
                     WRITE(IOUT10,2300) I,TIME,DELTAT,ITER,PNEW(INOD),
     1                                  K,INOD
                  END IF
               END DO
               GO TO 100
            END IF
         END IF
         IF (KEX .GT. 1) THEN
            J=KEX-1
            INOD=NSFNOD(I,J)
            IF (PNEW(INOD) .GE. 0.0D0) THEN
               SFEX(I)=J
               PNEW(INOD)=0.0D0
               DO K=J-1,1,-1
                  INOD=NSFNOD(I,K)
                  IF (PNEW(INOD) .GE. 0.0D0) THEN
                     SFFLAG(4)=SFFLAG(4) + 1
                     WRITE(IOUT10,2400) I,TIME,DELTAT,ITER,PNEW(INOD),
     1                                  K,INOD
                  END IF
               END DO
               IF (J .EQ. 1) THEN
                  SFFLAG(2)=SFFLAG(2) + 1
                  WRITE(IOUT10,2200) I,TIME,DELTAT,ITER
               END IF
               GO TO 100
            END IF
         END IF
  100 CONTINUE
      IF (DUPUIT.EQ.1) THEN
        DO I=1,NSF
         DO J=SFEX(I),NSFNUM(I)
            INOD=NSFNOD(I,J)
            PNEW(INOD)=0.0D0+Z(NSFNOD(I,SFEX(I)))-Z(INOD)
         END DO
        END DO
      END IF
C
C  check for convergence of seepage face exit points 
C
      CALL EXTCVG(NSF,NSFNOD,SFEX,SFEXIT,TIME,DELTAT,ITER,KSF)
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
 2300 FORMAT(  ' SFFLAG(3) AT SEEPAGE FACE ',I4,
     1         ' (TIME=',1PE8.2,', DELTAT=',1PE8.2,', ITER=',I4,')',
     2       /,11X,'NON-NEGATIVE PRESSURE HEAD OF ',1PE9.3,' AT NODE ',
     3         I3,' (NODE # ',I6,')',
     4       /,11X,'OCCURRED AT A POTENTIAL SEEPAGE FACE NODE DURING',
     5         ' EXIT POINT LOWERING')
 2400 FORMAT(  ' SFFLAG(4) AT SEEPAGE FACE ',I4,
     1         ' (TIME=',1PE8.2,', DELTAT=',1PE8.2,', ITER=',I4,')',
     2       /,11X,'NON-NEGATIVE PRESSURE HEAD OF ',1PE9.3,' AT NODE ',
     3         I3,' (NODE # ',I6,')',
     4       /,11X,'OCCURRED AT A POTENTIAL SEEPAGE FACE NODE DURING',
     5         ' EXIT POINT RAISING')
      END
