C
C**************************  PONDUPD ***********************************
C
C  update surface BCs for FLOW3D to reflect ponding situation as
C  calculated in SURF_ROUTE, and update PONDING flag
C
C***********************************************************************
C
      SUBROUTINE PONDUPD(NNOD,PONDING,IFATM,PONDH_MIN,DELTAT,PONDNOD,
     1                   ARENOD,ATMPOT,ATMACT,PNEW)
C
      IMPLICIT NONE
      INTEGER  NNOD 
      INTEGER  I
      INTEGER  IFATM(NNOD)
      LOGICAL  PONDING
      REAL*8   PONDH_MIN,DELTAT
      REAL*8   ZERO,ONE,DTR
      REAL*8   PONDNOD(NNOD),ARENOD(NNOD),ATMPOT(NNOD),ATMACT(NNOD)
      REAL*8   PNEW(NNOD)
      PARAMETER (ZERO = 0.0D+00, ONE = 1.0D+00)
      INCLUDE  'IOUNITS.H'
C    
      PONDING = .FALSE.
      DTR = ONE/DELTAT
      DO I = 1,NNOD
         IF(IFATM(I) .EQ. -1) GO TO 500
         IF (PONDNOD(I) .GE. PONDH_MIN) THEN
            IF (IFATM(I) .EQ. 1 .OR. IFATM(I) .EQ. 2) THEN
               PNEW(I) = PONDNOD(I)
               PONDING = .TRUE.
            ELSE IF (IFATM(I) .EQ. 0) THEN
               PONDING = .TRUE.
              ATMACT(I) = ATMPOT(I)+ PONDNOD(I)*ARENOD(I)*DTR
            END IF
         END IF

 500     CONTINUE
      END DO
C
      RETURN  
      END
