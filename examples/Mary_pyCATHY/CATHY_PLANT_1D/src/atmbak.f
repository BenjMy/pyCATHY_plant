C
C**************************  ATMBAK ************************************
C
C  special handling of atmospheric boundary conditions for
C  the case where ATMTIM(1) < TIME <= ATMTIM(2) during back-stepping
C
C***********************************************************************
C
      SUBROUTINE ATMBAK(NNOD,TIME,IFATM,AREA,ATMPOT,ATMACT,
     1                  ATMTIM,ATMINP,IETO,DELTAT,
     2                  ENKF,DSATM,QNEW,QTIMEP)
C
      IMPLICIT  NONE
      INCLUDE   'CATHY.H'
      INTEGER   I,J,IETO
      INTEGER   NNOD
      INTEGER   IFATM(NODMAX)
      LOGICAL   ENKF
      REAL*8    SLOPE,DSATM,ALPHA
      REAL*8    TIME,GASDEV,DELTAT
      REAL*8    AREA(*),ATMPOT(*),ATMACT(*),ATMTIM(*),ATMINP(3,*)
      REAL*8    NPER,DENPER
      REAL*8    QNEW,QTIMEP
cp      REAL*8    ENATMPOT(NODMAX,MAXNENS),ENATMACT(NODMAX,MAXNENS)
cp      REAL*8    ENATMOLD(NODMAX,MAXNENS)
      include   'RANDOM.H'
C
C  we don't need to do anything if first time step or if atmospheric
C  inputs are homogeneous in time
C
      IF (ATMTIM(1) .GE. ATMTIM(2)) GO TO 800
      IF (ENKF) THEN
         ALPHA=1.0d0-DELTAT/ATMTAU
c        ALPHA=0.0d0
         NPER=GASDEV(ISEED)
         QNEW=ALPHA*QTIMEP+
     &                 NPER*(1.0d0-ALPHA**2.0d0)**0.5d0
         DENPER=QNEW*(DLOG(1.0D0+DSATM**2))**0.5d0
     1                   -0.5d0*DLOG(1.0D0+DSATM**2)
      END IF
      DO I=1,NNOD
         SLOPE=(ATMINP(2,I) - ATMINP(1,I))/(ATMTIM(2) - ATMTIM(1))
         IF (IETO.NE.0) SLOPE=0.0d0
         ATMPOT(I)=(ATMINP(1,I) + SLOPE*(TIME - ATMTIM(1)))*AREA(I)
         IF (IFATM(I) .EQ. 0) ATMACT(I)=ATMPOT(I)
         IF (ENKF) THEN
            ATMPOT(I)=ATMPOT(I)*DEXP(DENPER)
            IF (IFATM(I) .EQ. 0) ATMACT(I)=ATMPOT(I)
         END IF
      END DO
C
  800 RETURN
      END
