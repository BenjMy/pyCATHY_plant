C
C************************ PERT_IC **************************************
C
C  computes the ensemble of the initial conditions and store the values
C  of the variables PNEW, POLD, PTNEW, PTOLD for each realization
C
C***********************************************************************
C
      SUBROUTINE PERT_IC(NENS,N,TETAF,PTIMEP,DSIC,ENPTIMEP,ENPNEW,
     1                   ENPOLD,ENPTNEW,ENPTOLD)
C
      IMPLICIT NONE
      INCLUDE 'CATHY.H'
      INTEGER I,J
      INTEGER NENS,N
cm    INTEGER COUNT,COUNT_RATE,COUNT_MAX
      REAL*8  GASDEV,DSIC,TETAF,AVEDENPER,VARDENPER
      REAL*8  PTIMEP(NMAX)
      REAL*8  TIMEP(NMAX),NEW(NMAX),TNEW(NMAX)
      REAL*8  NPER(MAXNENS),DENPER(NENS)
      REAL*8  ENPTIMEP(NMAX,NENS)
      REAL*8  ENPNEW(NMAX,NENS),ENPOLD(NMAX,NENS)
      REAL*8  ENPTNEW(NMAX,NENS),ENPTOLD(NMAX,NENS)
      INCLUDE 'RANDOM.H'
      INCLUDE 'IOUNITS.H'
C
C Computation of the initial condition ensemble
C (gaussian errors with mean 0 and standard deviation as
C an absolute pressure head value [L])
C
C     CALL SYSTEM_CLOCK(COUNT,COUNT_RATE,COUNT_MAX)
C     ISEED=-1*COUNT
cm    ISEED=0
cm    CALL RNSET(ISEED)
cm    print *,ISEED
      AVEDENPER=0.0
      WRITE(IOUT56,*) ISEED,'ISEED pert_ic'
      WRITE(IOUT56,*)
      WRITE(IOUT56,*)'Nominal perturbations on the initial conditions'
      WRITE(IOUT56,101)'mean value','stand dev'
      WRITE(IOUT56,102)AVEDENPER,DSIC
      WRITE(IOUT56,*)'Ensemble perturbations'
      WRITE(IOUT56,103) 'NRE','PERTURB'
      DO J=1,NENS
         NPER(J)=GASDEV(ISEED)
         DENPER(J)=NPER(J)*DSIC
         AVEDENPER=AVEDENPER+DENPER(J)
         WRITE(IOUT56,104)J,DENPER(J)
         DO I=1,N
            ENPTIMEP(I,J)=PTIMEP(I)+DENPER(J)
            ENPNEW(I,J)=ENPTIMEP(I,J)
            ENPOLD(I,J)=ENPNEW(I,J)
         END DO
      END DO
      AVEDENPER=AVEDENPER/NENS
      VARDENPER=0.0d0
      DO J=1,NENS
          VARDENPER=VARDENPER+(AVEDENPER-DENPER(J))**2
      END DO
      IF (NENS.EQ.1) THEN
          VARDENPER=0.0d0
      ELSE
          VARDENPER=VARDENPER/FLOAT(NENS-1)
      END IF
      WRITE(IOUT56,105)'mean', AVEDENPER
      WRITE(IOUT56,105)'var', VARDENPER
      WRITE(IOUT56,* )' '
      DO I=1,NENS
         DO J=1,N
            NEW(J)=ENPNEW(J,I)
            TIMEP(J)=ENPTIMEP(J,I)
         END DO
         CALL WEIGHT(N,TETAF,NEW,TIMEP,TNEW)
         DO J=1,N
            ENPTNEW(J,I)=TNEW(J)
            ENPTOLD(J,I)=ENPTNEW(J,I)
         END DO
      END DO
C
cm    print *,ISEED
cm    STOP
  101 format(2(2x,A15))
  102 format(2(2x,f15.8))
  103 format(2x,A7,2x,A15)
  104 format(2x,I7,2x,f15.8)
  105 format(2x,A7,2x,f15.8)
      RETURN
      END
