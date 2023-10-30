C
C************************** ENTIMUPD ************************************
C
C  update TIME and time step size for the next time step for all the 
C  the realizations of the ensemble  (ENKF=TRUE; EN_TOT=TRUE)
C
C***********************************************************************
C
      SUBROUTINE ENTIMUPD(TIMEP,TIME,ENDELTAT,DTMIN,DTMAX,TIMESTOP,
     1                  DTMAGA,DTMAGM,DTREDS,DTREDM,DTGMIN,DELTAT0,
     2                  NPRT,KPRT,TIMPRT,NOBS,ENKFCTR,ENKFT,ENKFTIM,
     3                  UPD,NENS,DELTAT,TMAX,ENPT)
C
      IMPLICIT  NONE
      INCLUDE   'CATHY.H'
      INTEGER   NPRT,KPRT,NOBS,ENKFCTR,ENKFT,NENS,J
      INTEGER   ENPT(MAXNENS)
      LOGICAL   DTGMIN,UPD
      REAL*8    TIMEP,TIME,DTMIN,DTMAX,TMAX,DELTAT,TIMESTOP
      REAL*8    DTMAGA,DTMAGM,DTREDS,DTREDM,DELTAT0
      REAL*8    TIMPRT(MAXPRT), ENKFTIM(MAXNUDT)
      REAL*8    ENDELTAT(MAXNENS)
C
C Calculate next TIMESTOP
C     
cp      write(111,*)' '
cp      write(111,*)'TIME,TIMEP,DELTAT,TIMESTOP'
cp      write(111,*)TIME,TIMEP,DELTAT,TIMESTOP
 
      TIMESTOP=TMAX
      IF (NPRT .GT. 0 .AND. KPRT .LE. NPRT) THEN
         IF (TIMESTOP- TIMPRT(KPRT) .GT. DTMAX) THEN
            TIMESTOP=TIMPRT(KPRT)
         END IF
      END IF
cp      write(111,*)'TIMESTOP'
cp      write(111,*)TIMESTOP
      IF ((NOBS .NE.0).AND.(ENKFCTR .LE. ENKFT)) THEN
         IF ( TIMESTOP.GE.ENKFTIM(ENKFCTR)) THEN
            TIMESTOP=ENKFTIM(ENKFCTR)
         ELSE IF (ENKFTIM(ENKFCTR)- TIMESTOP .LT. DTMAX  ) THEN
            ENKFTIM(ENKFCTR)=TIMESTOP
         END IF
      END IF
cp      write(111,*)'TIMESTOP'
cp      write(111,*)TIMESTOP
      DO J=1,NENS
         IF (UPD) THEN
            DELTAT = DELTAT0
         ELSE 
         DELTAT=ENDELTAT(ENPT(J))
         IF (DELTAT.GT.DTMAX) DELTAT=DTMAX
cp         DO J=2,NENS
cp            IF (ENDELTATP(ENPT(J)).GT.DELTAT) DELTAT=ENDELTATP(ENPT(J))
cp         END DO
         END IF
cp      WRITE(111,*) 'DELTAT'
cp      WRITE(111,*) DELTAT
         IF ((TIME+DELTAT) .GT. TIMESTOP  .OR.
     1              DABS(TIMESTOP-(TIME+DELTAT)).LT.
     2               0.5D0*(DELTAT*DTREDM-DTREDS)) THEN
                DELTAT=TIMESTOP - TIME
         END IF
cp      WRITE(111,*) DELTAT
         IF (DELTAT .LE. DTMIN) THEN
            DTGMIN=.FALSE.
         ELSE  
            DTGMIN=.TRUE.
         END IF
         ENDELTAT(ENPT(J))=DELTAT
         TIMEP = TIME
      END DO
cp      WRITE(111,*)'TIME'
cp      WRITE(111,*)TIME
      UPD = .FALSE.
C
      RETURN
      END
