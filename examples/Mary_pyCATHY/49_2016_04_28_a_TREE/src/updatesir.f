C
C**************************  UPDATESIR *********************************
C
C analysis step of the SIR (new algorithm, after ??? )
C
C***********************************************************************
C
      SUBROUTINE UPDATESIR(CTR,NENS,NOBS,IVGHU,DSMEAS,PNODI,
     1                     ENPNEW,ENKFNOD,ENKFVAL,ENQ_OUT_KKP1_SN_1,
     2                     ENQ_OUT_KKP1_SN_2,WSIR,RESAMP,NEFFMIN,ENPT,
     3                     DSMEASMAX,EN_ERT)
C
      IMPLICIT NONE
      INCLUDE 'CATHY.H'
      INCLUDE 'SOILCHAR.H'
      INCLUDE 'IOUNITS.H'
      LOGICAL RESAMP
      INTEGER I,J,IVGHU,K,CONT,CONTQ
      INTEGER CTR,NENS,NOBS,NEFFMIN
      INTEGER ENKFNOD(MAXNUDN,2),ENPT(MAXNENS)
      INTEGER CONT_ERT
      REAL*8  GASDEV,DSMEAS,DSMEAS_SIR,NEFF,DSMEASMAX
      REAL*8  SUMW,SUMPSI,SUMTHETA,SUMQ,D,SUMW1,SUM_ERT
      REAL*8  PSI,POR,SE,MUY,SIG2Y
      REAL*8  FVGSE,FHUSE,FBCSE
      REAL*8  PNODI(NMAX),NPER(NOBS)
      REAL*8  ENPNEW(NMAX,MAXNENS),WSIR(MAXNENS),WSIRN(MAXNENS)
      REAL*8  WSIRN1(MAXNENS)
      REAL*8  ENQ_OUT_KKP1_SN_1(MAXCEL,MAXNENS),ENQ_IN1
      REAL*8  ENQ_OUT_KKP1_SN_2(MAXCEL,MAXNENS)
c     REAL*8  ENQ_OUT_KKP1_SN_2(MAXCEL,MAXNENS)
c     REAL*8  ENH_POOL_KKP1_VEC(MAXRES,MAXNENS)
      REAL*8  ENKFVAL(MAXNUDT,NOBS),ENKFVAL1
      REAL*8  EN_ERT(MAXNUDN,MAXNENS)
cm    print *,CTR,NDIM,NENS,NOBS,DSMEAS 
C
C 
C
      DSMEAS_SIR=DSMEAS
      CONT=0
      write(IOUT56,*) 'UPDATE SIR NUMBER',CTR
      write(IOUT56,*)' ' 
 10   write(IOUT56,*)'DSMEAS_SIR=',DSMEAS_SIR
      CONTQ=0
      SUMW=0.0d0
      SUMW1=0.0d0
      write(IOUT56,101)'NRE','SUMPSI','SUMTHETA','SUMQ','SUM_ERT',
     1 'WSIRN'
      DO K=1,NENS
         J=ENPT(K)
c         WRITE(IOUT56,*)J,' REALIZATION'
         SUMPSI=0.0d0
         SUMTHETA=0.0d0
         SUMQ=0.0d0
         SUM_ERT=0.0d0
         CONT_ERT=0
         DO I=1,NOBS
            IF (ENKFNOD(I,2).EQ.1) THEN
c                SIG2Y=(DSMEAS_SIR*ENKFVAL(CTR,I))**2
                SIG2Y=DSMEAS_SIR**2
                SUMPSI=SUMPSI+(ENKFVAL(CTR,I)-
     &                ENPNEW(ENKFNOD(I,1),J))**2/
     &                SIG2Y
            ELSE IF (ENKFNOD(I,2).EQ.0) THEN
               PSI = ENPNEW(ENKFNOD(I,1),J)
               POR = PNODI(ENKFNOD(I,1))
cm
cm  Warning: the following assemblage stands only if no
cm  variability of the retention curve parameters is
cm  considered, i.e, only if DSRETC=0
cm
               IF (IVGHU .EQ. 0) THEN
                  SE = FVGSE(PSI)
                  D=(POR-VGRMC)*SE+VGRMC
               ELSE IF (IVGHU .EQ. 2) THEN
                  D=(ENKFVAL(CTR,I)/POR-HUSWR)/(1.0d0-HUSWR)
                  SE = FHUSE(PSI)
               ELSE IF (IVGHU .EQ. 4) THEN
                  D=(ENKFVAL(CTR,I)-BCRMC)/(POR-BCRMC)
                  SE = FBCSE(PSI)
               END IF
               SIG2Y=log(1.d0+(DSMEAS_SIR)**2)
               MUY=-0.5d0*SIG2Y
               SUMTHETA=SUMTHETA+(log(ENKFVAL(CTR,I)/D)-MUY)**2/SIG2Y
            ELSE IF (ENKFNOD(I,2).EQ.2) THEN
               IF (ENQ_OUT_KKP1_SN_1(ENKFNOD(I,1),J).NE.0.d0.AND.
     &             ENQ_OUT_KKP1_SN_2(ENKFNOD(I,1),J).NE.0.d0) THEN
                  WRITE(IOUT2,*)'Warning, discharge along 2 directions!'
                  STOP
               END IF
               IF(K.eq.1) write(IOUT2,*)'MEASURED STREAMFLOW= '
     1                                ,ENKFVAL(CTR,I)
               write(IOUT2,*)J,' ENQ_OUT=', 
     1                      ENQ_OUT_KKP1_SN_1(ENKFNOD(I,1),J)+
     +                      ENQ_OUT_KKP1_SN_2(ENKFNOD(I,1),J)
               IF (ENKFVAL(CTR,I) .LT. 1.0d-15) THEN
                  ENKFVAL1=ENKFVAL(CTR,I)+1.0d-15
               ELSE
                  ENKFVAL1=ENKFVAL(CTR,I)
               END IF
               IF (ENQ_OUT_KKP1_SN_1(ENKFNOD(I,1),J) +
     +             ENQ_OUT_KKP1_SN_2(ENKFNOD(I,1),J) .LT. 1.0d-15) THEN
                  ENQ_IN1=ENQ_OUT_KKP1_SN_1(ENKFNOD(I,1),J)+
     +                    ENQ_OUT_KKP1_SN_2(ENKFNOD(I,1),J)+1.0d-15
               ELSE
                  ENQ_IN1=ENQ_OUT_KKP1_SN_1(ENKFNOD(I,1),J)+
     +                    ENQ_OUT_KKP1_SN_2(ENKFNOD(I,1),J)
               END IF
               SIG2Y=log(1.d0+(DSMEAS_SIR)**2)
               MUY=-0.5d0*SIG2Y
               SUMQ=SUMQ+
     &                (log(ENKFVAL1/ENQ_IN1)-MUY)**2/SIG2Y
            ELSE IF (ENKFNOD(I,2).EQ.3) THEN

               CONT_ERT=CONT_ERT+1
               SIG2Y=(DSMEAS_SIR*ENKFVAL(CTR,I))**2

c MODIFICA - PROVA (Manoli)---------------------------------
c               ENKFVAL(CTR,I)=ENKFVAL(CTR,I)/1000
c               EN_ERT(CONT_ERT,J)=EN_ERT(CONT_ERT,J)/1000
c--------------------------------------------------------

               SUM_ERT=SUM_ERT+(ENKFVAL(CTR,I)-
     &                EN_ERT(CONT_ERT,J))**2/SIG2Y
            END IF
         END DO
         WSIRN(J)=WSIR(J)*exp(-0.5d0*(SUMPSI+SUMTHETA+SUMQ+SUM_ERT))
         WRITE(IOUT56,102)J,SUMPSI,SUMTHETA,SUMQ,SUM_ERT,WSIRN(J)
         SUMW=SUMW+WSIRN(J)
      END DO
      write(IOUT56,*)SUMW ,'SOMMA PESI (SUMW)'
      IF (SUMW .lt. 1.d-300) THEN
         IF (DSMEAS_SIR.lt.DSMEASMAX) THEN
            DSMEAS_SIR=DSMEAS_SIR*2
            CONT=CONT+1
            write(IOUT56,*) ' ' 
            write(IOUT56,*) 'Double DSMEAS ',CONT
            GO TO 10
         ELSE
            RESAMP=.FALSE.   
            WRITE(IOUT56,*)'SUMW = 0.0 , NO UPDATE'
            GO TO 20
         END IF
      END IF
      write(IOUT56,*)CONT,DSMEAS_SIR,'number of doubled DSMEAS and fina
     1l DSMEAS'
C
C     Normalize weights and calculate the effective number of realization NEFF
C 
      NEFF=0.0d0
      DO K=1,NENS
         I=ENPT(K)
         WSIR(I)=WSIRN(I)/SUMW
         NEFF=NEFF+WSIR(I)**2
      END DO
      NEFF=1.0d0/NEFF
      write(iout2,100)(WSIR(ENPT(J)),J=1,NENS)
      write(IOUT56,*)'SIR WEIGHTS AFTER UPDATE',CTR
      write(IOUT56,100)(WSIR(ENPT(J)),J=1,NENS)
c      write(IOUT56,*)
c      write(IOUT56,103)(ENPT(J),J=1,NENS)
      write(IOUT56,*)' ' 
      write(IOUT56,*)NEFF,'NEFF (estimate number of duplicates in the  
     1resample)'
      RESAMP=.FALSE.
      IF (NEFF.LE.NEFFMIN) RESAMP=.TRUE.
  20  CONTINUE  
      IF (RESAMP) THEN
          write(IOUT56,*)'RESAMPLING STEP '
      ELSE
          write(IOUT56,*)'NO RESAMPLING STEP '
      END IF
c
C  Update counter value
C
      CTR=CTR+1
C
      RETURN
100   FORMAT(5E15.7)
101   FORMAT(A4,5(2x,A13))
102   FORMAT(I4,5(2x,e13.5))
      END
