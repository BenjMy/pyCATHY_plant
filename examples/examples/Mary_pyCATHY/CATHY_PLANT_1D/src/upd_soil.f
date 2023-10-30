C
C************************ UPD_SOIL **********************************
C
C  computes the mean and cv of the the soil parameters of the
c  ensemble (Van Genuchten model,
C  saturated hydraulic conductivity, elastic storage and porosity)
C  after a SIR/ENKF with parameters updates. 
C  Output in file 'ensemble'
C
C***********************************************************************
C
      SUBROUTINE UPD_SOIL(NENS,N,NSTR,NZONE,NROW,NCOL,PERMX,PERMY,
     1                    PERMZ,ELSTOR,POROS,
     2                    ENRETC,ENKSX,ENKSY,ENKSZ,ENELSTOR,
     3                    ENPOROS,ENDTM_KSS1_SF_1,
     4                    ENDTM_WS1_SF_1,
     5                    ENDSRETC,ENDSKS,ENDSSTOR,
     6                    ENDSPOROS,ENDSSURF_KS,ENDSSURF_WS,ENPT,DAFLAG)
      IMPLICIT NONE
      INCLUDE 'CATHY.H'
      INCLUDE 'SOILCHAR.H'
      INCLUDE 'RIVERNETWORK.H'
      INCLUDE 'IOUNITS.H'

      INTEGER I,J,K,L,DAFLAG
      INTEGER NENS,N,NSTR,NZONE,NROW,NCOL,ENPT(NENS) 

      REAL*8 GASDEV
      REAL*8 PERMX(MAXSTR,NZONE),PERMY(MAXSTR,NZONE)
      REAL*8 PERMZ(MAXSTR,NZONE)
      REAL*8 POROS(MAXSTR,NZONE),ELSTOR(MAXSTR,NZONE)
      REAL*8 ENRETC(3,NENS)
      REAL*8 ENDSRETC(3)
      REAL*8 ENKSX(MAXSTR,MAXZON,MAXNENS)
      REAL*8 ENKSY(MAXSTR,MAXZON,MAXNENS)
      REAL*8 ENKSZ(MAXSTR,MAXZON,MAXNENS),ENDSKS(MAXSTR,MAXZON,3)
      REAL*8 ENPOROS(MAXSTR,MAXZON,MAXNENS),ENDSPOROS(MAXSTR,MAXZON)
      REAL*8 ENELSTOR(MAXSTR,MAXZON,MAXNENS),ENDSSTOR(MAXSTR,MAXZON)
      REAL*8 ENDSSURF_KS,ENDSSURF_WS
      REAL*8 ENDTM_KSS1_SF_1(COLMAX,ROWMAX,MAXNENS)
      REAL*8 ENDTM_WS1_SF_1(COLMAX,ROWMAX,MAXNENS)
      REAL*8 AVE_ENRETC(3),VAR_ENRETC(3)
      REAL*8 AVE_ENKSX(MAXSTR,NZONE),AVE_ENKSY(MAXSTR,NZONE)
      REAL*8 AVE_ENKSZ(MAXSTR,NZONE)
      REAL*8 AVE_DENPER,VAR_DENPER
      REAL*8 VAR_ENKSX(MAXSTR,NZONE),VAR_ENKSY(MAXSTR,NZONE)
      REAL*8 VAR_ENKSZ(MAXSTR,NZONE)
      REAL*8 AVE_ELSTOR(MAXSTR,NZONE),AVE_POROS(MAXSTR,NZONE)
      REAL*8 VAR_ELSTOR(MAXSTR,NZONE),VAR_POROS(MAXSTR,NZONE)
      
C 
C Computation of the Van Genuchten model parameters:
C mean and cv of the ensemble 
C
cm    print *,'Inserire un numero intero'
cm    read(5,*)ISEED
cm    CALL RNSET(ISEED)
cm    CALL RNGET(ISEED)
cm    ISEED=-1*ISEED
      write(IOUT56,*)'--- PARAMETERS UPDATE AFTER THE RESAMPLE ---'
      IF (DAFLAG .EQ.4)THEN
         write(IOUT56,*)'Nominal retention curve parameters and cv'
         write(IOUT56,103)' ','VGN','VGRMC','VGPSAT'
         write(IOUT56,105)'MEAN',VGN,VGRMC,VGPSAT
         write(IOUT56,105)'CV',ENDSRETC(1),ENDSRETC(2),ENDSRETC(3)
         write(IOUT56,*)' ' 
      END IF
      DO I=1,3
         AVE_ENRETC(i)=0
      END DO
      DO J=1,NENS
         AVE_ENRETC(1)=AVE_ENRETC(1)+ENRETC(1,ENPT(J))
         AVE_ENRETC(2)=AVE_ENRETC(2)+ENRETC(2,ENPT(J))
         AVE_ENRETC(3)=AVE_ENRETC(3)+ENRETC(3,ENPT(J))
      END DO
      DO J=1,3
         AVE_ENRETC(J)=AVE_ENRETC(J)/NENS
         VAR_ENRETC(J)=0.0d0
         DO I=1,NENS
            VAR_ENRETC(J)=VAR_ENRETC(J)+(ENRETC(J,ENPT(I))-
     1                                   AVE_ENRETC(J))**2
         END DO
         IF (NENS.GT.1) THEN
            VAR_ENRETC(J)=VAR_ENRETC(J)/FLOAT(NENS-1)
            ENDSRETC(J)=abs(sqrt(VAR_ENRETC(J))/AVE_ENRETC(J)) 
         ELSE
            VAR_ENRETC(J)=0.0d0
            ENDSRETC(J)=0.0d0
         END IF
      END DO
      VGN=AVE_ENRETC(1)
      VGRMC=AVE_ENRETC(2)
      VGPSAT=AVE_ENRETC(3)
      write(iout2,*)'Ensemble of retention curve parameters'
      do i=1,3
         write(iout2,10)(ENRETC(I,ENPT(J)),J=1,NENS)
      end do
      write(IOUT56,*)'Ensemble of retention curve parameters'
      write(IOUT56,103)'NRE','VGN','VGRMC','VGPSAT'
      do j=1,nens
         write(IOUT56,104) ENPT(J), (ENRETC(I,ENPT(J)),I=1,3)
      end do
      write(IOUT56,105) 'mean', (AVE_ENRETC(I),I=1,3)
      write(IOUT56,105) 'var', (VAR_ENRETC(I),I=1,3)
      write(IOUT56,105) 'cv', (ENDSRETC(I),I=1,3)

C
C Computation of the saturated hydraulic conductivity, porosity and
C elastic storage of the ensemble (only 
C zone-based heterogeneity is considered at this stage, plus anisotropy)
C
      DO J=1,NZONE
         DO K=1,NSTR
            AVE_ENKSX(K,J)=0.0d0
            AVE_ENKSY(K,J)=0.0d0
            AVE_ENKSZ(K,J)=0.0d0
            AVE_ELSTOR(K,J)=0.0d0
            AVE_POROS(K,J)=0.0d0
            VAR_ENKSX(K,J)=0.0d0
            VAR_ENKSY(K,J)=0.0d0
            VAR_ENKSZ(K,J)=0.0d0
            VAR_ELSTOR(K,J)=0.0d0
            VAR_POROS(K,J)=0.0d0
         END DO
c         write(iout2,*)'Ens. of hydraulic conductivity, elastic storage'
c         write(iout2,*)'and porosity perturbations for zone ',J
c         write(iout2,103)'NRE','KS','ELSTOR','POROS'
         DO I=1,NENS
            DO K=1,NSTR
               AVE_ENKSX(K,J)=AVE_ENKSX(K,J)+ENKSX(K,J,ENPT(I))
               AVE_ENKSY(K,J)=AVE_ENKSY(K,J)+ENKSY(K,J,ENPT(I))
               AVE_ENKSZ(K,J)=AVE_ENKSZ(K,J)+ENKSZ(K,J,ENPT(I))
               AVE_ELSTOR(K,J)=AVE_ELSTOR(K,J)+ENELSTOR(K,J,ENPT(I))
               AVE_POROS(K,J)=AVE_POROS(K,J)+ENPOROS(K,J,ENPT(I))
            END DO
         END DO
         DO K=1,NSTR
            AVE_ENKSX(K,J)=AVE_ENKSX(K,J)/FLOAT(NENS)
            AVE_ENKSY(K,J)=AVE_ENKSY(K,J)/FLOAT(NENS)
            AVE_ENKSZ(K,J)=AVE_ENKSZ(K,J)/FLOAT(NENS)
            AVE_ELSTOR(K,J)=AVE_ELSTOR(K,J)/FLOAT(NENS)
            AVE_POROS(K,J)=AVE_POROS(K,J)/FLOAT(NENS)
            DO I =1,NENS
               VAR_ENKSX(K,J)=VAR_ENKSX(K,J)+(ENKSX(K,J,ENPT(I))
     1     -AVE_ENKSX(K,J))**2     
               VAR_ENKSY(K,J)=VAR_ENKSY(K,J)+(ENKSY(K,J,ENPT(I))
     1     -AVE_ENKSY(K,J))**2     
               VAR_ENKSZ(K,J)=VAR_ENKSZ(K,J)+(ENKSZ(K,J,ENPT(I))
     1     -AVE_ENKSZ(K,J))**2     
               VAR_ELSTOR(K,J)=VAR_ELSTOR(K,J)+(ENELSTOR(K,J,ENPT(I))
     1     -AVE_ELSTOR(K,J))**2     
               VAR_POROS(K,J)=VAR_POROS(K,J)+(ENPOROS(K,J,ENPT(I))
     1     -AVE_POROS(K,J))**2     
            END DO
            IF (NENS.eq.1) then
               VAR_ENKSX(K,J)=0.0
               VAR_ENKSY(K,J)=0.0
               VAR_ENKSZ(K,J)=0.0
               VAR_ELSTOR(K,J)=0.0
               VAR_POROS(K,J)=0.0
            ELSE
               VAR_ENKSX(K,J)=VAR_ENKSX(K,J)/(NENS-1)
               VAR_ENKSY(K,J)=VAR_ENKSY(K,J)/(NENS-1)  
               VAR_ENKSZ(K,J)=VAR_ENKSZ(K,J)/(NENS-1) 
               VAR_ELSTOR(K,J)=VAR_ELSTOR(K,J)/(NENS-1)
               VAR_POROS(K,J)=VAR_POROS(K,J)/(NENS-1)
            END IF
         END DO
      END DO

      write(IOUT56,*) ' '
      write(IOUT56,*)'Soil param. (hydr. conductivity, elstor, poros)'
c      do K=1,NSTR
         K=1 
         do j=1,NZONE
            write(IOUT56,*) ' '
          IF (DAFLAG.EQ.4) THEN
            write(IOUT56,*) 'Input values and cv, str = ',K,' zone =', J
          write(IOUT56,108)' ', 'PERMX','PERMY','PERMZ','ELSTOR','POROS'
            write(IOUT56,110) 'MEAN', PERMX(K,J),PERMY(K,J),PERMZ(K,J),
     1          ELSTOR(K,J),POROS(K,J)
            write(IOUT56,110) 'CV', ENDSKS(K,J,1),ENDSKS(K,J,2),
     1          ENDSKS(K,J,3),ENDSSTOR(K,J),ENDSPOROS(K,J)
           END IF
           write(IOUT56,*)'Ensemble values, str=',K,' zone =',J
            write(IOUT56,108)'NRE','PERMX','PERMY','PERMZ'
            do L=1,nens
             I=ENPT(L)
             write(IOUT56,109) I,ENKSX(K,J,I),ENKSY(K,J,I),ENKSZ(K,J,I),
     1              ENELSTOR(K,J,I),ENPOROS(K,J,I)
            end do
            write(IOUT56,110) 'mean',AVE_ENKSX(K,J),AVE_ENKSY(K,J),
     1         AVE_ENKSZ(K,J),AVE_ELSTOR(K,J),AVE_POROS(K,J)
            write(IOUT56,110) 'var',VAR_ENKSX(K,J),VAR_ENKSY(K,J),
     1         VAR_ENKSZ(K,J),VAR_ELSTOR(K,J),VAR_POROS(K,J)
            write(IOUT56,110) 'cv',sqrt(VAR_ENKSX(K,J))/AVE_ENKSX(K,J),
     1                       sqrt(VAR_ENKSY(K,J))/AVE_ENKSY(K,J),
     2                       sqrt(VAR_ENKSZ(K,J))/AVE_ENKSZ(K,J),
     3                       SQRT(VAR_ELSTOR(K,J))/AVE_ELSTOR(K,J),
     4                       SQRT(VAR_POROS(K,J))/AVE_POROS(K,J)
         end do
c      END DO
      DO K=1,NSTR
         DO J=1,NZONE
            PERMX(K,J)=AVE_ENKSX(K,J)
            PERMY(K,J)=AVE_ENKSY(K,J)
            PERMZ(K,J)=AVE_ENKSZ(K,J)
            ELSTOR(K,J)=AVE_ELSTOR(K,J)
            POROS(K,J)=AVE_POROS(K,J)
            ENDSKS(K,J,1)=abs(SQRT(VAR_ENKSX(K,J))/AVE_ENKSX(K,J))
            ENDSKS(K,J,2)=ABS(SQRT(VAR_ENKSY(K,J))/AVE_ENKSY(K,J))
            ENDSKS(K,J,3)=ABS(SQRT(VAR_ENKSZ(K,J))/AVE_ENKSZ(K,J))
            ENDSSTOR(K,J)=ABS(SQRT(VAR_ELSTOR(K,J))/AVE_ELSTOR(K,J))
            ENDSPOROS(K,J)=ABS(SQRT(VAR_POROS(K,J))/AVE_POROS(K,J))
         END DO
      END DO
C
C
C Computation of the Gauckler-Strickler coefficient ensemble
C
      AVE_DENPER=0.0d0
      DO I=1,NENS
         AVE_DENPER=AVE_DENPER+ENDTM_KSS1_SF_1(1,1,ENPT(I))
      END DO
      write(iout2,*)'Ensemble of GS coefficient perturbations'
      write(IOUT56,*) ' '
      IF (DAFLAG.EQ.4) THEN
        write(IOUT56,*)'Nominal distribution of GS coeff. perturbations'
        write(IOUT56,101)'MEAN ', 1.0d0
        write(IOUT56,101)'CV  ', ENDSSURF_KS
      END IF
      AVE_DENPER=AVE_DENPER/NENS
      VAR_DENPER=0.0d0
      DO I=1,NENS
      VAR_DENPER=VAR_DENPER+(ENDTM_KSS1_SF_1(1,1,ENPT(I))-AVE_DENPER)**2
      END DO
      AVE_DENPER=AVE_DENPER/DTM_KSS1_SF_1(1,1)
      VAR_DENPER=VAR_DENPER/(DTM_KSS1_SF_1(1,1)**2)
       DO J=1,NCOL
          DO K=1,NROW
             DTM_KSS1_SF_1(J,K)=DTM_KSS1_SF_1(J,K)*AVE_DENPER
             DTM_KSS1_SF_2(J,K)=DTM_KSS1_SF_2(J,K)*AVE_DENPER
          END DO
      END DO
      IF (NENS.EQ.1) THEN
         VAR_DENPER=0.0d0
         ENDSSURF_KS=0.0d0
      ELSE
         VAR_DENPER=VAR_DENPER/(NENS-1)
         ENDSSURF_KS=abs(sqrt(VAR_DENPER)/AVE_DENPER)
      END IF
      write(IOUT56,*)'Ensemble distrib of GS coeff. perturbations'
      write(IOUT56,101)'MEAN ', AVE_DENPER
      write(IOUT56,101)'VAR  ', VAR_DENPER
      write(IOUT56,101)'CV   ', ENDSSURF_KS
C
C Computation of the water surface width coefficient ensemble
C
      AVE_DENPER=0.0d0
      DO I=1,NENS
         AVE_DENPER=AVE_DENPER+ENDTM_WS1_SF_1(1,1,ENPT(I))
      END DO
      write(IOUT56,*) ' '
      IF (DAFLAG.EQ.4) THEN
        write(IOUT56,*)'Nominal distribution of WS coeff. perturbations'
        write(IOUT56,101)'MEAN ', 1.0d0
        write(IOUT56,101)'CV  ', ENDSSURF_WS
      END IF
      AVE_DENPER=AVE_DENPER/NENS
      VAR_DENPER=0.0d0
      DO I=1,NENS
       VAR_DENPER=VAR_DENPER+(ENDTM_WS1_SF_1(1,1,ENPT(I))-AVE_DENPER)**2
      END DO
      AVE_DENPER=AVE_DENPER/DTM_WS1_SF_1(1,1)
      VAR_DENPER=VAR_DENPER/(DTM_WS1_SF_1(1,1)**2)
      DO J=1,NCOL
         DO K=1,NROW
            DTM_WS1_SF_1(J,K)=DTM_WS1_SF_1(J,K)*AVE_DENPER
            DTM_WS1_SF_2(J,K)=DTM_WS1_SF_2(J,K)*AVE_DENPER
         END DO
      END DO
      IF (NENS.EQ.1) THEN
         VAR_DENPER=0.0d0
         ENDSSURF_WS=0.0d0
      ELSE
         VAR_DENPER=VAR_DENPER/(NENS-1)
         ENDSSURF_WS=ABS(sqrt(VAR_DENPER)/AVE_DENPER)
      END IF
      write(IOUT56,*)'Ensemble distrib of WS coefficient perturbations'
      write(IOUT56,101)'MEAN ', AVE_DENPER
      write(IOUT56,101)'VAR  ', VAR_DENPER
      write(IOUT56,101)'CV   ', ENDSSURF_WS
cm    print *,ISEED
   10 format(10(1PE15.6))
  101 format(2x,A7,2x,f10.7)
c  102 format(4(2x,f12.6))
  103 format(2x,A7,3(2x,A12))
  104 format(2x,I7,3(2x,f12.6))
  105 format(2x,A7,3(2x,f12.6))
c  106 format(4(2x,A12))
c  107 format(4(2x,e12.6))
  108 format(2x,A7,5(2x,A12))
  109 format(2x,I7,5(2x,e12.6))
  110 format(2x,A7,5(2x,e12.6))
      RETURN
      END
