C
C**************************  RESAMPLE **********************************
C
C resample step of the SIR (new algorithm, after ??? )
C
C***********************************************************************
C
      SUBROUTINE RESAMPLE(NENS,N,NUMRES,WSIR,ENPNEW,ENQ_IN_KKP1_SN,
     1                    ENQ_OUT_KKP1_SN_1,ENQ_OUT_KKP1_SN_2,
     2                    ENH_POOL_KKP1_VEC,ENVOLUME_KKP1_SN,
     3                    ENPONDNOD,ELSTOR,POROS,DAFLAG,
     3                    PERMX,PERMY,PERMZ,
     3                    NNOD,ENPT,NCELNL,NZONE,NSTR,
     4                    ENRETC,ENKSX,ENKSY,ENKSZ,
     5                    NROW,NCOL,ENDTM_KSS1_SF_1,ENDTM_KSS1_SF_2,
     6                    ENDTM_WS1_SF_1,ENDTM_WS1_SF_2,
     7                    ENPOROS,ENELSTOR,ENDSRETC,ENDSKS,  
     8                    ENDSSTOR,ENDSPOROS,ENDSSURF_KS,ENDSSURF_WS) 

      IMPLICIT NONE
      INCLUDE 'CATHY.H'
      INCLUDE 'RANDOM.H'
      INCLUDE 'IOUNITS.H'
      INCLUDE 'SOILCHAR.H'
      INCLUDE 'RIVERNETWORK.H'

      INTEGER I,J,K,N,NUMRES,NENS,NNOD
      INTEGER ENPT(MAXNENS)
      INTEGER NCELNL,NZONE, NSTR,NROW,NCOL
      INTEGER NEWREP(NENS),CONTREP(NENS),NREP,CONT
      INTEGER SUCC,DAFLAG,ZON

      REAL*8 NPER,DENPER,GASDEV,DENPER1
      REAL*8 WSIR(MAXNENS),CUM(NENS),RAN1,U,U0
      REAL*8 ENPNEW(NMAX,NENS)
      REAL*8 ENQ_IN_KKP1_SN(MAXCEL,NENS)
      REAL*8 ENQ_OUT_KKP1_SN_1(MAXCEL,NENS)
      REAL*8 ENQ_OUT_KKP1_SN_2(MAXCEL,NENS)
      REAL*8 ENH_POOL_KKP1_VEC(MAXRES,NENS)
      REAL*8 ENVOLUME_KKP1_SN(MAXCEL,NENS)
      REAL*8 ENPONDNOD(NODMAX,NENS)
      REAL*8 ENPOROS(MAXSTR,MAXZON,MAXNENS)
      REAL*8 ENELSTOR(MAXSTR,MAXZON,MAXNENS)
      REAL*8 ENRETC(3,NENS),ENKSX(MAXSTR,MAXZON,MAXNENS)
      REAL*8 ENKSY(MAXSTR,MAXZON,MAXNENS),ENKSZ(MAXSTR,MAXZON,MAXNENS)
      REAL*8 ENDTM_KSS1_SF_1(COLMAX,ROWMAX,MAXNENS)
      REAL*8 ENDTM_KSS1_SF_2(COLMAX,ROWMAX,MAXNENS)
      REAL*8 ENDTM_WS1_SF_1(COLMAX,ROWMAX,MAXNENS)
      REAL*8 ENDTM_WS1_SF_2(COLMAX,ROWMAX,MAXNENS)
      REAL*8 PERMX(MAXSTR,NZONE),PERMY(MAXSTR,NZONE),PERMZ(MAXSTR,NZONE)
      REAL*8 ELSTOR(MAXSTR,NZONE),POROS(MAXSTR,NZONE) 
      REAL*8 ENDSRETC(3),OLENDSRETC(3)
      REAL*8 ENDSKS(MAXSTR,MAXZON,3),ENDSPOROS(MAXSTR,NZONE)
      REAL*8 OLENDSKS(MAXSTR,MAXZON,3),OLENDSPOROS(MAXSTR,NZONE)
      REAL*8 ENDSSTOR(MAXSTR,NZONE),OLENDSSTOR(MAXSTR,NZONE)
      REAL*8 ENDSSURF_KS,ENDSSURF_WS
      REAL*8 DSSTOR,DSKS,DSRETC,DSPOROS,DSSURF

      CUM(1)=WSIR(ENPT(1))
      DO I=2,NENS
         CUM(I)=CUM(I-1)+WSIR(ENPT(I))
      END DO
      I = 1
      write(IOUT56,*)'RESAMPLING STEP'
      U0 = RAN1(ISEED)
      write(IOUT56,*)'ISEED = ',ISEED,' U0 = ',U0 
      write(IOUT56,*)
      write(IOUT56,*)NENS,' REALIZATIONS (NENS)'
      U0 = U0*1.0d0/NENS
      NREP=0
      CONT=0
      DO J = 1,NENS
         U = U0 + (1.0d0/NENS)*(J-1)
         IF (U .LE. CUM(I)) THEN
            IF (J.EQ.1) THEN
               CONT=1
               NEWREP(CONT)=1
               CONTREP(CONT)=0
            END IF
            CONTREP(CONT)=CONTREP(CONT)+1
         ELSE
            CONT=CONT+1
            DO WHILE (U .GT. CUM(I))
               I = I+1
            END DO
            NEWREP(CONT)=I
            CONTREP(CONT)=1
         END IF
      END DO
      NREP=CONT
      WRITE(IOUT56,*)NREP,'NUMERO REALIZZAZIONI DUPLICATE'
      WRITE(IOUT56,*)'REALIZ, NUMERO COPIE'
      DO CONT=1,NREP
         WRITE(IOUT56,*)ENPT(NEWREP(CONT)),CONTREP(CONT)
      END DO
C
c parameters update
c
      IF (DAFLAG.EQ.4) THEN
         IF (NREP.EQ.1) THEN
c
c Only one realization is duplicated. The mean values of the ensemble 
c parameters correspond to the parameters associated to this
c realizations. The CVs are the same of the previous iteration.
c
            VGN= ENRETC(1,ENPT(NEWREP(1)))
            VGRMC= ENRETC(2,ENPT(NEWREP(1)))
            VGPSAT= ENRETC(3,ENPT(NEWREP(1)))
            DO K=1,NSTR
               DO J=1,NZONE
                  PERMX(K,J)=ENKSX(K,J,ENPT(NEWREP(1)))
                  PERMY(K,J)=ENKSY(K,J,ENPT(NEWREP(1)))
                  PERMZ(K,J)=ENKSZ(K,J,ENPT(NEWREP(1)))
                  ELSTOR(K,J)=ENELSTOR(K,J,ENPT(NEWREP(1)))
                  POROS(K,J)=ENPOROS(K,J,ENPT(NEWREP(1)))
               END DO
            END DO
            DO J=1,NCOL
               DO K=1,NROW
                DTM_KSS1_SF_1(J,K)=ENDTM_KSS1_SF_1(J,K,ENPT(NEWREP(1)))
                DTM_KSS1_SF_2(J,K)=ENDTM_KSS1_SF_2(J,K,ENPT(NEWREP(1)))
                DTM_WS1_SF_1(J,K)=ENDTM_WS1_SF_1(J,K,ENPT(NEWREP(1)))
                DTM_WS1_SF_2(J,K)=ENDTM_WS1_SF_2(J,K,ENPT(NEWREP(1)))
               END DO
            END DO
         ELSE
c
c If more than one realization is duplicated, the mean and the cv 
c of the ensemble parameters correspond to the weight mean and cv of 
c the parameters duplicated
c To avoid the ensemble degeneration, the cv will be the greatest 
c between the half of the old cv and the new cv.
c
            VGN= 0.0d0
            OLENDSRETC(1)=ENDSRETC(1)
            ENDSRETC(1)=0.0d0
            VGRMC= 0.0d0
            OLENDSRETC(2)=ENDSRETC(2) 
            ENDSRETC(2)=0.0d0
            VGPSAT= 0.0d0
            OLENDSRETC(3)=ENDSRETC(3) 
            ENDSRETC(3)=0.0d0
            DO K=1,NSTR
               DO J=1,NZONE
                  PERMX(K,J)=0.0d0
                  OLENDSKS(K,J,1)=ENDSKS(K,J,1)
                  ENDSKS(K,J,1)=0.0d0
                  PERMY(K,J)=0.0d0
                  OLENDSKS(K,J,2)=ENDSKS(K,J,2)
                  ENDSKS(K,J,2)=0.0d0
                  PERMZ(K,J)=0.0d0
                  OLENDSKS(K,J,3)=ENDSKS(K,J,3)
                  ENDSKS(K,J,3)=0.0d0
                  ELSTOR(K,J)=0.0d0
                  OLENDSSTOR(K,J)=ENDSSTOR(K,J)
                  ENDSSTOR(K,J)=0.0d0
                  POROS(K,J)=0.0d0
                  OLENDSPOROS(K,J)=ENDSPOROS(K,J)
                  ENDSPOROS(K,J)=0.0d0
               END DO
            END DO
            ENDSSURF_KS=0.0d0
            ENDSSURF_WS=0.0d0
            DO J=1,NCOL
               DO K=1,NROW
                 DTM_KSS1_SF_1(J,K)=0.0d0
                 DTM_KSS1_SF_2(J,K)=0.0d0
                 DTM_WS1_SF_1(J,K)=0.0d0
                 DTM_WS1_SF_2(J,K)=0.0d0
               END DO
            END DO
            DO I=1,NREP
              VGN= VGN+ENRETC(1,ENPT(NEWREP(I)))*CONTREP(I)/NENS
              VGRMC= VGRMC+ENRETC(2,ENPT(NEWREP(I)))*CONTREP(I)/NENS
              VGPSAT=VGPSAT+ ENRETC(3,ENPT(NEWREP(I)))*CONTREP(I)/NENS 
              DO K=1,NSTR
                 DO J=1,NZONE
                    PERMX(K,J)=PERMX(K,J)+ENKSX(K,J,ENPT(NEWREP(I)))*
     1                         CONTREP(I)/NENS
                    PERMY(K,J)=PERMY(K,J)+ENKSY(K,J,ENPT(NEWREP(I)))*
     1                         CONTREP(I)/NENS
                    PERMZ(K,J)=PERMZ(K,J)+ENKSZ(K,J,ENPT(NEWREP(I)))*
     1                         CONTREP(I)/NENS
                  ELSTOR(K,J)=ELSTOR(K,J)+ENELSTOR(K,J,ENPT(NEWREP(I)))*
     1                         CONTREP(I)/NENS
                    POROS(K,J)=POROS(K,J)+ENPOROS(K,J,ENPT(NEWREP(I)))*
     1                         CONTREP(I)/NENS
                 END DO
              END DO
              DO J=1,NCOL
                 DO K=1,NROW
                   DTM_KSS1_SF_1(J,K)=DTM_KSS1_SF_1(J,K)+
     1        ENDTM_KSS1_SF_1(J,K,ENPT(NEWREP(I)))*CONTREP(I)/NENS
                   DTM_KSS1_SF_2(J,K)=DTM_KSS1_SF_2(J,K)+
     1        ENDTM_KSS1_SF_2(J,K,ENPT(NEWREP(I)))*CONTREP(I)/NENS
                   DTM_WS1_SF_1(J,K)=DTM_WS1_SF_1(J,K)+
     1        ENDTM_WS1_SF_1(J,K,ENPT(NEWREP(I)))*CONTREP(I)/NENS
                   DTM_WS1_SF_2(J,K)=DTM_WS1_SF_2(J,K)+
     1        ENDTM_WS1_SF_2(J,K,ENPT(NEWREP(I)))*CONTREP(I)/NENS
                 END DO
              END DO
            END DO
            DO I=1,NREP
              ENDSRETC(1)=ENDSRETC(1)+(VGN-ENRETC(1,ENPT(NEWREP(I))))**2
     1                    *CONTREP(I)/(NENS-1)
              ENDSRETC(2)=ENDSRETC(2)+(VGRMC-ENRETC(2,ENPT(NEWREP(I))))
     1                    **2*CONTREP(I)/(NENS-1)
              ENDSRETC(3)=ENDSRETC(3)+(VGPSAT-ENRETC(3,ENPT(NEWREP(I))))
     1                    **2*CONTREP(I)/(NENS-1)
              DO K=1,NSTR
                 DO J=1,NZONE
                 ENDSKS(K,J,1)=ENDSKS(K,J,1)+(ENKSX(K,J,ENPT(NEWREP(I)))
     1                 -PERMX(K,J))**2*CONTREP(I)/(NENS-1)
                 ENDSKS(K,J,2)=ENDSKS(K,J,2)+(ENKSY(K,J,ENPT(NEWREP(I)))
     1                 -PERMY(K,J))**2*CONTREP(I)/(NENS-1)
                 ENDSKS(K,J,3)=ENDSKS(K,J,3)+(ENKSZ(K,J,ENPT(NEWREP(I)))
     1                 -PERMZ(K,J))**2*CONTREP(I)/(NENS-1)
             ENDSSTOR(K,J)=ENDSSTOR(K,J)+(ENELSTOR(K,J,ENPT(NEWREP(I)))
     1                 -ELSTOR(K,J))**2*CONTREP(I)/(NENS-1)
             ENDSPOROS(K,J)=ENDSPOROS(K,J)+(ENPOROS(K,J,ENPT(NEWREP(I)))
     1                 -POROS(K,J))**2*CONTREP(I)/(NENS-1)
                 END DO
              END DO
              ENDSSURF_KS=ENDSSURF_KS+(DTM_KSS1_SF_1(1,1)-
     1      ENDTM_KSS1_SF_1(1,1,ENPT(NEWREP(I))))**2*CONTREP(I)/(NENS-1)
              ENDSSURF_WS=ENDSSURF_WS+(DTM_WS1_SF_1(1,1)-
     1      ENDTM_WS1_SF_1(1,1,ENPT(NEWREP(I))))**2*CONTREP(I)/(NENS-1)
            END DO
            ENDSRETC(1)=max(abs(SQRT(ENDSRETC(1))/VGN),
     1                      OLENDSRETC(1)/2.0d0) 
            ENDSRETC(2)=max(abs(SQRT(ENDSRETC(2))/VGRMC),
     1                      OLENDSRETC(2)/2.0d0)  
            ENDSRETC(3)=max(abs(SQRT(ENDSRETC(3))/VGPSAT),
     1                      OLENDSRETC(3)/2.0d0)  
            DO K=1,NSTR
               DO J=1,NZONE
                 ENDSKS(K,J,1)=max(ABS(SQRT(ENDSKS(K,J,1))/PERMX(K,J)),
     1                         OLENDSKS(K,J,1)/2.0d0)
                 ENDSKS(K,J,2)=max(ABS(SQRT(ENDSKS(K,J,2))/PERMY(K,J)),
     1                         OLENDSKS(K,J,2)/2.0d0)
                 ENDSKS(K,J,3)=max(ABS(SQRT(ENDSKS(K,J,3))/PERMZ(K,J)),
     1                         OLENDSKS(K,J,3)/2.0d0)
                 ENDSSTOR(K,J)=max(ABS(SQRT(ENDSSTOR(K,J))/ELSTOR(K,J)),
     1                         OLENDSSTOR(K,J)/2.0d0)
                ENDSPOROS(K,J)=max(ABS(SQRT(ENDSPOROS(K,J))/POROS(K,J)),
     1                         OLENDSPOROS(K,J)/2.0d0)
               END DO
            END DO
            ENDSSURF_KS=ABS(SQRT(ENDSSURF_KS)/DTM_KSS1_SF_1(1,1))
            ENDSSURF_WS=ABS(SQRT(ENDSSURF_WS)/DTM_WS1_SF_1(1,1))
         END IF
      END IF
      CONT=1
      SUCC=1
      WRITE(IOUT56,*) 'OLD REP, NEW REP'
      DO J=1,NENS
         DO WHILE ((CONTREP(CONT) .EQ. 0).OR.(CONTREP(CONT)==1 .AND.
     1                                       (J.LE.NEWREP(CONT))))
             IF (CONT<NREP) THEN
                CONT=CONT+1
             ELSE
                CONTREP(CONT)=-1
             END IF
         END DO
         IF (J.NE.NEWREP(SUCC)) THEN
             I=NEWREP(CONT)
             CONTREP(CONT)=CONTREP(CONT)-1
             write(IOUT56,*)ENPT(J),ENPT(I)
         ELSE 
             I=J
             CONTREP(SUCC)=CONTREP(SUCC)-1
             write(IOUT56,*)ENPT(J),ENPT(I)
             IF (SUCC.LT.NREP) THEN
                SUCC=SUCC+1
             END IF
             GO TO 100
         END IF
         DO K=1,N
             ENPNEW(K,ENPT(J))=ENPNEW(K,ENPT(I))
         END DO
         DO K=1,MAXCEL
            ENQ_IN_KKP1_SN(K,ENPT(J))=ENQ_IN_KKP1_SN(K,ENPT(I))
            ENQ_OUT_KKP1_SN_1(K,ENPT(J))=ENQ_OUT_KKP1_SN_1(K,ENPT(I))
            ENQ_OUT_KKP1_SN_2(K,ENPT(J))=ENQ_OUT_KKP1_SN_2(K,ENPT(I))
            ENVOLUME_KKP1_SN(K,ENPT(J))=ENVOLUME_KKP1_SN(K,ENPT(I))
         END DO
         DO K=1,NNOD
            ENPONDNOD(K,ENPT(J))=ENPONDNOD(K,ENPT(I))
         END DO
         DO K=1,NUMRES
            ENH_POOL_KKP1_VEC(K,ENPT(J))=ENH_POOL_KKP1_VEC(K,ENPT(I))
         END DO
c
c if DAFLAG = 4 generate new soil parameters
c
         IF (DAFLAG.EQ.4) THEN
C 
C Perturbation of the Van Genuchten model parameters
C
            DSRETC=ENDSRETC(1) 
            NPER=GASDEV(ISEED)
            DENPER=NPER*(DLOG(1.0D0+DSRETC**2))**0.5
     1             -0.5*DLOG(1.0D0+DSRETC**2)
            ENRETC(1,ENPT(J))=VGN*DEXP(DENPER)
            DSRETC=ENDSRETC(2) 
            NPER=GASDEV(ISEED)
            DENPER=NPER*(DLOG(1.0D0+DSRETC**2))**0.5
     1             -0.5*DLOG(1.0D0+DSRETC**2)
            ENRETC(2,ENPT(J))=VGRMC*DEXP(DENPER)
            DSRETC=ENDSRETC(3) 
            NPER=GASDEV(ISEED)
            DENPER=NPER*(DLOG(1.0D0+DSRETC**2))**0.5
     1             -0.5*DLOG(1.0D0+DSRETC**2)
            ENRETC(3,ENPT(J))=VGPSAT*DEXP(DENPER)
C
C perturbation of the saturated hydraulic conductivity, porosity and
C elastic storage (only 
C zone-based heterogeneity is considered at this stage, plus anisotropy)
C            
            DO ZON=1,NZONE 
C hydraulic conductivity perturbation
               DSKS=ENDSKS(1,ZON,1)
               NPER=GASDEV(ISEED)
               DENPER=NPER*(DLOG(1.0D0+DSKS**2))**0.5
     1                -0.5*DLOG(1.0D0+DSKS**2)
               DSKS=ENDSKS(1,ZON,3)
               NPER=GASDEV(ISEED)
               DENPER1=NPER*(DLOG(1.0D0+DSKS**2))**0.5
     1                -0.5*DLOG(1.0D0+DSKS**2)
               DO K=1,NSTR
                  ENKSX(K,ZON,ENPT(J))=PERMX(K,ZON)*DEXP(DENPER)
                  ENKSY(K,ZON,ENPT(J))=PERMY(K,ZON)*DEXP(DENPER)
                  ENKSZ(K,ZON,ENPT(J))=PERMZ(K,ZON)*DEXP(DENPER1) 
               END DO
C elastic storage perturbation
               DSSTOR=ENDSSTOR(1,ZON)
               NPER=GASDEV(ISEED)
               DENPER=NPER*(DLOG(1.0D0+DSSTOR**2))**0.5
     2                 -0.5*DLOG(1.0D0+DSSTOR**2)
               DO K=1,NSTR
                  ENELSTOR(K,ZON,ENPT(J))=ELSTOR(K,ZON)*DEXP(DENPER)
               END DO
C porosity perturbation
               DSPOROS=ENDSPOROS(1,ZON)
               NPER=GASDEV(ISEED)
               DENPER=NPER*
     1                (DLOG(1.0D0+(DSPOROS)**2))**0.5
     2            -0.5*DLOG(1.0D0+(DSPOROS)**2)
               DO K=1,NSTR
                  ENPOROS(K,ZON,ENPT(J))=POROS(K,ZON)*DEXP(DENPER)
               END DO
            END DO
C Perturbation of the Gauckler-Strickler coefficient
            DSSURF=ENDSSURF_KS
            NPER=GASDEV(ISEED)
            DENPER=NPER*
     1             (DLOG(1.0D0+DSSURF**2))**0.5
     2              -0.5*DLOG(1.0D0+DSSURF**2)
            DO ZON=1,NCOL
               DO K=1,NROW
                  ENDTM_KSS1_SF_1(ZON,K,ENPT(J))=DTM_KSS1_SF_1(ZON,K)*
     1                  DEXP(DENPER)
                  ENDTM_KSS1_SF_2(ZON,K,ENPT(J))=DTM_KSS1_SF_2(ZON,K)*
     1                  DEXP(DENPER)
               END DO
            END DO
C Perturbation of the water surface width coefficient
            DSSURF=ENDSSURF_WS
            NPER=GASDEV(ISEED)
            DENPER=NPER*
     1             (DLOG(1.0D0+DSSURF**2))**0.5
     2              -0.5*DLOG(1.0D0+DSSURF**2)
            DO ZON=1,NCOL
               DO K=1,NROW
                  ENDTM_WS1_SF_1(ZON,K,ENPT(J))=DTM_WS1_SF_1(ZON,K)*
     1                  DEXP(DENPER)
                  ENDTM_WS1_SF_2(ZON,K,ENPT(J))=DTM_WS1_SF_2(ZON,K)*
     1                  DEXP(DENPER)
               END DO
            END DO
         END IF
 100     continue
         WSIR(ENPT(J)) = 1.0d0/DFLOAT(NENS) 
      END DO

      IF (DAFLAG .EQ. 4) THEN
c
c upd_soil generates an output with the new ensemble mean value 
c and CV of the soil parameters
c
         CALL UPD_SOIL(NENS,N,NSTR,NZONE,NROW,NCOL,PERMX,PERMY,
     1                 PERMZ,ELSTOR,POROS,
     2                 ENRETC,ENKSX,ENKSY,ENKSZ,ENELSTOR,
     3                 ENPOROS,ENDTM_KSS1_SF_1,
     4                 ENDTM_WS1_SF_1,
     5                 ENDSRETC,ENDSKS,ENDSSTOR,
     6                 ENDSPOROS,ENDSSURF_KS,ENDSSURF_WS,ENPT,DAFLAG)
      END IF
      RETURN
      END
