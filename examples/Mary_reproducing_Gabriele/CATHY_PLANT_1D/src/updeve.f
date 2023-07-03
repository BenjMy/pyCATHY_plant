C
C**************************  UPDEVE ************************************
C
C  analysis step of the EnKF algorithm
C
C***********************************************************************
C
      SUBROUTINE UPDEVE(CTR,N,NCELNL,NUMRES,NENS,NOBS,IVGHU,DAFLAG,
     1                  NZONE,NSTR,DSMEAS,ENKFNOD,ENKFVAL,DELTA_X,
     2                  PNODI,ENPOROS,ENRETC,ENKSX,ENKSY,ENKSZ,ENELSTOR,
     3                  NROW,NCOL,ENDTM_KSS1_SF_1,ENDTM_KSS1_SF_2,
     4                  ENDTM_WS1_SF_1,ENDTM_WS1_SF_2,
     5                  ENPNEW,ENQ_IN_KKP1_SN,ENQ_OUT_KKP1_SN_1,
     6                  ENQ_OUT_KKP1_SN_2,ENVOLUME_KKP1_SN,
     7                  ENH_POOL_KKP1_VEC,QOI_SN,ENPT)
C
      IMPLICIT NONE
      INCLUDE 'CATHY.H'
      INCLUDE 'SOILCHAR.H'
      INCLUDE 'IOUNITS.H'
      LOGICAL VERBOSE,UPDATE_RANDROT
      INTEGER I,J,K,L,DAFLAG,IVGHU,MODE,IND0
      INTEGER CTR,N,NCELNL,NENS,NOBS,NDIM,NUMRES,NSTR,NZONE,NROW,NCOL
      INTEGER QOI_SN(MAXCEL),ENPT(MAXNENS)
      INTEGER ENKFNOD(MAXNUDN,2)
      REAL*8  GASDEV,DSMEAS,TRUNCATION,DELTA_X
      REAL*8  PSI,POR,SE,FVGSE,FHUSE,FBCSE
      REAL*8  NPER(MAXNENS),DENPER(MAXNENS)
      REAL*8  PNODI(NMAX)
c      REAL*8  ENPNODI(NMAX,MAXNENS),ENSNODI(NMAX,MAXNENS)
      REAL*8  ENPOROS(MAXSTR,MAXZON,MAXNENS)
      REAL*8  ENELSTOR(MAXSTR,MAXZON,MAXNENS)
      REAL*8  ENRETC(3,MAXNENS),ENKSX(MAXSTR,MAXZON,MAXNENS)
      REAL*8  ENKSY(MAXSTR,MAXZON,MAXNENS),ENKSZ(MAXSTR,MAXZON,MAXNENS)
      REAL*8  ENDTM_KSS1_SF_1(COLMAX,ROWMAX,MAXNENS)
      REAL*8  ENDTM_KSS1_SF_2(COLMAX,ROWMAX,MAXNENS)
      REAL*8  ENDTM_WS1_SF_1(COLMAX,ROWMAX,MAXNENS)
      REAL*8  ENDTM_WS1_SF_2(COLMAX,ROWMAX,MAXNENS)
      REAL*8  ENPNEW(NMAX,MAXNENS),ENQ_IN_KKP1_SN(MAXCEL,MAXNENS)
      REAL*8  ENQ_OUT_KKP1_SN_1(MAXCEL,MAXNENS)
      REAL*8  ENQ_OUT_KKP1_SN_2(MAXCEL,MAXNENS)
      REAL*8  ENVOLUME_KKP1_SN(MAXCEL,MAXNENS)
      REAL*8  ENH_POOL_KKP1_VEC(MAXRES,MAXNENS)
      REAL*8  ENKFVAL(MAXNUDT,MAXNUDN)
      REAL*8  ONEN(NENS,NENS),IONEN(NENS,NENS)
      REAL*8, ALLOCATABLE :: D(:),INNOV(:)
      REAL*8, ALLOCATABLE :: HA(:,:),S(:,:),E(:,:),DD(:,:)
      REAL*8, ALLOCATABLE :: A(:,:),R(:,:),HAMEAN(:,:)
CM    REAL*8  HA(MAXNUDN,MAXNENS)
CM    REAL*8  A(MAXNDIM,MAXNENS),E(MAXNUDN,MAXNENS)
CM    REAL*8  S(MAXNUDN,MAXNENS)
CM    REAL*8  DD(MAXNUDN,MAXNENS)
      INCLUDE 'RANDOM.H'
      PARAMETER (MODE=23)
      PARAMETER (TRUNCATION=0.999d0)
cm    print *,CTR,NDIM,NENS,NOBS,DSMEAS 
C
C  Allocate and initialize various arrays
C
      IF (DAFLAG.EQ.1) NDIM=N+4*NCELNL+NUMRES
cp old update parameters:
cp      IF (DAFLAG.GT.1) NDIM=3*(1+N+NSTR*NZONE)+4*(NCELNL+NROW*NCOL)
cp new update parameters:
      IF (DAFLAG.GT.1) NDIM=N+3+5*NSTR*NZONE+4*(NCELNL+NROW*NCOL)
CM uncomment when lakes will be available again +NUMRES
      ALLOCATE (A(NDIM,NENS),E(NOBS,NENS),S(NOBS,NENS))
      ALLOCATE (DD(NOBS,NENS),R(NOBS,NOBS),D(NOBS),INNOV(NOBS))
      ALLOCATE (HA(NOBS,NENS),HAMEAN(NOBS,NENS))
      DO I=1,NOBS
         D(I)=0.0d0
         INNOV(I)=0.0d0
         DO J=1,NENS
            HA(I,J)=0.0d0
            HAMEAN(I,J)=0.0d0
            S(I,J)=0.0d0
            E(I,J)=0.0d0
            DD(I,J)=0.0d0
         END DO
         DO J=1,NOBS
            R(I,J)=0.0d0
         END DO
      END DO
C
C  I-1N and 1N matrix assemblage
C
      DO I=1,NENS
         DO J=1,NENS
            IF (I.EQ.J) THEN
               IONEN(I,J)=1.0d0-1.0D0/NENS
            ELSE
               IONEN(I,J)=-1.0D0/NENS
            END IF
            ONEN(I,J)=1.0D0/NENS
         END DO
      END DO
C
C  Assemble matrix HA for EnKF algorithm
C
      DO I=1,NOBS
         DO J=1,NENS
            IF (ENKFNOD(I,2).EQ.0) THEN
               PSI=ENPNEW(ENKFNOD(I,1),ENPT(J))
               IF (IVGHU .EQ. 0) THEN
                  SE = FVGSE(PSI)
               ELSE IF (IVGHU .EQ. 2) THEN
                  SE = FHUSE(PSI)
               ELSE IF (IVGHU .EQ. 4) THEN
                  SE = FBCSE(PSI)
               ELSE
                  write(iout2,*) 'Attenzione, correggere IVGHU!'
                  stop
               END IF
               HA(I,J)=SE
            ELSE IF (ENKFNOD(I,2).EQ.1) THEN
               HA(I,J)=ENPNEW(ENKFNOD(I,1),ENPT(J))/PSISTAR
            ELSE IF (ENKFNOD(I,2).EQ.2) THEN
               IF (ENQ_OUT_KKP1_SN_1(ENKFNOD(I,1),ENPT(J)).NE.0.d0.AND.
     &             ENQ_OUT_KKP1_SN_2(ENKFNOD(I,1),ENPT(J)).NE.0.d0) THEN
                  WRITE(IOUT2,*)'Warning, discharge along 2 directions!'
                  STOP
               ELSE
                  HA(I,J)=(ENQ_OUT_KKP1_SN_1(ENKFNOD(I,1),ENPT(J))+
     +                     ENQ_OUT_KKP1_SN_2(ENKFNOD(I,1),ENPT(J)))/
     /                     QSTAR
               END IF
            END IF
         END DO
cm       write(90,'(235e15.6)')(H(I,J),J=1,NDIM+NOBS)
      END DO
C
C  Measurement assignment and perturbations
C
      DO I=1,NOBS
         IF (ENKFNOD(I,2).EQ.0) THEN
            POR=PNODI(ENKFNOD(I,1))
            IF (IVGHU .EQ. 0) THEN
               D(I)=(ENKFVAL(CTR,I)-VGRMC)/(POR-VGRMC)
            ELSE IF (IVGHU .EQ. 2) THEN
               D(I)=(ENKFVAL(CTR,I)/POR-HUSWR)/(1.0d0-HUSWR)
            ELSE IF (IVGHU .EQ. 4) THEN
               D(I)=(ENKFVAL(CTR,I)-BCRMC)/(POR-BCRMC)
            ELSE
               write(iout2,*) 'Attenzione, correggere IVGHU!'
               stop
            END IF
         ELSE IF (ENKFNOD(I,2).EQ.1) THEN
            D(I)=ENKFVAL(CTR,I)/PSISTAR
         ELSE IF (ENKFNOD(I,2).EQ.2) THEN
            D(I)=ENKFVAL(CTR,I)/QSTAR
         END IF
      END DO
      DO I=1,NOBS
         DO J=1,NENS
            NPER(J)=GASDEV(ISEED)
            IF (ENKFNOD(I,2).EQ.1) THEN
CM                E(I,J)=NPER(J)*DSMEAS*D(I)
cp                DD(I,J)=D(I)+NPER(J)*DSMEAS*D(I)
                DD(I,J)=D(I)+NPER(J)*DSMEAS
            ELSE 
                DENPER(J)=NPER(J)*
     1                    (DLOG(1.0D0+DSMEAS**2))**0.5d0
     2        	           -0.5d0*DLOG(1.0D0+DSMEAS**2)
                DD(I,J)=D(I)*DEXP(DENPER(J))
                IF (ENKFNOD(I,2).EQ.0) THEN
                   IF (DD(I,J).GT.1.0d0) THEN
                      DD(I,J)=1.0d0
                   ELSE IF (DD(I,J).LT.0.0d0) THEN
                      DD(I,J)=0.0d0
                   END IF
                END IF
            END IF
            E(I,J)=DD(I,J)-D(I)
         END DO
      END DO
C
C  Computation of the covariance matrix R
C
      CALL DGEMM('N','T',NOBS,NOBS,NENS,1.0D0/DFLOAT(NENS-1),E,NOBS,
     1            E,NOBS,0.0d0,R,NOBS)
cm       print *,(dd(i,j),j=1,NENS)
cm       write(97,'(234e15.6)')(H(i,j),j=1,ndim)
cm       write(91,'(10e15.6)')(E(i,j),j=1,NENS)
cm       write(92,'(10e15.6)')(DD(i,j),j=1,NENS)
C
C  Initialization and assemblage of the system state matrix
C
      DO I=1,NDIM
         DO J=1,NENS
            A(I,J)=0.0d0
         END DO
      END DO
      DO J=1,NENS
         DO I=1,N
            A(I,J)=ENPNEW(I,ENPT(J))/PSISTAR
         END DO
         DO I=1,NCELNL
            A(N+I,J)=ENQ_IN_KKP1_SN(QOI_SN(I),ENPT(J))/QSTAR
            A(N+NCELNL+I,J)=ENQ_OUT_KKP1_SN_1(QOI_SN(I),ENPT(J))/QSTAR
            A(N+2*NCELNL+I,J)=ENQ_OUT_KKP1_SN_2(QOI_SN(I),ENPT(J))/QSTAR
            A(N+3*NCELNL+I,J)=ENVOLUME_KKP1_SN(QOI_SN(I),ENPT(J))/
     /                        (PSISTAR*DELTA_X**2.0d0)
         END DO
CM       DO I=1,NUMRES
CM          A(I+N+4*NCELNL,J)=ENH_POOL_KKP1_VEC(I,J)/PSISTAR
CM       END DO
C
C If we want to update also parameters, we need an augmented state
C
         IF (DAFLAG.GT.1) THEN
            A(N+4*NCELNL+1,J)=DLOG(ENRETC(1,ENPT(J))-1.0d0)
            A(N+4*NCELNL+2,J)=DLOG(ENRETC(2,ENPT(J)))
            A(N+4*NCELNL+3,J)=DLOG(DABS(ENRETC(3,ENPT(J))))
            L=0
            DO I=1,NSTR
               DO K=1,NZONE
                  L=L+1
                  A(N+4*NCELNL+3+L,J)=DLOG(ENKSX(I,K,ENPT(J)))
                  A(N+4*NCELNL+3+NSTR*NZONE+L,J)=
     1                                DLOG(ENKSY(I,K,ENPT(J)))
                  A(N+4*NCELNL+3+2*NSTR*NZONE+L,J)=
     1                                DLOG(ENKSZ(I,K,ENPT(J)))
cp new state augmentation
                  A(N+4*NCELNL+3+3*NSTR*NZONE+L,J)=
     1                                DLOG(ENPOROS(I,K,ENPT(J)))
                  A(N+4*NCELNL+3+4*NSTR*NZONE+L,J)=
     1                                DLOG(ENELSTOR(I,K,ENPT(J)))
               END DO
            END DO
cp old state augmentation
cp            DO I=1,N
cp               A(N+4*NCELNL+3+3*NSTR*NZONE+I,J)=DLOG(ENPNODI(I,ENPT(J)))
cp               A(N+4*NCELNL+3+3*NSTR*NZONE+N+I,J)=
cp     1                                          DLOG(ENSNODI(I,ENPT(J)))
cp            END DO
            L=0
cp            IND0=N+4*NCELNL+3+3*NSTR*NZONE+N+N
            IND0=N+4*NCELNL+3+5*NSTR*NZONE
            DO I=1,NCOL
               DO K=1,NROW
                  L=L+1
                  A(IND0+L,J)=DLOG(ENDTM_KSS1_SF_1(I,K,ENPT(J)))
                  A(IND0+NCOL*NROW+L,J)=
     1                        DLOG(ENDTM_KSS1_SF_2(I,K,ENPT(J)))
                  A(IND0+2*NCOL*NROW+L,J)=
     1                        DLOG(ENDTM_WS1_SF_1(I,K,ENPT(J)))
                  A(IND0+3*NCOL*NROW+L,J)=
     1                        DLOG(ENDTM_WS1_SF_2(I,K,ENPT(J)))
               END DO
            END DO
         END IF
      END DO      
C
C  Computation of matrix S=HA'
C
      CALL DGEMM('N','N',NOBS,NENS,NENS,1.0D0,HA,NOBS,
     1            IONEN,NENS,0.0D0,S,NOBS)
C
C  Computation of the mean innovation matrix
C
      CALL DGEMM('N','N',NOBS,NENS,NENS,1.0D0,HA,NOBS,ONEN,
     1            NENS,0.0D0,HAMEAN,NOBS)
      write(iout2,*) '--------------------------------------------'
      write(iout2,*) 'MEAN INNOVATION'
      DO I=1,NOBS
         INNOV(I)=D(I)-HAMEAN(I,1)
      END DO
CM    write(iout2,'(10e15.5)') (D(I),I=1,NOBS)
      VERBOSE=.TRUE.
      UPDATE_RANDROT=.TRUE.
      CALL ANALYSIS(A, R, E, S, DD, INNOV, NDIM, NENS, NOBS,
     &              VERBOSE, TRUNCATION, MODE, UPDATE_RANDROT)
CM    call analysis(A, E, S, D, NDIM, NENS, NOBS, VERBOSE)
C
C  Update of the pressure and discharge ensembles
C
      DO J=1,NENS
         DO I=1,N
            ENPNEW(I,ENPT(J))=A(I,J)*PSISTAR
         END DO
         DO I=1,NCELNL
            IF (A(I+N,J).ge.0.0d0) THEN
               ENQ_IN_KKP1_SN(QOI_SN(I),ENPT(J))=A(I+N,J)*QSTAR
            ELSE
               ENQ_IN_KKP1_SN(QOI_SN(I),ENPT(J))=0.0d0
            END IF
            IF (A(I+N+NCELNL,J).ge.0.0d0) THEN
               ENQ_OUT_KKP1_SN_1(QOI_SN(I),ENPT(J))=A(I+N+NCELNL,J)*
     1                                        QSTAR
            ELSE
               ENQ_OUT_KKP1_SN_1(QOI_SN(I),ENPT(J))=0.0d0
            END IF
            IF (A(I+N+2*NCELNL,J).ge.0.0d0) THEN
               ENQ_OUT_KKP1_SN_2(QOI_SN(I),ENPT(J))=A(I+N+2*NCELNL,J)*
     1                                        QSTAR
            ELSE
               ENQ_OUT_KKP1_SN_2(QOI_SN(I),ENPT(J))=0.0d0
            END IF
            IF (A(I+N+3*NCELNL,J).ge.0.0d0) THEN
               ENVOLUME_KKP1_SN(QOI_SN(I),ENPT(J))=A(I+N+3*NCELNL,J)*
     1                                       (PSISTAR*DELTA_X**2.0d0)
            ELSE
               ENVOLUME_KKP1_SN(QOI_SN(I),ENPT(J))=0.0d0
            END IF
         END DO
CM       DO I=1,NUMRES
CM          ENH_POOL_KKP1_VEC(I,J)=A(I+N+3*NCELNL,J)*PSISTAR
CM       END DO
         IF (DAFLAG.GT.1) THEN
            ENRETC(1,ENPT(J))=DEXP(A(N+4*NCELNL+1,J))+1.0d0
            ENRETC(2,ENPT(J))=DEXP(A(N+4*NCELNL+2,J))
            ENRETC(3,ENPT(J))=-1.d0*DEXP(A(N+4*NCELNL+3,J))
            L=0
            DO I=1,NSTR
               DO K=1,NZONE
                  L=L+1
                  ENKSX(I,K,ENPT(J))=DEXP(A(N+4*NCELNL+3+L,J))
                  ENKSY(I,K,ENPT(J))=
     1                    DEXP(A(N+4*NCELNL+3+NSTR*NZONE+L,J))
                  ENKSZ(I,K,ENPT(J))=
     1                    DEXP(A(N+4*NCELNL+3+2*NSTR*NZONE+L,J))
                  ENPOROS(I,K,ENPT(J))=
     1                    DEXP(A(N+4*NCELNL+3+3*NSTR*NZONE+L,J))
                  ENELSTOR(I,K,ENPT(J))=
     1                    DEXP(A(N+4*NCELNL+3+4*NSTR*NZONE+L,J))
               END DO
            END DO
cp            DO I=1,N
cp               ENPNODI(I,ENPT(J))=DEXP(A(N+4*NCELNL+3+3*NSTR*NZONE+I,J))
cp               ENSNODI(I,ENPT(J))=
cp     1                          DEXP(A(N+4*NCELNL+3+3*NSTR*NZONE+N+I,J))
cp            END DO
            L=0
cp            IND0=N+4*NCELNL+3+3*NSTR*NZONE+N+N
            IND0=N+4*NCELNL+3+5*NSTR*NZONE
            DO I=1,NCOL
               DO K=1,NROW
                  L=L+1
                  ENDTM_KSS1_SF_1(I,K,ENPT(J))=DEXP(A(IND0+L,J))
                  ENDTM_KSS1_SF_2(I,K,ENPT(J))=
     1                               DEXP(A(IND0+NCOL*NROW+L,J))
                  ENDTM_WS1_SF_1(I,K,ENPT(J))=
     1                               DEXP(A(IND0+2*NCOL*NROW+L,J))
                  ENDTM_WS1_SF_2(I,K,ENPT(J))=
     1                               DEXP(A(IND0+3*NCOL*NROW+L,J))
               END DO
            END DO
         END IF
      END DO      
C
C  Update counter value
C
      write(333,*)ctr
      write(111,*)ctr
      do i=1,nstr
         do j=1,nzone
            write(333,333)i,j,(enksx(i,j,k),k=1,nens)
            write(333,333)i,j,(enksy(i,j,k),k=1,nens)
            write(333,333)i,j,(enksz(i,j,k),k=1,nens)
            write(111,333)i,j,(enelstor(i,j,k),k=1,nens)
         end do
      end do
      CTR=CTR+1
C
      DEALLOCATE(A,E,S)
      DEALLOCATE(DD,R,D,INNOV)
      DEALLOCATE(HA,HAMEAN)
      RETURN
333   format(2(1i4),100(1pe15.6))
      END
