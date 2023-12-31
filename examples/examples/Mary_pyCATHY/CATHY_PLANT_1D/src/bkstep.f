C
C**************************  BKSTEP ************************************
C
C  reset variables & arrays for back-stepping
C
C***********************************************************************
C
      SUBROUTINE BKSTEP(IPRT1,N,NSTR,NDIR,NDIRC,NP,NQ,HSPATM,HTIATM,
     1                  IETO,HTIDIR,HTINEU,ITER,NITERT,KBACKT,KBACK,
     2                  NNOD,NSF,NUMRES,NCELNL,
     3                  CONTP,CONTQ,IFATM,IFATMP,NSFNUM,SFEX,SFEXIT,
     4                  SFEXP,SURF,PONDING,PONDP,DTGMIN,
     5                  TIME,TIMEP,DELTAT,DTMIN,DTREDM,DTREDS,TETAF,
     6                  OVFLNOD,OVFLP,PRESC,PTIM,PINP,Q,QTIM,QINP,
     7                  ATMTIM,ATMINP,ARENOD,
     8                  POLD,PTOLD,PNEW,PTNEW,PTIMEP,
     9                  ATMACT,ATMPOT,ATMOLD,
     A                  QPNEW,QPOLD,SFQ,SFQP,ANP,ANQ,ACONTP,ACONTQ,
     B                  NSFNOD,SFV,SFVNUM,SFVNOD,SFVTIM,
     C                  ENKF,DSATM,QNEW_1,QTIMEP_1,
     D                  PLANT_FLAG,HTIMETEO,METEO,METEOTIM,METEOINP,
     E                  NMETEODATA)
C
      IMPLICIT NONE
      INCLUDE 'CATHY.H'
      INTEGER  I,J
      INTEGER  IPRT1,N,NSTR
      INTEGER  HSPATM,HTIATM,IETO,HTIDIR,HTINEU
      INTEGER  ANP,ANQ
      INTEGER  NDIR(3),NDIRC(3),NP(3),NQ(3),ZERO(3),ZEROC(3)
      INTEGER  ITER,NITERT,KBACKT,KBACK
      INTEGER  NNOD,NSF,NUMRES,NCELNL
      INTEGER  NSFNOD(NSFMAX,NNSFMX)
      INTEGER  ACONTP(*),ACONTQ(*)
      INTEGER  CONTP(3,*),CONTQ(3,*)
      INTEGER  IFATM(*),IFATMP(*)
      INTEGER  NSFNUM(*)
      INTEGER  SFEX(*),SFEXIT(*),SFEXP(*)
      INTEGER  SFV(2),SFVNUM(2,NSFMAX),SFVNOD(2,NSFMAX,NNSFMX)
      INTEGER  ENIFATM(NODMAX,MAXNENS),ENIFATMP(NODMAX,MAXNENS)
      LOGICAL  ENKF,SURF,PONDING,PONDP,DTGMIN
      REAL*8   TIME,TIMEP,DELTAT,DTMIN,DTREDM,DTREDS
      REAL*8   TETAF,DSATM
      REAL*8   OVFLNOD(NNOD),OVFLP(NNOD)
      REAL*8   PRESC(*),PTIM(*),PINP(3,*)
      REAL*8   Q(*),QTIM(*),QINP(3,*)
      REAL*8   ATMTIM(*),ATMINP(3,*)
      REAL*8   ARENOD(*)
      REAL*8   POLD(*),PTOLD(*),PNEW(*),PTNEW(*),PTIMEP(*)
      REAL*8   ATMACT(*),ATMPOT(*),ATMOLD(*)
      REAL*8   QPNEW(*),QPOLD(*),SFQ(NSFMAX,*),SFQP(NSFMAX,*)
      REAL*8   SFVTIM(2)
      REAL*8   QNEW_1,QTIMEP_1
C  PLANT variables
      LOGICAL  PLANT_FLAG
      INTEGER  HTIMETEO,NMETEODATA
      REAL*8   METEO(*),METEOINP(3,*),METEOTIM(*)

      INCLUDE 'SOILCHAR.H'
      INCLUDE 'SURFWATER.H'
      INCLUDE 'IOUNITS.H'
C           
      CALL VCOPYR(N,PNEW,PTIMEP)   
      CALL VCOPYR(ANP,QPNEW,QPOLD)
      CALL VCOPYI(NSF,SFEX,SFEXP)
      CALL VCOPYI(NSF,SFEXIT,SFEXP)
      DO I=1,NSF 
         DO J=1,NSFNUM(I)
            SFQ(I,J)=SFQP(I,J)
         END DO
      END DO
      CALL VCOPYI(NNOD,IFATM,IFATMP)
      CALL VCOPYR(NNOD,ATMACT,ATMOLD)
      TIME=TIMEP
      DELTAT=DELTAT*DTREDM - DTREDS
      IF (DELTAT .LE. DTMIN) THEN
         DELTAT=DTMIN
         DTGMIN=.FALSE. 
      ELSE
         DTGMIN=.TRUE.
      END IF
      TIME=TIME + DELTAT
      KBACKT=KBACKT + 1
      KBACK=KBACK + 1
      ITER=1  
      NITERT=0
C  
C  we need to check whether PTIM(1) < TIME <= PTIM(2), in which
C  case we need special handling for the back-stepping case
C        
C     IF (NP .NE. 0) THEN
         IF (TIME .GT. PTIM(2)) THEN
            CALL BCNXT('BACK-STEP NATM,NSF DIR',IIN8,IOUT2,IOUT19,NNOD,
     1                 NPMAX,IPRT1,NDIR,NDIRC,NP,NSTR,HTIDIR,TIME,
     2                 PTIM,PINP,CONTP,PRESC,ANP,ACONTP)
         ELSE
            CALL BCBAK('BACK-STEP NATM,NSF DIR',IOUT19,
     1                 IPRT1,NP,TIME,PTIM,PINP,CONTP,PRESC,ANP,ACONTP)
         END IF
C     END IF
C        
C  we need to check whether QTIM(1) < TIME <= QTIM(2), in which
C  case we need special handling for the back-stepping case
C        
C     IF (NQ .NE. 0) THEN
      CALL INIT0I(3,ZERO)
      CALL INIT0I(3,ZEROC)
         IF (TIME .GT. QTIM(2)) THEN
            CALL BCNXT('BACK-STEP NATM,NSF NEU',IIN9,IOUT2,IOUT20,NNOD,
     1                 NQMAX,IPRT1,ZERO,ZEROC,NQ,NSTR,HTINEU,TIME,
     2                 QTIM,QINP,CONTQ,Q,ANQ,ACONTQ)
         ELSE
            CALL BCBAK('BACK-STEP NATM,NSF NEU',IOUT20,
     1                 IPRT1,NQ,TIME,QTIM,QINP,CONTQ,Q,ANQ,ACONTQ)
         END IF

C  PLANT -------------------------

      IF (PLANT_FLAG) THEN
         IF (TIME.GT.METEOTIM(2)) THEN
           CALL METEONXT('BACK-STEP METEO DATA',IIN62,IOUT2,
     1                   IPRT1,HTIMETEO,TIME,
     2                   METEOTIM,METEOINP,METEO,NMETEODATA)
         ELSE
           CALL METEOBAK('BACK-STEP METEO DATA',
     1                   IPRT1,TIME,METEOTIM,METEOINP,METEO,NMETEODATA)
         ENDIF
      ENDIF

C  -------------------------------

C     END IF
C           
C  we need to check whether ATMTIM(1) < TIME <= ATMTIM(2), in which
C  case we need special handling for the back-stepping case.
C  Note: SWITCH_OLD is the version of SWITCH for the uncoupled,
C  subsurface flow only version of the model. It should be superceded
C  by the new ADRSTN routine, but for the time being we keep the old
C  version of SWITCH active as well.
C           
         if (time.lt.sfvtim(1)) then
             call sfvbak(nsfmax,nnsfmx,nsf,nsfnum,nsfnod,
     1                   sfv,sfvnum,sfvnod)
         else
             call sfvrec(nsfmax,nnsfmx,nsf,nsfnum,nsfnod,
     1                   sfv,sfvnum,sfvnod)
         end if
      IF (TIME .GT. ATMTIM(2)) THEN
         CALL ATMNXT(NNOD,HSPATM,HTIATM,IETO,TIME,IFATM,ARENOD,
     1               ATMPOT,ATMACT,ATMTIM,ATMINP,DELTAT,
     2               ANP,ANQ,ACONTP,ACONTQ,NSF,NSFNUM,NSFNOD,
     3               ENKF,DSATM,QNEW_1,QTIMEP_1)
      ELSE
         CALL ATMBAK(NNOD,TIME,IFATM,ARENOD,ATMPOT,ATMACT,
     1               ATMTIM,ATMINP,IETO,DELTAT,
     2               ENKF,DSATM,QNEW_1,QTIMEP_1)
      END IF
      IF (.NOT. SURF) THEN
         CALL SWITCH_OLD(NNOD,IFATM,ATMACT,ATMPOT,PNEW)
      ELSE                
         CALL ADRSTN(NNOD,IFATM,ATMPOT,ATMACT,PNEW)
      END IF
C
C
      CALL VCOPYR(N,POLD,PNEW)
      CALL WEIGHT(N,TETAF,PNEW,PTIMEP,PTNEW)
      CALL VCOPYR(N,PTOLD,PTNEW)
C    
C  variables and arrays for SURF_ROUTE
C
      IF (SURF) THEN
         PONDING = PONDP
         CALL VCOPYR(NNOD,OVFLNOD,OVFLP)
         AK_MAX = AK_MAX_p
         CALL VCOPYR(MAXCEL,Q_IN_KK_SN,Q_IN_KK_SN_P)
         CALL VCOPYR(MAXCEL,Q_OUT_KK_SN_1,Q_OUT_KK_SN_1_P)
         CALL VCOPYR(MAXCEL,Q_OUT_KK_SN_2,Q_OUT_KK_SN_2_P)
         CALL VCOPYR(MAXCEL,VOLUME_KK_SN,VOLUME_KK_SN_P)
         CALL VCOPYR(NUMRES,H_POOL_KK_VEC,H_POOL_KK_VEC_P)
      END IF
C
      RETURN
      END
