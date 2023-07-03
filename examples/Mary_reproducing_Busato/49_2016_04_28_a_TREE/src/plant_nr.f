C
C ---------------------------  PLANT_NR  -----------------------------------
C
C  Calculate Leaf water potential PSILEAF (Newton - Raphson method)
C 
C --------------------------------------------------------------------------
C
      SUBROUTINE PLANT_NR(ITYPE,EPSW,VXYLEM,A1,A2,VPD,GXYLEM,KXYLEM,
     1                   DGXYLEM,DKXYLEM,LASTOMA,DLASTOMA,GSTOMA,
     2                   PSIR,DPSIR,PSILEAF,S1,S2,NRFUN,DNRFUN,
     3                   PSILMEAN,METEO,TIME,NP,DHCANO,LAD,LAIP,
     4                   GROWTH_FLAG,GROWFACTOR,SALE,SALT_FLAG,
     5                   SALT_FACTOR,COMP)


      IMPLICIT NONE

      INCLUDE 'CATHY.H'

C  Input variables
      INTEGER ITYPE,NP
      REAL*8  VPD,EPSW,VXYLEM,A1(DATALEAFMAX),A2(DATALEAFMAX)
      REAL*8  S1,S2,PSILMEAN,DHCANO,LAD,LAIP,GROWFACTOR
      REAL*8  SALT_FACTOR,SALE
      LOGICAL GROWTH_FLAG,SALT_FLAG
C  Input/Output variables
      REAL*8  GXYLEM,KXYLEM,DGXYLEM,DKXYLEM,LASTOMA,DLASTOMA,PSILEAF
      REAL*8  GXYLEM1,KXYLEM1,LASTOMA1
      REAL*8  GXYLEM2,KXYLEM2,LASTOMA2
      REAL*8  GSTOMA1(PLMAX,DATALEAFMAX),PSIR1
      REAL*8  GSTOMA2(PLMAX,DATALEAFMAX),PSIR2
      REAL*8  GSTOMA(PLMAX,DATALEAFMAX),DGSTOMA(DATALEAFMAX),PSIR,DPSIR
      REAL*8  E1,E2,DE1,DE2,NRFUN,DNRFUN,E11,E12,E22,E21,NRFUN1,NRFUN2
      REAL*8  COMP
C  Local variables 
      INTEGER ITER,I,J
      REAL*8  PSILEAFNEW,SCARTO,SCARTONEW,METEO(*),TIME
      REAL*8  E2_LEAF(DATALEAFMAX),DE2_LEAF(DATALEAFMAX)
      REAL*8  E2_LEAF1(DATALEAFMAX),E2_LEAF2(DATALEAFMAX)
      REAL*8  GS_NIGHT,II,X1,X2,NRFUN0

      INCLUDE 'IOUNITS.H' 
      INCLUDE 'PLANT.H'

C --------------------------------------------------------------------------

C Zeri iniziali Regula Falsi
C      WRITE(7894,*) TIME,'TIME'
      DO J=1,1001
         II=-(1002-J)


         IF(GROWTH_FLAG) THEN
         GXYLEM  = GXYLEM_MAX(ITYPE)*EXP(-(-II/D_GXYLEM(ITYPE))**
     1             C_GXYLEM(ITYPE))*LAIP
         KXYLEM  = GXYLEM/(HCANO(ITYPE)*GROWFACTOR)
         ELSE
         GXYLEM  = GXYLEM_MAX(ITYPE)*EXP(-(-II/D_GXYLEM(ITYPE))**
     1             C_GXYLEM(ITYPE))*LAI(ITYPE)
         KXYLEM  = GXYLEM/HCANO(ITYPE)
         ENDIF

         LASTOMA  = LA_MAX(ITYPE)*CCO2ATM/CCO2STAR*EXP(-LA_BETA(ITYPE)
     1              *(PSILMEAN-LA_PSILMAX(ITYPE))**2)
         E2 = 0
         DO I=1,DATA_LEAF(ITYPE)
c             GSTOMA(NP,I)= A1(I)/(A2(I)+SSTOMA*CCO2ATM)*
c     1                     (-1+SQRT(CCO2ATM/(ASTOMA*LASTOMA*VPD)))
c     2                     +GS_NIGHT
             GSTOMA(NP,I)=(-A1(I)*(A2(I)-CCO2ATM+2.D0*COMP)/(A2(I)+
     1                   CCO2ATM)**2.D0)+
     2                   SQRT(ASTOMA*VPD*LASTOMA*(A1(I)**2.D0)*
     3                   (CCO2ATM-COMP)*(A2(I)+COMP)*(A2(I)+
     4                   CCO2ATM-2.D0*ASTOMA*VPD*LASTOMA)**2.D0*
     5                   (A2(I)+CCO2ATM-ASTOMA*VPD*LASTOMA))/
     6                   (ASTOMA*VPD*LASTOMA*(A2(I)+CCO2ATM)**2
     7                   *(A2(I)+CCO2ATM-ASTOMA*VPD*LASTOMA))

             IF (GROWTH_FLAG) THEN 
             E2_LEAF(I)  = ASTOMA*GSTOMA(NP,I)*VPD*EPSW*LAD
     1                     *ACANO(ITYPE)*GROWFACTOR 
             ELSE
             E2_LEAF(I)  = ASTOMA*GSTOMA(NP,I)*VPD*EPSW*LAD
     1                     *ACANO(ITYPE)
             ENDIF

             E2 = E2 + E2_LEAF(I)
         ENDDO

         IF(GROWTH_FLAG) THEN
         PSIR  = (KXYLEM*VXYLEM*(II+HCANO(ITYPE)*GROWFACTOR)+S1)/
     1           (S2+KXYLEM*VXYLEM)
         E1     = -KXYLEM*VXYLEM*(II-PSIR+HCANO(ITYPE)*GROWFACTOR)
         ELSE
         PSIR  = (KXYLEM*VXYLEM*(II+HCANO(ITYPE))+S1)/
     1           (S2+KXYLEM*VXYLEM)
         E1     = -KXYLEM*VXYLEM*(II-PSIR+HCANO(ITYPE))
         ENDIF

         E2     = E2*DHCANO
         NRFUN  = E2-E1

         IF(NRFUN.LT.0.D0) X1=II
         IF(NRFUN.GT.0.D0) X2=II
c         WRITE(7894,*) II,NRFUN,E1,E2,KXYLEM,GROWFACTOR
      ENDDO

c      WRITE(7896,*) X1,X2

      SCARTO = 2.0D0*TOLLNR
      ITER      = 0
      PSILEAF   = PSILEAF0
      GS_NIGHT  = 0.018

      IF (NPLANT.LE.10) THEN
      WRITE(IOUT71,*) TIME,'TIME'
      WRITE(IOUT71,*) '     ITER        PSILEAF        x1          x2  
     1        SCARTO                FUN              FUN1          FUN2'
      ENDIF
C  ------------------- Start Newton - Raphson iteration --------------------  

      DO WHILE((SCARTO.GE.TOLLNR).AND.(ITER.LE.ITMAXNR))
         ITER    = ITER+1

         IF (GROWTH_FLAG) THEN
         GXYLEM1  = GXYLEM_MAX(ITYPE)*EXP(-(-X1/D_GXYLEM(ITYPE))**
     1              C_GXYLEM(ITYPE))*LAIP
         GXYLEM2  = GXYLEM_MAX(ITYPE)*EXP(-(-X2/D_GXYLEM(ITYPE))**
     1              C_GXYLEM(ITYPE))*LAIP
         KXYLEM1  = GXYLEM1/(HCANO(ITYPE)*GROWFACTOR)
         KXYLEM2  = GXYLEM2/(HCANO(ITYPE)*GROWFACTOR)
         ELSE
         GXYLEM1  = GXYLEM_MAX(ITYPE)*EXP(-(-X1/D_GXYLEM(ITYPE))**
     1             C_GXYLEM(ITYPE))*LAI(ITYPE)
         GXYLEM2  = GXYLEM_MAX(ITYPE)*EXP(-(-X2/D_GXYLEM(ITYPE))**
     1             C_GXYLEM(ITYPE))*LAI(ITYPE)
         KXYLEM1  = GXYLEM1/HCANO(ITYPE)
         KXYLEM2  = GXYLEM2/HCANO(ITYPE)
         ENDIF

         LASTOMA1  = LA_MAX(ITYPE)*CCO2ATM/CCO2STAR*EXP(-LA_BETA(ITYPE)
     1              *(PSILMEAN-LA_PSILMAX(ITYPE))**2)
         LASTOMA2  = LA_MAX(ITYPE)*CCO2ATM/CCO2STAR*EXP(-LA_BETA(ITYPE)
     1              *(PSILMEAN-LA_PSILMAX(ITYPE))**2)
         E21 = 0
         E22 = 0
         DO I=1,DATA_LEAF(ITYPE)
c             GSTOMA1(NP,I)= A1(I)/(A2(I)+SSTOMA*CCO2ATM)*
c     1                     (-1+SQRT(CCO2ATM/(ASTOMA*LASTOMA1*VPD)))
c     2                     +GS_NIGHT
c             GSTOMA2(NP,I)= A1(I)/(A2(I)+SSTOMA*CCO2ATM)*
c     1                     (-1+SQRT(CCO2ATM/(ASTOMA*LASTOMA2*VPD)))
c     2                     +GS_NIGHT
             GSTOMA1(NP,I)=(-A1(I)*(A2(I)-CCO2ATM+2.D0*COMP)/(A2(I)+
     1                   CCO2ATM)**2.D0)+
     2                   SQRT(ASTOMA*VPD*LASTOMA1*(A1(I)**2.D0)*
     3                   (CCO2ATM-COMP)*(A2(I)+COMP)*(A2(I)+
     4                   CCO2ATM-2.D0*ASTOMA*VPD*LASTOMA1)**2.D0*
     5                   (A2(I)+CCO2ATM-ASTOMA*VPD*LASTOMA1))/
     6                   (ASTOMA*VPD*LASTOMA1*(A2(I)+CCO2ATM)**2
     7                   *(A2(I)+CCO2ATM-ASTOMA*VPD*LASTOMA1))
             GSTOMA2(NP,I)=(-A1(I)*(A2(I)-CCO2ATM+2.D0*COMP)/(A2(I)+
     1                   CCO2ATM)**2.D0)+
     2                   SQRT(ASTOMA*VPD*LASTOMA2*(A1(I)**2.D0)*
     3                   (CCO2ATM-COMP)*(A2(I)+COMP)*(A2(I)+
     4                   CCO2ATM-2.D0*ASTOMA*VPD*LASTOMA2)**2.D0*
     5                   (A2(I)+CCO2ATM-ASTOMA*VPD*LASTOMA2))/
     6                   (ASTOMA*VPD*LASTOMA2*(A2(I)+CCO2ATM)**2
     7                   *(A2(I)+CCO2ATM-ASTOMA*VPD*LASTOMA2))

             IF (GROWTH_FLAG) THEN
             E2_LEAF1(I)  = ASTOMA*GSTOMA1(NP,I)*VPD*EPSW*LAD
     1                       *ACANO(ITYPE)*GROWFACTOR
             E2_LEAF2(I)  = ASTOMA*GSTOMA2(NP,I)*VPD*EPSW*LAD
     1                       *ACANO(ITYPE)*GROWFACTOR
             ELSE
             E2_LEAF1(I)  = ASTOMA*GSTOMA1(NP,I)*VPD*EPSW*LAD
     1                       *ACANO(ITYPE) 
             E2_LEAF2(I)  = ASTOMA*GSTOMA2(NP,I)*VPD*EPSW*LAD
     1                       *ACANO(ITYPE) 
             ENDIF

             E21 = E21 + E2_LEAF1(I)
             E22 = E22 + E2_LEAF2(I)
         ENDDO

         
         IF (GROWTH_FLAG) THEN
         PSIR1  = (KXYLEM1*VXYLEM*(X1+HCANO(ITYPE)*GROWFACTOR)+S1)/
     1           (S2+KXYLEM1*VXYLEM)
         PSIR2  = (KXYLEM2*VXYLEM*(X2+HCANO(ITYPE)*GROWFACTOR)+S1)/
     1           (S2+KXYLEM2*VXYLEM)
         E11     = -KXYLEM1*VXYLEM*(X1-PSIR1+HCANO(ITYPE)*GROWFACTOR)
         E12     = -KXYLEM2*VXYLEM*(X2-PSIR2+HCANO(ITYPE)*GROWFACTOR)
         ELSE 
         PSIR1  = (KXYLEM1*VXYLEM*(X1+HCANO(ITYPE))+S1)/
     1           (S2+KXYLEM1*VXYLEM)
         PSIR2  = (KXYLEM2*VXYLEM*(X2+HCANO(ITYPE))+S1)/
     1           (S2+KXYLEM2*VXYLEM)
         E11     = -KXYLEM1*VXYLEM*(X1-PSIR1+HCANO(ITYPE))
         E12     = -KXYLEM2*VXYLEM*(X2-PSIR2+HCANO(ITYPE))
         ENDIF

         E21     = E21*DHCANO
         E22     = E22*DHCANO
         NRFUN1  = E21-E11
         NRFUN2  = E22-E12

         PSILEAF = X1-NRFUN1*(X2-X1)/(NRFUN2-NRFUN1)
  
C  Update

         IF(GROWTH_FLAG) THEN
         GXYLEM  = GXYLEM_MAX(ITYPE)*EXP(-(-PSILEAF/D_GXYLEM(ITYPE))**
     1             C_GXYLEM(ITYPE))*LAIP
         KXYLEM  = GXYLEM/(HCANO(ITYPE)*GROWFACTOR)
         ELSE
         GXYLEM  = GXYLEM_MAX(ITYPE)*EXP(-(-PSILEAF/D_GXYLEM(ITYPE))**
     1             C_GXYLEM(ITYPE))*LAI(ITYPE)
         KXYLEM  = GXYLEM/HCANO(ITYPE)
         ENDIF
         LASTOMA  = LA_MAX(ITYPE)*CCO2ATM/CCO2STAR*EXP(-LA_BETA(ITYPE)
     1              *(PSILMEAN-LA_PSILMAX(ITYPE))**2.d0)
         E2 = 0
         DO I=1,DATA_LEAF(ITYPE)
c             GSTOMA(NP,I)= A1(I)/(A2(I)+SSTOMA*CCO2ATM)*
c     1                     (-1+SQRT(CCO2ATM/(ASTOMA*LASTOMA*VPD)))
c     2                     +GS_NIGHT
             GSTOMA(NP,I)=(-A1(I)*(A2(I)-CCO2ATM+2.D0*COMP)/(A2(I)+
     1                   CCO2ATM)**2.D0)+
     2                   SQRT(ASTOMA*VPD*LASTOMA*(A1(I)**2.D0)*
     3                   (CCO2ATM-COMP)*(A2(I)+COMP)*(A2(I)+
     4                   CCO2ATM-2.D0*ASTOMA*VPD*LASTOMA)**2.D0*
     5                   (A2(I)+CCO2ATM-ASTOMA*VPD*LASTOMA))/
     6                   (ASTOMA*VPD*LASTOMA*(A2(I)+CCO2ATM)**2
     7                   *(A2(I)+CCO2ATM-ASTOMA*VPD*LASTOMA))

             IF(GROWTH_FLAG) THEN
             E2_LEAF(I)  = ASTOMA*GSTOMA(NP,I)*VPD*EPSW*LAD
     1                     *ACANO(ITYPE)*GROWFACTOR 
             ELSE
             E2_LEAF(I)  = ASTOMA*GSTOMA(NP,I)*VPD*EPSW*LAD
     1                     *ACANO(ITYPE) 
             ENDIF

             E2 = E2 + E2_LEAF(I)
         ENDDO

         IF(GROWTH_FLAG) THEN
         PSIR  = (KXYLEM*VXYLEM*(PSILEAF+HCANO(ITYPE)*GROWFACTOR)+S1)/
     1           (S2+KXYLEM*VXYLEM)
         E1     = -KXYLEM*VXYLEM*(PSILEAF-PSIR+HCANO(ITYPE)*GROWFACTOR)
         ELSE
         PSIR  = (KXYLEM*VXYLEM*(PSILEAF+HCANO(ITYPE))+S1)/
     1           (S2+KXYLEM*VXYLEM)
         E1     = -KXYLEM*VXYLEM*(PSILEAF-PSIR+HCANO(ITYPE))
         ENDIF

         E2     = E2*DHCANO
         NRFUN  = E2-E1

         SCARTO = ABS(NRFUN*1e7)

         IF(NPLANT.LE.10) THEN
         WRITE(IOUT71,100) ITER,PSILEAF,X1,X2,SCARTO,tollnr,
     1                 NRFUN,NRFUN1,NRFUN2,GSTOMA(1,1),VPD,LASTOMA,
     2                 PSILMEAN,LA_PSILMAX(1)
         ENDIF

      IF(NRFUN*NRFUN1.LT.0) THEN
        X1=X1
        X2=PSILEAF
      ELSEIF(NRFUN*NRFUN2.LT.0) THEN
        X1=PSILEAF
        X2=X2
      ENDIF

        
      ENDDO

C Fine Newton Raphson (fine do while)

      IF (ITER.GT.ITMAXNR) THEN
c      WRITE(7892,*) TIME,'TIME'
      NRFUN=1.D0
      DO J=1,1001
         NRFUN0 = NRFUN
         II=-J
         
         IF (GROWTH_FLAG) THEN
         GXYLEM  = GXYLEM_MAX(ITYPE)*EXP(-(-II/D_GXYLEM(ITYPE))**
     1             C_GXYLEM(ITYPE))*LAIP
         KXYLEM  = GXYLEM/(HCANO(ITYPE)*GROWFACTOR)
         ELSE
         GXYLEM  = GXYLEM_MAX(ITYPE)*EXP(-(-II/D_GXYLEM(ITYPE))**
     1             C_GXYLEM(ITYPE))*LAI(ITYPE)
         KXYLEM  = GXYLEM/HCANO(ITYPE)
         ENDIF

         LASTOMA  = LA_MAX(ITYPE)*CCO2ATM/CCO2STAR*EXP(-LA_BETA(ITYPE)
     1              *(PSILMEAN-LA_PSILMAX(ITYPE))**2.d0)
         E2 = 0
         DO I=1,DATA_LEAF(ITYPE)
c             GSTOMA(NP,I)= A1(I)/(A2(I)+SSTOMA*CCO2ATM)*
c     1                     (-1+SQRT(CCO2ATM/(ASTOMA*LASTOMA*VPD)))
c     2                     +GS_NIGHT
             GSTOMA(NP,I)=(-A1(I)*(A2(I)-CCO2ATM+2.D0*COMP)/(A2(I)+
     1                   CCO2ATM)**2.D0)+
     2                   SQRT(ASTOMA*VPD*LASTOMA*(A1(I)**2.D0)*
     3                   (CCO2ATM-COMP)*(A2(I)+COMP)*(A2(I)+
     4                   CCO2ATM-2.D0*ASTOMA*VPD*LASTOMA)**2.D0*
     5                   (A2(I)+CCO2ATM-ASTOMA*VPD*LASTOMA))/
     6                   (ASTOMA*VPD*LASTOMA*(A2(I)+CCO2ATM)**2
     7                   *(A2(I)+CCO2ATM-ASTOMA*VPD*LASTOMA))

             IF (GROWTH_FLAG) THEN
             E2_LEAF(I)  = ASTOMA*GSTOMA(NP,I)*VPD*EPSW*LAD
     1                       *ACANO(ITYPE)*GROWFACTOR
             ELSE
             E2_LEAF(I)  = ASTOMA*GSTOMA(NP,I)*VPD*EPSW*LAD
     1                       *ACANO(ITYPE)
             ENDIF

             E2 = E2 + E2_LEAF(I)
         ENDDO

         IF(GROWTH_FLAG) THEN
         PSIR  = (KXYLEM*VXYLEM*(II+HCANO(ITYPE)*GROWFACTOR)+S1)/
     1           (S2+KXYLEM*VXYLEM)
         E1     = -KXYLEM*VXYLEM*(II-PSIR+HCANO(ITYPE)*GROWFACTOR)
         ELSE
         PSIR  = (KXYLEM*VXYLEM*(II+HCANO(ITYPE))+S1)/
     1           (S2+KXYLEM*VXYLEM)
         E1     = -KXYLEM*VXYLEM*(II-PSIR+HCANO(ITYPE))
         ENDIF

         E2     = E2*DHCANO
         NRFUN  = E2-E1
c         IF ((NRFUN.LT.NRFUN0).AND.(ABS(NRFUN).LT.TOLLNR)) THEN
         IF (NRFUN.LT.NRFUN0) THEN
            PSILEAF= II
c            write(7892,*)'ENTRA'
         ENDIF      

c         WRITE(7892,*) II,NRFUN,PSILEAF
      ENDDO

         IF (GROWTH_FLAG) THEN         
         GXYLEM  = GXYLEM_MAX(ITYPE)*EXP(-(-PSILEAF/D_GXYLEM(ITYPE))**
     1             C_GXYLEM(ITYPE))*LAIP
         KXYLEM  = GXYLEM/(HCANO(ITYPE)*GROWFACTOR)
         PSIR  = (KXYLEM*VXYLEM*(PSILEAF+HCANO(ITYPE)*GROWFACTOR)+S1)/
     1           (S2+KXYLEM*VXYLEM)
         ELSE
         GXYLEM  = GXYLEM_MAX(ITYPE)*EXP(-(-PSILEAF/D_GXYLEM(ITYPE))**
     1             C_GXYLEM(ITYPE))*LAI(ITYPE)
         KXYLEM  = GXYLEM/HCANO(ITYPE)
         PSIR  = (KXYLEM*VXYLEM*(PSILEAF+HCANO(ITYPE))+S1)/
     1           (S2+KXYLEM*VXYLEM)
         ENDIF

c      WRITE(7893,*) TIME,PSILEAF
      ENDIF







C ---------------------------------------------------------------------
C Se PSILEAF-PSIR < HCANO si impone PSILEAF-PSIR = HCANO, ovvero flusso nullo
 
c         IF (ABS(PSILEAF-PSIR).LE.HCANO(ITYPE)) THEN
c            PSILEAF = (S1-S2*HCANO(ITYPE))/S2
c            GXYLEM  = GXYLEM_MAX(ITYPE)*EXP(-(-PSILEAF/
c     1                D_GXYLEM(ITYPE))**C_GXYLEM(ITYPE))*LAI(ITYPE)
c            KXYLEM  = GXYLEM/HCANO(ITYPE)
c            LASTOMA = LA_MAX(ITYPE)*CCO2ATM/CCO2STAR*EXP(-LA_BETA(ITYPE)
c     1                *(PSILMEAN-LA_PSILMAX(ITYPE))**2)
c
c            DO I=1,DATA_LEAF(ITYPE) 
c               GSTOMA  = A1(I)/(A2(I)+SSTOMA*CCO2ATM)*
c     1                   (-1+SQRT(CCO2ATM/(ASTOMA*LASTOMA*VPD)))
c            ENDDO
c            PSIR    = (KXYLEM*VXYLEM*(PSILEAF+HCANO(ITYPE))
c     1                +S1)/(S2+KXYLEM*VXYLEM)
c         ENDIF 


C --------------------------------------------------------------------

C  Print

      IF(NPLANT.LE.10) THEN
      IF(SCARTO.LE.TOLLNR) THEN
         WRITE(IOUT71,200) ITER
         WRITE(IOUT2,200)  ITER
      ELSE
         WRITE(IOUT71,300) ITER
         WRITE(IOUT2,300)  ITER
      ENDIF
      ENDIF
c
C  Format

 100  FORMAT(i10,2x,16(e16.9,2x))
 200  FORMAT('Convergence achieved in',I4,' iterations')
 300  FORMAT('Convergence NOT achieved in',I4,' iterations')

C --------------------------------------------------------------------------

      RETURN
      END
