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
     5                   SALT_FACTOR)


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
      REAL*8  GSTOMA(PLMAX,DATALEAFMAX),DGSTOMA(DATALEAFMAX),PSIR,DPSIR
      REAL*8  E1,E2,DE1,DE2,NRFUN,DNRFUN
C  Local variables 
      INTEGER ITER,I,J
      REAL*8  PSILEAFNEW,SCARTO,SCARTONEW,METEO(*),TIME
      REAL*8  E2_LEAF(DATALEAFMAX),DE2_LEAF(DATALEAFMAX)
      REAL*8  GS_NIGHT,II

      INCLUDE 'IOUNITS.H' 
      INCLUDE 'PLANT.H'

C --------------------------------------------------------------------------

      SCARTONEW = 2.0D0*TOLLNR
      SCARTO    = SCARTONEW
      ITER      = 0
         PSILEAF   = PSILEAF0
      GS_NIGHT  = 0.018

      WRITE(IOUT71,*) TIME,'TIME'
      WRITE(IOUT71,*) '     ITER        PSILEAFNEW        PSILEAFOLD    
     1        SCARTO              PSIR           LASTOMA    GSTOMA    
     3       DGSTOMA     E1       E2     FUN       DFUN'

C  ------------------- Start Newton - Raphson iteration --------------------  

      DO WHILE((SCARTONEW.GE.TOLLNR).AND.(ITER.LE.ITMAXNR))
         ITER    = ITER+1

C  Calculate the xylem conductance and its derivative 
C  according to Daly et al., 2004.
C  Controlla se gxylem va moltiplicata o meno per il LAI

         IF(GROWTH_FLAG) THEN
         GXYLEM  = GXYLEM_MAX(ITYPE)*EXP(-(-PSILEAF/D_GXYLEM(ITYPE))**
     1             C_GXYLEM(ITYPE))*LAIP
         DGXYLEM = GXYLEM*(C_GXYLEM(ITYPE)/D_GXYLEM(ITYPE)*(-
     1             PSILEAF/D_GXYLEM(ITYPE))**(C_GXYLEM(ITYPE)-1))*
     2             LAIP
         KXYLEM  = GXYLEM/(HCANO(ITYPE)*GROWFACTOR)
         DKXYLEM = DGXYLEM/(HCANO(ITYPE)*GROWFACTOR)
         ELSE
         GXYLEM  = GXYLEM_MAX(ITYPE)*EXP(-(-PSILEAF/D_GXYLEM(ITYPE))**
     1             C_GXYLEM(ITYPE))*LAI(ITYPE)
         DGXYLEM = GXYLEM*(C_GXYLEM(ITYPE)/D_GXYLEM(ITYPE)*(-
     1             PSILEAF/D_GXYLEM(ITYPE))**(C_GXYLEM(ITYPE)-1))*
     2             LAI(ITYPE)
         KXYLEM  = GXYLEM/HCANO(ITYPE)
         DKXYLEM = DGXYLEM/HCANO(ITYPE)
         ENDIF

      
C  Calculate the stomatal conductance and its derivative 
C  according to optimization theories

         LASTOMA  = LA_MAX(ITYPE)*CCO2ATM/CCO2STAR*EXP(-LA_BETA(ITYPE)
     1              *(PSILMEAN-LA_PSILMAX(ITYPE))**2)
         DLASTOMA = LA_MAX(ITYPE)*CCO2ATM/CCO2STAR*EXP(-LA_BETA(ITYPE)
     1               *(PSILMEAN-LA_PSILMAX(ITYPE))**2)*(-2*
     2               LA_BETA(ITYPE)*(PSILMEAN-LA_PSILMAX(ITYPE)))

C  SALT -------------------------------------------
C      
C         IF(SALT_FLAG) THEN
C            LASTOMA = LASTOMA*SALT_FACTOR
C            DLASTOMA = DLASTOMA*SALT_FACTOR
C         ENDIF
C         
C  ------------------------------------------------

         E2 = 0
         DE2 = 0

         DO I=1,DATA_LEAF(ITYPE)

c         IF(METEO(3).LT.0.02)THEN
c             GSTOMA(NP,I)=0.018
c             DGSTOMA(I)=0
c         ELSE
             GSTOMA(NP,I)= A1(I)/(A2(I)+SSTOMA*CCO2ATM)*
     1                     (-1+SQRT(CCO2ATM/(ASTOMA*LASTOMA*VPD)))
     2                     +GS_NIGHT
             DGSTOMA(I) = -0.5*(A1(I)/(A2(I)+SSTOMA*CCO2ATM))*
     1                   (CCO2ATM/(ASTOMA*LASTOMA*VPD))**(-0.5)*CCO2ATM/
     2                   (ASTOMA*VPD*LASTOMA**2)*DLASTOMA
c         ENDIF

             IF (GROWTH_FLAG) THEN
               E2_LEAF(I)  = ASTOMA*GSTOMA(NP,I)*VPD*EPSW*LAD
     1                       *ACANO(ITYPE)*GROWFACTOR 
               DE2_LEAF(I) = ASTOMA*DGSTOMA(I)*VPD*EPSW*LAD
     1                       *ACANO(ITYPE)*GROWFACTOR
             ELSE
               E2_LEAF(I)  = ASTOMA*GSTOMA(NP,I)*VPD*EPSW*LAD
     1                       *ACANO(ITYPE) 
               DE2_LEAF(I) = ASTOMA*DGSTOMA(I)*VPD*EPSW*LAD
     1                       *ACANO(ITYPE)
             ENDIF

             E2 = E2 + E2_LEAF(I)
             DE2 = DE2 + DE2_LEAF(I)
         ENDDO

C  Calculate the water potential in the xylem and its derivative

         IF(GROWTH_FLAG) THEN
           PSIR  = (KXYLEM*VXYLEM*(PSILEAF+HCANO(ITYPE)*GROWFACTOR)+S1)/
     1             (S2+KXYLEM*VXYLEM)
           DPSIR = ((DKXYLEM*VXYLEM*(PSILEAF+HCANO(ITYPE)*GROWFACTOR)+
     1             KXYLEM*VXYLEM)*(S2+KXYLEM*VXYLEM)-(KXYLEM*VXYLEM*
     2             (PSILEAF+HCANO(ITYPE)*GROWFACTOR)+S1)*
     3             (DKXYLEM*VXYLEM))/(S2+KXYLEM*VXYLEM)**2
C          Calculate NRFUN and DNRFUN for Newton - Raphson iterations 
           E1    = -KXYLEM*VXYLEM*(PSILEAF-PSIR+HCANO(ITYPE)*GROWFACTOR)
           DE1   = -DKXYLEM*VXYLEM*(PSILEAF-PSIR+HCANO(ITYPE)*
     1              GROWFACTOR)-KXYLEM*VXYLEM*(1-DPSIR)
         ELSE
           PSIR  = (KXYLEM*VXYLEM*(PSILEAF+HCANO(ITYPE))+S1)/
     1             (S2+KXYLEM*VXYLEM)
           DPSIR =((DKXYLEM*VXYLEM*(PSILEAF+HCANO(ITYPE))+KXYLEM*VXYLEM)
     1            *(S2+KXYLEM*VXYLEM)-(KXYLEM*VXYLEM*(PSILEAF+HCANO
     2            (ITYPE))+S1)*(DKXYLEM*VXYLEM))/(S2+KXYLEM*VXYLEM)**2
C          Calculate NRFUN and DNRFUN for Newton - Raphson iterations 
           E1     = -KXYLEM*VXYLEM*(PSILEAF-PSIR+HCANO(ITYPE))
           DE1    = -DKXYLEM*VXYLEM*(PSILEAF-PSIR+HCANO(ITYPE))-KXYLEM*
     1               VXYLEM*(1-DPSIR)
         ENDIF

         E2     = E2*DHCANO
         DE2    = DE2*DHCANO
         NRFUN  = E2-E1
         DNRFUN = DE2-DE1

C  Calculate Psi_leaf and scarto

         PSILEAFNEW = PSILEAF -0.1*NRFUN/DNRFUN
         SCARTONEW  = ABS(PSILEAFNEW-PSILEAF)
  
C  Update

         SCARTO  = SCARTONEW    
         IF ((PSILEAF.LE.PSILMAX(ITYPE)).OR.(PSILEAF.GT.0)) THEN
c         IF (PSILEAF.LE.PSILMAX(ITYPE)) THEN
             ITER=999
         ENDIF
         WRITE(IOUT71,100) ITER,PSILEAFNEW,PSILEAF,SCARTONEW,PSIR,
     1                     LASTOMA,GSTOMA(1,10),DGSTOMA(10),E1,E2,NRFUN,
     2                     DNRFUN
         PSILEAF = PSILEAFNEW

        
      ENDDO

C Fine Newton Raphson (fine do while)


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IF(PSILEAF.LT.-250) THEN
      WRITE(7895,*) TIME,'TIME'
         DO J=1,500
           II=-j
           GXYLEM  = GXYLEM_MAX(ITYPE)*EXP(-(-II/D_GXYLEM(ITYPE))**
     1               C_GXYLEM(ITYPE))*LAI(ITYPE)
           DGXYLEM = GXYLEM*(C_GXYLEM(ITYPE)/D_GXYLEM(ITYPE)*(-
     1               II/D_GXYLEM(ITYPE))**(C_GXYLEM(ITYPE)-1))*
     2               LAI(ITYPE)
           KXYLEM  = GXYLEM/HCANO(ITYPE)
           DKXYLEM = DGXYLEM/HCANO(ITYPE)
           LASTOMA = LA_MAX(ITYPE)*CCO2ATM/CCO2STAR*EXP(-LA_BETA(ITYPE)
     1              *(PSILMEAN-LA_PSILMAX(ITYPE))**2)
           DLASTOMA = LA_MAX(ITYPE)*CCO2ATM/CCO2STAR*EXP(-LA_BETA(ITYPE)
     1               *(PSILMEAN-LA_PSILMAX(ITYPE))**2)*(-2*
     2               LA_BETA(ITYPE)*(PSILMEAN-LA_PSILMAX(ITYPE)))
         E2 = 0
         DE2 = 0
         DO I=1,DATA_LEAF(ITYPE)
             GSTOMA(NP,I)= A1(I)/(A2(I)+SSTOMA*CCO2ATM)*
     1                     (-1+SQRT(CCO2ATM/(ASTOMA*LASTOMA*VPD)))
     2                     +GS_NIGHT
             DGSTOMA(I) = -0.5*(A1(I)/(A2(I)+SSTOMA*CCO2ATM))*
     1                   (CCO2ATM/(ASTOMA*LASTOMA*VPD))**(-0.5)*CCO2ATM/
     2                   (ASTOMA*VPD*LASTOMA**2)*DLASTOMA
             E2_LEAF(I)  = ASTOMA*GSTOMA(NP,I)*VPD*EPSW*LAD
     1                     *ACANO(ITYPE) 
             DE2_LEAF(I) = ASTOMA*DGSTOMA(I)*VPD*EPSW*LAD
     1                     *ACANO(ITYPE)
             E2 = E2 + E2_LEAF(I)
             DE2 = DE2 + DE2_LEAF(I)
         ENDDO
         PSIR  = (KXYLEM*VXYLEM*(II+HCANO(ITYPE))+S1)/
     1           (S2+KXYLEM*VXYLEM)
         DPSIR =((DKXYLEM*VXYLEM*(II+HCANO(ITYPE))+KXYLEM*VXYLEM)
     1          *(S2+KXYLEM*VXYLEM)-(KXYLEM*VXYLEM*(II+HCANO
     2          (ITYPE))+S1)*(DKXYLEM*VXYLEM))/(S2+KXYLEM*VXYLEM)**2
         E1     = -KXYLEM*VXYLEM*(II-PSIR+HCANO(ITYPE))
         DE1    = -DKXYLEM*VXYLEM*(II-PSIR+HCANO(ITYPE))-KXYLEM*
     1             VXYLEM*(1-DPSIR)
         E2     = E2*DHCANO
         DE2    = DE2*DHCANO
         NRFUN  = E2-E1
         DNRFUN = DE2-DE1
         WRITE(7895,*) II,NRFUN
      ENDDO
      ENDIF

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC



C ---------------------------------------------------------------------
C Se PSILEAF-PSIR < HCANO si impone PSILEAF-PSIR = HCANO, ovvero flusso nullo
 
      IF (GROWTH_FLAG) THEN
         IF (ABS(PSILEAF-PSIR).LE.HCANO(ITYPE)*GROWFACTOR) THEN
            PSILEAF = (S1-S2*HCANO(ITYPE)*GROWFACTOR)/S2
            GXYLEM  = GXYLEM_MAX(ITYPE)*EXP(-(-PSILEAF/
     1                D_GXYLEM(ITYPE))**C_GXYLEM(ITYPE))*LAIP
            KXYLEM  = GXYLEM/(HCANO(ITYPE)*GROWFACTOR)
            LASTOMA = LA_MAX(ITYPE)*CCO2ATM/CCO2STAR*EXP(-LA_BETA(ITYPE)
     1                *(PSILMEAN-LA_PSILMAX(ITYPE))**2)

C  SALT -------------------------------------------
C      
C         IF(SALT_FLAG) THEN
C            LASTOMA = LASTOMA*SALT_FACTOR
C         ENDIF
C         
C  ------------------------------------------------

            DO I=1,DATA_LEAF(ITYPE) 
               GSTOMA  = A1(I)/(A2(I)+SSTOMA*CCO2ATM)*
     1                   (-1+SQRT(CCO2ATM/(ASTOMA*LASTOMA*VPD)))
            ENDDO
            PSIR    = (KXYLEM*VXYLEM*(PSILEAF+HCANO(ITYPE)*GROWFACTOR)
     1                +S1)/(S2+KXYLEM*VXYLEM)
         ENDIF
      ELSE
         IF (ABS(PSILEAF-PSIR).LE.HCANO(ITYPE)) THEN
            PSILEAF = (S1-S2*HCANO(ITYPE))/S2
            GXYLEM  = GXYLEM_MAX(ITYPE)*EXP(-(-PSILEAF/
     1                D_GXYLEM(ITYPE))**C_GXYLEM(ITYPE))*LAI(ITYPE)
            KXYLEM  = GXYLEM/HCANO(ITYPE)
            LASTOMA = LA_MAX(ITYPE)*CCO2ATM/CCO2STAR*EXP(-LA_BETA(ITYPE)
     1                *(PSILMEAN-LA_PSILMAX(ITYPE))**2)

C  SALT -------------------------------------------
C      
C         IF(SALT_FLAG) THEN
C            LASTOMA = LASTOMA*SALT_FACTOR
C         ENDIF
C
C  ------------------------------------------------

            DO I=1,DATA_LEAF(ITYPE) 
               GSTOMA  = A1(I)/(A2(I)+SSTOMA*CCO2ATM)*
     1                   (-1+SQRT(CCO2ATM/(ASTOMA*LASTOMA*VPD)))
            ENDDO
            PSIR    = (KXYLEM*VXYLEM*(PSILEAF+HCANO(ITYPE))
     1                +S1)/(S2+KXYLEM*VXYLEM)
         ENDIF 

      ENDIF

C --------------------------------------------------------------------
C Se PSI_LEAF > LA_PSILMAX oppure PSI_LEAF > 0 impone PSI_LEAF = LA_PSILMAX

      IF (ITER.EQ.999) THEN
         PSILEAF = PSILMAX(ITYPE)
         GXYLEM = 0 

         IF (GROWTH_FLAG) THEN
            KXYLEM   = GXYLEM/(HCANO(ITYPE)*GROWFACTOR)
            DGXYLEM  = GXYLEM*(C_GXYLEM(ITYPE)/D_GXYLEM(ITYPE)*(-
     1                 PSILEAF/D_GXYLEM(ITYPE))**(C_GXYLEM(ITYPE)-1))*
     2                 LAIP
            DKXYLEM  = DGXYLEM/(HCANO(ITYPE)*GROWFACTOR)
         ELSE
            KXYLEM   = GXYLEM/HCANO(ITYPE)
            DGXYLEM  = GXYLEM*(C_GXYLEM(ITYPE)/D_GXYLEM(ITYPE)*(-
     1                 PSILEAF/D_GXYLEM(ITYPE))**(C_GXYLEM(ITYPE)-1))*
     2                 LAI(ITYPE)
            DKXYLEM  = DGXYLEM/HCANO(ITYPE)
         ENDIF

         LASTOMA  = LA_MAX(ITYPE)*CCO2ATM/CCO2STAR*EXP(-LA_BETA(ITYPE)
     1              *(PSILMEAN-LA_PSILMAX(ITYPE))**2)
         DLASTOMA = LA_MAX(ITYPE)*CCO2ATM/CCO2STAR*EXP(-LA_BETA(ITYPE)
     1               *(PSILMEAN-LA_PSILMAX(ITYPE))**2)*(-2*
     2               LA_BETA(ITYPE)*(PSILMEAN-LA_PSILMAX(ITYPE)))

C  SALT -------------------------------------------
C      
C         IF(SALT_FLAG) THEN
C            LASTOMA = LASTOMA*SALT_FACTOR
C            DLASTOMA = DLASTOMA*SALT_FACTOR
C         ENDIF
C         
C  ------------------------------------------------

         E2 = 0
         DE2 = 0
         DO I=1,DATA_LEAF(ITYPE)
             GSTOMA(NP,I)= A1(I)/(A2(I)+SSTOMA*CCO2ATM)*
     1                     (-1+SQRT(CCO2ATM/(ASTOMA*LASTOMA*VPD)))
             DGSTOMA(I)=-0.5*(A1(I)/(A2(I)+SSTOMA*CCO2ATM))*
     1                  (CCO2ATM/(ASTOMA*LASTOMA*VPD))**(-0.5)*CCO2ATM/
     2                  (ASTOMA*VPD*LASTOMA**2)*DLASTOMA

             IF (GROWTH_FLAG) THEN
               E2_LEAF(I)  = ASTOMA*GSTOMA(NP,I)*VPD*EPSW*LAD
     1                       *ACANO(ITYPE)*GROWFACTOR 
               DE2_LEAF(I) = ASTOMA*DGSTOMA(I)*VPD*EPSW*LAD
     1                       *ACANO(ITYPE)*GROWFACTOR
             ELSE
               E2_LEAF(I)  = ASTOMA*GSTOMA(NP,I)*VPD*EPSW*LAD
     1                       *ACANO(ITYPE) 
               DE2_LEAF(I) = ASTOMA*DGSTOMA(I)*VPD*EPSW*LAD
     1                       *ACANO(ITYPE)
             ENDIF

             E2 = E2 + E2_LEAF(I)
             DE2 = DE2 + DE2_LEAF(I)
         ENDDO
         IF (GROWTH_FLAG) THEN
            PSIR  = (KXYLEM*VXYLEM*(PSILEAF+HCANO(ITYPE)*GROWFACTOR)+S1)
     1              /(S2+KXYLEM*VXYLEM)
            DPSIR = ((DKXYLEM*VXYLEM*(PSILEAF+HCANO(ITYPE)*GROWFACTOR)+
     1              KXYLEM*VXYLEM)*(S2+KXYLEM*VXYLEM)-(KXYLEM*VXYLEM*
     2              (PSILEAF+HCANO(ITYPE)*GROWFACTOR)+S1)*(DKXYLEM*
     3              VXYLEM))/(S2+KXYLEM*VXYLEM)**2
            E1    =-KXYLEM*VXYLEM*(PSILEAF-PSIR+HCANO(ITYPE)*GROWFACTOR)
            DE1   =-DKXYLEM*VXYLEM*(PSILEAF-PSIR+HCANO(ITYPE)*
     1              GROWFACTOR)-KXYLEM*VXYLEM*(1-DPSIR)
         ELSE
            PSIR  = (KXYLEM*VXYLEM*(PSILEAF+HCANO(ITYPE))+S1)/
     1              (S2+KXYLEM*VXYLEM)
            DPSIR = ((DKXYLEM*VXYLEM*(PSILEAF+HCANO(ITYPE))+KXYLEM*
     1              VXYLEM)*(S2+KXYLEM*VXYLEM)-(KXYLEM*VXYLEM*(PSILEAF+
     2              HCANO(ITYPE))+S1)*(DKXYLEM*VXYLEM))/
     3              (S2+KXYLEM*VXYLEM)**2
            E1    = -KXYLEM*VXYLEM*(PSILEAF-PSIR+HCANO(ITYPE))
            DE1   = -DKXYLEM*VXYLEM*(PSILEAF-PSIR+HCANO(ITYPE))-KXYLEM*
     1              VXYLEM*(1-DPSIR)
         ENDIF

         E2       = E2*DHCANO 
         DE2      = DE2*DHCANO
         NRFUN    = E2-E1
         DNRFUN   = DE2-DE1

      ENDIF 
C  -------------------- End Newton - Raphson iterations --------------------

C  Print

      IF(SCARTONEW.LE.TOLLNR) THEN
         WRITE(IOUT71,200) ITER
         WRITE(IOUT2,200)  ITER
      ELSE
         WRITE(IOUT71,300) ITER
         WRITE(IOUT2,300)  ITER
      ENDIF
c
C  Format

 100  FORMAT(i10,2x,11(e16.9,2x))
 200  FORMAT('Convergence achieved in',I4,' iterations')
 300  FORMAT('Convergence NOT achieved in',I4,' iterations')

C --------------------------------------------------------------------------

      RETURN
      END            
