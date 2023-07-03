C
C --------------------------  PLANT_LEAF  ----------------------------------
C
C  Calculate Leaf water potential PSILEAF (Newton - Raphson method)
C 
C --------------------------------------------------------------------------
C
      SUBROUTINE PLANT_LEAF(N,METEO,SUM1,SUM2,KRS,PNEW,VOLNOD,
     1                     PLNNOD,PLNODES,DELTAT,PSILEAF,LASTOMA,GSTOMA,
     2                    CCI,FC2,E1,E2,NRFUN,DNRFUN,GXYLEM,PSIR,VXYLEM,
     3                     KXYLEM,A1,A2,TIME,PSILMEAN,CUMSUMFC,
     4                     PLANT_FLAG,GROWTH_FLAG,GROWFACTOR,SALT_FLAG,
     5                     SALT_FACTOR,ROOTGROWTH,ABL_FLAG,TSNEW,
     6                     NMETEODATA,ARENOD)

      IMPLICIT NONE

      INCLUDE 'CATHY.H'

C  Input variables
      INTEGER N,PLNODES(PLMAX,PLNODMAX),PLNNOD(PLMAX),NMETEODATA
      REAL*8  DELTAT,KRS(PLMAX,PLNODMAX),SUM1(PLMAX),SUM2(PLMAX)
      REAL*8  PNEW(NMAX),VOLNOD(NMAX),METEO(*)
      REAL*8  TIME,PSILMEAN(PLMAX),GROWFACTOR(PLMAX),SALT3D(NMAX)
      LOGICAL PLANT_FLAG,GROWTH_FLAG,SALT_FLAG
C  Output variables
      REAL*8  GXYLEM(PLMAX),KXYLEM(PLMAX),GSTOMA(PLMAX,DATALEAFMAX)
      REAL*8  PSIR(PLMAX)
      REAL*8  LASTOMA(PLMAX),PSILEAF(PLMAX),E1(PLMAX),E2(PLMAX)
      REAL*8  VXYLEM(PLMAX),ARENOD(NODMAX)
      REAL*8  FC(PLMAX),CCI(PLMAX,DATALEAFMAX),CUMSUMFC(PLMAX)
      REAL*8  NRFUN(PLMAX),DNRFUN(PLMAX)
C  Local variables
      INTEGER I,J,ITYPE
      REAL*8  DGXYLEM(PLMAX),DKXYLEM(PLMAX),DE1(PLMAX),DE2(PLMAX)
      REAL*8  DGSTOMA(PLMAX),DPSIR(PLMAX),DLASTOMA(PLMAX),SIGNGS
      REAL*8  SATVAP,ACTVAP,VPD,EPSW,EPS,UF,FC_LEAF(DATALEAFMAX)    
      REAL*8  VCMAXTEMP,KCTEMP,KOTEMP,COMPTEMP,VCMAX,KC,KO
      REAL*8  A1(DATALEAFMAX),A2(DATALEAFMAX),LAD(PLMAX),DHCANO(PLMAX)
      REAL*8  SALE,SALT_FACTOR(PLMAX),ROOTGROWTH(PLMAX)

      REAL*8  FC2(PLMAX),FC2_LEAF(DATALEAFMAX),COMP(PLMAX)

      REAL*8  TSNEW
      LOGICAL ABL_FLAG

      INCLUDE 'IOUNITS.H' 
      INCLUDE 'PLANT.H'

C --------------------------------------------------------------------------

      DO I = 1,NPLANT
           FC(I)  = 0.D0
           FC2(I) = 0.D0
      ENDDO

C --------------------------------------------------------------------------
C  Calculate Vapour Pressure Deficit + Adjust parameters for temperature
C  Oss: VPD = [molar ratio] -> *0.01      

      IF((ABL_FLAG).AND.(ABL_MODEL.EQ.1)) THEN
         SATVAP = 0.611*EXP((17.502*TSNEW/(TSNEW+240.97)))      
         ACTVAP = 0.611*EXP((17.502*METEO(1)/(METEO(1)+240.97)))*
     1            METEO(2)/100
         VPD    = (SATVAP-ACTVAP)*0.01
      ELSE
         SATVAP = 0.611*EXP((17.502*METEO(1)/(METEO(1)+240.97)))      
         ACTVAP = SATVAP*METEO(2)/100
         VPD    = (SATVAP-ACTVAP)*0.01
      ENDIF

      IF(VPD.LT.0.d0) VPD = 1e-6
      
      EPSW   = 18/RHOW

      VCMAXTEMP = EXP(0.088*(METEO(1)-25))/(1+EXP(0.29*(METEO(1)-41)))
      KCTEMP    = EXP(0.074*(METEO(1)-25))
      KOTEMP    = EXP(0.018*(METEO(1)-25))
      COMPTEMP  = COATM/(2*EXP(-0.056*(METEO(1)-25)))

C ******************* START LOOP ON THE NUMBER OF PLANTS ************
  
      DO I = 1,NPLANT

         ITYPE = ITYPEP(I)  
C         VXYLEM(I) = AXYLEM(ITYPE)*HCANO(ITYPE)

      IF(ACANOFLAG.EQ.1) THEN
         ACANO(ITYPE) = ARENOD(INODP(I))
c         WRITE(9210,*) ACANO(ITYPE)
      ENDIF        

      IF(NMETEODATA.GE.5) THEN
         LAI(ITYPE) = METEO(5)
      ENDIF        

         IF ((PLANT_FLAG).AND.(.NOT.(GROWTH_FLAG))) THEN
            VXYLEM(I) = AXYLEM(ITYPE)*HCANO(ITYPE)
            LAD(I) = LAI(ITYPE)/HCANO(ITYPE)
            DHCANO(I) = HCANO(ITYPE)/DATA_LEAF(ITYPE)
         ELSEIF ((PLANT_FLAG).AND.(GROWTH_FLAG)) THEN
            LAD(I) = LAI2(I)/(HCANO(ITYPE))
c            LAD(I) = LAI2(I)/(HCANO(ITYPE)*GROWFACTOR(I))
            VXYLEM(I) = AXYLEM(ITYPE)*GROWFACTOR(I)*HCANO(ITYPE)*
     1                  GROWFACTOR(I)
            DHCANO(I) = HCANO(ITYPE)*GROWFACTOR(I)/DATA_LEAF(ITYPE)
         ENDIF


C         DHCANO(I) = HCANO(ITYPE)/DATA_LEAF(ITYPE)
         
C  Calculate A1,A2,FC and CUMSUMFC
     
         CALL PLANT_STOMA(ITYPE,VCMAXTEMP,KCTEMP,KOTEMP,COMPTEMP,CCI,
     1                  METEO,DELTAT,A1,A2,FC(I),CUMSUMFC(I),TIME,
     2                  I,LAD(I),DHCANO(I),FC_LEAF,LAI2(I),GROWTH_FLAG,
     3                  SALT_FLAG,SALT_FACTOR(I),COMP(I))

C  Newton-Raphson iterations for the calculation of Psi_leaf

c         WRITE(IOUT71,*) 'PLANT=',I

c         IF(SALT_FLAG) THEN
c            SALE = SALT(1,I)*(-0.165)+1
c         ELSE
            SALE = 1
c         ENDIF
         IF(PLNNOD(I).GT.1) THEN
         CALL PLANT_NR(ITYPE,EPSW,VXYLEM(I),A1,A2,VPD,GXYLEM(I),
     1                KXYLEM(I),DGXYLEM(I),DKXYLEM(I),LASTOMA(I),
     2                DLASTOMA(I),GSTOMA,PSIR(I),DPSIR(I),PSILEAF(I),
     3                SUM1(I),SUM2(I),NRFUN(I),DNRFUN(I),PSILMEAN(I),
     4                METEO,TIME,I,DHCANO(I),LAD(I),LAI2(I),GROWTH_FLAG,
     5                GROWFACTOR(I),SALE,SALT_FLAG,SALT_FACTOR(I),
     6                COMP(I))
         ENDIF      
C  Calculate internal C02 concentration (Fick)
         DO J=1,DATA_LEAF(ITYPE) 
            IF((GSTOMA(I,J)-0.001).GT.0) SIGNGS=1.0
            IF((GSTOMA(I,J)-0.001).EQ.0) SIGNGS=0.0
            IF((GSTOMA(I,J)-0.001).LT.0) SIGNGS=-1.0
            UF = (SIGNGS+1)/2   
c  LINEAR:
c            CCI(I,J)=CCO2ATM-(FC_LEAF(J)/(GSTOMA(I,J)))*UF
C  NON LINEAR:
            CCI(I,J)=CCO2ATM/2.D0+1/(2.D0*GSTOMA(I,J))*
     1               (-A1(J)-A2(J)*GSTOMA(I,J)+SQRT((A1(J)+
     2               (A2(J)-CCO2ATM)*GSTOMA(I,J))**2.D0+
     3               4.D0*GSTOMA(I,J)*(A1(J)*COMP(I)+A2(J)*
     4               CCO2ATM*GSTOMA(I,J)))) 

c
c  Calculate fc2 for crops (water stress) - PROVA ---------
c            IF(CCI(I,J).LT.350.0) CCI(I,J)=350.0

             FC2_LEAF(J)=GSTOMA(I,J)*(CCO2ATM-CCI(I,J))
c            FC2_LEAF(J)=GSTOMA(I,J)*A1(J)*CCO2ATM/
c     1                  (A1(J)+GSTOMA(I,J)*(A2(J)+SSTOMA*CCO2ATM))
             IF (GROWTH_FLAG) THEN
                FC2(I) = FC2(I)+FC2_LEAF(J)*LAD(I)/LAI2(I)
c                FC2(I) = FC2(I)+FC2_LEAF(J)*LAD(I)
             ELSE
                FC2(I) = FC2(I)+FC2_LEAF(J)*LAD(I)/LAI(ITYPE)
c                FC2(I) = FC2(I)+FC2_LEAF(J)*LAD(I)
             ENDIF
c----------------------------------------------------------            
         ENDDO
         IF(GROWTH_FLAG) THEN
            FC2(I)=FC2(I)*DHCANO(I)*LAI2(I)
         ELSE
            FC2(I)=FC2(I)*DHCANO(I)*LAI(ITYPE)
         ENDIF
      ENDDO

      
C ******************** END LOOP ON THE NUMBER OF PLANTS ************* 
c      WRITE(1970,*) TIME,(A1(I),A2(I),I=1,DATA_LEAF(1)),VPD
c      WRITE(1970,*) TIME,GSTOMA(1,15),CCI(1,15)/CCO2ATM,FC(1),FC2(1),
c     1               LAD(1)/LAI2(1),LAD(1),LAI2(1)

C Format ------------------------------------------------------------
 100  FORMAT(i10,2x,11(e16.9,2x))
 200  FORMAT('Convergence achieved in',I4,' iterations')
 300  FORMAT('Convergence NOT achieved in',I4,' iterations')
C -------------------------------------------------------------------

      RETURN
      END            
