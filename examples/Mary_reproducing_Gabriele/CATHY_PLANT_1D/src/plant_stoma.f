C
C --------------------------  PLANT_STOMA  ---------------------------------
C
C  Calculate A1, A2 and FC
C 
C --------------------------------------------------------------------------
C
      SUBROUTINE PLANT_STOMA(ITYPE,VCMAXTEMP,KCTEMP,KOTEMP,COMPTEMP,CCI,
     1                      METEO,DELTAT,A1,A2,FC,CUMSUMFC,TIME,NP,
     2                      LAD,DHCANO,FC_LEAF,LAIP,GROWTH_FLAG,
     3                      SALT_FLAG,SALT_FACTOR,COMP)

      IMPLICIT NONE

      INCLUDE 'CATHY.H'

C  Input variables
      INTEGER ITYPE,I,NP
      REAL*8  VCMAXTEMP,KCTEMP,KOTEMP,COMPTEMP,METEO(*),DELTAT,TIME
      REAL*8  LAD,DHCANO,LAIP,SALT_FACTOR
      LOGICAL GROWTH_FLAG,SALT_FLAG
C  Input/Output variables
      REAL*8  CCI(PLMAX,DATALEAFMAX),A1(DATALEAFMAX),A2(DATALEAFMAX)
      REAL*8  FC,CUMSUMFC
C  Local variables
      REAL*8  VCMAX,KC,KO,COMP
      REAL*8  A1_LIGHT(DATALEAFMAX),A2_LIGHT(DATALEAFMAX)
      REAL*8  A1_RUB(DATALEAFMAX),A2_RUB(DATALEAFMAX)
      REAL*8  FC_LIGHT(DATALEAFMAX),FC_RUB(DATALEAFMAX)
      REAL*8  FC_LEAF(DATALEAFMAX),PAR_ATT,PAR

      INCLUDE 'PLANT.H'

C --------------------------------------------------------------------------

C     Adjust parameters for temperature
      VCMAX = VCMAXTEMP*VCMAX25(ITYPE)
      KC    = KCTEMP*KC25(ITYPE)
      KO    = KOTEMP*KO25(ITYPE)
      COMP  = COMPTEMP/COMP25(ITYPE)


C     Start loop on the # of canopy layers
      DO I=1,DATA_LEAF(ITYPE)

C        Calculate PAR attenuation   
         CALL PLANT_PARATT(METEO(4),LAD,DHCANO,I,PAR_ATT)

C        Calculate A1, A2 and FC assuming RUBISCO limitation
         IF(SALT_FLAG) THEN
           A1_RUB(I) = VCMAX*SALT_FACTOR
         ELSE
           A1_RUB(I) = VCMAX
         ENDIF

         A2_RUB(I) = KC*(1+COATM/KO)
         FC_RUB(I) = A1_RUB(I)/(A2_RUB(I)+CCI(NP,I))*(CCI(NP,I)-COMP)
         IF (FC_RUB(I).LT.0.D0) FC_RUB(I)=0.D0

         PAR = METEO(3)
         IF(PAR.LT.1) PAR=1
C        Calculate A1, A2 and FC assuming LIGHT limitation
         IF(SALT_FLAG) THEN
            A1_LIGHT(I) = LIMIT(ITYPE)*PAR*PAR_ATT*
     1                    SALT_FACTOR
         ELSE
            A1_LIGHT(I) = LIMIT(ITYPE)*PAR*PAR_ATT
         ENDIF
         A2_LIGHT(I) = 2*COMP
         FC_LIGHT(I) = A1_LIGHT(I)/(A2_LIGHT(I)+CCI(NP,I))*
     1                 (CCI(NP,I)-COMP)
         IF (FC_LIGHT(I).LT.0.D0) FC_LIGHT(I)=0.D0

C        Calculate CO2 flux FC_LEAF for the i-layer: 
C        it is the minimum between FC_LIGHT and FC_RUB
         IF((FC_LIGHT(I).EQ.0).AND.(FC_RUB(I).EQ.0)) THEN
            FC_LEAF(I) = 0.D0
         ELSE
            FC_LEAF(I) = MIN(FC_LIGHT(I),FC_RUB(I))
         ENDIF

C        Calculate FC = sum FC_leaf --> somma dei flussi 
C        di carbonio relativi a ciascun layer
         IF (GROWTH_FLAG) THEN
            FC = FC+FC_LEAF(I)*LAD/LAIP
c            FC = FC+FC_LEAF(I)*LAD
         ELSE
            FC = FC+FC_LEAF(I)*LAD/LAI(ITYPE)
c            FC = FC+FC_LEAF(I)*LAD
         ENDIF

C        Select A1 and A2 values according to the kind 
C        of limitation (Light or Rubisco)
         IF (FC_LEAF(I).EQ.FC_LIGHT(I)) THEN
            A1(I) = A1_LIGHT(I)
            A2(I) = A2_LIGHT(I)
         ELSE
            A1(I) = A1_RUB(I)
            A2(I) = A2_RUB(I)
         ENDIF

C     End loop on the # of canopy layers      
      ENDDO

C     Calculate totalCO2flux FC and its cumulative for the 
C     current plant at the current time step

      IF(GROWTH_FLAG) THEN
         FC=FC*DHCANO*LAIP
      ELSE
         FC=FC*DHCANO*LAI(ITYPE)
      ENDIF

c      FC = FC*DHCANO

      IF(GROWTH_FLAG) THEN
c         CUMSUMFC = CUMSUMFC + FC*LAIP*DELTAT
         CUMSUMFC = CUMSUMFC + FC*DELTAT
      ELSE
c         CUMSUMFC = CUMSUMFC + FC*LAI(ITYPE)*DELTAT
         CUMSUMFC = CUMSUMFC + FC*DELTAT
      ENDIF      
C --------------------------------------------------------------------------
c      write(1970,*) TIME,FC,CUMSUMFC,LAD,LAI(ITYPE)
      RETURN
      END

