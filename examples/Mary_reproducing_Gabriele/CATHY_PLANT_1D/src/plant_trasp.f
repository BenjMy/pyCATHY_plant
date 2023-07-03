C
C --------------------------  PLANT_TRASP  ---------------------------------
C
C  Calculate the transpiration flux from each node 
C
C      QPLANT1 = temporary vector of nodal transpiration per plant
C      QPLANT  = nodal transpiration due to all plants
C      QTOT    = total transpiration calculated as sum of nodal values QPLANT
C      TRASP   = total transpiration calculated from PSIR
C 
C --------------------------------------------------------------------------
C
      SUBROUTINE PLANT_TRASP(N,PLNNOD,PLNODES,KRS,PSIR,PNEW,VOLNOD,
     1                      VXYLEM,PSILEAF,KXYLEM,QPLANT,QTOT,
     2                      TRASP,SUM1,SUM2,Z,GROWFACTOR,GROWTH_FLAG,
     3                      QTOT_HR,QTOT_TRASP)

      IMPLICIT NONE

      INCLUDE 'CATHY.H'

C  Input variables
      INTEGER N,PLNNOD(PLMAX),PLNODES(PLMAX,PLNODMAX)
      REAL*8  VOLNOD(NMAX),PNEW(NMAX),ZP,Z(NMAX),KXYLEM(PLMAX)
      REAL*8  KRS(PLMAX,PLNODMAX),PSIR(PLMAX),PSILEAF(PLMAX)
      REAL*8  VXYLEM(PLMAX),GROWFACTOR(PLMAX)
      LOGICAL GROWTH_FLAG
C  Output variables
      REAL*8  QPLANT(NMAX),QTOT(PLMAX),TRASP(PLMAX),QTOT_HR(PLMAX)
      REAL*8  QTOT_TRASP(PLMAX)
C  Local variables
      INTEGER I,J,PLNOD
      REAL*8  QPLANT1(NMAX),SUM1(PLMAX),SUM2(PLMAX),S1(PLMAX),S2(PLMAX)

      INCLUDE 'PLANT.H'

C --------------------------------------------------------------------------

      DO I=1,N
         QPLANT(I) = 0.D0
      ENDDO

      DO I=1,NPLANT
         QTOT(I) = 0.D0
c         S1(I) = 0.D0
c         S2(I) = 0.D0
         QTOT_HR(I)= 0.D0
         QTOT_TRASP(I)=0.D0
      ENDDO         

      DO I = 1,NPLANT

         DO J=1,N
            QPLANT1(J) = 0.D0
         ENDDO
         
         DO J = 1,PLNNOD(I)
            PLNOD = PLNODES(I,J) 
            ZP = ABS(Z(INODP(I))-Z(PLNOD))
            QPLANT1(PLNOD) = -KRS(I,J)*(PSIR(I)-PNEW(PLNOD)+
     1                       ZP)*VOLNOD(PLNOD)

C  Per SPEGNERE la REDISTRIBUTION --------------------------
C            IF (QPLANT1(PLNOD).LT.0.D0) QPLANT1(PLNOD)=0.D0
C  ---------------------------------------------------------

c  Calculate hYDRAULIC rEDISTRIBUTION plant i
            IF(QPLANT1(PLNOD).LT.0) THEN
               QTOT_HR(I)=QTOT_HR(I)+QPLANT1(PLNOD)
            ELSE 
               QTOT_TRASP(I)=QTOT_TRASP(I)+QPLANT1(PLNOD)
            ENDIF
C  ------------------------------------------

            QPLANT(PLNOD) = QPLANT(PLNOD)+QPLANT1(PLNOD)
            QTOT(I) = QTOT(I)+QPLANT1(PLNOD)
            
c            S1(I) = S1(I)+KRS(I,J)*VOLNOD(PLNOD)*
c     1                (PNEW(PLNOD)-ZP)
c            S2(I) = S2(I) +KRS(I,J)*VOLNOD(PLNOD)
         ENDDO
         IF (GROWTH_FLAG) THEN
            TRASP(I) = -KXYLEM(I)*VXYLEM(I)*(PSILEAF(I)-PSIR(I)+
     1                 HCANO(ITYPEP(I))*GROWFACTOR(I))
         ELSE
            TRASP(I) = -KXYLEM(I)*VXYLEM(I)*(PSILEAF(I)-PSIR(I)+
     1                 HCANO(ITYPEP(I)))
         ENDIF
      ENDDO

C --------------------------------------------------------------------------

      RETURN 
      END
