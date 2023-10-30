C
C --------------------------  PLANT_COND  ----------------------------------
C
C  Calculate Root - Soil conductances
C 
C --------------------------------------------------------------------------
C

      SUBROUTINE PLANT_COND(N,NT,TP,TETRA,KS,CKRW,VOLNOD,SUMRDF,
     1                      PLNNOD,PLNODES,PNEW,KRS,X,Y,Z,
     2                      SUM1,SUM2,SW,PNODI,IVGHU,SALT_FLAG,SALT3D,
     3                      ROOTGROWTH,GROWTH_FLAG,QPLANT,PSIR,TIME)

      IMPLICIT NONE 

      INCLUDE 'CATHY.H'
      INCLUDE 'SOILCHAR.H'

C  Input variables
      INTEGER N,NT,TP,TETRA(5,*),PLNNOD(PLMAX),PLNODES(PLMAX,PLNODMAX)
      INTEGER IVGHU
      REAL*8  KS(*),CKRW(NMAX),PNEW(NMAX),VOLNOD(NMAX)
      REAL*8  RDF,TIME
      REAL*8  SUMRDF(PLMAX),SW(NMAX),PNODI(NMAX)
      REAL*8  SALT3D(NMAX),ROOTGROWTH,PSIR(PLMAX)
      LOGICAL SALT_FLAG,GROWTH_FLAG
C  Output variables
      REAL*8  GRS,KRS(PLMAX,PLNODMAX)
      REAL*8  SUM1(PLMAX),SUM2(PLMAX)
C  Local variables
      INTEGER I,J,K,PLNOD
      INTEGER JNOD,JTYPE
      REAL*8  KSNOD(NMAX),KSOIL,GSOIL,GROOT,NRDF
      REAL*8  L,F1,F2,F,SWR,XP,ZP,YP,PLDIST
      REAL*8  X(NMAX),Y(NMAX),Z(NMAX),QPLANT(NMAX)
      REAL*8  AX,AY,AZ,BX,BY,BZ,XM,YM,ZM,PLANT_RDF,RAI,DISTSR
    
      INCLUDE 'PLANT.H'

C --------------------------------------------------------------------------

      
C     Inizialization
      NRDF  = 0.D0
      GSOIL = 0.D0
      GROOT = 0.D0
      GRS  = 0.D0

      DO J=1,NPLANT
         SUM1(J)      = 0.D0
         SUM2(J)      = 0.D0
         DO I=1,PLNNOD(J)
            KRS(J,I)  = 0.D0
         ENDDO
      ENDDO

C --------------------------------------------------------------------------

C     Start loop on the # of plants
      DO J = 1,NPLANT

        JNOD = INODP(J)
        JTYPE = ITYPEP(J)

C*****************************************
C        IF (GROWTH_FLAG) THEN
C           XM = XMAXP(JTYPE)
C           YM = YMAXP(JTYPE)
C           ZM = ZMAXP(JTYPE)*(ROOTGROWTH+1e-6)
C        ELSE
C           XM = XMAXP(JTYPE)
C           YM = YMAXP(JTYPE)
C           ZM = ZMAXP(JTYPE)
C        ENDIF
C*****************************************

C       Start loop on the # of nodes related to the current plant
        DO I = 1,PLNNOD(J)

C          Calculate the same parameters calculated in plant_root
           PLNOD = PLNODES(J,I)

C*****************************************
           ZP    = ABS(Z(INODP(J))-Z(PLNOD))
           XP    = ABS(X(INODP(J))-X(PLNOD))
           YP    = ABS(Y(INODP(J))-Y(PLNOD))
           PLDIST= SQRT(XP**2+YP**2+ZP**2)
C
C           IF (XP.LE.XM) THEN
C              AX = 1-XP/XM
C           ELSE
C              AX = 0.D0
C           ENDIF
C           IF (YP.LE.YM) THEN
C              AY = 1-YP/YM
C           ELSE
C              AY = 0.D0
C           ENDIF
C           IF (ZP.LE.ZM) THEN
C              AZ = 1-ZP/ZM
C           ELSE
C              AZ = 0.D0
C           ENDIF
C
C           BX = VRUX(JTYPE)/XM*
C     1          ABS(VRUX1(JTYPE)-XP)
C           BY = VRUY(JTYPE)/YM*
C     1          ABS(VRUY1(JTYPE)-YP)
C           BZ = VRUZ(JTYPE)/ZM*
C     1          ABS(VRUZ1(JTYPE)-ZP)
C
C           RDF = AX*AY*AZ*EXP(-(BX+BY+BZ))
C************************************************


           IF(GROWTH_FLAG) THEN
              RDF = PLANT_RDF(-ZP,JTYPE)*(ROOTGROWTH+1E-6)
c              RDF = PLANT_RDF(Z(PLNOD))*(ROOTGROWTH+1E-6)
           ELSE
              RDF = PLANT_RDF(-ZP,JTYPE)
c              RDF = PLANT_RDF(Z(PLNOD))
           ENDIF
C           RDF = RDF*(XMAXP(JTYPE)-SQRT(XP**2+YP**2))/XMAXP(JTYPE)
           IF (RDF.LT.0) RDF=0.D0
           CALL ELTNOD(N,NT,TP,TETRA,KS,KSNOD)

c           IF((PSIR(J)-PNEW(PLNOD)+ZP).LT.0.D0) THEN
              KSOIL = KSNOD(PLNOD)*CKRW(PLNOD)
c           ELSE
c              KSOIL = KSNOD(PLNOD)
c           ENDIF

C          Calculate RAI - Domec Site Only!!
C           IF(Z(PLNOD).GT.-0.3) THEN
C             RAI = 4.7
C           ELSEIF((Z(PLNOD).GE.-0.6).AND.(Z(PLNOD).LE.-0.3))THEN
C             RAI = 3.6
C           ELSE
C             RAI = 2.1 
C           ENDIF

c           PLDIST = SQRT(ZMAXP(JTYPE)*2*0.001/RAI)
c           DISTSR = SQRT(2*0.001/RDF)
           DISTSR =SQRT(1/3.14/RDF)
           PLDIST = DISTSR
c           PLDIST = DISTSR*ZP/ZMAXP(ITYPEP(J))
c           PLDIST = DISTSR*(1-RDF/7351.3)
c           IF(Z(PLNOD).GE.-0.2) PLDIST = ABS(Z(PLNOD))

C Fino al 12/06/13:
c           GSOIL = KSOIL*RDF*2*3.14*0.001/RAI
c Ora:
           GSOIL = KSOIL/PLDIST

c           GROOT = 0.0001*18/1e11/(2*3.14*0.001)*RDF*VOLNOD(PLNOD)
c           GROOT = GROOT_STAR(ITYPEP(J))*18/1e11/(2*3.14*0.001)*RDF*
c     1             VOLNOD(PLNOD)*EXP(-(PSIR(J)/(-120.d0))**2.d0)

C Fino al 12/06/13:
c           GROOT = GROOT_STAR(ITYPEP(J))*18/1e11/(2*3.14*0.001)*RDF*
c     1             VOLNOD(PLNOD)
c Ora:
           GROOT = GROOT_STAR(ITYPEP(J))


c           WRITE(3457,*) KSNOD(PLNOD),CKRW(PLNOD),KSOIL
C************************************************
C           NRDF  = RDF*VOLNOD(PLNOD)/SUMRDF(J)
C           L     = (VOLNOD(PLNOD)**(1/3))*SQRT(3.d0)
C
CC          Calculate Soil conductance GSOIL:          
C           GSOIL = GSOIL_STAR(ITYPEP(J))*NRDF/L
Cc          Calculate Root conductance GROOT:
C           GROOT = NRDF*GROOT_STAR(ITYPEP(J))
   
CC          Flux is limited to wetted roots using effective saturation
CC          Only Van Genuchten case is implemented -for the moment-
C           IF (IVGHU.EQ.0) THEN
C              SWR = VGRMC/PNODI(PLNOD)
C              F1 = (SW(PLNOD)-SWR)/(1-SWR)
C           ELSE
C              F1 = SW(PLNOD)    
C           ENDIF
C
C           F = F1

C          Calculate ROOT-SOIL conductance:
C          If GROOT = GSOIL = 0          ---> GRS=0 (No trasp flux from this node) 
           IF((GROOT.EQ.0).AND.(GSOIL.EQ.0)) THEN
              GRS = 0.D0
              KRS(J,I) = 0.D0
C          If GSOIL = 0 and GROOT not= 0 ---> GRS=GROOT
           ELSEIF ((GSOIL.EQ.0).AND.(GROOT.NE.0)) THEN
              GRS= GROOT
c              GRS=1E-7
c                If the node distance from the plant is =0 (i.e. the 
C                node coincides with the plant position node) Trasp flux 
C                from that node is put =0 (by imposing KRS=0)
                 IF (PLDIST.EQ.0) THEN
                    KRS(J,I) = 0.D0
                 ELSE
c Fino al 12/06/13:
c                    KRS(J,I) = GRS/PLDIST
c Ora:
                   KRS(J,I) = 2*3.14*0.001*RDF*GRS 

c                    KRS(J,I) = GRS*1E-7/(PLDIST*GRS+1E-7*DISTSR)
C*******************************************
c                    KRS(J,I) = GRS*F/PLDIST
C*******************************************
                 ENDIF
C          If GSOIL not= 0 and GROOT = 0 ---> GRS=GSOIL
           ELSEIF ((GSOIL.NE.0).AND.(GROOT.EQ.0)) THEN
              GRS= GSOIL
c              GRS=1E-7
c                If the node distance from the plant is =0 (i.e. the 
C                node coincides with the plant position node) Trasp flux 
C                from that node is put =0 (by imposing KRS=0)
                 IF (PLDIST.EQ.0) THEN
                    KRS(J,I) = 0.D0
                 ELSE
c Fino al 12/06/13:
c                    KRS(J,I) = GRS/PLDIST
c Ora:
                   KRS(J,I) = 2*3.14*0.001*RDF*GRS 
c                    KRS(J,I) = GRS*1E-7/(PLDIST*GRS+1E-7*DISTSR)
C*******************************************
c                    KRS(J,I) = GRS*F/PLDIST
C*******************************************
                 ENDIF

           ELSE
C          If GSOIL not= 0 and GROOT not= 0 ---> GRS is calculated as 
C          the series of GSOIL and GROOT
              GRS = GSOIL*GROOT/(GSOIL+GROOT)
c              GRS=1E-7
c                If the node distance from the plant is =0 (i.e. the 
C                node coincides with the plant position node) Trasp flux 
C                from that node is put =0 (by imposing KRS=0)
                 IF (PLDIST.EQ.0) THEN
                    KRS(J,I) = 0.D0
                 ELSE
c Fino al 12/06/13:
c                    KRS(J,I) = GRS/PLDIST
c Ora:
                   KRS(J,I) = 2*3.14*0.001*RDF*GRS 
c                    KRS(J,I) = GRS*1E-7/(PLDIST*GRS+1E-7*DISTSR)
C*******************************************
c                    KRS(J,I) = GRS*F/PLDIST
C*******************************************
                 ENDIF
            ENDIF
C         Calculate SUM1 and SUM2 (to be used in subroutine plant_nr)

          SUM1(J) = SUM1(J) + KRS(J,I)*VOLNOD(PLNOD)*
     1              (PNEW(PLNOD)-ZP)
          SUM2(J) = SUM2(J) + KRS(J,I)*VOLNOD(PLNOD)

c           WRITE(3456,*) TIME,GSOIL,GROOT,GRS,Z(PLNOD)

C       End loop on the # of nodes belonging to the current plant
        ENDDO

C     End loop on the # of plants       
      ENDDO

C --------------------------------------------------------------------------

      RETURN
 
      END 
