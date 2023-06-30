C
C************************  PLANT_ROOT  *********************************
C
C  Calculate the Root Distribution Function (According to Vrugt, 2001) 
C
C***********************************************************************
C

      SUBROUTINE PLANT_ROOT(N,X,Y,Z,RDFstamp,VOLNOD,PLNODES,
     1                     PLNNOD,SUMRDF,ROOTGROWTH,GROWTH_FLAG)

      IMPLICIT  NONE
      INCLUDE 'CATHY.H'
C  Input variables
      INTEGER  N
      REAL*8   X(NMAX),Y(NMAX),Z(NMAX),VOLNOD(NMAX),ROOTGROWTH
      LOGICAL  GROWTH_FLAG
C  Output variables
      INTEGER  PLNODES(PLMAX,PLNODMAX),PLNNOD(PLMAX)
      REAL*8   SUMRDF(PLMAX),RDFstamp(NMAX)
C  Local variables
      INTEGER  I,J,K,JNOD,JTYPE
      REAL*8   RDF1,XPLANT,YPLANT,ZP,PLD,AX,AY,AZ,BX,BY,BZ,XM,YM,ZM
      REAL*8   PLANT_RDF
      INCLUDE 'PLANT.H'

C ---------------------------------------------------------------------
C
C Variables:
C
C  JNOD          - Plant position node
C  JTYPE         - Plant type
C  PLNODES(J,K)  - Nodes belonging to the j-plant
C  PLNNOD(J)     - # of nodes belonging to the j-plant  
C  XPLANT        - [L]   - X-distance of the node from the plant
C  YPLANT        - [L]   - Y-distance of the node from the plant  
C  ZPLANT        - [L]   - Z-distanze of the node from the plant      
C  PLDIST        - [L]   - Distance of the node from the plant
C  RDF1          - [-]   - Root Distribution Function in the i-node for a single plant 
C  RDF(J,K)      - [-]   - Cumulative Root Distribution Function = sum of the RDF1 
C                          calculated for the NPLANT plants
C  RDFstamp      - [-]   - Total value of RDF for each node (it is used for printing vtk)
C  SUMRDF        - [L^3] - Integral of RDF over the domain occupied by roots
C  AX,AY,AZ      - Temporary scalars for the calculation of RDF1 
C  BX,BY,BZ      - Temporary scalars for the calculation of RDF1 
C  XMAX,YMAXP,   - [L]   - Max root length in the x,y and z directions,  
C  ZMAXP                   respectively (see Vrugt et al.,2001).
C  VRUX,VRUY,    - [-]   - Empirical parameters for the definition of different root
C  VRUZ                    distributions (see Vrugt et al.,2001).
C  VRUX1,VRUY1   - [L]   - Empirical parameters for the definition of different root
C  VRUZ1                   distributions (see Vrugt et al.,2001).
C
C ---------------------------------------------------------------------

      DO J = 1,NPLANT
         SUMRDF(J) = 0.D0
      ENDDO
      DO J=1,N
         RDFstamp(J) = 0.D0
      ENDDO

C     Loop on the number of plants
      DO J = 1,NPLANT

         K = 0
         JNOD  = INODP(J)
         JTYPE = ITYPEP(J)
 
c            IF(GROWTH_FLAG) THEN
c               XM = XMAXP(JTYPE)
c               YM = YMAXP(JTYPE)
c               ZM = ZMAXP(JTYPE)*(ROOTGROWTH+1e-6)
c            ELSE
               XM = XMAXP(JTYPE)
               YM = YMAXP(JTYPE)
               ZM = ZMAXP(JTYPE)
c            ENDIF
             
C        Loop on the 3D-grid nodes 
         DO I = 1,N
            ZP     = ABS(Z(JNOD)-Z(I))
            XPLANT = ABS(X(I)-X(JNOD))
            YPLANT = ABS(Y(I)-Y(JNOD))
            PLD    = SQRT(XPLANT**2+YPLANT**2+ZP**2)

C           If xplant,yplant or zplant are respectively bigger than 
C           xmax,ymax or zmax the Root Distribution Function is set 
C           equal to zero (by setting =0 the scalars AX,AY or AZ)
            IF (XPLANT.LE.XM) THEN
               AX = 1-XPLANT/XM
            ELSE
               AX = 0.D0
            ENDIF
            IF (YPLANT.LE.YM) THEN
               AY = 1-YPLANT/YM
            ELSE
               AY = 0.D0
            ENDIF
            IF (ZP.LE.ZM) THEN
               AZ = 1-ZP/ZM
            ELSE
               AZ = 0.D0
            ENDIF

            BX = VRUX(JTYPE)/XM*
     1           ABS(VRUX1(JTYPE)-XPLANT)
            BY = VRUY(JTYPE)/YM*
     1           ABS(VRUY1(JTYPE)-YPLANT)
            BZ = VRUZ(JTYPE)/ZM*
     1           ABS(VRUZ1(JTYPE)-ZP)

C           Calculate RDF1 (temporary RDF for the i-node, j-plant)            
            RDF1 = AX*AY*AZ*EXP(-(BX+BY+BZ))
C           Calculate RDF and SUMRDF only for the nodes occupied 
C           by roots, and fix the indices of these nodes in PLNODES
            IF(RDF1.NE.0) THEN
               K = K+1           
               PLNODES(J,K) = I 
               SUMRDF(J)    = SUMRDF(J)+RDF1*VOLNOD(I)            
               RDF1 = PLANT_RDF(-ZP,JTYPE)
c               RDF1 = PLANT_RDF(Z(I))
c               RDF1 = RDF1*(XMAXP(JTYPE)-SQRT(XPLANT**2+YPLANT**2))
c     1                /XMAXP(JTYPE)
               IF (RDF1.LT.0) RDF1=0.D0
            ENDIF
            


C TOGLI..Intanto si usa per la stampa in vtk della distribuzione di radici  
            RDFstamp(I) = RDFstamp(I)+RDF1                   
C ------------------------------------------------------------------------

         ENDDO 

         PLNNOD(J)=K

      ENDDO
c      WRITE(9090,*) ((ROOTGROWTH+1E-6), (ROOTGROWTH+1E-6)*ZMAXP(1))
c      write(9191,*) (PLNNOD(J),J=1,NPLANT)

      RETURN
      END
