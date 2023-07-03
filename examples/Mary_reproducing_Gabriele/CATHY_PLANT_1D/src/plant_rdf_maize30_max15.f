c----------------------------------------------------------------------
c
c FUNCTION PLANT_FSH
c
c Calculates the fraction of total dry matter allocated to shoots
c
c----------------------------------------------------------------------

      FUNCTION PLANT_RDF(DVS)

c----------------------------------------------------------------------
      INTEGER I,NDATA
      REAL*8 PLANT_RDF,DVS,X(4),Y(4)
c----------------------------------------------------------------------

      NDATA = 4

C DEPTH

      X(1)  = 0.0
      X(2)  = -0.1
      X(3)  = -0.2
      X(4)  = -0.3

C RDF      

      Y(1)  = 1.50*1e+4
      Y(2)  = 739 
      Y(3)  = 36
      Y(4)  = 1.8

      
      IF(DVS.GT.X(1)) THEN
              PLANT_RDF = Y(1)
      ELSEIF(DVS.LT.X(NDATA)) THEN
              PLANT_RDF = Y(NDATA)
      ENDIF

      DO I=1,NDATA-1
         IF((DVS.LE.X(I)).AND.(DVS.GT.X(I+1))) THEN
            PLANT_RDF =(Y(I+1)-Y(I))/(X(I+1)-X(I))*(DVS-X(I))+Y(I)
         ENDIF
      ENDDO

c----------------------------------------------------------------------
      RETURN
      END

