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
      REAL*8 PLANT_RDF,DVS,X(8),Y(8)
c----------------------------------------------------------------------

      NDATA = 8

C DEPTH

      X(1)  = 0.0
      X(2)  = -0.1
      X(3)  = -0.2
      X(4)  = -0.3
      X(5)  = -0.4
      X(6)  = -0.5
      X(7)  = -0.6
      X(8)  = -0.7

C RDF      

      Y(1)  = 6.00*1e+4
      Y(2)  = 2.36*1e+4 
      Y(3)  = 9.26*1e+3 
      Y(4)  = 3.64*1e+3 
      Y(5)  = 1.43*1e+3 
      Y(6)  = 5.62*1e+2 
      Y(7)  = 2.21*1e+2 
      Y(8)  = 8.68*1e+1 

      
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

