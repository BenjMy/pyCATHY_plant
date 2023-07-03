c----------------------------------------------------------------------
c
c FUNCTION PLANT_FST
c
c Calculates the fraction of total dry matter allocated to stems
c
c----------------------------------------------------------------------

      FUNCTION PLANT_FST(DVS)

c----------------------------------------------------------------------
      INTEGER I,NDATA
      REAL*8 PLANT_FST,DVS,X(8),Y(8)
c----------------------------------------------------------------------
      
      NDATA = 8

C DVS
      X(1) = 0.00
      X(2) = 0.10
      X(3) = 0.25
      X(4) = 0.50
      X(5) = 0.70
      X(6) = 0.95
      X(7) = 1.05
      X(8) = 2.50

C FSH      
      Y(1) = 0.35
      Y(2) = 0.35
      Y(3) = 0.30
      Y(4) = 0.50
      Y(5) = 0.85
      Y(6) = 1.00
      Y(7) = 0.00
      Y(8) = 0.00

      
      IF(DVS.LT.X(1)) THEN
              PLANT_FST = Y(1)
      ELSEIF(DVS.GT.X(NDATA)) THEN
              PLANT_FST = Y(NDATA)
      ENDIF

      DO I=1,NDATA-1
         IF((DVS.GE.X(I)).AND.(DVS.LT.X(I+1))) THEN
            PLANT_FST =(Y(I+1)-Y(I))/(X(I+1)-X(I))*(DVS-X(I))+Y(I)
         ENDIF
      ENDDO

c----------------------------------------------------------------------
      RETURN
      END

