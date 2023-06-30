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
      REAL*8 PLANT_FST,DVS,X(7),Y(7)
c----------------------------------------------------------------------
      
      NDATA = 7

C DVS
      X(1) = 0.00
      X(2) = 0.33
      X(3) = 0.88
      X(4) = 0.95
      X(5) = 1.10
      X(6) = 1.20
      X(7) = 2.00

C FSH      
      Y(1) = 0.38
      Y(2) = 0.38
      Y(3) = 0.85
      Y(4) = 0.85
      Y(5) = 0.40
      Y(6) = 0.00
      Y(7) = 0.00

      
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

