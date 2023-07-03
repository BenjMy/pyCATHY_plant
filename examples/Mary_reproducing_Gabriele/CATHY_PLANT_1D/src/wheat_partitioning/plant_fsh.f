c----------------------------------------------------------------------
c
c FUNCTION PLANT_FSH
c
c Calculates the fraction of total dry matter allocated to shoots
c
c----------------------------------------------------------------------

      FUNCTION PLANT_FSH(DVS)

c----------------------------------------------------------------------
      INTEGER I,NDATA
      REAL*8 PLANT_FSH,DVS,X(14),Y(14)
c----------------------------------------------------------------------

      NDATA = 14
C DVS
      X(1) = 0.00
      X(2) = 0.10
      X(3) = 0.20
      X(4) = 0.35
      X(5) = 0.40
      X(6) = 0.50
      X(7) = 0.60
      X(8) = 0.70
      X(9) = 0.80
      X(10) = 0.90
      X(11) = 1.00
      X(12) = 1.10
      X(13) = 1.20
      X(14) = 2.50

C FSH      
      Y(1) = 0.50
      Y(2) = 0.50
      Y(3) = 0.60
      Y(4) = 0.78
      Y(5) = 0.83
      Y(6) = 0.87
      Y(7) = 0.90
      Y(8) = 0.93
      Y(9) = 0.95
      Y(10) = 0.97
      Y(11) = 0.98
      Y(12) = 0.99
      Y(13) = 1.00
      Y(14) = 1.00

      
      IF(DVS.LT.X(1)) THEN
              PLANT_FSH = Y(1)
      ELSEIF(DVS.GT.X(NDATA)) THEN
              PLANT_FSH = Y(NDATA)
      ENDIF

      DO I=1,NDATA-1
         IF((DVS.GE.X(I)).AND.(DVS.LT.X(I+1))) THEN
            PLANT_FSH =(Y(I+1)-Y(I))/(X(I+1)-X(I))*(DVS-X(I))+Y(I)
         ENDIF
      ENDDO

c----------------------------------------------------------------------
      RETURN
      END

