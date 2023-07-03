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
      REAL*8 PLANT_FSH,DVS,X(12),Y(12)
c----------------------------------------------------------------------

      NDATA = 12
C DVS
      X(1)  = 0.00
      X(2)  = 0.10
      X(3)  = 0.20
      X(4)  = 0.30
      X(5)  = 0.40
      X(6)  = 0.50
      X(7)  = 0.60
      X(8)  = 0.70
      X(9)  = 0.80
      X(10) = 0.90
      X(11) = 1.00
      X(12) = 2.00

C FSH      
      Y(1)  = 0.60
      Y(2)  = 0.63
      Y(3)  = 0.66
      Y(4)  = 0.69
      Y(5)  = 0.73
      Y(6)  = 0.77
      Y(7)  = 0.81
      Y(8)  = 0.85
      Y(9)  = 0.90
      Y(10) = 0.94
      Y(11) = 1.00
      Y(12) = 1.00

      
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

