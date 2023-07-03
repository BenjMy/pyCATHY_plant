c----------------------------------------------------------------------
c
c FUNCTION PLANT_FSO
c
c Calculates the fraction of total dry matter allocated to storage
c organs
c
c----------------------------------------------------------------------

      FUNCTION PLANT_FSO(DVS)

c----------------------------------------------------------------------
      INTEGER I,NDATA
      REAL*8 PLANT_FSO,DVS,X(4),Y(4)
c----------------------------------------------------------------------
      
      NDATA = 4

C DVS
      X(1) = 0.00
      X(2) = 1.10
      X(3) = 1.20
      X(4) = 2.00

C FSO      
      Y(1) = 0.00
      Y(2) = 0.50
      Y(3) = 1.00
      Y(4) = 1.00

      
      IF(DVS.LT.X(1)) THEN
              PLANT_FSO = Y(1)
      ELSEIF(DVS.GT.X(NDATA)) THEN
              PLANT_FSO = Y(NDATA)
      ENDIF

      DO I=1,NDATA-1
         IF((DVS.GE.X(I)).AND.(DVS.LT.X(I+1))) THEN
            PLANT_FSO =(Y(I+1)-Y(I))/(X(I+1)-X(I))*(DVS-X(I))+Y(I)
         ENDIF
      ENDDO

c----------------------------------------------------------------------
      RETURN
      END

