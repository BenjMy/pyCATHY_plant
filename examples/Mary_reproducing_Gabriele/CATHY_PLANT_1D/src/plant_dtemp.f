c----------------------------------------------------------------------
c
c FUNCTION PLANT_DTEMP
c
c Calculates the daily increase temperature sum
c
c DVS here is the temperature!!
c----------------------------------------------------------------------

      FUNCTION PLANT_DTEMP(DVS)

c----------------------------------------------------------------------
      INTEGER I,NDATA
      REAL*8 PLANT_DTEMP,DVS,X(4),Y(4)
c----------------------------------------------------------------------

      NDATA = 4
C DVS
      X(1)  = 0.00
      X(2)  = 6.00
      X(3)  = 30.0
      X(4)  = 35.0

C FSH      
      Y(1)  = 0.00
      Y(2)  = 0.00
      Y(3)  = 24.0
      Y(4)  = 24.0

      
      IF(DVS.LT.X(1)) THEN
              PLANT_DTEMP = Y(1)
      ELSEIF(DVS.GT.X(NDATA)) THEN
              PLANT_DTEMP = Y(NDATA)
      ENDIF

      DO I=1,NDATA-1
         IF((DVS.GE.X(I)).AND.(DVS.LT.X(I+1))) THEN
            PLANT_DTEMP =(Y(I+1)-Y(I))/(X(I+1)-X(I))*(DVS-X(I))+Y(I)
         ENDIF
      ENDDO

c----------------------------------------------------------------------
      RETURN
      END

