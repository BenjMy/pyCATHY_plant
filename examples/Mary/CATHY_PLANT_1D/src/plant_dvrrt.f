c----------------------------------------------------------------------
c
c FUNCTION PLANT_DVRRT
c
c Calculates DeVELOPMENT RATE DVRRT
c
c DVS here is the temperature!!
c----------------------------------------------------------------------

      FUNCTION PLANT_DVRRT(DVS)

c----------------------------------------------------------------------
      INTEGER I,NDATA
      REAL*8 PLANT_DVRRT,DVS,X(3),Y(3)
c----------------------------------------------------------------------

      NDATA = 3
C DVS
      X(1)  = 0.00
      X(2)  = 10.0
      X(3)  = 30.0

C FSH      
      Y(1)  = 0.00
      Y(2)  = 0.00
      Y(3)  = 0.0471

      
      IF(DVS.LT.X(1)) THEN
              PLANT_DVRRT = Y(1)
      ELSEIF(DVS.GT.X(NDATA)) THEN
              PLANT_DVRRT = Y(NDATA)
      ENDIF

      DO I=1,NDATA-1
         IF((DVS.GE.X(I)).AND.(DVS.LT.X(I+1))) THEN
            PLANT_DVRRT =(Y(I+1)-Y(I))/(X(I+1)-X(I))*(DVS-X(I))+Y(I)
         ENDIF
      ENDDO

c----------------------------------------------------------------------
      RETURN
      END

