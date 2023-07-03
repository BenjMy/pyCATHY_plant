c----------------------------------------------------------------------
c
c FUNCTION PLANT_SLA
c
c Calculates SPECIFIC LEAF AREA
c
c----------------------------------------------------------------------

      FUNCTION PLANT_SLA(DVS)

c----------------------------------------------------------------------
      INTEGER I,NDATA
      REAL*8 PLANT_DVRRT,DVS,X(3),Y(3)
c----------------------------------------------------------------------

      NDATA = 3
C DVS
      X(1)  = 0.00
      X(2)  = 0.78
      X(3)  = 2.00

C FSH      
      Y(1)  = 0.026
      Y(2)  = 0.012
      Y(3)  = 0.012

      
      IF(DVS.LT.X(1)) THEN
              PLANT_SLA = Y(1)
      ELSEIF(DVS.GT.X(NDATA)) THEN
              PLANT_SLA = Y(NDATA)
      ENDIF

      DO I=1,NDATA-1
         IF((DVS.GE.X(I)).AND.(DVS.LT.X(I+1))) THEN
            PLANT_SLA =(Y(I+1)-Y(I))/(X(I+1)-X(I))*(DVS-X(I))+Y(I)
         ENDIF
      ENDDO

c----------------------------------------------------------------------
      RETURN
      END

