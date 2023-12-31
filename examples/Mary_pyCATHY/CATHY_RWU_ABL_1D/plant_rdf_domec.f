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
      REAL*8 PLANT_RDF,DVS,X(35),Y(35) 
c----------------------------------------------------------------------

      NDATA = 35

C DEPTH

      X(1)  = 0.0
      X(2)  = -0.14932200
c      X(1)  = -0.14932200
c      X(2)  = -0.18457833
      X(3)  = -0.21983467
      X(4)  = -0.24921500
      X(5)  = -0.30209967
      X(6)  = -0.34910833
      X(7)  = -0.40199300
      X(8)  = -0.46075367
      X(9)  = -0.51951433
      X(10) = -0.58415100
      X(11) = -0.64291167
      X(12) = -0.70167233
      X(13) = -0.76043300
      X(14) = -0.80156567
      X(15) = -0.84269800
      X(16) = -0.88970667
      X(17) = -0.93671500
      X(18) = -0.98959967
      X(19) = -1.04836033
      X(20) = -1.11299700
      X(21) = -1.18351000
      X(22) = -1.24814667
      X(23) = -1.30690733
      X(24) = -1.37742033
      X(25) = -1.41855267
      X(26) = -1.45968533
      X(27) = -1.48906567
      X(28) = -1.53607400
      X(29) = -1.58308267
      X(30) = -1.65359533
      X(31) = -1.71823233
      X(32) = -1.78874500
      X(33) = -1.86513400
      X(34) = -1.92977067
      X(35) = -1.99440733

C RDF      

      Y(1)  = 4*1e+4
      Y(2)  = 0.73512533*1e+4 
c      Y(1)  = 0.73512533*1e+4 
c      Y(2)  = 0.69986900*1e+4 
      Y(3)  = 0.65286033*1e+4 
      Y(4)  = 0.60585200*1e+4 
      Y(5)  = 0.57059533*1e+4 
      Y(6)  = 0.52946300*1e+4 
      Y(7)  = 0.52358700*1e+4 
      Y(8)  = 0.51183467*1e+4 
      Y(9)  = 0.50595867*1e+4 
      Y(10) = 0.51183467*1e+4 
      Y(11) = 0.51771100*1e+4 
      Y(12) = 0.50008267*1e+4 
      Y(13) = 0.48245433*1e+4 
      Y(14) = 0.44719800*1e+4 
      Y(15) = 0.40606567*1e+4 
      Y(16) = 0.37668533*1e+4 
      Y(17) = 0.32967667*1e+4 
      Y(18) = 0.29442033*1e+4 
      Y(19) = 0.28266800*1e+4 
      Y(20) = 0.28266800*1e+4 
      Y(21) = 0.31204833*1e+4 
      Y(22) = 0.32967667*1e+4 
      Y(23) = 0.32967667*1e+4 
      Y(24) = 0.30617233*1e+4 
      Y(25) = 0.25916400*1e+4 
      Y(26) = 0.21215533*1e+4 
      Y(27) = 0.16514667*1e+4 
      Y(28) = 0.12401433*1e+4 
      Y(29) = 0.08288167*1e+4 
      Y(30) = 0.07700567*1e+4 
      Y(31) = 0.09463400*1e+4 
      Y(32) = 0.09463400*1e+4 
      Y(33) = 0.08875800*1e+4 
      Y(34) = 0.07700567*1e+4 
      Y(35) = 0.07112967*1e+4 

      
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

