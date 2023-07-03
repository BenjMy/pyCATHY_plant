       SUBROUTINE PLANT_PARATT(ZEN,LAD,DHCANO,J,PAR_ATT)

       IMPLICIT NONE

       INCLUDE 'CATHY.H'
       INCLUDE 'PLANT.H'

       INTEGER  I,J
       REAL*8   X,CLUMP,LAD,ZLEAF,CUMSUM,XN1,XD1,EX1
       REAL*8   PAR_ATT,DHCANO,ZEN
c Reference: Plants and Microclimate: a quantitative approach to 
C            environmental plant physiology - H.J. Jones 1992

C x = 1        : Spherical leaf angle distribution
C x = INFINITY : Horizontal leaf angle distribution
C x = 0        : Vertical leaf angle distribution           
c       X = 0
       X = 1e+20
       CLUMP = 0.85 
       ZLEAF = J* DHCANO
       CUMSUM = LAD*DHCANO*J
       XN1 = SQRT(X*X+(TAN(ZEN))**2)
       XD1 = X+1.774*(X+1.182)**(-0.733)
       EX1 = XN1/XD1
       PAR_ATT = EXP(-EX1*CUMSUM*CLUMP)

       RETURN
       END
