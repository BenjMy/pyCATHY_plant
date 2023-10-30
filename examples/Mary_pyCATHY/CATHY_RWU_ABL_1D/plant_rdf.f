c----------------------------------------------------------------------
c
c FUNCTION PLANT_RDF
c
c Calculates the Root Distribution Function
c
c----------------------------------------------------------------------

      FUNCTION PLANT_RDF(ZZ,JTYPE)

      IMPLICIT NONE
      INCLUDE 'CATHY.H'
      INCLUDE 'PLANT.H'

c----------------------------------------------------------------------
      INTEGER I,JTYPE
      REAL*8 PLANT_RDF,ZZ
c----------------------------------------------------------------------

      IF(ZZ.GT.ZRDF(1)) THEN
              PLANT_RDF = RDFVAL(1,JTYPE)
      ELSEIF(ZZ.LT.ZRDF(NDATA_RDF)) THEN
              PLANT_RDF = RDFVAL(NDATA_RDF,JTYPE)
      ENDIF

      DO I=1,NDATA_RDF-1
         IF((ZZ.LE.ZRDF(I)).AND.(ZZ.GT.ZRDF(I+1))) THEN
            PLANT_RDF =(RDFVAL(I+1,JTYPE)-RDFVAL(I,JTYPE))/
     1                 (ZRDF(I+1)-ZRDF(I))*
     2                 (ZZ-ZRDF(I))+RDFVAL(I,JTYPE)
         ENDIF
      ENDDO

c----------------------------------------------------------------------
      RETURN
      END

