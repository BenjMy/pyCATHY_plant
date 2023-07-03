C
C**************************  MCOPYR ************************************
C
C  copy real matrix IMAT2 into real matrix IMAT1
C
C***********************************************************************
C
      SUBROUTINE MCOPYR(ROW,COL,ROWMAX,COLMAX,IMAT1,IMAT2)
C
      IMPLICIT  NONE
      INTEGER   I,J
      INTEGER   ROW,COL,ROWMAX,COLMAX
      REAL*8    IMAT1(ROWMAX,COLMAX),IMAT2(ROWMAX,COLMAX)
C
      DO I=1,ROW
	   DO J=1,COL
            IMAT1(I,J)=IMAT2(I,J)
	   END DO
      END DO
C
      RETURN
      END

