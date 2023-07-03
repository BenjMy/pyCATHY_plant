C
C**************************  MCOPYI ************************************
C
C  copy integer matrix IMAT2 into integer matrix IMAT1
C
C***********************************************************************
C
      SUBROUTINE MCOPYI(ROW,COL,ROWMAX,COLMAX,IMAT1,IMAT2)
C
      IMPLICIT  NONE
      INTEGER   I,J
      INTEGER   ROW,COL,ROWMAX,COLMAX
      INTEGER   IMAT1(ROWMAX,COLMAX),IMAT2(ROWMAX,COLMAX)
C
      DO I=1,ROW
	   DO J=1,COL
            IMAT1(I,J)=IMAT2(I,J)
	   END DO
      END DO
C
      RETURN
      END

