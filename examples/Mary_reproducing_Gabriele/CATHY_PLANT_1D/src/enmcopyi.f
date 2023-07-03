C
C************************** ENMCOPYI ************************************
C
C  copy integer matrix IMAT2 into integer matrix IMAT1
C
C***********************************************************************
C
      SUBROUTINE ENMCOPYI(ROW,COL,ROWMAX,COLMAX,IMAT1,IMAT2,ENFLAG)
C
      IMPLICIT  NONE
      INTEGER   I,J
      INTEGER   ROW,COL,ROWMAX,COLMAX
      LOGICAL   ENFLAG(COLMAX)
      INTEGER   IMAT1(ROWMAX,COLMAX),IMAT2(ROWMAX,COLMAX)
C
      DO I=1,ROW
	   DO J=1,COL
            IF(ENFLAG(J)) THEN
               IMAT1(I,J)=IMAT2(I,J)
            END IF
	   END DO
      END DO
C
      RETURN
      END

