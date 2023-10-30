C
C**************************  RDNDBC ************************************
C
C  read boundary conditions nodes and if necessary duplicate along the 
C  vertical
C
C***********************************************************************
C
      SUBROUTINE RDNDBC(IUNIT,NB2C,NBFC,NBC,NSTR,BCNOD,NNOD)
C
      IMPLICIT  NONE
      INTEGER   I,J,K1,K2,NODIN2,NODINF
      INTEGER   IUNIT,NSTR,NNOD
      INTEGER   NB2C(3),NBFC(3),NBC(3)
      INTEGER   BCNOD(3,*)
C
      READ(IUNIT,*,END=100) NODIN2,NODINF
      NB2C(3) = NODIN2
      NBFC(3) = NODINF
      NBC(3) = NODIN2*(NSTR+1) + NODINF
      IF (NBC(3).EQ.0.0D0) GO TO 100
      IF (NODIN2 .NE. 0) THEN
         READ(IUNIT,*) (BCNOD(3,I),I=1,NODIN2)
         DO I=1,NSTR
            DO J=1,NODIN2
               BCNOD(3,I*NODIN2+J)=BCNOD(3,J)+I*NNOD 
            END DO
         END DO
         IF (NODINF .NE. 0) THEN
            K1=NODIN2*(NSTR+1)+1
            K2=NODIN2*(NSTR+1)+NODINF
            READ(IUNIT,*) (BCNOD(3,I),I=K1,K2)
         END IF
      ELSE
         READ(IUNIT,*) (BCNOD(3,I),I=1,NBC(3))
      END IF
C
 100  RETURN
      END
