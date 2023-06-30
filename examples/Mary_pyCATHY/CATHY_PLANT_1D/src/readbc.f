C
C**************************  READBC ************************************
C
C  read boundary conditions and if necessary duplicate along the 
C  vertical
C
C***********************************************************************
C
      SUBROUTINE READBC(IUNIT,NB2C,NBFC,NBC,NSTR,BCINP)
C
      IMPLICIT  NONE
      INTEGER   I,J,K1,K2
      INTEGER   IUNIT,NB2C,NBFC,NBC,NSTR
      REAL*8    BCINP(3,*)
C
      IF (NBC.EQ.0.0D0) GO TO 100
      IF (NB2C .NE. 0) THEN
         READ(IUNIT,*) (BCINP(3,I),I=1,NB2C)
         DO I=1,NSTR
            DO J=1,NB2C
               BCINP(3,I*NB2C+J)=BCINP(3,J)
            END DO
         END DO
         IF (NBFC .NE. 0) THEN
            K1=NB2C*(NSTR+1)+1
            K2=NB2C*(NSTR+1)+NBFC
            READ(IUNIT,*) (BCINP(3,I),I=K1,K2)
         END IF
      ELSE
         READ(IUNIT,*) (BCINP(3,I),I=1,NBC)
      END IF
C
 100  RETURN
      END
