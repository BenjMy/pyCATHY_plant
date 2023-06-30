      SUBROUTINE PLANT_SALT(NNOD,SALT3D)

      IMPLICIT NONE

      INCLUDE 'CATHY.H'
      INCLUDE 'PLANT.H'

      INTEGER I,J,K,A,STRATO,NNOD
      REAL*8  SALT3D(NMAX) 

      DO I=1,NNOD
         SALT3D(I)=SALT(1,I)
         WRITE(1133,*) I,SALT(1,I),SALT3D(I)
      ENDDO

      STRATO=0

      DO I=1,NSTRSALT
        DO J=1,STRSALT(I)
           STRATO=STRATO+1
            DO K=1,NNOD
               A=K+STRATO*NNOD
               SALT3D(A)=SALT(I,K)
               WRITE(1133,*) A,SALT3D(A)
            ENDDO
        ENDDO
      ENDDO

  
      RETURN
      END
