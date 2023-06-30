C
C     Subroutine for calculation of the ensemble average
C
      SUBROUTINE AVRG(NMAX,MMAX,N,M,MATR,AVE)

      INTEGER I,J
      INTEGER N,M,NMAX,MMAX
      REAL*8 MATR(NMAX,MMAX)
      REAL*8 AVE(NMAX)

      DO I=1,N
         AVE(I)=0.0d0
         DO J=1,M
            AVE(I)=AVE(I)+MATR(I,J)
         END DO
         AVE(I)=AVE(I)/M
      END DO
      RETURN
      END
