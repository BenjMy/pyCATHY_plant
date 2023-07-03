C
C     Subroutine for calculation of the ensemble average, weighted by 
C     WSIR
C
      SUBROUTINE ENAVRG(NMAX,MMAX,N,M,MATR,AVE,WSIR,RESAMP,ENPT)

      INTEGER I,J
      INTEGER N,M,NMAX,MMAX
      INTEGER ENPT(MMAX)
      LOGICAL RESAMP
      REAL*8 MATR(NMAX,MMAX)
      REAL*8 AVE(NMAX),WSIR(MMAX)
C
      IF (.not.RESAMP) THEN
         DO I=1,N
            AVE(I)=0.0d0
            DO J=1,M
               AVE(I)=AVE(I)+MATR(I,ENPT(J))*WSIR(ENPT(J))
            END DO
            AVE(I)=AVE(I)
         END DO
      ELSE
         DO I=1,N
            AVE(I)=0.0d0
            DO J=1,M
               AVE(I)=AVE(I)+MATR(I,ENPT(J))
            END DO
            AVE(I)=AVE(I)/FLOAT(M)
         END DO
      END IF
      RETURN
      END
