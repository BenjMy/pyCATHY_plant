C
C     Subroutine for calculation of the ensemble average, weighted by 
C     WSIR
C
      SUBROUTINE AVRG3D(NMAX,MMAX,LMAX,N,M,L,MATR,AVE,WSIR,RESAMP,ENPT)

      INTEGER I,J,K
      INTEGER N,M,L,NMAX,MMAX,LMAX
      INTEGER ENPT(LMAX)
      LOGICAL RESAMP
      REAL*8  MATR(NMAX,MMAX,LMAX)
      REAL*8  AVE(NMAX,MMAX),WSIR(LMAX)
C
      IF (.not.RESAMP) THEN
         DO I=1,N
            DO J=1,M
               AVE(I,J)=0.0d0
               DO K=1,L
                  AVE(I,J)=AVE(I,J)+MATR(I,J,ENPT(K))*WSIR(ENPT(K))
               END DO
               AVE(I,J)=AVE(I,J)
            END DO
         END DO
      ELSE
         DO I=1,N
            DO J=1,M
               AVE(I,J)=0.0d0
               DO K=1,L
                  AVE(I,J)=AVE(I,J)+MATR(I,J,ENPT(K))
               END DO
               AVE(I,J)=AVE(I,J)/FLOAT(L)
            END DO
         END DO
      END IF
      RETURN
      END
