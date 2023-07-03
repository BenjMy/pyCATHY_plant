C
C**************************  METEOBAK **********************************
C
C  special handling of time variable meteorological data
C  for the case where BCTIM(1) < TIME <= BCTIM(2)
C  during back-stepping 
C
C***********************************************************************
C
      SUBROUTINE METEOBAK(BCTYPE,
     1                 IPRT1,TIME,BCTIM,BCINP,BC,NMETEODATA)
C
      IMPLICIT     NONE
      INTEGER      I,NMETEODATA
      INTEGER      IOUNI2,IPRT1
      REAL*8       SLOPE,TIM21R,TIMW1
      REAL*8       TIME
      REAL*8       BCTIM(*),BCINP(3,*),BC(*)
      CHARACTER*22 BCTYPE
C      INCLUDE 'PLANT.H'

C
      IF (BCTIM(1) .GE. BCTIM(2)) GO TO 800
C     TIM21R=1.0D0/(BCTIM(2) - BCTIM(1))
C     TIMW1=TIME - BCTIM(1)
      DO I=1,NMETEODATA
C        SLOPE=(BCINP(2,I) - BCINP(1,I))*TIM21R
C        BC(I)=BCINP(1,I) + SLOPE*TIMW1
         BC(I)=BCINP(1,I)
      END DO
c      IF (IPRT1 .GE. 1) THEN
c         WRITE(IOUNI2,1060) BCTYPE,TIME
c         WRITE(IOUNI2,*) '1=Temp 2=RH 3=Radiation'
c         WRITE(1990,*) time,(BC(I),I=1,NMETEODATA)
c         WRITE(IOUNI2,1010) (I,BC(I),I=1,NMETEODATA)
c      END IF
C
  800 RETURN
 1010 FORMAT((1X,3(I5,1X,1PE13.5)))
 1060 FORMAT(2X,A22,' BCs VALUES AT TIME: ',1PE12.4,/,
     1       1X,3('    #  NODE  BC VALUE    '))
      END
