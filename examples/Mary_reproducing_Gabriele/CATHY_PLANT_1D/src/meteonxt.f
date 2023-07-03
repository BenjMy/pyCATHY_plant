C
C**************************  METEONXT  *********************************
C
C  input (if necessary) and interpolate time variable meteorological
C  data
C
C***********************************************************************
C
      SUBROUTINE METEONXT(BCTYPE,IUNIT,IOUNI1,
     1                 IPRT1,HTIBC,TIME,
     1                 BCTIM,BCINP,BC,NMETEODATA)
C
      IMPLICIT     NONE
      INTEGER      I,NMETEODATA
      INTEGER      IUNIT,IOUNI1,IOUNI2,IPRT1,HTIBC
      REAL*8       TIMEIN,SLOPE,TIM32R,TIMW1
      REAL*8       TIME
      REAL*8       BCTIM(*),BCINP(3,*),BC(*)
      CHARACTER*22 BCTYPE
C      INCLUDE 'PLANT.H'
C     INCLUDE 'CATHY.H'
C
      IF (HTIBC .NE. 0) GO TO 800
  200 IF (TIME .LE. BCTIM(3)) GO TO 300
      BCTIM(1)=BCTIM(2)
      BCTIM(2)=BCTIM(3)
      DO I=1,NMETEODATA
         BCINP(1,I)=BCINP(2,I)
         BCINP(2,I)=BCINP(3,I)
      END DO
      IF (IPRT1 .GE. 1) THEN
         WRITE(IOUNI1,1050) BCTYPE,BCTIM(2)
         WRITE(IOUNI1,1055) BCTYPE,BCTIM(2)
         WRITE(IOUNI1,*)'1=Temp 2=RH 3=Radiation'
         WRITE(IOUNI1,1010) (I,BCINP(2,I),I=1,NMETEODATA)
      END IF
      READ(IUNIT,*,END=700) TIMEIN
      BCTIM(3)=TIMEIN
      READ(IUNIT,*)(BCINP(3,I),I=1,NMETEODATA)
      GO TO 200
  300 IF (BCTIM(3) .GT. BCTIM(2)) THEN
C        TIM32R=1.0D0/(BCTIM(3) - BCTIM(2))
C        TIMW1=TIME - BCTIM(2)
         DO I=1,NMETEODATA
C           SLOPE=(BCINP(3,I) - BCINP(2,I))*TIM32R
C           BC(I)=BCINP(2,I) + SLOPE*TIMW1
            BC(I)=BCINP(2,I)
         END DO
      ELSE
         DO I=1,NMETEODATA
            BC(I)=BCINP(3,I)
         END DO
      END IF
c      IF (IPRT1 .GE. 1) THEN
c         WRITE(IOUNI2,1060) BCTYPE,TIME
c         WRITE(IOUNI2,*)'1=Temp 2=RH 3=Radiation'
c         WRITE(1990,*) time,(BC(I),I=1,NMETEODATA)
c         WRITE(IOUNI2,1010) (I,BC(I),I=1,NMETEODATA)
c      END IF
      GO TO 800
C
  700 HTIBC=1
      GO TO 300
C
  800 RETURN
 1010 FORMAT((1X,3(I5,1X,1PE13.5)))
 1050 FORMAT(/,2X,A22,' BCs INPUT AT TIME: ',1PE12.4)
 1055 FORMAT(2X,A22,' BCs  INPUT AT TIME: ',1PE12.4,/,
     1       1X,3('    #  NODE  BC VALUE    '))
 1060 FORMAT(2X,A22,' BCs VALUES AT TIME: ',1PE12.4,/,
     1       1X,3('    #  NODE  BC VALUE    '))
      END
