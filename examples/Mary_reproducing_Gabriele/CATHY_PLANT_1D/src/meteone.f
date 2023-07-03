C
C**************************  METEONE  **********************************
C
C  read and initialize time variable meteorological data
C  METEO(1)  = Air Temperature [degC]
C  METEO(2)  = Relative Humidity [%]
C  METEO(3)  = PAR
C  METEO(4)  = Zenith Angle
C  METEO(5)  = LAI
C  METEO(6)  = Net Radiation [W/m2]
C  METEO(7)  = Shortwave incoming radiation [W/m2]
C  METEO(8)  = Wind speed [m/s]
C  METEO(9)  = G = Soil Heat Flux [W/m2]
C  METEO(10) = Sensible Heat Flux (uso solo per stampare DATO agli stessi tempi del modello)
C  METEO(11) = ..used to be the Friction Velocity...not needed anymore! :)
C  METEO(12) = Atmospheric Pressure
C***********************************************************************
C
      SUBROUTINE METEONE(BCTYPE,IUNIT,IOUNI1,
     1                   IPRT1,HTIBC,TIME,
     2                   DELTAT,BCTIM,BCINP,BC,NMETEODATA)
c      SUBROUTINE BCONE(BCTYPE,IUNIT,IOUNI1,IOUNI2,NNOD,
c     1                 NBCMAX,IPRT1,NB2C,NBFC,NBC,NSTR,HTIBC,TIME,
c     2                 DELTAT,BCTIM,BCINP,BCNOD,BC,ANBC,ABCNOD)
C
      IMPLICIT     NONE
      INTEGER      I,J,NMETEODATA
      INTEGER      IUNIT,IOUNI1,IOUNI2,IPRT1,HTIBC
      REAL*8       TIMEIN,SLOPE,TIM32R,TIMW1
      REAL*8       TIME,DELTAT
      REAL*8       BCTIM(*),BCINP(3,*),BC(*)
      CHARACTER*22 BCTYPE
c      INCLUDE 'PLANT.H'
C     INCLUDE 'CATHY.H'
C
      HTIBC=0
      BCTIM(1)=0.0D0
      BCTIM(2)=0.0D0
      DO J=1,3
         DO I=1,NMETEODATA
            BCINP(J,I)=0.0D0
         END DO
      END DO
      READ(IUNIT,*) TIMEIN
      IF (DELTAT .GE. 1.0D+15) TIMEIN=0.0D0
      BCTIM(3)=TIMEIN
      READ(IUNIT,*)(BCINP(3,I),I=1,NMETEODATA)
      IF (IPRT1 .GE. 1) THEN
         WRITE(IOUNI1,1000) BCTYPE
         WRITE(IOUNI1,*) '1=Temp 2=RH 3=Radiation'
         WRITE(IOUNI1,1010) (I,BCINP(3,I),I=1,NMETEODATA)
      END IF
      IF (DELTAT .GE. 1.0D+15) GO TO 300
  200 IF (TIME .LE. BCTIM(3)) GO TO 300
      BCTIM(1)=BCTIM(2)
      BCTIM(2)=BCTIM(3)
      DO I=1,NMETEODATA 
            BCINP(1,I)=BCINP(2,I)
            BCINP(2,I)=BCINP(3,I)
      END DO
      IF (IPRT1 .GE. 1) THEN
         WRITE(IOUNI1,1050) BCTYPE,BCTIM(2)
c         WRITE(IOUNI2,1055) BCTYPE,BCTIM(2)
c         WRITE(IOUNI2,*) '1=Temp 2=RH 3=Radiation'
c         WRITE(1990,*) time,(BCINP(2,I),I=1,NMETEODATA)
c         WRITE(IOUNI2,1010) (I,BCINP(2,I),I=1,NMETEODATA)
      END IF
      READ(IUNIT,*,END=700) TIMEIN
      BCTIM(3)=TIMEIN
      READ(IUNIT,*)(BCINP(3,I),I=1,NMETEODATA)
      GO TO 200
  300 IF (BCTIM(3) .GT. BCTIM(2)) THEN
C        TIM32R=1.0D0/(BCTIM(3) - BCTIM(2))
C        TIMW1=TIME - BCTIM(2)
C           SLOPE=(BCINP(3,I) - BCINP(2,I))*TIM32R
C           BC(I)=BCINP(2,I) + SLOPE*TIMW1
            DO I=1,NMETEODATA
            BC(I)=BCINP(2,I)
            END DO
      ELSE
         DO I=1,NMETEODATA
            BC(I)=BCINP(3,I)
         END DO
      END IF
      IF (IPRT1 .GE. 1) THEN
         WRITE(IOUNI2,1060) BCTYPE,TIME
         WRITE(IOUNI2,1010) (I,BC(I),I=1,NMETEODATA)
      END IF
      GO TO 800
C
  700 HTIBC=1
      GO TO 300
C
  800 RETURN
 1000 FORMAT(/,2X,A22,' BCs AT BEGINNING OF SIMULATION',/,
     1       1X,3('    #  NODE  BC VALUE    '))
 1010 FORMAT((4(I5,1X,1PE13.5)))
 1050 FORMAT(/,2X,A22,' BCs INPUT AT TIME: ',1PE12.4)
 1055 FORMAT(2X,A22,' BCs  INPUT AT TIME: ',1PE12.4,/,
     1       1X,3('    #  NODE  BC VALUE    '))
 1060 FORMAT(2X,A22,' BCs VALUES AT TIME: ',1PE12.4,/,
     1       1X,3('    #  NODE  BC VALUE    '))
      END
