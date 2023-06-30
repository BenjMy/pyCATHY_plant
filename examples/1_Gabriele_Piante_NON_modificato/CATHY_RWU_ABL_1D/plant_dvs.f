      SUBROUTINE PLANT_DVS(TIME,DELTAT,DAVTMP,ROOTGROWTH)

C --------------------------------------------------------------------------  
C     Calculate phenological DeVelopment Stage (DVS)
C --------------------------------------------------------------------------  

      IMPLICIT NONE 

      INCLUDE 'CATHY.H'
      INCLUDE 'PLANT.H'
      INCLUDE 'IOUNITS.H'

      REAL*8 TIME,DELTAT,DAVTMP
      REAL*8 ERRSH,FSHLV,FSHST,FSHSO,DTSUM
      REAL*8 PLANT_DTEMP,PLANT_DVRRT,PLANT_DVRVT
      REAL*8 PLANT_FSH,PLANT_FLV,PLANT_FST,PLANT_FSO
      REAL*8 FRT_MAX,ROOTGROWTH


c FUNCTION DVRVT(DAVTMP)  ------------------------------
c
c calculate the development rate DVR as a function of 
c daily average temperature (DAVTMP) for the period 
c from emergence till flowering (pre-anthesis, DVRVT) 
c and from flowering till maturity (post-anthesis, DVRRT)


c           DTSUM = PLANT_DTEMP(DAVTMP)
          IF(DAVTMP.LT.0) THEN
                  DVRVT =0.0
                  DVRRT =0.0
          ELSE 
                  DVRVT = PLANT_DVRVT(DAVTMP)
                  DVRRT = PLANT_DVRRT(DAVTMP)
          ENDIF
c------------------------------------------------------

c phenological development stage (DVS) of the crop 
c d-1 ---> s-1 (LDAY)

          IF(DVS-1.LT.0) THEN
c                  DVR = DTSUM/935/LDAY
                  DVR = DVRVT/LDAY
          ELSE
c                  DVR = DTSUM/920/LDAY
                  DVR = DVRRT/LDAY
          ENDIF

c          IF(DAVTMP.LT.0)DVR = 0.0


          DVS = DVS + DVR*DELTAT

C DRY MATTER PARTITIONING
C --------------------------------------------------------------------------  
  
          FSH = PLANT_FSH(DVS)
          FRT = 1-FSH
          FRT_MAX = 1-PLANT_FSH(0)
          ROOTGROWTH = 1-FRT/FRT_MAX

          FLV = PLANT_FLV(DVS)
          FST = PLANT_FST(DVS)
          FSO = PLANT_FSO(DVS)

c for plot
          FSHLV = FSH*FLV
          FSHST = FSH*FST
          FSHSO = FSH*FSO

c check
          ERRSH = ABS(FLV+FST+FSO-1)
c         IF(ERRSH.GT.1e-6) THEN
c            WRITE(IOUT2,*) 'ERROR in DRY MATTER PARTITIONING'
c            WRITE(IOUT2,*) 'ERRSH > 1e-6'
c            STOP
c         ENDIF         

c PLOT
         WRITE(1122,*) TIME,DVS,FSH,FRT,FSHLV,FSHST,FSHSO,ERRSH
C          WRITE(1122,*) TIME, DAVTMP,DTSUM,DVR,DVS
C          write(1919,*) TIME,FRT,FRT_MAX,ROOTGROWTH
C **************************************************************************

      RETURN
 
      END 
