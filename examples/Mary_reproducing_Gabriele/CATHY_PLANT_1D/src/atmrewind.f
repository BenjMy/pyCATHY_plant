C
C**************************  ATREWIND ************************************
C
C  rewinds the file containing the atmospheric boundary conditions 
C  and reads it untill ATMTIM(2)<TIME<= ATMTIM(3)
C
C***********************************************************************
C
      SUBROUTINE ATMREWIND(TIME,ATMTIM,ATMINP,NNOD,HTIATM)
C
      IMPLICIT  NONE
      INCLUDE   'CATHY.H'
      INCLUDE   'IOUNITS.H'
      INTEGER   I
      INTEGER   NNOD,HTIATM, HSPATM_1,IETO_1
      REAL*8    TIME
      REAL*8    ATMTIM(3),ATMINP(3,NNOD)
C
      IF(HTIATM .NE.0) GO TO 800
      REWIND(IIN6)
      READ(IIN6,*,END=700) HSPATM_1,IETO_1
      READ(IIN6,*,END=700) ATMTIM(3)
      READ(IIN6,*) ATMINP(3,1)
      ATMTIM(2)=ATMTIM(3)
      DO WHILE( TIME.GT.ATMTIM(3))
         ATMTIM(2)=ATMTIM(3)
         ATMINP(2,1)=ATMINP(3,1)
         READ(IIN6,*,END=700) ATMTIM(3)
         READ(IIN6,*)ATMINP(3,1)
      END DO
  300 ATMTIM(1)=ATMTIM(2)
      DO I=1,NNOD
         ATMINP(3,I)=ATMINP(3,1)
         ATMINP(2,I)=ATMINP(2,1)
         ATMINP(1,I)=ATMINP(2,1)
      END DO
      GO TO 800
  700 HTIATM=1
      GO TO 300
  800 RETURN
      END
