C
C************************** DETOUTQ *********************************
C 
C  OUTPUT OF THE DISCHARGE AT THE OUTLET OF THE BASIN AND IN THE OTHER
C  CELLS SUCH AS SELECTED IN PARM INPUT FILE 
C
C***********************************************************************
C
      SUBROUTINE DETOUTQ(NCELNL,NUM_QOUT,TIME,ID_QOUT,QOI_SN,
     1     Q_IN_KKP1_SN,Q_OUT_KKP1_SN_1,Q_OUT_KKP1_SN_2)

      IMPLICIT NONE   
      INCLUDE 'CATHY.H'
      INTEGER I,NUM_QOUT,NCELNL
      INTEGER ID_QOUT(MAXQOUT),QOI_SN(MAXCEL)
      REAL*8  TIME
      REAL*8  Q_OUTPUT(MAXQOUT)
      REAL*8  Q_IN_KKP1_SN(MAXCEL)
      REAL*8  Q_OUT_KKP1_SN_1(MAXCEL)
      REAL*8  Q_OUT_KKP1_SN_2(MAXCEL)

     
      INCLUDE 'IOUNITS.H'

C      WRITING STATEMENT DIFFERENT IF THE CONSIDERED CELL IS THE OUTLET OR NOT.
C      FOR THE INTERNAL CELLS OF THE BASIN WE'RE 'OBLIGED' TO CONSDER 
C      THE Q_IN_KKP1_SN IN CASE THERE IS DISPERSION IN THE FLUXES        

      DO I=1,NUM_QOUT
             Q_OUTPUT(I)=Q_IN_KKP1_SN(ID_QOUT(I))
      END DO
      IF (Q_OUT_KKP1_SN_1(QOI_SN(NCELNL)).NE.0.d0.AND.
     &    Q_OUT_KKP1_SN_2(QOI_SN(NCELNL)).NE.0.d0) THEN
          WRITE(IOUT2,*)'Warning, discharge along 2 directions!'
          STOP
      ELSE
          WRITE(IOUT41,1530)TIME,Q_OUT_KKP1_SN_1(QOI_SN(NCELNL))+
     1                      Q_OUT_KKP1_SN_2(QOI_SN(NCELNL)),
     2                      (Q_OUTPUT(I),I=1,NUM_QOUT)
      END IF
      RETURN
 1530 FORMAT(1PE16.8,5000E16.6E3)
      END
