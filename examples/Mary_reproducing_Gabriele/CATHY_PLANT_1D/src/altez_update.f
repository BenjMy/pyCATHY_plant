C
C************************** ALTEZ_UPDATE *************************************
C This subroutines calculates the water heights after an update, using the
C updated variable VOLUME_KKP1_SN
C***********************************************************************
C

      SUBROUTINE ALTEZ_UPDATE(TIPO_R,DELTA_X,DELTA_Y,
     4                  VOLUME_KKP1_SN,H_WATER_KKP1_SN,
     6                  H_POOL_KKP1_VEC,NUM_R,ELEVATION) 

      IMPLICIT NONE
      INCLUDE 'CATHY.H'
      INCLUDE 'IOUNITS.H'

      INTEGER INDEX_C,AA
      INTEGER TIPO_R(*),NUM_R(*)
      REAL*8  H_WATER_KKP1
      REAL*8  DELTA_X,DELTA_Y,VOLUME_KKP1
      REAL*8  H_WATER_KKP1_SN(*)
      REAL*8  ELEVATION(*)
      REAL*8  H_POOL_KKP1_VEC(*)
      REAL*8  VOLUME_KKP1_SN(*)
C
      DO INDEX_C=1,MAXCEL
 
       IF(TIPO_R(INDEX_C).EQ.0) THEN
          VOLUME_KKP1=VOLUME_KKP1_SN(INDEX_C) 
C   CALCULATION OF WATER HEIGHTS
          IF(VOLUME_KKP1 .GE. 0.D0) THEN
             H_WATER_KKP1=VOLUME_KKP1/(DELTA_X*DELTA_Y)
          ELSE
             VOLUME_KKP1_SN(INDEX_C)=0.0D0
             H_WATER_KKP1=0.0D0
          END IF
          H_WATER_KKP1_SN(INDEX_C)=H_WATER_KKP1
       ELSE IF(TIPO_R(INDEX_C).LT.1000) THEN
          H_WATER_KKP1_SN(INDEX_C)=H_POOL_KKP1_VEC(NUM_R(INDEX_C))-
     &                            ELEVATION(INDEX_C) 
          IF(H_WATER_KKP1_SN(INDEX_C).LT.0) THEN
              H_WATER_KKP1_SN(INDEX_C)=0.0D0
          END IF
       ELSE  
          AA=NUM_R(TIPO_R(INDEX_C)-1000)
          H_WATER_KKP1=H_POOL_KKP1_VEC(AA)-
     &                 ELEVATION(INDEX_C)
          IF(H_WATER_KKP1.GT.0.D0) THEN
             H_WATER_KKP1_SN(INDEX_C)=H_WATER_KKP1
          ELSE 
            H_WATER_KKP1_SN(INDEX_C)=0.0D0
          END IF
        END IF
      END DO
      RETURN
      END
