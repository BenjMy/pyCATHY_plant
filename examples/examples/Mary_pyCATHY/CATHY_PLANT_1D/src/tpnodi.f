C
C**************************  TPNODI ************************************
C
C  ripartisco ai nodi le quantita' che sono date per 
C  elementi (POROS, ELSTOR) e che devono andare ai nodi
C
C***********************************************************************
C
      SUBROUTINE TPNODI(N,NT,NTRI,TETRA,PNODI,SNODI,TP,IPEAT,
     1                  POROS,ELSTOR,PORE,INDE,DEF)
C
      IMPLICIT NONE
      INCLUDE 'CATHY.H'
      INTEGER  K,II,IR,ISTR,MTYPE,INOD,ISTOP,IPEAT
      INTEGER  N,NT,NTRI
      INTEGER  TETRA(5,*),TP(*)
      REAL*8   PNODI(*),SNODI(*),POROS(MAXSTR,*),ELSTOR(MAXSTR,*)
      REAL*8   PORE(*),INDE(*),DEF(*)
      INCLUDE 'IOUNITS.H'
C
      DO K=1,NT
         MTYPE=TETRA(5,K)
         ISTR=1+K/(NTRI*3)
         IR=MOD(K,NTRI*3)
         IF (IR .EQ. 0) ISTR=ISTR-1
         DO II=1,4
            INOD=TETRA(II,K)
            PNODI(INOD)=PNODI(INOD) + POROS(ISTR,MTYPE)
            SNODI(INOD)=SNODI(INOD) + ELSTOR(ISTR,MTYPE)
            TP(INOD)=TP(INOD)+1
         END DO
      END DO
      ISTOP=0
      DO K=1,N
         IF (TP(K) .EQ. 0) THEN
            WRITE(IOUT2,1000) K
            ISTOP=1
         ELSE
            PNODI(K)=PNODI(K)/TP(K)
            SNODI(K)=SNODI(K)/TP(K)
            IF (IPEAT .EQ. 1) THEN
               PORE(K)=PNODI(K)
               INDE(K)=PNODI(K)/(1.0D0-PNODI(K))
               DEF(K)=1.0D0 
            ELSE
               PORE(K)=0.0D0
               INDE(K)=0.0D0
               DEF(K)=0.0D0
            END IF
         END IF
      END DO
      IF (ISTOP .NE. 0) THEN
         CALL CLOSIO
         STOP
      END IF
C
      RETURN
 1000 FORMAT(/,' ERROR IN MESH CONSTRUCTION:    NODE ',I7,
     1         ' IS NOT CONNECTED TO ANY ELEMENTS')
      END
