C
C**************************  SWITCH ************************************
C
C  switching control of atmospheric boundary conditions
C
C***********************************************************************
C
      SUBROUTINE SWITCH(NNOD,IFATM,PONDING,TIME,DELTAT,PONDH_MIN,
     1                  ARENOD,PONDNOD,ATMPOT,ATMACT,PNEW,OVFLNOD)
C
      IMPLICIT  NONE
      INCLUDE  'CATHY.H'
      INTEGER   I
      INTEGER   NNOD
      INTEGER   IFATM(NNOD)
      LOGICAL   PONDING
      REAL*8    PL,ATMDIF
      REAL*8    TIME,DELTAT,PONDH_MIN
      REAL*8    ZERO
      REAL*8    ARENOD(NNOD),PONDNOD(NNOD),ATMPOT(NNOD),ATMACT(NNOD)
      REAL*8    PNEW(NNOD),OVFLNOD(NNOD)
      INCLUDE  'IOUNITS.H'
      INCLUDE  'SOILCHAR.H'
      PARAMETER (ZERO = 0.0)
C
      ponding = .false.

cm    write(99,*)'nodo atmpot       atmact       pondnod      pnew',
cm   1'         pl            ifatm'
      DO I=1,NNOD
         IF (IFATM(I) .EQ. -1) THEN
cc            PONDING = .FALSE.
            OVFLNOD(I) = ZERO
            GO TO 500
         END IF
C
         ATMDIF = ATMPOT(I) - ATMACT(I)
         PL = PONDNOD(I) + (ATMDIF * DELTAT / ARENOD(I))
cm       if (i.eq.264) then
cm          write(666,*) 'prima',time,atmpot(i),atmact(i),pl,pnew(i),
cm   1                    ifatm(i)
cm       end if
cm       write(99,1000) i,atmpot(i),atmact(i),
cm   1                              pondnod(i),pnew(i),pl,
cm   2                              ifatm(i)
cm 1000 FORMAT(i4,5e13.5,i4)
         IF (IFATM(I) .EQ. 2) THEN
            IF (ATMPOT(I) .GE. ZERO) THEN
               IF (ATMACT(I) .GE. ZERO) THEN
C
C  ponding, rainfall, infiltration 
C  (case A)
C
                  IF (PL .GE. PONDH_MIN) THEN
                     PONDING = .TRUE.
                     IFATM(I) = 2
                     PNEW(I) = PL
                     OVFLNOD(I) = ATMDIF
                     GO TO 500
                  END IF
                  IF (PL .GE. ZERO .AND. PL .LT. PONDH_MIN) THEN
cc                     PONDING = .FALSE.
                     IFATM(I) = 1
C  need to decide how to treat this case: set PNEW to 0, to PL, or
C  leave it as is??
C+                   PNEW(I) = ZERO
C+                   PNEW(I) = PL
                     OVFLNOD(I) = ATMDIF
                     GO TO 500
                  END IF
                  IF (PL .GT. PMIN .AND. PL .LT. ZERO) THEN
cc                     PONDING = .FALSE.
                     IFATM(I) = 0
                     ATMACT(I) = ATMPOT(I)
C  need to decide how to treat this case: set PNEW to 0 or leave it
C  as is??
C+                   PNEW(I) = ZERO
                     OVFLNOD(I) = -PONDNOD(I) * ARENOD(I) / DELTAT
                     GO TO 500
                  END IF
                  IF (PL .LE. PMIN) THEN
cc                     PONDING = .FALSE.
                     IFATM(I) = 0
                     ATMACT(I) = ATMPOT(I)
C  need to decide how to treat this case: set PNEW to 0 or leave it
C  as is??
C+                   PNEW(I) = ZERO
                     OVFLNOD(I) = -PONDNOD(I) * ARENOD(I) / DELTAT
                     GO TO 500
                  END IF
               ELSE
C
C  ponding, rainfall, exfiltration (==> return flow or seepage) 
C  (case B)
C
                  IF (PL .GE. PONDH_MIN) THEN
                     PONDING = .TRUE.
                     IFATM(I) = 2
                     PNEW(I) = PL
                     OVFLNOD(I) = ATMDIF
                     GO TO 500
                  END IF
                  IF (PL .GE. ZERO .AND. PL .LT. PONDH_MIN) THEN
cc                     PONDING = .FALSE.
                     IFATM(I) = 1
C  need to decide how to treat this case: set PNEW to 0, to PL, or
C  leave it as is??
C+                   PNEW(I) = ZERO
C+                   PNEW(I) = PL
                     OVFLNOD(I) = ATMDIF
                     GO TO 500
                  END IF
                  IF (PL .LT. ZERO) THEN
cc                     PONDING = .FALSE.
                     IFATM(I) = 1
C  need to decide how to treat this case: set PNEW to 0 or leave it
C  as is??
C+                   PNEW(I) = ZERO
                     OVFLNOD(I) = ATMDIF
cxcx                 WRITE(IOUT9,9000) I,PL,TIME,DELTAT
                     GO TO 500
                  END IF
               END IF
            ELSE
               IF (ATMACT(I) .GE. ZERO) THEN
C
C  ponding, evaporation, infiltration (typical scenario: surface is
C  saturated and rainfall has just switched to evaporation)
C  (case C)
C
                  IF (PL .GE. PONDH_MIN) THEN
                     PONDING = .TRUE.
                     IFATM(I) = 2
                     PNEW(I) = PL
                     OVFLNOD(I) = ATMDIF
                     GO TO 500
                  END IF
                  IF (PL .GE. ZERO .AND. PL .LT. PONDH_MIN) THEN
cc                     PONDING = .FALSE.
                     IFATM(I) = 1
C  need to decide how to treat this case: set PNEW to 0, to PL, or
C  leave it as is??
                     PNEW(I) = ZERO
C+                   PNEW(I) = PL
                     OVFLNOD(I) = ATMDIF
                     GO TO 500
                  END IF
                  IF (PL .GT. PMIN .AND. PL .LT. ZERO) THEN
cc                     PONDING = .FALSE.
                     IFATM(I) = 0
                     ATMACT(I) = ATMPOT(I)
C  need to decide how to treat this case: set PNEW to 0 or leave it
C  as is??
C+                   PNEW(I) = ZERO
                     PNEW(I) = ZERO
                     OVFLNOD(I) = -PONDNOD(I) * ARENOD(I) / DELTAT
                     GO TO 500
                  END IF
                  IF (PL .LE. PMIN) THEN
cc                     PONDING = .FALSE.
                     IFATM(I) = 0
                     ATMACT(I) = ATMPOT(I)
C  need to decide how to treat this case: set PNEW to 0 or leave it
C  as is??
C+                   PNEW(I) = ZERO
                     OVFLNOD(I) = -PONDNOD(I) * ARENOD(I) / DELTAT
                     GO TO 500
                  END IF
               ELSE
C
C  ponding, evaporation, exfiltration (==> return flow or seepage) 
C  (case D)
C
                  IF (PL .GE. PONDH_MIN) THEN
                     PONDING = .TRUE.
                     IFATM(I) = 2
                     PNEW(I) = PL
                     OVFLNOD(I) = ATMDIF
                     GO TO 500
                  END IF
                  IF (PL .GE. ZERO .AND. PL .LT. PONDH_MIN) THEN
cc                     PONDING = .FALSE.
                     IFATM(I) = 1
C  need to decide how to treat this case: set PNEW to 0, to PL, or
C  leave it as is??
                     PNEW(I) = ZERO
C+                   PNEW(I) = PL
                     OVFLNOD(I) = ATMDIF
                     GO TO 500
                  END IF
                  IF (PL .GT. PMIN .AND. PL .LT. ZERO) THEN
cc                     PONDING = .FALSE.
                     IFATM(I) = 0
                     ATMACT(I) = ATMPOT(I)
C  need to decide how to treat this case: set PNEW to 0 or leave it
C  as is??
C+                   PNEW(I) = ZERO
                     PNEW(I) = ZERO
                     OVFLNOD(I) = -PONDNOD(I) * ARENOD(I) / DELTAT
                     GO TO 500
                  END IF
                  IF (PL .LE. PMIN) THEN
cc                     PONDING = .FALSE.
                     IFATM(I) = 0
                     ATMACT(I) = ATMPOT(I)
C  need to decide how to treat this case: set PNEW to 0 or leave it
C  as is??
C+                   PNEW(I) = ZERO
                     OVFLNOD(I) = -PONDNOD(I) * ARENOD(I) / DELTAT
                     GO TO 500
                  END IF
               END IF
            END IF
         END IF
C
         IF (IFATM(I) .EQ. 1) THEN
            IF (ATMPOT(I) .GE. ZERO) THEN
               IF (ATMACT(I) .GE. ZERO) THEN
C
C  saturated (but not ponded) or 'air dry', rainfall, infiltration 
C  (this is correct even for 'air dry' case, though this scenario
C  should not occur because it is pre-empted by the call to ADRSTN)
C  (case E)
C
                  IF (PL .GE. PONDH_MIN) THEN
                     PONDING = .TRUE.
                     IFATM(I) = 2
                     PNEW(I) = PL
                     OVFLNOD(I) = ATMDIF
                     GO TO 500
                  END IF
                  IF (PL .GE. ZERO .AND. PL .LT. PONDH_MIN) THEN
cc                     PONDING = .FALSE.
                     IFATM(I) = 1
C  need to decide how to treat this case: set PNEW to 0, to PL, or
C  leave it as is??
C+                   PNEW(I) = ZERO
C+                   PNEW(I) = PL
                     OVFLNOD(I) = ATMDIF
                     GO TO 500
                  END IF
                  IF (PL .GT. PMIN .AND. PL .LT. ZERO) THEN
cc                     PONDING = .FALSE.
                     IFATM(I) = 0
                     ATMACT(I) = ATMPOT(I)
C  need to decide how to treat this case: set PNEW to 0 or leave it
C  as is??
C+                   PNEW(I) = ZERO
                     OVFLNOD(I) = -PONDNOD(I) * ARENOD(I) / DELTAT
                     GO TO 500
                  END IF
                  IF (PL .LE. PMIN) THEN
cc                     PONDING = .FALSE.
                     IFATM(I) = 0
                     ATMACT(I) = ATMPOT(I)
C  need to decide how to treat this case: set PNEW to 0 or leave it
C  as is??
C+                   PNEW(I) = ZERO
                     OVFLNOD(I) = -PONDNOD(I) * ARENOD(I) / DELTAT
                     GO TO 500
                  END IF
               ELSE
C
C  saturated (but not ponded) or 'air dry', rainfall, exfiltration
C  (==> return flow or seepage) 
C  (this is correct even for 'air dry' case, though this scenario
C  should not occur because it is pre-empted by the call to ADRSTN)
C  (case F)
C
                  IF (PL .GE. PONDH_MIN) THEN
                     PONDING = .TRUE.
                     IFATM(I) = 2
                     PNEW(I) = PL
                     OVFLNOD(I) = ATMDIF
                     GO TO 500
                  END IF
                  IF (PL .GE. ZERO .AND. PL .LT. PONDH_MIN) THEN
cc                     PONDING = .FALSE.
                     IFATM(I) = 1
C  need to decide how to treat this case: set PNEW to 0, to PL, or
C  leave it as is??
C+                   PNEW(I) = ZERO
C+                   PNEW(I) = PL
                     OVFLNOD(I) = ATMDIF
                     GO TO 500
                  END IF
                  IF (PL .LT. ZERO) THEN
cc                     PONDING = .FALSE.
                     IFATM(I) = 1
C  need to decide how to treat this case: set PNEW to 0 or leave it
C  as is??
C+                   PNEW(I) = ZERO
                     OVFLNOD(I) = ATMDIF
cxcx                 WRITE(IOUT9,9020) I,PL,TIME,DELTAT
                     GO TO 500
                  END IF
               END IF
            ELSE
               IF (ATMACT(I) .GE. ZERO) THEN
C
C  saturated (but not ponded) or 'air dry', evaporation, infiltration
C  (case G)
C
                  IF (PNEW(I) .LE. PMIN) THEN
                     if (atmact(i) .lt. atmpot(i)) then
                        ifatm(i)=0
                        atmact(i)=atmpot(i)
                        pnew(i)=pmin
                        ovflnod(i)=zero
                        go to 500
                     else 
                        go to 500
                     end if
cxcx                 WRITE(IOUT9,9040) I,ATMACT(I),TIME,DELTAT
                  END IF
                  IF (PL .GE. PONDH_MIN) THEN
                     PONDING = .TRUE.
                     IFATM(I) = 2
                     PNEW(I) = PL
                     OVFLNOD(I) = ATMDIF
                     GO TO 500
                  END IF
                  IF (PL .GE. ZERO .AND. PL .LT. PONDH_MIN) THEN
cc                     PONDING = .FALSE.
                     IFATM(I) = 1
C  need to decide how to treat this case: set PNEW to 0, to PL, or
C  leave it as is??
C+                   PNEW(I) = ZERO
C+                   PNEW(I) = PL
                     OVFLNOD(I) = ATMDIF
                     GO TO 500
                  END IF
                  IF (PL .GT. PMIN .AND. PL .LT. ZERO) THEN
cc                     PONDING = .FALSE.
                     IFATM(I) = 0
                     ATMACT(I) = ATMPOT(I)
C  need to decide how to treat this case: set PNEW to 0 or leave it
C  as is??
C+                   PNEW(I) = ZERO
                     OVFLNOD(I) = -PONDNOD(I) * ARENOD(I) / DELTAT
                     GO TO 500
                  END IF
                  IF (PL .LE. PMIN) THEN
cc                     PONDING = .FALSE.
                     IFATM(I) = 0
                     ATMACT(I) = ATMPOT(I)
C  need to decide how to treat this case: set PNEW to 0 or leave it
C  as is??
C+                   PNEW(I) = ZERO
                     OVFLNOD(I) = -PONDNOD(I) * ARENOD(I) / DELTAT
                     GO TO 500
                  END IF
               ELSE
C
C  saturated (but not ponded) or 'air dry', evaporation, exfiltration
C  (==> return flow or seepage for the case 'saturated but not ponded') 
C  (case H)
C
                  IF (PNEW(I) .LE. PMIN) THEN
                     IF (ATMACT(I) .LT. ATMPOT(I)) THEN
cc                        PONDING = .FALSE.
                        IFATM(I) = 0
                        ATMACT(I) = ATMPOT(I)
                        PNEW(I) = PMIN
                        OVFLNOD(I) = ZERO
                        GO TO 500
                     else
                        go to 500
                     END IF
                  END IF
                  IF (PL .GE. PONDH_MIN) THEN
                     PONDING = .TRUE.
                     IFATM(I) = 2
                     PNEW(I) = PL
                     OVFLNOD(I) = ATMDIF
                     GO TO 500
                  END IF
                  IF (PL .GE. ZERO .AND. PL .LT. PONDH_MIN) THEN
cc                     PONDING = .FALSE.
                     IFATM(I) = 1
C  need to decide how to treat this case: set PNEW to 0, to PL, or
C  leave it as is??
C+                   PNEW(I) = ZERO
C+                   PNEW(I) = PL
                     OVFLNOD(I) = ATMDIF
                     GO TO 500
                  END IF
                  IF (PL .GT. PMIN .AND. PL .LT. ZERO) THEN
cc                     PONDING = .FALSE.
                     IFATM(I) = 0
                     ATMACT(I) = ATMPOT(I)
C  need to decide how to treat this case: set PNEW to 0 or leave it
C  as is??
C+                   PNEW(I) = ZERO
                     OVFLNOD(I) = -PONDNOD(I) * ARENOD(I) / DELTAT
                     GO TO 500
                  END IF
                  IF (PL .LE. PMIN) THEN
cc                     PONDING = .FALSE.
                     IFATM(I) = 0
                     ATMACT(I) = ATMPOT(I)
C  need to decide how to treat this case: set PNEW to 0 or leave it
C  as is??
C+                   PNEW(I) = ZERO
                     OVFLNOD(I) = -PONDNOD(I) * ARENOD(I) / DELTAT
                     GO TO 500
                  END IF
               END IF
            END IF
         END IF
C
         IF (IFATM(I) .EQ. 0) THEN
            IF (ATMPOT(I) .GE. ZERO) THEN
C
C  unsaturated (Neumann BC) with rainfall
C  (case I)
C
               IF (PNEW(I) .GE. PONDH_MIN) THEN
                  PONDING = .TRUE.
                  IFATM(I) = 2
                  OVFLNOD(I) = (PNEW(I)-PONDNOD(I)) * ARENOD(I) / DELTAT
                  GO TO 500
               END IF
               IF (PNEW(I) .GE. ZERO .AND. PNEW(I) .LT. PONDH_MIN) THEN
cc                  PONDING = .FALSE.
                  IFATM(I) = 1
C  need to decide how to treat this case: set PNEW to 0 or leave it
C  as is??
C+                PNEW(I) = ZERO
                  OVFLNOD(I) = (PNEW(I)-PONDNOD(I)) * ARENOD(I) / DELTAT
                  GO TO 500
               END IF
               IF (PNEW(I) .GT. PMIN .AND. PNEW(I) .LT. ZERO) THEN
cc                  PONDING = .FALSE.
                  IFATM(I) = 0
                  OVFLNOD(I) = -PONDNOD(I) * ARENOD(I) / DELTAT
                  GO TO 500
               END IF
               IF (PNEW(I) .LE. PMIN) THEN
cc                  PONDING = .FALSE.
                  IFATM(I) = 1
                  PNEW(I) = PMIN
                  OVFLNOD(I) = -PONDNOD(I) * ARENOD(I) / DELTAT
cxcx              WRITE(IOUT9,9060) I,ATMPOT(I),TIME,DELTAT
                  GO TO 500
               END IF
            ELSE
C
C  unsaturated (Neumann BC) with evaporation
C  (case J)
C
               IF (PNEW(I) .GE. PONDH_MIN) THEN
                  PONDING = .TRUE.
                  IFATM(I) = 2
                  OVFLNOD(I) = (PNEW(I)-PONDNOD(I)) * ARENOD(I) / DELTAT
                  GO TO 500
               END IF
               IF (PNEW(I) .GE. ZERO .AND. PNEW(I) .LT. PONDH_MIN) THEN
cc                  PONDING = .FALSE.
                  IFATM(I) = 1
C  need to decide how to treat this case: set PNEW to 0 or leave it
C  as is??
C+                PNEW(I) = ZERO
                  OVFLNOD(I) = (PNEW(I)-PONDNOD(I)) * ARENOD(I) / DELTAT
                  GO TO 500
               END IF
               IF (PNEW(I) .GT. PMIN .AND. PNEW(I) .LT. ZERO) THEN
cc                  PONDING = .FALSE.
                  IFATM(I) = 0
                  OVFLNOD(I) = -PONDNOD(I) * ARENOD(I) / DELTAT
                  GO TO 500
               END IF
               IF (PNEW(I) .LE. PMIN) THEN
cc                  PONDING = .FALSE.
                  IFATM(I) = 1
                  PNEW(I) = PMIN
                  OVFLNOD(I) = -PONDNOD(I) * ARENOD(I) / DELTAT
                  GO TO 500
               END IF
            END IF
         END IF
  500    CONTINUE
cm       if (i.eq.264) then
cm          write(666,*) 'dopo ',time,atmpot(i),atmact(i),pl,pnew(i),
cm   1                    ifatm(i)
cm       end if
cd         write(99,'(i2,2x,7f10.5,i4)') i,atmpot(i),atmact(i),atmdif,
cd     1                              pondnod(i),pl,pnew(i),ovflnod(i),
cd     2                              ifatm(i)
      END DO
C
      RETURN
 9000 FORMAT(  'SWITCH anomaly: unsaturated node ',I6,
     1         ' (PL = ',1PE12.5,')',
     2       /,'  generated where IFATM was 2 and under ',
     3         'rainfall + exfiltration conditions ',
     4       /,'  (ie, return flow/seepage)',
     5       /,'  (TIME=',1PE8.2,', DELTAT=',1PE8.2,')')
 9020 FORMAT(  'SWITCH anomaly: unsaturated node ',I6,
     1         ' (PL = ',1PE12.5,')',
     2       /,'  generated where IFATM was 1 and under ',
     3         'rainfall + exfiltration conditions ',
     4       /,'  (ie, return flow/seepage)',
     5       /,'  (TIME=',1PE8.2,', DELTAT=',1PE8.2,')')
 9040 FORMAT(  'SWITCH anomaly: air dry node ',I6,' with ',
     1         'Dirichlet BC yields infiltration flux ',
     2       /,'  (ATMACT = ',1PE12.5,')',
     3       /,'  (TIME=',1PE8.2,', DELTAT=',1PE8.2,')')
 9060 FORMAT(  'SWITCH anomaly: Neumann node ',I6,' under ',
     1         'rainfall conditions becomes air dry ',
     2       /,'  (ATMPOT = ',1PE12.5,')',
     3       /,'  (TIME=',1PE8.2,', DELTAT=',1PE8.2,')')
      END
