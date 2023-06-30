      SUBROUTINE TIM(TIME,ICOD)
      INTEGER ICOD
      REAL TIME
c
      IF (ICOD.EQ.1) THEN
c
C  SETTING INIZIALE
c
         ITIME = mclock()
         TIME = FLOAT(ITIME)/100.0
      ELSE
c
C  RILEVAZIONE PERIODICA
c
         ITIME = mclock()
         TIME = FLOAT(ITIME)/100.0 - TIME
      ENDIF
c
      RETURN
      END
