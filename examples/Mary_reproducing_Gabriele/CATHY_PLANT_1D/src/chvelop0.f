C
C**************************  CHVELOP0 **********************************
C
C  calculate soil moisture characteristics needed for storage and
C  velocity calculations and for output (peat soil case)
C
C***********************************************************************
C
      SUBROUTINE CHVELOP0(NLKP,NNOD,N,NT,NTRI,IVGHU,TETRA,TP,
     1                    PTIMEP,Z,INDE,INDE0,PORE,
     2                    SNODI,PNODI,SW,CKRW,SENODI,CKRWE,SEELT)
C
      IMPLICIT NONE
      INCLUDE 'CATHY.H'
      INTEGER  I,NS,ITERV,IVMAX
      INTEGER  NLKP,NNOD,N,NT,NTRI,IVGHU
      INTEGER  TETRA(5,*),TP(*)
      REAL*8   FVGKR,FVGSE,FBCKR,FBCSE,FCINDE,CORNR
      REAL*8   PSI,SE,D,SCARTO,SCAROLD,TOL,CBETA
      REAL*8   PTIMEP(*),Z(*),INDE(*),INDE0(*),PORE(*)
      REAL*8   SNODI(*),PNODI(*),SW(*),CKRW(*)
      REAL*8   SENODI(*),CKRWE(*),SEELT(*)
      INCLUDE 'SOILCHAR.H'
C
      IF (IVGHU .EQ. -1) THEN
         CALL MOISTABVELPT(NLKP,NT,NTRI,TETRA,PTIMEP,SEELT,CKRWE)
         CALL ELTNOD(N,NT,TP,TETRA,CKRWE,CKRW)
         CALL ELTNOD(N,NT,TP,TETRA,SEELT,SENODI)
         DO I=1,N
            NS       = I-(NNOD*INT((I-1)/NNOD))
            D        = Z(NS) - Z(I)
            CBETA    = CBETA0+CANG*D
            IF (CBETA .GT. 1.0D0) CBETA=1.0D0
            PSI      = PTIMEP(I)
            SE       = SENODI(I)
            IVMAX=100
            TOL = 1.0D-3
            ITERV = 0
            SCARTO  =  2.0d0*TOL
C***********************************************************************
C Ciclo per l'inizializzazione dell'indice dei vuoti (porosita')
C***********************************************************************
            DO WHILE (SCARTO.GT.TOL .AND. ITERV.LT.IVMAX)
               ITERV = ITERV +1 
               IF (ITERV .EQ. IVMAX) WRITE(500,77) I
               VGPNOT(I)= (PORE(I)-VGRMC)/PORE(I)
               SW(I)    = VGPNOT(I)*SE + VGRMC/PORE(I)
cm             IF (SW(I) .EQ. 1.0D0) THEN
               INDE0(I) = INDE(I)-CORNR(INDE(I),SW(I),D,PSI,PNODI(I),
     1                                  SNODI(I))
cm             ELSE
cm             INDE0(I)= FCINDE(INDE(I),SW(I),D,PSI,PNODI(I),SNODI(I),
cm   1                          PORE(I))
cm             END IF
               SCAROLD = SCARTO
               SCARTO  = ABS(INDE0(I)-INDE(I))
C            WRITE(500,'(4e12.4)')    INDE0(I),INDE(I),SW(I),PSI
C            WRITE(500,'(I5,3e12.4)') ITERV,INDE0(I),INDE(I),SCARTO/SCAROLD
               INDE(I) = INDE0(I)
               PORE(I)= INDE(I)/(1.0D0+INDE(I))
            END DO
	 END DO 
      ELSE IF (IVGHU .EQ. 0) THEN
         DO I=1,N
            NS       = I-(NNOD*INT((I-1)/NNOD))
            D        = Z(NS) - Z(I)
            CBETA    = CBETA0+CANG*D
            IF (CBETA .GT. 1.0D0) CBETA=1.0D0
            PSI      = PTIMEP(I)
            SE       = FVGSE(PSI)
            IVMAX=100
            TOL = 1.0D-3
            ITERV = 0
            SCARTO  =  2.0d0*TOL
C***********************************************************************
C Ciclo per l'inizializzazione dell'indice dei vuoti (porosita')
C***********************************************************************
            DO WHILE (SCARTO.GT.TOL .AND. ITERV.LT.IVMAX)
               ITERV = ITERV +1 
               IF (ITERV .EQ. IVMAX) WRITE(500,77) I
               VGPNOT(I)= (PORE(I)-VGRMC)/PORE(I)
               SW(I)    = VGPNOT(I)*SE + VGRMC/PORE(I)
               CKRW(I)  = FVGKR(PSI,SE)
cm             IF (SW(I) .EQ. 1.0D0) THEN
               INDE0(I) = INDE(I)-CORNR(INDE(I),SW(I),D,PSI,PNODI(I),
     1                                  SNODI(I))
cm             ELSE
cm             INDE0(I)= FCINDE(INDE(I),SW(I),D,PSI,PNODI(I),SNODI(I),
cm   1                          PORE(I))
cm             END IF
               SCAROLD = SCARTO
               SCARTO  = ABS(INDE0(I)-INDE(I))
C            WRITE(500,'(4e12.4)')INDE0(I),INDE(I),SW(I),PSI
C            WRITE(500,'(I5,3e12.4)') ITERV,INDE0(I),INDE(I),SCARTO/SCAROLD
               INDE(I) = INDE0(I)
               PORE(I)= INDE(I)/(1.0D0+INDE(I))
            END DO 
         END DO
      ELSE IF (IVGHU .EQ. 4) THEN
         DO I=1,N
            NS       = I-(NNOD*INT((I-1)/NNOD))
            D        = Z(NS) - Z(I)
            CBETA    = CBETA0+CANG*D
            IF (CBETA .GT. 1.0D0) CBETA=1.0D0
            PSI      = PTIMEP(I)
            SE       = FBCSE(PSI)
            IVMAX=100
            TOL = 1.0D-3
            ITERV = 0
            SCARTO  =  2.0d0*TOL
C***********************************************************************
C Ciclo per l'inizializzazione dell'indice dei vuoti (porosita')
C***********************************************************************
            DO WHILE (SCARTO.GT.TOL .AND. ITERV.LT.IVMAX)
               ITERV = ITERV +1 
               IF (ITERV .EQ. IVMAX) WRITE(500,77) I
               SW(I)    = BCPORM(I)*SE + BCRMC/PORE(I)
               CKRW(I)  = FBCKR(PSI)
cm             IF (SW(I) .EQ. 1.0D0) THEN
               INDE0(I) = INDE(I)-CORNR(INDE(I),SW(I),D,PSI,PNODI(I),
     1                                  SNODI(I))
cm             ELSE
cm             INDE0(I)= FCINDE(INDE(I),SW(I),D,PSI,PNODI(I),SNODI(I),
cm   1                          PORE(I))
cm             END IF
               SCAROLD = SCARTO
               SCARTO  = ABS(INDE0(I)-INDE(I))
C            WRITE(500,'(4e12.4)')INDE0(I),INDE(I),SW(I),PSI
C            WRITE(500,'(I5,3e12.4)') ITERV,INDE0(I),INDE(I),SCARTO/SCAROLD
               INDE(I) = INDE0(I)
               PORE(I)= INDE(I)/(1.0D0+INDE(I))
            END DO 
cm          INDE0(I) = FCINDE(INDE(I),SW(I),D,PSI,PNODI(I),SNODI(I),
cm   1                        PORE(I))
cm          SCARTO   = ABS(INDE0(I)-INDE(I))
cm          IF (SCARTO .GT. 1.0D-09) THEN
cm             INDE(I) = INDE0(I)
cm             GO TO 200
cm          END IF
         END DO 
      END IF
C
      RETURN
  77  FORMAT('La convergenza  non e''stata raggiunta nel nodo:',i6)
      END
