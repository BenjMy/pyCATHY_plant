C
C ----------------------------  ABLMODEL -----------------------------------
C   BETA     = 0.2-0.4 ..we assumed 0.2
C   LAPSE    = adiabatic lapse rate = 9.8 degK/km = 0.0098 degK/m 
C   LATHEAT  = latent heat of vaporization =2.25e+6 J/kg  
C   CP       = Specific Heat at constant pressure 
C            = 4183 J/kg/degK
C   RHO_V    = Vapor density = 0.7464 kg/m^3
C   RHO_A    = Air density = 1.2 kg/m^3
C   METEO(6) = Net Radiation [W/m^2]
C   METEO(7) = Shortwave incoming radiation [W/m^2]
C   METEO(8) = Wind speed [m/s]
C   RNS      = Net Shortwave Radiation [W/m^2]
C   RNL      = Net Longwave Radiation [W/m^2]
C   SIGMA    = Stefan-Boltzmann constant = 4.903E-9 MJ/degK^4/m^2/day 
C                                        = 5.67E-8 W/degK^4/m^2
C   ALBEDO   = 0.23
C   EPS_S    = Surface Emissivity [-] = 1 according to Brunt, 1952
C   EPS_A    = Atmospheric Emissivity [-]
C   A, B     = Parameters for the calculation of EPS_A
C              A = 0.62  and B = 0.05
C   SATVAP   = Saturation Vapor Pressure [kPa] 
C   ACTVAP   = Actual Vapor Pressure [mbar --> *10]
C   ZM  [m]  = Height of wind measurements
C   ZH  [m]  = Height of humidity measurements
C   D   [m]  = Zero plane displacement height 
C   ZOM [m]  = Roughness length governing momentum transfer
C   ZOH [m]  = Roughness length governing heat and vapor transfer
C   RA       = Aerodynamic Resistance [s/m]
C
C

C --------------------------------------------------------------------------
C
      SUBROUTINE ABLMODEL(TRASP,METEO,ZABL,DELTAT,TIME,ATMPOT,ZLCLOLD,
     1                     ZABLOLD,ATMINP,TSNEW,TAIR,WAIR)

      IMPLICIT NONE

      INCLUDE 'CATHY.H'
      INTEGER I,ITER
      REAL*8  TRASP(PLMAX),METEO(*),DELTAT,TIME,ZABLOLD
      REAL*8  ZABL,ZLCL,H,ATMPOT(*),ZLCLOLD
      REAL*8  ATMINP(3,*),FLAG,MODEL_FLAG,ST_FACTOR
      REAL*8  A,B,SATVAP,ACTVAP,EPS_A,ZM,D
      REAL*8  ZOM,ZOH,GA,RA,TOLL,RNS,TSOLD,TSNEW
      REAL*8  H1,H2,RNL,TS,RNET,ZOS,X,TS_OK,G,LE,ST,DRNET,FUN,DFUN
      REAL*8  ZLCL1,MIXRATIO,PLCL,TLCL,R,MA,GRAV,TAIR,DTAIR,TFA,TFA0
      REAL*8  DZABL,WFA,WFA0,DWAIR,WAIR,WSAT,TAREAL,RHREAL


      INCLUDE 'IOUNITS.H'
      INCLUDE 'PLANT.H'

C --------------------------------------------------------------------------

      IF (ABL_MODEL.GE.1) THEN
          TAREAL  = METEO(1)
          RHREAL  = METEO(2)
      ENDIF
      IF (ABL_MODEL.GT.1)THEN
          METEO(1)= TAIR-273.15
              SATVAP = 0.6108*EXP(17.27*METEO(1)/(METEO(1)+237.3))
              WSAT   = 0.622*SATVAP/(METEO(12)-SATVAP)
          METEO(2)= WAIR/WSAT*100
          write(1985,*) METEO(1),METEO(2),RHREAL,wair,satvap,wsat
      ENDIF

      ST_FACTOR=0.05
      
      IF (ABL_MODEL.GE.1) THEN
      
      SATVAP  = 0.6108*EXP(17.27*METEO(1)/(METEO(1)+237.3))
      ACTVAP  = SATVAP*METEO(2)/100*10
      EPS_A   = 1.24*((ACTVAP/(METEO(1)+273.15))**(1.d0/7.d0))
c      EPS_A   = (1.24*(1.d0+0.22*METEO(11)))*((ACTVAP/(METEO(1)+273.15))
c     1          **(1.d0/7.d0))
          
      D       = 2.d0/3.d0*HCANO(1)
C      D       = HCANO(1)*(1.D0-((0.32*METEO(8)/METEO(11))**(1.12))*
C     1          EXP(-0.197*METEO(8)/METEO(11)))
      ZOM     = 0.123*HCANO(1)
      ZOH     = 0.1*ZOM
      RA      = LOG((ZM_ABL-D)/ZOM)*LOG((ZH-D)/ZOH)/0.41**2.d0/METEO(8)
      GA      = 1.d0/RA
      RNS     = (1.d0-ALBEDO)*METEO(7)
      IF(METEO(7).LT.0) RNS=0
      write(1985,*)'ra...',d,hcano(1),(zm_abl-d)/zom,(zh-d)/zoh

        TSOLD   = METEO(1)
        ITER    = 0
        TOLL    = 1E-2

        WRITE(7890,*) TIME,'TIME',tsold,'TSOLD'

C       Ciclo NR per trovare TS:

        DO WHILE((TOLL.GT.1E-6).AND.(ITER.LE.30))
          ITER  = ITER + 1          

          RNET  = RNS-SIGMA*(EPS_S*(TSOLD+273.15)**4.d0-EPS_A*(METEO(1)+
     1            273.15)**4.d0)

          DRNET = -SIGMA*EPS_S*4.d0*(TSOLD+273.15)**3.d0

          H     = GA*(TSOLD-METEO(1))

c          LE   = LATHEAT*TRASP(1)/ACANO(1)*RHO_V
          LE    = LATHEAT*TRASP(1)/ACANO(1)/18*1000*1000
  
          G     = METEO(9)
      
          ST    = ST_FACTOR*RNET

          FUN   = H-1.d0/(CP*RHO_A)*(RNET-LE-G-ST)

          DFUN  = GA-1.d0/(CP*RHO_A)*(DRNET-ST_FACTOR*DRNET)

          TSNEW = TSOLD-FUN/DFUN


          TOLL = ABS(TSNEW-TSOLD)
          WRITE(7890,*) ITER,TSNEW,TSOLD,TOLL,fun,dfun,rnet,drnet,ga
          IF(ITER.GT.30) WRITE(7890,*) 'Convergence NOT achieved'
          TSOLD = TSNEW
         ENDDO

c          IF (ABS(TSNEW-METEO(1)).GT.5) THEN
c             TSNEW=METEO(1)
c             WRITE(7890,*) 'Impone TSNEW=TAIR'
c          ENDIF

C      if(meteo(7).le.1e-1) then
C        tsnew=meteo(1)
C      endif

C     Trovata TS calcolo i flussi
      H2 = GA*(TSNEW-METEO(1))

      RNL = SIGMA*(EPS_S*(TSNEW+273.15)**4.d0-EPS_A*(METEO(1)+
     1       273.15)**4.d0)
     
      RNET = RNS - RNL

      H    = 1/(CP*RHO_A)*(RNET-LE-G-ST)
C      if(meteo(7).le.1e-1) then
C        h=0.d0
C      endif

C --------------------------------------
C Calculation of LCL according to Lawrence:

c      ZLCL = (20+METEO(1)/5.d0)*(100-METEO(2))

C Calculation of LCL according to Stull:

      MIXRATIO = 0.622*SATVAP*METEO(2)/100/
     1           (METEO(12)-SATVAP*METEO(2)/100)

      TLCL = 2840/(3.5*LOG(METEO(1)+273.15)-LOG(METEO(12)*MIXRATIO/
     1       (0.622+MIXRATIO))-7.108)+55

      PLCL = METEO(12)*(TLCL/(METEO(1)+273.15))**3.5

      R  = 8314
      MA = 29
      GRAV = 9.8

      ZLCL = R*(METEO(1)+273.15)/(MA*GRAV)*LOG(METEO(12)/PLCL) 
C --------------------------------------

C     Calcola ZABL:
      DZABL=H*(1+BETA)/(LAPSE*ZABL)
c      IF(H.GT.0) ZABL = ZABL + H*(1+BETA)/(LAPSE*ZABL)*DELTAT
      IF(H.GT.0) ZABL=ZABL+DZABL*DELTAT
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Calcolo Tair (in degK..nella stampa e' in degC)
      
      TFA0  = 288.d0
      TFA   = TFA0+(LAPSE*ZABL)
      DTAIR = H+(TFA-TAIR)*DZABL
      IF(H.GT.0) TAIR=TAIR+DTAIR*DELTAT/ZABL

C     Calcola WAIR [kg/kg]
C     LAPSEW = [kg/kg/m]
      WFA0  = 0.00785 
c      LAPSEW= -0.00145*10e-3
C      LAPSEW= -0.0005*10e-3
      WFA   = WFA0+(LAPSEW*ZABL)
      DWAIR = LE/(RHO_A*LATHEAT)*18/1000+(WFA-WAIR)*DZABL
      IF(H.GT.0)      WAIR=WAIR+DWAIR*DELTAT/ZABL
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      WSAT = 0.622*SATVAP/(METEO(12)-SATVAP)

      IF(RNET.LT.0) THEN
         ZABL = ZABL0
         TAIR = TAREAL+273.15
         WAIR = RHREAL/100*WSAT
      WRITE(1985,*)'WAIR=',WAIR,RHREAL,TAIR,ZABL
      ENDIF



      IF (ZABL.GE.ZLCL) THEN
          ZABL = ZLCL   
          IF(ZLCLOLD.NE.ZABLOLD) THEN
              FLAG = 1
          ELSE
             FLAG = 0
          ENDIF
      ELSE
          FLAG = 0
      ENDIF
      ENDIF




C     Modello base..-----------------------------------------------

      IF (ABL_MODEL.EQ.0) THEN

C Calculate Lifting Condensation level  ZLCL
        ZLCL = (20+METEO(1)/5)*(100-METEO(2))

C Calculate Atmospheric boundary Layer ZABL
        G    = METEO(9)
        ST   = ST_FACTOR*METEO(6)
c        LE   = TRASP(1)/ACANO(ITYPEP(1))*RHO_V*LATHEAT
        LE    = LATHEAT*TRASP(1)/ACANO(1)/18*1000*1000
        H    =(METEO(6)-LE-G-ST)/RHO_A/CP 

        IF(H.GT.0)  ZABL = ZABL+H*(1+BETA)/(LAPSE*ZABL)*DELTAT      

        IF (METEO(6).LT.0) ZABL = ZABL0

        IF (ZABL.GE.ZLCL) THEN
          ZABL = ZLCL   
          IF(ZLCLOLD.NE.ZABLOLD) THEN
              FLAG = 1
          ELSE
             FLAG = 0
          ENDIF
        ELSE
          FLAG = 0
        ENDIF
      
      ENDIF


      IF (ATMPOT(5).GT.0) THEN
        WRITE(IOUT73,*) TIME,RNS,RNL,METEO(6),RNET,METEO(1),TSNEW,H,
     1                  METEO(10),LE,H2,TAIR-273.15,WAIR/WSAT*100,
     2                  METEO(2),WSAT,METEO(12),SATVAP,TAREAL,RHREAL
      ELSE
        WRITE(IOUT73,*) TIME,RNS,RNL,METEO(6),RNET,METEO(1),TSNEW,H,
     1                  METEO(10),LE,H2,TAIR-273.15,WAIR/WSAT*100,
     2                  METEO(2),WSAT,METEO(12),SATVAP,TAREAL,RHREAL
      ENDIF


      IF (ATMPOT(5).GT.0) THEN
        WRITE(IOUT74,*) TIME,FLAG,ATMPOT(5),ZABL,ZLCL,METEO(2),
     1                  METEO(9)
      ELSE
        WRITE(IOUT74,*) TIME,FLAG,0,ZABL,ZLCL,METEO(2),METEO(9)
      ENDIF
     
     

c      WRITE(IOUT73,*) TIME,DELTAT,ATMPOT(5),TRASP(1),H,ZABL,ZLCL
c      WRITE(IOUT73,*) TIME,ATMPOT(5),ZABL

      ZLCLOLD =ZLCL
      ZABLOLD = ZABL

c      IF (ZABL.GT.ZLCL) WRITE(IOUT74,*) TIME,ZABL
C --------------------------------------------------------------------------

      RETURN 
      END
