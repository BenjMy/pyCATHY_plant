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
     1                     ZABLOLD,ATMINP,TSNEW)

      IMPLICIT NONE

      INCLUDE 'CATHY.H'
      INTEGER I,ITER
      REAL*8  TRASP(PLMAX),METEO(*),DELTAT,TIME,ZABLOLD
      REAL*8  ZABL,ZLCL,BETA,LAPSE,LATHEAT,H,ATMPOT(*),ZLCLOLD
      REAL*8  ATMINP(3,*),FLAG,MODEL_FLAG,ST_FACTOR
      REAL*8  ALBEDO,A,B,EPS_S,SIGMA,SATVAP,ACTVAP,EPS_A,ZM,ZH,D
      REAL*8  ZOM,ZOH,GA,RA,TOLL,RNS,TSOLD,TSNEW,RHO_A,RHO_V,CP
      REAL*8  H1,H2,RNL,TS,RNET,ZOS,X,TS_OK,G,LE,ST,DRNET,FUN,DFUN
      REAL*8  ZLCL1,MIXRATIO,PLCL,TLCL,R,MA,GRAV
      REAL*8  LOBU,XSTAB,STABH,H_temp,dlobu,dstabh,dga,zita,dzita
      REAL*8  dh_temp,ga1,TSOLDOLD,rnet1,gaprova,ii,X1,X2,TT
      REAL*8  ST2,RNET2,H_temp2,LOBU2,ZITA2,XSTAB2,STABH2,FUN2,GA2
      REAL*8  H3,ST3,RNET3,H_temp3,LOBU3,ZITA3,XSTAB3,STABH3,FUN3,GA3
      REAL*8  REY
      LOGICAL PROVA

      INCLUDE 'IOUNITS.H'
      INCLUDE 'PLANT.H'

C --------------------------------------------------------------------------
     
      ST_FACTOR=0.15
c      ST_FACTOR=0.d0
      
      IF (ABL_MODEL.EQ.1) THEN
      
      SATVAP  = 0.6108*EXP(17.27*METEO(1)/(METEO(1)+237.3))
      ACTVAP  = SATVAP*METEO(2)/100*10
      EPS_A   = 1.24*((ACTVAP/(METEO(1)+273.15))**(1.d0/7.d0))
      D       = 2.d0/3.d0*HCANO(1)
      ZOM     = 0.123*HCANO(1)
      ZOH     = 0.1*ZOM
      RA      = LOG((ZM_ABL-D)/ZOM)*LOG((ZH-D)/ZOH)/0.41**2.d0/METEO(8)
      GA1      = 1.d0/RA
      RNS     = (1.d0-ALBEDO)*METEO(7)
      R  = 8314
      MA = 29
      GRAV = 9.8
      rey = meteo(11)*zom/1.4e-5
      if (rey.lt.20000) meteo(11)=1.4e-5*20000/zom

      IF(METEO(7).LT.0) RNS=0






C Trova punti iniziali per regula falsi 
        X1 = METEO(1)
        X2=METEO(1)
        PROVA=.FALSE.
        DO I=1,201
          II = I-101
          TT=METEO(1)+11
          RNET  = RNS-SIGMA*(EPS_S*(TT+273.15)**4.d0-EPS_A*(METEO(1)+
     1            273.15)**4.d0) 
          LE    = LATHEAT*TRASP(1)/ACANO(1)/18*1000*1000
          G     = METEO(9)
          ST    = ST_FACTOR*RNET
          H_temp= 1.d0/(CP*RHO_A)*(RNET-LE-G-ST)  
          LOBU  =-(METEO(11)**3.D0)*(METEO(1)+273.15)/(0.41*GRAV*H_temp)
          ZITA  = (ZH-D)/LOBU
          XSTAB = (1.d0-15.d0*ZITA)**(1.D0/4.D0)
          IF (ZITA.LT.0.d0) THEN
             STABH = 2.d0*LOG((1+XSTAB**2.d0)/2.d0)
          ELSEIF(ZITA.GT.1.D0) THEN
             STABH = -5.d0
          ELSE
             STABH = -5.d0*ZITA
          ENDIF
          GA   = METEO(11)*0.41/(LOG((ZH-D)/ZOH)-STABH)
          H     = GA*(TT-METEO(1))
          FUN   = H-1.d0/(CP*RHO_A)*(RNET-LE-G-ST)
          IF (FUN.LT.0.D0) X1=TT
          IF (FUN.GT.0.D0) X2=TT
      ENDDO 
          WRITE(7893,*) TIME, METEO(1), 'X1=',X1,'  X2=',X2
c        IF((X1.EQ.METEO(1)).OR.(X2.EQ.METEO(1))) THEN
c           TSNEW = METEO(1)
c           GO TO 200   
c        ENDIF    

        ITER    = 0
        TOLL    = 1E-2

        WRITE(7890,*) TIME,'TIME'
        WRITE(7890,*) 'ITER TS_NEW TS_OLD SCARTO FUN DFUN'
        DO WHILE((TOLL.GT.1E-6).AND.(ITER.LE.30))
          ITER  = ITER + 1          

          RNET  = RNS-SIGMA*(EPS_S*(X1+273.15)**4.d0-EPS_A*(METEO(1)+
     1            273.15)**4.d0)
          RNET2 = RNS-SIGMA*(EPS_S*(X2+273.15)**4.d0-EPS_A*(METEO(1)+
     1            273.15)**4.d0)

          LE    = LATHEAT*TRASP(1)/ACANO(1)/18*1000*1000
          G     = METEO(9)
          ST    = ST_FACTOR*RNET
          ST2    = ST_FACTOR*RNET2

          H_temp= 1.d0/(CP*RHO_A)*(RNET-LE-G-ST)
          H_temp2= 1.d0/(CP*RHO_A)*(RNET2-LE-G-ST2)


          LOBU  =-(METEO(11)**3.D0)*(METEO(1)+273.15)/(0.41*GRAV*H_temp)
          LOBU2=-(METEO(11)**3.D0)*(METEO(1)+273.15)/(0.41*GRAV*H_temp2)

          ZITA  = (ZH-D)/LOBU
          ZITA2  = (ZH-D)/LOBU2

          
          XSTAB = (1.d0-15.d0*ZITA)**(1.D0/4.D0)
          XSTAB2 = (1.d0-15.d0*ZITA2)**(1.D0/4.D0)

          IF (ZITA.LT.0.d0) THEN
             STABH = 2.d0*LOG((1+XSTAB**2.d0)/2.d0)
          ELSEIF(ZITA.GT.1.D0) THEN
             STABH = -5.d0
          ELSE
             STABH = -5.d0*ZITA
          ENDIF
          IF (ZITA2.LT.0.d0) THEN
             STABH2 = 2.d0*LOG((1+XSTAB2**2.d0)/2.d0)
          ELSEIF(ZITA2.GT.1.D0) THEN
             STABH2 = -5.d0
          ELSE
             STABH2 = -5.d0*ZITA2
          ENDIF

          GA   = METEO(11)*0.41/(LOG((ZH-D)/ZOH)-STABH)
          GA2   = METEO(11)*0.41/(LOG((ZH-D)/ZOH)-STABH2)
    
 
          H     = GA*(X1-METEO(1))
          H2     = GA*(X2-METEO(1))

          FUN   = H-1.d0/(CP*RHO_A)*(RNET-LE-G-ST)
          FUN2   = H2-1.d0/(CP*RHO_A)*(RNET2-LE-G-ST)

cccccccccccccccccccccccccccccccccccccccccccccccccccccc
          TSNEW = X1-FUN*(X2-X1)/(FUN2-FUN)
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
        
          RNET3 = RNS-SIGMA*(EPS_S*(TSNEW+273.15)**4.d0-EPS_A*(METEO(1)+
     1            273.15)**4.d0)

          ST3    = ST_FACTOR*RNET3

          H_temp3= 1.d0/(CP*RHO_A)*(RNET3-LE-G-ST3)


          LOBU3=-(METEO(11)**3.D0)*(METEO(1)+273.15)/(0.41*GRAV*H_temp3)
          ZITA3  = (ZH-D)/LOBU3

          
          XSTAB3 = (1.d0-15.d0*ZITA3)**(1.D0/4.D0)

          IF (ZITA3.LT.0.d0) THEN
             STABH3 = 2.d0*LOG((1+XSTAB3**2.d0)/2.d0)
          ELSEIF(ZITA3.GT.1.D0) THEN
             STABH3 = -5.d0
          ELSE
             STABH3 = -5.d0*ZITA3
          ENDIF

          GA3   = METEO(11)*0.41/(LOG((ZH-D)/ZOH)-STABH3)
 
          H3     = GA3*(tsnew-METEO(1))

          FUN3   = H3-1.d0/(CP*RHO_A)*(RNET3-LE-G-ST3)


          TOLL = ABS(FUN3)
          WRITE(7890,*) ITER,PROVA,TSNEW,METEO(1),TOLL,
     1                  ZITA,METEO(11),REY,GA
          IF(ITER.GT.30) WRITE(7890,*) 'Convergence NOT achieved'
c          IF((GA.LT.0).AND.(ITER.LT.30)) WRITE(7890,*) 'NEGATIVE GA'
          IF(FUN*FUN3.LT.0.D0) THEN
            X1=X1
            X2=TSNEW
          ELSE
            X1=TSNEW
            X2=X2
          ENDIF


         ENDDO


c ---------------------------- 
200      IF((ITER.GT.30)) THEN
         WRITE(7892,*) TIME,'TIME'
         DO I=1,201,1
          II=I-101
          RNET1  = RNS-SIGMA*(EPS_S*(II+273.15)**4.d0-EPS_A*(METEO(1)+
     1            273.15)**4.d0)
          h1 = 1.d0/cp/rho_a*(rnet1-le-g-0.15*rnet1)
          LOBU  =-(METEO(11)**3.D0)*(METEO(1)+273.15)/(0.41*GRAV*H1)

          ZITA  = (ZH-D)/LOBU
          
          XSTAB = (1.d0-15.d0*ZITA)**(1.D0/4.D0)

          IF (ZITA.LT.0.d0) THEN
             STABH = 2.d0*LOG((1+XSTAB**2.d0)/2.d0)
c          ELSEIF(ZITA.GT.1.D0) THEN
c             STABH = -5.d0
          ELSE
             STABH = -5.d0*ZITA

          ENDIF
          GAPROVA   = METEO(11)*0.41/(LOG((ZH-D)/ZOH)-STABH)
          H1=GAPROVA*(II-METEO(1))
          FUN   = H1-1.d0/(CP*RHO_A)*(RNET1-LE-G-0.15*rnet1)
          write(7892,*) ii,lobu,ZITA,xstab,GAPROVA,h1,fun       
          ENDDO
      ENDIF
c ---------------------------- 

    
      H1 = GA*(TSNEW-METEO(1))

      RNL = SIGMA*(EPS_S*(TSNEW+273.15)**4.d0-EPS_A*(METEO(1)+
     1       273.15)**4.d0)
     
      RNET = RNS - RNL

      H    = 1/(CP*RHO_A)*(RNET-LE-G-ST)

c      IF (ITER.GT.30) THEN
c         H=0.D0
c      ENDIF

C --------------------------------------
C Calculation of LCL according to Lawrence:

c      ZLCL = (20+METEO(1)/5.d0)*(100-METEO(2))

C Calculation of LCL according to Stull:

      MIXRATIO = 0.622*SATVAP*METEO(2)/100/
     1           (METEO(12)-SATVAP*METEO(2)/100)

      TLCL = 2840/(3.6*LOG(METEO(1)+273.15)-LOG(METEO(12)*MIXRATIO/
     1       (0.622+MIXRATIO))-7.108)+55

      PLCL = METEO(12)*(TLCL/(METEO(1)+273.15))**3.5


      ZLCL = R*(METEO(1)+263.15)/(MA*GRAV)*LOG10(METEO(12)/PLCL) 
C --------------------------------------

      IF(H.GT.0) ZABL = ZABL + H*(1+BETA)/(LAPSE*ZABL)*DELTAT


      IF(RNET.LT.0) ZABL=ZABL0

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
     1                  METEO(10),LE
      ELSE
        WRITE(IOUT73,*) TIME,RNS,RNL,METEO(6),RNET,METEO(1),TSNEW,H,
     1                  METEO(10),LE
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
