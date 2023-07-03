C
C --------------------------  PLANT_GROWTH  --------------------------------
C 
C   Subroutine for the simulation of crop growth, modified from the model
C   SUCROS (SUCROS1_97 V1.0)
C
C --------------------------------------------------------------------------  
C
      SUBROUTINE PLANT_GROWTH(TIME,DELTAT,DAVTMP,FC,RLAI)


C 1.1 INTRODUCTION
C --------------------------------------------------------------------------  
C    The model for potential growth (example: spring wheat) 
C    is described in:
C      
C    J. Goudriaan & H.H. van Laar, 1994. Modelling Potential Crop     
C       Growth Processes. Textbook with Exercises. Kluwer Academic    
C       Publishers, Dordrecht, The Netherlands, 238 pp.               
C
C    The model is modified and coupled with CATHY Root Water Uptake
C    (August 2012 - Manoli & Bonetti)
C
C --------------------------------------------------------------------------  


      IMPLICIT NONE 

      INCLUDE 'CATHY.H'
      INCLUDE 'PLANT.H'
      INCLUDE 'IOUNITS.H'

      INTEGER I,J,K,ITYPE
      REAL*8 TIME,DELTAT,FC(PLMAX)
      REAL*8 DAVTMP
      REAL*8 GPHOT,MAINT,MAINTS,MNDVS,TEFF
      REAL*8 ASRQ,TRANSL
      REAL*8 GTW,GRT,GLV,GST,GSO
      REAL*8 DLV,RWLVG
      REAL*8 TAI,REAI,RLAI,GLAI,DLAI,RDR,RDRDV,RDRSH
      REAL*8 AA,BB,CC,DD,DTEFF
      REAL*8 SLATAB,PLANT_SLA

C --------------------------------------------------------------------------  

      DO I=1,NPLANT
         ITYPE = ITYPEP(I)  



C CARBOHYDRATE PRODUCTION
C --------------------------------------------------------------------------  

c conversion of CO2 to carbohydrates assuming that:
c CO2 + H2O --> CH2O + O2
c
c MW(CO2)  = 44 g/mol
c MW(CH2O) = 30 g/mol

c FC = [mumol/m2/s] 
c GPHOT daily total gross CH2O assimilation of the crop [gCH2O/m2/d]

      
          GPHOT = FC(I)*1e-6*44


C MAINTENANCE RESPIRATION
C --------------------------------------------------------------------------  

          MAINTS = MAINLV(ITYPE)*WLVG(I)+MAINST(ITYPE)*WST(I)+
     1             MAINRT(ITYPE)*WRT(I)+MAINSO(ITYPE)*WSO(I)  
          TEFF = Q10(ITYPE)**((DAVTMP-TREF(ITYPE))/10)

          IF(WLV(I).NE.0) THEN
                  MNDVS = WLVG(I)/WLV(I)
          ELSE
                  MNDVS = WLVG(I)
          ENDIF

          MAINT = MAINTS*TEFF*MNDVS



C GROWTH RATES OF PLANT ORGANS AND TRANSLOCATION
C --------------------------------------------------------------------------  

c ASRQR = [g CH2O / g DM]

           ASRQ = FSH*(ASRQLV(ITYPE)*FLV+ASRQST(ITYPE)*FST
     1            +ASRQSO(ITYPE)*FSO)+ASRQRT(ITYPE)*FRT

           IF(DVS-1.LT.0) THEN
                   TRANSL = 0
           ELSE
                   TRANSL = WST(I)*DVR*FRTRL(ITYPE)
           ENDIF

           GTW = (GPHOT-MAINT+CONVL(ITYPE)*TRANSL*CFST(ITYPE)*30./12.)
     1           /ASRQ

           GRT = FRT*GTW 
           GLV = FLV*FSH*GTW
           GST = FST*FSH*GTW-TRANSL
           GSO = FSO*FSH*GTW


C 1.8 LEAF AND EAR DEVELOPMENT
C --------------------------------------------------------------------------  

C           TAI = 0.5*EAI(I) + LAI2(I)

c--------------------------------
c GLAI
c---------------------------------
c computes daily increase of leaf 
c area index LAI2
c GLAI [m2 leaf / m2 ground / d]           
           
c maturation stage
           SLATAB = PLANT_SLA(DVS)
           IF(SLA(ITYPE).EQ.0) THEN
              GLAI = SLATAB*GLV
           ELSE
              GLAI = SLA(ITYPE)*GLV
           ENDIF
c juvenile stage
           IF((DVS.LT.0.3).AND.(LAI2(I).LT.0.75)) THEN
                   CC = 0.0
                   DD = DAVTMP - TBASE(ITYPE)
                   DTEFF = MAX(CC,DD)
                   GLAI =(LAI2(I)*(EXP(RGRL(ITYPE)*DTEFF*DELTAT)-1))
     1                    /DELTAT
           ENDIF
c---------------------------------

c RLAI & REAI

c RLAI
           IF(DVS-1.LT.0) THEN
                   RDRDV = 0.0
           ELSE
                   AA = 0.1
                   RDRDV = DVR/(MAX(AA,2-DVS))*FRDR(ITYPE)
           ENDIF

           BB = 0.03*(LAI2(I)-LAICR(ITYPE))/LAICR(ITYPE)

           IF((BB.GE.0).AND.(BB.LE.0.03)) THEN
                   RDRSH = BB
           ELSEIF(BB.LT.0) THEN
                   RDRSH = 0
           ELSEIF(BB.GT.0.03) THEN
                   RDRSH = 0.03
           ENDIF

c RDRDV = [s-1]
c RDRSH = [d-1] ---> ... / LDAY = [s-1] 
           RDR = MAX(RDRDV,RDRSH/LDAY)

c --------------------------------
c REAI
c---------------------------------
c calculates ear area index

           IF (DVS.LT.0.8) REAI = 0.
           IF ((DVS.GE.0.8).AND.(EAI(I).EQ.0.)) THEN
                   REAI = (EAR(ITYPE)*TADRW(I))/DELTAT
           ELSE
                   REAI = 0.
           ENDIF
           IF (DVS.GE.1.3) REAI = -RDRDV * EAI(I)
c----------------------------------
c LAI & EAI

           DLAI = LAI2(I) * RDR
           RLAI = GLAI - DLAI

           LAI2(I) = LAI2(I) + RLAI*DELTAT
           EAI(I) = EAI(I) + REAI*DELTAT 

           IF(LAI2(I).NE.0) THEN
                  DLV = WLVG(I)*DLAI/LAI2(I)
           ELSE
                  DLV = WLVG(I)*DLAI
           ENDIF
                 
           RWLVG = GLV - DLV

C 1.9 DRY MATTER PRODUCTION
C --------------------------------------------------------------------------  
           IF(DVS.LE.2.D0) THEN
           WRT(I)  = WRT(I) + GRT*DELTAT
           WLVG(I) = WLVG(I) + RWLVG*DELTAT
           WLVD(I) = WLVD(I) + DLV*DELTAT
           WST(I)  = WST(I) + GST*DELTAT
           WSO(I)  = WSO(I) + GSO*DELTAT

c Total above ground dry matter TADRM = [g DM / m2]
           WLV(I) = WLVG(I) + WLVD(I)
           TADRW(I) = WLV(I) + WST(I) + WSO(I)

c Total dry matter TDRW = [g DM / m2]           
           TDRW(I) = TADRW(I) + WRT(I)
           ENDIF


C **************************************************************************

      ENDDO


C **************************************************************************
c        Sistema stampa crescita piante ---
C         WRITE(IOUT72,*) TIME,TADRW(50),TADRW(51),LAI2(50),
C     1                   LAI2(51),TDRW(50),TDRW(51),WRT(50),
C     2                   WRT(51),WLVG(50),WLVG(51),WLVD(50),
C     3                   WLVD(51),WLV(50),WLV(51),WST(50),
C     4        i           WST(51),WSO(50),WSO(51)

         IF(NPLANT.LE.10) THEN
         WRITE(IOUT72,*) TIME,(WRT(I),WLVG(I),WLVD(I),WLV(I),WST(I),
     1                   WSO(I),TADRW(I),TDRW(I),LAI2(I),FC(I)*LAI2(I),
     2                   i=1,NPLANT)
         ELSE
c griglia senza fossi                 
c         WRITE(IOUT72,*) TIME,WSO(311),WSO(91),WSO(127),
c     1                   TADRW(311),TADRW(91),TADRW(127),
c     2                   LAI2(311),LAI2(91),LAI2(127)
c griglia con fossi
         WRITE(IOUT72,*) TIME,WSO(635),WSO(503),WSO(667),
     1                   TADRW(635),TADRW(503),TADRW(667),
     2                   LAI2(635),LAI2(503),LAI2(667)
         ENDIF        
c         WRITE(IOUT72,*) TIME,(WSO(I),TADRW(I),TDRW(I),i=1,NPLANT)

c Write WSO(i) in [tonn/ha]: ---------------------------------------------- 
c         WRITE(IOUT72,*) TIME,(WSO(I)/100,i=1,NPLANT)
c--------------------------------------------------------------------------

C         WRITE(1111,*) TIME,(LAI2(I),i=1,NPLANT)

      RETURN
 
      END 
