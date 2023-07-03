C
C**************************  DATIN  ************************************
C
C  read in some of the input data (other data is read in subroutines
C  BCONE, BCNXT, ATMONE, ATMNXT, NUDONE, NUDNXT, EFFONE, EFFNXT, RAST_INPUT,
C  FORCE_FLOWDIR, and FORCE_HG)
C
C***********************************************************************
C
      SUBROUTINE DATIN(ISIMGR,IVERT,ISP,WTPOSITION,BASE,ZRATIO, 
     1                 TRIANG,X,Y,Z,PERMX,PERMY,PERMZ,ELSTOR,
     2                 POROS,CONTR,NODVP,ID_QOUT,
     3                 TIMPRT,PONDNOD,PTIMEP,TETAF,
     4                 DELTAT,DTMIN,DTMAX,TMAX,DTMAGA,DTMAGM,
     5                 DTREDS,DTREDM,DTOUT,ITUNS,ITUNS1,ITUNS2,
     6                 TOLUNS,TOLSWI,ERNLMX,ITMXCG,TOLCG,
     7                 KSLOPE,LUMP,IPEAT,IVGHU,NLKP,
     8                 IOPT,ISOLV,IPRT1,IPRT,IPOND,INDP,NNOD,NTRI,NSTR,
     9                 NZONE,N1,NR,NUMVP,NUM_QOUT,NPRT,N,NT,
     A                 ISFONE,ISFCVG,DUPUIT,
     B                 L2NORM,NLRELX,OMEGA,
     C                 PONDH_MIN,
     D                 FL3D,SURF,DEM,GRID,NROW,NCOL,
     E                 NCELNL,NCELL,DOSTEP,NCELL_COARSE,
     F                 DEM_MAP,ZONE,LAKES_MAP,INDEX,INDEX_WITH_LAKES,
     G                 CELL,CELLCOL,CELLROW,TP2D,NODI,
     H                 TIPO_R,RESERVR,N_HA,CELLCOL_WL,CELLROW_WL,
     I                 BASE_MAP,DEPTH,ELTRIA,
     J                 DAFLAG,ERT_FLAG,NENS,NOBS,DSRETC,
     K                 DSATM,DSKS,DSSTOR,DSPOROS,DSSURF,DSIC,DSMEAS,
     L                 ENKFT,ENKFTIM,ENKFNOD,ENKFVAL,ENKF,NENSMIN,
     M                 NEFFMIN,DSMEASMAX,DSKSZ,PLANT_FLAG,NMETEODATA,
     N                 PLANT_PRINT,GROWTH_FLAG,SALT_FLAG,ABL_FLAG)
C
      IMPLICIT  NONE
      INCLUDE  'CATHY.H'
      INTEGER   IROW,ICOL,INDEX_C,ILAKE
      INTEGER   I,J,K,FLAGASC,PP
      INTEGER   ISIMGR,IVERT,ISP
      INTEGER   ITUNS,ITUNS1,ITUNS2,ITMXCG,KSLOPE,LUMP,IPEAT
      INTEGER   IVGHU,IOPT,ISOLV,IPRT1,IPRT,IPOND,INDP,NZONE
      INTEGER   NNOD,NTRI,NSTR
      INTEGER   NLKP
      INTEGER   N1,NR,NUMVP,NPRT,N,NT,ISFONE,ISFCVG,DUPUIT
      INTEGER   NUM_QOUT,L2NORM,NLRELX
      INTEGER   NROW,NCOL,NCELNL,NCELL,NUM_TOT_R
      INTEGER   DOSTEP,NCELL_COARSE
      INTEGER   DAFLAG,ERT_FLAG,NENS,NOBS,ENKFT,NENSMIN,NEFFMIN
      INTEGER   TRIANG(4,*)
      INTEGER   CONTR(*),NODVP(*),ID_QOUT(*)
      INTEGER   CELLCOL(*),CELLROW(*),TP2D(*)
      INTEGER   CELLCOL_WL(*),CELLROW_WL(*)
      INTEGER   LAKES_MAP(ROWMAX,COLMAX)
      INTEGER   ZONE(ROWMAX,COLMAX)
      INTEGER   INDEX(ROWMAX,COLMAX),INDEX_WITH_LAKES(ROWMAX,COLMAX)
      INTEGER   NODI(ROWMAX+1,COLMAX+1),CELL(5,*)
      INTEGER   RESERVR(*)
      INTEGER   TIPO_R(MAXCEL)
      INTEGER   N_HA(*)
      INTEGER   ENKFNOD(MAXNUDN,2)
      LOGICAL   GRID,DEM,FL3D,SURF,ENKF,PLANT_FLAG
      LOGICAL   GROWTH_FLAG,SALT_FLAG,ABL_FLAG
      REAL*8    SUMZ
      REAL*8    WTPOSITION,BASE,TETAF,DELTAT,DTMIN,DTMAX,TMAX
      REAL*8    DTMAGA,DTMAGM,DTREDS,DTREDM,DTOUT
      REAL*8    TOLUNS,TOLSWI,ERNLMX,TOLCG
      REAL*8    OMEGA
      REAL*8    PONDH_MIN,DSMEASMAX,DSKSZ
      REAL*8    DSRETC,DSKS,DSSTOR,DSPOROS,DSSURF,DSIC,DSATM
      REAL*8    ZRATIO(*),X(*),Y(*),Z(*),DEPTH(NODMAX)
      REAL*8    ELTRIA(NTRMAX),DSMEAS(MAXNUDT) 
      REAL*8    PERMX(MAXSTR,*),PERMY(MAXSTR,*),PERMZ(MAXSTR,*)
      REAL*8    ELSTOR(MAXSTR,*),POROS(MAXSTR,*)
      REAL*8    TIMPRT(*),PONDNOD(*),PTIMEP(*)
      REAL*8    DEM_MAP(ROWMAX,COLMAX),BASE_MAP(ROWMAX,COLMAX)
      REAL*8    ENKFTIM(MAXNUDT),ENKFVAL(MAXNUDT,MAXNUDN)

      INTEGER   NMETEODATA,PLANT_PRINT    
      INCLUDE  'IOUNITS.H'
      INCLUDE  'SOILCHAR.H'
      INCLUDE  'SURFWATER.H'
      INCLUDE  'RIVERNETWORK.H'
      INCLUDE  'RANDOM.H'
      INCLUDE  'PLANT.H'


C
C  unit IIN1 input 
C
      READ(IIN1,*) IPRT1,DAFLAG
      READ(IIN1,*) ISIMGR,PONDH_MIN
      IF (ISIMGR .EQ. 0) THEN
         FL3D=.TRUE.
         SURF=.FALSE.
         DEM=.FALSE.
         GRID=.TRUE.
      ELSE IF (ISIMGR .EQ. 1) THEN
         FL3D=.TRUE.
         SURF=.FALSE.
         DEM=.TRUE.
         GRID=.FALSE.
      ELSE IF (ISIMGR .EQ. 2) THEN
         FL3D=.TRUE.
         SURF=.TRUE.
         DEM=.TRUE.
         GRID=.FALSE.
      ELSE IF (ISIMGR .EQ. 3) THEN
         FL3D=.FALSE.
         SURF=.TRUE.
         DEM=.TRUE.
         GRID=.FALSE.   
      ELSE
         WRITE(IOUT2,1900) ISIMGR
         CALL CLOSIO
         STOP
      END IF
      READ(IIN1,*) KSLOPE,TOLKSL
      READ(IIN1,*) PKRL,PKRR,PSEL,PSER
      READ(IIN1,*) PDSE1L,PDSE1R,PDSE2L,PDSE2R
      READ(IIN1,*) ISFONE,ISFCVG,DUPUIT
      READ(IIN1,*) TETAF,LUMP,IOPT
      READ(IIN1,*) NLRELX,OMEGA
      READ(IIN1,*) L2NORM,TOLUNS,TOLSWI,ERNLMX
      READ(IIN1,*) ITUNS,ITUNS1,ITUNS2
      READ(IIN1,*) ISOLV,ITMXCG,TOLCG
      READ(IIN1,*) DELTAT,DTMIN,DTMAX,TMAX
      READ(IIN1,*) DTMAGA,DTMAGM,DTREDS,DTREDM
      READ(IIN1,*) IPRT,NPRT,(TIMPRT(I),I=1,NPRT)
      READ(IIN1,*) NUMVP,(NODVP(I),I=1,NUMVP)
      READ(IIN1,*) NR
      IF (NR .NE. 0) READ(IIN1,*) (CONTR(I),I=1,NR)
      READ(IIN1,*) NUM_QOUT,(ID_QOUT(I),I=1,NUM_QOUT)
      if (isimgr.le.1) pondh_min=1.0d+10
      WRITE(IOUT2,1000) ISIMGR,PONDH_MIN,DAFLAG
      WRITE(IOUT2,1010) IPRT1,IPRT,NPRT,NUMVP,NR
      WRITE(IOUT2,1015) KSLOPE,TOLKSL
      IF (KSLOPE .EQ. 3  .OR.  KSLOPE .EQ. 4) WRITE(IOUT2,1040)
     1                                      PKRL,  PKRR,  PSEL,  PSER,
     2                                      PDSE1L,PDSE1R,PDSE2L,PDSE2R
      WRITE(IOUT2,1045) ISFONE,ISFCVG
      WRITE(IOUT2,1004) TETAF,LUMP,IOPT
      WRITE(IOUT2,1050) NLRELX
      IF (NLRELX .EQ. 1) WRITE(IOUT2,1055) OMEGA
      WRITE(IOUT2,1060) L2NORM,TOLUNS,TOLSWI,ERNLMX
      WRITE(IOUT2,1025) ITUNS,ITUNS1,ITUNS2
      WRITE(IOUT2,1005) ISOLV,ITMXCG,TOLCG
      WRITE(IOUT2,1030) DELTAT,DTMIN,DTMAX,TMAX
      WRITE(IOUT2,1035) DTMAGA,DTMAGM,DTREDS,DTREDM
C
C  unit IIN2 input; if (grid) we just read the grid file; if (dem) we
C                   build up the surface mesh starting from the DEM
C                   of the basin.
C  
      IF (GRID) THEN 
         READ(IIN2,*) NZONE,NSTR,N1
         READ(IIN2,*) NNOD,NTRI
         WRITE(IOUT2,1020) NNOD,NTRI,NZONE,NSTR,N1
         READ(IIN2,*) IVERT,ISP,BASE
         WRITE(IOUT2,1350) IVERT,ISP,BASE
         READ(IIN2,*) (ZRATIO(I),I=1,NSTR)
         SUMZ=0.0D0
         DO I=1,NSTR
            SUMZ=SUMZ + ZRATIO(I)
         END DO
         WRITE(IOUT2,1340) (I,ZRATIO(I),I=1,NSTR)
         IF (DABS(SUMZ - 1.0D0) .GT. 1.0D-14) THEN
            WRITE(IOUT2,1345) 
            WRITE(IOUT2,*) SUMZ
            CALL CLOSIO
            STOP
         END IF
         IF (ISP .NE. 0) THEN
            READ(IIN2,*) (Z(I),I=1,NNOD)
         ELSE
            READ(IIN2,*) Z(1)
            DO I=2,NNOD
               Z(I)=Z(1)
            END DO
         END IF
         IF (IPRT1 .GE .1) WRITE(IOUT2,1360) (I,Z(I),I=1,NNOD)
         READ(IIN2,*) ((TRIANG(I,K),I=1,4),K=1,NTRI)
         READ(IIN2,*) (X(K),Y(K),K=1,NNOD)
        
         
      ELSE IF (DEM) THEN
        
         CALL RAST_INPUT_DEM(IIN10,NROW,NCOL,NORTH,SOUTH,EAST,WEST,
     1        DEM_MAP)
    
         CALL RAST_INPUT_LZ(IIN21,NROW,NCOL,NORTH,SOUTH,EAST,WEST,
     1        ZONE)
        
         READ(IIN11,*) DELTA_X,DELTA_Y
         WRITE(IOUT40,*) 'DELTA_X=',DELTA_X,'DELTA_Y=',DELTA_Y
         READ(IIN11,*) FACTOR
         WRITE(IOUT40,*)'FACTOR=',FACTOR
         READ(IIN11,*) DOSTEP
         WRITE(IOUT40,*) 'DOSTEP=',DOSTEP                  
        
         IF (FL3D) THEN
            READ(IIN11,*) NZONE,NSTR,N1 
            READ(IIN11,*) IVERT,ISP,BASE
            IF (IVERT.EQ.3) THEN
               CALL RAST_INPUT_DEM(IIN60,NROW,NCOL,NORTH,SOUTH,EAST,WEST
     1              ,BASE_MAP)
               CALL TRIANGOLI(NROW,NCOL,DELTA_X,DELTA_Y,WEST,SOUTH,
     1              BASE_MAP,ZONE,FACTOR,NNOD,NTRI,DOSTEP,
     2              NCELL_COARSE,NODI,TRIANG,TP2D,X,Y,DEPTH,ELTRIA,
     3              CELL)
            END IF
         END IF
         
c
c  reading of lakes_map
c              
         CALL RAST_INPUT_LZ(IIN20,NROW,NCOL,NORTH,SOUTH,EAST,
     1        WEST,LAKES_MAP)
c
c  cells numbering with and without lakes
c
        CALL INDEX_DEM(ROWMAX,MAXCEL,NROW,NCOL,NCELNL,NCELL,
     1       ZONE,INDEX,CELLCOL,CELLROW,
     2       LAKES_MAP,INDEX_WITH_LAKES,
     3       CELLCOL_WL,CELLROW_WL)
        WRITE(IOUT40,*) 'NCELNL=',NCELNL
        WRITE(IOUT40,*) 'NCELL=',NCELL
c
c  from cells to triangles: construction of traingles, numbering of 
c  nodes ed elements, assignment of x y z coordinates, assignment of elevation
c  values to triangles
c
         IF (FL3D) THEN
            CALL ASSIGN_DEM(NROW,NCOL,DEM_MAP,INDEX,INDEX_WITH_LAKES,
     1           LAKES_MAP,FACTOR,ELEVATION,ELEVATION_WITH_LAKES)
            
            CALL TRIANGOLI(NROW,NCOL,DELTA_X,DELTA_Y,WEST,SOUTH,
     1           DEM_MAP,ZONE,FACTOR,NNOD,NTRI,DOSTEP,
     2           NCELL_COARSE,NODI,TRIANG,TP2D,X,Y,Z,ELTRIA,
     3           CELL)        
c
c     reading of grid parameters
c
            WRITE(IOUT40,1020) NNOD,NTRI,NZONE,NSTR,N1
            WRITE(IOUT40,1350) IVERT,ISP,BASE
            READ(IIN11,*) (ZRATIO(I),I=1,NSTR)
            SUMZ=0.0D0
            DO I=1,NSTR
               SUMZ=SUMZ + ZRATIO(I)
            END DO
            WRITE(IOUT40,1340) (I,ZRATIO(I),I=1,NSTR)
            IF (DABS(SUMZ - 1.0D0) .GT. 1.0D-14) THEN
               WRITE(IOUT2,1345)
               CALL CLOSIO
               STOP
            END IF
c     IF (IPRT1 .GE .1) WRITE(IOUT40,1360) (I,Z(I),I=1,NNOD)
         END IF
      END IF
      N=NNOD*(NSTR + 1)
      NT=3*NTRI*NSTR

      IF (SURF) THEN
         read(IIN17,*) num_tot_r
         write(IOUT40,*) 'numero serbatoi=',num_tot_r
         if (num_tot_r.ne.0) then
            write(IOUT40,*) 'index_c  tipo_r(index_c) reservr(index_c)'
            do j=1,ncelnl
               read(IIN17,*) i,tipo_r(i),reservr(i)
            end do
            write(IOUT40,*) 'livelli iniziali nei serbatoi'
c            do j=1,num_tot_r
c               read(IIN18,*) i,h_pool_kkp1_vec(i)
c               h_pool_kkp1_vec(i)=0
c               write(IOUT40,*) i,h_pool_kkp1_vec(i)
c            end do
         else
            do index_c=1,maxcel
               tipo_r(index_c)=0
               reservr(index_c)=0
            end do
         end if
         if (num_tot_r.ne.0) then
            write(iout40,*) 'caratteristiche',num_tot_r,' serbatoi'
            do i=1,num_tot_r
               read(iin18,*)ilake,n_hA(ilake)
               write(iout40,*) ilake,n_hA(ilake)
               pp=n_hA(ilake)
               write(iout40,*) 'n_ha=',pp
               do k=1,pp
                  read(iin18,*)hres(ilake,k),Ares(ilake,k)
                  write(iout40,*) hres(ilake,k),Ares(ilake,k)   
               end do
               h_fondo(ilake)=hres(ilake,2)
               write(iout40,*) 'h_fondo=',h_fondo(ilake)
               read(iin18,*) h_soglia(ilake)
               write(iout40,*) 'h_soglia=',h_soglia(ilake)
               read(iin18,*) h_pool_kkp1_vec(i)
ccc            h_pool_kkp1_vec(i)=h_fondo(ilake)
               write(iout40,*) 'liv.iniz.=',h_pool_kkp1_vec(i)
            end do
         end if
                 
         
         READ(IIN23,*)
         DO I=1,NCELL
            READ(IIN23,*) QOI_SN(I)
         END DO
    
         CALL RAST_INPUT_REAL(IIN25,NROW,NCOL,NORTH,SOUTH,EAST,WEST,
     1        DTM_W_1)

         CALL RAST_INPUT_REAL(IIN26,NROW,NCOL,NORTH,SOUTH,EAST,WEST,
     1        DTM_W_2)
        
         CALL RAST_INPUT_INT(IIN27,NROW,NCOL,NORTH,SOUTH,EAST,WEST,
     1        DTM_P_OUTFLOW_1)
        
         CALL RAST_INPUT_INT(IIN28,NROW,NCOL,NORTH,SOUTH,EAST,WEST,
     1        DTM_P_OUTFLOW_2)

         CALL RAST_INPUT_REAL(IIN29,NROW,NCOL,NORTH,SOUTH,EAST,WEST,
     1        DTM_LOCAL_SLOPE_1)
         
         CALL RAST_INPUT_REAL(IIN30,NROW,NCOL,NORTH,SOUTH,EAST,WEST,
     1        DTM_LOCAL_SLOPE_2)

         CALL RAST_INPUT_REAL(IIN31,NROW,NCOL,NORTH,SOUTH,EAST,WEST,
     1        DTM_EPL_1)
         
         CALL RAST_INPUT_REAL(IIN32,NROW,NCOL,NORTH,SOUTH,EAST,WEST,
     1        DTM_EPL_2)
               
         CALL RAST_INPUT_REAL(IIN33,NROW,NCOL,NORTH,SOUTH,EAST,WEST,
     1        DTM_KSS1_SF_1)
         
         CALL RAST_INPUT_REAL(IIN34,NROW,NCOL,NORTH,SOUTH,EAST,WEST,
     1        DTM_KSS1_SF_2)

         CALL RAST_INPUT_REAL(IIN35,NROW,NCOL,NORTH,SOUTH,EAST,WEST,
     1        DTM_WS1_SF_1)
         
         CALL RAST_INPUT_REAL(IIN36,NROW,NCOL,NORTH,SOUTH,EAST,WEST,
     1        DTM_WS1_SF_2)

         CALL RAST_INPUT_REAL(IIN37,NROW,NCOL,NORTH,SOUTH,EAST,WEST,
     1        DTM_B1_SF)

         CALL RAST_INPUT_REAL(IIN38,NROW,NCOL,NORTH,SOUTH,EAST,WEST,
     1        DTM_Y1_SF)  

         CALL RAST_INPUT_REAL(IIN39,NROW,NCOL,NORTH,SOUTH,EAST,WEST,
     1        DTM_NRC)
      END IF
C     
      IF (FL3D) THEN
C
C     unit IIN5 input
C
         READ(IIN5,*) INDP,IPOND
         IF (.NOT. SURF) IPOND = 0
         IF ((INDP .EQ. 3).OR.(INDP .EQ. 4)) THEN
            READ(IIN5,*) WTPOSITION
         END IF
         IF (INDP .EQ. 0) THEN
            READ(IIN5,*) PTIMEP(1)
            WRITE(IOUT2,1070) PTIMEP(1)
            DO K=2,N
               PTIMEP(K)=PTIMEP(1)
            END DO
         ELSE IF (INDP .EQ. 1) THEN
            READ(IIN5,*) (PTIMEP(K),K=1,N)
         END IF 
         IF (IPOND .EQ. 1) THEN
            READ(IIN5,*) PONDNOD(1)
            WRITE(IOUT2,1071) PONDNOD(1)   
            IF (PONDNOD(1) .GT. 0.0D0) PTIMEP(1) = PONDNOD(1)
            DO K=2,NNOD
               PONDNOD(K)=PONDNOD(1)   
               IF (PONDNOD(K) .GT. 0.0D0) PTIMEP(K) = PONDNOD(K)
            END DO
         ELSE IF (IPOND .EQ. 2) THEN   
            READ(IIN5,*) (PONDNOD(K),K=1,NNOD)   
            DO K=1,NNOD
               IF (PONDNOD(K) .GT. 0.0D0) PTIMEP(K) = PONDNOD(K)
            END DO
            IF (INDP .EQ. 0 .AND. IPRT1 .GE. 1) THEN
               WRITE(IOUT2,1072) (K,PONDNOD(K),K=1,NNOD)
            END IF
         END IF           
         IF (INDP .EQ. 1 .AND. IPRT1 .GE. 1) THEN
            IF (IPOND.EQ.0) THEN
               WRITE(IOUT2,1090) (K,PTIMEP(K),K=1,N)
            ELSE
               WRITE(IOUT2,1091) (K,PTIMEP(K),K=1,N)
            END IF
         END IF
C
C     unit IIN4 input
C
         READ(IIN4,*) PMIN
         READ(IIN4,*) IPEAT
         READ(IIN4,*) CBETA0,CANG
c     peat soil deformation is not yet supported for the Newton
c     scheme
         IF (IOPT .NE. 1  .AND.  IPEAT .EQ. 1) THEN
            WRITE(IOUT2,1910)
            CALL CLOSIO
            STOP
         END IF
c     peat soil deformation is not yet supported for the EnKF and
c     SIR DA schemes
         IF (DAFLAG .GT. 0  .AND.  IPEAT .EQ. 1) THEN
            WRITE(IOUT2,1920)
            CALL CLOSIO
            STOP
         END IF
         READ(IIN4,*) IVGHU
c     peat soil deformation is not yet supported for the chord slope
c     and tangent slope differentiation of moisture curves
         IF (KSLOPE .NE. 0  .AND.  IPEAT .EQ. 1) THEN
            WRITE(IOUT2,1930)
            CALL CLOSIO
            STOP
         END IF
c     peat soil deformation is not yet supported for the extended
c     van Genuchten or Huyakorn moisture curves (unless they are
c     given in lookup table form)
         IF ( (IPEAT .EQ. 1)  .AND.
     1        (IVGHU .EQ. 1 .OR. IVGHU .EQ. 2 .OR. IVGHU .EQ. 3) ) THEN
            WRITE(IOUT2,1940)
            CALL CLOSIO
            STOP
         END IF
c     moisture curve lookup table option is not yet supported for
c     the EnKF and SIR DA schemes
         IF (DAFLAG .GT. 0  .AND.  IVGHU .EQ. -1) THEN
            WRITE(IOUT2,1950)
            CALL CLOSIO
            STOP
         END IF
c     for IVGHU = -1 the following "VG"/"HU"/"BC" parameters are not
c     needed but are read in anyway. Arbitrary values can be assigned.
         READ(IIN4,*) VGN,VGRMC,VGPSAT
         READ(IIN4,*) HUALFA,HUBETA,HUGAMA,HUPSIA,HUSWR
         READ(IIN4,*) HUN
         READ(IIN4,*) HUA,HUB
         READ(IIN4,*) BCBETA,BCRMC,BCPSAT
C
         IF (IVGHU .EQ. -1) THEN
            READ(IIN16,*) NLKP
            IF (NLKP .LT. 3) THEN
               WRITE(IOUT2,1960)
               CALL CLOSIO
               STOP
            END IF
            DO I=1,NSTR
               DO J=1,NZONE
                  DO K=1,NLKP
                     READ(IIN16,*) PCAP(I,J,K),SATC(I,J,K),KRWC(I,J,K)
                  END DO
               END DO
            END DO
            FLAGASC=0
            DO I=1,NSTR
               DO J=1,NZONE
                  DO K=2,NLKP
                     IF (PCAP(I,J,K) .LE. PCAP(I,J,K-1)) FLAGASC=1
                  END DO
               END DO
            END DO
            IF (FLAGASC .EQ. 1) THEN
               WRITE(IOUT2,1970)
               CALL CLOSIO
               STOP
            END IF
         END IF
C
         WRITE(IOUT2,1215) PMIN
         WRITE(IOUT2,1216) IPEAT
         IF (IPEAT .EQ. 1) THEN
            WRITE(IOUT2,1217) CBETA0,CANG
         END IF
         WRITE(IOUT2,1220) IVGHU
         IF (IVGHU .EQ. 0 .OR. IVGHU .EQ. 1) THEN
            WRITE(IOUT2,1230) VGN,VGRMC,VGPSAT
         ELSE IF (IVGHU .EQ. 2 .OR. IVGHU .EQ. 3) THEN
            WRITE(IOUT2,1240) HUALFA,HUBETA,HUGAMA,HUPSIA,HUSWR
            IF (IVGHU .EQ. 2) THEN
               WRITE(IOUT2,1250) HUN
            ELSE
               WRITE(IOUT2,1260) HUA,HUB
            END IF
         ELSE IF (IVGHU .EQ. 4) THEN
            WRITE(IOUT2,1235) BCBETA,BCRMC,BCPSAT
         END IF
c
         WRITE(IOUT2,1100)
         DO I=1,NSTR
            DO J=1,NZONE
              READ(IIN4,*) PERMX(I,J),PERMY(I,J),PERMZ(I,J),ELSTOR(I,J),
     1              POROS(I,J)
            IF(IPRT1.GE.2) WRITE(IOUT2,1110) I,J,PERMX(I,J), PERMY(I,J),
     1              PERMZ(I,J),ELSTOR(I,J),POROS(I,J)
            END DO
         END DO
      END IF
      
C
C  unit IIN40 input
C
      IF (DAFLAG.NE.0) THEN
         ENKF=.TRUE.
         READ(IIN40,*) NENS,NOBS,NENSMIN
         READ(IIN40,*) DSRETC,DSKS,DSSTOR,DSPOROS,DSSURF,DSIC,DSATM
         READ(IIN40,*) DSKSZ
         READ(IIN40,*) NEFFMIN,DSMEASMAX
         READ(IIN40,*) ERT_FLAG
         DSRETC=DSRETC/1.0D+2
         DSKS=DSKS/1.0D+2
         DSKSZ=DSKSZ/1.0D+2
         DSSTOR=DSSTOR/1.0D+2
         DSPOROS=DSPOROS/1.0D+2
         DSSURF=DSSURF/1.0D+2
         DSIC=DSIC/1.0D+2
         DSATM=DSATM/1.0D+2
         DSMEASMAX=DSMEASMAX/1.0D+2
         READ(IIN40,*) PSISTAR,QSTAR,ATMTAU,DTOUT
         READ(IIN40,*) ISEED
c         write(IOUT56,*) ISEED,'ISEED datin'
         READ(IIN40,*) ENKFT
         READ(IIN40,*) (ENKFTIM(I),I=1,ENKFT)
         IF (NOBS.EQ.0) GO TO 666
         DO I=1,NOBS
            READ(IIN40,*) (ENKFNOD(I,J),J=1,2)
         END DO
         DO I=1,ENKFT
            READ(IIN40,*) DSMEAS(I)
            READ(IIN40,*) (ENKFVAL(I,J),J=1,NOBS)
            DSMEAS(I)=DSMEAS(I)/1.0d+2
         END DO
      ELSE
         ENKF=.FALSE.
      END IF

C PLANT ------------------------------------
c ------------------------------------------
C SISTEMARE FORMATI PER LA STAMPA SUL risul!
c ------------------------------------------

      READ(IIN61,*) NPLANT, NPTYPE
      IF (IPRT1.GE.0) THEN
         WRITE(IOUT2,*) 'PLANT PARAMETERS'
         WRITE(IOUT2,'(A8,I6,2X,A8,I6)') 'NPLANT =',NPLANT,
     1'NPTYPE =',NPTYPE
      ENDIF
      IF (NPLANT.eq.0) THEN
         PLANT_FLAG=.FALSE.
      ELSE
c ------------------------------------------ PLANT_FLAG
         PLANT_FLAG=.TRUE.
         READ(IIN61,*) ACANOFLAG
         READ(IIN61,*) PLANT_PRINT
         READ(IIN61,*) NMETEODATA
         IF (IPRT1.GE.0) THEN
            WRITE(IOUT2,'(A12,I1)') 'NMETEODATA =',NMETEODATA
         ENDIF
         READ(IIN61,*) RHOW
         IF (IPRT1.GE.0) THEN
            WRITE(IOUT2,*) 'RHOW =',RHOW
         ENDIF
         READ(IIN61,*) COATM,CCO2STAR,CCO2ATM
         IF (IPRT1.GE.0) THEN
            WRITE(IOUT2,*)'COATM,CCO2STAR,CCO2ATM=',COATM,CCO2STAR,
     1      CCO2ATM
         ENDIF
         READ(IIN61,*) SSTOMA,ASTOMA
         IF (IPRT1.GE.0) THEN
            WRITE(IOUT2,*)'SSTOMA,ASTOMA =',SSTOMA,ASTOMA
         ENDIF
         READ(IIN61,*) TOLLNR,ITMAXNR,PSILEAF0
         IF (IPRT1.GE.0) THEN
            WRITE(IOUT2,*)'TOLLNR,ITMAXNR,PSILEAF0=',TOLLNR,ITMAXNR,
     1      PSILEAF0
         ENDIF
         READ(IIN61,*) PSTEP,NPMED
         IF (IPRT1.GE.0) THEN
            WRITE(IOUT2,*)'PSTEP,NPMED =',PSTEP,NPMED
         ENDIF
         READ(IIN61,*)(INODP(I),ITYPEP(I),I=1,NPLANT)
         IF (IPRT1.GE.0) THEN
            WRITE(IOUT2,*)'INODP,ITYPEP =',(INODP(I),ITYPEP(I),
     1      I=1,NPLANT)
         ENDIF
         READ(IIN61,*)(XMAXP(I),I=1,NPTYPE)
         READ(IIN61,*)(YMAXP(I),I=1,NPTYPE)
         READ(IIN61,*)(ZMAXP(I),I=1,NPTYPE)
         IF (IPRT1.GE.0) THEN
            WRITE(IOUT2,*)'XMAXP,YMAXP,ZMAXP =',(XMAXP(I),YMAXP(I),
     1      ZMAXP(I),I=1,NPTYPE)
         ENDIF
         READ(IIN61,*)(VRUX(I),I=1,NPTYPE)
         READ(IIN61,*)(VRUY(I),I=1,NPTYPE)
         READ(IIN61,*)(VRUZ(I),I=1,NPTYPE)
         IF (IPRT1.GE.0) THEN
            WRITE(IOUT2,*)'VRUX,VRUY,VRUZ =',(VRUX(I),VRUY(I),
     1      VRUZ(I),I=1,NPTYPE)
         ENDIF
         READ(IIN61,*)(VRUX1(I),I=1,NPTYPE)
         READ(IIN61,*)(VRUY1(I),I=1,NPTYPE)
         READ(IIN61,*)(VRUZ1(I),I=1,NPTYPE)
         IF (IPRT1.GE.0) THEN
            WRITE(IOUT2,*)'VRUX1,VRUY1,VRUZ1 =',(VRUX1(I),VRUY1(I),
     1      VRUZ1(I),I=1,NPTYPE)
         ENDIF
         READ(IIN61,*)(LAI(I),I=1,NPTYPE)
         IF (IPRT1.GE.0) THEN
            WRITE(IOUT2,*)'LAI =',(LAI(I),I=1,NPTYPE)
         ENDIF
         READ(IIN61,*)(GXYLEM_MAX(I),I=1,NPTYPE)
         READ(IIN61,*)(C_GXYLEM(I),I=1,NPTYPE)
         READ(IIN61,*)(D_GXYLEM(I),I=1,NPTYPE)
         IF (IPRT1.GE.0) THEN
            WRITE(IOUT2,*)'GXYLEM_MAX,C_GXYLEM,D_GXYLEM =',
     1      (GXYLEM_MAX(I),C_GXYLEM(I),D_GXYLEM(I),I=1,NPTYPE)
         ENDIF
         READ(IIN61,*)(VCMAX25(I),I=1,NPTYPE)
         READ(IIN61,*)(KC25(I),I=1,NPTYPE)
         READ(IIN61,*)(KO25(I),I=1,NPTYPE)
         READ(IIN61,*)(COMP25(I),I=1,NPTYPE)
         IF (IPRT1.GE.0) THEN
            WRITE(IOUT2,*)'VCMAX25,KC25,KO25,COMP25 =',(VCMAX25(I),
     1      KC25(I),KO25(I),COMP25(I),I=1,NPTYPE)
         ENDIF
         READ(IIN61,*)(LA_MAX(I),I=1,NPTYPE)
         READ(IIN61,*)(LA_BETA(I),I=1,NPTYPE)
         READ(IIN61,*)(LA_PSILMAX(I),I=1,NPTYPE)
         IF (IPRT1.GE.0) THEN
            WRITE(IOUT2,*)'LA_MAX,LA_BETA,LA_PSILMAX =',(LA_MAX(I),
     1      LA_BETA(I),LA_PSILMAX(I),I=1,NPTYPE)
         ENDIF
         READ(IIN61,*)(PSILMAX(I),I=1,NPTYPE)
         READ(IIN61,*)(LIMIT(I),I=1,NPTYPE)
         IF (IPRT1.GE.0) THEN
            WRITE(IOUT2,*)'PSILMAX,LIMIT=',(PSILMAX(I),LIMIT(I),
     1      I=1,NPTYPE)
         ENDIF
         READ(IIN61,*)(HCANO(I),I=1,NPTYPE)
         READ(IIN61,*)(DATA_LEAF(I),I=1,NPTYPE)
         READ(IIN61,*)(ACANO(I),I=1,NPTYPE)
         READ(IIN61,*)(AXYLEM(I),I=1,NPTYPE)
         IF (IPRT1.GE.0) THEN
            WRITE(IOUT2,*)'HCANO,DATALEAF,ACANO,AXYLEM =',(HCANO(I),
     1      DATA_LEAF(I),ACANO(I),AXYLEM(I),I=1,NPTYPE)
         ENDIF
         READ(IIN61,*)(GROOT_STAR(I),I=1,NPTYPE)
         READ(IIN61,*)(GSOIL_STAR(I),I=1,NPTYPE)
         IF (IPRT1.GE.0) THEN
            WRITE(IOUT2,*)'GROOT_STAR,GSOIL_STAR =',(GROOT_STAR(I),
     1      GSOIL_STAR(I),I=1,NPTYPE)
         ENDIF
         READ(IIN61,*)(SALT_TOX(I),I=1,NPTYPE)
         IF (IPRT1.GE.0) THEN
            WRITE(IOUT2,*)'SALT_TOX=',(SALT_TOX(I),
     1      I=1,NPTYPE)
         ENDIF
         READ(IIN61,*) NDATA_RDF
         DO I=1,NDATA_RDF
            READ(IIN61,*) ZRDF(I),(RDFVAL(I,J),J=1,NPTYPE)
         ENDDO
c      ENDIF
  
C     Read Soil Salinity

c  Old plant_salt input file
c      READ(IIN63,*) NSTRSALT
c      IF (NSTRSALT.eq.0) THEN
c         SALT_FLAG=.FALSE.
c      ELSE
c         SALT_FLAG=.TRUE.
c         DO I=1,NSTRSALT
c            READ(IIN63,*) STRSALT(I)
c            DO J=1,NNOD
c               READ(IIN63,*) SALT(I,J)
c            ENDDO
c         ENDDO
c      ENDIF

c  New plant_salt input file
      READ(IIN63,*) NSALTPLANT
      IF(NSALTPLANT.EQ.0) THEN
         SALT_FLAG=.FALSE.
      ELSE
         SALT_FLAG=.TRUE.
         READ(IIN63,*) SALT_C,SALT_S
         DO I=1,NSALTPLANT
            READ(IIN63,*) SALINITY(I)
         ENDDO
      ENDIF


C     Read Plant_Growth parameters  

      READ(IIN64,*) GFLAG
        IF (GFLAG.EQ.1) THEN
              GROWTH_FLAG = .TRUE.
              READ(IIN64,*) LDAY
              READ(IIN64,*) GROWIN  
              READ(IIN64,*) NGROWDAY
              READ(IIN64,*) IDVS 
              READ(IIN64,*)(ILAI(I),I=1,NPTYPE) 
              READ(IIN64,*)(WLVI(I),I=1,NPTYPE)    
              READ(IIN64,*)(WSTI(I),I=1,NPTYPE)  
              READ(IIN64,*)(WRTI(I),I=1,NPTYPE)  
              READ(IIN64,*)(WLVDI(I),I=1,NPTYPE)  
              READ(IIN64,*)(WSOI(I),I=1,NPTYPE)  
              READ(IIN64,*)(IEAI(I),I=1,NPTYPE)  
              READ(IIN64,*)(CFLV(I),I=1,NPTYPE)  
              READ(IIN64,*)(CFST(I),I=1,NPTYPE)  
              READ(IIN64,*)(CFRT(I),I=1,NPTYPE)  
              READ(IIN64,*)(CFSO(I),I=1,NPTYPE)  
              READ(IIN64,*)(Q10(I),I=1,NPTYPE)  
              READ(IIN64,*)(TREF(I),I=1,NPTYPE)  
              READ(IIN64,*)(MAINLV(I),I=1,NPTYPE)  
              READ(IIN64,*)(MAINST(I),I=1,NPTYPE)  
              READ(IIN64,*)(MAINRT(I),I=1,NPTYPE)  
              READ(IIN64,*)(MAINSO(I),I=1,NPTYPE)  
              READ(IIN64,*)(ASRQRT(I),I=1,NPTYPE)  
              READ(IIN64,*)(ASRQLV(I),I=1,NPTYPE)  
              READ(IIN64,*)(ASRQST(I),I=1,NPTYPE)  
              READ(IIN64,*)(ASRQSO(I),I=1,NPTYPE)  
              READ(IIN64,*)(FRTRL(I),I=1,NPTYPE)  
              READ(IIN64,*)(CONVL(I),I=1,NPTYPE)  
              READ(IIN64,*)(FRDR(I),I=1,NPTYPE)  
              READ(IIN64,*)(LAICR(I),I=1,NPTYPE)  
              READ(IIN64,*)(RGRL(I),I=1,NPTYPE)  
              READ(IIN64,*)(SLA(I),I=1,NPTYPE)  
              READ(IIN64,*)(EAR(I),I=1,NPTYPE)  
              READ(IIN64,*)(TBASE(I),I=1,NPTYPE)  
        ELSE
              GROWTH_FLAG = .FALSE.
        ENDIF
C     Read ABL model input data
        READ(IIN65,*) FLAG_ABL
        IF(FLAG_ABL.EQ.0) THEN
          ABL_FLAG = .FALSE.
        ELSE
          ABL_FLAG = .TRUE.
          READ(IIN65,*) ABL_MODEL
          READ(IIN65,*) ZABL0
          READ(IIN65,*) BETA
          READ(IIN65,*) LAPSE
          READ(IIN65,*) LATHEAT
          READ(IIN65,*) CP
          READ(IIN65,*) RHO_A
          READ(IIN65,*) RHO_V
          WRITE(7777,*) ABL_MODEL,ZABL0,BETA,LAPSE,LATHEAT,CP,RHO_A,
     1                  RHO_V
          IF (ABL_MODEL.GE.1) THEN
             READ(IIN65,*) ALBEDO
             READ(IIN65,*) EPS_S
             READ(IIN65,*) SIGMA
             READ(IIN65,*) ZM_ABL
             READ(IIN65,*) ZH
          WRITE(7777,*) ALBEDO,EPS_S,SIGMA,ZM_ABL,ZH
          ENDIF  
          IF(ABL_MODEL.GT.1) READ(IIN65,*)LAPSEW
        ENDIF
      ENDIF
c ------------------------------------------ PLANT_FLAG


666   CONTINUE
C
      WRITE(IOUT2,1300) N,NT
C

      RETURN
C
C  format statements
C
 1000 FORMAT(/,5X,'ISIMGR (0 FLOW3D only w/ grid input, ',
     1       /,5X,'        1 FLOW3D only w/ DEM input, ',
     2       /,5X,'        2 FLOW3D and SURF_ROUTE w/ DEM) = ',I6,
     3       /,5X,'PONDH_MIN (MIN. PONDING HEAD)           = ',1PE15.5,
     4       /,5X,'DAFLAG                                  = ',I6)
 1004 FORMAT(  5X,'TETAF  (EG: 1 BACKWARD EULER, 0.5 C-N)  = ',1PE15.5,
     1       /,5X,'LUMP   (MASS LUMPING IF NONZERO)        = ',I6,
     2       /,5X,'IOPT   (1 PICARD, 2 NEWTON)             = ',I6)
 1005 FORMAT(/,5X,'ISOLV  (-5 BiCGSTAB w/ diag precond, ',
     1       /,5X,'        -4 BiCGSTAB without precond, ',
     2       /,5X,'        -3 TFQMR   w/ diag precond, ',
     3       /,5X,'        -2 TFQMR   without precond, ',
     4       /,5X,'        -1 TFQMR   w/ K^-1  precond, ',
     5       /,5X,'         0 BiCGSTAB w/ K^-1  precond, ',
     6       /,5X,'         1 GRAMRB (min residual), ',
     7       /,5X,'         2 GCRK(5) (ORTHOMIN), ',
     8       /,5X,'         3 NONSYM (direct solver))      = ',I6,
     9       /,5X,'ITMXCG (MAX ITER FOR CG LINEAR SOLVERS) = ',I6,
     A       /,5X,'TOLCG  (TOLER. FOR CG LINEAR SOLVERS)   = ',1PE15.5)
 1010 FORMAT(  5X,'IPRT1  (FOR OUTPUT OF INPUT DATA)       = ',I6,
     1       /,5X,'IPRT   (FOR ELEM VEL & DET NODAL OUTPUT)= ',I6,
     2       /,5X,'NPRT   (# OF TIME VALUES FOR DET OUTPUT)= ',I6,
     3       /,5X,'NUMVP  (# OF SURF NODES FOR VP OUTPUT)  = ',I6,
     4       /,5X,'NR     (# OF NODES FOR PARTIAL OUTPUT)  = ',I6)
 1015 FORMAT(/,5X,'KSLOPE (0 ANA, 1 KSL/ANA, 2 KSL/C-DIFF,',
     1       /,5X,'        3 LOC KSL/ANA, 4 LOC TAN-SLOPE) = ',I6,
     2       /,5X,'TOLKSL (TOLERANCE FOR CHORD SLOPE)      = ',1PE15.5)
 1020 FORMAT(/,5X,'NNOD   (# OF NODES IN 2-D MESH)         = ',I6,
     1       /,5X,'NTRI   (# OF TRIANGLES IN 2-D MESH)     = ',I6,
     2       /,5X,'NZONE  (NUMERO ZONE (MATERIAL TYPES))   = ',I6,
     3       /,5X,'NSTR   (NUMERO STRATI)                  = ',I6,
     4       /,5X,'N1     (NUM. MAX CONTATTI NODALI)       = ',I6)
 1025 FORMAT(  5X,'ITUNS  (MAX NONLINEAR ITER / TIME STEP) = ',I6,
     1       /,5X,'ITUNS1 (DELTAT INCREASE THRESHOLD)      = ',I6,
     2       /,5X,'ITUNS2 (DELTAT DECREASE THRESHOLD)      = ',I6)
 1030 FORMAT(/,5X,'DELTAT (INITIAL TIME STEP SIZE)         = ',1PE15.5,
     1       /,5X,'DTMIN  (MINIMUM TIME STEP SIZE)         = ',1PE15.5,
     2       /,5X,'DTMAX  (MAXIMUM TIME STEP SIZE)         = ',1PE15.5,
     3       /,5X,'TMAX   (TIME AT END OF SIMULATION)      = ',1PE15.5)
 1035 FORMAT(  5X,'DTMAGA (MAG. FACTOR FOR DELTAT, ADD.)   = ',1PE15.5,
     1       /,5X,'DTMAGM (MAG. FACTOR FOR DELTAT, MULT.)  = ',1PE15.5,
     2       /,5X,'DTREDS (RED. FACTOR FOR DELTAT, SUB.)   = ',1PE15.5,
     3       /,5X,'DTREDM (RED. FACTOR FOR DELTAT, MULT.)  = ',1PE15.5)
 1036 FORMAT(  5X,'DTOUT  (TIME STEP FOR OUTPUTS)          = ',1PE15.5)
 1040 FORMAT(  5X,'PKRL   (LEFT  PrHead ENDPT FOR d(Kr)/dP)= ',1PE15.5,
     1       /,5X,'PKRR   (RIGHT PrHead ENDPT FOR d(Kr)/dP)= ',1PE15.5,
     2       /,5X,'PSEL   (LEFT  PrHead ENDPT FOR d(Se)/dP)= ',1PE15.5,
     3       /,5X,'PSER   (RIGHT PrHead ENDPT FOR d(Se)/dP)= ',1PE15.5,
     4       /,5X,'PDSE1L (LEFT  P FOR dd(Se)/dPP, RANGE 1)= ',1PE15.5,
     5       /,5X,'PDSE1R (RIGHT P FOR dd(Se)/dPP, RANGE 1)= ',1PE15.5,
     6       /,5X,'PDSE2L (LEFT  P FOR dd(Se)/dPP, RANGE 2)= ',1PE15.5,
     7       /,5X,'PDSE2R (RIGHT P FOR dd(Se)/dPP, RANGE 2)= ',1PE15.5)
 1045 FORMAT(/,5X,'ISFONE (0 ALL-NODE SFEX UPD, 1 1-NODE)  = ',I6,
     1       /,5X,'ISFCVG (0 NL CONVG W/O SFEX, 1 W/ SFEX) = ',I6)
 1050 FORMAT(  5X,'NLRELX (0 NORELX,1 CONS RELX,2 VAR RELX)= ',I6)
 1055 FORMAT(  5X,'OMEGA  (RELAXATION PARAMETER, CONSTANT) = ',1PE15.5)
 1060 FORMAT(  5X,'L2NORM (0 INFINITY NORM, ELSE L2 NORM)  = ',I6,
     1       /,5X,'TOLUNS (TOLERANCE FOR NONLINEAR ITER)   = ',1PE15.5,
     2       /,5X,'TOLSWI (TOLERANCE FOR BC SWITCHING)     = ',1PE15.5,
     3       /,5X,'ERNLMX (MAX ALLOWABLE CVG OR RESID ERR) = ',1PE15.5)
 1070 FORMAT(/,5X,'CONSTANT INITIAL PRESSURE HEAD          = ',1PE15.5)
 1071 FORMAT(/,5X,'CONSTANT INITIAL PONDING HEAD           = ',1PE15.5)
 1072 FORMAT(/,1X,'INITIAL PONDING HEAD',/,(4(I6,2X,1PE11.3)))
 1090 FORMAT(/,1X,'INITIAL PRESSURE HEAD',/,(4(I6,2X,1PE11.3)))
 1091 FORMAT(/,1X,'INITIAL PONDING+PRESSURE HEAD',/,(4(I6,2X,1PE11.3)))
 1100 FORMAT(//,5X,
     1   'SATURATED HYDRAULIC CONDUCTIVITY, SPECIFIC STORAGE, AND ',
     2   'POROSITY VALUES',/,
     3   1X,' LAYER MAT.TYPE  X-PERM       Y-PERM       Z-PERM',
     4      '       STORAGE      POROSITY')
 1110 FORMAT(1X,I4,I8,2X,5(1PE13.5))
C1120 FORMAT(/,5X,'NDIR (# OF NON-ATM, NON-SF DIR NODES 2D)= ',I6)  
C1150 FORMAT(/,5X,'NON-ATM, NON-SF DIR NODES IN 2-D MESH')
C1155 FORMAT(/,5X,'NATM, NSF FIXED DIR NODES IN 3-D MESH')
C1160 FORMAT(1X,10I7)
C1170 FORMAT(/,5X,'NDIRC(# OF FIXED NATM,NSF DIR NODES 3-D)= ',I6)  
C1175 FORMAT(/,5X,'NP   (TOTAL # OF NATM,NSF DIR NODES 3-D)= ',I6)
C1180 FORMAT(/,5X,'NQ   (# OF NON-ATM, NON-SF NEU NODES 3D)= ',I6)  
C1200 FORMAT(/,5X,'NON-ATM, NON-SF NEU NODES IN 3-D MESH')
 1215 FORMAT(/,5X,'PMIN  (AIR DRY PRESSURE HEAD VALUE)     = ',1PE15.5) 
 1216 FORMAT(  5X,'IPEAT (flag for peat deformation)       = ',I6)
 1217 FORMAT(  5X,'CBETA0                                  = ',1PE15.5,
     1       /,5X,'CANG                                    = ',1PE15.5)
 1220 FORMAT(  5X,'IVGHU (-1 moisture curve lookup table, ',
     1       /,5X,'       0 van Genuchten, ',
     2       /,5X,'       1 extended van Genuchten, ',
     3       /,5X,'       2 Huyakorn with Kr=Se**n, ',
     4       /,5X,'       3 Huyakorn with log_10 Kr(Se), ',
     5       /,5X,'       4 Brooks-Corey)                  = ',I6)
 1230 FORMAT(  5X,'VGN                                     = ',1PE15.5,
     1       /,5X,'VGRMC                                   = ',1PE15.5,
     2       /,5X,'VGPSAT                                  = ',1PE15.5)
 1235 FORMAT(  5X,'BCBETA                                  = ',1PE15.5,
     1       /,5X,'BCRMC                                   = ',1PE15.5,
     2       /,5X,'BCPSAT                                  = ',1PE15.5)
 1240 FORMAT(  5X,'HUALFA                                  = ',1PE15.5,
     1       /,5X,'HUBETA                                  = ',1PE15.5,
     2       /,5X,'HUGAMA                                  = ',1PE15.5,
     3       /,5X,'HUPSIA                                  = ',1PE15.5,
     4       /,5X,'HUSWR                                   = ',1PE15.5)
 1250 FORMAT(  5X,'HUN                                     = ',1PE15.5)
 1260 FORMAT(  5X,'HUA                                     = ',1PE15.5,
     1       /,5X,'HUB                                     = ',1PE15.5)
 1300 FORMAT(/,5X,'N     (# OF NODES IN 3-D MESH)          = ',I8,
     1       /,5X,'NT    (# OF TETRAHEDRA IN 3-D MESH)     = ',I8)
 1340 FORMAT(  5X,'LAYER ',I3,5X,'ZRATIO = ',1PE13.5)
 1345 FORMAT(//,' INPUT ERROR : ZRATIO VALUES MUST SUM TO 1.0')
 1350 FORMAT(/,5X,'IVERT  (TYPE OF VERTICAL DISCRETIZATION)= ',I6,
     1       /,5X,'ISP    (0 FLAT SURFACE, ELSE NOT FLAT)  = ',I6,
     2       /,5X,'BASE   (THICKNESS OR BASE OF 3-D MESH)  = ',1PE15.5)
 1360 FORMAT(/,5X,' SURFACE ELEVATION VALUES'/(4(I5,1X,1PE12.4)))
C1400 FORMAT(/,5X,'NSF   (# OF SEEPAGE FACES)              = ',I6)
C1410 FORMAT(  5X,'NUMBER OF NODES ON SEEPAGE FACE ',I6,'  = ',I6)
C1420 FORMAT(  5X,'NODE #''S : ',10I6)
C1500 FORMAT(//,' INPUT ERROR : ELEVATION VALUES NOT IN DESCENDING',
C    1          ' ORDER ON SEEPAGE FACE ',I6,
C    2       /,4X,'NODES',I4,' (NODE #',I6,') AND',I4,' (NODE #',I6,')')
C1600 format(I5,I5,I5/I5,I5,I5/I5,I5,I5) 
 1900 FORMAT(/,5X,'INCORRECT INPUT VALUE FOR ISIMGR: ',I6)
 1910 FORMAT(/,5X,'INVALID IOPT AND IPEAT COMBINATION (peat soil ',
     1            'deformation is not yet ',
     2       /,5X,'supported for the Newton scheme)')
 1920 FORMAT(/,5X,'INVALID DAFLAG AND IPEAT COMBINATION (peat soil ',
     1            'deformation is not yet ',
     2       /,5X,'supported for EnKF and SIR DA schemes)')
 1930 FORMAT(/,5X,'INVALID KSLOPE AND IPEAT COMBINATION (peat soil ',
     1            'deformation is not yet ',
     2       /,5X,'supported for chord slope and tangent slope ',
     3            'differentiation of moisture curves)')
 1940 FORMAT(/,5X,'INVALID IVGHU AND IPEAT COMBINATION (peat soil ',
     1            'deformation is not yet ',
     2       /,5X,'supported for the extended van Genuchten or ',
     3            'Huyakorn moisture curves ',
     4       /,5X,'(unless they are given in lookup table form))')
 1950 FORMAT(/,5X,'INVALID DAFLAG AND IVGHU COMBINATION (moisture ',
     1            'curve lookup table option ',
     2       /,5X,'is not yet supported for EnKF and SIR DA schemes)')
 1960 FORMAT(/,5X,'NLKP must be at least 3')
 1970 FORMAT(/,5X,'PCAP values must be in ascending order')
      END
