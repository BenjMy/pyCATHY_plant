C
C************************** SURF_ROUTE *********************************
C
C water surface routing procedure:
C takes OVFCEL as input from FLOW3D and returns `routed' PONDNOD
C as BC for FLOW3D
C 
C***********************************************************************
C
      subroutine surf_route(ncell,nnod,nrow,ncol,ntri,dostep,numres,
     1                      ncell_coarse,ncelnl,cell,
     2                      indcel,indcelwl,tipo_r,reservr,dem_map,
     3                      lakes_map,cellcol,cellrow,
     4                      celtype,cells_r,tp2d,triang,
     5                      time,deltat,arenod,OVFLNOD,PONDNOD,
     6                      ovflcel,pondcel,n_hA,
     7                      NSURF,NSURFT,NSURFT_TB)
C
      implicit none
      include 'CATHY.H'
      integer  ncell,nnod,nrow,ncol,ntri,dostep,numres
      integer  ncell_coarse,ncelnl
      INTEGER  NSURF,NSURFT,NSURFT_TB
      integer  i,j,k
      integer  cell(5,*),indcel(rowmax,ncol),indcelwl(rowmax,ncol)
      integer  tipo_r(*),reservr(*)
      integer  lakes_map(rowmax,ncol)
      integer  cellcol(*),cellrow(*)
      integer  celtype(*),cells_r(*)
      integer  tp2d(*),triang(4,*)
      integer  n_hA(*)
      real*8   time,deltat
      real*8   dtsurf
      real*8   arenod(*)
      real*8   OVFLNOD(*),PONDNOD(*)
      real*8   ovflcel(*),pondcel(*),dem_map(rowmax,ncol)
      INCLUDE 'SURFWATER.H'
      INCLUDE 'IOUNITS.H'
      INCLUDE 'RIVERNETWORK.H'
C
c
c  from node to cell
c    
      call nod_cell(ncell,nrow,ncol,dostep,ncell_coarse,
     1              nnod,cell,dem_map,indcelwl,cellcoarse,OVFLNOD,
     2              ovflcel,arenod,delta_x,delta_y)

c   trasferimento informazioni al surf_route (cambia la numerazione
c   delle celle se ci sono laghi!)
c     
      call transfer_f3d_surf(nrow,ncol,indcel,indcelwl,
     1                       surface_water_sn,ovflcel,
     2                       reservr,lakes_map)
c    
c  surface routing
c  
      cu_max=ak_max*deltat
      CUTRGT = 1.0d0
      if(cu_max .gt. CUTRGT) then
         dtsurf = CUTRGT/ak_max
         nsurf = IDINT((deltat/dtsurf))+1
         dtsurf = deltat/nsurf
      else
         dtsurf = deltat
         nsurf=1
      end if
      NSURFT=NSURFT + NSURF
      NSURFT_TB=NSURFT_TB + NSURF


      
      DO I=1,NSURF
         CALL ROUTE(NROW,NCELL,TIPO_R,RESERVR,
     1        CELLS_R,N_HA,TIME,DTSURF)
C     CALCULATION OF SURFACE WATER HEIGHTS
               CALL ALTEZZE(NCELNL,CELTYPE,DELTA_X,
     1             DELTA_Y,Q_IN_KK_SN,Q_IN_KKP1_SN,
     2             Q_OUT_KK_SN_1,Q_OUT_KK_SN_2,
     3             Q_OUT_KKP1_SN_1,Q_OUT_KKP1_SN_2,
     4             VOLUME_KK_SN,VOLUME_KKP1_SN,SURFACE_WATER_SN,
     5             H_WATER_KKP1_SN,DTSURF,
     6             H_POOL_KKP1_VEC,RESERVR,ELEVATION)

C INTERNAL UPDATES OF SURFACE ROUTING VARIABLES (I.E., FLOW RATE, WATER HEIGHTS)
            IF(NSURF.GT.1 .AND. I.LT.NSURF) THEN
               CALL VCOPYR(MAXCEL,Q_IN_KK_SN,Q_IN_KKP1_SN)
               CALL INIT0R(MAXCEL,Q_IN_KKP1_SN)
               CALL VCOPYR(MAXCEL,Q_OUT_KK_SN_1,Q_OUT_KKP1_SN_1)
               CALL INIT0R(MAXCEL,Q_OUT_KKP1_SN_1)
               CALL VCOPYR(MAXCEL,Q_OUT_KK_SN_2,Q_OUT_KKP1_SN_2)
               CALL INIT0R(MAXCEL,Q_OUT_KKP1_SN_2)
               CALL VCOPYR(MAXCEL,VOLUME_KK_SN,VOLUME_KKP1_SN)
               CALL INIT0R(MAXCEL,VOLUME_KKP1_SN)
C     CALL VCOPYR(NUMRES,H_POOL_KK_VEC,H_POOL_KKP1_VEC)       
            END IF
         END DO
                        
c
c trasferimento informazioni per il flow3d
c
      call transfer_surf_f3d(nrow,ncol,indcel,indcelwl,
     1                       ncell,h_water_kkp1_sn,
     2                       pondcel,h_pool_kkp1_vec,
     3                       lakes_map,elevation_with_lakes)
c
c  calcola PONDNOD da PONDCEL passando per i triangoli
c
      call cell_nod(ncell,nrow,ncol,nnod,ntri,dostep,
     1              tp2d,triang,indcelwl,dem_map,
     2              pondcel,cellcoarse,PONDNOD)
c
      return
      end

