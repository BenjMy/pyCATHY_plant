C
C**************************  COMMON BLOCK SOILCHAR *********************
C            
C  real scalars and real arrays for the PMIN parameter, for unsaturated
C  soil characteristics, for chord and tangent slope formulas, for
C  peat soil deformation parameters, and for the moisture curve
C  lookup table option
C            
C***********************************************************************
C
      REAL*8  PMIN
      REAL*8  VGN,VGM,VGRMC,VGPSAT
      REAL*8  VGN1,VGM52,VGNR,VGMM1,VGMR,VGPSN
      REAL*8  HUN,HUA,HUB,HUALFA,HUBETA,HUGAMA,HUPSIA,HUSWR
      REAL*8  HUSWR1,HUALB,HUGAM1,HUGB,HUGB1,HU1BET,HUGAM2
      REAL*8  HUN1,HU2A,HUB2A,HUAB,HULN10
      REAL*8  BCBETA,BCRMC,BCPSAT
      REAL*8  BCB1,BCB2,BCBB1P,BCBPS,BC23B,BC23BP,BC33B
      REAL*8  TOLKSL,PKRL,PKRR,PSEL,PSER
      REAL*8  PDSE1L,PDSE1R,PDSE2L,PDSE2R,DKRTAN
      REAL*8  CBETA0,CANG
      REAL*8  VGPNOT(NMAX),BCPORM(NMAX)
      REAL*8  DSETAN(NMAX),DDSE1T(NMAX),DDSE2T(NMAX)
      REAL*8  PCAP(MAXSTR,MAXZON,MAXLKP),SATC(MAXSTR,MAXZON,MAXLKP)
      REAL*8  KRWC(MAXSTR,MAXZON,MAXLKP)
C
C  common block for PMIN
C
      COMMON  /COMPMIN/PMIN
C
C  common blocks for van Genuchten, extended van Genuchten,
C  Huyakorn, and Brooks-Corey moisture curve parameters
C
      COMMON  /VG/VGPNOT,VGN,VGM,VGRMC,VGPSAT,
     1            VGN1,VGM52,VGNR,VGMM1,VGMR,VGPSN
      COMMON  /HU/HUN,HUA,HUB,HUALFA,HUBETA,HUGAMA,HUPSIA,HUSWR,
     1            HUSWR1,HUALB,HUGAM1,HUGB,HUGB1,HU1BET,HUGAM2,
     2            HUN1,HU2A,HUB2A,HUAB,HULN10
      COMMON  /BC/BCPORM,BCBETA,BCRMC,BCPSAT,
     1            BCB1,BCB2,BCBB1P,BCBPS,BC23B,BC23BP,BC33B
C
C  common block for chord and tangent slope formulas
C
      COMMON  /CHORD/TOLKSL,PKRL,PKRR,PSEL,PSER,
     1               PDSE1L,PDSE1R,PDSE2L,PDSE2R,DKRTAN,
     2               DSETAN,DDSE1T,DDSE2T
C
C  common block for peat soil deformation parameters
C
      COMMON  /CAMPO/CBETA0,CANG
C
C  common block for moisture curve lookup table parameters
C
      COMMON  /TBLKUP/PCAP,SATC,KRWC

