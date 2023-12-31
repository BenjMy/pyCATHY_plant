C
C********************  RAST_INPUT_DLZ ***************************************
C
C  reads a real raster map in GRASS (ASCII) format from unit IIN according 
C  to matrix convenction. This subroutine is used to read the dem and base map 
C  raster input file
C
C*****************************************************************************
C
      SUBROUTINE RAST_INPUT_DEM(IIN,NROW,NCOL,NORTH,SOUTH,EAST,WEST,
     1     RAST)


      IMPLICIT NONE
      INCLUDE 'CATHY.H'
     
      INTEGER  IIN,NROW,NCOL
      INTEGER  IROW,ICOL
      REAL*8 NORTH,SOUTH,EAST,WEST
      CHARACTER RDWR*80,LINE*80
      REAl*8 RAST(ROWMAX,COLMAX)

c
c  reads the first four lines (north,south,east,west)
c
      READ(IIN,'(A80)') LINE
      WRITE(RDWR,'(A74)') LINE(8:80)
      READ(RDWR,*) NORTH
      READ(IIN,'(A80)') LINE
      WRITE(RDWR,'(A74)') LINE(8:80)
      READ(RDWR,*) SOUTH 
      READ(IIN,'(A80)') LINE
      WRITE(RDWR,'(A74)') LINE(8:80)
      READ(RDWR,*) EAST 
      READ(IIN,'(A80)') LINE
      WRITE(RDWR,'(A74)') LINE(8:80)
      READ(RDWR,*) WEST 
c
c  reads in NROW and NCOL
c
      READ(IIN,'(A80)') LINE
      WRITE(RDWR,'(A75)') LINE(7:80)
      READ(RDWR,*) NROW
      READ(IIN,'(A80)') LINE
      WRITE(RDWR,'(A75)') LINE(7:80)
      READ(RDWR,*) NCOL
c  reads the raster map 
      DO IROW=1,NROW
         READ(IIN,*) (RAST(IROW,ICOL),ICOL=1,NCOL)
      END DO

      RETURN
      END
