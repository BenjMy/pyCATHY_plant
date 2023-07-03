C
C**************************  RAST_OUTPUT ********************************
C
C  outputs a raster map in GRASS (ASCII) format to unit IOUT
C
C***********************************************************************
C
      subroutine rast_output(iout,rowmax,nrow,ncol,
     1                       north,south,east,west,rast)


      implicit none
     
      integer  iout,rowmax,nrow,ncol
      integer  irow,icol
      integer  rast(rowmax,ncol)
C     real*8 rast(rowmax,ncol)
      real*8 north,south,east,west
      character rdwr*80,line*80
      integer  out_of_dem
      parameter (out_of_dem=0)

c
c  reads the first four lines (north,south,east,west)
c
      write(iout,'(a7,f20.8)') 'north: ',north
      write(iout,'(a7,f20.8)') 'south: ',south
      write(iout,'(a7,f20.8)') 'east:  ',east
      write(iout,'(a7,f20.8)') 'west:  ',west
c
c  reads in NROW and NCOL
c
      write(iout,'(a6,I5)') 'rows: ',nrow
      write(iout,'(a6,I5)') 'cols: ',ncol
c
c  reads the raster map 
c
      do irow=1,nrow
         do icol=1,ncol-1
            if(rast(irow,icol).eq.out_of_dem) then
               write(iout,'(1x,i10,1x,$)') rast(irow,icol)
            else
               write(iout,'(1x,i10,1x,$)') rast(irow,icol)
            end if
         end do
         if(rast(irow,icol).eq.out_of_dem) then
            write(iout,'(1x,i10)') rast(irow,ncol)
         else
            write(iout,'(1x,i10)') rast(irow,ncol)
         end if
      end do
c
      return
      end
