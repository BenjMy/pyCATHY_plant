         subroutine vtkris1(ntria, nnode,iunit,triang,u,x,z)
         implicit none
         integer nnode, ntria, iunit,  i, j, numb
         integer triang(4,ntria)
         real*8 u(nnode)
         real*8 x(nnode), z(nnode)
         character*15 nome

         write(nome,'(i3,a4)') iunit, '.vtk'
         open(iunit, file=nome)
        
         write(iunit,78) 
 78      FORMAT('# vtk DataFile Version 2.0',
     1     /,'2D Unstructured Grid of Linear Triangles', 
     2     /,'ASCII')

        write(iunit,79) nnode
 79     FORMAT('DATASET UNSTRUCTURED_GRID', 
     1    /, 'POINTS',1x, i8, ' float')

         do i=1,nnode
             write(iunit,100) x(i), z(i), 0.d0 
         end do
 100     format(4x, 3(1pe16.4))

         numb=ntria*4
         write(iunit, 80) ntria,numb
 80      format('CELLS',1x,i8,1x,i8)
         j=3
         do i=1,ntria
           write(iunit,101) j, triang(1,i)-1, triang(2,i)-1, 
     1                    triang(3,i)-1
         end do
 101     format(i1, 3x,i8,3x,i8,3x,i8)

         write(iunit,81)  ntria
 81      format('CELL_TYPES',i8)
         do i=1,ntria
          write(iunit,*) 5
         end do

c         write(iunit,82)  ntria
         write(iunit,82)  nnode
c 82      format('CELL_DATA', 1x,i8, 
 82      format('POINT_DATA', 1x,i8, 
     1    /, 'SCALARS atmbc float',
     2    /, 'LOOKUP_TABLE default' )
c         do i=1,ntria
         do i=1,nnode
          write(iunit,*) u(i)
         end do

 
         close(iunit)
         return

         end
         
