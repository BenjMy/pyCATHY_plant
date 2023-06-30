      program create_nansf

      implicit none
      integer n
      parameter(n=23000)
      integer i 
      integer nnod,nnod3d,nn(n),counter,counter1,counterneu,counterneu1
      real*8  x(n),y(n),z(n),depthMESH,max_x_MESH,max_y_MESH,depthWT
      real*8  zero

C Valori da inserire di volta in volta:
C Profondità della tavola d'acqua dalla superficie

      counter=0
      counter1=0
      counterneu=0
      counterneu1=0
      zero=0.0

C Open input file

      open(21,file='../../output/xyz',status='old')
      open(22,file='nansfdirbc',status='unknown')
      open(23,file='nansfneubc',status='unknown')
      open(24,file='input',status='unknown')

C Read input file
C nnod = numero nodi mesh 2D
C nnod3d = numero nodi mesh 3D

      read(21,*) nnod, nnod3d
      do i=1,nnod3d
         read(21,*) nn(i),x(i),y(i),z(i)
      end do

      read(24,*) depthWT
      read(24,*) depthMESH
      read(24,*) max_x_MESH
      read(24,*) max_y_MESH

C ** Write nansfdirbc and nansfneubc **
C
C Conto nodi per nansfDIRbc e li scrivo nel rispettivo file
     
      write(22,*) 0.0, 'TIME'
      do i=1,nnod3d
         if ((z(i).le.depthWT).and.(z(i).gt.depthMESH)) then
            if ((x(i).eq.zero).or.(x(i).eq.max_x_MESH).or.
     1(y(i).eq.zero).or.(y(i).eq.max_y_MESH)) then
                 counter=counter+1 
            end if
         end if
      end do
      write(22,100) 0, counter
 100  format(1x,i1,2x,i5,5x,'NDIR',3x,'NDIRC')
      do i=1,nnod3d
         if ((z(i).le.depthWT).and.(z(i).gt.depthMESH)) then
            if ((x(i).eq.zero).or.(x(i).eq.max_x_MESH).or.
     1(y(i).eq.zero).or.(y(i).eq.max_y_MESH)) then
                 write(22,*) nn(i)
            end if
         end if
      end do

      do i=1,nnod3d
         if ((z(i).le.depthWT).and.(z(i).gt.depthMESH)) then
            if ((x(i).eq.zero).or.(x(i).eq.max_x_MESH).or.
     1(y(i).eq.zero).or.(y(i).eq.max_y_MESH)) then
C                 write(22,201) nn(i),x(i),y(i),z(i),-z(i)+depthWT
                  write(22,200) -z(i)+depthWT
            end if
         end if
      end do      
C 201  format(1x,i3,f6.2,f6.2,f6.2,f6.2) 
  200  format(1x,f6.2)
      write(22,*) 

      write(22,*) 1.0E+8, 'TIME'
      write(22,100) 0, counter
      do i=1,nnod3d
         if ((z(i).le.depthWT).and.(z(i).gt.depthMESH)) then
            if ((x(i).eq.zero).or.(x(i).eq.max_x_MESH).or.
     1(y(i).eq.zero).or.(y(i).eq.max_y_MESH)) then
                 write(22,*) nn(i)
            end if
         end if
      end do

      do i=1,nnod3d
         if ((z(i).le.depthWT).and.(z(i).gt.depthMESH)) then
            if ((x(i).eq.zero).or.(x(i).eq.max_x_MESH).or.
     1(y(i).eq.zero).or.(y(i).eq.max_y_MESH)) then
                  write(22,200) -z(i)+depthWT
            end if
         end if
      end do      
 
C Conto nodi per nansfNEUbc e li scrivo nel rispettivo file
      write(23,*) 0.0, 'TIME'
      do i=1,nnod3d
         if (z(i).gt.depthWT) then
            if ((x(i).eq.zero).or.(x(i).eq.max_x_MESH).or.
     1(y(i).eq.zero).or.(y(i).eq.max_y_MESH)) then
                 counterneu=counterneu+1
            end if
         elseif (z(i).eq.depthMESH) then
                 counterneu=counterneu+1
         end if                 
      end do
      write(23,101) 0, counterneu
 101  format(1x,i1,2x,i5,5x,'ZERO',3x,'NQ')

      do i=1,nnod3d
         if (z(i).gt.depthWT) then
            if ((x(i).eq.zero).or.(x(i).eq.max_x_MESH).or.
     1(y(i).eq.zero).or.(y(i).eq.max_y_MESH)) then
                 write(23,*) nn(i)
            end if
         elseif (z(i).eq.depthMESH) then
                 write(23,*) nn(i)
         end if                 
      end do

      do i=1,nnod3d
         if (z(i).gt.depthWT) then
            if ((x(i).eq.zero).or.(x(i).eq.max_x_MESH).or.
     1(y(i).eq.zero).or.(y(i).eq.max_y_MESH)) then
                 write(23,201) 0.0
            end if
         elseif (z(i).eq.depthMESH) then
                 write(23,201) 0.0
         end if                 
      end do
  201 format(1x,f6.2)
   
      write(23,*) 
      write(23,*) 1.0E+8, 'TIME'
      write(23,101) 0, counterneu
      do i=1,nnod3d
         if (z(i).gt.depthWT) then
            if ((x(i).eq.zero).or.(x(i).eq.max_x_MESH).or.
     1(y(i).eq.zero).or.(y(i).eq.max_y_MESH)) then
                 write(23,*) nn(i)
            end if
         elseif (z(i).eq.depthMESH) then
                 write(23,*) nn(i)
         end if                 
      end do

      do i=1,nnod3d
         if (z(i).gt.depthWT) then
            if ((x(i).eq.zero).or.(x(i).eq.max_x_MESH).or.
     1(y(i).eq.zero).or.(y(i).eq.max_y_MESH)) then
                 write(23,201) 0.0
            end if
         elseif (z(i).eq.depthMESH) then
                 write(23,201) 0.0
         end if                 
      end do


C Chiuso i file aperti

      close(21)
      close(22)
      close(23)
      close(24)

      end
