      program create_atmbc

      implicit none
      integer n
      parameter(n=23000)
      integer i,j,k
      integer nnod,nnod3d,nn(n),HSPATM,IETO,counter,time_0,time_new
      integer sec_h,nn2d(n),sec_g,giorni_irrigazione,sec_7h
      integer giorni_stop,giorni_irrigazione_2
      real*8  x(n),y(n),z(n),x_min,x_max,y_min,y_max,flux
      real*8  zero,t,areanod2d(n),totarea,time

C Valori da inserire di volta in volta:
C HSPATM e IETO sono parametri del CATHY (si veda la guida per il significato)
C time_h = secondi in un'ora

      HSPATM=0
      IETO=1
      zero=0.0
      time=25200
      counter=0
      totarea=0.0
      sec_h=3600
      sec_g=24*3600
C      sec_7h=25200

C Open input file

      open(21,file='../../output/xyz',status='old')
      open(22,file='atmbc',status='unknown')
      open(23,file='input_atmbc',status='unknown')
      open(24,file='../../output/areanod',status='old')
      open(25,file='area_e_nodo',status='unknown')

C Read input file
C nnod = numero nodi mesh 2D
C nnod3d = numero nodi mesh 3D

      read(21,*) nnod, nnod3d
      do i=1,nnod3d
         read(21,*) nn(i),x(i),y(i),z(i)
      end do

      read(23,*) x_min
      read(23,*) x_max
      read(23,*) y_min
      read(23,*) y_max
      read(23,*) flux
      read(23,300) giorni_irrigazione
      read(23,300) giorni_stop
      read(23,300) giorni_irrigazione_2
 300  format(i3)

      do j=1,nnod
         read(24,*) nn2d(j),areanod2d(j)
      end do


C ** Write atmbc **
C
C Conto nodi per atmbc e li scrivo nel rispettivo file

      if (giorni_stop .ne. 0 .and. giorni_irrigazione_2 .ne. 0) then

          do i=1,nnod3d
             if (z(i).eq.zero) then
                if (((x(i).gt.x_min).and.(x(i).lt.x_max)).and.
     1    ((y(i).gt.y_min).and.(y(i).lt.y_max))) then
                counter=counter+1
                    if (nn(i).eq.nn2d(i)) then
                       write(25,*) nn(i),nn2d(i),areanod2d(i)
                       totarea=totarea+areanod2d(i)
                    end if
                end if
             end if
          end do
          write(*,*) counter
          write(*,*) totarea

          write(22,100) HSPATM,IETO
 100      format(1x,i1,1x,i1,3x,'HSPATM',2x,'IETO')
          write(22,200) zero
 200      format(1x,f16.2,3x,'TIME')
          do i=1,nnod3d
             if (z(i).eq.zero) then
                write(22,*) zero
             end if
          end do

          do k=1,giorni_irrigazione
             write(22,200) time
             do i=1,nnod3d
                if (z(i).eq.zero) then
                   if (((x(i).gt.x_min).and.(x(i).lt.x_max)).and.
     1    ((y(i).gt.y_min).and.(y(i).lt.y_max))) then
                      write(22,*) ((flux)/totarea)
                   else 
                      write(22,*) zero
                   end if
                end if
             end do 
          
             t=time+5*sec_h+1
C          write(*,*) t
             write(22,*) 
             write(22,200) t                
             do i=1,nnod3d
                if (z(i).eq.zero) then
                   write(22,*) zero
                end if
             end do
             write(22,*)      
             time=time+sec_g
             t=t+sec_g
          end do 

C Aggiungo lo stop per 5 giorni e quindi la nuova irrigazione di 10 giorni 

          write(*,*) time
          write(*,*) t
          time=time+sec_g*giorni_stop
          write(22,200) time
          do i=1,nnod3d
             if (z(i).eq.zero) then
                if (((x(i).gt.x_min).and.(x(i).lt.x_max)).and.
     1    ((y(i).gt.y_min).and.(y(i).lt.y_max))) then
                   write(22,*) ((flux)/totarea)
                else 
                   write(22,*) zero
                end if
             end if
          end do 
          write(22,*) 
          t=t+sec_g*giorni_stop
          write(22,200) t                
          do i=1,nnod3d
             if (z(i).eq.zero) then
                write(22,*) zero
             end if
          end do
          write(22,*)      
          time=time+sec_g
          t=t+sec_g
      
          do k=1,giorni_irrigazione_2-1
             write(22,200) time
             do i=1,nnod3d
                if (z(i).eq.zero) then
                   if (((x(i).gt.x_min).and.(x(i).lt.x_max)).and.
     1    ((y(i).gt.y_min).and.(y(i).lt.y_max))) then
                      write(22,*) ((flux)/totarea)
                   else 
                      write(22,*) zero
                   end if
                end if
             end do 
          
             t=time+5*sec_h+1
C          write(*,*) t
             write(22,*) 
             write(22,200) t                
             do i=1,nnod3d
                if (z(i).eq.zero) then
                   write(22,*) zero
                end if
             end do
             write(22,*)      
             time=time+sec_g
             t=t+sec_g
          end do 
          write(*,*) 'atmbc con sospensione'
C Atmbc nel caso in cui non vi sia la sospensione

      elseif (giorni_stop .eq. 0 .and. giorni_irrigazione_2 .eq. 0) then

          do i=1,nnod3d
             if (z(i).eq.zero) then
                if (((x(i).gt.x_min).and.(x(i).lt.x_max)).and.
     1    ((y(i).gt.y_min).and.(y(i).lt.y_max))) then
                counter=counter+1
                    if (nn(i).eq.nn2d(i)) then
                       write(25,*) nn(i),nn2d(i),areanod2d(i)
                       totarea=totarea+areanod2d(i)
                    end if
                end if
             end if
          end do
          write(*,*) counter
          write(*,*) totarea

          write(22,100) HSPATM,IETO
          write(22,200) zero
          do i=1,nnod3d
             if (z(i).eq.zero) then
                write(22,*) zero
             end if
          end do

          do k=1,giorni_irrigazione
             write(22,200) time
             do i=1,nnod3d
                if (z(i).eq.zero) then
                   if (((x(i).gt.x_min).and.(x(i).lt.x_max)).and.
     1    ((y(i).gt.y_min).and.(y(i).lt.y_max))) then
                      write(22,*) ((flux)/totarea)
                   else 
                      write(22,*) zero
                   end if
                end if
             end do 
          
             t=time+5*sec_h+1
C          write(*,*) t
             write(22,*) 
             write(22,200) t                
             do i=1,nnod3d
                if (z(i).eq.zero) then
                   write(22,*) zero
                end if
             end do
             write(22,*)      
             time=time+sec_g
             t=t+sec_g
          end do 
          write(*,*) 'atmbc senza sospensione'

      end if
      
C Chiuso i file aperti

      close(21)
      close(22)
      close(23)
      close(24)
      close(25)

      end
