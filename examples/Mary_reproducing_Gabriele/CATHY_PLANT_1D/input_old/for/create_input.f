      program input
      implicit none
      integer   i,j,k,n,nmax,ndir
      parameter (nmax=2000000)
      integer   nodidir(nmax),nsup,nnod,nel
      real*8    time(nmax),temp(nmax),rh(nmax),rad(nmax),zen(nmax)
      real*8    lai(nmax),rain(nmax),wtlevel(nmax),x(nmax),y(nmax)
      real*8    z(nmax),ic1(7),ic2(7),rn(nmax),sw(nmax),wind(nmax)
      real*8    g(nmax),hs(nmax),cloud(nmax),patm(nmax),ustar(nmax)

      open(14,file='grid3d')
      open(23,file='nansfdirbc')      

C     Read INPUT data
      read(14,*) nsup,nnod,nel
      do i=1,nel
         read(14,*)
      enddo
      do i=1,nnod
         read(14,*) x(i),y(i),z(i)
      enddo

C Se ho Dirichlet sotto la falda e no flux sopra uso questa ***


      write(23,*) 0.00, 'TIME'
      ndir= 0
      do j=1,nnod
         if((x(j).eq.0).or.(x(j).eq.5).or.(y(j).eq.0).or.(y(j).eq.5))
     1   then
             if (z(j).lt.-1.5)then
               ndir =ndir+1
               nodidir(ndir)=j
             endif
          endif
       enddo

      write(23,*) '0', ndir
      do j=1,ndir
            write(23,*) nodidir(j)
      enddo
      do j=1,ndir
             write(23,*) -z(nodidir(j))-1.5
      enddo

      write(23,*) 1e+20, 'TIME'
      write(23,*) '0', ndir
      do j=1,ndir
            write(23,*) nodidir(j)
      enddo
      do j=1,ndir
             write(23,*) -z(nodidir(j))-1.5
      enddo
      
     




 100  format (e16.9,2x,e16.9,2x,e16.9,2x,e16.9,2x,e16.9,2x,e16.9,2x,
     1        e16.9,2x,e16.9,2x,e16.9,2x,e16.9,2x,e16.9,2x,e16.9)

      stop
      end
      
 
