      program meteo_input

      implicit none
      integer n
      parameter(n=20000)
      integer i 
      integer a,nt,nnod,nnod3,nel
      real*8  time(n),rh(n),par(n),temp(n),zen(n),atmbc(n)
      real*8  x(n),y(n),z(n),WTdepth,lai(n)
      character*13 date(n)

      nt=17556
      WTdepth=2.0

C     Open Input file
      open(21,file='time.txt',status='old')
      open(22,file='Temperatura.txt',status='old')
      open(23,file='RH.txt')
      open(24,file='PAR.txt')
      open(25,file='ZEN.txt')
c      open(26,file='RainEvap')
      open(27,file='LAI.txt')
      open(28,file='grid3d')

C     Open Output files
      open(31,file='plant_meteo')
      open(32,file='atmbc')
      open(33,file='nansfdirbc')

C     Read Input file
      do i=1,nt
         read(21,*) time(i)
         read(22,*) temp(i)
         read(23,*) rh(i)
         read(24,*) par(i)
         read(25,*) zen(i)
c         read(26,*) atmbc(i)
         read(27,*) lai(i)
      enddo

c     Read grid3d
      read(28,*) nnod,nnod3,nel
      do i=1,nel
         read(28,*) 
      enddo
      a=0
      do i=1,nnod3
         read(28,*) x(i),y(i),z(i)
         if ((x(i).eq.0).or.(x(i).eq.5).or.(y(i).eq.0).or.
     1       (y(i).eq.5))then
         a=a+1
         endif
      enddo

C     Write plant_meteo file
      do i=1,nt
         write(31,*) time(i), 'TIME'
         write(31,*)temp(i),rh(i),par(i),zen(i),lai(i)
      enddo

C     Write atmbc file
      write(32,*) '1  1  HSPATM IETO'
         write(32,*) 0.D0, 'TIME'
         write(32,*) 0.D0
         write(32,*) 1e+20, 'TIME'
         write(32,*) 0.D0
c      do i=1,nt
c         write(32,*) time(i), 'TIME'
c         write(32,*) atmbc(i)
c      enddo

C     Write dirbc 
      write(33,*) 0.0, 'TIME'
      write(33,*) '0', a
      do i=1,nnod3
         if ((x(i).eq.0).or.(x(i).eq.5).or.(y(i).eq.0).or.
     1       (y(i).eq.5))then
         write(33,*) i
         endif
      enddo
      do i=1,nnod3
         if ((x(i).eq.0).or.(x(i).eq.5).or.(y(i).eq.0).or.
     1       (y(i).eq.5))then
         write(33,*) -z(i)-WTdepth
         endif
      enddo
         
      write(33,*) 2e+20, 'TIME'
      write(33,*) '0', a
      do i=1,nnod3
         if ((x(i).eq.0).or.(x(i).eq.5).or.(y(i).eq.0).or.
     1       (y(i).eq.5))then
         write(33,*) i
         endif
      enddo
      do i=1,nnod3
         if ((x(i).eq.0).or.(x(i).eq.5).or.(y(i).eq.0).or.
     1       (y(i).eq.5))then
         write(33,*) -z(i)-WTdepth
         endif
      enddo
      




      stop
      end
