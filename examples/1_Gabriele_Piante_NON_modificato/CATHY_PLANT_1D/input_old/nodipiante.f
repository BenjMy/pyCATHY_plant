      program nodipiante
      implicit none
      integer i

      open(21,file='nodipiante')
      do i=1,14471
         write(21,*) i, 1
      enddo

      stop
      end
