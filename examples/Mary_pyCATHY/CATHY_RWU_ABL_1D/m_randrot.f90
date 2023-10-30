module m_randrot
contains
subroutine randrot(Q,nrens)
   implicit none
   integer, intent(in)  :: nrens
   real*8,  intent(out) :: Q(nrens,nrens)

   real*8, dimension(nrens,nrens) ::  A, B
   real*8  sigma(nrens), work(10*nrens)
   real*8, parameter :: pi=3.14159253589d0
   integer ierr
   real*8 meanB

   call random_number(B)
   call random_number(A)
   Q = dsqrt(-2.0d0*dlog(A+tiny(A))) * dcos(2.0d0*pi*B)

!$OMP CRITICAL
! QR factorization
   call dgeqrf(nrens, nrens, Q, nrens, sigma, work, 10*nrens, ierr )
   if (ierr /= 0) print *, 'randrot: dgeqrf ierr=',ierr

! Construction of Q
   call dorgqr(nrens, nrens, nrens, Q, nrens, sigma, work, 10*nrens, ierr )
   if (ierr /= 0) print *, 'randrot: dorgqr ierr=',ierr
!$OMP END CRITICAL


end subroutine randrot
end module m_randrot

