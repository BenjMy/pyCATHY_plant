module m_multa
contains
subroutine multa(A, X, ndim, nrens, iblkmax)
implicit none
integer, intent(in) :: ndim
integer, intent(in) :: nrens
integer, intent(in) :: iblkmax
real*8, intent(in)    :: X(nrens,nrens)
real*8, intent(inout) :: A(ndim,nrens)
real*8 v(iblkmax,nrens)  ! Automatic work array

integer ia,ib
do ia = 1,ndim,iblkmax
  ib = min(ia+iblkmax-1,ndim)
  v(1:ib-ia+1,1:nrens) = A(ia:ib,1:nrens)
  call dgemm('n','n', ib-ia+1, nrens, nrens, &
              1.0d0, v(1,1), iblkmax, &
              X(1,1), nrens, &
              0.0d0, A(ia,1), ndim)
enddo
end subroutine multa
end module m_multa
