!--------------------------------------------------------------------
!  SPLBRK -- make a spline with a break (C0 only) at specified locn
!
subroutine splbrk (IOPT, N, NBRK, X, Y, B, C, D)
  use precision_mod, only: fp
  !
  !  Spline the whole interval, and then respline a subinterval, saving
  !  the endpt coeffs from the original spline.
  !
  !  Result is a spline that is C2 everywhere except C0 only at one
  !  interior break point.
  !
  !  IOPT=0  use SPLAAN for dy/dx=0 at LHS bc
  !  IOPT=1  use SPLINE for standard bc
  !
  !  N -- no. of data pts
  !  NBRK -- break point
  !
  !============
  implicit none
  integer n,nbrk,iopt
  !============
  real(fp) :: bsave,csave,dsave
  !============
  real(fp) :: X(N),Y(N),B(N),C(N),D(N)
  !
  if(iopt.eq.0) then
     call splaan(n, x, y, b, c, d)
  else
     call spline(n, x, y, b, c, d)
  end if
  !
  bsave=b(nbrk)
  csave=c(nbrk)
  dsave=d(nbrk)
  !
  call nspline(nbrk, x, y, b, c, d)
  !
  b(nbrk)=bsave
  c(nbrk)=csave
  d(nbrk)=dsave
  !
  return
end subroutine splbrk
