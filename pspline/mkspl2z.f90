subroutine mkspl2z(fun,x,nx,th,nth,fspl,nf2,wk,inwk, &
     ilinx,ilinth,ier)
  use precision_mod, only: fp
  !
  !  create a bicubic periodic spline with knots at grid points and
  !  function values from the callable function `fun' passed.
  !
  !  use bpspline to setup the spline coefficients
  !
  !============
  implicit none
  integer nth,nf2,inwk,nx,ix,ith
  !============
  real(fp) :: fun,zdummy
  !============
  external fun                      ! passed real function(x,th)
  real(fp) :: x(nx)                        ! x coordinate array
  real(fp) :: th(nth)                      ! th coordinate array
  !
  !  output:
  !
  real(fp) :: fspl(4,nf2,nth)              ! function data / spline coeff array
  real(fp) :: wk(inwk)                     ! workspace -- at least 5*max(nx,nth)
  !
  integer ilinx                     ! output =1 if x(1...nx) evenly spaced
  integer ilinth                    ! output =1 if th(1..nth) evenly spaced
  !
  integer ier                       ! completion code from bpspline 0=OK
  !
  !----------------------------
  !
  ier=0
  if(nf2.lt.nx) then
     write(6,'('' ?mkspl2p -- array dim error, nf2.lt.nx'')')
     ier=1
  end if
  !
  do ix=1,nx
     do ith=1,nth
        fspl(1,ix,ith)=fun(x(ix),th(ith))
     end do
  end do
  !
  call mkbicubw(x,nx,th,nth,fspl,nf2, &
       0,zdummy,0,zdummy,-1,zdummy,-1,zdummy, &
       wk,inwk,ilinx,ilinth,ier)
  !
  return
end subroutine mkspl2z
