subroutine mkspl2pb(fun,x,nx,th,nth,fspl,nf3, &
     ibcxmin,bcxmin,ibcxmax,bcxmax, &
     wk,inwk,ilinx,ilinth,ier)
  use psp_precision_mod, only: fp
  !
  !  create a bicubic periodic spline with knots at grid points and
  !  function values from the callable function `fun' passed.
  !
  !  use bpsplinb to setup the spline coefficients
  !
  !  periodic boundary condition for theta grid;
  !  boundary condition data may be specified at x(1) and x(nx)
  !     ibcxmin=0 ==> "not a knot" boundary condition, cf cubspl.f90
  !     ibcxmin=1 ==> match slope, bcxmin(ith) gives df/dx at x(1),th(ith).
  !     ibcxmin=2 ==> match 2nd derivative, given at x(1),th(ith) by bcxmin(ith)
  !
  !     ibcxmax,bcxmax have analogous interpretation -- at x(nx)
  !
  !============
  implicit none
  integer nth,nf3,inwk,nx,ix,ith
  !============
  real(fp) :: fun
  !============
  external fun                      ! passed real function(x,th)
  real(fp) :: x(nx)                        ! x coordinate array
  real(fp) :: th(nth)                      ! th coordinate array
  !
  real(fp) :: fspl(4,4,nf3,nth)            ! function data / spline coeff array
  real(fp) :: wk(inwk)                     ! workspace for bpsplinb
  !
  integer ibcxmin                   ! boundary condition indicator
  real(fp) :: bcxmin(nth)                  ! boundary condition data
  integer ibcxmax                   ! boundary condition indicator
  real(fp) :: bcxmax(nth)                  ! boundary condition data
  !
  integer ilinx                     ! output =1 if x(...) evenly spaced
  integer ilinth                    ! output =1 if th(...) evenly spaced
  !
  integer ier                       ! completion code from bpspline
  !
  !----------------------------
  !
  ier=0
  if(nf3.lt.nx) then
     write(6,'('' ?mkspl2pb -- array dim error, nf3.lt.nx'')')
     ier=1
  end if
  if(inwk.lt.5*max(nx,nth)) then
     write(6,'('' ?mkspl2pb -- array dim error, inwk too small'')')
     ier=2
  end if
  !
  do ix=1,nx
     do ith=1,nth
        fspl(1,1,ix,ith)=fun(x(ix),th(ith))
     end do
  end do
  !
  call bpsplinb(x,nx,th,nth,fspl,nx,ibcxmin,bcxmin,ibcxmax,bcxmax, &
       wk,inwk,ilinx,ilinth,ier)
  !
  return
end subroutine mkspl2pb
