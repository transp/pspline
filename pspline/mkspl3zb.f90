subroutine mkspl3zb(fun,x,nx,th,nth,ph,nph,fspl,nf2,nf3, &
     ibcxmin,bcxmin,ibcxmax,bcxmax,nb1, &
     wk,inwk,ilinx,ilinth,ilinph,ier)
  use psp_precision_mod, only: fp
  !
  !  create a tricubic biperiodic spline with knots at grid points and
  !  function values from the callable function `fun' passed.
  !
  !  use tpsplinb to setup the spline coefficients
  !
  !  periodic boundary condition for theta & phi grids.
  !  boundary condition data may be specified at x(1) and x(nx) for each
  !  theta & phi:
  !     ibcxmin=0 ==> "not a knot" boundary condition, cf cubspl.f90
  !     ibcxmin=1 ==> match slope, bcxmin(ith,iph) gives df/dx at
  !         x(1),th(ith),ph(iph)
  !     ibcxmin=2 ==> match 2nd derivative, given at x(1),th(ith),ph(iph)
  !          by bcxmin(ith,iph)
  !
  !     ibcxmax,bcxmax have analogous interpretation -- at x(nx).
  !
  !============
  implicit none
  integer nth,nph,nf2,nf3,nb1,inwk,nx,iph,ith,ix
  !============
  real(fp) :: fun,zdummy
  !============
  external fun                      ! passed real function(x,th)
  real(fp) :: x(nx)                        ! x coordinate array
  real(fp) :: th(nth)                      ! th coordinate array
  real(fp) :: ph(nph)                      ! ph coordinate array
  !
  real(fp) :: fspl(8,nf2,nf3,nph)          ! function data / spline coeff array
  real(fp) :: wk(inwk)                     ! workspace for bpsplinb
  !
  integer ibcxmin                   ! boundary condition indicator
  real(fp) :: bcxmin(nb1,nph)              ! boundary condition data
  integer ibcxmax                   ! boundary condition indicator
  real(fp) :: bcxmax(nb1,nph)              ! boundary condition data
  !
  integer ilinx                     ! output =1 if x(...) evenly spaced
  integer ilinth                    ! output =1 if th(...) evenly spaced
  integer ilinph                    ! output =1 if ph(...) evenly spaced
  !
  integer ier                       ! completion code from bpspline
  !
  !----------------------------
  !
  ier=0
  if(nf2.lt.nx) then
     write(6,'('' ?mkspl3pb -- array dim error, nf2 .lt. nx'')')
     ier=1
  end if
  if(nf3.lt.nth) then
     write(6,'('' ?mkspl3pb -- array dim error, nf3 .lt. nth'')')
     ier=2
  end if
  if(nb1.lt.nth) then
     write(6,'('' ?mkspl3pb -- array dim error, nb1 .lt. nth'')')
     ier=3
  end if

  if(ier.ne.0) return
  !
  do iph=1,nph
     do ith=1,nth
        do ix=1,nx
           fspl(1,ix,ith,iph)=fun(x(ix),th(ith),ph(iph))
        end do
     end do
  end do
  !
  call mktricubw(x,nx,th,nth,ph,nph,fspl,nf2,nf3, &
       ibcxmin,bcxmin,ibcxmax,bcxmax,nb1, &
       -1,zdummy,-1,zdummy,max(nx,nth,nph), &
       -1,zdummy,-1,zdummy,max(nx,nth,nph), &
       wk,inwk,ilinx,ilinth,ilinph,ier)
  !
  return
end subroutine mkspl3zb
