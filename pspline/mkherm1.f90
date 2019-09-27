subroutine mkherm1(fun,x,nx,fherm)
  use precision_mod, only: fp
  !
  !  create a data set for Hermite interpolation, from evaluation of
  !  function and derivatives.  function of 1 indep. coordinate.
  !
  !  input:
  !
  !============
  implicit none
  integer nx,ix
  !============
  real(fp) :: fun,dfdx
  !============
  external fun                      ! passed real function(x,dfdx)
  real(fp) :: x(nx)                        ! x coordinate array
  !
  !  the passed function fun must have the interface:
  !
  !        real function <name>(x,dfdx)
  !  where x is input, the function returns the function value and the
  !  argument dfdx returns as output the function derivative at x.
  !
  !  output:
  !
  real(fp) :: fherm(0:1,nx)                ! function data & derivative
  !
  !  fherm(0,j) = function value f at x(j)
  !  fherm(1,j) = derivative df/dx at x(j)
  !
  !----------------------------
  !
  do ix=1,nx
     fherm(0,ix)=fun(x(ix),dfdx)
     fherm(1,ix)=dfdx
  end do
  !
  return
end subroutine mkherm1
