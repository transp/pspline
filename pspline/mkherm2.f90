subroutine mkherm2(fun,x,nx,y,ny,fherm)
  use psp_precision_mod, only: fp
  !
  !  create a data set for Hermite interpolation, from evaluation of
  !  function and derivatives.  function of 2 indep. coordinates.
  !
  !  input:
  !
  !============
  implicit none
  integer ny,nx,ix,iy
  !============
  real(fp) :: fun,dfdx,dfdy,d2fdxdy
  !============
  external fun                      ! passed real function(x,y)
  real(fp) :: x(nx)                        ! x coordinate array
  real(fp) :: y(ny)                        ! y coordinate array
  !
  !  the passed function fun must have the interface:
  !
  !        real function <name>(x,y,dfdx,dfdy,d2fdxdy)
  !  where x,y are input, the function returns the function value,
  !  and the arguments dfdx and dfdy return as output the function
  !  derivative at the point (x,y).
  !
  !  output:
  !
  real(fp) :: fherm(0:3,nx,ny)             ! function data & derivatives
  !
  !  fherm(0,i,j) = function value f at x(i),y(j)
  !  fherm(1,i,j) = derivative df/dx at x(i),y(j)
  !  fherm(2,i,j) = derivative df/dy at x(i),y(j)
  !  fherm(3,i,j) = derivative d2f/dxdy at x(i),y(j)
  !
  !----------------------------
  !
  do ix=1,nx
     do iy=1,ny
        fherm(0,ix,iy)=fun(x(ix),y(iy),dfdx,dfdy,d2fdxdy)
        fherm(1,ix,iy)=dfdx
        fherm(2,ix,iy)=dfdy
        fherm(3,ix,iy)=d2fdxdy
     end do
  end do
  !
  return
end subroutine mkherm2
