subroutine dnherm1(x,nx,fherm,ilinx,ier)
  use psp_precision_mod, only: fp
  !
  !  create a data set for Hermite interpolation, based on simple
  !  numerical differentiation using the given grid.
  !
  !  1d routine
  !
  !  input:
  !
  !============
  implicit none
  integer ier,ilinx,ix,ixp,ixm
  !============
  real(fp) :: ztol,zd
  !============
  integer nx                        ! array dimensions
  real(fp) :: x(nx)                        ! x coordinate array
  real(fp) :: fherm(0:1,nx)                ! data/Hermite array
  !
  !  fherm(0,i) = function value f at x(i)       **on input**
  !
  !  fherm(1,i) = derivative df/dx at x(i)       **on output**
  !
  ! addl output:
  !  ilinx=1 if x is "evenly spaced" ier=0 if no errors
  !
  !  ** x must be strict ascending **
  !
  !----------------------------
  !
  ztol=1.0E-3_fp
  call splinck(x,nx,ilinx,ztol,ier)
  if(ier.ne.0) then
     write(6,*) '?dnherm1:  x axis not strict ascending.'
     return
  end if
  !
  do ix=1,nx
     !
     !  x div. diffs in vicinity
     !
     ixp=min(nx,ix+1)
     ixm=max(1,ix-1)
     zd=(fherm(0,ixp)-fherm(0,ixm))/(x(ixp)-x(ixm))
     !
     fherm(1,ix)=zd
  end do
  !
  return
end subroutine dnherm1
