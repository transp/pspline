subroutine dnherm2(x,nx,y,ny,fherm,nf2,ilinx,iliny,ier)
  use psp_precision_mod, only: fp
  !
  !  create a data set for Hermite interpolation, based on simple
  !  numerical differentiation using the given grid.  The grid should
  !  be "approximately" evenly spaced.
  !
  !  2d routine
  !
  !  input:
  !
  !  nf2.gt.nx expected!
  !
  !============
  implicit none
  integer iliny,ier,ilinx,ierx,iery,iy,iyp,iym,ix,ixp,ixm
  !============
  real(fp) :: zd
  !============
  integer nx,ny,nf2                 ! array dimensions
  real(fp) :: x(nx)                        ! x coordinate array
  real(fp) :: y(ny)                        ! y coordinate array
  real(fp) :: fherm(0:3,nf2,ny)            ! data/Hermite array
  !
  !  fherm(0,i,j) = function value f at x(i),y(j)  **on input**
  !
  !  fherm(1,i,j) = derivative df/dx at x(i),y(j)  **on output**
  !  fherm(2,i,j) = derivative df/dy at x(i),y(j)  **on output**
  !  fherm(3,i,j) = derivative d2f/dxdy at x(i),y(j)  **on output**
  !
  !  addl output:
  !    ilinx=1 if x axis is evenly spaced
  !    iliny=1 if y axis is evenly spaced
  !    ier=0 if no error:
  !      x, y must both be strict ascending
  !      nf2.ge.nx is required.
  !
  !----------------------------
  !
  ier=0
  !
  call splinck(x,nx,ilinx,1.0E-3_fp,ierx)
  if(ierx.ne.0) ier=2
  !
  if(ier.eq.2) then
     write(6,'('' ?dnherm2:  x axis not strict ascending'')')
  end if
  !
  call splinck(y,ny,iliny,1.0E-3_fp,iery)
  if(iery.ne.0) ier=3
  !
  if(ier.eq.3) then
     write(6,'('' ?dnherm2:  y axis not strict ascending'')')
  end if
  !
  if(nf2.lt.nx) then
     ier=4
     write(6,*) '?dnherm2:  fherm array dimension too small.'
  end if
  !
  if(ier.ne.0) return
  !
  do iy=1,ny
     !
     iyp=min(ny,iy+1)
     iym=max(1,iy-1)
     !
     do ix=1,nx
        !
        ixp=min(nx,ix+1)
        ixm=max(1,ix-1)
        !
        !  x diffs in vicinity
        !
        zd=(fherm(0,ixp,iy)-fherm(0,ixm,iy))/(x(ixp)-x(ixm))
        !
        fherm(1,ix,iy)=zd
        !
        !  y diffs in vicinity
        !
        zd=(fherm(0,ix,iyp)-fherm(0,ix,iym))/(y(iyp)-y(iym))
        !
        fherm(2,ix,iy)=zd
        !
        !  xy cross derivs (except at corners, iedge=2)
        !
        fherm(3,ix,iy)=(fherm(0,ixp,iyp)-fherm(0,ixm,iyp) &
             -fherm(0,ixp,iym)+fherm(0,ixm,iym))/ &
             ((x(ixp)-x(ixm))*(y(iyp)-y(iym)))
        !
     end do
  end do
  !
  return
end subroutine dnherm2
