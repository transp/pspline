subroutine dnherm3(x,nx,y,ny,z,nz,fherm,nf2,nf3, &
     ilinx,iliny,ilinz,ier)
  use precision_mod, only: fp
  !
  !  create a data set for Hermite interpolation, based on simple
  !  numerical differentiation using the given grid.  The grid should
  !  be "approximately" evenly spaced.
  !
  !  3d routine
  !
  !  input:
  !
  !  nf2.get.nx and nf3.ge.ny  expected!
  !
  !============
  implicit none
  integer iliny,ilinz,ier,ilinx,ierx,iery,ierz,iz,izp,izm,iy
  integer iyp,iym,ix,ixp,ixm
  !============
  real(fp) :: zd
  !============
  integer nx,ny,nz,nf2,nf3          ! array dimensions
  real(fp) :: x(nx)                        ! x coordinate array
  real(fp) :: y(ny)                        ! y coordinate array
  real(fp) :: z(nz)                        ! z coordinate array
  real(fp) :: fherm(0:7,nf2,nf3,nz)        ! data/Hermite array
  !
  !  fherm(0,i,j,k) = function value f at x(i),y(j),z(k)  **on input**
  !
  !  fherm(1,i,j,k) = derivative df/dx at x(i),y(j),z(k)  **on output**
  !  fherm(2,i,j,k) = derivative df/dy at x(i),y(j),z(k)  **on output**
  !  fherm(3,i,j,k) = derivative df/dz at x(i),y(j),z(k)  **on output**
  !  fherm(4,i,j,k) = derivative d2f/dxdy at x(i),y(j),z(k)  **on output**
  !  fherm(5,i,j,k) = derivative d2f/dxdz at x(i),y(j),z(k)  **on output**
  !  fherm(6,i,j,k) = derivative d2f/dydz at x(i),y(j),z(k)  **on output**
  !  fherm(7,i,j,k) = derivative d3f/dxdydz at x(i),y(j),z(k)  **on output**
  !
  !  additional outputs:
  !
  !    ilinx = 1  if x is evenly spaced (approximately)
  !    iliny = 1  if y is evenly spaced (approximately)
  !    ilinz = 1  if z is evenly spaced (approximately)
  !
  !    ier = 0 if there are no errors
  !
  !  note possible errors:  x y or z NOT strict ascending
  !  nf2.lt.nx   .or.   nf3.lt.ny
  !
  !----------------------------
  !
  !
  !  error check
  !
  ier=0
  !
  call splinck(x,nx,ilinx,1.0E-3_fp,ierx)
  if(ierx.ne.0) ier=2
  !
  if(ier.eq.2) then
     write(6,'('' ?dnherm3:  x axis not strict ascending'')')
  end if
  !
  call splinck(y,ny,iliny,1.0E-3_fp,iery)
  if(iery.ne.0) ier=3
  !
  if(ier.eq.3) then
     write(6,'('' ?dnherm3:  y axis not strict ascending'')')
  end if
  !
  call splinck(z,nz,ilinz,1.0E-3_fp,ierz)
  if(ierz.ne.0) ier=4
  !
  if(ier.eq.4) then
     write(6,'('' ?dnherm3:  z axis not strict ascending'')')
  end if
  !
  if(nf2.lt.nx) then
     ier=5
     write(6,*) '?dnherm3:  fherm (x) array dimension too small.'
  end if
  !
  if(nf3.lt.ny) then
     ier=6
     write(6,*) '?dnherm3:  fherm (y) array dimension too small.'
  end if
  !
  if(ier.ne.0) return
  !
  !  error check OK
  !
  do iz=1,nz
     !
     izp=min(nz,iz+1)
     izm=max(1,iz-1)
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
           zd=(fherm(0,ixp,iy,iz)-fherm(0,ixm,iy,iz))/ &
                (x(ixp)-x(ixm))
           !
           fherm(1,ix,iy,iz)=zd
           !
           !  y diffs in vicinity
           !
           zd=(fherm(0,ix,iyp,iz)-fherm(0,ix,iym,iz))/ &
                (y(iyp)-y(iym))
           !
           fherm(2,ix,iy,iz)=zd
           !
           !  z diffs in vicinity
           !
           zd=(fherm(0,ix,iy,izp)-fherm(0,ix,iy,izm))/ &
                (z(izp)-z(izm))
           !
           fherm(3,ix,iy,iz)=zd
           !
           !  xy cross derivs
           !
           fherm(4,ix,iy,iz)= &
                (fherm(0,ixp,iyp,iz)-fherm(0,ixm,iyp,iz) &
                -fherm(0,ixp,iym,iz)+fherm(0,ixm,iym,iz))/ &
                ((x(ixp)-x(ixm))*(y(iyp)-y(iym)))
           !
           !  xz cross derivs
           !
           fherm(5,ix,iy,iz)= &
                (fherm(0,ixp,iy,izp)-fherm(0,ixm,iy,izp) &
                -fherm(0,ixp,iy,izm)+fherm(0,ixm,iy,izm))/ &
                ((x(ixp)-x(ixm))*(z(izp)-z(izm)))
           !
           !  yz cross derivs
           !
           fherm(6,ix,iy,iz)= &
                (fherm(0,ix,iyp,izp)-fherm(0,ix,iym,izp) &
                -fherm(0,ix,iyp,izm)+fherm(0,ix,iym,izm))/ &
                ((y(iyp)-y(iym))*(z(izp)-z(izm)))
           !
           !  xyz cross deriv
           !
           fherm(7,ix,iy,iz)= &
                ((fherm(0,ixp,iyp,izp)-fherm(0,ixp,iym,izp) &
                -fherm(0,ixp,iyp,izm)+fherm(0,ixp,iym,izm))- &
                (fherm(0,ixm,iyp,izp)-fherm(0,ixm,iym,izp) &
                -fherm(0,ixm,iyp,izm)+fherm(0,ixm,iym,izm)))/ &
                ((x(ixp)-x(ixm))*(y(iyp)-y(iym))*(z(izp)-z(izm)))
           !
        end do
     end do
  end do
  !
  return
end subroutine dnherm3
