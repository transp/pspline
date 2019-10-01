subroutine mktricubw(x,nx,y,ny,z,nz,f,nf2,nf3, &
     ibcxmin,bcxmin,ibcxmax,bcxmax,inb1x, &
     ibcymin,bcymin,ibcymax,bcymax,inb1y, &
     ibczmin,bczmin,ibczmax,bczmax,inb1z, &
     wk,nwk,ilinx,iliny,ilinz,ier)
  use precision_mod, only: fp
  !
  !  setup a tricubic spline; store coefficients in compatct form
  !  (as per suggestion of L. Zakharov, PPPL, Feb. 1999)
  !  8 coeffs per (x,y,z) grid point:
  !          f,fxx,fyy,fzz,fxxyy,fxxzz,fyyzz,fxxyyzz
  !
  !  input:
  !============
  implicit none
  integer itest,iadfp,isiz1,iadfw,inwk
  !============
  real(fp) :: hxlast,hylast,hzlast
  !============
  integer nx                        ! length of x vector
  integer ny                        ! length of y vector
  integer nz                        ! length of z vector
  real(fp) :: x(nx)                        ! x vector, strict ascending
  real(fp) :: y(ny)                        ! y vector, strict ascending
  real(fp) :: z(nz)                        ! z vector, strict ascending
  !
  integer nf2                       ! 2nd dim. of f array, nf2.ge.nx
  integer nf3                       ! 3rd dim. of f array, nf3.ge.ny
  !
  !  input/output:
  !
  real(fp) :: f(8,nf2,nf3,nz)              ! data and spline coefficients
  !
  !  on input:  f(1,i,j,k) = f(x(i),y(j),z(k))
  !  on output:  f(1,i,j,k) unchanged
  !              f(2,i,j,k) = d2f/dx2(x(i),y(j),z(k))
  !              f(3,i,j,k) = d2f/dy2(x(i),y(j),z(k))
  !              f(4,i,j,k) = d2f/dz2(x(i),y(j),z(k))
  !              f(5,i,j,k) = d4f/dx2dy2(x(i),y(j),z(k))
  !              f(6,i,j,k) = d4f/dx2dz2(x(i),y(j),z(k))
  !              f(7,i,j,k) = d4f/dy2dz2(x(i),y(j),z(k))
  !              f(8,i,j,k) = d6f/dx2dy2dz2(x(i),y(j),z(k))
  !
  !  there is a rather Hermite like interpolation formula to go with
  !  this-- see evtricub.f90.  Also the bicubic formula is given in
  !  mkbicubw.f90; the tricubic formula is precisely analogous.
  !
  !  boundary condition data
  !  inputs:
  integer inb1x                     ! 1st dim of xmin & xmax bc arrays
  integer inb1y                     ! 1st dim of ymin & ymax bc arrays
  integer inb1z                     ! 1st dim of zmin & zmax bc arrays
  !
  integer ibcxmin,ibcxmax           ! BC type flag @xmin, xmax
  integer ibcymin,ibcymax           ! BC type flag @ymin, ymax
  integer ibczmin,ibczmax           ! BC type flag @zmin, zmax
  !
  real(fp) :: bcxmin(inb1x,nz),bcxmax(inb1x,nz) ! xmin & xmax BC data, ny x nz
  real(fp) :: bcymin(inb1y,nz),bcymax(inb1y,nz) ! ymin & ymax BC data, nx x nz
  real(fp) :: bczmin(inb1z,ny),bczmax(inb1z,ny) ! zmin & zmax BC data, nx x ny
  !
  !  where BC data is not required, dummy scalars may be passed.
  !  the ibc* flags determine whether BC data isneeded.
  !
  !  BC data:  bcxmin & bcxmax:  BC vs. y,z @xmin,xmax
  !            bcymin & bcymax:  BC vs. x,z @ymin,ymax
  !            bczmin & bczmax:  BC vs. x,y @zmin,zmax
  !
  !   ibcxmin -- indicator for boundary condition at xmin=x(1):
  !    bcxmin(...) -- boundary condition data
  !     =-1 -- use periodic boundary condition
  !     =0 -- use "not a knot"
  !     =1 -- match slope, specified at x(1),y(iy),z(iz) by bcxmin(iy,iz)
  !     =2 -- match 2nd derivative, specified at x(1),y(iy),z(iz)
  !           by bcxmin(iy,iz
  !     =3 -- boundary condition is slope=0 (df/dx=0) at x(1), all y(j)
  !     =4 -- boundary condition is d2f/dx2=0 at x(1), all y(j)
  !     =5 -- match 1st derivative to 1st divided difference
  !     =6 -- match 2nd derivative to 2nd divided difference
  !     =7 -- match 3rd derivative to 3rd divided difference
  !           (for more detailed definition of BCs 5-7, see the
  !           comments of subroutine mkspline)
  !   ***NOTE bcxmin(...) referenced ONLY if ibcxmin=1 or ibcxmin=2
  !
  !   ibcxmax -- indicator for boundary condition at x(nx):
  !    bcxmax(...) -- boundary condition data
  !     (interpretation as with ibcxmin, bcxmin)
  !     NOTE:  if ibcxmin=-1 then the periodic BC applies on both sides
  !            and ibcxmax, bcxmax are ignored.
  !   inb1x -- 1st dimension of bcxmin, bcxmax: if ibcxmin or ibcxmax .gt. 0
  !            this must be .ge. ny.
  !
  !   interpretation of ibcymin,bcymin,ibcymax,bcymax,inb1y
  !     is same as with ibcxmin,...
  !
  !   interpretation of ibczmin,bczmin,ibczmax,bczmax,inb1z
  !     is same as with ibcxmin,...
  !
  !   the explicit bdy condition arrays are referenced only if the
  !     corresponding "ibc" flag values are set to 1 or 2.
  !
  !  workspace:
  integer nwk                       ! size of workspace
  real(fp) :: wk(nwk)                      ! workspace array
  !
  !  this version requires a very large workspace, nwk.ge.80*nx*ny*nz
  !  so as to be able to use tcspline to calculate the spline coefficients.
  !
  !  output:
  integer ilinx                     ! x vector equal spacing flag
  integer iliny                     ! y vector equal spacing flag
  integer ilinz                     ! z vector equal spacing flag
  !
  !   ilinx -- =1 on output if x(nx) pts are nearly evenly spaced (tol=1e-3)
  !   iliny -- =1 on output if y(ny) evenly spaced (tol=1e-3)
  !   ilinz -- =1 on output if z(nz) evenly spaced (tol=1e-3)
  !
  integer ier                       ! exit code
  !   ier -- completion code, 0 for normal
  !
  !-----------------------------------------------------
  !
  itest=80*nx*ny*nz
  if(nwk.lt.itest) then
     write(6,9901) nwk,itest
9901 format(' ?mktricubw:  workspace too small:'/ &
          '  user supplied:  nwk=',i7,'; need at least:  ',i7/ &
          '  nwk = at least 21*nx*ny is required.')
     ier=1
     return
  end if
  !
  iadfp=1
  isiz1=64*nx*ny*nz
  iadfw=iadfp+isiz1
  inwk = nwk-isiz1
  !
  call mktricop(f,nf2,nf3,wk(iadfp),nx,ny,nz)
  !
  !  evaluate 4x4x4 continuous tricubic spline
  !
  call tcspline(x,nx,y,ny,z,nz,wk(iadfp),nx,ny, &
       ibcxmin,bcxmin,ibcxmax,bcxmax,inb1x, &
       ibcymin,bcymin,ibcymax,bcymax,inb1y, &
       ibczmin,bczmin,ibczmax,bczmax,inb1z, &
       wk(iadfw),inwk,ilinx,iliny,ilinz,ier)
  !
  !  convert to 8-coefficient form
  !
  hxlast=x(nx)-x(nx-1)
  hylast=y(ny)-y(ny-1)
  hzlast=z(nz)-z(nz-1)
  call mktricon(f,nf2,nf3,wk(iadfp),nx,ny,nz,hxlast,hylast,hzlast)
  !
  return
end subroutine mktricubw
!----------------------------------------------------------------
!  mktricop -- copy spline function input data
!
subroutine mktricop(fin,nf2,nf3,fwk,nx,ny,nz)
  use precision_mod, only: fp
  !
  !============
  implicit none
  integer nf3,nx,ny,nz,nf2,iz,iy,ix
  !============
  real(fp) :: fin(8,nf2,nf3,nz)
  real(fp) :: fwk(4,4,4,nx,ny,nz)
  !
  do iz=1,nz
     do iy=1,ny
        do ix=1,nx
           fwk(1,1,1,ix,iy,iz)=fin(1,ix,iy,iz)
        end do
     end do
  end do
  !
  return
end subroutine mktricop
!----------------------------------------------------------------
!  mktricon -- create compact spline representation from 4x4
!             (bcspline) representation
!
subroutine mktricon(fin,nf2,nf3,fwk,nx,ny,nz, &
     hxlast,hylast,hzlast)
  use precision_mod, only: fp
  !
  !============
  implicit none
  integer nf3,nx,ny,nz,nf2,iz,iy,ix,iflag,ixuse,iyuse,izuse,j
  !============
  real(fp) :: hxlast,hylast,hzlast,dxuse,dyuse,dzuse
  !============
  real(fp) :: fin(8,nf2,nf3,nz)
  real(fp) :: fwk(4,4,4,nx,ny,nz)
  !-----------------------------------------------------
  !  local arrays
  !
  integer iselect(10)
  real(fp) :: zvalues(10)
  !
  iselect = 0
  iselect(1) = -1
  !
  !-----------------------------------------------------
  !
  do iz=1,nz
     do iy=1,ny
        do ix=1,nx
           !
           !  copy derivatives from result.  Special treatment needed for end zones
           !
           iflag=0
           dxuse=0.0_fp
           dyuse=0.0_fp
           dzuse=0.0_fp
           ixuse=ix
           iyuse=iy
           izuse=iz
           !
           if(ix.eq.nx) then
              iflag=1
              dxuse=hxlast
              ixuse=ix-1
           end if
           if(iy.eq.ny) then
              iflag=1
              dyuse=hylast
              iyuse=iy-1
           end if
           if(iz.eq.nz) then
              iflag=1
              dzuse=hzlast
              izuse=iz-1
           end if
           !
           if(iflag.eq.1) then
              call tcspevfn(iselect,1,1,zvalues, &
                   ixuse,iyuse,izuse,dxuse,dyuse,dzuse, &
                   fwk,nx,ny,nz)
              do j=2,8
                 fin(j,ix,iy,iz)=zvalues(j)
              end do
           else
              fin(2,ix,iy,iz)=2.0_fp*fwk(3,1,1,ix,iy,iz)
              fin(3,ix,iy,iz)=2.0_fp*fwk(1,3,1,ix,iy,iz)
              fin(4,ix,iy,iz)=2.0_fp*fwk(1,1,3,ix,iy,iz)
              fin(5,ix,iy,iz)=4.0_fp*fwk(3,3,1,ix,iy,iz)
              fin(6,ix,iy,iz)=4.0_fp*fwk(3,1,3,ix,iy,iz)
              fin(7,ix,iy,iz)=4.0_fp*fwk(1,3,3,ix,iy,iz)
              fin(8,ix,iy,iz)=8.0_fp*fwk(3,3,3,ix,iy,iz)
           end if
           !
        end do
     end do
  end do
  !
  return
end subroutine mktricon
