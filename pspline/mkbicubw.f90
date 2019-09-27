subroutine mkbicubw(x,nx,y,ny,f,nf2, &
     ibcxmin,bcxmin,ibcxmax,bcxmax, &
     ibcymin,bcymin,ibcymax,bcymax, &
     wk,nwk, &
     ilinx,iliny,ier)
  use precision_mod, only: fp
  !
  !  setup a bicubic spline; store coefficients in compact form
  !  (as per suggestion of L. Zakharov, PPPL, Feb. 1999)
  !     4 coeffs per grid point:  f,fxx,fyy,fxxyy
  !
  !
  !  input:
  !============
  implicit none
  integer itest,iadfp,isiz1,iadfw,inwk
  !============
  real(fp) :: hxlast,hylast
  !============
  integer nx                        ! length of x vector
  integer ny                        ! length of y vector
  real(fp) :: x(nx)                        ! x vector, strict ascending
  real(fp) :: y(ny)                        ! y vector, strict ascending
  !
  integer nf2                       ! 2nd dimension of f, nf2.ge.nx
  !  input/output:
  !
  real(fp) :: f(4,nf2,ny)                  ! data & spline coefficients
  !
  !  on input:  f(1,i,j) = f(x(i),y(j))
  !  on output:  f(1,i,j) unchanged
  !              f(2,i,j) = d2f/dx2(x(i),y(j))
  !              f(3,i,j) = d2f/dy2(x(i),y(j))
  !              f(4,i,j) = d4f/dx2dy2(x(i),y(j))
  !
  !  and the interpolation formula for (x,y) in (x(i),x(i+1))^(y(j),y(j+1))
  !  is:
  !        hx = x(i+1)-x(i)   hy = y(j+1)-y(j)
  !        dxp= (x-x(i))/hx   dxm= 1-dxp     dxp,dxm in (0,1)
  !        dyp= (x-x(i))/hx   dym= 1-dyp     dyp,dym in (0,1)
  !        dx3p = dxp**3-dxp  dx3m = dxm**3-dxm     dxp3,dxm3 in (0,1)
  !
  !   finterp = dxm*(dym*f(1,i,j)+dyp*f(1,i,j+1))
  !            +dxp*(dym*f(1,i+1,j)+dyp*f(1,i+1,j+1))
  !     +1/6*hx**2*
  !            dx3m*(dym*f(2,i,j)+dyp*f(2,i,j+1))
  !           +dx3p*(dym*f(2,i+1,j)+dyp*f(2,i+1,j+1))
  !     +1/6*hy**2*
  !            dxm*(dy3m*f(3,i,j)+dy3p*f(3,i,j+1))
  !           +dxp*(dy3m*f(3,i+1,j)+dy3p*f(3,i+1,j+1))
  !     +1/36*hx**2*hy**2*
  !            dx3m*(dym*f(4,i,j)+dyp*f(4,i,j+1))
  !           +dx3p*(dym*f(4,i+1,j)+dyp*f(4,i+1,j+1))
  !
  !  where the f(2:4,*,*) are cleverly evaluated to assure
  !  (a) finterp is continuous and twice differentiable across all
  !      grid cell boundaries, and
  !  (b) all boundary conditions are satisfied.
  !
  !  the boundary conditions specification options are:
  !  inputs:
  !
  integer ibcxmin                   ! bc flag for x=xmin
  real(fp) :: bcxmin(ny)                   ! bc data vs. y at x=xmin
  integer ibcxmax                   ! bc flag for x=xmax
  real(fp) :: bcxmax(ny)                   ! bc data vs. y at x=xmax
  !
  integer ibcymin                   ! bc flag for y=ymin
  real(fp) :: bcymin(nx)                   ! bc data vs. x at y=ymin
  integer ibcymax                   ! bc flag for y=ymax
  real(fp) :: bcymax(nx)                   ! bc data vs. x at y=ymax
  !
  !  with interpretation:
  !   ibcxmin -- indicator for boundary condition at x(1):
  !    bcxmin(...) -- boundary condition data
  !     =-1 -- periodic boundary condition
  !     =0 -- use "not a knot"
  !     =1 -- match slope, specified at x(1),th(ith) by bcxmin(ith)
  !     =2 -- match 2nd derivative, specified at x(1),th(ith) by bcxmin(ith)
  !     =3 -- boundary condition is slope=0 (df/dx=0) at x(1), all th(j)
  !     =4 -- boundary condition is d2f/dx2=0 at x(1), all th(j)
  !     =5 -- match 1st derivative to 1st divided difference
  !     =6 -- match 2nd derivative to 2nd divided difference
  !     =7 -- match 3rd derivative to 3rd divided difference
  !           (for more detailed definition of BCs 5-7, see the
  !           comments of subroutine mkspline)
  !   NOTE bcxmin(...) referenced ONLY if ibcxmin=1 or ibcxmin=2
  !
  !   ibcxmax -- indicator for boundary condition at x(nx):
  !    bcxmax(...) -- boundary condition data
  !     (interpretation as with ibcxmin, bcxmin)
  !   NOTE:  if ibcxmin=-1, ibcxmax is ignored! ...and the BC is periodic.
  !
  !  and analogous interpretation for ibcymin,bcymin,ibcymax,bcymax
  !  (df/dy or d2f/dy2 boundary conditions at y=ymin and y=ymax).
  !
  !  *** workspace ***
  !  This FORTRAN-77 routine requires a workspace of
  !    nwk = (AT LEAST) 21*nx*ny
  !  *** for dynamic allocation of this workspace use subroutine mkbicub
  !    which has the same arguments as mkbicubw, except that the workspace
  !    does not have to be provided -- mkbicub is FORTRAN-90.
  !
  !  input:
  integer nwk                       ! size of workspace
  !  modified on output:
  real(fp) :: wk(nwk)                      ! workspace
  !
  !  and output arguments
  !
  integer ilinx                     ! =1: x grid is "nearly" equally spaced
  integer iliny                     ! =1: y grid is "nearly" equally spaced
  !
  !  ilinx and iliny are set to zero if corresponding grids are not equally
  !  spaced
  !
  !  and completion code
  !
  integer ier                       ! =0 on exit if there is no error.
  !
  !  if there is an error, ier is set and a message is output on unit 6.
  !  these are considered programming errors in the calling routine.
  !
  !  possible errors:
  !    x(...) not strict ascending
  !    y(...) not strict ascending
  !    nx .lt. 4
  !    ny .lt. 4
  !    invalid boundary condition flag
  !    workspace too small
  !
  !-----------------------------------------------------
  !  local arrays
  !
  integer iselect(10)
  integer zvalues(10)
  !
  data iselect/-1,0,0,0,0,0,0,0,0,0/
  !
  !-----------------------------------------------------
  !
  itest=21*nx*ny
  if(nwk.lt.itest) then
     write(6,9901) nwk,itest
9901 format(' ?mkbicubw:  workspace too small:'/ &
          '  user supplied:  nwk=',i7,'; need at least:  ',i7/ &
          '  nwk = at least 21*nx*ny is required.')
     ier=1
     return
  end if
  !
  iadfp=1
  isiz1=16*nx*ny
  iadfw=iadfp+isiz1
  inwk = nwk-isiz1
  !
  call mkbicop(f,nf2,wk(iadfp),nx,ny)
  !
  !  evaluate 4x4 continuous bicubic spline
  !
  call bcspline(x,nx,y,ny,wk(iadfp),nx, &
       ibcxmin,bcxmin,ibcxmax,bcxmax, &
       ibcymin,bcymin,ibcymax,bcymax, &
       wk(iadfw),inwk, &
       ilinx,iliny,ier)
  !
  if(ier.ne.0) return
  !
  !  convert to 4-coefficient (2nd/4th partials form)
  !
  hxlast=x(nx)-x(nx-1)
  hylast=y(ny)-y(ny-1)
  call mkbicon(f,nf2,wk(iadfp),nx,ny,hxlast,hylast)
  !
  return
end subroutine mkbicubw
!
!----------------------------------------------------------------
!  mkbicop -- copy spline function input data
!
subroutine mkbicop(fin,nf2,fwk,nx,ny)
  use precision_mod, only: fp
  !
  !============
  implicit none
  integer nx,ny,nf2,iy,ix
  !============
  real(fp) :: fin(4,nf2,ny)
  real(fp) :: fwk(4,4,nx,ny)
  !
  do iy=1,ny
     do ix=1,nx
        fwk(1,1,ix,iy)=fin(1,ix,iy)
     end do
  end do
  !
  return
end subroutine mkbicop
!----------------------------------------------------------------
!  mkbicon -- create compact spline representation from 4x4
!             (bcspline) representation
!
subroutine mkbicon(fin,nf2,fwk,nx,ny,hxlast,hylast)
  use precision_mod, only: fp
  !
  !============
  implicit none
  integer nx,ny,nf2,iy,ix,iflag,ixuse,iyuse,j
  !============
  real(fp) :: hxlast,hylast,dxuse,dyuse
  !============
  real(fp) :: fin(4,nf2,ny)
  real(fp) :: fwk(4,4,nx,ny)
  !
  !-----------------------------------------------------
  !  local arrays
  !
  integer iselect(10)
  real(fp) :: zvalues(10)
  !
  data iselect/-1,0,0,0,0,0,0,0,0,0/
  !
  !-----------------------------------------------------
  !
  do iy=1,ny
     do ix=1,nx
        !
        !  copy derivatives from result.  Special treatment needed for end zones
        !
        iflag=0
        dxuse=0.0_fp
        dyuse=0.0_fp
        ixuse=ix
        iyuse=iy
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
        !
        if(iflag.eq.1) then
           call bcspevfn(iselect,1,1,zvalues, &
                ixuse,iyuse,dxuse,dyuse, &
                fwk,nx,ny)
           do j=2,4
              fin(j,ix,iy)=zvalues(j)
           end do
        else
           fin(2,ix,iy)=2.0_fp*fwk(3,1,ix,iy)
           fin(3,ix,iy)=2.0_fp*fwk(1,3,ix,iy)
           fin(4,ix,iy)=4.0_fp*fwk(3,3,ix,iy)
        end if
        !
     end do
  end do
  !
  return
end subroutine mkbicon
