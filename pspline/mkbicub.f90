subroutine mkbicub(x,nx,y,ny,f,nf2, &
     ibcxmin,bcxmin,ibcxmax,bcxmax, &
     ibcymin,bcymin,ibcymax,bcymax, &
     ilinx,iliny,ier)
  use precision_mod, only: fp
  !
  !  setup bicubic spline, dynamic allocation of workspace
  !  fortran-90 fixed form source
  !
  !  ** see mkbicubw.f90 **
  !  arguments are identical, except the workspace argument is
  !  omitted (created here)
  !
  !  --NOTE-- dmc 22 Feb 2004 -- rewrite for direct calculation of
  !  coefficients, to avoid large transient use of memory.  mkbicubw is
  !  no longer called (but reference to mkbicubw comments is still correct).
  !
  !  abbreviated description of arguments, mkbicubw has details:
  !
  !
  !  input:
  implicit none
  integer nx                        ! length of x vector
  integer ny                        ! length of y vector
  real(fp) :: x(nx)                        ! x vector, strict ascending
  real(fp) :: y(ny)                        ! y vector, strict ascending
  !
  integer nf2                       ! 2nd dimension of f, nf2.ge.nx
  !  input/output:
  real(fp) :: f(4,nf2,ny)                  ! data & spline coefficients
  !
  !  input bdy condition data:
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
  !  output linear grid flags and completion code (ier=0 is normal):
  !
  integer ilinx                     ! =1: x grid is "nearly" equally spaced
  integer iliny                     ! =1: y grid is "nearly" equally spaced
  !
  integer ier                       ! =0 on exit if there is no error.
  !
  !-----------------------
  integer ierx,iery
  !
  real(fp), dimension(:,:), allocatable :: fwk
  real(fp) :: zbcmin,zbcmax
  integer ix,iy,ibcmin,ibcmax
  !
  real(fp), dimension(:,:,:), allocatable :: fcorr
  integer iflg2
  real(fp) :: zdiff(2),hy
  !
  !-----------------------
  !
  !  see if 2nd pass is needed due to inhomogeneous d/dy bdy cond.
  !
  iflg2=0
  if(ibcymin.ne.-1) then
     if((ibcymin.eq.1).or.(ibcymin.eq.2)) then
        do ix=1,nx
           if (bcymin(ix).ne.0.0_fp) iflg2=1
        end do
     end if
     if((ibcymax.eq.1).or.(ibcymax.eq.2)) then
        do ix=1,nx
           if (bcymax(ix).ne.0.0_fp) iflg2=1
        end do
     end if
  end if
  !
  !  check boundary condition specifications
  !
  ier=0
  !
  call ibc_ck(ibcxmin,'bcspline','xmin',-1,7,ier)
  if(ibcxmin.ge.0) call ibc_ck(ibcxmax,'bcspline','xmax',0,7,ier)
  call ibc_ck(ibcymin,'bcspline','ymin',-1,7,ier)
  if(ibcymin.ge.0) call ibc_ck(ibcymax,'bcspline','ymax',0,7,ier)
  !
  !  check ilinx & x vector
  !
  call splinck(x,nx,ilinx,1.0E-3_fp,ierx)
  if(ierx.ne.0) ier=2
  !
  if(ier.eq.2) then
     write(6,'('' ?bcspline:  x axis not strict ascending'')')
  end if
  !
  !  check iliny & y vector
  !
  call splinck(y,ny,iliny,1.0E-3_fp,iery)
  if(iery.ne.0) ier=3
  !
  if(ier.eq.3) then
     write(6,'('' ?bcspline:  y axis not strict ascending'')')
  end if
  !
  if(ier.ne.0) return
  !
  !------------------------------------
  allocate(fwk(2,max(nx,ny)))
  !
  !  evaluate fxx (spline in x direction)
  !
  zbcmin=0
  zbcmax=0
  do iy=1,ny
     fwk(1,1:nx) = f(1,1:nx,iy)
     if((ibcxmin.eq.1).or.(ibcxmin.eq.2)) zbcmin=bcxmin(iy)
     if((ibcxmax.eq.1).or.(ibcxmax.eq.2)) zbcmax=bcxmax(iy)
     call mkspline(x,nx,fwk, &
          ibcxmin,zbcmin,ibcxmax,zbcmax,ilinx,ier)
     if(ier.ne.0) return
     f(2,1:nx,iy)=fwk(2,1:nx)
  end do
  !
  !  evaluate fyy (spline in y direction)
  !  use homogeneous boundary condition; correction done later if necessary
  !
  zbcmin=0
  zbcmax=0
  ibcmin=ibcymin
  ibcmax=ibcymax
  do ix=1,nx
     fwk(1,1:ny) = f(1,ix,1:ny)
     if(iflg2.eq.1) then
        if((ibcymin.eq.1).or.(ibcymin.eq.2)) ibcmin=0
        if((ibcymax.eq.1).or.(ibcymax.eq.2)) ibcmax=0
     end if
     call mkspline(y,ny,fwk, &
          ibcmin,zbcmin,ibcmax,zbcmax,iliny,ier)
     if(ier.ne.0) return
     f(3,ix,1:ny)=fwk(2,1:ny)
  end do
  !
  !  evaluate fxxyy (spline fxx in y direction; BC simplified; avg
  !  d2(d2f/dx2)/dy2 and d2(df2/dy2)/dx2
  !
  zbcmin=0
  zbcmax=0
  ibcmin=ibcymin
  ibcmax=ibcymax
  do ix=1,nx
     fwk(1,1:ny) = f(2,ix,1:ny)
     if(iflg2.eq.1) then
        if((ibcymin.eq.1).or.(ibcymin.eq.2)) ibcmin=0
        if((ibcymax.eq.1).or.(ibcymax.eq.2)) ibcmax=0
     end if
     call mkspline(y,ny,fwk, &
          ibcmin,zbcmin,ibcmax,zbcmax,iliny,ier)
     if(ier.ne.0) return
     f(4,ix,1:ny)= fwk(2,1:ny)
  end do
  !
  if(iflg2.eq.1) then
     allocate(fcorr(2,nx,ny))
     !
     !  correct for inhomogeneous y boundary condition
     !
     do ix=1,nx
        !  the desired inhomogenous BC is the difference btw the
        !  requested derivative (1st or 2nd) and the current value

        zdiff(1)=0.0_fp
        if(ibcymin.eq.1) then
           hy=y(2)-y(1)
           zdiff(1)=(f(1,ix,2)-f(1,ix,1))/hy + &
                hy*(-2*f(3,ix,1)-f(3,ix,2))/6
           zdiff(1)=bcymin(ix)-zdiff(1)
        else if(ibcymin.eq.2) then
           zdiff(1)=bcymin(ix)-f(3,ix,1)
        end if

        zdiff(2)=0.0_fp
        if(ibcymax.eq.1) then
           hy=y(ny)-y(ny-1)
           zdiff(2)=(f(1,ix,ny)-f(1,ix,ny-1))/hy + &
                hy*(2*f(3,ix,ny)+f(3,ix,ny-1))/6
           zdiff(2)=bcymax(ix)-zdiff(2)
        else if(ibcymax.eq.2) then
           zdiff(2)=bcymax(ix)-f(3,ix,ny)
        end if
        !
        fwk(1,1:ny)=0.0_fp  ! values are zero; only BC is not
        call mkspline(y,ny,fwk,ibcymin,zdiff(1),ibcymax,zdiff(2), &
             iliny,ier)
        if(ier.ne.0) return
        fcorr(1,ix,1:ny)=fwk(2,1:ny)  ! fyy-correction
     end do
     !
     zbcmin=0
     zbcmax=0
     do iy=1,ny
        fwk(1,1:nx)=fcorr(1,1:nx,iy)
        call mkspline(x,nx,fwk,ibcxmin,zbcmin,ibcxmax,zbcmax, &
             ilinx,ier)
        if(ier.ne.0) return
        fcorr(2,1:nx,iy)=fwk(2,1:nx)  ! fxxyy-correction
     end do
     !
     f(3:4,1:nx,1:ny)=f(3:4,1:nx,1:ny)+fcorr(1:2,1:nx,1:ny)
     !
     deallocate(fcorr)
  end if
  !
  !  correction spline -- f=fxx=zero; fyy & fxxyy are affected
  !
  deallocate(fwk)
  !------------------------------------
  !
  !  thats all
  !
  return
end subroutine mkbicub
