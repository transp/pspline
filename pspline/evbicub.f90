subroutine evbicub(xget,yget,x,nx,y,ny,ilinx,iliny, &
     f,inf2,ict,fval,ier)
  use precision_mod, only: fp
  !
  !  evaluate a 2d cubic Spline interpolant on a rectilinear
  !  grid -- this is C2 in both directions.
  !
  !  this subroutine calls two subroutines:
  !     herm2xy  -- find cell containing (xget,yget)
  !     fvbicub  -- evaluate interpolant function and (optionally) derivatives
  !
  !  input arguments:
  !  ================
  !
  !============
  implicit none
  integer inf2
  !============
  integer nx,ny                     ! grid sizes
  real(fp) :: xget,yget                    ! target of this interpolation
  real(fp) :: x(nx)                        ! ordered x grid
  real(fp) :: y(ny)                        ! ordered y grid
  integer ilinx                     ! ilinx=1 => assume x evenly spaced
  integer iliny                     ! iliny=1 => assume y evenly spaced
  !
  real(fp) :: f(0:3,inf2,ny)               ! function data
  !
  !       f 2nd dimension inf2 must be .ge. nx
  !       contents of f:
  !
  !  f(0,i,j) = f @ x(i),y(j)
  !  f(1,i,j) = d2f/dx2 @ x(i),y(j)
  !  f(2,i,j) = d2f/dy2 @ x(i),y(j)
  !  f(3,i,j) = d4f/dx2dy2 @ x(i),y(j)
  !
  !      (these are spline coefficients selected for continuous 2-
  !      diffentiability, see mkbicub[w].f90)
  !
  integer ict(6)                    ! code specifying output desired
  !
  !  ict(1)=1 -- return f  (0, don't)
  !  ict(2)=1 -- return df/dx  (0, don't)
  !  ict(3)=1 -- return df/dy  (0, don't)
  !  ict(4)=1 -- return d2f/dx2  (0, don't)
  !  ict(5)=1 -- return d2f/dy2  (0, don't)
  !  ict(6)=1 -- return d2f/dxdy (0, don't)
  !                   the number of non zero values ict(1:6)
  !                   determines the number of outputs...
  !
  !  new dmc December 2005 -- access to higher derivatives (even if not
  !  continuous-- but can only go up to 3rd derivatives on any one coordinate.
  !     if ict(1)=3 -- want 3rd derivatives
  !          ict(2)=1 for d3f/dx3
  !          ict(3)=1 for d3f/dx2dy
  !          ict(4)=1 for d3f/dxdy2
  !          ict(5)=1 for d3f/dy3
  !               number of non-zero values ict(2:5) gives no. of outputs
  !     if ict(1)=4 -- want 4th derivatives
  !          ict(2)=1 for d4f/dx3dy
  !          ict(3)=1 for d4f/dx2dy2
  !          ict(4)=1 for d4f/dxdy3
  !               number of non-zero values ict(2:4) gives no. of outputs
  !     if ict(1)=5 -- want 5th derivatives
  !          ict(2)=1 for d5f/dx3dy2
  !          ict(3)=1 for d5f/dx2dy3
  !               number of non-zero values ict(2:3) gives no. of outputs
  !     if ict(1)=6 -- want 6th derivatives
  !          d6f/dx3dy3 -- one value is returned.
  !
  ! output arguments:
  ! =================
  !
  real(fp) :: fval(*)                      ! output data
  integer ier                       ! error code =0 ==> no error
  !
  !  fval(1) receives the first output (depends on ict(...) spec)
  !  fval(2) receives the second output (depends on ict(...) spec)
  !  fval(3) receives the third output (depends on ict(...) spec)
  !  fval(4) receives the fourth output (depends on ict(...) spec)
  !  fval(5) receives the fourth output (depends on ict(...) spec)
  !  fval(6) receives the fourth output (depends on ict(...) spec)
  !
  !  examples:
  !    on input ict = [1,1,1,0,0,1]
  !   on output fval= [f,df/dx,df/dy,d2f/dxdy], elements 5 & 6 not referenced.
  !
  !    on input ict = [1,0,0,0,0,0]
  !   on output fval= [f] ... elements 2 -- 6 never referenced.
  !
  !    on input ict = [0,0,0,1,1,0]
  !   on output fval= [d2f/dx2,d2f/dy2] ... elements 3 -- 6 never referenced.
  !
  !    on input ict = [0,0,1,0,0,0]
  !   on output fval= [df/dy] ... elements 2 -- 6 never referenced.
  !
  !  ier -- completion code:  0 means OK
  !-------------------
  !  local:
  !
  integer i(1),j(1)                       ! cell indices
  !
  !  normalized displacement from (x(i),y(j)) corner of cell.
  !    xparam=0 @x(i)  xparam=1 @x(i+1)
  !    yparam=0 @y(j)  yparam=1 @y(j+1)
  !
  real(fp) :: xparam(1),yparam(1)
  !
  !  cell dimensions and
  !  inverse cell dimensions hxi = 1/(x(i+1)-x(i)), hyi = 1/(y(j+1)-y(j))
  !
  real(fp) :: hx(1),hy(1)
  real(fp) :: hxi(1),hyi(1)
  !
  !  0 .le. xparam .le. 1
  !  0 .le. yparam .le. 1
  !
  !  ** the interface is very similar to herm2ev.f90; can use herm2xy **
  !---------------------------------------------------------------------
  !
  i(1)=0
  j(1)=0
  call herm2xy(xget,yget,x,nx,y,ny,ilinx,iliny, &
       i(1),j(1),xparam(1),yparam(1),hx(1),hxi(1),hy(1),hyi(1),ier)
  if(ier.ne.0) return
  !
  call fvbicub(ict,1,1, &
       fval,i,j,xparam,yparam,hx,hxi,hy,hyi,f,inf2,ny)
  !
  return
end subroutine evbicub
!---------------------------------------------------------------------
!  evaluate C1 cubic Hermite function interpolation -- 2d fcn
!   --vectorized-- dmc 10 Feb 1999
!
subroutine fvbicub(ict,ivec,ivecd, &
     fval,ii,jj,xparam,yparam,hx,hxi,hy,hyi, &
     fin,inf2,ny)
  use precision_mod, only: fp
  !
  !============
  implicit none
  integer ny,inf2,iadr,i,j
  !============
  real(fp) :: z36th,xp,xpi,xp2,xpi2,cx,cxi,hx2,yp,ypi,yp2,ypi2,cy
  real(fp) :: cyi,hy2,cxd,cxdi,cyd,cydi
  !============
  integer ict(6)                    ! requested output control
  integer ivec                      ! vector length
  integer ivecd                     ! vector dimension (1st dim of fval)
  !
  integer ii(ivec),jj(ivec)         ! target cells (i,j)
  real(fp) :: xparam(ivec),yparam(ivec)
  ! normalized displacements from (i,j) corners
  !
  real(fp) :: hx(ivec),hy(ivec)            ! grid spacing, and
  real(fp) :: hxi(ivec),hyi(ivec)          ! inverse grid spacing 1/(x(i+1)-x(i))
  ! & 1/(y(j+1)-y(j))
  !
  real(fp) :: fin(0:3,inf2,ny)             ! interpolant data (cf "evbicub")
  !
  real(fp) :: fval(ivecd,*)                ! output returned
  !
  !  for detailed description of fin, ict and fval see subroutine
  !  evbicub comments.  Note ict is not vectorized; the same output
  !  is expected to be returned for all input vector data points.
  !
  !  note that the index inputs ii,jj and parameter inputs
  !     xparam,yparam,hx,hxi,hy,hyi are vectorized, and the
  !     output array fval has a vector ** 1st dimension ** whose
  !     size must be given as a separate argument
  !
  !  to use this routine in scalar mode, pass in ivec=ivecd=1
  !
  !---------------
  !  Spline evaluation consists of a "mixing" of the interpolant
  !  data using the linear functionals xparam, xpi = 1-xparam,
  !  yparam, ypi = 1-yparam, and the cubic functionals
  !  xparam**3-xparam, xpi**3-xpi, yparam**3-yparam, ypi**3-ypi ...
  !  and their derivatives as needed.
  !
  integer v
  real(fp) :: sum
  !
  real(fp), parameter :: sixth = 0.166666666666666667_fp
  !
  !---------------
  !   ...in x direction
  !
  z36th=sixth*sixth
  iadr=0
  !
  if(ict(1).le.2) then
     !
     !  get desired values:
     !
     if(ict(1).eq.1) then
        !
        !  function value:
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           j=jj(v)
           !
           !  in x direction...
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           xp2=xp*xp
           xpi2=xpi*xpi
           !
           cx=xp*(xp2-1.0_fp)
           cxi=xpi*(xpi2-1.0_fp)
           hx2=hx(v)*hx(v)
           !
           !   ...and in y direction
           !
           yp=yparam(v)
           ypi=1.0_fp-yp
           yp2=yp*yp
           ypi2=ypi*ypi
           !
           cy=yp*(yp2-1.0_fp)
           cyi=ypi*(ypi2-1.0_fp)
           hy2=hy(v)*hy(v)
           !
           sum=xpi*(ypi*fin(0,i,j)  +yp*fin(0,i,j+1))+ &
                xp*(ypi*fin(0,i+1,j)+yp*fin(0,i+1,j+1))
           !
           sum=sum+sixth*hx2*( &
                cxi*(ypi*fin(1,i,j)  +yp*fin(1,i,j+1))+ &
                cx*(ypi*fin(1,i+1,j)+yp*fin(1,i+1,j+1)))
           !
           sum=sum+sixth*hy2*( &
                xpi*(cyi*fin(2,i,j)  +cy*fin(2,i,j+1))+ &
                xp*(cyi*fin(2,i+1,j)+cy*fin(2,i+1,j+1)))
           !
           sum=sum+z36th*hx2*hy2*( &
                cxi*(cyi*fin(3,i,j)  +cy*fin(3,i,j+1))+ &
                cx*(cyi*fin(3,i+1,j)+cy*fin(3,i+1,j+1)))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(2).eq.1) then
        !
        !  df/dx:
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           j=jj(v)
           !
           !  in x direction...
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           xp2=xp*xp
           xpi2=xpi*xpi

           cxd=3.0_fp*xp2-1.0_fp
           cxdi=-3.0_fp*xpi2+1.0_fp
           !
           !   ...and in y direction
           !
           yp=yparam(v)
           ypi=1.0_fp-yp
           yp2=yp*yp
           ypi2=ypi*ypi
           !
           cy=yp*(yp2-1.0_fp)
           cyi=ypi*(ypi2-1.0_fp)
           hy2=hy(v)*hy(v)
           !
           sum=hxi(v)*( &
                -(ypi*fin(0,i,j)  +yp*fin(0,i,j+1)) &
                +(ypi*fin(0,i+1,j)+yp*fin(0,i+1,j+1)))
           !
           sum=sum+sixth*hx(v)*( &
                cxdi*(ypi*fin(1,i,j)  +yp*fin(1,i,j+1))+ &
                cxd*(ypi*fin(1,i+1,j)+yp*fin(1,i+1,j+1)))
           !
           sum=sum+sixth*hxi(v)*hy2*( &
                -(cyi*fin(2,i,j)  +cy*fin(2,i,j+1)) &
                +(cyi*fin(2,i+1,j)+cy*fin(2,i+1,j+1)))
           !
           sum=sum+z36th*hx(v)*hy2*( &
                cxdi*(cyi*fin(3,i,j)  +cy*fin(3,i,j+1))+ &
                cxd*(cyi*fin(3,i+1,j)+cy*fin(3,i+1,j+1)))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(3).eq.1) then
        !
        !  df/dy:
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           j=jj(v)
           !
           !  in x direction...
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           xp2=xp*xp
           xpi2=xpi*xpi
           !
           cx=xp*(xp2-1.0_fp)
           cxi=xpi*(xpi2-1.0_fp)
           hx2=hx(v)*hx(v)
           !
           !   ...and in y direction
           !
           yp=yparam(v)
           ypi=1.0_fp-yp
           yp2=yp*yp
           ypi2=ypi*ypi

           cyd=3.0_fp*yp2-1.0_fp
           cydi=-3.0_fp*ypi2+1.0_fp
           !
           sum=hyi(v)*( &
                xpi*(-fin(0,i,j)  +fin(0,i,j+1))+ &
                xp*(-fin(0,i+1,j)+fin(0,i+1,j+1)))
           !
           sum=sum+sixth*hx2*hyi(v)*( &
                cxi*(-fin(1,i,j)  +fin(1,i,j+1))+ &
                cx*(-fin(1,i+1,j)+fin(1,i+1,j+1)))
           !
           sum=sum+sixth*hy(v)*( &
                xpi*(cydi*fin(2,i,j)  +cyd*fin(2,i,j+1))+ &
                xp*(cydi*fin(2,i+1,j)+cyd*fin(2,i+1,j+1)))
           !
           sum=sum+z36th*hx2*hy(v)*( &
                cxi*(cydi*fin(3,i,j)  +cyd*fin(3,i,j+1))+ &
                cx*(cydi*fin(3,i+1,j)+cyd*fin(3,i+1,j+1)))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(4).eq.1) then
        !
        !  d2f/dx2:
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           j=jj(v)
           !
           !  in x direction...
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           !
           !   ...and in y direction
           !
           yp=yparam(v)
           ypi=1.0_fp-yp
           yp2=yp*yp
           ypi2=ypi*ypi
           !
           cy=yp*(yp2-1.0_fp)
           cyi=ypi*(ypi2-1.0_fp)
           hy2=hy(v)*hy(v)
           !
           sum=( &
                xpi*(ypi*fin(1,i,j)  +yp*fin(1,i,j+1))+ &
                xp*(ypi*fin(1,i+1,j)+yp*fin(1,i+1,j+1)))
           !
           sum=sum+sixth*hy2*( &
                xpi*(cyi*fin(3,i,j)  +cy*fin(3,i,j+1))+ &
                xp*(cyi*fin(3,i+1,j)+cy*fin(3,i+1,j+1)))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(5).eq.1) then
        !
        !  d2f/dy2:
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           j=jj(v)
           !
           !  in x direction...
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           xp2=xp*xp
           xpi2=xpi*xpi
           !
           cx=xp*(xp2-1.0_fp)
           cxi=xpi*(xpi2-1.0_fp)
           hx2=hx(v)*hx(v)
           !
           !   ...and in y direction
           !
           yp=yparam(v)
           ypi=1.0_fp-yp
           !
           sum=( &
                xpi*(ypi*fin(2,i,j)  +yp*fin(2,i,j+1))+ &
                xp*(ypi*fin(2,i+1,j)+yp*fin(2,i+1,j+1)))
           !
           sum=sum+sixth*hx2*( &
                cxi*(ypi*fin(3,i,j)  +yp*fin(3,i,j+1))+ &
                cx*(ypi*fin(3,i+1,j)+yp*fin(3,i+1,j+1)))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(6).eq.1) then
        !
        !  d2f/dxdy:
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           j=jj(v)
           !
           !  in x direction...
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           xp2=xp*xp
           xpi2=xpi*xpi

           cxd=3.0_fp*xp2-1.0_fp
           cxdi=-3.0_fp*xpi2+1.0_fp
           !
           !   ...and in y direction
           !
           yp=yparam(v)
           ypi=1.0_fp-yp
           yp2=yp*yp
           ypi2=ypi*ypi

           cyd=3.0_fp*yp2-1.0_fp
           cydi=-3.0_fp*ypi2+1.0_fp
           !
           sum=hxi(v)*hyi(v)*( &
                fin(0,i,j)  -fin(0,i,j+1) &
                -fin(0,i+1,j)+fin(0,i+1,j+1))
           !
           sum=sum+sixth*hx(v)*hyi(v)*( &
                cxdi*(-fin(1,i,j)  +fin(1,i,j+1))+ &
                cxd*(-fin(1,i+1,j)+fin(1,i+1,j+1)))
           !
           sum=sum+sixth*hxi(v)*hy(v)*( &
                -(cydi*fin(2,i,j)  +cyd*fin(2,i,j+1)) &
                +(cydi*fin(2,i+1,j)+cyd*fin(2,i+1,j+1)))
           !
           sum=sum+z36th*hx(v)*hy(v)*( &
                cxdi*(cydi*fin(3,i,j)  +cyd*fin(3,i,j+1))+ &
                cxd*(cydi*fin(3,i+1,j)+cyd*fin(3,i+1,j+1)))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     !-------------------------------------------------
     !
  else if(ict(1).eq.3) then
     if(ict(2).eq.1) then
        !  evaluate d3f/dx3 (not continuous)
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           j=jj(v)
           yp=yparam(v)
           ypi=1.0_fp-yp
           yp2=yp*yp
           ypi2=ypi*ypi
           cy=yp*(yp2-1.0_fp)
           cyi=ypi*(ypi2-1.0_fp)
           hy2=hy(v)*hy(v)
           sum=hxi(v)*( &
                -(ypi*fin(1,i,j)  +yp*fin(1,i,j+1)) &
                +(ypi*fin(1,i+1,j)+yp*fin(1,i+1,j+1)))
           !
           sum=sum+sixth*hy2*hxi(v)*( &
                -(cyi*fin(3,i,j)  +cy*fin(3,i,j+1)) &
                +(cyi*fin(3,i+1,j)+cy*fin(3,i+1,j+1)))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(3).eq.1) then
        !  evaluate d3f/dx2dy
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           j=jj(v)
           xp=xparam(v)
           xpi=1.0_fp-xp
           yp=yparam(v)
           ypi=1.0_fp-yp
           yp2=yp*yp
           ypi2=ypi*ypi
           cyd=3.0_fp*yp2-1.0_fp
           cydi=-3.0_fp*ypi2+1.0_fp
           !
           sum=hyi(v)*( &
                xpi*(-fin(1,i,j)  +fin(1,i,j+1))+ &
                xp*(-fin(1,i+1,j) +fin(1,i+1,j+1)))
           !
           sum=sum+sixth*hy(v)*( &
                xpi*(cydi*fin(3,i,j) +cyd*fin(3,i,j+1))+ &
                xp*(cydi*fin(3,i+1,j)+cyd*fin(3,i+1,j+1)))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(4).eq.1) then
        !  evaluate d3f/dxdy2
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           j=jj(v)
           xp=xparam(v)
           xpi=1.0_fp-xp
           xp2=xp*xp
           xpi2=xpi*xpi
           cxd=3.0_fp*xp2-1.0_fp
           cxdi=-3.0_fp*xpi2+1.0_fp
           yp=yparam(v)
           ypi=1.0_fp-yp
           !
           sum=hxi(v)*( &
                -(ypi*fin(2,i,j)  +yp*fin(2,i,j+1)) &
                +(ypi*fin(2,i+1,j)+yp*fin(2,i+1,j+1)))
           !
           sum=sum+sixth*hx(v)*( &
                cxdi*(ypi*fin(3,i,j)  +yp*fin(3,i,j+1))+ &
                cxd*(ypi*fin(3,i+1,j)+yp*fin(3,i+1,j+1)))
           !
           fval(v,iadr)=sum
        end do
     end if

     if(ict(5).eq.1) then
        !  evaluate d3f/dy3 (not continuous)
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           j=jj(v)
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           xp2=xp*xp
           xpi2=xpi*xpi
           !
           cx=xp*(xp2-1.0_fp)
           cxi=xpi*(xpi2-1.0_fp)
           hx2=hx(v)*hx(v)
           !
           sum=hyi(v)*( &
                xpi*(-fin(2,i,j)  +fin(2,i,j+1))+ &
                xp*(-fin(2,i+1,j) +fin(2,i+1,j+1)))
           !
           sum=sum+sixth*hx2*hyi(v)*( &
                cxi*(-fin(3,i,j)  +fin(3,i,j+1))+ &
                cx*(-fin(3,i+1,j) +fin(3,i+1,j+1)))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     !-----------------------------------
     !  access to 4th derivatives
     !
  else if(ict(1).eq.4) then
     if(ict(2).eq.1) then
        !  evaluate d4f/dx3dy (not continuous)
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           j=jj(v)
           yp=yparam(v)
           ypi=1.0_fp-yp
           yp2=yp*yp
           ypi2=ypi*ypi
           cyd=3.0_fp*yp2-1.0_fp
           cydi=-3.0_fp*ypi2+1.0_fp
           !
           sum=hxi(v)*hyi(v)*( &
                +( fin(1,i,j)  -fin(1,i,j+1)) &
                +(-fin(1,i+1,j)+fin(1,i+1,j+1)))
           !
           sum=sum+sixth*hy(v)*hxi(v)*( &
                -(cydi*fin(3,i,j)  +cyd*fin(3,i,j+1)) &
                +(cydi*fin(3,i+1,j)+cyd*fin(3,i+1,j+1)))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(3).eq.1) then
        !  evaluate d4f/dx2dy2
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           j=jj(v)
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           yp=yparam(v)
           ypi=1.0_fp-yp
           !
           sum=xpi*(ypi*fin(3,i,j)  +yp*fin(3,i,j+1))+ &
                xp*(ypi*fin(3,i+1,j)+yp*fin(3,i+1,j+1))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(4).eq.1) then
        !  evaluate d4f/dxdy3 (not continuous)
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           j=jj(v)
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           xp2=xp*xp
           xpi2=xpi*xpi
           !
           cxd=3.0_fp*xp2-1.0_fp
           cxdi=-3.0_fp*xpi2+1.0_fp
           !
           sum=hyi(v)*hxi(v)*( &
                +( fin(2,i,j)  -fin(2,i,j+1)) &
                +(-fin(2,i+1,j)+fin(2,i+1,j+1)))
           !
           sum=sum+sixth*hx(v)*hyi(v)*( &
                cxdi*(-fin(3,i,j)  +fin(3,i,j+1))+ &
                cxd*(-fin(3,i+1,j) +fin(3,i+1,j+1)))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     !-----------------------------------
     !  access to 5th derivatives
     !
  else if(ict(1).eq.5) then
     if(ict(2).eq.1) then
        !  evaluate d5f/dx3dy2 (not continuous)
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           j=jj(v)
           !
           yp=yparam(v)
           ypi=1.0_fp-yp
           !
           sum=hxi(v)*( &
                -(ypi*fin(3,i,j)  +yp*fin(3,i,j+1)) &
                +(ypi*fin(3,i+1,j)+yp*fin(3,i+1,j+1)))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(3).eq.1) then
        !  evaluate d5f/dx2dy3 (not continuous)
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           j=jj(v)
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           !
           sum=hyi(v)*( &
                xpi*(-fin(3,i,j)  +fin(3,i,j+1))+ &
                xp*(-fin(3,i+1,j)+fin(3,i+1,j+1)))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     !-----------------------------------
     !  access to 6th derivatives
     !
  else if(ict(1).eq.6) then
     !  evaluate d6f/dx3dy3 (not continuous)
     iadr=iadr+1
     do v=1,ivec
        i=ii(v)
        j=jj(v)
        sum=hxi(v)*hyi(v)*( &
             +( fin(3,i,j)  -fin(3,i,j+1)) &
             +(-fin(3,i+1,j)+fin(3,i+1,j+1)))
        fval(v,iadr)=sum
     end do
  end if
  !
  return
end subroutine fvbicub
!---------------------------------------------------------------------
!  evaluate C1 cubic Hermite function interpolation -- 2d fcn
!   --vectorized-- dmc 10 Feb 1999
!    --optimized for VARIATION along x axis ONLY--
!
subroutine fvbicubx(ict,ivec,ivecd, &
     fval,ii,jj,xparam,yparam,hx,hxi,hy,hyi, &
     fin,inf2,ny)
  use precision_mod, only: fp
  !
  !============
  implicit none
  integer ny,inf2,iadr,j,i
  !============
  real(fp) :: z36th,yp,ypi,yp2,ypi2,cy,cyi,hy2,xp,xpi,xp2,xpi2,cx
  real(fp) :: cxi,hx2,cxd,cxdi,cyd,cydi
  !============
  integer ict(6)                    ! requested output control
  integer ivec                      ! vector length
  integer ivecd                     ! vector dimension (1st dim of fval)
  !
  integer ii(ivec),jj               ! target cells (i,j)
  real(fp) :: xparam(ivec),yparam
  ! normalized displacements from (i,j) corners
  !
  real(fp) :: hx(ivec),hy                  ! grid spacing, and
  real(fp) :: hxi(ivec),hyi                ! inverse grid spacing 1/(x(i+1)-x(i))
  ! & 1/(y(j+1)-y(j))
  !
  real(fp) :: fin(0:3,inf2,ny)             ! interpolant data (cf "evbicub")
  !
  real(fp) :: fval(ivecd,*)                ! output returned
  !
  !  for detailed description of fin, ict and fval see subroutine
  !  evbicub comments.  Note ict is not vectorized; the same output
  !  is expected to be returned for all input vector data points.
  !
  !  note that the index inputs ii,jj and parameter inputs
  !     xparam,yparam,hx,hxi,hy,hyi are vectorized, and the
  !     output array fval has a vector ** 1st dimension ** whose
  !     size must be given as a separate argument
  !
  !  to use this routine in scalar mode, pass in ivec=ivecd=1
  !
  !---------------
  !  Spline evaluation consists of a "mixing" of the interpolant
  !  data using the linear functionals xparam, xpi = 1-xparam,
  !  yparam, ypi = 1-yparam, and the cubic functionals
  !  xparam**3-xparam, xpi**3-xpi, yparam**3-yparam, ypi**3-ypi ...
  !  and their derivatives as needed.
  !
  integer v
  real(fp) :: sum
  !
  real(fp), parameter :: sixth = 0.166666666666666667_fp
  !
  !---------------
  !   ...in x direction
  !
  z36th=sixth*sixth
  iadr=0
  !
  if(ict(1).le.2) then
     !
     !  get desired values:
     !
     if(ict(1).eq.1) then
        !
        !  function value:
        !
        j=jj
        !
        !   ...and in y direction
        !
        yp=yparam
        ypi=1.0_fp-yp
        yp2=yp*yp
        ypi2=ypi*ypi
        !
        cy=yp*(yp2-1.0_fp)
        cyi=ypi*(ypi2-1.0_fp)
        hy2=hy*hy
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           !
           !  in x direction...
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           xp2=xp*xp
           xpi2=xpi*xpi
           !
           cx=xp*(xp2-1.0_fp)
           cxi=xpi*(xpi2-1.0_fp)
           hx2=hx(v)*hx(v)
           !
           sum=xpi*(ypi*fin(0,i,j)  +yp*fin(0,i,j+1))+ &
                xp*(ypi*fin(0,i+1,j)+yp*fin(0,i+1,j+1))
           !
           sum=sum+sixth*hx2*( &
                cxi*(ypi*fin(1,i,j)  +yp*fin(1,i,j+1))+ &
                cx*(ypi*fin(1,i+1,j)+yp*fin(1,i+1,j+1)))
           !
           sum=sum+sixth*hy2*( &
                xpi*(cyi*fin(2,i,j)  +cy*fin(2,i,j+1))+ &
                xp*(cyi*fin(2,i+1,j)+cy*fin(2,i+1,j+1)))
           !
           sum=sum+z36th*hx2*hy2*( &
                cxi*(cyi*fin(3,i,j)  +cy*fin(3,i,j+1))+ &
                cx*(cyi*fin(3,i+1,j)+cy*fin(3,i+1,j+1)))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(2).eq.1) then
        !
        !  df/dx:
        !
        j=jj
        !
        !   ...and in y direction
        !
        yp=yparam
        ypi=1.0_fp-yp
        yp2=yp*yp
        ypi2=ypi*ypi
        !
        cy=yp*(yp2-1.0_fp)
        cyi=ypi*(ypi2-1.0_fp)
        hy2=hy*hy
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           !
           !  in x direction...
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           xp2=xp*xp
           xpi2=xpi*xpi

           cxd=3.0_fp*xp2-1.0_fp
           cxdi=-3.0_fp*xpi2+1.0_fp
           !
           sum=hxi(v)*( &
                -(ypi*fin(0,i,j)  +yp*fin(0,i,j+1)) &
                +(ypi*fin(0,i+1,j)+yp*fin(0,i+1,j+1)))
           !
           sum=sum+sixth*hx(v)*( &
                cxdi*(ypi*fin(1,i,j)  +yp*fin(1,i,j+1))+ &
                cxd*(ypi*fin(1,i+1,j)+yp*fin(1,i+1,j+1)))
           !
           sum=sum+sixth*hxi(v)*hy2*( &
                -(cyi*fin(2,i,j)  +cy*fin(2,i,j+1)) &
                +(cyi*fin(2,i+1,j)+cy*fin(2,i+1,j+1)))
           !
           sum=sum+z36th*hx(v)*hy2*( &
                cxdi*(cyi*fin(3,i,j)  +cy*fin(3,i,j+1))+ &
                cxd*(cyi*fin(3,i+1,j)+cy*fin(3,i+1,j+1)))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(3).eq.1) then
        !
        !  df/dy:
        !
        j=jj
        !
        !   ...and in y direction
        !
        yp=yparam
        ypi=1.0_fp-yp
        yp2=yp*yp
        ypi2=ypi*ypi

        cyd=3.0_fp*yp2-1.0_fp
        cydi=-3.0_fp*ypi2+1.0_fp
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           !
           !  in x direction...
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           xp2=xp*xp
           xpi2=xpi*xpi
           !
           cx=xp*(xp2-1.0_fp)
           cxi=xpi*(xpi2-1.0_fp)
           hx2=hx(v)*hx(v)
           !
           sum=hyi*( &
                xpi*(-fin(0,i,j)  +fin(0,i,j+1))+ &
                xp*(-fin(0,i+1,j)+fin(0,i+1,j+1)))
           !
           sum=sum+sixth*hx2*hyi*( &
                cxi*(-fin(1,i,j)  +fin(1,i,j+1))+ &
                cx*(-fin(1,i+1,j)+fin(1,i+1,j+1)))
           !
           sum=sum+sixth*hy*( &
                xpi*(cydi*fin(2,i,j)  +cyd*fin(2,i,j+1))+ &
                xp*(cydi*fin(2,i+1,j)+cyd*fin(2,i+1,j+1)))
           !
           sum=sum+z36th*hx2*hy*( &
                cxi*(cydi*fin(3,i,j)  +cyd*fin(3,i,j+1))+ &
                cx*(cydi*fin(3,i+1,j)+cyd*fin(3,i+1,j+1)))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(4).eq.1) then
        !
        !  d2f/dx2:
        !
        j=jj
        !
        !   ...and in y direction
        !
        yp=yparam
        ypi=1.0_fp-yp
        yp2=yp*yp
        ypi2=ypi*ypi
        !
        cy=yp*(yp2-1.0_fp)
        cyi=ypi*(ypi2-1.0_fp)
        hy2=hy*hy
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           !
           !  in x direction...
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           !
           sum=( &
                xpi*(ypi*fin(1,i,j)  +yp*fin(1,i,j+1))+ &
                xp*(ypi*fin(1,i+1,j)+yp*fin(1,i+1,j+1)))
           !
           sum=sum+sixth*hy2*( &
                xpi*(cyi*fin(3,i,j)  +cy*fin(3,i,j+1))+ &
                xp*(cyi*fin(3,i+1,j)+cy*fin(3,i+1,j+1)))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(5).eq.1) then
        !
        !  d2f/dy2:
        !
        j=jj
        !
        !   ...and in y direction
        !
        yp=yparam
        ypi=1.0_fp-yp
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           !
           !  in x direction...
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           xp2=xp*xp
           xpi2=xpi*xpi
           !
           cx=xp*(xp2-1.0_fp)
           cxi=xpi*(xpi2-1.0_fp)
           hx2=hx(v)*hx(v)
           !
           sum=( &
                xpi*(ypi*fin(2,i,j)  +yp*fin(2,i,j+1))+ &
                xp*(ypi*fin(2,i+1,j)+yp*fin(2,i+1,j+1)))
           !
           sum=sum+sixth*hx2*( &
                cxi*(ypi*fin(3,i,j)  +yp*fin(3,i,j+1))+ &
                cx*(ypi*fin(3,i+1,j)+yp*fin(3,i+1,j+1)))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(6).eq.1) then
        !
        !  d2f/dxdy:
        !
        j=jj
        !
        !   ...and in y direction
        !
        yp=yparam
        ypi=1.0_fp-yp
        yp2=yp*yp
        ypi2=ypi*ypi

        cyd=3.0_fp*yp2-1.0_fp
        cydi=-3.0_fp*ypi2+1.0_fp
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           !
           !  in x direction...
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           xp2=xp*xp
           xpi2=xpi*xpi

           cxd=3.0_fp*xp2-1.0_fp
           cxdi=-3.0_fp*xpi2+1.0_fp
           !
           sum=hxi(v)*hyi*( &
                fin(0,i,j)  -fin(0,i,j+1) &
                -fin(0,i+1,j)+fin(0,i+1,j+1))
           !
           sum=sum+sixth*hx(v)*hyi*( &
                cxdi*(-fin(1,i,j)  +fin(1,i,j+1))+ &
                cxd*(-fin(1,i+1,j)+fin(1,i+1,j+1)))
           !
           sum=sum+sixth*hxi(v)*hy*( &
                -(cydi*fin(2,i,j)  +cyd*fin(2,i,j+1)) &
                +(cydi*fin(2,i+1,j)+cyd*fin(2,i+1,j+1)))
           !
           sum=sum+z36th*hx(v)*hy*( &
                cxdi*(cydi*fin(3,i,j)  +cyd*fin(3,i,j+1))+ &
                cxd*(cydi*fin(3,i+1,j)+cyd*fin(3,i+1,j+1)))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     !-------------------------------------------------
     !
  else if(ict(1).eq.3) then
     if(ict(2).eq.1) then
        !  evaluate d3f/dx3 (not continuous)
        j=jj
        yp=yparam
        ypi=1.0_fp-yp
        yp2=yp*yp
        ypi2=ypi*ypi
        cy=yp*(yp2-1.0_fp)
        cyi=ypi*(ypi2-1.0_fp)
        hy2=hy*hy
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           sum=hxi(v)*( &
                -(ypi*fin(1,i,j)  +yp*fin(1,i,j+1)) &
                +(ypi*fin(1,i+1,j)+yp*fin(1,i+1,j+1)))
           !
           sum=sum+sixth*hy2*hxi(v)*( &
                -(cyi*fin(3,i,j)  +cy*fin(3,i,j+1)) &
                +(cyi*fin(3,i+1,j)+cy*fin(3,i+1,j+1)))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(3).eq.1) then
        !  evaluate d3f/dx2dy
        j=jj
        !
        yp=yparam
        ypi=1.0_fp-yp
        yp2=yp*yp
        ypi2=ypi*ypi
        cyd=3.0_fp*yp2-1.0_fp
        cydi=-3.0_fp*ypi2+1.0_fp
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           xp=xparam(v)
           xpi=1.0_fp-xp
           !
           sum=hyi*( &
                xpi*(-fin(1,i,j)  +fin(1,i,j+1))+ &
                xp*(-fin(1,i+1,j) +fin(1,i+1,j+1)))
           !
           sum=sum+sixth*hy*( &
                xpi*(cydi*fin(3,i,j) +cyd*fin(3,i,j+1))+ &
                xp*(cydi*fin(3,i+1,j)+cyd*fin(3,i+1,j+1)))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(4).eq.1) then
        !  evaluate d3f/dxdy2
        j=jj
        yp=yparam
        ypi=1.0_fp-yp
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           xp=xparam(v)
           xpi=1.0_fp-xp
           xp2=xp*xp
           xpi2=xpi*xpi
           cxd=3.0_fp*xp2-1.0_fp
           cxdi=-3.0_fp*xpi2+1.0_fp
           !
           sum=hxi(v)*( &
                -(ypi*fin(2,i,j)  +yp*fin(2,i,j+1)) &
                +(ypi*fin(2,i+1,j)+yp*fin(2,i+1,j+1)))
           !
           sum=sum+sixth*hx(v)*( &
                cxdi*(ypi*fin(3,i,j)  +yp*fin(3,i,j+1))+ &
                cxd*(ypi*fin(3,i+1,j)+yp*fin(3,i+1,j+1)))
           !
           fval(v,iadr)=sum
        end do
     end if

     if(ict(5).eq.1) then
        !  evaluate d3f/dy3 (not continuous)
        iadr=iadr+1
        j=jj
        do v=1,ivec
           i=ii(v)
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           xp2=xp*xp
           xpi2=xpi*xpi
           !
           cx=xp*(xp2-1.0_fp)
           cxi=xpi*(xpi2-1.0_fp)
           hx2=hx(v)*hx(v)
           !
           sum=hyi*( &
                xpi*(-fin(2,i,j)  +fin(2,i,j+1))+ &
                xp*(-fin(2,i+1,j) +fin(2,i+1,j+1)))
           !
           sum=sum+sixth*hx2*hyi*( &
                cxi*(-fin(3,i,j)  +fin(3,i,j+1))+ &
                cx*(-fin(3,i+1,j) +fin(3,i+1,j+1)))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     !-----------------------------------
     !  access to 4th derivatives
     !
  else if(ict(1).eq.4) then
     if(ict(2).eq.1) then
        !  evaluate d4f/dx3dy (not continuous)
        iadr=iadr+1
        j=jj
        !
        yp=yparam
        ypi=1.0_fp-yp
        yp2=yp*yp
        ypi2=ypi*ypi
        cyd=3.0_fp*yp2-1.0_fp
        cydi=-3.0_fp*ypi2+1.0_fp
        !
        do v=1,ivec
           i=ii(v)
           !
           sum=hxi(v)*hyi*( &
                +( fin(1,i,j)  -fin(1,i,j+1)) &
                +(-fin(1,i+1,j)+fin(1,i+1,j+1)))
           !
           sum=sum+sixth*hy*hxi(v)*( &
                -(cydi*fin(3,i,j)  +cyd*fin(3,i,j+1)) &
                +(cydi*fin(3,i+1,j)+cyd*fin(3,i+1,j+1)))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(3).eq.1) then
        !  evaluate d4f/dx2dy2
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           j=jj
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           yp=yparam
           ypi=1.0_fp-yp
           !
           sum=xpi*(ypi*fin(3,i,j)  +yp*fin(3,i,j+1))+ &
                xp*(ypi*fin(3,i+1,j)+yp*fin(3,i+1,j+1))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(4).eq.1) then
        !  evaluate d4f/dxdy3 (not continuous)
        j=jj
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           xp2=xp*xp
           xpi2=xpi*xpi
           !
           cxd=3.0_fp*xp2-1.0_fp
           cxdi=-3.0_fp*xpi2+1.0_fp
           !
           sum=hyi*hxi(v)*( &
                +( fin(2,i,j)  -fin(2,i,j+1)) &
                +(-fin(2,i+1,j)+fin(2,i+1,j+1)))
           !
           sum=sum+sixth*hx(v)*hyi*( &
                cxdi*(-fin(3,i,j)  +fin(3,i,j+1))+ &
                cxd*(-fin(3,i+1,j) +fin(3,i+1,j+1)))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     !-----------------------------------
     !  access to 5th derivatives
     !
  else if(ict(1).eq.5) then
     if(ict(2).eq.1) then
        !  evaluate d5f/dx3dy2 (not continuous)
        j=jj
        !
        yp=yparam
        ypi=1.0_fp-yp
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           !
           sum=hxi(v)*( &
                -(ypi*fin(3,i,j)  +yp*fin(3,i,j+1)) &
                +(ypi*fin(3,i+1,j)+yp*fin(3,i+1,j+1)))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(3).eq.1) then
        !  evaluate d5f/dx2dy3 (not continuous)
        j=jj
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           !
           sum=hyi*( &
                xpi*(-fin(3,i,j)  +fin(3,i,j+1))+ &
                xp*(-fin(3,i+1,j)+fin(3,i+1,j+1)))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     !-----------------------------------
     !  access to 6th derivatives
     !
  else if(ict(1).eq.6) then
     !  evaluate d6f/dx3dy3 (not continuous)
     iadr=iadr+1
     j=jj
     do v=1,ivec
        i=ii(v)
        sum=hxi(v)*hyi*( &
             +( fin(3,i,j)  -fin(3,i,j+1)) &
             +(-fin(3,i+1,j)+fin(3,i+1,j+1)))
        fval(v,iadr)=sum
     end do
  end if
  !
  return
end subroutine fvbicubx
