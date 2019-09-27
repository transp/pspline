subroutine evspline(xget,x,nx,ilinx,f,ict,fval,ier)
  use precision_mod, only: fp
  !
  !  use mkspline to set up spline coefficients...
  !
  !  evaluate a 1d cubic Spline interpolant -- this is C2
  !
  !  this subroutine calls two subroutines:
  !     herm1x  -- find cell containing (xget)
  !     fvspline  -- evaluate interpolant function and (optionally) derivatives
  !
  !  input arguments:
  !  ================
  !
  implicit none
  integer nx                        ! grid size
  real(fp) :: xget                         ! target of this interpolation
  real(fp) :: x(nx)                        ! ordered x grid
  integer ilinx                     ! ilinx=1 => assume x evenly spaced
  !
  real(fp) :: f(0:1,nx)                    ! function data
  !
  !  f(0,i) = f @ x(i)
  !  f(1,i) = d2f/dx2 @ x(i)
  !
  !      (these are spline coefficients selected for continuous 2-
  !      diffentiability, see mkspline.f90)
  !
  integer ict(3)                    ! code specifying output desired
  !
  !  ict(1)=1 -- return f  (0, don't)
  !  ict(2)=1 -- return df/dx  (0, don't)
  !  ict(3)=1 -- return d2f/dx2  (0, don't)
  !
  !        set ict(1)=3 to get d3f/dx3 (only)
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
  !
  !  examples:
  !    on input ict = [1,1,1]
  !   on output fval= [f,df/dx,d2f/dx2]
  !
  !    on input ict = [1,0,0]
  !   on output fval= [f] ... elements 2 -- 3 never referenced.
  !
  !    on input ict = [0,0,1]
  !   on output fval= [d2f/dx2] ... elements 2 -- 3 never referenced.
  !
  !  ier -- completion code:  0 means OK
  !-------------------
  !  local:
  !
  integer i(1)                         ! cell indices
  !
  !  normalized displacement from x(i) grid point.
  !    xparam=0 @x(i)  xparam=1 @x(i+1)
  !
  real(fp) :: xparam(1)
  !
  !  cell dimensions and
  !  inverse cell dimensions hxi = 1/(x(i+1)-x(i))
  !
  real(fp) :: hx(1)
  real(fp) :: hxi(1)
  !
  !  0 .le. xparam .le. 1
  !
  !  ** the interface is very similar to herm2ev.f90; can use herm2xy **
  !---------------------------------------------------------------------
  !
  i(1)=0
  call herm1x(xget,x,nx,ilinx,i(1),xparam(1),hx(1),hxi(1),ier)
  if(ier.ne.0) return
  !
  call fvspline(ict,1,1,fval,i,xparam,hx,hxi,f,nx)
  !
  return
end subroutine evspline
!---------------------------------------------------------------------
!  evaluate C1 cubic Hermite function interpolation -- 2d fcn
!   --vectorized-- dmc 10 Feb 1999
!
subroutine fvspline(ict,ivec,ivecd, &
     fval,ii,xparam,hx,hxi,fin,nx)
  use precision_mod, only: fp
  !
  !  use mkspline to set up spline coefficients...
  !
  !============
  implicit none
  integer nx,iadr,i
  !============
  real(fp) :: xp,xpi,xp2,xpi2,cx,cxi,hx2,cxd,cxdi
  !============
  integer ict(3)                    ! requested output control
  integer ivec                      ! vector length
  integer ivecd                     ! vector dimension (1st dim of fval)
  !
  integer ii(ivec)                  ! target cells (i)
  real(fp) :: xparam(ivec)
  ! normalized displacements from (i) corners
  !
  real(fp) :: hx(ivec)                     ! grid spacing, and
  real(fp) :: hxi(ivec)                    ! inverse grid spacing 1/(x(i+1)-x(i))
  !
  real(fp) :: fin(0:1,nx)                  ! interpolant data (cf "evspline")
  !
  real(fp) :: fval(ivecd,*)                ! output returned
  !
  !  for detailed description of fin, ict and fval see subroutine
  !  evspline comments.  Note ict is not vectorized; the same output
  !  is expected to be returned for all input vector data points.
  !
  !  note that the index inputs ii and parameter inputs
  !     xparam,hx,hxi, are vectorized, and the
  !     output array fval has a vector ** 1st dimension ** whose
  !     size must be given as a separate argument
  !
  !  to use this routine in scalar mode, pass in ivec=ivecd=1
  !
  !---------------
  !  Spline evaluation consists of a "mixing" of the interpolant
  !  data using the linear functionals xparam, xpi = 1-xparam,
  !  xparam**3-xparam, xpi**3-xpi ...
  !  and their derivatives as needed.
  !
  integer v
  !
  real(fp) :: sum
  real(fp) :: sixth
  !
  data sixth/0.166666666666666667_fp/
  !
  !---------------
  !   ...in x direction
  !
  iadr=0
  !
  if(ict(1).le.2) then
     !
     !  normal call
     !
     if(ict(1).eq.1) then
        !
        !  get desired values:
        !
        iadr = iadr + 1
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
           !  function value:
           !
           sum=xpi*fin(0,i) +xp*fin(0,i+1)
           sum=sum+sixth*hx2*(cxi*fin(1,i) +cx*fin(1,i+1))
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
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           xp2=xp*xp
           xpi2=xpi*xpi
           !
           cxd=3.0_fp*xp2-1.0_fp
           cxdi=-3.0_fp*xpi2+1.0_fp

           sum=hxi(v)*(fin(0,i+1)-fin(0,i))
           sum=sum+sixth*hx(v)*(cxdi*fin(1,i) +cxd*fin(1,i+1))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(3).eq.1) then
        !
        !  d2f/dx2:
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           !
           sum=xpi*fin(1,i) +xp*fin(1,i+1)
           fval(v,iadr)=sum
        end do
     end if
     !
  else
     !
     !  return fxxx = d3f/dx3
     !
     iadr=1
     do v=1,ivec
        i=ii(v)
        fval(v,iadr)=hxi(v)*(fin(1,i+1)-fin(1,i))
     end do
     !
  end if
  !
  return
end subroutine fvspline
