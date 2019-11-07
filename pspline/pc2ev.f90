subroutine pc2ev(xget,yget,x,nx,y,ny,ilinx,iliny, &
     f,inf2,ict,fval,ier)
  use psp_precision_mod, only: fp
  !
  !  evaluate a piecewise bilinear interpolant on a rectilinear
  !  grid -- this is C1 in both directions.
  !
  !  this subroutine calls two subroutines:
  !     herm2xy  -- find cell containing (xget,yget)
  !     pc2fcn -- evaluate interpolant function and (optionally) derivatives
  !
  !  input arguments:
  !  ================
  !
  !============
  implicit none
  integer ny,inf2,nx
  !============
  real(fp) :: xget,yget                    ! target of this interpolation
  real(fp) :: x(nx)                        ! ordered x grid
  real(fp) :: y(ny)                        ! ordered y grid
  integer ilinx                     ! ilinx=1 => assume x evenly spaced
  integer iliny                     ! iliny=1 => assume y evenly spaced
  !
  real(fp) :: f(inf2,ny)                   ! function data
  !
  !       f 2nd dimension inf2 must be .ge. nx
  !       contents of f:
  !
  !  f(i,j) = f @ x(i),y(j)
  !
  integer ict(4)                    ! code specifying output desired
  !
  !  ict(1)=1 -- return f  (0, don't)
  !  ict(2)=1 -- return df/dx  (0, don't)
  !  ict(3)=1 -- return df/dy  (0, don't)
  !  ict(4)=1 -- return d2f/dxdy (0, don't)
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
  !
  !  examples:
  !    on input ict = [1,1,1,1]
  !   on output fval= [f,df/dx,df/dy,d2f/dxdy]
  !
  !    on input ict = [1,0,0,0]
  !   on output fval= [f] ... elements 2 & 3 & 4 never referenced
  !
  !    on input ict = [0,1,1,0]
  !   on output fval= [df/dx,df/dy] ... element 3 & 4 never referenced
  !
  !    on input ict = [0,0,1,0]
  !   on output fval= [df/dy] ... elements 2 & 3 & 4 never referenced.
  !
  !  ier -- completion code:  0 means OK
  !-------------------
  !  local:
  !
  integer, dimension(1) :: i,j                       ! cell indices
  !
  !  normalized displacement from (x(i),y(j)) corner of cell.
  !    xparam=0 @x(i)  xparam=1 @x(i+1)
  !    yparam=0 @y(j)  yparam=1 @y(j+1)
  !
  real(fp), dimension(1) :: xparam,yparam
  !
  !  cell dimensions and
  !  inverse cell dimensions hxi = 1/(x(i+1)-x(i)), hyi = 1/(y(j+1)-y(j))
  !
  real(fp), dimension(1) :: hx,hy
  real(fp), dimension(1) :: hxi,hyi
  !
  !  0 .le. xparam .le. 1
  !  0 .le. yparam .le. 1
  !
  !---------------------------------------------------------------------
  !
  call herm2xy(xget,yget,x,nx,y,ny,ilinx,iliny, &
       i(1),j(1),xparam(1),yparam(1), &
       hx(1),hxi(1),hy(1),hyi(1),ier)
  if(ier.ne.0) return
  !
  call pc2fcn(ict,1,1, &
       fval,i,j,xparam,yparam,hx,hxi,hy,hyi,f,inf2,ny)
  !
  return
end subroutine pc2ev
!---------------------------------------------------------------------
!  evaluate piecewise bilinear function interpolation -- 2d fcn
!   --vectorized-- dmc 10 Feb 1999
!
subroutine pc2fcn(ict,ivec,ivecd, &
     fval,ii,jj,xparam,yparam,hx,hxi,hy,hyi, &
     fin,inf2,ny)
  use psp_precision_mod, only: fp
  !
  !============
  implicit none
  integer ny,inf2,i,j,iadr
  !============
  real(fp) :: xp,xpi,yp,ypi
  !============
  integer ict(4)                    ! requested output control
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
  real(fp) :: fin(inf2,ny)                 ! interpolant data (cf "pc2ev")
  !
  real(fp) :: fval(ivecd,*)                ! output returned
  !
  !  for detailed description of fin, ict and fval see subroutine
  !  pc2ev comments.  Note ict is not vectorized; the same output
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
  !
  real(fp) :: sum
  integer v
  !
  !   ...in x direction
  !
  do v=1,ivec
     i=ii(v)
     j=jj(v)
     !
     xp=xparam(v)
     xpi=1.0_fp-xp
     !
     !   ...in y direction
     !
     yp=yparam(v)
     ypi=1.0_fp-yp
     !
     iadr=0
     !
     !  get desired values:
     !
     if(ict(1).eq.1) then
        !
        !  function value:
        !
        iadr=iadr+1
        sum=ypi*(xpi*fin(i,j)+xp*fin(i+1,j)) &
             + yp*(xpi*fin(i,j+1)+xp*fin(i+1,j+1))
        !
        fval(v,iadr)=sum
     end if
     !
     if(ict(2).eq.1) then
        !
        !  df/dx:
        !
        iadr=iadr+1
        !
        sum=ypi*(fin(i+1,j)-fin(i,j)) &
             + yp*(fin(i+1,j+1)-fin(i,j+1))
        fval(v,iadr)=sum*hxi(v)
        !
     end if
     !
     if(ict(3).eq.1) then
        !
        !  df/dy:
        !
        iadr=iadr+1
        !
        sum=xpi*(fin(i,j+1)-fin(i,j)) &
             + xp*(fin(i+1,j+1)-fin(i+1,j))
        fval(v,iadr)=sum*hyi(v)
     end if
     !
     if(ict(4).eq.1) then
        !
        !  d2f/dxdy:
        !
        iadr=iadr+1
        !
        sum=fin(i+1,j+1)-fin(i,j+1)-fin(i+1,j)+fin(i,j)
        fval(v,iadr)=sum*hxi(v)*hyi(v)
     end if
     !
  end do                             ! vector loop
  !
  return
end subroutine pc2fcn
!---------------------------------------------------------------------
!  evaluate piecewise bilinear function interpolation -- 2d fcn
!   --vectorized-- dmc 10 Feb 1999
!    --optimized for VARIATION along x axis ONLY--
!
subroutine pc2fcnx(ict,ivec,ivecd, &
     fval,ii,jj,xparam,yparam,hx,hxi,hy,hyi, &
     fin,inf2,ny)
  use psp_precision_mod, only: fp
  !
  !============
  implicit none
  integer ny,inf2,j,i,iadr
  !============
  real(fp) :: yp,ypi,xp,xpi
  !============
  integer ict(4)                    ! requested output control
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
  real(fp) :: fin(inf2,ny)                 ! interpolant data (cf "pc2ev")
  !
  real(fp) :: fval(ivecd,*)                ! output returned
  !
  !  for detailed description of fin, ict and fval see subroutine
  !  pc2ev comments.  Note ict is not vectorized; the same output
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
  !
  real(fp) :: sum
  integer v
  !
  !   ...in x direction
  !
  j=jj
  !
  !   ...in y direction
  !
  yp=yparam
  ypi=1.0_fp-yp
  !
  do v=1,ivec
     i=ii(v)
     !
     xp=xparam(v)
     xpi=1.0_fp-xp
     !
     iadr=0
     !
     !  get desired values:
     !
     if(ict(1).eq.1) then
        !
        !  function value:
        !
        iadr=iadr+1
        sum=ypi*(xpi*fin(i,j)+xp*fin(i+1,j)) &
             + yp*(xpi*fin(i,j+1)+xp*fin(i+1,j+1))
        !
        fval(v,iadr)=sum
     end if
     !
     if(ict(2).eq.1) then
        !
        !  df/dx:
        !
        iadr=iadr+1
        !
        sum=ypi*(fin(i+1,j)-fin(i,j)) &
             + yp*(fin(i+1,j+1)-fin(i,j+1))
        fval(v,iadr)=sum*hxi(v)
        !
     end if
     !
     if(ict(3).eq.1) then
        !
        !  df/dy:
        !
        iadr=iadr+1
        !
        sum=xpi*(fin(i,j+1)-fin(i,j)) &
             + xp*(fin(i+1,j+1)-fin(i+1,j))
        fval(v,iadr)=sum*hyi
     end if
     !
     if(ict(4).eq.1) then
        !
        !  d2f/dxdy:
        !
        iadr=iadr+1
        !
        sum=fin(i+1,j+1)-fin(i,j+1)-fin(i+1,j)+fin(i,j)
        fval(v,iadr)=sum*hxi(v)*hyi
     end if
     !
  end do                             ! vector loop
  !
  return
end subroutine pc2fcnx
