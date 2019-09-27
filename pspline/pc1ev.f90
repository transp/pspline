subroutine pc1ev(xget,x,nx,ilinx,f,ict,fval,ier)
  use precision_mod, only: fp
  !
  !  evaluate a 1d piecewise linear interpolant -- this is C0.
  !    ...a derivative can be evaluated but it is not continuous.
  !
  !  this subroutine calls two subroutines:
  !     herm1x   -- find cell containing (xget,yget)
  !     pc1fcn -- evaluate interpolant function and (optionally) derivatives
  !
  !  input arguments:
  !  ================
  !
  !============
  implicit none
  integer nx
  !============
  real(fp) :: xget                         ! target of this interpolation
  real(fp) :: x(nx)                        ! ordered x grid
  integer ilinx                     ! ilinx=1 => assume x evenly spaced
  !
  real(fp) :: f(nx)                        ! function data
  !
  !       contents of f:
  !
  !  f(i) = f @ x(i)
  !
  integer ict(2)                    ! code specifying output desired
  !
  !  ict(1)=1 -- return f  (0, don't)
  !  ict(2)=1 -- return df/dx  (0, don't)
  !
  ! output arguments:
  ! =================
  !
  real(fp) :: fval(*)                      ! output data
  integer ier                       ! error code =0 ==> no error
  !
  !  fval(1) receives the first output (depends on ict(...) spec)
  !  fval(2) receives the second output (depends on ict(...) spec)
  !
  !  examples:
  !    on input ict = [1,1]
  !   on output fval= [f,df/dx]
  !
  !    on input ict = [1,0]
  !   on output fval= [f] ... element 2 never referenced
  !
  !    on input ict = [0,1]
  !   on output fval= [df/dx] ... element 2 never referenced
  !
  !  ier -- completion code:  0 means OK
  !-------------------
  !  local:
  !
  integer, dimension(1) :: i=0                         ! cell index
  !
  !  normalized displacement from (x(i)) corner of cell.
  !    xparam=0 @x(i)  xparam=1 @x(i+1)
  !
  real(fp), dimension(1) :: xparam
  !
  !  cell dimensions and
  !  inverse cell dimensions hxi = 1/(x(i+1)-x(i))
  !
  real(fp), dimension(1) :: hx
  real(fp), dimension(1) :: hxi
  !
  !  0 .le. xparam .le. 1
  !
  !---------------------------------------------------------------------
  !
  call herm1x(xget,x,nx,ilinx,i(1),xparam(1),hx(1),hxi(1),ier)
  if(ier.ne.0) return
  !
  call pc1fcn(ict,1,1,fval,i,xparam,hx,hxi,f,nx)
  !
  return
end subroutine pc1ev
!---------------------------------------------------------------------
!  evaluate C0 piecewise linear function interpolation -- 1d fcn
!   --vectorized-- dmc 10 Feb 1999
!
subroutine pc1fcn(ict,ivec,ivecd, &
     fval,ii,xparam,hx,hxi,fin,nx)
  use precision_mod, only: fp
  !
  !  input:
  !
  !============
  implicit none
  integer nx,i,iadr
  !============
  real(fp) :: xp,xpi
  !============
  integer ict(2)                    ! requested output control
  integer ivec                      ! vector length
  integer ivecd                     ! vector dimension (1st dim of fval)
  !
  integer ii(ivec)                  ! target cells
  real(fp) :: xparam(ivec)
  ! normalized displacements in cells
  !
  real(fp) :: hx(ivec)                     ! grid spacing, and
  real(fp) :: hxi(ivec)                    ! inverse grid spacing 1/(x(i+1)-x(i))
  !
  real(fp) :: fin(nx)                      ! the data
  !
  !  output:
  !
  real(fp) :: fval(ivecd,*)                ! interpolation results
  !
  !  for detailed description of fin, ict and fval see subroutine
  !  pc1ev comments.  Note ict is not vectorized -- the same output
  !  is expected to be returned for all input vector data points.
  !
  !  note that the index inputs ii, and parameter inputs
  !     xparam,hx,hxi,are vectorized, and the
  !
  !     output array fval has a vector ** 1st dimension ** whose
  !     size must be given as a separate argument; ivecd.ge.ivec
  !     expected!
  !
  !  to use this routine in scalar mode, pass in ivec=ivecd=1
  !
  !---------------
  !
  integer v
  !
  do v=1,ivec
     i=ii(v)
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
        fval(v,iadr)=xpi*fin(i)+xp*fin(i+1)
        !
     end if
     !
     if(ict(2).eq.1) then
        !
        !  df/dx:
        !
        iadr=iadr+1
        !
        fval(v,iadr)=(fin(i+1)-fin(i))*hxi(v)
     end if
     !
  end do                             ! vector loop
  !
  return
end subroutine pc1fcn
