subroutine vecintrp2d(ict,ivec,xvec,yvec,ivd,fval, &
     nx,xpkg,ny,ypkg,jspline,fspl,icoeff,ixdim,iydim, &
     iwarn,ier)
  use precision_mod, only: fp
  !
  !
  !  vectorized hybrid spline evaluation routine -- 2d
  !  1.  call vectorized zone lookup routine
  !  2.  call vectorized hybrid spline evaluation routine
  !
  !--------------------------
  !  input:
  implicit none
  integer ict(6)                    ! selector:
  !        ict(1)=1 for f      (don't evaluate if ict(1)=0)
  !        ict(2)=1 for df/dx  (don't evaluate if ict(2)=0)
  !        ict(3)=1 for df/dy  (don't evaluate if ict(3)=0)
  !        ict(4)=1 for d2f/dx2 (don't evaluate if ict(4)=0)
  !        ict(5)=1 for d2f/dy2 (don't evaluate if ict(5)=0)
  !        ict(6)=1 for d2f/dxdy (don't evaluate if ict(6)=0)
  !
  !        ict(1)=3 followed by ict(2:5) = 1 or 0 allow 3rd derivatives to
  !          be selected:  fxxx fxxy fxyy fyyy
  !
  !        ict(1)=4 followed by ict(2:4) allow 4th derivatives to
  !          be selected:  fxxxy fxxyy fxyyy; fxxxx=fyyyy=0
  !
  !        ict(1)=5 followed by ict(2:4) allow 4th derivatives to
  !          be selected:  fxxxyy fxxyyy
  !
  !        ict(1)=6 specifies 6th derivative:  fxxxyyy (step fcn)
  !
  !     in hybrid spline evaluation, any derivative along a dimension
  !     with zonal interpolation (jspline(i)=-1) gives zero;
  !
  !     piecewise linear and Hermite interpolation give zero if a 2nd or
  !     higher derivative is sought along any dimension.
  !
  integer ivec                      ! vector dimensioning
  !
  !    ivec-- number of vector pts (spline values to look up)
  !
  !  list of (x,y) pairs:
  !
  real(fp) :: xvec(ivec)                   ! x-locations at which to evaluate
  real(fp) :: yvec(ivec)                   ! y-locations at which to evaluate
  !
  integer ivd                       ! 1st dimension of output array
  !
  !    ivd -- 1st dimension of fval, .ge.ivec
  !
  ! output:
  real(fp) :: fval(ivd,*)                  ! output array
  !
  !  fval(1:ivec,1) -- values as per 1st non-zero ict(...) element
  !  fval(1:ivec,2) -- values as per 2nd non-zero ict(...) element
  !   --etc--
  !
  ! input:
  integer nx,ny                     ! dimension of spline grids
  real(fp) :: xpkg(nx,4)                   ! x grid "package" (cf genxpkg)
  real(fp) :: ypkg(ny,4)                   ! y grid "package" (cf genxpkg)
  !
  integer :: jspline(2)             ! interpolation method on each dim:
  !        jspline(1) for x; jspline(2) for y
  !        -1: zonal step fcn; 0: pcwise linear; 1: Hermite; 2: compact spline
  !
  integer icoeff                    ! #coeffs per data node
  integer ixdim                     ! x dim for fspl
  integer iydim                     ! y dim for fspl
  real(fp) :: fspl(icoeff,ixdim,iydim)     ! hybrid spline coefficients
  !
  !  ixdim=nx-1 for zonal step function along x dimension; o.w. expect ixdim=nx
  !  iydim=ny-1 for zonal step function along y dimension; o.w. expect iydim=ny
  !
  ! output:
  ! condition codes, 0 = normal return
  integer iwarn                     ! =1 if an x value was out of range
  integer ier                       ! =1 if argument error detected
  !
  !---------------------------------------------------------------
  !  local variables and arrays
  !
  integer :: iwarn1,iwarn2
  !
  integer ix(ivec)                  ! zone indices {j}
  real(fp) :: dxn(ivec)                    ! normalized displacements w/in zones
  real(fp) :: hx(ivec)                     ! h(j) vector
  real(fp) :: hxi(ivec)                    ! 1/h(j) vector
  !
  integer iy(ivec)                  ! zone indices {j}
  real(fp) :: dyn(ivec)                    ! normalized displacements w/in zones
  real(fp) :: hy(ivec)                     ! h(j) vector
  real(fp) :: hyi(ivec)                    ! 1/h(j) vector
  !
  !---------------------------------------------------------------
  !
  !  error checks
  !
  ier=0
  !
  if(nx.lt.2) then
     write(6,*) ' ?vecintrp2d:  nx.lt.2:  nx = ',nx
     ier=1
  end if
  !
  if(ny.lt.2) then
     write(6,*) ' ?vecintrp2d:  ny.lt.2:  ny = ',ny
     ier=1
  end if
  !
  call vecin2d_argchk('vecintrp2d',jspline, &
       icoeff,nx,ny,ixdim,iydim,ier)

  if(ivec.le.0) then
     write(6,*) ' ?vecintrp2d:  vector dimension .le. 0:  ivec = ', &
          ivec
     ier=1
  end if
  !
  if(ivd.lt.ivec) then
     write(6,*) &
          ' ?vecintrp2d:  output vector dimension less than input ', &
          'vector dimension.'
     write(6,*) ' ivec=',ivec,' ivd=',ivd
     ier=1
  end if
  !
  if(ier.ne.0) return
  !
  !  vectorized lookup
  !
  ix=0
  iy=0
  call xlookup(ivec,xvec,nx,xpkg,2,ix,dxn,hx,hxi,iwarn1)
  call xlookup(ivec,yvec,ny,ypkg,2,iy,dyn,hy,hyi,iwarn2)
  iwarn=iwarn1+iwarn2
  !
  !  vectorized evaluation
  !
  call fvintrp2d(ict,ivec,ivd,fval,ix,iy,dxn,dyn, &
       hx,hxi,hy,hyi,jspline,fspl,icoeff,ixdim,iydim)
  !
  return
end subroutine vecintrp2d
