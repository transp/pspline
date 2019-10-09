subroutine vecintrp3d(ict,ivec,xvec,yvec,zvec,ivd,fval, &
     nx,xpkg,ny,ypkg,nz,zpkg,jspline,fspl,icoeff,ixdim,iydim,izdim, &
     iwarn,ier)
  use precision_mod, only: fp
  !
  !
  !  vectorized hybrid spline evaluation routine -- 3d
  !  1.  call vectorized zone lookup routine
  !  2.  call vectorized hybrid spline evaluation routine
  !
  !--------------------------
  !  input:
  !============
  implicit none
  !============
  integer ict(10)                    ! selector:
  !        ict(1)=1 for f      (don't evaluate if ict(1)=0)
  !        ict(2)=1 for df/dx  (don't evaluate if ict(2)=0)
  !        ict(3)=1 for df/dy  (don't evaluate if ict(3)=0)
  !        ict(4)=1 for df/dy  (don't evaluate if ict(4)=0)
  !        ict(5)=1 for d2f/dx2 (don't evaluate if ict(5)=0)
  !        ict(6)=1 for d2f/dy2 (don't evaluate if ict(6)=0)
  !        ict(7)=1 for d2f/dz2 (don't evaluate if ict(7)=0)
  !        ict(8)=1 for d2f/dxdy (don't evaluate if ict(8)=0)
  !        ict(9)=1 for d2f/dxdz (don't evaluate if ict(9)=0)
  !        ict(10)=1 for d2f/dydz (don't evaluate if ict(10)=0)
  !
  !     higher derivatives can be selected by setting ict(1)=3,-3,4,-4,...
  !     see fvtricub.f90 comments.
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
  !  list of (x,y,z) triples:
  !
  real(fp) :: xvec(ivec)                   ! x-locations at which to evaluate
  real(fp) :: yvec(ivec)                   ! y-locations at which to evaluate
  real(fp) :: zvec(ivec)                   ! z-locations at which to evaluate
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
  integer nx,ny,nz                  ! dimension of spline grids
  real(fp) :: xpkg(nx,4)                   ! x grid "package" (cf genxpkg)
  real(fp) :: ypkg(ny,4)                   ! y grid "package" (cf genxpkg)
  real(fp) :: zpkg(nz,4)                   ! z grid "package" (cf genxpkg)
  !
  integer :: jspline(3)             ! interpolation method on each dim:
  !        jspline(1) for x; jspline(2) for y; jspline(3) for z
  !        -1: zonal step fcn; 0: pcwise linear; 1: Hermite; 2: compact spline
  !
  integer icoeff                    ! #coeffs per data node
  integer ixdim                     ! x dim for fspl
  integer iydim                     ! y dim for fspl
  integer izdim                     ! z dim for fspl
  real(fp) :: fspl(icoeff,ixdim,iydim,izdim)     ! hybrid spline coefficients
  !
  !  ixdim=nx-1 for zonal step function along x dimension; o.w. expect ixdim=nx
  !  iydim=ny-1 for zonal step function along y dimension; o.w. expect iydim=ny
  !  izdim=nz-1 for zonal step function along z dimension; o.w. expect izdim=nz
  !
  ! output:
  ! condition codes, 0 = normal return
  integer iwarn                     ! =1 if an x value was out of range
  integer ier                       ! =1 if argument error detected
  !
  !---------------------------------------------------------------
  !  local arrays
  !
  integer :: iwarn1,iwarn2,iwarn3
  !
  integer, dimension(:), allocatable :: ix,iy,iz
  real(fp), dimension(:), allocatable :: dxn,dyn,dzn
  real(fp), dimension(:), allocatable :: hx,hxi,hy,hyi,hz,hzi
  !
  !---------------------------------------------------------------
  !
  !  error checks
  !
  ier=0
  !
  if(nx.lt.2) then
     write(6,*) ' ?vecintrp3d:  nx.lt.2:  nx = ',nx
     ier=1
  end if
  !
  if(ny.lt.2) then
     write(6,*) ' ?vecintrp3d:  ny.lt.2:  ny = ',ny
     ier=1
  end if
  !
  if(nz.lt.2) then
     write(6,*) ' ?vecintrp3d:  nz.lt.2:  nz = ',nz
     ier=1
  end if
  !
  call vecin3d_argchk('vecintrp3d',jspline, &
       icoeff,nx,ny,nz,ixdim,iydim,izdim,ier)
  !
  if(ivec.le.0) then
     write(6,*) ' ?vecintrp3d:  vector dimension .le. 0:  ivec = ', &
          ivec
     ier=1
  end if
  !
  if(ivd.lt.ivec) then
     write(6,*) &
          ' ?vecintrp3d:  output vector dimension less than input ', &
          'vector dimension.'
     write(6,*) ' ivec=',ivec,' ivd=',ivd
     ier=1
  end if
  !
  if(ier.ne.0) return
  !
  allocate(ix(ivec), iy(ivec), iz(ivec), &
       dxn(ivec), dyn(ivec), dzn(ivec), &
       hx(ivec),  hy(ivec),  hz(ivec), &
       hxi(ivec), hyi(ivec), hzi(ivec), stat=ier)
  !
  if(ier.ne.0) then
     write(6,*) &
          ' ?vecintrp3d: memory allocation failure.'
     ier=99
  end if
  !
  if(ier.ne.0) return
  !
  !  vectorized lookup
  !
  ix=0; iy=0; iz=0
  call xlookup(ivec,xvec,nx,xpkg,2,ix,dxn,hx,hxi,iwarn1)
  call xlookup(ivec,yvec,ny,ypkg,2,iy,dyn,hy,hyi,iwarn2)
  call xlookup(ivec,zvec,nz,zpkg,2,iz,dzn,hz,hzi,iwarn3)
  iwarn=max(iwarn1,iwarn2,iwarn3)
  !
  !  vectorized evaluation
  !
  call fvintrp3d(ict,ivec,ivd,fval,ix,iy,iz,dxn,dyn,dzn, &
       hx,hxi,hy,hyi,hz,hzi,jspline,fspl,icoeff,ixdim,iydim,izdim)
  !
  deallocate(ix,iy,iz,dxn,dyn,dzn,hx,hy,hz,hxi,hyi,hzi)
  !
  return
end subroutine vecintrp3d
