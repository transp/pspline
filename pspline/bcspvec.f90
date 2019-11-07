subroutine bcspvec(ict,ivec,xvec,yvec,ivd,fval, &
     nx,xpkg,ny,ypkg,fspl,inf3, &
     iwarn,ier)
  use psp_precision_mod, only: fp
  !
  !  vectorized spline evaluation routine -- 2d spline
  !  1.  call vectorized zone lookup routine
  !  2.  call vectorized spline evaluation routine
  !
  !--------------------------
  !  input:
  !============
  implicit none
  integer iwarn1,iwarn2
  !============
  integer ict(6)                    ! selector:
  !        ict(1)=1 for f      (don't evaluate if ict(1)=0)
  !        ict(2)=1 for df/dx  (don't evaluate if ict(2)=0)
  !        ict(3)=1 for df/dy  (don't evaluate if ict(3)=0)
  !        ict(4)=1 for d2f/dx2 (don't evaluate if ict(4)=0)
  !        ict(5)=1 for d2f/dy2 (don't evaluate if ict(5)=0)
  !        ict(6)=1 for d2f/dxdy (don't evaluate if ict(6)=0)
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
  integer inf3                      ! fspl 3rd array dimension, .ge.nx
  real(fp) :: fspl(4,4,inf3,ny)            ! (non-compact) spline coefficients
  !
  ! output:
  ! condition codes, 0 = normal return
  integer iwarn                     ! =1 if an x value was out of range
  integer ier                       ! =1 if argument error detected
  !
  !---------------------------------------------------------------
  !  local arrays
  !
  integer ix(ivec)                  ! x zone indices
  real(fp) :: dxv(ivec)                    ! x displacements w/in zones
  integer iy(ivec)                  ! y zone indices
  real(fp) :: dyv(ivec)                    ! y displacements w/in zones
  !
  !---------------------------------------------------------------
  !
  !  error checks
  !
  ier=0
  !
  if(nx.lt.2) then
     write(6,*) ' ?bcspvec:  nx.lt.2:  nx = ',nx
     ier=1
  end if
  !
  if(ny.lt.2) then
     write(6,*) ' ?bcspvec:  ny.lt.2:  ny = ',ny
     ier=1
  end if
  !
  if(ivec.le.0) then
     write(6,*) ' ?bcspvec:  vector dimension .le. 0:  ivec = ', &
          ivec
     ier=1
  end if
  !
  if(ivd.lt.ivec) then
     write(6,*) &
          ' ?bcspvec:  output vector dimension less than input ', &
          'vector dimension.'
     write(6,*) ' ivec=',ivec,' ivd=',ivd
     ier=1
  end if
  !
  if(ier.ne.0) return
  !
  !  vectorized lookups
  !
  ix=0
  iy=0
  call xlookup(ivec,xvec,nx,xpkg,1,ix,dxv,dxv,dxv,iwarn1)
  call xlookup(ivec,yvec,ny,ypkg,1,iy,dyv,dyv,dyv,iwarn2)
  iwarn=max(iwarn1,iwarn2)
  !
  !  vectorized evaluation
  !
  call bcspevfn(ict,ivec,ivd,fval,ix,iy,dxv,dyv,fspl,inf3,ny)
  !
  return
end subroutine bcspvec
