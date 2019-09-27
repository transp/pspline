subroutine vecherm2(ict,ivec,xvec,yvec,ivd,fval, &
     nx,xpkg,ny,ypkg,fspl,inf2, &
     iwarn,ier)
  use precision_mod, only: fp
  !
  !  vectorized hermite evaluation routine -- 2d
  !  1.  call vectorized zone lookup routine
  !  2.  call vectorized hermite evaluation routine
  !
  !--------------------------
  !  input:
  !============
  implicit none
  integer iwarn1,iwarn2
  !============
  integer ict(4)                    ! selector:
  !        ict(1)=1 for f      (don't evaluate if ict(1)=0)
  !        ict(2)=1 for df/dx  (don't evaluate if ict(2)=0)
  !        ict(3)=1 for df/dy  (don't evaluate if ict(3)=0)
  !        ict(4)=1 for d2f/dxdy (don't evaluate if ict(4)=0)
  !
  integer ivec                      ! vector dimensioning
  !
  !    ivec-- number of vector pts (hermite values to look up)
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
  integer nx,ny                     ! dimension of hermite grids
  real(fp) :: xpkg(nx,4)                   ! x grid "package" (cf genxpkg)
  real(fp) :: ypkg(ny,4)                   ! y grid "package" (cf genxpkg)
  integer inf2                      ! fspl 3rd array dimension, .ge.nx
  real(fp) :: fspl(0:3,inf2,ny)            ! Hermite coefficients, cf herm2ev
  !
  ! output:
  ! condition codes, 0 = normal return
  integer iwarn                     ! =1 if an x value was out of range
  integer ier                       ! =1 if argument error detected
  !
  !---------------------------------------------------------------
  !  local arrays
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
     write(6,*) ' ?vecherm2:  nx.lt.2:  nx = ',nx
     ier=1
  end if
  !
  if(ny.lt.2) then
     write(6,*) ' ?vecherm2:  ny.lt.2:  ny = ',ny
     ier=1
  end if
  !
  if(ivec.le.0) then
     write(6,*) ' ?vecherm2:  vector dimension .le. 0:  ivec = ', &
          ivec
     ier=1
  end if
  !
  if(ivd.lt.ivec) then
     write(6,*) &
          ' ?vecherm2:  output vector dimension less than input ', &
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
  call herm2fcn(ict,ivec,ivd,fval,ix,iy,dxn,dyn, &
       hx,hxi,hy,hyi,fspl,inf2,ny)
  !
  return
end subroutine vecherm2
