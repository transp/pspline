subroutine vectricub(ict,ivec,xvec,yvec,zvec,ivd,fval, &
     nx,xpkg,ny,ypkg,nz,zpkg,fspl,inf4,inf5, &
     iwarn,ier)
  use precision_mod, only: fp
  !
  !  vectorized spline evaluation routine -- 3d *compact* spline
  !  1.  call vectorized zone lookup routine
  !  2.  call vectorized spline evaluation routine
  !
  !--------------------------
  !  input:
  implicit none
  integer iwarn1,iwarn2,iwarn3
  real(fp) :: stat
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
  integer inf4                      ! fspl 4th array dimension, .ge.nx
  integer inf5                      ! fspl 5th array dimension, .ge.ny
  real(fp) :: fspl(8,inf4,inf5,nz)         ! (compact) spline coefficients
  !
  ! output:
  ! condition codes, 0 = normal return
  integer iwarn                     ! =1 if an x value was out of range
  integer ier                       ! =1 if argument error detected
  !
  !---------------------------------------------------------------
  !  local arrays
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
     write(6,*) ' ?vectricub:  nx.lt.2:  nx = ',nx
     ier=1
  end if
  !
  if(ny.lt.2) then
     write(6,*) ' ?vectricub:  ny.lt.2:  ny = ',ny
     ier=1
  end if
  !
  if(nz.lt.2) then
     write(6,*) ' ?vectricub:  nz.lt.2:  nz = ',nz
     ier=1
  end if
  !
  if(ivec.le.0) then
     write(6,*) ' ?vectricub:  vector dimension .le. 0:  ivec = ', &
          ivec
     ier=1
  end if
  !
  if(ivd.lt.ivec) then
     write(6,*) &
          ' ?vectricub:  output vector dimension less than input ', &
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
          ' ?vectricub: memory allocation failure.'
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
  call fvtricub(ict,ivec,ivd,fval,ix,iy,iz,dxn,dyn,dzn, &
       hx,hxi,hy,hyi,hz,hzi,fspl,inf4,inf5,nz)
  !
  deallocate(ix,iy,iz,dxn,dyn,dzn,hx,hy,hz,hxi,hyi,hzi)
  !
  return
end subroutine vectricub
