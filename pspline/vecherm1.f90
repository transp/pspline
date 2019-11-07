subroutine vecherm1(ict,ivec,xvec,ivd,fval,nx,xpkg,fspl, &
     iwarn,ier)
  use psp_precision_mod, only: fp
  !
  !  vectorized hermite evaluation routine -- 1d
  !  1.  call vectorized zone lookup routine
  !  2.  call vectorized hermite evaluation routine
  !
  !--------------------------
  !  input:
  implicit none
  integer ict(2)                    ! selector:
  !        ict(1)=1 for f      (don't evaluate if ict(1)=0)
  !        ict(2)=1 for df/dx  (don't evaluate if ict(2)=0)
  !
  integer ivec                      ! vector dimensioning
  !
  !    ivec-- number of vector pts (hermite values to look up)
  !
  real(fp) :: xvec(ivec)                   ! x-locations at which to evaluate
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
  integer nx                        ! dimension of hermite x grid
  real(fp) :: xpkg(nx,4)                   ! x grid "package" (cf genxpkg)
  real(fp) :: fspl(2,nx)                   ! Hermite coefficients (cf herm1ev)
  !
  ! output:
  ! condition codes, 0 = normal return
  integer iwarn                     ! =1 if an x value was out of range
  integer ier                       ! =1 if argument error detected
  !
  !---------------------------------------------------------------
  !  local arrays
  !
  integer iv(ivec)                  ! zone indices {j}
  real(fp) :: dxn(ivec)                    ! normalized displacements w/in zones
  real(fp) :: h(ivec)                      ! h(j) vector
  real(fp) :: hi(ivec)                     ! 1/h(j) vector
  !
  !---------------------------------------------------------------
  !
  !  error checks
  !
  ier=0
  if(nx.lt.2) then
     write(6,*) ' ?vecherm1:  nx.lt.2:  nx = ',nx
     ier=1
  end if
  !
  if(ivec.le.0) then
     write(6,*) ' ?vecherm1:  vector dimension .le. 0:  ivec = ', &
          ivec
     ier=1
  end if
  !
  if(ivd.lt.ivec) then
     write(6,*) &
          ' ?vecherm1:  output vector dimension less than input ', &
          'vector dimension.'
     write(6,*) ' ivec=',ivec,' ivd=',ivd
     ier=1
  end if
  !
  if(ier.ne.0) return
  !
  !  vectorized lookup
  !
  iv=0
  call xlookup(ivec,xvec,nx,xpkg,2,iv,dxn,h,hi,iwarn)
  !
  !  vectorized evaluation
  !
  call herm1fcn(ict,ivec,ivd,fval,iv,dxn,h,hi,fspl,nx)
  !
  return
end subroutine vecherm1
