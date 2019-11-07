subroutine spvec(ict,ivec,xvec,ivd,fval,nx,xpkg,fspl,iwarn,ier)
  use psp_precision_mod, only: fp
  !
  !  vectorized spline evaluation routine -- 1d spline
  !  1.  call vectorized zone lookup routine
  !  2.  call vectorized spline evaluation routine
  !
  !--------------------------
  !  input:
  implicit none
  integer ict(3)                    ! selector:
  !        ict(1)=1 for f      (don't evaluate if ict(1)=0)
  !        ict(2)=1 for df/dx  (don't evaluate if ict(2)=0)
  !        ict(3)=1 for d2f/dx2 (don't evaluate if ict(3)=0)
  !
  integer ivec                      ! vector dimensioning
  !
  !    ivec-- number of vector pts (spline values to look up)
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
  integer nx                        ! dimension of spline x grid
  real(fp) :: xpkg(nx,4)                   ! x grid "package" (cf genxpkg)
  real(fp) :: fspl(4,nx)                   ! (non-compact) spline coefficients
  !
  ! output:
  ! condition codes, 0 = normal return
  integer iwarn                     ! =1 if an x value was out of range
  integer ier                       ! =1 if argument error detected
  !
  !---------------------------------------------------------------
  !  local arrays
  !
  integer iv(ivec)                  ! zone indices
  real(fp) :: dxv(ivec)                    ! displacements w/in zones
  !
  !---------------------------------------------------------------
  !
  !  error checks
  !
  ier=0
  if(nx.lt.2) then
     write(6,*) ' ?spvec:  nx.lt.2:  nx = ',nx
     ier=1
  end if
  !
  if(ivec.le.0) then
     write(6,*) ' ?spvec:  vector dimension .le. 0:  ivec = ',ivec
     ier=1
  end if
  !
  if(ivd.lt.ivec) then
     write(6,*) &
          ' ?spvec:  output vector dimension less than input ', &
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
  call xlookup(ivec,xvec,nx,xpkg,1,iv,dxv,dxv,dxv,iwarn)
  !
  !  vectorized evaluation
  !
  call cspevfn(ict,ivec,ivd,fval,iv,dxv,fspl,nx)
  !
  return
end subroutine spvec
