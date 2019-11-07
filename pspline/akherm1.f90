subroutine akherm1(x,nx,fherm,ilinx,ier)
  use psp_precision_mod, only: fp
  !
  !  create a data set for Hermite interpolation, based on Akima's method
  !  [Hiroshi Akima, Communications of the ACM, Jan 1974, Vol. 17 No. 1]
  !
  !  1d routine -- default boundary condition (based on divided differences)
  !
  !  input:
  !
  !============
  implicit none
  integer ier,ilinx
  !============
  integer nx                        ! array dimensions
  real(fp) :: x(nx)                        ! x coordinate array
  real(fp) :: fherm(0:1,nx)                ! data/Hermite array
  !
  !  fherm(0,i) = function value f at x(i)       **on input**
  !
  !  fherm(1,i) = derivative df/dx at x(i)       **on output**
  !
  ! addl output:
  !  ilinx=1 if x is "evenly spaced" ier=0 if no errors
  !
  !  ** x must be strict ascending **
  !
  ! work by calling akherm1p; no periodic boundary condition
  !
  call akherm1p(x,nx,fherm,ilinx,0,ier)
  !
  return
end subroutine akherm1
!----------------------------
subroutine akherm1p(x,nx,fherm,ilinx,ipx,ier)
  use psp_precision_mod, only: fp
  !
  !  Akima subroutine with boundary condition option:
  !
  !  =>ipx=1 for periodic boundary condition
  !
  !  =>ipx=2 for user specified boundary condition:
  !    in which case fherm(1,1) and fherm(1,nx) must contain
  !    the derivatives fx(x(1)) and fx(x(nx)), respectively, and these
  !    will not be modified on output.
  !
  !  =>ipx=0: default boundary conditions
  !
  !  other arguments as with akherm1.
  !
  !  input:
  !
  !============
  implicit none
  integer ier,ilinx,ierbc,ix,ixm2,ixm1,ixmm2,ixmm1,ixp2,ixp1
  integer ixpp2,ixpp1
  !============
  real(fp) :: ztol,cxp,cxm,cxpp,cxmm,cxtrap0,cxtrap1
  !============
  integer nx                        ! array dimensions
  real(fp) :: x(nx)                        ! x coordinate array
  real(fp) :: fherm(0:1,nx)                ! data/Hermite array
  integer ipx                       ! =1:  f' periodic in x
  !
  !----------------------------
  !
  real(fp) :: wx(2)
  !
  !  error checks...
  !
  ztol=1.0E-3_fp
  ier=0
  !
  call splinck(x,nx,ilinx,ztol,ier)
  if(ier.ne.0) then
     write(6,*) '?akherm1:  x axis not strict ascending.'
     ier=ier+1
  end if
  !
  ierbc=0
  call ibc_ck(ipx,'akherm1','Bdy Cond',0,2,ierbc)
  ier=ier+ierbc
  if(ier.gt.0) return
  !
  !  deal with boundary region.  The boundary derivatives are set
  !  and "ghost" points are provided...
  !
  !  all paths need the 1st div. diffs at the bdy...
  !
  cxp=(fherm(0,2)-fherm(0,1))/(x(2)-x(1))
  cxm=(fherm(0,nx)-fherm(0,nx-1))/(x(nx)-x(nx-1))
  !
  if(ipx.eq.1) then
     !
     !  periodic BC
     !
     !  LHS/RHS
     !
     if(nx.gt.2) then
        cxpp=(fherm(0,3)-fherm(0,2))/(x(3)-x(2))
        cxmm=(fherm(0,nx-1)-fherm(0,nx-2))/(x(nx-1)-x(nx-2))
        !
        call akherm0(cxmm,cxm,cxp,cxpp,wx,fherm(1,1))
        fherm(1,nx)=fherm(1,1)
     else
        ! nx=2
        fherm(1,1)=cxp  ! =cxm
        fherm(1,nx)=cxm ! =cxp
     end if
     !
     cxtrap0=cxm
     cxtrap1=cxp
     !
  else if(ipx.eq.0) then
     !
     !  default BC -- standard numeric extrapolation
     !
     if(nx.gt.2) then
        cxpp=(fherm(0,3)-fherm(0,2))/(x(3)-x(2))
        fherm(1,1)=1.5_fp*cxp-0.5_fp*cxpp
        !
        cxmm=(fherm(0,nx-1)-fherm(0,nx-2))/(x(nx-1)-x(nx-2))
        fherm(1,nx)=1.5_fp*cxm-0.5_fp*cxmm
        !
     else
        ! nx=2
        fherm(1,1)=cxp  ! =cxm
        fherm(1,nx)=cxm ! =cxp
     end if
     !
     !  extrapolate to slope to ghost points just past bdy...
     !
     cxtrap0=2.0_fp*fherm(1,1)-cxp
     cxtrap1=2.0_fp*fherm(1,nx)-cxm
     !
  else
     !
     !  BC supplied by user.  Also use this for extrapolation...
     !  extrapolate to slope to ghost points just past bdy...
     !
     cxtrap0=2.0_fp*fherm(1,1)-cxp
     cxtrap1=2.0_fp*fherm(1,nx)-cxm
     !
  end if
  !
  ! NOTE: this loop is inactive if nx=2
  !
  do ix=2,nx-1
     !
     !  x div. diffs in vicinity
     !
     ixm2=ix
     ixm1=ixm2-1
     ixmm2=ixm1
     ixmm1=ixm1-1
     !
     ixp2=ix+1
     ixp1=ix
     ixpp2=ix+2
     ixpp1=ixp2
     !
     if(ix.eq.2) then
        cxmm=cxtrap0
     else
        cxmm=(fherm(0,ixmm2)-fherm(0,ixmm1))/(x(ixmm2)-x(ixmm1))
     end if
     !
     if(ix.eq.nx-1) then
        cxpp=cxtrap1
     else
        cxpp=(fherm(0,ixpp2)-fherm(0,ixpp1))/(x(ixpp2)-x(ixpp1))
     end if
     !
     cxm=(fherm(0,ixm2)-fherm(0,ixm1))/(x(ixm2)-x(ixm1))
     cxp=(fherm(0,ixp2)-fherm(0,ixp1))/(x(ixp2)-x(ixp1))
     !
     call akherm0(cxmm,cxm,cxp,cxpp,wx,fherm(1,ix))
     !
  end do
  !
  return
end subroutine akherm1p
!--------------------------------------
subroutine akherm0(cxmm,cxm,cxp,cxpp,wx,slp)
  use psp_precision_mod, only: fp
  !
  !  basic akima formula for 1st derivative at pt p:
  !
  !     cxmm = numerical slope 2 zones to left
  !     cxm  = numerical slope 1 zone to left
  !     cxp  = numerical slope 1 zone to right
  !     cxpp = numerical slope 2 zones to right
  !
  !  return slope at point p and weighting factors
  !
  implicit none
  real(fp) :: cxmm,cxm,cxp,cxpp            ! nearby slopes (in)
  real(fp) :: wx(2)                        ! weights (out)
  real(fp) :: slp                          ! Akima nominal slope (out)
  !
  wx(2)=abs(cxm-cxmm)
  wx(1)=abs(cxp-cxpp)
  if(wx(1)+wx(2).eq.0.0_fp) then
     wx(1)=1.0_fp
     wx(2)=1.0_fp
  end if
  !
  !  the derivative -- weighted average of neighbouring slopes
  !
  slp=(wx(1)*cxm+wx(2)*cxp)/(wx(1)+wx(2))
  !
  return
end subroutine akherm0
