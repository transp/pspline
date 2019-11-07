subroutine akherm2(x,nx,y,ny,fherm,nf2,ilinx,iliny,ier)
  use psp_precision_mod, only: fp
  !
  !  create a data set for Hermite interpolation, based on Akima's method
  !  [Hiroshi Akima, Communications of the ACM, Jan 1974, Vol. 17 No. 1]
  !
  !  input:
  !
  !============
  implicit none
  integer ilinx,iliny,ier,nf2
  !============
  integer nx,ny                            ! array dimensions
  real(fp) :: x(nx)                        ! x coordinate array
  real(fp) :: y(ny)                        ! y coordinate array
  real(fp) :: fherm(0:3,nf2,ny)            ! data/Hermite array
  !
  !  fherm(0,i,j) = function value f at x(i),y(j)  **on input**
  !
  !  fherm(1,i,j) = derivative df/dx at x(i),y(j)  **on output**
  !  fherm(2,i,j) = derivative df/dy at x(i),y(j)  **on output**
  !  fherm(3,i,j) = derivative d2f/dxdy at x(i),y(j)  **on output**
  !
  !  addl output:
  !    ilinx=1 if x axis is evenly spaced
  !    iliny=1 if y axis is evenly spaced
  !    ier=0 if no error:
  !      x, y must both be strict ascending
  !      nf2.ge.nx is required.
  !
  !  a default boundary condition is used, based on divided differences
  !  in the edge region.  For more control of BC, use akherm2p...
  !
  call akherm2p(x,nx,y,ny,fherm,nf2,ilinx,iliny,0,0,ier)
  !
  return
end subroutine akherm2
!----------------------------
subroutine akherm2p(x,nx,y,ny,fherm,nf2,ilinx,iliny,ipx,ipy,ier)
  use psp_precision_mod, only: fp
  !
  !  create a data set for Hermite interpolation, based on Akima's method
  !  [Hiroshi Akima, Communications of the ACM, Jan 1974, Vol. 17 No. 1]
  !
  !  with independently settable boundary condition options:
  !     ipx or ipy  =0  -- default, boundary conditions from divided diffs
  !     ipx or ipy  =1  -- periodic boundary condition
  !     ipx or ipy  =2  -- user supplied df/dx or df/dy
  !  input:
  !
  !============
  implicit none
  integer ilinx,iliny,ier,nf2,ierx,iery,ierbc,iy,ix,icx,icy
  integer incx,incy,ix1,iy1,iym2,iym1,iymm2,iymm1,iyp2,iyp1
  integer iypp2,iypp1,ixm2,ixm1,ixmm2,ixmm1,iflagx,ixp2,ixp1
  integer ixpp2,ixpp1,iflagy
  !============
  real(fp) :: cxp,cxm,cxpp,cxmm,cxtrap0,cxtrapn,cyp,cym,cypp,cymm
  real(fp) :: cytrap0,cytrapn,zansr,cxm2,cxp2
  !============
  integer nx,ny                            ! array dimensions
  real(fp) :: x(nx)                        ! x coordinate array
  real(fp) :: y(ny)                        ! y coordinate array
  real(fp) :: fherm(0:3,nf2,ny)            ! data/Hermite array
  !
  integer ipx                       ! =1 if df/dx periodic in x
  integer ipy                       ! =1 if df/dy periodic in y
  !
  !  fherm(0,1:nx,1:ny) supplied by user; this routine computes the
  !  rest of the elements, but note:
  !
  !    if ipx=2:  fherm(1,1,1:ny) & fherm(1,nx,1:ny) = INPUT df/dx BCs
  !    if ipy=2:  fherm(2,1:nx,1) & fherm(2,1:nx,ny) = INPUT df/dy BCs
  !
  !  on output, at all grid points (i,j) covering (1:nx,1:ny):
  !  fherm1(1,i,j) -- df/dx at the grid point
  !  fherm1(2,i,j) -- df/dy at the grid point
  !  fherm1(3,i,j) -- d2f/dxdy at the grid point
  !
  !---------------------
  !  local...
  !
  real(fp) :: wx(2),wy(2),e(2,2)
  !
  real(fp), dimension(:,:), allocatable ::  ftmp
  !
  real(fp) :: xx(0:nx+1)
  real(fp) :: yy(0:ny+1)
  !
  !---------------------
  !
  !  error checks
  !
  ier=0
  !
  call splinck(x,nx,ilinx,1.0E-3_fp,ierx)
  if(ierx.ne.0) ier=ier+1
  !
  if(ierx.ne.0) then
     write(6,'('' ?akherm2:  x axis not strict ascending'')')
  end if
  !
  call splinck(y,ny,iliny,1.0E-3_fp,iery)
  if(iery.ne.0) ier=ier+1
  !
  if(iery.ne.0) then
     write(6,'('' ?akherm2:  y axis not strict ascending'')')
  end if
  !
  if(nf2.lt.nx) then
     ier=ier+1
     write(6,*) '?akherm2:  fherm array dimension too small.'
  end if
  !
  ierbc=0
  call ibc_ck(ipx,'akherm2','X Bdy Cond',0,2,ierbc)
  ier=ier+ierbc
  !
  ierbc=0
  call ibc_ck(ipy,'akherm2','Y Bdy Cond',0,2,ierbc)
  ier=ier+ierbc
  !
  if(ier.ne.0) return
  !
  !---------------------------------------
  !
  !  get a  temporary array for f -- will extend out by 1 zone in
  !  each direction so that numerical derivative evaluation can be
  !  done without a lot of special case logic near the edges...
  !
  allocate(ftmp(0:nx+1,0:ny+1))
  !
  do iy=1,ny
     do ix=1,nx
        ftmp(ix,iy)=fherm(0,ix,iy)
     end do
  end do
  !
  !  also create expanded axes grids...
  !
  xx(1:nx)=x
  yy(1:ny)=y
  xx(0)=2*x(1)-x(2)
  xx(nx+1)=2*x(nx)-x(nx-1)
  yy(0)=2*y(1)-y(2)
  yy(ny+1)=2*y(ny)-y(ny-1)
  !
  !---------------------------------------
  !
  !  handle boundary conditions and create rows of extrapolated points
  !  in ftmp.  first do ftmp(0,1:ny), ftmp(nx+1,1:ny),
  !                then ftmp(1:nx,0), ftmp(1:nx,ny+1),
  !                then ... fill in the corners
  !
  !  also, for ipx.le.1 fill in the bdy fherm(1,*) values;
  !        for ipy.le.1 fill in the bdy fherm(2,*) values
  !
  !  x bc's
  !
  do iy=1,ny
     !
     cxp=(ftmp(2,iy)-ftmp(1,iy))/(xx(2)-xx(1))
     cxm=(ftmp(nx,iy)-ftmp(nx-1,iy))/(xx(nx)-xx(nx-1))
     !
     if(ipx.eq.1) then
        !
        !  periodic BC
        !
        if(nx.gt.2) then
           cxpp=(ftmp(3,iy)-ftmp(2,iy))/(xx(3)-xx(2))
           cxmm=(ftmp(nx-1,iy)-ftmp(nx-2,iy))/(xx(nx-1)-xx(nx-2))
           !
           call akherm0(cxmm,cxm,cxp,cxpp,wx,fherm(1,1,iy))
           fherm(1,nx,iy)=fherm(1,1,iy)
        else
           fherm(1,1,iy) = cxp  ! =cxm, nx=2
           fherm(1,nx,iy) = fherm(1,1,iy)
        end if
        !
        cxtrap0=cxm
        cxtrapn=cxp
        !
     else if(ipx.eq.0) then
        !
        !  default BC -- standard numeric extrapolation
        !
        if(nx.gt.2) then
           cxpp=(ftmp(3,iy)-ftmp(2,iy))/(xx(3)-xx(2))
           fherm(1,1,iy)=1.5_fp*cxp-0.5_fp*cxpp
           !
           cxmm=(ftmp(nx-1,iy)-ftmp(nx-2,iy))/(xx(nx-1)-xx(nx-2))
           fherm(1,nx,iy)=1.5_fp*cxm-0.5_fp*cxmm
           !
        else
           fherm(1,1,iy) = cxp  ! =cxm, nx=2
           fherm(1,nx,iy) = fherm(1,1,iy)
        end if
        !
        !  extrapolate to slope to ghost points just past bdy...
        !
        cxtrap0=2.0_fp*fherm(1,1,iy)-cxp
        cxtrapn=2.0_fp*fherm(1,nx,iy)-cxm
        !
     else
        !
        !  BC supplied by user.  Also use this for extrapolation...
        !
        cxtrap0=2.0_fp*fherm(1,1,iy)-cxp
        cxtrapn=2.0_fp*fherm(1,nx,iy)-cxm
        !
     end if
     !
     ftmp(0,iy)=ftmp(1,iy)-cxtrap0*(xx(1)-xx(0))
     ftmp(nx+1,iy)=ftmp(nx,iy)+cxtrapn*(xx(nx+1)-xx(nx))
     !
  end do
  !
  !  y bc's
  !
  do ix=1,nx
     !
     cyp=(ftmp(ix,2)-ftmp(ix,1))/(yy(2)-yy(1))
     cym=(ftmp(ix,ny)-ftmp(ix,ny-1))/(yy(ny)-yy(ny-1))
     !
     if(ipy.eq.1) then
        !
        !  periodic BC
        !
        if(ny.gt.2) then
           cypp=(ftmp(ix,3)-ftmp(ix,2))/(yy(3)-yy(2))
           cymm=(ftmp(ix,ny-1)-ftmp(ix,ny-2))/(yy(ny-1)-yy(ny-2))
           !
           call akherm0(cymm,cym,cyp,cypp,wy,fherm(2,ix,1))
           fherm(2,ix,ny)=fherm(2,ix,1)
           !
        else
           fherm(2,ix,1) = cyp  ! =cym, ny=2
           fherm(2,ix,ny)=fherm(2,ix,1)
        end if
        !
        cytrap0=cym
        cytrapn=cyp
        !
     else if(ipy.eq.0) then
        !
        !  default BC -- standard numeric extrapolation
        !
        if(ny.gt.2) then
           cypp=(ftmp(ix,3)-ftmp(ix,2))/(yy(3)-yy(2))
           fherm(2,ix,1)=1.5_fp*cyp-0.5_fp*cypp
           !
           cymm=(ftmp(ix,ny-1)-ftmp(ix,ny-2))/(yy(ny-1)-yy(ny-2))
           fherm(2,ix,ny)=1.5_fp*cym-0.5_fp*cymm
           !
        else
           fherm(2,ix,1) = cyp  ! =cym, ny=2
           fherm(2,ix,ny)=fherm(2,ix,1)
        end if
        !
        !  extrapolate to slope to ghost points just past bdy...
        !
        cytrap0=2.0_fp*fherm(2,ix,1)-cyp
        cytrapn=2.0_fp*fherm(2,ix,ny)-cym
        !
     else
        !
        !  BC supplied by user.  Also use this for extrapolation...
        !
        cytrap0=2.0_fp*fherm(2,ix,1)-cyp
        cytrapn=2.0_fp*fherm(2,ix,ny)-cym
        !
     end if
     !
     ftmp(ix,0)=ftmp(ix,1)-cytrap0*(yy(1)-yy(0))
     ftmp(ix,ny+1)=ftmp(ix,ny)+cytrapn*(yy(ny+1)-yy(ny))
     !
  end do
  !
  !  and do something for the corners...
  !
  do ix=0,1
     do iy=0,1
        icx=ix*(nx+1)
        icy=iy*(ny+1)
        incx=1-2*ix
        incy=1-2*iy

        ix1=icx+incx
        iy1=icy+incy

        ftmp(icx,icy)=ftmp(icx,iy1)+ftmp(ix1,icy)-ftmp(ix1,iy1)
        !
     end do
  end do
  !
  !----------------------------------------------------------------
  !  OK, now ready to compute all the interior coefficients and the
  !  rest of the edge coefficients as well...
  !
  do iy=1,ny
     iym2=iy
     iym1=iym2-1
     !
     iymm2=iy-1
     iymm1=iymm2-1
     !
     iyp2=iy+1
     iyp1=iyp2-1
     !
     iypp2=iy+2
     iypp1=iypp2-1
     !
     do ix=1,nx
        !
        !  x div. diffs in vicinity
        !
        ixm2=ix
        ixm1=ixm2-1
        !
        ixmm2=ix-1
        ixmm1=ixmm2-1
        !
        iflagx=0
        cxm=(ftmp(ixm2,iy)-ftmp(ixm1,iy))/(xx(ixm2)-xx(ixm1))
        if(ix.gt.1) then
           cxmm=(ftmp(ixmm2,iy)-ftmp(ixmm1,iy))/ &
                (xx(ixmm2)-xx(ixmm1))
        else
           if(ipx.eq.1) then
              cxmm=(ftmp(nx-1,iy)-ftmp(nx-2,iy))/ &
                   (xx(nx-1)-xx(nx-2))
           else
              iflagx=1
           end if
        end if
        !
        ixp2=ix+1
        ixp1=ixp2-1
        !
        ixpp2=ix+2
        ixpp1=ixpp2-1
        !
        cxp=(ftmp(ixp2,iy)-ftmp(ixp1,iy))/(xx(ixp2)-xx(ixp1))
        if(ix.lt.nx) then
           cxpp=(ftmp(ixpp2,iy)-ftmp(ixpp1,iy))/ &
                (xx(ixpp2)-xx(ixpp1))
        else
           if(ipx.eq.1) then
              cxpp=(ftmp(3,iy)-ftmp(2,iy))/(xx(3)-xx(2))
           else
              cxpp=cxp+(cxm-cxmm)
           end if
        end if
        !
        if(iflagx.eq.1) then
           cxmm=cxm+(cxp-cxpp)
        end if
        !
        !  Akima weightings + df/dx for interior pts
        !
        call akherm0(cxmm,cxm,cxp,cxpp,wx,zansr)
        if((ix.gt.1).and.(ix.lt.nx)) fherm(1,ix,iy)=zansr
        !
        !  y div. diffs in vicinity
        !
        iflagy=0
        cym=(ftmp(ix,iym2)-ftmp(ix,iym1))/(yy(iym2)-yy(iym1))
        if(iy.gt.1) then
           cymm=(ftmp(ix,iymm2)-ftmp(ix,iymm1))/ &
                (yy(iymm2)-yy(iymm1))
        else
           if(ipy.eq.1) then
              cymm=(ftmp(ix,ny-1)-ftmp(ix,ny-2))/ &
                   (yy(ny-1)-yy(ny-2))
           else
              iflagy=1
           end if
        end if
        !
        cyp=(ftmp(ix,iyp2)-ftmp(ix,iyp1))/(yy(iyp2)-yy(iyp1))
        if(iy.lt.ny) then
           cypp=(ftmp(ix,iypp2)-ftmp(ix,iypp1))/ &
                (yy(iypp2)-yy(iypp1))
        else
           if(ipy.eq.1) then
              cypp=(ftmp(ix,3)-ftmp(ix,2))/(yy(3)-yy(2))
           else
              cypp=cyp+(cym-cymm)
           end if
        end if
        !
        if(iflagy.eq.1) then
           cymm=cym+(cyp-cypp)
        end if
        !
        !  Akima weightings + df/dy for interior pts
        !
        call akherm0(cymm,cym,cyp,cypp,wy,zansr)
        if((iy.gt.1).and.(iy.lt.ny)) fherm(2,ix,iy)=zansr
        !
        !  cross derivatives (2nd order divided differences)
        !
        cxm2=(ftmp(ixm2,iym1)-ftmp(ixm1,iym1))/ &
             (xx(ixm2)-xx(ixm1))
        e(1,1)=(cxm-cxm2)/(yy(iym2)-yy(iym1))
        !
        cxm2=(ftmp(ixm2,iyp2)-ftmp(ixm1,iyp2))/ &
             (xx(ixm2)-xx(ixm1))
        e(1,2)=(cxm2-cxm)/(yy(iyp2)-yy(iyp1))
        !
        cxp2=(ftmp(ixp2,iym1)-ftmp(ixp1,iym1))/ &
             (xx(ixp2)-xx(ixp1))
        e(2,1)=(cxp-cxp2)/(yy(iym2)-yy(iym1))
        !
        cxp2=(ftmp(ixp2,iyp2)-ftmp(ixp1,iyp2))/ &
             (xx(ixp2)-xx(ixp1))
        e(2,2)=(cxp2-cxp)/(yy(iyp2)-yy(iyp1))
        !
        !  the values
        !
        fherm(3,ix,iy)=(wx(1)*(wy(1)*e(1,1)+wy(2)*e(1,2))+ &
             wx(2)*(wy(1)*e(2,1)+wy(2)*e(2,2)))/ &
             ((wx(1)+wx(2))*(wy(1)+wy(2)))
        !
     end do
  end do
  !
  deallocate (ftmp)
  return
end subroutine akherm2p
