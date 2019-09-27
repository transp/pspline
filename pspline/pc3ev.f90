subroutine pc3ev(xget,yget,zget,x,nx,y,ny,z,nz, &
     ilinx,iliny,ilinz, &
     f,inf2,inf3,ict,fval,ier)
  use precision_mod, only: fp
  !
  !  evaluate a trilinear interpolant on a rectilinear grid
  !  derivatives are available, but, not continuous across grid planes.
  !
  !  this subroutine calls two subroutines:
  !     herm3xyz  -- find cell containing (xget,yget,zget)
  !     pc3fcn -- evaluate interpolant function and (optionally) derivatives
  !
  !  input arguments:
  !  ================
  !
  !============
  implicit none
  integer ny,nz,inf2,inf3,nx
  !============
  real(fp) :: xget,yget,zget               ! target of this interpolation
  real(fp) :: x(nx)                        ! ordered x grid
  real(fp) :: y(ny)                        ! ordered y grid
  real(fp) :: z(nz)                        ! ordered z grid
  integer ilinx                     ! ilinx=1 => assume x evenly spaced
  integer iliny                     ! iliny=1 => assume y evenly spaced
  integer ilinz                     ! ilinz=1 => assume z evenly spaced
  !
  real(fp) :: f(inf2,inf3,nz)              ! function data
  !
  !       f 2nd dimension inf2 must be .ge. nx; 3rd dim inf3 .ge. ny
  !       contents of f:
  !
  !  f(i,j,k) = f @ x(i),y(j),z(k)
  !
  integer ict(8)                    ! code specifying output desired
  !
  !  ict(1)=1 -- return f  (0, don't)
  !  ict(2)=1 -- return df/dx  (0, don't)
  !  ict(3)=1 -- return df/dy  (0, don't)
  !  ict(4)=1 -- return df/dz  (0, don't)
  !  ict(5)=1 -- return d2f/dxdy  (0, don't)
  !  ict(6)=1 -- return d2f/dxdz  (0, don't)
  !  ict(7)=1 -- return d2f/dydz  (0, don't)
  !  ict(8)=1 -- return d3f/dxdydz  (0, don't)
  !
  ! output arguments:
  ! =================
  !
  real(fp) :: fval(*)                      ! output data
  integer ier                       ! error code =0 ==> no error
  !
  !  fval(1) receives the first output (depends on ict(...) spec)
  !  fval(2) receives the second output (depends on ict(...) spec)
  !  fval(3) receives the third output (depends on ict(...) spec)
  !  fval(4) receives the 4th output (depends on ict(...) spec)
  !  fval(5-8) receive 5th thru 8th outputs (if required by ict(...) spec)
  !
  !  examples:
  !    on input ict = [1,1,1,1,0,0,0,0]
  !   on output fval= [f,df/dx,df/dy,df/dz]
  !
  !    on input ict = [1,0,0,0,0,0,0,0]
  !   on output fval= [f] ... elements 2-8 never referenced
  !
  !    on input ict = [0,1,1,0,0,0,0,0]
  !   on output fval= [df/dx,df/dy] ... elements 3-8 never referenced
  !
  !    on input ict = [0,0,0,0,1,0,0,0]
  !   on output fval= [d2f/dxdy] ... elements 2-8 never referenced.
  !
  !  ier -- completion code:  0 means OK
  !-------------------
  !  local:
  !
  integer, dimension(1) :: i,j,k  ! cell indices
  !
  !  normalized displacement from (x(i),y(j)) corner of cell.
  !    xparam=0 @x(i)  xparam=1 @x(i+1)
  !    yparam=0 @y(j)  yparam=1 @y(j+1)
  !    zparam=0 @z(k)  zparam=1 @z(k+1)
  !
  real(fp), dimension(1) :: xparam,yparam,zparam
  !
  !  cell dimensions and
  !  inverse cell dimensions hxi = 1/(x(i+1)-x(i)), hyi = 1/(y(j+1)-y(j))
  !
  real(fp), dimension(1) :: hx,hy,hz
  real(fp), dimension(1) :: hxi,hyi,hzi
  !
  !  0 .le. xparam .le. 1
  !  0 .le. yparam .le. 1
  !  0 .le. zparam .le. 1
  !
  !---------------------------------------------------------------------
  !
  call herm3xyz(xget,yget,zget,x,nx,y,ny,z,nz,ilinx,iliny,ilinz, &
       i(1),j(1),k(1),xparam(1),yparam(1),zparam(1), &
       hx(1),hxi(1),hy(1),hyi(1),hz(1),hzi(1),ier)
  if(ier.ne.0) return
  !
  call pc3fcn(ict,1,1, &
       fval,i,j,k,xparam,yparam,zparam, &
       hx,hxi,hy,hyi,hz,hzi, &
       f,inf2,inf3,nz)
  !
  return
end subroutine pc3ev
!---------------------------------------------------------------------
!  evaluate trilinear function interpolation -- 3d fcn
!   --vectorized-- dmc 10 Feb 1999
!
subroutine pc3fcn(ict,ivec,ivecd, &
     fval,ii,jj,kk,xparam,yparam,zparam, &
     hx,hxi,hy,hyi,hz,hzi, &
     fin,inf2,inf3,nz)
  use precision_mod, only: fp
  !
  !============
  implicit none
  integer inf3,nz,inf2,i,j,k,iadr
  !============
  real(fp) :: xp,xpi,yp,ypi,zp,zpi
  !============
  integer ict(8)                    ! requested output control
  integer ivec                      ! vector length
  integer ivecd                     ! vector dimension (1st dim of fval)
  !
  integer ii(ivec),jj(ivec),kk(ivec) ! target cells (i,j,k)
  real(fp) :: xparam(ivec),yparam(ivec),zparam(ivec)
  ! normalized displacements from (i,j,k) corners
  !
  real(fp) :: hx(ivec),hy(ivec),hz(ivec)   ! grid spacing, and
  real(fp) :: hxi(ivec),hyi(ivec),hzi(ivec) ! inverse grid spacing
  ! 1/(x(i+1)-x(i)) & 1/(y(j+1)-y(j)) & 1/(z(k+1)-z(i))
  !
  real(fp) :: fin(inf2,inf3,nz)            ! interpolant data (cf "pc3ev")
  !
  real(fp) :: fval(ivecd,*)                ! output returned
  !
  !  for detailed description of fin, ict and fval see subroutine pc3ev
  !  comments.  Note ict is not vectorized; the same output
  !  is expected to be returned for all input vector data points.
  !
  !  note that the index inputs ii,jj,kk and parameter inputs
  !     xparam,yparam,zparam,hx,hxi,hy,hyi,hz,hzi are vectorized, and the
  !     output array fval has a vector ** 1st dimension ** whose
  !     size must be given as a separate argument
  !
  !  to use this routine in scalar mode, pass in ivec=ivecd=1
  !
  !---------------
  !
  real(fp) :: sum
  integer v
  !
  do v=1,ivec
     i=ii(v)
     j=jj(v)
     k=kk(v)
     !
     !   ...in x direction
     !
     xp=xparam(v)
     xpi=1.0_fp-xp
     !
     !   ...in y direction
     !
     yp=yparam(v)
     ypi=1.0_fp-yp
     !
     !   ...in z direction
     !
     zp=zparam(v)
     zpi=1.0_fp-zp
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
        sum=zpi*( &
             xpi*(ypi*fin(i,j,k)  +yp*fin(i,j+1,k))+ &
             xp*(ypi*fin(i+1,j,k)+yp*fin(i+1,j+1,k))) &
             +zp*( &
             xpi*(ypi*fin(i,j,k+1)  +yp*fin(i,j+1,k+1))+ &
             xp*(ypi*fin(i+1,j,k+1)+yp*fin(i+1,j+1,k+1)))
        !
        fval(v,iadr)=sum
     end if
     !
     if(ict(2).eq.1) then
        !
        !  df/dx:
        !
        iadr=iadr+1
        sum=zpi*( &
             -(ypi*fin(i,j,k)  +yp*fin(i,j+1,k)) &
             +(ypi*fin(i+1,j,k)+yp*fin(i+1,j+1,k))) &
             +zp*( &
             -(ypi*fin(i,j,k+1)  +yp*fin(i,j+1,k+1)) &
             +(ypi*fin(i+1,j,k+1)+yp*fin(i+1,j+1,k+1)))
        !
        fval(v,iadr)=sum*hxi(v)
     end if
     !
     if(ict(3).eq.1) then
        !
        !  df/dy:
        !
        iadr=iadr+1
        sum=zpi*( &
             xpi*(-fin(i,j,k)  +fin(i,j+1,k))+ &
             xp*(-fin(i+1,j,k)+fin(i+1,j+1,k))) &
             +zp*( &
             xpi*(-fin(i,j,k+1)  +fin(i,j+1,k+1))+ &
             xp*(-fin(i+1,j,k+1)+fin(i+1,j+1,k+1)))
        !
        fval(v,iadr)=sum*hyi(v)
     end if
     !
     if(ict(4).eq.1) then
        !
        !  df/dz:
        !
        iadr=iadr+1
        sum=   -( &
             xpi*(ypi*fin(i,j,k)  +yp*fin(i,j+1,k))+ &
             xp*(ypi*fin(i+1,j,k)+yp*fin(i+1,j+1,k))) &
             +( &
             xpi*(ypi*fin(i,j,k+1)  +yp*fin(i,j+1,k+1))+ &
             xp*(ypi*fin(i+1,j,k+1)+yp*fin(i+1,j+1,k+1)))
        !
        fval(v,iadr)=sum*hzi(v)
     end if
     !
     if(ict(5).eq.1) then
        !
        !  d2f/dxdy:
        !
        iadr=iadr+1
        sum=zpi*( &
             -(-fin(i,j,k)  +fin(i,j+1,k)) &
             +(-fin(i+1,j,k)+fin(i+1,j+1,k))) &
             +zp*( &
             -(-fin(i,j,k+1)  +fin(i,j+1,k+1)) &
             +(-fin(i+1,j,k+1)+fin(i+1,j+1,k+1)))
        !
        fval(v,iadr)=sum*hxi(v)*hyi(v)
     end if
     !
     if(ict(6).eq.1) then
        !
        !  d2f/dxdz:
        !
        iadr=iadr+1
        sum=  -( &
             -(ypi*fin(i,j,k)  +yp*fin(i,j+1,k)) &
             +(ypi*fin(i+1,j,k)+yp*fin(i+1,j+1,k))) &
             +( &
             -(ypi*fin(i,j,k+1)  +yp*fin(i,j+1,k+1)) &
             +(ypi*fin(i+1,j,k+1)+yp*fin(i+1,j+1,k+1)))
        !
        fval(v,iadr)=sum*hxi(v)*hzi(v)
     end if
     !
     if(ict(7).eq.1) then
        !
        !  d2f/dydz:
        !
        iadr=iadr+1
        sum=  -( &
             xpi*(-fin(i,j,k)  +fin(i,j+1,k))+ &
             xp*(-fin(i+1,j,k)+fin(i+1,j+1,k))) &
             +( &
             xpi*(-fin(i,j,k+1)  +fin(i,j+1,k+1))+ &
             xp*(-fin(i+1,j,k+1)+fin(i+1,j+1,k+1)))
        !
        fval(v,iadr)=sum*hyi(v)*hzi(v)
     end if
     !
     if(ict(8).eq.1) then
        !
        !  d3f/dxdydz:
        !
        iadr=iadr+1
        sum=  -( &
             -(-fin(i,j,k)  +fin(i,j+1,k)) &
             +(-fin(i+1,j,k)+fin(i+1,j+1,k))) &
             +( &
             -(-fin(i,j,k+1)  +fin(i,j+1,k+1)) &
             +(-fin(i+1,j,k+1)+fin(i+1,j+1,k+1)))
        !
        fval(v,iadr)=sum*hxi(v)*hyi(v)*hzi(v)
     end if
     !
  end do                             ! vector loop
  !
  return
end subroutine pc3fcn
!---------------------------------------------------------------------
!  evaluate trilinear function interpolation -- 3d fcn
!   --vectorized-- dmc 10 Feb 1999
!    --optimized for VARIATION along x axis ONLY--
!
subroutine pc3fcnx(ict,ivec,ivecd, &
     fval,ii,jj,kk,xparam,yparam,zparam, &
     hx,hxi,hy,hyi,hz,hzi, &
     fin,inf2,inf3,nz)
  use precision_mod, only: fp
  !
  !============
  implicit none
  integer inf3,nz,inf2,j,k,i,iadr
  !============
  real(fp) :: yp,ypi,zp,zpi,xp,xpi
  !============
  integer ict(8)                    ! requested output control
  integer ivec                      ! vector length
  integer ivecd                     ! vector dimension (1st dim of fval)
  !
  integer ii(ivec),jj,kk ! target cells (i,j,k)
  real(fp) :: xparam(ivec),yparam,zparam
  ! normalized displacements from (i,j,k) corners
  !
  real(fp) :: hx(ivec),hy,hz               ! grid spacing, and
  real(fp) :: hxi(ivec),hyi,hzi            ! inverse grid spacing
  ! 1/(x(i+1)-x(i)) & 1/(y(j+1)-y(j)) & 1/(z(k+1)-z(i))
  !
  real(fp) :: fin(inf2,inf3,nz)            ! interpolant data (cf "pc3ev")
  !
  real(fp) :: fval(ivecd,*)                ! output returned
  !
  !  for detailed description of fin, ict and fval see subroutine pc3ev
  !  comments.  Note ict is not vectorized; the same output
  !  is expected to be returned for all input vector data points.
  !
  !  note that the index inputs ii,jj,kk and parameter inputs
  !     xparam,yparam,zparam,hx,hxi,hy,hyi,hz,hzi are vectorized, and the
  !     output array fval has a vector ** 1st dimension ** whose
  !     size must be given as a separate argument
  !
  !  to use this routine in scalar mode, pass in ivec=ivecd=1
  !
  !---------------
  !
  real(fp) :: sum
  integer v
  !
  j=jj
  k=kk
  !
  !   ...in y direction
  !
  yp=yparam
  ypi=1.0_fp-yp
  !
  !   ...in z direction
  !
  zp=zparam
  zpi=1.0_fp-zp

  do v=1,ivec
     i=ii(v)
     !
     !   ...in x direction
     !
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
        sum=zpi*( &
             xpi*(ypi*fin(i,j,k)  +yp*fin(i,j+1,k))+ &
             xp*(ypi*fin(i+1,j,k)+yp*fin(i+1,j+1,k))) &
             +zp*( &
             xpi*(ypi*fin(i,j,k+1)  +yp*fin(i,j+1,k+1))+ &
             xp*(ypi*fin(i+1,j,k+1)+yp*fin(i+1,j+1,k+1)))
        !
        fval(v,iadr)=sum
     end if
     !
     if(ict(2).eq.1) then
        !
        !  df/dx:
        !
        iadr=iadr+1
        sum=zpi*( &
             -(ypi*fin(i,j,k)  +yp*fin(i,j+1,k)) &
             +(ypi*fin(i+1,j,k)+yp*fin(i+1,j+1,k))) &
             +zp*( &
             -(ypi*fin(i,j,k+1)  +yp*fin(i,j+1,k+1)) &
             +(ypi*fin(i+1,j,k+1)+yp*fin(i+1,j+1,k+1)))
        !
        fval(v,iadr)=sum*hxi(v)
     end if
     !
     if(ict(3).eq.1) then
        !
        !  df/dy:
        !
        iadr=iadr+1
        sum=zpi*( &
             xpi*(-fin(i,j,k)  +fin(i,j+1,k))+ &
             xp*(-fin(i+1,j,k)+fin(i+1,j+1,k))) &
             +zp*( &
             xpi*(-fin(i,j,k+1)  +fin(i,j+1,k+1))+ &
             xp*(-fin(i+1,j,k+1)+fin(i+1,j+1,k+1)))
        !
        fval(v,iadr)=sum*hyi
     end if
     !
     if(ict(4).eq.1) then
        !
        !  df/dz:
        !
        iadr=iadr+1
        sum=   -( &
             xpi*(ypi*fin(i,j,k)  +yp*fin(i,j+1,k))+ &
             xp*(ypi*fin(i+1,j,k)+yp*fin(i+1,j+1,k))) &
             +( &
             xpi*(ypi*fin(i,j,k+1)  +yp*fin(i,j+1,k+1))+ &
             xp*(ypi*fin(i+1,j,k+1)+yp*fin(i+1,j+1,k+1)))
        !
        fval(v,iadr)=sum*hzi
     end if
     !
     if(ict(5).eq.1) then
        !
        !  d2f/dxdy:
        !
        iadr=iadr+1
        sum=zpi*( &
             -(-fin(i,j,k)  +fin(i,j+1,k)) &
             +(-fin(i+1,j,k)+fin(i+1,j+1,k))) &
             +zp*( &
             -(-fin(i,j,k+1)  +fin(i,j+1,k+1)) &
             +(-fin(i+1,j,k+1)+fin(i+1,j+1,k+1)))
        !
        fval(v,iadr)=sum*hxi(v)*hyi
     end if
     !
     if(ict(6).eq.1) then
        !
        !  d2f/dxdz:
        !
        iadr=iadr+1
        sum=  -( &
             -(ypi*fin(i,j,k)  +yp*fin(i,j+1,k)) &
             +(ypi*fin(i+1,j,k)+yp*fin(i+1,j+1,k))) &
             +( &
             -(ypi*fin(i,j,k+1)  +yp*fin(i,j+1,k+1)) &
             +(ypi*fin(i+1,j,k+1)+yp*fin(i+1,j+1,k+1)))
        !
        fval(v,iadr)=sum*hxi(v)*hzi
     end if
     !
     if(ict(7).eq.1) then
        !
        !  d2f/dydz:
        !
        iadr=iadr+1
        sum=  -( &
             xpi*(-fin(i,j,k)  +fin(i,j+1,k))+ &
             xp*(-fin(i+1,j,k)+fin(i+1,j+1,k))) &
             +( &
             xpi*(-fin(i,j,k+1)  +fin(i,j+1,k+1))+ &
             xp*(-fin(i+1,j,k+1)+fin(i+1,j+1,k+1)))
        !
        fval(v,iadr)=sum*hyi*hzi
     end if
     !
     if(ict(8).eq.1) then
        !
        !  d3f/dxdydz:
        !
        iadr=iadr+1
        sum=  -( &
             -(-fin(i,j,k)  +fin(i,j+1,k)) &
             +(-fin(i+1,j,k)+fin(i+1,j+1,k))) &
             +( &
             -(-fin(i,j,k+1)  +fin(i,j+1,k+1)) &
             +(-fin(i+1,j,k+1)+fin(i+1,j+1,k+1)))
        !
        fval(v,iadr)=sum*hxi(v)*hyi*hzi
     end if
     !
  end do                             ! vector loop
  !
  return
end subroutine pc3fcnx
