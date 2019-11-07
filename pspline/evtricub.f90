subroutine evtricub(xget,yget,zget,x,nx,y,ny,z,nz, &
     ilinx,iliny,ilinz, &
     f,inf2,inf3,ict,fval,ier)
  use psp_precision_mod, only: fp
  !
  !  use mktricub to set up spline coefficients...
  !
  !  evaluate a 3d cubic Spline interpolant on a rectilinear
  !  grid -- this is C2 in all directions.
  !
  !  this subroutine calls two subroutines:
  !     herm3xyz  -- find cell containing (xget,yget,zget)
  !     fvtricub  -- evaluate the spline function (w/derivatives if req.)
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
  real(fp) :: f(0:7,inf2,inf3,nz)          ! function data
  !
  !       f 2nd dimension inf2 must be .ge. nx; 3rd dim inf3 .ge. ny
  !       contents of f:
  !
  !  f(0,i,j,k) = f @ x(i),y(j),z(k)
  !  f(1,i,j,k) = d2f/dx2 @ x(i),y(j),z(k)
  !  f(2,i,j,k) = d2f/dy2 @ x(i),y(j),z(k)
  !  f(3,i,j,k) = d2f/dz2 @ x(i),y(j),z(k)
  !  f(4,i,j,k) = d4f/dx2dy2 @ x(i),y(j),z(k)
  !  f(5,i,j,k) = d4f/dx2dz2 @ x(i),y(j),z(k)
  !  f(6,i,j,k) = d4f/dy2dz2 @ x(i),y(j),z(k)
  !  f(7,i,j,k) = d6f/dx2dy2dz2 @ x(i),y(j),z(k)
  !
  integer ict(10)                   ! code specifying output desired
  !
  !  ict(1)=1 -- return f  (0, don't)
  !  ict(2)=1 -- return df/dx  (0, don't)
  !  ict(3)=1 -- return df/dy  (0, don't)
  !  ict(4)=1 -- return df/dz  (0, don't)
  !  ict(5)=1 -- return d2f/dx2  (0, don't)
  !  ict(6)=1 -- return d2f/dy2  (0, don't)
  !  ict(7)=1 -- return d2f/dz2  (0, don't)
  !  ict(8)=1 -- return d2f/dxdy (0, don't)
  !  ict(9)=1 -- return d2f/dxdz (0, don't)
  !  ict(10)=1-- return d2f/dydz (0, don't)
  !
  !  (new dmc Dec 2005 -- higher derivatives available)
  !    ict(1)=3 --> 3rd derivative, .le.2 diff. in any coordinate
  !      ict(2:8) select: fxxy fxxz fxyy fxyz fxzz fyyz fyzz
  !      ->note ict(1)=3, ict(5)=1 gives fxyz = d3f/dxdydz
  !    ict(1)=-3 --> 3rd derivative, 3 in one coordinate
  !      ict(2:4) select: fxxx fyyy fzzz
  !    ict(1)=4 --> 3rd derivative, .le.2 diff. in any coordinate
  !      ict(2:7) select: fxxyy fxxyz fxxzz fxyyz fxyzz fyyzz
  !    ict(1)=-4 --> 3rd derivative, 3 in one coordinate
  !      ict(2:7) select: fxxxy fxxxz fxyyy fxzzz fyyyz fyzzz
  !    ict(1)=5 --> 3rd derivative, .le.2 diff. in any coordinate
  !      ict(2:4) select: fxxyyz fxxyzz fxyyzz
  !    ict(1)=-5 --> 3rd derivative, 3 in one coordinate
  !      ict(2:10) select:  fxxxyy fxxxyz fxxxzz fxxyyy fxxzzz
  !                         fxyyyz fxyzzz fyyyzz fzzzyy
  !    ict(1)=6 --> 3rd derivative, .le.2 diff. in any coordinate
  !      fxxyyzz
  !    ict(1)=-6 --> 3rd derivative, 3 in one coordinate
  !      ict(2:10) select: fxxxyyy fxxxyyz fxxxyzz fxxxyyz
  !                        fxxyyyz fxxyzzz fxyyyzz fxyyzzz fyyyzzz
  !    ict(1)=-7 --> 7th derivative
  !      ict(2:7) select: fxxxyyyz fxxxyyzz fxxxyzzz
  !                       fxxyyyzz fxxyyzzz fxyyyzzz
  !    ict(1)=-8 --> 8th derivative
  !      ict(2:4) select: fxxxyyyzz fxxxyyzzz fxxyyyzzz
  !    ict(1)=-9 --> 9th derivative:  fxxxyyyzzz
  !
  !
  ! output arguments:
  ! =================
  !
  real(fp) :: fval(*)                     ! output data
  integer ier                       ! error code =0 ==> no error
  !
  !  fval(1) receives the first output (depends on ict(...) spec)
  !  fval(2) receives the second output (depends on ict(...) spec)
  !  fval(3) receives the third output (depends on ict(...) spec)
  !  fval(4) receives the 4th output (depends on ict(...) spec)
  !  fval(5-10) receive 5th thru 10th outputs (if required by ict(...) spec)
  !
  !  examples:
  !    on input ict = [1,1,1,1,0,0,0,0,0,0,0]
  !   on output fval= [f,df/dx,df/dy,df/dz]
  !
  !    on input ict = [1,0,0,0,0,0,0,0,0,0,0]
  !   on output fval= [f] ... elements 2-10 never referenced
  !
  !    on input ict = [0,1,1,0,0,0,0,0,0,0,0]
  !   on output fval= [df/dx,df/dy] ... elements 3-10 never referenced
  !
  !    on input ict = [0,0,0,0,1,0,0,0,0,0,0]
  !   on output fval= [d2f/dx2] ... elements 2-10 never referenced.
  !
  !  ier -- completion code:  0 means OK
  !-------------------
  !  local:
  !
  integer i(1),j(1),k(1)                     ! cell indices
  !
  !  normalized displacement from (x(i),y(j),z(k)) corner of cell.
  !    xparam=0 @x(i)  xparam=1 @x(i+1)
  !    yparam=0 @y(j)  yparam=1 @y(j+1)
  !    zparam=0 @z(k)  zparam=1 @z(k+1)
  !
  real(fp) :: xparam(1),yparam(1),zparam(1)
  !
  !  cell dimensions and
  !  inverse cell dimensions hxi = 1/(x(i+1)-x(i)), hyi = 1/(y(j+1)-y(j))
  !
  real(fp) :: hx(1),hy(1),hz(1)
  real(fp) :: hxi(1),hyi(1),hzi(1)
  !
  !  0 .le. xparam .le. 1
  !  0 .le. yparam .le. 1
  !  0 .le. zparam .le. 1
  !
  !---------------------------------------------------------------------
  !  use lookup routine as in Hermite interpolation
  !
  i(1)=0
  j(1)=0
  k(1)=0
  call herm3xyz(xget,yget,zget,x,nx,y,ny,z,nz,ilinx,iliny,ilinz, &
       i(1),j(1),k(1),xparam(1),yparam(1),zparam(1), &
       hx(1),hxi(1),hy(1),hyi(1),hz(1),hzi(1),ier)
  if(ier.ne.0) return
  !
  call fvtricub(ict,1,1, &
       fval,i,j,k,xparam,yparam,zparam, &
       hx,hxi,hy,hyi,hz,hzi, &
       f,inf2,inf3,nz)
  !
  return
end subroutine evtricub
!---------------------------------------------------------------------
!  evaluate C1 cubic Hermite function interpolation -- 3d fcn
!   --vectorized-- dmc 10 Feb 1999
!
subroutine fvtricub(ict,ivec,ivecd, &
     fval,ii,jj,kk,xparam,yparam,zparam, &
     hx,hxi,hy,hyi,hz,hzi, &
     fin,inf2,inf3,nz)
  use psp_precision_mod, only: fp
  !
  !  use mktricub to set up spline coefficients...
  !
  !============
  implicit none
  integer inf3,nz,inf2,iadr,i,j,k
  !============
  real(fp) :: z36th,z216th,xp,xpi,xp2,xpi2,cx,cxi,hx2,yp,ypi,yp2
  real(fp) :: ypi2,cy,cyi,hy2,zp,zpi,zp2,zpi2,cz,czi,hz2,cxd,cxdi
  real(fp) :: cyd,cydi,czd,czdi
  !============
  integer ict(10)                   ! requested output control
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
  real(fp) :: fin(0:7,inf2,inf3,nz)        ! interpolant data (cf "evtricub")
  !
  real(fp) :: fval(ivecd,*)               ! output returned
  !
  !  for detailed description of fin, ict and fval see subroutine evtricub
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
  integer v
  !
  real(fp) :: sum
  real(fp), parameter :: sixth = 0.166666666666666667_fp
  !
  !---------------
  !
  z36th=sixth*sixth
  z216th=sixth*sixth*sixth
  !
  iadr=0
  if(abs(ict(1)).le.2) then
     !
     !  0, 1st, 2nd derivatives...
     !
     !  get desired values:
     !
     if(ict(1).eq.1) then
        !
        !  function value...
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           j=jj(v)
           k=kk(v)
           !
           !   ...in x direction
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           xp2=xp*xp
           xpi2=xpi*xpi
           !
           cx=xp*(xp2-1.0_fp)
           cxi=xpi*(xpi2-1.0_fp)
           hx2=hx(v)*hx(v)
           !
           !   ...and in y direction
           !
           yp=yparam(v)
           ypi=1.0_fp-yp
           yp2=yp*yp
           ypi2=ypi*ypi
           !
           cy=yp*(yp2-1.0_fp)
           cyi=ypi*(ypi2-1.0_fp)
           hy2=hy(v)*hy(v)
           !
           !   ...and in z direction
           !
           zp=zparam(v)
           zpi=1.0_fp-zp
           zp2=zp*zp
           zpi2=zpi*zpi
           !
           cz=zp*(zp2-1.0_fp)
           czi=zpi*(zpi2-1.0_fp)
           hz2=hz(v)*hz(v)
           !
           sum=( &
                zpi*( &
                xpi*(ypi*fin(0,i,j,k)  +yp*fin(0,i,j+1,k))+ &
                xp*(ypi*fin(0,i+1,j,k)+yp*fin(0,i+1,j+1,k))) &
                +zp*( &
                xpi*(ypi*fin(0,i,j,k+1)  +yp*fin(0,i,j+1,k+1))+ &
                xp*(ypi*fin(0,i+1,j,k+1)+yp*fin(0,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hx2*( &
                zpi*( &
                cxi*(ypi*fin(1,i,j,k)  +yp*fin(1,i,j+1,k))+ &
                cx*(ypi*fin(1,i+1,j,k)+yp*fin(1,i+1,j+1,k))) &
                +zp*( &
                cxi*(ypi*fin(1,i,j,k+1)  +yp*fin(1,i,j+1,k+1))+ &
                cx*(ypi*fin(1,i+1,j,k+1)+yp*fin(1,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hy2*( &
                zpi*( &
                xpi*(cyi*fin(2,i,j,k)  +cy*fin(2,i,j+1,k))+ &
                xp*(cyi*fin(2,i+1,j,k)+cy*fin(2,i+1,j+1,k))) &
                +zp*( &
                xpi*(cyi*fin(2,i,j,k+1)  +cy*fin(2,i,j+1,k+1))+ &
                xp*(cyi*fin(2,i+1,j,k+1)+cy*fin(2,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hz2*( &
                czi*( &
                xpi*(ypi*fin(3,i,j,k)  +yp*fin(3,i,j+1,k))+ &
                xp*(ypi*fin(3,i+1,j,k)+yp*fin(3,i+1,j+1,k))) &
                +cz*( &
                xpi*(ypi*fin(3,i,j,k+1)  +yp*fin(3,i,j+1,k+1))+ &
                xp*(ypi*fin(3,i+1,j,k+1)+yp*fin(3,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hx2*hy2*( &
                zpi*( &
                cxi*(cyi*fin(4,i,j,k)  +cy*fin(4,i,j+1,k))+ &
                cx*(cyi*fin(4,i+1,j,k)+cy*fin(4,i+1,j+1,k))) &
                +zp*( &
                cxi*(cyi*fin(4,i,j,k+1)  +cy*fin(4,i,j+1,k+1))+ &
                cx*(cyi*fin(4,i+1,j,k+1)+cy*fin(4,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hx2*hz2*( &
                czi*( &
                cxi*(ypi*fin(5,i,j,k)  +yp*fin(5,i,j+1,k))+ &
                cx*(ypi*fin(5,i+1,j,k)+yp*fin(5,i+1,j+1,k))) &
                +cz*( &
                cxi*(ypi*fin(5,i,j,k+1)  +yp*fin(5,i,j+1,k+1))+ &
                cx*(ypi*fin(5,i+1,j,k+1)+yp*fin(5,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hy2*hz2*( &
                czi*( &
                xpi*(cyi*fin(6,i,j,k)  +cy*fin(6,i,j+1,k))+ &
                xp*(cyi*fin(6,i+1,j,k)+cy*fin(6,i+1,j+1,k))) &
                +cz*( &
                xpi*(cyi*fin(6,i,j,k+1)  +cy*fin(6,i,j+1,k+1))+ &
                xp*(cyi*fin(6,i+1,j,k+1)+cy*fin(6,i+1,j+1,k+1))))
           !
           sum=sum+z216th*hx2*hy2*hz2*( &
                czi*( &
                cxi*(cyi*fin(7,i,j,k)  +cy*fin(7,i,j+1,k))+ &
                cx*(cyi*fin(7,i+1,j,k)+cy*fin(7,i+1,j+1,k))) &
                +cz*( &
                cxi*(cyi*fin(7,i,j,k+1)  +cy*fin(7,i,j+1,k+1))+ &
                cx*(cyi*fin(7,i+1,j,k+1)+cy*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(2).eq.1) then
        !
        !  df/dx:
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           j=jj(v)
           k=kk(v)
           !
           !   ...in x direction
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           xp2=xp*xp
           xpi2=xpi*xpi

           cxd=3.0_fp*xp2-1.0_fp
           cxdi=-3.0_fp*xpi2+1.0_fp
           !
           !   ...and in y direction
           !
           yp=yparam(v)
           ypi=1.0_fp-yp
           yp2=yp*yp
           ypi2=ypi*ypi
           !
           cy=yp*(yp2-1.0_fp)
           cyi=ypi*(ypi2-1.0_fp)
           hy2=hy(v)*hy(v)
           !
           !   ...and in z direction
           !
           zp=zparam(v)
           zpi=1.0_fp-zp
           zp2=zp*zp
           zpi2=zpi*zpi
           !
           cz=zp*(zp2-1.0_fp)
           czi=zpi*(zpi2-1.0_fp)
           hz2=hz(v)*hz(v)
           !
           sum=hxi(v)*( &
                zpi*( &
                -(ypi*fin(0,i,j,k)  +yp*fin(0,i,j+1,k)) &
                +(ypi*fin(0,i+1,j,k)+yp*fin(0,i+1,j+1,k))) &
                +zp*( &
                -(ypi*fin(0,i,j,k+1)  +yp*fin(0,i,j+1,k+1)) &
                +(ypi*fin(0,i+1,j,k+1)+yp*fin(0,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hx(v)*( &
                zpi*( &
                cxdi*(ypi*fin(1,i,j,k)  +yp*fin(1,i,j+1,k))+ &
                cxd*(ypi*fin(1,i+1,j,k)+yp*fin(1,i+1,j+1,k))) &
                +zp*( &
                cxdi*(ypi*fin(1,i,j,k+1)  +yp*fin(1,i,j+1,k+1))+ &
                cxd*(ypi*fin(1,i+1,j,k+1)+yp*fin(1,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hxi(v)*hy2*( &
                zpi*( &
                -(cyi*fin(2,i,j,k)  +cy*fin(2,i,j+1,k)) &
                +(cyi*fin(2,i+1,j,k)+cy*fin(2,i+1,j+1,k))) &
                +zp*( &
                -(cyi*fin(2,i,j,k+1)  +cy*fin(2,i,j+1,k+1)) &
                +(cyi*fin(2,i+1,j,k+1)+cy*fin(2,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hxi(v)*hz2*( &
                czi*( &
                -(ypi*fin(3,i,j,k)  +yp*fin(3,i,j+1,k)) &
                +(ypi*fin(3,i+1,j,k)+yp*fin(3,i+1,j+1,k))) &
                +cz*( &
                -(ypi*fin(3,i,j,k+1)  +yp*fin(3,i,j+1,k+1)) &
                +(ypi*fin(3,i+1,j,k+1)+yp*fin(3,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hx(v)*hy2*( &
                zpi*( &
                cxdi*(cyi*fin(4,i,j,k)  +cy*fin(4,i,j+1,k))+ &
                cxd*(cyi*fin(4,i+1,j,k)+cy*fin(4,i+1,j+1,k))) &
                +zp*( &
                cxdi*(cyi*fin(4,i,j,k+1)  +cy*fin(4,i,j+1,k+1))+ &
                cxd*(cyi*fin(4,i+1,j,k+1)+cy*fin(4,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hx(v)*hz2*( &
                czi*( &
                cxdi*(ypi*fin(5,i,j,k)  +yp*fin(5,i,j+1,k))+ &
                cxd*(ypi*fin(5,i+1,j,k)+yp*fin(5,i+1,j+1,k))) &
                +cz*( &
                cxdi*(ypi*fin(5,i,j,k+1)  +yp*fin(5,i,j+1,k+1))+ &
                cxd*(ypi*fin(5,i+1,j,k+1)+yp*fin(5,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hxi(v)*hy2*hz2*( &
                czi*( &
                -(cyi*fin(6,i,j,k)  +cy*fin(6,i,j+1,k)) &
                +(cyi*fin(6,i+1,j,k)+cy*fin(6,i+1,j+1,k))) &
                +cz*( &
                -(cyi*fin(6,i,j,k+1)  +cy*fin(6,i,j+1,k+1)) &
                +(cyi*fin(6,i+1,j,k+1)+cy*fin(6,i+1,j+1,k+1))))
           !
           sum=sum+z216th*hx(v)*hy2*hz2*( &
                czi*( &
                cxdi*(cyi*fin(7,i,j,k)  +cy*fin(7,i,j+1,k))+ &
                cxd*(cyi*fin(7,i+1,j,k)+cy*fin(7,i+1,j+1,k))) &
                +cz*( &
                cxdi*(cyi*fin(7,i,j,k+1)  +cy*fin(7,i,j+1,k+1))+ &
                cxd*(cyi*fin(7,i+1,j,k+1)+cy*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(3).eq.1) then
        !
        !  df/dy:
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           j=jj(v)
           k=kk(v)
           !
           !   ...in x direction
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           xp2=xp*xp
           xpi2=xpi*xpi
           !
           cx=xp*(xp2-1.0_fp)
           cxi=xpi*(xpi2-1.0_fp)
           hx2=hx(v)*hx(v)
           !
           !   ...and in y direction
           !
           yp=yparam(v)
           ypi=1.0_fp-yp
           yp2=yp*yp
           ypi2=ypi*ypi

           cyd=3.0_fp*yp2-1.0_fp
           cydi=-3.0_fp*ypi2+1.0_fp
           !
           !   ...and in z direction
           !
           zp=zparam(v)
           zpi=1.0_fp-zp
           zp2=zp*zp
           zpi2=zpi*zpi
           !
           cz=zp*(zp2-1.0_fp)
           czi=zpi*(zpi2-1.0_fp)
           hz2=hz(v)*hz(v)
           !
           sum=hyi(v)*( &
                zpi*( &
                xpi*(-fin(0,i,j,k)  +fin(0,i,j+1,k))+ &
                xp*(-fin(0,i+1,j,k)+fin(0,i+1,j+1,k))) &
                +zp*( &
                xpi*(-fin(0,i,j,k+1)  +fin(0,i,j+1,k+1))+ &
                xp*(-fin(0,i+1,j,k+1)+fin(0,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hyi(v)*hx2*( &
                zpi*( &
                cxi*(-fin(1,i,j,k)  +fin(1,i,j+1,k))+ &
                cx*(-fin(1,i+1,j,k)+fin(1,i+1,j+1,k))) &
                +zp*( &
                cxi*(-fin(1,i,j,k+1)  +fin(1,i,j+1,k+1))+ &
                cx*(-fin(1,i+1,j,k+1)+fin(1,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hy(v)*( &
                zpi*( &
                xpi*(cydi*fin(2,i,j,k)  +cyd*fin(2,i,j+1,k))+ &
                xp*(cydi*fin(2,i+1,j,k)+cyd*fin(2,i+1,j+1,k))) &
                +zp*( &
                xpi*(cydi*fin(2,i,j,k+1)  +cyd*fin(2,i,j+1,k+1))+ &
                xp*(cydi*fin(2,i+1,j,k+1)+cyd*fin(2,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hyi(v)*hz2*( &
                czi*( &
                xpi*(-fin(3,i,j,k)  +fin(3,i,j+1,k))+ &
                xp*(-fin(3,i+1,j,k)+fin(3,i+1,j+1,k))) &
                +cz*( &
                xpi*(-fin(3,i,j,k+1)  +fin(3,i,j+1,k+1))+ &
                xp*(-fin(3,i+1,j,k+1)+fin(3,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hx2*hy(v)*( &
                zpi*( &
                cxi*(cydi*fin(4,i,j,k)  +cyd*fin(4,i,j+1,k))+ &
                cx*(cydi*fin(4,i+1,j,k)+cyd*fin(4,i+1,j+1,k))) &
                +zp*( &
                cxi*(cydi*fin(4,i,j,k+1)  +cyd*fin(4,i,j+1,k+1))+ &
                cx*(cydi*fin(4,i+1,j,k+1)+cyd*fin(4,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hyi(v)*hx2*hz2*( &
                czi*( &
                cxi*(-fin(5,i,j,k)  +fin(5,i,j+1,k))+ &
                cx*(-fin(5,i+1,j,k)+fin(5,i+1,j+1,k))) &
                +cz*( &
                cxi*(-fin(5,i,j,k+1)  +fin(5,i,j+1,k+1))+ &
                cx*(-fin(5,i+1,j,k+1)+fin(5,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hy(v)*hz2*( &
                czi*( &
                xpi*(cydi*fin(6,i,j,k)  +cyd*fin(6,i,j+1,k))+ &
                xp*(cydi*fin(6,i+1,j,k)+cyd*fin(6,i+1,j+1,k))) &
                +cz*( &
                xpi*(cydi*fin(6,i,j,k+1)  +cyd*fin(6,i,j+1,k+1))+ &
                xp*(cydi*fin(6,i+1,j,k+1)+cyd*fin(6,i+1,j+1,k+1))))
           !
           sum=sum+z216th*hx2*hy(v)*hz2*( &
                czi*( &
                cxi*(cydi*fin(7,i,j,k)  +cyd*fin(7,i,j+1,k))+ &
                cx*(cydi*fin(7,i+1,j,k)+cyd*fin(7,i+1,j+1,k))) &
                +cz*( &
                cxi*(cydi*fin(7,i,j,k+1)  +cyd*fin(7,i,j+1,k+1))+ &
                cx*(cydi*fin(7,i+1,j,k+1)+cyd*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(4).eq.1) then
        !
        !  df/dz:
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           j=jj(v)
           k=kk(v)
           !
           !   ...in x direction
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           xp2=xp*xp
           xpi2=xpi*xpi
           !
           cx=xp*(xp2-1.0_fp)
           cxi=xpi*(xpi2-1.0_fp)
           hx2=hx(v)*hx(v)
           !
           !   ...and in y direction
           !
           yp=yparam(v)
           ypi=1.0_fp-yp
           yp2=yp*yp
           ypi2=ypi*ypi
           !
           cy=yp*(yp2-1.0_fp)
           cyi=ypi*(ypi2-1.0_fp)
           hy2=hy(v)*hy(v)
           !
           !   ...and in z direction
           !
           zp=zparam(v)
           zpi=1.0_fp-zp
           zp2=zp*zp
           zpi2=zpi*zpi

           czd=3.0_fp*zp2-1.0_fp
           czdi=-3.0_fp*zpi2+1.0_fp
           !
           sum=hzi(v)*( &
                -( &
                xpi*(ypi*fin(0,i,j,k)  +yp*fin(0,i,j+1,k))+ &
                xp*(ypi*fin(0,i+1,j,k)+yp*fin(0,i+1,j+1,k))) &
                +( &
                xpi*(ypi*fin(0,i,j,k+1)  +yp*fin(0,i,j+1,k+1))+ &
                xp*(ypi*fin(0,i+1,j,k+1)+yp*fin(0,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hx2*hzi(v)*( &
                -( &
                cxi*(ypi*fin(1,i,j,k)  +yp*fin(1,i,j+1,k))+ &
                cx*(ypi*fin(1,i+1,j,k)+yp*fin(1,i+1,j+1,k))) &
                +( &
                cxi*(ypi*fin(1,i,j,k+1)  +yp*fin(1,i,j+1,k+1))+ &
                cx*(ypi*fin(1,i+1,j,k+1)+yp*fin(1,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hy2*hzi(v)*( &
                -( &
                xpi*(cyi*fin(2,i,j,k)  +cy*fin(2,i,j+1,k))+ &
                xp*(cyi*fin(2,i+1,j,k)+cy*fin(2,i+1,j+1,k))) &
                +( &
                xpi*(cyi*fin(2,i,j,k+1)  +cy*fin(2,i,j+1,k+1))+ &
                xp*(cyi*fin(2,i+1,j,k+1)+cy*fin(2,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hz(v)*( &
                czdi*( &
                xpi*(ypi*fin(3,i,j,k)  +yp*fin(3,i,j+1,k))+ &
                xp*(ypi*fin(3,i+1,j,k)+yp*fin(3,i+1,j+1,k))) &
                +czd*( &
                xpi*(ypi*fin(3,i,j,k+1)  +yp*fin(3,i,j+1,k+1))+ &
                xp*(ypi*fin(3,i+1,j,k+1)+yp*fin(3,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hx2*hy2*hzi(v)*( &
                -( &
                cxi*(cyi*fin(4,i,j,k)  +cy*fin(4,i,j+1,k))+ &
                cx*(cyi*fin(4,i+1,j,k)+cy*fin(4,i+1,j+1,k))) &
                +( &
                cxi*(cyi*fin(4,i,j,k+1)  +cy*fin(4,i,j+1,k+1))+ &
                cx*(cyi*fin(4,i+1,j,k+1)+cy*fin(4,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hx2*hz(v)*( &
                czdi*( &
                cxi*(ypi*fin(5,i,j,k)  +yp*fin(5,i,j+1,k))+ &
                cx*(ypi*fin(5,i+1,j,k)+yp*fin(5,i+1,j+1,k))) &
                +czd*( &
                cxi*(ypi*fin(5,i,j,k+1)  +yp*fin(5,i,j+1,k+1))+ &
                cx*(ypi*fin(5,i+1,j,k+1)+yp*fin(5,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hy2*hz(v)*( &
                czdi*( &
                xpi*(cyi*fin(6,i,j,k)  +cy*fin(6,i,j+1,k))+ &
                xp*(cyi*fin(6,i+1,j,k)+cy*fin(6,i+1,j+1,k))) &
                +czd*( &
                xpi*(cyi*fin(6,i,j,k+1)  +cy*fin(6,i,j+1,k+1))+ &
                xp*(cyi*fin(6,i+1,j,k+1)+cy*fin(6,i+1,j+1,k+1))))
           !
           sum=sum+z216th*hx2*hy2*hz(v)*( &
                czdi*( &
                cxi*(cyi*fin(7,i,j,k)  +cy*fin(7,i,j+1,k))+ &
                cx*(cyi*fin(7,i+1,j,k)+cy*fin(7,i+1,j+1,k))) &
                +czd*( &
                cxi*(cyi*fin(7,i,j,k+1)  +cy*fin(7,i,j+1,k+1))+ &
                cx*(cyi*fin(7,i+1,j,k+1)+cy*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(5).eq.1) then
        !
        !  d2f/dx2:
        !
        iadr=iadr+1
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
           !   ...and in y direction
           !
           yp=yparam(v)
           ypi=1.0_fp-yp
           yp2=yp*yp
           ypi2=ypi*ypi
           !
           cy=yp*(yp2-1.0_fp)
           cyi=ypi*(ypi2-1.0_fp)
           hy2=hy(v)*hy(v)
           !
           !   ...and in z direction
           !
           zp=zparam(v)
           zpi=1.0_fp-zp
           zp2=zp*zp
           zpi2=zpi*zpi
           !
           cz=zp*(zp2-1.0_fp)
           czi=zpi*(zpi2-1.0_fp)
           hz2=hz(v)*hz(v)
           !
           sum=( &
                zpi*( &
                xpi*(ypi*fin(1,i,j,k)  +yp*fin(1,i,j+1,k))+ &
                xp*(ypi*fin(1,i+1,j,k)+yp*fin(1,i+1,j+1,k))) &
                +zp*( &
                xpi*(ypi*fin(1,i,j,k+1)  +yp*fin(1,i,j+1,k+1))+ &
                xp*(ypi*fin(1,i+1,j,k+1)+yp*fin(1,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hy2*( &
                zpi*( &
                xpi*(cyi*fin(4,i,j,k)  +cy*fin(4,i,j+1,k))+ &
                xp*(cyi*fin(4,i+1,j,k)+cy*fin(4,i+1,j+1,k))) &
                +zp*( &
                xpi*(cyi*fin(4,i,j,k+1)  +cy*fin(4,i,j+1,k+1))+ &
                xp*(cyi*fin(4,i+1,j,k+1)+cy*fin(4,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hz2*( &
                czi*( &
                xpi*(ypi*fin(5,i,j,k)  +yp*fin(5,i,j+1,k))+ &
                xp*(ypi*fin(5,i+1,j,k)+yp*fin(5,i+1,j+1,k))) &
                +cz*( &
                xpi*(ypi*fin(5,i,j,k+1)  +yp*fin(5,i,j+1,k+1))+ &
                xp*(ypi*fin(5,i+1,j,k+1)+yp*fin(5,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hy2*hz2*( &
                czi*( &
                xpi*(cyi*fin(7,i,j,k)  +cy*fin(7,i,j+1,k))+ &
                xp*(cyi*fin(7,i+1,j,k)+cy*fin(7,i+1,j+1,k))) &
                +cz*( &
                xpi*(cyi*fin(7,i,j,k+1)  +cy*fin(7,i,j+1,k+1))+ &
                xp*(cyi*fin(7,i+1,j,k+1)+cy*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(6).eq.1) then
        !
        !  d2f/dy2:
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           j=jj(v)
           k=kk(v)
           !
           !   ...in x direction
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           xp2=xp*xp
           xpi2=xpi*xpi
           !
           cx=xp*(xp2-1.0_fp)
           cxi=xpi*(xpi2-1.0_fp)
           hx2=hx(v)*hx(v)
           !
           !   ...and in y direction
           !
           yp=yparam(v)
           ypi=1.0_fp-yp
           !
           !   ...and in z direction
           !
           zp=zparam(v)
           zpi=1.0_fp-zp
           zp2=zp*zp
           zpi2=zpi*zpi
           !
           cz=zp*(zp2-1.0_fp)
           czi=zpi*(zpi2-1.0_fp)
           hz2=hz(v)*hz(v)
           !
           sum=( &
                zpi*( &
                xpi*(ypi*fin(2,i,j,k)  +yp*fin(2,i,j+1,k))+ &
                xp*(ypi*fin(2,i+1,j,k)+yp*fin(2,i+1,j+1,k))) &
                +zp*( &
                xpi*(ypi*fin(2,i,j,k+1)  +yp*fin(2,i,j+1,k+1))+ &
                xp*(ypi*fin(2,i+1,j,k+1)+yp*fin(2,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hx2*( &
                zpi*( &
                cxi*(ypi*fin(4,i,j,k)  +yp*fin(4,i,j+1,k))+ &
                cx*(ypi*fin(4,i+1,j,k)+yp*fin(4,i+1,j+1,k))) &
                +zp*( &
                cxi*(ypi*fin(4,i,j,k+1)  +yp*fin(4,i,j+1,k+1))+ &
                cx*(ypi*fin(4,i+1,j,k+1)+yp*fin(4,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hz2*( &
                czi*( &
                xpi*(ypi*fin(6,i,j,k)  +yp*fin(6,i,j+1,k))+ &
                xp*(ypi*fin(6,i+1,j,k)+yp*fin(6,i+1,j+1,k))) &
                +cz*( &
                xpi*(ypi*fin(6,i,j,k+1)  +yp*fin(6,i,j+1,k+1))+ &
                xp*(ypi*fin(6,i+1,j,k+1)+yp*fin(6,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hx2*hz2*( &
                czi*( &
                cxi*(ypi*fin(7,i,j,k)  +yp*fin(7,i,j+1,k))+ &
                cx*(ypi*fin(7,i+1,j,k)+yp*fin(7,i+1,j+1,k))) &
                +cz*( &
                cxi*(ypi*fin(7,i,j,k+1)  +yp*fin(7,i,j+1,k+1))+ &
                cx*(ypi*fin(7,i+1,j,k+1)+yp*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(7).eq.1) then
        !
        !  d2f/dz2:
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           j=jj(v)
           k=kk(v)
           !
           !   ...in x direction
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           xp2=xp*xp
           xpi2=xpi*xpi
           !
           cx=xp*(xp2-1.0_fp)
           cxi=xpi*(xpi2-1.0_fp)
           hx2=hx(v)*hx(v)
           !
           !   ...and in y direction
           !
           yp=yparam(v)
           ypi=1.0_fp-yp
           yp2=yp*yp
           ypi2=ypi*ypi
           !
           cy=yp*(yp2-1.0_fp)
           cyi=ypi*(ypi2-1.0_fp)
           hy2=hy(v)*hy(v)
           !
           !   ...and in z direction
           !
           zp=zparam(v)
           zpi=1.0_fp-zp
           !
           sum=( &
                zpi*( &
                xpi*(ypi*fin(3,i,j,k)  +yp*fin(3,i,j+1,k))+ &
                xp*(ypi*fin(3,i+1,j,k)+yp*fin(3,i+1,j+1,k))) &
                +zp*( &
                xpi*(ypi*fin(3,i,j,k+1)  +yp*fin(3,i,j+1,k+1))+ &
                xp*(ypi*fin(3,i+1,j,k+1)+yp*fin(3,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hx2*( &
                zpi*( &
                cxi*(ypi*fin(5,i,j,k)  +yp*fin(5,i,j+1,k))+ &
                cx*(ypi*fin(5,i+1,j,k)+yp*fin(5,i+1,j+1,k))) &
                +zp*( &
                cxi*(ypi*fin(5,i,j,k+1)  +yp*fin(5,i,j+1,k+1))+ &
                cx*(ypi*fin(5,i+1,j,k+1)+yp*fin(5,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hy2*( &
                zpi*( &
                xpi*(cyi*fin(6,i,j,k)  +cy*fin(6,i,j+1,k))+ &
                xp*(cyi*fin(6,i+1,j,k)+cy*fin(6,i+1,j+1,k))) &
                +zp*( &
                xpi*(cyi*fin(6,i,j,k+1)  +cy*fin(6,i,j+1,k+1))+ &
                xp*(cyi*fin(6,i+1,j,k+1)+cy*fin(6,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hx2*hy2*( &
                zpi*( &
                cxi*(cyi*fin(7,i,j,k)  +cy*fin(7,i,j+1,k))+ &
                cx*(cyi*fin(7,i+1,j,k)+cy*fin(7,i+1,j+1,k))) &
                +zp*( &
                cxi*(cyi*fin(7,i,j,k+1)  +cy*fin(7,i,j+1,k+1))+ &
                cx*(cyi*fin(7,i+1,j,k+1)+cy*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(8).eq.1) then
        !
        !  d2f/dxdy:
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           j=jj(v)
           k=kk(v)
           !
           !   ...in x direction
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           xp2=xp*xp
           xpi2=xpi*xpi

           cxd=3.0_fp*xp2-1.0_fp
           cxdi=-3.0_fp*xpi2+1.0_fp
           !
           !   ...and in y direction
           !
           yp=yparam(v)
           ypi=1.0_fp-yp
           yp2=yp*yp
           ypi2=ypi*ypi

           cyd=3.0_fp*yp2-1.0_fp
           cydi=-3.0_fp*ypi2+1.0_fp
           !
           !   ...and in z direction
           !
           zp=zparam(v)
           zpi=1.0_fp-zp
           zp2=zp*zp
           zpi2=zpi*zpi
           !
           cz=zp*(zp2-1.0_fp)
           czi=zpi*(zpi2-1.0_fp)
           hz2=hz(v)*hz(v)
           !
           sum=hxi(v)*hyi(v)*( &
                zpi*( &
                (fin(0,i,j,k)  -fin(0,i,j+1,k))- &
                (fin(0,i+1,j,k)-fin(0,i+1,j+1,k))) &
                +zp*( &
                (fin(0,i,j,k+1)  -fin(0,i,j+1,k+1))- &
                (fin(0,i+1,j,k+1)-fin(0,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hyi(v)*hx(v)*( &
                zpi*( &
                cxdi*(-fin(1,i,j,k)  +fin(1,i,j+1,k))+ &
                cxd*(-fin(1,i+1,j,k)+fin(1,i+1,j+1,k))) &
                +zp*( &
                cxdi*(-fin(1,i,j,k+1)  +fin(1,i,j+1,k+1))+ &
                cxd*(-fin(1,i+1,j,k+1)+fin(1,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hxi(v)*hy(v)*( &
                zpi*( &
                -(cydi*fin(2,i,j,k)  +cyd*fin(2,i,j+1,k)) &
                +(cydi*fin(2,i+1,j,k)+cyd*fin(2,i+1,j+1,k))) &
                +zp*( &
                -(cydi*fin(2,i,j,k+1)  +cyd*fin(2,i,j+1,k+1)) &
                +(cydi*fin(2,i+1,j,k+1)+cyd*fin(2,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hxi(v)*hyi(v)*hz2*( &
                czi*( &
                (fin(3,i,j,k)  -fin(3,i,j+1,k))- &
                (fin(3,i+1,j,k)-fin(3,i+1,j+1,k))) &
                +cz*( &
                (fin(3,i,j,k+1)  -fin(3,i,j+1,k+1))- &
                (fin(3,i+1,j,k+1)-fin(3,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hx(v)*hy(v)*( &
                zpi*( &
                cxdi*(cydi*fin(4,i,j,k)  +cyd*fin(4,i,j+1,k))+ &
                cxd*(cydi*fin(4,i+1,j,k)+cyd*fin(4,i+1,j+1,k))) &
                +zp*( &
                cxdi*(cydi*fin(4,i,j,k+1)  +cyd*fin(4,i,j+1,k+1))+ &
                cxd*(cydi*fin(4,i+1,j,k+1)+cyd*fin(4,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hyi(v)*hx(v)*hz2*( &
                czi*( &
                cxdi*(-fin(5,i,j,k)  +fin(5,i,j+1,k))+ &
                cxd*(-fin(5,i+1,j,k)+fin(5,i+1,j+1,k))) &
                +cz*( &
                cxdi*(-fin(5,i,j,k+1)  +fin(5,i,j+1,k+1))+ &
                cxd*(-fin(5,i+1,j,k+1)+fin(5,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hxi(v)*hy(v)*hz2*( &
                czi*( &
                -(cydi*fin(6,i,j,k)  +cyd*fin(6,i,j+1,k)) &
                +(cydi*fin(6,i+1,j,k)+cyd*fin(6,i+1,j+1,k))) &
                +cz*( &
                -(cydi*fin(6,i,j,k+1)  +cyd*fin(6,i,j+1,k+1)) &
                +(cydi*fin(6,i+1,j,k+1)+cyd*fin(6,i+1,j+1,k+1))))
           !
           sum=sum+z216th*hx(v)*hy(v)*hz2*( &
                czi*( &
                cxdi*(cydi*fin(7,i,j,k)  +cyd*fin(7,i,j+1,k))+ &
                cxd*(cydi*fin(7,i+1,j,k)+cyd*fin(7,i+1,j+1,k))) &
                +cz*( &
                cxdi*(cydi*fin(7,i,j,k+1)  +cyd*fin(7,i,j+1,k+1))+ &
                cxd*(cydi*fin(7,i+1,j,k+1)+cyd*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(9).eq.1) then
        !
        !  d2f/dxdz:
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           j=jj(v)
           k=kk(v)
           !
           !   ...in x direction
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           xp2=xp*xp
           xpi2=xpi*xpi

           cxd=3.0_fp*xp2-1.0_fp
           cxdi=-3.0_fp*xpi2+1.0_fp
           !
           !   ...and in y direction
           !
           yp=yparam(v)
           ypi=1.0_fp-yp
           yp2=yp*yp
           ypi2=ypi*ypi
           !
           cy=yp*(yp2-1.0_fp)
           cyi=ypi*(ypi2-1.0_fp)
           hy2=hy(v)*hy(v)
           !
           !   ...and in z direction
           !
           zp=zparam(v)
           zpi=1.0_fp-zp
           zp2=zp*zp
           zpi2=zpi*zpi

           czd=3.0_fp*zp2-1.0_fp
           czdi=-3.0_fp*zpi2+1.0_fp
           !
           sum=hxi(v)*hzi(v)*( &
                ( &
                (ypi*fin(0,i,j,k)  +yp*fin(0,i,j+1,k)) - &
                (ypi*fin(0,i+1,j,k)+yp*fin(0,i+1,j+1,k))) &
                -( &
                (ypi*fin(0,i,j,k+1)  +yp*fin(0,i,j+1,k+1)) - &
                (ypi*fin(0,i+1,j,k+1)+yp*fin(0,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hx(v)*hzi(v)*( &
                -( &
                cxdi*(ypi*fin(1,i,j,k)  +yp*fin(1,i,j+1,k))+ &
                cxd*(ypi*fin(1,i+1,j,k)+yp*fin(1,i+1,j+1,k))) &
                +( &
                cxdi*(ypi*fin(1,i,j,k+1)  +yp*fin(1,i,j+1,k+1))+ &
                cxd*(ypi*fin(1,i+1,j,k+1)+yp*fin(1,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hxi(v)*hy2*hzi(v)*( &
                ( &
                (cyi*fin(2,i,j,k)  +cy*fin(2,i,j+1,k)) - &
                (cyi*fin(2,i+1,j,k)+cy*fin(2,i+1,j+1,k))) &
                -( &
                (cyi*fin(2,i,j,k+1)  +cy*fin(2,i,j+1,k+1)) - &
                (cyi*fin(2,i+1,j,k+1)+cy*fin(2,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hxi(v)*hz(v)*( &
                czdi*( &
                -(ypi*fin(3,i,j,k)  +yp*fin(3,i,j+1,k)) &
                +(ypi*fin(3,i+1,j,k)+yp*fin(3,i+1,j+1,k))) &
                +czd*( &
                -(ypi*fin(3,i,j,k+1)  +yp*fin(3,i,j+1,k+1)) &
                +(ypi*fin(3,i+1,j,k+1)+yp*fin(3,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hx(v)*hy2*hzi(v)*( &
                -( &
                cxdi*(cyi*fin(4,i,j,k)  +cy*fin(4,i,j+1,k))+ &
                cxd*(cyi*fin(4,i+1,j,k)+cy*fin(4,i+1,j+1,k))) &
                +( &
                cxdi*(cyi*fin(4,i,j,k+1)  +cy*fin(4,i,j+1,k+1))+ &
                cxd*(cyi*fin(4,i+1,j,k+1)+cy*fin(4,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hx(v)*hz(v)*( &
                czdi*( &
                cxdi*(ypi*fin(5,i,j,k)  +yp*fin(5,i,j+1,k))+ &
                cxd*(ypi*fin(5,i+1,j,k)+yp*fin(5,i+1,j+1,k))) &
                +czd*( &
                cxdi*(ypi*fin(5,i,j,k+1)  +yp*fin(5,i,j+1,k+1))+ &
                cxd*(ypi*fin(5,i+1,j,k+1)+yp*fin(5,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hxi(v)*hy2*hz(v)*( &
                czdi*( &
                -(cyi*fin(6,i,j,k)  +cy*fin(6,i,j+1,k)) &
                +(cyi*fin(6,i+1,j,k)+cy*fin(6,i+1,j+1,k))) &
                +czd*( &
                -(cyi*fin(6,i,j,k+1)  +cy*fin(6,i,j+1,k+1)) &
                +(cyi*fin(6,i+1,j,k+1)+cy*fin(6,i+1,j+1,k+1))))
           !
           sum=sum+z216th*hx(v)*hy2*hz(v)*( &
                czdi*( &
                cxdi*(cyi*fin(7,i,j,k)  +cy*fin(7,i,j+1,k))+ &
                cxd*(cyi*fin(7,i+1,j,k)+cy*fin(7,i+1,j+1,k))) &
                +czd*( &
                cxdi*(cyi*fin(7,i,j,k+1)  +cy*fin(7,i,j+1,k+1))+ &
                cxd*(cyi*fin(7,i+1,j,k+1)+cy*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(10).eq.1) then
        !
        !  d2f/dydz:
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           j=jj(v)
           k=kk(v)
           !
           !   ...in x direction
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           xp2=xp*xp
           xpi2=xpi*xpi
           !
           cx=xp*(xp2-1.0_fp)
           cxi=xpi*(xpi2-1.0_fp)
           hx2=hx(v)*hx(v)
           !
           !   ...and in y direction
           !
           yp=yparam(v)
           ypi=1.0_fp-yp
           yp2=yp*yp
           ypi2=ypi*ypi

           cyd=3.0_fp*yp2-1.0_fp
           cydi=-3.0_fp*ypi2+1.0_fp
           !
           !   ...and in z direction
           !
           zp=zparam(v)
           zpi=1.0_fp-zp
           zp2=zp*zp
           zpi2=zpi*zpi

           czd=3.0_fp*zp2-1.0_fp
           czdi=-3.0_fp*zpi2+1.0_fp
           !
           sum=hyi(v)*hzi(v)*( &
                ( &
                xpi*(fin(0,i,j,k)  -fin(0,i,j+1,k))+ &
                xp*(fin(0,i+1,j,k)-fin(0,i+1,j+1,k))) &
                -( &
                xpi*(fin(0,i,j,k+1)  -fin(0,i,j+1,k+1))+ &
                xp*(fin(0,i+1,j,k+1)-fin(0,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hyi(v)*hx2*hzi(v)*( &
                ( &
                cxi*(fin(1,i,j,k)  -fin(1,i,j+1,k))+ &
                cx*(fin(1,i+1,j,k)-fin(1,i+1,j+1,k))) &
                -( &
                cxi*(fin(1,i,j,k+1)  -fin(1,i,j+1,k+1))+ &
                cx*(fin(1,i+1,j,k+1)-fin(1,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hy(v)*hzi(v)*( &
                -( &
                xpi*(cydi*fin(2,i,j,k)  +cyd*fin(2,i,j+1,k))+ &
                xp*(cydi*fin(2,i+1,j,k)+cyd*fin(2,i+1,j+1,k))) &
                +( &
                xpi*(cydi*fin(2,i,j,k+1)  +cyd*fin(2,i,j+1,k+1))+ &
                xp*(cydi*fin(2,i+1,j,k+1)+cyd*fin(2,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hyi(v)*hz(v)*( &
                czdi*( &
                xpi*(-fin(3,i,j,k)  +fin(3,i,j+1,k))+ &
                xp*(-fin(3,i+1,j,k)+fin(3,i+1,j+1,k))) &
                +czd*( &
                xpi*(-fin(3,i,j,k+1)  +fin(3,i,j+1,k+1))+ &
                xp*(-fin(3,i+1,j,k+1)+fin(3,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hx2*hy(v)*hzi(v)*( &
                -( &
                cxi*(cydi*fin(4,i,j,k)  +cyd*fin(4,i,j+1,k))+ &
                cx*(cydi*fin(4,i+1,j,k)+cyd*fin(4,i+1,j+1,k))) &
                +( &
                cxi*(cydi*fin(4,i,j,k+1)  +cyd*fin(4,i,j+1,k+1))+ &
                cx*(cydi*fin(4,i+1,j,k+1)+cyd*fin(4,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hyi(v)*hx2*hz(v)*( &
                czdi*( &
                cxi*(-fin(5,i,j,k)  +fin(5,i,j+1,k))+ &
                cx*(-fin(5,i+1,j,k)+fin(5,i+1,j+1,k))) &
                +czd*( &
                cxi*(-fin(5,i,j,k+1)  +fin(5,i,j+1,k+1))+ &
                cx*(-fin(5,i+1,j,k+1)+fin(5,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hy(v)*hz(v)*( &
                czdi*( &
                xpi*(cydi*fin(6,i,j,k)  +cyd*fin(6,i,j+1,k))+ &
                xp*(cydi*fin(6,i+1,j,k)+cyd*fin(6,i+1,j+1,k))) &
                +czd*( &
                xpi*(cydi*fin(6,i,j,k+1)  +cyd*fin(6,i,j+1,k+1))+ &
                xp*(cydi*fin(6,i+1,j,k+1)+cyd*fin(6,i+1,j+1,k+1))))
           !
           sum=sum+z216th*hx2*hy(v)*hz(v)*( &
                czdi*( &
                cxi*(cydi*fin(7,i,j,k)  +cyd*fin(7,i,j+1,k))+ &
                cx*(cydi*fin(7,i+1,j,k)+cyd*fin(7,i+1,j+1,k))) &
                +czd*( &
                cxi*(cydi*fin(7,i,j,k+1)  +cyd*fin(7,i,j+1,k+1))+ &
                cx*(cydi*fin(7,i+1,j,k+1)+cyd*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     !----------------------------------
     !  3rd derivatives (.le.2 in each coordinate)
     !
  else if(ict(1).eq.3) then
     if(ict(2).eq.1) then
        !                               ! fxxy
        iadr=iadr+1
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
           !   ...and in y direction
           !
           yp=yparam(v)
           ypi=1.0_fp-yp
           yp2=yp*yp
           ypi2=ypi*ypi
           !
           cyd=3.0_fp*yp2-1.0_fp
           cydi=-3.0_fp*ypi2+1.0_fp
           !
           !   ...and in z direction
           !
           zp=zparam(v)
           zpi=1.0_fp-zp
           zp2=zp*zp
           zpi2=zpi*zpi
           !
           cz=zp*(zp2-1.0_fp)
           czi=zpi*(zpi2-1.0_fp)
           hz2=hz(v)*hz(v)
           !
           sum=hyi(v)*( &
                zpi*( &
                xpi*(-fin(1,i,j,k)  +fin(1,i,j+1,k))+ &
                xp*( -fin(1,i+1,j,k)+fin(1,i+1,j+1,k))) &
                +zp*( &
                xpi*(-fin(1,i,j,k+1)  +fin(1,i,j+1,k+1))+ &
                xp*( -fin(1,i+1,j,k+1)+fin(1,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hy(v)*( &
                zpi*( &
                xpi*(cydi*fin(4,i,j,k)  +cyd*fin(4,i,j+1,k))+ &
                xp*(cydi*fin(4,i+1,j,k)+cyd*fin(4,i+1,j+1,k))) &
                +zp*( &
                xpi*(cydi*fin(4,i,j,k+1)  +cyd*fin(4,i,j+1,k+1))+ &
                xp*(cydi*fin(4,i+1,j,k+1)+cyd*fin(4,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hz2*hyi(v)*( &
                czi*( &
                xpi*(-fin(5,i,j,k)  +fin(5,i,j+1,k))+ &
                xp*( -fin(5,i+1,j,k)+fin(5,i+1,j+1,k))) &
                +cz*( &
                xpi*(-fin(5,i,j,k+1)  +fin(5,i,j+1,k+1))+ &
                xp*( -fin(5,i+1,j,k+1)+fin(5,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hy(v)*hz2*( &
                czi*( &
                xpi*(cydi*fin(7,i,j,k)  +cyd*fin(7,i,j+1,k))+ &
                xp*(cydi*fin(7,i+1,j,k)+cyd*fin(7,i+1,j+1,k))) &
                +cz*( &
                xpi*(cydi*fin(7,i,j,k+1)  +cyd*fin(7,i,j+1,k+1))+ &
                xp*(cydi*fin(7,i+1,j,k+1)+cyd*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(3).eq.1) then
        !                               ! fxxz
        iadr=iadr+1
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
           !   ...and in y direction
           !
           yp=yparam(v)
           ypi=1.0_fp-yp
           yp2=yp*yp
           ypi2=ypi*ypi
           !
           cy=yp*(yp2-1.0_fp)
           cyi=ypi*(ypi2-1.0_fp)
           hy2=hy(v)*hy(v)
           !
           !   ...and in z direction
           !
           zp=zparam(v)
           zpi=1.0_fp-zp
           zp2=zp*zp
           zpi2=zpi*zpi
           !
           czd=3.0_fp*zp2-1.0_fp
           czdi=-3.0_fp*zpi2+1.0_fp
           !
           sum=hzi(v)*( &
                -( &
                xpi*(ypi*fin(1,i,j,k)  +yp*fin(1,i,j+1,k))+ &
                xp*(ypi*fin(1,i+1,j,k)+yp*fin(1,i+1,j+1,k))) &
                +( &
                xpi*(ypi*fin(1,i,j,k+1)  +yp*fin(1,i,j+1,k+1))+ &
                xp*(ypi*fin(1,i+1,j,k+1)+yp*fin(1,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hy2*hzi(v)*( &
                -( &
                xpi*(cyi*fin(4,i,j,k)  +cy*fin(4,i,j+1,k))+ &
                xp*(cyi*fin(4,i+1,j,k)+cy*fin(4,i+1,j+1,k))) &
                +( &
                xpi*(cyi*fin(4,i,j,k+1)  +cy*fin(4,i,j+1,k+1))+ &
                xp*(cyi*fin(4,i+1,j,k+1)+cy*fin(4,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hz(v)*( &
                czdi*( &
                xpi*(ypi*fin(5,i,j,k)  +yp*fin(5,i,j+1,k))+ &
                xp*(ypi*fin(5,i+1,j,k)+yp*fin(5,i+1,j+1,k))) &
                +czd*( &
                xpi*(ypi*fin(5,i,j,k+1)  +yp*fin(5,i,j+1,k+1))+ &
                xp*(ypi*fin(5,i+1,j,k+1)+yp*fin(5,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hy2*hz(v)*( &
                czdi*( &
                xpi*(cyi*fin(7,i,j,k)  +cy*fin(7,i,j+1,k))+ &
                xp*(cyi*fin(7,i+1,j,k)+cy*fin(7,i+1,j+1,k))) &
                +czd*( &
                xpi*(cyi*fin(7,i,j,k+1)  +cy*fin(7,i,j+1,k+1))+ &
                xp*(cyi*fin(7,i+1,j,k+1)+cy*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(4).eq.1) then
        !                               ! fxyy
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           j=jj(v)
           k=kk(v)
           !
           !   ...in x direction
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           xp2=xp*xp
           xpi2=xpi*xpi
           !
           cxd=3.0_fp*xp2-1.0_fp
           cxdi=-3.0_fp*xpi2+1.0_fp
           !
           !   ...and in y direction
           !
           yp=yparam(v)
           ypi=1.0_fp-yp
           !
           !   ...and in z direction
           !
           zp=zparam(v)
           zpi=1.0_fp-zp
           zp2=zp*zp
           zpi2=zpi*zpi
           !
           cz=zp*(zp2-1.0_fp)
           czi=zpi*(zpi2-1.0_fp)
           hz2=hz(v)*hz(v)
           !
           sum=hxi(v)*( &
                zpi*( &
                -(ypi*fin(2,i,j,k)  +yp*fin(2,i,j+1,k)) &
                +(ypi*fin(2,i+1,j,k)+yp*fin(2,i+1,j+1,k))) &
                +zp*( &
                -(ypi*fin(2,i,j,k+1)  +yp*fin(2,i,j+1,k+1)) &
                +(ypi*fin(2,i+1,j,k+1)+yp*fin(2,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hx(v)*( &
                zpi*( &
                cxdi*(ypi*fin(4,i,j,k)  +yp*fin(4,i,j+1,k))+ &
                cxd*(ypi*fin(4,i+1,j,k)+yp*fin(4,i+1,j+1,k))) &
                +zp*( &
                cxdi*(ypi*fin(4,i,j,k+1)  +yp*fin(4,i,j+1,k+1))+ &
                cxd*(ypi*fin(4,i+1,j,k+1)+yp*fin(4,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hz2*hxi(v)*( &
                czi*( &
                -(ypi*fin(6,i,j,k)  +yp*fin(6,i,j+1,k)) &
                +(ypi*fin(6,i+1,j,k)+yp*fin(6,i+1,j+1,k))) &
                +cz*( &
                -(ypi*fin(6,i,j,k+1)  +yp*fin(6,i,j+1,k+1)) &
                +(ypi*fin(6,i+1,j,k+1)+yp*fin(6,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hx(v)*hz2*( &
                czi*( &
                cxdi*(ypi*fin(7,i,j,k)  +yp*fin(7,i,j+1,k))+ &
                cxd*(ypi*fin(7,i+1,j,k)+yp*fin(7,i+1,j+1,k))) &
                +cz*( &
                cxdi*(ypi*fin(7,i,j,k+1)  +yp*fin(7,i,j+1,k+1))+ &
                cxd*(ypi*fin(7,i+1,j,k+1)+yp*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(5).eq.1) then
        !                               ! fxyz
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           j=jj(v)
           k=kk(v)
           !
           !   ...in x direction
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           xp2=xp*xp
           xpi2=xpi*xpi
           !
           cxd=3.0_fp*xp2-1.0_fp
           cxdi=-3.0_fp*xpi2+1.0_fp
           !
           !   ...and in y direction
           !
           yp=yparam(v)
           ypi=1.0_fp-yp
           yp2=yp*yp
           ypi2=ypi*ypi
           !
           cyd=3.0_fp*yp2-1.0_fp
           cydi=-3.0_fp*ypi2+1.0_fp
           !
           !   ...and in z direction
           !
           zp=zparam(v)
           zpi=1.0_fp-zp
           zp2=zp*zp
           zpi2=zpi*zpi
           !
           czd=3.0_fp*zp2-1.0_fp
           czdi=-3.0_fp*zpi2+1.0_fp
           !
           sum=hxi(v)*hyi(v)*hzi(v)*( &
                -( &
                (fin(0,i,j,k)  -fin(0,i,j+1,k))- &
                (fin(0,i+1,j,k)-fin(0,i+1,j+1,k))) &
                +( &
                (fin(0,i,j,k+1)  -fin(0,i,j+1,k+1))- &
                (fin(0,i+1,j,k+1)-fin(0,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hyi(v)*hx(v)*hzi(v)*( &
                -( &
                cxdi*(-fin(1,i,j,k)  +fin(1,i,j+1,k))+ &
                cxd*(-fin(1,i+1,j,k)+fin(1,i+1,j+1,k))) &
                +( &
                cxdi*(-fin(1,i,j,k+1)  +fin(1,i,j+1,k+1))+ &
                cxd*(-fin(1,i+1,j,k+1)+fin(1,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hxi(v)*hy(v)*hzi(v)*( &
                -( &
                -(cydi*fin(2,i,j,k)  +cyd*fin(2,i,j+1,k)) &
                +(cydi*fin(2,i+1,j,k)+cyd*fin(2,i+1,j+1,k))) &
                +( &
                -(cydi*fin(2,i,j,k+1)  +cyd*fin(2,i,j+1,k+1)) &
                +(cydi*fin(2,i+1,j,k+1)+cyd*fin(2,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hxi(v)*hyi(v)*hz(v)*( &
                czdi*( &
                (fin(3,i,j,k)  -fin(3,i,j+1,k))- &
                (fin(3,i+1,j,k)-fin(3,i+1,j+1,k))) &
                +czd*( &
                (fin(3,i,j,k+1)  -fin(3,i,j+1,k+1))- &
                (fin(3,i+1,j,k+1)-fin(3,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hx(v)*hy(v)*hzi(v)*( &
                -( &
                cxdi*(cydi*fin(4,i,j,k)  +cyd*fin(4,i,j+1,k))+ &
                cxd*(cydi*fin(4,i+1,j,k)+cyd*fin(4,i+1,j+1,k))) &
                +( &
                cxdi*(cydi*fin(4,i,j,k+1)  +cyd*fin(4,i,j+1,k+1))+ &
                cxd*(cydi*fin(4,i+1,j,k+1)+cyd*fin(4,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hyi(v)*hx(v)*hz(v)*( &
                czdi*( &
                cxdi*(-fin(5,i,j,k)  +fin(5,i,j+1,k))+ &
                cxd*(-fin(5,i+1,j,k)+fin(5,i+1,j+1,k))) &
                +czd*( &
                cxdi*(-fin(5,i,j,k+1)  +fin(5,i,j+1,k+1))+ &
                cxd*(-fin(5,i+1,j,k+1)+fin(5,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hxi(v)*hy(v)*hz(v)*( &
                czdi*( &
                -(cydi*fin(6,i,j,k)  +cyd*fin(6,i,j+1,k)) &
                +(cydi*fin(6,i+1,j,k)+cyd*fin(6,i+1,j+1,k))) &
                +czd*( &
                -(cydi*fin(6,i,j,k+1)  +cyd*fin(6,i,j+1,k+1)) &
                +(cydi*fin(6,i+1,j,k+1)+cyd*fin(6,i+1,j+1,k+1))))
           !
           sum=sum+z216th*hx(v)*hy(v)*hz(v)*( &
                czdi*( &
                cxdi*(cydi*fin(7,i,j,k)  +cyd*fin(7,i,j+1,k))+ &
                cxd*(cydi*fin(7,i+1,j,k)+cyd*fin(7,i+1,j+1,k))) &
                +czd*( &
                cxdi*(cydi*fin(7,i,j,k+1)  +cyd*fin(7,i,j+1,k+1))+ &
                cxd*(cydi*fin(7,i+1,j,k+1)+cyd*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(6).eq.1) then
        !                               ! fxzz
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           j=jj(v)
           k=kk(v)
           !
           !   ...in x direction
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           xp2=xp*xp
           xpi2=xpi*xpi
           !
           cxd=3.0_fp*xp2-1.0_fp
           cxdi=-3.0_fp*xpi2+1.0_fp
           !
           !   ...and in y direction
           !
           yp=yparam(v)
           ypi=1.0_fp-yp
           yp2=yp*yp
           ypi2=ypi*ypi
           !
           cy=yp*(yp2-1.0_fp)
           cyi=ypi*(ypi2-1.0_fp)
           hy2=hy(v)*hy(v)

           !
           !   ...and in z direction
           !
           zp=zparam(v)
           zpi=1.0_fp-zp
           !
           sum=hxi(v)*( &
                zpi*( &
                -(ypi*fin(3,i,j,k)  +yp*fin(3,i,j+1,k)) &
                +(ypi*fin(3,i+1,j,k)+yp*fin(3,i+1,j+1,k))) &
                +zp*( &
                -(ypi*fin(3,i,j,k+1)  +yp*fin(3,i,j+1,k+1)) &
                +(ypi*fin(3,i+1,j,k+1)+yp*fin(3,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hx(v)*( &
                zpi*( &
                cxdi*(ypi*fin(5,i,j,k)  +yp*fin(5,i,j+1,k))+ &
                cxd*(ypi*fin(5,i+1,j,k)+yp*fin(5,i+1,j+1,k))) &
                +zp*( &
                cxdi*(ypi*fin(5,i,j,k+1)  +yp*fin(5,i,j+1,k+1))+ &
                cxd*(ypi*fin(5,i+1,j,k+1)+yp*fin(5,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hy2*hxi(v)*( &
                zpi*( &
                -(cyi*fin(6,i,j,k)  +cy*fin(6,i,j+1,k)) &
                +(cyi*fin(6,i+1,j,k)+cy*fin(6,i+1,j+1,k))) &
                +zp*( &
                -(cyi*fin(6,i,j,k+1)  +cy*fin(6,i,j+1,k+1)) &
                +(cyi*fin(6,i+1,j,k+1)+cy*fin(6,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hx(v)*hy2*( &
                zpi*( &
                cxdi*(cyi*fin(7,i,j,k)  +cy*fin(7,i,j+1,k))+ &
                cxd*(cyi*fin(7,i+1,j,k)+cy*fin(7,i+1,j+1,k))) &
                +zp*( &
                cxdi*(cyi*fin(7,i,j,k+1)  +cy*fin(7,i,j+1,k+1))+ &
                cxd*(cyi*fin(7,i+1,j,k+1)+cy*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(7).eq.1) then
        !                               ! fyyz
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           j=jj(v)
           k=kk(v)
           !
           !   ...in x direction
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           xp2=xp*xp
           xpi2=xpi*xpi
           !
           cx=xp*(xp2-1.0_fp)
           cxi=xpi*(xpi2-1.0_fp)
           hx2=hx(v)*hx(v)
           !
           !   ...and in y direction
           !
           yp=yparam(v)
           ypi=1.0_fp-yp
           !
           !   ...and in z direction
           !
           zp=zparam(v)
           zpi=1.0_fp-zp
           zp2=zp*zp
           zpi2=zpi*zpi
           !
           czd=3.0_fp*zp2-1.0_fp
           czdi=-3.0_fp*zpi2+1.0_fp
           !
           sum=hzi(v)*( &
                -( &
                xpi*(ypi*fin(2,i,j,k)  +yp*fin(2,i,j+1,k))+ &
                xp*(ypi*fin(2,i+1,j,k)+yp*fin(2,i+1,j+1,k))) &
                +( &
                xpi*(ypi*fin(2,i,j,k+1)  +yp*fin(2,i,j+1,k+1))+ &
                xp*(ypi*fin(2,i+1,j,k+1)+yp*fin(2,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hx2*hzi(v)*( &
                -( &
                cxi*(ypi*fin(4,i,j,k)  +yp*fin(4,i,j+1,k))+ &
                cx*(ypi*fin(4,i+1,j,k)+yp*fin(4,i+1,j+1,k))) &
                +( &
                cxi*(ypi*fin(4,i,j,k+1)  +yp*fin(4,i,j+1,k+1))+ &
                cx*(ypi*fin(4,i+1,j,k+1)+yp*fin(4,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hz(v)*( &
                czdi*( &
                xpi*(ypi*fin(6,i,j,k)  +yp*fin(6,i,j+1,k))+ &
                xp*(ypi*fin(6,i+1,j,k)+yp*fin(6,i+1,j+1,k))) &
                +czd*( &
                xpi*(ypi*fin(6,i,j,k+1)  +yp*fin(6,i,j+1,k+1))+ &
                xp*(ypi*fin(6,i+1,j,k+1)+yp*fin(6,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hx2*hz(v)*( &
                czdi*( &
                cxi*(ypi*fin(7,i,j,k)  +yp*fin(7,i,j+1,k))+ &
                cx*(ypi*fin(7,i+1,j,k)+yp*fin(7,i+1,j+1,k))) &
                +czd*( &
                cxi*(ypi*fin(7,i,j,k+1)  +yp*fin(7,i,j+1,k+1))+ &
                cx*(ypi*fin(7,i+1,j,k+1)+yp*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(8).eq.1) then
        !                               ! fyzz
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           j=jj(v)
           k=kk(v)
           !
           !   ...in x direction
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           xp2=xp*xp
           xpi2=xpi*xpi
           !
           cx=xp*(xp2-1.0_fp)
           cxi=xpi*(xpi2-1.0_fp)
           hx2=hx(v)*hx(v)
           !
           !   ...and in y direction
           !
           yp=yparam(v)
           ypi=1.0_fp-yp
           yp2=yp*yp
           ypi2=ypi*ypi
           !
           cyd=3.0_fp*yp2-1.0_fp
           cydi=-3.0_fp*ypi2+1.0_fp
           !
           !   ...and in z direction
           !
           zp=zparam(v)
           zpi=1.0_fp-zp
           !
           sum=hyi(v)*( &
                zpi*( &
                xpi*(-fin(3,i,j,k)  +fin(3,i,j+1,k))+ &
                xp*( -fin(3,i+1,j,k)+fin(3,i+1,j+1,k))) &
                +zp*( &
                xpi*(-fin(3,i,j,k+1)  +fin(3,i,j+1,k+1))+ &
                xp*( -fin(3,i+1,j,k+1)+fin(3,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hx2*hyi(v)*( &
                zpi*( &
                cxi*(-fin(5,i,j,k)  +fin(5,i,j+1,k))+ &
                cx*( -fin(5,i+1,j,k)+fin(5,i+1,j+1,k))) &
                +zp*( &
                cxi*(-fin(5,i,j,k+1)  +fin(5,i,j+1,k+1))+ &
                cx*( -fin(5,i+1,j,k+1)+fin(5,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hy(v)*( &
                zpi*( &
                xpi*(cydi*fin(6,i,j,k) +cyd*fin(6,i,j+1,k))+ &
                xp*(cydi*fin(6,i+1,j,k)+cyd*fin(6,i+1,j+1,k))) &
                +zp*( &
                xpi*(cydi*fin(6,i,j,k+1) +cyd*fin(6,i,j+1,k+1))+ &
                xp*(cydi*fin(6,i+1,j,k+1)+cyd*fin(6,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hx2*hy(v)*( &
                zpi*( &
                cxi*(cydi*fin(7,i,j,k) +cyd*fin(7,i,j+1,k))+ &
                cx*(cydi*fin(7,i+1,j,k)+cyd*fin(7,i+1,j+1,k))) &
                +zp*( &
                cxi*(cydi*fin(7,i,j,k+1) +cyd*fin(7,i,j+1,k+1))+ &
                cx*(cydi*fin(7,i+1,j,k+1)+cyd*fin(7,i+1,j+1,k+1))))
           fval(v,iadr)=sum
        end do
     end if
     !
     !----------------------------------
     !  3rd derivatives (3 in each coordinate)
     !
  else if(ict(1).eq.-3) then
     if(ict(2).eq.1) then
        !                               ! fxxx
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           j=jj(v)
           k=kk(v)
           !
           !   ...and in y direction
           !
           yp=yparam(v)
           ypi=1.0_fp-yp
           yp2=yp*yp
           ypi2=ypi*ypi
           !
           cy=yp*(yp2-1.0_fp)
           cyi=ypi*(ypi2-1.0_fp)
           hy2=hy(v)*hy(v)
           !
           !   ...and in z direction
           !
           zp=zparam(v)
           zpi=1.0_fp-zp
           zp2=zp*zp
           zpi2=zpi*zpi
           !
           cz=zp*(zp2-1.0_fp)
           czi=zpi*(zpi2-1.0_fp)
           hz2=hz(v)*hz(v)
           !
           sum=hxi(v)*( &
                zpi*( &
                -(ypi*fin(1,i,j,k)  +yp*fin(1,i,j+1,k)) &
                +(ypi*fin(1,i+1,j,k)+yp*fin(1,i+1,j+1,k))) &
                +zp*( &
                -(ypi*fin(1,i,j,k+1)  +yp*fin(1,i,j+1,k+1)) &
                +(ypi*fin(1,i+1,j,k+1)+yp*fin(1,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hy2*hxi(v)*( &
                zpi*( &
                -(cyi*fin(4,i,j,k)  +cy*fin(4,i,j+1,k)) &
                +(cyi*fin(4,i+1,j,k)+cy*fin(4,i+1,j+1,k))) &
                +zp*( &
                -(cyi*fin(4,i,j,k+1)  +cy*fin(4,i,j+1,k+1)) &
                +(cyi*fin(4,i+1,j,k+1)+cy*fin(4,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hz2*hxi(v)*( &
                czi*( &
                -(ypi*fin(5,i,j,k)  +yp*fin(5,i,j+1,k)) &
                +(ypi*fin(5,i+1,j,k)+yp*fin(5,i+1,j+1,k))) &
                +cz*( &
                -(ypi*fin(5,i,j,k+1)  +yp*fin(5,i,j+1,k+1)) &
                +(ypi*fin(5,i+1,j,k+1)+yp*fin(5,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hy2*hz2*hxi(v)*( &
                czi*( &
                -(cyi*fin(7,i,j,k)  +cy*fin(7,i,j+1,k)) &
                +(cyi*fin(7,i+1,j,k)+cy*fin(7,i+1,j+1,k))) &
                +cz*( &
                -(cyi*fin(7,i,j,k+1)  +cy*fin(7,i,j+1,k+1)) &
                +(cyi*fin(7,i+1,j,k+1)+cy*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(3).eq.1) then
        !                               ! fyyy
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           j=jj(v)
           k=kk(v)
           !
           !   ...in x direction
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           xp2=xp*xp
           xpi2=xpi*xpi
           !
           cx=xp*(xp2-1.0_fp)
           cxi=xpi*(xpi2-1.0_fp)
           hx2=hx(v)*hx(v)
           !
           !   ...and in z direction
           !
           zp=zparam(v)
           zpi=1.0_fp-zp
           zp2=zp*zp
           zpi2=zpi*zpi
           !
           cz=zp*(zp2-1.0_fp)
           czi=zpi*(zpi2-1.0_fp)
           hz2=hz(v)*hz(v)
           !
           sum=hyi(v)*( &
                zpi*( &
                xpi*(-fin(2,i,j,k)  +fin(2,i,j+1,k))+ &
                xp*( -fin(2,i+1,j,k)+fin(2,i+1,j+1,k))) &
                +zp*( &
                xpi*(-fin(2,i,j,k+1)  +fin(2,i,j+1,k+1))+ &
                xp*( -fin(2,i+1,j,k+1)+fin(2,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hx2*hyi(v)*( &
                zpi*( &
                cxi*(-fin(4,i,j,k)  +fin(4,i,j+1,k))+ &
                cx*( -fin(4,i+1,j,k)+fin(4,i+1,j+1,k))) &
                +zp*( &
                cxi*(-fin(4,i,j,k+1)  +fin(4,i,j+1,k+1))+ &
                cx*( -fin(4,i+1,j,k+1)+fin(4,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hz2*hyi(v)*( &
                czi*( &
                xpi*(-fin(6,i,j,k)  +fin(6,i,j+1,k))+ &
                xp*( -fin(6,i+1,j,k)+fin(6,i+1,j+1,k))) &
                +cz*( &
                xpi*(-fin(6,i,j,k+1)  +fin(6,i,j+1,k+1))+ &
                xp*( -fin(6,i+1,j,k+1)+fin(6,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hx2*hz2*hyi(v)*( &
                czi*( &
                cxi*(-fin(7,i,j,k)  +fin(7,i,j+1,k))+ &
                cx*( -fin(7,i+1,j,k)+fin(7,i+1,j+1,k))) &
                +cz*( &
                cxi*(-fin(7,i,j,k+1)  +fin(7,i,j+1,k+1))+ &
                cx*( -fin(7,i+1,j,k+1)+fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(4).eq.1) then
        !                               ! fzzz
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           j=jj(v)
           k=kk(v)
           !
           !   ...in x direction
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           xp2=xp*xp
           xpi2=xpi*xpi
           !
           cx=xp*(xp2-1.0_fp)
           cxi=xpi*(xpi2-1.0_fp)
           hx2=hx(v)*hx(v)
           !
           !   ...and in y direction
           !
           yp=yparam(v)
           ypi=1.0_fp-yp
           yp2=yp*yp
           ypi2=ypi*ypi
           !
           cy=yp*(yp2-1.0_fp)
           cyi=ypi*(ypi2-1.0_fp)
           hy2=hy(v)*hy(v)
           !
           sum=hzi(v)*( &
                -( &
                xpi*(ypi*fin(3,i,j,k)  +yp*fin(3,i,j+1,k))+ &
                xp*(ypi*fin(3,i+1,j,k)+yp*fin(3,i+1,j+1,k))) &
                +( &
                xpi*(ypi*fin(3,i,j,k+1)  +yp*fin(3,i,j+1,k+1))+ &
                xp*(ypi*fin(3,i+1,j,k+1)+yp*fin(3,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hx2*hzi(v)*( &
                -( &
                cxi*(ypi*fin(5,i,j,k)  +yp*fin(5,i,j+1,k))+ &
                cx*(ypi*fin(5,i+1,j,k)+yp*fin(5,i+1,j+1,k))) &
                +( &
                cxi*(ypi*fin(5,i,j,k+1)  +yp*fin(5,i,j+1,k+1))+ &
                cx*(ypi*fin(5,i+1,j,k+1)+yp*fin(5,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hy2*hzi(v)*( &
                -( &
                xpi*(cyi*fin(6,i,j,k)  +cy*fin(6,i,j+1,k))+ &
                xp*(cyi*fin(6,i+1,j,k)+cy*fin(6,i+1,j+1,k))) &
                +( &
                xpi*(cyi*fin(6,i,j,k+1)  +cy*fin(6,i,j+1,k+1))+ &
                xp*(cyi*fin(6,i+1,j,k+1)+cy*fin(6,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hx2*hy2*hzi(v)*( &
                -( &
                cxi*(cyi*fin(7,i,j,k)  +cy*fin(7,i,j+1,k))+ &
                cx*(cyi*fin(7,i+1,j,k)+cy*fin(7,i+1,j+1,k))) &
                +( &
                cxi*(cyi*fin(7,i,j,k+1)  +cy*fin(7,i,j+1,k+1))+ &
                cx*(cyi*fin(7,i+1,j,k+1)+cy*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     !----------------------------------
     !  4th derivatives (.le.2 in each coordinate)
     !
  else if(ict(1).eq.4) then
     if(ict(2).eq.1) then
        !                               ! fxxyy
        iadr=iadr+1
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
           !   ...and in y direction
           !
           yp=yparam(v)
           ypi=1.0_fp-yp
           !
           !   ...and in z direction
           !
           zp=zparam(v)
           zpi=1.0_fp-zp
           zp2=zp*zp
           zpi2=zpi*zpi
           !
           cz=zp*(zp2-1.0_fp)
           czi=zpi*(zpi2-1.0_fp)
           hz2=hz(v)*hz(v)
           !
           sum=( &
                zpi*( &
                xpi*(ypi*fin(4,i,j,k)  +yp*fin(4,i,j+1,k))+ &
                xp*( ypi*fin(4,i+1,j,k)+yp*fin(4,i+1,j+1,k))) &
                +zp*( &
                xpi*(ypi*fin(4,i,j,k+1) +yp*fin(4,i,j+1,k+1))+ &
                xp*(ypi*fin(4,i+1,j,k+1)+yp*fin(4,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hz2*( &
                czi*( &
                xpi*(ypi*fin(7,i,j,k) +yp*fin(7,i,j+1,k))+ &
                xp*(ypi*fin(7,i+1,j,k)+yp*fin(7,i+1,j+1,k))) &
                +cz*( &
                xpi*(ypi*fin(7,i,j,k+1) +yp*fin(7,i,j+1,k+1))+ &
                xp*(ypi*fin(7,i+1,j,k+1)+yp*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(3).eq.1) then
        !                               ! fxxyz
        iadr=iadr+1
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
           !   ...and in y direction
           !
           yp=yparam(v)
           ypi=1.0_fp-yp
           yp2=yp*yp
           ypi2=ypi*ypi
           !
           cyd=3.0_fp*yp2-1.0_fp
           cydi=-3.0_fp*ypi2+1.0_fp
           !
           !   ...and in z direction
           !
           zp=zparam(v)
           zpi=1.0_fp-zp
           zp2=zp*zp
           zpi2=zpi*zpi
           !
           czd=3.0_fp*zp2-1.0_fp
           czdi=-3.0_fp*zpi2+1.0_fp
           !
           sum=hyi(v)*hzi(v)*( &
                -( &
                xpi*(-fin(1,i,j,k)  +fin(1,i,j+1,k))+ &
                xp*( -fin(1,i+1,j,k)+fin(1,i+1,j+1,k))) &
                +( &
                xpi*(-fin(1,i,j,k+1)  +fin(1,i,j+1,k+1))+ &
                xp*( -fin(1,i+1,j,k+1)+fin(1,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hy(v)*hzi(v)*( &
                -( &
                xpi*(cydi*fin(4,i,j,k)  +cyd*fin(4,i,j+1,k))+ &
                xp*(cydi*fin(4,i+1,j,k)+cyd*fin(4,i+1,j+1,k))) &
                +( &
                xpi*(cydi*fin(4,i,j,k+1)  +cyd*fin(4,i,j+1,k+1))+ &
                xp*(cydi*fin(4,i+1,j,k+1)+cyd*fin(4,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hz(v)*hyi(v)*( &
                czdi*( &
                xpi*(-fin(5,i,j,k)  +fin(5,i,j+1,k))+ &
                xp*( -fin(5,i+1,j,k)+fin(5,i+1,j+1,k))) &
                +czd*( &
                xpi*(-fin(5,i,j,k+1)  +fin(5,i,j+1,k+1))+ &
                xp*( -fin(5,i+1,j,k+1)+fin(5,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hy(v)*hz(v)*( &
                czdi*( &
                xpi*(cydi*fin(7,i,j,k)  +cyd*fin(7,i,j+1,k))+ &
                xp*(cydi*fin(7,i+1,j,k)+cyd*fin(7,i+1,j+1,k))) &
                +czd*( &
                xpi*(cydi*fin(7,i,j,k+1)  +cyd*fin(7,i,j+1,k+1))+ &
                xp*(cydi*fin(7,i+1,j,k+1)+cyd*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(4).eq.1) then
        !                               ! fxxzz
        iadr=iadr+1
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
           !   ...and in y direction
           !
           yp=yparam(v)
           ypi=1.0_fp-yp
           yp2=yp*yp
           ypi2=ypi*ypi
           !
           cy=yp*(yp2-1.0_fp)
           cyi=ypi*(ypi2-1.0_fp)
           hy2=hy(v)*hy(v)
           !
           !   ...and in z direction
           !
           zp=zparam(v)
           zpi=1.0_fp-zp
           !
           sum=( &
                zpi*( &
                xpi*(ypi*fin(5,i,j,k)  +yp*fin(5,i,j+1,k))+ &
                xp*(ypi*fin(5,i+1,j,k)+yp*fin(5,i+1,j+1,k))) &
                +zp*( &
                xpi*(ypi*fin(5,i,j,k+1)  +yp*fin(5,i,j+1,k+1))+ &
                xp*(ypi*fin(5,i+1,j,k+1)+yp*fin(5,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hy2*( &
                zpi*( &
                xpi*(cyi*fin(7,i,j,k)  +cy*fin(7,i,j+1,k))+ &
                xp*(cyi*fin(7,i+1,j,k)+cy*fin(7,i+1,j+1,k))) &
                +zp*( &
                xpi*(cyi*fin(7,i,j,k+1)  +cy*fin(7,i,j+1,k+1))+ &
                xp*(cyi*fin(7,i+1,j,k+1)+cy*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(5).eq.1) then
        !                               ! fxyyz
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           j=jj(v)
           k=kk(v)
           !
           !   ...in x direction
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           xp2=xp*xp
           xpi2=xpi*xpi
           !
           cxd=3.0_fp*xp2-1.0_fp
           cxdi=-3.0_fp*xpi2+1.0_fp
           !
           !   ...and in y direction
           !
           yp=yparam(v)
           ypi=1.0_fp-yp
           !
           !   ...and in z direction
           !
           zp=zparam(v)
           zpi=1.0_fp-zp
           zp2=zp*zp
           zpi2=zpi*zpi

           czd=3.0_fp*zp2-1.0_fp
           czdi=-3.0_fp*zpi2+1.0_fp
           !
           sum=hxi(v)*hzi(v)*( &
                -( &
                -(ypi*fin(2,i,j,k)  +yp*fin(2,i,j+1,k)) &
                +(ypi*fin(2,i+1,j,k)+yp*fin(2,i+1,j+1,k))) &
                +( &
                -(ypi*fin(2,i,j,k+1)  +yp*fin(2,i,j+1,k+1)) &
                +(ypi*fin(2,i+1,j,k+1)+yp*fin(2,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hx(v)*hzi(v)*( &
                -( &
                cxdi*(ypi*fin(4,i,j,k)  +yp*fin(4,i,j+1,k))+ &
                cxd*(ypi*fin(4,i+1,j,k)+yp*fin(4,i+1,j+1,k))) &
                +( &
                cxdi*(ypi*fin(4,i,j,k+1)  +yp*fin(4,i,j+1,k+1))+ &
                cxd*(ypi*fin(4,i+1,j,k+1)+yp*fin(4,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hz(v)*hxi(v)*( &
                czdi*( &
                -(ypi*fin(6,i,j,k)  +yp*fin(6,i,j+1,k)) &
                +(ypi*fin(6,i+1,j,k)+yp*fin(6,i+1,j+1,k))) &
                +czd*( &
                -(ypi*fin(6,i,j,k+1)  +yp*fin(6,i,j+1,k+1)) &
                +(ypi*fin(6,i+1,j,k+1)+yp*fin(6,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hx(v)*hz(v)*( &
                czdi*( &
                cxdi*(ypi*fin(7,i,j,k)  +yp*fin(7,i,j+1,k))+ &
                cxd*(ypi*fin(7,i+1,j,k)+yp*fin(7,i+1,j+1,k))) &
                +czd*( &
                cxdi*(ypi*fin(7,i,j,k+1)  +yp*fin(7,i,j+1,k+1))+ &
                cxd*(ypi*fin(7,i+1,j,k+1)+yp*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(6).eq.1) then
        !                               ! fxyzz
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           j=jj(v)
           k=kk(v)
           !
           !   ...in x direction
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           xp2=xp*xp
           xpi2=xpi*xpi

           cxd=3.0_fp*xp2-1.0_fp
           cxdi=-3.0_fp*xpi2+1.0_fp
           !
           !   ...and in y direction
           !
           yp=yparam(v)
           ypi=1.0_fp-yp
           yp2=yp*yp
           ypi2=ypi*ypi

           cyd=3.0_fp*yp2-1.0_fp
           cydi=-3.0_fp*ypi2+1.0_fp
           !
           !   ...and in z direction
           !
           zp=zparam(v)
           zpi=1.0_fp-zp
           !
           sum=hyi(v)*hxi(v)*( &
                zpi*( &
                ( +fin(3,i,j,k)  -fin(3,i,j+1,k)) &
                +( -fin(3,i+1,j,k)+fin(3,i+1,j+1,k))) &
                +zp*( &
                ( +fin(3,i,j,k+1)  -fin(3,i,j+1,k+1)) &
                +( -fin(3,i+1,j,k+1)+fin(3,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hx(v)*hyi(v)*( &
                zpi*( &
                cxdi*(-fin(5,i,j,k)  +fin(5,i,j+1,k))+ &
                cxd*( -fin(5,i+1,j,k)+fin(5,i+1,j+1,k))) &
                +zp*( &
                cxdi*(-fin(5,i,j,k+1)  +fin(5,i,j+1,k+1))+ &
                cxd*( -fin(5,i+1,j,k+1)+fin(5,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hy(v)*hxi(v)*( &
                zpi*( &
                -(cydi*fin(6,i,j,k) +cyd*fin(6,i,j+1,k)) &
                +(cydi*fin(6,i+1,j,k)+cyd*fin(6,i+1,j+1,k))) &
                +zp*( &
                -(cydi*fin(6,i,j,k+1) +cyd*fin(6,i,j+1,k+1)) &
                +(cydi*fin(6,i+1,j,k+1)+cyd*fin(6,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hx(v)*hy(v)*( &
                zpi*( &
                cxdi*(cydi*fin(7,i,j,k) +cyd*fin(7,i,j+1,k))+ &
                cxd*(cydi*fin(7,i+1,j,k)+cyd*fin(7,i+1,j+1,k))) &
                +zp*( &
                cxdi*(cydi*fin(7,i,j,k+1) +cyd*fin(7,i,j+1,k+1))+ &
                cxd*(cydi*fin(7,i+1,j,k+1)+cyd*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(7).eq.1) then
        !                               ! fyyzz
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           j=jj(v)
           k=kk(v)
           !
           !   ...in x direction
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           xp2=xp*xp
           xpi2=xpi*xpi
           !
           cx=xp*(xp2-1.0_fp)
           cxi=xpi*(xpi2-1.0_fp)
           hx2=hx(v)*hx(v)
           !
           !   ...and in y direction
           !
           yp=yparam(v)
           ypi=1.0_fp-yp
           !
           !   ...and in z direction
           !
           zp=zparam(v)
           zpi=1.0_fp-zp
           !
           sum=( &
                zpi*( &
                xpi*(ypi*fin(6,i,j,k)  +yp*fin(6,i,j+1,k))+ &
                xp*( ypi*fin(6,i+1,j,k)+yp*fin(6,i+1,j+1,k))) &
                +zp*( &
                xpi*(ypi*fin(6,i,j,k+1) +yp*fin(6,i,j+1,k+1))+ &
                xp*(ypi*fin(6,i+1,j,k+1)+yp*fin(6,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hx2*( &
                zpi*( &
                cxi*(ypi*fin(7,i,j,k)  +yp*fin(7,i,j+1,k))+ &
                cx*( ypi*fin(7,i+1,j,k)+yp*fin(7,i+1,j+1,k))) &
                +zp*( &
                cxi*(ypi*fin(7,i,j,k+1)  +yp*fin(7,i,j+1,k+1))+ &
                cx*( ypi*fin(7,i+1,j,k+1)+yp*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     !----------------------------------
     !  4th derivatives (3 in a coordinate)
     !
  else if(ict(1).eq.-4) then
     if(ict(2).eq.1) then
        !                               ! fxxxy
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           j=jj(v)
           k=kk(v)
           !
           !   ...and in y direction
           !
           yp=yparam(v)
           ypi=1.0_fp-yp
           yp2=yp*yp
           ypi2=ypi*ypi

           cyd=3.0_fp*yp2-1.0_fp
           cydi=-3.0_fp*ypi2+1.0_fp
           !
           !   ...and in z direction
           !
           zp=zparam(v)
           zpi=1.0_fp-zp
           zp2=zp*zp
           zpi2=zpi*zpi
           !
           cz=zp*(zp2-1.0_fp)
           czi=zpi*(zpi2-1.0_fp)
           hz2=hz(v)*hz(v)
           !
           sum=hyi(v)*hxi(v)*( &
                zpi*( &
                (  fin(1,i,j,k)  -fin(1,i,j+1,k))+ &
                ( -fin(1,i+1,j,k)+fin(1,i+1,j+1,k))) &
                +zp*( &
                (  fin(1,i,j,k+1)  -fin(1,i,j+1,k+1))+ &
                ( -fin(1,i+1,j,k+1)+fin(1,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hy(v)*hxi(v)*( &
                zpi*( &
                -(cydi*fin(4,i,j,k)  +cyd*fin(4,i,j+1,k))+ &
                (cydi*fin(4,i+1,j,k)+cyd*fin(4,i+1,j+1,k))) &
                +zp*( &
                -(cydi*fin(4,i,j,k+1)  +cyd*fin(4,i,j+1,k+1))+ &
                (cydi*fin(4,i+1,j,k+1)+cyd*fin(4,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hz2*hyi(v)*hxi(v)*( &
                czi*( &
                (  fin(5,i,j,k)  -fin(5,i,j+1,k))+ &
                ( -fin(5,i+1,j,k)+fin(5,i+1,j+1,k))) &
                +cz*( &
                (  fin(5,i,j,k+1)  -fin(5,i,j+1,k+1))+ &
                ( -fin(5,i+1,j,k+1)+fin(5,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hy(v)*hz2*hxi(v)*( &
                czi*( &
                -(cydi*fin(7,i,j,k)  +cyd*fin(7,i,j+1,k))+ &
                (cydi*fin(7,i+1,j,k)+cyd*fin(7,i+1,j+1,k))) &
                +cz*( &
                -(cydi*fin(7,i,j,k+1)  +cyd*fin(7,i,j+1,k+1))+ &
                (cydi*fin(7,i+1,j,k+1)+cyd*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(3).eq.1) then
        !                               ! fxxxz
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           j=jj(v)
           k=kk(v)
           !
           !   ...and in y direction
           !
           yp=yparam(v)
           ypi=1.0_fp-yp
           yp2=yp*yp
           ypi2=ypi*ypi
           !
           cy=yp*(yp2-1.0_fp)
           cyi=ypi*(ypi2-1.0_fp)
           hy2=hy(v)*hy(v)
           !
           !   ...and in z direction
           !
           zp=zparam(v)
           zpi=1.0_fp-zp
           zp2=zp*zp
           zpi2=zpi*zpi

           czd=3.0_fp*zp2-1.0_fp
           czdi=-3.0_fp*zpi2+1.0_fp
           !
           sum=hzi(v)*hxi(v)*( &
                ( &
                +(ypi*fin(1,i,j,k)  +yp*fin(1,i,j+1,k)) &
                -(ypi*fin(1,i+1,j,k)+yp*fin(1,i+1,j+1,k))) &
                +( &
                -(ypi*fin(1,i,j,k+1)  +yp*fin(1,i,j+1,k+1)) &
                +(ypi*fin(1,i+1,j,k+1)+yp*fin(1,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hy2*hzi(v)*hxi(v)*( &
                ( &
                +(cyi*fin(4,i,j,k)  +cy*fin(4,i,j+1,k)) &
                -(cyi*fin(4,i+1,j,k)+cy*fin(4,i+1,j+1,k))) &
                +( &
                -(cyi*fin(4,i,j,k+1)  +cy*fin(4,i,j+1,k+1)) &
                +(cyi*fin(4,i+1,j,k+1)+cy*fin(4,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hz(v)*hxi(v)*( &
                czdi*( &
                -(ypi*fin(5,i,j,k)  +yp*fin(5,i,j+1,k)) &
                +(ypi*fin(5,i+1,j,k)+yp*fin(5,i+1,j+1,k))) &
                +czd*( &
                -(ypi*fin(5,i,j,k+1)  +yp*fin(5,i,j+1,k+1)) &
                +(ypi*fin(5,i+1,j,k+1)+yp*fin(5,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hy2*hz(v)*hxi(v)*( &
                czdi*( &
                -(cyi*fin(7,i,j,k)  +cy*fin(7,i,j+1,k)) &
                +(cyi*fin(7,i+1,j,k)+cy*fin(7,i+1,j+1,k))) &
                +czd*( &
                -(cyi*fin(7,i,j,k+1)  +cy*fin(7,i,j+1,k+1)) &
                +(cyi*fin(7,i+1,j,k+1)+cy*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(4).eq.1) then
        !                               ! fxyyy
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           j=jj(v)
           k=kk(v)
           !
           !   ...in x direction
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           xp2=xp*xp
           xpi2=xpi*xpi

           cxd=3.0_fp*xp2-1.0_fp
           cxdi=-3.0_fp*xpi2+1.0_fp
           !
           !   ...and in z direction
           !
           zp=zparam(v)
           zpi=1.0_fp-zp
           zp2=zp*zp
           zpi2=zpi*zpi
           !
           cz=zp*(zp2-1.0_fp)
           czi=zpi*(zpi2-1.0_fp)
           hz2=hz(v)*hz(v)
           !
           sum=hxi(v)*hyi(v)*( &
                zpi*( &
                ( fin(2,i,j,k)  -fin(2,i,j+1,k)) &
                +(-fin(2,i+1,j,k)+fin(2,i+1,j+1,k))) &
                +zp*( &
                ( fin(2,i,j,k+1)  -fin(2,i,j+1,k+1)) &
                +(-fin(2,i+1,j,k+1)+fin(2,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hx(v)*hyi(v)*( &
                zpi*( &
                cxdi*(-fin(4,i,j,k)  +fin(4,i,j+1,k))+ &
                cxd*( -fin(4,i+1,j,k)+fin(4,i+1,j+1,k))) &
                +zp*( &
                cxdi*(-fin(4,i,j,k+1)  +fin(4,i,j+1,k+1))+ &
                cxd*( -fin(4,i+1,j,k+1)+fin(4,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hz2*hxi(v)*hyi(v)*( &
                czi*( &
                ( fin(6,i,j,k)  -fin(6,i,j+1,k)) &
                +(-fin(6,i+1,j,k)+fin(6,i+1,j+1,k))) &
                +cz*( &
                ( fin(6,i,j,k+1)  -fin(6,i,j+1,k+1)) &
                +(-fin(6,i+1,j,k+1)+fin(6,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hx(v)*hz2*hyi(v)*( &
                czi*( &
                cxdi*(-fin(7,i,j,k)  +fin(7,i,j+1,k))+ &
                cxd*( -fin(7,i+1,j,k)+fin(7,i+1,j+1,k))) &
                +cz*( &
                cxdi*(-fin(7,i,j,k+1)  +fin(7,i,j+1,k+1))+ &
                cxd*( -fin(7,i+1,j,k+1)+fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(5).eq.1) then
        !                               ! fxzzz
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           j=jj(v)
           k=kk(v)
           !
           !   ...in x direction
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           xp2=xp*xp
           xpi2=xpi*xpi

           cxd=3.0_fp*xp2-1.0_fp
           cxdi=-3.0_fp*xpi2+1.0_fp
           !
           !   ...and in y direction
           !
           yp=yparam(v)
           ypi=1.0_fp-yp
           yp2=yp*yp
           ypi2=ypi*ypi
           !
           cy=yp*(yp2-1.0_fp)
           cyi=ypi*(ypi2-1.0_fp)
           hy2=hy(v)*hy(v)
           !
           sum=hxi(v)*hzi(v)*( &
                -( &
                -(ypi*fin(3,i,j,k)  +yp*fin(3,i,j+1,k)) &
                +(ypi*fin(3,i+1,j,k)+yp*fin(3,i+1,j+1,k))) &
                +( &
                -(ypi*fin(3,i,j,k+1)  +yp*fin(3,i,j+1,k+1)) &
                +(ypi*fin(3,i+1,j,k+1)+yp*fin(3,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hx(v)*hzi(v)*( &
                -( &
                cxdi*(ypi*fin(5,i,j,k)  +yp*fin(5,i,j+1,k))+ &
                cxd*(ypi*fin(5,i+1,j,k)+yp*fin(5,i+1,j+1,k))) &
                +( &
                cxdi*(ypi*fin(5,i,j,k+1)  +yp*fin(5,i,j+1,k+1))+ &
                cxd*(ypi*fin(5,i+1,j,k+1)+yp*fin(5,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hy2*hxi(v)*hzi(v)*( &
                -( &
                -(cyi*fin(6,i,j,k)  +cy*fin(6,i,j+1,k)) &
                +(cyi*fin(6,i+1,j,k)+cy*fin(6,i+1,j+1,k))) &
                +( &
                -(cyi*fin(6,i,j,k+1)  +cy*fin(6,i,j+1,k+1)) &
                +(cyi*fin(6,i+1,j,k+1)+cy*fin(6,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hx(v)*hy2*hzi(v)*( &
                -( &
                cxdi*(cyi*fin(7,i,j,k)  +cy*fin(7,i,j+1,k))+ &
                cxd*(cyi*fin(7,i+1,j,k)+cy*fin(7,i+1,j+1,k))) &
                +( &
                cxdi*(cyi*fin(7,i,j,k+1)  +cy*fin(7,i,j+1,k+1))+ &
                cxd*(cyi*fin(7,i+1,j,k+1)+cy*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(6).eq.1) then
        !                               ! fyyyz
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           j=jj(v)
           k=kk(v)
           !
           !   ...in x direction
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           xp2=xp*xp
           xpi2=xpi*xpi
           !
           cx=xp*(xp2-1.0_fp)
           cxi=xpi*(xpi2-1.0_fp)
           hx2=hx(v)*hx(v)
           !
           !   ...and in z direction
           !
           zp=zparam(v)
           zpi=1.0_fp-zp
           zp2=zp*zp
           zpi2=zpi*zpi

           czd=3.0_fp*zp2-1.0_fp
           czdi=-3.0_fp*zpi2+1.0_fp
           !
           sum=hyi(v)*hzi(v)*( &
                -( &
                xpi*(-fin(2,i,j,k)  +fin(2,i,j+1,k))+ &
                xp*( -fin(2,i+1,j,k)+fin(2,i+1,j+1,k))) &
                +( &
                xpi*(-fin(2,i,j,k+1)  +fin(2,i,j+1,k+1))+ &
                xp*( -fin(2,i+1,j,k+1)+fin(2,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hx2*hyi(v)*hzi(v)*( &
                -( &
                cxi*(-fin(4,i,j,k)  +fin(4,i,j+1,k))+ &
                cx*( -fin(4,i+1,j,k)+fin(4,i+1,j+1,k))) &
                +( &
                cxi*(-fin(4,i,j,k+1)  +fin(4,i,j+1,k+1))+ &
                cx*( -fin(4,i+1,j,k+1)+fin(4,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hz(v)*hyi(v)*( &
                czdi*( &
                xpi*(-fin(6,i,j,k)  +fin(6,i,j+1,k))+ &
                xp*( -fin(6,i+1,j,k)+fin(6,i+1,j+1,k))) &
                +czd*( &
                xpi*(-fin(6,i,j,k+1)  +fin(6,i,j+1,k+1))+ &
                xp*( -fin(6,i+1,j,k+1)+fin(6,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hx2*hz(v)*hyi(v)*( &
                czdi*( &
                cxi*(-fin(7,i,j,k)  +fin(7,i,j+1,k))+ &
                cx*( -fin(7,i+1,j,k)+fin(7,i+1,j+1,k))) &
                +czd*( &
                cxi*(-fin(7,i,j,k+1)  +fin(7,i,j+1,k+1))+ &
                cx*( -fin(7,i+1,j,k+1)+fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(7).eq.1) then
        !                               ! fyzzz
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           j=jj(v)
           k=kk(v)
           !
           !   ...in x direction
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           xp2=xp*xp
           xpi2=xpi*xpi
           !
           cx=xp*(xp2-1.0_fp)
           cxi=xpi*(xpi2-1.0_fp)
           hx2=hx(v)*hx(v)
           !
           !   ...and in y direction
           !
           yp=yparam(v)
           ypi=1.0_fp-yp
           yp2=yp*yp
           ypi2=ypi*ypi

           cyd=3.0_fp*yp2-1.0_fp
           cydi=-3.0_fp*ypi2+1.0_fp
           !
           sum=hyi(v)*hzi(v)*( &
                -( &
                xpi*(-fin(3,i,j,k)  +fin(3,i,j+1,k))+ &
                xp*( -fin(3,i+1,j,k)+fin(3,i+1,j+1,k))) &
                +( &
                xpi*(-fin(3,i,j,k+1)  +fin(3,i,j+1,k+1))+ &
                xp*( -fin(3,i+1,j,k+1)+fin(3,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hx2*hyi(v)*hzi(v)*( &
                -( &
                cxi*(-fin(5,i,j,k)  +fin(5,i,j+1,k))+ &
                cx*( -fin(5,i+1,j,k)+fin(5,i+1,j+1,k))) &
                +( &
                cxi*(-fin(5,i,j,k+1)  +fin(5,i,j+1,k+1))+ &
                cx*( -fin(5,i+1,j,k+1)+fin(5,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hy(v)*hzi(v)*( &
                -( &
                xpi*(cydi*fin(6,i,j,k) +cyd*fin(6,i,j+1,k))+ &
                xp*(cydi*fin(6,i+1,j,k)+cyd*fin(6,i+1,j+1,k))) &
                +( &
                xpi*(cydi*fin(6,i,j,k+1) +cyd*fin(6,i,j+1,k+1))+ &
                xp*(cydi*fin(6,i+1,j,k+1)+cyd*fin(6,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hx2*hy(v)*hzi(v)*( &
                -( &
                cxi*(cydi*fin(7,i,j,k) +cyd*fin(7,i,j+1,k))+ &
                cx*(cydi*fin(7,i+1,j,k)+cyd*fin(7,i+1,j+1,k))) &
                +( &
                cxi*(cydi*fin(7,i,j,k+1) +cyd*fin(7,i,j+1,k+1))+ &
                cx*(cydi*fin(7,i+1,j,k+1)+cyd*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     !----------------------------------
     !  5th derivatives (.le.2 in each coordinate)
     !
  else if(ict(1).eq.5) then
     if(ict(2).eq.1) then
        !                               ! fxxyyz
        iadr=iadr+1
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
           !   ...and in y direction
           !
           yp=yparam(v)
           ypi=1.0_fp-yp
           !
           !   ...and in z direction
           !
           zp=zparam(v)
           zpi=1.0_fp-zp
           zp2=zp*zp
           zpi2=zpi*zpi

           czd=3.0_fp*zp2-1.0_fp
           czdi=-3.0_fp*zpi2+1.0_fp
           !
           sum=hzi(v)*( &
                -( &
                xpi*(ypi*fin(4,i,j,k)  +yp*fin(4,i,j+1,k))+ &
                xp*( ypi*fin(4,i+1,j,k)+yp*fin(4,i+1,j+1,k))) &
                +( &
                xpi*(ypi*fin(4,i,j,k+1) +yp*fin(4,i,j+1,k+1))+ &
                xp*(ypi*fin(4,i+1,j,k+1)+yp*fin(4,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hz(v)*( &
                czdi*( &
                xpi*(ypi*fin(7,i,j,k) +yp*fin(7,i,j+1,k))+ &
                xp*(ypi*fin(7,i+1,j,k)+yp*fin(7,i+1,j+1,k))) &
                +czd*( &
                xpi*(ypi*fin(7,i,j,k+1) +yp*fin(7,i,j+1,k+1))+ &
                xp*(ypi*fin(7,i+1,j,k+1)+yp*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(3).eq.1) then
        !                               ! fxxyzz
        iadr=iadr+1
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
           !   ...and in y direction
           !
           yp=yparam(v)
           ypi=1.0_fp-yp
           yp2=yp*yp
           ypi2=ypi*ypi

           cyd=3.0_fp*yp2-1.0_fp
           cydi=-3.0_fp*ypi2+1.0_fp
           !
           !   ...and in z direction
           !
           zp=zparam(v)
           zpi=1.0_fp-zp
           !
           sum=hyi(v)*( &
                zpi*( &
                xpi*(-fin(5,i,j,k)  +fin(5,i,j+1,k))+ &
                xp*( -fin(5,i+1,j,k)+fin(5,i+1,j+1,k))) &
                +zp*( &
                xpi*(-fin(5,i,j,k+1)  +fin(5,i,j+1,k+1))+ &
                xp*( -fin(5,i+1,j,k+1)+fin(5,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hy(v)*( &
                zpi*( &
                xpi*(cydi*fin(7,i,j,k)  +cyd*fin(7,i,j+1,k))+ &
                xp*( cydi*fin(7,i+1,j,k)+cyd*fin(7,i+1,j+1,k))) &
                +zp*( &
                xpi*(cydi*fin(7,i,j,k+1)  +cyd*fin(7,i,j+1,k+1))+ &
                xp*( cydi*fin(7,i+1,j,k+1)+cyd*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(4).eq.1) then
        !                               ! fxyyzz
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           j=jj(v)
           k=kk(v)
           !
           !   ...in x direction
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           xp2=xp*xp
           xpi2=xpi*xpi

           cxd=3.0_fp*xp2-1.0_fp
           cxdi=-3.0_fp*xpi2+1.0_fp
           !
           !   ...and in y direction
           !
           yp=yparam(v)
           ypi=1.0_fp-yp
           !
           !   ...and in z direction
           !
           zp=zparam(v)
           zpi=1.0_fp-zp
           !
           sum=hxi(v)*( &
                zpi*( &
                -(ypi*fin(6,i,j,k)  +yp*fin(6,i,j+1,k)) &
                +(ypi*fin(6,i+1,j,k)+yp*fin(6,i+1,j+1,k))) &
                +zp*( &
                -(ypi*fin(6,i,j,k+1)  +yp*fin(6,i,j+1,k+1)) &
                +(ypi*fin(6,i+1,j,k+1)+yp*fin(6,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hx(v)*( &
                zpi*( &
                cxdi*(ypi*fin(7,i,j,k)  +yp*fin(7,i,j+1,k))+ &
                cxd*(ypi*fin(7,i+1,j,k) +yp*fin(7,i+1,j+1,k))) &
                +zp*( &
                cxdi*(ypi*fin(7,i,j,k+1)  +yp*fin(7,i,j+1,k+1))+ &
                cxd*( ypi*fin(7,i+1,j,k+1)+yp*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     !----------------------------------
     !  5th derivatives (3 in a coordinate)
     !
  else if(ict(1).eq.-5) then
     if(ict(2).eq.1) then
        !                               ! fxxxyy
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           j=jj(v)
           k=kk(v)
           !
           !   ...and in y direction
           !
           yp=yparam(v)
           ypi=1.0_fp-yp
           !
           !   ...and in z direction
           !
           zp=zparam(v)
           zpi=1.0_fp-zp
           zp2=zp*zp
           zpi2=zpi*zpi
           !
           cz=zp*(zp2-1.0_fp)
           czi=zpi*(zpi2-1.0_fp)
           hz2=hz(v)*hz(v)
           !
           sum=hxi(v)*( &
                zpi*( &
                -(ypi*fin(4,i,j,k)  +yp*fin(4,i,j+1,k)) &
                +( ypi*fin(4,i+1,j,k)+yp*fin(4,i+1,j+1,k))) &
                +zp*( &
                -(ypi*fin(4,i,j,k+1) +yp*fin(4,i,j+1,k+1)) &
                +(ypi*fin(4,i+1,j,k+1)+yp*fin(4,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hz2*hxi(v)*( &
                czi*( &
                -(ypi*fin(7,i,j,k) +yp*fin(7,i,j+1,k)) &
                +(ypi*fin(7,i+1,j,k)+yp*fin(7,i+1,j+1,k))) &
                +cz*( &
                -(ypi*fin(7,i,j,k+1) +yp*fin(7,i,j+1,k+1)) &
                +(ypi*fin(7,i+1,j,k+1)+yp*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(3).eq.1) then
        !                               ! fxxxyz
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           j=jj(v)
           k=kk(v)
           !
           !   ...and in y direction
           !
           yp=yparam(v)
           ypi=1.0_fp-yp
           yp2=yp*yp
           ypi2=ypi*ypi

           cyd=3.0_fp*yp2-1.0_fp
           cydi=-3.0_fp*ypi2+1.0_fp
           !
           !   ...and in z direction
           !
           zp=zparam(v)
           zpi=1.0_fp-zp
           zp2=zp*zp
           zpi2=zpi*zpi

           czd=3.0_fp*zp2-1.0_fp
           czdi=-3.0_fp*zpi2+1.0_fp
           !
           sum=hyi(v)*hzi(v)*hxi(v)*( &
                -( &
                -(-fin(1,i,j,k)  +fin(1,i,j+1,k)) &
                +( -fin(1,i+1,j,k)+fin(1,i+1,j+1,k))) &
                +( &
                -(-fin(1,i,j,k+1)  +fin(1,i,j+1,k+1)) &
                +( -fin(1,i+1,j,k+1)+fin(1,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hy(v)*hzi(v)*hxi(v)*( &
                -( &
                -(cydi*fin(4,i,j,k)  +cyd*fin(4,i,j+1,k)) &
                +(cydi*fin(4,i+1,j,k)+cyd*fin(4,i+1,j+1,k))) &
                +( &
                -(cydi*fin(4,i,j,k+1)  +cyd*fin(4,i,j+1,k+1)) &
                +(cydi*fin(4,i+1,j,k+1)+cyd*fin(4,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hz(v)*hyi(v)*hxi(v)*( &
                czdi*( &
                -(-fin(5,i,j,k)  +fin(5,i,j+1,k)) &
                +( -fin(5,i+1,j,k)+fin(5,i+1,j+1,k))) &
                +czd*( &
                -(-fin(5,i,j,k+1)  +fin(5,i,j+1,k+1)) &
                +( -fin(5,i+1,j,k+1)+fin(5,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hy(v)*hz(v)*hxi(v)*( &
                czdi*( &
                -(cydi*fin(7,i,j,k)  +cyd*fin(7,i,j+1,k)) &
                +(cydi*fin(7,i+1,j,k)+cyd*fin(7,i+1,j+1,k))) &
                +czd*( &
                -(cydi*fin(7,i,j,k+1)  +cyd*fin(7,i,j+1,k+1)) &
                +(cydi*fin(7,i+1,j,k+1)+cyd*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(4).eq.1) then
        !                               ! fxxxzz
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           j=jj(v)
           k=kk(v)
           !
           !   ...and in y direction
           !
           yp=yparam(v)
           ypi=1.0_fp-yp
           yp2=yp*yp
           ypi2=ypi*ypi
           !
           cy=yp*(yp2-1.0_fp)
           cyi=ypi*(ypi2-1.0_fp)
           hy2=hy(v)*hy(v)
           !
           !   ...and in z direction
           !
           zp=zparam(v)
           zpi=1.0_fp-zp
           !
           sum=hxi(v)*( &
                zpi*( &
                -(ypi*fin(5,i,j,k)  +yp*fin(5,i,j+1,k)) &
                +(ypi*fin(5,i+1,j,k)+yp*fin(5,i+1,j+1,k))) &
                +zp*( &
                -(ypi*fin(5,i,j,k+1)  +yp*fin(5,i,j+1,k+1)) &
                +(ypi*fin(5,i+1,j,k+1)+yp*fin(5,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hy2*hxi(v)*( &
                zpi*( &
                -(cyi*fin(7,i,j,k)  +cy*fin(7,i,j+1,k)) &
                +(cyi*fin(7,i+1,j,k)+cy*fin(7,i+1,j+1,k))) &
                +zp*( &
                -(cyi*fin(7,i,j,k+1)  +cy*fin(7,i,j+1,k+1)) &
                +(cyi*fin(7,i+1,j,k+1)+cy*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(5).eq.1) then
        !                               ! fxxyyy
        iadr=iadr+1
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
           !   ...and in z direction
           !
           zp=zparam(v)
           zpi=1.0_fp-zp
           zp2=zp*zp
           zpi2=zpi*zpi
           !
           cz=zp*(zp2-1.0_fp)
           czi=zpi*(zpi2-1.0_fp)
           hz2=hz(v)*hz(v)

           !
           sum=hyi(v)*( &
                zpi*( &
                xpi*(-fin(4,i,j,k)  +fin(4,i,j+1,k))+ &
                xp*( -fin(4,i+1,j,k)+fin(4,i+1,j+1,k))) &
                +zp*( &
                xpi*(-fin(4,i,j,k+1)  +fin(4,i,j+1,k+1))+ &
                xp*( -fin(4,i+1,j,k+1)+fin(4,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hz2*hyi(v)*( &
                czi*( &
                xpi*(-fin(7,i,j,k)  +fin(7,i,j+1,k))+ &
                xp*( -fin(7,i+1,j,k)+fin(7,i+1,j+1,k))) &
                +cz*( &
                xpi*(-fin(7,i,j,k+1)  +fin(7,i,j+1,k+1))+ &
                xp*( -fin(7,i+1,j,k+1)+fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(6).eq.1) then
        !                               ! fxxzzz
        iadr=iadr+1
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
           !   ...and in y direction
           !
           yp=yparam(v)
           ypi=1.0_fp-yp
           yp2=yp*yp
           ypi2=ypi*ypi
           !
           cy=yp*(yp2-1.0_fp)
           cyi=ypi*(ypi2-1.0_fp)
           hy2=hy(v)*hy(v)
           !
           sum=hzi(v)*( &
                -( &
                xpi*(ypi*fin(5,i,j,k)  +yp*fin(5,i,j+1,k))+ &
                xp*(ypi*fin(5,i+1,j,k)+yp*fin(5,i+1,j+1,k))) &
                +( &
                xpi*(ypi*fin(5,i,j,k+1)  +yp*fin(5,i,j+1,k+1))+ &
                xp*(ypi*fin(5,i+1,j,k+1)+yp*fin(5,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hy2*hzi(v)*( &
                -( &
                xpi*(cyi*fin(7,i,j,k)  +cy*fin(7,i,j+1,k))+ &
                xp*(cyi*fin(7,i+1,j,k)+cy*fin(7,i+1,j+1,k))) &
                +( &
                xpi*(cyi*fin(7,i,j,k+1)  +cy*fin(7,i,j+1,k+1))+ &
                xp*(cyi*fin(7,i+1,j,k+1)+cy*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(7).eq.1) then
        !                               ! fxyyyz
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           j=jj(v)
           k=kk(v)
           !
           !   ...in x direction
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           xp2=xp*xp
           xpi2=xpi*xpi

           cxd=3.0_fp*xp2-1.0_fp
           cxdi=-3.0_fp*xpi2+1.0_fp
           !
           !   ...and in z direction
           !
           zp=zparam(v)
           zpi=1.0_fp-zp
           zp2=zp*zp
           zpi2=zpi*zpi

           czd=3.0_fp*zp2-1.0_fp
           czdi=-3.0_fp*zpi2+1.0_fp
           !
           sum=hxi(v)*hzi(v)*hyi(v)*( &
                -( &
                ( fin(2,i,j,k)  -fin(2,i,j+1,k)) &
                +(-fin(2,i+1,j,k)+fin(2,i+1,j+1,k))) &
                +( &
                ( fin(2,i,j,k+1)  -fin(2,i,j+1,k+1)) &
                +(-fin(2,i+1,j,k+1)+fin(2,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hx(v)*hzi(v)*hyi(v)*( &
                -( &
                cxdi*(-fin(4,i,j,k)  +fin(4,i,j+1,k))+ &
                cxd*(-fin(4,i+1,j,k) +fin(4,i+1,j+1,k))) &
                +( &
                cxdi*(-fin(4,i,j,k+1)  +fin(4,i,j+1,k+1))+ &
                cxd*(-fin(4,i+1,j,k+1) +fin(4,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hz(v)*hxi(v)*hyi(v)*( &
                czdi*( &
                ( fin(6,i,j,k)  -fin(6,i,j+1,k)) &
                +(-fin(6,i+1,j,k)+fin(6,i+1,j+1,k))) &
                +czd*( &
                ( fin(6,i,j,k+1)  -fin(6,i,j+1,k+1)) &
                +(-fin(6,i+1,j,k+1)+fin(6,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hx(v)*hz(v)*hyi(v)*( &
                czdi*( &
                cxdi*(-fin(7,i,j,k)  +fin(7,i,j+1,k))+ &
                cxd*(-fin(7,i+1,j,k) +fin(7,i+1,j+1,k))) &
                +czd*( &
                cxdi*(-fin(7,i,j,k+1)  +fin(7,i,j+1,k+1))+ &
                cxd*(-fin(7,i+1,j,k+1) +fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(8).eq.1) then
        !                               ! fxyzzz
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           j=jj(v)
           k=kk(v)
           !
           !   ...in x direction
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           xp2=xp*xp
           xpi2=xpi*xpi

           cxd=3.0_fp*xp2-1.0_fp
           cxdi=-3.0_fp*xpi2+1.0_fp
           !
           !   ...and in y direction
           !
           yp=yparam(v)
           ypi=1.0_fp-yp
           yp2=yp*yp
           ypi2=ypi*ypi

           cyd=3.0_fp*yp2-1.0_fp
           cydi=-3.0_fp*ypi2+1.0_fp
           !
           sum=hyi(v)*hxi(v)*hzi(v)*( &
                -( &
                ( +fin(3,i,j,k)  -fin(3,i,j+1,k)) &
                +( -fin(3,i+1,j,k)+fin(3,i+1,j+1,k))) &
                +( &
                ( +fin(3,i,j,k+1)  -fin(3,i,j+1,k+1)) &
                +( -fin(3,i+1,j,k+1)+fin(3,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hx(v)*hyi(v)*hzi(v)*( &
                -( &
                cxdi*(-fin(5,i,j,k)  +fin(5,i,j+1,k))+ &
                cxd*( -fin(5,i+1,j,k)+fin(5,i+1,j+1,k))) &
                +( &
                cxdi*(-fin(5,i,j,k+1)  +fin(5,i,j+1,k+1))+ &
                cxd*( -fin(5,i+1,j,k+1)+fin(5,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hy(v)*hxi(v)*hzi(v)*( &
                -( &
                -(cydi*fin(6,i,j,k) +cyd*fin(6,i,j+1,k)) &
                +(cydi*fin(6,i+1,j,k)+cyd*fin(6,i+1,j+1,k))) &
                +( &
                -(cydi*fin(6,i,j,k+1) +cyd*fin(6,i,j+1,k+1)) &
                +(cydi*fin(6,i+1,j,k+1)+cyd*fin(6,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hx(v)*hy(v)*hzi(v)*( &
                -( &
                cxdi*(cydi*fin(7,i,j,k) +cyd*fin(7,i,j+1,k))+ &
                cxd*(cydi*fin(7,i+1,j,k)+cyd*fin(7,i+1,j+1,k))) &
                +( &
                cxdi*(cydi*fin(7,i,j,k+1) +cyd*fin(7,i,j+1,k+1))+ &
                cxd*(cydi*fin(7,i+1,j,k+1)+cyd*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(9).eq.1) then
        !                               ! fyyyzz
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           j=jj(v)
           k=kk(v)
           !
           !   ...in x direction
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           xp2=xp*xp
           xpi2=xpi*xpi
           !
           cx=xp*(xp2-1.0_fp)
           cxi=xpi*(xpi2-1.0_fp)
           hx2=hx(v)*hx(v)
           !
           !   ...and in z direction
           !
           zp=zparam(v)
           zpi=1.0_fp-zp
           !
           sum=hyi(v)*( &
                zpi*( &
                xpi*(-fin(6,i,j,k)  +fin(6,i,j+1,k))+ &
                xp*( -fin(6,i+1,j,k)+fin(6,i+1,j+1,k))) &
                +zp*( &
                xpi*(-fin(6,i,j,k+1)  +fin(6,i,j+1,k+1))+ &
                xp*( -fin(6,i+1,j,k+1)+fin(6,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hx2*hyi(v)*( &
                zpi*( &
                cxi*(-fin(7,i,j,k)  +fin(7,i,j+1,k))+ &
                cx*( -fin(7,i+1,j,k)+fin(7,i+1,j+1,k))) &
                +zp*( &
                cxi*(-fin(7,i,j,k+1)  +fin(7,i,j+1,k+1))+ &
                cx*( -fin(7,i+1,j,k+1)+fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(10).eq.1) then
        !                               ! fyyzzz
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           j=jj(v)
           k=kk(v)
           !
           !   ...in x direction
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           xp2=xp*xp
           xpi2=xpi*xpi
           !
           cx=xp*(xp2-1.0_fp)
           cxi=xpi*(xpi2-1.0_fp)
           hx2=hx(v)*hx(v)
           !
           !   ...and in y direction
           !
           yp=yparam(v)
           ypi=1.0_fp-yp
           !
           sum=hzi(v)*( &
                -( &
                xpi*(ypi*fin(6,i,j,k)  +yp*fin(6,i,j+1,k))+ &
                xp*( ypi*fin(6,i+1,j,k)+yp*fin(6,i+1,j+1,k))) &
                +( &
                xpi*(ypi*fin(6,i,j,k+1) +yp*fin(6,i,j+1,k+1))+ &
                xp*(ypi*fin(6,i+1,j,k+1)+yp*fin(6,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hx2*hzi(v)*( &
                -( &
                cxi*(ypi*fin(7,i,j,k)  +yp*fin(7,i,j+1,k))+ &
                cx*( ypi*fin(7,i+1,j,k)+yp*fin(7,i+1,j+1,k))) &
                +( &
                cxi*(ypi*fin(7,i,j,k+1)  +yp*fin(7,i,j+1,k+1))+ &
                cx*( ypi*fin(7,i+1,j,k+1)+yp*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     !----------------------------------
     !  6th derivatives (2 in each coordinate)
     !
  else if(ict(1).eq.6) then
     !                               ! fxxyyzz
     iadr=iadr+1
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
        !   ...and in y direction
        !
        yp=yparam(v)
        ypi=1.0_fp-yp
        !
        !   ...and in z direction
        !
        zp=zparam(v)
        zpi=1.0_fp-zp
        !
        sum=( &
             zpi*( &
             xpi*(ypi*fin(7,i,j,k) +yp*fin(7,i,j+1,k))+ &
             xp*(ypi*fin(7,i+1,j,k)+yp*fin(7,i+1,j+1,k))) &
             +zp*( &
             xpi*(ypi*fin(7,i,j,k+1) +yp*fin(7,i,j+1,k+1))+ &
             xp*(ypi*fin(7,i+1,j,k+1)+yp*fin(7,i+1,j+1,k+1))))
        !
        fval(v,iadr)=sum
     end do
  end if
  !
  !----------------------------------
  !  6th derivatives (3 in a coordinate)
  !
  if(ict(1).eq.-6) then
     if(ict(2).eq.1) then
        !                               ! fxxxyyy
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           j=jj(v)
           k=kk(v)
           !
           !   ...and in z direction
           !
           zp=zparam(v)
           zpi=1.0_fp-zp
           zp2=zp*zp
           zpi2=zpi*zpi
           !
           cz=zp*(zp2-1.0_fp)
           czi=zpi*(zpi2-1.0_fp)
           hz2=hz(v)*hz(v)
           !
           sum=hyi(v)*hxi(v)*( &
                zpi*( &
                ( fin(4,i,j,k)  -fin(4,i,j+1,k)) &
                +(-fin(4,i+1,j,k)+fin(4,i+1,j+1,k))) &
                +zp*( &
                ( fin(4,i,j,k+1)  -fin(4,i,j+1,k+1)) &
                +(-fin(4,i+1,j,k+1)+fin(4,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hz2*hyi(v)*hxi(v)*( &
                czi*( &
                ( fin(7,i,j,k)  -fin(7,i,j+1,k)) &
                +(-fin(7,i+1,j,k)+fin(7,i+1,j+1,k))) &
                +cz*( &
                ( fin(7,i,j,k+1)  -fin(7,i,j+1,k+1)) &
                +(-fin(7,i+1,j,k+1)+fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(3).eq.1) then
        !                               ! fxxxyyz
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           j=jj(v)
           k=kk(v)
           !
           !   ...and in y direction
           !
           yp=yparam(v)
           ypi=1.0_fp-yp
           !
           !   ...and in z direction
           !
           zp=zparam(v)
           zpi=1.0_fp-zp
           zp2=zp*zp
           zpi2=zpi*zpi

           czd=3.0_fp*zp2-1.0_fp
           czdi=-3.0_fp*zpi2+1.0_fp
           !
           sum=hxi(v)*hzi(v)*( &
                -( &
                -(ypi*fin(4,i,j,k)  +yp*fin(4,i,j+1,k)) &
                +( ypi*fin(4,i+1,j,k)+yp*fin(4,i+1,j+1,k))) &
                +( &
                -(ypi*fin(4,i,j,k+1) +yp*fin(4,i,j+1,k+1)) &
                +(ypi*fin(4,i+1,j,k+1)+yp*fin(4,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hz(v)*hxi(v)*( &
                czdi*( &
                -(ypi*fin(7,i,j,k) +yp*fin(7,i,j+1,k)) &
                +(ypi*fin(7,i+1,j,k)+yp*fin(7,i+1,j+1,k))) &
                +czd*( &
                -(ypi*fin(7,i,j,k+1) +yp*fin(7,i,j+1,k+1)) &
                +(ypi*fin(7,i+1,j,k+1)+yp*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(4).eq.1) then
        !                               ! fxxxyzz
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           j=jj(v)
           k=kk(v)
           !
           !   ...and in y direction
           !
           yp=yparam(v)
           ypi=1.0_fp-yp
           yp2=yp*yp
           ypi2=ypi*ypi

           cyd=3.0_fp*yp2-1.0_fp
           cydi=-3.0_fp*ypi2+1.0_fp
           !
           !   ...and in z direction
           !
           zp=zparam(v)
           zpi=1.0_fp-zp
           !
           sum=hxi(v)*hyi(v)*( &
                zpi*( &
                ( fin(5,i,j,k)  -fin(5,i,j+1,k)) &
                +(-fin(5,i+1,j,k)+fin(5,i+1,j+1,k))) &
                +zp*( &
                ( fin(5,i,j,k+1)  -fin(5,i,j+1,k+1)) &
                +(-fin(5,i+1,j,k+1)+fin(5,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hy(v)*hxi(v)*( &
                zpi*( &
                -(cydi*fin(7,i,j,k)  +cyd*fin(7,i,j+1,k)) &
                +(cydi*fin(7,i+1,j,k)+cyd*fin(7,i+1,j+1,k))) &
                +zp*( &
                -(cydi*fin(7,i,j,k+1)  +cyd*fin(7,i,j+1,k+1)) &
                +(cydi*fin(7,i+1,j,k+1)+cyd*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(5).eq.1) then
        !                               ! fxxxzzz
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           j=jj(v)
           k=kk(v)
           !
           !   ...and in y direction
           !
           yp=yparam(v)
           ypi=1.0_fp-yp
           yp2=yp*yp
           ypi2=ypi*ypi
           !
           cy=yp*(yp2-1.0_fp)
           cyi=ypi*(ypi2-1.0_fp)
           hy2=hy(v)*hy(v)
           !
           sum=hxi(v)*hzi(v)*( &
                -( &
                -(ypi*fin(5,i,j,k)  +yp*fin(5,i,j+1,k)) &
                +(ypi*fin(5,i+1,j,k)+yp*fin(5,i+1,j+1,k))) &
                +( &
                -(ypi*fin(5,i,j,k+1)  +yp*fin(5,i,j+1,k+1)) &
                +(ypi*fin(5,i+1,j,k+1)+yp*fin(5,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hy2*hxi(v)*hzi(v)*( &
                -( &
                -(cyi*fin(7,i,j,k)  +cy*fin(7,i,j+1,k)) &
                +(cyi*fin(7,i+1,j,k)+cy*fin(7,i+1,j+1,k))) &
                +( &
                -(cyi*fin(7,i,j,k+1)  +cy*fin(7,i,j+1,k+1)) &
                +(cyi*fin(7,i+1,j,k+1)+cy*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(6).eq.1) then
        !                               ! fxxyyyz
        iadr=iadr+1
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
           !   ...and in z direction
           !
           zp=zparam(v)
           zpi=1.0_fp-zp
           zp2=zp*zp
           zpi2=zpi*zpi

           czd=3.0_fp*zp2-1.0_fp
           czdi=-3.0_fp*zpi2+1.0_fp
           !
           sum=hzi(v)*hyi(v)*( &
                -( &
                xpi*(-fin(4,i,j,k)  +fin(4,i,j+1,k))+ &
                xp*(-fin(4,i+1,j,k) +fin(4,i+1,j+1,k))) &
                +( &
                xpi*(-fin(4,i,j,k+1)  +fin(4,i,j+1,k+1))+ &
                xp*(-fin(4,i+1,j,k+1) +fin(4,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hz(v)*hyi(v)*( &
                czdi*( &
                xpi*(-fin(7,i,j,k)  +fin(7,i,j+1,k))+ &
                xp*(-fin(7,i+1,j,k) +fin(7,i+1,j+1,k))) &
                +czd*( &
                xpi*(-fin(7,i,j,k+1)  +fin(7,i,j+1,k+1))+ &
                xp*(-fin(7,i+1,j,k+1) +fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(7).eq.1) then
        !                               ! fxxyzzz
        iadr=iadr+1
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
           !   ...and in y direction
           !
           yp=yparam(v)
           ypi=1.0_fp-yp
           yp2=yp*yp
           ypi2=ypi*ypi

           cyd=3.0_fp*yp2-1.0_fp
           cydi=-3.0_fp*ypi2+1.0_fp
           !
           sum=hyi(v)*hzi(v)*( &
                -( &
                xpi*(-fin(5,i,j,k)  +fin(5,i,j+1,k))+ &
                xp*( -fin(5,i+1,j,k)+fin(5,i+1,j+1,k))) &
                +( &
                xpi*(-fin(5,i,j,k+1)  +fin(5,i,j+1,k+1))+ &
                xp*( -fin(5,i+1,j,k+1)+fin(5,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hy(v)*hzi(v)*( &
                -( &
                xpi*(cydi*fin(7,i,j,k)  +cyd*fin(7,i,j+1,k))+ &
                xp*( cydi*fin(7,i+1,j,k)+cyd*fin(7,i+1,j+1,k))) &
                +( &
                xpi*(cydi*fin(7,i,j,k+1)  +cyd*fin(7,i,j+1,k+1))+ &
                xp*( cydi*fin(7,i+1,j,k+1)+cyd*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(8).eq.1) then
        !                               ! fxyyyzz
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           j=jj(v)
           k=kk(v)
           !
           !   ...in x direction
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           xp2=xp*xp
           xpi2=xpi*xpi

           cxd=3.0_fp*xp2-1.0_fp
           cxdi=-3.0_fp*xpi2+1.0_fp
           !
           !   ...and in z direction
           !
           zp=zparam(v)
           zpi=1.0_fp-zp
           !
           sum=hxi(v)*hyi(v)*( &
                zpi*( &
                ( fin(6,i,j,k)  -fin(6,i,j+1,k)) &
                +(-fin(6,i+1,j,k)+fin(6,i+1,j+1,k))) &
                +zp*( &
                ( fin(6,i,j,k+1)  -fin(6,i,j+1,k+1)) &
                +(-fin(6,i+1,j,k+1)+fin(6,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hx(v)*hyi(v)*( &
                zpi*( &
                cxdi*(-fin(7,i,j,k)  +fin(7,i,j+1,k))+ &
                cxd*(-fin(7,i+1,j,k) +fin(7,i+1,j+1,k))) &
                +zp*( &
                cxdi*(-fin(7,i,j,k+1)  +fin(7,i,j+1,k+1))+ &
                cxd*(-fin(7,i+1,j,k+1) +fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(9).eq.1) then
        !                               ! fxyyzzz
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           j=jj(v)
           k=kk(v)
           !
           !   ...in x direction
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           xp2=xp*xp
           xpi2=xpi*xpi

           cxd=3.0_fp*xp2-1.0_fp
           cxdi=-3.0_fp*xpi2+1.0_fp
           !
           !   ...and in y direction
           !
           yp=yparam(v)
           ypi=1.0_fp-yp
           !
           sum=hxi(v)*hzi(v)*( &
                -( &
                -(ypi*fin(6,i,j,k)  +yp*fin(6,i,j+1,k)) &
                +(ypi*fin(6,i+1,j,k)+yp*fin(6,i+1,j+1,k))) &
                +( &
                -(ypi*fin(6,i,j,k+1)  +yp*fin(6,i,j+1,k+1)) &
                +(ypi*fin(6,i+1,j,k+1)+yp*fin(6,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hx(v)*hzi(v)*( &
                -( &
                cxdi*(ypi*fin(7,i,j,k)  +yp*fin(7,i,j+1,k))+ &
                cxd*(ypi*fin(7,i+1,j,k) +yp*fin(7,i+1,j+1,k))) &
                +( &
                cxdi*(ypi*fin(7,i,j,k+1)  +yp*fin(7,i,j+1,k+1))+ &
                cxd*( ypi*fin(7,i+1,j,k+1)+yp*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(10).eq.1) then
        !                               ! fyyyzzz
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           j=jj(v)
           k=kk(v)
           !
           !   ...in x direction
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           xp2=xp*xp
           xpi2=xpi*xpi
           !
           cx=xp*(xp2-1.0_fp)
           cxi=xpi*(xpi2-1.0_fp)
           hx2=hx(v)*hx(v)
           !
           sum=hyi(v)*hzi(v)*( &
                -( &
                xpi*(-fin(6,i,j,k)  +fin(6,i,j+1,k))+ &
                xp*( -fin(6,i+1,j,k)+fin(6,i+1,j+1,k))) &
                +( &
                xpi*(-fin(6,i,j,k+1)  +fin(6,i,j+1,k+1))+ &
                xp*( -fin(6,i+1,j,k+1)+fin(6,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hx2*hyi(v)*hzi(v)*( &
                -( &
                cxi*(-fin(7,i,j,k)  +fin(7,i,j+1,k))+ &
                cx*( -fin(7,i+1,j,k)+fin(7,i+1,j+1,k))) &
                +( &
                cxi*(-fin(7,i,j,k+1)  +fin(7,i,j+1,k+1))+ &
                cx*( -fin(7,i+1,j,k+1)+fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     !----------------------------------
     !  7th derivatives
     !
  else if(abs(ict(1)).eq.7) then
     if(ict(2).eq.1) then
        !                               ! fxxxyyyz
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           j=jj(v)
           k=kk(v)
           !
           !   ...and in z direction
           !
           zp=zparam(v)
           zpi=1.0_fp-zp
           zp2=zp*zp
           zpi2=zpi*zpi

           czd=3.0_fp*zp2-1.0_fp
           czdi=-3.0_fp*zpi2+1.0_fp
           !
           sum=hyi(v)*hxi(v)*hzi(v)*( &
                -( &
                ( fin(4,i,j,k)  -fin(4,i,j+1,k)) &
                +(-fin(4,i+1,j,k)+fin(4,i+1,j+1,k))) &
                +( &
                ( fin(4,i,j,k+1)  -fin(4,i,j+1,k+1)) &
                +(-fin(4,i+1,j,k+1)+fin(4,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hz(v)*hyi(v)*hxi(v)*( &
                czdi*( &
                ( fin(7,i,j,k)  -fin(7,i,j+1,k)) &
                +(-fin(7,i+1,j,k)+fin(7,i+1,j+1,k))) &
                +czd*( &
                ( fin(7,i,j,k+1)  -fin(7,i,j+1,k+1)) &
                +(-fin(7,i+1,j,k+1)+fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(3).eq.1) then
        !                               ! fxxxyyzz
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           j=jj(v)
           k=kk(v)
           !
           !   ...and in y direction
           !
           yp=yparam(v)
           ypi=1.0_fp-yp
           !
           !   ...and in z direction
           !
           zp=zparam(v)
           zpi=1.0_fp-zp
           !
           sum=hxi(v)*( &
                zpi*( &
                -(ypi*fin(7,i,j,k) +yp*fin(7,i,j+1,k)) &
                +(ypi*fin(7,i+1,j,k)+yp*fin(7,i+1,j+1,k))) &
                +zp*( &
                -(ypi*fin(7,i,j,k+1) +yp*fin(7,i,j+1,k+1)) &
                +(ypi*fin(7,i+1,j,k+1)+yp*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(4).eq.1) then
        !                               ! fxxxyzzz
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           j=jj(v)
           k=kk(v)
           !
           !   ...and in y direction
           !
           yp=yparam(v)
           ypi=1.0_fp-yp
           yp2=yp*yp
           ypi2=ypi*ypi

           cyd=3.0_fp*yp2-1.0_fp
           cydi=-3.0_fp*ypi2+1.0_fp
           !
           sum=hxi(v)*hyi(v)*hzi(v)*( &
                -( &
                ( fin(5,i,j,k)  -fin(5,i,j+1,k)) &
                +(-fin(5,i+1,j,k)+fin(5,i+1,j+1,k))) &
                +( &
                ( fin(5,i,j,k+1)  -fin(5,i,j+1,k+1)) &
                +(-fin(5,i+1,j,k+1)+fin(5,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hy(v)*hxi(v)*hzi(v)*( &
                -( &
                -(cydi*fin(7,i,j,k)  +cyd*fin(7,i,j+1,k)) &
                +(cydi*fin(7,i+1,j,k)+cyd*fin(7,i+1,j+1,k))) &
                +( &
                -(cydi*fin(7,i,j,k+1)  +cyd*fin(7,i,j+1,k+1)) &
                +(cydi*fin(7,i+1,j,k+1)+cyd*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(5).eq.1) then
        !                               ! fxxyyyzz
        iadr=iadr+1
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
           !   ...and in z direction
           !
           zp=zparam(v)
           zpi=1.0_fp-zp
           !
           sum=hyi(v)*( &
                zpi*( &
                xpi*(-fin(7,i,j,k)  +fin(7,i,j+1,k))+ &
                xp*( -fin(7,i+1,j,k)+fin(7,i+1,j+1,k))) &
                +zp*( &
                xpi*(-fin(7,i,j,k+1)  +fin(7,i,j+1,k+1))+ &
                xp*( -fin(7,i+1,j,k+1)+fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(6).eq.1) then
        !                               ! fxxyyzzz
        iadr=iadr+1
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
           !   ...and in y direction
           !
           yp=yparam(v)
           ypi=1.0_fp-yp
           !
           sum=hzi(v)*( &
                -( &
                xpi*(ypi*fin(7,i,j,k) +yp*fin(7,i,j+1,k))+ &
                xp*(ypi*fin(7,i+1,j,k)+yp*fin(7,i+1,j+1,k))) &
                +( &
                xpi*(ypi*fin(7,i,j,k+1) +yp*fin(7,i,j+1,k+1))+ &
                xp*(ypi*fin(7,i+1,j,k+1)+yp*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(7).eq.1) then
        !                               ! fxyyyzzz
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           j=jj(v)
           k=kk(v)
           !
           !   ...in x direction
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           xp2=xp*xp
           xpi2=xpi*xpi

           cxd=3.0_fp*xp2-1.0_fp
           cxdi=-3.0_fp*xpi2+1.0_fp
           !
           sum=hxi(v)*hyi(v)*hzi(v)*( &
                -( &
                ( fin(6,i,j,k)  -fin(6,i,j+1,k)) &
                +(-fin(6,i+1,j,k)+fin(6,i+1,j+1,k))) &
                +( &
                ( fin(6,i,j,k+1)  -fin(6,i,j+1,k+1)) &
                +(-fin(6,i+1,j,k+1)+fin(6,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hx(v)*hyi(v)*hzi(v)*( &
                -( &
                cxdi*(-fin(7,i,j,k)  +fin(7,i,j+1,k))+ &
                cxd*(-fin(7,i+1,j,k) +fin(7,i+1,j+1,k))) &
                +( &
                cxdi*(-fin(7,i,j,k+1)  +fin(7,i,j+1,k+1))+ &
                cxd*(-fin(7,i+1,j,k+1) +fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     !----------------------------------
     !  8th derivatives
     !
  else if(abs(ict(1)).eq.8) then
     if(ict(2).eq.1) then
        !                               ! fxxxyyyzz
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           j=jj(v)
           k=kk(v)
           !
           !   ...and in z direction
           !
           zp=zparam(v)
           zpi=1.0_fp-zp
           !
           sum=hyi(v)*hxi(v)*( &
                zpi*( &
                ( fin(7,i,j,k)  -fin(7,i,j+1,k)) &
                +(-fin(7,i+1,j,k)+fin(7,i+1,j+1,k))) &
                +zp*( &
                ( fin(7,i,j,k+1)  -fin(7,i,j+1,k+1)) &
                +(-fin(7,i+1,j,k+1)+fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(3).eq.1) then
        !                               ! fxxxyyzzz
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           j=jj(v)
           k=kk(v)
           !
           !   ...and in y direction
           !
           yp=yparam(v)
           ypi=1.0_fp-yp
           !
           sum=hxi(v)*hzi(v)*( &
                -( &
                -(ypi*fin(7,i,j,k) +yp*fin(7,i,j+1,k)) &
                +(ypi*fin(7,i+1,j,k)+yp*fin(7,i+1,j+1,k))) &
                +( &
                -(ypi*fin(7,i,j,k+1) +yp*fin(7,i,j+1,k+1)) &
                +(ypi*fin(7,i+1,j,k+1)+yp*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(4).eq.1) then
        !                               ! fxxyyyzzz
        iadr=iadr+1
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
           sum=hyi(v)*hzi(v)*( &
                -( &
                xpi*(-fin(7,i,j,k)  +fin(7,i,j+1,k))+ &
                xp*( -fin(7,i+1,j,k)+fin(7,i+1,j+1,k))) &
                +( &
                xpi*(-fin(7,i,j,k+1)  +fin(7,i,j+1,k+1))+ &
                xp*( -fin(7,i+1,j,k+1)+fin(7,i+1,j+1,k+1))))
           fval(v,iadr)=sum
           !
        end do
     end if
     !
     !----------------------------------
     !  9th derivative
     !
  else if(abs(ict(1)).eq.9) then
     !                               ! fxxxyyyzzz
     iadr=iadr+1
     do v=1,ivec
        i=ii(v)
        j=jj(v)
        k=kk(v)
        !
        sum=hyi(v)*hxi(v)*hzi(v)*( &
             -( &
             ( fin(7,i,j,k)  -fin(7,i,j+1,k)) &
             +(-fin(7,i+1,j,k)+fin(7,i+1,j+1,k))) &
             +( &
             ( fin(7,i,j,k+1)  -fin(7,i,j+1,k+1)) &
             +(-fin(7,i+1,j,k+1)+fin(7,i+1,j+1,k+1))))
        !
        fval(v,iadr)=sum
     end do
  end if
  !
  return
end subroutine fvtricub
!---------------------------------------------------------------------
!  evaluate C1 cubic Hermite function interpolation -- 3d fcn
!   --vectorized-- dmc 10 Feb 1999
!    --optimized for VARIATION along x axis ONLY--
!
subroutine fvtricubx(ict,ivec,ivecd, &
     fval,ii,jj,kk,xparam,yparam,zparam, &
     hx,hxi,hy,hyi,hz,hzi, &
     fin,inf2,inf3,nz)
  use psp_precision_mod, only: fp
  !
  !  use mktricub to set up spline coefficients...
  !
  !============
  implicit none
  integer inf3,nz,inf2,iadr,j,k,i
  !============
  real(fp) :: z36th,z216th,yp,ypi,yp2,ypi2,cy,cyi,hy2,zp,zpi,zp2
  real(fp) :: zpi2,cz,czi,hz2,xp,xpi,xp2,xpi2,cx,cxi,hx2,cxd,cxdi
  real(fp) :: cyd,cydi,czd,czdi
  !============
  integer ict(10)                   ! requested output control
  integer ivec                      ! vector length
  integer ivecd                     ! vector dimension (1st dim of fval)
  !
  integer ii(ivec),jj,kk            ! target cells (i,j,k)
  real(fp) :: xparam(ivec),yparam,zparam
  ! normalized displacements from (i,j,k) corners
  !
  real(fp) :: hx(ivec),hy,hz               ! grid spacing, and
  real(fp) :: hxi(ivec),hyi,hzi            ! inverse grid spacing
  ! 1/(x(i+1)-x(i)) & 1/(y(j+1)-y(j)) & 1/(z(k+1)-z(i))
  !
  real(fp) :: fin(0:7,inf2,inf3,nz)        ! interpolant data (cf "evtricub")
  !
  real(fp) :: fval(ivecd,*)               ! output returned
  !
  !  for detailed description of fin, ict and fval see subroutine evtricub
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
  integer v
  !
  real(fp) :: sum
  real(fp), parameter :: sixth = 0.166666666666666667_fp
  !
  !---------------
  !
  z36th=sixth*sixth
  z216th=sixth*sixth*sixth
  !
  iadr=0
  if(abs(ict(1)).le.2) then
     !
     !  0, 1st, 2nd derivatives...
     !
     !  get desired values:
     !
     if(ict(1).eq.1) then
        !
        !  function value...
        !
        j=jj
        k=kk
        !
        !   ...and in y direction
        !
        yp=yparam
        ypi=1.0_fp-yp
        yp2=yp*yp
        ypi2=ypi*ypi
        !
        cy=yp*(yp2-1.0_fp)
        cyi=ypi*(ypi2-1.0_fp)
        hy2=hy*hy
        !
        !   ...and in z direction
        !
        zp=zparam
        zpi=1.0_fp-zp
        zp2=zp*zp
        zpi2=zpi*zpi
        !
        cz=zp*(zp2-1.0_fp)
        czi=zpi*(zpi2-1.0_fp)
        hz2=hz*hz

        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           !
           !   ...in x direction
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           xp2=xp*xp
           xpi2=xpi*xpi
           !
           cx=xp*(xp2-1.0_fp)
           cxi=xpi*(xpi2-1.0_fp)
           hx2=hx(v)*hx(v)
           !
           sum=( &
                zpi*( &
                xpi*(ypi*fin(0,i,j,k)  +yp*fin(0,i,j+1,k))+ &
                xp*(ypi*fin(0,i+1,j,k)+yp*fin(0,i+1,j+1,k))) &
                +zp*( &
                xpi*(ypi*fin(0,i,j,k+1)  +yp*fin(0,i,j+1,k+1))+ &
                xp*(ypi*fin(0,i+1,j,k+1)+yp*fin(0,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hx2*( &
                zpi*( &
                cxi*(ypi*fin(1,i,j,k)  +yp*fin(1,i,j+1,k))+ &
                cx*(ypi*fin(1,i+1,j,k)+yp*fin(1,i+1,j+1,k))) &
                +zp*( &
                cxi*(ypi*fin(1,i,j,k+1)  +yp*fin(1,i,j+1,k+1))+ &
                cx*(ypi*fin(1,i+1,j,k+1)+yp*fin(1,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hy2*( &
                zpi*( &
                xpi*(cyi*fin(2,i,j,k)  +cy*fin(2,i,j+1,k))+ &
                xp*(cyi*fin(2,i+1,j,k)+cy*fin(2,i+1,j+1,k))) &
                +zp*( &
                xpi*(cyi*fin(2,i,j,k+1)  +cy*fin(2,i,j+1,k+1))+ &
                xp*(cyi*fin(2,i+1,j,k+1)+cy*fin(2,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hz2*( &
                czi*( &
                xpi*(ypi*fin(3,i,j,k)  +yp*fin(3,i,j+1,k))+ &
                xp*(ypi*fin(3,i+1,j,k)+yp*fin(3,i+1,j+1,k))) &
                +cz*( &
                xpi*(ypi*fin(3,i,j,k+1)  +yp*fin(3,i,j+1,k+1))+ &
                xp*(ypi*fin(3,i+1,j,k+1)+yp*fin(3,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hx2*hy2*( &
                zpi*( &
                cxi*(cyi*fin(4,i,j,k)  +cy*fin(4,i,j+1,k))+ &
                cx*(cyi*fin(4,i+1,j,k)+cy*fin(4,i+1,j+1,k))) &
                +zp*( &
                cxi*(cyi*fin(4,i,j,k+1)  +cy*fin(4,i,j+1,k+1))+ &
                cx*(cyi*fin(4,i+1,j,k+1)+cy*fin(4,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hx2*hz2*( &
                czi*( &
                cxi*(ypi*fin(5,i,j,k)  +yp*fin(5,i,j+1,k))+ &
                cx*(ypi*fin(5,i+1,j,k)+yp*fin(5,i+1,j+1,k))) &
                +cz*( &
                cxi*(ypi*fin(5,i,j,k+1)  +yp*fin(5,i,j+1,k+1))+ &
                cx*(ypi*fin(5,i+1,j,k+1)+yp*fin(5,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hy2*hz2*( &
                czi*( &
                xpi*(cyi*fin(6,i,j,k)  +cy*fin(6,i,j+1,k))+ &
                xp*(cyi*fin(6,i+1,j,k)+cy*fin(6,i+1,j+1,k))) &
                +cz*( &
                xpi*(cyi*fin(6,i,j,k+1)  +cy*fin(6,i,j+1,k+1))+ &
                xp*(cyi*fin(6,i+1,j,k+1)+cy*fin(6,i+1,j+1,k+1))))
           !
           sum=sum+z216th*hx2*hy2*hz2*( &
                czi*( &
                cxi*(cyi*fin(7,i,j,k)  +cy*fin(7,i,j+1,k))+ &
                cx*(cyi*fin(7,i+1,j,k)+cy*fin(7,i+1,j+1,k))) &
                +cz*( &
                cxi*(cyi*fin(7,i,j,k+1)  +cy*fin(7,i,j+1,k+1))+ &
                cx*(cyi*fin(7,i+1,j,k+1)+cy*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(2).eq.1) then
        !
        !  df/dx:
        !
        j=jj
        k=kk
        !
        !   ...and in y direction
        !
        yp=yparam
        ypi=1.0_fp-yp
        yp2=yp*yp
        ypi2=ypi*ypi
        !
        cy=yp*(yp2-1.0_fp)
        cyi=ypi*(ypi2-1.0_fp)
        hy2=hy*hy
        !
        !   ...and in z direction
        !
        zp=zparam
        zpi=1.0_fp-zp
        zp2=zp*zp
        zpi2=zpi*zpi
        !
        cz=zp*(zp2-1.0_fp)
        czi=zpi*(zpi2-1.0_fp)
        hz2=hz*hz

        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           !
           !   ...in x direction
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           xp2=xp*xp
           xpi2=xpi*xpi

           cxd=3.0_fp*xp2-1.0_fp
           cxdi=-3.0_fp*xpi2+1.0_fp
           !
           sum=hxi(v)*( &
                zpi*( &
                -(ypi*fin(0,i,j,k)  +yp*fin(0,i,j+1,k)) &
                +(ypi*fin(0,i+1,j,k)+yp*fin(0,i+1,j+1,k))) &
                +zp*( &
                -(ypi*fin(0,i,j,k+1)  +yp*fin(0,i,j+1,k+1)) &
                +(ypi*fin(0,i+1,j,k+1)+yp*fin(0,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hx(v)*( &
                zpi*( &
                cxdi*(ypi*fin(1,i,j,k)  +yp*fin(1,i,j+1,k))+ &
                cxd*(ypi*fin(1,i+1,j,k)+yp*fin(1,i+1,j+1,k))) &
                +zp*( &
                cxdi*(ypi*fin(1,i,j,k+1)  +yp*fin(1,i,j+1,k+1))+ &
                cxd*(ypi*fin(1,i+1,j,k+1)+yp*fin(1,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hxi(v)*hy2*( &
                zpi*( &
                -(cyi*fin(2,i,j,k)  +cy*fin(2,i,j+1,k)) &
                +(cyi*fin(2,i+1,j,k)+cy*fin(2,i+1,j+1,k))) &
                +zp*( &
                -(cyi*fin(2,i,j,k+1)  +cy*fin(2,i,j+1,k+1)) &
                +(cyi*fin(2,i+1,j,k+1)+cy*fin(2,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hxi(v)*hz2*( &
                czi*( &
                -(ypi*fin(3,i,j,k)  +yp*fin(3,i,j+1,k)) &
                +(ypi*fin(3,i+1,j,k)+yp*fin(3,i+1,j+1,k))) &
                +cz*( &
                -(ypi*fin(3,i,j,k+1)  +yp*fin(3,i,j+1,k+1)) &
                +(ypi*fin(3,i+1,j,k+1)+yp*fin(3,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hx(v)*hy2*( &
                zpi*( &
                cxdi*(cyi*fin(4,i,j,k)  +cy*fin(4,i,j+1,k))+ &
                cxd*(cyi*fin(4,i+1,j,k)+cy*fin(4,i+1,j+1,k))) &
                +zp*( &
                cxdi*(cyi*fin(4,i,j,k+1)  +cy*fin(4,i,j+1,k+1))+ &
                cxd*(cyi*fin(4,i+1,j,k+1)+cy*fin(4,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hx(v)*hz2*( &
                czi*( &
                cxdi*(ypi*fin(5,i,j,k)  +yp*fin(5,i,j+1,k))+ &
                cxd*(ypi*fin(5,i+1,j,k)+yp*fin(5,i+1,j+1,k))) &
                +cz*( &
                cxdi*(ypi*fin(5,i,j,k+1)  +yp*fin(5,i,j+1,k+1))+ &
                cxd*(ypi*fin(5,i+1,j,k+1)+yp*fin(5,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hxi(v)*hy2*hz2*( &
                czi*( &
                -(cyi*fin(6,i,j,k)  +cy*fin(6,i,j+1,k)) &
                +(cyi*fin(6,i+1,j,k)+cy*fin(6,i+1,j+1,k))) &
                +cz*( &
                -(cyi*fin(6,i,j,k+1)  +cy*fin(6,i,j+1,k+1)) &
                +(cyi*fin(6,i+1,j,k+1)+cy*fin(6,i+1,j+1,k+1))))
           !
           sum=sum+z216th*hx(v)*hy2*hz2*( &
                czi*( &
                cxdi*(cyi*fin(7,i,j,k)  +cy*fin(7,i,j+1,k))+ &
                cxd*(cyi*fin(7,i+1,j,k)+cy*fin(7,i+1,j+1,k))) &
                +cz*( &
                cxdi*(cyi*fin(7,i,j,k+1)  +cy*fin(7,i,j+1,k+1))+ &
                cxd*(cyi*fin(7,i+1,j,k+1)+cy*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(3).eq.1) then
        !
        !  df/dy:
        !
        j=jj
        k=kk
        !
        !   ...and in y direction
        !
        yp=yparam
        ypi=1.0_fp-yp
        yp2=yp*yp
        ypi2=ypi*ypi

        cyd=3.0_fp*yp2-1.0_fp
        cydi=-3.0_fp*ypi2+1.0_fp
        !
        !   ...and in z direction
        !
        zp=zparam
        zpi=1.0_fp-zp
        zp2=zp*zp
        zpi2=zpi*zpi
        !
        cz=zp*(zp2-1.0_fp)
        czi=zpi*(zpi2-1.0_fp)
        hz2=hz*hz
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           !
           !   ...in x direction
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           xp2=xp*xp
           xpi2=xpi*xpi
           !
           cx=xp*(xp2-1.0_fp)
           cxi=xpi*(xpi2-1.0_fp)
           hx2=hx(v)*hx(v)
           !
           sum=hyi*( &
                zpi*( &
                xpi*(-fin(0,i,j,k)  +fin(0,i,j+1,k))+ &
                xp*(-fin(0,i+1,j,k)+fin(0,i+1,j+1,k))) &
                +zp*( &
                xpi*(-fin(0,i,j,k+1)  +fin(0,i,j+1,k+1))+ &
                xp*(-fin(0,i+1,j,k+1)+fin(0,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hyi*hx2*( &
                zpi*( &
                cxi*(-fin(1,i,j,k)  +fin(1,i,j+1,k))+ &
                cx*(-fin(1,i+1,j,k)+fin(1,i+1,j+1,k))) &
                +zp*( &
                cxi*(-fin(1,i,j,k+1)  +fin(1,i,j+1,k+1))+ &
                cx*(-fin(1,i+1,j,k+1)+fin(1,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hy*( &
                zpi*( &
                xpi*(cydi*fin(2,i,j,k)  +cyd*fin(2,i,j+1,k))+ &
                xp*(cydi*fin(2,i+1,j,k)+cyd*fin(2,i+1,j+1,k))) &
                +zp*( &
                xpi*(cydi*fin(2,i,j,k+1)  +cyd*fin(2,i,j+1,k+1))+ &
                xp*(cydi*fin(2,i+1,j,k+1)+cyd*fin(2,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hyi*hz2*( &
                czi*( &
                xpi*(-fin(3,i,j,k)  +fin(3,i,j+1,k))+ &
                xp*(-fin(3,i+1,j,k)+fin(3,i+1,j+1,k))) &
                +cz*( &
                xpi*(-fin(3,i,j,k+1)  +fin(3,i,j+1,k+1))+ &
                xp*(-fin(3,i+1,j,k+1)+fin(3,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hx2*hy*( &
                zpi*( &
                cxi*(cydi*fin(4,i,j,k)  +cyd*fin(4,i,j+1,k))+ &
                cx*(cydi*fin(4,i+1,j,k)+cyd*fin(4,i+1,j+1,k))) &
                +zp*( &
                cxi*(cydi*fin(4,i,j,k+1)  +cyd*fin(4,i,j+1,k+1))+ &
                cx*(cydi*fin(4,i+1,j,k+1)+cyd*fin(4,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hyi*hx2*hz2*( &
                czi*( &
                cxi*(-fin(5,i,j,k)  +fin(5,i,j+1,k))+ &
                cx*(-fin(5,i+1,j,k)+fin(5,i+1,j+1,k))) &
                +cz*( &
                cxi*(-fin(5,i,j,k+1)  +fin(5,i,j+1,k+1))+ &
                cx*(-fin(5,i+1,j,k+1)+fin(5,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hy*hz2*( &
                czi*( &
                xpi*(cydi*fin(6,i,j,k)  +cyd*fin(6,i,j+1,k))+ &
                xp*(cydi*fin(6,i+1,j,k)+cyd*fin(6,i+1,j+1,k))) &
                +cz*( &
                xpi*(cydi*fin(6,i,j,k+1)  +cyd*fin(6,i,j+1,k+1))+ &
                xp*(cydi*fin(6,i+1,j,k+1)+cyd*fin(6,i+1,j+1,k+1))))
           !
           sum=sum+z216th*hx2*hy*hz2*( &
                czi*( &
                cxi*(cydi*fin(7,i,j,k)  +cyd*fin(7,i,j+1,k))+ &
                cx*(cydi*fin(7,i+1,j,k)+cyd*fin(7,i+1,j+1,k))) &
                +cz*( &
                cxi*(cydi*fin(7,i,j,k+1)  +cyd*fin(7,i,j+1,k+1))+ &
                cx*(cydi*fin(7,i+1,j,k+1)+cyd*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(4).eq.1) then
        !
        !  df/dz:
        !
        j=jj
        k=kk
        !
        !   ...and in y direction
        !
        yp=yparam
        ypi=1.0_fp-yp
        yp2=yp*yp
        ypi2=ypi*ypi
        !
        cy=yp*(yp2-1.0_fp)
        cyi=ypi*(ypi2-1.0_fp)
        hy2=hy*hy
        !
        !   ...and in z direction
        !
        zp=zparam
        zpi=1.0_fp-zp
        zp2=zp*zp
        zpi2=zpi*zpi

        czd=3.0_fp*zp2-1.0_fp
        czdi=-3.0_fp*zpi2+1.0_fp
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           !
           !   ...in x direction
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           xp2=xp*xp
           xpi2=xpi*xpi
           !
           cx=xp*(xp2-1.0_fp)
           cxi=xpi*(xpi2-1.0_fp)
           hx2=hx(v)*hx(v)
           !
           sum=hzi*( &
                -( &
                xpi*(ypi*fin(0,i,j,k)  +yp*fin(0,i,j+1,k))+ &
                xp*(ypi*fin(0,i+1,j,k)+yp*fin(0,i+1,j+1,k))) &
                +( &
                xpi*(ypi*fin(0,i,j,k+1)  +yp*fin(0,i,j+1,k+1))+ &
                xp*(ypi*fin(0,i+1,j,k+1)+yp*fin(0,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hx2*hzi*( &
                -( &
                cxi*(ypi*fin(1,i,j,k)  +yp*fin(1,i,j+1,k))+ &
                cx*(ypi*fin(1,i+1,j,k)+yp*fin(1,i+1,j+1,k))) &
                +( &
                cxi*(ypi*fin(1,i,j,k+1)  +yp*fin(1,i,j+1,k+1))+ &
                cx*(ypi*fin(1,i+1,j,k+1)+yp*fin(1,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hy2*hzi*( &
                -( &
                xpi*(cyi*fin(2,i,j,k)  +cy*fin(2,i,j+1,k))+ &
                xp*(cyi*fin(2,i+1,j,k)+cy*fin(2,i+1,j+1,k))) &
                +( &
                xpi*(cyi*fin(2,i,j,k+1)  +cy*fin(2,i,j+1,k+1))+ &
                xp*(cyi*fin(2,i+1,j,k+1)+cy*fin(2,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hz*( &
                czdi*( &
                xpi*(ypi*fin(3,i,j,k)  +yp*fin(3,i,j+1,k))+ &
                xp*(ypi*fin(3,i+1,j,k)+yp*fin(3,i+1,j+1,k))) &
                +czd*( &
                xpi*(ypi*fin(3,i,j,k+1)  +yp*fin(3,i,j+1,k+1))+ &
                xp*(ypi*fin(3,i+1,j,k+1)+yp*fin(3,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hx2*hy2*hzi*( &
                -( &
                cxi*(cyi*fin(4,i,j,k)  +cy*fin(4,i,j+1,k))+ &
                cx*(cyi*fin(4,i+1,j,k)+cy*fin(4,i+1,j+1,k))) &
                +( &
                cxi*(cyi*fin(4,i,j,k+1)  +cy*fin(4,i,j+1,k+1))+ &
                cx*(cyi*fin(4,i+1,j,k+1)+cy*fin(4,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hx2*hz*( &
                czdi*( &
                cxi*(ypi*fin(5,i,j,k)  +yp*fin(5,i,j+1,k))+ &
                cx*(ypi*fin(5,i+1,j,k)+yp*fin(5,i+1,j+1,k))) &
                +czd*( &
                cxi*(ypi*fin(5,i,j,k+1)  +yp*fin(5,i,j+1,k+1))+ &
                cx*(ypi*fin(5,i+1,j,k+1)+yp*fin(5,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hy2*hz*( &
                czdi*( &
                xpi*(cyi*fin(6,i,j,k)  +cy*fin(6,i,j+1,k))+ &
                xp*(cyi*fin(6,i+1,j,k)+cy*fin(6,i+1,j+1,k))) &
                +czd*( &
                xpi*(cyi*fin(6,i,j,k+1)  +cy*fin(6,i,j+1,k+1))+ &
                xp*(cyi*fin(6,i+1,j,k+1)+cy*fin(6,i+1,j+1,k+1))))
           !
           sum=sum+z216th*hx2*hy2*hz*( &
                czdi*( &
                cxi*(cyi*fin(7,i,j,k)  +cy*fin(7,i,j+1,k))+ &
                cx*(cyi*fin(7,i+1,j,k)+cy*fin(7,i+1,j+1,k))) &
                +czd*( &
                cxi*(cyi*fin(7,i,j,k+1)  +cy*fin(7,i,j+1,k+1))+ &
                cx*(cyi*fin(7,i+1,j,k+1)+cy*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(5).eq.1) then
        !
        !  d2f/dx2:
        !
        j=jj
        k=kk
        !
        !   ...and in y direction
        !
        yp=yparam
        ypi=1.0_fp-yp
        yp2=yp*yp
        ypi2=ypi*ypi
        !
        cy=yp*(yp2-1.0_fp)
        cyi=ypi*(ypi2-1.0_fp)
        hy2=hy*hy
        !
        !   ...and in z direction
        !
        zp=zparam
        zpi=1.0_fp-zp
        zp2=zp*zp
        zpi2=zpi*zpi
        !
        cz=zp*(zp2-1.0_fp)
        czi=zpi*(zpi2-1.0_fp)
        hz2=hz*hz
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           !
           !   ...in x direction
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           !
           sum=( &
                zpi*( &
                xpi*(ypi*fin(1,i,j,k)  +yp*fin(1,i,j+1,k))+ &
                xp*(ypi*fin(1,i+1,j,k)+yp*fin(1,i+1,j+1,k))) &
                +zp*( &
                xpi*(ypi*fin(1,i,j,k+1)  +yp*fin(1,i,j+1,k+1))+ &
                xp*(ypi*fin(1,i+1,j,k+1)+yp*fin(1,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hy2*( &
                zpi*( &
                xpi*(cyi*fin(4,i,j,k)  +cy*fin(4,i,j+1,k))+ &
                xp*(cyi*fin(4,i+1,j,k)+cy*fin(4,i+1,j+1,k))) &
                +zp*( &
                xpi*(cyi*fin(4,i,j,k+1)  +cy*fin(4,i,j+1,k+1))+ &
                xp*(cyi*fin(4,i+1,j,k+1)+cy*fin(4,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hz2*( &
                czi*( &
                xpi*(ypi*fin(5,i,j,k)  +yp*fin(5,i,j+1,k))+ &
                xp*(ypi*fin(5,i+1,j,k)+yp*fin(5,i+1,j+1,k))) &
                +cz*( &
                xpi*(ypi*fin(5,i,j,k+1)  +yp*fin(5,i,j+1,k+1))+ &
                xp*(ypi*fin(5,i+1,j,k+1)+yp*fin(5,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hy2*hz2*( &
                czi*( &
                xpi*(cyi*fin(7,i,j,k)  +cy*fin(7,i,j+1,k))+ &
                xp*(cyi*fin(7,i+1,j,k)+cy*fin(7,i+1,j+1,k))) &
                +cz*( &
                xpi*(cyi*fin(7,i,j,k+1)  +cy*fin(7,i,j+1,k+1))+ &
                xp*(cyi*fin(7,i+1,j,k+1)+cy*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(6).eq.1) then
        !
        !  d2f/dy2:
        !
        j=jj
        k=kk
        !
        !   ...and in y direction
        !
        yp=yparam
        ypi=1.0_fp-yp
        !
        !   ...and in z direction
        !
        zp=zparam
        zpi=1.0_fp-zp
        zp2=zp*zp
        zpi2=zpi*zpi
        !
        cz=zp*(zp2-1.0_fp)
        czi=zpi*(zpi2-1.0_fp)
        hz2=hz*hz
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           !
           !   ...in x direction
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           xp2=xp*xp
           xpi2=xpi*xpi
           !
           cx=xp*(xp2-1.0_fp)
           cxi=xpi*(xpi2-1.0_fp)
           hx2=hx(v)*hx(v)
           !
           sum=( &
                zpi*( &
                xpi*(ypi*fin(2,i,j,k)  +yp*fin(2,i,j+1,k))+ &
                xp*(ypi*fin(2,i+1,j,k)+yp*fin(2,i+1,j+1,k))) &
                +zp*( &
                xpi*(ypi*fin(2,i,j,k+1)  +yp*fin(2,i,j+1,k+1))+ &
                xp*(ypi*fin(2,i+1,j,k+1)+yp*fin(2,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hx2*( &
                zpi*( &
                cxi*(ypi*fin(4,i,j,k)  +yp*fin(4,i,j+1,k))+ &
                cx*(ypi*fin(4,i+1,j,k)+yp*fin(4,i+1,j+1,k))) &
                +zp*( &
                cxi*(ypi*fin(4,i,j,k+1)  +yp*fin(4,i,j+1,k+1))+ &
                cx*(ypi*fin(4,i+1,j,k+1)+yp*fin(4,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hz2*( &
                czi*( &
                xpi*(ypi*fin(6,i,j,k)  +yp*fin(6,i,j+1,k))+ &
                xp*(ypi*fin(6,i+1,j,k)+yp*fin(6,i+1,j+1,k))) &
                +cz*( &
                xpi*(ypi*fin(6,i,j,k+1)  +yp*fin(6,i,j+1,k+1))+ &
                xp*(ypi*fin(6,i+1,j,k+1)+yp*fin(6,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hx2*hz2*( &
                czi*( &
                cxi*(ypi*fin(7,i,j,k)  +yp*fin(7,i,j+1,k))+ &
                cx*(ypi*fin(7,i+1,j,k)+yp*fin(7,i+1,j+1,k))) &
                +cz*( &
                cxi*(ypi*fin(7,i,j,k+1)  +yp*fin(7,i,j+1,k+1))+ &
                cx*(ypi*fin(7,i+1,j,k+1)+yp*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(7).eq.1) then
        !
        !  d2f/dz2:
        !
        j=jj
        k=kk
        !
        !   ...and in y direction
        !
        yp=yparam
        ypi=1.0_fp-yp
        yp2=yp*yp
        ypi2=ypi*ypi
        !
        cy=yp*(yp2-1.0_fp)
        cyi=ypi*(ypi2-1.0_fp)
        hy2=hy*hy
        !
        !   ...and in z direction
        !
        zp=zparam
        zpi=1.0_fp-zp

        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           !
           !   ...in x direction
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           xp2=xp*xp
           xpi2=xpi*xpi
           !
           cx=xp*(xp2-1.0_fp)
           cxi=xpi*(xpi2-1.0_fp)
           hx2=hx(v)*hx(v)
           !
           sum=( &
                zpi*( &
                xpi*(ypi*fin(3,i,j,k)  +yp*fin(3,i,j+1,k))+ &
                xp*(ypi*fin(3,i+1,j,k)+yp*fin(3,i+1,j+1,k))) &
                +zp*( &
                xpi*(ypi*fin(3,i,j,k+1)  +yp*fin(3,i,j+1,k+1))+ &
                xp*(ypi*fin(3,i+1,j,k+1)+yp*fin(3,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hx2*( &
                zpi*( &
                cxi*(ypi*fin(5,i,j,k)  +yp*fin(5,i,j+1,k))+ &
                cx*(ypi*fin(5,i+1,j,k)+yp*fin(5,i+1,j+1,k))) &
                +zp*( &
                cxi*(ypi*fin(5,i,j,k+1)  +yp*fin(5,i,j+1,k+1))+ &
                cx*(ypi*fin(5,i+1,j,k+1)+yp*fin(5,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hy2*( &
                zpi*( &
                xpi*(cyi*fin(6,i,j,k)  +cy*fin(6,i,j+1,k))+ &
                xp*(cyi*fin(6,i+1,j,k)+cy*fin(6,i+1,j+1,k))) &
                +zp*( &
                xpi*(cyi*fin(6,i,j,k+1)  +cy*fin(6,i,j+1,k+1))+ &
                xp*(cyi*fin(6,i+1,j,k+1)+cy*fin(6,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hx2*hy2*( &
                zpi*( &
                cxi*(cyi*fin(7,i,j,k)  +cy*fin(7,i,j+1,k))+ &
                cx*(cyi*fin(7,i+1,j,k)+cy*fin(7,i+1,j+1,k))) &
                +zp*( &
                cxi*(cyi*fin(7,i,j,k+1)  +cy*fin(7,i,j+1,k+1))+ &
                cx*(cyi*fin(7,i+1,j,k+1)+cy*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(8).eq.1) then
        !
        !  d2f/dxdy:
        !
        j=jj
        k=kk
        !
        !   ...and in y direction
        !
        yp=yparam
        ypi=1.0_fp-yp
        yp2=yp*yp
        ypi2=ypi*ypi

        cyd=3.0_fp*yp2-1.0_fp
        cydi=-3.0_fp*ypi2+1.0_fp
        !
        !   ...and in z direction
        !
        zp=zparam
        zpi=1.0_fp-zp
        zp2=zp*zp
        zpi2=zpi*zpi
        !
        cz=zp*(zp2-1.0_fp)
        czi=zpi*(zpi2-1.0_fp)
        hz2=hz*hz
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           !
           !   ...in x direction
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           xp2=xp*xp
           xpi2=xpi*xpi

           cxd=3.0_fp*xp2-1.0_fp
           cxdi=-3.0_fp*xpi2+1.0_fp
           !
           sum=hxi(v)*hyi*( &
                zpi*( &
                (fin(0,i,j,k)  -fin(0,i,j+1,k))- &
                (fin(0,i+1,j,k)-fin(0,i+1,j+1,k))) &
                +zp*( &
                (fin(0,i,j,k+1)  -fin(0,i,j+1,k+1))- &
                (fin(0,i+1,j,k+1)-fin(0,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hyi*hx(v)*( &
                zpi*( &
                cxdi*(-fin(1,i,j,k)  +fin(1,i,j+1,k))+ &
                cxd*(-fin(1,i+1,j,k)+fin(1,i+1,j+1,k))) &
                +zp*( &
                cxdi*(-fin(1,i,j,k+1)  +fin(1,i,j+1,k+1))+ &
                cxd*(-fin(1,i+1,j,k+1)+fin(1,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hxi(v)*hy*( &
                zpi*( &
                -(cydi*fin(2,i,j,k)  +cyd*fin(2,i,j+1,k)) &
                +(cydi*fin(2,i+1,j,k)+cyd*fin(2,i+1,j+1,k))) &
                +zp*( &
                -(cydi*fin(2,i,j,k+1)  +cyd*fin(2,i,j+1,k+1)) &
                +(cydi*fin(2,i+1,j,k+1)+cyd*fin(2,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hxi(v)*hyi*hz2*( &
                czi*( &
                (fin(3,i,j,k)  -fin(3,i,j+1,k))- &
                (fin(3,i+1,j,k)-fin(3,i+1,j+1,k))) &
                +cz*( &
                (fin(3,i,j,k+1)  -fin(3,i,j+1,k+1))- &
                (fin(3,i+1,j,k+1)-fin(3,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hx(v)*hy*( &
                zpi*( &
                cxdi*(cydi*fin(4,i,j,k)  +cyd*fin(4,i,j+1,k))+ &
                cxd*(cydi*fin(4,i+1,j,k)+cyd*fin(4,i+1,j+1,k))) &
                +zp*( &
                cxdi*(cydi*fin(4,i,j,k+1)  +cyd*fin(4,i,j+1,k+1))+ &
                cxd*(cydi*fin(4,i+1,j,k+1)+cyd*fin(4,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hyi*hx(v)*hz2*( &
                czi*( &
                cxdi*(-fin(5,i,j,k)  +fin(5,i,j+1,k))+ &
                cxd*(-fin(5,i+1,j,k)+fin(5,i+1,j+1,k))) &
                +cz*( &
                cxdi*(-fin(5,i,j,k+1)  +fin(5,i,j+1,k+1))+ &
                cxd*(-fin(5,i+1,j,k+1)+fin(5,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hxi(v)*hy*hz2*( &
                czi*( &
                -(cydi*fin(6,i,j,k)  +cyd*fin(6,i,j+1,k)) &
                +(cydi*fin(6,i+1,j,k)+cyd*fin(6,i+1,j+1,k))) &
                +cz*( &
                -(cydi*fin(6,i,j,k+1)  +cyd*fin(6,i,j+1,k+1)) &
                +(cydi*fin(6,i+1,j,k+1)+cyd*fin(6,i+1,j+1,k+1))))
           !
           sum=sum+z216th*hx(v)*hy*hz2*( &
                czi*( &
                cxdi*(cydi*fin(7,i,j,k)  +cyd*fin(7,i,j+1,k))+ &
                cxd*(cydi*fin(7,i+1,j,k)+cyd*fin(7,i+1,j+1,k))) &
                +cz*( &
                cxdi*(cydi*fin(7,i,j,k+1)  +cyd*fin(7,i,j+1,k+1))+ &
                cxd*(cydi*fin(7,i+1,j,k+1)+cyd*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(9).eq.1) then
        !
        !  d2f/dxdz:
        !
        j=jj
        k=kk
        !
        !   ...and in y direction
        !
        yp=yparam
        ypi=1.0_fp-yp
        yp2=yp*yp
        ypi2=ypi*ypi
        !
        cy=yp*(yp2-1.0_fp)
        cyi=ypi*(ypi2-1.0_fp)
        hy2=hy*hy
        !
        !   ...and in z direction
        !
        zp=zparam
        zpi=1.0_fp-zp
        zp2=zp*zp
        zpi2=zpi*zpi

        czd=3.0_fp*zp2-1.0_fp
        czdi=-3.0_fp*zpi2+1.0_fp
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           !
           !   ...in x direction
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           xp2=xp*xp
           xpi2=xpi*xpi

           cxd=3.0_fp*xp2-1.0_fp
           cxdi=-3.0_fp*xpi2+1.0_fp
           !
           sum=hxi(v)*hzi*( &
                ( &
                (ypi*fin(0,i,j,k)  +yp*fin(0,i,j+1,k)) - &
                (ypi*fin(0,i+1,j,k)+yp*fin(0,i+1,j+1,k))) &
                -( &
                (ypi*fin(0,i,j,k+1)  +yp*fin(0,i,j+1,k+1)) - &
                (ypi*fin(0,i+1,j,k+1)+yp*fin(0,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hx(v)*hzi*( &
                -( &
                cxdi*(ypi*fin(1,i,j,k)  +yp*fin(1,i,j+1,k))+ &
                cxd*(ypi*fin(1,i+1,j,k)+yp*fin(1,i+1,j+1,k))) &
                +( &
                cxdi*(ypi*fin(1,i,j,k+1)  +yp*fin(1,i,j+1,k+1))+ &
                cxd*(ypi*fin(1,i+1,j,k+1)+yp*fin(1,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hxi(v)*hy2*hzi*( &
                ( &
                (cyi*fin(2,i,j,k)  +cy*fin(2,i,j+1,k)) - &
                (cyi*fin(2,i+1,j,k)+cy*fin(2,i+1,j+1,k))) &
                -( &
                (cyi*fin(2,i,j,k+1)  +cy*fin(2,i,j+1,k+1)) - &
                (cyi*fin(2,i+1,j,k+1)+cy*fin(2,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hxi(v)*hz*( &
                czdi*( &
                -(ypi*fin(3,i,j,k)  +yp*fin(3,i,j+1,k)) &
                +(ypi*fin(3,i+1,j,k)+yp*fin(3,i+1,j+1,k))) &
                +czd*( &
                -(ypi*fin(3,i,j,k+1)  +yp*fin(3,i,j+1,k+1)) &
                +(ypi*fin(3,i+1,j,k+1)+yp*fin(3,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hx(v)*hy2*hzi*( &
                -( &
                cxdi*(cyi*fin(4,i,j,k)  +cy*fin(4,i,j+1,k))+ &
                cxd*(cyi*fin(4,i+1,j,k)+cy*fin(4,i+1,j+1,k))) &
                +( &
                cxdi*(cyi*fin(4,i,j,k+1)  +cy*fin(4,i,j+1,k+1))+ &
                cxd*(cyi*fin(4,i+1,j,k+1)+cy*fin(4,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hx(v)*hz*( &
                czdi*( &
                cxdi*(ypi*fin(5,i,j,k)  +yp*fin(5,i,j+1,k))+ &
                cxd*(ypi*fin(5,i+1,j,k)+yp*fin(5,i+1,j+1,k))) &
                +czd*( &
                cxdi*(ypi*fin(5,i,j,k+1)  +yp*fin(5,i,j+1,k+1))+ &
                cxd*(ypi*fin(5,i+1,j,k+1)+yp*fin(5,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hxi(v)*hy2*hz*( &
                czdi*( &
                -(cyi*fin(6,i,j,k)  +cy*fin(6,i,j+1,k)) &
                +(cyi*fin(6,i+1,j,k)+cy*fin(6,i+1,j+1,k))) &
                +czd*( &
                -(cyi*fin(6,i,j,k+1)  +cy*fin(6,i,j+1,k+1)) &
                +(cyi*fin(6,i+1,j,k+1)+cy*fin(6,i+1,j+1,k+1))))
           !
           sum=sum+z216th*hx(v)*hy2*hz*( &
                czdi*( &
                cxdi*(cyi*fin(7,i,j,k)  +cy*fin(7,i,j+1,k))+ &
                cxd*(cyi*fin(7,i+1,j,k)+cy*fin(7,i+1,j+1,k))) &
                +czd*( &
                cxdi*(cyi*fin(7,i,j,k+1)  +cy*fin(7,i,j+1,k+1))+ &
                cxd*(cyi*fin(7,i+1,j,k+1)+cy*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(10).eq.1) then
        !
        !  d2f/dydz:
        !
        j=jj
        k=kk
        !
        !   ...and in y direction
        !
        yp=yparam
        ypi=1.0_fp-yp
        yp2=yp*yp
        ypi2=ypi*ypi

        cyd=3.0_fp*yp2-1.0_fp
        cydi=-3.0_fp*ypi2+1.0_fp
        !
        !   ...and in z direction
        !
        zp=zparam
        zpi=1.0_fp-zp
        zp2=zp*zp
        zpi2=zpi*zpi

        czd=3.0_fp*zp2-1.0_fp
        czdi=-3.0_fp*zpi2+1.0_fp
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           !
           !   ...in x direction
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           xp2=xp*xp
           xpi2=xpi*xpi
           !
           cx=xp*(xp2-1.0_fp)
           cxi=xpi*(xpi2-1.0_fp)
           hx2=hx(v)*hx(v)
           !
           sum=hyi*hzi*( &
                ( &
                xpi*(fin(0,i,j,k)  -fin(0,i,j+1,k))+ &
                xp*(fin(0,i+1,j,k)-fin(0,i+1,j+1,k))) &
                -( &
                xpi*(fin(0,i,j,k+1)  -fin(0,i,j+1,k+1))+ &
                xp*(fin(0,i+1,j,k+1)-fin(0,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hyi*hx2*hzi*( &
                ( &
                cxi*(fin(1,i,j,k)  -fin(1,i,j+1,k))+ &
                cx*(fin(1,i+1,j,k)-fin(1,i+1,j+1,k))) &
                -( &
                cxi*(fin(1,i,j,k+1)  -fin(1,i,j+1,k+1))+ &
                cx*(fin(1,i+1,j,k+1)-fin(1,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hy*hzi*( &
                -( &
                xpi*(cydi*fin(2,i,j,k)  +cyd*fin(2,i,j+1,k))+ &
                xp*(cydi*fin(2,i+1,j,k)+cyd*fin(2,i+1,j+1,k))) &
                +( &
                xpi*(cydi*fin(2,i,j,k+1)  +cyd*fin(2,i,j+1,k+1))+ &
                xp*(cydi*fin(2,i+1,j,k+1)+cyd*fin(2,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hyi*hz*( &
                czdi*( &
                xpi*(-fin(3,i,j,k)  +fin(3,i,j+1,k))+ &
                xp*(-fin(3,i+1,j,k)+fin(3,i+1,j+1,k))) &
                +czd*( &
                xpi*(-fin(3,i,j,k+1)  +fin(3,i,j+1,k+1))+ &
                xp*(-fin(3,i+1,j,k+1)+fin(3,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hx2*hy*hzi*( &
                -( &
                cxi*(cydi*fin(4,i,j,k)  +cyd*fin(4,i,j+1,k))+ &
                cx*(cydi*fin(4,i+1,j,k)+cyd*fin(4,i+1,j+1,k))) &
                +( &
                cxi*(cydi*fin(4,i,j,k+1)  +cyd*fin(4,i,j+1,k+1))+ &
                cx*(cydi*fin(4,i+1,j,k+1)+cyd*fin(4,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hyi*hx2*hz*( &
                czdi*( &
                cxi*(-fin(5,i,j,k)  +fin(5,i,j+1,k))+ &
                cx*(-fin(5,i+1,j,k)+fin(5,i+1,j+1,k))) &
                +czd*( &
                cxi*(-fin(5,i,j,k+1)  +fin(5,i,j+1,k+1))+ &
                cx*(-fin(5,i+1,j,k+1)+fin(5,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hy*hz*( &
                czdi*( &
                xpi*(cydi*fin(6,i,j,k)  +cyd*fin(6,i,j+1,k))+ &
                xp*(cydi*fin(6,i+1,j,k)+cyd*fin(6,i+1,j+1,k))) &
                +czd*( &
                xpi*(cydi*fin(6,i,j,k+1)  +cyd*fin(6,i,j+1,k+1))+ &
                xp*(cydi*fin(6,i+1,j,k+1)+cyd*fin(6,i+1,j+1,k+1))))
           !
           sum=sum+z216th*hx2*hy*hz*( &
                czdi*( &
                cxi*(cydi*fin(7,i,j,k)  +cyd*fin(7,i,j+1,k))+ &
                cx*(cydi*fin(7,i+1,j,k)+cyd*fin(7,i+1,j+1,k))) &
                +czd*( &
                cxi*(cydi*fin(7,i,j,k+1)  +cyd*fin(7,i,j+1,k+1))+ &
                cx*(cydi*fin(7,i+1,j,k+1)+cyd*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     !----------------------------------
     !  3rd derivatives (.le.2 in each coordinate)
     !
  else if(ict(1).eq.3) then
     if(ict(2).eq.1) then
        !                               ! fxxy
        j=jj
        k=kk
        !
        !   ...and in y direction
        !
        yp=yparam
        ypi=1.0_fp-yp
        yp2=yp*yp
        ypi2=ypi*ypi
        !
        cyd=3.0_fp*yp2-1.0_fp
        cydi=-3.0_fp*ypi2+1.0_fp
        !
        !   ...and in z direction
        !
        zp=zparam
        zpi=1.0_fp-zp
        zp2=zp*zp
        zpi2=zpi*zpi
        !
        cz=zp*(zp2-1.0_fp)
        czi=zpi*(zpi2-1.0_fp)
        hz2=hz*hz
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           !
           !   ...in x direction
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           !
           sum=hyi*( &
                zpi*( &
                xpi*(-fin(1,i,j,k)  +fin(1,i,j+1,k))+ &
                xp*( -fin(1,i+1,j,k)+fin(1,i+1,j+1,k))) &
                +zp*( &
                xpi*(-fin(1,i,j,k+1)  +fin(1,i,j+1,k+1))+ &
                xp*( -fin(1,i+1,j,k+1)+fin(1,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hy*( &
                zpi*( &
                xpi*(cydi*fin(4,i,j,k)  +cyd*fin(4,i,j+1,k))+ &
                xp*(cydi*fin(4,i+1,j,k)+cyd*fin(4,i+1,j+1,k))) &
                +zp*( &
                xpi*(cydi*fin(4,i,j,k+1)  +cyd*fin(4,i,j+1,k+1))+ &
                xp*(cydi*fin(4,i+1,j,k+1)+cyd*fin(4,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hz2*hyi*( &
                czi*( &
                xpi*(-fin(5,i,j,k)  +fin(5,i,j+1,k))+ &
                xp*( -fin(5,i+1,j,k)+fin(5,i+1,j+1,k))) &
                +cz*( &
                xpi*(-fin(5,i,j,k+1)  +fin(5,i,j+1,k+1))+ &
                xp*( -fin(5,i+1,j,k+1)+fin(5,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hy*hz2*( &
                czi*( &
                xpi*(cydi*fin(7,i,j,k)  +cyd*fin(7,i,j+1,k))+ &
                xp*(cydi*fin(7,i+1,j,k)+cyd*fin(7,i+1,j+1,k))) &
                +cz*( &
                xpi*(cydi*fin(7,i,j,k+1)  +cyd*fin(7,i,j+1,k+1))+ &
                xp*(cydi*fin(7,i+1,j,k+1)+cyd*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(3).eq.1) then
        !                               ! fxxz
        j=jj
        k=kk
        !
        !   ...and in y direction
        !
        yp=yparam
        ypi=1.0_fp-yp
        yp2=yp*yp
        ypi2=ypi*ypi
        !
        cy=yp*(yp2-1.0_fp)
        cyi=ypi*(ypi2-1.0_fp)
        hy2=hy*hy
        !
        !   ...and in z direction
        !
        zp=zparam
        zpi=1.0_fp-zp
        zp2=zp*zp
        zpi2=zpi*zpi
        !
        czd=3.0_fp*zp2-1.0_fp
        czdi=-3.0_fp*zpi2+1.0_fp
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           !
           !   ...in x direction
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           !
           sum=hzi*( &
                -( &
                xpi*(ypi*fin(1,i,j,k)  +yp*fin(1,i,j+1,k))+ &
                xp*(ypi*fin(1,i+1,j,k)+yp*fin(1,i+1,j+1,k))) &
                +( &
                xpi*(ypi*fin(1,i,j,k+1)  +yp*fin(1,i,j+1,k+1))+ &
                xp*(ypi*fin(1,i+1,j,k+1)+yp*fin(1,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hy2*hzi*( &
                -( &
                xpi*(cyi*fin(4,i,j,k)  +cy*fin(4,i,j+1,k))+ &
                xp*(cyi*fin(4,i+1,j,k)+cy*fin(4,i+1,j+1,k))) &
                +( &
                xpi*(cyi*fin(4,i,j,k+1)  +cy*fin(4,i,j+1,k+1))+ &
                xp*(cyi*fin(4,i+1,j,k+1)+cy*fin(4,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hz*( &
                czdi*( &
                xpi*(ypi*fin(5,i,j,k)  +yp*fin(5,i,j+1,k))+ &
                xp*(ypi*fin(5,i+1,j,k)+yp*fin(5,i+1,j+1,k))) &
                +czd*( &
                xpi*(ypi*fin(5,i,j,k+1)  +yp*fin(5,i,j+1,k+1))+ &
                xp*(ypi*fin(5,i+1,j,k+1)+yp*fin(5,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hy2*hz*( &
                czdi*( &
                xpi*(cyi*fin(7,i,j,k)  +cy*fin(7,i,j+1,k))+ &
                xp*(cyi*fin(7,i+1,j,k)+cy*fin(7,i+1,j+1,k))) &
                +czd*( &
                xpi*(cyi*fin(7,i,j,k+1)  +cy*fin(7,i,j+1,k+1))+ &
                xp*(cyi*fin(7,i+1,j,k+1)+cy*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(4).eq.1) then
        !                               ! fxyy
        j=jj
        k=kk
        !
        !   ...and in y direction
        !
        yp=yparam
        ypi=1.0_fp-yp
        !
        !   ...and in z direction
        !
        zp=zparam
        zpi=1.0_fp-zp
        zp2=zp*zp
        zpi2=zpi*zpi
        !
        cz=zp*(zp2-1.0_fp)
        czi=zpi*(zpi2-1.0_fp)
        hz2=hz*hz
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           !
           !   ...in x direction
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           xp2=xp*xp
           xpi2=xpi*xpi
           !
           cxd=3.0_fp*xp2-1.0_fp
           cxdi=-3.0_fp*xpi2+1.0_fp
           !
           sum=hxi(v)*( &
                zpi*( &
                -(ypi*fin(2,i,j,k)  +yp*fin(2,i,j+1,k)) &
                +(ypi*fin(2,i+1,j,k)+yp*fin(2,i+1,j+1,k))) &
                +zp*( &
                -(ypi*fin(2,i,j,k+1)  +yp*fin(2,i,j+1,k+1)) &
                +(ypi*fin(2,i+1,j,k+1)+yp*fin(2,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hx(v)*( &
                zpi*( &
                cxdi*(ypi*fin(4,i,j,k)  +yp*fin(4,i,j+1,k))+ &
                cxd*(ypi*fin(4,i+1,j,k)+yp*fin(4,i+1,j+1,k))) &
                +zp*( &
                cxdi*(ypi*fin(4,i,j,k+1)  +yp*fin(4,i,j+1,k+1))+ &
                cxd*(ypi*fin(4,i+1,j,k+1)+yp*fin(4,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hz2*hxi(v)*( &
                czi*( &
                -(ypi*fin(6,i,j,k)  +yp*fin(6,i,j+1,k)) &
                +(ypi*fin(6,i+1,j,k)+yp*fin(6,i+1,j+1,k))) &
                +cz*( &
                -(ypi*fin(6,i,j,k+1)  +yp*fin(6,i,j+1,k+1)) &
                +(ypi*fin(6,i+1,j,k+1)+yp*fin(6,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hx(v)*hz2*( &
                czi*( &
                cxdi*(ypi*fin(7,i,j,k)  +yp*fin(7,i,j+1,k))+ &
                cxd*(ypi*fin(7,i+1,j,k)+yp*fin(7,i+1,j+1,k))) &
                +cz*( &
                cxdi*(ypi*fin(7,i,j,k+1)  +yp*fin(7,i,j+1,k+1))+ &
                cxd*(ypi*fin(7,i+1,j,k+1)+yp*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(5).eq.1) then
        !                               ! fxyz
        j=jj
        k=kk
        !
        !   ...and in y direction
        !
        yp=yparam
        ypi=1.0_fp-yp
        yp2=yp*yp
        ypi2=ypi*ypi
        !
        cyd=3.0_fp*yp2-1.0_fp
        cydi=-3.0_fp*ypi2+1.0_fp
        !
        !   ...and in z direction
        !
        zp=zparam
        zpi=1.0_fp-zp
        zp2=zp*zp
        zpi2=zpi*zpi
        !
        czd=3.0_fp*zp2-1.0_fp
        czdi=-3.0_fp*zpi2+1.0_fp
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           !
           !   ...in x direction
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           xp2=xp*xp
           xpi2=xpi*xpi
           !
           cxd=3.0_fp*xp2-1.0_fp
           cxdi=-3.0_fp*xpi2+1.0_fp
           !
           sum=hxi(v)*hyi*hzi*( &
                -( &
                (fin(0,i,j,k)  -fin(0,i,j+1,k))- &
                (fin(0,i+1,j,k)-fin(0,i+1,j+1,k))) &
                +( &
                (fin(0,i,j,k+1)  -fin(0,i,j+1,k+1))- &
                (fin(0,i+1,j,k+1)-fin(0,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hyi*hx(v)*hzi*( &
                -( &
                cxdi*(-fin(1,i,j,k)  +fin(1,i,j+1,k))+ &
                cxd*(-fin(1,i+1,j,k)+fin(1,i+1,j+1,k))) &
                +( &
                cxdi*(-fin(1,i,j,k+1)  +fin(1,i,j+1,k+1))+ &
                cxd*(-fin(1,i+1,j,k+1)+fin(1,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hxi(v)*hy*hzi*( &
                -( &
                -(cydi*fin(2,i,j,k)  +cyd*fin(2,i,j+1,k)) &
                +(cydi*fin(2,i+1,j,k)+cyd*fin(2,i+1,j+1,k))) &
                +( &
                -(cydi*fin(2,i,j,k+1)  +cyd*fin(2,i,j+1,k+1)) &
                +(cydi*fin(2,i+1,j,k+1)+cyd*fin(2,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hxi(v)*hyi*hz*( &
                czdi*( &
                (fin(3,i,j,k)  -fin(3,i,j+1,k))- &
                (fin(3,i+1,j,k)-fin(3,i+1,j+1,k))) &
                +czd*( &
                (fin(3,i,j,k+1)  -fin(3,i,j+1,k+1))- &
                (fin(3,i+1,j,k+1)-fin(3,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hx(v)*hy*hzi*( &
                -( &
                cxdi*(cydi*fin(4,i,j,k)  +cyd*fin(4,i,j+1,k))+ &
                cxd*(cydi*fin(4,i+1,j,k)+cyd*fin(4,i+1,j+1,k))) &
                +( &
                cxdi*(cydi*fin(4,i,j,k+1)  +cyd*fin(4,i,j+1,k+1))+ &
                cxd*(cydi*fin(4,i+1,j,k+1)+cyd*fin(4,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hyi*hx(v)*hz*( &
                czdi*( &
                cxdi*(-fin(5,i,j,k)  +fin(5,i,j+1,k))+ &
                cxd*(-fin(5,i+1,j,k)+fin(5,i+1,j+1,k))) &
                +czd*( &
                cxdi*(-fin(5,i,j,k+1)  +fin(5,i,j+1,k+1))+ &
                cxd*(-fin(5,i+1,j,k+1)+fin(5,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hxi(v)*hy*hz*( &
                czdi*( &
                -(cydi*fin(6,i,j,k)  +cyd*fin(6,i,j+1,k)) &
                +(cydi*fin(6,i+1,j,k)+cyd*fin(6,i+1,j+1,k))) &
                +czd*( &
                -(cydi*fin(6,i,j,k+1)  +cyd*fin(6,i,j+1,k+1)) &
                +(cydi*fin(6,i+1,j,k+1)+cyd*fin(6,i+1,j+1,k+1))))
           !
           sum=sum+z216th*hx(v)*hy*hz*( &
                czdi*( &
                cxdi*(cydi*fin(7,i,j,k)  +cyd*fin(7,i,j+1,k))+ &
                cxd*(cydi*fin(7,i+1,j,k)+cyd*fin(7,i+1,j+1,k))) &
                +czd*( &
                cxdi*(cydi*fin(7,i,j,k+1)  +cyd*fin(7,i,j+1,k+1))+ &
                cxd*(cydi*fin(7,i+1,j,k+1)+cyd*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(6).eq.1) then
        !                               ! fxzz
        j=jj
        k=kk
        !
        !   ...and in y direction
        !
        yp=yparam
        ypi=1.0_fp-yp
        yp2=yp*yp
        ypi2=ypi*ypi
        !
        cy=yp*(yp2-1.0_fp)
        cyi=ypi*(ypi2-1.0_fp)
        hy2=hy*hy

        !
        !   ...and in z direction
        !
        zp=zparam
        zpi=1.0_fp-zp
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           !
           !   ...in x direction
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           xp2=xp*xp
           xpi2=xpi*xpi
           !
           cxd=3.0_fp*xp2-1.0_fp
           cxdi=-3.0_fp*xpi2+1.0_fp
           !
           sum=hxi(v)*( &
                zpi*( &
                -(ypi*fin(3,i,j,k)  +yp*fin(3,i,j+1,k)) &
                +(ypi*fin(3,i+1,j,k)+yp*fin(3,i+1,j+1,k))) &
                +zp*( &
                -(ypi*fin(3,i,j,k+1)  +yp*fin(3,i,j+1,k+1)) &
                +(ypi*fin(3,i+1,j,k+1)+yp*fin(3,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hx(v)*( &
                zpi*( &
                cxdi*(ypi*fin(5,i,j,k)  +yp*fin(5,i,j+1,k))+ &
                cxd*(ypi*fin(5,i+1,j,k)+yp*fin(5,i+1,j+1,k))) &
                +zp*( &
                cxdi*(ypi*fin(5,i,j,k+1)  +yp*fin(5,i,j+1,k+1))+ &
                cxd*(ypi*fin(5,i+1,j,k+1)+yp*fin(5,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hy2*hxi(v)*( &
                zpi*( &
                -(cyi*fin(6,i,j,k)  +cy*fin(6,i,j+1,k)) &
                +(cyi*fin(6,i+1,j,k)+cy*fin(6,i+1,j+1,k))) &
                +zp*( &
                -(cyi*fin(6,i,j,k+1)  +cy*fin(6,i,j+1,k+1)) &
                +(cyi*fin(6,i+1,j,k+1)+cy*fin(6,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hx(v)*hy2*( &
                zpi*( &
                cxdi*(cyi*fin(7,i,j,k)  +cy*fin(7,i,j+1,k))+ &
                cxd*(cyi*fin(7,i+1,j,k)+cy*fin(7,i+1,j+1,k))) &
                +zp*( &
                cxdi*(cyi*fin(7,i,j,k+1)  +cy*fin(7,i,j+1,k+1))+ &
                cxd*(cyi*fin(7,i+1,j,k+1)+cy*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(7).eq.1) then
        !                               ! fyyz
        j=jj
        k=kk
        !
        !
        !   ...and in y direction
        !
        yp=yparam
        ypi=1.0_fp-yp
        !
        !   ...and in z direction
        !
        zp=zparam
        zpi=1.0_fp-zp
        zp2=zp*zp
        zpi2=zpi*zpi
        !
        czd=3.0_fp*zp2-1.0_fp
        czdi=-3.0_fp*zpi2+1.0_fp
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           !
           !   ...in x direction
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           xp2=xp*xp
           xpi2=xpi*xpi
           !
           cx=xp*(xp2-1.0_fp)
           cxi=xpi*(xpi2-1.0_fp)
           hx2=hx(v)*hx(v)
           !
           sum=hzi*( &
                -( &
                xpi*(ypi*fin(2,i,j,k)  +yp*fin(2,i,j+1,k))+ &
                xp*(ypi*fin(2,i+1,j,k)+yp*fin(2,i+1,j+1,k))) &
                +( &
                xpi*(ypi*fin(2,i,j,k+1)  +yp*fin(2,i,j+1,k+1))+ &
                xp*(ypi*fin(2,i+1,j,k+1)+yp*fin(2,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hx2*hzi*( &
                -( &
                cxi*(ypi*fin(4,i,j,k)  +yp*fin(4,i,j+1,k))+ &
                cx*(ypi*fin(4,i+1,j,k)+yp*fin(4,i+1,j+1,k))) &
                +( &
                cxi*(ypi*fin(4,i,j,k+1)  +yp*fin(4,i,j+1,k+1))+ &
                cx*(ypi*fin(4,i+1,j,k+1)+yp*fin(4,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hz*( &
                czdi*( &
                xpi*(ypi*fin(6,i,j,k)  +yp*fin(6,i,j+1,k))+ &
                xp*(ypi*fin(6,i+1,j,k)+yp*fin(6,i+1,j+1,k))) &
                +czd*( &
                xpi*(ypi*fin(6,i,j,k+1)  +yp*fin(6,i,j+1,k+1))+ &
                xp*(ypi*fin(6,i+1,j,k+1)+yp*fin(6,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hx2*hz*( &
                czdi*( &
                cxi*(ypi*fin(7,i,j,k)  +yp*fin(7,i,j+1,k))+ &
                cx*(ypi*fin(7,i+1,j,k)+yp*fin(7,i+1,j+1,k))) &
                +czd*( &
                cxi*(ypi*fin(7,i,j,k+1)  +yp*fin(7,i,j+1,k+1))+ &
                cx*(ypi*fin(7,i+1,j,k+1)+yp*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(8).eq.1) then
        !                               ! fyzz
        j=jj
        k=kk
        !
        !   ...and in y direction
        !
        yp=yparam
        ypi=1.0_fp-yp
        yp2=yp*yp
        ypi2=ypi*ypi
        !
        cyd=3.0_fp*yp2-1.0_fp
        cydi=-3.0_fp*ypi2+1.0_fp
        !
        !   ...and in z direction
        !
        zp=zparam
        zpi=1.0_fp-zp
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           !
           !   ...in x direction
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           xp2=xp*xp
           xpi2=xpi*xpi
           !
           cx=xp*(xp2-1.0_fp)
           cxi=xpi*(xpi2-1.0_fp)
           hx2=hx(v)*hx(v)
           !
           sum=hyi*( &
                zpi*( &
                xpi*(-fin(3,i,j,k)  +fin(3,i,j+1,k))+ &
                xp*( -fin(3,i+1,j,k)+fin(3,i+1,j+1,k))) &
                +zp*( &
                xpi*(-fin(3,i,j,k+1)  +fin(3,i,j+1,k+1))+ &
                xp*( -fin(3,i+1,j,k+1)+fin(3,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hx2*hyi*( &
                zpi*( &
                cxi*(-fin(5,i,j,k)  +fin(5,i,j+1,k))+ &
                cx*( -fin(5,i+1,j,k)+fin(5,i+1,j+1,k))) &
                +zp*( &
                cxi*(-fin(5,i,j,k+1)  +fin(5,i,j+1,k+1))+ &
                cx*( -fin(5,i+1,j,k+1)+fin(5,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hy*( &
                zpi*( &
                xpi*(cydi*fin(6,i,j,k) +cyd*fin(6,i,j+1,k))+ &
                xp*(cydi*fin(6,i+1,j,k)+cyd*fin(6,i+1,j+1,k))) &
                +zp*( &
                xpi*(cydi*fin(6,i,j,k+1) +cyd*fin(6,i,j+1,k+1))+ &
                xp*(cydi*fin(6,i+1,j,k+1)+cyd*fin(6,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hx2*hy*( &
                zpi*( &
                cxi*(cydi*fin(7,i,j,k) +cyd*fin(7,i,j+1,k))+ &
                cx*(cydi*fin(7,i+1,j,k)+cyd*fin(7,i+1,j+1,k))) &
                +zp*( &
                cxi*(cydi*fin(7,i,j,k+1) +cyd*fin(7,i,j+1,k+1))+ &
                cx*(cydi*fin(7,i+1,j,k+1)+cyd*fin(7,i+1,j+1,k+1))))
           fval(v,iadr)=sum
        end do
     end if
     !
     !----------------------------------
     !  3rd derivatives (3 in each coordinate)
     !
  else if(ict(1).eq.-3) then
     if(ict(2).eq.1) then
        !                               ! fxxx
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           j=jj
           k=kk
           !
           !   ...and in y direction
           !
           yp=yparam
           ypi=1.0_fp-yp
           yp2=yp*yp
           ypi2=ypi*ypi
           !
           cy=yp*(yp2-1.0_fp)
           cyi=ypi*(ypi2-1.0_fp)
           hy2=hy*hy
           !
           !   ...and in z direction
           !
           zp=zparam
           zpi=1.0_fp-zp
           zp2=zp*zp
           zpi2=zpi*zpi
           !
           cz=zp*(zp2-1.0_fp)
           czi=zpi*(zpi2-1.0_fp)
           hz2=hz*hz
           !
           sum=hxi(v)*( &
                zpi*( &
                -(ypi*fin(1,i,j,k)  +yp*fin(1,i,j+1,k)) &
                +(ypi*fin(1,i+1,j,k)+yp*fin(1,i+1,j+1,k))) &
                +zp*( &
                -(ypi*fin(1,i,j,k+1)  +yp*fin(1,i,j+1,k+1)) &
                +(ypi*fin(1,i+1,j,k+1)+yp*fin(1,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hy2*hxi(v)*( &
                zpi*( &
                -(cyi*fin(4,i,j,k)  +cy*fin(4,i,j+1,k)) &
                +(cyi*fin(4,i+1,j,k)+cy*fin(4,i+1,j+1,k))) &
                +zp*( &
                -(cyi*fin(4,i,j,k+1)  +cy*fin(4,i,j+1,k+1)) &
                +(cyi*fin(4,i+1,j,k+1)+cy*fin(4,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hz2*hxi(v)*( &
                czi*( &
                -(ypi*fin(5,i,j,k)  +yp*fin(5,i,j+1,k)) &
                +(ypi*fin(5,i+1,j,k)+yp*fin(5,i+1,j+1,k))) &
                +cz*( &
                -(ypi*fin(5,i,j,k+1)  +yp*fin(5,i,j+1,k+1)) &
                +(ypi*fin(5,i+1,j,k+1)+yp*fin(5,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hy2*hz2*hxi(v)*( &
                czi*( &
                -(cyi*fin(7,i,j,k)  +cy*fin(7,i,j+1,k)) &
                +(cyi*fin(7,i+1,j,k)+cy*fin(7,i+1,j+1,k))) &
                +cz*( &
                -(cyi*fin(7,i,j,k+1)  +cy*fin(7,i,j+1,k+1)) &
                +(cyi*fin(7,i+1,j,k+1)+cy*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(3).eq.1) then
        !                               ! fyyy
        j=jj
        k=kk
        !
        !   ...and in z direction
        !
        zp=zparam
        zpi=1.0_fp-zp
        zp2=zp*zp
        zpi2=zpi*zpi
        !
        cz=zp*(zp2-1.0_fp)
        czi=zpi*(zpi2-1.0_fp)
        hz2=hz*hz
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           !
           !   ...in x direction
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           xp2=xp*xp
           xpi2=xpi*xpi
           !
           cx=xp*(xp2-1.0_fp)
           cxi=xpi*(xpi2-1.0_fp)
           hx2=hx(v)*hx(v)
           !
           sum=hyi*( &
                zpi*( &
                xpi*(-fin(2,i,j,k)  +fin(2,i,j+1,k))+ &
                xp*( -fin(2,i+1,j,k)+fin(2,i+1,j+1,k))) &
                +zp*( &
                xpi*(-fin(2,i,j,k+1)  +fin(2,i,j+1,k+1))+ &
                xp*( -fin(2,i+1,j,k+1)+fin(2,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hx2*hyi*( &
                zpi*( &
                cxi*(-fin(4,i,j,k)  +fin(4,i,j+1,k))+ &
                cx*( -fin(4,i+1,j,k)+fin(4,i+1,j+1,k))) &
                +zp*( &
                cxi*(-fin(4,i,j,k+1)  +fin(4,i,j+1,k+1))+ &
                cx*( -fin(4,i+1,j,k+1)+fin(4,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hz2*hyi*( &
                czi*( &
                xpi*(-fin(6,i,j,k)  +fin(6,i,j+1,k))+ &
                xp*( -fin(6,i+1,j,k)+fin(6,i+1,j+1,k))) &
                +cz*( &
                xpi*(-fin(6,i,j,k+1)  +fin(6,i,j+1,k+1))+ &
                xp*( -fin(6,i+1,j,k+1)+fin(6,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hx2*hz2*hyi*( &
                czi*( &
                cxi*(-fin(7,i,j,k)  +fin(7,i,j+1,k))+ &
                cx*( -fin(7,i+1,j,k)+fin(7,i+1,j+1,k))) &
                +cz*( &
                cxi*(-fin(7,i,j,k+1)  +fin(7,i,j+1,k+1))+ &
                cx*( -fin(7,i+1,j,k+1)+fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(4).eq.1) then
        !                               ! fzzz
        j=jj
        k=kk
        !
        !   ...and in y direction
        !
        yp=yparam
        ypi=1.0_fp-yp
        yp2=yp*yp
        ypi2=ypi*ypi
        !
        cy=yp*(yp2-1.0_fp)
        cyi=ypi*(ypi2-1.0_fp)
        hy2=hy*hy
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           !
           !   ...in x direction
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           xp2=xp*xp
           xpi2=xpi*xpi
           !
           cx=xp*(xp2-1.0_fp)
           cxi=xpi*(xpi2-1.0_fp)
           hx2=hx(v)*hx(v)
           !
           sum=hzi*( &
                -( &
                xpi*(ypi*fin(3,i,j,k)  +yp*fin(3,i,j+1,k))+ &
                xp*(ypi*fin(3,i+1,j,k)+yp*fin(3,i+1,j+1,k))) &
                +( &
                xpi*(ypi*fin(3,i,j,k+1)  +yp*fin(3,i,j+1,k+1))+ &
                xp*(ypi*fin(3,i+1,j,k+1)+yp*fin(3,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hx2*hzi*( &
                -( &
                cxi*(ypi*fin(5,i,j,k)  +yp*fin(5,i,j+1,k))+ &
                cx*(ypi*fin(5,i+1,j,k)+yp*fin(5,i+1,j+1,k))) &
                +( &
                cxi*(ypi*fin(5,i,j,k+1)  +yp*fin(5,i,j+1,k+1))+ &
                cx*(ypi*fin(5,i+1,j,k+1)+yp*fin(5,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hy2*hzi*( &
                -( &
                xpi*(cyi*fin(6,i,j,k)  +cy*fin(6,i,j+1,k))+ &
                xp*(cyi*fin(6,i+1,j,k)+cy*fin(6,i+1,j+1,k))) &
                +( &
                xpi*(cyi*fin(6,i,j,k+1)  +cy*fin(6,i,j+1,k+1))+ &
                xp*(cyi*fin(6,i+1,j,k+1)+cy*fin(6,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hx2*hy2*hzi*( &
                -( &
                cxi*(cyi*fin(7,i,j,k)  +cy*fin(7,i,j+1,k))+ &
                cx*(cyi*fin(7,i+1,j,k)+cy*fin(7,i+1,j+1,k))) &
                +( &
                cxi*(cyi*fin(7,i,j,k+1)  +cy*fin(7,i,j+1,k+1))+ &
                cx*(cyi*fin(7,i+1,j,k+1)+cy*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     !----------------------------------
     !  4th derivatives (.le.2 in each coordinate)
     !
  else if(ict(1).eq.4) then
     if(ict(2).eq.1) then
        !                               ! fxxyy
        j=jj
        k=kk
        !
        !   ...and in y direction
        !
        yp=yparam
        ypi=1.0_fp-yp
        !
        !   ...and in z direction
        !
        zp=zparam
        zpi=1.0_fp-zp
        zp2=zp*zp
        zpi2=zpi*zpi
        !
        cz=zp*(zp2-1.0_fp)
        czi=zpi*(zpi2-1.0_fp)
        hz2=hz*hz
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           !
           !   ...in x direction
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           !
           sum=( &
                zpi*( &
                xpi*(ypi*fin(4,i,j,k)  +yp*fin(4,i,j+1,k))+ &
                xp*( ypi*fin(4,i+1,j,k)+yp*fin(4,i+1,j+1,k))) &
                +zp*( &
                xpi*(ypi*fin(4,i,j,k+1) +yp*fin(4,i,j+1,k+1))+ &
                xp*(ypi*fin(4,i+1,j,k+1)+yp*fin(4,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hz2*( &
                czi*( &
                xpi*(ypi*fin(7,i,j,k) +yp*fin(7,i,j+1,k))+ &
                xp*(ypi*fin(7,i+1,j,k)+yp*fin(7,i+1,j+1,k))) &
                +cz*( &
                xpi*(ypi*fin(7,i,j,k+1) +yp*fin(7,i,j+1,k+1))+ &
                xp*(ypi*fin(7,i+1,j,k+1)+yp*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(3).eq.1) then
        !                               ! fxxyz
        j=jj
        k=kk
        !
        !   ...and in y direction
        !
        yp=yparam
        ypi=1.0_fp-yp
        yp2=yp*yp
        ypi2=ypi*ypi
        !
        cyd=3.0_fp*yp2-1.0_fp
        cydi=-3.0_fp*ypi2+1.0_fp
        !
        !   ...and in z direction
        !
        zp=zparam
        zpi=1.0_fp-zp
        zp2=zp*zp
        zpi2=zpi*zpi
        !
        czd=3.0_fp*zp2-1.0_fp
        czdi=-3.0_fp*zpi2+1.0_fp
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           !
           !   ...in x direction
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           !
           sum=hyi*hzi*( &
                -( &
                xpi*(-fin(1,i,j,k)  +fin(1,i,j+1,k))+ &
                xp*( -fin(1,i+1,j,k)+fin(1,i+1,j+1,k))) &
                +( &
                xpi*(-fin(1,i,j,k+1)  +fin(1,i,j+1,k+1))+ &
                xp*( -fin(1,i+1,j,k+1)+fin(1,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hy*hzi*( &
                -( &
                xpi*(cydi*fin(4,i,j,k)  +cyd*fin(4,i,j+1,k))+ &
                xp*(cydi*fin(4,i+1,j,k)+cyd*fin(4,i+1,j+1,k))) &
                +( &
                xpi*(cydi*fin(4,i,j,k+1)  +cyd*fin(4,i,j+1,k+1))+ &
                xp*(cydi*fin(4,i+1,j,k+1)+cyd*fin(4,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hz*hyi*( &
                czdi*( &
                xpi*(-fin(5,i,j,k)  +fin(5,i,j+1,k))+ &
                xp*( -fin(5,i+1,j,k)+fin(5,i+1,j+1,k))) &
                +czd*( &
                xpi*(-fin(5,i,j,k+1)  +fin(5,i,j+1,k+1))+ &
                xp*( -fin(5,i+1,j,k+1)+fin(5,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hy*hz*( &
                czdi*( &
                xpi*(cydi*fin(7,i,j,k)  +cyd*fin(7,i,j+1,k))+ &
                xp*(cydi*fin(7,i+1,j,k)+cyd*fin(7,i+1,j+1,k))) &
                +czd*( &
                xpi*(cydi*fin(7,i,j,k+1)  +cyd*fin(7,i,j+1,k+1))+ &
                xp*(cydi*fin(7,i+1,j,k+1)+cyd*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(4).eq.1) then
        !                               ! fxxzz
        j=jj
        k=kk
        !
        !   ...and in y direction
        !
        yp=yparam
        ypi=1.0_fp-yp
        yp2=yp*yp
        ypi2=ypi*ypi
        !
        cy=yp*(yp2-1.0_fp)
        cyi=ypi*(ypi2-1.0_fp)
        hy2=hy*hy
        !
        !   ...and in z direction
        !
        zp=zparam
        zpi=1.0_fp-zp
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           !
           !   ...in x direction
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           !
           sum=( &
                zpi*( &
                xpi*(ypi*fin(5,i,j,k)  +yp*fin(5,i,j+1,k))+ &
                xp*(ypi*fin(5,i+1,j,k)+yp*fin(5,i+1,j+1,k))) &
                +zp*( &
                xpi*(ypi*fin(5,i,j,k+1)  +yp*fin(5,i,j+1,k+1))+ &
                xp*(ypi*fin(5,i+1,j,k+1)+yp*fin(5,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hy2*( &
                zpi*( &
                xpi*(cyi*fin(7,i,j,k)  +cy*fin(7,i,j+1,k))+ &
                xp*(cyi*fin(7,i+1,j,k)+cy*fin(7,i+1,j+1,k))) &
                +zp*( &
                xpi*(cyi*fin(7,i,j,k+1)  +cy*fin(7,i,j+1,k+1))+ &
                xp*(cyi*fin(7,i+1,j,k+1)+cy*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(5).eq.1) then
        !                               ! fxyyz
        j=jj
        k=kk
        !
        !   ...and in y direction
        !
        yp=yparam
        ypi=1.0_fp-yp
        !
        !   ...and in z direction
        !
        zp=zparam
        zpi=1.0_fp-zp
        zp2=zp*zp
        zpi2=zpi*zpi

        czd=3.0_fp*zp2-1.0_fp
        czdi=-3.0_fp*zpi2+1.0_fp
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           !
           !   ...in x direction
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           xp2=xp*xp
           xpi2=xpi*xpi
           !
           cxd=3.0_fp*xp2-1.0_fp
           cxdi=-3.0_fp*xpi2+1.0_fp
           !
           sum=hxi(v)*hzi*( &
                -( &
                -(ypi*fin(2,i,j,k)  +yp*fin(2,i,j+1,k)) &
                +(ypi*fin(2,i+1,j,k)+yp*fin(2,i+1,j+1,k))) &
                +( &
                -(ypi*fin(2,i,j,k+1)  +yp*fin(2,i,j+1,k+1)) &
                +(ypi*fin(2,i+1,j,k+1)+yp*fin(2,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hx(v)*hzi*( &
                -( &
                cxdi*(ypi*fin(4,i,j,k)  +yp*fin(4,i,j+1,k))+ &
                cxd*(ypi*fin(4,i+1,j,k)+yp*fin(4,i+1,j+1,k))) &
                +( &
                cxdi*(ypi*fin(4,i,j,k+1)  +yp*fin(4,i,j+1,k+1))+ &
                cxd*(ypi*fin(4,i+1,j,k+1)+yp*fin(4,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hz*hxi(v)*( &
                czdi*( &
                -(ypi*fin(6,i,j,k)  +yp*fin(6,i,j+1,k)) &
                +(ypi*fin(6,i+1,j,k)+yp*fin(6,i+1,j+1,k))) &
                +czd*( &
                -(ypi*fin(6,i,j,k+1)  +yp*fin(6,i,j+1,k+1)) &
                +(ypi*fin(6,i+1,j,k+1)+yp*fin(6,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hx(v)*hz*( &
                czdi*( &
                cxdi*(ypi*fin(7,i,j,k)  +yp*fin(7,i,j+1,k))+ &
                cxd*(ypi*fin(7,i+1,j,k)+yp*fin(7,i+1,j+1,k))) &
                +czd*( &
                cxdi*(ypi*fin(7,i,j,k+1)  +yp*fin(7,i,j+1,k+1))+ &
                cxd*(ypi*fin(7,i+1,j,k+1)+yp*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(6).eq.1) then
        !                               ! fxyzz
        j=jj
        k=kk
        !
        !   ...and in y direction
        !
        yp=yparam
        ypi=1.0_fp-yp
        yp2=yp*yp
        ypi2=ypi*ypi

        cyd=3.0_fp*yp2-1.0_fp
        cydi=-3.0_fp*ypi2+1.0_fp
        !
        !   ...and in z direction
        !
        zp=zparam
        zpi=1.0_fp-zp
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           !
           !   ...in x direction
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           xp2=xp*xp
           xpi2=xpi*xpi

           cxd=3.0_fp*xp2-1.0_fp
           cxdi=-3.0_fp*xpi2+1.0_fp
           !
           sum=hyi*hxi(v)*( &
                zpi*( &
                ( +fin(3,i,j,k)  -fin(3,i,j+1,k)) &
                +( -fin(3,i+1,j,k)+fin(3,i+1,j+1,k))) &
                +zp*( &
                ( +fin(3,i,j,k+1)  -fin(3,i,j+1,k+1)) &
                +( -fin(3,i+1,j,k+1)+fin(3,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hx(v)*hyi*( &
                zpi*( &
                cxdi*(-fin(5,i,j,k)  +fin(5,i,j+1,k))+ &
                cxd*( -fin(5,i+1,j,k)+fin(5,i+1,j+1,k))) &
                +zp*( &
                cxdi*(-fin(5,i,j,k+1)  +fin(5,i,j+1,k+1))+ &
                cxd*( -fin(5,i+1,j,k+1)+fin(5,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hy*hxi(v)*( &
                zpi*( &
                -(cydi*fin(6,i,j,k) +cyd*fin(6,i,j+1,k)) &
                +(cydi*fin(6,i+1,j,k)+cyd*fin(6,i+1,j+1,k))) &
                +zp*( &
                -(cydi*fin(6,i,j,k+1) +cyd*fin(6,i,j+1,k+1)) &
                +(cydi*fin(6,i+1,j,k+1)+cyd*fin(6,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hx(v)*hy*( &
                zpi*( &
                cxdi*(cydi*fin(7,i,j,k) +cyd*fin(7,i,j+1,k))+ &
                cxd*(cydi*fin(7,i+1,j,k)+cyd*fin(7,i+1,j+1,k))) &
                +zp*( &
                cxdi*(cydi*fin(7,i,j,k+1) +cyd*fin(7,i,j+1,k+1))+ &
                cxd*(cydi*fin(7,i+1,j,k+1)+cyd*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(7).eq.1) then
        !                               ! fyyzz
        j=jj
        k=kk
        !
        !   ...and in y direction
        !
        yp=yparam
        ypi=1.0_fp-yp
        !
        !   ...and in z direction
        !
        zp=zparam
        zpi=1.0_fp-zp
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           !
           !   ...in x direction
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           xp2=xp*xp
           xpi2=xpi*xpi
           !
           cx=xp*(xp2-1.0_fp)
           cxi=xpi*(xpi2-1.0_fp)
           hx2=hx(v)*hx(v)
           !
           sum=( &
                zpi*( &
                xpi*(ypi*fin(6,i,j,k)  +yp*fin(6,i,j+1,k))+ &
                xp*( ypi*fin(6,i+1,j,k)+yp*fin(6,i+1,j+1,k))) &
                +zp*( &
                xpi*(ypi*fin(6,i,j,k+1) +yp*fin(6,i,j+1,k+1))+ &
                xp*(ypi*fin(6,i+1,j,k+1)+yp*fin(6,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hx2*( &
                zpi*( &
                cxi*(ypi*fin(7,i,j,k)  +yp*fin(7,i,j+1,k))+ &
                cx*( ypi*fin(7,i+1,j,k)+yp*fin(7,i+1,j+1,k))) &
                +zp*( &
                cxi*(ypi*fin(7,i,j,k+1)  +yp*fin(7,i,j+1,k+1))+ &
                cx*( ypi*fin(7,i+1,j,k+1)+yp*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     !----------------------------------
     !  4th derivatives (3 in a coordinate)
     !
  else if(ict(1).eq.-4) then
     if(ict(2).eq.1) then
        !                               ! fxxxy
        j=jj
        k=kk
        !
        !   ...and in y direction
        !
        yp=yparam
        ypi=1.0_fp-yp
        yp2=yp*yp
        ypi2=ypi*ypi

        cyd=3.0_fp*yp2-1.0_fp
        cydi=-3.0_fp*ypi2+1.0_fp
        !
        !   ...and in z direction
        !
        zp=zparam
        zpi=1.0_fp-zp
        zp2=zp*zp
        zpi2=zpi*zpi
        !
        cz=zp*(zp2-1.0_fp)
        czi=zpi*(zpi2-1.0_fp)
        hz2=hz*hz
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           !
           sum=hyi*hxi(v)*( &
                zpi*( &
                (  fin(1,i,j,k)  -fin(1,i,j+1,k))+ &
                ( -fin(1,i+1,j,k)+fin(1,i+1,j+1,k))) &
                +zp*( &
                (  fin(1,i,j,k+1)  -fin(1,i,j+1,k+1))+ &
                ( -fin(1,i+1,j,k+1)+fin(1,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hy*hxi(v)*( &
                zpi*( &
                -(cydi*fin(4,i,j,k)  +cyd*fin(4,i,j+1,k))+ &
                (cydi*fin(4,i+1,j,k)+cyd*fin(4,i+1,j+1,k))) &
                +zp*( &
                -(cydi*fin(4,i,j,k+1)  +cyd*fin(4,i,j+1,k+1))+ &
                (cydi*fin(4,i+1,j,k+1)+cyd*fin(4,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hz2*hyi*hxi(v)*( &
                czi*( &
                (  fin(5,i,j,k)  -fin(5,i,j+1,k))+ &
                ( -fin(5,i+1,j,k)+fin(5,i+1,j+1,k))) &
                +cz*( &
                (  fin(5,i,j,k+1)  -fin(5,i,j+1,k+1))+ &
                ( -fin(5,i+1,j,k+1)+fin(5,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hy*hz2*hxi(v)*( &
                czi*( &
                -(cydi*fin(7,i,j,k)  +cyd*fin(7,i,j+1,k))+ &
                (cydi*fin(7,i+1,j,k)+cyd*fin(7,i+1,j+1,k))) &
                +cz*( &
                -(cydi*fin(7,i,j,k+1)  +cyd*fin(7,i,j+1,k+1))+ &
                (cydi*fin(7,i+1,j,k+1)+cyd*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(3).eq.1) then
        !                               ! fxxxz
        j=jj
        k=kk
        !
        !   ...and in y direction
        !
        yp=yparam
        ypi=1.0_fp-yp
        yp2=yp*yp
        ypi2=ypi*ypi
        !
        cy=yp*(yp2-1.0_fp)
        cyi=ypi*(ypi2-1.0_fp)
        hy2=hy*hy
        !
        !   ...and in z direction
        !
        zp=zparam
        zpi=1.0_fp-zp
        zp2=zp*zp
        zpi2=zpi*zpi

        czd=3.0_fp*zp2-1.0_fp
        czdi=-3.0_fp*zpi2+1.0_fp
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           !
           sum=hzi*hxi(v)*( &
                ( &
                +(ypi*fin(1,i,j,k)  +yp*fin(1,i,j+1,k)) &
                -(ypi*fin(1,i+1,j,k)+yp*fin(1,i+1,j+1,k))) &
                +( &
                -(ypi*fin(1,i,j,k+1)  +yp*fin(1,i,j+1,k+1)) &
                +(ypi*fin(1,i+1,j,k+1)+yp*fin(1,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hy2*hzi*hxi(v)*( &
                ( &
                +(cyi*fin(4,i,j,k)  +cy*fin(4,i,j+1,k)) &
                -(cyi*fin(4,i+1,j,k)+cy*fin(4,i+1,j+1,k))) &
                +( &
                -(cyi*fin(4,i,j,k+1)  +cy*fin(4,i,j+1,k+1)) &
                +(cyi*fin(4,i+1,j,k+1)+cy*fin(4,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hz*hxi(v)*( &
                czdi*( &
                -(ypi*fin(5,i,j,k)  +yp*fin(5,i,j+1,k)) &
                +(ypi*fin(5,i+1,j,k)+yp*fin(5,i+1,j+1,k))) &
                +czd*( &
                -(ypi*fin(5,i,j,k+1)  +yp*fin(5,i,j+1,k+1)) &
                +(ypi*fin(5,i+1,j,k+1)+yp*fin(5,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hy2*hz*hxi(v)*( &
                czdi*( &
                -(cyi*fin(7,i,j,k)  +cy*fin(7,i,j+1,k)) &
                +(cyi*fin(7,i+1,j,k)+cy*fin(7,i+1,j+1,k))) &
                +czd*( &
                -(cyi*fin(7,i,j,k+1)  +cy*fin(7,i,j+1,k+1)) &
                +(cyi*fin(7,i+1,j,k+1)+cy*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(4).eq.1) then
        !                               ! fxyyy
        j=jj
        k=kk
        !
        !   ...and in z direction
        !
        zp=zparam
        zpi=1.0_fp-zp
        zp2=zp*zp
        zpi2=zpi*zpi
        !
        cz=zp*(zp2-1.0_fp)
        czi=zpi*(zpi2-1.0_fp)
        hz2=hz*hz
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           !
           !   ...in x direction
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           xp2=xp*xp
           xpi2=xpi*xpi

           cxd=3.0_fp*xp2-1.0_fp
           cxdi=-3.0_fp*xpi2+1.0_fp
           !
           sum=hxi(v)*hyi*( &
                zpi*( &
                ( fin(2,i,j,k)  -fin(2,i,j+1,k)) &
                +(-fin(2,i+1,j,k)+fin(2,i+1,j+1,k))) &
                +zp*( &
                ( fin(2,i,j,k+1)  -fin(2,i,j+1,k+1)) &
                +(-fin(2,i+1,j,k+1)+fin(2,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hx(v)*hyi*( &
                zpi*( &
                cxdi*(-fin(4,i,j,k)  +fin(4,i,j+1,k))+ &
                cxd*( -fin(4,i+1,j,k)+fin(4,i+1,j+1,k))) &
                +zp*( &
                cxdi*(-fin(4,i,j,k+1)  +fin(4,i,j+1,k+1))+ &
                cxd*( -fin(4,i+1,j,k+1)+fin(4,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hz2*hxi(v)*hyi*( &
                czi*( &
                ( fin(6,i,j,k)  -fin(6,i,j+1,k)) &
                +(-fin(6,i+1,j,k)+fin(6,i+1,j+1,k))) &
                +cz*( &
                ( fin(6,i,j,k+1)  -fin(6,i,j+1,k+1)) &
                +(-fin(6,i+1,j,k+1)+fin(6,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hx(v)*hz2*hyi*( &
                czi*( &
                cxdi*(-fin(7,i,j,k)  +fin(7,i,j+1,k))+ &
                cxd*( -fin(7,i+1,j,k)+fin(7,i+1,j+1,k))) &
                +cz*( &
                cxdi*(-fin(7,i,j,k+1)  +fin(7,i,j+1,k+1))+ &
                cxd*( -fin(7,i+1,j,k+1)+fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(5).eq.1) then
        !                               ! fxzzz
        j=jj
        k=kk
        !
        !   ...and in y direction
        !
        yp=yparam
        ypi=1.0_fp-yp
        yp2=yp*yp
        ypi2=ypi*ypi
        !
        cy=yp*(yp2-1.0_fp)
        cyi=ypi*(ypi2-1.0_fp)
        hy2=hy*hy
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           !
           !   ...in x direction
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           xp2=xp*xp
           xpi2=xpi*xpi

           cxd=3.0_fp*xp2-1.0_fp
           cxdi=-3.0_fp*xpi2+1.0_fp
           !
           sum=hxi(v)*hzi*( &
                -( &
                -(ypi*fin(3,i,j,k)  +yp*fin(3,i,j+1,k)) &
                +(ypi*fin(3,i+1,j,k)+yp*fin(3,i+1,j+1,k))) &
                +( &
                -(ypi*fin(3,i,j,k+1)  +yp*fin(3,i,j+1,k+1)) &
                +(ypi*fin(3,i+1,j,k+1)+yp*fin(3,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hx(v)*hzi*( &
                -( &
                cxdi*(ypi*fin(5,i,j,k)  +yp*fin(5,i,j+1,k))+ &
                cxd*(ypi*fin(5,i+1,j,k)+yp*fin(5,i+1,j+1,k))) &
                +( &
                cxdi*(ypi*fin(5,i,j,k+1)  +yp*fin(5,i,j+1,k+1))+ &
                cxd*(ypi*fin(5,i+1,j,k+1)+yp*fin(5,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hy2*hxi(v)*hzi*( &
                -( &
                -(cyi*fin(6,i,j,k)  +cy*fin(6,i,j+1,k)) &
                +(cyi*fin(6,i+1,j,k)+cy*fin(6,i+1,j+1,k))) &
                +( &
                -(cyi*fin(6,i,j,k+1)  +cy*fin(6,i,j+1,k+1)) &
                +(cyi*fin(6,i+1,j,k+1)+cy*fin(6,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hx(v)*hy2*hzi*( &
                -( &
                cxdi*(cyi*fin(7,i,j,k)  +cy*fin(7,i,j+1,k))+ &
                cxd*(cyi*fin(7,i+1,j,k)+cy*fin(7,i+1,j+1,k))) &
                +( &
                cxdi*(cyi*fin(7,i,j,k+1)  +cy*fin(7,i,j+1,k+1))+ &
                cxd*(cyi*fin(7,i+1,j,k+1)+cy*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(6).eq.1) then
        !                               ! fyyyz
        j=jj
        k=kk
        !
        !   ...and in z direction
        !
        zp=zparam
        zpi=1.0_fp-zp
        zp2=zp*zp
        zpi2=zpi*zpi

        czd=3.0_fp*zp2-1.0_fp
        czdi=-3.0_fp*zpi2+1.0_fp
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           !
           !   ...in x direction
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           xp2=xp*xp
           xpi2=xpi*xpi
           !
           cx=xp*(xp2-1.0_fp)
           cxi=xpi*(xpi2-1.0_fp)
           hx2=hx(v)*hx(v)
           !
           sum=hyi*hzi*( &
                -( &
                xpi*(-fin(2,i,j,k)  +fin(2,i,j+1,k))+ &
                xp*( -fin(2,i+1,j,k)+fin(2,i+1,j+1,k))) &
                +( &
                xpi*(-fin(2,i,j,k+1)  +fin(2,i,j+1,k+1))+ &
                xp*( -fin(2,i+1,j,k+1)+fin(2,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hx2*hyi*hzi*( &
                -( &
                cxi*(-fin(4,i,j,k)  +fin(4,i,j+1,k))+ &
                cx*( -fin(4,i+1,j,k)+fin(4,i+1,j+1,k))) &
                +( &
                cxi*(-fin(4,i,j,k+1)  +fin(4,i,j+1,k+1))+ &
                cx*( -fin(4,i+1,j,k+1)+fin(4,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hz*hyi*( &
                czdi*( &
                xpi*(-fin(6,i,j,k)  +fin(6,i,j+1,k))+ &
                xp*( -fin(6,i+1,j,k)+fin(6,i+1,j+1,k))) &
                +czd*( &
                xpi*(-fin(6,i,j,k+1)  +fin(6,i,j+1,k+1))+ &
                xp*( -fin(6,i+1,j,k+1)+fin(6,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hx2*hz*hyi*( &
                czdi*( &
                cxi*(-fin(7,i,j,k)  +fin(7,i,j+1,k))+ &
                cx*( -fin(7,i+1,j,k)+fin(7,i+1,j+1,k))) &
                +czd*( &
                cxi*(-fin(7,i,j,k+1)  +fin(7,i,j+1,k+1))+ &
                cx*( -fin(7,i+1,j,k+1)+fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(7).eq.1) then
        !                               ! fyzzz
        j=jj
        k=kk
        !
        !   ...and in y direction
        !
        yp=yparam
        ypi=1.0_fp-yp
        yp2=yp*yp
        ypi2=ypi*ypi

        cyd=3.0_fp*yp2-1.0_fp
        cydi=-3.0_fp*ypi2+1.0_fp
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           !
           !   ...in x direction
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           xp2=xp*xp
           xpi2=xpi*xpi
           !
           cx=xp*(xp2-1.0_fp)
           cxi=xpi*(xpi2-1.0_fp)
           hx2=hx(v)*hx(v)
           !
           sum=hyi*hzi*( &
                -( &
                xpi*(-fin(3,i,j,k)  +fin(3,i,j+1,k))+ &
                xp*( -fin(3,i+1,j,k)+fin(3,i+1,j+1,k))) &
                +( &
                xpi*(-fin(3,i,j,k+1)  +fin(3,i,j+1,k+1))+ &
                xp*( -fin(3,i+1,j,k+1)+fin(3,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hx2*hyi*hzi*( &
                -( &
                cxi*(-fin(5,i,j,k)  +fin(5,i,j+1,k))+ &
                cx*( -fin(5,i+1,j,k)+fin(5,i+1,j+1,k))) &
                +( &
                cxi*(-fin(5,i,j,k+1)  +fin(5,i,j+1,k+1))+ &
                cx*( -fin(5,i+1,j,k+1)+fin(5,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hy*hzi*( &
                -( &
                xpi*(cydi*fin(6,i,j,k) +cyd*fin(6,i,j+1,k))+ &
                xp*(cydi*fin(6,i+1,j,k)+cyd*fin(6,i+1,j+1,k))) &
                +( &
                xpi*(cydi*fin(6,i,j,k+1) +cyd*fin(6,i,j+1,k+1))+ &
                xp*(cydi*fin(6,i+1,j,k+1)+cyd*fin(6,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hx2*hy*hzi*( &
                -( &
                cxi*(cydi*fin(7,i,j,k) +cyd*fin(7,i,j+1,k))+ &
                cx*(cydi*fin(7,i+1,j,k)+cyd*fin(7,i+1,j+1,k))) &
                +( &
                cxi*(cydi*fin(7,i,j,k+1) +cyd*fin(7,i,j+1,k+1))+ &
                cx*(cydi*fin(7,i+1,j,k+1)+cyd*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     !----------------------------------
     !  5th derivatives (.le.2 in each coordinate)
     !
  else if(ict(1).eq.5) then
     if(ict(2).eq.1) then
        !                               ! fxxyyz
        j=jj
        k=kk
        !
        !   ...and in y direction
        !
        yp=yparam
        ypi=1.0_fp-yp
        !
        !   ...and in z direction
        !
        zp=zparam
        zpi=1.0_fp-zp
        zp2=zp*zp
        zpi2=zpi*zpi

        czd=3.0_fp*zp2-1.0_fp
        czdi=-3.0_fp*zpi2+1.0_fp
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           !
           !   ...in x direction
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           !
           sum=hzi*( &
                -( &
                xpi*(ypi*fin(4,i,j,k)  +yp*fin(4,i,j+1,k))+ &
                xp*( ypi*fin(4,i+1,j,k)+yp*fin(4,i+1,j+1,k))) &
                +( &
                xpi*(ypi*fin(4,i,j,k+1) +yp*fin(4,i,j+1,k+1))+ &
                xp*(ypi*fin(4,i+1,j,k+1)+yp*fin(4,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hz*( &
                czdi*( &
                xpi*(ypi*fin(7,i,j,k) +yp*fin(7,i,j+1,k))+ &
                xp*(ypi*fin(7,i+1,j,k)+yp*fin(7,i+1,j+1,k))) &
                +czd*( &
                xpi*(ypi*fin(7,i,j,k+1) +yp*fin(7,i,j+1,k+1))+ &
                xp*(ypi*fin(7,i+1,j,k+1)+yp*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(3).eq.1) then
        !                               ! fxxyzz
        j=jj
        k=kk
        !
        !   ...and in y direction
        !
        yp=yparam
        ypi=1.0_fp-yp
        yp2=yp*yp
        ypi2=ypi*ypi

        cyd=3.0_fp*yp2-1.0_fp
        cydi=-3.0_fp*ypi2+1.0_fp
        !
        !   ...and in z direction
        !
        zp=zparam
        zpi=1.0_fp-zp
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           !
           !   ...in x direction
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           !
           sum=hyi*( &
                zpi*( &
                xpi*(-fin(5,i,j,k)  +fin(5,i,j+1,k))+ &
                xp*( -fin(5,i+1,j,k)+fin(5,i+1,j+1,k))) &
                +zp*( &
                xpi*(-fin(5,i,j,k+1)  +fin(5,i,j+1,k+1))+ &
                xp*( -fin(5,i+1,j,k+1)+fin(5,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hy*( &
                zpi*( &
                xpi*(cydi*fin(7,i,j,k)  +cyd*fin(7,i,j+1,k))+ &
                xp*( cydi*fin(7,i+1,j,k)+cyd*fin(7,i+1,j+1,k))) &
                +zp*( &
                xpi*(cydi*fin(7,i,j,k+1)  +cyd*fin(7,i,j+1,k+1))+ &
                xp*( cydi*fin(7,i+1,j,k+1)+cyd*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(4).eq.1) then
        !                               ! fxyyzz
        j=jj
        k=kk
        !
        !   ...and in y direction
        !
        yp=yparam
        ypi=1.0_fp-yp
        !
        !   ...and in z direction
        !
        zp=zparam
        zpi=1.0_fp-zp
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           !
           !   ...in x direction
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           xp2=xp*xp
           xpi2=xpi*xpi

           cxd=3.0_fp*xp2-1.0_fp
           cxdi=-3.0_fp*xpi2+1.0_fp
           !
           sum=hxi(v)*( &
                zpi*( &
                -(ypi*fin(6,i,j,k)  +yp*fin(6,i,j+1,k)) &
                +(ypi*fin(6,i+1,j,k)+yp*fin(6,i+1,j+1,k))) &
                +zp*( &
                -(ypi*fin(6,i,j,k+1)  +yp*fin(6,i,j+1,k+1)) &
                +(ypi*fin(6,i+1,j,k+1)+yp*fin(6,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hx(v)*( &
                zpi*( &
                cxdi*(ypi*fin(7,i,j,k)  +yp*fin(7,i,j+1,k))+ &
                cxd*(ypi*fin(7,i+1,j,k) +yp*fin(7,i+1,j+1,k))) &
                +zp*( &
                cxdi*(ypi*fin(7,i,j,k+1)  +yp*fin(7,i,j+1,k+1))+ &
                cxd*( ypi*fin(7,i+1,j,k+1)+yp*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     !----------------------------------
     !  5th derivatives (3 in a coordinate)
     !
  else if(ict(1).eq.-5) then
     if(ict(2).eq.1) then
        !                               ! fxxxyy
        j=jj
        k=kk
        !
        !   ...and in y direction
        !
        yp=yparam
        ypi=1.0_fp-yp
        !
        !   ...and in z direction
        !
        zp=zparam
        zpi=1.0_fp-zp
        zp2=zp*zp
        zpi2=zpi*zpi
        !
        cz=zp*(zp2-1.0_fp)
        czi=zpi*(zpi2-1.0_fp)
        hz2=hz*hz
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           !
           sum=hxi(v)*( &
                zpi*( &
                -(ypi*fin(4,i,j,k)  +yp*fin(4,i,j+1,k)) &
                +( ypi*fin(4,i+1,j,k)+yp*fin(4,i+1,j+1,k))) &
                +zp*( &
                -(ypi*fin(4,i,j,k+1) +yp*fin(4,i,j+1,k+1)) &
                +(ypi*fin(4,i+1,j,k+1)+yp*fin(4,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hz2*hxi(v)*( &
                czi*( &
                -(ypi*fin(7,i,j,k) +yp*fin(7,i,j+1,k)) &
                +(ypi*fin(7,i+1,j,k)+yp*fin(7,i+1,j+1,k))) &
                +cz*( &
                -(ypi*fin(7,i,j,k+1) +yp*fin(7,i,j+1,k+1)) &
                +(ypi*fin(7,i+1,j,k+1)+yp*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(3).eq.1) then
        !                               ! fxxxyz
        j=jj
        k=kk
        !
        !   ...and in y direction
        !
        yp=yparam
        ypi=1.0_fp-yp
        yp2=yp*yp
        ypi2=ypi*ypi

        cyd=3.0_fp*yp2-1.0_fp
        cydi=-3.0_fp*ypi2+1.0_fp
        !
        !   ...and in z direction
        !
        zp=zparam
        zpi=1.0_fp-zp
        zp2=zp*zp
        zpi2=zpi*zpi

        czd=3.0_fp*zp2-1.0_fp
        czdi=-3.0_fp*zpi2+1.0_fp
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           !
           sum=hyi*hzi*hxi(v)*( &
                -( &
                -(-fin(1,i,j,k)  +fin(1,i,j+1,k)) &
                +( -fin(1,i+1,j,k)+fin(1,i+1,j+1,k))) &
                +( &
                -(-fin(1,i,j,k+1)  +fin(1,i,j+1,k+1)) &
                +( -fin(1,i+1,j,k+1)+fin(1,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hy*hzi*hxi(v)*( &
                -( &
                -(cydi*fin(4,i,j,k)  +cyd*fin(4,i,j+1,k)) &
                +(cydi*fin(4,i+1,j,k)+cyd*fin(4,i+1,j+1,k))) &
                +( &
                -(cydi*fin(4,i,j,k+1)  +cyd*fin(4,i,j+1,k+1)) &
                +(cydi*fin(4,i+1,j,k+1)+cyd*fin(4,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hz*hyi*hxi(v)*( &
                czdi*( &
                -(-fin(5,i,j,k)  +fin(5,i,j+1,k)) &
                +( -fin(5,i+1,j,k)+fin(5,i+1,j+1,k))) &
                +czd*( &
                -(-fin(5,i,j,k+1)  +fin(5,i,j+1,k+1)) &
                +( -fin(5,i+1,j,k+1)+fin(5,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hy*hz*hxi(v)*( &
                czdi*( &
                -(cydi*fin(7,i,j,k)  +cyd*fin(7,i,j+1,k)) &
                +(cydi*fin(7,i+1,j,k)+cyd*fin(7,i+1,j+1,k))) &
                +czd*( &
                -(cydi*fin(7,i,j,k+1)  +cyd*fin(7,i,j+1,k+1)) &
                +(cydi*fin(7,i+1,j,k+1)+cyd*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(4).eq.1) then
        !                               ! fxxxzz
        j=jj
        k=kk
        !
        !   ...and in y direction
        !
        yp=yparam
        ypi=1.0_fp-yp
        yp2=yp*yp
        ypi2=ypi*ypi
        !
        cy=yp*(yp2-1.0_fp)
        cyi=ypi*(ypi2-1.0_fp)
        hy2=hy*hy
        !
        !   ...and in z direction
        !
        zp=zparam
        zpi=1.0_fp-zp
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           !
           sum=hxi(v)*( &
                zpi*( &
                -(ypi*fin(5,i,j,k)  +yp*fin(5,i,j+1,k)) &
                +(ypi*fin(5,i+1,j,k)+yp*fin(5,i+1,j+1,k))) &
                +zp*( &
                -(ypi*fin(5,i,j,k+1)  +yp*fin(5,i,j+1,k+1)) &
                +(ypi*fin(5,i+1,j,k+1)+yp*fin(5,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hy2*hxi(v)*( &
                zpi*( &
                -(cyi*fin(7,i,j,k)  +cy*fin(7,i,j+1,k)) &
                +(cyi*fin(7,i+1,j,k)+cy*fin(7,i+1,j+1,k))) &
                +zp*( &
                -(cyi*fin(7,i,j,k+1)  +cy*fin(7,i,j+1,k+1)) &
                +(cyi*fin(7,i+1,j,k+1)+cy*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(5).eq.1) then
        !                               ! fxxyyy
        j=jj
        k=kk
        !
        !   ...and in z direction
        !
        zp=zparam
        zpi=1.0_fp-zp
        zp2=zp*zp
        zpi2=zpi*zpi
        !
        cz=zp*(zp2-1.0_fp)
        czi=zpi*(zpi2-1.0_fp)
        hz2=hz*hz
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           !
           !   ...in x direction
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           !
           sum=hyi*( &
                zpi*( &
                xpi*(-fin(4,i,j,k)  +fin(4,i,j+1,k))+ &
                xp*( -fin(4,i+1,j,k)+fin(4,i+1,j+1,k))) &
                +zp*( &
                xpi*(-fin(4,i,j,k+1)  +fin(4,i,j+1,k+1))+ &
                xp*( -fin(4,i+1,j,k+1)+fin(4,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hz2*hyi*( &
                czi*( &
                xpi*(-fin(7,i,j,k)  +fin(7,i,j+1,k))+ &
                xp*( -fin(7,i+1,j,k)+fin(7,i+1,j+1,k))) &
                +cz*( &
                xpi*(-fin(7,i,j,k+1)  +fin(7,i,j+1,k+1))+ &
                xp*( -fin(7,i+1,j,k+1)+fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(6).eq.1) then
        !                               ! fxxzzz
        j=jj
        k=kk
        !
        !   ...and in y direction
        !
        yp=yparam
        ypi=1.0_fp-yp
        yp2=yp*yp
        ypi2=ypi*ypi
        !
        cy=yp*(yp2-1.0_fp)
        cyi=ypi*(ypi2-1.0_fp)
        hy2=hy*hy
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           !
           !   ...in x direction
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           !
           sum=hzi*( &
                -( &
                xpi*(ypi*fin(5,i,j,k)  +yp*fin(5,i,j+1,k))+ &
                xp*(ypi*fin(5,i+1,j,k)+yp*fin(5,i+1,j+1,k))) &
                +( &
                xpi*(ypi*fin(5,i,j,k+1)  +yp*fin(5,i,j+1,k+1))+ &
                xp*(ypi*fin(5,i+1,j,k+1)+yp*fin(5,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hy2*hzi*( &
                -( &
                xpi*(cyi*fin(7,i,j,k)  +cy*fin(7,i,j+1,k))+ &
                xp*(cyi*fin(7,i+1,j,k)+cy*fin(7,i+1,j+1,k))) &
                +( &
                xpi*(cyi*fin(7,i,j,k+1)  +cy*fin(7,i,j+1,k+1))+ &
                xp*(cyi*fin(7,i+1,j,k+1)+cy*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(7).eq.1) then
        !                               ! fxyyyz
        j=jj
        k=kk
        !
        !   ...and in z direction
        !
        zp=zparam
        zpi=1.0_fp-zp
        zp2=zp*zp
        zpi2=zpi*zpi

        czd=3.0_fp*zp2-1.0_fp
        czdi=-3.0_fp*zpi2+1.0_fp
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           !
           !   ...in x direction
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           xp2=xp*xp
           xpi2=xpi*xpi

           cxd=3.0_fp*xp2-1.0_fp
           cxdi=-3.0_fp*xpi2+1.0_fp
           !
           sum=hxi(v)*hzi*hyi*( &
                -( &
                ( fin(2,i,j,k)  -fin(2,i,j+1,k)) &
                +(-fin(2,i+1,j,k)+fin(2,i+1,j+1,k))) &
                +( &
                ( fin(2,i,j,k+1)  -fin(2,i,j+1,k+1)) &
                +(-fin(2,i+1,j,k+1)+fin(2,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hx(v)*hzi*hyi*( &
                -( &
                cxdi*(-fin(4,i,j,k)  +fin(4,i,j+1,k))+ &
                cxd*(-fin(4,i+1,j,k) +fin(4,i+1,j+1,k))) &
                +( &
                cxdi*(-fin(4,i,j,k+1)  +fin(4,i,j+1,k+1))+ &
                cxd*(-fin(4,i+1,j,k+1) +fin(4,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hz*hxi(v)*hyi*( &
                czdi*( &
                ( fin(6,i,j,k)  -fin(6,i,j+1,k)) &
                +(-fin(6,i+1,j,k)+fin(6,i+1,j+1,k))) &
                +czd*( &
                ( fin(6,i,j,k+1)  -fin(6,i,j+1,k+1)) &
                +(-fin(6,i+1,j,k+1)+fin(6,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hx(v)*hz*hyi*( &
                czdi*( &
                cxdi*(-fin(7,i,j,k)  +fin(7,i,j+1,k))+ &
                cxd*(-fin(7,i+1,j,k) +fin(7,i+1,j+1,k))) &
                +czd*( &
                cxdi*(-fin(7,i,j,k+1)  +fin(7,i,j+1,k+1))+ &
                cxd*(-fin(7,i+1,j,k+1) +fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(8).eq.1) then
        !                               ! fxyzzz
        j=jj
        k=kk
        !
        !   ...and in y direction
        !
        yp=yparam
        ypi=1.0_fp-yp
        yp2=yp*yp
        ypi2=ypi*ypi

        cyd=3.0_fp*yp2-1.0_fp
        cydi=-3.0_fp*ypi2+1.0_fp
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           !
           !   ...in x direction
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           xp2=xp*xp
           xpi2=xpi*xpi

           cxd=3.0_fp*xp2-1.0_fp
           cxdi=-3.0_fp*xpi2+1.0_fp
           !
           sum=hyi*hxi(v)*hzi*( &
                -( &
                ( +fin(3,i,j,k)  -fin(3,i,j+1,k)) &
                +( -fin(3,i+1,j,k)+fin(3,i+1,j+1,k))) &
                +( &
                ( +fin(3,i,j,k+1)  -fin(3,i,j+1,k+1)) &
                +( -fin(3,i+1,j,k+1)+fin(3,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hx(v)*hyi*hzi*( &
                -( &
                cxdi*(-fin(5,i,j,k)  +fin(5,i,j+1,k))+ &
                cxd*( -fin(5,i+1,j,k)+fin(5,i+1,j+1,k))) &
                +( &
                cxdi*(-fin(5,i,j,k+1)  +fin(5,i,j+1,k+1))+ &
                cxd*( -fin(5,i+1,j,k+1)+fin(5,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hy*hxi(v)*hzi*( &
                -( &
                -(cydi*fin(6,i,j,k) +cyd*fin(6,i,j+1,k)) &
                +(cydi*fin(6,i+1,j,k)+cyd*fin(6,i+1,j+1,k))) &
                +( &
                -(cydi*fin(6,i,j,k+1) +cyd*fin(6,i,j+1,k+1)) &
                +(cydi*fin(6,i+1,j,k+1)+cyd*fin(6,i+1,j+1,k+1))))
           !
           sum=sum+z36th*hx(v)*hy*hzi*( &
                -( &
                cxdi*(cydi*fin(7,i,j,k) +cyd*fin(7,i,j+1,k))+ &
                cxd*(cydi*fin(7,i+1,j,k)+cyd*fin(7,i+1,j+1,k))) &
                +( &
                cxdi*(cydi*fin(7,i,j,k+1) +cyd*fin(7,i,j+1,k+1))+ &
                cxd*(cydi*fin(7,i+1,j,k+1)+cyd*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(9).eq.1) then
        !                               ! fyyyzz
        j=jj
        k=kk
        !
        !   ...and in z direction
        !
        zp=zparam
        zpi=1.0_fp-zp
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           !
           !   ...in x direction
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           xp2=xp*xp
           xpi2=xpi*xpi
           !
           cx=xp*(xp2-1.0_fp)
           cxi=xpi*(xpi2-1.0_fp)
           hx2=hx(v)*hx(v)
           !
           sum=hyi*( &
                zpi*( &
                xpi*(-fin(6,i,j,k)  +fin(6,i,j+1,k))+ &
                xp*( -fin(6,i+1,j,k)+fin(6,i+1,j+1,k))) &
                +zp*( &
                xpi*(-fin(6,i,j,k+1)  +fin(6,i,j+1,k+1))+ &
                xp*( -fin(6,i+1,j,k+1)+fin(6,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hx2*hyi*( &
                zpi*( &
                cxi*(-fin(7,i,j,k)  +fin(7,i,j+1,k))+ &
                cx*( -fin(7,i+1,j,k)+fin(7,i+1,j+1,k))) &
                +zp*( &
                cxi*(-fin(7,i,j,k+1)  +fin(7,i,j+1,k+1))+ &
                cx*( -fin(7,i+1,j,k+1)+fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(10).eq.1) then
        !                               ! fyyzzz
        j=jj
        k=kk
        !
        !   ...and in y direction
        !
        yp=yparam
        ypi=1.0_fp-yp
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           !
           !   ...in x direction
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           xp2=xp*xp
           xpi2=xpi*xpi
           !
           cx=xp*(xp2-1.0_fp)
           cxi=xpi*(xpi2-1.0_fp)
           hx2=hx(v)*hx(v)
           !
           sum=hzi*( &
                -( &
                xpi*(ypi*fin(6,i,j,k)  +yp*fin(6,i,j+1,k))+ &
                xp*( ypi*fin(6,i+1,j,k)+yp*fin(6,i+1,j+1,k))) &
                +( &
                xpi*(ypi*fin(6,i,j,k+1) +yp*fin(6,i,j+1,k+1))+ &
                xp*(ypi*fin(6,i+1,j,k+1)+yp*fin(6,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hx2*hzi*( &
                -( &
                cxi*(ypi*fin(7,i,j,k)  +yp*fin(7,i,j+1,k))+ &
                cx*( ypi*fin(7,i+1,j,k)+yp*fin(7,i+1,j+1,k))) &
                +( &
                cxi*(ypi*fin(7,i,j,k+1)  +yp*fin(7,i,j+1,k+1))+ &
                cx*( ypi*fin(7,i+1,j,k+1)+yp*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     !----------------------------------
     !  6th derivatives (2 in each coordinate)
     !
  else if(ict(1).eq.6) then
     !                               ! fxxyyzz
     j=jj
     k=kk
     !
     !   ...and in y direction
     !
     yp=yparam
     ypi=1.0_fp-yp
     !
     !   ...and in z direction
     !
     zp=zparam
     zpi=1.0_fp-zp
     !
     iadr=iadr+1
     do v=1,ivec
        i=ii(v)
        !
        !   ...in x direction
        !
        xp=xparam(v)
        xpi=1.0_fp-xp
        !
        sum=( &
             zpi*( &
             xpi*(ypi*fin(7,i,j,k) +yp*fin(7,i,j+1,k))+ &
             xp*(ypi*fin(7,i+1,j,k)+yp*fin(7,i+1,j+1,k))) &
             +zp*( &
             xpi*(ypi*fin(7,i,j,k+1) +yp*fin(7,i,j+1,k+1))+ &
             xp*(ypi*fin(7,i+1,j,k+1)+yp*fin(7,i+1,j+1,k+1))))
        !
        fval(v,iadr)=sum
     end do
  end if
  !
  !----------------------------------
  !  6th derivatives (3 in a coordinate)
  !
  if(ict(1).eq.-6) then
     if(ict(2).eq.1) then
        !                               ! fxxxyyy
        j=jj
        k=kk
        !
        !   ...and in z direction
        !
        zp=zparam
        zpi=1.0_fp-zp
        zp2=zp*zp
        zpi2=zpi*zpi
        !
        cz=zp*(zp2-1.0_fp)
        czi=zpi*(zpi2-1.0_fp)
        hz2=hz*hz
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           !
           sum=hyi*hxi(v)*( &
                zpi*( &
                ( fin(4,i,j,k)  -fin(4,i,j+1,k)) &
                +(-fin(4,i+1,j,k)+fin(4,i+1,j+1,k))) &
                +zp*( &
                ( fin(4,i,j,k+1)  -fin(4,i,j+1,k+1)) &
                +(-fin(4,i+1,j,k+1)+fin(4,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hz2*hyi*hxi(v)*( &
                czi*( &
                ( fin(7,i,j,k)  -fin(7,i,j+1,k)) &
                +(-fin(7,i+1,j,k)+fin(7,i+1,j+1,k))) &
                +cz*( &
                ( fin(7,i,j,k+1)  -fin(7,i,j+1,k+1)) &
                +(-fin(7,i+1,j,k+1)+fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(3).eq.1) then
        !                               ! fxxxyyz
        j=jj
        k=kk
        !
        !   ...and in y direction
        !
        yp=yparam
        ypi=1.0_fp-yp
        !
        !   ...and in z direction
        !
        zp=zparam
        zpi=1.0_fp-zp
        zp2=zp*zp
        zpi2=zpi*zpi

        czd=3.0_fp*zp2-1.0_fp
        czdi=-3.0_fp*zpi2+1.0_fp
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           !
           sum=hxi(v)*hzi*( &
                -( &
                -(ypi*fin(4,i,j,k)  +yp*fin(4,i,j+1,k)) &
                +( ypi*fin(4,i+1,j,k)+yp*fin(4,i+1,j+1,k))) &
                +( &
                -(ypi*fin(4,i,j,k+1) +yp*fin(4,i,j+1,k+1)) &
                +(ypi*fin(4,i+1,j,k+1)+yp*fin(4,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hz*hxi(v)*( &
                czdi*( &
                -(ypi*fin(7,i,j,k) +yp*fin(7,i,j+1,k)) &
                +(ypi*fin(7,i+1,j,k)+yp*fin(7,i+1,j+1,k))) &
                +czd*( &
                -(ypi*fin(7,i,j,k+1) +yp*fin(7,i,j+1,k+1)) &
                +(ypi*fin(7,i+1,j,k+1)+yp*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(4).eq.1) then
        !                               ! fxxxyzz
        j=jj
        k=kk
        !
        !   ...and in y direction
        !
        yp=yparam
        ypi=1.0_fp-yp
        yp2=yp*yp
        ypi2=ypi*ypi

        cyd=3.0_fp*yp2-1.0_fp
        cydi=-3.0_fp*ypi2+1.0_fp
        !
        !   ...and in z direction
        !
        zp=zparam
        zpi=1.0_fp-zp
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           !
           sum=hxi(v)*hyi*( &
                zpi*( &
                ( fin(5,i,j,k)  -fin(5,i,j+1,k)) &
                +(-fin(5,i+1,j,k)+fin(5,i+1,j+1,k))) &
                +zp*( &
                ( fin(5,i,j,k+1)  -fin(5,i,j+1,k+1)) &
                +(-fin(5,i+1,j,k+1)+fin(5,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hy*hxi(v)*( &
                zpi*( &
                -(cydi*fin(7,i,j,k)  +cyd*fin(7,i,j+1,k)) &
                +(cydi*fin(7,i+1,j,k)+cyd*fin(7,i+1,j+1,k))) &
                +zp*( &
                -(cydi*fin(7,i,j,k+1)  +cyd*fin(7,i,j+1,k+1)) &
                +(cydi*fin(7,i+1,j,k+1)+cyd*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(5).eq.1) then
        !                               ! fxxxzzz
        j=jj
        k=kk
        !
        !   ...and in y direction
        !
        yp=yparam
        ypi=1.0_fp-yp
        yp2=yp*yp
        ypi2=ypi*ypi
        !
        cy=yp*(yp2-1.0_fp)
        cyi=ypi*(ypi2-1.0_fp)
        hy2=hy*hy
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           !
           sum=hxi(v)*hzi*( &
                -( &
                -(ypi*fin(5,i,j,k)  +yp*fin(5,i,j+1,k)) &
                +(ypi*fin(5,i+1,j,k)+yp*fin(5,i+1,j+1,k))) &
                +( &
                -(ypi*fin(5,i,j,k+1)  +yp*fin(5,i,j+1,k+1)) &
                +(ypi*fin(5,i+1,j,k+1)+yp*fin(5,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hy2*hxi(v)*hzi*( &
                -( &
                -(cyi*fin(7,i,j,k)  +cy*fin(7,i,j+1,k)) &
                +(cyi*fin(7,i+1,j,k)+cy*fin(7,i+1,j+1,k))) &
                +( &
                -(cyi*fin(7,i,j,k+1)  +cy*fin(7,i,j+1,k+1)) &
                +(cyi*fin(7,i+1,j,k+1)+cy*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(6).eq.1) then
        !                               ! fxxyyyz
        j=jj
        k=kk
        !
        !   ...and in z direction
        !
        zp=zparam
        zpi=1.0_fp-zp
        zp2=zp*zp
        zpi2=zpi*zpi

        czd=3.0_fp*zp2-1.0_fp
        czdi=-3.0_fp*zpi2+1.0_fp
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           !
           !   ...in x direction
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           !
           sum=hzi*hyi*( &
                -( &
                xpi*(-fin(4,i,j,k)  +fin(4,i,j+1,k))+ &
                xp*(-fin(4,i+1,j,k) +fin(4,i+1,j+1,k))) &
                +( &
                xpi*(-fin(4,i,j,k+1)  +fin(4,i,j+1,k+1))+ &
                xp*(-fin(4,i+1,j,k+1) +fin(4,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hz*hyi*( &
                czdi*( &
                xpi*(-fin(7,i,j,k)  +fin(7,i,j+1,k))+ &
                xp*(-fin(7,i+1,j,k) +fin(7,i+1,j+1,k))) &
                +czd*( &
                xpi*(-fin(7,i,j,k+1)  +fin(7,i,j+1,k+1))+ &
                xp*(-fin(7,i+1,j,k+1) +fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(7).eq.1) then
        !                               ! fxxyzzz
        j=jj
        k=kk
        !
        !   ...and in y direction
        !
        yp=yparam
        ypi=1.0_fp-yp
        yp2=yp*yp
        ypi2=ypi*ypi

        cyd=3.0_fp*yp2-1.0_fp
        cydi=-3.0_fp*ypi2+1.0_fp
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           !
           !   ...in x direction
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           !
           sum=hyi*hzi*( &
                -( &
                xpi*(-fin(5,i,j,k)  +fin(5,i,j+1,k))+ &
                xp*( -fin(5,i+1,j,k)+fin(5,i+1,j+1,k))) &
                +( &
                xpi*(-fin(5,i,j,k+1)  +fin(5,i,j+1,k+1))+ &
                xp*( -fin(5,i+1,j,k+1)+fin(5,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hy*hzi*( &
                -( &
                xpi*(cydi*fin(7,i,j,k)  +cyd*fin(7,i,j+1,k))+ &
                xp*( cydi*fin(7,i+1,j,k)+cyd*fin(7,i+1,j+1,k))) &
                +( &
                xpi*(cydi*fin(7,i,j,k+1)  +cyd*fin(7,i,j+1,k+1))+ &
                xp*( cydi*fin(7,i+1,j,k+1)+cyd*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(8).eq.1) then
        !                               ! fxyyyzz
        j=jj
        k=kk
        !
        !   ...and in z direction
        !
        zp=zparam
        zpi=1.0_fp-zp
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           !
           !   ...in x direction
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           xp2=xp*xp
           xpi2=xpi*xpi

           cxd=3.0_fp*xp2-1.0_fp
           cxdi=-3.0_fp*xpi2+1.0_fp
           !
           sum=hxi(v)*hyi*( &
                zpi*( &
                ( fin(6,i,j,k)  -fin(6,i,j+1,k)) &
                +(-fin(6,i+1,j,k)+fin(6,i+1,j+1,k))) &
                +zp*( &
                ( fin(6,i,j,k+1)  -fin(6,i,j+1,k+1)) &
                +(-fin(6,i+1,j,k+1)+fin(6,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hx(v)*hyi*( &
                zpi*( &
                cxdi*(-fin(7,i,j,k)  +fin(7,i,j+1,k))+ &
                cxd*(-fin(7,i+1,j,k) +fin(7,i+1,j+1,k))) &
                +zp*( &
                cxdi*(-fin(7,i,j,k+1)  +fin(7,i,j+1,k+1))+ &
                cxd*(-fin(7,i+1,j,k+1) +fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(9).eq.1) then
        !                               ! fxyyzzz
        j=jj
        k=kk
        !
        !   ...and in y direction
        !
        yp=yparam
        ypi=1.0_fp-yp
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           !
           !   ...in x direction
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           xp2=xp*xp
           xpi2=xpi*xpi

           cxd=3.0_fp*xp2-1.0_fp
           cxdi=-3.0_fp*xpi2+1.0_fp
           !
           sum=hxi(v)*hzi*( &
                -( &
                -(ypi*fin(6,i,j,k)  +yp*fin(6,i,j+1,k)) &
                +(ypi*fin(6,i+1,j,k)+yp*fin(6,i+1,j+1,k))) &
                +( &
                -(ypi*fin(6,i,j,k+1)  +yp*fin(6,i,j+1,k+1)) &
                +(ypi*fin(6,i+1,j,k+1)+yp*fin(6,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hx(v)*hzi*( &
                -( &
                cxdi*(ypi*fin(7,i,j,k)  +yp*fin(7,i,j+1,k))+ &
                cxd*(ypi*fin(7,i+1,j,k) +yp*fin(7,i+1,j+1,k))) &
                +( &
                cxdi*(ypi*fin(7,i,j,k+1)  +yp*fin(7,i,j+1,k+1))+ &
                cxd*( ypi*fin(7,i+1,j,k+1)+yp*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(10).eq.1) then
        !                               ! fyyyzzz
        j=jj
        k=kk
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           !
           !   ...in x direction
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           xp2=xp*xp
           xpi2=xpi*xpi
           !
           cx=xp*(xp2-1.0_fp)
           cxi=xpi*(xpi2-1.0_fp)
           hx2=hx(v)*hx(v)
           !
           sum=hyi*hzi*( &
                -( &
                xpi*(-fin(6,i,j,k)  +fin(6,i,j+1,k))+ &
                xp*( -fin(6,i+1,j,k)+fin(6,i+1,j+1,k))) &
                +( &
                xpi*(-fin(6,i,j,k+1)  +fin(6,i,j+1,k+1))+ &
                xp*( -fin(6,i+1,j,k+1)+fin(6,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hx2*hyi*hzi*( &
                -( &
                cxi*(-fin(7,i,j,k)  +fin(7,i,j+1,k))+ &
                cx*( -fin(7,i+1,j,k)+fin(7,i+1,j+1,k))) &
                +( &
                cxi*(-fin(7,i,j,k+1)  +fin(7,i,j+1,k+1))+ &
                cx*( -fin(7,i+1,j,k+1)+fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     !----------------------------------
     !  7th derivatives
     !
  else if(abs(ict(1)).eq.7) then
     if(ict(2).eq.1) then
        !                               ! fxxxyyyz
        j=jj
        k=kk
        !
        !   ...and in z direction
        !
        zp=zparam
        zpi=1.0_fp-zp
        zp2=zp*zp
        zpi2=zpi*zpi

        czd=3.0_fp*zp2-1.0_fp
        czdi=-3.0_fp*zpi2+1.0_fp
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           !
           sum=hyi*hxi(v)*hzi*( &
                -( &
                ( fin(4,i,j,k)  -fin(4,i,j+1,k)) &
                +(-fin(4,i+1,j,k)+fin(4,i+1,j+1,k))) &
                +( &
                ( fin(4,i,j,k+1)  -fin(4,i,j+1,k+1)) &
                +(-fin(4,i+1,j,k+1)+fin(4,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hz*hyi*hxi(v)*( &
                czdi*( &
                ( fin(7,i,j,k)  -fin(7,i,j+1,k)) &
                +(-fin(7,i+1,j,k)+fin(7,i+1,j+1,k))) &
                +czd*( &
                ( fin(7,i,j,k+1)  -fin(7,i,j+1,k+1)) &
                +(-fin(7,i+1,j,k+1)+fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(3).eq.1) then
        !                               ! fxxxyyzz
        j=jj
        k=kk
        !
        !   ...and in y direction
        !
        yp=yparam
        ypi=1.0_fp-yp
        !
        !   ...and in z direction
        !
        zp=zparam
        zpi=1.0_fp-zp
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           !
           sum=hxi(v)*( &
                zpi*( &
                -(ypi*fin(7,i,j,k) +yp*fin(7,i,j+1,k)) &
                +(ypi*fin(7,i+1,j,k)+yp*fin(7,i+1,j+1,k))) &
                +zp*( &
                -(ypi*fin(7,i,j,k+1) +yp*fin(7,i,j+1,k+1)) &
                +(ypi*fin(7,i+1,j,k+1)+yp*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(4).eq.1) then
        !                               ! fxxxyzzz
        j=jj
        k=kk
        !
        !   ...and in y direction
        !
        yp=yparam
        ypi=1.0_fp-yp
        yp2=yp*yp
        ypi2=ypi*ypi

        cyd=3.0_fp*yp2-1.0_fp
        cydi=-3.0_fp*ypi2+1.0_fp
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           !
           sum=hxi(v)*hyi*hzi*( &
                -( &
                ( fin(5,i,j,k)  -fin(5,i,j+1,k)) &
                +(-fin(5,i+1,j,k)+fin(5,i+1,j+1,k))) &
                +( &
                ( fin(5,i,j,k+1)  -fin(5,i,j+1,k+1)) &
                +(-fin(5,i+1,j,k+1)+fin(5,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hy*hxi(v)*hzi*( &
                -( &
                -(cydi*fin(7,i,j,k)  +cyd*fin(7,i,j+1,k)) &
                +(cydi*fin(7,i+1,j,k)+cyd*fin(7,i+1,j+1,k))) &
                +( &
                -(cydi*fin(7,i,j,k+1)  +cyd*fin(7,i,j+1,k+1)) &
                +(cydi*fin(7,i+1,j,k+1)+cyd*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(5).eq.1) then
        !                               ! fxxyyyzz
        j=jj
        k=kk
        !
        !   ...and in z direction
        !
        zp=zparam
        zpi=1.0_fp-zp
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           !
           !   ...in x direction
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           !
           sum=hyi*( &
                zpi*( &
                xpi*(-fin(7,i,j,k)  +fin(7,i,j+1,k))+ &
                xp*( -fin(7,i+1,j,k)+fin(7,i+1,j+1,k))) &
                +zp*( &
                xpi*(-fin(7,i,j,k+1)  +fin(7,i,j+1,k+1))+ &
                xp*( -fin(7,i+1,j,k+1)+fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(6).eq.1) then
        !                               ! fxxyyzzz
        j=jj
        k=kk
        !
        !   ...and in y direction
        !
        yp=yparam
        ypi=1.0_fp-yp
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           !
           !   ...in x direction
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           !
           sum=hzi*( &
                -( &
                xpi*(ypi*fin(7,i,j,k) +yp*fin(7,i,j+1,k))+ &
                xp*(ypi*fin(7,i+1,j,k)+yp*fin(7,i+1,j+1,k))) &
                +( &
                xpi*(ypi*fin(7,i,j,k+1) +yp*fin(7,i,j+1,k+1))+ &
                xp*(ypi*fin(7,i+1,j,k+1)+yp*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(7).eq.1) then
        !                               ! fxyyyzzz
        j=jj
        k=kk
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           !
           !   ...in x direction
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           xp2=xp*xp
           xpi2=xpi*xpi

           cxd=3.0_fp*xp2-1.0_fp
           cxdi=-3.0_fp*xpi2+1.0_fp
           !
           sum=hxi(v)*hyi*hzi*( &
                -( &
                ( fin(6,i,j,k)  -fin(6,i,j+1,k)) &
                +(-fin(6,i+1,j,k)+fin(6,i+1,j+1,k))) &
                +( &
                ( fin(6,i,j,k+1)  -fin(6,i,j+1,k+1)) &
                +(-fin(6,i+1,j,k+1)+fin(6,i+1,j+1,k+1))))
           !
           sum=sum+sixth*hx(v)*hyi*hzi*( &
                -( &
                cxdi*(-fin(7,i,j,k)  +fin(7,i,j+1,k))+ &
                cxd*(-fin(7,i+1,j,k) +fin(7,i+1,j+1,k))) &
                +( &
                cxdi*(-fin(7,i,j,k+1)  +fin(7,i,j+1,k+1))+ &
                cxd*(-fin(7,i+1,j,k+1) +fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     !----------------------------------
     !  8th derivatives
     !
  else if(abs(ict(1)).eq.8) then
     if(ict(2).eq.1) then
        !                               ! fxxxyyyzz
        j=jj
        k=kk
        !
        !   ...and in z direction
        !
        zp=zparam
        zpi=1.0_fp-zp
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           !
           sum=hyi*hxi(v)*( &
                zpi*( &
                ( fin(7,i,j,k)  -fin(7,i,j+1,k)) &
                +(-fin(7,i+1,j,k)+fin(7,i+1,j+1,k))) &
                +zp*( &
                ( fin(7,i,j,k+1)  -fin(7,i,j+1,k+1)) &
                +(-fin(7,i+1,j,k+1)+fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(3).eq.1) then
        !                               ! fxxxyyzzz
        j=jj
        k=kk
        !
        !   ...and in y direction
        !
        yp=yparam
        ypi=1.0_fp-yp
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           !
           sum=hxi(v)*hzi*( &
                -( &
                -(ypi*fin(7,i,j,k) +yp*fin(7,i,j+1,k)) &
                +(ypi*fin(7,i+1,j,k)+yp*fin(7,i+1,j+1,k))) &
                +( &
                -(ypi*fin(7,i,j,k+1) +yp*fin(7,i,j+1,k+1)) &
                +(ypi*fin(7,i+1,j,k+1)+yp*fin(7,i+1,j+1,k+1))))
           !
           fval(v,iadr)=sum
        end do
     end if
     !
     if(ict(4).eq.1) then
        !                               ! fxxyyyzzz
        j=jj
        k=kk
        !
        iadr=iadr+1
        do v=1,ivec
           i=ii(v)
           !
           !   ...in x direction
           !
           xp=xparam(v)
           xpi=1.0_fp-xp
           !
           sum=hyi*hzi*( &
                -( &
                xpi*(-fin(7,i,j,k)  +fin(7,i,j+1,k))+ &
                xp*( -fin(7,i+1,j,k)+fin(7,i+1,j+1,k))) &
                +( &
                xpi*(-fin(7,i,j,k+1)  +fin(7,i,j+1,k+1))+ &
                xp*( -fin(7,i+1,j,k+1)+fin(7,i+1,j+1,k+1))))
           fval(v,iadr)=sum
           !
        end do
     end if
     !
     !----------------------------------
     !  9th derivative
     !
  else if(abs(ict(1)).eq.9) then
     !                               ! fxxxyyyzzz
     j=jj
     k=kk
     !
     iadr=iadr+1
     do v=1,ivec
        i=ii(v)
        !
        sum=hyi*hxi(v)*hzi*( &
             -( &
             ( fin(7,i,j,k)  -fin(7,i,j+1,k)) &
             +(-fin(7,i+1,j,k)+fin(7,i+1,j+1,k))) &
             +( &
             ( fin(7,i,j,k+1)  -fin(7,i,j+1,k+1)) &
             +(-fin(7,i+1,j,k+1)+fin(7,i+1,j+1,k+1))))
        !
        fval(v,iadr)=sum
     end do
  end if
  !
  return
end subroutine fvtricubx
