subroutine herm3ev(xget,yget,zget,x,nx,y,ny,z,nz, &
     ilinx,iliny,ilinz, &
     f,inf2,inf3,ict,fval,ier)
  use psp_precision_mod, only: fp
  !
  !  evaluate a 3d cubic Hermite interpolant on a rectilinear
  !  grid -- this is C1 in all directions.
  !
  !  this subroutine calls two subroutines:
  !     herm3xyz  -- find cell containing (xget,yget,zget)
  !     herm3fcn -- evaluate interpolant function and (optionally) derivatives
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
  !  f(1,i,j,k) = df/dx @ x(i),y(j),z(k)
  !  f(2,i,j,k) = df/dy @ x(i),y(j),z(k)
  !  f(3,i,j,k) = df/dz @ x(i),y(j),z(k)
  !  f(4,i,j,k) = d2f/dxdy @ x(i),y(j),z(k)
  !  f(5,i,j,k) = d2f/dxdz @ x(i),y(j),z(k)
  !  f(6,i,j,k) = d2f/dydz @ x(i),y(j),z(k)
  !  f(7,i,j,k) = d3f/dxdydz @ x(i),y(j),z(k)
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
  integer, dimension(1) :: i,j,k                     ! cell indices
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
  call herm3fcn(ict,1,1, &
       fval,i,j,k,xparam,yparam,zparam, &
       hx,hxi,hy,hyi,hz,hzi, &
       f,inf2,inf3,nz)
  !
  return
end subroutine herm3ev
!---------------------------------------------------------------------
!  herm3xyz -- look up x-y-z zone
!
!  this is the "first part" of herm3ev, see comments, above.
!
subroutine herm3xyz(xget,yget,zget,x,nx,y,ny,z,nz, &
     ilinx,iliny,ilinz, &
     i,j,k,xparam,yparam,zparam, &
     hx,hxi,hy,hyi,hz,hzi,ier)
  use psp_precision_mod, only: fp
  !
  !  input of herm3xyz
  !  ================
  !
  !============
  implicit none
  integer nxm,nym,nzm,ii,jj,kk
  !============
  real(fp) :: zxget,zyget,zzget,zxtol,zytol,zztol
  !============
  integer nx,ny,nz                  ! coord. grid dimensions
  !
  real(fp) :: xget,yget,zget               ! target point
  real(fp) :: x(nx),y(ny),z(nz)            ! indep. coords. (ascending order)
  !
  integer ilinx                     ! =1:  x evenly spaced
  integer iliny                     ! =1:  y evenly spaced
  integer ilinz                     ! =1:  z evenly spaced
  !
  !  output of herm3xyz
  !  =================
  integer i,j,k                     ! index to cell containing target pt
  !          on exit:  1.le.i.le.nx-1   1.le.j.le.ny-1  1.le.k.le.nz-1
  !
  !  normalized position w/in (i,j) cell (btw 0 and 1):
  !
  real(fp) :: xparam                       ! (xget-x(i))/(x(i+1)-x(i))
  real(fp) :: yparam                       ! (yget-y(j))/(y(j+1)-y(j))
  real(fp) :: zparam                       ! (zget-z(k))/(z(k+1)-z(k))
  !
  !  grid spacing
  !
  real(fp) :: hx                           ! hx = x(i+1)-x(i)
  real(fp) :: hy                           ! hy = y(j+1)-y(j)
  real(fp) :: hz                           ! hz = z(k+1)-z(k)
  !
  !  inverse grid spacing:
  !
  real(fp) :: hxi                          ! 1/hx = 1/(x(i+1)-x(i))
  real(fp) :: hyi                          ! 1/hy = 1/(y(j+1)-y(j))
  real(fp) :: hzi                          ! 1/hz = 1/(z(k+1)-z(k))
  !
  integer ier                       ! return ier.ne.0 on error
  !
  !------------------------------------
  !
  ier=0
  !
  !  range check
  !
  zxget=xget
  zyget=yget
  zzget=zget
  if((xget.lt.x(1)).or.(xget.gt.x(nx))) then
     zxtol=4.0E-7_fp*max(abs(x(1)),abs(x(nx)))
     if((xget.lt.x(1)-zxtol).or.(xget.gt.x(nx)+zxtol)) then
        ier=1
        write(6,1001) xget,x(1),x(nx)
1001    format(' ?herm3ev:  xget=',1pe11.4,' out of range ', &
             1pe11.4,' to ',1pe11.4)
     else
        if((xget.lt.x(1)-0.5_fp*zxtol).or. &
             (xget.gt.x(nx)+0.5_fp*zxtol)) &
             write(6,1011) xget,x(1),x(nx)
1011    format(' %herm3ev:  xget=',1pe15.8,' beyond range ', &
             1pe15.8,' to ',1pe15.8,' (fixup applied)')
        if(xget.lt.x(1)) then
           zxget=x(1)
        else
           zxget=x(nx)
        end if
     end if
  end if
  if((yget.lt.y(1)).or.(yget.gt.y(ny))) then
     zytol=4.0E-7_fp*max(abs(y(1)),abs(y(ny)))
     if((yget.lt.y(1)-zytol).or.(yget.gt.y(ny)+zytol)) then
        ier=1
        write(6,1002) yget,y(1),y(ny)
1002    format(' ?herm3ev:  yget=',1pe11.4,' out of range ', &
             1pe11.4,' to ',1pe11.4)
     else
        if((yget.lt.y(1)-0.5_fp*zytol).or. &
             (yget.gt.y(ny)+0.5_fp*zytol)) &
             write(6,1012) yget,y(1),y(ny)
1012    format(' %herm3ev:  yget=',1pe15.8,' beyond range ', &
             1pe15.8,' to ',1pe15.8,' (fixup applied)')
        if(yget.lt.y(1)) then
           zyget=y(1)
        else
           zyget=y(ny)
        end if
     end if
  end if
  if((zget.lt.z(1)).or.(zget.gt.z(nz))) then
     zztol=4.0E-7_fp*max(abs(z(1)),abs(z(nz)))
     if((zget.lt.z(1)-zztol).or.(zget.gt.z(nz)+zztol)) then
        ier=1
        write(6,1003) zget,z(1),z(nz)
1003    format(' ?herm3ev:  zget=',1pe11.4,' out of range ', &
             1pe11.4,' to ',1pe11.4)
     else
        if((zget.lt.z(1)-0.5_fp*zztol).or. &
             (zget.gt.z(nz)+0.5_fp*zztol)) &
             write(6,1013) zget,z(1),z(nz)
1013    format(' %herm3ev:  zget=',1pe15.8,' beyond range ', &
             1pe15.8,' to ',1pe15.8,' (fixup applied)')
        if(zget.lt.z(1)) then
           zzget=z(1)
        else
           zzget=z(nz)
        end if
     end if
  end if
  if(ier.ne.0) return
  !
  !  now find interval in which target point lies..
  !
  nxm=nx-1
  nym=ny-1
  nzm=nz-1
  !
  if(ilinx.eq.1) then
     ii=1+nxm*(zxget-x(1))/(x(nx)-x(1))
     i=min(nxm, ii)
     if(zxget.lt.x(i)) then
        i=i-1
     else if(zxget.gt.x(i+1)) then
        i=i+1
     end if
  else
     if((1.le.i).and.(i.lt.nxm)) then
        if((x(i).le.zxget).and.(zxget.le.x(i+1))) then
           continue  ! already have the zone
        else
           call zonfind(x,nx,zxget,i)
        end if
     else
        i=nx/2
        call zonfind(x,nx,zxget,i)
     end if
  end if
  !
  if(iliny.eq.1) then
     jj=1+nym*(zyget-y(1))/(y(ny)-y(1))
     j=min(nym, jj)
     if(zyget.lt.y(j)) then
        j=j-1
     else if(zyget.gt.y(j+1)) then
        j=j+1
     end if
  else
     if((1.le.j).and.(j.lt.nym)) then
        if((y(j).le.zyget).and.(zyget.le.y(j+1))) then
           continue  ! already have the zone
        else
           call zonfind(y,ny,zyget,j)
        end if
     else
        j=ny/2
        call zonfind(y,ny,zyget,j)
     end if
  end if
  !
  if(ilinz.eq.1) then
     kk=1+nzm*(zzget-z(1))/(z(nz)-z(1))
     k=min(nzm, kk)
     if(zzget.lt.z(k)) then
        k=k-1
     else if(zzget.gt.z(k+1)) then
        k=k+1
     end if
  else
     if((1.le.k).and.(k.lt.nzm)) then
        if((z(k).le.zzget).and.(zzget.le.z(k+1))) then
           continue  ! already have the zone
        else
           call zonfind(z,nz,zzget,k)
        end if
     else
        k=nz/2
        call zonfind(z,nz,zzget,k)
     end if
  end if
  !
  hx=(x(i+1)-x(i))
  hy=(y(j+1)-y(j))
  hz=(z(k+1)-z(k))
  !
  hxi=1.0_fp/hx
  hyi=1.0_fp/hy
  hzi=1.0_fp/hz
  !
  xparam=(zxget-x(i))*hxi
  yparam=(zyget-y(j))*hyi
  zparam=(zzget-z(k))*hzi
  !
  return
end subroutine herm3xyz
!---------------------------------------------------------------------
!  evaluate C1 cubic Hermite function interpolation -- 3d fcn
!   --vectorized-- dmc 10 Feb 1999
!
subroutine herm3fcn(ict,ivec,ivecd, &
     fval,ii,jj,kk,xparam,yparam,zparam, &
     hx,hxi,hy,hyi,hz,hzi, &
     fin,inf2,inf3,nz)
  use psp_precision_mod, only: fp
  !
  !============
  implicit none
  integer inf3,nz,inf2,i,j,k,iadr
  !============
  real(fp) :: xp,xpi,xp2,xpi2,ax,axbar,bx,bxbar,yp,ypi,yp2,ypi2,ay
  real(fp) :: aybar,by,bybar,zp,zpi,zp2,zpi2,az,azbar,bz,bzbar,axp
  real(fp) :: axbarp,bxp,bxbarp,ayp,aybarp,byp,bybarp,azp,azbarp,bzp
  real(fp) :: bzbarp
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
  real(fp) :: fin(0:7,inf2,inf3,nz)        ! interpolant data (cf "herm3ev")
  !
  real(fp) :: fval(ivecd,*)                ! output returned
  !
  !  for detailed description of fin, ict and fval see subroutine herm3ev
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
  !  Hermite cubic basis functions
  !  -->for function value matching
  !     a(0)=0 a(1)=1        a'(0)=0 a'(1)=0
  !   abar(0)=1 abar(1)=0  abar'(0)=0 abar'(1)=0
  !
  !   a(x)=-2*x**3 + 3*x**2    = x*x*(-2.0*x+3.0)
  !   abar(x)=1-a(x)
  !   a'(x)=-abar'(x)          = 6.0*x*(1.0-x)
  !
  !  -->for derivative matching
  !     b(0)=0 b(1)=0          b'(0)=0 b'(1)=1
  !   bbar(0)=0 bbar(1)=0  bbar'(0)=1 bbar'(1)=0
  !
  !     b(x)=x**3-x**2         b'(x)=3*x**2-2*x
  !     bbar(x)=x**3-2*x**2+x  bbar'(x)=3*x**2-4*x+1
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
     xp2=xp*xp
     xpi2=xpi*xpi
     ax=xp2*(3.0_fp-2.0_fp*xp)
     axbar=1.0_fp-ax
     bx=-xp2*xpi
     bxbar=xpi2*xp
     !
     !   ...in y direction
     !
     yp=yparam(v)
     ypi=1.0_fp-yp
     yp2=yp*yp
     ypi2=ypi*ypi
     ay=yp2*(3.0_fp-2.0_fp*yp)
     aybar=1.0_fp-ay
     by=-yp2*ypi
     bybar=ypi2*yp
     !
     !   ...in z direction
     !
     zp=zparam(v)
     zpi=1.0_fp-zp
     zp2=zp*zp
     zpi2=zpi*zpi
     az=zp2*(3.0_fp-2.0_fp*zp)
     azbar=1.0_fp-az
     bz=-zp2*zpi
     bzbar=zpi2*zp
     !
     iadr=0
     !
     !  derivatives:
     !
     axp=6.0_fp*xp*xpi
     axbarp=-axp
     bxp=xp*(3.0_fp*xp-2.0_fp)
     bxbarp=xpi*(3.0_fp*xpi-2.0_fp)
     !
     ayp=6.0_fp*yp*ypi
     aybarp=-ayp
     byp=yp*(3.0_fp*yp-2.0_fp)
     bybarp=ypi*(3.0_fp*ypi-2.0_fp)
     !
     azp=6.0_fp*zp*zpi
     azbarp=-azp
     bzp=zp*(3.0_fp*zp-2.0_fp)
     bzbarp=zpi*(3.0_fp*zpi-2.0_fp)
     !
     !  get desired values:
     !
     if(ict(1).eq.1) then
        !
        !  function value:
        !
        iadr=iadr+1
        sum=azbar*( &
             axbar*(aybar*fin(0,i,j,k)  +ay*fin(0,i,j+1,k))+ &
             ax*(aybar*fin(0,i+1,j,k)+ay*fin(0,i+1,j+1,k))) &
             +  az*( &
             axbar*(aybar*fin(0,i,j,k+1)  +ay*fin(0,i,j+1,k+1))+ &
             ax*(aybar*fin(0,i+1,j,k+1)+ay*fin(0,i+1,j+1,k+1)))
        !
        sum=sum+hx(v)*( &
             azbar*( &
             bxbar*(aybar*fin(1,i,j,k)  +ay*fin(1,i,j+1,k))+ &
             bx*(aybar*fin(1,i+1,j,k)+ay*fin(1,i+1,j+1,k))) &
             + az*( &
             bxbar*(aybar*fin(1,i,j,k+1)  +ay*fin(1,i,j+1,k+1))+ &
             bx*(aybar*fin(1,i+1,j,k+1)+ay*fin(1,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hy(v)*( &
             azbar*( &
             axbar*(bybar*fin(2,i,j,k)  +by*fin(2,i,j+1,k))+ &
             ax*(bybar*fin(2,i+1,j,k)+by*fin(2,i+1,j+1,k))) &
             + az*( &
             axbar*(bybar*fin(2,i,j,k+1)  +by*fin(2,i,j+1,k+1))+ &
             ax*(bybar*fin(2,i+1,j,k+1)+by*fin(2,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hz(v)*( &
             bzbar*( &
             axbar*(aybar*fin(3,i,j,k)  +ay*fin(3,i,j+1,k))+ &
             ax*(aybar*fin(3,i+1,j,k)+ay*fin(3,i+1,j+1,k))) &
             + bz*( &
             axbar*(aybar*fin(3,i,j,k+1)  +ay*fin(3,i,j+1,k+1))+ &
             ax*(aybar*fin(3,i+1,j,k+1)+ay*fin(3,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hx(v)*hy(v)*( &
             azbar*( &
             bxbar*(bybar*fin(4,i,j,k)  +by*fin(4,i,j+1,k))+ &
             bx*(bybar*fin(4,i+1,j,k)+by*fin(4,i+1,j+1,k))) &
             + az*( &
             bxbar*(bybar*fin(4,i,j,k+1)  +by*fin(4,i,j+1,k+1))+ &
             bx*(bybar*fin(4,i+1,j,k+1)+by*fin(4,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hx(v)*hz(v)*( &
             bzbar*( &
             bxbar*(aybar*fin(5,i,j,k)  +ay*fin(5,i,j+1,k))+ &
             bx*(aybar*fin(5,i+1,j,k)+ay*fin(5,i+1,j+1,k))) &
             + bz*( &
             bxbar*(aybar*fin(5,i,j,k+1)  +ay*fin(5,i,j+1,k+1))+ &
             bx*(aybar*fin(5,i+1,j,k+1)+ay*fin(5,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hy(v)*hz(v)*( &
             bzbar*( &
             axbar*(bybar*fin(6,i,j,k)  +by*fin(6,i,j+1,k))+ &
             ax*(bybar*fin(6,i+1,j,k)+by*fin(6,i+1,j+1,k))) &
             + bz*( &
             axbar*(bybar*fin(6,i,j,k+1)  +by*fin(6,i,j+1,k+1))+ &
             ax*(bybar*fin(6,i+1,j,k+1)+by*fin(6,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hx(v)*hy(v)*hz(v)*( &
             bzbar*( &
             bxbar*(bybar*fin(7,i,j,k)  +by*fin(7,i,j+1,k))+ &
             bx*(bybar*fin(7,i+1,j,k)+by*fin(7,i+1,j+1,k))) &
             + bz*( &
             bxbar*(bybar*fin(7,i,j,k+1)  +by*fin(7,i,j+1,k+1))+ &
             bx*(bybar*fin(7,i+1,j,k+1)+by*fin(7,i+1,j+1,k+1))) &
             )
        !
        fval(v,iadr)=sum
     end if
     !
     if(ict(2).eq.1) then
        !
        !     df/dx:
        !
        iadr=iadr+1
        !
        sum=hxi(v)*( &
             azbar*( &
             axbarp*(aybar*fin(0,i,j,k)  +ay*fin(0,i,j+1,k))+ &
             axp*(aybar*fin(0,i+1,j,k)+ay*fin(0,i+1,j+1,k))) &
             +  az*( &
             axbarp*(aybar*fin(0,i,j,k+1)  +ay*fin(0,i,j+1,k+1))+ &
             axp*(aybar*fin(0,i+1,j,k+1)+ay*fin(0,i+1,j+1,k+1))) &
             )
        !
        sum=sum+( &
             azbar*( &
             bxbarp*(aybar*fin(1,i,j,k)  +ay*fin(1,i,j+1,k))+ &
             bxp*(aybar*fin(1,i+1,j,k)+ay*fin(1,i+1,j+1,k))) &
             + az*( &
             bxbarp*(aybar*fin(1,i,j,k+1)  +ay*fin(1,i,j+1,k+1))+ &
             bxp*(aybar*fin(1,i+1,j,k+1)+ay*fin(1,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hxi(v)*hy(v)*( &
             azbar*( &
             axbarp*(bybar*fin(2,i,j,k)  +by*fin(2,i,j+1,k))+ &
             axp*(bybar*fin(2,i+1,j,k)+by*fin(2,i+1,j+1,k))) &
             + az*( &
             axbarp*(bybar*fin(2,i,j,k+1)  +by*fin(2,i,j+1,k+1))+ &
             axp*(bybar*fin(2,i+1,j,k+1)+by*fin(2,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hxi(v)*hz(v)*( &
             bzbar*( &
             axbarp*(aybar*fin(3,i,j,k)  +ay*fin(3,i,j+1,k))+ &
             axp*(aybar*fin(3,i+1,j,k)+ay*fin(3,i+1,j+1,k))) &
             + bz*( &
             axbarp*(aybar*fin(3,i,j,k+1)  +ay*fin(3,i,j+1,k+1))+ &
             axp*(aybar*fin(3,i+1,j,k+1)+ay*fin(3,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hy(v)*( &
             azbar*( &
             bxbarp*(bybar*fin(4,i,j,k)  +by*fin(4,i,j+1,k))+ &
             bxp*(bybar*fin(4,i+1,j,k)+by*fin(4,i+1,j+1,k))) &
             + az*( &
             bxbarp*(bybar*fin(4,i,j,k+1)  +by*fin(4,i,j+1,k+1))+ &
             bxp*(bybar*fin(4,i+1,j,k+1)+by*fin(4,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hz(v)*( &
             bzbar*( &
             bxbarp*(aybar*fin(5,i,j,k)  +ay*fin(5,i,j+1,k))+ &
             bxp*(aybar*fin(5,i+1,j,k)+ay*fin(5,i+1,j+1,k))) &
             + bz*( &
             bxbarp*(aybar*fin(5,i,j,k+1)  +ay*fin(5,i,j+1,k+1))+ &
             bxp*(aybar*fin(5,i+1,j,k+1)+ay*fin(5,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hxi(v)*hy(v)*hz(v)*( &
             bzbar*( &
             axbarp*(bybar*fin(6,i,j,k)  +by*fin(6,i,j+1,k))+ &
             axp*(bybar*fin(6,i+1,j,k)+by*fin(6,i+1,j+1,k))) &
             + bz*( &
             axbarp*(bybar*fin(6,i,j,k+1)  +by*fin(6,i,j+1,k+1))+ &
             axp*(bybar*fin(6,i+1,j,k+1)+by*fin(6,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hy(v)*hz(v)*( &
             bzbar*( &
             bxbarp*(bybar*fin(7,i,j,k)  +by*fin(7,i,j+1,k))+ &
             bxp*(bybar*fin(7,i+1,j,k)+by*fin(7,i+1,j+1,k))) &
             + bz*( &
             bxbarp*(bybar*fin(7,i,j,k+1)  +by*fin(7,i,j+1,k+1))+ &
             bxp*(bybar*fin(7,i+1,j,k+1)+by*fin(7,i+1,j+1,k+1))) &
             )
        !
        fval(v,iadr)=sum
     end if
     !
     if(ict(3).eq.1) then
        !
        !  df/dy:
        !
        iadr=iadr+1
        !
        sum=hyi(v)*( &
             azbar*( &
             axbar*(aybarp*fin(0,i,j,k)  +ayp*fin(0,i,j+1,k))+ &
             ax*(aybarp*fin(0,i+1,j,k)+ayp*fin(0,i+1,j+1,k))) &
             +  az*( &
             axbar*(aybarp*fin(0,i,j,k+1)  +ayp*fin(0,i,j+1,k+1))+ &
             ax*(aybarp*fin(0,i+1,j,k+1)+ayp*fin(0,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hyi(v)*hx(v)*( &
             azbar*( &
             bxbar*(aybarp*fin(1,i,j,k)  +ayp*fin(1,i,j+1,k))+ &
             bx*(aybarp*fin(1,i+1,j,k)+ayp*fin(1,i+1,j+1,k))) &
             + az*( &
             bxbar*(aybarp*fin(1,i,j,k+1)  +ayp*fin(1,i,j+1,k+1))+ &
             bx*(aybarp*fin(1,i+1,j,k+1)+ayp*fin(1,i+1,j+1,k+1))) &
             )
        !
        sum=sum+( &
             azbar*( &
             axbar*(bybarp*fin(2,i,j,k)  +byp*fin(2,i,j+1,k))+ &
             ax*(bybarp*fin(2,i+1,j,k)+byp*fin(2,i+1,j+1,k))) &
             + az*( &
             axbar*(bybarp*fin(2,i,j,k+1)  +byp*fin(2,i,j+1,k+1))+ &
             ax*(bybarp*fin(2,i+1,j,k+1)+byp*fin(2,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hyi(v)*hz(v)*( &
             bzbar*( &
             axbar*(aybarp*fin(3,i,j,k)  +ayp*fin(3,i,j+1,k))+ &
             ax*(aybarp*fin(3,i+1,j,k)+ayp*fin(3,i+1,j+1,k))) &
             + bz*( &
             axbar*(aybarp*fin(3,i,j,k+1)  +ayp*fin(3,i,j+1,k+1))+ &
             ax*(aybarp*fin(3,i+1,j,k+1)+ayp*fin(3,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hx(v)*( &
             azbar*( &
             bxbar*(bybarp*fin(4,i,j,k)  +byp*fin(4,i,j+1,k))+ &
             bx*(bybarp*fin(4,i+1,j,k)+byp*fin(4,i+1,j+1,k))) &
             + az*( &
             bxbar*(bybarp*fin(4,i,j,k+1)  +byp*fin(4,i,j+1,k+1))+ &
             bx*(bybarp*fin(4,i+1,j,k+1)+byp*fin(4,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hx(v)*hyi(v)*hz(v)*( &
             bzbar*( &
             bxbar*(aybarp*fin(5,i,j,k)  +ayp*fin(5,i,j+1,k))+ &
             bx*(aybarp*fin(5,i+1,j,k)+ayp*fin(5,i+1,j+1,k))) &
             + bz*( &
             bxbar*(aybarp*fin(5,i,j,k+1)  +ayp*fin(5,i,j+1,k+1))+ &
             bx*(aybarp*fin(5,i+1,j,k+1)+ayp*fin(5,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hz(v)*( &
             bzbar*( &
             axbar*(bybarp*fin(6,i,j,k)  +byp*fin(6,i,j+1,k))+ &
             ax*(bybarp*fin(6,i+1,j,k)+byp*fin(6,i+1,j+1,k))) &
             + bz*( &
             axbar*(bybarp*fin(6,i,j,k+1)  +byp*fin(6,i,j+1,k+1))+ &
             ax*(bybarp*fin(6,i+1,j,k+1)+byp*fin(6,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hx(v)*hz(v)*( &
             bzbar*( &
             bxbar*(bybarp*fin(7,i,j,k)  +byp*fin(7,i,j+1,k))+ &
             bx*(bybarp*fin(7,i+1,j,k)+byp*fin(7,i+1,j+1,k))) &
             + bz*( &
             bxbar*(bybarp*fin(7,i,j,k+1)  +byp*fin(7,i,j+1,k+1))+ &
             bx*(bybarp*fin(7,i+1,j,k+1)+byp*fin(7,i+1,j+1,k+1))) &
             )
        !
        fval(v,iadr)=sum
     end if
     !
     if(ict(4).eq.1) then
        !
        !  df/dz:
        !
        iadr=iadr+1
        !
        sum=hzi(v)*( &
             azbarp*( &
             axbar*(aybar*fin(0,i,j,k)  +ay*fin(0,i,j+1,k))+ &
             ax*(aybar*fin(0,i+1,j,k)+ay*fin(0,i+1,j+1,k))) &
             +  azp*( &
             axbar*(aybar*fin(0,i,j,k+1)  +ay*fin(0,i,j+1,k+1))+ &
             ax*(aybar*fin(0,i+1,j,k+1)+ay*fin(0,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hzi(v)*hx(v)*( &
             azbarp*( &
             bxbar*(aybar*fin(1,i,j,k)  +ay*fin(1,i,j+1,k))+ &
             bx*(aybar*fin(1,i+1,j,k)+ay*fin(1,i+1,j+1,k))) &
             + azp*( &
             bxbar*(aybar*fin(1,i,j,k+1)  +ay*fin(1,i,j+1,k+1))+ &
             bx*(aybar*fin(1,i+1,j,k+1)+ay*fin(1,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hzi(v)*hy(v)*( &
             azbarp*( &
             axbar*(bybar*fin(2,i,j,k)  +by*fin(2,i,j+1,k))+ &
             ax*(bybar*fin(2,i+1,j,k)+by*fin(2,i+1,j+1,k))) &
             + azp*( &
             axbar*(bybar*fin(2,i,j,k+1)  +by*fin(2,i,j+1,k+1))+ &
             ax*(bybar*fin(2,i+1,j,k+1)+by*fin(2,i+1,j+1,k+1))) &
             )
        !
        sum=sum+( &
             bzbarp*( &
             axbar*(aybar*fin(3,i,j,k)  +ay*fin(3,i,j+1,k))+ &
             ax*(aybar*fin(3,i+1,j,k)+ay*fin(3,i+1,j+1,k))) &
             + bzp*( &
             axbar*(aybar*fin(3,i,j,k+1)  +ay*fin(3,i,j+1,k+1))+ &
             ax*(aybar*fin(3,i+1,j,k+1)+ay*fin(3,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hzi(v)*hx(v)*hy(v)*( &
             azbarp*( &
             bxbar*(bybar*fin(4,i,j,k)  +by*fin(4,i,j+1,k))+ &
             bx*(bybar*fin(4,i+1,j,k)+by*fin(4,i+1,j+1,k))) &
             + azp*( &
             bxbar*(bybar*fin(4,i,j,k+1)  +by*fin(4,i,j+1,k+1))+ &
             bx*(bybar*fin(4,i+1,j,k+1)+by*fin(4,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hx(v)*( &
             bzbarp*( &
             bxbar*(aybar*fin(5,i,j,k)  +ay*fin(5,i,j+1,k))+ &
             bx*(aybar*fin(5,i+1,j,k)+ay*fin(5,i+1,j+1,k))) &
             + bzp*( &
             bxbar*(aybar*fin(5,i,j,k+1)  +ay*fin(5,i,j+1,k+1))+ &
             bx*(aybar*fin(5,i+1,j,k+1)+ay*fin(5,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hy(v)*( &
             bzbarp*( &
             axbar*(bybar*fin(6,i,j,k)  +by*fin(6,i,j+1,k))+ &
             ax*(bybar*fin(6,i+1,j,k)+by*fin(6,i+1,j+1,k))) &
             + bzp*( &
             axbar*(bybar*fin(6,i,j,k+1)  +by*fin(6,i,j+1,k+1))+ &
             ax*(bybar*fin(6,i+1,j,k+1)+by*fin(6,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hx(v)*hy(v)*( &
             bzbarp*( &
             bxbar*(bybar*fin(7,i,j,k)  +by*fin(7,i,j+1,k))+ &
             bx*(bybar*fin(7,i+1,j,k)+by*fin(7,i+1,j+1,k))) &
             + bzp*( &
             bxbar*(bybar*fin(7,i,j,k+1)  +by*fin(7,i,j+1,k+1))+ &
             bx*(bybar*fin(7,i+1,j,k+1)+by*fin(7,i+1,j+1,k+1))) &
             )
        !
        fval(v,iadr)=sum
     end if
     !
     if(ict(5).eq.1) then
        !
        !  d2f/dxdy:
        !
        iadr=iadr+1
        !
        sum=hxi(v)*hyi(v)*( &
             azbar*( &
             axbarp*(aybarp*fin(0,i,j,k)  +ayp*fin(0,i,j+1,k))+ &
             axp*(aybarp*fin(0,i+1,j,k)+ayp*fin(0,i+1,j+1,k))) &
             +  az*( &
             axbarp*(aybarp*fin(0,i,j,k+1)  +ayp*fin(0,i,j+1,k+1))+ &
             axp*(aybarp*fin(0,i+1,j,k+1)+ayp*fin(0,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hyi(v)*( &
             azbar*( &
             bxbarp*(aybarp*fin(1,i,j,k)  +ayp*fin(1,i,j+1,k))+ &
             bxp*(aybarp*fin(1,i+1,j,k)+ayp*fin(1,i+1,j+1,k))) &
             + az*( &
             bxbarp*(aybarp*fin(1,i,j,k+1)  +ayp*fin(1,i,j+1,k+1))+ &
             bxp*(aybarp*fin(1,i+1,j,k+1)+ayp*fin(1,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hxi(v)*( &
             azbar*( &
             axbarp*(bybarp*fin(2,i,j,k)  +byp*fin(2,i,j+1,k))+ &
             axp*(bybarp*fin(2,i+1,j,k)+byp*fin(2,i+1,j+1,k))) &
             + az*( &
             axbarp*(bybarp*fin(2,i,j,k+1)  +byp*fin(2,i,j+1,k+1))+ &
             axp*(bybarp*fin(2,i+1,j,k+1)+byp*fin(2,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hxi(v)*hyi(v)*hz(v)*( &
             bzbar*( &
             axbarp*(aybarp*fin(3,i,j,k)  +ayp*fin(3,i,j+1,k))+ &
             axp*(aybarp*fin(3,i+1,j,k)+ayp*fin(3,i+1,j+1,k))) &
             + bz*( &
             axbarp*(aybarp*fin(3,i,j,k+1)  +ayp*fin(3,i,j+1,k+1))+ &
             axp*(aybarp*fin(3,i+1,j,k+1)+ayp*fin(3,i+1,j+1,k+1))) &
             )
        !
        sum=sum+( &
             azbar*( &
             bxbarp*(bybarp*fin(4,i,j,k)  +byp*fin(4,i,j+1,k))+ &
             bxp*(bybarp*fin(4,i+1,j,k)+byp*fin(4,i+1,j+1,k))) &
             + az*( &
             bxbarp*(bybarp*fin(4,i,j,k+1)  +byp*fin(4,i,j+1,k+1))+ &
             bxp*(bybarp*fin(4,i+1,j,k+1)+byp*fin(4,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hyi(v)*hz(v)*( &
             bzbar*( &
             bxbarp*(aybarp*fin(5,i,j,k)  +ayp*fin(5,i,j+1,k))+ &
             bxp*(aybarp*fin(5,i+1,j,k)+ayp*fin(5,i+1,j+1,k))) &
             + bz*( &
             bxbarp*(aybarp*fin(5,i,j,k+1)  +ayp*fin(5,i,j+1,k+1))+ &
             bxp*(aybarp*fin(5,i+1,j,k+1)+ayp*fin(5,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hxi(v)*hz(v)*( &
             bzbar*( &
             axbarp*(bybarp*fin(6,i,j,k)  +byp*fin(6,i,j+1,k))+ &
             axp*(bybarp*fin(6,i+1,j,k)+byp*fin(6,i+1,j+1,k))) &
             + bz*( &
             axbarp*(bybarp*fin(6,i,j,k+1)  +byp*fin(6,i,j+1,k+1))+ &
             axp*(bybarp*fin(6,i+1,j,k+1)+byp*fin(6,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hz(v)*( &
             bzbar*( &
             bxbarp*(bybarp*fin(7,i,j,k)  +byp*fin(7,i,j+1,k))+ &
             bxp*(bybarp*fin(7,i+1,j,k)+byp*fin(7,i+1,j+1,k))) &
             + bz*( &
             bxbarp*(bybarp*fin(7,i,j,k+1)  +byp*fin(7,i,j+1,k+1))+ &
             bxp*(bybarp*fin(7,i+1,j,k+1)+byp*fin(7,i+1,j+1,k+1))) &
             )
        !
        fval(v,iadr)=sum
     end if
     !
     if(ict(6).eq.1) then
        !
        !  d2f/dxdz:
        !
        iadr=iadr+1
        !
        sum=hxi(v)*hzi(v)*( &
             azbarp*( &
             axbarp*(aybar*fin(0,i,j,k)  +ay*fin(0,i,j+1,k))+ &
             axp*(aybar*fin(0,i+1,j,k)+ay*fin(0,i+1,j+1,k))) &
             +  azp*( &
             axbarp*(aybar*fin(0,i,j,k+1)  +ay*fin(0,i,j+1,k+1))+ &
             axp*(aybar*fin(0,i+1,j,k+1)+ay*fin(0,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hzi(v)*( &
             azbarp*( &
             bxbarp*(aybar*fin(1,i,j,k)  +ay*fin(1,i,j+1,k))+ &
             bxp*(aybar*fin(1,i+1,j,k)+ay*fin(1,i+1,j+1,k))) &
             + azp*( &
             bxbarp*(aybar*fin(1,i,j,k+1)  +ay*fin(1,i,j+1,k+1))+ &
             bxp*(aybar*fin(1,i+1,j,k+1)+ay*fin(1,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hxi(v)*hy(v)*hzi(v)*( &
             azbarp*( &
             axbarp*(bybar*fin(2,i,j,k)  +by*fin(2,i,j+1,k))+ &
             axp*(bybar*fin(2,i+1,j,k)+by*fin(2,i+1,j+1,k))) &
             + azp*( &
             axbarp*(bybar*fin(2,i,j,k+1)  +by*fin(2,i,j+1,k+1))+ &
             axp*(bybar*fin(2,i+1,j,k+1)+by*fin(2,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hxi(v)*( &
             bzbarp*( &
             axbarp*(aybar*fin(3,i,j,k)  +ay*fin(3,i,j+1,k))+ &
             axp*(aybar*fin(3,i+1,j,k)+ay*fin(3,i+1,j+1,k))) &
             + bzp*( &
             axbarp*(aybar*fin(3,i,j,k+1)  +ay*fin(3,i,j+1,k+1))+ &
             axp*(aybar*fin(3,i+1,j,k+1)+ay*fin(3,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hy(v)*hzi(v)*( &
             azbarp*( &
             bxbarp*(bybar*fin(4,i,j,k)  +by*fin(4,i,j+1,k))+ &
             bxp*(bybar*fin(4,i+1,j,k)+by*fin(4,i+1,j+1,k))) &
             + azp*( &
             bxbarp*(bybar*fin(4,i,j,k+1)  +by*fin(4,i,j+1,k+1))+ &
             bxp*(bybar*fin(4,i+1,j,k+1)+by*fin(4,i+1,j+1,k+1))) &
             )
        !
        sum=sum+( &
             bzbarp*( &
             bxbarp*(aybar*fin(5,i,j,k)  +ay*fin(5,i,j+1,k))+ &
             bxp*(aybar*fin(5,i+1,j,k)+ay*fin(5,i+1,j+1,k))) &
             + bzp*( &
             bxbarp*(aybar*fin(5,i,j,k+1)  +ay*fin(5,i,j+1,k+1))+ &
             bxp*(aybar*fin(5,i+1,j,k+1)+ay*fin(5,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hxi(v)*hy(v)*( &
             bzbarp*( &
             axbarp*(bybar*fin(6,i,j,k)  +by*fin(6,i,j+1,k))+ &
             axp*(bybar*fin(6,i+1,j,k)+by*fin(6,i+1,j+1,k))) &
             + bzp*( &
             axbarp*(bybar*fin(6,i,j,k+1)  +by*fin(6,i,j+1,k+1))+ &
             axp*(bybar*fin(6,i+1,j,k+1)+by*fin(6,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hy(v)*( &
             bzbarp*( &
             bxbarp*(bybar*fin(7,i,j,k)  +by*fin(7,i,j+1,k))+ &
             bxp*(bybar*fin(7,i+1,j,k)+by*fin(7,i+1,j+1,k))) &
             + bzp*( &
             bxbarp*(bybar*fin(7,i,j,k+1)  +by*fin(7,i,j+1,k+1))+ &
             bxp*(bybar*fin(7,i+1,j,k+1)+by*fin(7,i+1,j+1,k+1))) &
             )
        !
        fval(v,iadr)=sum
     end if
     !
     if(ict(7).eq.1) then
        !
        !  d2f/dydz:
        !
        iadr=iadr+1
        !
        sum=hyi(v)*hzi(v)*( &
             azbarp*( &
             axbar*(aybarp*fin(0,i,j,k)  +ayp*fin(0,i,j+1,k))+ &
             ax*(aybarp*fin(0,i+1,j,k)+ayp*fin(0,i+1,j+1,k))) &
             +  azp*( &
             axbar*(aybarp*fin(0,i,j,k+1)  +ayp*fin(0,i,j+1,k+1))+ &
             ax*(aybarp*fin(0,i+1,j,k+1)+ayp*fin(0,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hyi(v)*hzi(v)*hx(v)*( &
             azbarp*( &
             bxbar*(aybarp*fin(1,i,j,k)  +ayp*fin(1,i,j+1,k))+ &
             bx*(aybarp*fin(1,i+1,j,k)+ayp*fin(1,i+1,j+1,k))) &
             + azp*( &
             bxbar*(aybarp*fin(1,i,j,k+1)  +ayp*fin(1,i,j+1,k+1))+ &
             bx*(aybarp*fin(1,i+1,j,k+1)+ayp*fin(1,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hzi(v)*( &
             azbarp*( &
             axbar*(bybarp*fin(2,i,j,k)  +byp*fin(2,i,j+1,k))+ &
             ax*(bybarp*fin(2,i+1,j,k)+byp*fin(2,i+1,j+1,k))) &
             + azp*( &
             axbar*(bybarp*fin(2,i,j,k+1)  +byp*fin(2,i,j+1,k+1))+ &
             ax*(bybarp*fin(2,i+1,j,k+1)+byp*fin(2,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hyi(v)*( &
             bzbarp*( &
             axbar*(aybarp*fin(3,i,j,k)  +ayp*fin(3,i,j+1,k))+ &
             ax*(aybarp*fin(3,i+1,j,k)+ayp*fin(3,i+1,j+1,k))) &
             + bzp*( &
             axbar*(aybarp*fin(3,i,j,k+1)  +ayp*fin(3,i,j+1,k+1))+ &
             ax*(aybarp*fin(3,i+1,j,k+1)+ayp*fin(3,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hx(v)*hzi(v)*( &
             azbarp*( &
             bxbar*(bybarp*fin(4,i,j,k)  +byp*fin(4,i,j+1,k))+ &
             bx*(bybarp*fin(4,i+1,j,k)+byp*fin(4,i+1,j+1,k))) &
             + azp*( &
             bxbar*(bybarp*fin(4,i,j,k+1)  +byp*fin(4,i,j+1,k+1))+ &
             bx*(bybarp*fin(4,i+1,j,k+1)+byp*fin(4,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hx(v)*hyi(v)*( &
             bzbarp*( &
             bxbar*(aybarp*fin(5,i,j,k)  +ayp*fin(5,i,j+1,k))+ &
             bx*(aybarp*fin(5,i+1,j,k)+ayp*fin(5,i+1,j+1,k))) &
             + bzp*( &
             bxbar*(aybarp*fin(5,i,j,k+1)  +ayp*fin(5,i,j+1,k+1))+ &
             bx*(aybarp*fin(5,i+1,j,k+1)+ayp*fin(5,i+1,j+1,k+1))) &
             )
        !
        sum=sum+( &
             bzbarp*( &
             axbar*(bybarp*fin(6,i,j,k)  +byp*fin(6,i,j+1,k))+ &
             ax*(bybarp*fin(6,i+1,j,k)+byp*fin(6,i+1,j+1,k))) &
             + bzp*( &
             axbar*(bybarp*fin(6,i,j,k+1)  +byp*fin(6,i,j+1,k+1))+ &
             ax*(bybarp*fin(6,i+1,j,k+1)+byp*fin(6,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hx(v)*( &
             bzbarp*( &
             bxbar*(bybarp*fin(7,i,j,k)  +byp*fin(7,i,j+1,k))+ &
             bx*(bybarp*fin(7,i+1,j,k)+byp*fin(7,i+1,j+1,k))) &
             + bzp*( &
             bxbar*(bybarp*fin(7,i,j,k+1)  +byp*fin(7,i,j+1,k+1))+ &
             bx*(bybarp*fin(7,i+1,j,k+1)+byp*fin(7,i+1,j+1,k+1))) &
             )
        !
        fval(v,iadr)=sum
     end if
     !
     if(ict(8).eq.1) then
        !
        !  d3f/dxdydz:
        !
        iadr=iadr+1
        !
        sum=hxi(v)*hyi(v)*hzi(v)*( &
             azbarp*( &
             axbarp*(aybarp*fin(0,i,j,k)  +ayp*fin(0,i,j+1,k))+ &
             axp*(aybarp*fin(0,i+1,j,k)+ayp*fin(0,i+1,j+1,k))) &
             +  azp*( &
             axbarp*(aybarp*fin(0,i,j,k+1)  +ayp*fin(0,i,j+1,k+1))+ &
             axp*(aybarp*fin(0,i+1,j,k+1)+ayp*fin(0,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hyi(v)*hzi(v)*( &
             azbarp*( &
             bxbarp*(aybarp*fin(1,i,j,k)  +ayp*fin(1,i,j+1,k))+ &
             bxp*(aybarp*fin(1,i+1,j,k)+ayp*fin(1,i+1,j+1,k))) &
             + azp*( &
             bxbarp*(aybarp*fin(1,i,j,k+1)  +ayp*fin(1,i,j+1,k+1))+ &
             bxp*(aybarp*fin(1,i+1,j,k+1)+ayp*fin(1,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hxi(v)*hzi(v)*( &
             azbarp*( &
             axbarp*(bybarp*fin(2,i,j,k)  +byp*fin(2,i,j+1,k))+ &
             axp*(bybarp*fin(2,i+1,j,k)+byp*fin(2,i+1,j+1,k))) &
             + azp*( &
             axbarp*(bybarp*fin(2,i,j,k+1)  +byp*fin(2,i,j+1,k+1))+ &
             axp*(bybarp*fin(2,i+1,j,k+1)+byp*fin(2,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hxi(v)*hyi(v)*( &
             bzbarp*( &
             axbarp*(aybarp*fin(3,i,j,k)  +ayp*fin(3,i,j+1,k))+ &
             axp*(aybarp*fin(3,i+1,j,k)+ayp*fin(3,i+1,j+1,k))) &
             + bzp*( &
             axbarp*(aybarp*fin(3,i,j,k+1)  +ayp*fin(3,i,j+1,k+1))+ &
             axp*(aybarp*fin(3,i+1,j,k+1)+ayp*fin(3,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hzi(v)*( &
             azbarp*( &
             bxbarp*(bybarp*fin(4,i,j,k)  +byp*fin(4,i,j+1,k))+ &
             bxp*(bybarp*fin(4,i+1,j,k)+byp*fin(4,i+1,j+1,k))) &
             + azp*( &
             bxbarp*(bybarp*fin(4,i,j,k+1)  +byp*fin(4,i,j+1,k+1))+ &
             bxp*(bybarp*fin(4,i+1,j,k+1)+byp*fin(4,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hyi(v)*( &
             bzbarp*( &
             bxbarp*(aybarp*fin(5,i,j,k)  +ayp*fin(5,i,j+1,k))+ &
             bxp*(aybarp*fin(5,i+1,j,k)+ayp*fin(5,i+1,j+1,k))) &
             + bzp*( &
             bxbarp*(aybarp*fin(5,i,j,k+1)  +ayp*fin(5,i,j+1,k+1))+ &
             bxp*(aybarp*fin(5,i+1,j,k+1)+ayp*fin(5,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hxi(v)*( &
             bzbarp*( &
             axbarp*(bybarp*fin(6,i,j,k)  +byp*fin(6,i,j+1,k))+ &
             axp*(bybarp*fin(6,i+1,j,k)+byp*fin(6,i+1,j+1,k))) &
             + bzp*( &
             axbarp*(bybarp*fin(6,i,j,k+1)  +byp*fin(6,i,j+1,k+1))+ &
             axp*(bybarp*fin(6,i+1,j,k+1)+byp*fin(6,i+1,j+1,k+1))) &
             )
        !
        sum=sum+( &
             bzbarp*( &
             bxbarp*(bybarp*fin(7,i,j,k)  +byp*fin(7,i,j+1,k))+ &
             bxp*(bybarp*fin(7,i+1,j,k)+byp*fin(7,i+1,j+1,k))) &
             + bzp*( &
             bxbarp*(bybarp*fin(7,i,j,k+1)  +byp*fin(7,i,j+1,k+1))+ &
             bxp*(bybarp*fin(7,i+1,j,k+1)+byp*fin(7,i+1,j+1,k+1))) &
             )
        !
        fval(v,iadr)=sum
     end if
     !
  end do                             ! vector loop
  !
  return
end subroutine herm3fcn
!---------------------------------------------------------------------
!  evaluate C1 cubic Hermite function interpolation -- 3d fcn
!   --vectorized-- dmc 10 Feb 1999
!    --optimized for VARIATION along x axis ONLY--
!
subroutine herm3fcnx(ict,ivec,ivecd, &
     fval,ii,jj,kk,xparam,yparam,zparam, &
     hx,hxi,hy,hyi,hz,hzi, &
     fin,inf2,inf3,nz)
  use psp_precision_mod, only: fp
  !
  !============
  implicit none
  integer inf3,nz,inf2,j,k,i,iadr
  !============
  real(fp) :: yp,ypi,yp2,ypi2,ay,aybar,by,bybar,zp,zpi,zp2,zpi2,az
  real(fp) :: azbar,bz,bzbar,ayp,aybarp,byp,bybarp,azp,azbarp,bzp
  real(fp) :: bzbarp,xp,xpi,xp2,xpi2,ax,axbar,bx,bxbar,axp,axbarp
  real(fp) :: bxp,bxbarp
  !============
  integer ict(8)                    ! requested output control
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
  real(fp) :: fin(0:7,inf2,inf3,nz)        ! interpolant data (cf "herm3ev")
  !
  real(fp) :: fval(ivecd,*)                ! output returned
  !
  !  for detailed description of fin, ict and fval see subroutine herm3ev
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
  !  Hermite cubic basis functions
  !  -->for function value matching
  !     a(0)=0 a(1)=1        a'(0)=0 a'(1)=0
  !   abar(0)=1 abar(1)=0  abar'(0)=0 abar'(1)=0
  !
  !   a(x)=-2*x**3 + 3*x**2    = x*x*(-2.0*x+3.0)
  !   abar(x)=1-a(x)
  !   a'(x)=-abar'(x)          = 6.0*x*(1.0-x)
  !
  !  -->for derivative matching
  !     b(0)=0 b(1)=0          b'(0)=0 b'(1)=1
  !   bbar(0)=0 bbar(1)=0  bbar'(0)=1 bbar'(1)=0
  !
  !     b(x)=x**3-x**2         b'(x)=3*x**2-2*x
  !     bbar(x)=x**3-2*x**2+x  bbar'(x)=3*x**2-4*x+1
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
  yp2=yp*yp
  ypi2=ypi*ypi
  ay=yp2*(3.0_fp-2.0_fp*yp)
  aybar=1.0_fp-ay
  by=-yp2*ypi
  bybar=ypi2*yp
  !
  !   ...in z direction
  !
  zp=zparam
  zpi=1.0_fp-zp
  zp2=zp*zp
  zpi2=zpi*zpi
  az=zp2*(3.0_fp-2.0_fp*zp)
  azbar=1.0_fp-az
  bz=-zp2*zpi
  bzbar=zpi2*zp
  !
  !  derivatives
  !
  ayp=6.0_fp*yp*ypi
  aybarp=-ayp
  byp=yp*(3.0_fp*yp-2.0_fp)
  bybarp=ypi*(3.0_fp*ypi-2.0_fp)
  !
  azp=6.0_fp*zp*zpi
  azbarp=-azp
  bzp=zp*(3.0_fp*zp-2.0_fp)
  bzbarp=zpi*(3.0_fp*zpi-2.0_fp)

  do v=1,ivec
     i=ii(v)
     !
     !   ...in x direction
     !
     xp=xparam(v)
     xpi=1.0_fp-xp
     xp2=xp*xp
     xpi2=xpi*xpi
     ax=xp2*(3.0_fp-2.0_fp*xp)
     axbar=1.0_fp-ax
     bx=-xp2*xpi
     bxbar=xpi2*xp
     !
     iadr=0
     !
     !  derivatives:
     !
     axp=6.0_fp*xp*xpi
     axbarp=-axp
     bxp=xp*(3.0_fp*xp-2.0_fp)
     bxbarp=xpi*(3.0_fp*xpi-2.0_fp)
     !
     !  get desired values:
     !
     if(ict(1).eq.1) then
        !
        !  function value:
        !
        iadr=iadr+1
        sum=azbar*( &
             axbar*(aybar*fin(0,i,j,k)  +ay*fin(0,i,j+1,k))+ &
             ax*(aybar*fin(0,i+1,j,k)+ay*fin(0,i+1,j+1,k))) &
             +  az*( &
             axbar*(aybar*fin(0,i,j,k+1)  +ay*fin(0,i,j+1,k+1))+ &
             ax*(aybar*fin(0,i+1,j,k+1)+ay*fin(0,i+1,j+1,k+1)))
        !
        sum=sum+hx(v)*( &
             azbar*( &
             bxbar*(aybar*fin(1,i,j,k)  +ay*fin(1,i,j+1,k))+ &
             bx*(aybar*fin(1,i+1,j,k)+ay*fin(1,i+1,j+1,k))) &
             + az*( &
             bxbar*(aybar*fin(1,i,j,k+1)  +ay*fin(1,i,j+1,k+1))+ &
             bx*(aybar*fin(1,i+1,j,k+1)+ay*fin(1,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hy*( &
             azbar*( &
             axbar*(bybar*fin(2,i,j,k)  +by*fin(2,i,j+1,k))+ &
             ax*(bybar*fin(2,i+1,j,k)+by*fin(2,i+1,j+1,k))) &
             + az*( &
             axbar*(bybar*fin(2,i,j,k+1)  +by*fin(2,i,j+1,k+1))+ &
             ax*(bybar*fin(2,i+1,j,k+1)+by*fin(2,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hz*( &
             bzbar*( &
             axbar*(aybar*fin(3,i,j,k)  +ay*fin(3,i,j+1,k))+ &
             ax*(aybar*fin(3,i+1,j,k)+ay*fin(3,i+1,j+1,k))) &
             + bz*( &
             axbar*(aybar*fin(3,i,j,k+1)  +ay*fin(3,i,j+1,k+1))+ &
             ax*(aybar*fin(3,i+1,j,k+1)+ay*fin(3,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hx(v)*hy*( &
             azbar*( &
             bxbar*(bybar*fin(4,i,j,k)  +by*fin(4,i,j+1,k))+ &
             bx*(bybar*fin(4,i+1,j,k)+by*fin(4,i+1,j+1,k))) &
             + az*( &
             bxbar*(bybar*fin(4,i,j,k+1)  +by*fin(4,i,j+1,k+1))+ &
             bx*(bybar*fin(4,i+1,j,k+1)+by*fin(4,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hx(v)*hz*( &
             bzbar*( &
             bxbar*(aybar*fin(5,i,j,k)  +ay*fin(5,i,j+1,k))+ &
             bx*(aybar*fin(5,i+1,j,k)+ay*fin(5,i+1,j+1,k))) &
             + bz*( &
             bxbar*(aybar*fin(5,i,j,k+1)  +ay*fin(5,i,j+1,k+1))+ &
             bx*(aybar*fin(5,i+1,j,k+1)+ay*fin(5,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hy*hz*( &
             bzbar*( &
             axbar*(bybar*fin(6,i,j,k)  +by*fin(6,i,j+1,k))+ &
             ax*(bybar*fin(6,i+1,j,k)+by*fin(6,i+1,j+1,k))) &
             + bz*( &
             axbar*(bybar*fin(6,i,j,k+1)  +by*fin(6,i,j+1,k+1))+ &
             ax*(bybar*fin(6,i+1,j,k+1)+by*fin(6,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hx(v)*hy*hz*( &
             bzbar*( &
             bxbar*(bybar*fin(7,i,j,k)  +by*fin(7,i,j+1,k))+ &
             bx*(bybar*fin(7,i+1,j,k)+by*fin(7,i+1,j+1,k))) &
             + bz*( &
             bxbar*(bybar*fin(7,i,j,k+1)  +by*fin(7,i,j+1,k+1))+ &
             bx*(bybar*fin(7,i+1,j,k+1)+by*fin(7,i+1,j+1,k+1))) &
             )
        !
        fval(v,iadr)=sum
     end if
     !
     if(ict(2).eq.1) then
        !
        !     df/dx:
        !
        iadr=iadr+1
        !
        sum=hxi(v)*( &
             azbar*( &
             axbarp*(aybar*fin(0,i,j,k)  +ay*fin(0,i,j+1,k))+ &
             axp*(aybar*fin(0,i+1,j,k)+ay*fin(0,i+1,j+1,k))) &
             +  az*( &
             axbarp*(aybar*fin(0,i,j,k+1)  +ay*fin(0,i,j+1,k+1))+ &
             axp*(aybar*fin(0,i+1,j,k+1)+ay*fin(0,i+1,j+1,k+1))) &
             )
        !
        sum=sum+( &
             azbar*( &
             bxbarp*(aybar*fin(1,i,j,k)  +ay*fin(1,i,j+1,k))+ &
             bxp*(aybar*fin(1,i+1,j,k)+ay*fin(1,i+1,j+1,k))) &
             + az*( &
             bxbarp*(aybar*fin(1,i,j,k+1)  +ay*fin(1,i,j+1,k+1))+ &
             bxp*(aybar*fin(1,i+1,j,k+1)+ay*fin(1,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hxi(v)*hy*( &
             azbar*( &
             axbarp*(bybar*fin(2,i,j,k)  +by*fin(2,i,j+1,k))+ &
             axp*(bybar*fin(2,i+1,j,k)+by*fin(2,i+1,j+1,k))) &
             + az*( &
             axbarp*(bybar*fin(2,i,j,k+1)  +by*fin(2,i,j+1,k+1))+ &
             axp*(bybar*fin(2,i+1,j,k+1)+by*fin(2,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hxi(v)*hz*( &
             bzbar*( &
             axbarp*(aybar*fin(3,i,j,k)  +ay*fin(3,i,j+1,k))+ &
             axp*(aybar*fin(3,i+1,j,k)+ay*fin(3,i+1,j+1,k))) &
             + bz*( &
             axbarp*(aybar*fin(3,i,j,k+1)  +ay*fin(3,i,j+1,k+1))+ &
             axp*(aybar*fin(3,i+1,j,k+1)+ay*fin(3,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hy*( &
             azbar*( &
             bxbarp*(bybar*fin(4,i,j,k)  +by*fin(4,i,j+1,k))+ &
             bxp*(bybar*fin(4,i+1,j,k)+by*fin(4,i+1,j+1,k))) &
             + az*( &
             bxbarp*(bybar*fin(4,i,j,k+1)  +by*fin(4,i,j+1,k+1))+ &
             bxp*(bybar*fin(4,i+1,j,k+1)+by*fin(4,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hz*( &
             bzbar*( &
             bxbarp*(aybar*fin(5,i,j,k)  +ay*fin(5,i,j+1,k))+ &
             bxp*(aybar*fin(5,i+1,j,k)+ay*fin(5,i+1,j+1,k))) &
             + bz*( &
             bxbarp*(aybar*fin(5,i,j,k+1)  +ay*fin(5,i,j+1,k+1))+ &
             bxp*(aybar*fin(5,i+1,j,k+1)+ay*fin(5,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hxi(v)*hy*hz*( &
             bzbar*( &
             axbarp*(bybar*fin(6,i,j,k)  +by*fin(6,i,j+1,k))+ &
             axp*(bybar*fin(6,i+1,j,k)+by*fin(6,i+1,j+1,k))) &
             + bz*( &
             axbarp*(bybar*fin(6,i,j,k+1)  +by*fin(6,i,j+1,k+1))+ &
             axp*(bybar*fin(6,i+1,j,k+1)+by*fin(6,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hy*hz*( &
             bzbar*( &
             bxbarp*(bybar*fin(7,i,j,k)  +by*fin(7,i,j+1,k))+ &
             bxp*(bybar*fin(7,i+1,j,k)+by*fin(7,i+1,j+1,k))) &
             + bz*( &
             bxbarp*(bybar*fin(7,i,j,k+1)  +by*fin(7,i,j+1,k+1))+ &
             bxp*(bybar*fin(7,i+1,j,k+1)+by*fin(7,i+1,j+1,k+1))) &
             )
        !
        fval(v,iadr)=sum
     end if
     !
     if(ict(3).eq.1) then
        !
        !  df/dy:
        !
        iadr=iadr+1
        !
        sum=hyi*( &
             azbar*( &
             axbar*(aybarp*fin(0,i,j,k)  +ayp*fin(0,i,j+1,k))+ &
             ax*(aybarp*fin(0,i+1,j,k)+ayp*fin(0,i+1,j+1,k))) &
             +  az*( &
             axbar*(aybarp*fin(0,i,j,k+1)  +ayp*fin(0,i,j+1,k+1))+ &
             ax*(aybarp*fin(0,i+1,j,k+1)+ayp*fin(0,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hyi*hx(v)*( &
             azbar*( &
             bxbar*(aybarp*fin(1,i,j,k)  +ayp*fin(1,i,j+1,k))+ &
             bx*(aybarp*fin(1,i+1,j,k)+ayp*fin(1,i+1,j+1,k))) &
             + az*( &
             bxbar*(aybarp*fin(1,i,j,k+1)  +ayp*fin(1,i,j+1,k+1))+ &
             bx*(aybarp*fin(1,i+1,j,k+1)+ayp*fin(1,i+1,j+1,k+1))) &
             )
        !
        sum=sum+( &
             azbar*( &
             axbar*(bybarp*fin(2,i,j,k)  +byp*fin(2,i,j+1,k))+ &
             ax*(bybarp*fin(2,i+1,j,k)+byp*fin(2,i+1,j+1,k))) &
             + az*( &
             axbar*(bybarp*fin(2,i,j,k+1)  +byp*fin(2,i,j+1,k+1))+ &
             ax*(bybarp*fin(2,i+1,j,k+1)+byp*fin(2,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hyi*hz*( &
             bzbar*( &
             axbar*(aybarp*fin(3,i,j,k)  +ayp*fin(3,i,j+1,k))+ &
             ax*(aybarp*fin(3,i+1,j,k)+ayp*fin(3,i+1,j+1,k))) &
             + bz*( &
             axbar*(aybarp*fin(3,i,j,k+1)  +ayp*fin(3,i,j+1,k+1))+ &
             ax*(aybarp*fin(3,i+1,j,k+1)+ayp*fin(3,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hx(v)*( &
             azbar*( &
             bxbar*(bybarp*fin(4,i,j,k)  +byp*fin(4,i,j+1,k))+ &
             bx*(bybarp*fin(4,i+1,j,k)+byp*fin(4,i+1,j+1,k))) &
             + az*( &
             bxbar*(bybarp*fin(4,i,j,k+1)  +byp*fin(4,i,j+1,k+1))+ &
             bx*(bybarp*fin(4,i+1,j,k+1)+byp*fin(4,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hx(v)*hyi*hz*( &
             bzbar*( &
             bxbar*(aybarp*fin(5,i,j,k)  +ayp*fin(5,i,j+1,k))+ &
             bx*(aybarp*fin(5,i+1,j,k)+ayp*fin(5,i+1,j+1,k))) &
             + bz*( &
             bxbar*(aybarp*fin(5,i,j,k+1)  +ayp*fin(5,i,j+1,k+1))+ &
             bx*(aybarp*fin(5,i+1,j,k+1)+ayp*fin(5,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hz*( &
             bzbar*( &
             axbar*(bybarp*fin(6,i,j,k)  +byp*fin(6,i,j+1,k))+ &
             ax*(bybarp*fin(6,i+1,j,k)+byp*fin(6,i+1,j+1,k))) &
             + bz*( &
             axbar*(bybarp*fin(6,i,j,k+1)  +byp*fin(6,i,j+1,k+1))+ &
             ax*(bybarp*fin(6,i+1,j,k+1)+byp*fin(6,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hx(v)*hz*( &
             bzbar*( &
             bxbar*(bybarp*fin(7,i,j,k)  +byp*fin(7,i,j+1,k))+ &
             bx*(bybarp*fin(7,i+1,j,k)+byp*fin(7,i+1,j+1,k))) &
             + bz*( &
             bxbar*(bybarp*fin(7,i,j,k+1)  +byp*fin(7,i,j+1,k+1))+ &
             bx*(bybarp*fin(7,i+1,j,k+1)+byp*fin(7,i+1,j+1,k+1))) &
             )
        !
        fval(v,iadr)=sum
     end if
     !
     if(ict(4).eq.1) then
        !
        !  df/dz:
        !
        iadr=iadr+1
        !
        sum=hzi*( &
             azbarp*( &
             axbar*(aybar*fin(0,i,j,k)  +ay*fin(0,i,j+1,k))+ &
             ax*(aybar*fin(0,i+1,j,k)+ay*fin(0,i+1,j+1,k))) &
             +  azp*( &
             axbar*(aybar*fin(0,i,j,k+1)  +ay*fin(0,i,j+1,k+1))+ &
             ax*(aybar*fin(0,i+1,j,k+1)+ay*fin(0,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hzi*hx(v)*( &
             azbarp*( &
             bxbar*(aybar*fin(1,i,j,k)  +ay*fin(1,i,j+1,k))+ &
             bx*(aybar*fin(1,i+1,j,k)+ay*fin(1,i+1,j+1,k))) &
             + azp*( &
             bxbar*(aybar*fin(1,i,j,k+1)  +ay*fin(1,i,j+1,k+1))+ &
             bx*(aybar*fin(1,i+1,j,k+1)+ay*fin(1,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hzi*hy*( &
             azbarp*( &
             axbar*(bybar*fin(2,i,j,k)  +by*fin(2,i,j+1,k))+ &
             ax*(bybar*fin(2,i+1,j,k)+by*fin(2,i+1,j+1,k))) &
             + azp*( &
             axbar*(bybar*fin(2,i,j,k+1)  +by*fin(2,i,j+1,k+1))+ &
             ax*(bybar*fin(2,i+1,j,k+1)+by*fin(2,i+1,j+1,k+1))) &
             )
        !
        sum=sum+( &
             bzbarp*( &
             axbar*(aybar*fin(3,i,j,k)  +ay*fin(3,i,j+1,k))+ &
             ax*(aybar*fin(3,i+1,j,k)+ay*fin(3,i+1,j+1,k))) &
             + bzp*( &
             axbar*(aybar*fin(3,i,j,k+1)  +ay*fin(3,i,j+1,k+1))+ &
             ax*(aybar*fin(3,i+1,j,k+1)+ay*fin(3,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hzi*hx(v)*hy*( &
             azbarp*( &
             bxbar*(bybar*fin(4,i,j,k)  +by*fin(4,i,j+1,k))+ &
             bx*(bybar*fin(4,i+1,j,k)+by*fin(4,i+1,j+1,k))) &
             + azp*( &
             bxbar*(bybar*fin(4,i,j,k+1)  +by*fin(4,i,j+1,k+1))+ &
             bx*(bybar*fin(4,i+1,j,k+1)+by*fin(4,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hx(v)*( &
             bzbarp*( &
             bxbar*(aybar*fin(5,i,j,k)  +ay*fin(5,i,j+1,k))+ &
             bx*(aybar*fin(5,i+1,j,k)+ay*fin(5,i+1,j+1,k))) &
             + bzp*( &
             bxbar*(aybar*fin(5,i,j,k+1)  +ay*fin(5,i,j+1,k+1))+ &
             bx*(aybar*fin(5,i+1,j,k+1)+ay*fin(5,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hy*( &
             bzbarp*( &
             axbar*(bybar*fin(6,i,j,k)  +by*fin(6,i,j+1,k))+ &
             ax*(bybar*fin(6,i+1,j,k)+by*fin(6,i+1,j+1,k))) &
             + bzp*( &
             axbar*(bybar*fin(6,i,j,k+1)  +by*fin(6,i,j+1,k+1))+ &
             ax*(bybar*fin(6,i+1,j,k+1)+by*fin(6,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hx(v)*hy*( &
             bzbarp*( &
             bxbar*(bybar*fin(7,i,j,k)  +by*fin(7,i,j+1,k))+ &
             bx*(bybar*fin(7,i+1,j,k)+by*fin(7,i+1,j+1,k))) &
             + bzp*( &
             bxbar*(bybar*fin(7,i,j,k+1)  +by*fin(7,i,j+1,k+1))+ &
             bx*(bybar*fin(7,i+1,j,k+1)+by*fin(7,i+1,j+1,k+1))) &
             )
        !
        fval(v,iadr)=sum
     end if
     !
     if(ict(5).eq.1) then
        !
        !  d2f/dxdy:
        !
        iadr=iadr+1
        !
        sum=hxi(v)*hyi*( &
             azbar*( &
             axbarp*(aybarp*fin(0,i,j,k)  +ayp*fin(0,i,j+1,k))+ &
             axp*(aybarp*fin(0,i+1,j,k)+ayp*fin(0,i+1,j+1,k))) &
             +  az*( &
             axbarp*(aybarp*fin(0,i,j,k+1)  +ayp*fin(0,i,j+1,k+1))+ &
             axp*(aybarp*fin(0,i+1,j,k+1)+ayp*fin(0,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hyi*( &
             azbar*( &
             bxbarp*(aybarp*fin(1,i,j,k)  +ayp*fin(1,i,j+1,k))+ &
             bxp*(aybarp*fin(1,i+1,j,k)+ayp*fin(1,i+1,j+1,k))) &
             + az*( &
             bxbarp*(aybarp*fin(1,i,j,k+1)  +ayp*fin(1,i,j+1,k+1))+ &
             bxp*(aybarp*fin(1,i+1,j,k+1)+ayp*fin(1,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hxi(v)*( &
             azbar*( &
             axbarp*(bybarp*fin(2,i,j,k)  +byp*fin(2,i,j+1,k))+ &
             axp*(bybarp*fin(2,i+1,j,k)+byp*fin(2,i+1,j+1,k))) &
             + az*( &
             axbarp*(bybarp*fin(2,i,j,k+1)  +byp*fin(2,i,j+1,k+1))+ &
             axp*(bybarp*fin(2,i+1,j,k+1)+byp*fin(2,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hxi(v)*hyi*hz*( &
             bzbar*( &
             axbarp*(aybarp*fin(3,i,j,k)  +ayp*fin(3,i,j+1,k))+ &
             axp*(aybarp*fin(3,i+1,j,k)+ayp*fin(3,i+1,j+1,k))) &
             + bz*( &
             axbarp*(aybarp*fin(3,i,j,k+1)  +ayp*fin(3,i,j+1,k+1))+ &
             axp*(aybarp*fin(3,i+1,j,k+1)+ayp*fin(3,i+1,j+1,k+1))) &
             )
        !
        sum=sum+( &
             azbar*( &
             bxbarp*(bybarp*fin(4,i,j,k)  +byp*fin(4,i,j+1,k))+ &
             bxp*(bybarp*fin(4,i+1,j,k)+byp*fin(4,i+1,j+1,k))) &
             + az*( &
             bxbarp*(bybarp*fin(4,i,j,k+1)  +byp*fin(4,i,j+1,k+1))+ &
             bxp*(bybarp*fin(4,i+1,j,k+1)+byp*fin(4,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hyi*hz*( &
             bzbar*( &
             bxbarp*(aybarp*fin(5,i,j,k)  +ayp*fin(5,i,j+1,k))+ &
             bxp*(aybarp*fin(5,i+1,j,k)+ayp*fin(5,i+1,j+1,k))) &
             + bz*( &
             bxbarp*(aybarp*fin(5,i,j,k+1)  +ayp*fin(5,i,j+1,k+1))+ &
             bxp*(aybarp*fin(5,i+1,j,k+1)+ayp*fin(5,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hxi(v)*hz*( &
             bzbar*( &
             axbarp*(bybarp*fin(6,i,j,k)  +byp*fin(6,i,j+1,k))+ &
             axp*(bybarp*fin(6,i+1,j,k)+byp*fin(6,i+1,j+1,k))) &
             + bz*( &
             axbarp*(bybarp*fin(6,i,j,k+1)  +byp*fin(6,i,j+1,k+1))+ &
             axp*(bybarp*fin(6,i+1,j,k+1)+byp*fin(6,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hz*( &
             bzbar*( &
             bxbarp*(bybarp*fin(7,i,j,k)  +byp*fin(7,i,j+1,k))+ &
             bxp*(bybarp*fin(7,i+1,j,k)+byp*fin(7,i+1,j+1,k))) &
             + bz*( &
             bxbarp*(bybarp*fin(7,i,j,k+1)  +byp*fin(7,i,j+1,k+1))+ &
             bxp*(bybarp*fin(7,i+1,j,k+1)+byp*fin(7,i+1,j+1,k+1))) &
             )
        !
        fval(v,iadr)=sum
     end if
     !
     if(ict(6).eq.1) then
        !
        !  d2f/dxdz:
        !
        iadr=iadr+1
        !
        sum=hxi(v)*hzi*( &
             azbarp*( &
             axbarp*(aybar*fin(0,i,j,k)  +ay*fin(0,i,j+1,k))+ &
             axp*(aybar*fin(0,i+1,j,k)+ay*fin(0,i+1,j+1,k))) &
             +  azp*( &
             axbarp*(aybar*fin(0,i,j,k+1)  +ay*fin(0,i,j+1,k+1))+ &
             axp*(aybar*fin(0,i+1,j,k+1)+ay*fin(0,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hzi*( &
             azbarp*( &
             bxbarp*(aybar*fin(1,i,j,k)  +ay*fin(1,i,j+1,k))+ &
             bxp*(aybar*fin(1,i+1,j,k)+ay*fin(1,i+1,j+1,k))) &
             + azp*( &
             bxbarp*(aybar*fin(1,i,j,k+1)  +ay*fin(1,i,j+1,k+1))+ &
             bxp*(aybar*fin(1,i+1,j,k+1)+ay*fin(1,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hxi(v)*hy*hzi*( &
             azbarp*( &
             axbarp*(bybar*fin(2,i,j,k)  +by*fin(2,i,j+1,k))+ &
             axp*(bybar*fin(2,i+1,j,k)+by*fin(2,i+1,j+1,k))) &
             + azp*( &
             axbarp*(bybar*fin(2,i,j,k+1)  +by*fin(2,i,j+1,k+1))+ &
             axp*(bybar*fin(2,i+1,j,k+1)+by*fin(2,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hxi(v)*( &
             bzbarp*( &
             axbarp*(aybar*fin(3,i,j,k)  +ay*fin(3,i,j+1,k))+ &
             axp*(aybar*fin(3,i+1,j,k)+ay*fin(3,i+1,j+1,k))) &
             + bzp*( &
             axbarp*(aybar*fin(3,i,j,k+1)  +ay*fin(3,i,j+1,k+1))+ &
             axp*(aybar*fin(3,i+1,j,k+1)+ay*fin(3,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hy*hzi*( &
             azbarp*( &
             bxbarp*(bybar*fin(4,i,j,k)  +by*fin(4,i,j+1,k))+ &
             bxp*(bybar*fin(4,i+1,j,k)+by*fin(4,i+1,j+1,k))) &
             + azp*( &
             bxbarp*(bybar*fin(4,i,j,k+1)  +by*fin(4,i,j+1,k+1))+ &
             bxp*(bybar*fin(4,i+1,j,k+1)+by*fin(4,i+1,j+1,k+1))) &
             )
        !
        sum=sum+( &
             bzbarp*( &
             bxbarp*(aybar*fin(5,i,j,k)  +ay*fin(5,i,j+1,k))+ &
             bxp*(aybar*fin(5,i+1,j,k)+ay*fin(5,i+1,j+1,k))) &
             + bzp*( &
             bxbarp*(aybar*fin(5,i,j,k+1)  +ay*fin(5,i,j+1,k+1))+ &
             bxp*(aybar*fin(5,i+1,j,k+1)+ay*fin(5,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hxi(v)*hy*( &
             bzbarp*( &
             axbarp*(bybar*fin(6,i,j,k)  +by*fin(6,i,j+1,k))+ &
             axp*(bybar*fin(6,i+1,j,k)+by*fin(6,i+1,j+1,k))) &
             + bzp*( &
             axbarp*(bybar*fin(6,i,j,k+1)  +by*fin(6,i,j+1,k+1))+ &
             axp*(bybar*fin(6,i+1,j,k+1)+by*fin(6,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hy*( &
             bzbarp*( &
             bxbarp*(bybar*fin(7,i,j,k)  +by*fin(7,i,j+1,k))+ &
             bxp*(bybar*fin(7,i+1,j,k)+by*fin(7,i+1,j+1,k))) &
             + bzp*( &
             bxbarp*(bybar*fin(7,i,j,k+1)  +by*fin(7,i,j+1,k+1))+ &
             bxp*(bybar*fin(7,i+1,j,k+1)+by*fin(7,i+1,j+1,k+1))) &
             )
        !
        fval(v,iadr)=sum
     end if
     !
     if(ict(7).eq.1) then
        !
        !  d2f/dydz:
        !
        iadr=iadr+1
        !
        sum=hyi*hzi*( &
             azbarp*( &
             axbar*(aybarp*fin(0,i,j,k)  +ayp*fin(0,i,j+1,k))+ &
             ax*(aybarp*fin(0,i+1,j,k)+ayp*fin(0,i+1,j+1,k))) &
             +  azp*( &
             axbar*(aybarp*fin(0,i,j,k+1)  +ayp*fin(0,i,j+1,k+1))+ &
             ax*(aybarp*fin(0,i+1,j,k+1)+ayp*fin(0,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hyi*hzi*hx(v)*( &
             azbarp*( &
             bxbar*(aybarp*fin(1,i,j,k)  +ayp*fin(1,i,j+1,k))+ &
             bx*(aybarp*fin(1,i+1,j,k)+ayp*fin(1,i+1,j+1,k))) &
             + azp*( &
             bxbar*(aybarp*fin(1,i,j,k+1)  +ayp*fin(1,i,j+1,k+1))+ &
             bx*(aybarp*fin(1,i+1,j,k+1)+ayp*fin(1,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hzi*( &
             azbarp*( &
             axbar*(bybarp*fin(2,i,j,k)  +byp*fin(2,i,j+1,k))+ &
             ax*(bybarp*fin(2,i+1,j,k)+byp*fin(2,i+1,j+1,k))) &
             + azp*( &
             axbar*(bybarp*fin(2,i,j,k+1)  +byp*fin(2,i,j+1,k+1))+ &
             ax*(bybarp*fin(2,i+1,j,k+1)+byp*fin(2,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hyi*( &
             bzbarp*( &
             axbar*(aybarp*fin(3,i,j,k)  +ayp*fin(3,i,j+1,k))+ &
             ax*(aybarp*fin(3,i+1,j,k)+ayp*fin(3,i+1,j+1,k))) &
             + bzp*( &
             axbar*(aybarp*fin(3,i,j,k+1)  +ayp*fin(3,i,j+1,k+1))+ &
             ax*(aybarp*fin(3,i+1,j,k+1)+ayp*fin(3,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hx(v)*hzi*( &
             azbarp*( &
             bxbar*(bybarp*fin(4,i,j,k)  +byp*fin(4,i,j+1,k))+ &
             bx*(bybarp*fin(4,i+1,j,k)+byp*fin(4,i+1,j+1,k))) &
             + azp*( &
             bxbar*(bybarp*fin(4,i,j,k+1)  +byp*fin(4,i,j+1,k+1))+ &
             bx*(bybarp*fin(4,i+1,j,k+1)+byp*fin(4,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hx(v)*hyi*( &
             bzbarp*( &
             bxbar*(aybarp*fin(5,i,j,k)  +ayp*fin(5,i,j+1,k))+ &
             bx*(aybarp*fin(5,i+1,j,k)+ayp*fin(5,i+1,j+1,k))) &
             + bzp*( &
             bxbar*(aybarp*fin(5,i,j,k+1)  +ayp*fin(5,i,j+1,k+1))+ &
             bx*(aybarp*fin(5,i+1,j,k+1)+ayp*fin(5,i+1,j+1,k+1))) &
             )
        !
        sum=sum+( &
             bzbarp*( &
             axbar*(bybarp*fin(6,i,j,k)  +byp*fin(6,i,j+1,k))+ &
             ax*(bybarp*fin(6,i+1,j,k)+byp*fin(6,i+1,j+1,k))) &
             + bzp*( &
             axbar*(bybarp*fin(6,i,j,k+1)  +byp*fin(6,i,j+1,k+1))+ &
             ax*(bybarp*fin(6,i+1,j,k+1)+byp*fin(6,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hx(v)*( &
             bzbarp*( &
             bxbar*(bybarp*fin(7,i,j,k)  +byp*fin(7,i,j+1,k))+ &
             bx*(bybarp*fin(7,i+1,j,k)+byp*fin(7,i+1,j+1,k))) &
             + bzp*( &
             bxbar*(bybarp*fin(7,i,j,k+1)  +byp*fin(7,i,j+1,k+1))+ &
             bx*(bybarp*fin(7,i+1,j,k+1)+byp*fin(7,i+1,j+1,k+1))) &
             )
        !
        fval(v,iadr)=sum
     end if
     !
     if(ict(8).eq.1) then
        !
        !  d3f/dxdydz:
        !
        iadr=iadr+1
        !
        sum=hxi(v)*hyi*hzi*( &
             azbarp*( &
             axbarp*(aybarp*fin(0,i,j,k)  +ayp*fin(0,i,j+1,k))+ &
             axp*(aybarp*fin(0,i+1,j,k)+ayp*fin(0,i+1,j+1,k))) &
             +  azp*( &
             axbarp*(aybarp*fin(0,i,j,k+1)  +ayp*fin(0,i,j+1,k+1))+ &
             axp*(aybarp*fin(0,i+1,j,k+1)+ayp*fin(0,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hyi*hzi*( &
             azbarp*( &
             bxbarp*(aybarp*fin(1,i,j,k)  +ayp*fin(1,i,j+1,k))+ &
             bxp*(aybarp*fin(1,i+1,j,k)+ayp*fin(1,i+1,j+1,k))) &
             + azp*( &
             bxbarp*(aybarp*fin(1,i,j,k+1)  +ayp*fin(1,i,j+1,k+1))+ &
             bxp*(aybarp*fin(1,i+1,j,k+1)+ayp*fin(1,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hxi(v)*hzi*( &
             azbarp*( &
             axbarp*(bybarp*fin(2,i,j,k)  +byp*fin(2,i,j+1,k))+ &
             axp*(bybarp*fin(2,i+1,j,k)+byp*fin(2,i+1,j+1,k))) &
             + azp*( &
             axbarp*(bybarp*fin(2,i,j,k+1)  +byp*fin(2,i,j+1,k+1))+ &
             axp*(bybarp*fin(2,i+1,j,k+1)+byp*fin(2,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hxi(v)*hyi*( &
             bzbarp*( &
             axbarp*(aybarp*fin(3,i,j,k)  +ayp*fin(3,i,j+1,k))+ &
             axp*(aybarp*fin(3,i+1,j,k)+ayp*fin(3,i+1,j+1,k))) &
             + bzp*( &
             axbarp*(aybarp*fin(3,i,j,k+1)  +ayp*fin(3,i,j+1,k+1))+ &
             axp*(aybarp*fin(3,i+1,j,k+1)+ayp*fin(3,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hzi*( &
             azbarp*( &
             bxbarp*(bybarp*fin(4,i,j,k)  +byp*fin(4,i,j+1,k))+ &
             bxp*(bybarp*fin(4,i+1,j,k)+byp*fin(4,i+1,j+1,k))) &
             + azp*( &
             bxbarp*(bybarp*fin(4,i,j,k+1)  +byp*fin(4,i,j+1,k+1))+ &
             bxp*(bybarp*fin(4,i+1,j,k+1)+byp*fin(4,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hyi*( &
             bzbarp*( &
             bxbarp*(aybarp*fin(5,i,j,k)  +ayp*fin(5,i,j+1,k))+ &
             bxp*(aybarp*fin(5,i+1,j,k)+ayp*fin(5,i+1,j+1,k))) &
             + bzp*( &
             bxbarp*(aybarp*fin(5,i,j,k+1)  +ayp*fin(5,i,j+1,k+1))+ &
             bxp*(aybarp*fin(5,i+1,j,k+1)+ayp*fin(5,i+1,j+1,k+1))) &
             )
        !
        sum=sum+hxi(v)*( &
             bzbarp*( &
             axbarp*(bybarp*fin(6,i,j,k)  +byp*fin(6,i,j+1,k))+ &
             axp*(bybarp*fin(6,i+1,j,k)+byp*fin(6,i+1,j+1,k))) &
             + bzp*( &
             axbarp*(bybarp*fin(6,i,j,k+1)  +byp*fin(6,i,j+1,k+1))+ &
             axp*(bybarp*fin(6,i+1,j,k+1)+byp*fin(6,i+1,j+1,k+1))) &
             )
        !
        sum=sum+( &
             bzbarp*( &
             bxbarp*(bybarp*fin(7,i,j,k)  +byp*fin(7,i,j+1,k))+ &
             bxp*(bybarp*fin(7,i+1,j,k)+byp*fin(7,i+1,j+1,k))) &
             + bzp*( &
             bxbarp*(bybarp*fin(7,i,j,k+1)  +byp*fin(7,i,j+1,k+1))+ &
             bxp*(bybarp*fin(7,i+1,j,k+1)+byp*fin(7,i+1,j+1,k+1))) &
             )
        !
        fval(v,iadr)=sum
     end if
     !
  end do                             ! vector loop
  !
  return
end subroutine herm3fcnx
