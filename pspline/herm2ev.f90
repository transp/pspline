subroutine herm2ev(xget,yget,x,nx,y,ny,ilinx,iliny, &
     f,inf2,ict,fval,ier)
  use precision_mod, only: fp
  !
  !  evaluate a 2d cubic Hermite interpolant on a rectilinear
  !  grid -- this is C1 in both directions.
  !
  !  this subroutine calls two subroutines:
  !     herm2xy  -- find cell containing (xget,yget)
  !     herm2fcn -- evaluate interpolant function and (optionally) derivatives
  !
  !  input arguments:
  !  ================
  !
  !============
  implicit none
  integer ny,inf2,nx
  !============
  real(fp) :: xget,yget                    ! target of this interpolation
  real(fp) :: x(nx)                        ! ordered x grid
  real(fp) :: y(ny)                        ! ordered y grid
  integer ilinx                     ! ilinx=1 => assume x evenly spaced
  integer iliny                     ! iliny=1 => assume y evenly spaced
  !
  real(fp) :: f(0:3,inf2,ny)               ! function data
  !
  !       f 2nd dimension inf2 must be .ge. nx
  !       contents of f:
  !
  !  f(0,i,j) = f @ x(i),y(j)
  !  f(1,i,j) = df/dx @ x(i),y(j)
  !  f(2,i,j) = df/dy @ x(i),y(j)
  !  f(3,i,j) = d2f/dxdy @ x(i),y(j)
  !
  integer ict(4)                    ! code specifying output desired
  !
  !  ict(1)=1 -- return f  (0, don't)
  !  ict(2)=1 -- return df/dx  (0, don't)
  !  ict(3)=1 -- return df/dy  (0, don't)
  !  ict(4)=1 -- return d2f/dxdy (0, don't)
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
  !  fval(4) receives the fourth output (depends on ict(...) spec)
  !
  !  examples:
  !    on input ict = [1,1,1,1]
  !   on output fval= [f,df/dx,df/dy,d2f/dxdy]
  !
  !    on input ict = [1,0,0,0]
  !   on output fval= [f] ... elements 2 & 3 & 4 never referenced
  !
  !    on input ict = [0,1,1,0]
  !   on output fval= [df/dx,df/dy] ... element 3 & 4 never referenced
  !
  !    on input ict = [0,0,1,0]
  !   on output fval= [df/dy] ... elements 2 & 3 & 4 never referenced.
  !
  !  ier -- completion code:  0 means OK
  !-------------------
  !  local:
  !
  integer, dimension(1) :: i,j                       ! cell indices
  !
  !  normalized displacement from (x(i),y(j)) corner of cell.
  !    xparam=0 @x(i)  xparam=1 @x(i+1)
  !    yparam=0 @y(j)  yparam=1 @y(j+1)
  !
  real(fp), dimension(1) :: xparam,yparam
  !
  !  cell dimensions and
  !  inverse cell dimensions hxi = 1/(x(i+1)-x(i)), hyi = 1/(y(j+1)-y(j))
  !
  real(fp), dimension(1) :: hx,hy
  real(fp), dimension(1) :: hxi,hyi
  !
  !  0 .le. xparam .le. 1
  !  0 .le. yparam .le. 1
  !
  !---------------------------------------------------------------------
  !
  call herm2xy(xget,yget,x,nx,y,ny,ilinx,iliny, &
       i(1),j(1),xparam(1),yparam(1), &
       hx(1),hxi(1),hy(1),hyi(1),ier)
  if(ier.ne.0) return
  !
  call herm2fcn(ict,1,1, &
       fval,i,j,xparam,yparam,hx,hxi,hy,hyi,f,inf2,ny)
  !
  return
end subroutine herm2ev
!---------------------------------------------------------------------
!  herm2xy -- look up x-y zone
!
!  this is the "first part" of herm2ev, see comments, above.
!
subroutine herm2xy(xget,yget,x,nx,y,ny,ilinx,iliny, &
     i,j,xparam,yparam,hx,hxi,hy,hyi,ier)
  use precision_mod, only: fp
  !
  !  input of herm2xy
  !  ================
  !
  !============
  implicit none
  integer nxm,nym,ii,jj
  !============
  real(fp) :: zxget,zyget,zxtol,zytol
  !============
  integer nx,ny                     ! array dimensions
  !
  real(fp) :: xget,yget                    ! target point
  real(fp) :: x(nx),y(ny)                  ! indep. coords., strict ascending
  !
  integer ilinx                     ! =1:  x evenly spaced
  integer iliny                     ! =1:  y evenly spaced
  !
  !  output of herm2xy
  !  =================
  integer i,j                       ! index to cell containing target pt
  !          on exit:  1.le.i.le.nx-1   1.le.j.le.ny-1
  !
  !  normalized position w/in (i,j) cell (btw 0 and 1):
  !
  real(fp) :: xparam                       ! (xget-x(i))/(x(i+1)-x(i))
  real(fp) :: yparam                       ! (yget-y(j))/(y(j+1)-y(j))
  !
  !  grid spacing
  !
  real(fp) :: hx                           ! hx = x(i+1)-x(i)
  real(fp) :: hy                           ! hy = y(j+1)-y(j)
  !
  !  inverse grid spacing:
  !
  real(fp) :: hxi                          ! 1/hx = 1/(x(i+1)-x(i))
  real(fp) :: hyi                          ! 1/hy = 1/(y(j+1)-y(j))
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
  if((xget.lt.x(1)).or.(xget.gt.x(nx))) then
     zxtol=4.0E-7_fp*max(abs(x(1)),abs(x(nx)))
     if((xget.lt.x(1)-zxtol).or.(xget.gt.x(nx)+zxtol)) then
        ier=1
        write(6,1001) xget,x(1),x(nx)
1001    format(' ?herm2ev:  xget=',1pe11.4,' out of range ', &
             1pe11.4,' to ',1pe11.4)
     else
        if((xget.lt.x(1)-0.5_fp*zxtol).or. &
             (xget.gt.x(nx)+0.5_fp*zxtol)) &
             write(6,1011) xget,x(1),x(nx)
1011    format(' %herm2ev:  xget=',1pe15.8,' beyond range ', &
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
1002    format(' ?herm2ev:  yget=',1pe11.4,' out of range ', &
             1pe11.4,' to ',1pe11.4)
     else
        if((yget.lt.y(1)-0.5_fp*zytol).or. &
             (yget.gt.y(ny)+0.5_fp*zytol)) &
             write(6,1012) yget,y(1),y(ny)
1012    format(' %herm2ev:  yget=',1pe15.8,' beyond range ', &
             1pe15.8,' to ',1pe15.8,' (fixup applied)')
        if(yget.lt.y(1)) then
           zyget=y(1)
        else
           zyget=y(ny)
        end if
     end if
  end if
  if(ier.ne.0) return
  !
  !  now find interval in which target point lies..
  !
  nxm=nx-1
  nym=ny-1
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
  hx=(x(i+1)-x(i))
  hy=(y(j+1)-y(j))
  !
  hxi=1.0_fp/hx
  hyi=1.0_fp/hy
  !
  xparam=(zxget-x(i))*hxi
  yparam=(zyget-y(j))*hyi
  !
  return
end subroutine herm2xy
!---------------------------------------------------------------------
!  evaluate C1 cubic Hermite function interpolation -- 2d fcn
!   --vectorized-- dmc 10 Feb 1999
!
subroutine herm2fcn(ict,ivec,ivecd, &
     fval,ii,jj,xparam,yparam,hx,hxi,hy,hyi, &
     fin,inf2,ny)
  use precision_mod, only: fp
  !
  !============
  implicit none
  integer ny,inf2,i,j,iadr
  !============
  real(fp) :: xp,xpi,xp2,xpi2,ax,axbar,bx,bxbar,yp,ypi,yp2,ypi2,ay
  real(fp) :: aybar,by,bybar,axp,axbarp,bxp,bxbarp,ayp,aybarp,byp
  real(fp) :: bybarp
  !============
  integer ict(4)                    ! requested output control
  integer ivec                      ! vector length
  integer ivecd                     ! vector dimension (1st dim of fval)
  !
  integer ii(ivec),jj(ivec)         ! target cells (i,j)
  real(fp) :: xparam(ivec),yparam(ivec)
  ! normalized displacements from (i,j) corners
  !
  real(fp) :: hx(ivec),hy(ivec)            ! grid spacing, and
  real(fp) :: hxi(ivec),hyi(ivec)          ! inverse grid spacing 1/(x(i+1)-x(i))
  ! & 1/(y(j+1)-y(j))
  !
  real(fp) :: fin(0:3,inf2,ny)             ! interpolant data (cf "herm2ev")
  !
  real(fp) :: fval(ivecd,*)                ! output returned
  !
  !  for detailed description of fin, ict and fval see subroutine
  !  herm2ev comments.  Note ict is not vectorized; the same output
  !  is expected to be returned for all input vector data points.
  !
  !  note that the index inputs ii,jj and parameter inputs
  !     xparam,yparam,hx,hxi,hy,hyi are vectorized, and the
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
     !   ...derivatives...
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
     iadr=0
     !
     !  get desired values:
     !
     if(ict(1).eq.1) then
        !
        !  function value:
        !
        iadr=iadr+1
        sum=axbar*(aybar*fin(0,i,j)  +ay*fin(0,i,j+1))+ &
             ax*(aybar*fin(0,i+1,j)+ay*fin(0,i+1,j+1))
        !
        sum=sum+hx(v)*( &
             bxbar*(aybar*fin(1,i,j)  +ay*fin(1,i,j+1))+ &
             bx*(aybar*fin(1,i+1,j)+ay*fin(1,i+1,j+1)))
        !
        sum=sum+hy(v)*( &
             axbar*(bybar*fin(2,i,j)  +by*fin(2,i,j+1))+ &
             ax*(bybar*fin(2,i+1,j)+by*fin(2,i+1,j+1)))
        !
        sum=sum+hx(v)*hy(v)*( &
             bxbar*(bybar*fin(3,i,j)  +by*fin(3,i,j+1))+ &
             bx*(bybar*fin(3,i+1,j)+by*fin(3,i+1,j+1)))
        !
        fval(v,iadr)=sum
     end if
     !
     if(ict(2).eq.1) then
        !
        !  df/dx:
        !
        iadr=iadr+1
        !
        sum=hxi(v)*( &
             axbarp*(aybar*fin(0,i,j)  +ay*fin(0,i,j+1))+ &
             axp*(aybar*fin(0,i+1,j)+ay*fin(0,i+1,j+1)))
        !
        sum=sum+ &
             bxbarp*(aybar*fin(1,i,j)  +ay*fin(1,i,j+1))+ &
             bxp*(aybar*fin(1,i+1,j)+ay*fin(1,i+1,j+1))
        !
        sum=sum+hxi(v)*hy(v)*( &
             axbarp*(bybar*fin(2,i,j)  +by*fin(2,i,j+1))+ &
             axp*(bybar*fin(2,i+1,j)+by*fin(2,i+1,j+1)))
        !
        sum=sum+hy(v)*( &
             bxbarp*(bybar*fin(3,i,j)  +by*fin(3,i,j+1))+ &
             bxp*(bybar*fin(3,i+1,j)+by*fin(3,i+1,j+1)))
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
             axbar*(aybarp*fin(0,i,j)  +ayp*fin(0,i,j+1))+ &
             ax*(aybarp*fin(0,i+1,j)+ayp*fin(0,i+1,j+1)))
        !
        sum=sum+hx(v)*hyi(v)*( &
             bxbar*(aybarp*fin(1,i,j)  +ayp*fin(1,i,j+1))+ &
             bx*(aybarp*fin(1,i+1,j)+ayp*fin(1,i+1,j+1)))
        !
        sum=sum+ &
             axbar*(bybarp*fin(2,i,j)  +byp*fin(2,i,j+1))+ &
             ax*(bybarp*fin(2,i+1,j)+byp*fin(2,i+1,j+1))
        !
        sum=sum+hx(v)*( &
             bxbar*(bybarp*fin(3,i,j)  +byp*fin(3,i,j+1))+ &
             bx*(bybarp*fin(3,i+1,j)+byp*fin(3,i+1,j+1)))
        !
        fval(v,iadr)=sum
     end if
     !
     if(ict(4).eq.1) then
        !
        !  d2f/dxdy:
        !
        iadr=iadr+1
        !
        sum=hxi(v)*hyi(v)*( &
             axbarp*(aybarp*fin(0,i,j)  +ayp*fin(0,i,j+1))+ &
             axp*(aybarp*fin(0,i+1,j)+ayp*fin(0,i+1,j+1)))
        !
        sum=sum+hyi(v)*( &
             bxbarp*(aybarp*fin(1,i,j)  +ayp*fin(1,i,j+1))+ &
             bxp*(aybarp*fin(1,i+1,j)+ayp*fin(1,i+1,j+1)))
        !
        sum=sum+hxi(v)*( &
             axbarp*(bybarp*fin(2,i,j)  +byp*fin(2,i,j+1))+ &
             axp*(bybarp*fin(2,i+1,j)+byp*fin(2,i+1,j+1)))
        !
        sum=sum+ &
             bxbarp*(bybarp*fin(3,i,j)  +byp*fin(3,i,j+1))+ &
             bxp*(bybarp*fin(3,i+1,j)+byp*fin(3,i+1,j+1))
        !
        fval(v,iadr)=sum
     end if
     !
  end do                             ! vector loop
  !
  return
end subroutine herm2fcn
!---------------------------------------------------------------------
!  evaluate C1 cubic Hermite function interpolation -- 2d fcn
!   --vectorized-- dmc 10 Feb 1999
!    --optimized for VARIATION along x axis ONLY--
!
subroutine herm2fcnx(ict,ivec,ivecd, &
     fval,ii,jj,xparam,yparam,hx,hxi,hy,hyi, &
     fin,inf2,ny)
  use precision_mod, only: fp
  !
  !============
  implicit none
  integer ny,inf2,j,i,iadr
  !============
  real(fp) :: yp,ypi,yp2,ypi2,ay,aybar,by,bybar,ayp,aybarp,byp
  real(fp) :: bybarp,xp,xpi,xp2,xpi2,ax,axbar,bx,bxbar,axp,axbarp
  real(fp) :: bxp,bxbarp
  !============
  integer ict(4)                    ! requested output control
  integer ivec                      ! vector length
  integer ivecd                     ! vector dimension (1st dim of fval)
  !
  integer ii(ivec),jj               ! target cells (i,j)
  real(fp) :: xparam(ivec),yparam
  ! normalized displacements from (i,j) corners
  !
  real(fp) :: hx(ivec),hy                  ! grid spacing, and
  real(fp) :: hxi(ivec),hyi                ! inverse grid spacing 1/(x(i+1)-x(i))
  ! & 1/(y(j+1)-y(j))
  !
  real(fp) :: fin(0:3,inf2,ny)             ! interpolant data (cf "herm2ev")
  !
  real(fp) :: fval(ivecd,*)                ! output returned
  !
  !  for detailed description of fin, ict and fval see subroutine
  !  herm2ev comments.  Note ict is not vectorized; the same output
  !  is expected to be returned for all input vector data points.
  !
  !  note that the index inputs ii,jj and parameter inputs
  !     xparam,yparam,hx,hxi,hy,hyi are vectorized, and the
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
  !   ...derivatives...
  !
  ayp=6.0_fp*yp*ypi
  aybarp=-ayp
  byp=yp*(3.0_fp*yp-2.0_fp)
  bybarp=ypi*(3.0_fp*ypi-2.0_fp)
  !
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
     !   ...derivatives...
     !
     axp=6.0_fp*xp*xpi
     axbarp=-axp
     bxp=xp*(3.0_fp*xp-2.0_fp)
     bxbarp=xpi*(3.0_fp*xpi-2.0_fp)
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
        sum=axbar*(aybar*fin(0,i,j)  +ay*fin(0,i,j+1))+ &
             ax*(aybar*fin(0,i+1,j)+ay*fin(0,i+1,j+1))
        !
        sum=sum+hx(v)*( &
             bxbar*(aybar*fin(1,i,j)  +ay*fin(1,i,j+1))+ &
             bx*(aybar*fin(1,i+1,j)+ay*fin(1,i+1,j+1)))
        !
        sum=sum+hy*( &
             axbar*(bybar*fin(2,i,j)  +by*fin(2,i,j+1))+ &
             ax*(bybar*fin(2,i+1,j)+by*fin(2,i+1,j+1)))
        !
        sum=sum+hx(v)*hy*( &
             bxbar*(bybar*fin(3,i,j)  +by*fin(3,i,j+1))+ &
             bx*(bybar*fin(3,i+1,j)+by*fin(3,i+1,j+1)))
        !
        fval(v,iadr)=sum
     end if
     !
     if(ict(2).eq.1) then
        !
        !  df/dx:
        !
        iadr=iadr+1
        !
        sum=hxi(v)*( &
             axbarp*(aybar*fin(0,i,j)  +ay*fin(0,i,j+1))+ &
             axp*(aybar*fin(0,i+1,j)+ay*fin(0,i+1,j+1)))
        !
        sum=sum+ &
             bxbarp*(aybar*fin(1,i,j)  +ay*fin(1,i,j+1))+ &
             bxp*(aybar*fin(1,i+1,j)+ay*fin(1,i+1,j+1))
        !
        sum=sum+hxi(v)*hy*( &
             axbarp*(bybar*fin(2,i,j)  +by*fin(2,i,j+1))+ &
             axp*(bybar*fin(2,i+1,j)+by*fin(2,i+1,j+1)))
        !
        sum=sum+hy*( &
             bxbarp*(bybar*fin(3,i,j)  +by*fin(3,i,j+1))+ &
             bxp*(bybar*fin(3,i+1,j)+by*fin(3,i+1,j+1)))
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
             axbar*(aybarp*fin(0,i,j)  +ayp*fin(0,i,j+1))+ &
             ax*(aybarp*fin(0,i+1,j)+ayp*fin(0,i+1,j+1)))
        !
        sum=sum+hx(v)*hyi*( &
             bxbar*(aybarp*fin(1,i,j)  +ayp*fin(1,i,j+1))+ &
             bx*(aybarp*fin(1,i+1,j)+ayp*fin(1,i+1,j+1)))
        !
        sum=sum+ &
             axbar*(bybarp*fin(2,i,j)  +byp*fin(2,i,j+1))+ &
             ax*(bybarp*fin(2,i+1,j)+byp*fin(2,i+1,j+1))
        !
        sum=sum+hx(v)*( &
             bxbar*(bybarp*fin(3,i,j)  +byp*fin(3,i,j+1))+ &
             bx*(bybarp*fin(3,i+1,j)+byp*fin(3,i+1,j+1)))
        !
        fval(v,iadr)=sum
     end if
     !
     if(ict(4).eq.1) then
        !
        !  d2f/dxdy:
        !
        iadr=iadr+1
        !
        sum=hxi(v)*hyi*( &
             axbarp*(aybarp*fin(0,i,j)  +ayp*fin(0,i,j+1))+ &
             axp*(aybarp*fin(0,i+1,j)+ayp*fin(0,i+1,j+1)))
        !
        sum=sum+hyi*( &
             bxbarp*(aybarp*fin(1,i,j)  +ayp*fin(1,i,j+1))+ &
             bxp*(aybarp*fin(1,i+1,j)+ayp*fin(1,i+1,j+1)))
        !
        sum=sum+hxi(v)*( &
             axbarp*(bybarp*fin(2,i,j)  +byp*fin(2,i,j+1))+ &
             axp*(bybarp*fin(2,i+1,j)+byp*fin(2,i+1,j+1)))
        !
        sum=sum+ &
             bxbarp*(bybarp*fin(3,i,j)  +byp*fin(3,i,j+1))+ &
             bxp*(bybarp*fin(3,i+1,j)+byp*fin(3,i+1,j+1))
        !
        fval(v,iadr)=sum
     end if
     !
  end do                             ! vector loop
  !
  return
end subroutine herm2fcnx
