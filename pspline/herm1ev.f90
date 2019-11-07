subroutine herm1ev(xget,x,nx,ilinx,f,ict,fval,ier)
  use psp_precision_mod, only: fp
  !
  !  evaluate a 1d cubic Hermite interpolant -- this is C1.
  !
  !  this subroutine calls two subroutines:
  !     herm1x   -- find cell containing (xget,yget)
  !     herm1fcn -- evaluate interpolant function and (optionally) derivatives
  !
  !  input arguments:
  !  ================
  !
  !============
  implicit none
  integer nx
  !============
  real(fp) :: xget                         ! target of this interpolation
  real(fp) :: x(nx)                        ! ordered x grid
  integer ilinx                     ! ilinx=1 => assume x evenly spaced
  !
  real(fp) :: f(0:1,nx)                    ! function data
  !
  !       contents of f:
  !
  !  f(0,i) = f @ x(i)
  !  f(1,i) = df/dx @ x(i)
  !
  integer ict(2)                    ! code specifying output desired
  !
  !  ict(1)=1 -- return f  (0, don't)
  !  ict(2)=1 -- return df/dx  (0, don't)
  !
  ! output arguments:
  ! =================
  !
  real(fp) :: fval(*)                      ! output data
  integer ier                       ! error code =0 ==> no error
  !
  !  fval(1) receives the first output (depends on ict(...) spec)
  !  fval(2) receives the second output (depends on ict(...) spec)
  !
  !  examples:
  !    on input ict = [1,1]
  !   on output fval= [f,df/dx]
  !
  !    on input ict = [1,0]
  !   on output fval= [f] ... element 2 never referenced
  !
  !    on input ict = [0,1]
  !   on output fval= [df/dx] ... element 2 never referenced
  !
  !  ier -- completion code:  0 means OK
  !-------------------
  !  local:
  !
  integer, dimension(1) :: i   ! cell index
  !
  !  normalized displacement from (x(i)) corner of cell.
  !    xparam=0 @x(i)  xparam=1 @x(i+1)
  !
  real(fp), dimension(1) :: xparam
  !
  !  cell dimensions and
  !  inverse cell dimensions hxi = 1/(x(i+1)-x(i))
  !
  real(fp), dimension(1) :: hx
  real(fp), dimension(1) :: hxi
  !
  !  0 .le. xparam .le. 1
  !
  !---------------------------------------------------------------------
  !
  call herm1x(xget,x,nx,ilinx,i(1),xparam(1),hx(1),hxi(1),ier)
  if(ier.ne.0) return
  !
  call herm1fcn(ict,1,1,fval,i,xparam,hx,hxi,f,nx)
  !
  return
end subroutine herm1ev
!---------------------------------------------------------------------
!  herm1xy -- look up x zone
!
!  this is the "first part" of herm1ev, see comments, above.
!
subroutine herm1x(xget,x,nx,ilinx,i,xparam,hx,hxi,ier)
  use psp_precision_mod, only: fp
  !
  !  input of herm1x
  !  ===============
  !
  !============
  implicit none
  integer nxm,ii
  !============
  real(fp) :: zxget,zxtol
  !============
  integer nx                        ! x array dimension
  !
  real(fp) :: xget                         ! target point
  real(fp) :: x(nx)                        ! indep. coordinate, strict ascending
  !
  integer ilinx                     ! =1:  x evenly spaced
  !
  !  output of herm1x
  !  =================
  integer i                         ! index to cell containing target pt
  !          on exit:  1.le.i.le.nx-1
  !
  !  normalized position w/in (i) cell (always btw 0 and 1):
  !
  real(fp) :: xparam                       ! (xget-x(i))/(x(i+1)-x(i))
  !
  !  grid spacing
  !
  real(fp) :: hx                           ! hx = x(i+1)-x(i)
  !
  !  inverse grid spacing:
  !
  real(fp) :: hxi                          ! 1/hx = 1/(x(i+1)-x(i))
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
  if((xget.lt.x(1)).or.(xget.gt.x(nx))) then
     zxtol=4.0E-7_fp*max(abs(x(1)),abs(x(nx)))
     if((xget.lt.x(1)-zxtol).or.(xget.gt.x(nx)+zxtol)) then
        ier=1
        write(6,1001) xget,x(1),x(nx)
1001    format(' ?herm1ev:  xget=',1pe11.4,' out of range ', &
             1pe11.4,' to ',1pe11.4)
     else
        if((xget.lt.x(1)-0.5_fp*zxtol).or. &
             (xget.gt.x(nx)+0.5_fp*zxtol)) &
             write(6,1011) xget,x(1),x(nx)
1011    format(' %herm1ev:  xget=',1pe15.8,' beyond range ', &
             1pe15.8,' to ',1pe15.8,' (fixup applied)')
        if(xget.lt.x(1)) then
           zxget=x(1)
        else
           zxget=x(nx)
        end if
     end if
  end if
  if(ier.ne.0) return
  !
  !  now find interval in which target point lies..
  !
  nxm=nx-1
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
  hx=(x(i+1)-x(i))
  !
  hxi=1.0_fp/hx
  !
  xparam=(zxget-x(i))*hxi
  !
  return
end subroutine herm1x
!---------------------------------------------------------------------
!  evaluate C1 cubic Hermite function interpolation -- 1d fcn
!   --vectorized-- dmc 10 Feb 1999
!
subroutine herm1fcn(ict,ivec,ivecd, &
     fval,ii,xparam,hx,hxi,fin,nx)
  use psp_precision_mod, only: fp
  !
  !  input:
  !
  !============
  implicit none
  integer nx,i,iadr
  !============
  real(fp) :: xp,xpi,xp2,xpi2,ax,axbar,bx,bxbar,axp,axbarp,bxp,bxbarp
  !============
  integer ict(2)                    ! requested output control
  integer ivec                      ! vector length
  integer ivecd                     ! vector dimension (1st dim of fval)
  !
  integer ii(ivec)                  ! target cells
  real(fp) :: xparam(ivec)
  ! normalized displacements in cells
  !
  real(fp) :: hx(ivec)                     ! grid spacing, and
  real(fp) :: hxi(ivec)                    ! inverse grid spacing 1/(x(i+1)-x(i))
  !
  real(fp) :: fin(0:1,nx)                  ! Hermite dataset
  !
  !  output:
  !
  real(fp) :: fval(ivecd,*)                ! interpolation results
  !
  !  for detailed description of fin, ict and fval see subroutine
  !  herm1ev comments.  Note ict is not vectorized -- the same output
  !  is expected to be returned for all input vector data points.
  !
  !  note that the index inputs ii,jj and parameter inputs
  !     xparam,yparam,hx,hxi,hy,hyi are vectorized, and the
  !
  !     output array fval has a vector ** 1st dimension ** whose
  !     size must be given as a separate argument; ivecd.ge.ivec
  !     expected!
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
  !   ...in x direction
  !
  real(fp) :: sum
  integer v
  !
  do v=1,ivec
     i=ii(v)
     xp=xparam(v)
     xpi=1.0_fp-xp
     xp2=xp*xp
     xpi2=xpi*xpi
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
        !
        !  Hermite basis functions
        !
        ax=xp2*(3.0_fp-2.0_fp*xp)
        axbar=1.0_fp-ax
        bx=-xp2*xpi
        bxbar=xpi2*xp
        !
        sum=axbar*fin(0,i) + ax*fin(0,i+1)
        !
        sum=sum+hx(v)*(bxbar*fin(1,i) + bx*fin(1,i+1))
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
        axp=6.0_fp*xp*xpi
        axbarp=-axp
        bxp=xp*(3.0_fp*xp-2.0_fp)
        bxbarp=xpi*(3.0_fp*xpi-2.0_fp)
        !
        sum=hxi(v)*(axbarp*fin(0,i) +axp*fin(0,i+1))
        !
        sum=sum+ bxbarp*fin(1,i) + bxp*fin(1,i+1)
        !
        fval(v,iadr)=sum
     end if
     !
  end do                             ! vector loop
  !
  return
end subroutine herm1fcn
