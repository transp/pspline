!     tcspeval -- eval tricubic spline function and/or derivatives
!
subroutine tcspeval(xget,yget,zget,iselect,fval, &
     x,nx,y,ny,z,nz,ilinx,iliny,ilinz,f,inf4,inf5, &
     ier)
  use psp_precision_mod, only: fp
  !
  implicit none
  integer iselect(10)
  integer ilinx,iliny,ilinz,nx,ny,nz,inf4,inf5,ier
  !
  real(fp) :: xget,yget,zget
  real(fp) :: fval(*)
  real(fp) :: x(nx),y(ny),z(nz),f(4,4,4,inf4,inf5,nz)
  !
  !  modification -- dmc 11 Jan 1999 -- remove SAVE stmts; break routine
  !    into these parts:
  !
  !    tcspevxyz -- find grid cell of target pt.
  !    tcspevfn -- evaluate function using output of tcspevxyz
  !
  !    in cases where multiple functions are defined on the same grid,
  !    time can be saved by using tcspevxyz once and then tcspevfn
  !    multiple times.
  !
  !  input:
  !     (xget,yget,zget)   location where interpolated value is desired
  !                   x(1).le.xget.le.x(nx) expected
  !                   y(1).le.yget.le.y(ny) expected
  !                   z(1).le.zget.le.z(nz) expected
  !
  !     iselect       select desired output
  !
  !                     iselect(1)=1 -- want function value (f) itself
  !                     iselect(2)=1 -- want  df/dx
  !                     iselect(3)=1 -- want  df/dy
  !                     iselect(4)=1 -- want  df/dz
  !                     iselect(5)=1 -- want  d2f/dx2
  !                     iselect(6)=1 -- want  d2f/dy2
  !                     iselect(7)=1 -- want  d2f/dz2
  !                     iselect(8)=1 -- want  d2f/dxdy
  !                     iselect(9)=1 -- want  d2f/dxdz
  !                     iselect(10)=1 -- want  d2f/dydz
  !
  !
  !              example:  iselect(1)=iselect(2)=iselect(3)=iselect(4)=1
  !                            f, df/dx, df/dy, and df/dz all evaluated
  !                        iselect(5)=iselect(6)=iselect(7)=0
  !                        iselect(8)=iselect(9)=iselect(10)=0
  !                            2nd derivatives not evaluated.
  !
  !  (new dmc Dec 2005 -- higher derivatives available)
  !    iselect(1)=3 --> 3rd derivative, .le.2 diff. in any coordinate
  !      iselect(2:8) select: fxxy fxxz fxyy fxyz fxzz fyyz fyzz
  !      ->note iselect(1)=3, iselect(5)=1 gives fxyz = d3f/dxdydz
  !    iselect(1)=-3 --> 3rd derivative, 3 in one coordinate
  !      iselect(2:4) select: fxxx fyyy fzzz
  !    iselect(1)=4 --> 3rd derivative, .le.2 diff. in any coordinate
  !      iselect(2:7) select: fxxyy fxxyz fxxzz fxyyz fxyzz fyyzz
  !    iselect(1)=-4 --> 3rd derivative, 3 in one coordinate
  !      iselect(2:7) select: fxxxy fxxxz fxyyy fxzzz fyyyz fyzzz
  !    iselect(1)=5 --> 3rd derivative, .le.2 diff. in any coordinate
  !      iselect(2:4) select: fxxyyz fxxyzz fxyyzz
  !    iselect(1)=-5 --> 3rd derivative, 3 in one coordinate
  !      iselect(2:10) select:  fxxxyy fxxxyz fxxxzz fxxyyy fxxzzz
  !                             fxyyyz fxyzzz fyyyzz fzzzyy
  !    iselect(1)=6 --> 3rd derivative, .le.2 diff. in any coordinate
  !      fxxyyzz
  !    iselect(1)=-6 --> 3rd derivative, 3 in one coordinate
  !      iselect(2:10) select: fxxxyyy fxxxyyz fxxxyzz fxxxyyz
  !                            fxxyyyz fxxyzzz fxyyyzz fxyyzzz fyyyzzz
  !    iselect(1)=-7 --> 7th derivative
  !      iselect(2:7) select: fxxxyyyz fxxxyyzz fxxxyzzz
  !                           fxxyyyzz fxxyyzzz fxyyyzzz
  !    iselect(1)=-8 --> 8th derivative
  !      iselect(2:4) select: fxxxyyyzz fxxxyyzzz fxxyyyzzz
  !    iselect(1)=-9 --> 9th derivative:  fxxxyyyzzz
  !
  !-------
  !
  !     x(1...nx)     independent coordinate x, strict ascending
  !     y(1...ny)     independent coordinate y, strict ascending
  !     z(1...nz)     independent coordinate y, strict ascending
  !
  !     ilinx  --  =1: flag that x is linearly spaced
  !
  !                   see fval (output) description.
  !
  !     x(1...nx)     independent coordinate x, strict ascending
  !     y(1...ny)     independent coordinate y, strict ascending
  !     z(1...nz)     independent coordinate y, strict ascending
  !
  !     ilinx  --  =1: flag that x is linearly spaced (avoid search for speed)
  !     iliny  --  =1: flag that y is linearly spaced (avoid search for speed)
  !     ilinz  --  =1: flat that z is linearly spaced (avoid search for speed)
  !
  !  **CAUTION** actual even spacing of x, y, z is NOT CHECKED HERE!
  !
  !
  !     f             the function values (at grid points) and spline coefs
  !
  !  evaluation formula:  for point x btw x(i) and x(i+1), dx=x-x(i)
  !                             and y btw y(j) and y(j+1), dy=y-y(j),
  !                             and z btw z(k) and z(k+1), dz=z-z(k)
  !
  !  do m=1,4
  !   p(m) =
  !    f(1,1,m,i,j,k)+dx*f(2,1,m,i,j,k)+dx**2*f(3,1,m,i,j,k)+dx**3*f(4,1,m,i,j,k)
  !   +dy*(
  !   f(1,2,m,i,j,k)+dx*f(2,2,m,i,j,k)+dx**2*f(3,2,m,i,j,k)+dx**3*f(4,2,m,i,j,k))
  !   +dy**2*(
  !   f(1,3,m,i,j,k)+dx*f(2,3,m,i,j,k)+dx**2*f(3,3,m,i,j,k)+dx**3*f(4,3,m,i,j,k))
  !   +dy**3*(
  !   f(1,4,m,i,j,k)+dx*f(2,4,m,i,j,k)+dx**2*f(3,4,m,i,j,k)+dx**3*f(4,4,m,i,j,k))
  !  end do
  !  answer = p(1)+dz*p(2)+dz**2*p(3)+dz**3*p(4)
  !
  !      where d2=dy**2 and d3=dy**3.
  !
  !  nb dmc Feb 1999 -- p loops unrolled, by hand, to aid vector compilers
  !
  !  output:
  !      up to 10 elements of fval, ordered as follows:
  !        fval(1)=function value or lowest order derivative requested
  !        fval(2)=next order derivative
  !             etc
  !        the ordering is a subset of the sequence given under the "iselect"
  !        description; the first M elements of fval are used, where M = the
  !        number of non-zero elements of iselect.
  !
  !      ier = 0 -- successful completion; = 1 -- an error occurred.
  !
  !-------------------------------------------------------------------
  !  local
  !
  integer :: i(1),j(1),k(1)
  !
  real(fp) :: dx(1),dy(1),dz(1)
  !
  !--------------------------
  !
  i(1)=0
  j(1)=0
  k(1)=0
  !
  call tcspevxyz(xget,yget,zget,x,nx,y,ny,z,nz,ilinx,iliny,ilinz, &
       i(1),j(1),k(1),dx(1),dy(1),dz(1),ier)
  if(ier.ne.0) return
  !
  call tcspevfn(iselect,1,1,fval,i,j,k,dx,dy,dz,f,inf4,inf5,nz)
  !
  return
end subroutine tcspeval
!
!-------------------------------------------------------------------------
!  tcspevxyz -- look up x-y zone
!
!  this is the "first part" of tcspeval, see comments, above.
!
subroutine tcspevxyz(xget,yget,zget,x,nx,y,ny,z,nz, &
     ilinx,iliny,ilinz, &
     i,j,k,dx,dy,dz,ier)
  use psp_precision_mod, only: fp
  !
  !============
  implicit none
  integer nxm,nym,nzm,ii,jj,kk
  !============
  real(fp) :: zxget,zyget,zzget,zxtol,zytol,zztol
  !============
  integer nx,ny,nz                  ! array dimensions
  !
  real(fp) :: xget,yget,zget               ! target point
  real(fp) :: x(nx),y(ny),z(nz)            ! indep. coords.
  !
  integer ilinx                     ! =1:  assume x evenly spaced
  integer iliny                     ! =1:  assume y evenly spaced
  integer ilinz                     ! =1:  assume z evenly spaced
  !
  !  output of tcspevxyz
  !
  integer i,j,k                     ! index to cell containing target pt
  real(fp) :: dx,dy,dz                     ! displacement of target pt w/in cell
  ! dx=x-x(i)  dy=y-y(j)  dz=z-z(k)
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
  !
  if((xget.lt.x(1)).or.(xget.gt.x(nx))) then
     zxtol=4.0E-7_fp*max(abs(x(1)),abs(x(nx)))
     if((xget.lt.x(1)-zxtol).or.(xget.gt.x(nx)+zxtol)) then
        ier=1
        write(6,1001) xget,x(1),x(nx)
1001    format(' ?tcspeval:  xget=',1pe11.4,' out of range ', &
             1pe11.4,' to ',1pe11.4)
     else
        if((xget.lt.x(1)-0.5_fp*zxtol).or. &
             (xget.gt.x(nx)+0.5_fp*zxtol)) &
             write(6,1011) xget,x(1),x(nx)
1011    format(' %tcspeval:  xget=',1pe15.8,' beyond range ', &
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
1002    format(' ?tcspeval:  yget=',1pe11.4,' out of range ', &
             1pe11.4,' to ',1pe11.4)
     else
        if((yget.lt.y(1)-0.5_fp*zytol).or. &
             (yget.gt.y(ny)+0.5_fp*zytol)) &
             write(6,1012) yget,y(1),y(ny)
1012    format(' %tcspeval:  yget=',1pe15.8,' beyond range ', &
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
1003    format(' ?tcspeval:  zget=',1pe11.4,' out of range ', &
             1pe11.4,' to ',1pe11.4)
     else
        if((zget.lt.z(1)-0.5_fp*zztol).or. &
             (zget.gt.z(nz)+0.5_fp*zztol)) &
             write(6,1013) zget,z(1),z(nz)
1013    format(' %tcspeval:  zget=',1pe15.8,' beyond range ', &
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
  dx=zxget-x(i)
  dy=zyget-y(j)
  dz=zzget-z(k)
  !
  return
end subroutine tcspevxyz
!------------------------------------------------------------------------
!  tcspevfn -- OK now evaluate the tricubic spline
!     evaluate at a vector set of target locations as specified by
!     input vectors (iv,jv,kv), (dxv,dyv,dzv)
!
subroutine tcspevfn(ict,ivec,ivd,fval,iv,jv,kv,dxv,dyv,dzv, &
     f,inf4,inf5,nz)
  use psp_precision_mod, only: fp
  !
  !  input:
  !============
  implicit none
  integer nz,iaval,i,j,k
  !============
  real(fp) :: dx,dy,dz,p1,p2,p3,p4
  !============
  integer ict(10)                   ! selector:
  !        ict(1)=1 for f      (don't evaluate f if ict(1)=0)
  !        ict(2)=1 for df/dx   ""
  !        ict(3)=1 for df/dy   ""
  !        ict(4)=1 for df/dz   ""
  !        ict(5)=1 for d2f/dx2
  !        ict(6)=1 for d2f/dy2
  !        ict(7)=1 for d2f/dz2
  !        ict(8)=1 for d2f/dxdy
  !        ict(9)=1 for d2f/dxdz
  !        ict(10)=1 for d2f/dydz
  !
  integer ivec,ivd                  ! vector dimensioning
  !
  !    ivec-- number of vector pts (spline values to look up)
  !    ivd -- 1st dimension of fval, .ge. ivec
  !
  ! output:
  real(fp) :: fval(ivd,*)                     ! output array
  !
  !  for vector elements v,  (iv(v),jv(v),kv(v),dxv(v),dyv(v),dzv(v))
  !
  !  fval(v,1) = first item requested by ict(...),
  !  fval(v,2) = 2nd item requested,  ...etc...
  !
  !  input:
  integer iv(ivec),jv(ivec),kv(ivec) ! grid cell indices
  real(fp) :: dxv(ivec),dyv(ivec),dzv(ivec) ! displacements w/in cell
  !
  integer inf4                      ! 4th dimension of f, .le.nx
  integer inf5                      ! 5th dimension of f, .le.ny
  real(fp) :: f(4,4,4,inf4,inf5,nz)        ! tricubic fcn spline coeffs array
  !
  !  usage example:
  !
  !  1.  for each element (xx(v),yy(v),zz(v)) in a vector of (x,y,z)
  !    triples, find the x,y,z zone indices and displacements with
  !    to the "lower left corner" of the zone; store these in vectors
  !    iv,jv,kv and dxv,dyv,dzv
  !
  !  2.  set ict(1)=0, ict(2)=1, ict(3)=1, ict(4)=1 & the rest zero --
  !      to get only the 1st derivatives.
  !
  !  3.  ivec is the length of the vector; ivd is the 1st dimension
  !      of the array fval to receive the output
  !
  !      real fval(ivd,10)
  !      real xv(ivd),yv(ivd),zv(ivd)
  !      integer iv(ivd),jv(ivd),kv(ivd)
  !      real dxv(ivd),dyv(ivd),dzv(ivd)
  !      integer ict(10)
  !
  !      real fspline(4,4,4,nx,ny,nz)  ! spline coeffs
  !      data ict/0,1,1,1,0,0,0,0,0,0/    ! this call:  want 1st
  !                               ! derivatives only ... these will
  !                               ! be output to
  !                               ! fval(*,1) fval(*,2) fval(*,3)
  !      ...
  !      do iv=1,ivec
  !        ...                    ! find indices and displacements
  !      end do
  !      call tcspevfn(ict,ivec,ivd,fval,iv,jv,kv,dxv,dyv,dzv, &
  !                    fspline,nx,ny,nz)
  !
  !-------------------------------------------------------------------
  !
  !  local --
  !
  !cc      real p(4)  use p1,p2,p3,p4 now
  !
  integer v                         ! vector index
  !
  !  OK can now do evaluations
  !
  iaval=0  ! fval addressing
  !
  if(abs(ict(1)).le.2) then
     if((ict(1).gt.0).or.(ict(1).eq.-1)) then
        !  evaluate f
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           k=kv(v)
           dx=dxv(v)
           dy=dyv(v)
           dz=dzv(v)
           p1= &
                f(1,1,1,i,j,k)+dy*(f(1,2,1,i,j,k)+dy*(f(1,3,1,i,j,k)+ &
                dy*f(1,4,1,i,j,k))) &
                +dx*(f(2,1,1,i,j,k)+dy*(f(2,2,1,i,j,k)+dy*(f(2,3,1,i,j,k)+ &
                dy*f(2,4,1,i,j,k))) &
                +dx*(f(3,1,1,i,j,k)+dy*(f(3,2,1,i,j,k)+dy*(f(3,3,1,i,j,k)+ &
                dy*f(3,4,1,i,j,k))) &
                +dx*(f(4,1,1,i,j,k)+dy*(f(4,2,1,i,j,k)+dy*(f(4,3,1,i,j,k)+ &
                dy*f(4,4,1,i,j,k))) &
                )))
           p2= &
                f(1,1,2,i,j,k)+dy*(f(1,2,2,i,j,k)+dy*(f(1,3,2,i,j,k)+ &
                dy*f(1,4,2,i,j,k))) &
                +dx*(f(2,1,2,i,j,k)+dy*(f(2,2,2,i,j,k)+dy*(f(2,3,2,i,j,k)+ &
                dy*f(2,4,2,i,j,k))) &
                +dx*(f(3,1,2,i,j,k)+dy*(f(3,2,2,i,j,k)+dy*(f(3,3,2,i,j,k)+ &
                dy*f(3,4,2,i,j,k))) &
                +dx*(f(4,1,2,i,j,k)+dy*(f(4,2,2,i,j,k)+dy*(f(4,3,2,i,j,k)+ &
                dy*f(4,4,2,i,j,k))) &
                )))
           p3= &
                f(1,1,3,i,j,k)+dy*(f(1,2,3,i,j,k)+dy*(f(1,3,3,i,j,k)+ &
                dy*f(1,4,3,i,j,k))) &
                +dx*(f(2,1,3,i,j,k)+dy*(f(2,2,3,i,j,k)+dy*(f(2,3,3,i,j,k)+ &
                dy*f(2,4,3,i,j,k))) &
                +dx*(f(3,1,3,i,j,k)+dy*(f(3,2,3,i,j,k)+dy*(f(3,3,3,i,j,k)+ &
                dy*f(3,4,3,i,j,k))) &
                +dx*(f(4,1,3,i,j,k)+dy*(f(4,2,3,i,j,k)+dy*(f(4,3,3,i,j,k)+ &
                dy*f(4,4,3,i,j,k))) &
                )))
           p4= &
                f(1,1,4,i,j,k)+dy*(f(1,2,4,i,j,k)+dy*(f(1,3,4,i,j,k)+ &
                dy*f(1,4,4,i,j,k))) &
                +dx*(f(2,1,4,i,j,k)+dy*(f(2,2,4,i,j,k)+dy*(f(2,3,4,i,j,k)+ &
                dy*f(2,4,4,i,j,k))) &
                +dx*(f(3,1,4,i,j,k)+dy*(f(3,2,4,i,j,k)+dy*(f(3,3,4,i,j,k)+ &
                dy*f(3,4,4,i,j,k))) &
                +dx*(f(4,1,4,i,j,k)+dy*(f(4,2,4,i,j,k)+dy*(f(4,3,4,i,j,k)+ &
                dy*f(4,4,4,i,j,k))) &
                )))
           fval(v,iaval)=p1+dz*(p2+dz*(p3+dz*p4))
        end do
     end if
     !
     if((ict(2).gt.0).and.(ict(1).ne.-1)) then
        !  evaluate df/dx
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           k=kv(v)
           dx=dxv(v)
           dy=dyv(v)
           dz=dzv(v)
           p1= &
                f(2,1,1,i,j,k)+dy*(f(2,2,1,i,j,k)+dy*(f(2,3,1,i,j,k)+ &
                dy*f(2,4,1,i,j,k))) &
                +2._fp*dx* &
                (f(3,1,1,i,j,k)+dy*(f(3,2,1,i,j,k)+dy*(f(3,3,1,i,j,k)+ &
                dy*f(3,4,1,i,j,k))) &
                +1.5_fp*dx* &
                (f(4,1,1,i,j,k)+dy*(f(4,2,1,i,j,k)+dy*(f(4,3,1,i,j,k)+ &
                dy*f(4,4,1,i,j,k))) &
                ))
           p2= &
                f(2,1,2,i,j,k)+dy*(f(2,2,2,i,j,k)+dy*(f(2,3,2,i,j,k)+ &
                dy*f(2,4,2,i,j,k))) &
                +2._fp*dx* &
                (f(3,1,2,i,j,k)+dy*(f(3,2,2,i,j,k)+dy*(f(3,3,2,i,j,k)+ &
                dy*f(3,4,2,i,j,k))) &
                +1.5_fp*dx* &
                (f(4,1,2,i,j,k)+dy*(f(4,2,2,i,j,k)+dy*(f(4,3,2,i,j,k)+ &
                dy*f(4,4,2,i,j,k))) &
                ))
           p3= &
                f(2,1,3,i,j,k)+dy*(f(2,2,3,i,j,k)+dy*(f(2,3,3,i,j,k)+ &
                dy*f(2,4,3,i,j,k))) &
                +2._fp*dx* &
                (f(3,1,3,i,j,k)+dy*(f(3,2,3,i,j,k)+dy*(f(3,3,3,i,j,k)+ &
                dy*f(3,4,3,i,j,k))) &
                +1.5_fp*dx* &
                (f(4,1,3,i,j,k)+dy*(f(4,2,3,i,j,k)+dy*(f(4,3,3,i,j,k)+ &
                dy*f(4,4,3,i,j,k))) &
                ))
           p4= &
                f(2,1,4,i,j,k)+dy*(f(2,2,4,i,j,k)+dy*(f(2,3,4,i,j,k)+ &
                dy*f(2,4,4,i,j,k))) &
                +2._fp*dx* &
                (f(3,1,4,i,j,k)+dy*(f(3,2,4,i,j,k)+dy*(f(3,3,4,i,j,k)+ &
                dy*f(3,4,4,i,j,k))) &
                +1.5_fp*dx* &
                (f(4,1,4,i,j,k)+dy*(f(4,2,4,i,j,k)+dy*(f(4,3,4,i,j,k)+ &
                dy*f(4,4,4,i,j,k))) &
                ))
           fval(v,iaval)=p1+dz*(p2+dz*(p3+dz*p4))
        end do
     end if
     !
     if((ict(3).gt.0).and.(ict(1).ne.-1)) then
        !  evaluate df/dy
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           k=kv(v)
           dx=dxv(v)
           dy=dyv(v)
           dz=dzv(v)
           p1= &
                f(1,2,1,i,j,k)+ &
                dy*(2.0_fp*f(1,3,1,i,j,k)+dy*3.0_fp*f(1,4,1,i,j,k)) &
                +dx*(f(2,2,1,i,j,k)+ &
                dy*(2.0_fp*f(2,3,1,i,j,k)+dy*3.0_fp*f(2,4,1,i,j,k)) &
                +dx*(f(3,2,1,i,j,k)+ &
                dy*(2.0_fp*f(3,3,1,i,j,k)+dy*3.0_fp*f(3,4,1,i,j,k)) &
                +dx*(f(4,2,1,i,j,k)+ &
                dy*(2.0_fp*f(4,3,1,i,j,k)+dy*3.0_fp*f(4,4,1,i,j,k)) &
                )))
           p2= &
                f(1,2,2,i,j,k)+ &
                dy*(2.0_fp*f(1,3,2,i,j,k)+dy*3.0_fp*f(1,4,2,i,j,k)) &
                +dx*(f(2,2,2,i,j,k)+ &
                dy*(2.0_fp*f(2,3,2,i,j,k)+dy*3.0_fp*f(2,4,2,i,j,k)) &
                +dx*(f(3,2,2,i,j,k)+ &
                dy*(2.0_fp*f(3,3,2,i,j,k)+dy*3.0_fp*f(3,4,2,i,j,k)) &
                +dx*(f(4,2,2,i,j,k)+ &
                dy*(2.0_fp*f(4,3,2,i,j,k)+dy*3.0_fp*f(4,4,2,i,j,k)) &
                )))
           p3= &
                f(1,2,3,i,j,k)+ &
                dy*(2.0_fp*f(1,3,3,i,j,k)+dy*3.0_fp*f(1,4,3,i,j,k)) &
                +dx*(f(2,2,3,i,j,k)+ &
                dy*(2.0_fp*f(2,3,3,i,j,k)+dy*3.0_fp*f(2,4,3,i,j,k)) &
                +dx*(f(3,2,3,i,j,k)+ &
                dy*(2.0_fp*f(3,3,3,i,j,k)+dy*3.0_fp*f(3,4,3,i,j,k)) &
                +dx*(f(4,2,3,i,j,k)+ &
                dy*(2.0_fp*f(4,3,3,i,j,k)+dy*3.0_fp*f(4,4,3,i,j,k)) &
                )))
           p4= &
                f(1,2,4,i,j,k)+ &
                dy*(2.0_fp*f(1,3,4,i,j,k)+dy*3.0_fp*f(1,4,4,i,j,k)) &
                +dx*(f(2,2,4,i,j,k)+ &
                dy*(2.0_fp*f(2,3,4,i,j,k)+dy*3.0_fp*f(2,4,4,i,j,k)) &
                +dx*(f(3,2,4,i,j,k)+ &
                dy*(2.0_fp*f(3,3,4,i,j,k)+dy*3.0_fp*f(3,4,4,i,j,k)) &
                +dx*(f(4,2,4,i,j,k)+ &
                dy*(2.0_fp*f(4,3,4,i,j,k)+dy*3.0_fp*f(4,4,4,i,j,k)) &
                )))
           fval(v,iaval)=p1+dz*(p2+dz*(p3+dz*p4))
        end do
     end if
     !
     if((ict(4).gt.0).and.(ict(1).ne.-1)) then
        !  evaluate df/dz
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           k=kv(v)
           dx=dxv(v)
           dy=dyv(v)
           dz=dzv(v)
           p2= &
                f(1,1,2,i,j,k)+dy*(f(1,2,2,i,j,k)+dy*(f(1,3,2,i,j,k)+ &
                dy*f(1,4,2,i,j,k))) &
                +dx*(f(2,1,2,i,j,k)+dy*(f(2,2,2,i,j,k)+dy*(f(2,3,2,i,j,k)+ &
                dy*f(2,4,2,i,j,k))) &
                +dx*(f(3,1,2,i,j,k)+dy*(f(3,2,2,i,j,k)+dy*(f(3,3,2,i,j,k)+ &
                dy*f(3,4,2,i,j,k))) &
                +dx*(f(4,1,2,i,j,k)+dy*(f(4,2,2,i,j,k)+dy*(f(4,3,2,i,j,k)+ &
                dy*f(4,4,2,i,j,k))) &
                )))
           p3= &
                f(1,1,3,i,j,k)+dy*(f(1,2,3,i,j,k)+dy*(f(1,3,3,i,j,k)+ &
                dy*f(1,4,3,i,j,k))) &
                +dx*(f(2,1,3,i,j,k)+dy*(f(2,2,3,i,j,k)+dy*(f(2,3,3,i,j,k)+ &
                dy*f(2,4,3,i,j,k))) &
                +dx*(f(3,1,3,i,j,k)+dy*(f(3,2,3,i,j,k)+dy*(f(3,3,3,i,j,k)+ &
                dy*f(3,4,3,i,j,k))) &
                +dx*(f(4,1,3,i,j,k)+dy*(f(4,2,3,i,j,k)+dy*(f(4,3,3,i,j,k)+ &
                dy*f(4,4,3,i,j,k))) &
                )))
           p4= &
                f(1,1,4,i,j,k)+dy*(f(1,2,4,i,j,k)+dy*(f(1,3,4,i,j,k)+ &
                dy*f(1,4,4,i,j,k))) &
                +dx*(f(2,1,4,i,j,k)+dy*(f(2,2,4,i,j,k)+dy*(f(2,3,4,i,j,k)+ &
                dy*f(2,4,4,i,j,k))) &
                +dx*(f(3,1,4,i,j,k)+dy*(f(3,2,4,i,j,k)+dy*(f(3,3,4,i,j,k)+ &
                dy*f(3,4,4,i,j,k))) &
                +dx*(f(4,1,4,i,j,k)+dy*(f(4,2,4,i,j,k)+dy*(f(4,3,4,i,j,k)+ &
                dy*f(4,4,4,i,j,k))) &
                )))
           fval(v,iaval)=p2+dz*(2.0_fp*p3+dz*3.0_fp*p4)
        end do
     end if
     !
     if((ict(5).gt.0).or.(ict(1).eq.-1)) then
        !  evaluate d2f/dx2
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           k=kv(v)
           dx=dxv(v)
           dy=dyv(v)
           dz=dzv(v)
           p1= 2.0_fp* &
                (f(3,1,1,i,j,k)+dy*(f(3,2,1,i,j,k)+dy*(f(3,3,1,i,j,k)+ &
                dy*f(3,4,1,i,j,k)))) &
                +6.0_fp*dx* &
                (f(4,1,1,i,j,k)+dy*(f(4,2,1,i,j,k)+dy*(f(4,3,1,i,j,k)+ &
                dy*f(4,4,1,i,j,k))))
           p2= 2.0_fp* &
                (f(3,1,2,i,j,k)+dy*(f(3,2,2,i,j,k)+dy*(f(3,3,2,i,j,k)+ &
                dy*f(3,4,2,i,j,k)))) &
                +6.0_fp*dx* &
                (f(4,1,2,i,j,k)+dy*(f(4,2,2,i,j,k)+dy*(f(4,3,2,i,j,k)+ &
                dy*f(4,4,2,i,j,k))))
           p3= 2.0_fp* &
                (f(3,1,3,i,j,k)+dy*(f(3,2,3,i,j,k)+dy*(f(3,3,3,i,j,k)+ &
                dy*f(3,4,3,i,j,k)))) &
                +6.0_fp*dx* &
                (f(4,1,3,i,j,k)+dy*(f(4,2,3,i,j,k)+dy*(f(4,3,3,i,j,k)+ &
                dy*f(4,4,3,i,j,k))))
           p4= 2.0_fp* &
                (f(3,1,4,i,j,k)+dy*(f(3,2,4,i,j,k)+dy*(f(3,3,4,i,j,k)+ &
                dy*f(3,4,4,i,j,k)))) &
                +6.0_fp*dx* &
                (f(4,1,4,i,j,k)+dy*(f(4,2,4,i,j,k)+dy*(f(4,3,4,i,j,k)+ &
                dy*f(4,4,4,i,j,k))))
           fval(v,iaval)=p1+dz*(p2+dz*(p3+dz*p4))
        end do
     end if
     !
     if((ict(6).gt.0).or.(ict(1).eq.-1)) then
        !  evaluate d2f/dy2
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           k=kv(v)
           dx=dxv(v)
           dy=dyv(v)
           dz=dzv(v)
           p1= &
                2.0_fp*f(1,3,1,i,j,k)+6.0_fp*dy*f(1,4,1,i,j,k) &
                +dx*(2.0_fp*f(2,3,1,i,j,k)+6.0_fp*dy*f(2,4,1,i,j,k) &
                +dx*(2.0_fp*f(3,3,1,i,j,k)+6.0_fp*dy*f(3,4,1,i,j,k) &
                +dx*(2.0_fp*f(4,3,1,i,j,k)+6.0_fp*dy*f(4,4,1,i,j,k))))
           p2= &
                2.0_fp*f(1,3,2,i,j,k)+6.0_fp*dy*f(1,4,2,i,j,k) &
                +dx*(2.0_fp*f(2,3,2,i,j,k)+6.0_fp*dy*f(2,4,2,i,j,k) &
                +dx*(2.0_fp*f(3,3,2,i,j,k)+6.0_fp*dy*f(3,4,2,i,j,k) &
                +dx*(2.0_fp*f(4,3,2,i,j,k)+6.0_fp*dy*f(4,4,2,i,j,k))))
           p3= &
                2.0_fp*f(1,3,3,i,j,k)+6.0_fp*dy*f(1,4,3,i,j,k) &
                +dx*(2.0_fp*f(2,3,3,i,j,k)+6.0_fp*dy*f(2,4,3,i,j,k) &
                +dx*(2.0_fp*f(3,3,3,i,j,k)+6.0_fp*dy*f(3,4,3,i,j,k) &
                +dx*(2.0_fp*f(4,3,3,i,j,k)+6.0_fp*dy*f(4,4,3,i,j,k))))
           p4= &
                2.0_fp*f(1,3,4,i,j,k)+6.0_fp*dy*f(1,4,4,i,j,k) &
                +dx*(2.0_fp*f(2,3,4,i,j,k)+6.0_fp*dy*f(2,4,4,i,j,k) &
                +dx*(2.0_fp*f(3,3,4,i,j,k)+6.0_fp*dy*f(3,4,4,i,j,k) &
                +dx*(2.0_fp*f(4,3,4,i,j,k)+6.0_fp*dy*f(4,4,4,i,j,k))))
           fval(v,iaval)=p1+dz*(p2+dz*(p3+dz*p4))
        end do
     end if
     !
     if((ict(7).gt.0).or.(ict(1).eq.-1)) then
        !  evaluate df2/dz2
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           k=kv(v)
           dx=dxv(v)
           dy=dyv(v)
           dz=dzv(v)
           p3= &
                f(1,1,3,i,j,k)+dy*(f(1,2,3,i,j,k)+dy*(f(1,3,3,i,j,k)+ &
                dy*f(1,4,3,i,j,k))) &
                +dx*(f(2,1,3,i,j,k)+dy*(f(2,2,3,i,j,k)+dy*(f(2,3,3,i,j,k)+ &
                dy*f(2,4,3,i,j,k))) &
                +dx*(f(3,1,3,i,j,k)+dy*(f(3,2,3,i,j,k)+dy*(f(3,3,3,i,j,k)+ &
                dy*f(3,4,3,i,j,k))) &
                +dx*(f(4,1,3,i,j,k)+dy*(f(4,2,3,i,j,k)+dy*(f(4,3,3,i,j,k)+ &
                dy*f(4,4,3,i,j,k))) &
                )))
           p4= &
                f(1,1,4,i,j,k)+dy*(f(1,2,4,i,j,k)+dy*(f(1,3,4,i,j,k)+ &
                dy*f(1,4,4,i,j,k))) &
                +dx*(f(2,1,4,i,j,k)+dy*(f(2,2,4,i,j,k)+dy*(f(2,3,4,i,j,k)+ &
                dy*f(2,4,4,i,j,k))) &
                +dx*(f(3,1,4,i,j,k)+dy*(f(3,2,4,i,j,k)+dy*(f(3,3,4,i,j,k)+ &
                dy*f(3,4,4,i,j,k))) &
                +dx*(f(4,1,4,i,j,k)+dy*(f(4,2,4,i,j,k)+dy*(f(4,3,4,i,j,k)+ &
                dy*f(4,4,4,i,j,k))) &
                )))
           fval(v,iaval)=2.0_fp*p3+6.0_fp*dz*p4
        end do
     end if
     !
     if((ict(8).gt.0).and.(ict(1).ne.-1)) then
        !  evaluate d2f/dxdy
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           k=kv(v)
           dx=dxv(v)
           dy=dyv(v)
           dz=dzv(v)
           p1= &
                f(2,2,1,i,j,k)+ &
                dy*(2.0_fp*f(2,3,1,i,j,k)+dy*3.0_fp*f(2,4,1,i,j,k)) &
                +2._fp*dx*( &
                f(3,2,1,i,j,k)+ &
                dy*(2.0_fp*f(3,3,1,i,j,k)+dy*3.0_fp*f(3,4,1,i,j,k)) &
                +1.5_fp*dx*( &
                f(4,2,1,i,j,k)+ &
                dy*(2.0_fp*f(4,3,1,i,j,k)+dy*3.0_fp*f(4,4,1,i,j,k)) &
                ))
           p2= &
                f(2,2,2,i,j,k)+ &
                dy*(2.0_fp*f(2,3,2,i,j,k)+dy*3.0_fp*f(2,4,2,i,j,k)) &
                +2._fp*dx*( &
                f(3,2,2,i,j,k)+ &
                dy*(2.0_fp*f(3,3,2,i,j,k)+dy*3.0_fp*f(3,4,2,i,j,k)) &
                +1.5_fp*dx*( &
                f(4,2,2,i,j,k)+ &
                dy*(2.0_fp*f(4,3,2,i,j,k)+dy*3.0_fp*f(4,4,2,i,j,k)) &
                ))
           p3= &
                f(2,2,3,i,j,k)+ &
                dy*(2.0_fp*f(2,3,3,i,j,k)+dy*3.0_fp*f(2,4,3,i,j,k)) &
                +2._fp*dx*( &
                f(3,2,3,i,j,k)+ &
                dy*(2.0_fp*f(3,3,3,i,j,k)+dy*3.0_fp*f(3,4,3,i,j,k)) &
                +1.5_fp*dx*( &
                f(4,2,3,i,j,k)+ &
                dy*(2.0_fp*f(4,3,3,i,j,k)+dy*3.0_fp*f(4,4,3,i,j,k)) &
                ))
           p4= &
                f(2,2,4,i,j,k)+ &
                dy*(2.0_fp*f(2,3,4,i,j,k)+dy*3.0_fp*f(2,4,4,i,j,k)) &
                +2._fp*dx*( &
                f(3,2,4,i,j,k)+ &
                dy*(2.0_fp*f(3,3,4,i,j,k)+dy*3.0_fp*f(3,4,4,i,j,k)) &
                +1.5_fp*dx*( &
                f(4,2,4,i,j,k)+ &
                dy*(2.0_fp*f(4,3,4,i,j,k)+dy*3.0_fp*f(4,4,4,i,j,k)) &
                ))
           fval(v,iaval)=p1+dz*(p2+dz*(p3+dz*p4))
        end do
     end if
     !
     if((ict(9).gt.0).and.(ict(1).ne.-1)) then
        !  evaluate d2f/dxdz
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           k=kv(v)
           dx=dxv(v)
           dy=dyv(v)
           dz=dzv(v)
           p2= &
                f(2,1,2,i,j,k)+dy*(f(2,2,2,i,j,k)+dy*(f(2,3,2,i,j,k)+ &
                dy*f(2,4,2,i,j,k))) &
                +2._fp*dx* &
                (f(3,1,2,i,j,k)+dy*(f(3,2,2,i,j,k)+dy*(f(3,3,2,i,j,k)+ &
                dy*f(3,4,2,i,j,k))) &
                +1.5_fp*dx* &
                (f(4,1,2,i,j,k)+dy*(f(4,2,2,i,j,k)+dy*(f(4,3,2,i,j,k)+ &
                dy*f(4,4,2,i,j,k))) &
                ))
           p3= &
                f(2,1,3,i,j,k)+dy*(f(2,2,3,i,j,k)+dy*(f(2,3,3,i,j,k)+ &
                dy*f(2,4,3,i,j,k))) &
                +2._fp*dx* &
                (f(3,1,3,i,j,k)+dy*(f(3,2,3,i,j,k)+dy*(f(3,3,3,i,j,k)+ &
                dy*f(3,4,3,i,j,k))) &
                +1.5_fp*dx* &
                (f(4,1,3,i,j,k)+dy*(f(4,2,3,i,j,k)+dy*(f(4,3,3,i,j,k)+ &
                dy*f(4,4,3,i,j,k))) &
                ))
           p4= &
                f(2,1,4,i,j,k)+dy*(f(2,2,4,i,j,k)+dy*(f(2,3,4,i,j,k)+ &
                dy*f(2,4,4,i,j,k))) &
                +2._fp*dx* &
                (f(3,1,4,i,j,k)+dy*(f(3,2,4,i,j,k)+dy*(f(3,3,4,i,j,k)+ &
                dy*f(3,4,4,i,j,k))) &
                +1.5_fp*dx* &
                (f(4,1,4,i,j,k)+dy*(f(4,2,4,i,j,k)+dy*(f(4,3,4,i,j,k)+ &
                dy*f(4,4,4,i,j,k))) &
                ))
           fval(v,iaval)=p2+dz*(2.0_fp*p3+dz*3.0_fp*p4)
        end do
     end if
     !
     if((ict(10).gt.0).and.(ict(1).ne.-1)) then
        !  evaluate d2f/dydz
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           k=kv(v)
           dx=dxv(v)
           dy=dyv(v)
           dz=dzv(v)
           p2= &
                f(1,2,2,i,j,k)+ &
                dy*(2.0_fp*f(1,3,2,i,j,k)+dy*3.0_fp*f(1,4,2,i,j,k)) &
                +dx*( &
                f(2,2,2,i,j,k)+ &
                dy*(2.0_fp*f(2,3,2,i,j,k)+dy*3.0_fp*f(2,4,2,i,j,k)) &
                +dx*( &
                f(3,2,2,i,j,k)+ &
                dy*(2.0_fp*f(3,3,2,i,j,k)+dy*3.0_fp*f(3,4,2,i,j,k)) &
                +dx*( &
                f(4,2,2,i,j,k)+ &
                dy*(2.0_fp*f(4,3,2,i,j,k)+dy*3.0_fp*f(4,4,2,i,j,k)) &
                )))
           p3= &
                f(1,2,3,i,j,k)+ &
                dy*(2.0_fp*f(1,3,3,i,j,k)+dy*3.0_fp*f(1,4,3,i,j,k)) &
                +dx*( &
                f(2,2,3,i,j,k)+ &
                dy*(2.0_fp*f(2,3,3,i,j,k)+dy*3.0_fp*f(2,4,3,i,j,k)) &
                +dx*( &
                f(3,2,3,i,j,k)+ &
                dy*(2.0_fp*f(3,3,3,i,j,k)+dy*3.0_fp*f(3,4,3,i,j,k)) &
                +dx*( &
                f(4,2,3,i,j,k)+ &
                dy*(2.0_fp*f(4,3,3,i,j,k)+dy*3.0_fp*f(4,4,3,i,j,k)) &
                )))
           p4= &
                f(1,2,4,i,j,k)+ &
                dy*(2.0_fp*f(1,3,4,i,j,k)+dy*3.0_fp*f(1,4,4,i,j,k)) &
                +dx*( &
                f(2,2,4,i,j,k)+ &
                dy*(2.0_fp*f(2,3,4,i,j,k)+dy*3.0_fp*f(2,4,4,i,j,k)) &
                +dx*( &
                f(3,2,4,i,j,k)+ &
                dy*(2.0_fp*f(3,3,4,i,j,k)+dy*3.0_fp*f(3,4,4,i,j,k)) &
                +dx*( &
                f(4,2,4,i,j,k)+ &
                dy*(2.0_fp*f(4,3,4,i,j,k)+dy*3.0_fp*f(4,4,4,i,j,k)) &
                )))
           fval(v,iaval)=p2+dz*(2.0_fp*p3+dz*3.0_fp*p4)
        end do
     end if
     !
     if(ict(1).eq.-1) then
        !  evaluate d4f/dx2dy2
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           k=kv(v)
           dx=dxv(v)
           dy=dyv(v)
           dz=dzv(v)
           p1= &
                4.0_fp*f(3,3,1,i,j,k)+12.0_fp*dy*f(3,4,1,i,j,k) &
                +dx*(12.0_fp*f(4,3,1,i,j,k)+36.0_fp*dy*f(4,4,1,i,j,k))
           p2= &
                4.0_fp*f(3,3,2,i,j,k)+12.0_fp*dy*f(3,4,2,i,j,k) &
                +dx*(12.0_fp*f(4,3,2,i,j,k)+36.0_fp*dy*f(4,4,2,i,j,k))
           p3= &
                4.0_fp*f(3,3,3,i,j,k)+12.0_fp*dy*f(3,4,3,i,j,k) &
                +dx*(12.0_fp*f(4,3,3,i,j,k)+36.0_fp*dy*f(4,4,3,i,j,k))
           p4= &
                4.0_fp*f(3,3,4,i,j,k)+12.0_fp*dy*f(3,4,4,i,j,k) &
                +dx*(12.0_fp*f(4,3,4,i,j,k)+36.0_fp*dy*f(4,4,4,i,j,k))
           fval(v,iaval)=p1+dz*(p2+dz*(p3+dz*p4))
        end do
     end if
     !
     if(ict(1).eq.-1) then
        !  evaluate d4f/dx2dz2
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           k=kv(v)
           dx=dxv(v)
           dy=dyv(v)
           dz=dzv(v)
           p3= 2.0_fp* &
                (f(3,1,3,i,j,k)+dy*(f(3,2,3,i,j,k)+dy*(f(3,3,3,i,j,k)+ &
                dy*f(3,4,3,i,j,k)))) &
                +6.0_fp*dx* &
                (f(4,1,3,i,j,k)+dy*(f(4,2,3,i,j,k)+dy*(f(4,3,3,i,j,k)+ &
                dy*f(4,4,3,i,j,k))))
           p4= 2.0_fp* &
                (f(3,1,4,i,j,k)+dy*(f(3,2,4,i,j,k)+dy*(f(3,3,4,i,j,k)+ &
                dy*f(3,4,4,i,j,k)))) &
                +6.0_fp*dx* &
                (f(4,1,4,i,j,k)+dy*(f(4,2,4,i,j,k)+dy*(f(4,3,4,i,j,k)+ &
                dy*f(4,4,4,i,j,k))))
           fval(v,iaval)=2.0_fp*p3+6.0_fp*dz*p4
        end do
     end if
     !
     if(ict(1).eq.-1) then
        !  evaluate d4f/dy2dz2
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           k=kv(v)
           dx=dxv(v)
           dy=dyv(v)
           dz=dzv(v)
           p3= &
                2.0_fp*f(1,3,3,i,j,k)+6.0_fp*dy*f(1,4,3,i,j,k) &
                +dx*(2.0_fp*f(2,3,3,i,j,k)+6.0_fp*dy*f(2,4,3,i,j,k) &
                +dx*(2.0_fp*f(3,3,3,i,j,k)+6.0_fp*dy*f(3,4,3,i,j,k) &
                +dx*(2.0_fp*f(4,3,3,i,j,k)+6.0_fp*dy*f(4,4,3,i,j,k))))
           p4= &
                2.0_fp*f(1,3,4,i,j,k)+6.0_fp*dy*f(1,4,4,i,j,k) &
                +dx*(2.0_fp*f(2,3,4,i,j,k)+6.0_fp*dy*f(2,4,4,i,j,k) &
                +dx*(2.0_fp*f(3,3,4,i,j,k)+6.0_fp*dy*f(3,4,4,i,j,k) &
                +dx*(2.0_fp*f(4,3,4,i,j,k)+6.0_fp*dy*f(4,4,4,i,j,k))))
           fval(v,iaval)=2.0_fp*p3+6.0_fp*dz*p4
        end do
     end if
     !
     if(ict(1).eq.-1) then
        !  evaluate d6f/dx2dy2dz2
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           k=kv(v)
           dx=dxv(v)
           dy=dyv(v)
           dz=dzv(v)
           p3= &
                4.0_fp*f(3,3,3,i,j,k)+12.0_fp*dy*f(3,4,3,i,j,k) &
                +dx*(12.0_fp*f(4,3,3,i,j,k)+36.0_fp*dy*f(4,4,3,i,j,k))
           p4= &
                4.0_fp*f(3,3,4,i,j,k)+12.0_fp*dy*f(3,4,4,i,j,k) &
                +dx*(12.0_fp*f(4,3,4,i,j,k)+36.0_fp*dy*f(4,4,4,i,j,k))
           fval(v,iaval)=2.0_fp*p3+6.0_fp*dz*p4
        end do
     end if
     !
  end if
  !
  !----------------------------------
  !  3rd derivatives (.le.2 in each coordinate)
  !
  if(ict(1).eq.3) then
     if(ict(2).eq.1) then
        !                               ! fxxy
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           k=kv(v)
           dx=dxv(v)
           dy=dyv(v)
           dz=dzv(v)
           p1= 2.0_fp* &
                (f(3,2,1,i,j,k)+2.0_fp*dy*(f(3,3,1,i,j,k)+ &
                1.5_fp*dy*f(3,4,1,i,j,k))) &
                +6.0_fp*dx* &
                (f(4,2,1,i,j,k)+2.0_fp*dy*(f(4,3,1,i,j,k)+ &
                1.5_fp*dy*f(4,4,1,i,j,k)))
           p2= 2.0_fp* &
                (f(3,2,2,i,j,k)+2.0_fp*dy*(f(3,3,2,i,j,k)+ &
                1.5_fp*dy*f(3,4,2,i,j,k))) &
                +6.0_fp*dx* &
                (f(4,2,2,i,j,k)+2.0_fp*dy*(f(4,3,2,i,j,k)+ &
                1.5_fp*dy*f(4,4,2,i,j,k)))
           p3= 2.0_fp* &
                (f(3,2,3,i,j,k)+2.0_fp*dy*(f(3,3,3,i,j,k)+ &
                1.5_fp*dy*f(3,4,3,i,j,k))) &
                +6.0_fp*dx* &
                (f(4,2,3,i,j,k)+2.0_fp*dy*(f(4,3,3,i,j,k)+ &
                1.5_fp*dy*f(4,4,3,i,j,k)))
           p4= 2.0_fp* &
                (f(3,2,4,i,j,k)+2.0_fp*dy*(f(3,3,4,i,j,k)+ &
                1.5_fp*dy*f(3,4,4,i,j,k))) &
                +6.0_fp*dx* &
                (f(4,2,4,i,j,k)+2.0_fp*dy*(f(4,3,4,i,j,k)+ &
                1.5_fp*dy*f(4,4,4,i,j,k)))
           !
           fval(v,iaval)=p1+dz*(p2+dz*(p3+dz*p4))
        end do
     end if
     !
     if(ict(3).eq.1) then
        !                               ! fxxz
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           k=kv(v)
           dx=dxv(v)
           dy=dyv(v)
           dz=dzv(v)
           !
           p2= 2.0_fp* &
                (f(3,1,2,i,j,k)+dy*(f(3,2,2,i,j,k)+dy*(f(3,3,2,i,j,k)+ &
                dy*f(3,4,2,i,j,k)))) &
                +6.0_fp*dx* &
                (f(4,1,2,i,j,k)+dy*(f(4,2,2,i,j,k)+dy*(f(4,3,2,i,j,k)+ &
                dy*f(4,4,2,i,j,k))))
           p3= 2.0_fp* &
                (f(3,1,3,i,j,k)+dy*(f(3,2,3,i,j,k)+dy*(f(3,3,3,i,j,k)+ &
                dy*f(3,4,3,i,j,k)))) &
                +6.0_fp*dx* &
                (f(4,1,3,i,j,k)+dy*(f(4,2,3,i,j,k)+dy*(f(4,3,3,i,j,k)+ &
                dy*f(4,4,3,i,j,k))))
           p4= 2.0_fp* &
                (f(3,1,4,i,j,k)+dy*(f(3,2,4,i,j,k)+dy*(f(3,3,4,i,j,k)+ &
                dy*f(3,4,4,i,j,k)))) &
                +6.0_fp*dx* &
                (f(4,1,4,i,j,k)+dy*(f(4,2,4,i,j,k)+dy*(f(4,3,4,i,j,k)+ &
                dy*f(4,4,4,i,j,k))))
           fval(v,iaval)=p2+dz*(2.0_fp*p3+dz*3.0_fp*p4)
        end do
     end if
     !
     if(ict(4).eq.1) then
        !                               ! fxyy
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           k=kv(v)
           dx=dxv(v)
           dy=dyv(v)
           dz=dzv(v)
           !
           p1= &
                2.0_fp*f(2,3,1,i,j,k)+6.0_fp*dy*f(2,4,1,i,j,k) &
                +2.0_fp*dx*( &
                2.0_fp*f(3,3,1,i,j,k)+6.0_fp*dy*f(3,4,1,i,j,k) &
                +1.5_fp*dx*( &
                2.0_fp*f(4,3,1,i,j,k)+6.0_fp*dy*f(4,4,1,i,j,k)))
           p2= &
                2.0_fp*f(2,3,2,i,j,k)+6.0_fp*dy*f(2,4,2,i,j,k) &
                +2.0_fp*dx*( &
                2.0_fp*f(3,3,2,i,j,k)+6.0_fp*dy*f(3,4,2,i,j,k) &
                +1.5_fp*dx*( &
                2.0_fp*f(4,3,2,i,j,k)+6.0_fp*dy*f(4,4,2,i,j,k)))
           p3= &
                2.0_fp*f(2,3,3,i,j,k)+6.0_fp*dy*f(2,4,3,i,j,k) &
                +2.0_fp*dx*( &
                2.0_fp*f(3,3,3,i,j,k)+6.0_fp*dy*f(3,4,3,i,j,k) &
                +1.5_fp*dx*( &
                2.0_fp*f(4,3,3,i,j,k)+6.0_fp*dy*f(4,4,3,i,j,k)))
           p4= &
                2.0_fp*f(2,3,4,i,j,k)+6.0_fp*dy*f(2,4,4,i,j,k) &
                +2.0_fp*dx*( &
                2.0_fp*f(3,3,4,i,j,k)+6.0_fp*dy*f(3,4,4,i,j,k) &
                +1.5_fp*dx*( &
                2.0_fp*f(4,3,4,i,j,k)+6.0_fp*dy*f(4,4,4,i,j,k)))
           fval(v,iaval)=p1+dz*(p2+dz*(p3+dz*p4))
        end do
     end if
     !
     if(ict(5).eq.1) then
        !                               ! fxyz
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           k=kv(v)
           dx=dxv(v)
           dy=dyv(v)
           dz=dzv(v)
           !
           p2= &
                f(2,2,2,i,j,k)+ &
                dy*(2.0_fp*f(2,3,2,i,j,k)+dy*3.0_fp*f(2,4,2,i,j,k)) &
                +2._fp*dx*( &
                f(3,2,2,i,j,k)+ &
                dy*(2.0_fp*f(3,3,2,i,j,k)+dy*3.0_fp*f(3,4,2,i,j,k)) &
                +1.5_fp*dx*( &
                f(4,2,2,i,j,k)+ &
                dy*(2.0_fp*f(4,3,2,i,j,k)+dy*3.0_fp*f(4,4,2,i,j,k)) &
                ))
           p3= &
                f(2,2,3,i,j,k)+ &
                dy*(2.0_fp*f(2,3,3,i,j,k)+dy*3.0_fp*f(2,4,3,i,j,k)) &
                +2._fp*dx*( &
                f(3,2,3,i,j,k)+ &
                dy*(2.0_fp*f(3,3,3,i,j,k)+dy*3.0_fp*f(3,4,3,i,j,k)) &
                +1.5_fp*dx*( &
                f(4,2,3,i,j,k)+ &
                dy*(2.0_fp*f(4,3,3,i,j,k)+dy*3.0_fp*f(4,4,3,i,j,k)) &
                ))
           p4= &
                f(2,2,4,i,j,k)+ &
                dy*(2.0_fp*f(2,3,4,i,j,k)+dy*3.0_fp*f(2,4,4,i,j,k)) &
                +2._fp*dx*( &
                f(3,2,4,i,j,k)+ &
                dy*(2.0_fp*f(3,3,4,i,j,k)+dy*3.0_fp*f(3,4,4,i,j,k)) &
                +1.5_fp*dx*( &
                f(4,2,4,i,j,k)+ &
                dy*(2.0_fp*f(4,3,4,i,j,k)+dy*3.0_fp*f(4,4,4,i,j,k)) &
                ))
           fval(v,iaval)=p2+dz*(2.0_fp*p3+dz*3.0_fp*p4)
        end do
     end if
     !
     if(ict(6).eq.1) then
        !                               ! fxzz
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           k=kv(v)
           dx=dxv(v)
           dy=dyv(v)
           dz=dzv(v)
           !
           p3= &
                f(2,1,3,i,j,k)+dy*(f(2,2,3,i,j,k)+dy*(f(2,3,3,i,j,k)+ &
                dy*f(2,4,3,i,j,k))) &
                +2._fp*dx* &
                (f(3,1,3,i,j,k)+dy*(f(3,2,3,i,j,k)+dy*(f(3,3,3,i,j,k)+ &
                dy*f(3,4,3,i,j,k))) &
                +1.5_fp*dx* &
                (f(4,1,3,i,j,k)+dy*(f(4,2,3,i,j,k)+dy*(f(4,3,3,i,j,k)+ &
                dy*f(4,4,3,i,j,k))) &
                ))
           p4= &
                f(2,1,4,i,j,k)+dy*(f(2,2,4,i,j,k)+dy*(f(2,3,4,i,j,k)+ &
                dy*f(2,4,4,i,j,k))) &
                +2._fp*dx* &
                (f(3,1,4,i,j,k)+dy*(f(3,2,4,i,j,k)+dy*(f(3,3,4,i,j,k)+ &
                dy*f(3,4,4,i,j,k))) &
                +1.5_fp*dx* &
                (f(4,1,4,i,j,k)+dy*(f(4,2,4,i,j,k)+dy*(f(4,3,4,i,j,k)+ &
                dy*f(4,4,4,i,j,k))) &
                ))
           fval(v,iaval)=2.0_fp*p3+6.0_fp*dz*p4
        end do
     end if
     !
     if(ict(7).eq.1) then
        !                               ! fyyz
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           k=kv(v)
           dx=dxv(v)
           dy=dyv(v)
           dz=dzv(v)
           !
           p2= &
                2.0_fp*f(1,3,2,i,j,k)+6.0_fp*dy*f(1,4,2,i,j,k) &
                +dx*(2.0_fp*f(2,3,2,i,j,k)+6.0_fp*dy*f(2,4,2,i,j,k) &
                +dx*(2.0_fp*f(3,3,2,i,j,k)+6.0_fp*dy*f(3,4,2,i,j,k) &
                +dx*(2.0_fp*f(4,3,2,i,j,k)+6.0_fp*dy*f(4,4,2,i,j,k))))
           p3= &
                2.0_fp*f(1,3,3,i,j,k)+6.0_fp*dy*f(1,4,3,i,j,k) &
                +dx*(2.0_fp*f(2,3,3,i,j,k)+6.0_fp*dy*f(2,4,3,i,j,k) &
                +dx*(2.0_fp*f(3,3,3,i,j,k)+6.0_fp*dy*f(3,4,3,i,j,k) &
                +dx*(2.0_fp*f(4,3,3,i,j,k)+6.0_fp*dy*f(4,4,3,i,j,k))))
           p4= &
                2.0_fp*f(1,3,4,i,j,k)+6.0_fp*dy*f(1,4,4,i,j,k) &
                +dx*(2.0_fp*f(2,3,4,i,j,k)+6.0_fp*dy*f(2,4,4,i,j,k) &
                +dx*(2.0_fp*f(3,3,4,i,j,k)+6.0_fp*dy*f(3,4,4,i,j,k) &
                +dx*(2.0_fp*f(4,3,4,i,j,k)+6.0_fp*dy*f(4,4,4,i,j,k))))
           fval(v,iaval)=p2+dz*(2.0_fp*p3+dz*3.0_fp*p4)
        end do
     end if
     !
     if(ict(8).eq.1) then
        !                               ! fyzz
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           k=kv(v)
           dx=dxv(v)
           dy=dyv(v)
           dz=dzv(v)
           !
           p3= &
                f(1,2,3,i,j,k)+ &
                dy*(2.0_fp*f(1,3,3,i,j,k)+dy*3.0_fp*f(1,4,3,i,j,k)) &
                +dx*( &
                f(2,2,3,i,j,k)+ &
                dy*(2.0_fp*f(2,3,3,i,j,k)+dy*3.0_fp*f(2,4,3,i,j,k)) &
                +dx*( &
                f(3,2,3,i,j,k)+ &
                dy*(2.0_fp*f(3,3,3,i,j,k)+dy*3.0_fp*f(3,4,3,i,j,k)) &
                +dx*( &
                f(4,2,3,i,j,k)+ &
                dy*(2.0_fp*f(4,3,3,i,j,k)+dy*3.0_fp*f(4,4,3,i,j,k)) &
                )))
           p4= &
                f(1,2,4,i,j,k)+ &
                dy*(2.0_fp*f(1,3,4,i,j,k)+dy*3.0_fp*f(1,4,4,i,j,k)) &
                +dx*( &
                f(2,2,4,i,j,k)+ &
                dy*(2.0_fp*f(2,3,4,i,j,k)+dy*3.0_fp*f(2,4,4,i,j,k)) &
                +dx*( &
                f(3,2,4,i,j,k)+ &
                dy*(2.0_fp*f(3,3,4,i,j,k)+dy*3.0_fp*f(3,4,4,i,j,k)) &
                +dx*( &
                f(4,2,4,i,j,k)+ &
                dy*(2.0_fp*f(4,3,4,i,j,k)+dy*3.0_fp*f(4,4,4,i,j,k)) &
                )))
           fval(v,iaval)=2.0_fp*p3+6.0_fp*dz*p4
        end do
     end if
     !
  end if
  !
  !----------------------------------
  !  3rd derivatives (3 in each coordinate)
  !
  if(ict(1).eq.-3) then
     if(ict(2).eq.1) then
        !                               ! fxxx
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           k=kv(v)
           dx=dxv(v)
           dy=dyv(v)
           dz=dzv(v)
           !
           p1= 6.0_fp* &
                (f(4,1,1,i,j,k)+dy*(f(4,2,1,i,j,k)+dy*(f(4,3,1,i,j,k)+ &
                dy*f(4,4,1,i,j,k))))
           p2= 6.0_fp* &
                (f(4,1,2,i,j,k)+dy*(f(4,2,2,i,j,k)+dy*(f(4,3,2,i,j,k)+ &
                dy*f(4,4,2,i,j,k))))
           p3= 6.0_fp* &
                (f(4,1,3,i,j,k)+dy*(f(4,2,3,i,j,k)+dy*(f(4,3,3,i,j,k)+ &
                dy*f(4,4,3,i,j,k))))
           p4= 6.0_fp* &
                (f(4,1,4,i,j,k)+dy*(f(4,2,4,i,j,k)+dy*(f(4,3,4,i,j,k)+ &
                dy*f(4,4,4,i,j,k))))
           fval(v,iaval)=p1+dz*(p2+dz*(p3+dz*p4))
        end do
     end if
     !
     if(ict(3).eq.1) then
        !                               ! fyyy
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           k=kv(v)
           dx=dxv(v)
           dy=dyv(v)
           dz=dzv(v)
           !
           p1= &
                6.0_fp*(f(1,4,1,i,j,k) &
                +dx*(f(2,4,1,i,j,k) &
                +dx*(f(3,4,1,i,j,k) &
                +dx*f(4,4,1,i,j,k))))
           p2= &
                6.0_fp*(f(1,4,2,i,j,k) &
                +dx*(f(2,4,2,i,j,k) &
                +dx*(f(3,4,2,i,j,k) &
                +dx*f(4,4,2,i,j,k))))
           p3= &
                6.0_fp*(f(1,4,3,i,j,k) &
                +dx*(f(2,4,3,i,j,k) &
                +dx*(f(3,4,3,i,j,k) &
                +dx*f(4,4,3,i,j,k))))
           p4= &
                6.0_fp*(f(1,4,4,i,j,k) &
                +dx*(f(2,4,4,i,j,k) &
                +dx*(f(3,4,4,i,j,k) &
                +dx*f(4,4,4,i,j,k))))
           fval(v,iaval)=p1+dz*(p2+dz*(p3+dz*p4))
        end do
     end if
     !
     if(ict(4).eq.1) then
        !                               ! fzzz
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           k=kv(v)
           dx=dxv(v)
           dy=dyv(v)
           dz=dzv(v)
           !
           p4= &
                f(1,1,4,i,j,k)+dy*(f(1,2,4,i,j,k)+dy*(f(1,3,4,i,j,k)+ &
                dy*f(1,4,4,i,j,k))) &
                +dx*(f(2,1,4,i,j,k)+dy*(f(2,2,4,i,j,k)+dy*(f(2,3,4,i,j,k)+ &
                dy*f(2,4,4,i,j,k))) &
                +dx*(f(3,1,4,i,j,k)+dy*(f(3,2,4,i,j,k)+dy*(f(3,3,4,i,j,k)+ &
                dy*f(3,4,4,i,j,k))) &
                +dx*(f(4,1,4,i,j,k)+dy*(f(4,2,4,i,j,k)+dy*(f(4,3,4,i,j,k)+ &
                dy*f(4,4,4,i,j,k))) &
                )))
           fval(v,iaval)=6.0_fp*p4
        end do
     end if
     !
  end if
  !
  !----------------------------------
  !  4th derivatives (.le.2 in each coordinate)
  !
  if(ict(1).eq.4) then
     if(ict(2).eq.1) then
        !                               ! fxxyy
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           k=kv(v)
           dx=dxv(v)
           dy=dyv(v)
           dz=dzv(v)
           !
           p1= &
                4.0_fp*f(3,3,1,i,j,k)+12.0_fp*dy*f(3,4,1,i,j,k) &
                +12.0_fp*dx*( &
                f(4,3,1,i,j,k)+3.0_fp*dy*f(4,4,1,i,j,k))
           p2= &
                4.0_fp*f(3,3,2,i,j,k)+12.0_fp*dy*f(3,4,2,i,j,k) &
                +12.0_fp*dx*( &
                f(4,3,2,i,j,k)+3.0_fp*dy*f(4,4,2,i,j,k))
           p3= &
                4.0_fp*f(3,3,3,i,j,k)+12.0_fp*dy*f(3,4,3,i,j,k) &
                +12.0_fp*dx*( &
                f(4,3,3,i,j,k)+3.0_fp*dy*f(4,4,3,i,j,k))
           p4= &
                4.0_fp*f(3,3,4,i,j,k)+12.0_fp*dy*f(3,4,4,i,j,k) &
                +12.0_fp*dx*( &
                f(4,3,4,i,j,k)+3.0_fp*dy*f(4,4,4,i,j,k))
           !
           fval(v,iaval)=p1+dz*(p2+dz*(p3+dz*p4))
        end do
     end if
     !
     if(ict(3).eq.1) then
        !                               ! fxxyz
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           k=kv(v)
           dx=dxv(v)
           dy=dyv(v)
           dz=dzv(v)
           !
           p2= 2.0_fp* &
                (f(3,2,2,i,j,k)+2.0_fp*dy*(f(3,3,2,i,j,k)+ &
                1.5_fp*dy*f(3,4,2,i,j,k))) &
                +6.0_fp*dx* &
                (f(4,2,2,i,j,k)+2.0_fp*dy*(f(4,3,2,i,j,k)+ &
                1.5_fp*dy*f(4,4,2,i,j,k)))
           p3= 2.0_fp* &
                (f(3,2,3,i,j,k)+2.0_fp*dy*(f(3,3,3,i,j,k)+ &
                1.5_fp*dy*f(3,4,3,i,j,k))) &
                +6.0_fp*dx* &
                (f(4,2,3,i,j,k)+2.0_fp*dy*(f(4,3,3,i,j,k)+ &
                1.5_fp*dy*f(4,4,3,i,j,k)))
           p4= 2.0_fp* &
                (f(3,2,4,i,j,k)+2.0_fp*dy*(f(3,3,4,i,j,k)+ &
                1.5_fp*dy*f(3,4,4,i,j,k))) &
                +6.0_fp*dx* &
                (f(4,2,4,i,j,k)+2.0_fp*dy*(f(4,3,4,i,j,k)+ &
                1.5_fp*dy*f(4,4,4,i,j,k)))
           !
           fval(v,iaval)=p2+dz*(2.0_fp*p3+dz*3.0_fp*p4)
        end do
     end if
     !
     if(ict(4).eq.1) then
        !                               ! fxxzz
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           k=kv(v)
           dx=dxv(v)
           dy=dyv(v)
           dz=dzv(v)
           !
           p3= 2.0_fp* &
                (f(3,1,3,i,j,k)+dy*(f(3,2,3,i,j,k)+dy*(f(3,3,3,i,j,k)+ &
                dy*f(3,4,3,i,j,k)))) &
                +6.0_fp*dx* &
                (f(4,1,3,i,j,k)+dy*(f(4,2,3,i,j,k)+dy*(f(4,3,3,i,j,k)+ &
                dy*f(4,4,3,i,j,k))))
           p4= 2.0_fp* &
                (f(3,1,4,i,j,k)+dy*(f(3,2,4,i,j,k)+dy*(f(3,3,4,i,j,k)+ &
                dy*f(3,4,4,i,j,k)))) &
                +6.0_fp*dx* &
                (f(4,1,4,i,j,k)+dy*(f(4,2,4,i,j,k)+dy*(f(4,3,4,i,j,k)+ &
                dy*f(4,4,4,i,j,k))))
           !
           fval(v,iaval)=2.0_fp*p3+6.0_fp*dz*p4
        end do
     end if
     !
     if(ict(5).eq.1) then
        !                               ! fxyyz
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           k=kv(v)
           dx=dxv(v)
           dy=dyv(v)
           dz=dzv(v)
           !
           p2= &
                2.0_fp*f(2,3,2,i,j,k)+6.0_fp*dy*f(2,4,2,i,j,k) &
                +2.0_fp*dx*( &
                2.0_fp*f(3,3,2,i,j,k)+6.0_fp*dy*f(3,4,2,i,j,k) &
                +1.5_fp*dx*( &
                2.0_fp*f(4,3,2,i,j,k)+6.0_fp*dy*f(4,4,2,i,j,k)))
           p3= &
                2.0_fp*f(2,3,3,i,j,k)+6.0_fp*dy*f(2,4,3,i,j,k) &
                +2.0_fp*dx*( &
                2.0_fp*f(3,3,3,i,j,k)+6.0_fp*dy*f(3,4,3,i,j,k) &
                +1.5_fp*dx*( &
                2.0_fp*f(4,3,3,i,j,k)+6.0_fp*dy*f(4,4,3,i,j,k)))
           p4= &
                2.0_fp*f(2,3,4,i,j,k)+6.0_fp*dy*f(2,4,4,i,j,k) &
                +2.0_fp*dx*( &
                2.0_fp*f(3,3,4,i,j,k)+6.0_fp*dy*f(3,4,4,i,j,k) &
                +1.5_fp*dx*( &
                2.0_fp*f(4,3,4,i,j,k)+6.0_fp*dy*f(4,4,4,i,j,k)))
           !
           fval(v,iaval)=p2+dz*(2.0_fp*p3+dz*3.0_fp*p4)
        end do
     end if
     !
     if(ict(6).eq.1) then
        !                               ! fxyzz
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           k=kv(v)
           dx=dxv(v)
           dy=dyv(v)
           dz=dzv(v)
           !
           p3= &
                f(2,2,3,i,j,k)+ &
                dy*(2.0_fp*f(2,3,3,i,j,k)+dy*3.0_fp*f(2,4,3,i,j,k)) &
                +2._fp*dx*( &
                f(3,2,3,i,j,k)+ &
                dy*(2.0_fp*f(3,3,3,i,j,k)+dy*3.0_fp*f(3,4,3,i,j,k)) &
                +1.5_fp*dx*( &
                f(4,2,3,i,j,k)+ &
                dy*(2.0_fp*f(4,3,3,i,j,k)+dy*3.0_fp*f(4,4,3,i,j,k)) &
                ))
           p4= &
                f(2,2,4,i,j,k)+ &
                dy*(2.0_fp*f(2,3,4,i,j,k)+dy*3.0_fp*f(2,4,4,i,j,k)) &
                +2._fp*dx*( &
                f(3,2,4,i,j,k)+ &
                dy*(2.0_fp*f(3,3,4,i,j,k)+dy*3.0_fp*f(3,4,4,i,j,k)) &
                +1.5_fp*dx*( &
                f(4,2,4,i,j,k)+ &
                dy*(2.0_fp*f(4,3,4,i,j,k)+dy*3.0_fp*f(4,4,4,i,j,k)) &
                ))
           fval(v,iaval)=2.0_fp*p3+6.0_fp*dz*p4
        end do
     end if
     !
     if(ict(7).eq.1) then
        !                               ! fyyzz
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           k=kv(v)
           dx=dxv(v)
           dy=dyv(v)
           dz=dzv(v)
           !
           p3= &
                2.0_fp*f(1,3,3,i,j,k)+6.0_fp*dy*f(1,4,3,i,j,k) &
                +dx*(2.0_fp*f(2,3,3,i,j,k)+6.0_fp*dy*f(2,4,3,i,j,k) &
                +dx*(2.0_fp*f(3,3,3,i,j,k)+6.0_fp*dy*f(3,4,3,i,j,k) &
                +dx*(2.0_fp*f(4,3,3,i,j,k)+6.0_fp*dy*f(4,4,3,i,j,k))))
           p4= &
                2.0_fp*f(1,3,4,i,j,k)+6.0_fp*dy*f(1,4,4,i,j,k) &
                +dx*(2.0_fp*f(2,3,4,i,j,k)+6.0_fp*dy*f(2,4,4,i,j,k) &
                +dx*(2.0_fp*f(3,3,4,i,j,k)+6.0_fp*dy*f(3,4,4,i,j,k) &
                +dx*(2.0_fp*f(4,3,4,i,j,k)+6.0_fp*dy*f(4,4,4,i,j,k))))
           !
           fval(v,iaval)=2.0_fp*p3+6.0_fp*dz*p4
        end do
     end if
     !
  end if
  !
  !
  !----------------------------------
  !  4th derivatives (3 in a coordinate)
  !
  if(ict(1).eq.-4) then
     if(ict(2).eq.1) then
        !                               ! fxxxy
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           k=kv(v)
           dx=dxv(v)
           dy=dyv(v)
           dz=dzv(v)
           !
           p1= 6.0_fp* &
                (f(4,2,1,i,j,k)+2.0_fp*dy*(f(4,3,1,i,j,k)+ &
                1.5_fp*dy*f(4,4,1,i,j,k)))
           p2= 6.0_fp* &
                (f(4,2,2,i,j,k)+2.0_fp*dy*(f(4,3,2,i,j,k)+ &
                1.5_fp*dy*f(4,4,2,i,j,k)))
           p3= 6.0_fp* &
                (f(4,2,3,i,j,k)+2.0_fp*dy*(f(4,3,3,i,j,k)+ &
                1.5_fp*dy*f(4,4,3,i,j,k)))
           p4= 6.0_fp* &
                (f(4,2,4,i,j,k)+2.0_fp*dy*(f(4,3,4,i,j,k)+ &
                1.5_fp*dy*f(4,4,4,i,j,k)))
           !
           fval(v,iaval)=p1+dz*(p2+dz*(p3+dz*p4))
        end do
     end if
     !
     if(ict(3).eq.1) then
        !                               ! fxxxz
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           k=kv(v)
           dx=dxv(v)
           dy=dyv(v)
           dz=dzv(v)
           !
           p2= 6.0_fp* &
                (f(4,1,2,i,j,k)+dy*(f(4,2,2,i,j,k)+dy*(f(4,3,2,i,j,k)+ &
                dy*f(4,4,2,i,j,k))))
           p3= 6.0_fp* &
                (f(4,1,3,i,j,k)+dy*(f(4,2,3,i,j,k)+dy*(f(4,3,3,i,j,k)+ &
                dy*f(4,4,3,i,j,k))))
           p4= 6.0_fp* &
                (f(4,1,4,i,j,k)+dy*(f(4,2,4,i,j,k)+dy*(f(4,3,4,i,j,k)+ &
                dy*f(4,4,4,i,j,k))))
           !
           fval(v,iaval)=p2+dz*(2.0_fp*p3+dz*3.0_fp*p4)
        end do
     end if
     !
     if(ict(4).eq.1) then
        !                               ! fxyyy
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           k=kv(v)
           dx=dxv(v)
           dy=dyv(v)
           dz=dzv(v)
           !
           p1= &
                6.0_fp*(f(2,4,1,i,j,k) &
                +2.0_fp*dx*(f(3,4,1,i,j,k) &
                +1.5_fp*dx*f(4,4,1,i,j,k)))
           p2= &
                6.0_fp*(f(2,4,2,i,j,k) &
                +2.0_fp*dx*(f(3,4,2,i,j,k) &
                +1.5_fp*dx*f(4,4,2,i,j,k)))
           p3= &
                6.0_fp*(f(2,4,3,i,j,k) &
                +2.0_fp*dx*(f(3,4,3,i,j,k) &
                +1.5_fp*dx*f(4,4,3,i,j,k)))
           p4= &
                6.0_fp*(f(2,4,4,i,j,k) &
                +2.0_fp*dx*(f(3,4,4,i,j,k) &
                +1.5_fp*dx*f(4,4,4,i,j,k)))
           !
           fval(v,iaval)=p1+dz*(p2+dz*(p3+dz*p4))
        end do
     end if
     !
     if(ict(5).eq.1) then
        !                               ! fxzzz
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           k=kv(v)
           dx=dxv(v)
           dy=dyv(v)
           dz=dzv(v)
           !
           p4= &
                f(2,1,4,i,j,k)+dy*(f(2,2,4,i,j,k)+dy*(f(2,3,4,i,j,k)+ &
                dy*f(2,4,4,i,j,k))) &
                +2._fp*dx* &
                (f(3,1,4,i,j,k)+dy*(f(3,2,4,i,j,k)+dy*(f(3,3,4,i,j,k)+ &
                dy*f(3,4,4,i,j,k))) &
                +1.5_fp*dx* &
                (f(4,1,4,i,j,k)+dy*(f(4,2,4,i,j,k)+dy*(f(4,3,4,i,j,k)+ &
                dy*f(4,4,4,i,j,k))) &
                ))
           fval(v,iaval)=6.0_fp*p4
        end do
     end if
     !
     if(ict(6).eq.1) then
        !                               ! fyyyz
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           k=kv(v)
           dx=dxv(v)
           dy=dyv(v)
           dz=dzv(v)
           !
           p2= &
                6.0_fp*(f(1,4,2,i,j,k) &
                +dx*(f(2,4,2,i,j,k) &
                +dx*(f(3,4,2,i,j,k) &
                +dx*f(4,4,2,i,j,k))))
           p3= &
                6.0_fp*(f(1,4,3,i,j,k) &
                +dx*(f(2,4,3,i,j,k) &
                +dx*(f(3,4,3,i,j,k) &
                +dx*f(4,4,3,i,j,k))))
           p4= &
                6.0_fp*(f(1,4,4,i,j,k) &
                +dx*(f(2,4,4,i,j,k) &
                +dx*(f(3,4,4,i,j,k) &
                +dx*f(4,4,4,i,j,k))))
           fval(v,iaval)=p2+dz*(2.0_fp*p3+dz*3.0_fp*p4)
           !
        end do
     end if
     !
     if(ict(7).eq.1) then
        !                               ! fyzzz
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           k=kv(v)
           dx=dxv(v)
           dy=dyv(v)
           dz=dzv(v)
           !
           p4= &
                f(1,2,4,i,j,k)+ &
                dy*(2.0_fp*f(1,3,4,i,j,k)+dy*3.0_fp*f(1,4,4,i,j,k)) &
                +dx*( &
                f(2,2,4,i,j,k)+ &
                dy*(2.0_fp*f(2,3,4,i,j,k)+dy*3.0_fp*f(2,4,4,i,j,k)) &
                +dx*( &
                f(3,2,4,i,j,k)+ &
                dy*(2.0_fp*f(3,3,4,i,j,k)+dy*3.0_fp*f(3,4,4,i,j,k)) &
                +dx*( &
                f(4,2,4,i,j,k)+ &
                dy*(2.0_fp*f(4,3,4,i,j,k)+dy*3.0_fp*f(4,4,4,i,j,k)) &
                )))
           fval(v,iaval)=6.0_fp*p4
        end do
     end if
     !
  end if
  !
  !----------------------------------
  !  5th derivatives (.le.2 in each coordinate)
  !
  if(ict(1).eq.5) then
     if(ict(2).eq.1) then
        !                               ! fxxyyz
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           k=kv(v)
           dx=dxv(v)
           dy=dyv(v)
           dz=dzv(v)
           !
           p2= &
                4.0_fp*f(3,3,2,i,j,k)+12.0_fp*dy*f(3,4,2,i,j,k) &
                +12.0_fp*dx*( &
                f(4,3,2,i,j,k)+3.0_fp*dy*f(4,4,2,i,j,k))
           p3= &
                4.0_fp*f(3,3,3,i,j,k)+12.0_fp*dy*f(3,4,3,i,j,k) &
                +12.0_fp*dx*( &
                f(4,3,3,i,j,k)+3.0_fp*dy*f(4,4,3,i,j,k))
           p4= &
                4.0_fp*f(3,3,4,i,j,k)+12.0_fp*dy*f(3,4,4,i,j,k) &
                +12.0_fp*dx*( &
                f(4,3,4,i,j,k)+3.0_fp*dy*f(4,4,4,i,j,k))
           !
           fval(v,iaval)=p2+dz*(2.0_fp*p3+dz*3.0_fp*p4)
        end do
     end if
     !
     if(ict(3).eq.1) then
        !                               ! fxxyzz
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           k=kv(v)
           dx=dxv(v)
           dy=dyv(v)
           dz=dzv(v)
           !
           p3= 2.0_fp* &
                (f(3,2,3,i,j,k)+2.0_fp*dy*(f(3,3,3,i,j,k)+ &
                1.5_fp*dy*f(3,4,3,i,j,k))) &
                +6.0_fp*dx* &
                (f(4,2,3,i,j,k)+2.0_fp*dy*(f(4,3,3,i,j,k)+ &
                1.5_fp*dy*f(4,4,3,i,j,k)))
           p4= 2.0_fp* &
                (f(3,2,4,i,j,k)+2.0_fp*dy*(f(3,3,4,i,j,k)+ &
                1.5_fp*dy*f(3,4,4,i,j,k))) &
                +6.0_fp*dx* &
                (f(4,2,4,i,j,k)+2.0_fp*dy*(f(4,3,4,i,j,k)+ &
                1.5_fp*dy*f(4,4,4,i,j,k)))
           !
           fval(v,iaval)=2.0_fp*p3+6.0_fp*dz*p4
        end do
     end if
     !
     if(ict(4).eq.1) then
        !                               ! fxyyzz
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           k=kv(v)
           dx=dxv(v)
           dy=dyv(v)
           dz=dzv(v)
           !
           p3= &
                2.0_fp*f(2,3,3,i,j,k)+6.0_fp*dy*f(2,4,3,i,j,k) &
                +2.0_fp*dx*( &
                2.0_fp*f(3,3,3,i,j,k)+6.0_fp*dy*f(3,4,3,i,j,k) &
                +1.5_fp*dx*( &
                2.0_fp*f(4,3,3,i,j,k)+6.0_fp*dy*f(4,4,3,i,j,k)))
           p4= &
                2.0_fp*f(2,3,4,i,j,k)+6.0_fp*dy*f(2,4,4,i,j,k) &
                +2.0_fp*dx*( &
                2.0_fp*f(3,3,4,i,j,k)+6.0_fp*dy*f(3,4,4,i,j,k) &
                +1.5_fp*dx*( &
                2.0_fp*f(4,3,4,i,j,k)+6.0_fp*dy*f(4,4,4,i,j,k)))
           !
           fval(v,iaval)=2.0_fp*p3+6.0_fp*dz*p4
        end do
     end if
     !
  end if
  !
  !----------------------------------
  !  5th derivatives (3 in a coordinate)
  !
  if(ict(1).eq.-5) then
     if(ict(2).eq.1) then
        !                               ! fxxxyy
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           k=kv(v)
           dx=dxv(v)
           dy=dyv(v)
           dz=dzv(v)
           !
           p1=12.0_fp*( &
                f(4,3,1,i,j,k)+3.0_fp*dy*f(4,4,1,i,j,k))
           p2=12.0_fp*( &
                f(4,3,2,i,j,k)+3.0_fp*dy*f(4,4,2,i,j,k))
           p3=12.0_fp*( &
                f(4,3,3,i,j,k)+3.0_fp*dy*f(4,4,3,i,j,k))
           p4=12.0_fp*( &
                f(4,3,4,i,j,k)+3.0_fp*dy*f(4,4,4,i,j,k))
           !
           fval(v,iaval)=p1+dz*(p2+dz*(p3+dz*p4))
        end do
     end if
     !
     if(ict(3).eq.1) then
        !                               ! fxxxyz
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           k=kv(v)
           dx=dxv(v)
           dy=dyv(v)
           dz=dzv(v)
           !
           p2= 6.0_fp* &
                (f(4,2,2,i,j,k)+2.0_fp*dy*(f(4,3,2,i,j,k)+ &
                1.5_fp*dy*f(4,4,2,i,j,k)))
           p3= 6.0_fp* &
                (f(4,2,3,i,j,k)+2.0_fp*dy*(f(4,3,3,i,j,k)+ &
                1.5_fp*dy*f(4,4,3,i,j,k)))
           p4= 6.0_fp* &
                (f(4,2,4,i,j,k)+2.0_fp*dy*(f(4,3,4,i,j,k)+ &
                1.5_fp*dy*f(4,4,4,i,j,k)))
           !
           fval(v,iaval)=p2+dz*(2.0_fp*p3+dz*3.0_fp*p4)
        end do
     end if
     !
     if(ict(4).eq.1) then
        !                               ! fxxxzz
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           k=kv(v)
           dx=dxv(v)
           dy=dyv(v)
           dz=dzv(v)
           !
           p3= 6.0_fp* &
                (f(4,1,3,i,j,k)+dy*(f(4,2,3,i,j,k)+dy*(f(4,3,3,i,j,k)+ &
                dy*f(4,4,3,i,j,k))))
           p4= 6.0_fp* &
                (f(4,1,4,i,j,k)+dy*(f(4,2,4,i,j,k)+dy*(f(4,3,4,i,j,k)+ &
                dy*f(4,4,4,i,j,k))))
           !
           fval(v,iaval)=2.0_fp*p3+6.0_fp*dz*p4
        end do
     end if
     !
     if(ict(5).eq.1) then
        !                               ! fxxyyy
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           k=kv(v)
           dx=dxv(v)
           dy=dyv(v)
           dz=dzv(v)
           !
           p1=12.0_fp*(f(3,4,1,i,j,k)+3.0_fp*dx*f(4,4,1,i,j,k))
           p2=12.0_fp*(f(3,4,2,i,j,k)+3.0_fp*dx*f(4,4,2,i,j,k))
           p3=12.0_fp*(f(3,4,3,i,j,k)+3.0_fp*dx*f(4,4,3,i,j,k))
           p4=12.0_fp*(f(3,4,4,i,j,k)+3.0_fp*dx*f(4,4,4,i,j,k))
           !
           fval(v,iaval)=p1+dz*(p2+dz*(p3+dz*p4))
        end do
     end if
     !
     if(ict(6).eq.1) then
        !                               ! fxxzzz
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           k=kv(v)
           dx=dxv(v)
           dy=dyv(v)
           dz=dzv(v)
           !
           p4= 2.0_fp* &
                (f(3,1,4,i,j,k)+dy*(f(3,2,4,i,j,k)+dy*(f(3,3,4,i,j,k)+ &
                dy*f(3,4,4,i,j,k)))) &
                +6.0_fp*dx* &
                (f(4,1,4,i,j,k)+dy*(f(4,2,4,i,j,k)+dy*(f(4,3,4,i,j,k)+ &
                dy*f(4,4,4,i,j,k))))
           fval(v,iaval)=6.0_fp*p4
        end do
     end if
     !
     if(ict(7).eq.1) then
        !                               ! fxyyyz
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           k=kv(v)
           dx=dxv(v)
           dy=dyv(v)
           dz=dzv(v)
           !
           p2= &
                6.0_fp*(f(2,4,2,i,j,k) &
                +2.0_fp*dx*(f(3,4,2,i,j,k) &
                +1.5_fp*dx*f(4,4,2,i,j,k)))
           p3= &
                6.0_fp*(f(2,4,3,i,j,k) &
                +2.0_fp*dx*(f(3,4,3,i,j,k) &
                +1.5_fp*dx*f(4,4,3,i,j,k)))
           p4= &
                6.0_fp*(f(2,4,4,i,j,k) &
                +2.0_fp*dx*(f(3,4,4,i,j,k) &
                +1.5_fp*dx*f(4,4,4,i,j,k)))
           !
           fval(v,iaval)=p2+dz*(2.0_fp*p3+dz*3.0_fp*p4)
        end do
     end if
     !
     if(ict(8).eq.1) then
        !                               ! fxyzzz
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           k=kv(v)
           dx=dxv(v)
           dy=dyv(v)
           dz=dzv(v)
           !
           p4= &
                f(2,2,4,i,j,k)+ &
                dy*(2.0_fp*f(2,3,4,i,j,k)+dy*3.0_fp*f(2,4,4,i,j,k)) &
                +2._fp*dx*( &
                f(3,2,4,i,j,k)+ &
                dy*(2.0_fp*f(3,3,4,i,j,k)+dy*3.0_fp*f(3,4,4,i,j,k)) &
                +1.5_fp*dx*( &
                f(4,2,4,i,j,k)+ &
                dy*(2.0_fp*f(4,3,4,i,j,k)+dy*3.0_fp*f(4,4,4,i,j,k)) &
                ))
           fval(v,iaval)=6.0_fp*p4
        end do
     end if
     !
     if(ict(9).eq.1) then
        !                               ! fyyyzz
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           k=kv(v)
           dx=dxv(v)
           dy=dyv(v)
           dz=dzv(v)
           !
           p3= &
                6.0_fp*(f(1,4,3,i,j,k) &
                +dx*(f(2,4,3,i,j,k) &
                +dx*(f(3,4,3,i,j,k) &
                +dx*f(4,4,3,i,j,k))))
           p4= &
                6.0_fp*(f(1,4,4,i,j,k) &
                +dx*(f(2,4,4,i,j,k) &
                +dx*(f(3,4,4,i,j,k) &
                +dx*f(4,4,4,i,j,k))))
           fval(v,iaval)=2.0_fp*p3+6.0_fp*dz*p4
        end do
     end if
     !
     if(ict(10).eq.1) then
        !                               ! fyyzzz
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           k=kv(v)
           dx=dxv(v)
           dy=dyv(v)
           dz=dzv(v)
           !
           p4= &
                2.0_fp*f(1,3,4,i,j,k)+6.0_fp*dy*f(1,4,4,i,j,k) &
                +dx*(2.0_fp*f(2,3,4,i,j,k)+6.0_fp*dy*f(2,4,4,i,j,k) &
                +dx*(2.0_fp*f(3,3,4,i,j,k)+6.0_fp*dy*f(3,4,4,i,j,k) &
                +dx*(2.0_fp*f(4,3,4,i,j,k)+6.0_fp*dy*f(4,4,4,i,j,k))))
           fval(v,iaval)=6.0_fp*p4
        end do
     end if
     !
  end if
  !
  !----------------------------------
  !  6th derivatives (2 in each coordinate)
  !
  if(ict(1).eq.6) then
     !                               ! fxxyyzz
     iaval=iaval+1
     do v=1,ivec
        i=iv(v)
        j=jv(v)
        k=kv(v)
        dx=dxv(v)
        dy=dyv(v)
        dz=dzv(v)
        !
        p3= &
             4.0_fp*f(3,3,3,i,j,k)+12.0_fp*dy*f(3,4,3,i,j,k) &
             +12.0_fp*dx*( &
             f(4,3,3,i,j,k)+3.0_fp*dy*f(4,4,3,i,j,k))
        p4= &
             4.0_fp*f(3,3,4,i,j,k)+12.0_fp*dy*f(3,4,4,i,j,k) &
             +12.0_fp*dx*( &
             f(4,3,4,i,j,k)+3.0_fp*dy*f(4,4,4,i,j,k))
        !
        fval(v,iaval)=2.0_fp*p3+6.0_fp*dz*p4
     end do
  end if
  !
  !----------------------------------
  !  6th derivatives (3 in a coordinate)
  !
  if(ict(1).eq.-6) then
     if(ict(2).eq.1) then
        !                               ! fxxxyyy
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           k=kv(v)
           dx=dxv(v)
           dy=dyv(v)
           dz=dzv(v)
           !
           p1=36.0_fp*f(4,4,1,i,j,k)
           p2=36.0_fp*f(4,4,2,i,j,k)
           p3=36.0_fp*f(4,4,3,i,j,k)
           p4=36.0_fp*f(4,4,4,i,j,k)
           !
           fval(v,iaval)=p1+dz*(p2+dz*(p3+dz*p4))
        end do
     end if
     !
     if(ict(3).eq.1) then
        !                               ! fxxxyyz
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           k=kv(v)
           dx=dxv(v)
           dy=dyv(v)
           dz=dzv(v)
           !
           p2=12.0_fp*( &
                f(4,3,2,i,j,k)+3.0_fp*dy*f(4,4,2,i,j,k))
           p3=12.0_fp*( &
                f(4,3,3,i,j,k)+3.0_fp*dy*f(4,4,3,i,j,k))
           p4=12.0_fp*( &
                f(4,3,4,i,j,k)+3.0_fp*dy*f(4,4,4,i,j,k))
           !
           fval(v,iaval)=p2+dz*(2.0_fp*p3+dz*3.0_fp*p4)
        end do
     end if
     !
     if(ict(4).eq.1) then
        !                               ! fxxxyzz
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           k=kv(v)
           dx=dxv(v)
           dy=dyv(v)
           dz=dzv(v)
           !
           p3= 6.0_fp* &
                (f(4,2,3,i,j,k)+2.0_fp*dy*(f(4,3,3,i,j,k)+ &
                1.5_fp*dy*f(4,4,3,i,j,k)))
           p4= 6.0_fp* &
                (f(4,2,4,i,j,k)+2.0_fp*dy*(f(4,3,4,i,j,k)+ &
                1.5_fp*dy*f(4,4,4,i,j,k)))
           !
           fval(v,iaval)=2.0_fp*p3+6.0_fp*dz*p4
        end do
     end if
     !
     if(ict(5).eq.1) then
        !                               ! fxxxzzz
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           k=kv(v)
           dx=dxv(v)
           dy=dyv(v)
           dz=dzv(v)
           !
           p4= 6.0_fp* &
                (f(4,1,4,i,j,k)+dy*(f(4,2,4,i,j,k)+dy*(f(4,3,4,i,j,k)+ &
                dy*f(4,4,4,i,j,k))))
           fval(v,iaval)=6.0_fp*p4
        end do
     end if
     !
     if(ict(6).eq.1) then
        !                               ! fxxyyyz
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           k=kv(v)
           dx=dxv(v)
           dy=dyv(v)
           dz=dzv(v)
           !
           p2=12.0_fp*(f(3,4,2,i,j,k)+3.0_fp*dx*f(4,4,2,i,j,k))
           p3=12.0_fp*(f(3,4,3,i,j,k)+3.0_fp*dx*f(4,4,3,i,j,k))
           p4=12.0_fp*(f(3,4,4,i,j,k)+3.0_fp*dx*f(4,4,4,i,j,k))
           !
           fval(v,iaval)=p2+dz*(2.0_fp*p3+dz*3.0_fp*p4)
        end do
     end if
     !
     if(ict(7).eq.1) then
        !                               ! fxxyzzz
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           k=kv(v)
           dx=dxv(v)
           dy=dyv(v)
           dz=dzv(v)
           !
           p4= 2.0_fp* &
                (f(3,2,4,i,j,k)+2.0_fp*dy*(f(3,3,4,i,j,k)+ &
                1.5_fp*dy*f(3,4,4,i,j,k))) &
                +6.0_fp*dx* &
                (f(4,2,4,i,j,k)+2.0_fp*dy*(f(4,3,4,i,j,k)+ &
                1.5_fp*dy*f(4,4,4,i,j,k)))
           !
           fval(v,iaval)=6.0_fp*p4
        end do
     end if
     !
     if(ict(8).eq.1) then
        !                               ! fxyyyzz
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           k=kv(v)
           dx=dxv(v)
           dy=dyv(v)
           dz=dzv(v)
           !
           p3= &
                6.0_fp*(f(2,4,3,i,j,k) &
                +2.0_fp*dx*(f(3,4,3,i,j,k) &
                +1.5_fp*dx*f(4,4,3,i,j,k)))
           p4= &
                6.0_fp*(f(2,4,4,i,j,k) &
                +2.0_fp*dx*(f(3,4,4,i,j,k) &
                +1.5_fp*dx*f(4,4,4,i,j,k)))
           !
           fval(v,iaval)=2.0_fp*p3+6.0_fp*dz*p4
        end do
     end if
     !
     if(ict(9).eq.1) then
        !                               ! fxyyzzz
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           k=kv(v)
           dx=dxv(v)
           dy=dyv(v)
           dz=dzv(v)
           !
           p4= &
                2.0_fp*f(2,3,4,i,j,k)+6.0_fp*dy*f(2,4,4,i,j,k) &
                +2.0_fp*dx*( &
                2.0_fp*f(3,3,4,i,j,k)+6.0_fp*dy*f(3,4,4,i,j,k) &
                +1.5_fp*dx*( &
                2.0_fp*f(4,3,4,i,j,k)+6.0_fp*dy*f(4,4,4,i,j,k)))
           !
           fval(v,iaval)=6.0_fp*p4
        end do
     end if
     !
     if(ict(10).eq.1) then
        !                               ! fyyyzzz
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           k=kv(v)
           dx=dxv(v)
           dy=dyv(v)
           dz=dzv(v)
           !
           p4= &
                6.0_fp*(f(1,4,4,i,j,k) &
                +dx*(f(2,4,4,i,j,k) &
                +dx*(f(3,4,4,i,j,k) &
                +dx*f(4,4,4,i,j,k))))
           fval(v,iaval)=6.0_fp*p4
        end do
     end if
     !
  end if
  !
  !----------------------------------
  !  7th derivatives
  !
  if(abs(ict(1)).eq.7) then
     if(ict(2).eq.1) then
        !                               ! fxxxyyyz
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           k=kv(v)
           dz=dzv(v)
           !
           p2=36.0_fp*f(4,4,2,i,j,k)
           p3=36.0_fp*f(4,4,3,i,j,k)
           p4=36.0_fp*f(4,4,4,i,j,k)
           fval(v,iaval)=p2+dz*(2.0_fp*p3+dz*3.0_fp*p4)
        end do
     end if
     !
     if(ict(3).eq.1) then
        !                               ! fxxxyyzz
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           k=kv(v)
           dx=dxv(v)
           dy=dyv(v)
           dz=dzv(v)
           !
           p3=12.0_fp*( &
                f(4,3,3,i,j,k)+3.0_fp*dy*f(4,4,3,i,j,k))
           p4=12.0_fp*( &
                f(4,3,4,i,j,k)+3.0_fp*dy*f(4,4,4,i,j,k))
           fval(v,iaval)=2.0_fp*p3+6.0_fp*dz*p4
        end do
     end if
     !
     if(ict(4).eq.1) then
        !                               ! fxxxyzzz
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           k=kv(v)
           dx=dxv(v)
           dy=dyv(v)
           dz=dzv(v)
           !
           p4= 6.0_fp* &
                (f(4,2,4,i,j,k)+2.0_fp*dy*(f(4,3,4,i,j,k)+ &
                1.5_fp*dy*f(4,4,4,i,j,k)))
           !
           fval(v,iaval)=6.0_fp*p4
        end do
     end if
     !
     if(ict(5).eq.1) then
        !                               ! fxxyyyzz
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           k=kv(v)
           dx=dxv(v)
           dy=dyv(v)
           dz=dzv(v)
           !
           p3=12.0_fp*(f(3,4,3,i,j,k)+3.0_fp*dx*f(4,4,3,i,j,k))
           p4=12.0_fp*(f(3,4,4,i,j,k)+3.0_fp*dx*f(4,4,4,i,j,k))
           !
           fval(v,iaval)=2.0_fp*p3+6.0_fp*dz*p4
        end do
     end if
     !
     if(ict(6).eq.1) then
        !                               ! fxxyyzzz
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           k=kv(v)
           dx=dxv(v)
           dy=dyv(v)
           dz=dzv(v)
           !
           p4= &
                4.0_fp*f(3,3,4,i,j,k)+12.0_fp*dy*f(3,4,4,i,j,k) &
                +12.0_fp*dx*( &
                f(4,3,4,i,j,k)+3.0_fp*dy*f(4,4,4,i,j,k))
           !
           fval(v,iaval)=6.0_fp*p4
        end do
     end if
     !
     if(ict(7).eq.1) then
        !                               ! fxyyyzzz
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           k=kv(v)
           dx=dxv(v)
           dy=dyv(v)
           dz=dzv(v)
           !
           p4= &
                6.0_fp*(f(2,4,4,i,j,k) &
                +2.0_fp*dx*(f(3,4,4,i,j,k) &
                +1.5_fp*dx*f(4,4,4,i,j,k)))
           !
           fval(v,iaval)=6.0_fp*p4
        end do
     end if
     !
  end if
  !
  !----------------------------------
  !  8th derivatives
  !
  if(abs(ict(1)).eq.8) then
     if(ict(2).eq.1) then
        !                               ! fxxxyyyzz
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           k=kv(v)
           dz=dzv(v)
           !
           p3=36.0_fp*f(4,4,3,i,j,k)
           p4=36.0_fp*f(4,4,4,i,j,k)
           fval(v,iaval)=2.0_fp*p3+6.0_fp*dz*p4
        end do
     end if
     !
     if(ict(3).eq.1) then
        !                               ! fxxxyyzzz
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           k=kv(v)
           dy=dyv(v)
           !
           p4=12.0_fp*( &
                f(4,3,4,i,j,k)+3.0_fp*dy*f(4,4,4,i,j,k))
           fval(v,iaval)=6.0_fp*p4
        end do
     end if
     !
     if(ict(4).eq.1) then
        !                               ! fxxyyyzzz
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           k=kv(v)
           dx=dxv(v)
           !
           p4=12.0_fp*(f(3,4,4,i,j,k)+3.0_fp*dx*f(4,4,4,i,j,k))
           fval(v,iaval)=6.0_fp*p4
        end do
     end if
     !
  end if
  !
  !----------------------------------
  !  9th derivative
  !
  if(abs(ict(1)).eq.9) then
     !                               ! fxxxyyyzzz
     iaval=iaval+1
     do v=1,ivec
        i=iv(v)
        j=jv(v)
        k=kv(v)
        !
        p4=36.0_fp*f(4,4,4,i,j,k)
        fval(v,iaval)=6.0_fp*p4
     end do
  end if
  !
  return
end subroutine tcspevfn
