subroutine mkintrp3d(x,nx,y,ny,z,nz,jspline, &
     f,icoeff,ixdim,iydim,izdim, &
     ibcxmin,bcxmin,ibcxmax,bcxmax, &
     ibcymin,bcymin,ibcymax,bcymax, &
     ibczmin,bczmin,ibczmax,bczmax, &
     ier)
  use precision_mod, only: fp
  !
  !  setup a tricubic spline, or tricubic Hermite, or hybrid linear/zonal 2d or
  !  1d with 1d or 2d cubic or Hermite spline interpolation
  !
  !
  !  input:
  implicit none
  integer nx                        ! length of x vector
  integer ny                        ! length of y vector
  integer nz                        ! length of z vector
  real(fp) :: x(nx)                        ! x vector, strict ascending
  real(fp) :: y(ny)                        ! y vector, strict ascending
  real(fp) :: z(nz)                        ! z vector, strict ascending
  !
  integer :: jspline(3)             ! interpolation method control
  !        (1) -- 1st dimension: -1: zone, 0: pclin, 1: Hermite, 2: spline
  !        (2) -- 2nd dimension: -1: zone, 0: pclin, 1: Hermite, 2: spline
  !        (3) -- 3rd dimension: -1: zone, 0: pclin, 1: Hermite, 2: spline
  !
  !    Standard interpolation-- all jspline values match
  !      e.g. jspline(1)=jspline(2)=jspline(3)=2 for tricubic spline
  !
  !    Hybrid interpolation-- not all jspline values are the same.
  !
  !    RESTRICTION: if any jspline(...) element has value 1, none can have
  !    value 2.  I.e. Spline and Hermite interpolation cannot currently be mixed.
  !    This restriction exists because of technical issues in the
  !    implementation (it could be removed in principle but the work to do
  !    this has not been scheduled).
  !
  !  coefficient buffer dimensions
  integer :: icoeff                 ! #coefficients per data point
  integer :: ixdim                  ! nx; nx-1 if jspline(1)==-1
  integer :: iydim                  ! ny; ny-1 if jspline(2)==-1
  integer :: izdim                  ! nz; nz-1 if jspline(3)==-1
  !  input/output:
  real(fp) :: f(icoeff,ixdim,iydim,izdim)  ! data and spline coefficients
  !
  !  boundary condition data
  !    Contiguous storage is assumed-- 1st dimension size must match
  !      actual use
  !
  !
  integer ibcxmin,ibcxmax           ! BC type flag @xmin, xmax
  integer ibcymin,ibcymax           ! BC type flag @ymin, ymax
  integer ibczmin,ibczmax           ! BC type flag @zmin, zmax
  !
  real(fp) :: bcxmin(iydim,izdim),bcxmax(iydim,izdim) ! xmin & xmax BC data
  real(fp) :: bcymin(ixdim,izdim),bcymax(ixdim,izdim) ! ymin & ymax BC data
  real(fp) :: bczmin(ixdim,iydim),bczmax(ixdim,iydim) ! zmin & zmax BC data
  !
  !  where BC data is not required, dummy scalars may be passed.
  !  the ibc* flags determine whether BC data isneeded.
  !
  !  BC data not required for zonal or piecewise linear interpolation
  !  for Hermite interpolation ibc* values from set {-1,0,1} are accepted.
  !
  !  BC data:  bcxmin & bcxmax:  BC vs. y,z @xmin,xmax
  !            bcymin & bcymax:  BC vs. x,z @ymin,ymax
  !            bczmin & bczmax:  BC vs. x,y @zmin,zmax
  !
  !   ibcxmin -- indicator for boundary condition at xmin=x(1):
  !    bcxmin(...) -- boundary condition data
  !     =-1 -- use periodic boundary condition
  !     =0 -- use "not a knot"
  !     =1 -- match slope, specified at x(1),y(iy),z(iz) by bcxmin(iy,iz)
  !     =2 -- match 2nd derivative, specified at x(1),y(iy),z(iz)
  !           by bcxmin(iy,iz
  !     =3 -- boundary condition is slope=0 (df/dx=0) at x(1), all y(j)
  !     =4 -- boundary condition is d2f/dx2=0 at x(1), all y(j)
  !     =5 -- df/dx BC from 1st divided difference
  !     =6 -- d2f/dx2 BC from 2nd divided difference (parabolic fit)
  !     =7 -- d3f/dx3 BC from 3rd divided difference (cubic fit)
  !   ***NOTE bcxmin(...) referenced ONLY if ibcxmin=1 or ibcxmin=2
  !
  !   ibcxmax -- indicator for boundary condition at x(nx):
  !    bcxmax(...) -- boundary condition data
  !     (interpretation as with ibcxmin, bcxmin)
  !     NOTE:  if ibcxmin=-1 then the periodic BC applies on both sides
  !            and ibcxmax, bcxmax are ignored.
  !
  !   interpretation of ibcymin,bcymin,ibcymax,bcymax
  !     is same as with ibcxmin,...
  !
  !   interpretation of ibczmin,bczmin,ibczmax,bczmax
  !     is same as with ibcxmin,...
  !
  !   the explicit bdy condition arrays are referenced only if the
  !     corresponding "ibc" flag values are set > 0.
  !
  integer ier                       ! exit code
  !   ier -- completion code, 0 for normal
  !
  !-----------------------
  integer :: kspline
  integer :: ii,jj,imul,imin,imax,ickx,icky,ickz,inum
  integer :: idum1,idum2,idum3,idimcu
  integer :: ipx,ipy,ipz
  real(fp) :: ztol = 0.0001_fp
  logical :: ifound(-1:2)

  real(fp), dimension(:,:), allocatable :: wk2
  real(fp), dimension(:,:,:), allocatable :: wk3
  !-----------------------
  !
  ier=0
  !
  imin=3
  imax=-2
  imul=1
  inum=0
  idimcu=0
  ifound = .FALSE.

  do ii=1,3
    imin=min(imin,jspline(ii))
    imax=max(imax,jspline(ii))
    if(jspline(ii).gt.0) then
      idimcu=idimcu+1
      imul=imul*2
    end if
    if(.not.ifound(jspline(ii))) then
      ifound(jspline(ii)) = .TRUE.
      inum = inum+1
    end if
  end do

  if((imin.lt.-1).or.(imax.gt.2)) then
    ier = 1
    write(6,*) ' ?mkintrp3d: spline type control out of range -1 to 2: ', jspline
  end if
  if(ier.ne.0) return

  if(inum.eq.1) then
    kspline=imin   ! same interp type on all dimensions
  else
    kspline=-99    ! hybrid
    if(ifound(1).and.ifound(2)) then
      ier = 1
      write(6,*) &
           ' ?mkintrp3d: spline/Hermite hybrid not supported (', &
           jspline,')'
    end if
  end if
  if(ier.ne.0) return
  !
  if(imul.ne.icoeff) then
    write(6,*) &
         ' ?coeff dimension inconsistency for spline type codes ', &
         jspline
    write(6,*) ' in mkintrp3d: expected: ',imul,' got: ',icoeff
    ier=1
    return
  end if
  !
  !
  !  check dimensioning consistency
  !
  if(jspline(1).eq.-1) then
    ickx=nx-1
  else
    ickx=nx
  end if
  !
  if(jspline(2).eq.-1) then
    icky=ny-1
  else
    icky=ny
  end if
  !
  if(jspline(3).eq.-1) then
    ickz=nz-1
  else
    ickz=nz
  end if
  !
  if((ickx.ne.ixdim).or.(icky.ne.iydim).or.(ickz.ne.izdim)) then
    write(6,*) &
         ' ?mkintrp2d: dimensioning inconsistent with '// &
         'interpolation controls: ',jspline
    write(6,*) '  expected: ',ickx,icky,ickz, &
         '; got: ',ixdim,iydim,izdim
    ier=1
    return
  end if
  !
  if(jspline(1).le.0) then
    call splinck(x,nx,idum1,ztol,ier)
    if(ier.ne.0) then
      write(6,*) ' ?mkintrp2d: x axis not strict ascending.'
      return
    end if
  end if
  !
  if(jspline(2).le.0) then
    call splinck(y,ny,idum1,ztol,ier)
    if(ier.ne.0) then
      write(6,*) ' ?mkintrp2d: y axis not strict ascending.'
      return
    end if
  end if
  !
  if(jspline(3).le.0) then
    call splinck(z,nz,idum1,ztol,ier)
    if(ier.ne.0) then
      write(6,*) ' ?mkintrp2d: z axis not strict ascending.'
      return
    end if
  end if
  !
  !  if no work to be done: exit now
  if(imul.eq.1) return
  !
  !  check Hermite BCs if necessary
  !
  if(jspline(1).eq.1) then
    if((min(ibcxmin,ibcxmax).lt.-1).or. &
         (max(ibcxmin,ibcxmax).gt.1)) then
      write(6,*) ' ?mkintrp2d: Bdy Cond code out of range for'
      write(6,*) '  Hermite interpolation; (-1:1) allowed, '// &
           'found: ',ibcxmin,ibcxmax
      ier=1
      return
    end if
    ipx=0
    if(ibcxmin.eq.-1) then
      ipx=1
    else if((ibcxmin.eq.1).or.(ibcxmax.eq.1)) then
      ipx=2
    end if
  end if
  !
  if(jspline(2).eq.1) then
    if((min(ibcymin,ibcymax).lt.-1).or. &
         (max(ibcymin,ibcymax).gt.1)) then
      write(6,*) ' ?mkintrp2d: Bdy Cond code out of range for'
      write(6,*) '  Hermite interpolation; (-1:1) allowed, '// &
           'found: ',ibcymin,ibcymax
      ier=1
      return
    end if
    ipy=0
    if(ibcymin.eq.-1) then
      ipy=1
    else if((ibcymin.eq.1).or.(ibcymax.eq.1)) then
      ipy=2
    end if
  end if
  !
  if(jspline(3).eq.1) then
    if((min(ibczmin,ibczmax).lt.-1).or. &
         (max(ibczmin,ibczmax).gt.1)) then
      write(6,*) ' ?mkintrp2d: Bdy Cond code out of range for'
      write(6,*) '  Hermite interpolation; (-1:1) allowed, '// &
           'found: ',ibczmin,ibczmax
      ier=1
      return
    end if
    ipz=0
    if(ibczmin.eq.-1) then
      ipz=1
    else if((ibczmin.eq.1).or.(ibczmax.eq.1)) then
      ipz=2
    end if
  end if
  !
  if(kspline.eq.1) then
    ! tricubic Hermite

    ! put the BCs inside the function data at the right locations...
    call util_bcherm3(f, ixdim, iydim, izdim, &
         ibcxmin,ibcxmax, ibcymin,ibcymax, ibczmin,ibczmax, &
         bcxmin, bcxmax,  bcymin, bcymax,  bczmin, bczmax, &
         x, y, z)

    call akherm3p(x,ixdim, y,iydim, z,izdim, f,ixdim,iydim, &
         idum1,idum2,idum3, ipx,ipy,ipz, ier)

  else if(kspline.eq.2) then
    ! tricubic Spline

    call mktricub(x,nx,y,ny,z,nz, &
         f,ixdim,iydim, &
         ibcxmin,bcxmin,ibcxmax,bcxmax,iydim, &
         ibcymin,bcymin,ibcymax,bcymax,ixdim, &
         ibczmin,bczmin,ibczmax,bczmax,ixdim, &
         idum1,idum2,idum3, ier)

  else
     ! Hybrid
    if(idimcu.eq.1) then
      ! cubic along 1 dimension; other dims are step or pclin

      if(jspline(1).gt.0) then
        ! cubic in x direction
        do jj=1,izdim
          do ii=1,iydim
            if(jspline(1).eq.1) then
              call util_bcherm1(f(1,1,ii,jj), ixdim, &
                   ibcxmin, ibcxmax, &
                   bcxmin(ii,jj), bcxmax(ii,jj), x)
              call akherm1p(x,ixdim,f(1,1,ii,jj),idum1,ipx, &
                   ier)

            else if(jspline(1).eq.2) then
              call mkspline(x,ixdim,f(1,1,ii,jj), &
                   ibcxmin,bcxmin(ii,jj), &
                   ibcxmax,bcxmax(ii,jj), &
                   idum1,ier)
            end if
            if(ier.ne.0) exit
          end do
          if(ier.ne.0) exit
        end do

      else if(jspline(2).gt.0) then
        ! cubic in y direction
        allocate(wk2(2,ny))

        do jj=1,izdim
          do ii=1,ixdim
            wk2(1,1:iydim) = f(1,ii,1:iydim,jj)
            wk2(2,1:iydim) = 0.0_fp
            if(jspline(2).eq.1) then
              call util_bcherm1(wk2, iydim, &
                   ibcymin, ibcymax, &
                   bcymin(ii,jj), bcymax(ii,jj), y)
              call akherm1p(y,iydim,wk2,idum2,ipy,ier)

            else if(jspline(2).eq.2) then
              call mkspline(y,iydim,wk2, &
                   ibcymin,bcymin(ii,jj), &
                   ibcymax,bcymax(ii,jj), &
                   idum2,ier)

            end if
            if(ier.ne.0) exit
            f(1:2,ii,1:iydim,jj) = wk2(1:2,1:iydim)
          end do
          if(ier.ne.0) exit
        end do

        deallocate(wk2)

      else
        ! cubic in z direction
        allocate(wk2(2,nz))

        do jj=1,iydim
          do ii=1,ixdim
            wk2(1,1:izdim) = f(1,ii,jj,1:izdim)
            wk2(2,1:izdim) = 0.0_fp
            if(jspline(3).eq.1) then
              call util_bcherm1(wk2, izdim, &
                   ibczmin, ibczmax, &
                   bczmin(ii,jj), bczmax(ii,jj), z)
              call akherm1p(z,izdim,wk2,idum2,ipz,ier)

            else if(jspline(3).eq.2) then
              call mkspline(z,izdim,wk2, &
                   ibczmin,bczmin(ii,jj), &
                   ibczmax,bczmax(ii,jj), &
                   idum2,ier)

            end if
            if(ier.ne.0) exit
            f(1:2,ii,jj,1:izdim) = wk2(1:2,1:izdim)
          end do
          if(ier.ne.0) exit
        end do

        deallocate(wk2)

      end if
    else
      ! cubic along 2 dimensions
      if(jspline(3).le.0) then
        do ii=1,izdim
          if(jspline(1).eq.1) then
            call util_bcherm2(f(1,1,1,ii),ixdim,iydim, &
                 ibcxmin, ibcxmax, ibcymin, ibcymax, &
                 bcxmin(1,ii), bcxmax(1,ii), &
                 bcymin(1,ii), bcymax(1,ii), &
                 x,y)
            call akherm2p(x,ixdim,y,iydim,f(1,1,1,ii),ixdim, &
                 idum1,idum2,ipx,ipy,ier)
          else
            call mkbicub(x,ixdim,y,iydim,f(1,1,1,ii),ixdim, &
                 ibcxmin, bcxmin(1,ii), ibcxmax, bcxmax(1,ii), &
                 ibcymin, bcymin(1,ii), ibcymax, bcymax(1,ii), &
                 idum1,idum2,ier)
          end if
          if(ier.ne.0) exit
        end do

      else if(jspline(2).le.0) then
        allocate(wk3(4,nx,nz))

        do ii=1,iydim
          wk3 = f(1:4,1:nx,ii,1:nz)
          if(jspline(1).eq.1) then
            call util_bcherm2(wk3,ixdim,izdim, &
                 ibcxmin, ibcxmax, ibczmin, ibczmax, &
                 bcxmin(ii,1:nz), bcxmax(ii,1:nz), &
                 bczmin(1,ii), bczmax(1,ii), x, z)
            call akherm2p(x,ixdim,z,izdim,wk3,ixdim, &
                 idum1,idum2,ipx,ipz,ier)
          else
            call mkbicub(x,ixdim,z,izdim,wk3,ixdim, &
                 ibcxmin,bcxmin(ii,1:nz), &
                 ibcxmax,bcxmax(ii,1:nz), &
                 ibczmin,bczmin(1,ii), ibczmax,bczmax(1,ii), &
                 idum1,idum2,ier)
          end if
          if(ier.ne.0) exit
          f(1:4,1:nx,ii,1:nz) = wk3
        end do

        deallocate(wk3)

      else
        allocate(wk3(4,ny,nz))

        do ii=1,ixdim
          wk3 = f(1:4,ii,1:ny,1:nz)
          if(jspline(2).eq.1) then
            call util_bcherm2(wk3,iydim,izdim, &
                 ibcymin, ibcymax, ibczmin, ibczmax, &
                 bcymin(ii,1:nz), bcymax(ii,1:nz), &
                 bczmin(ii,1:ny), bczmax(ii,1:ny), y, z)
            call akherm2p(y,iydim,z,izdim,wk3,iydim, &
                 idum1,idum2,ipy,ipz,ier)
          else
            call mkbicub(y,iydim,z,izdim,wk3,iydim, &
                 ibcymin,bcymin(ii,1:nz), &
                 ibcymax,bcymax(ii,1:nz), &
                 ibczmin,bczmin(ii,1:ny), &
                 ibczmax,bczmax(ii,1:ny), &
                 idum1,idum2,ier)
          end if
          if(ier.ne.0) exit
          f(1:4,ii,1:ny,1:nz) = wk3
        end do

        deallocate(wk3)

      end if
    end if

  end if

  return
end subroutine mkintrp3d
