#ifdef _NETCDF

!  DMC modifications -- spl_name optional argument.
!  a single file can contain more than one spline (or not).  If it does
!  contain more than one spline, then each one must have a distinct name,
!  which must be specified in the call to tell which one to read this time.

!  fullsv -- ezspline_save now has a "full_save" option which means, save
!  the spline including coefficients.  ezspline_load is modified to look
!  for and read fullsv data if it is present; if not, the spline coefficients
!  are recalculated.

!
! 1-D
!

subroutine EZspline_load1(spline_o, filename, ier, spl_name)
  use psp_precision_mod, only: fp
  use EZspline_obj
  use netcdf
  implicit none
  type(EZspline1),            intent(out) :: spline_o
  character(len=*),           intent(in)  :: filename
  integer,                    intent(out) :: ier
  character(len=*), optional, intent(in)  :: spl_name  ! specify name of spline to be read
                                                       ! in case file contains more than one
  integer :: ncid
  integer :: dimid_n1, dimid_n2, dimid_xnpkg
  integer :: n1
  integer :: varid
  integer, dimension(2) :: bcs1
  real(fp), dimension(:), allocatable :: f

  character(len=32) :: zpre
  logical :: fullsv

  if(present(spl_name)) then
    zpre = spl_name
  else
    zpre = ' '
  end if

  ier = nf90_open(filename,NF90_NOWRITE,ncid)
  if(ier.ne.NF90_NOERR .or. ncid.eq.0) then
    ier = 20
    return
  end if

  ! check if n1 is present; check if spline object has correct rank
  !                         n2 should NOT be present.

  ier = nf90_inq_dimid(ncid,trim(zpre)//'n1',dimid_n1)
  if(ier.ne.NF90_NOERR) then
    ier = 21
    return
  end if
  ier = nf90_inq_dimid(ncid,trim(zpre)//'n2',dimid_n2)
  if(ier.eq.NF90_NOERR) then !n2 should NOT be present
    ier = 22
    return
  end if

  ! check if fullsv was in effect -- if so, x1pkg will be found.

  fullsv = .false.
  ier = nf90_inq_dimid(ncid,trim(zpre)//'xnpkg',dimid_xnpkg)
  if(ier.eq.NF90_NOERR) then
    fullsv = .true.
  end if

  if(spline_o%nguard.ne.123456789) call EZspline_preInit(spline_o)

  if(ezspline_allocated(spline_o)) call EZspline_free1(spline_o,ier)

  ier = nf90_inquire_dimension(ncid,dimid_n1,len=n1)
  ier = nf90_inq_varid(ncid,trim(zpre)//'isLinear',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_get_var(ncid,varid,spline_o%isLinear)
  if(ier.ne.NF90_NOERR) spline_o%isLinear = 0

  if(spline_o%isLinear.eq.1) then
    call EZlinear_init1(spline_o,n1,ier)
  else
    bcs1 = (/0, 0/)
    call EZspline_init1(spline_o,n1,bcs1,ier)
  end if
  if(ier.ne.0) return

  ier = nf90_inq_varid(ncid,trim(zpre)//'klookup1',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_get_var(ncid,varid,spline_o%klookup1)
  if(ier.ne.NF90_NOERR) spline_o%klookup1=-3  ! the old default

  if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'isHermite',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_get_var(ncid,varid,spline_o%isHermite)
  if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'isReady',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_get_var(ncid,varid,spline_o%isReady)
  if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'ibctype1',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_get_var(ncid,varid,spline_o%ibctype1)
  if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'x1',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_get_var(ncid,varid,spline_o%x1)
  if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'bcval1min',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_get_var(ncid,varid,spline_o%bcval1min)
  if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'bcval1max',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_get_var(ncid,varid,spline_o%bcval1max)
  if(ier.ne.NF90_NOERR) then
    ier = 23
    return
  end if

  if(.not.fullsv) then
    allocate(f(n1))
    if(spline_o%isReady==1) then
      ier = nf90_inq_varid(ncid,trim(zpre)//'f',varid)
      if(ier.eq.NF90_NOERR) ier = nf90_get_var(ncid,varid,f)
      if(ier.ne.NF90_NOERR) ier = 24
    else
      ier = 24
    end if
  else
    ier = nf90_inq_varid(ncid,trim(zpre)//'x1min',varid)
    if(ier.eq.NF90_NOERR) ier = nf90_get_var(ncid,varid,spline_o%x1min)
    if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'x1max',varid)
    if(ier.eq.NF90_NOERR) ier = nf90_get_var(ncid,varid,spline_o%x1max)
    if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'ilin1',varid)
    if(ier.eq.NF90_NOERR) ier = nf90_get_var(ncid,varid,spline_o%ilin1)
    if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'x1pkg',varid)
    if(ier.eq.NF90_NOERR) ier = nf90_get_var(ncid,varid,spline_o%x1pkg)
    if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'fspl',varid)
    if(ier.eq.NF90_NOERR) ier = nf90_get_var(ncid,varid,spline_o%fspl)
    if(ier.ne.NF90_NOERR) ier = 24
  end if
  if(ier.ne.0) return

  ier = nf90_close(ncid)
  if(ier.ne.NF90_NOERR) ier = 25

  if(.not.fullsv) then
    call EZspline_setup1_x(spline_o, f, ier)
    deallocate(f)
    if(ier.ne.0) ier = 26
  end if

  return
end subroutine EZspline_load1



!
! 2-D
!

subroutine EZspline_load2(spline_o, filename, ier, spl_name)
  use psp_precision_mod, only: fp
  use EZspline_obj
  use netcdf
  implicit none
  type(EZspline2),            intent(out) :: spline_o
  character(len=*),           intent(in)  :: filename
  integer,                    intent(out) :: ier
  character(len=*), optional, intent(in)  :: spl_name  ! specify name of spline to be read
                                                       ! in case file contains more than one

  integer :: ncid
  integer :: dimid_n1, dimid_n2, dimid_n3, dimid_xnpkg
  integer :: dimid_in1, dimid_in2
  integer :: n1, n2
  integer :: in1, in2
  integer :: varid
  integer, dimension(2) :: bcs1, bcs2
  integer, dimension(2) :: hspline
  real(fp), dimension(:,:), allocatable :: f

  character(len=32) zpre
  logical :: fullsv

  if(present(spl_name)) then
    zpre = spl_name
  else
    zpre = ' '
  end if

  ier = nf90_open(filename,NF90_NOWRITE,ncid)
  if(ier.ne.NF90_NOERR .or. ncid.eq.0) then
    ier = 20
    return
  end if

  ! check if n1 is present; check if spline object has correct rank
  !                         n2 should be present but not n3...

  ier = nf90_inq_dimid(ncid,trim(zpre)//'n1',dimid_n1)
  if(ier.ne.NF90_NOERR) then
    ier = 21
    return
  end if
  ier = nf90_inq_dimid(ncid,trim(zpre)//'n2',dimid_n2)
  if(ier.ne.NF90_NOERR) then
    ier = 21
    return
  end if
  ier = nf90_inq_dimid(ncid,trim(zpre)//'n3',dimid_n3)
  if(ier.eq.NF90_NOERR) then !n3 should NOT be present
    ier = 22
    return
  end if

  ! check if fullsv was in effect -- if so, x1pkg will be found.

  fullsv = .false.
  ier = nf90_inq_dimid(ncid,trim(zpre)//'xnpkg',dimid_xnpkg)
  if(ier.eq.NF90_NOERR) then
    fullsv = .true.
  end if

  if(spline_o%nguard /= 123456789) call EZspline_preInit(spline_o)

  if(ezspline_allocated(spline_o)) call EZspline_free2(spline_o,ier)

  ier = nf90_inquire_dimension(ncid,dimid_n1,len=n1)
  ier = nf90_inquire_dimension(ncid,dimid_n2,len=n2)

  ier = nf90_inq_varid(ncid,trim(zpre)//'isLinear',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_get_var(ncid,varid,spline_o%isLinear)
  if(ier.ne.NF90_NOERR) spline_o%isLinear = 0

  ier = nf90_inq_varid(ncid,trim(zpre)//'isHybrid',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_get_var(ncid,varid,spline_o%isHybrid)
  if(ier.ne.NF90_NOERR) spline_o%isHybrid = 0
  spline_o%hspline = 0

  if(spline_o%isLinear.eq.1) then
    call EZlinear_init2(spline_o,n1,n2,ier)
  else if(spline_o%isHybrid.eq.1) then
    ier = nf90_inq_varid(ncid,trim(zpre)//'hspline',varid)
    ier = nf90_get_var(ncid,varid,hspline)
    call EZhybrid_init2_x(spline_o,n1,n2,hspline,ier)
  else
    bcs1 = (/0, 0/); bcs2 = (/0, 0/)
    call EZspline_init2(spline_o,n1,n2,bcs1,bcs2,ier)
  end if
  if(ier.ne.0) return

  ier = nf90_inq_varid(ncid,trim(zpre)//'klookup1',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_get_var(ncid,varid,spline_o%klookup1)
  if(ier.ne.NF90_NOERR) spline_o%klookup1=-3  ! the old default
  ier = nf90_inq_varid(ncid,trim(zpre)//'klookup2',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_get_var(ncid,varid,spline_o%klookup2)
  if(ier.ne.NF90_NOERR) spline_o%klookup2=-3  ! the old default

  if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'isHermite',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_get_var(ncid,varid,spline_o%isHermite)
  if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'isReady',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_get_var(ncid,varid,spline_o%isReady)
  if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'ibctype1',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_get_var(ncid,varid,spline_o%ibctype1)
  if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'ibctype2',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_get_var(ncid,varid,spline_o%ibctype2)
  if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'x1',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_get_var(ncid,varid,spline_o%x1)
  if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'x2',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_get_var(ncid,varid,spline_o%x2)
  if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'bcval1min',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_get_var(ncid,varid,spline_o%bcval1min)
  if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'bcval1max',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_get_var(ncid,varid,spline_o%bcval1max)
  if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'bcval2min',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_get_var(ncid,varid,spline_o%bcval2min)
  if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'bcval2max',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_get_var(ncid,varid,spline_o%bcval2max)
  if(ier.ne.NF90_NOERR) then
    ier = 23
    return
  end if

  if(.not.fullsv) then
    ier = nf90_inq_dimid(ncid,trim(zpre)//'in1',dimid_in1)
    ier = nf90_inquire_dimension(ncid,dimid_in1,len=in1)
    ier = nf90_inq_dimid(ncid,trim(zpre)//'in2',dimid_in2)
    ier = nf90_inquire_dimension(ncid,dimid_in2,len=in2)
    allocate(f(in1,in2))
    if(spline_o%isReady==1) then
      ier = nf90_inq_varid(ncid,trim(zpre)//'f',varid)
      if(ier.eq.NF90_NOERR) ier = nf90_get_var(ncid,varid,f)
      if(ier.ne.NF90_NOERR) ier = 24
    else
      ier = 24
    end if
  else
    ier = nf90_inq_varid(ncid,trim(zpre)//'x1min',varid)
    if(ier.eq.NF90_NOERR) ier = nf90_get_var(ncid,varid,spline_o%x1min)
    if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'x1max',varid)
    if(ier.eq.NF90_NOERR) ier = nf90_get_var(ncid,varid,spline_o%x1max)
    if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'ilin1',varid)
    if(ier.eq.NF90_NOERR) ier = nf90_get_var(ncid,varid,spline_o%ilin1)
    if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'x2min',varid)
    if(ier.eq.NF90_NOERR) ier = nf90_get_var(ncid,varid,spline_o%x2min)
    if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'x2max',varid)
    if(ier.eq.NF90_NOERR) ier = nf90_get_var(ncid,varid,spline_o%x2max)
    if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'ilin2',varid)
    if(ier.eq.NF90_NOERR) ier = nf90_get_var(ncid,varid,spline_o%ilin2)
    if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'x1pkg',varid)
    if(ier.eq.NF90_NOERR) ier = nf90_get_var(ncid,varid,spline_o%x1pkg)
    if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'x2pkg',varid)
    if(ier.eq.NF90_NOERR) ier = nf90_get_var(ncid,varid,spline_o%x2pkg)
    if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'fspl',varid)
    if(ier.eq.NF90_NOERR) ier = nf90_get_var(ncid,varid,spline_o%fspl)
    if(ier.ne.NF90_NOERR) ier = 24
  end if
  if(ier.ne.0) return

  ier = nf90_close(ncid)
  if(ier.ne.NF90_NOERR) ier = 25

  if(.not.fullsv) then
    call EZspline_setup2_x(spline_o, f, ier)
    deallocate(f)
    if(ier.ne.0) ier = 26
  end if

  return
end subroutine EZspline_load2


!
! 3-D
!

subroutine EZspline_load3(spline_o, filename, ier, spl_name)
  use psp_precision_mod, only: fp
  use EZspline_obj
  use netcdf
  implicit none
  type(EZspline3),            intent(out) :: spline_o
  character(len=*),           intent(in)  :: filename
  integer,                    intent(out) :: ier
  character(len=*), optional, intent(in)  :: spl_name  ! specify name of spline to be read
                                                       ! in case file contains more than one

  integer :: ncid
  integer :: dimid_n1, dimid_n2, dimid_n3, dimid_xnpkg
  integer :: dimid_in1, dimid_in2, dimid_in3
  integer :: n1, n2, n3
  integer :: in1, in2, in3
  integer :: varid
  integer, dimension(2) :: bcs1, bcs2, bcs3
  integer, dimension(3) :: hspline
  real(fp), dimension(:,:,:), allocatable :: f

  character(len=32) :: zpre
  logical :: fullsv

  if(present(spl_name)) then
    zpre = spl_name
  else
    zpre = ' '
  end if

  ier = nf90_open(filename,NF90_NOWRITE,ncid)
  if(ier.ne.NF90_NOERR .or. ncid.eq.0) then
    ier = 20
    return
  end if

  ! check if n1 is present; check if spline object has correct rank
  !                         n2 and n3 should also be present.

  ier = nf90_inq_dimid(ncid,trim(zpre)//'n1',dimid_n1)
  if(ier.ne.NF90_NOERR) then
    ier = 21
    return
  end if
  ier = nf90_inq_dimid(ncid,trim(zpre)//'n2',dimid_n2)
  if(ier.ne.NF90_NOERR) then
    ier = 21
    return
  end if
  ier = nf90_inq_dimid(ncid,trim(zpre)//'n3',dimid_n3)
  if(ier.ne.NF90_NOERR) then
    ier = 21
    return
  end if

  ! check if fullsv was in effect -- if so, x1pkg will be found.

  fullsv = .false.
  ier = nf90_inq_dimid(ncid,trim(zpre)//'xnpkg',dimid_xnpkg)
  if(ier.eq.NF90_NOERR) then
    fullsv = .true.
  end if

  if(spline_o%nguard /= 123456789) call EZspline_preInit(spline_o)

  if(ezspline_allocated(spline_o)) call EZspline_free3(spline_o, ier)

  ier = nf90_inquire_dimension(ncid,dimid_n1,len=n1)
  ier = nf90_inquire_dimension(ncid,dimid_n2,len=n2)
  ier = nf90_inquire_dimension(ncid,dimid_n3,len=n3)

  ier = nf90_inq_varid(ncid,trim(zpre)//'isLinear',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_get_var(ncid,varid,spline_o%isLinear)
  if(ier.ne.NF90_NOERR) spline_o%isLinear = 0

  ier = nf90_inq_varid(ncid,trim(zpre)//'isHybrid',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_get_var(ncid,varid,spline_o%isHybrid)
  if(ier.ne.NF90_NOERR) spline_o%isHybrid = 0
  spline_o%hspline = 0

  if(spline_o%isLinear.eq.1) then
    call EZlinear_init3(spline_o,n1,n2,n3,ier)
  else if(spline_o%isHybrid.eq.1) then
    ier = nf90_inq_varid(ncid,trim(zpre)//'hspline',varid)
    ier = nf90_get_var(ncid,varid,hspline)
    call EZhybrid_init3_x(spline_o,n1,n2,n3,hspline,ier)
  else
    bcs1 = (/0, 0/); bcs2 = (/0, 0/); bcs3 = (/0, 0/)
    call EZspline_init3(spline_o,n1,n2,n3,bcs1,bcs2,bcs3,ier)
  end if
  if(ier.ne.0) return

  ier = nf90_inq_varid(ncid,trim(zpre)//'klookup1',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_get_var(ncid,varid,spline_o%klookup1)
  if(ier.ne.NF90_NOERR) spline_o%klookup1=-3  ! the old default
  ier = nf90_inq_varid(ncid,trim(zpre)//'klookup2',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_get_var(ncid,varid,spline_o%klookup2)
  if(ier.ne.NF90_NOERR) spline_o%klookup2=-3  ! the old default
  ier = nf90_inq_varid(ncid,trim(zpre)//'klookup3',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_get_var(ncid,varid,spline_o%klookup3)
  if(ier.ne.NF90_NOERR) spline_o%klookup3=-3  ! the old default

  if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'isHermite',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_get_var(ncid,varid,spline_o%isHermite)
  if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'isReady',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_get_var(ncid,varid,spline_o%isReady)
  if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'ibctype1',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_get_var(ncid,varid,spline_o%ibctype1)
  if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'ibctype2',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_get_var(ncid,varid,spline_o%ibctype2)
  if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'ibctype3',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_get_var(ncid,varid,spline_o%ibctype3)
  if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'x1',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_get_var(ncid,varid,spline_o%x1)
  if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'x2',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_get_var(ncid,varid,spline_o%x2)
  if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'x3',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_get_var(ncid,varid,spline_o%x3)
  if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'bcval1min',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_get_var(ncid,varid,spline_o%bcval1min)
  if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'bcval1max',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_get_var(ncid,varid,spline_o%bcval1max)
  if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'bcval2min',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_get_var(ncid,varid,spline_o%bcval2min)
  if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'bcval2max',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_get_var(ncid,varid,spline_o%bcval2max)
  if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'bcval3min',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_get_var(ncid,varid,spline_o%bcval3min)
  if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'bcval3max',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_get_var(ncid,varid,spline_o%bcval3max)
  if(ier.ne.NF90_NOERR) then
    ier = 23
    return
  end if

  if(.not.fullsv) then
    ier = nf90_inq_dimid(ncid,trim(zpre)//'in1',dimid_in1)
    ier = nf90_inquire_dimension(ncid,dimid_in1,len=in1)
    ier = nf90_inq_dimid(ncid,trim(zpre)//'in2',dimid_in2)
    ier = nf90_inquire_dimension(ncid,dimid_in2,len=in2)
    ier = nf90_inq_dimid(ncid,trim(zpre)//'in2',dimid_in3)
    ier = nf90_inquire_dimension(ncid,dimid_in2,len=in3)
    allocate(f(in1,in2,in3))
    if(spline_o%isReady==1) then
      ier = nf90_inq_varid(ncid,trim(zpre)//'f',varid)
      if(ier.eq.NF90_NOERR) ier = nf90_get_var(ncid,varid,f)
      if(ier.ne.NF90_NOERR) ier = 24
    else
      ier = 24
    end if
  else
    ier = nf90_inq_varid(ncid,trim(zpre)//'x1min',varid)
    if(ier.eq.NF90_NOERR) ier = nf90_get_var(ncid,varid,spline_o%x1min)
    if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'x1max',varid)
    if(ier.eq.NF90_NOERR) ier = nf90_get_var(ncid,varid,spline_o%x1max)
    if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'ilin1',varid)
    if(ier.eq.NF90_NOERR) ier = nf90_get_var(ncid,varid,spline_o%ilin1)
    if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'x2min',varid)
    if(ier.eq.NF90_NOERR) ier = nf90_get_var(ncid,varid,spline_o%x2min)
    if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'x2max',varid)
    if(ier.eq.NF90_NOERR) ier = nf90_get_var(ncid,varid,spline_o%x2max)
    if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'ilin2',varid)
    if(ier.eq.NF90_NOERR) ier = nf90_get_var(ncid,varid,spline_o%ilin2)
    if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'x3min',varid)
    if(ier.eq.NF90_NOERR) ier = nf90_get_var(ncid,varid,spline_o%x3min)
    if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'x3max',varid)
    if(ier.eq.NF90_NOERR) ier = nf90_get_var(ncid,varid,spline_o%x3max)
    if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'ilin3',varid)
    if(ier.eq.NF90_NOERR) ier = nf90_get_var(ncid,varid,spline_o%ilin3)
    if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'x1pkg',varid)
    if(ier.eq.NF90_NOERR) ier = nf90_get_var(ncid,varid,spline_o%x1pkg)
    if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'x2pkg',varid)
    if(ier.eq.NF90_NOERR) ier = nf90_get_var(ncid,varid,spline_o%x2pkg)
    if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'x3pkg',varid)
    if(ier.eq.NF90_NOERR) ier = nf90_get_var(ncid,varid,spline_o%x3pkg)
    if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'fspl',varid)
    if(ier.eq.NF90_NOERR) ier = nf90_get_var(ncid,varid,spline_o%fspl)
    if(ier.ne.NF90_NOERR) ier = 24
  end if
  if(ier.ne.0) return

  ier = nf90_close(ncid)
  if(ier.ne.NF90_NOERR) ier = 25

  if(.not.fullsv) then
    call EZspline_setup3_x(spline_o, f, ier)
    deallocate(f)
    if(ier.ne.0) ier = 26
  end if

  return
end subroutine EZspline_load3

#endif
