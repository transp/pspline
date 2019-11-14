#ifdef _NETCDF

!
!  MOD DMC Apr 2007 -- support Hybrid spline objects.
!   these objects can have "zonal step function" dimensions 1 less than
!   the grid dimension; also, the number of spline coefficients varies
!   depending on how many dimensions are using cubic interpolation.  So,
!   local variables in0, in1, in2, in3 are defined to distinguish these
!   dimensions from the grid dimensions spline_o%n1, etc.

!
! 1-D
!

subroutine EZspline_save1(spline_o, filename, ier, spl_name, fullsave)
  use psp_precision_mod, only: fp
  use EZspline_obj
  use netcdf
  implicit none
  type(EZspline1),  intent(in)  :: spline_o
  character(len=*), intent(in)  :: filename
  integer,          intent(out) :: ier
  !  if SPL_NAME is set, APPEND to file instead of creating new; prepend
  !  name to all NetCDF data items.  This allows one file to contain
  !  multiple items.
  character(len=*), intent(in), optional :: spl_name
  !  if FULLSAVE is set .TRUE., save derived coefficients along with data
  !  this saves recalculating the coefficients when file is read
  logical, intent(in), optional :: fullsave

  integer :: ncid, in0, in1, ix2
  integer :: varid
  integer :: ndims
  integer :: dimid_n1, dimid_n2, dimid_n3
  integer :: dimid_nbcs
  integer :: dimid_x2
  integer :: dimid_in0, dimid_in1
  integer, dimension(3) :: dimids

  character(len=32) :: zpre
  logical :: fullsv, imodify

  ier = 0
  if(spline_o%isReady.ne.1) then
    ier = 91
    return
  end if

  if(present(spl_name)) then
    call ezspline_spl_name_chk(spl_name,ier)
    if(ier.ne.0) return
    zpre=spl_name
  else
    zpre=' '
  end if

  if(present(fullsave)) then
    fullsv=fullsave
  else
    fullsv=.FALSE.
  end if

  if(fullsv) then
    in0 = size(spline_o%fspl,1)
    in1 = size(spline_o%fspl,2)
    ix2 = size(spline_o%x1pkg,2)
  end if

  if(zpre.eq.' ') then
    ier = nf90_create(filename,IOR(NF90_CLOBBER,NF90_64BIT_OFFSET),ncid)
    imodify=.FALSE.
  else
    ier = nf90_open(filename,NF90_WRITE,ncid)
    if(ier.eq.NF90_NOERR) then
      ier = nf90_redef(ncid)  ! start in "define mode".
    else
      ier = nf90_create(filename,IOR(NF90_CLOBBER,NF90_64BIT_OFFSET),ncid)
    endif
    imodify=.TRUE.
  end if
  if(ier.ne.NF90_NOERR .or. ncid.eq.0) then
    ier = 92
    return
  end if

  ier = NF90_NOERR
  if(imodify) then
    ier = nf90_inq_dimid(ncid,trim(zpre)//'n1',dimid_n1)
    ier = nf90_inq_dimid(ncid,trim(zpre)//'nbcs',dimid_nbcs)
    if(fullsv) then
      ier = nf90_inq_dimid(ncid,trim(zpre)//'xnpkg',dimid_x2)
      ier = nf90_inq_dimid(ncid,trim(zpre)//'in0',dimid_in0)
      ier = nf90_inq_dimid(ncid,trim(zpre)//'in1',dimid_in1)
    end if
  end if
  if(ier.ne.NF90_NOERR .or. .not.imodify) then
    ier = nf90_def_dim(ncid,trim(zpre)//'n1',spline_o%n1,dimid_n1)
    ier = nf90_def_dim(ncid,trim(zpre)//'nbcs',2,dimid_nbcs)
    if(fullsv) then
      ier = nf90_def_dim(ncid,trim(zpre)//'xnpkg',ix2,dimid_x2)
      ier = nf90_def_dim(ncid,trim(zpre)//'in0',in0,dimid_in0)
      ier = nf90_def_dim(ncid,trim(zpre)//'in1',in1,dimid_in1)
    end if
  end if
  if(ier.ne.NF90_NOERR) then
    ier = 93
    return
  end if

  ndims = 0
  dimids = (/0, 0, 0/)
  call ezspline_defVar(ncid,trim(zpre)//'klookup1',ndims,dimids,NF90_INT,imodify,varid,ier)
  call ezspline_defVar(ncid,trim(zpre)//'isHermite',ndims,dimids,NF90_INT,imodify,varid,ier)
  call ezspline_defVar(ncid,trim(zpre)//'isLinear',ndims,dimids,NF90_INT,imodify,varid,ier)
  call ezspline_defVar(ncid,trim(zpre)//'isReady',ndims,dimids,NF90_INT,imodify,varid,ier)
  ndims = 1
  dimids = (/dimid_nbcs, 0, 0/)
  call ezspline_defVar(ncid, trim(zpre)//'ibctype1',ndims,dimids,NF90_INT,imodify,varid,ier)
  dimids = (/dimid_n1, 0, 0/)
  call ezspline_defVar(ncid,trim(zpre)//'x1',ndims,dimids,NF90_DOUBLE,imodify,varid,ier)
  ndims = 0
  dimids = (/0, 0, 0/)
  call ezspline_defVar(ncid,trim(zpre)//'bcval1min',ndims,dimids,NF90_DOUBLE,imodify,varid,ier)
  call ezspline_defVar(ncid,trim(zpre)//'bcval1max',ndims,dimids,NF90_DOUBLE,imodify,varid,ier)
  !
  ! save only the function values at the nodes
  !
  if(.not.fullsv) then
    ndims = 1
    dimids = (/dimid_n1, 0, 0/)
    call ezspline_defVar(ncid,trim(zpre)//'f',ndims,dimids,NF90_DOUBLE,imodify,varid,ier)
  else
    ndims = 0
    dimids = (/0, 0, 0/) ! scalar
    call ezspline_defVar(ncid,trim(zpre)//'x1min',ndims,dimids,NF90_DOUBLE,imodify,varid,ier)
    call ezspline_defVar(ncid,trim(zpre)//'x1max',ndims,dimids,NF90_DOUBLE,imodify,varid,ier)
    call ezspline_defVar(ncid,trim(zpre)//'ilin1',ndims,dimids,NF90_INT,imodify,varid,ier)
    ndims = 2
    dimids = (/dimid_n1, dimid_x2, 0/)
    call ezspline_defVar(ncid,trim(zpre)//'x1pkg',ndims,dimids,NF90_DOUBLE,imodify,varid,ier)
    ndims = 2
    dimids = (/dimid_in0, dimid_in1, 0/)
    call ezspline_defVar(ncid,trim(zpre)//'fspl',ndims,dimids,NF90_DOUBLE,imodify,varid,ier)
  end if
  if(ier.eq.NF90_NOERR) ier = nf90_enddef(ncid)
  if(ier.ne.NF90_NOERR) then
    ier = 94
    return
  end if

  ier = nf90_inq_varid(ncid,trim(zpre)//'klookup1',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid,spline_o%klookup1)
  if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'isHermite',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid,spline_o%isHermite)
  if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'isLinear',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid,spline_o%isLinear)
  if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'isReady',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid,spline_o%isReady)
  if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'ibctype1',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid,spline_o%ibctype1)
  if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'x1',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid,spline_o%x1)
  if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'bcval1min',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid,spline_o%bcval1min)
  if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'bcval1max',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid,spline_o%bcval1max)
  if(ier.ne.NF90_NOERR) then
    ier = 95
    return
  end if

  if(spline_o%isReady == 1) then
    if(.not.fullsv) then
      ier = nf90_inq_varid(ncid,trim(zpre)//'f',varid)
      if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid,spline_o%fspl(1,:))
    else
      ier = nf90_inq_varid(ncid,trim(zpre)//'x1min',varid)
      if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid,spline_o%x1min)
      if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'x1max',varid)
      if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid,spline_o%x1max)
      if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'ilin1',varid)
      if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid,spline_o%ilin1)
      if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'x1pkg',varid)
      if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid,spline_o%x1pkg)
      if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'fspl',varid)
      if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid,spline_o%fspl)
    end if
    if(ier.ne.NF90_NOERR) ier = 96
  else
    ier = 96
  end if
  if(ier.ne.0) return

  ier = nf90_close(ncid)
  if(ier.ne.NF90_NOERR) ier = 97

  return
end subroutine EZspline_save1


!! 
!! 2-D
!!

subroutine EZspline_save2(spline_o, filename, ier, spl_name, fullsave)
  use psp_precision_mod, only: fp
  use EZspline_obj
  use netcdf
  implicit none
  type(EZspline2),  intent(in)  :: spline_o
  character(len=*), intent(in)  :: filename
  integer,          intent(out) :: ier
  !  if SPL_NAME is set, APPEND to file instead of creating new; prepend
  !  name to all NetCDF data items.  This allows one file to contain
  !  multiple items.
  character(len=*), intent(in), optional :: spl_name  ! optional spline "name"
  !  if FULLSAVE is set .TRUE., save derived coefficients along with data
  !  this saves recalculating the coefficients when file is read (Mar 2006)
  logical, intent(in), optional :: fullsave

  integer :: ncid, in0, in1, in2, ix2
  integer :: varid
  integer :: ndims
  integer :: dimid_n1, dimid_n2, dimid_n3
  integer :: dimid_nbcs
  integer :: dimid_x2
  integer :: dimid_in0, dimid_in1, dimid_in2
  integer :: dimid_hs
  integer, dimension(3) :: dimids

  character(len=32) :: zpre
  logical :: fullsv, imodify

  ier = 0
  if(spline_o%isReady.ne.1) then
    ier = 91
    return
  end if

  if(present(spl_name)) then
    call ezspline_spl_name_chk(spl_name,ier)
    if(ier.ne.0) return
    zpre=spl_name
  else
    zpre=' '
  end if

  if(present(fullsave)) then
    fullsv=fullsave
  else
    fullsv=.FALSE.
  end if

  in0 = size(spline_o%fspl,1)
  in1 = size(spline_o%fspl,2)
  in2 = size(spline_o%fspl,3)

  if(fullsv) then
    ix2 = size(spline_o%x1pkg,2)
  end if

  if(zpre.eq.' ') then
    ier = nf90_create(filename,IOR(NF90_CLOBBER,NF90_64BIT_OFFSET),ncid)
    imodify=.FALSE.
  else
    ier = nf90_open(filename,NF90_WRITE,ncid)
    if(ier.eq.NF90_NOERR) then
      ier = nf90_redef(ncid)  ! start in "define mode".
    else
      ier = nf90_create(filename,IOR(NF90_CLOBBER,NF90_64BIT_OFFSET),ncid)
    endif
    imodify=.TRUE.
  end if
  if(ier.ne.NF90_NOERR .or. ncid.eq.0) then
    ier = 92
    return
  end if

  ier = NF90_NOERR
  if(imodify) then
    ier = nf90_inq_dimid(ncid,trim(zpre)//'n1',dimid_n1)
    ier = nf90_inq_dimid(ncid,trim(zpre)//'n2',dimid_n2)
    ier = nf90_inq_dimid(ncid,trim(zpre)//'nbcs',dimid_nbcs)
    ier = nf90_inq_dimid(ncid,trim(zpre)//'in0',dimid_in0)
    ier = nf90_inq_dimid(ncid,trim(zpre)//'in1',dimid_in1)
    ier = nf90_inq_dimid(ncid,trim(zpre)//'in2',dimid_in2)
    ier = nf90_inq_dimid(ncid,trim(zpre)//'hs',dimid_hs)
    if(fullsv) then
      ier = nf90_inq_dimid(ncid,trim(zpre)//'xnpkg',dimid_x2)
    end if
  end if
  if(ier.ne.NF90_NOERR .or. .not.imodify) then
    ier = nf90_def_dim(ncid,trim(zpre)//'n1',spline_o%n1,dimid_n1)
    ier = nf90_def_dim(ncid,trim(zpre)//'n2',spline_o%n2,dimid_n2)
    ier = nf90_def_dim(ncid,trim(zpre)//'nbcs',2,dimid_nbcs)
    ier = nf90_def_dim(ncid,trim(zpre)//'in0',in0,dimid_in0)
    ier = nf90_def_dim(ncid,trim(zpre)//'in1',in1,dimid_in1)
    ier = nf90_def_dim(ncid,trim(zpre)//'in2',in2,dimid_in2)
    ier = nf90_def_dim(ncid,trim(zpre)//'hs',2,dimid_hs)
    if(fullsv) then
      ier = nf90_def_dim(ncid,trim(zpre)//'xnpkg',ix2,dimid_x2)
    end if
  end if
  if(ier.ne.NF90_NOERR) then
    ier = 93
    return
  end if

  ndims = 0
  dimids = (/0, 0, 0/)
  call ezspline_defVar(ncid,trim(zpre)//'klookup1',ndims,dimids,NF90_INT,imodify,varid,ier)
  call ezspline_defVar(ncid,trim(zpre)//'klookup2',ndims,dimids,NF90_INT,imodify,varid,ier)
  call ezspline_defVar(ncid,trim(zpre)//'isHermite',ndims,dimids,NF90_INT,imodify,varid,ier)
  call ezspline_defVar(ncid,trim(zpre)//'isLinear',ndims,dimids,NF90_INT,imodify,varid,ier)
  call ezspline_defVar(ncid,trim(zpre)//'isHybrid',ndims,dimids,NF90_INT,imodify,varid,ier)
  call ezspline_defVar(ncid,trim(zpre)//'isReady',ndims,dimids,NF90_INT,imodify,varid,ier)
  ndims = 1
  dimids = (/dimid_hs, 0, 0/)
  call ezspline_defvar(ncid,trim(zpre)//'hspline',ndims,dimids,NF90_INT,imodify,varid,ier)
  dimids = (/dimid_nbcs, 0, 0/)
  call ezspline_defVar(ncid,trim(zpre)//'ibctype1',ndims,dimids,NF90_INT,imodify,varid,ier)
  call ezspline_defVar(ncid,trim(zpre)//'ibctype2',ndims,dimids,NF90_INT,imodify,varid,ier)
  dimids = (/dimid_n1, 0, 0/)
  call ezspline_defVar(ncid,trim(zpre)//'x1',ndims,dimids,NF90_DOUBLE,imodify,varid,ier)
  dimids = (/dimid_n2, 0, 0/)
  call ezspline_defVar(ncid,trim(zpre)//'x2',ndims,dimids,NF90_DOUBLE,imodify,varid,ier)
  dimids = (/dimid_in2, 0, 0/)
  call ezspline_defVar(ncid,trim(zpre)//'bcval1min',ndims,dimids,NF90_DOUBLE,imodify,varid,ier)
  call ezspline_defVar(ncid,trim(zpre)//'bcval1max',ndims,dimids,NF90_DOUBLE,imodify,varid,ier)
  dimids = (/dimid_in1, 0, 0/)
  call ezspline_defVar(ncid,trim(zpre)//'bcval2min',ndims,dimids,NF90_DOUBLE,imodify,varid,ier)
  call ezspline_defVar(ncid,trim(zpre)//'bcval2max',ndims,dimids,NF90_DOUBLE,imodify,varid,ier)
  !
  ! save only the function values at the nodes
  !
  if(.not.fullsv) then
    ndims = 2
    dimids = (/dimid_in1, dimid_in2, 0/)
    call ezspline_defVar(ncid,trim(zpre)//'f',ndims,dimids,NF90_DOUBLE,imodify,varid,ier)
  else
    ndims = 0
    dimids = (/0, 0, 0/) ! scalar
    call ezspline_defVar(ncid,trim(zpre)//'x1min',ndims,dimids,NF90_DOUBLE,imodify,varid,ier)
    call ezspline_defVar(ncid,trim(zpre)//'x1max',ndims,dimids,NF90_DOUBLE,imodify,varid,ier)
    call ezspline_defVar(ncid,trim(zpre)//'ilin1',ndims,dimids,NF90_INT,imodify,varid,ier)
    call ezspline_defVar(ncid,trim(zpre)//'x2min',ndims,dimids,NF90_DOUBLE,imodify,varid,ier)
    call ezspline_defVar(ncid,trim(zpre)//'x2max',ndims,dimids,NF90_DOUBLE,imodify,varid,ier)
    call ezspline_defVar(ncid,trim(zpre)//'ilin2',ndims,dimids,NF90_INT,imodify,varid,ier)
    ndims = 2
    dimids = (/dimid_n1, dimid_x2, 0/)
    call ezspline_defVar(ncid,trim(zpre)//'x1pkg',ndims,dimids,NF90_DOUBLE,imodify,varid,ier)
    dimids = (/dimid_n2, dimid_x2, 0/)
    call ezspline_defVar(ncid,trim(zpre)//'x2pkg',ndims,dimids,NF90_DOUBLE,imodify,varid,ier)
    ndims = 3
    dimids = (/dimid_in0, dimid_in1, dimid_in2/)
    call ezspline_defVar(ncid,trim(zpre)//'fspl',ndims,dimids,NF90_DOUBLE,imodify,varid,ier)
  end if
  if(ier.eq.NF90_NOERR) ier = nf90_enddef(ncid)
  if(ier.ne.NF90_NOERR) then
    ier = 94
    return
  end if

  ier = nf90_inq_varid(ncid,trim(zpre)//'klookup1',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid,spline_o%klookup1)
  if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'klookup2',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid,spline_o%klookup2)
  if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'isHermite',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid,spline_o%isHermite)
  if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'isLinear',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid,spline_o%isLinear)
  if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'isHybrid',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid,spline_o%isHybrid)
  if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'isReady',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid,spline_o%isReady)
  if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'hspline',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid,spline_o%hspline)
  if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'ibctype1',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid,spline_o%ibctype1)
  if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'ibctype2',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid,spline_o%ibctype2)
  if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'x1',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid,spline_o%x1)
  if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'x2',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid,spline_o%x2)
  if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'bcval1min',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid,spline_o%bcval1min)
  if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'bcval1max',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid,spline_o%bcval1max)
  if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'bcval2min',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid,spline_o%bcval2min)
  if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'bcval2max',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid,spline_o%bcval2max)
  if(ier.ne.NF90_NOERR) then
    ier = 95
    return
  end if

  if(spline_o%isReady == 1) then
    if(.not.fullsv) then
      ier = nf90_inq_varid(ncid,trim(zpre)//'f',varid)
      if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid,spline_o%fspl(1,:,:))
    else
      ier = nf90_inq_varid(ncid,trim(zpre)//'x1min',varid)
      if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid,spline_o%x1min)
      if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'x1max',varid)
      if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid,spline_o%x1max)
      if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'ilin1',varid)
      if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid,spline_o%ilin1)
      if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'x2min',varid)
      if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid,spline_o%x2min)
      if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'x2max',varid)
      if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid,spline_o%x2max)
      if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'ilin2',varid)
      if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid,spline_o%ilin2)
      if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'x1pkg',varid)
      if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid,spline_o%x1pkg)
      if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'x2pkg',varid)
      if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid,spline_o%x2pkg)
      if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'fspl',varid)
      if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid,spline_o%fspl)
    end if
    if(ier.ne.NF90_NOERR) ier = 96
  else
    ier = 96
  end if
  if(ier.ne.0) return

  ier = nf90_close(ncid)
  if(ier.ne.0) ier = 97

  return
end subroutine EZspline_save2


!!! 
!!! 3-D
!!!

subroutine EZspline_save3(spline_o, filename, ier, spl_name, fullsave)
  use psp_precision_mod, only: fp
  use EZspline_obj
  use netcdf
  implicit none
  type(EZspline3),  intent(in)  :: spline_o
  character(len=*), intent(in)  :: filename
  integer,          intent(out) :: ier
  !  if SPL_NAME is set, APPEND to file instead of creating new; prepend
  !  name to all NetCDF data items.  This allows one file to contain
  !  multiple items.
  character(len=*), intent(in), optional :: spl_name  ! optional spline "name"
  !  if FULLSAVE is set .TRUE., save derived coefficients along with data
  !  this saves recalculating the coefficients when file is read (Mar 2006)
  logical, intent(in), optional :: fullsave

  integer :: ncid, in0, in1, in2, in3, ix2
  integer :: varid
  integer :: ndims
  integer :: dimid_n1, dimid_n2, dimid_n3
  integer :: dimid_nbcs
  integer :: dimid_x2
  integer :: dimid_in0, dimid_in1, dimid_in2, dimid_in3
  integer :: dimid_hs
  integer, dimension(3) :: dimids

  character(len=32) zpre
  logical :: fullsv, imodify

  ier = 0
  if(spline_o%isReady.ne.1) then
    ier = 91
    return
  end if

  if(present(spl_name)) then
    call ezspline_spl_name_chk(spl_name,ier)
    if(ier.ne.0) return
    zpre=spl_name
  else
    zpre=' '
  end if

  if(present(fullsave)) then
    fullsv=fullsave
  else
    fullsv=.FALSE.
  end if

  in0 = size(spline_o%fspl,1)
  in1 = size(spline_o%fspl,2)
  in2 = size(spline_o%fspl,3)
  in3 = size(spline_o%fspl,4)

  if(fullsv) then
    ix2 = size(spline_o%x1pkg,2)
  end if

  if(zpre.eq.' ') then
    ier = nf90_create(filename,IOR(NF90_CLOBBER,NF90_64BIT_OFFSET),ncid)
    imodify=.FALSE.
  else
    ier = nf90_open(filename,NF90_WRITE,ncid)
    if(ier.eq.NF90_NOERR) then
      ier = nf90_redef(ncid)  ! start in "define mode".
    else
      ier = nf90_create(filename,IOR(NF90_CLOBBER,NF90_64BIT_OFFSET),ncid)
    endif
    imodify=.TRUE.
  end if
  if(ier.ne.NF90_NOERR .or. ncid.eq.0) then
    ier = 92
    return
  end if

  ier = NF90_NOERR
  if(imodify) then
    ier = nf90_inq_dimid(ncid,trim(zpre)//'n1',dimid_n1)
    ier = nf90_inq_dimid(ncid,trim(zpre)//'n2',dimid_n2)
    ier = nf90_inq_dimid(ncid,trim(zpre)//'n3',dimid_n3)
    ier = nf90_inq_dimid(ncid,trim(zpre)//'nbcs',dimid_nbcs)
    ier = nf90_inq_dimid(ncid,trim(zpre)//'in0',dimid_in0)
    ier = nf90_inq_dimid(ncid,trim(zpre)//'in1',dimid_in1)
    ier = nf90_inq_dimid(ncid,trim(zpre)//'in2',dimid_in2)
    ier = nf90_inq_dimid(ncid,trim(zpre)//'in3',dimid_in3)
    ier = nf90_inq_dimid(ncid,trim(zpre)//'hs',dimid_hs)
    if(fullsv) then
      ier = nf90_inq_dimid(ncid,trim(zpre)//'xnpkg',dimid_x2)
    end if
  end if
  if(ier.ne.NF90_NOERR .or. .not.imodify) then
    ier = nf90_def_dim(ncid,trim(zpre)//'n1',spline_o%n1,dimid_n1)
    ier = nf90_def_dim(ncid,trim(zpre)//'n2',spline_o%n2,dimid_n2)
    ier = nf90_def_dim(ncid,trim(zpre)//'n3',spline_o%n3,dimid_n3)
    ier = nf90_def_dim(ncid,trim(zpre)//'nbcs',2,dimid_nbcs)
    ier = nf90_def_dim(ncid,trim(zpre)//'in0',in0*in1,dimid_in0)
    ier = nf90_def_dim(ncid,trim(zpre)//'in1',in1,dimid_in1)
    ier = nf90_def_dim(ncid,trim(zpre)//'in2',in2,dimid_in2)
    ier = nf90_def_dim(ncid,trim(zpre)//'in3',in3,dimid_in3)
    ier = nf90_def_dim(ncid,trim(zpre)//'hs',3,dimid_hs)
    if(fullsv) then
      ier = nf90_def_dim(ncid,trim(zpre)//'xnpkg',ix2,dimid_x2)
    end if
  end if
  if(ier.ne.NF90_NOERR) then
    ier = 93
    return
  end if

  ndims = 0
  dimids = (/0, 0, 0/)
  call ezspline_defVar(ncid,trim(zpre)//'klookup1',ndims,dimids,NF90_INT,imodify,varid,ier)
  call ezspline_defVar(ncid,trim(zpre)//'klookup2',ndims,dimids,NF90_INT,imodify,varid,ier)
  call ezspline_defVar(ncid,trim(zpre)//'klookup3',ndims,dimids,NF90_INT,imodify,varid,ier)
  call ezspline_defVar(ncid,trim(zpre)//'isHermite',ndims,dimids,NF90_INT,imodify,varid,ier)
  call ezspline_defVar(ncid,trim(zpre)//'isLinear',ndims,dimids,NF90_INT,imodify,varid,ier)
  call ezspline_defVar(ncid,trim(zpre)//'isHybrid',ndims,dimids,NF90_INT,imodify,varid,ier)
  call ezspline_defVar(ncid,trim(zpre)//'isReady',ndims,dimids,NF90_INT,imodify,varid,ier)
  ndims = 1
  dimids = (/dimid_hs, 0, 0/)
  call ezspline_defvar(ncid,trim(zpre)//'hspline',ndims,dimids,NF90_INT,imodify,varid,ier)
  dimids = (/dimid_nbcs, 0, 0/)
  call ezspline_defVar(ncid,trim(zpre)//'ibctype1',ndims,dimids,NF90_INT,imodify,varid,ier)
  call ezspline_defVar(ncid,trim(zpre)//'ibctype2',ndims,dimids,NF90_INT,imodify,varid,ier)
  call ezspline_defVar(ncid,trim(zpre)//'ibctype3',ndims,dimids,NF90_INT,imodify,varid,ier)
  dimids = (/dimid_n1, 0, 0/)
  call ezspline_defVar(ncid,trim(zpre)//'x1',ndims,dimids,NF90_DOUBLE,imodify,varid,ier)
  dimids = (/dimid_n2, 0, 0/)
  call ezspline_defVar(ncid,trim(zpre)//'x2',ndims,dimids,NF90_DOUBLE,imodify,varid,ier)
  dimids = (/dimid_n3, 0, 0/)
  call ezspline_defVar(ncid,trim(zpre)//'x3',ndims,dimids,NF90_DOUBLE,imodify,varid,ier)
  ndims = 2
  dimids = (/dimid_in2, dimid_in3, 0/)
  call ezspline_defVar(ncid,trim(zpre)//'bcval1min',ndims,dimids,NF90_DOUBLE,imodify,varid,ier)
  call ezspline_defVar(ncid,trim(zpre)//'bcval1max',ndims,dimids,NF90_DOUBLE,imodify,varid,ier)
  dimids = (/dimid_in1, dimid_in3, 0/)
  call ezspline_defVar(ncid,trim(zpre)//'bcval2min',ndims,dimids,NF90_DOUBLE,imodify,varid,ier)
  call ezspline_defVar(ncid,trim(zpre)//'bcval2max',ndims,dimids,NF90_DOUBLE,imodify,varid,ier)
  dimids = (/dimid_in1, dimid_in2, 0/)
  call ezspline_defVar(ncid,trim(zpre)//'bcval3min',ndims,dimids,NF90_DOUBLE,imodify,varid,ier)
  call ezspline_defVar(ncid,trim(zpre)//'bcval3max',ndims,dimids,NF90_DOUBLE,imodify,varid,ier)
  !
  ! save only the function values at the nodes
  !
  if(.not.fullsv) then
    ndims = 3
    dimids = (/dimid_in1, dimid_in2, dimid_in3/)
    call ezspline_defVar(ncid,trim(zpre)//'f',ndims,dimids,NF90_DOUBLE,imodify,varid,ier)
  else
    ndims = 0
    dimids = (/0, 0, 0/) ! scalar
    call ezspline_defVar(ncid,trim(zpre)//'x1min',ndims,dimids,NF90_DOUBLE,imodify,varid,ier)
    call ezspline_defVar(ncid,trim(zpre)//'x1max',ndims,dimids,NF90_DOUBLE,imodify,varid,ier)
    call ezspline_defVar(ncid,trim(zpre)//'ilin1',ndims,dimids,NF90_INT,imodify,varid,ier)
    call ezspline_defVar(ncid,trim(zpre)//'x2min',ndims,dimids,NF90_DOUBLE,imodify,varid,ier)
    call ezspline_defVar(ncid,trim(zpre)//'x2max',ndims,dimids,NF90_DOUBLE,imodify,varid,ier)
    call ezspline_defVar(ncid,trim(zpre)//'ilin2',ndims,dimids,NF90_INT,imodify,varid,ier)
    call ezspline_defVar(ncid,trim(zpre)//'x3min',ndims,dimids,NF90_DOUBLE,imodify,varid,ier)
    call ezspline_defVar(ncid,trim(zpre)//'x3max',ndims,dimids,NF90_DOUBLE,imodify,varid,ier)
    call ezspline_defVar(ncid,trim(zpre)//'ilin3',ndims,dimids,NF90_INT,imodify,varid,ier)
    ndims = 2
    dimids = (/dimid_n1, dimid_x2, 0/)
    call ezspline_defVar(ncid,trim(zpre)//'x1pkg',ndims,dimids,NF90_DOUBLE,imodify,varid,ier)
    dimids = (/dimid_n2, dimid_x2, 0/)
    call ezspline_defVar(ncid,trim(zpre)//'x2pkg',ndims,dimids,NF90_DOUBLE,imodify,varid,ier)
    dimids = (/dimid_n3, dimid_x2, 0/)
    call ezspline_defVar(ncid,trim(zpre)//'x3pkg',ndims,dimids,NF90_DOUBLE,imodify,varid,ier)
    ndims = 3
    dimids = (/dimid_in0, dimid_in2, dimid_in3/)
    call ezspline_defVar(ncid,trim(zpre)//'fspl',ndims,dimids,NF90_DOUBLE,imodify,varid,ier)
  end if
  if(ier.eq.NF90_NOERR) ier = nf90_enddef(ncid)
  if(ier.ne.NF90_NOERR) then
    ier = 94
    return
  end if

  ier = nf90_inq_varid(ncid,trim(zpre)//'klookup1',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid,spline_o%klookup1)
  if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'klookup2',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid,spline_o%klookup2)
  if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'klookup3',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid,spline_o%klookup3)
  if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'isHermite',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid,spline_o%isHermite)
  if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'isLinear',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid,spline_o%isLinear)
  if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'isHybrid',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid,spline_o%isHybrid)
  if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'isReady',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid,spline_o%isReady)
  if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'hspline',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid,spline_o%hspline)
  if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'ibctype1',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid,spline_o%ibctype1)
  if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'ibctype2',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid,spline_o%ibctype2)
  if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'ibctype3',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid,spline_o%ibctype3)
  if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'x1',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid,spline_o%x1)
  if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'x2',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid,spline_o%x2)
  if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'x3',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid,spline_o%x3)
  if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'bcval1min',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid,spline_o%bcval1min)
  if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'bcval1max',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid,spline_o%bcval1max)
  if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'bcval2min',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid,spline_o%bcval2min)
  if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'bcval2max',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid,spline_o%bcval2max)
  if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'bcval3min',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid,spline_o%bcval3min)
  if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'bcval3max',varid)
  if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid,spline_o%bcval3max)
  if(ier.ne.NF90_NOERR) then
    ier=95
    return
  end if

  if(spline_o%isReady == 1) then
    if(.not.fullsv) then
      ier = nf90_inq_varid(ncid,trim(zpre)//'f',varid)
      if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid,spline_o%fspl(1,:,:,:))
    else
      ier = nf90_inq_varid(ncid,trim(zpre)//'x1min',varid)
      if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid,spline_o%x1min)
      if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'x1max',varid)
      if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid,spline_o%x1max)
      if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'ilin1',varid)
      if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid,spline_o%ilin1)
      if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'x2min',varid)
      if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid,spline_o%x2min)
      if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'x2max',varid)
      if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid,spline_o%x2max)
      if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'ilin2',varid)
      if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid,spline_o%ilin2)
      if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'x3min',varid)
      if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid,spline_o%x3min)
      if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'x3max',varid)
      if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid,spline_o%x3max)
      if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'ilin3',varid)
      if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid,spline_o%ilin3)
      if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'x1pkg',varid)
      if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid,spline_o%x1pkg)
      if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'x2pkg',varid)
      if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid,spline_o%x2pkg)
      if(ier.eq.NF90_NOERR) ier = nf90_inq_varid(ncid,trim(zpre)//'x3pkg',varid)
      if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid,spline_o%x3pkg)
      if(ier.eq.NF90_NOERR) &
           call ezspline_cdfput3(ncid,trim(zpre)//'fspl',spline_o%fspl,in0*in1,in2,in3,ier)
    end if
    if(ier.ne.NF90_NOERR) ier = 96
  else
    ier = 96
  end if
  if(ier.ne.0) return

  ier = nf90_close(ncid)
  if(ier.ne.NF90_NOERR) ier = 97

  return
end subroutine EZspline_save3


subroutine ezspline_spl_name_chk(spl_name,ier)
  character(len=*), intent(in) :: spl_name
  integer, intent(out) :: ier
  integer :: ic,ilen,indx
  character(len=63) :: zlegal = &
       '0123456789_abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'

  ier=0

  ilen = len(trim(spl_name))

  if(ilen.le.0) then
    ier=50
  else if(ilen.gt.20) then
    ier=51
  else
    indx=index(zlegal,spl_name(1:1))
    if(indx.le.10) then
      ier=52
    else
      do ic=2,ilen
        indx=index(zlegal,spl_name(ic:ic))
        if(indx.le.0) ier=52
      end do
    end if
  end if

end subroutine ezspline_spl_name_chk


subroutine ezspline_defVar(ncid, name, ndims, dimids, itype, lmodify, ivarid, ier)
  use netcdf
  implicit none

  ! inquire about variable with logic specific to ezspline_save

  integer, intent(in)               :: ncid    ! NetCDF handle
  character(len=*), intent(in)      :: name    ! item name to define
  integer, intent(in)               :: ndims   ! dimensioning information
  integer, dimension(3), intent(in) :: dimids  ! dimensioning information
  integer, intent(in)               :: itype   ! data type
  logical, intent(in)               :: lmodify ! T if item may already exist
  integer, intent(out)              :: ivarid  ! NetCDF variable id
  integer, intent(inout)            :: ier

  integer               :: itype_loc
  integer               :: ndims_loc
  integer               :: id1,id2,ii
  character(len=128)    :: name_loc
  integer, dimension(3) :: dimids_loc
  logical               :: lmodify_loc

  if(.not.lmodify) lmodify_loc = .false.
  if(lmodify) then
    ier = nf90_inq_varid(ncid,name,ivarid)
    if(ier.ne.NF90_NOERR) then
      lmodify_loc = .false.
    else
      lmodify_loc = .true.
      ! item DOES exist; require dimension and type match...
      ier = nf90_inquire_variable(ncid,ivarid,name_loc,xtype=itype_loc,ndims=ndims_loc)
      if(itype.ne.itype_loc) ier = 51
      if(ndims.ne.ndims_loc) then
        ier = 52
      else
        do ii = 1, ndims
          if(dimids(ii).ne.dimids_loc(ii)) ier = 53
        end do
      end if
    end if
  endif

  if(.not.lmodify_loc) then
    if(ndims.eq.0) then
      ier = nf90_def_var(ncid,name,itype,ivarid)
    else
      ier = nf90_def_var(ncid,name,itype,dimids(1:ndims),ivarid)
    end if
    if(ier.ne.NF90_NOERR) ier=54
  end if

  return
end subroutine ezspline_defVar
#endif
