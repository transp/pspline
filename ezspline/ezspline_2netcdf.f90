#ifdef _NETCDF

subroutine EZspline_2NetCDF_array3(n1, n2, n3, x1, x2, x3, f, filename, ier)
  use psp_precision_mod, only: fp
  use EZspline_obj   
  use netcdf
  implicit none
  integer,                    intent(in)  :: n1, n2, n3
  real(fp), dimension(:),     intent(in)  :: x1, x2, x3
  real(fp), dimension(:,:,:), intent(in)  :: f
  character(len=*),           intent(in)  :: filename
  integer,                    intent(out) :: ier

  integer, dimension(3) :: dimid_n
  integer               :: varid_x1, varid_x2, varid_x3
  integer               :: varid_f
  integer               :: ncid

  ier = nf90_create(filename,IOR(NF90_CLOBBER,NF90_64BIT_OFFSET),ncid)
  if(ier.ne.NF90_NOERR .or. ncid.eq.0) then
    ier = 33
    return
  end if

  ier = nf90_def_dim(ncid,'n1',n1,dimid_n(1))
  if(ier.eq.NF90_NOERR) ier = nf90_def_dim(ncid,'n2',n2,dimid_n(2))
  if(ier.eq.NF90_NOERR) ier = nf90_def_dim(ncid,'n3',n3,dimid_n(3))
  if(ier.eq.NF90_NOERR) ier = nf90_def_var(ncid,'x1',NF90_DOUBLE,dimid_n(1),varid_x1)
  if(ier.eq.NF90_NOERR) ier = nf90_def_var(ncid,'x2',NF90_DOUBLE,dimid_n(2),varid_x2)
  if(ier.eq.NF90_NOERR) ier = nf90_def_var(ncid,'x3',NF90_DOUBLE,dimid_n(3),varid_x3)
  if(ier.eq.NF90_NOERR) ier = nf90_def_var(ncid,'f',NF90_DOUBLE,dimid_n,varid_f)
  if(ier.eq.NF90_NOERR) ier = nf90_enddef(ncid)
  if(ier.ne.NF90_NOERR) then
    ier = 34
    return
  end if

  ier = nf90_put_var(ncid,varid_x1,x1)
  if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid_x2,x2)
  if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid_x3,x3)
  if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid_f,f)
  if(ier.ne.NF90_NOERR) then
    ier = 35
    return
  end if

  ier = nf90_close(ncid)
  if(ier.ne.NF90_NOERR) ier = 36

  return
end subroutine EZspline_2NetCDF_array3


subroutine EZspline_2NetCDF_array2(n1, n2, x1, x2, f, filename, ier)
  use psp_precision_mod, only: fp
  use EZspline_obj   
  use netcdf
  implicit none
  integer,                  intent(in)  :: n1, n2
  real(fp), dimension(:),   intent(in)  :: x1, x2
  real(fp), dimension(:,:), intent(in)  :: f
  character(len=*),         intent(in)  :: filename
  integer,                  intent(out) :: ier

  integer, dimension(2) :: dimid_n
  integer               :: varid_x1, varid_x2
  integer               :: varid_f
  integer               :: ncid

  ier = nf90_create(filename,IOR(NF90_CLOBBER,NF90_64BIT_OFFSET),ncid)
  if(ier.ne.NF90_NOERR .or. ncid.eq.0) then
    ier = 33
    return
  end if

  ier = nf90_def_dim(ncid,'n1',n1,dimid_n(1))
  if(ier.eq.NF90_NOERR) ier = nf90_def_dim(ncid,'n2',n2,dimid_n(2))
  if(ier.eq.NF90_NOERR) ier = nf90_def_var(ncid,'x1',NF90_DOUBLE,dimid_n(1),varid_x1)
  if(ier.eq.NF90_NOERR) ier = nf90_def_var(ncid,'x2',NF90_DOUBLE,dimid_n(2),varid_x2)
  if(ier.eq.NF90_NOERR) ier = nf90_def_var(ncid,'f',NF90_DOUBLE,dimid_n,varid_f)
  if(ier.eq.NF90_NOERR) ier = nf90_enddef(ncid)
  if(ier.ne.NF90_NOERR) then
    ier = 34
    return
  end if

  ier = nf90_put_var(ncid,varid_x1,x1)
  if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid_x2,x2)
  if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid_f,f)
  if(ier.ne.NF90_NOERR) then
    ier = 35
    return
  end if

  ier = nf90_close(ncid)
  if(ier.ne.NF90_NOERR) ier = 36

  return
end subroutine EZspline_2NetCDF_array2


subroutine EZspline_2NetCDF_array1(n1, x1, f, filename, ier)
  use psp_precision_mod, only: fp
  use EZspline_obj   
  use netcdf
  implicit none
  integer,                intent(in)  :: n1
  real(fp), dimension(:), intent(in)  :: x1
  real(fp), dimension(:), intent(in)  :: f
  character(len=*),       intent(in)  :: filename
  integer,                intent(out) :: ier

  integer :: dimid_n1
  integer :: varid_x1
  integer :: varid_f
  integer :: ncid

  ier = nf90_create(filename,IOR(NF90_CLOBBER,NF90_64BIT_OFFSET),ncid)
  if(ier.ne.NF90_NOERR .or. ncid.eq.0) then
    ier = 33
    return
  end if

  ier = nf90_def_dim(ncid,'n1',n1,dimid_n1)
  if(ier.eq.NF90_NOERR) ier = nf90_def_var(ncid,'x1',NF90_DOUBLE,dimid_n1,varid_x1)
  if(ier.eq.NF90_NOERR) ier = nf90_def_var(ncid,'f',NF90_DOUBLE,dimid_n1,varid_f)
  if(ier.eq.NF90_NOERR) ier = nf90_enddef(ncid)
  if(ier.ne.NF90_NOERR) then
    ier = 34
    return
  end if

  ier = nf90_put_var(ncid,varid_x1,x1)
  if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid_f,f)
  if(ier.ne.NF90_NOERR) then
    ier = 35
    return
  end if

  ier = nf90_close(ncid)
  if(ier.ne.NF90_NOERR) ier = 36

  return
end subroutine EZspline_2NetCDF_array1


subroutine EZspline_2NetCDF_cloud3(n, x1, x2, x3, f, filename, ier)
  use psp_precision_mod, only: fp
  use EZspline_obj   
  use netcdf
  implicit none
  integer,                intent(in)  :: n
  real(fp), dimension(:), intent(in)  :: x1, x2, x3, f
  character(len=*),       intent(in)  :: filename
  integer,                intent(out) :: ier

  integer :: dimid_n
  integer :: varid_x1, varid_x2, varid_x3
  integer :: varid_f
  integer :: ncid

  ier = nf90_create(filename,IOR(NF90_CLOBBER,NF90_64BIT_OFFSET),ncid)
  if(ier.ne.NF90_NOERR .or. ncid.eq.0) then
    ier = 33
    return
  end if

  ier = nf90_def_dim(ncid,'n',n,dimid_n)
  if(ier.eq.NF90_NOERR) ier = nf90_def_var(ncid,'x1',NF90_DOUBLE,dimid_n,varid_x1)
  if(ier.eq.NF90_NOERR) ier = nf90_def_var(ncid,'x2',NF90_DOUBLE,dimid_n,varid_x2)
  if(ier.eq.NF90_NOERR) ier = nf90_def_var(ncid,'x3',NF90_DOUBLE,dimid_n,varid_x3)
  if(ier.eq.NF90_NOERR) ier = nf90_def_var(ncid,'f',NF90_DOUBLE,dimid_n,varid_f)
  if(ier.eq.NF90_NOERR) ier = nf90_enddef(ncid)
  if(ier.ne.NF90_NOERR) then
    ier = 34
    return
  end if

  ier = nf90_put_var(ncid,varid_x1,x1)
  if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid_x2,x2)
  if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid_x3,x3)
  if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid_f,f)
  if(ier.ne.NF90_NOERR) then
    ier = 35
    return
  end if

  ier = nf90_close(ncid)
  if(ier.ne.NF90_NOERR) ier = 36

  return
end subroutine EZspline_2NetCDF_cloud3

subroutine EZspline_2NetCDF_cloud2(n, x1, x2, f, filename, ier)
  use psp_precision_mod, only: fp
  use EZspline_obj   
  use netcdf
  implicit none
  integer,                intent(in)  :: n
  real(fp), dimension(:), intent(in)  :: x1, x2
  real(fp), dimension(:), intent(in)  :: f
  character(len=*),       intent(in)  :: filename
  integer,                intent(out) :: ier

  integer :: dimid_n
  integer :: varid_x1, varid_x2
  integer :: varid_f
  integer :: ncid

  ier = nf90_create(filename,IOR(NF90_CLOBBER,NF90_64BIT_OFFSET),ncid)
  if(ier.ne.NF90_NOERR .or. ncid.eq.0) then
    ier = 33
    return
  end if

  ier = nf90_def_dim(ncid,'n',n,dimid_n)
  if(ier.eq.NF90_NOERR) ier = nf90_def_var(ncid,'x1',NF90_DOUBLE,dimid_n,varid_x1)
  if(ier.eq.NF90_NOERR) ier = nf90_def_var(ncid,'x2',NF90_DOUBLE,dimid_n,varid_x2)
  if(ier.eq.NF90_NOERR) ier = nf90_def_var(ncid,'f',NF90_DOUBLE,dimid_n,varid_f)
  if(ier.eq.NF90_NOERR) ier = nf90_enddef(ncid)
  if(ier.ne.NF90_NOERR) then
    ier = 34
    return
  end if

  ier = nf90_put_var(ncid,varid_x1,x1)
  if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid_x2,x2)
  if(ier.eq.NF90_NOERR) ier = nf90_put_var(ncid,varid_f,f)
  if(ier.ne.NF90_NOERR) then
    ier = 35
    return
  end if

  ier = nf90_close(ncid)
  if(ier.ne.NF90_NOERR) ier = 36

  return
end subroutine EZspline_2NetCDF_cloud2

#endif
