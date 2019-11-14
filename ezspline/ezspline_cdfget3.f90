#ifdef _NETCDF
subroutine ezspline_cdfget3(ncid,zname,fspl,i1,i2,i3,ier)
  use psp_precision_mod, only: fp
  use netcdf
  implicit NONE

  ! (due to ezspline rank limitation) read 4d object as 3d object
  ! use f77 interface style, should prevent unnecessary array copy.

  integer,                       intent(in)  :: ncid  ! opened NetCDF file
  character(len=*),              intent(in)  :: zname ! name to us for writing
  integer,                       intent(in)  :: i1    ! spline data dimensions
  integer,                       intent(in)  :: i2    ! spline data dimensions
  integer,                       intent(in)  :: i3    ! spline data dimensions
  real(fp), dimension(i1,i2,i3), intent(out) :: fspl  ! spline data and coefficients
  integer,                       intent(out) :: ier   ! status returned from NetCDF
  integer :: varid

  ier = nf90_inq_varid(ncid,trim(zname),varid)
  if(ier.eq.NF90_NOERR) ier = nf90_get_var(ncid,varid,fspl)

  return
end subroutine ezspline_cdfget3
#endif
