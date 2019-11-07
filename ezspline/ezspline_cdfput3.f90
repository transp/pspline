#ifdef _EZCDF
subroutine ezspline_cdfput3(ncid,zname,fspl,idim1,idim2,idim3,ifail)
  use psp_precision_mod, only: fp
  use EZcdf
  implicit NONE

  ! (due to ezspline rank limitation) write 4d object as 3d object
  ! use f77 interface style, should prevent unnecessary array copy.

  integer, intent(in) :: ncid             ! opened NetCDF file
  character(len=*),intent(in) :: zname       ! name to us for writing
  integer, intent(in) :: idim1,idim2,idim3   ! spline data dimensions
  real(fp), intent(in) :: fspl(idim1,idim2,idim3)  ! spline data & coefficients
  integer, intent(out) :: ifail           ! status retruned from NetCDF

  call cdfPutVar(ncid, trim(zname), fspl, ifail)

  return
end subroutine ezspline_cdfput3
#endif
