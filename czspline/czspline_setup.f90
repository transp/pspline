!---------------------------------------------------------------------------
! This code was developed at Tech-X (www.txcorp.com). It is free for any one
! to use but comes with no warranty whatsoever. Use at your own risk. 
! Thanks for reporting bugs to pletzer@txcorp.com. 
!---------------------------------------------------------------------------
! Compute spline coefficients

#include "czspline_handle_size.h"

subroutine czspline_setup1(handle, n1, f, ier)
  use precision_mod, only: fp
  use ezspline_obj
  use ezspline
  use czspline_pointer_types
  implicit none
  integer, intent(inout) :: handle(_ARRSZ)
  integer, intent(in) :: n1
  real(fp), intent(in) :: f(n1)
  integer, intent(out) :: ier
  type(czspline1) :: self
  self = transfer(handle, self)
  call ezspline_setup(self % ptr, f, ier)
end subroutine czspline_setup1

subroutine czspline_setup2(handle, n1, n2, f, ier)
  use precision_mod, only: fp
  use ezspline_obj
  use ezspline
  use czspline_pointer_types
  implicit none
  integer, intent(inout) :: handle(_ARRSZ)
  integer, intent(in) :: n1, n2
  real(fp), intent(in) :: f(n1, n2)
  integer, intent(out) :: ier
  type(czspline2) :: self
  self = transfer(handle, self)
  call ezspline_setup(self % ptr, f, ier)
end subroutine czspline_setup2

subroutine czspline_setup3(handle, n1, n2, n3, f, ier)
  use precision_mod, only: fp
  use ezspline_obj
  use ezspline
  use czspline_pointer_types
  implicit none
  integer, intent(inout) :: handle(_ARRSZ)
  integer, intent(in) :: n1, n2, n3
  real(fp), intent(in) :: f(n1, n2, n3)
  integer, intent(out) :: ier
  type(czspline3) :: self
  self = transfer(handle, self)
  call ezspline_setup(self % ptr, f, ier)
end subroutine czspline_setup3
