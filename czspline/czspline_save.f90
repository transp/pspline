!---------------------------------------------------------------------------
! This code was developed at Tech-X (www.txcorp.com). It is free for any one
! to use but comes with no warranty whatsoever. Use at your own risk. 
! Thanks for reporting bugs to pletzer@txcorp.com. 
!---------------------------------------------------------------------------
! Save object in file

#include "czspline_handle_size.h"

subroutine czspline_save1(handle, filename, ier)
  use ezspline_obj
  use ezspline
  use czspline_pointer_types
  implicit none
  integer, intent(inout) :: handle(_ARRSZ)
  character(len=*), intent(in) :: filename
  integer, intent(out) :: ier
  type(czspline1) :: self
  self = transfer(handle, self)
  call ezspline_save(self % ptr, filename, ier)
end subroutine czspline_save1

subroutine czspline_save2(handle, filename, ier)
  use ezspline_obj
  use ezspline
  use czspline_pointer_types
  implicit none
  integer, intent(inout) :: handle(_ARRSZ)
  character(len=*), intent(in) :: filename
  integer, intent(out) :: ier
  type(czspline2) :: self
  self = transfer(handle, self)
  call ezspline_save(self % ptr, filename, ier)
end subroutine czspline_save2

subroutine czspline_save3(handle, filename, ier)
  use ezspline_obj
  use ezspline
  use czspline_pointer_types
  implicit none
  integer, intent(inout) :: handle(_ARRSZ)
  character(len=*), intent(in) :: filename
  integer, intent(out) :: ier
  type(czspline3) :: self
  self = transfer(handle, self)
  call ezspline_save(self % ptr, filename, ier)
end subroutine czspline_save3
