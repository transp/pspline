!---------------------------------------------------------------------------
! This code was developed at Tech-X (www.txcorp.com). It is free for any one
! to use but comes with no warranty whatsoever. Use at your own risk. 
! Thanks for reporting bugs to pletzer@txcorp.com. 
!---------------------------------------------------------------------------
! Load object from file

#include "czspline_handle_size.h"

subroutine czspline_load1(handle, filename, ier)
  use ezspline_obj
  use ezspline
  use czspline_pointer_types
  implicit none
  integer, intent(inout) :: handle(_ARRSZ)
  character(len=*), intent(in) :: filename
  integer, intent(out) :: ier
  type(czspline1) :: self
  allocate(self % ptr)
  write(6,*) "before ezspline_load" ; flush(6)
  call ezspline_load(self % ptr, filename, ier)
  write(6,*) "after ezspline_load" ; flush(6)
  handle = 0
  handle = transfer(self, handle)
end subroutine czspline_load1

subroutine czspline_load2(handle, filename, ier)
  use ezspline_obj
  use ezspline
  use czspline_pointer_types
  implicit none
  integer, intent(inout) :: handle(_ARRSZ)
  character(len=*), intent(in) :: filename
  integer, intent(out) :: ier
  type(czspline2) :: self
  allocate(self % ptr)
  call ezspline_load(self % ptr, filename, ier, ' ')
  handle = 0
  handle = transfer(self, handle)
end subroutine czspline_load2

subroutine czspline_load3(handle, filename, ier)
  use ezspline_obj
  use ezspline
  use czspline_pointer_types
  implicit none
  integer, intent(inout) :: handle(_ARRSZ)
  character(len=*), intent(in) :: filename
  integer, intent(out) :: ier
  type(czspline3) :: self
  allocate(self % ptr)
  call ezspline_load(self % ptr, filename, ier)
  handle = 0
  handle = transfer(self, handle)
end subroutine czspline_load3
