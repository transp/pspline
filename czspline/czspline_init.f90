!---------------------------------------------------------------------------
! This code was developed at Tech-X (www.txcorp.com). It is free for any one
! to use but comes with no warranty whatsoever. Use at your own risk. 
! Thanks for reporting bugs to pletzer@txcorp.com. 
!---------------------------------------------------------------------------
! Initialization of ezspline

#include "czspline_handle_size.h"

subroutine czspline_init1(handle, n1, BCS1, ier)
  use ezspline_obj
  use ezspline
  use czspline_pointer_types
  implicit none
  integer, intent(inout) :: handle(_ARRSZ)
  integer, intent(in) :: n1
  integer, intent(in) :: BCS1(2)
  integer, intent(out) :: ier
  type(czspline1) :: self
  allocate(self % ptr)
  call ezspline_init(self % ptr, n1, BCS1, ier)
  handle = 0
  handle = transfer(self, handle)
end subroutine czspline_init1

subroutine czspline_init2(handle, n1, n2, BCS1, BCS2, ier)
  use ezspline_obj
  use ezspline
  use czspline_pointer_types
  implicit none
  integer, intent(inout) :: handle(_ARRSZ)
  integer, intent(in) :: n1, n2
  integer, intent(in) :: BCS1(2), BCS2(2)
  integer, intent(out) :: ier
  type(czspline2) :: self
  allocate(self % ptr)
  call ezspline_init(self % ptr, n1, n2, BCS1, BCS2, ier)
  handle = 0
  handle = transfer(self, handle)
end subroutine czspline_init2

subroutine czspline_init3(handle, n1, n2, n3, BCS1, BCS2, BCS3, ier)
  use ezspline_obj
  use ezspline
  use czspline_pointer_types
  implicit none
  integer, intent(inout) :: handle(_ARRSZ)
  integer, intent(in) :: n1, n2, n3
  integer, intent(in) :: BCS1(2), BCS2(2), BCS3(2)
  integer, intent(out) :: ier
  type(czspline3) :: self
  allocate(self % ptr)
  call ezspline_init(self % ptr, n1, n2, n3, BCS1, BCS2, BCS3, ier)
  handle = 0
  handle = transfer(self, handle)
end subroutine czspline_init3
