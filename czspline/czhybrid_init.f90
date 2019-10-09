!---------------------------------------------------------------------------
! This code was developed at Tech-X (www.txcorp.com). It is free for any one
! to use but comes with no warranty whatsoever. Use at your own risk. 
! Thanks for reporting bugs to pletzer@txcorp.com. 
!---------------------------------------------------------------------------
! Initialization of ezspline

#include "czspline_handle_size.h"

subroutine czhybrid_init2(handle, n1, n2, hspline, BCS1, BCS2, ier)
  use ezspline_obj
  use ezspline
  use czspline_pointer_types
  implicit none
  integer, intent(inout) :: handle(_ARRSZ)
  integer, intent(in) :: n1, n2
  integer, intent(in) :: hspline(2)
  integer, intent(in) :: BCS1(2), BCS2(2)
  integer, intent(out) :: ier
  type(czspline2) :: self
  allocate(self % ptr)
  call ezhybrid_init(self % ptr, n1, n2, hspline, ier, BCS1, BCS2)
  handle = 0
  handle = transfer(self, handle)
end subroutine czhybrid_init2

subroutine czhybrid_init3(handle, n1, n2, n3, hspline, BCS1, BCS2, BCS3, ier)
  use ezspline_obj
  use ezspline
  use czspline_pointer_types
  implicit none
  integer, intent(inout) :: handle(_ARRSZ)
  integer, intent(in) :: n1, n2, n3
  integer, intent(in) :: hspline(3)
  integer, intent(in) :: BCS1(2), BCS2(2), BCS3(2)
  integer, intent(out) :: ier
  type(czspline3) :: self
  allocate(self % ptr)
  call ezhybrid_init(self % ptr, n1, n2, n3, hspline, ier, BCS1, BCS2, BCS3)
  handle = 0
  handle = transfer(self, handle)
end subroutine czhybrid_init3
