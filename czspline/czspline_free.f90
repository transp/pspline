!---------------------------------------------------------------------------
! This code was developed at Tech-X (www.txcorp.com). It is free for any one
! to use but comes with no warranty whatsoever. Use at your own risk. 
! Thanks for reporting bugs to pletzer@txcorp.com. 
!---------------------------------------------------------------------------
! Finalization of ezspline
  
#include "czspline_handle_size.h"

subroutine czspline_free1(handle, ier)
  use ezspline_obj
  use ezspline
  use czspline_pointer_types
  implicit none
  integer, intent(inout) :: handle(_ARRSZ)
  integer, intent(out) :: ier
  type(czspline1) :: self
  self = transfer(handle, self)
  call ezspline_free(self %ptr, ier)
  deallocate(self %ptr)
  handle = 0
end subroutine czspline_free1

subroutine czspline_free2(handle, ier)
  use ezspline_obj
  use ezspline
  use czspline_pointer_types
  implicit none
  integer, intent(inout) :: handle(_ARRSZ)
  integer, intent(out) :: ier
  type(czspline2) :: self
  self = transfer(handle, self)
  call ezspline_free(self %ptr, ier)
  deallocate(self %ptr)
  handle = 0
end subroutine czspline_free2

subroutine czspline_free3(handle, ier)
  use ezspline_obj
  use ezspline
  use czspline_pointer_types
  implicit none
  integer, intent(inout) :: handle(_ARRSZ)
  integer, intent(out) :: ier
  type(czspline3) :: self
  self = transfer(handle, self)
  call ezspline_free(self %ptr, ier)
  deallocate(self %ptr)
  handle = 0
end subroutine czspline_free3

