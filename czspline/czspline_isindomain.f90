!---------------------------------------------------------------------------
! This code was developed at Tech-X (www.txcorp.com). It is free for any one
! to use but comes with no warranty whatsoever. Use at your own risk. 
! Thanks for reporting bugs to pletzer@txcorp.com. 
!---------------------------------------------------------------------------
! Return error if any of the target points lie outside the domain

#include "czspline_handle_size.h"

!point
subroutine czspline_isindomain1(handle, p1, ier)
  use precision_mod, only: fp
  use ezspline_obj
  use ezspline
  use czspline_pointer_types
  implicit none
  integer, intent(inout) :: handle(_ARRSZ)
  real(fp), intent(in) :: p1
  integer, intent(out) :: ier
  type(czspline1) :: self
  self = transfer(handle, self)
  call ezspline_isindomain(self % ptr, p1, ier)
end subroutine czspline_isindomain1

!cloud
subroutine czspline_isindomain1_cloud(handle, k, p1, ier)
  use precision_mod, only: fp
  use ezspline_obj
  use ezspline
  use czspline_pointer_types
  implicit none
  integer, intent(inout) :: handle(_ARRSZ)
  integer, intent(in) :: k
  real(fp), intent(in) :: p1(k)
  integer, intent(out) :: ier
  type(czspline1) :: self
  self = transfer(handle, self)
  call ezspline_isindomain(self % ptr, k, p1, ier)
end subroutine czspline_isindomain1_cloud

!array
subroutine czspline_isindomain1_array(handle, k1, p1, ier)
  use precision_mod, only: fp
  use ezspline_obj
  use ezspline
  use czspline_pointer_types
  implicit none
  integer, intent(inout) :: handle(_ARRSZ)
  integer, intent(in) :: k1
  real(fp), intent(in) :: p1(k1)
  integer, intent(out) :: ier
  type(czspline1) :: self
  self = transfer(handle, self)
  call ezspline_isindomain(self % ptr, k1, p1, ier)
end subroutine czspline_isindomain1_array

!point
subroutine czspline_isindomain2(handle, p1, p2, ier)
  use ezspline_obj
  use ezspline
  use czspline_pointer_types
  implicit none
  integer, intent(inout) :: handle(_ARRSZ)
  real(fp), intent(in) :: p1, p2
  integer, intent(out) :: ier
  type(czspline2) :: self
  self = transfer(handle, self)
  call ezspline_isindomain(self % ptr, p1, p2, ier)
end subroutine czspline_isindomain2

!cloud
subroutine czspline_isindomain2_cloud(handle, k, p1, p2, ier)
  use precision_mod, only: fp
  use ezspline_obj
  use ezspline
  use czspline_pointer_types
  implicit none
  integer, intent(inout) :: handle(_ARRSZ)
  integer, intent(in) :: k
  real(fp), intent(in) :: p1(k), p2(k)
  integer, intent(out) :: ier
  type(czspline2) :: self
  self = transfer(handle, self)
  call ezspline_isindomain(self % ptr, k, p1, p2, ier)
end subroutine czspline_isindomain2_cloud

!array
subroutine czspline_isindomain2_array(handle, k1, k2, p1, p2, ier)
  use precision_mod, only: fp
  use ezspline_obj
  use ezspline
  use czspline_pointer_types
  implicit none
  integer, intent(inout) :: handle(_ARRSZ)
  integer, intent(in) :: k1, k2
  real(fp), intent(in)  :: p1(k1), p2(k2)
  integer, intent(out) :: ier
  type(czspline2) :: self
  self = transfer(handle, self)
  call ezspline_isindomain(self % ptr, k1, k2, p1, p2, ier)
end subroutine czspline_isindomain2_array

!point
subroutine czspline_isindomain3(handle, p1, p2, p3, ier)
  use precision_mod, only: fp
  use ezspline_obj
  use ezspline
  use czspline_pointer_types
  implicit none
  integer, intent(inout) :: handle(_ARRSZ)
  real(fp), intent(in) :: p1, p2, p3
  integer, intent(out) :: ier
  type(czspline3) :: self
  self = transfer(handle, self)
  call ezspline_isindomain(self % ptr, p1, p2, p3, ier)
end subroutine czspline_isindomain3

!cloud
subroutine czspline_isindomain3_cloud(handle, k, p1, p2, p3, ier)
  use precision_mod, only: fp
  use ezspline_obj
  use ezspline
  use czspline_pointer_types
  implicit none
  integer, intent(inout) :: handle(_ARRSZ)
  integer, intent(in) :: k
  real(fp), intent(in) :: p1(k), p2(k), p3(k)
  integer, intent(out) :: ier
  type(czspline3) :: self
  self = transfer(handle, self)
  call ezspline_isindomain(self % ptr, k, p1, p2, p3, ier)
end subroutine czspline_isindomain3_cloud

!array
subroutine czspline_isindomain3_array(handle, k1, k2, k3, p1, p2, p3, ier)
  use precision_mod, only: fp
  use ezspline_obj
  use ezspline
  use czspline_pointer_types
  implicit none
  integer, intent(inout) :: handle(_ARRSZ)
  integer, intent(in) :: k1, k2, k3
  real(fp), intent(in) :: p1(k1), p2(k2), p3(k3)
  integer, intent(out) :: ier
  type(czspline3) :: self
  self = transfer(handle, self)
  call ezspline_isindomain(self % ptr, k1, k2, k3, p1, p2, p3, ier)
end subroutine czspline_isindomain3_array
