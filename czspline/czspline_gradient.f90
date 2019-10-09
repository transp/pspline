!---------------------------------------------------------------------------
! This code was developed at Tech-X (www.txcorp.com). It is free for any one
! to use but comes with no warranty whatsoever. Use at your own risk. 
! Thanks for reporting bugs to pletzer@txcorp.com. 
!---------------------------------------------------------------------------
! Compute gradients
  
#include "czspline_handle_size.h"

! point
subroutine czspline_gradient1(handle, p1, df, ier)
  use precision_mod, only: fp
  use ezspline_obj
  use ezspline
  use czspline_pointer_types
  implicit none
  integer, intent(inout) :: handle(_ARRSZ)
  real(fp), intent(in) :: p1
  real(fp), intent(out) :: df
  integer, intent(out) :: ier
  type(czspline1) :: self
  self = transfer(handle, self)
  call ezspline_gradient(self % ptr, p1, df, ier)
end subroutine czspline_gradient1

! cloud
subroutine czspline_gradient1_cloud(handle, k, p1, df, ier)
  use precision_mod, only: fp
  use ezspline_obj
  use ezspline
  use czspline_pointer_types
  implicit none
  integer, intent(inout) :: handle(_ARRSZ)
  integer, intent(in)  :: k
  real(fp), intent(in) :: p1(k)
  real(fp), intent(out) :: df(k)
  integer, intent(out) :: ier
  type(czspline1) :: self
  self = transfer(handle, self)
  call ezspline_gradient(self % ptr, k, p1, df, ier)
end subroutine czspline_gradient1_cloud

! array
subroutine czspline_gradient1_array(handle, k1, p1, df, ier)
  use precision_mod, only: fp
  use ezspline_obj
  use ezspline
  use czspline_pointer_types
  implicit none
  integer, intent(inout) :: handle(_ARRSZ)
  integer, intent(in)  :: k1
  real(fp), intent(in) :: p1(k1)
  real(fp), intent(out) :: df(k1)
  integer, intent(out) :: ier
  type(czspline1) :: self
  self = transfer(handle, self)
  call ezspline_gradient(self % ptr, k1, p1, df, ier)
end subroutine czspline_gradient1_array

! point
subroutine czspline_gradient2(handle, p1, p2, df, ier)
  use precision_mod, only: fp
  use ezspline_obj
  use ezspline
  use czspline_pointer_types
  implicit none
  integer, intent(inout) :: handle(_ARRSZ)
  real(fp), intent(in) :: p1, p2
  real(fp), intent(out) :: df(2)
  integer, intent(out) :: ier
  type(czspline2) :: self
  self = transfer(handle, self)
  call ezspline_gradient(self % ptr, p1, p2, df, ier)
end subroutine czspline_gradient2

! cloud
subroutine czspline_gradient2_cloud(handle, k, p1, p2, df, ier)
  use precision_mod, only: fp
  use ezspline_obj
  use ezspline
  use czspline_pointer_types
  implicit none
  integer, intent(inout) :: handle(_ARRSZ)
  integer, intent(in)  :: k
  real(fp), intent(in) :: p1(k), p2(k)
  real(fp), intent(out) :: df(k, 2)
  integer, intent(out) :: ier
  type(czspline2) :: self
  self = transfer(handle, self)
  call ezspline_gradient(self % ptr, k, p1, p2, df, ier)
end subroutine czspline_gradient2_cloud

! array
subroutine czspline_gradient2_array(handle, k1, k2, p1, p2, df, ier)
  use precision_mod, only: fp
  use ezspline_obj
  use ezspline
  use czspline_pointer_types
  implicit none
  integer, intent(inout) :: handle(_ARRSZ)
  integer, intent(in)  :: k1, k2
  real(fp), intent(in) :: p1(k1), p2(k2)
  real(fp), intent(out) :: df(k1, k2, 2)
  integer, intent(out) :: ier
  type(czspline2) :: self
  self = transfer(handle, self)
  call ezspline_gradient(self % ptr, k1, k2, p1, p2, df, ier)
end subroutine czspline_gradient2_array

! point
subroutine czspline_gradient3(handle, p1, p2, p3, df, ier)
  use precision_mod, only: fp
  use ezspline_obj
  use ezspline
  use czspline_pointer_types
  implicit none
  integer, intent(inout) :: handle(_ARRSZ)
  real(fp), intent(in) :: p1, p2, p3
  real(fp), intent(out) :: df(3)
  integer, intent(out) :: ier
  type(czspline3) :: self
  self = transfer(handle, self)
  call ezspline_gradient(self % ptr, p1, p2, p3, df, ier)
end subroutine czspline_gradient3

! cloud
subroutine czspline_gradient3_cloud(handle, k, p1, p2, p3, df, ier)
  use precision_mod, only: fp
  use ezspline_obj
  use ezspline
  use czspline_pointer_types
  implicit none
  integer, intent(inout) :: handle(_ARRSZ)
  integer, intent(in)  :: k
  real(fp), intent(in) :: p1(k), p2(k), p3(k)
  real(fp), intent(out) :: df(k, 3)
  integer, intent(out) :: ier
  type(czspline3) :: self
  self = transfer(handle, self)
  call ezspline_gradient(self % ptr, k, p1, p2, p3, df, ier)
end subroutine czspline_gradient3_cloud

! array
subroutine czspline_gradient3_array(handle, k1, k2, k3, p1, p2, p3, df, ier)
  use precision_mod, only: fp
  use ezspline_obj
  use ezspline
  use czspline_pointer_types
  implicit none
  integer, intent(inout) :: handle(_ARRSZ)
  integer, intent(in)  :: k1, k2, k3
  real(fp), intent(in) :: p1(k1), p2(k2), p3(k3)
  real(fp), intent(out) :: df(k1, k2, k3, 3)
  integer, intent(out) :: ier
  type(czspline3) :: self
  self = transfer(handle, self)
  call ezspline_gradient(self % ptr, k1, k2, k3, p1, p2, p3, df, ier)
end subroutine czspline_gradient3_array


