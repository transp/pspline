!---------------------------------------------------------------------------
! This code was developed at Tech-X (www.txcorp.com). It is free for any one
! to use but comes with no warranty whatsoever. Use at your own risk. 
! Thanks for reporting bugs to pletzer@txcorp.com. 
!---------------------------------------------------------------------------
! Compute derivatives

#include "czspline_handle_size.h"

! point
subroutine czspline_derivative1(handle, i1, p1, f, ier)
  use precision_mod, only: fp
  use ezspline_obj
  use ezspline
  use czspline_pointer_types
  implicit none
  integer, intent(inout) :: handle(_ARRSZ)
  integer, intent(in) :: i1
  real(fp), intent(in) :: p1
  real(fp), intent(out) :: f
  integer, intent(out) :: ier
  type(czspline1) :: self
  self = transfer(handle, self)
  call ezspline_derivative(self % ptr, i1, p1, f, ier)
end subroutine czspline_derivative1

! cloud
subroutine czspline_derivative1_cloud(handle, i1, k, p1, f, ier)
  use precision_mod, only: fp
  use ezspline_obj
  use ezspline
  use czspline_pointer_types
  implicit none
  integer, intent(inout) :: handle(_ARRSZ)
  integer, intent(in) :: i1
  integer, intent(in) :: k
  real(fp), intent(in) :: p1(k)
  real(fp), intent(out) :: f(k)
  integer, intent(out) :: ier
  type(czspline1) :: self
  self = transfer(handle, self)
  call ezspline_derivative(self % ptr, i1, k, p1, f, ier)
end subroutine czspline_derivative1_cloud

!array
subroutine czspline_derivative1_array(handle, i1, k1, p1, f, ier)
  use precision_mod, only: fp
  use ezspline_obj
  use ezspline
  use czspline_pointer_types
  implicit none
  integer, intent(inout) :: handle(_ARRSZ)
  integer, intent(in) :: i1
  integer, intent(in) :: k1
  real(fp), intent(in) :: p1(k1)
  real(fp), intent(out) :: f(k1)
  integer, intent(out) :: ier
  type(czspline1) :: self
  self = transfer(handle, self)
  call ezspline_derivative(self % ptr, i1, k1, p1, f, ier)
end subroutine czspline_derivative1_array

! point
subroutine czspline_derivative2(handle, i1, i2, p1, p2, f, ier)
  use precision_mod, only: fp
  use ezspline_obj
  use ezspline
  use czspline_pointer_types
  implicit none
  integer, intent(inout) :: handle(_ARRSZ)
  integer, intent(in) :: i1, i2
  real(fp), intent(in) :: p1, p2
  real(fp), intent(out) :: f
  integer, intent(out) :: ier
  type(czspline2) :: self
  self = transfer(handle, self)
  call ezspline_derivative(self % ptr, i1, i2, p1, p2, f, ier)
end subroutine czspline_derivative2

! cloud
subroutine czspline_derivative2_cloud(handle, i1, i2, k, p1, p2, f, ier)
  use precision_mod, only: fp
  use ezspline_obj
  use ezspline
  use czspline_pointer_types
  implicit none
  integer, intent(inout) :: handle(_ARRSZ)
  integer, intent(in) :: i1, i2
  integer, intent(in) :: k
  real(fp), intent(in) :: p1(k), p2(k)
  real(fp), intent(out) :: f(k)
  integer, intent(out) :: ier
  type(czspline2) :: self
  self = transfer(handle, self)
  call ezspline_derivative(self % ptr, i1, i2, k, p1, p2, f, ier)
end subroutine czspline_derivative2_cloud

!array
subroutine czspline_derivative2_array(handle, i1, i2, k1, k2, p1, p2, f, ier)
  use precision_mod, only: fp
  use ezspline_obj
  use ezspline
  use czspline_pointer_types
  implicit none
  integer, intent(inout) :: handle(_ARRSZ)
  integer, intent(in) :: i1, i2
  integer, intent(in) :: k1, k2
  real(fp), intent(in) :: p1(k1), p2(k2)
  real(fp), intent(out) :: f(k1, k2)
  integer, intent(out) :: ier
  type(czspline2) :: self
  self = transfer(handle, self)
  call ezspline_derivative(self % ptr, i1, i2, k1, k2, p1, p2, f, ier)
end subroutine czspline_derivative2_array

! point
subroutine czspline_derivative3(handle, i1, i2, i3, p1, p2, p3, f, ier)
  use precision_mod, only: fp
  use ezspline_obj
  use ezspline
  use czspline_pointer_types
  implicit none
  integer, intent(inout) :: handle(_ARRSZ)
  integer, intent(in) :: i1, i2, i3
  real(fp), intent(in) :: p1, p2, p3
  real(fp), intent(out) :: f
  integer, intent(out) :: ier
  type(czspline3) :: self
  self = transfer(handle, self)
  call ezspline_derivative(self % ptr, i1, i2, i3, p1, p2, p3, f, ier)
end subroutine czspline_derivative3

! cloud
subroutine czspline_derivative3_cloud(handle, i1, i2, i3, k, p1, p2, p3, f, ier)
  use precision_mod, only: fp
  use ezspline_obj
  use ezspline
  use czspline_pointer_types
  implicit none
  integer, intent(inout) :: handle(_ARRSZ)
  integer, intent(in) :: i1, i2, i3
  integer, intent(in) :: k
  real(fp), intent(in) :: p1(k), p2(k), p3(k)
  real(fp), intent(out) :: f(k)
  integer, intent(out) :: ier
  type(czspline3) :: self
  self = transfer(handle, self)
  call ezspline_derivative(self % ptr, i1, i2, i3, k, p1, p2, p3, f, ier)
end subroutine czspline_derivative3_cloud

!array
subroutine czspline_derivative3_array(handle, i1, i2, i3, k1, k2, k3, p1, p2, p3, f, ier)
  use precision_mod, only: fp
  use ezspline_obj
  use ezspline
  use czspline_pointer_types
  implicit none
  integer, intent(inout) :: handle(_ARRSZ)
  integer, intent(in) :: i1, i2, i3
  integer, intent(in) :: k1, k2, k3
  real(fp), intent(in) :: p1(k1), p2(k2), p3(k3)
  real(fp), intent(out) :: f(k1, k2, k3)
  integer, intent(out) :: ier
  type(czspline3) :: self
  self = transfer(handle, self)
  call ezspline_derivative(self % ptr, i1, i2, i3, k1, k2, k3, p1, p2, p3, f, ier)
end subroutine czspline_derivative3_array
