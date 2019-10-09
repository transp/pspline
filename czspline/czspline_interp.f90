!---------------------------------------------------------------------------
! This code was developed at Tech-X (www.txcorp.com). It is free for any one
! to use but comes with no warranty whatsoever. Use at your own risk. 
! Thanks for reporting bugs to pletzer@txcorp.com. 
!---------------------------------------------------------------------------
! Interpolate

#include "czspline_handle_size.h"

! point interpolation
subroutine czspline_interp1(handle, p1, f, ier)
  use precision_mod, only: fp
  use ezspline_obj
  use ezspline
  use czspline_pointer_types
  implicit none
  integer, intent(inout) :: handle(_ARRSZ)
  real(fp), intent(in) :: p1
  real(fp), intent(out) :: f
  integer, intent(out) :: ier
  type(czspline1) :: self
  self = transfer(handle, self)
  call ezspline_interp(self % ptr, p1, f, ier)
end subroutine czspline_interp1

! cloud interpolation
subroutine czspline_interp1_cloud(handle, k, p1, f, ier)
  use precision_mod, only: fp
  use ezspline_obj
  use ezspline
  use czspline_pointer_types
  implicit none
  integer, intent(inout) :: handle(_ARRSZ)
  integer, intent(in)  :: k
  real(fp), intent(in) :: p1(k)
  real(fp), intent(out) :: f(k)
  integer, intent(out) :: ier
  type(czspline1) :: self
  self = transfer(handle, self)
  call ezspline_interp(self % ptr, k, p1, f, ier)
end subroutine czspline_interp1_cloud

! array interpolation
subroutine czspline_interp1_array(handle, k1, p1, f, ier)
  use precision_mod, only: fp
  use ezspline_obj
  use ezspline
  use czspline_pointer_types
  implicit none
  integer, intent(inout) :: handle(_ARRSZ)
  integer, intent(in)  :: k1 
  real(fp), intent(in) :: p1(k1)
  real(fp), intent(out) :: f(k1)
  integer, intent(out) :: ier
  type(czspline1) :: self
  self = transfer(handle, self)
  call ezspline_interp(self % ptr, k1, p1, f, ier)
end subroutine czspline_interp1_array

! point interpolation
subroutine czspline_interp2(handle, p1, p2, f, ier)
  use precision_mod, only: fp
  use ezspline_obj
  use ezspline
  use czspline_pointer_types
  implicit none
  integer, intent(inout) :: handle(_ARRSZ)
  real(fp), intent(in) :: p1, p2
  real(fp), intent(out) :: f
  integer, intent(out) :: ier
  type(czspline2) :: self
  self = transfer(handle, self)
  call ezspline_interp(self % ptr, p1, p2, f, ier)
end subroutine czspline_interp2

! cloud interpolation
subroutine czspline_interp2_cloud(handle, k, p1, p2, f, ier)
  use precision_mod, only: fp
  use ezspline_obj
  use ezspline
  use czspline_pointer_types
  implicit none
  integer, intent(inout) :: handle(_ARRSZ)
  integer, intent(in)  :: k
  real(fp), intent(in) :: p1(k), p2(k)
  real(fp), intent(out) :: f(k)
  integer, intent(out) :: ier
  type(czspline2) :: self
  self = transfer(handle, self)
  call ezspline_interp(self % ptr, k, p1, p2, f, ier)
end subroutine czspline_interp2_cloud

! array interpolation
subroutine czspline_interp2_array(handle, k1, k2, p1, p2, f, ier)
  use precision_mod, only: fp
  use ezspline_obj
  use ezspline
  use czspline_pointer_types
  implicit none
  integer, intent(inout) :: handle(_ARRSZ)
  integer, intent(in)  :: k1, k2
  real(fp), intent(in) :: p1(k1), p2(k2)
  real(fp), intent(out) :: f(k1, k2)
  integer, intent(out) :: ier
  type(czspline2) :: self
  self = transfer(handle, self)
  call ezspline_interp(self % ptr, k1, k2, p1, p2, f, ier)
end subroutine czspline_interp2_array

! point interpolation
subroutine czspline_interp3(handle, p1, p2, p3, f, ier)
  use precision_mod, only: fp
  use ezspline_obj
  use ezspline
  use czspline_pointer_types
  implicit none
  integer, intent(inout) :: handle(_ARRSZ)
  real(fp), intent(in) :: p1, p2, p3
  real(fp), intent(out) :: f
  integer, intent(out) :: ier
  type(czspline3) :: self
  self = transfer(handle, self)
  call ezspline_interp(self % ptr, p1, p2, p3, f, ier)
end subroutine czspline_interp3

! cloud interpolation
subroutine czspline_interp3_cloud(handle, k, p1, p2, p3, f, ier)
  use precision_mod, only: fp
  use ezspline_obj
  use ezspline
  use czspline_pointer_types
  implicit none
  integer, intent(inout) :: handle(_ARRSZ)
  integer, intent(in)  :: k
  real(fp), intent(in) :: p1(k), p2(k), p3(k)
  real(fp), intent(out) :: f(k)
  integer, intent(out) :: ier
  type(czspline3) :: self
  self = transfer(handle, self)
  call ezspline_interp(self % ptr, k, p1, p2, p3, f, ier)
end subroutine czspline_interp3_cloud

! array interpolation
subroutine czspline_interp3_array(handle, k1, k2, k3, p1, p2, p3, f, ier)
  use precision_mod, only: fp
  use ezspline_obj
  use ezspline
  use czspline_pointer_types
  implicit none
  integer, intent(inout) :: handle(_ARRSZ)
  integer, intent(in)  :: k1, k2, k3
  real(fp), intent(in) :: p1(k1), p2(k2), p3(k3)
  real(fp), intent(out) :: f(k1, k2 ,k3)
  integer, intent(out) :: ier
  type(czspline3) :: self
  self = transfer(handle, self)
  call ezspline_interp(self % ptr, k1, k2, k3, p1, p2, p3, f, ier)
end subroutine czspline_interp3_array
