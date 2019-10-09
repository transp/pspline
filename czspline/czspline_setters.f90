!---------------------------------------------------------------------------
! This code was developed at Tech-X (www.txcorp.com). It is free for any one
! to use but comes with no warranty whatsoever. Use at your own risk. 
! Thanks for reporting bugs to pletzer@txcorp.com. 
!---------------------------------------------------------------------------
! Member setters

#include "czspline_handle_size.h"

! x1, x2, x3
subroutine czspline_set_axes1(handle, k1, x1, ier)
  use precision_mod, only: fp
  use ezspline_obj
  use ezspline
  use czspline_pointer_types
  implicit none
  integer, intent(inout) :: handle(_ARRSZ)
  integer, intent(in)  :: k1
  real(fp), intent(in) :: x1(k1)
  integer, intent(out) :: ier
  type(czspline1) :: self
  self = transfer(handle, self)
  ier = 0
  self % ptr % x1 = x1
end subroutine czspline_set_axes1

subroutine czspline_set_axes2(handle, k1, k2, x1, x2, ier)
  use precision_mod, only: fp
  use ezspline_obj
  use ezspline
  use czspline_pointer_types
  implicit none
  integer, intent(inout) :: handle(_ARRSZ)
  integer, intent(in)  :: k1, k2
  real(fp), intent(in) :: x1(k1), x2(k2)
  integer, intent(out) :: ier
  type(czspline2) :: self
  self = transfer(handle, self)
  ier = 0
  self % ptr % x1 = x1
  self % ptr % x2 = x2
end subroutine czspline_set_axes2

subroutine czspline_set_axes3(handle, k1, k2, k3, x1, x2, x3, ier)
  use precision_mod, only: fp
  use ezspline_obj
  use ezspline
  use czspline_pointer_types
  implicit none
  integer, intent(inout) :: handle(_ARRSZ)
  integer, intent(in)  :: k1, k2, k3
  real(fp), intent(in) :: x1(k1), x2(k2), x3(k3)
  integer, intent(out) :: ier
  type(czspline3) :: self
  self = transfer(handle, self)
  ier = 0
  self % ptr % x1 = x1
  self % ptr % x2 = x2
  self % ptr % x3 = x3
end subroutine czspline_set_axes3

! isHermite
subroutine czspline_set_ishermite1(handle, flag, ier)
  use ezspline_obj
  use ezspline
  use czspline_pointer_types
  implicit none
  integer, intent(inout) :: handle(_ARRSZ)
  integer, intent(in)  :: flag
  integer, intent(out) :: ier
  type(czspline1) :: self
  self = transfer(handle, self)
  ier = 0
  self % ptr % isHermite = flag
end subroutine czspline_set_ishermite1

subroutine czspline_set_ishermite2(handle, flag, ier)
  use ezspline_obj
  use ezspline
  use czspline_pointer_types
  implicit none
  integer, intent(inout) :: handle(_ARRSZ)
  integer, intent(in)  :: flag
  integer, intent(out) :: ier
  type(czspline2) :: self
  self = transfer(handle, self)
  ier = 0
  self % ptr % isHermite = flag
end subroutine czspline_set_ishermite2

subroutine czspline_set_ishermite3(handle, flag, ier)
  use ezspline_obj
  use ezspline
  use czspline_pointer_types
  implicit none
  integer, intent(inout) :: handle(_ARRSZ)
  integer, intent(in)  :: flag
  integer, intent(out) :: ier
  type(czspline3) :: self
  self = transfer(handle, self)
  ier = 0
  self % ptr % isHermite = flag
end subroutine czspline_set_ishermite3

! Boundary types
subroutine czspline_set_bctypes1(handle, bctype1, ier)
  use ezspline_obj
  use ezspline
  use czspline_pointer_types
  implicit none
  integer, intent(inout) :: handle(_ARRSZ)
  integer, intent(in)    :: bctype1(2)
  integer, intent(out) :: ier
  type(czspline1) :: self
  self = transfer(handle, self)
  ier = 0
  self % ptr % ibctype1 = bctype1
end subroutine czspline_set_bctypes1

subroutine czspline_set_bctypes2(handle, bctype1, bctype2, ier)
  use ezspline_obj
  use ezspline
  use czspline_pointer_types
  implicit none
  integer, intent(inout) :: handle(_ARRSZ)
  integer, intent(in)    :: bctype1(2), bctype2(2)
  integer, intent(out) :: ier
  type(czspline2) :: self
  self = transfer(handle, self)
  ier = 0
  self % ptr % ibctype1 = bctype1
  self % ptr % ibctype2 = bctype2
end subroutine czspline_set_bctypes2

subroutine czspline_set_bctypes3(handle, bctype1, bctype2, bctype3, ier)
  use ezspline_obj
  use ezspline
  use czspline_pointer_types
  implicit none
  integer, intent(inout) :: handle(_ARRSZ)
  integer, intent(in)    :: bctype1(2), bctype2(2), bctype3(2)
  integer, intent(out) :: ier
  type(czspline3) :: self
  self = transfer(handle, self)
  ier = 0
  self % ptr % ibctype1 = bctype1
  self % ptr % ibctype2 = bctype2
  self % ptr % ibctype3 = bctype3
end subroutine czspline_set_bctypes3

! Boundary conditions bcval1min, bcval1max, ...
subroutine czspline_set_bcvals1(handle, bcval1, ier)
  use precision_mod, only: fp
  use ezspline_obj
  use ezspline
  use czspline_pointer_types
  implicit none
  integer, intent(inout) :: handle(_ARRSZ)
  real(fp), intent(in)   :: bcval1(2)
  integer, intent(out) :: ier
  type(czspline1) :: self
  self = transfer(handle, self)
  ier = 0
  self % ptr % bcval1min = bcval1(1)
  self % ptr % bcval1max = bcval1(2)
end subroutine czspline_set_bcvals1

subroutine czspline_set_bcvals2(handle, bcval1, bcval2, ier)
  use precision_mod, only: fp
  use ezspline_obj
  use ezspline
  use czspline_pointer_types
  implicit none
  integer, intent(inout) :: handle(_ARRSZ)
  real(fp), intent(in)   :: bcval1(2), bcval2(2)
  integer, intent(out) :: ier
  type(czspline2) :: self
  self = transfer(handle, self)
  ier = 0
  self % ptr % bcval1min = bcval1(1)
  self % ptr % bcval1max = bcval1(2)
  self % ptr % bcval2min = bcval2(1)
  self % ptr % bcval2max = bcval2(2)
end subroutine czspline_set_bcvals2

subroutine czspline_set_bcvals3(handle, bcval1, bcval2, bcval3, ier)
  use precision_mod, only: fp
  use ezspline_obj
  use ezspline
  use czspline_pointer_types
  implicit none
  integer, intent(inout) :: handle(_ARRSZ)
  real(fp), intent(in)   :: bcval1(2), bcval2(2), bcval3(2)
  integer, intent(out) :: ier
  type(czspline3) :: self
  self = transfer(handle, self)
  ier = 0
  self % ptr % bcval1min = bcval1(1)
  self % ptr % bcval1max = bcval1(2)
  self % ptr % bcval2min = bcval2(1)
  self % ptr % bcval2max = bcval2(2)
  self % ptr % bcval3min = bcval3(1)
  self % ptr % bcval3max = bcval3(2)
end subroutine czspline_set_bcvals3
