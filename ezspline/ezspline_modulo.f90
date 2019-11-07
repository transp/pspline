!
! map point into (xmin, xmax) cell when boundary conditions are periodic.

subroutine EZspline_modulo1(spline_o, p1, ier)
  use psp_precision_mod, only: fp
  use EZspline_obj
  implicit none
  type(EZspline1) spline_o
  real(fp) :: p1 ! the location
  integer, intent(out) :: ier

  ier = 0
  if(spline_o%ibctype1(1)==-1 .OR. spline_o%ibctype1(2)==-1) then
     p1 = MOD( p1, spline_o%x1max - spline_o%x1min )
  end if

end subroutine EZspline_modulo1

subroutine EZspline_modulo_array1(spline_o, k1, p1, ier)
  use psp_precision_mod, only: fp
  use EZspline_obj
  implicit none
  type(EZspline1) spline_o
  integer, intent(in) :: k1 
  real(fp) :: p1(k1) ! the points
  integer, intent(out) :: ier

  ier = 0
  if(spline_o%ibctype1(1)==-1 .OR. spline_o%ibctype1(2)==-1) then
     p1(1:k1) = MOD( p1(1:k1), spline_o%x1max - spline_o%x1min )
  end if

end subroutine EZspline_modulo_array1


subroutine EZspline_modulo2(spline_o, p1, p2, ier)
  use psp_precision_mod, only: fp
  use EZspline_obj
  implicit none
  type(EZspline2) spline_o
  real(fp) :: p1, p2 ! the location
  integer, intent(out) :: ier

  ier = 0
  if(spline_o%ibctype1(1)==-1 .OR. spline_o%ibctype1(2)==-1) then
     p1 = MOD( p1, spline_o%x1max - spline_o%x1min )
  end if

  if(spline_o%ibctype2(1)==-1 .OR. spline_o%ibctype2(2)==-1) then
     p2 = MOD( p2, spline_o%x2max - spline_o%x2min )
  end if

end subroutine EZspline_modulo2

subroutine EZspline_modulo_array2(spline_o, k1, k2, p1, p2, ier)
  use psp_precision_mod, only: fp
  use EZspline_obj
  implicit none
  type(EZspline2) spline_o
  integer, intent(in) :: k1, k2
  real(fp) :: p1(k1), p2(k2) ! the points
  integer, intent(out) :: ier

  ier = 0
  if(spline_o%ibctype1(1)==-1 .OR. spline_o%ibctype1(2)==-1) then
     p1(1:k1) = MOD( p1(1:k1), spline_o%x1max - spline_o%x1min )
  end if

  if(spline_o%ibctype2(1)==-1 .OR. spline_o%ibctype2(2)==-1) then
     p2(1:k2) = MOD( p2(1:k2), spline_o%x2max - spline_o%x2min )
  end if

end subroutine EZspline_modulo_array2

subroutine EZspline_modulo_cloud2(spline_o, k, p1, p2, ier)
  use psp_precision_mod, only: fp
  use EZspline_obj
  implicit none
  type(EZspline2) spline_o
  integer, intent(in) :: k
  real(fp) :: p1(k), p2(k) ! the points
  integer, intent(out) :: ier

  ier = 0
  if(spline_o%ibctype1(1)==-1 .OR. spline_o%ibctype1(2)==-1) then
     p1(1:k) = MOD( p1(1:k), spline_o%x1max - spline_o%x1min )
  end if

  if(spline_o%ibctype2(1)==-1 .OR. spline_o%ibctype2(2)==-1) then
     p2(1:k) = MOD( p2(1:k), spline_o%x2max - spline_o%x2min )
  end if

end subroutine EZspline_modulo_cloud2


subroutine EZspline_modulo3(spline_o, p1, p2, p3, ier)
  use psp_precision_mod, only: fp
  use EZspline_obj
  implicit none
  type(EZspline3) spline_o
  real(fp) :: p1, p2, p3 ! the location
  integer, intent(out) :: ier

  ier = 0
  if(spline_o%ibctype1(1)==-1 .OR. spline_o%ibctype1(2)==-1) then
     p1 = MOD( p1, spline_o%x1max - spline_o%x1min )
  end if

  if(spline_o%ibctype2(1)==-1 .OR. spline_o%ibctype2(2)==-1) then
     p2 = MOD( p2, spline_o%x2max - spline_o%x2min )
  end if

  if(spline_o%ibctype3(1)==-1 .OR. spline_o%ibctype3(2)==-1) then
     p3 = MOD( p3, spline_o%x3max - spline_o%x3min )
  end if

end subroutine EZspline_modulo3

subroutine EZspline_modulo_array3(spline_o, k1, k2, k3, p1, p2, p3, ier)
  use psp_precision_mod, only: fp
  use EZspline_obj
  implicit none
  type(EZspline3) spline_o
  integer, intent(in) :: k1, k2, k3
  real(fp) :: p1(k1), p2(k2), p3(k3) ! the points
  integer, intent(out) :: ier

  ier = 0
  if(spline_o%ibctype1(1)==-1 .OR. spline_o%ibctype1(2)==-1) then
     p1(1:k1) = MOD( p1(1:k1), spline_o%x1max - spline_o%x1min )
  end if

  if(spline_o%ibctype2(1)==-1 .OR. spline_o%ibctype2(2)==-1) then
     p2(1:k2) = MOD( p2(1:k2), spline_o%x2max - spline_o%x2min )
  end if

  if(spline_o%ibctype3(1)==-1 .OR. spline_o%ibctype3(2)==-1) then
     p3(1:k3) = MOD( p3(1:k3), spline_o%x3max - spline_o%x3min )
  end if

end subroutine EZspline_modulo_array3

subroutine EZspline_modulo_cloud3(spline_o, k, p1, p2, p3, ier)
  use psp_precision_mod, only: fp
  use EZspline_obj
  implicit none
  type(EZspline3) spline_o
  integer, intent(in) :: k
  real(fp) :: p1(k), p2(k), p3(k) ! the points
  integer, intent(out) :: ier

  ier = 0
  if(spline_o%ibctype1(1)==-1 .OR. spline_o%ibctype1(2)==-1) then
     p1(1:k) = MOD( p1(1:k), spline_o%x1max - spline_o%x1min )
  end if

  if(spline_o%ibctype2(1)==-1 .OR. spline_o%ibctype2(2)==-1) then
     p2(1:k) = MOD( p2(1:k), spline_o%x2max - spline_o%x2min )
  end if

  if(spline_o%ibctype3(1)==-1 .OR. spline_o%ibctype3(2)==-1) then
     p3(1:k) = MOD( p3(1:k), spline_o%x3max - spline_o%x3min )
  end if

end subroutine EZspline_modulo_cloud3
