!
! 1-D
!


subroutine EZspline_isInDomain1(spline_o, p1, ier)
  use psp_precision_mod, only: fp
  use EZspline_obj
  implicit none
  type(EZspline1) :: spline_o
  real(fp), intent(in) :: p1

  ! ier:
  ! 6=out of interval p1 < min(x1)
  ! 7=out of interval p1 > max(x1)
  integer, intent(out) :: ier

  real(fp), parameter :: small = 1.e-10_fp
  real(fp) :: eps1

  ier = 0

  eps1 = small*(spline_o%x1(spline_o%n1) - spline_o%x1(1))


  if(p1 < spline_o%x1(          1)-eps1) then
     ier = 6
     return
  end if
  if(p1 > spline_o%x1(spline_o%n1)+eps1) then
     ier = 7
     return
  end if

end subroutine EZspline_isInDomain1

subroutine EZspline_isInDomain1_array(spline_o, k1, p1, ier)
  use psp_precision_mod, only: fp
  use EZspline_obj
  type(EZspline1) :: spline_o
  integer, intent(in) :: k1
  real(fp), intent(in) :: p1(k1)
  integer, intent(out) :: ier

  real(fp), parameter :: small = 1.e-10_fp
  real(fp) :: eps1

  ier = 0

  eps1 = small*(spline_o%x1(spline_o%n1) - spline_o%x1(1))

  if(minval(p1(1:k1)) < spline_o%x1(          1)-eps1) then
     ier = 6
     return
  end if
  if(maxval(p1(1:k1)) > spline_o%x1(spline_o%n1)+eps1) then
     ier = 7
     return
  end if

end subroutine EZspline_isInDomain1_array



!!
!! 2-D
!!


subroutine EZspline_isInDomain2(spline_o, p1, p2, ier)
  use psp_precision_mod, only: fp
  use EZspline_obj
  implicit none
  type(EZspline2) :: spline_o
  real(fp), intent(in) :: p1, p2

  ! ier:
  ! 6=out of interval p1 < min(x1)
  ! 7=out of interval p1 > max(x1)
  ! 8=out of interval p2 < min(x2)
  ! 9=out of interval p2 > max(x2)
  integer, intent(out) :: ier

  real(fp), parameter :: small = 1.e-10_fp
  real(fp) :: eps1, eps2

  ier = 0

  eps1 = small*(spline_o%x1(spline_o%n1) - spline_o%x1(1))
  eps2 = small*(spline_o%x2(spline_o%n2) - spline_o%x2(1))


  if(p1 < spline_o%x1(          1)-eps1) then
     ier = 6
     return
  end if
  if(p1 > spline_o%x1(spline_o%n1)+eps1) then
     ier = 7
     return
  end if
  if(p2 < spline_o%x2(          1)-eps2) then
     ier = 8
     return
  end if
  if(p2 > spline_o%x2(spline_o%n2)+eps2) then
     ier = 9
     return
  end if

end subroutine EZspline_isInDomain2

subroutine EZspline_isInDomain2_array(spline_o, k1, k2, p1, p2, ier)
  use psp_precision_mod, only: fp
  use EZspline_obj
  type(EZspline2) :: spline_o
  integer, intent(in) :: k1, k2
  real(fp), intent(in) :: p1(k1), p2(k2)
  integer, intent(out) :: ier

  real(fp), parameter :: small = 1.e-10_fp
  real(fp) :: eps1, eps2

  ier = 0

  eps1 = small*(spline_o%x1(spline_o%n1) - spline_o%x1(1))
  eps2 = small*(spline_o%x2(spline_o%n2) - spline_o%x2(1))

  if(minval(p1(1:k1)) < spline_o%x1(          1)-eps1) then
     ier = 6
     return
  end if
  if(maxval(p1(1:k1)) > spline_o%x1(spline_o%n1)+eps1) then
     ier = 7
     return
  end if
  if(minval(p2(1:k2)) < spline_o%x2(          1)-eps2) then
     ier = 8
     return
  end if
  if(maxval(p2(1:k2)) > spline_o%x2(spline_o%n2)+eps2) then
     ier = 9
     return
  end if

end subroutine EZspline_isInDomain2_array

subroutine EZspline_isInDomain2_cloud(spline_o, k, p1, p2, ier)
  use psp_precision_mod, only: fp
  use EZspline_obj
  type(EZspline2) :: spline_o
  integer, intent(in) :: k
  real(fp), intent(in) :: p1(k), p2(k)
  integer, intent(out) :: ier

  real(fp), parameter :: small = 1.e-10_fp
  real(fp) :: eps1, eps2

  ier = 0

  eps1 = small*(spline_o%x1(spline_o%n1) - spline_o%x1(1))
  eps2 = small*(spline_o%x2(spline_o%n2) - spline_o%x2(1))

  if(minval(p1(1:k)) < spline_o%x1(          1)-eps1) then
     ier = 6
     return
  end if
  if(maxval(p1(1:k)) > spline_o%x1(spline_o%n1)+eps1) then
     ier = 7
     return
  end if
  if(minval(p2(1:k)) < spline_o%x2(          1)-eps2) then
     ier = 8
     return
  end if
  if(maxval(p2(1:k)) > spline_o%x2(spline_o%n2)+eps2) then
     ier = 9
     return
  end if

end subroutine EZspline_isInDomain2_cloud



!!!
!!! 3-D
!!!


subroutine EZspline_isInDomain3(spline_o, p1, p2, p3, ier)
  use psp_precision_mod, only: fp
  use EZspline_obj
  implicit none
  type(EZspline3) :: spline_o
  real(fp), intent(in) :: p1, p2, p3

  ! ier:
  ! 6=out of interval p1 < min(x1)
  ! 7=out of interval p1 > max(x1)
  ! 8=out of interval p2 < min(x2)
  ! 9=out of interval p2 > max(x2)
  ! 10=out of interval p3 < min(x3)
  ! 11=out of interval p3 > max(x3)
  integer, intent(out) :: ier

  real(fp), parameter :: small = 1.e-10_fp
  real(fp) :: eps1, eps2, eps3

  ier = 0

  eps1 = small*(spline_o%x1(spline_o%n1) - spline_o%x1(1))
  eps2 = small*(spline_o%x2(spline_o%n2) - spline_o%x2(1))
  eps3 = small*(spline_o%x3(spline_o%n3) - spline_o%x3(1))


  if(p1 < spline_o%x1(          1)-eps1) then
     ier = 6
     return
  end if
  if(p1 > spline_o%x1(spline_o%n1)+eps1) then
     ier = 7
     return
  end if
  if(p2 < spline_o%x2(          1)-eps2) then
     ier = 8
     return
  end if
  if(p2 > spline_o%x2(spline_o%n2)+eps2) then
     ier = 9
     return
  end if
  if(p3 < spline_o%x3(          1)-eps3) then
     ier = 10
     return
  end if
  if(p3 > spline_o%x3(spline_o%n3)+eps3) then
     ier = 11
     return
  end if

end subroutine EZspline_isInDomain3

subroutine EZspline_isInDomain3_array(spline_o, k1, k2, k3, p1, p2, p3, ier)
  use psp_precision_mod, only: fp
  use EZspline_obj
  type(EZspline3) :: spline_o
  integer, intent(in) :: k1, k2, k3
  real(fp), intent(in) :: p1(k1), p2(k2), p3(k3)
  integer, intent(out) :: ier

  real(fp), parameter :: small = 1.e-10_fp
  real(fp) :: eps1, eps2, eps3

  ier = 0

  eps1 = small*(spline_o%x1(spline_o%n1) - spline_o%x1(1))
  eps2 = small*(spline_o%x2(spline_o%n2) - spline_o%x2(1))
  eps3 = small*(spline_o%x3(spline_o%n3) - spline_o%x3(1))

  if(minval(p1(1:k1)) < spline_o%x1(          1)-eps1) then
     ier = 6
     return
  end if
  if(maxval(p1(1:k1)) > spline_o%x1(spline_o%n1)+eps1) then
     ier = 7
     return
  end if
  if(minval(p2(1:k2)) < spline_o%x2(          1)-eps2) then
     ier = 8
     return
  end if
  if(maxval(p2(1:k2)) > spline_o%x2(spline_o%n2)+eps2) then
     ier = 9
     return
  end if
  if(minval(p3(1:k3)) < spline_o%x3(          1)-eps3) then
     ier = 10
     return
  end if
  if(maxval(p3(1:k3)) > spline_o%x3(spline_o%n3)+eps3) then
     ier = 11
     return
  end if

end subroutine EZspline_isInDomain3_array

subroutine EZspline_isInDomain3_cloud(spline_o, k, p1, p2, p3, ier)
  use psp_precision_mod, only: fp
  use EZspline_obj
  type(EZspline3) :: spline_o
  integer, intent(in) :: k
  real(fp), intent(in) :: p1(k), p2(k), p3(k)
  integer, intent(out) :: ier

  real(fp), parameter :: small = 1.e-10_fp
  real(fp) :: eps1, eps2, eps3

  ier = 0

  eps1 = small*(spline_o%x1(spline_o%n1) - spline_o%x1(1))
  eps2 = small*(spline_o%x2(spline_o%n2) - spline_o%x2(1))
  eps3 = small*(spline_o%x3(spline_o%n3) - spline_o%x3(1))

  if(minval(p1(1:k)) < spline_o%x1(          1)-eps1) then
     ier = 6
     return
  end if
  if(maxval(p1(1:k)) > spline_o%x1(spline_o%n1)+eps1) then
     ier = 7
     return
  end if
  if(minval(p2(1:k)) < spline_o%x2(          1)-eps2) then
     ier = 8
     return
  end if
  if(maxval(p2(1:k)) > spline_o%x2(spline_o%n2)+eps2) then
     ier = 9
     return
  end if
  if(minval(p3(1:k)) < spline_o%x3(          1)-eps3) then
     ier = 10
     return
  end if
  if(maxval(p3(1:k)) > spline_o%x3(spline_o%n3)+eps3) then
     ier = 11
     return
  end if
end subroutine EZspline_isInDomain3_cloud
