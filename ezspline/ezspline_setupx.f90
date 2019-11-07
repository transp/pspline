subroutine ezspline_setup1_x(spline_o, f, ier)
  use psp_precision_mod, only: fp

  !  setup call (set array size from object)

  use EZspline_obj
  use EZspline
  implicit NONE

  type(EZspline1) :: spline_o
  real(fp) :: f(size(spline_o%fspl,2))
  integer, intent(out) :: ier

  call ezspline_setup(spline_o, f, ier)

end subroutine ezspline_setup1_x

subroutine ezspline_setup2_x(spline_o, f, ier)
  use psp_precision_mod, only: fp

  !  setup call (set array size from object)

  use EZspline_obj
  use EZspline
  implicit NONE

  type(EZspline2) :: spline_o
  real(fp) :: f(size(spline_o%fspl,2),size(spline_o%fspl,3))
  integer, intent(out) :: ier

  call ezspline_setup(spline_o, f, ier)

end subroutine ezspline_setup2_x

subroutine ezspline_setup3_x(spline_o, f, ier)
  use psp_precision_mod, only: fp

  !  setup call (set array size from object)

  use EZspline_obj
  use EZspline
  implicit NONE

  type(EZspline3) :: spline_o
  real(fp) :: f(size(spline_o%fspl,2),size(spline_o%fspl,3), &
       size(spline_o%fspl,4))
  integer, intent(out) :: ier

  call ezspline_setup(spline_o, f, ier)

end subroutine ezspline_setup3_x
