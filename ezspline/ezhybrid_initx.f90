!  EZhybrid routines callable without a module interface...

subroutine EZhybrid_init2_x(spline_o, n1, n2, hspline, ier)

  !  init call (set array size from object)

  use EZspline_obj
  use EZspline
  implicit NONE

  type(EZspline2) spline_o
  integer, intent(in) :: n1, n2
  integer, intent(in) :: hspline(2)
  integer, intent(out) :: ier

  call EZhybrid_init(spline_o, n1, n2, hspline, ier)

end subroutine EZhybrid_init2_x

subroutine EZhybrid_init3_x(spline_o, n1, n2, n3, hspline, ier)

  !  init call (set array size from object)

  use EZspline_obj
  use EZspline
  implicit NONE

  type(EZspline3) spline_o
  integer, intent(in) :: n1, n2, n3
  integer, intent(in) :: hspline(3)
  integer, intent(out) :: ier

  call EZhybrid_init(spline_o, n1, n2, n3, hspline, ier)

end subroutine EZhybrid_init3_x
