subroutine EZspline_isGridRegular1(spline_o, ier)
  use precision_mod, only: fp
  use EZspline_type
  use EZspline_obj
  implicit none
  type(EZspline1) :: spline_o
  ! ier:
  ! 14=x1 grid is not strictly increasing
  integer, intent(out) :: ier
  integer i

  ier = 0
  do i=1, spline_o%n1-1
     if(spline_o%x1(i+1) <= spline_o%x1(i)) then
        ier = 14
        return
     end if
  end do

end subroutine EZspline_isGridRegular1

subroutine EZspline_isGridRegular2(spline_o, ier)
  use precision_mod, only: fp
  use EZspline_obj
  implicit none
  type(EZspline2) :: spline_o
  ! ier:
  ! 14=x1 grid is not strictly increasing
  ! 15=x2 grid is not strictly increasing
  integer, intent(out) :: ier
  integer i

  ier = 0
  do i=1, spline_o%n1-1
     if(spline_o%x1(i+1) <= spline_o%x1(i)) then
        ier = 14
        return
     end if
  end do
  do i=1, spline_o%n2-1
     if(spline_o%x2(i+1) <= spline_o%x2(i)) then
        ier = 15
        return
     end if
  end do

end subroutine EZspline_isGridRegular2


subroutine EZspline_isGridRegular3(spline_o, ier)
  use precision_mod, only: fp
  use EZspline_obj
  implicit none
  type(EZspline3) :: spline_o
  ! ier:
  ! 14=x1 grid is not strictly increasing
  ! 15=x2 grid is not strictly increasing
  ! 16=x3 grid is not strictly increasing
  integer, intent(out) :: ier
  integer i

  ier = 0
  do i=1, spline_o%n1-1
     if(spline_o%x1(i+1) <= spline_o%x1(i)) then
        ier = 14
        return
     end if
  end do
  do i=1, spline_o%n2-1
     if(spline_o%x2(i+1) <= spline_o%x2(i)) then
        ier = 15
        return
     end if
  end do
  do i=1, spline_o%n3-1
     if(spline_o%x3(i+1) <= spline_o%x3(i)) then
        ier = 16
        return
     end if
  end do

end subroutine EZspline_isGridRegular3
