subroutine EZspline_free1(spline_o, ier)
  use psp_precision_mod, only: fp
  use EZspline_obj
  implicit none
  type(EZspline1) spline_o
  ! ier:
  ! 101= warning, spline object was never allocated
  integer, intent(out) :: ier
  integer ifail

  ier = 0
  if(.not.EZspline_allocated(spline_o)) ier=101

  deallocate(spline_o%x1, stat=ifail)
  deallocate(spline_o%fspl, stat=ifail)
  deallocate(spline_o%x1pkg, stat=ifail)

  call EZspline_preInit(spline_o)

  spline_o%n1 = 0

  return
end subroutine EZspline_free1


subroutine EZspline_free2(spline_o, ier)
  use psp_precision_mod, only: fp
  use EZspline_obj
  implicit none
  type(EZspline2) spline_o
  ! ier:
  ! 101= warning, spline object was never allocated
  integer, intent(out) :: ier
  integer ifail

  ier = 0
  if(.not.EZspline_allocated(spline_o)) ier=101

  deallocate(spline_o%x1, stat=ifail)
  deallocate(spline_o%x2, stat=ifail)
  deallocate(spline_o%fspl, stat=ifail)
  deallocate(spline_o%bcval1min, stat=ifail)
  deallocate(spline_o%bcval1max, stat=ifail)
  deallocate(spline_o%bcval2min, stat=ifail)
  deallocate(spline_o%bcval2max, stat=ifail)
  deallocate(spline_o%x1pkg, stat=ifail)
  deallocate(spline_o%x2pkg, stat=ifail)

  call EZspline_preInit(spline_o)

  spline_o%n1 = 0
  spline_o%n2 = 0

  return
end subroutine EZspline_free2


subroutine EZspline_free3(spline_o, ier)
  use psp_precision_mod, only: fp
  use EZspline_obj
  implicit none
  type(EZspline3) spline_o
  ! ier:
  ! 101= warning, spline object was never allocated
  integer, intent(out) :: ier
  integer ifail

  ier = 0
  if(.not.EZspline_allocated(spline_o)) ier=101

  deallocate(spline_o%x1, stat=ifail)
  deallocate(spline_o%x2, stat=ifail)
  deallocate(spline_o%x3, stat=ifail)
  deallocate(spline_o%fspl, stat=ifail)
  deallocate(spline_o%bcval1min, stat=ifail)
  deallocate(spline_o%bcval1max, stat=ifail)
  deallocate(spline_o%bcval2min, stat=ifail)
  deallocate(spline_o%bcval2max, stat=ifail)
  deallocate(spline_o%bcval3min, stat=ifail)
  deallocate(spline_o%bcval3max, stat=ifail)
  deallocate(spline_o%x1pkg, stat=ifail)
  deallocate(spline_o%x2pkg, stat=ifail)
  deallocate(spline_o%x3pkg, stat=ifail)

  call EZspline_preInit(spline_o)

  spline_o%n1 = 0
  spline_o%n2 = 0
  spline_o%n3 = 0

  return
end subroutine EZspline_free3
