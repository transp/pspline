!  utility routines for pspltsub

subroutine psp_tolsum(refval,tolval,sum)
  use precision_mod, only: fp

  !  difference will be detectable
  !  refval =~ 1.0; tolval =~1.0e-6

  implicit none
  real(fp) :: refval, tolval  ! input
  real(fp) :: sum             ! output

  sum = refval + tolval*tolval

  return
end subroutine psp_tolsum
