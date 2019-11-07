!  pspline library test routine
!
program pspltest
  use psp_precision_mod, only: fp
  real(fp) :: zdum = 0.0_fp
  call pspltsub(' ',zdum)
  stop
end program pspltest
