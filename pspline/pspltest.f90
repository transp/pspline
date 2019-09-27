!  pspline library test routine
!
program pspltest
  use precision_mod, only: fp
  real(fp) :: zdum = 0.0_fp
  call pspltsub(' ',zdum)
  stop
end program pspltest
