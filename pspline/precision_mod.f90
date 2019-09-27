module precision_mod
  use iso_c_binding, only: c_int, c_float, c_double
  integer, parameter :: ip = c_int
#ifdef _SINGLE_PRECISION
  integer, parameter :: fp = c_float
#else
  integer, parameter :: fp = c_double
#endif
end module precision_mod
