subroutine gridpc1(x_newgrid,nx_new,f_new,nx,xpkg,fspl, &
     iwarn,ier)
  use precision_mod, only: fp
  !
  !  regrid a piecewise linear function f defined vs. x as in xpkg
  !  to a new grid, given by x_newgrid.
  !
  !  set warning flag if the range x_newgrid exceeds the range of the
  !  original xpkg.
  !
  !  (xpkg -- see genxpkg subroutine)
  !
  !  input:
  !
  !============
  implicit none
  integer nx_new
  !============
  real(fp) :: x_newgrid(nx_new)            ! new grid
  !
  !  output:
  !
  real(fp) :: f_new(nx_new)                ! f evaluated on this grid
  !
  !  input:
  !
  integer nx                        ! size of old grid
  real(fp) :: xpkg(nx,4)                   ! old grid "package"
  real(fp) :: fspl(nx)                     ! the function data
  !
  !  output:
  !  condition codes, =0 for normal exit
  !
  integer iwarn                     ! =1 if new grid points out of range
  integer ier                       ! =1 if there is an argument error
  !
  !--------------------------------------------
  !  local
  !
  integer ict(2)
  !
  ict(1) = 1
  ict(2) = 0
  !
  !--------------------------------------------
  !
  call vecpc1(ict,nx_new,x_newgrid,nx_new,f_new,nx,xpkg,fspl, &
       iwarn,ier)
  !
  return
end subroutine gridpc1
