subroutine gridpc2( &
     x_newgrid,nx_new,y_newgrid,ny_new,f_new,if1, &
     nx,xpkg,ny,ypkg,fspl,inf1,iwarn,ier)
  use precision_mod, only: fp
  !
  !  regrid a piecewise linear function f defined vs. x,y as in xpkg,ypkg
  !  to a new grid, given by x_newgrid, y_newgrid
  !
  !  set warning flag if the range x_newgrid, y_newgrid
  !  exceeds the range of the original xpkg,ypkg.
  !
  !  (xpkg,ypkg -- axis data, see genxpkg subroutine)
  !
  !  input:
  !
  !============
  implicit none
  integer ny_new,nx_new,iy
  !============
  real(fp) :: x_newgrid(nx_new)            ! new x grid
  real(fp) :: y_newgrid(ny_new)            ! new y grid
  !
  !  output:
  !
  integer if1                       ! 1st dimension of f_new
  real(fp) :: f_new(if1,ny_new)            ! f evaluated on this grid
  !
  !  input:
  !
  integer nx                        ! size of old grid
  real(fp) :: xpkg(nx,4)                   ! old grid "package"
  integer ny                        ! size of old grid
  real(fp) :: ypkg(ny,4)                   ! old grid "package"
  !
  integer inf1                      ! array dimension
  real(fp) :: fspl(inf1,ny)                ! f array
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
  real(fp) :: ytmp(nx_new)
  integer ict(4)
  !
  data ict/1,0,0,0/
  !
  !--------------------------------------------
  !
  do iy=1,ny_new
     ytmp=y_newgrid(iy)
     call vecpc2(ict,nx_new,x_newgrid,ytmp,nx_new,f_new(1,iy), &
          nx,xpkg,ny,ypkg,fspl,inf1,iwarn,ier)
  end do
  !
  return
end subroutine gridpc2
