subroutine mkspline(x,nx,fspl,ibcxmin,bcxmin,ibcxmax,bcxmax, &
     ilinx,ier)
  use precision_mod, only: fp
  !
  !  make a 2-coefficient 1d spline
  !
  !  only 2 coefficients, the data and its 2nd derivative, are needed to
  !  fully specify a spline.  See e.g. Numerical Recipies in Fortran-77
  !  (2nd edition) chapter 3, section on cubic splines.
  !
  !  input:
  !============
  implicit none
  integer i,inwk
  !============
  real(fp) :: bcxmax
  !============
  integer nx                        ! no. of data points
  real(fp) :: x(nx)                        ! x axis data, strict ascending order
  !
  !  input/output:
  real(fp) :: fspl(2,nx)                   ! f(1,*):  data in; f(2,*):  coeffs out
  !
  !     f(1,j) = f(x(j))  on input (unchanged on output)
  !     f(2,j) = f''(x(j)) (of interpolating spline) (on output).
  !
  !  ...boundary conditions...
  !
  !  input:
  !
  integer ibcxmin                   ! b.c. flag @ x=xmin=x(1)
  real(fp) :: bcxmin                       ! b.c. data @xmin
  !
  integer ibcxmax                   ! b.c. flag @ x=xmax=x(nx)
  real(fp) :: bcxcmax                      ! b.c. data @xmax
  !
  !  ibcxmin=-1 -- periodic boundary condition
  !                (bcxmin,ibcxmax,bcxmax are ignored)
  !
  !                the output spline s satisfies
  !                s'(x(1))=s'(x(nx)) ..and.. s''(x(1))=s''(x(nx))
  !
  !  if non-periodic boundary conditions are used, then the xmin and xmax
  !  boundary conditions can be specified independently:
  !
  !  ibcxmin (ibcxmax) = 0 -- this specifies a "not a knot" boundary
  !                condition, see "cubsplb.f90".  This is a common way
  !                for inferring a "good" spline boundary condition
  !                automatically from data in the vicinity of the
  !                boundary.  ... bcxmin (bcxmax) are ignored.
  !
  !  ibcxmin (ibcxmax) = 1 -- boundary condition is to have s'(x(1))
  !                ( s'(x(nx)) ) match the passed value bcxmin (bcxmax).
  !
  !  ibcxmin (ibcxmax) = 2 -- boundary condition is to have s''(x(1))
  !                ( s''(x(nx)) ) match the passed value bcxmin (bcxmax).
  !
  !  ibcxmin (ibcxmax) = 3 -- boundary condition is to have s'(x(1))=0.0
  !                ( s'(x(nx))=0.0 )
  !
  !  ibcxmin (ibcxmax) = 4 -- boundary condition is to have s''(x(1))=0.0
  !                ( s''(x(nx))=0.0 )
  !
  !  ibcxmin (ibcxmax) = 5 -- boundary condition is to have s'(x(1))
  !                ( s'(x(nx)) ) match the 1st divided difference
  !                e.g. at x(1):  d(1)/h(1), where
  !                           d(j)=f(1,j+1)-f(1,j)
  !                           h(j)=x(j+1)-x(j)
  !
  !  ibcxmin (ibcxmax) = 6 -- BC is to have s''(x(1)) ( s''(x(nx)) )
  !                match the 2nd divided difference
  !                e.g. at x(1):
  !                     e(1) = [d(2)/h(2) - d(1)/h(1)]/(0.5*(h(1)+h(2)))
  !
  !  ibcxmin (ibcxmax) = 7 -- BC is to have s'''(x(1)) ( s'''(x(nx)) )
  !                match the 3rd divided difference
  !                e.g. at x(1): [e(2)-e(1)]/(0.33333*(h(1)+h(2)+h(3)))
  !
  !  output:
  !
  integer ilinx                     ! =1: hint, x axis is ~evenly spaced
  !
  !  let dx[avg] = (x(nx)-x(1))/(nx-1)
  !  let dx[j] = x(j+1)-x(j), for all j satisfying 1.le.j.lt.nx
  !
  !  if for all such j, abs(dx[j]-dx[avg]).le.(1.0e-3*dx[avg]) then
  !  ilinx=1 is returned, indicating the data is (at least nearly)
  !  evenly spaced.  Even spacing is useful, for speed of zone lookup,
  !  when evaluating a spline.
  !
  !  if the even spacing condition is not satisfied, ilinx=2 is returned.
  !
  integer ier                       ! exit code, 0=OK
  !
  !  an error code is returned if the x axis is not strict ascending,
  !  or if nx.lt.4, or if an invalid boundary condition specification was
  !  input.
  !
  !------------------------------------
  !
  !  this routine calls traditional 4 coefficient spline software, and
  !  translates the result to 2 coefficient form.
  !
  !  this could be done more efficiently but we decided out of conservatism
  !  to use the traditional software.
  !
  !------------------------------------
  !  workspaces -- f90 dynamically allocated
  !
  real(fp), dimension(:,:), allocatable :: fspl4 ! traditional 4-spline
  real(fp), dimension(:), allocatable :: wk ! cspline workspace
  !
  !------------------------------------
  !
  allocate(fspl4(4,nx),wk(nx))
  !
  !  make the traditional call
  !
  do i=1,nx
     fspl4(1,i)=fspl(1,i)
     fspl(2,i)=0.0_fp  ! for now
  end do
  !
  inwk=nx
  !
  !  boundary conditions imposed by cspline...
  !
  call cspline(x,nx,fspl4,ibcxmin,bcxmin,ibcxmax,bcxmax, &
       wk,inwk,ilinx,ier)
  !
  if(ier.eq.0) then
     !
     !  copy the output -- careful of end point.
     !
     do i=1,nx-1
        fspl(2,i)=2.0_fp*fspl4(3,i)
     end do
     fspl(2,nx)=2.0_fp*fspl4(3,nx-1) + &
          (x(nx)-x(nx-1))*6.0_fp*fspl4(4,nx-1)
  end if
  !
  deallocate(fspl4,wk)
  !
  return
end subroutine mkspline
