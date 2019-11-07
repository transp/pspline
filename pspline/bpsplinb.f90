!  bpsplinb -- dmc 30 May 1996
!
!  set up coefficients for bicubic spline with following BC's:
!  * LHS and RHS BC under user control (see comments)
!  * derivatives periodic in second coordinate (use pspline.f90)
!
subroutine bpsplinb(x,inx,th,inth,fspl,inf3, &
     ibcxmin,bcxmin,ibcxmax,bcxmax, &
     wk,nwk,ilinx,ilinth,ier)
  use psp_precision_mod, only: fp
  !
  implicit none
  integer inx,inth,inf3,nwk
  real(fp) :: x(inx),th(inth),fspl(4,4,inf3,inth),wk(nwk)
  integer ibcxmin,ibcxmax
  real(fp) :: bcxmin(inth),bcxmax(inth)
  !
  !  input:
  !    x(1...inx) -- abscissae, first dimension of data
  !   th(1...inth) -- abscissae, second (periodic) dimension of data
  !   fspl(1,1,1..inx,1..inth) -- function values
  !   inf3 -- fspl dimensioning, inf3.ge.inx required.
  !
  !  boundary conditions input:
  !   ibcxmin -- indicator for boundary condition at x(1):
  !    bcxmin(...) -- boundary condition data
  !     =0 -- use "not a knot", bcxmin(...) ignored
  !     =1 -- match slope, specified at x(1),th(ith) by bcxmin(ith)
  !     =2 -- match 2nd derivative, specified at x(1),th(ith) by bcxmin(ith)
  !     =3 -- boundary condition is slope=0 (df/dx=0) at x(1), all th(j)
  !     =4 -- boundary condition is d2f/dx2=0 at x(1), all th(j)
  !     =5 -- match 1st derivative to 1st divided difference
  !     =6 -- match 2nd derivative to 2nd divided difference
  !     =7 -- match 3rd derivative to 3rd divided difference
  !           (for more detailed definition of BCs 5-7, see the
  !           comments of subroutine mkspline)
  !   NOTE bcxmin(...) referenced ONLY if ibcxmin=1 or ibcxmin=2
  !
  !   ibcxmax -- indicator for boundary condition at x(nx):
  !    bcxmax(...) -- boundary condition data
  !     (interpolation as with ibcxmin, bcxmin)
  !
  !  output:
  !   fspl(*,*,1..inx,1..inth) -- bicubic spline coeffs (4x4)
  !   ...fspl(1,1,*,*) is not replaced.
  !
  integer ilinx,ilinth,ier
  !
  !   ilinx -- =1 on output if x(inx) pts are nearly evenly spaced (tol=1e-3)
  !   ilinth-- =1 on output if th(inth) evenly spaced (tol=1e-3)
  !
  !   ier -- completion code, 0 for normal
  !
  !  workspace:
  !   wk -- must be at least 5*max(inx,inth) large
  !  nwk -- size of workspace
  !
  !---------------------------------
  !  compute bicubic spline of 2d function, given values at the grid
  !  grid crossing points, f(1,1,i,j)=f(x(i),th(j)).
  !
  !  on evaluation:  for point x btw x(i) and x(i+1), dx=x-x(i)
  !                       and th btw th(j) and th(j+1), dt=th-th(j),
  !
  !      spline =
  !        f(1,1,i,j) + dx*f(2,1,i,j) + dx**2*f(3,1,i,j) + dx**3*f(4,1,i,j)
  !   +dt*(f(1,2,i,j) + dx*f(2,2,i,j) + dx**2*f(3,2,i,j) + dx**3*f(4,2,i,j))
  !   +d2*(f(1,3,i,j) + dx*f(2,3,i,j) + dx**2*f(3,3,i,j) + dx**3*f(4,3,i,j))
  !   +d3*(f(1,4,i,j) + dx*f(2,4,i,j) + dx**2*f(3,4,i,j) + dx**3*f(4,4,i,j))
  !
  !      where d2=dt**2 and d3=dt**3.
  !
  !---------------------------------
  !
  integer ierx,ierth
  real(fp) :: zdum(1)
  !
  !---------------------------------
  !
  ier=0
  if(nwk.lt.5*max(inx,inth)) then
     write(6,'('' ?bpsplinb:  workspace too small.'')')
     ier=1
  end if
  if(inx.lt.2) then
     write(6,'('' ?bpsplinb:  at least 2 x points required.'')')
     ier=1
  end if
  if(inth.lt.2) then
     write(6,'('' ?bpsplinb:  need at least 2 theta points.'')')
     ier=1
  end if
  !
  call ibc_ck(ibcxmin,'bcspline','xmin',0,7,ier)
  call ibc_ck(ibcxmax,'bcspline','xmax',0,7,ier)
  !
  !  check ilinx & x vector
  !
  call splinck(x,inx,ilinx,1.0E-3_fp,ierx)
  if(ierx.ne.0) ier=2
  !
  if(ier.eq.2) then
     write(6,'('' ?bpsplinb:  x axis not strict ascending'')')
  end if
  !
  !  check ilinth & th vector
  !
  call splinck(th,inth,ilinth,1.0E-3_fp,ierth)
  if(ierth.ne.0) ier=3
  !
  if(ier.eq.3) then
     write(6,'('' ?bpsplinb:  th axis not strict ascending'')')
  end if
  !
  if(ier.ne.0) return
  !
  !------------------------------------
  zdum=0.0_fp
  !
  call bcspline(x,inx,th,inth,fspl,inf3, &
       ibcxmin,bcxmin,ibcxmax,bcxmax, &
       -1,zdum,-1,zdum, &
       wk,nwk,ilinx,ilinth,ier)
  !
end subroutine bpsplinb
