!  bpspline -- dmc 30 May 1996
!
!  set up coefficients for bicubic spline with following BC's:
!  * LHS and RHS BCs -- d3f/dx3 from 3rd divided difference, 1st coordinate
!  * derivatives periodic in second coordinate
!
!  see similar routine, bpsplinb, for more control over x boundary
!  conditions...
!
subroutine bpspline(x,inx,th,inth,fspl,inf3, &
     wk,nwk,ilinx,ilinth,ier)
  use precision_mod, only: fp
  !
  implicit none
  integer inx,inth,inf3,nwk
  real(fp) :: x(inx),th(inth),fspl(4,4,inf3,inth),wk(nwk)
  integer ilinx,ilinth
  integer ier
  !
  !  input:
  !    x(1...inx) -- abscissae, first dimension of data
  !   th(1...inth) -- abscissae, second (periodic) dimension of data
  !   fspl(1,1,1..inx,1..inth) -- function values
  !   inf3 -- fspl dimensioning, inf3.ge.inx required.
  !
  !  output:
  !   fspl(*,*,1..inx,1..inth) -- bicubic spline coeffs (4x4)
  !   ...fspl(1,1,*,*) is not replaced.
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
  integer ierx,ierth
  real(fp) :: zdum(1)
  !---------------------------------
  !
  ier=0
  if(nwk.lt.5*max(inx,inth)) then
     write(6,'('' ?bpspline:  workspace too small.'')')
     ier=1
  end if
  if(inx.lt.2) then
     write(6,'('' ?bpspline:  at least 2 x points required.'')')
     ier=1
  end if
  if(inth.lt.2) then
     write(6,'('' ?bpspline:  need at least 2 theta points.'')')
     ier=1
  end if
  !
  !  check ilinx & x vector
  !
  call splinck(x,inx,ilinx,1.0E-3_fp,ierx)
  if(ierx.ne.0) ier=2
  !
  if(ier.eq.2) then
     write(6,'('' ?bpspline:  x axis not strict ascending'')')
  end if
  !
  !  check ilinth & th vector
  !
  call splinck(th,inth,ilinth,1.0E-3_fp,ierth)
  if(ierth.ne.0) ier=3
  !
  if(ier.eq.3) then
     write(6,'('' ?bpspline:  th axis not strict ascending'')')
  end if
  !
  if(ier.ne.0) return
  !
  !------------------------------------
  !
  zdum=0.0_fp
  !
  !  not-a-knot BCs for x, periodic for theta
  !
  call bcspline(x,inx,th,inth,fspl,inf3, &
       7,zdum,7,zdum, &
       -1,zdum,-1,zdum, &
       wk,nwk,ilinx,ilinth,ier)
  !
  return
end subroutine bpspline
