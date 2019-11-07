!  tpsplinb -- dmc 20 Jan 1999
!
!  set up coefficients for bicubic spline with following BC's:
!  * 1st (x) dimension BCs under user control
!  * derivatives periodic in second & third coordinates
!
subroutine tpsplinb(x,inx,th,inth,ph,inph,fspl,inf4,inf5, &
     ibcxmin,bcxmin,ibcxmax,bcxmax,inb1, &
     wk,nwk,ilinx,ilinth,ilinph,ier)
  use psp_precision_mod, only: fp
  !
  implicit none
  integer inx,inth,inph,inf4,inf5,nwk,inb1
  real(fp) :: x(inx),th(inth),ph(inph)
  real(fp) :: fspl(4,4,4,inf4,inf5,inph),wk(nwk)
  integer ibcxmin,ibcxmax
  real(fp) :: bcxmin(inb1,inph),bcxmax(inb1,inph)
  integer ilinx,ilinth,ilinph,ier
  !
  !  input:
  !    x(1...inx) -- abscissae, first dimension of data
  !   th(1...inth) -- abscissae, second (periodic) dimension of data
  !   ph(1...inph) -- abscissae, third (periodic) dimension of data
  !   fspl(1,1,1,1..inx,1..inth,1..inph) -- function values
  !   inf4 -- fspl dimensioning, inf4.ge.inx required.
  !   inf5 -- fspl dimensioning, inf5.ge.inth required.
  !
  !  boundary conditions input:
  !   ibcxmin -- indicator for boundary condition at x(1):
  !    bcxmin(...) -- boundary condition data
  !     =0 -- use "not a knot", bcxmin(...) ignored
  !     =1 -- match slope, specified at x(1),th(ith),ph(iph) by bcxmin(ith,iph)
  !     =2 -- match 2nd derivative, specified at x(1),th(ith),ph(iph)
  !           by bcxmin(ith,iph
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
  !   fspl(*,*,*,1..inx,1..inth,1..inph) -- bicubic spline coeffs (4x4)
  !   ...fspl(1,1,1,*,*,*) is not replaced.
  !
  !   ilinx -- =1 on output if x(inx) pts are nearly evenly spaced (tol=1e-3)
  !   ilinth-- =1 on output if th(inth) evenly spaced (tol=1e-3)
  !   ilinph-- =1 on output if ph(inph) evenly spaced (tol=1e-3)
  !
  !   ier -- completion code, 0 for normal
  !
  !  workspace:
  !   wk -- must be at least 5*max(inx,inth,inph) large
  !  nwk -- size of workspace
  !
  !---------------------------------
  !  compute tricubic spline of 3d function, given values at the
  !  grid crossing points, f(1,1,1,i,j,k)=f(x(i),th(j),ph(k)).
  !
  !  on evaluation:  for point x btw x(i) and x(i+1), dx=x-x(i)
  !                       and th btw th(j) and th(j+1), dt=th-th(j),
  !                       and ph btw ph(k) and ph(k+1), dp=ph-ph(k),
  !
  !      spline =
  !        f(1,1,1,i,j,k)+dx*f(2,1,1,i,j,k)+dx2*f(3,1,1,i,j,k)+dx3*f(4,1,1,i,j,k)
  !  +dt*(f(1,2,1,i,j,k)+dx*f(2,2,1,i,j,k)+dx2*f(3,2,1,i,j,k)+dx3*f(4,2,1,i,j,k))
  ! +dt2*(f(1,3,1,i,j,k)+dx*f(2,3,1,i,j,k)+dx2*f(3,3,1,i,j,k)+dx3*f(4,3,1,i,j,k))
  ! +dt3*(f(1,4,1,i,j,k)+dx*f(2,4,1,i,j,k)+dx2*f(3,4,1,i,j,k)+dx3*f(4,4,1,i,j,k))
  !        +dp*(
  !        f(1,1,2,i,j,k)+dx*f(2,1,2,i,j,k)+dx2*f(3,1,2,i,j,k)+dx3*f(4,1,2,i,j,k)
  !  +dt*(f(1,2,2,i,j,k)+dx*f(2,2,2,i,j,k)+dx2*f(3,2,2,i,j,k)+dx3*f(4,2,2,i,j,k))
  ! +dt2*(f(1,3,2,i,j,k)+dx*f(2,3,2,i,j,k)+dx2*f(3,3,2,i,j,k)+dx3*f(4,3,2,i,j,k))
  ! +dt3*(f(1,4,2,i,j,k)+dx*f(2,4,2,i,j,k)+dx2*f(3,4,2,i,j,k)+dx3*f(4,4,2,i,j,k)))
  !        +dp2*(
  !        f(1,1,3,i,j,k)+dx*f(2,1,3,i,j,k)+dx2*f(3,1,3,i,j,k)+dx3*f(4,1,3,i,j,k)
  !  +dt*(f(1,2,3,i,j,k)+dx*f(2,2,3,i,j,k)+dx2*f(3,2,3,i,j,k)+dx3*f(4,2,3,i,j,k))
  ! +dt2*(f(1,3,3,i,j,k)+dx*f(2,3,3,i,j,k)+dx2*f(3,3,3,i,j,k)+dx3*f(4,3,3,i,j,k))
  ! +dt3*(f(1,4,3,i,j,k)+dx*f(2,4,3,i,j,k)+dx2*f(3,4,3,i,j,k)+dx3*f(4,4,3,i,j,k)))
  !        +dp3*(
  !        f(1,1,4,i,j,k)+dx*f(2,1,4,i,j,k)+dx2*f(3,1,4,i,j,k)+dx3*f(4,1,4,i,j,k)
  !  +dt*(f(1,2,4,i,j,k)+dx*f(2,2,4,i,j,k)+dx2*f(3,2,4,i,j,k)+dx3*f(4,2,4,i,j,k))
  ! +dt2*(f(1,3,4,i,j,k)+dx*f(2,3,4,i,j,k)+dx2*f(3,3,4,i,j,k)+dx3*f(4,3,4,i,j,k))
  ! +dt3*(f(1,4,4,i,j,k)+dx*f(2,4,4,i,j,k)+dx2*f(3,4,4,i,j,k)+dx3*f(4,4,4,i,j,k)))
  !
  !      where dx2=dx**2 and dx3=dx**3.
  !      where dt2=dt**2 and dt3=dt**3.
  !      where dp2=dp**2 and dp3=dp**3.
  !
  !---------------------------------
  integer ierx,ierth,ierph
  real(fp) :: zdumx(inx)
  !---------------------------------
  !
  ier=0
  if(nwk.lt.5*max(inx,inth,inph)) then
     write(6,'('' ?tpsplinb:  workspace too small.'')')
     ier=1
  end if
  if(inx.lt.2) then
     write(6,'('' ?tpsplinb:  at least 2 x points required.'')')
     ier=1
  end if
  if(inth.lt.2) then
     write(6,'('' ?tpsplinb:  need at least 2 theta points.'')')
     ier=1
  end if
  if(inph.lt.2) then
     write(6,'('' ?tpsplinb:  need at least 2 phi points.'')')
     ier=1
  end if
  !
  if((ibcxmin.eq.1).or.(ibcxmax.eq.1).or.(ibcxmin.eq.2).or. &
       (ibcxmax.eq.2)) then
     if(inb1.lt.inth) then
        ier=1
        write(6, &
             '('' ?tpsplinb:  1st dim of bcxmin/max arrays .lt. inth'')')
     end if
  end if
  !
  call ibc_ck(ibcxmin,'tpsplinb','xmin',0,7,ier)
  call ibc_ck(ibcxmax,'tpsplinb','xmax',0,7,ier)
  !
  !  check ilinx & x vector
  !
  call splinck(x,inx,ilinx,1.0E-3_fp,ierx)
  if(ierx.ne.0) ier=2
  !
  if(ier.eq.2) then
     write(6,'('' ?tpsplinb:  x axis not strict ascending'')')
  end if
  !
  !  check ilinth & th vector
  !
  call splinck(th,inth,ilinth,1.0E-3_fp,ierth)
  if(ierth.ne.0) ier=3
  !
  if(ier.eq.3) then
     write(6,'('' ?tpsplinb:  theta axis not strict ascending'')')
  end if
  !
  !  check ilinth & th vector
  !
  call splinck(ph,inph,ilinph,1.0E-3_fp,ierph)
  if(ierph.ne.0) ier=4
  !
  if(ier.eq.4) then
     write(6,'('' ?tpsplinb:  phi axis not strict ascending'')')
  end if
  !
  if(ier.ne.0) return
  !
  !------------------------------------
  !
  call tcspline(x,inx,th,inth,ph,inph,fspl,inf4,inf5, &
       ibcxmin,bcxmin,ibcxmax,bcxmax,inb1, &
       -1,zdumx,-1,zdumx,inx, &
       -1,zdumx,-1,zdumx,inx, &
       wk,nwk,ilinx,ilinth,ilinph,ier)
  !
  return
end subroutine tpsplinb
