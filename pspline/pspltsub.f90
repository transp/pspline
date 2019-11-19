!  pspline library test routine
!
subroutine pspltsub(filename,zctrl)
  use psp_precision_mod, only: fp
  !
  !  write output to file
  !
  implicit none
  character(len=*), intent(in) :: filename  ! write output here; ' ' for stdio
  real(fp), intent(in) :: zctrl              ! control (not yet used...)
  !
  !  if filename ends in ".tmp" delete when done.
  !
  common/bc_ctrl/ nbc
  common/pspltest_io/ m
  !
  integer :: ilen,m,nbc
  logical :: tmpflag
  !-----------------------------------
  !
  tmpflag=.FALSE.
  if(filename.eq.' ') then
     m=6
  else
     m=99
     ilen=len(trim(filename))
     if(filename(max(1,(ilen-3)):ilen).eq.'.tmp') tmpflag = .TRUE.
     !
     open(unit=m,file=filename,status='unknown')
  end if
  !
  !-----------------------------------
  !  set nbc=0 to use "not a knot" instead of explicit 1st deriv bc
  !  set nbc=1 to use explicit condition based on analytic expression
  !  (more accurate)
  !
  nbc=1
  !
  write(m,1000)
1000 format(/' *** 1d spline tests ***'/ &
       '     f(x)=2+sin(x)'/ &
       '     x in the interval [0,2pi]')
  call pspltest1(zctrl)
  write(m,1001)
1001 format(/' *** 2d spline tests ***'/ &
       '     f(x,theta)=exp(2x-1)*(2+sin(theta))'/ &
       '     x in [0,1],  theta in [0,2pi]'/)
  call pspltest2(zctrl)
  write(m,1002)
1002 format(/' *** 3d spline tests ***'/ &
       '     f(x,theta,phi)=exp(2x-1)*(2+sin(theta))*(3+sin(phi))'/ &
       '     x in [0,1],  theta in [0,2pi],  phi in [0,2pi]'/)
  call pspltest3(zctrl)
  !
  !-----------------------------------
  !  dmc Feb 2010:  add tests to verify extension of pspline to
  !  small grids: nx=2, nx=3
  !
  call smalltest(zctrl)
  !
  !-----------------------------------
  !  close file
  !
  if(tmpflag) then
     close(unit=m,status='delete')
  else
     close(unit=m)
  end if
  !
  !-----------------------------------
  !
  return
end subroutine pspltsub
!
!------------------------------------------------
!
subroutine pspltest1(zctrl)
  use psp_precision_mod, only: fp
  !
  !============
  implicit none
  integer m,is,inum,ix
  !============
  real(fp) :: pi2
  !============
  real(fp), intent(in) :: zctrl
  !
  common/pspltest_io/ m
  !
  real(fp) :: x(80)
  real(fp) :: fs(320),fsp(320),xpkg(320)
  real(fp) :: fh(160),fs2(160)
  !
  real(fp) :: xtest(1000),ftest(1000),zdum(1000)
  real(fp) :: testa1(1000),testa2(1000),testa3(1000),testa4(1000)
  real(fp) :: testa5(1000)
  !
  real(fp) :: zcos(80),z2sin(80)
  real(fp) :: zero
  !
  integer isize(4)
  !
  isize(1) = 10
  isize(2) = 20
  isize(3) = 40
  isize(4) = 80
  !
  pi2 = 6.28318530718_fp
  zero = 0.0_fp
  !
  !  test splining of function:  f(x)=2+sin(x)
  !
  !-------------------------------------------
  !
  call tset(1000,xtest,ftest,zdum,zero-0.1_fp,pi2+0.1_fp)
  !
  do is=1,4
     inum=isize(is)
     call tset(inum,x,z2sin,zcos,zero,pi2)
     do ix=1,inum
        zcos(ix)=cos(x(ix))         ! df/dx
     end do
     !
     write(m,*) ' '
     call dotest1(inum,x,z2sin,zcos,fs,fsp,fh,fs2,1000, &
          xtest,ftest,xpkg,testa1,testa2,testa3,testa4,testa5,zdum)
  end do
  !
  return
end subroutine pspltest1
!------------------------------------------------
!
subroutine dotest1(ns,x,f,fd,fspl,fspp,fherm,fs2,nt,xt,ft,xpkg, &
     testa1,testa2,testa3,testa4,testa5,wk)
  use psp_precision_mod, only: fp
  !
  !
  !============
  implicit none
  integer nt,ns,m,nbc,i,ierg,ilinx,ier,iwarn
  !============
  real(fp) :: zdum,fmin,fmax,sdif,hdif,pdif,s2dif,pldif,sdifr,hdifr
  real(fp) :: pdifr,s2difr,pldifr,hermv,difabs,splinv,zlinv
  !============
  common/pspltest_io/ m
  !
  !  interpolant
  !
  real(fp) :: x(ns)                        ! interpolant grid
  real(fp) :: f(ns)                        ! fcn values at grid pts
  real(fp) :: fd(ns)                       ! derivatives at grid pts
  real(fp) :: fspl(4,ns)                   ! spline coeff. array (calc. here)
  real(fp) :: fspp(4,ns)                   ! spline coeff. array (calc. here)
  real(fp) :: fherm(0:1,ns)                ! hermite array (calc. here)
  real(fp) :: fs2(0:1,ns)                  ! compact spline coeff. array calc here
  !
  real(fp) :: xpkg(ns,4)                   ! x "package"
  !
  real(fp) :: wk(nt)                       ! workspace
  !
  real(fp) :: xt(nt)                       ! test grid
  real(fp) :: ft(nt)                       ! fcn values at test grid pts
  !
  real(fp) :: testa1(nt),testa2(nt),testa3(nt),testa4(nt),testa5(nt)
  !
  integer ict(3)
  !
  common/bc_ctrl/ nbc
  !
  ict = 0
  ict(1) = 1
  !
  !-------------------
  !
  do i=1,ns
     fherm(0,i)=f(i)
     fherm(1,i)=fd(i)
     fspl(1,i)=f(i)
     fspp(1,i)=f(i)
     fs2(0,i)=f(i)
  end do
  !
  call genxpkg(ns,x,xpkg,1,1,1,4.0E-7_fp,1,ierg)
  if(ierg.ne.0) write(m,*) ' ??dotest1:  genxpkg:  ierg=',ierg
  !
  !  spline setup calls
  !
  call cspline(x,ns,fspl,nbc,fd(1),nbc,fd(ns),wk,nt,ilinx,ier) ! 1st deriv bc
  call cspline(x,ns,fspp,-1,zdum,-1,zdum,wk,nt,ilinx,ier) ! periodic
  !
  call mkspline(x,ns,fs2,nbc,fd(1),nbc,fd(ns),ilinx,ier)
  !
  call akherm1p(x,ns,fherm,ilinx,1,ier)
  !
  !  vectorize spline/hermite evaluation calls
  !
  call spgrid(xt,nt,testa1,ns,xpkg,fspl,iwarn,ier)
  if(iwarn.ne.0) write(m,*) ' ?dotest1:  spgrid(1):  iwarn=',iwarn
  if(ier.ne.0) write(m,*) ' ?dotest1:  spgrid(1):  ier=',ier
  !
  call spgrid(xt,nt,testa2,ns,xpkg,fspp,iwarn,ier)
  if(iwarn.ne.0) write(m,*) ' ?dotest1:  spgrid(2):  iwarn=',iwarn
  if(ier.ne.0) write(m,*) ' ?dotest1:  spgrid(2):  ier=',ier
  !
  call gridspline(xt,nt,testa3,ns,xpkg,fs2,iwarn,ier)
  if(iwarn.ne.0) write(m,*) ' ?dotest1:  gridspline:  iwarn=', &
       iwarn
  if(ier.ne.0) write(m,*) ' ?dotest1:  gridspline:  ier=',ier
  !
  call gridherm1(xt,nt,testa4,ns,xpkg,fherm,iwarn,ier)
  if(iwarn.ne.0) write(m,*) ' ?dotest1:  gridherm1:  iwarn=', &
       iwarn
  if(ier.ne.0) write(m,*) ' ?dotest1:  gridherm1:  ier=',ier
  !
  call gridpc1(xt,nt,testa5,ns,xpkg,f,iwarn,ier)
  if(iwarn.ne.0) write(m,*) ' ?dotest1:  gridpc1:  iwarn=', &
       iwarn
  if(ier.ne.0) write(m,*) ' ?dotest1:  gridpc1:  ier=',ier
  !
  fmin=1.0E30_fp
  fmax=-fmin
  sdif=0.0_fp
  hdif=0.0_fp
  pdif=0.0_fp
  s2dif=0.0_fp
  pldif=0.0_fp
  sdifr=0.0_fp
  hdifr=0.0_fp
  pdifr=0.0_fp
  s2difr=0.0_fp
  pldifr=0.0_fp
  !
  ier=0
  !
  do i=1,nt
     fmin=min(ft(i),fmin)
     fmax=max(ft(i),fmax)
     !
     hermv=testa4(i)
     !xx         call herm1ev(xt(i),x,ns,ilinx,fherm,ict,hermv,ier)
     !xx         if(ier.ne.0) write(m,'('' ?? ier.ne.0 in dotest1 (herm1ev)'')')
     !
     difabs=abs(hermv-ft(i))
     hdif=max(hdif,difabs)
     hdifr=max(hdifr,difabs/ft(i))
     !
     splinv=testa1(i)
     !xx         call cspeval(xt(i),ict,splinv,x,ns,ilinx,fspl,ier)
     !xx         if(ier.ne.0) write(m,'('' ?? ier.ne.0 in dotest1 (cspeval)'')')
     !
     difabs=abs(splinv-ft(i))
     sdif=max(sdif,difabs)
     sdifr=max(sdifr,difabs/ft(i))
     !
     splinv=testa2(i)
     !xx         call cspeval(xt(i),ict,splinv,x,ns,ilinx,fspp,ier)
     !xx         if(ier.ne.0) write(m,'('' ?? ier.ne.0 in dotest1 (cspeval)'')')
     !
     difabs=abs(splinv-ft(i))
     pdif=max(pdif,difabs)
     pdifr=max(pdifr,difabs/ft(i))
     !
     splinv=testa3(i)
     !xx         call evspline(xt(i),x,ns,ilinx,fs2,ict,splinv,ier)
     !xx         if(ier.ne.0) write(m,'('' ?? ier.ne.0 in dotest1 (cspeval)'')')
     !
     difabs=abs(splinv-ft(i))
     s2dif=max(s2dif,difabs)
     s2difr=max(s2difr,difabs/ft(i))
     !
     zlinv=testa5(i)
     difabs=abs(zlinv-ft(i))
     pldif=max(pldif,difabs)
     pldifr=max(pldifr,difabs/ft(i))
     !
  end do
  !
  write(m,1000) '1d periodic spline',ns,fmin,fmax,pdif,pdifr
  write(m,1000) '1d spline w/bdy cond.',ns,fmin,fmax,sdif,sdifr
  write(m,1000) '1d compact spline w/bdy cond.',ns,fmin,fmax, &
       s2dif,s2difr
  write(m,1000) '1d Hermite interpolation',ns,fmin,fmax,hdif,hdifr
  write(m,1000) '1d piecewise linear',ns,fmin,fmax,pldif,pldifr
1000 format(/2x,a,' setup on ',i3,' point grid'/ &
       '   fmin=',1pe11.4,' fmax=',1pe11.4,' max(|diff|)=',1pe11.4/ &
       '     max(|diff|/f) = ',1pe11.4)
  !
  return
end subroutine dotest1
!
!------------------------------------------------
!
subroutine pspltest2(zctrl)
  use psp_precision_mod, only: fp
  !
  !============
  implicit none
  integer m
  !============
  real(fp) :: pi2
  !============
  real(fp), intent(in) :: zctrl
  !
  common/pspltest_io/ m
  !
  !  test various bicubic spline routines
  !  for fitting the function
  !
  !     f(x,th)=exp(2*x-1)*(2+sin(th))
  !
  !  on 3 grid sizes; compare accuracy on test grid.
  !
  !  grid ranges:  x in [0,1], th in [0,2pi]
  !
  real(fp) :: x1(10),x2(20),x4(40),ex1(10),ex2(20),ex4(40)
  real(fp) :: t1(10),t2(20),t4(40),st1(10),st2(20),st4(40), &
       ct1(10),ct2(20),ct4(40)
  real(fp) :: bcx1(40),bcx2(40),bcth1(40),bcth2(40)
  !
  real(fp) :: f1(4,4,10,10),f2(4,4,20,20),f4(4,4,40,40),fh(6400)
  real(fp) :: flin(1600)
  !
  real(fp) :: xtest(200),extest(200),ttest(200),stest(200),ctest(200)
  real(fp) :: zero,one
  !
  pi2 = 6.28318530718_fp
  zero = 0.0_fp
  one = 1.0_fp
  !
  !---------------------------
  !
  call xset(10,x1,ex1,zero,one)
  call xset(20,x2,ex2,zero,one)
  call xset(40,x4,ex4,zero,one)
  call xset(200,xtest,extest,zero,one)
  !
  call tset(10,t1,st1,ct1,zero,pi2)
  call tset(20,t2,st2,ct2,zero,pi2)
  call tset(40,t4,st4,ct4,zero,pi2)
  call tset(200,ttest,stest,ctest,zero,pi2)
  !
  call ffset(10,ex1,st1,f1)
  call ffset(20,ex2,st2,f2)
  call ffset(40,ex4,st4,f4)
  !
  call dotest2(x1,ex1,10,t1,st1,ct1,10,f1,fh,flin, &
       bcx1,bcx2,bcth1,bcth2, &
       xtest,extest,ttest,stest,200)
  !
  call dotest2(x2,ex2,20,t2,st2,ct2,20,f2,fh,flin, &
       bcx1,bcx2,bcth1,bcth2, &
       xtest,extest,ttest,stest,200)
  !
  call dotest2(x4,ex4,40,t4,st4,ct4,40,f4,fh,flin, &
       bcx1,bcx2,bcth1,bcth2, &
       xtest,extest,ttest,stest,200)
  !
  return
end subroutine pspltest2
!-----------------------------------------------------------
subroutine xset(nx,x,ex,xmin,xmax)
  use psp_precision_mod, only: fp
  !
  !============
  implicit none
  integer nx,ix
  !============
  real(fp) :: xmin,xmax
  !============
  real(fp) :: x(nx)
  real(fp) :: ex(nx)
  !
  do ix=1,nx
     x(ix)=xmin + real(ix-1,fp)*(xmax-xmin)/real(nx-1,fp)
     ex(ix)=exp(2.0_fp*x(ix)-1.0_fp)
  end do
  !
  return
end subroutine xset
!-----------------------------------------------------------
subroutine tset(nth,th,sth,cth,thmin,thmax)
  use psp_precision_mod, only: fp
  !
  !============
  implicit none
  integer nth,ith
  !============
  real(fp) :: thmin,thmax
  !============
  real(fp) :: th(nth)
  real(fp) :: sth(nth)
  real(fp) :: cth(nth)
  !
  do ith=1,nth
     th(ith)=thmin + real(ith-1,fp)*(thmax-thmin)/real(nth-1,fp)
     sth(ith)=2.0_fp+sin(th(ith))
     cth(ith)=cos(th(ith))
  end do
  !
  return
end subroutine tset
!-----------------------------------------------------------
subroutine tset3(nph,ph,sph,cph,phmin,phmax)
  use psp_precision_mod, only: fp
  !
  !============
  implicit none
  integer nph,iph
  !============
  real(fp) :: phmin,phmax
  !============
  real(fp) :: ph(nph)
  real(fp) :: sph(nph)
  real(fp) :: cph(nph)
  !
  do iph=1,nph
     ph(iph)=phmin + real(iph-1,fp)*(phmax-phmin)/real(nph-1,fp)
     sph(iph)=3.0_fp+sin(ph(iph))
     cph(iph)=cos(ph(iph))
  end do
  !
  return
end subroutine tset3
!-----------------------------------------------------------
subroutine ffset(num,xf,tf,f)
  use psp_precision_mod, only: fp
  !
  !============
  implicit none
  integer num,j,i
  !============
  real(fp) :: xf(num),tf(num),f(4,4,num,num)
  !
  do j=1,num
     do i=1,num
        f(1,1,i,j)=xf(i)*tf(j)
     end do
  end do
  !
  return
end subroutine ffset
!-----------------------------------------------------------
subroutine fset3(num,xf,tf,pf,f)
  use psp_precision_mod, only: fp
  !
  !============
  implicit none
  integer num,k,j,i
  !============
  real(fp) :: xf(num),tf(num),pf(num),f(4,4,4,num,num,num)
  !
  do k=1,num
     do j=1,num
        do i=1,num
           f(1,1,1,i,j,k)=xf(i)*tf(j)*pf(k)
        end do
     end do
  end do
  !
  return
end subroutine fset3
!----------------------------------------------------------
subroutine bset(fx,nx,fth,nth,bcx1,bcx2,bcth1,bcth2)
  use psp_precision_mod, only: fp
  !
  !============
  implicit none
  integer nth,nx,ith,ix
  !============
  real(fp) :: fx(nx)                       ! x factor of test fcn
  real(fp) :: fth(nth)                     ! th factor of test fcn
  real(fp) :: bcx1(nth),bcx2(nth)          ! df/dx bdyy conds at x(1),x(nx)
  real(fp) :: bcth1(nx),bcth2(nx)          ! df/dth bdy conds at th(1),th(nth)
  !
  !  df/dx = 2*exp(2x-1)*(2+sin(th)) = 2*f
  !
  do ith=1,nth
     bcx1(ith)=2.0_fp*fx(1)*fth(ith)   ! df/dx @ x(1)
     bcx2(ith)=2.0_fp*fx(nx)*fth(ith)  ! df/dx @ x(nx)
  end do
  !
  !  df/dth = exp(2x-1)*cos(th)
  !  cos(0)=cos(2pi)=1
  !
  do ix=1,nx
     bcth1(ix)=fx(ix)               ! df/dth @ th=0 = th(1)
     bcth2(ix)=fx(ix)               ! df/dth @ th=2pi = th(nth)
  end do
  !
  return
end subroutine bset
!---------------------------------------------------------
subroutine dotest2(x,fx,nx,th,fth,dfth,nth,f,fh,flin, &
     bcx1,bcx2,bcth1,bcth2, &
     xtest,fxtest,thtest,fthtest,ntest)
  use psp_precision_mod, only: fp
  !
  !============
  implicit none
  integer nth,ntest,nx,m,nbc,ith,ix,ilinx,ilinth,ier
  !============
  common/pspltest_io/ m
  !
  !  test spline of f; f(i,j)=fx(i)*fth(j) on spline grid
  !                    f(i,j)=fxtest(i)*fthtest(j) on test grid
  !
  !  f is exp(2x-1)*(2+sin(th))  df/dx = 2*f
  !
  real(fp) :: x(nx),fx(nx)                 ! x & fx vectors (already set)
  real(fp) :: th(nth),fth(nth),dfth(nth)   ! th & fth already set
  !
  real(fp) :: f(4,4,nx,nth)                ! spline, f(1,1,*,*) already set
  real(fp) :: fh(0:3,nx,nth)               ! hermite array
  real(fp) :: flin(nx,nth)                 ! piecewise linear array
  !
  real(fp) :: bcx1(nth),bcx2(nth)          ! bcs: dfdx vs. th @ x(1), x(nx).
  real(fp) :: bcth1(nx),bcth2(nx)          ! bcs: dfdth vs. x @ th(1), th(nth).
  !
  real(fp) :: xtest(ntest),fxtest(ntest)   ! x test grid & fx
  real(fp) :: thtest(ntest),fthtest(ntest) ! th test grid & fth
  !
  real(fp) :: wk(64000)                    ! workspace
  !
  common/bc_ctrl/ nbc
  !
  !--------------
  !
  if(ntest*ntest.gt.64000) then
     write(m,*) ' ?dotest2:  ntest*ntest exceeds workspace dim.'
     return
  end if
  !
  !  set up hermite array
  !
  do ith=1,nth
     do ix=1,nx
        flin(ix,ith)=f(1,1,ix,ith)         ! f
        fh(0,ix,ith)=f(1,1,ix,ith)         ! f
        fh(1,ix,ith)=2.0_fp*f(1,1,ix,ith)     ! df/dx
        fh(2,ix,ith)=fx(ix)*dfth(ith)      ! df/dy
        fh(3,ix,ith)=2.0_fp*fx(ix)*dfth(ith)  ! d2f/dxdy
     end do
  end do
  !
  call akherm2p(x,nx,th,nth,fh,nx,ilinx,ilinth,2,1,ier)
  !
  call compare('hermite',x,nx,th,nth,f,fh,flin,ilinx,ilinth, &
       xtest,fxtest,thtest,fthtest,ntest,wk)
  !
  !  set bdy conds
  !
  call bset(fx,nx,fth,nth,bcx1,bcx2,bcth1,bcth2)
  !
  call bpspline(x,nx,th,nth,f,nx,wk,15000,ilinx,ilinth,ier)
  if(ier.ne.0) then
     write(m,*) ' ?? error in pspltest:  dotest2(bpspline)'
  end if
  !
  call compare('bpspline',x,nx,th,nth,f,fh,flin,ilinx,ilinth, &
       xtest,fxtest,thtest,fthtest,ntest,wk)
  !
  !
  call bpsplinb(x,nx,th,nth,f,nx, &
       nbc,bcx1,nbc,bcx2, &
       wk,15000,ilinx,ilinth,ier)
  if(ier.ne.0) then
     write(m,*) ' ?? error in pspltest:  dotest2(bpsplinb)'
  end if
  !
  call compare('bpsplinb',x,nx,th,nth,f,fh,flin,ilinx,ilinth, &
       xtest,fxtest,thtest,fthtest,ntest,wk)
  !
  !
  call bcspline(x,nx,th,nth,f,nx, &
       nbc,bcx1,nbc,bcx2, &
       nbc,bcth1,nbc,bcth2, &
       wk,15000,ilinx,ilinth,ier)
  if(ier.ne.0) then
     write(m,*) ' ?? error in pspltest:  dotest2(bcspline)'
     return
  end if
  !
  call compare('bcspline',x,nx,th,nth,f,fh,flin,ilinx,ilinth, &
       xtest,fxtest,thtest,fthtest,ntest,wk)
  !
  !  compact spline representation (as per L. Zakharov)
  !
  do ith=1,nth
     do ix=1,nx
        fh(0,ix,ith)=f(1,1,ix,ith)         ! f
     end do
  end do
  !
  call mkbicub(x,nx,th,nth,fh,nx, &
       nbc,bcx1,nbc,bcx2, &
       nbc,bcth1,nbc,bcth2, &
       ilinx,ilinth,ier)
  !
  call compare('mkbicub',x,nx,th,nth,f,fh,flin,ilinx,ilinth, &
       xtest,fxtest,thtest,fthtest,ntest,wk)
  !
  !  piecewise linear...
  !
  call compare('piecewise linear',x,nx,th,nth,f, &
       fh,flin,ilinx,ilinth, &
       xtest,fxtest,thtest,fthtest,ntest,wk)
  !
  return
end subroutine dotest2
!--------------------------------
!
subroutine compare(slbl,x,nx,th,nth,f,fh,fl,ilinx,ilinth, &
     xtest,fxtest,thtest,fthtest,ntest,wk)
  use psp_precision_mod, only: fp
  !
  !============
  implicit none
  integer nth,ntest,nx,m,icycle,iherm,ier,ijk,iwarn,j,i
  !============
  real(fp) :: fmin,fmax,fdif,fdifr,zdum,zth,zx,ff,fs
#ifdef _TIMER
  real(fp) :: ztime1,ztime2,zctime
#endif
  !============
  common/pspltest_io/ m
  !
  character(len=*) slbl                ! spline coeff routine:  label
  real(fp) :: x(nx),th(nth)                ! indep. coords.
  real(fp) :: f(4,4,nx,nth)                ! spline data
  real(fp) :: fh(0:3,nx,nth)               ! hermite data
  real(fp) :: fl(nx,nth)                   ! piecewise linear data
  integer ilinx,ilinth              ! even spacing flags
  !
  !  test data grid & data:
  !
  real(fp) :: xtest(ntest),fxtest(ntest),thtest(ntest),fthtest(ntest)
  real(fp) :: wk(ntest,ntest)
  !
  !-------------
  !  select spline fcn eval only (no derivatives)
  !
  integer isel(10)
  !
  real(fp) :: fget(10)
  !
  real(fp) :: xpkg(200,1),thpkg(200,1)
  !-------------
  !
  isel = 0
  isel(1) = 1
  !
  icycle=20
  !
  iherm=0
  if(slbl.eq.'hermite') then
     write(m,*) ' '
     iherm=1
  end if
  if(slbl.eq.'mkbicub') then
     iherm=2
  end if
  if(slbl.eq.'piecewise linear') then
     iherm=3
  end if
  !
  fmin=1.0E30_fp
  fmax=-1.0E30_fp
  fdif=0.0_fp
  fdifr=0.0_fp
  !
  call genxpkg(nx,x,xpkg,0,1,0,zdum,1,ier)
  call genxpkg(nth,th,thpkg,1,1,0,zdum,1,ier)
  !
  if(iherm.eq.0) then
#ifdef _TIMER
     call cpu_time(ztime1)
#endif
     do ijk=1,icycle
        call bcspgrid(xtest,ntest,thtest,ntest,wk,ntest, &
             nx,xpkg,nth,thpkg,f,nx,iwarn,ier)
     end do
#ifdef _TIMER
     call cpu_time(ztime2)
#endif
  else if(iherm.eq.1) then
#ifdef _TIMER
     call cpu_time(ztime1)
#endif
     do ijk=1,icycle
        call gridherm2(xtest,ntest,thtest,ntest,wk,ntest, &
             nx,xpkg,nth,thpkg,fh,nx,iwarn,ier)
     end do
#ifdef _TIMER
     call cpu_time(ztime2)
#endif
  else if(iherm.eq.2) then
#ifdef _TIMER
     call cpu_time(ztime1)
#endif
     do ijk=1,icycle
        call gridbicub(xtest,ntest,thtest,ntest,wk,ntest, &
             nx,xpkg,nth,thpkg,fh,nx,iwarn,ier)
     end do
#ifdef _TIMER
     call cpu_time(ztime2)
#endif
  else
#ifdef _TIMER
     call cpu_time(ztime1)
#endif
     do ijk=1,icycle
        call gridpc2(xtest,ntest,thtest,ntest,wk,ntest, &
             nx,xpkg,nth,thpkg,fl,nx,iwarn,ier)
     end do
#ifdef _TIMER
     call cpu_time(ztime2)
#endif
  end if
  !
  do j=1,ntest
     zth=thtest(j)
     do i=1,ntest
        zx=xtest(i)
        !
        ff=fxtest(i)*fthtest(j)
        fmin=min(fmin,ff)
        fmax=max(fmax,ff)
        !
        if(iherm.eq.0) then
           fget(1)=wk(i,j)
           !xx               call bcspeval(zx,zth, &
           !xx                 isel,fget,x,nx,th,nth,ilinx,ilinth, &
           !xx                 f,nx,ier)
        else if(iherm.eq.1) then
           fget(1)=wk(i,j)
           !xx               call herm2ev(zx,zth,x,nx,th,nth,ilinx,ilinth, &
           !xx                  f,nx,isel,fget,ier)
        else if(iherm.eq.2) then
           fget(1)=wk(i,j)
           !xx               call evbicub(zx,zth,x,nx,th,nth,ilinx,ilinth, &
           !xx                  f,nx,isel,fget,ier)
        else
           fget(1)=wk(i,j)
        end if
        !xx            if(ier.ne.0) then
        !xx               write(m,*) ' ??compare ('//slbl//') ier.ne.0 exit'
        !xx               return
        !xx            end if
        !
        fs=fget(1)
        !
        fdif=max(fdif,abs(ff-fs))
        fdifr=max(fdifr,abs(ff-fs)/(0.5_fp*(ff+fs)))
        wk(i,j)=ff-fs
        !
     end do
  end do
  !
#ifdef _TIMER
  zctime=ztime2-ztime1
#endif
  write(m,1000) slbl,nx,nth,fmin,fmax,fdif,fdifr
1000 format(2x,a,' setup on ',i3,' x ',i3,' grid'/ &
       '   fmin=',1pe11.4,' fmax=',1pe11.4,' max(|diff|)=',1pe11.4/ &
       '     max(|diff|/f) = ',1pe11.4)
#ifdef _TIMER
  write(m,1001) icycle,ntest,ntest,zctime
#else
  write(m,1002) icycle,ntest,ntest
#endif
1001 format(2x,i3,' x ',i3,' x ',i3,' evaluations, cpu = ', &
       1pe11.4,' (s)')
1002 format(2x,i3,' x ',i3,' x ',i3,' evaluations')
  !
  return
end subroutine compare
!------------------------
subroutine pspltest3(zctrl)
  use psp_precision_mod, only: fp
  !
  !============
  implicit none
  integer m
  !============
  real(fp) :: pi2
  !============
  real(fp), intent(in) :: zctrl
  !
  common/pspltest_io/ m
  !
  !  test various bicubic spline routines
  !  for fitting the function
  !
  !     f(x,th,ph)=exp(2*x-1)*(2+sin(th))*(3+sin(ph))
  !
  !  on 3 grid sizes; compare accuracy on test grid.
  !
  !  grid ranges:  x in [0,1], th in [0,2pi], ph in [0,2pi].
  !
  real(fp) :: x1(10),x2(20),x4(40),ex1(10),ex2(20),ex4(40)
  real(fp) :: t1(10),t2(20),t4(40),st1(10),st2(20),st4(40), &
       ct1(10),ct2(20),ct4(40)
  real(fp) :: p1(10),p2(20),p4(40),sp1(10),sp2(20),sp4(40), &
       cp1(10),cp2(20),cp4(40)
  real(fp) :: bcx1(1600),bcx2(1600)
  real(fp) :: bcth1(1600),bcth2(1600)
  real(fp) :: bcph1(1600),bcph2(1600)
  !
  real(fp) :: f1(4,4,4,10,10,10),f2(4,4,4,20,20,20),f4(4,4,4,40,40,40)
  real(fp) :: fh(8,40,40,40),flin(40,40,40)
  !
  real(fp) :: xtest(100),extest(100),ttest(100),stest(100),zdum(100)
  real(fp) :: phtest(100),sptest(100)
  real(fp) :: zero,one
  !
  pi2 = 6.28318530718_fp
  zero = 0.0_fp
  one = 1.0_fp
  !
  !---------------------------
  !
  call xset(10,x1,ex1,zero,one)
  call xset(20,x2,ex2,zero,one)
  call xset(40,x4,ex4,zero,one)
  call xset(100,xtest,extest,zero,one)
  !
  call tset(10,t1,st1,ct1,zero,pi2)
  call tset(20,t2,st2,ct2,zero,pi2)
  call tset(40,t4,st4,ct4,zero,pi2)
  call tset(100,ttest,stest,zdum,zero,pi2)
  !
  call tset3(10,p1,sp1,cp1,zero,pi2)
  call tset3(20,p2,sp2,cp2,zero,pi2)
  call tset3(40,p4,sp4,cp4,zero,pi2)
  call tset3(100,phtest,sptest,zdum,zero,pi2)
  !
  call fset3(10,ex1,st1,sp1,f1)
  call fset3(20,ex2,st2,sp2,f2)
  call fset3(40,ex4,st4,sp4,f4)
  !
  call dotest3(x1,ex1,10,t1,st1,ct1,10,p1,sp1,cp1,10,f1,fh,flin, &
       bcx1,bcx2,bcth1,bcth2,bcph1,bcph2, &
       xtest,extest,ttest,stest,phtest,sptest,100)
  !
  call dotest3(x2,ex2,20,t2,st2,ct2,20,p2,sp2,cp2,20,f2,fh,flin, &
       bcx1,bcx2,bcth1,bcth2,bcph1,bcph2, &
       xtest,extest,ttest,stest,phtest,sptest,100)
  !
  call dotest3(x4,ex4,40,t4,st4,ct4,40,p4,sp4,cp4,40,f4,fh,flin, &
       bcx1,bcx2,bcth1,bcth2,bcph1,bcph2, &
       xtest,extest,ttest,stest,phtest,sptest,100)
  !
  return
end subroutine pspltest3
!---------------------------------------------------------
subroutine dotest3(x,fx,nx,th,fth,dfth,nth,ph,fph,dfph,nph, &
     f,fh,flin, &
     bcx1,bcx2,bcth1,bcth2,bcph1,bcph2, &
     xtest,fxtest,thtest,fthtest,phtest,fphtest,ntest)
  use psp_precision_mod, only: fp
  !
  !
  !============
  implicit none
  integer nth,nph,ntest,nx,m,nbc,inwk,itot,iph,ith,ix,ilinx
  integer ilinth,ilinph,ier
  !============
  real(fp) :: ztime1,ztime2,zdiff,zdiffr,zbc
  !============
  common/pspltest_io/ m
  !
  !  test spline of f; f(i,j,k)=fx(i)*fth(j)*fph(k) on spline grid
  !                    f(i,j,k)=fxtest(i)*fthtest(j)*fphtest(k) on test grid
  !
  !  f is exp(2x-1)*(2+sin(th))*(3+sin(ph)) -- derivatives for BCs evaluated
  !  here...
  !
  !  df/dx = 2*f
  !
  real(fp) :: x(nx),fx(nx)                 ! x & fx vectors (already set)
  real(fp) :: th(nth),fth(nth),dfth(nth)   ! th & fth & fth' already set
  real(fp) :: ph(nph),fph(nph),dfph(nth)   ! ph & fph & fph' already set
  !
  real(fp) :: f(4,4,4,nx,nth,nph)          ! spline, f(1,1,1,*,*,*) already set
  real(fp) :: fh(0:7,nx,nth,nph)           ! hermite array
  real(fp) :: flin(nx,nth,nph)             ! function data only -- array
  !
  real(fp) :: bcx1(nth,nph),bcx2(nth,nph)  ! bcs: dfdx vs. th,ph @ x(1), x(nx).
  real(fp) :: bcth1(nx,nph),bcth2(nx,nph)  ! bcs: dfdth vs. x,ph @ th(1), th(nth).
  real(fp) :: bcph1(nx,nth),bcph2(nx,nth)  ! bcs: dfdph vs. th,ph @ ph(1), ph(nph)
  !
  real(fp) :: xtest(ntest),fxtest(ntest)   ! x test grid & fx
  real(fp) :: thtest(ntest),fthtest(ntest) ! th test grid & fth
  real(fp) :: phtest(ntest),fphtest(ntest) ! ph test grid & fph
  !
  real(fp) :: wk(80*40*40*40)
  !
  common/bc_ctrl/ nbc
  !
  !--------------
  real(fp) :: zsave(20,20,20)
  real(fp) :: zvals(10)
  integer iselect(10)
  iselect = 0
  iselect(1:4) = 1
  !--------------
  !
  inwk=80*40*40*40
  !
  itot=ntest*ntest*ntest
  write(m,999) itot
999 format(/ &
       ' %dotest3:  4 x ',i7, &
       ' evaluations in progress -- be patient.'/)
  !
  do iph=1,nph
     do ith=1,nth
        do ix=1,nx
           flin(ix,ith,iph)=f(1,1,1,ix,ith,iph)
           fh(0,ix,ith,iph)=f(1,1,1,ix,ith,iph)
           fh(1,ix,ith,iph)=2.0_fp*f(1,1,1,ix,ith,iph)
           fh(2,ix,ith,iph)=fx(ix)*dfth(ith)*fph(iph)
           fh(3,ix,ith,iph)=fx(ix)*fth(ith)*dfph(iph)
           fh(4,ix,ith,iph)=2.0_fp*fx(ix)*dfth(ith)*fph(iph)
           fh(5,ix,ith,iph)=2.0_fp*fx(ix)*fth(ith)*dfph(iph)
           fh(6,ix,ith,iph)=fx(ix)*dfth(ith)*dfph(iph)
           fh(7,ix,ith,iph)=2.0_fp*fx(ix)*dfth(ith)*dfph(iph)
        end do
     end do
  end do
  !
  call akherm3p(x,nx,th,nth,ph,nph,fh,nx,nth, &
       ilinx,ilinth,ilinph,2,1,1,ier)
  !
  if(ier.ne.0) then
     write(m,*) ' ?? error in pspltest:  dotest3(akherm3p)'
  end if
  !
  call compare3('hermite',x,nx,th,nth,ph,nph,f,fh,flin, &
       ilinx,ilinth,ilinph, &
       xtest,fxtest,thtest,fthtest,phtest,fphtest,ntest)
  !
  !  set bdy conds
  !
  call bset3(fx,nx,fth,nth,fph,nph, &
       bcx1,bcx2,bcth1,bcth2,bcph1,bcph2)
  !
  call tpspline(x,nx,th,nth,ph,nph,f,nx,nth, &
       wk,inwk,ilinx,ilinth,ilinph,ier)
  if(ier.ne.0) then
     write(m,*) ' ?? error in pspltest:  dotest3(tpspline)'
  end if
  !
  call compare3('tpspline',x,nx,th,nth,ph,nph,f,fh,flin, &
       ilinx,ilinth,ilinph, &
       xtest,fxtest,thtest,fthtest,phtest,fphtest,ntest)
  !
  !
  call tpsplinb(x,nx,th,nth,ph,nph,f,nx,nth, &
       nbc,bcx1,nbc,bcx2,nth, &
       wk,inwk,ilinx,ilinth,ilinph,ier)
  if(ier.ne.0) then
     write(m,*) ' ?? error in pspltest:  dotest3(tpsplinb)'
  end if
  !
  call compare3('tpsplinb',x,nx,th,nth,ph,nph,f,fh,flin, &
       ilinx,ilinth,ilinph, &
       xtest,fxtest,thtest,fthtest,phtest,fphtest,ntest)
  !
  if(max(nx,nth,nph).le.20) then
     do iph=1,nph
        do ith=1,nth
           do ix=1,nx
              zsave(ix,ith,iph)=f(1,1,1,ix,ith,iph)
           end do
        end do
     end do
  end if
  !
#ifdef _TIMER
  call cpu_time(ztime1)
#endif
  call tcspline(x,nx,th,nth,ph,nph,f,nx,nth, &
       nbc,bcx1,nbc,bcx2,nth, &
       nbc,bcth1,nbc,bcth2,nx, &
       nbc,bcph1,nbc,bcph2,nx, &
       wk,inwk,ilinx,ilinth,ilinph,ier)
  if(ier.ne.0) then
     write(m,*) ' ?? error in pspltest:  dotest3(tcspline)'
     return
  end if
#ifdef _TIMER
  call cpu_time(ztime2)
  write(m,7706) nx,ztime2-ztime1
7706 format(' ==> tcspline setup (nx=',i3,') cpu time (s):',1pe11.4, ' <== ')
#else
  write(m,7706) nx
7706 format(' ==> tcspline setup (nx=',i3,') <== ')
#endif
  !
  if(max(nx,nth,nph).le.20) then
     do iph=1,nph
        do ith=1,nth
           do ix=1,nx
              call tcspeval(x(ix),th(ith),ph(iph), &
                   iselect, zvals, &
                   x,nx,th,nth,ph,nph, 1, 1, 1, f,nx,nth, ier)
              zdiff=abs(zsave(ix,ith,iph)-zvals(1))
              zdiffr=zdiff/zvals(1)
              if(zdiffr.gt.(1.0E-6_fp)) then
                 write(m,7701) ix,ith,iph,zsave(ix,ith,iph),zvals(1)
              end if
           end do
        end do
     end do
     !
7701 format(' ix=',i2,' ith=',i2,' iph=',i2,' ** f changed:', &
          2(1x,1pe13.6))
     !
     do iph=1,nph
        do ith=1,nth
           do ix=1,nx,nx-1
              if(ix.eq.1) then
                 zbc=bcx1(ith,iph)
              else
                 zbc=bcx2(ith,iph)
              end if
              call tcspeval(x(ix),th(ith),ph(iph), iselect, zvals, &
                   x,nx,th,nth,ph,nph, 1, 1, 1, f,nx,nth, ier)
              if(nbc.eq.1) then
                 if(abs(zbc-zvals(2)).gt.1.0E-4_fp) &
                      write(m,7702) ix,ith,iph,zbc,zvals(2)
              end if
           end do
        end do
     end do
     !
7702 format(' df/dx BC check @ix=',i2,',ith=',i2,',iph=',i2,': ', &
          2(1x,1pe13.6))
     !
     do iph=1,nph
        do ix=1,nx
           do ith=1,nth,nth-1
              if(ith.eq.1) then
                 zbc=bcth1(ix,iph)
              else
                 zbc=bcth2(ix,iph)
              end if
              call tcspeval(x(ix),th(ith),ph(iph), iselect, zvals, &
                   x,nx,th,nth,ph,nph, 1, 1, 1, f,nx,nth, ier)
              if(nbc.eq.1) then
                 if(abs(zbc-zvals(3)).gt.1.0E-4_fp) &
                      write(m,7703) ix,ith,iph,zbc,zvals(3)
              end if
           end do
        end do
     end do
     !
7703 format(' df/dth BC check @ix=',i2,',ith=',i2,',iph=',i2,': ', &
          2(1x,1pe13.6))
     !
     do iph=1,nph,nph-1
        do ith=1,nth
           do ix=1,nx
              if(iph.eq.1) then
                 zbc=bcph1(ix,ith)
              else
                 zbc=bcph2(ix,ith)
              end if
              call tcspeval(x(ix),th(ith),ph(iph), iselect, zvals, &
                   x,nx,th,nth,ph,nph, 1, 1, 1, f,nx,nth, ier)
              if(nbc.eq.1) then
                 if(abs(zbc-zvals(4)).gt.1.0E-4_fp) &
                      write(m,7704) ix,ith,iph,zbc,zvals(4)
              end if
           end do
        end do
     end do
     !
7704 format(' df/dph BC check @ix=',i2,',ith=',i2,',iph=',i2,': ', &
          2(1x,1pe13.6))
     !
  end if
  !
  call compare3('tcspline',x,nx,th,nth,ph,nph,f,fh,flin, &
       ilinx,ilinth,ilinph, &
       xtest,fxtest,thtest,fthtest,phtest,fphtest,ntest)
  !
  do iph=1,nph
     do ith=1,nth
        do ix=1,nx
           fh(0,ix,ith,iph)=f(1,1,1,ix,ith,iph)
        end do
     end do
  end do
  !
#ifdef _TIMER
  call cpu_time(ztime1)
#endif
  call mktricub(x,nx,th,nth,ph,nph,fh,nx,nth, &
       nbc,bcx1,nbc,bcx2,nth, &
       nbc,bcth1,nbc,bcth2,nx, &
       nbc,bcph1,nbc,bcph2,nx, &
       ilinx,ilinth,ilinph,ier)
  if(ier.ne.0) then
     write(m,*) ' ?? error in pspltest:  dotest3(mktricub)'
     return
  end if
#ifdef _TIMER
  call cpu_time(ztime2)
  write(m,7707) nx,ztime2-ztime1
7707 format(' ==> mktricub setup (nx=',i3,') cpu time (s):',1pe11.4, ' <== ')
#else
  write(m,7707) nx
7707 format(' ==> mktricub setup (nx=',i3,') <== ')
#endif
  !
  call compare3('mktricub', &
       x,nx,th,nth,ph,nph,f,fh,flin, &
       ilinx,ilinth,ilinph, &
       xtest,fxtest,thtest,fthtest,phtest,fphtest,ntest)
  !
  call compare3('piecewise linear', &
       x,nx,th,nth,ph,nph,f,fh,flin, &
       ilinx,ilinth,ilinph, &
       xtest,fxtest,thtest,fthtest,phtest,fphtest,ntest)
  !
  return
end subroutine dotest3
!----------------------------------------------------------
subroutine bset3(fx,nx,fth,nth,fph,nph, &
     bcx1,bcx2,bcth1,bcth2,bcph1,bcph2)
  use psp_precision_mod, only: fp
  !
  !============
  implicit none
  integer nth,nph,nx,iph,ith,ix
  !============
  real(fp) :: fx(nx)                       ! x factor of test fcn
  real(fp) :: fth(nth)                     ! th factor of test fcn
  real(fp) :: fph(nph)                     ! ph factor of test fcn
  real(fp) :: bcx1(nth,nph),bcx2(nth,nph)  ! df/dx bdy conds at x(1),x(nx)
  real(fp) :: bcth1(nx,nph),bcth2(nx,nph)  ! df/dth bdy conds at th(1),th(nth)
  real(fp) :: bcph1(nx,nth),bcph2(nx,nth)  ! df/dph bdy conds at ph(1),ph(nph)
  !
  !  df/dx = 2*exp(2x-1)*(2+sin(th))*(3+sin(ph)) = 2*f
  !
  do iph=1,nph
     do ith=1,nth
        bcx1(ith,iph)=2.0_fp*fx(1)*fth(ith)*fph(iph) ! df/dx @ x(1)
        bcx2(ith,iph)=2.0_fp*fx(nx)*fth(ith)*fph(iph) ! df/dx @ x(nx)
     end do
  end do
  !
  !  df/dth = exp(2x-1)*cos(th)*(3+sin(ph))
  !  cos(0)=cos(2pi)=1
  !
  do iph=1,nph
     do ix=1,nx
        bcth1(ix,iph)=fx(ix)*fph(iph)    ! df/dth @ th=0 = th(1)
        bcth2(ix,iph)=fx(ix)*fph(iph)    ! df/dth @ th=2pi = th(nth)
     end do
  end do
  !
  !  df/dph = exp(2x-1)*(2+sin(th))*cos(ph)
  !  cos(0)=cos(2pi)=1
  !
  do ith=1,nth
     do ix=1,nx
        bcph1(ix,ith)=fx(ix)*fth(ith)    ! df/dph @ ph=0
        bcph2(ix,ith)=fx(ix)*fth(ith)    ! df/dph @ ph=2pi
     end do
  end do
  !
  return
end subroutine bset3
!--------------------------------
!
subroutine compare3(slbl,x,nx,th,nth,ph,nph,f,fh,flin, &
     ilinx,ilinth,ilinph, &
     xtest,fxtest,thtest,fthtest,phtest,fphtest,ntest)
  use psp_precision_mod, only: fp
  !
  !
  !============
  implicit none
  integer nth,nph,ntest,nx,m,iherm,ier,k,j,i,iwarn
  !============
  real(fp) :: fmin,fmax,fdif,fdifr,zdum,zph,zth,zx,ff,fs
#ifdef _TIMER
  real(fp) :: ztime1,ztime2,zctime
#endif
  !============
  common/pspltest_io/ m
  !
  character(len=*) slbl                ! spline coeff routine:  label
  real(fp) :: x(nx),th(nth),ph(nph)        ! indep. coords.
  real(fp) :: f(4,4,4,nx,nth,nph)          ! spline data
  real(fp) :: fh(0:7,nx,nth,nph)           ! hermite data
  real(fp) :: flin(nx,nth,nph)             ! pc linear data
  integer ilinx,ilinth,ilinph       ! even spacing flags
  !
  !  test data grid & data:
  !
  real(fp) :: xtest(ntest),fxtest(ntest)
  real(fp) :: thtest(ntest),fthtest(ntest)
  real(fp) :: phtest(ntest),fphtest(ntest)
  !
  !-------------
  real(fp) :: xpkg(ntest,4),thpkg(ntest,4),phpkg(ntest,4)
  real(fp) :: thvec(ntest),phvec(ntest),fvec(ntest,1)
  !-------------
  !  select spline fcn eval only (no derivatives)
  !
  integer isel(10)
  !
  real(fp) :: fget(10)
  !
  !-------------
  !
  isel = 0
  isel(1) = 1
  !
  iherm=0
  if(slbl.eq.'hermite') iherm=1
  if(slbl.eq.'mktricub') iherm=2
  if(slbl.eq.'piecewise linear') iherm=3
  !
  fmin=1.0E30_fp
  fmax=-1.0E30_fp
  fdif=0.0_fp
  fdifr=0.0_fp
  !
  call genxpkg(nx,x,xpkg,0,1,0,zdum,1,ier)
  call genxpkg(nth,th,thpkg,1,1,0,zdum,1,ier)
  call genxpkg(nph,ph,phpkg,1,1,0,zdum,1,ier)
  !
#ifdef _TIMER
  call cpu_time(ztime1)
#endif
  do k=1,ntest
     zph=phtest(k)
     do j=1,ntest
        phvec(j)=zph
     end do
     do j=1,ntest
        zth=thtest(j)
        do i=1,ntest
           thvec(i)=zth
        end do
        if(iherm.eq.0) then
           call tcspvec(isel,ntest,xtest,thvec,phvec, &
                ntest,fvec,nx,xpkg,nth,thpkg,nph,phpkg, &
                f,nx,nth,iwarn,ier)
        else if(iherm.eq.1) then
           call vecherm3(isel,ntest,xtest,thvec,phvec, &
                ntest,fvec,nx,xpkg,nth,thpkg,nph,phpkg, &
                fh,nx,nth,iwarn,ier)
        else if(iherm.eq.2) then
           call vectricub(isel,ntest,xtest,thvec,phvec, &
                ntest,fvec,nx,xpkg,nth,thpkg,nph,phpkg, &
                fh,nx,nth,iwarn,ier)
        else if(iherm.eq.3) then
           call vecpc3(isel,ntest,xtest,thvec,phvec, &
                ntest,fvec,nx,xpkg,nth,thpkg,nph,phpkg, &
                flin,nx,nth,iwarn,ier)
        end if
        do i=1,ntest
           zx=xtest(i)
           !
           ff=fxtest(i)*fthtest(j)*fphtest(k)
           fmin=min(fmin,ff)
           fmax=max(fmax,ff)
           !
           if(iherm.eq.0) then
              !xx                  call tcspeval(zx,zth,zph, &
              !xx                    isel,fget,x,nx,th,nth,ph,nph, &
              !xx                    ilinx,ilinth,ilinph, &
              !xx                    f,nx,nth,ier)
              fget(1)=fvec(i,1)
           else if(iherm.eq.1) then
              !xx                  call herm3ev(zx,zth,zph,x,nx,th,nth,ph,nph, &
              !xx                    ilinx,ilinth,ilinph, &
              !xx                    fh,nx,nth,isel,fget,ier)
              fget(1)=fvec(i,1)
           else if(iherm.eq.2) then
              !xx                  call evtricub(zx,zth,zph,x,nx,th,nth,ph,nph, &
              !xx                    ilinx,ilinth,ilinph, &
              !xx                    fh,nx,nth,isel,fget,ier)
              fget(1)=fvec(i,1)
           else
              fget(1)=fvec(i,1)
           end if
           !xx               if(ier.ne.0) then
           !xx                  write(m,*) ' ??compare3 ('//slbl//') ier.ne.0 exit'
           !xx                  return
           !xx               end if
           !
           fs=fget(1)
           !
           fdif=max(fdif,abs(ff-fs))
           fdifr=max(fdifr,abs(ff-fs)/(0.5_fp*(ff+fs)))
           !
        end do
     end do
  end do
#ifdef _TIMER
  call cpu_time(ztime2)
#endif
  !
  write(m,1000) slbl,nx,nth,nph,fmin,fmax,fdif,fdifr
1000 format(2x,a,' setup on ',i3,' x ',i3,' x ',i3,' grid'/ &
       '   fmin=',1pe11.4,' fmax=',1pe11.4,' max(|diff|)=',1pe11.4/ &
       '     max(|diff|/f) = ',1pe11.4)
  !
#ifdef _TIMER
  zctime=ztime2-ztime1
  write(m,1001) ntest,ntest,ntest,zctime
#else
  write(m,1002) ntest,ntest,ntest
#endif
1001 format('  ',i3,' x ',i3,' x ',i3,' evaluations, cpu = ',1pe11.4,' (s)')
1002 format('  ',i3,' x ',i3,' x ',i3,' evaluations')
  !
  return
end subroutine compare3
!================================
!  smalltest added: DMC Feb 2010
!  test splines with very small grids: nx=2, nx=3

subroutine smalltest(zctrl)
  use psp_precision_mod, only: fp
  implicit none
  common/pspltest_io/ m
  integer :: m

  !  1d spline & Hermite test

  real(fp) :: zctrl  ! control info (not yet used)
  real(fp) :: x1(2),xpkg(2,4)
  real(fp) :: xx1(3),xxpkg(3,4)
  real(fp) :: xxx1(4),xxxpkg(4,4)
  real(fp) :: s1da(2,2),s1db(2,2),h1da(2,2)
  real(fp) :: s1dc(2,3),s1dd(2,4),szp40
  real(fp) :: test,zbc1,zbc2
  real(fp) :: feval0(5),feval1(5),feval2(5)
  real(fp) :: feval0p(5),feval1p(5),feval2p(5)
  integer :: ier,idum,iwarn,ii,jj

  real(fp), dimension(5) :: xeval = &
       (/ 0.0_fp, 0.25_fp, 0.50_fp, 0.75_fp, 1.00_fp/)
  integer, dimension(3) :: ict01 = (/ 1, 0, 0 /)
  integer, dimension(3) :: ict11 = (/ 0, 1, 0 /)
  integer, dimension(3) :: ict21 = (/ 0, 0, 1 /)

  real(fp), parameter :: ZERO = 0.0_fp
  real(fp), parameter :: ZP40 = 0.40_fp
  real(fp), parameter :: HALF = 0.5_fp
  real(fp), parameter :: ONE = 1.0_fp
  real(fp), parameter :: TWO = 2.0_fp

  real(fp) :: ztol = 1.0E-6_fp
  real(fp) :: ztest

  !---------------------------------
  !  1d splines & hermites
  !  test data with zero at node points and zero, then non-zero BCs

  call psp_tolsum(ONE,ztol,ztest)
  if(ztest.ne.ONE) ztol=ztol*ztol  ! get finer precision tol: r8

  write(m,*) ' '
  write(m,*) '  ...small grid spline tests... '

  !  2 pt grid
  x1(1) = ZERO
  x1(2) = ONE

  !  3 pt grid
  xx1(1) = ZERO
  xx1(2) = ZP40
  xx1(3) = ONE

  !  4 pt grid
  xxx1(1) = ZERO
  xxx1(2) = ZP40
  xxx1(3) = HALF
  xxx1(4) = ONE

  !  2 pt grid
  call genxpkg(2,x1,xpkg,0,1,0,ZERO,0,ier)
  call ckerr('genxpkg#1')
  if(ier.ne.0) return

  !  3 pt grid
  call genxpkg(3,xx1,xxpkg,0,1,0,ZERO,0,ier)
  call ckerr('genxpkg#2')
  if(ier.ne.0) return

  !  4 pt grid
  call genxpkg(4,xxx1,xxxpkg,0,1,0,ZERO,0,ier)
  call ckerr('genxpkg#3')
  if(ier.ne.0) return

  !  all zero spline mini-test

  s1da = ZERO
  ! 1d spline, periodic BC  -->  all zero
  call mkspline(x1,2,s1da,-1,ZERO,0,ZERO,idum,ier) ! periodic -> all zero
  call ckerr('mkspline#1')
  if(ier.ne.0) return

  test=maxval(s1da)
  if(test.ne.ZERO) then
     write(m,*) &
          ' ? short spline: expected ZERO result not obtained.'
  end if

  !  loop over two simple data tests.
  !  from the nx=2 spline results, nx=3 and nx=4 splines are constructed
  !  that should match exactly...

  do jj=1,2
     ! input datasets, all zero, then, not all zero...
     if(jj.eq.1) then
        !
        !  all ZERO data but non-zero boundary conditions
        !
        s1da=ZERO
        s1db=ZERO
        s1dc=ZERO
        s1dd=ZERO
        h1da=ZERO
     else
        !
        !  non-zero data, non-zero boundary conditions
        !
        s1da(1,1)=0.25_fp
        s1da(1,2)=0.75_fp
        s1db(1,1:2)=s1da(1,1:2)
        h1da(1,1:2)=s1da(1,1:2)
     end if

     !  3 pt copy
     s1dc(1,1)=s1da(1,1)
     s1dc(1,3)=s1da(1,2)

     !  4 pt copy
     s1dd(1,1)=s1da(1,1)
     s1dd(1,4)=s1da(1,2)

     ! 1d spline, non-zero f' boundary conditions

     call mkspline(x1,2,s1da,1,ONE,1,TWO,idum,ier)
     call ckerr('mkspline#2')
     if(ier.ne.0) return

     !  evaluate f, f', f''
     !  augment data for nx=3 and nx=4 splines

     call vecspline(ict01,1,ZP40,1,szp40,2,xpkg,s1da,iwarn,ier)
     call ckerr('vecspline#0')
     if(ier.ne.0) return

     s1dc(1,2) = szp40
     s1dd(1,2) = szp40

     call vecspline(ict01,5,xeval,5,feval0,2,xpkg,s1da,iwarn,ier)
     call ckerr('vecspline#1')
     if(ier.ne.0) return

     s1dd(1,3) = feval0(3)

     call vecspline(ict11,5,xeval,5,feval1,2,xpkg,s1da,iwarn,ier)
     call ckerr('vecspline#2')
     if(ier.ne.0) return

     call vecspline(ict21,5,xeval,5,feval2,2,xpkg,s1da,iwarn,ier)
     call ckerr('vecspline#3')
     if(ier.ne.0) return

     !  print results

     write(m,*) ' '
     write(m,*) " nx=2 1d spline, f' boundary conditions set:"
     write(m,*) &
          "  #    x            f(x)         f'(x)        f''(x)"
     do ii=1,5
        write(m,'(3x,i1,2x,4(1x,1pe12.5))') ii,xeval(ii), &
             feval0(ii),feval1(ii),feval2(ii)
     end do

     ! save...
     feval0p = feval0
     feval1p = feval1
     feval2p = feval2
     ! stash 2nd derivative BC values...
     zbc1=feval2(1)
     zbc2=feval2(5)
     ! 1d hermite with same f' BCs

     !  Hermite test
     h1da(2,1)=ONE
     h1da(2,2)=TWO
     call akherm1p(x1,2,h1da,idum,2,ier)
     call ckerr('akherm1p#2')
     if(ier.ne.0) return

     !  Hermite evaluation
     call vecherm1(ict01,5,xeval,5,feval0,2,xpkg,h1da,iwarn,ier)
     call ckerr('vecherm1#1')
     if(ier.ne.0) return

     call vecherm1(ict11,5,xeval,5,feval1,2,xpkg,h1da,iwarn,ier)
     call ckerr('vecherm1#1')
     if(ier.ne.0) return

     !  print results
     write(m,*) ' '
     write(m,*) " 1d hermite, f' boundary conditions (BCs) set:"
     write(m,*) "  #    x            f(x)         f'(x)"
     do ii=1,5
        write(m,'(3x,i1,2x,3(1x,1pe12.5))') ii,xeval(ii), &
             feval0(ii),feval1(ii)
     end do

     !  compare
     call ckdiff(2)

     !  nx=3 spline test

     call mkspline(xx1,3,s1dc,1,ONE,1,TWO,idum,ier)
     call ckerr('mkspline')
     if(ier.ne.0) return

     !  evaluate f, f', f''

     call vecspline(ict01,5,xeval,5,feval0,3,xxpkg,s1dc,iwarn, &
          ier)
     call ckerr('vecspline')
     if(ier.ne.0) return

     call vecspline(ict11,5,xeval,5,feval1,3,xxpkg,s1dc,iwarn, &
          ier)
     call ckerr('vecspline')
     if(ier.ne.0) return

     call vecspline(ict21,5,xeval,5,feval2,3,xxpkg,s1dc,iwarn, &
          ier)
     call ckerr('vecspline')
     if(ier.ne.0) return

     !  print results

     write(m,*) ' '
     write(m,*) " nx=3 1d spline, f' LHS BC, f' RHS BC"
     write(m,*) &
          "  #    x            f(x)         f'(x)        f''(x)"
     do ii=1,5
        write(m,'(3x,i1,2x,4(1x,1pe12.5))') ii,xeval(ii), &
             feval0(ii),feval1(ii),feval2(ii)
     end do

     !  compare
     call ckdiff(3)

     !  nx=4 spline test

     call mkspline(xxx1,4,s1dd,1,ONE,1,TWO,idum,ier)
     call ckerr('mkspline')
     if(ier.ne.0) return

     !  evaluate f, f', f''

     call vecspline(ict01,5,xeval,5,feval0,4,xxxpkg,s1dd,iwarn, &
          ier)
     call ckerr('vecspline')
     if(ier.ne.0) return

     call vecspline(ict11,5,xeval,5,feval1,4,xxxpkg,s1dd,iwarn, &
          ier)
     call ckerr('vecspline')
     if(ier.ne.0) return

     call vecspline(ict21,5,xeval,5,feval2,4,xxxpkg,s1dd,iwarn, &
          ier)
     call ckerr('vecspline')
     if(ier.ne.0) return

     !  print results

     write(m,*) ' '
     write(m,*) " nx=4 1d spline, f' LHS BC, f' RHS BC"
     write(m,*) &
          "  #    x            f(x)         f'(x)        f''(x)"
     do ii=1,5
        write(m,'(3x,i1,2x,4(1x,1pe12.5))') ii,xeval(ii), &
             feval0(ii),feval1(ii),feval2(ii)
     end do

     !  compare
     call ckdiff(3)
     ! 1d spline, f' LHS, f'' RHS

     call mkspline(x1,2,s1db,1,ONE,2,zbc2,idum,ier)
     call ckerr('mkspline#3')
     if(ier.ne.0) return

     !  evaluate f, f', f''

     call vecspline(ict01,5,xeval,5,feval0,2,xpkg,s1db,iwarn,ier)
     call ckerr('vecspline#4')
     if(ier.ne.0) return

     call vecspline(ict11,5,xeval,5,feval1,2,xpkg,s1db,iwarn,ier)
     call ckerr('vecspline#5')
     if(ier.ne.0) return

     call vecspline(ict21,5,xeval,5,feval2,2,xpkg,s1db,iwarn,ier)
     call ckerr('vecspline#6')
     if(ier.ne.0) return

     !  print results

     write(m,*) ' '
     write(m,*) " nx=2 1d spline, f' LHS BC, f'' RHS BC"
     write(m,*) &
          "  #    x            f(x)         f'(x)        f''(x)"
     do ii=1,5
        write(m,'(3x,i1,2x,4(1x,1pe12.5))') ii,xeval(ii), &
             feval0(ii),feval1(ii),feval2(ii)
     end do

     !  compare
     call ckdiff(3)

     !  nx=3 spline test

     !  nx=3 spline test

     call mkspline(xx1,3,s1dc,1,ONE,2,zbc2,idum,ier)
     call ckerr('mkspline')
     if(ier.ne.0) return

     !  evaluate f, f', f''

     call vecspline(ict01,5,xeval,5,feval0,3,xxpkg,s1dc,iwarn, &
          ier)
     call ckerr('vecspline')
     if(ier.ne.0) return

     call vecspline(ict11,5,xeval,5,feval1,3,xxpkg,s1dc,iwarn, &
          ier)
     call ckerr('vecspline')
     if(ier.ne.0) return

     call vecspline(ict21,5,xeval,5,feval2,3,xxpkg,s1dc,iwarn, &
          ier)
     call ckerr('vecspline')
     if(ier.ne.0) return

     !  print results

     write(m,*) ' '
     write(m,*) " nx=3 1d spline, f' LHS BC, f'' RHS BC"
     write(m,*) &
          "  #    x            f(x)         f'(x)        f''(x)"
     do ii=1,5
        write(m,'(3x,i1,2x,4(1x,1pe12.5))') ii,xeval(ii), &
             feval0(ii),feval1(ii),feval2(ii)
     end do

     !  compare
     call ckdiff(3)

     !  nx=4 spline test

     call mkspline(xxx1,4,s1dd,1,ONE,2,zbc2,idum,ier)
     call ckerr('mkspline')
     if(ier.ne.0) return

     !  evaluate f, f', f''

     call vecspline(ict01,5,xeval,5,feval0,4,xxxpkg,s1dd,iwarn, &
          ier)
     call ckerr('vecspline')
     if(ier.ne.0) return

     call vecspline(ict11,5,xeval,5,feval1,4,xxxpkg,s1dd,iwarn, &
          ier)
     call ckerr('vecspline')
     if(ier.ne.0) return

     call vecspline(ict21,5,xeval,5,feval2,4,xxxpkg,s1dd,iwarn, &
          ier)
     call ckerr('vecspline')
     if(ier.ne.0) return

     !  print results

     write(m,*) ' '
     write(m,*) " nx=4 1d spline, f' LHS BC, f'' RHS BC"
     write(m,*) &
          "  #    x            f(x)         f'(x)        f''(x)"
     do ii=1,5
        write(m,'(3x,i1,2x,4(1x,1pe12.5))') ii,xeval(ii), &
             feval0(ii),feval1(ii),feval2(ii)
     end do

     !  compare
     call ckdiff(3)

     ! 1d spline, f'' LHS, f' RHS

     call mkspline(x1,2,s1db,2,zbc1,1,TWO,idum,ier)
     call ckerr('mkspline#4')
     if(ier.ne.0) return

     !  evaluate f, f', f''

     call vecspline(ict01,5,xeval,5,feval0,2,xpkg,s1db,iwarn,ier)
     call ckerr('vecspline#7')
     if(ier.ne.0) return

     call vecspline(ict11,5,xeval,5,feval1,2,xpkg,s1db,iwarn,ier)
     call ckerr('vecspline#8')
     if(ier.ne.0) return

     call vecspline(ict21,5,xeval,5,feval2,2,xpkg,s1db,iwarn,ier)
     call ckerr('vecspline#9')
     if(ier.ne.0) return

     !  print results

     write(m,*) ' '
     write(m,*) " nx=2 1d spline, f'' LHS BC, f' RHS BC"
     write(m,*) &
          "  #    x            f(x)         f'(x)        f''(x)"
     do ii=1,5
        write(m,'(3x,i1,2x,4(1x,1pe12.5))') ii,xeval(ii), &
             feval0(ii),feval1(ii),feval2(ii)
     end do

     !  compare
     call ckdiff(3)

     !  nx=3 spline test

     !  nx=3 spline test

     call mkspline(xx1,3,s1dc,2,zbc1,1,TWO,idum,ier)
     call ckerr('mkspline')
     if(ier.ne.0) return

     !  evaluate f, f', f''

     call vecspline(ict01,5,xeval,5,feval0,3,xxpkg,s1dc,iwarn, &
          ier)
     call ckerr('vecspline')
     if(ier.ne.0) return

     call vecspline(ict11,5,xeval,5,feval1,3,xxpkg,s1dc,iwarn, &
          ier)
     call ckerr('vecspline')
     if(ier.ne.0) return

     call vecspline(ict21,5,xeval,5,feval2,3,xxpkg,s1dc,iwarn, &
          ier)
     call ckerr('vecspline')
     if(ier.ne.0) return

     !  print results

     write(m,*) ' '
     write(m,*) " nx=3 1d spline, f'' LHS BC, f' RHS BC"
     write(m,*) &
          "  #    x            f(x)         f'(x)        f''(x)"
     do ii=1,5
        write(m,'(3x,i1,2x,4(1x,1pe12.5))') ii,xeval(ii), &
             feval0(ii),feval1(ii),feval2(ii)
     end do

     !  compare
     call ckdiff(3)

     !  nx=4 spline test

     call mkspline(xxx1,4,s1dd,2,zbc1,1,TWO,idum,ier)
     call ckerr('mkspline')
     if(ier.ne.0) return

     !  evaluate f, f', f''

     call vecspline(ict01,5,xeval,5,feval0,4,xxxpkg,s1dd,iwarn, &
          ier)
     call ckerr('vecspline')
     if(ier.ne.0) return

     call vecspline(ict11,5,xeval,5,feval1,4,xxxpkg,s1dd,iwarn, &
          ier)
     call ckerr('vecspline')
     if(ier.ne.0) return

     call vecspline(ict21,5,xeval,5,feval2,4,xxxpkg,s1dd,iwarn, &
          ier)
     call ckerr('vecspline')
     if(ier.ne.0) return

     !  print results

     write(m,*) ' '
     write(m,*) " nx=4 1d spline, f'' LHS BC, f' RHS BC"
     write(m,*) &
          "  #    x            f(x)         f'(x)        f''(x)"
     do ii=1,5
        write(m,'(3x,i1,2x,4(1x,1pe12.5))') ii,xeval(ii), &
             feval0(ii),feval1(ii),feval2(ii)
     end do

     !  compare
     call ckdiff(3)
     ! 1d spline, f'' LHS, f'' RHS

     call mkspline(x1,2,s1db,2,zbc1,2,zbc2,idum,ier)
     call ckerr('mkspline#5')
     if(ier.ne.0) return

     !  evaluate f, f', f''

     call vecspline(ict01,5,xeval,5,feval0,2,xpkg,s1db,iwarn,ier)
     call ckerr('vecspline#10')
     if(ier.ne.0) return

     call vecspline(ict11,5,xeval,5,feval1,2,xpkg,s1db,iwarn,ier)
     call ckerr('vecspline#11')
     if(ier.ne.0) return

     call vecspline(ict21,5,xeval,5,feval2,2,xpkg,s1db,iwarn,ier)
     call ckerr('vecspline#12')
     if(ier.ne.0) return

     !  print results

     write(m,*) ' '
     write(m,*) " nx=2 1d spline, f'' LHS BC, f'' RHS BC"
     write(m,*) &
          "  #    x            f(x)         f'(x)        f''(x)"
     do ii=1,5
        write(m,'(3x,i1,2x,4(1x,1pe12.5))') ii,xeval(ii), &
             feval0(ii),feval1(ii),feval2(ii)
     end do

     !  compare
     call ckdiff(3)

     !  nx=3 spline test

     call mkspline(xx1,3,s1dc,2,zbc1,2,zbc2,idum,ier)
     call ckerr('mkspline')
     if(ier.ne.0) return

     !  evaluate f, f', f''

     call vecspline(ict01,5,xeval,5,feval0,3,xxpkg,s1dc,iwarn, &
          ier)
     call ckerr('vecspline')
     if(ier.ne.0) return

     call vecspline(ict11,5,xeval,5,feval1,3,xxpkg,s1dc,iwarn, &
          ier)
     call ckerr('vecspline')
     if(ier.ne.0) return

     call vecspline(ict21,5,xeval,5,feval2,3,xxpkg,s1dc,iwarn, &
          ier)
     call ckerr('vecspline')
     if(ier.ne.0) return

     !  print results

     write(m,*) ' '
     write(m,*) " nx=3 1d spline, f'' LHS BC, f'' RHS BC"
     write(m,*) &
          "  #    x            f(x)         f'(x)        f''(x)"
     do ii=1,5
        write(m,'(3x,i1,2x,4(1x,1pe12.5))') ii,xeval(ii), &
             feval0(ii),feval1(ii),feval2(ii)
     end do

     !  compare
     call ckdiff(3)

     !  nx=4 spline test

     call mkspline(xxx1,4,s1dd,2,zbc1,2,zbc2,idum,ier)
     call ckerr('mkspline')
     if(ier.ne.0) return

     !  evaluate f, f', f''

     call vecspline(ict01,5,xeval,5,feval0,4,xxxpkg,s1dd,iwarn, &
          ier)
     call ckerr('vecspline')
     if(ier.ne.0) return

     call vecspline(ict11,5,xeval,5,feval1,4,xxxpkg,s1dd,iwarn, &
          ier)
     call ckerr('vecspline')
     if(ier.ne.0) return

     call vecspline(ict21,5,xeval,5,feval2,4,xxxpkg,s1dd,iwarn, &
          ier)
     call ckerr('vecspline')
     if(ier.ne.0) return

     !  print results

     write(m,*) ' '
     write(m,*) " nx=4 1d spline, f'' LHS BC, f'' RHS BC"
     write(m,*) &
          "  #    x            f(x)         f'(x)        f''(x)"
     do ii=1,5
        write(m,'(3x,i1,2x,4(1x,1pe12.5))') ii,xeval(ii), &
             feval0(ii),feval1(ii),feval2(ii)
     end do

     !  compare
     call ckdiff(3)

     !  nx=4 "not a knot" spline test

     call mkspline(xxx1,4,s1dd,0,ZERO,0,ZERO,idum,ier)
     call ckerr('mkspline')
     if(ier.ne.0) return

     !  evaluate f, f', f''

     call vecspline(ict01,5,xeval,5,feval0,4,xxxpkg,s1dd,iwarn, &
          ier)
     call ckerr('vecspline')
     if(ier.ne.0) return

     call vecspline(ict11,5,xeval,5,feval1,4,xxxpkg,s1dd,iwarn, &
          ier)
     call ckerr('vecspline')
     if(ier.ne.0) return

     call vecspline(ict21,5,xeval,5,feval2,4,xxxpkg,s1dd,iwarn, &
          ier)
     call ckerr('vecspline')
     if(ier.ne.0) return

     !  print results

     write(m,*) ' '
     write(m,*) " nx=4 1d spline, not-a-knot BC, both sides."
     write(m,*) &
          "  #    x            f(x)         f'(x)        f''(x)"
     do ii=1,5
        write(m,'(3x,i1,2x,4(1x,1pe12.5))') ii,xeval(ii), &
             feval0(ii),feval1(ii),feval2(ii)
     end do

     !  compare
     call ckdiff(3)
  end do

  ! nx=3 & nx=4 periodic spline comparison

  write(m,*) ' '
  write(m,*) ' ---------------------------------'
  write(m,*) ' nx=3 & nx=4 periodic spline test.'

  s1dc(1,1)=ONE
  s1dc(1,2)=TWO
  s1dc(1,3)=TWO*TWO

  s1dd(1,1)=ONE
  s1dd(1,2)=TWO
  s1dd(1,4)=TWO*TWO

  !  nx=3 periodic

  call mkspline(xx1,3,s1dc,-1,ZERO,-1,ZERO,idum,ier)
  call ckerr('mkspline')
  if(ier.ne.0) return

  !  evaluate f, f', f''

  call vecspline(ict01,5,xeval,5,feval0,3,xxpkg,s1dc,iwarn, &
       ier)
  call ckerr('vecspline')
  if(ier.ne.0) return

  s1dd(1,3)=feval0(3)  ! copy result for nx=4 test data

  call vecspline(ict11,5,xeval,5,feval1,3,xxpkg,s1dc,iwarn, &
       ier)
  call ckerr('vecspline')
  if(ier.ne.0) return

  call vecspline(ict21,5,xeval,5,feval2,3,xxpkg,s1dc,iwarn, &
       ier)
  call ckerr('vecspline')
  if(ier.ne.0) return

  !  print results

  write(m,*) ' '
  write(m,*) " nx=3 periodic spline "
  write(m,*) &
       "  #    x            f(x)         f'(x)        f''(x)"
  do ii=1,5
     write(m,'(3x,i1,2x,4(1x,1pe12.5))') ii,xeval(ii), &
          feval0(ii),feval1(ii),feval2(ii)
  end do

  ! save...
  feval0p = feval0
  feval1p = feval1
  feval2p = feval2

  ! nx=4 periodic

  call mkspline(xxx1,4,s1dd,-1,ZERO,-1,ZERO,idum,ier)
  call ckerr('mkspline')
  if(ier.ne.0) return

  !  evaluate f, f', f''

  call vecspline(ict01,5,xeval,5,feval0,4,xxxpkg,s1dd,iwarn, &
       ier)
  call ckerr('vecspline')
  if(ier.ne.0) return

  call vecspline(ict11,5,xeval,5,feval1,4,xxxpkg,s1dd,iwarn, &
       ier)
  call ckerr('vecspline')
  if(ier.ne.0) return

  call vecspline(ict21,5,xeval,5,feval2,4,xxxpkg,s1dd,iwarn, &
       ier)
  call ckerr('vecspline')
  if(ier.ne.0) return

  !  print results

  write(m,*) ' '
  write(m,*) " nx=4 periodic spline "
  write(m,*) &
       "  #    x            f(x)         f'(x)        f''(x)"
  do ii=1,5
     write(m,'(3x,i1,2x,4(1x,1pe12.5))') ii,xeval(ii), &
          feval0(ii),feval1(ii),feval2(ii)
  end do

  !  compare
  call ckdiff(3)

  !  nx=3 two sided not-a-knot

  write(m,*) ' '
  write(m,*) ' ---------------------------------'
  write(m,*) ' nx=3 & nx=4 two sided not-a-knot.'

  call mkspline(xx1,3,s1dc,0,ZERO,0,ZERO,idum,ier)
  call ckerr('mkspline')
  if(ier.ne.0) return

  !  evaluate f, f', f''

  call vecspline(ict01,5,xeval,5,feval0,3,xxpkg,s1dc,iwarn, &
       ier)
  call ckerr('vecspline')
  if(ier.ne.0) return

  s1dd(1,3)=feval0(3)  ! copy result for nx=4 test data

  call vecspline(ict11,5,xeval,5,feval1,3,xxpkg,s1dc,iwarn, &
       ier)
  call ckerr('vecspline')
  if(ier.ne.0) return

  call vecspline(ict21,5,xeval,5,feval2,3,xxpkg,s1dc,iwarn, &
       ier)
  call ckerr('vecspline')
  if(ier.ne.0) return

  !  print results

  write(m,*) ' '
  write(m,*) " nx=3 double not-a-knot spline "
  write(m,*) &
       "  #    x            f(x)         f'(x)        f''(x)"
  do ii=1,5
     write(m,'(3x,i1,2x,4(1x,1pe12.5))') ii,xeval(ii), &
          feval0(ii),feval1(ii),feval2(ii)
  end do

  ! save...
  feval0p = feval0
  feval1p = feval1
  feval2p = feval2

  ! nx=4 comparison spline

  call mkspline(xxx1,4,s1dd,0,ZERO,0,ZERO,idum,ier)
  call ckerr('mkspline')
  if(ier.ne.0) return

  !  evaluate f, f', f''

  call vecspline(ict01,5,xeval,5,feval0,4,xxxpkg,s1dd,iwarn, &
       ier)
  call ckerr('vecspline')
  if(ier.ne.0) return

  call vecspline(ict11,5,xeval,5,feval1,4,xxxpkg,s1dd,iwarn, &
       ier)
  call ckerr('vecspline')
  if(ier.ne.0) return

  call vecspline(ict21,5,xeval,5,feval2,4,xxxpkg,s1dd,iwarn, &
       ier)
  call ckerr('vecspline')
  if(ier.ne.0) return

  !  print results

  write(m,*) ' '
  write(m,*) " nx=4 double not-a-knot"
  write(m,*) &
       "  #    x            f(x)         f'(x)        f''(x)"
  do ii=1,5
     write(m,'(3x,i1,2x,4(1x,1pe12.5))') ii,xeval(ii), &
          feval0(ii),feval1(ii),feval2(ii)
  end do

  !  compare
  call ckdiff(3)

  !  nx=3 RHS not-a-knot, LHS f' BC

  write(m,*) ' '
  write(m,*) ' ---------------------------------'
  write(m,*) " nx=3 RHS not-a-knot & LHS f' BC..."

  call mkspline(xx1,3,s1dc,1,ONE,0,ZERO,idum,ier)
  call ckerr('mkspline')
  if(ier.ne.0) return

  !  evaluate f, f', f''

  call vecspline(ict01,5,xeval,5,feval0,3,xxpkg,s1dc,iwarn, &
       ier)
  call ckerr('vecspline')
  if(ier.ne.0) return

  s1dd(1,3)=feval0(3)  ! copy result for nx=4 test data

  call vecspline(ict11,5,xeval,5,feval1,3,xxpkg,s1dc,iwarn, &
       ier)
  call ckerr('vecspline')
  if(ier.ne.0) return

  call vecspline(ict21,5,xeval,5,feval2,3,xxpkg,s1dc,iwarn, &
       ier)
  call ckerr('vecspline')
  if(ier.ne.0) return

  !  print results

  write(m,*) ' '
  write(m,*) " nx=3 RHS not-a-knot LHS f' BC"
  write(m,*) &
       "  #    x            f(x)         f'(x)        f''(x)"
  do ii=1,5
     write(m,'(3x,i1,2x,4(1x,1pe12.5))') ii,xeval(ii), &
          feval0(ii),feval1(ii),feval2(ii)
  end do

  ! save...
  feval0p = feval0
  feval1p = feval1
  feval2p = feval2

  zbc1=feval2(1)
  zbc2=feval2(5)

  ! nx=4 comparison spline

  call mkspline(xxx1,4,s1dd,2,zbc1,2,zbc2,idum,ier)
  call ckerr('mkspline')
  if(ier.ne.0) return

  !  evaluate f, f', f''

  call vecspline(ict01,5,xeval,5,feval0,4,xxxpkg,s1dd,iwarn, &
       ier)
  call ckerr('vecspline')
  if(ier.ne.0) return

  call vecspline(ict11,5,xeval,5,feval1,4,xxxpkg,s1dd,iwarn, &
       ier)
  call ckerr('vecspline')
  if(ier.ne.0) return

  call vecspline(ict21,5,xeval,5,feval2,4,xxxpkg,s1dd,iwarn, &
       ier)
  call ckerr('vecspline')
  if(ier.ne.0) return

  !  print results

  write(m,*) ' '
  write(m,*) " nx=4 comparison spline"
  write(m,*) &
       "  #    x            f(x)         f'(x)        f''(x)"
  do ii=1,5
     write(m,'(3x,i1,2x,4(1x,1pe12.5))') ii,xeval(ii), &
          feval0(ii),feval1(ii),feval2(ii)
  end do

  !  compare
  call ckdiff(3)

  !  nx=3 RHS not-a-knot, LHS f'' BC

  write(m,*) ' '
  write(m,*) ' ---------------------------------'
  write(m,*) " nx=3 RHS not-a-knot & LHS f'' BC..."

  call mkspline(xx1,3,s1dc,2,TWO*TWO*TWO,0,ZERO,idum,ier)
  call ckerr('mkspline')
  if(ier.ne.0) return

  !  evaluate f, f', f''

  call vecspline(ict01,5,xeval,5,feval0,3,xxpkg,s1dc,iwarn, &
       ier)
  call ckerr('vecspline')
  if(ier.ne.0) return

  s1dd(1,3)=feval0(3)  ! copy result for nx=4 test data

  call vecspline(ict11,5,xeval,5,feval1,3,xxpkg,s1dc,iwarn, &
       ier)
  call ckerr('vecspline')
  if(ier.ne.0) return

  call vecspline(ict21,5,xeval,5,feval2,3,xxpkg,s1dc,iwarn, &
       ier)
  call ckerr('vecspline')
  if(ier.ne.0) return

  !  print results

  write(m,*) ' '
  write(m,*) " nx=3 RHS not-a-knot LHS f'' BC"
  write(m,*) &
       "  #    x            f(x)         f'(x)        f''(x)"
  do ii=1,5
     write(m,'(3x,i1,2x,4(1x,1pe12.5))') ii,xeval(ii), &
          feval0(ii),feval1(ii),feval2(ii)
  end do

  ! save...
  feval0p = feval0
  feval1p = feval1
  feval2p = feval2

  zbc1=feval2(1)
  zbc2=feval2(5)

  ! nx=4 comparison spline

  call mkspline(xxx1,4,s1dd,2,zbc1,2,zbc2,idum,ier)
  call ckerr('mkspline')
  if(ier.ne.0) return

  !  evaluate f, f', f''

  call vecspline(ict01,5,xeval,5,feval0,4,xxxpkg,s1dd,iwarn, &
       ier)
  call ckerr('vecspline')
  if(ier.ne.0) return

  call vecspline(ict11,5,xeval,5,feval1,4,xxxpkg,s1dd,iwarn, &
       ier)
  call ckerr('vecspline')
  if(ier.ne.0) return

  call vecspline(ict21,5,xeval,5,feval2,4,xxxpkg,s1dd,iwarn, &
       ier)
  call ckerr('vecspline')
  if(ier.ne.0) return

  !  print results

  write(m,*) ' '
  write(m,*) " nx=4 comparison spline"
  write(m,*) &
       "  #    x            f(x)         f'(x)        f''(x)"
  do ii=1,5
     write(m,'(3x,i1,2x,4(1x,1pe12.5))') ii,xeval(ii), &
          feval0(ii),feval1(ii),feval2(ii)
  end do

  !  compare
  call ckdiff(3)


  !  nx=3 LHS not-a-knot, RHS f' BC

  write(m,*) ' '
  write(m,*) ' ---------------------------------'
  write(m,*) " nx=3 LHS not-a-knot & RHS f' BC..."

  call mkspline(xx1,3,s1dc,0,ZERO,1,ONE,idum,ier)
  call ckerr('mkspline')
  if(ier.ne.0) return

  !  evaluate f, f', f''

  call vecspline(ict01,5,xeval,5,feval0,3,xxpkg,s1dc,iwarn, &
       ier)
  call ckerr('vecspline')
  if(ier.ne.0) return

  s1dd(1,3)=feval0(3)  ! copy result for nx=4 test data

  call vecspline(ict11,5,xeval,5,feval1,3,xxpkg,s1dc,iwarn, &
       ier)
  call ckerr('vecspline')
  if(ier.ne.0) return

  call vecspline(ict21,5,xeval,5,feval2,3,xxpkg,s1dc,iwarn, &
       ier)
  call ckerr('vecspline')
  if(ier.ne.0) return

  !  print results

  write(m,*) ' '
  write(m,*) " nx=3 LHS not-a-knot RHS f' BC"
  write(m,*) &
       "  #    x            f(x)         f'(x)        f''(x)"
  do ii=1,5
     write(m,'(3x,i1,2x,4(1x,1pe12.5))') ii,xeval(ii), &
          feval0(ii),feval1(ii),feval2(ii)
  end do

  ! save...
  feval0p = feval0
  feval1p = feval1
  feval2p = feval2

  zbc1=feval2(1)
  zbc2=feval2(5)

  ! nx=4 comparison spline

  call mkspline(xxx1,4,s1dd,2,zbc1,2,zbc2,idum,ier)
  call ckerr('mkspline')
  if(ier.ne.0) return

  !  evaluate f, f', f''

  call vecspline(ict01,5,xeval,5,feval0,4,xxxpkg,s1dd,iwarn, &
       ier)
  call ckerr('vecspline')
  if(ier.ne.0) return

  call vecspline(ict11,5,xeval,5,feval1,4,xxxpkg,s1dd,iwarn, &
       ier)
  call ckerr('vecspline')
  if(ier.ne.0) return

  call vecspline(ict21,5,xeval,5,feval2,4,xxxpkg,s1dd,iwarn, &
       ier)
  call ckerr('vecspline')
  if(ier.ne.0) return

  !  print results

  write(m,*) ' '
  write(m,*) " nx=4 comparison spline"
  write(m,*) &
       "  #    x            f(x)         f'(x)        f''(x)"
  do ii=1,5
     write(m,'(3x,i1,2x,4(1x,1pe12.5))') ii,xeval(ii), &
          feval0(ii),feval1(ii),feval2(ii)
  end do

  !  compare
  call ckdiff(3)

  !  nx=3 LHS not-a-knot, RHS f'' BC

  write(m,*) ' '
  write(m,*) ' ---------------------------------'
  write(m,*) " nx=3 LHS not-a-knot & RHS f'' BC..."

  call mkspline(xx1,3,s1dc,0,ZERO,2,TWO*TWO*TWO,idum,ier)
  call ckerr('mkspline')
  if(ier.ne.0) return

  !  evaluate f, f', f''

  call vecspline(ict01,5,xeval,5,feval0,3,xxpkg,s1dc,iwarn, &
       ier)
  call ckerr('vecspline')
  if(ier.ne.0) return

  s1dd(1,3)=feval0(3)  ! copy result for nx=4 test data

  call vecspline(ict11,5,xeval,5,feval1,3,xxpkg,s1dc,iwarn, &
       ier)
  call ckerr('vecspline')
  if(ier.ne.0) return

  call vecspline(ict21,5,xeval,5,feval2,3,xxpkg,s1dc,iwarn, &
       ier)
  call ckerr('vecspline')
  if(ier.ne.0) return

  !  print results

  write(m,*) ' '
  write(m,*) " nx=3 LHS not-a-knot RHS f'' BC"
  write(m,*) &
       "  #    x            f(x)         f'(x)        f''(x)"
  do ii=1,5
     write(m,'(3x,i1,2x,4(1x,1pe12.5))') ii,xeval(ii), &
          feval0(ii),feval1(ii),feval2(ii)
  end do

  ! save...
  feval0p = feval0
  feval1p = feval1
  feval2p = feval2

  zbc1=feval2(1)
  zbc2=feval2(5)

  ! nx=4 comparison spline

  call mkspline(xxx1,4,s1dd,2,zbc1,2,zbc2,idum,ier)
  call ckerr('mkspline')
  if(ier.ne.0) return

  !  evaluate f, f', f''

  call vecspline(ict01,5,xeval,5,feval0,4,xxxpkg,s1dd,iwarn, &
       ier)
  call ckerr('vecspline')
  if(ier.ne.0) return

  call vecspline(ict11,5,xeval,5,feval1,4,xxxpkg,s1dd,iwarn, &
       ier)
  call ckerr('vecspline')
  if(ier.ne.0) return

  call vecspline(ict21,5,xeval,5,feval2,4,xxxpkg,s1dd,iwarn, &
       ier)
  call ckerr('vecspline')
  if(ier.ne.0) return

  !  print results

  write(m,*) ' '
  write(m,*) " nx=4 comparison spline"
  write(m,*) &
       "  #    x            f(x)         f'(x)        f''(x)"
  do ii=1,5
     write(m,'(3x,i1,2x,4(1x,1pe12.5))') ii,xeval(ii), &
          feval0(ii),feval1(ii),feval2(ii)
  end do

  !  compare
  call ckdiff(3)


CONTAINS

  subroutine ckerr(slbl)
    character(len=*), intent(in) :: slbl

    if(ier.ne.0) then
       write(m,*) ' ? smalltest call returned error: '//trim(slbl)
    end if
  end subroutine ckerr

  subroutine ckdiff(iordp1)

    integer, intent(in) :: iordp1  ! =2: ck f & f'; =3: ck f,f',f''

    real(fp) :: zdiff

    zdiff = maxval(abs(feval0-feval0p))
    if(zdiff.gt.ztol) then
       write(m,99) " *** f value differences exceed tolerance:", &
            ztol,zdiff
    end if

    zdiff = maxval(abs(feval1-feval1p))/10
    if(zdiff.gt.ztol) then
       write(m,99) " *** f' differences exceed tolerance:", &
            ztol,zdiff
    end if

    if(iordp1.lt.3) return

    zdiff = maxval(abs(feval2-feval2p))/100
    if(zdiff.gt.ztol) then
       write(m,99) " *** f'' differences exceed tolerance:", &
            ztol,zdiff
    end if

99  format(1x,a,2(1x,1pe12.5))

  end subroutine ckdiff

end subroutine smalltest
