!  tcspline -- dmc 20 Jan 1999
!
!  set up coefficients for bicubic spline with following BC's:
!  * LHS and RHS handled as in cubspl.f90 for 1st coordinate
!  * derivatives periodic in second coordinate (use pspline.f90)
!
! workspace:
!  if phi bdy cond. is periodic, not-a-knot, df/dphi = 0 everywhere,
!  or d2f/dphi2 = 0 everywhere, then the phi boundary condition is
!  "linear" and a workspace of size at least:
!
!     nwk = 20*inx*inth + 10*max(inx,inth,inph)
!
!  will suffice.
!
!  if the phi bdy cond. involves specification of df/dphi .ne. 0 or
!  d2f/dphi .ne. 0 at any (x,theta) grid point, then, the phi boundary
!  condition is "non-linear", a correction step is needed, and a workspace
!  of size at least:
!
!     nwk = 16*inx*inth*inph
!
!  is required.
!
subroutine tcspline(x,inx,th,inth,ph,inph,fspl,inf4,inf5, &
     ibcxmin,bcxmin,ibcxmax,bcxmax,inb1x, &
     ibcthmin,bcthmin,ibcthmax,bcthmax,inb1th, &
     ibcphmin,bcphmin,ibcphmax,bcphmax,inb1ph, &
     wk,nwk,ilinx,ilinth,ilinph,ier)
  use psp_precision_mod, only: fp
  !
  !============
  implicit none
  integer inth,inph,inf4,inf5,ibcxmin,ibcxmax,inb1x,ibcthmin
  integer ibcthmax,inb1th,ibcphmin,ibcphmax,inb1ph,nwk,ilinx
  integer ilinth,ilinph,ier,inx,iflg,ith,ix,itest,ierx,ierth
  integer ierph,iaspl2,iabcx1,iabcx2,iabcth1,iabcth2,iawk,inwk
  integer iph,inpho,ic1,ic2,ibcphmina,ibcphmaxa,iabcph1,iabcph2
  integer iaccoef,i,iskip1,iskip2,ia1,ia2
  !============
  real(fp) :: xo2,xo6
  !============
  real(fp) :: x(inx),th(inth),ph(inph)
  real(fp) :: fspl(4,4,4,inf4,inf5,inph),wk(nwk)
  real(fp) :: bcxmin(inb1x,*),bcxmax(inb1x,*) ! inth x inph defined (if used)
  real(fp) :: bcthmin(inb1th,*),bcthmax(inb1th,*) ! inx x inph defined (if used)
  real(fp) :: bcphmin(inb1ph,*),bcphmax(inb1ph,*) ! inx x inth defined (if used)
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
  !
  !   bc data at xmin, xmax  vs.  theta,phi
  !   bc data at thmin, thmax  vs.  x,phi
  !   bc data at phmin, phmax  vs.  x,theta
  !
  !   ibcxmin -- indicator for boundary condition at x(1):
  !    bcxmin(...) -- boundary condition data
  !     =-1 -- use periodic boundary condition
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
  !     (interpretation as with ibcxmin, bcxmin)
  !     NOTE:  if ibcxmin=-1 then the periodic BC applies on both sides
  !            and ibcxmax is ignored.
  !   inb1x -- 1st dimension of bcxmin, bcxmax: if ibcxmin or ibcxmax .gt. 0
  !            this must be .ge. inth:
  !
  !   interpretation of ibcthmin,bcthmin,ibcthmax,bcthmax,inb1th
  !     is same as with ibcxmin,...
  !
  !   interpretation of ibcphmin,bcphmin,ibcphmax,bcphmax,inb1ph
  !     is same as with ibcxmin,...
  !
  !   the explicit bdy condition arrays are referenced only if the
  !     corresponding "ibc" flag values are set to 1 or 2.
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
  !   wk -- must be at least 5*max(inx,inth,inph) large -- or more, see
  !         comments, above.
  !  nwk -- size of workspace
  !
  !---------------------------------
  !  ** in what follows, f is an abbreviation for fspl **
  !
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
  integer iselect1(10)
  integer iselect2(10)
  !
  real(fp) :: z0,z1,ztol
  real(fp) :: zcur(1)
  !
  z0 = 0.0_fp
  z1 = 1.0_fp
  ztol = 1.0E-3_fp
  !
  !---------------------------------
  !
  ier=0
  !
  iflg=0
  !
  !  check phi bdy condition "linearity"
  !
  if(ibcphmin.ne.-1) then
     if((ibcphmin.eq.1).or.(ibcphmin.eq.2)) then
        do ith=1,inth
           do ix=1,inx
              if(bcphmin(ix,ith).ne.z0) iflg=1
           end do
        end do
     end if
     if((ibcphmax.eq.1).or.(ibcphmax.eq.2)) then
        do ith=1,inth
           do ix=1,inx
              if(bcphmax(ix,ith).ne.z0) iflg=1
           end do
        end do
     end if
  end if
  !
  itest=10*max(inx,inth,inph) + 20*inx*inth
  if(iflg.eq.1) then
     itest=16*inx*inth*inph
  end if
  !
  if(nwk.lt.itest) then
     write(6,9901) nwk,itest
     ier=1
9901 format(' ?tcspline:  workspace too small.'/ &
          '  user supplied nwk=',i7,'; need at least: ',i7/ &
          '  If no explicit df/dph boundary condition is set,'/ &
          '  nwk = 20*inx*inth + 10*max(inx,inth,inph) can be used.'/ &
          '  If an explicit df/dph or d2f/dph2 boundary condition'/ &
          '  is set, nwk=16*inx*inth*inph is required.')
  end if
  if(inx.lt.2) then
     write(6,'('' ?tcspline:  at least 2 x points required.'')')
     ier=1
  end if
  if(inth.lt.2) then
     write(6,'('' ?tcspline:  need at least 2 theta points.'')')
     ier=1
  end if
  if(inph.lt.2) then
     write(6,'('' ?tcspline:  need at least 2 phi points.'')')
     ier=1
  end if
  !
  if((ibcxmin.eq.1).or.(ibcxmax.eq.1).or.(ibcxmin.eq.2).or. &
       (ibcxmax.eq.2)) then
     if(inb1x.lt.inth) then
        ier=1
        write(6, &
             '('' ?tcspline:  1st dim of bcxmin/max arrays .lt. inth'')')
     end if
  end if
  !
  if((ibcthmin.eq.1).or.(ibcthmax.eq.1).or.(ibcthmin.eq.2).or. &
       (ibcthmax.eq.2)) then
     if(inb1th.lt.inx) then
        ier=1
        write(6, &
             '('' ?tcspline:  1st dim of bcthmin/max arrays .lt. inx'')')
     end if
  end if
  !
  if((ibcphmin.eq.1).or.(ibcphmax.eq.1).or.(ibcphmin.eq.2).or. &
       (ibcphmax.eq.2)) then
     if(inb1ph.lt.inx) then
        ier=1
        write(6, &
             '('' ?tcspline:  1st dim of bphmin/max arrays .lt. inx'')')
     end if
  end if
  !
  call ibc_ck(ibcxmin,'tcspline','xmin',-1,7,ier)
  if(ibcxmin.ge.0) call ibc_ck(ibcxmax,'tcspline','xmax',0,7,ier)
  !
  call ibc_ck(ibcthmin,'tcspline','thmin',-1,7,ier)
  if(ibcthmin.ge.0) call ibc_ck(ibcthmax,'tcspline','thmax',0,7,ier)
  !
  call ibc_ck(ibcphmin,'tcspline','phmin',-1,7,ier)
  if(ibcphmax.ge.0) call ibc_ck(ibcphmax,'tcspline','phmax',0,7,ier)
  !
  !  check ilinx & x vector
  !
  call splinck(x,inx,ilinx,ztol,ierx)
  if(ierx.ne.0) ier=2
  !
  if(ier.eq.2) then
     write(6,'('' ?tcspline:  x axis not strict ascending'')')
  end if
  !
  !  check ilinth & th vector
  !
  call splinck(th,inth,ilinth,ztol,ierth)
  if(ierth.ne.0) ier=3
  !
  if(ier.eq.3) then
     write(6,'('' ?tcspline:  theta axis not strict ascending'')')
  end if
  !
  !  check ilinth & th vector
  !
  call splinck(ph,inph,ilinph,ztol,ierph)
  if(ierph.ne.0) ier=4
  !
  if(ier.eq.4) then
     write(6,'('' ?tcspline:  phi axis not strict ascending'')')
  end if
  !
  if(ier.ne.0) return
  !
  !------------------------------------
  !
  !  part 1.  compute (x,theta) spline coeffs via an intermediate
  !  routine that call bcspline
  !
  !  workspace addresses
  !
  iaspl2=1
  iabcx1=iaspl2+16*inx*inth
  iabcx2=iabcx1+inth
  iabcth1=iabcx2+inth
  iabcth2=iabcth1+inx
  iawk=iabcth2+inx
  inwk=nwk-iawk+1
  !
  do iph=1,inph
     !
     !  copy bc data
     !
     do ix=1,inx
        wk(iabcth1+ix-1)=0.0_fp
        wk(iabcth2+ix-1)=0.0_fp
        if((ibcthmin.eq.1).or.(ibcthmin.eq.2)) then
           wk(iabcth1+ix-1)=bcthmin(ix,iph)
        end if
        if((ibcthmin.ne.-1).and. &
             ((ibcthmax.eq.1).or.(ibcthmax.eq.2))) then
           wk(iabcth2+ix-1)=bcthmax(ix,iph)
        end if
     end do
     do ith=1,inth
        wk(iabcx1+ith-1)=0.0_fp
        wk(iabcx2+ith-1)=0.0_fp
        if((ibcxmin.eq.1).or.(ibcxmin.eq.2)) then
           wk(iabcx1+ith-1)=bcxmin(ith,iph)
        end if
        if((ibcxmin.ne.-1).and. &
             ((ibcxmax.eq.1).or.(ibcxmax.eq.2))) then
           wk(iabcx2+ith-1)=bcxmax(ith,iph)
        end if
     end do
     !
     !  call 2d spline intermediary routine
     !
     call tcsp23(x,inx,th,inth,fspl(1,1,1,1,1,iph),inf4, &
          ibcxmin,wk(iabcx1),ibcxmax,wk(iabcx2), &
          ibcthmin,wk(iabcth1),ibcthmax,wk(iabcth2), &
          wk(iaspl2),wk(iawk),inwk,ilinx,ilinth,ier)
     !
     if(ier.ne.0) then
        write(6,*) ' ?tcspline:  error in 2d spline, exiting.'
        return
     end if
     !
  end do
  !
  !  ok now fspl(*,*,1,*,*,*) have been evaluated and C2 in (x,theta)
  !  now need to extend to coeffs in phi direction.
  !
  xo2=0.5_fp
  xo6=1.0_fp/6.0_fp
  !
  !  spline each (x,th) coeff in the phi direction
  !
  inpho=4*(inph-1)
  do ith=1,inth-1
     do ix=1,inx-1
        !
        do ic1=1,4
           do ic2=1,4
              !
              !  copy coeff. ordinates in
              !
              do iph=1,inph
                 wk(4*(iph-1)+1)=fspl(ic1,ic2,1,ix,ith,iph)
              end do
              !
              !  use linear BC on this first pass; will correct later if
              !  necessary
              !
              wk(2)=0.0_fp
              wk(3)=0.0_fp
              wk(inpho+2)=0.0_fp
              wk(inpho+3)=0.0_fp
              !
              ibcphmina=ibcphmin
              ibcphmaxa=ibcphmax
              if(iflg.eq.1) then
                 if((ibcphmin.eq.1).or.(ibcphmin.eq.2)) ibcphmina=0
                 if((ibcphmax.eq.1).or.(ibcphmax.eq.2)) ibcphmaxa=0
              end if
              !
              call v_spline(ibcphmina,ibcphmaxa,inph,ph,wk, &
                   wk(4*inph+1))
              !
              !  copy coeffs out
              !
              do iph=1,inph-1
                 fspl(ic1,ic2,2,ix,ith,iph)=wk(4*(iph-1)+2)
                 fspl(ic1,ic2,3,ix,ith,iph)=wk(4*(iph-1)+3)*xo2
                 fspl(ic1,ic2,4,ix,ith,iph)=wk(4*(iph-1)+4)*xo6
              end do
              !
           end do                    ! ic2
        end do                       ! ic1
        !
     end do                          ! ix
  end do                             ! ith
  !
  !  if there are "non-linear" BCs requiring correction...
  !
  !  at each (x(ix),th(ith)) get the d/dph BC's right while preserving C2
  !  everywhere...
  !
  if(iflg.eq.1) then
     !
     !  first get BC correction numbers
     !
     iabcph1=1
     iabcph2=iabcph1+inx*inth
     iaccoef=iabcph2+inx*inth
     iawk=iaccoef+12*inx*inth*inph
     inwk=nwk-iawk+1
     !
     do i=1,10
        iselect1(i)=0
        iselect2(i)=0
     end do
     !
     !  note because iflg=1, we know at least one of ibcphmin/max = 1 or 2
     !
     iskip1=0
     if(ibcphmin.eq.1) then
        iselect1(4)=1               ! df/dph
     else if(ibcphmin.eq.2) then
        iselect1(7)=1               ! d2f/dph2
     else
        iskip1=1
     end if
     !
     iskip2=0
     if(ibcphmax.eq.1) then
        iselect2(4)=1               ! df/dph
     else if(ibcphmax.eq.2) then
        iselect2(7)=1               ! d2f/dph2
     else
        iskip2=1
     end if
     !
     ia1=iabcph1-1
     ia2=iabcph2-1
     do ith=1,inth
        do ix=1,inx
           ia1=ia1+1
           ia2=ia2+1
           !
           if(iskip1.eq.0) then
              call tcspeval(x(ix),th(ith),ph(1),iselect1, zcur, &
                   x,inx,th,inth,ph,inph,ilinx,ilinth,ilinph, &
                   fspl,inf4,inf5,ier)
              if(ier.ne.0) then
                 write(6,*) ' ?? tcspline:  error in tcspeval call'
                 return
              end if
              wk(ia1)=bcphmin(ix,ith)-zcur(1) ! correction needed
           else
              wk(ia1)=z0
           end if
           !
           if(iskip2.eq.0) then
              call tcspeval(x(ix),th(ith),ph(inph),iselect2, zcur, &
                   x,inx,th,inth,ph,inph,ilinx,ilinth,ilinph, &
                   fspl,inf4,inf5,ier)
              if(ier.ne.0) then
                 write(6,*) ' ?? tcspline:  error in tcspeval call'
                 return
              end if
              wk(ia2)=bcphmax(ix,ith)-zcur(1) ! correction needed
           else
              wk(ia2)=z0
           end if
        end do
     end do
     !
     call tcspcorr(x,inx,th,inth,ph,inph,fspl,inf4,inf5, &
          ibcxmin,ibcxmax,ibcthmin,ibcthmax, &
          ibcphmin,wk(iabcph1),ibcphmax,wk(iabcph2), &
          wk(iaccoef),wk(iawk),inwk)
     !
  end if
  !
  return
end subroutine tcspline
!-----------------------------------------
subroutine tcspcorr(x,inx,th,inth,ph,inph,fspl,inf4,inf5, &
     ibcxmin,ibcxmax,ibcthmin,ibcthmax, &
     ibcphmin,bcph1,ibcphmax,bcph2,ccorr,wk,nwk)
  use psp_precision_mod, only: fp
  !
  !  intermediary routine for tcspline:
  !  do correction needed to get C2 3d spline with phi bdy conditions
  !  matched.
  !
  !  all input unless noted:
  !
  !============
  implicit none
  integer inth,inph,inf4,inf5,nwk,inx,iawk2,inpho,ith,ix,iph
  integer inxo,icph,intho,icx,icth
  !============
  real(fp) :: xo2,xo6,z0,zfac
  !============
  real(fp) :: x(inx)                       ! x axis
  real(fp) :: th(inth)                     ! th axis
  real(fp) :: ph(inph)                     ! ph axis
  !
  real(fp) :: fspl(4,4,4,inf4,inf5,inph)   ! spline coeffs -- adjusted
  !
  integer ibcxmin,ibcxmax           ! x BC flags
  integer ibcthmin,ibcthmax         ! th BC flags
  integer ibcphmin,ibcphmax         ! ph BC flags
  !
  real(fp) :: bcph1(inx,inth)              ! ph BC correction array @ phi(1)
  real(fp) :: bcph2(inx,inth)              ! ph BC correction array @ phi(inph)
  !
  !  workspaces:
  !
  real(fp) :: ccorr(3,4,inx,inth,inph)     ! correction coefficients (partial)
  !
  real(fp) :: wk(nwk)                      ! workspace
  !
  !---------------------------
  !
  xo2=0.5_fp
  xo6=1.0_fp/6.0_fp
  !
  !  1.  splines in phi -- fcns are zero everywhere but have non-zero BCs
  !
  if(nwk.lt.10*max(inx,inth,inph)) then
     write(6,*) ' ?? programming error in tcspcorr (tcspline)'
     return
  end if
  !
  z0=0.0_fp
  iawk2=4*inph+1
  !
  inpho=4*(inph-1)
  do ith=1,inth
     do ix=1,inx
        !
        do iph=1,inph
           wk(4*(iph-1)+1)=z0
        end do
        !
        !  set BC for this 1d spline correction
        !
        if(ibcphmin.eq.1) then
           wk(2)=bcph1(ix,ith)
        else if(ibcphmin.eq.2) then
           wk(3)=bcph1(ix,ith)
        end if
        !
        if(ibcphmax.eq.1) then
           wk(inpho+2)=bcph2(ix,ith)
        else if(ibcphmax.eq.2) then
           wk(inpho+3)=bcph2(ix,ith)
        end if
        !
        call v_spline(ibcphmin,ibcphmax,inph,ph,wk,wk(iawk2))
        !
        !  copy non-zero coeffs out to ccorr
        !
        do iph=1,inph-1
           ccorr(1,1,ix,ith,iph)=wk(4*(iph-1)+2)
           ccorr(2,1,ix,ith,iph)=wk(4*(iph-1)+3)*xo2
           ccorr(3,1,ix,ith,iph)=wk(4*(iph-1)+4)*xo6
        end do
        !
     end do
  end do
  !
  !  2. spline the coeffs in x -- use ibcx flags & zero for derivative
  !  bc if necessary
  !
  iawk2=4*inx+1
  !
  inxo=4*(inx-1)
  do iph=1,inph-1
     do ith=1,inth
        !
        do icph=1,3
           !
           do ix=1,inx
              wk(4*(ix-1)+1)=ccorr(icph,1,ix,ith,iph)
           end do
           !
           !  zero BC:  correction spline
           !
           wk(2)=0.0_fp
           wk(3)=0.0_fp
           wk(inxo+2)=0.0_fp
           wk(inxo+3)=0.0_fp
           !
           call v_spline(ibcxmin,ibcxmax,inx,x,wk,wk(iawk2))
           !
           do ix=1,inx-1
              ccorr(icph,2,ix,ith,iph)=wk(4*(ix-1)+2)
              ccorr(icph,3,ix,ith,iph)=wk(4*(ix-1)+3)*xo2
              ccorr(icph,4,ix,ith,iph)=wk(4*(ix-1)+4)*xo6
           end do
           !
        end do
        !
     end do
  end do
  !
  !  3.  spline all the ccorr coefs in th -- use ibcth flags & zero for
  !      derivative correction BC if necessary
  !
  !      add the results into fspl
  !
  iawk2=4*inth+1
  !
  intho=4*(inth-1)
  do iph=1,inph-1
     do ix=1,inx-1
        !
        do icx=1,4
           do icph=1,3
              !
              do ith=1,inth
                 wk(4*(ith-1)+1)=ccorr(icph,icx,ix,ith,iph)
              end do
              !
              !  zero BC:  correction spline
              !
              wk(2)=0.0_fp
              wk(3)=0.0_fp
              wk(intho+2)=0.0_fp
              wk(intho+3)=0.0_fp
              !
              call v_spline(ibcthmin,ibcthmax,inth,th,wk, &
                   wk(iawk2))
              !
              do ith=1,inth-1
                 do icth=1,4
                    zfac=1.0_fp
                    if(icth.eq.3) zfac=xo2
                    if(icth.eq.4) zfac=xo6
                    fspl(icx,icth,icph+1,ix,ith,iph)= &
                         fspl(icx,icth,icph+1,ix,ith,iph)+ &
                         wk(4*(ith-1)+icth)*zfac
                 end do
              end do
              !
           end do
        end do
        !
     end do
  end do
  !
  return
end subroutine tcspcorr
!-----------------------------------------
subroutine tcsp23(x,inx,th,inth,fspl,inf4, &
     ibcxmin,bcxmin,ibcxmax,bcxmax, &
     ibcthmin,bcthmin,ibcthmax,bcthmax, &
     fspl2,wk,nwk,ilinx,ilinth,ier)
  use psp_precision_mod, only: fp
  !
  !  intermediary routines
  !  call bcspline from tcspline loop
  !  to set up 2d splines in each phi plane
  !
  !============
  implicit none
  integer inth,inf4,ibcxmin,ibcxmax,ibcthmin,ibcthmax,nwk,ilinx
  integer ilinth,ier,inx,ith,ix,j,i
  !============
  real(fp) :: x(inx)                       ! x axis
  real(fp) :: th(inth)                     ! th axis
  !
  real(fp) :: fspl(4,4,4,inf4,inth)        ! fspl array at one phi pt.
  real(fp) :: fspl2(4,4,inx,inth)          ! temp fspl array for bcspline
  !
  real(fp) :: bcxmin(inth),bcxmax(inth)    ! d/dx BC's @ x(1),x(inx), th(*)
  real(fp) :: bcthmin(inx),bcthmax(inx)    ! d/dth BC's @ th(1),th(inth), x(*)
  !
  real(fp) :: wk(nwk)
  !
  !--------------------
  !
  !  1.  copy spline data in
  !
  do ith=1,inth
     do ix=1,inx
        fspl2(1,1,ix,ith)=fspl(1,1,1,ix,ith)
     end do
  end do
  !
  !  2.  compute the 2d spline
  !
  call bcspline(x,inx,th,inth,fspl2,inx, &
       ibcxmin,bcxmin,ibcxmax,bcxmax, &
       ibcthmin,bcthmin,ibcthmax,bcthmax, &
       wk,nwk,ilinx,ilinth,ier)
  if(ier.ne.0) return
  !
  !  3.  copy spline coeff results out
  !
  do ith=1,inth-1
     do ix=1,inx-1
        do j=1,4
           do i=1,4
              fspl(i,j,1,ix,ith)=fspl2(i,j,ix,ith)
           end do
        end do
     end do
  end do
  !
  return
end subroutine tcsp23
