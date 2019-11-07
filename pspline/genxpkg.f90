subroutine genxpkg(nx,x,xpkg,iper,imsg,itol,ztol,ialg,ier)
  use psp_precision_mod, only: fp
  !
  !  from an x axis assemble a "package":
  !
  !  MOD DMC Feb 2010: handle very small grids: nx=2, nx=3...
  !  interchanged meaning of xpkg(1,4) and xpkg(3,4);
  !  xpkg(3,4) only set if nx.ge.3;
  !  xpkg(4,4) only set if nx.ge.4.
  !
  !  there are corresponding changes in xlookup: simplified code lookup
  !  code for the cases nx=2 and nx=3.
  !
  !     xpkg(j,1) = x(j), j = 1 to nx    ... nx.ge.2
  !
  !     if(abs(ialg).ne.3) then...
  !       for j=1:nx-1
  !       xpkg(j,2) = h(j) = x(j+1)-x(j), j=1 to nx-1
  !     else
  !       for j=1:nx-1
  !       xpkg(j,2) = index location, with linear offset, of
  !                   xpkg(1,1)+<h>*(j-1) in the original x(1:nx)
  !            (piecewise linear indexing function)
  !     end if
  !     xpkg(nx,2) = <h> = (x(nx)-x(1))/(nx-1)
  !
  !     xpkg(j,3) = 1/h(j)
  !     xpkg(nx,3) = 1/<h>
  !
  !     xpkg(1,4) = +/- tolerance epsilon for out-of-range warnings/errors
  !                 +if message is to be written on out of range errors
  !                 -if no message to be written.  In either case, a
  !                        warning flag is set.
  !
  !     xpkg(2,4) = 1.0 if x is *periodic*, else 0.0
  !
  !  only set if nx.ge.3:
  !     xpkg(3,4) = 0.0 if x is *evenly spaced* (to ztol or 5.e-7 rel) else:
  !                   = 1.0 => use (1/h) newton's method like search algorithm
  !                   = 2.0 => use binary search algorithm
  !                   = 3.0 => use piecewise-linear indexing function...
  !
  !  only set if nx.ge.4:
  !     xpkg(4,4) = 0.0 -- do not use past search result as start point
  !                        for next search;
  !               = 1.0 -- do use past search result as start point for
  !                        next search.
  !
  !  tolerance epsilon means:
  !     if xtest is within epsilon*max(abs(x(1)),abs(x(nx))) outside the
  !     range [x(1),x(nx)] treat xtest as being at the endpoint x(1) or
  !     x(nx) whichever is closer.
  !
  ! input arguments:
  !
  !============
  implicit none
  integer ialgu,iabs,ix,ixp,itest,i
  !============
  real(fp) :: ztolr,ztola,zh,xtest
  !============
  integer nx                        ! length of x, .ge.4
  real(fp) :: x(nx)                        ! x axis vector, strict ascending
  !
  integer iper                      ! =1 if x is periodic
  integer imsg                      ! =1 for range error messages
  integer itol                      ! =1 to specify tolerance, else default
  !
  ! default tolerance is 5.0e-7
  !
  real(fp) :: ztol                         ! range tolerance, if itol=1
  !
  ! lookup algorithm selection for uneven grids:
  !
  integer ialg                      ! = +/- 1:  <1/h(j)> "Newton" method
  !                                       ! = +/- 2:  binary search
  !                                       ! = +/- 3:  indexing function
  !
  !       to use past search result as init. cond. for next search, give
  !       ialg .lt. 0; otherwise give ialg .gt. 0.
  !
  ! output arguments:
  !
  real(fp) :: xpkg(nx,4)                   ! xpkg, as described above
  integer ier                       ! completion code, 0=OK
  !
  !------------------------------------------------
  !
  if(nx.lt.2) then
     write(6,*) ' %genxpkg:  nx.ge.2 required!'
     ier=1
     return
  else
     ier=0
  end if
  !
  ialgu=ialg
  if(ialgu.eq.0) ialgu=3
  if(iabs(ialgu).gt.3) ialgu=3
  !
  !  get tolerance for end point range check & even spacing check
  !
  if(itol.eq.1) then
     ztolr=abs(ztol)
  else
     ztolr=5.0E-7_fp
  end if
  !
  ztola=max(abs(x(1)),abs(x(nx)))*ztolr
  !
  !  assume even spacing for now...
  !
  if(nx.ge.3) then
     xpkg(3,4)=0.0_fp
  end if
  !
  !  mark if x axis is a periodic coordinate
  !
  if(iper.eq.1) then
     xpkg(2,4)=1.0_fp
  else
     xpkg(2,4)=0.0_fp
  end if
  !
  !  store tolerance parameter
  !
  xpkg(1,4)=ztola
  !
  !  mark if messages are to be written if range lookup errors occur
  !
  if(imsg.eq.1) then
     continue                       ! xpkg(1,4) left .gt. 0.0
  else
     xpkg(1,4)=-xpkg(1,4)
  end if
  !
  !  OK check linearity and spacing
  !
  ier=0
  !
  xpkg(nx,2)=(x(nx)-x(1))/(nx-1)    ! average spacing
  !
  do ix=1,nx
     xpkg(ix,1)=x(ix)
     if((ier.eq.0).and.(ix.lt.nx)) then
        if(x(ix+1).le.x(ix)) then
           ier=1
           write(6,*) ' %genxpkg:  x axis not strict ascending!'
        else
           zh=x(ix+1)-x(ix)
           !
           !  checking even spacing now...
           !
           if(nx.ge.3) then
              if(abs(zh-xpkg(nx,2)).gt.ztola) xpkg(3,4)=1.0_fp
           end if
           !
           xpkg(ix,2)=zh
           xpkg(ix,3)=1.0_fp/zh
           !
        end if
     end if
  end do
  !
  if(ier.ne.0) return
  !
  !  OK, store inverse average spacing
  !
  xpkg(nx,3)=1.0_fp/xpkg(nx,2)
  !
  !  if even spacing is detected, redefine x axis slightly, for
  !  improved regularity of behaviour
  !
  if(nx.ge.3) then
     !
     if(xpkg(3,4).eq.0.0_fp) then
        do ix=1,nx-2
           ixp=ix+1
           xpkg(ixp,1)=xpkg(1,1)+ix*xpkg(nx,2)
        end do
        if(nx.gt.3) then
           if(ialgu.lt.0) then
              xpkg(4,4)=1.0_fp ! check init guess
           else
              xpkg(4,4)=0.0_fp
           end if
        end if
     end if
     !
     !  if uneven spacing is detected, must use an algorithm...
     !
     if(xpkg(3,4).ne.0.0_fp) then
        xpkg(3,4)=abs(ialgu)
        if(nx.gt.3) then
           if(ialgu.lt.0) then
              xpkg(4,4)=1.0_fp ! check init guess
           else
              xpkg(4,4)=0.0_fp
           end if
        end if
        !
        if(abs(ialgu).eq.3) then
           !
           !  construct a piecewise linear indexing function
           !
           xpkg(1,2)=1.0_fp
           xtest=xpkg(1,1)
           itest=1
           do i=2,nx-1
              xtest=xtest+xpkg(nx,2) ! x1 + (i-1)*<h>
10            continue
              if((xpkg(itest,1).le.xtest).and. &
                   (xtest.le.xpkg(itest+1,1))) then
                 xpkg(i,2)=itest+(xtest-xpkg(itest,1))/ &
                      (xpkg(itest+1,1)-xpkg(itest,1))
              else
                 itest=itest+1
                 go to 10
              end if
           end do
           !
           !  (implicitly, xpkg(nx,2)=nx*1.0; but we leave <h> in xpkg(nx,2)
           !   instead...)
           !
        end if
     end if
     !
  end if  ! nx.ge.3
  !
  return
end subroutine genxpkg
