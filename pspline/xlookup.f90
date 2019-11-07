subroutine xlookup(ivec,xvec,nx,xpkg,imode,iv,dxn,hv,hiv,iwarn)

  !  vector lookup routine
  !
  !   given a set of x points xvec(...) and an x grid (xpkg, nx x pts)
  !   return the vector of indices iv(...) and displacements dxv(...)
  !   within the indexed zones, corresponding to each x point.
  !
  !   if any of the x points in the vector are out of range the warning
  !   flag iwarn is set.
  !
  !  MOD DMC Feb 2010: changes related to supporting nx=2 and nx=3 small grids.
  !  Changes are consistent with Feb 2010 changes in genxpkg:
  !    meanings of xpkg(1,4) and xpkg(3,4) interchanged.
  !
  !    if nx.eq.2:  xpkg(3,4) and xpkg(4,4) never referenced;
  !    if nx.eq.3:  xpkg(4,4) never referenced.
  !
  !  input:
  use psp_precision_mod, only: fp
  implicit none
  integer inum,istat,ilin,ialg,iper,imsg,init_guess,iprev,i
  integer init,iprob,isrch,i_sign,inc
  real(kind=fp) :: ztola,period,hav,havi,hloci,hloc,zdelta,xfac
  REAL(kind=fp) :: zindx0,zdindx,zindex
  !============
  integer ivec                      ! size of vector
  real(kind=fp) :: xvec(ivec)       ! x points to lookup on xpkg grid
  integer nx                        ! size of grid
  real(kind=fp) :: xpkg(nx,4)       ! grid data
  integer imode                     ! output control flag
  !  imode=1:  return indices iv(...) and un-normalized displacements dxn(...)
  !            ignore hv and hiv
  !                 dxn(j)=xvec(j)-xpkg(iv(j),1)
  !
  !  imode=2   return indices iv(...) and *normalized* displacements dxn(...)
  !            and hv(...) and hiv(...)
  !                 dxn(j)=(xvec(j)-xpkg(iv(j),1))*hiv(j)
  !
  !  output:
  integer iv(ivec)                  ! index into grid for each xvec(j)
  !  note: old values of iv(...) may be used as start point for grid
  !  searches, depending on xpkg controls
  real(fp) :: dxn(ivec)        ! displacement w/in zone (see imode)
  real(fp) :: hv(ivec)         ! zone width (if imode=2)
  real(fp) :: hiv(ivec)        ! inverse zone width (if imode=2)
  integer iwarn                     ! =0: OK; =n:  n points out of range

  !  xpkg is a "structure" constructed by subroutine genxpkg, which
  !  contains the x grid, xpkg(1:nx,1), spacing and error handling
  !  information -- see genxpkg.f90

  real(fp), dimension(:), allocatable :: xuse
  logical, dimension(:), allocatable :: iok
  integer, dimension(:), allocatable :: imina,imaxa

  if(nx.lt.2) then
     iwarn=1
     write(6,*) ' ?? xlookup: nx.lt.2, nx=',nx
     go to 1100
  end if

  inum=ivec
  allocate(xuse(inum),stat=istat)
  if(istat.ne.0) then
     iwarn=1
     write(6,*) ' ?? xlookup "xuse" vector allocation failure!'
     go to 1000
  end if

  if(nx.eq.2) then
     ilin=1
  end if

  if(nx.gt.2) then
     if(xpkg(3,4).eq.0.0_fp) then
        ilin=1              ! evenly spaced grid
     else
        ilin=0
        if(xpkg(3,4).gt.2.5_fp) then
           ialg=3
        else if(xpkg(3,4).gt.1.5_fp) then
           ialg=2
        else
           ialg=1
        end if
     end if
  end if

  if(xpkg(2,4).ne.0.0_fp) then
     iper=1                         ! periodic grid
  else
     iper=0
  end if

  ztola=abs(xpkg(1,4))              ! tolerance for range checking

  if(xpkg(1,4).ge.0.0_fp) then
     imsg=1                         ! write message on range error
  else
     imsg=0
  end if

  init_guess=0
  if(nx.gt.3) then
     if(xpkg(4,4).gt.0.0_fp) then
        init_guess=1
        iprev=min((nx-1),max(1,iv(ivec)))
     else
        init_guess=0
     end if
  end if

  iwarn=0

  !---------------------
  !  range check
  !
  if(iper.eq.0) then
     !
     !  check min/max with tolerance
     !
     do i=1,ivec
        if(xvec(i).lt.xpkg(1,1)) then
           xuse(i)=xpkg(1,1)
           if((xpkg(1,1)-xvec(i)).gt.ztola) iwarn=iwarn+1
        else if(xvec(i).gt.xpkg(nx,1)) then
           xuse(i)=xpkg(nx,1)
           if((xvec(i)-xpkg(nx,1)).gt.ztola) iwarn=iwarn+1
        else
           xuse(i)=xvec(i)
        end if
     end do

  else

     ! normalize to interval
     period=xpkg(nx,1)-xpkg(1,1)
     do i=1,ivec
        if((xvec(i).lt.xpkg(1,1)).or.(xvec(i).gt.xpkg(nx,1))) then
           xuse(i)=mod(xvec(i)-xpkg(1,1),period)
           if(xuse(i).lt.0.0_fp) xuse(i)=xuse(i)+period
           xuse(i)=xuse(i)+xpkg(1,1)
           xuse(i)=max(xpkg(1,1),min(xpkg(nx,1),xuse(i)))
        else
           xuse(i)=xvec(i)
        end if
     end do

  end if

  if((imsg.eq.1).and.(iwarn.gt.0)) then
     write(6,*) ' %xlookup:  ',iwarn,' points not in range: ', &
          xpkg(1,1),' to ',xpkg(nx,1)
  end if

  !---------------------
  !  actual lookup -- initially, assume even spacing
  !
  !   initial index guess:  1 + <1/h>*(x-x0) ... then refine
  !
  hav=xpkg(nx,2)
  havi=xpkg(nx,3)
  if(ilin.eq.1) then
     !
     !  faster lookup OK:  even spacing
     !
     if(init_guess.eq.0) then
        !
        !  even spacing lookup, no initial guess from previous iteration
        !
        do i=1,ivec
           iv(i)=1+havi*(xuse(i)-xpkg(1,1))
           iv(i)=max(1,min((nx-1),iv(i)))
           if(imode.eq.1) then
              dxn(i)=(xuse(i)-xpkg(iv(i),1))
           else
              dxn(i)=(xuse(i)-xpkg(iv(i),1))*havi
              hiv(i)=havi
              hv(i)=hav
           end if
        end do

     else
        !
        !  even spacing lookup, do use initial guess from previous iteration
        !
        do i=1,ivec
           if((xpkg(iprev,1).le.xuse(i)).and. &
                (xuse(i).le.xpkg(iprev+1,1))) then
              iv(i)=iprev
           else
              iv(i)=1+havi*(xuse(i)-xpkg(1,1))
              iv(i)=max(1,min((nx-1),iv(i)))
              iprev=iv(i)
           end if
           if(imode.eq.1) then
              dxn(i)=(xuse(i)-xpkg(iv(i),1))
           else
              dxn(i)=(xuse(i)-xpkg(iv(i),1))*havi
              hiv(i)=havi
              hv(i)=hav
           end if
        end do
     end if
     go to 1000
  end if

4 continue

  init=-1
  allocate(iok(inum),stat=istat)
  if(istat.ne.0) then
     iwarn=1
     write(6,*) ' ?? xlookup "iok" vector allocation failure!'
     go to 1000
  end if

  if(ialg.lt.3) then
     allocate(imina(inum),stat=istat)
     if(istat.ne.0) then
        iwarn=1
        write(6,*) ' ?? xlookup "imina" vector allocation failure!'
        go to 1000
     end if

     allocate(imaxa(inum),stat=istat)
     if(istat.ne.0) then
        iwarn=1
        write(6,*) ' ?? xlookup "imaxa" vector allocation failure!'
        go to 1000
     end if
  end if

5 continue                          ! re-entry -- hit problem cases...

  iprob=0                           ! count "problem" cases...
  init=init+1

  if(init_guess.eq.0) then
     go to (100,200,300) ialg
  else
     go to (150,250,350) ialg
  end if
  !-----------------------------------------------------------
  !  Newton like algorithm:  use local spacing to estimate
  !   step to next zone guess; don't use prior guess
  !
100 continue
  do i=1,ivec
     !
     !  first iteration
     !         
     if(init.eq.0) then
        iok(i)=.FALSE.
        imina(i)=1
        imaxa(i)=nx-1
        iv(i)=1+havi*(xuse(i)-xpkg(1,1))
        iv(i)=max(imina(i),min(imaxa(i),iv(i)))
        if(xuse(i).lt.xpkg(iv(i),1)) then
           isrch=0
           i_sign=-1
           imaxa(i)=max(1,(iv(i)-1))
        else if(xuse(i).gt.xpkg(iv(i)+1,1)) then
           isrch=1
           i_sign=+1
           imina(i)=min((nx-1),(iv(i)+1))
        else
           iok(i)=.TRUE.
        end if
     end if
     !
     !  second iteration
     !
     if(.not.iok(i)) then
        hloci=xpkg(iv(i),3)
        hloc=xpkg(iv(i),2)
        zdelta=(xuse(i)-xpkg(iv(i)+isrch,1))
        if(i_sign*zdelta.le.hloc) then
           inc=i_sign
        else
           inc=zdelta*hloci
           inc=inc+i_sign
        end if
        iv(i)=iv(i)+inc
        iv(i)=max(imina(i),min(imaxa(i),iv(i)))
        if(xuse(i).lt.xpkg(iv(i),1)) then
           isrch=0
           i_sign=-1
           imaxa(i)=max(1,(iv(i)-1))
        else if(xuse(i).gt.xpkg(iv(i)+1,1)) then
           isrch=1
           i_sign=+1
           imina(i)=min((nx-1),(iv(i)+1))
        else
           iok(i)=.TRUE.
        end if
        !
        !  third iteration
        !
        if(.not.iok(i)) then
           hloci=xpkg(iv(i),3)
           hloc=xpkg(iv(i),2)
           zdelta=(xuse(i)-xpkg(iv(i)+isrch,1))
           if(i_sign*zdelta.le.hloc) then
              inc=i_sign
           else
              inc=zdelta*hloci
              inc=inc+i_sign
           end if
           iv(i)=iv(i)+inc
           iv(i)=max(imina(i),min(imaxa(i),iv(i)))
           if(xuse(i).lt.xpkg(iv(i),1)) then
              isrch=0
              i_sign=-1
              imaxa(i)=max(1,(iv(i)-1))
           else if(xuse(i).gt.xpkg(iv(i)+1,1)) then
              isrch=1
              i_sign=+1
              imina(i)=min((nx-1),(iv(i)+1))
           else
              iok(i)=.TRUE.
           end if
           !
           !  fourth iteration
           !
           if(.not.iok(i)) then
              hloci=xpkg(iv(i),3)
              hloc=xpkg(iv(i),2)
              zdelta=(xuse(i)-xpkg(iv(i)+isrch,1))
              if(i_sign*zdelta.le.hloc) then
                 inc=i_sign
              else
                 inc=zdelta*hloci
                 inc=inc+i_sign
              end if
              iv(i)=iv(i)+inc
              iv(i)=max(imina(i),min(imaxa(i),iv(i)))
              if(xuse(i).lt.xpkg(iv(i),1)) then
                 isrch=0
                 i_sign=-1
                 imaxa(i)=max(1,(iv(i)-1))
              else if(xuse(i).gt.xpkg(iv(i)+1,1)) then
                 isrch=1
                 i_sign=+1
                 imina(i)=min((nx-1),(iv(i)+1))
              else
                 iok(i)=.TRUE.
              end if
              !
              !  fifth iteration
              !
              if(.not.iok(i)) then
                 hloci=xpkg(iv(i),3)
                 hloc=xpkg(iv(i),2)
                 zdelta=(xuse(i)-xpkg(iv(i)+isrch,1))
                 if(i_sign*zdelta.le.hloc) then
                    inc=i_sign
                 else
                    inc=zdelta*hloci
                    inc=inc+i_sign
                 end if
                 iv(i)=iv(i)+inc
                 iv(i)=max(imina(i),min(imaxa(i),iv(i)))
                 if(xuse(i).lt.xpkg(iv(i),1)) then
                    isrch=0
                    i_sign=-1
                    imaxa(i)=max(1,(iv(i)-1))
                 else if(xuse(i).gt.xpkg(iv(i)+1,1)) then
                    isrch=1
                    i_sign=+1
                    imina(i)=min((nx-1),(iv(i)+1))
                 else
                    iok(i)=.TRUE.
                 end if
                 if(.not.iok(i)) iprob=iprob+1
                 !
                 !  end chain of iteration if-then-else blocks
                 !
              end if
           end if
        end if
     end if
     !
     !  end of loop
     !
  end do

  go to 500

  !-----------------------------------------------------------
  !  Newton like algorithm:  use local spacing to estimate
  !   step to next zone guess; DO use prior guess
  !
150 continue

  do i=1,ivec
     !
     !  first iteration
     !
     if(init.eq.0) then
        if((xpkg(iprev,1).le.xuse(i)).and. &
             (xuse(i).le.xpkg(iprev+1,1))) then
           iok(i)=.TRUE.
           iv(i)=iprev
        else
           iok(i)=.FALSE.
           imina(i)=1
           imaxa(i)=nx-1
           iv(i)=1+havi*(xuse(i)-xpkg(1,1))
           iv(i)=max(imina(i),min(imaxa(i),iv(i)))
           if(xuse(i).lt.xpkg(iv(i),1)) then
              isrch=0
              i_sign=-1
              imaxa(i)=max(1,(iv(i)-1))
           else if(xuse(i).gt.xpkg(iv(i)+1,1)) then
              isrch=1
              i_sign=+1
              imina(i)=min((nx-1),(iv(i)+1))
           else
              iok(i)=.TRUE.
           end if
        end if
     end if
     !
     !  second iteration
     !
     if(.not.iok(i)) then
        hloci=xpkg(iv(i),3)
        hloc=xpkg(iv(i),2)
        zdelta=(xuse(i)-xpkg(iv(i)+isrch,1))
        if(i_sign*zdelta.le.hloc) then
           inc=i_sign
        else
           inc=zdelta*hloci
           inc=inc+i_sign
        end if
        iv(i)=iv(i)+inc
        iv(i)=max(imina(i),min(imaxa(i),iv(i)))
        if(xuse(i).lt.xpkg(iv(i),1)) then
           isrch=0
           i_sign=-1
           imaxa(i)=max(1,(iv(i)-1))
        else if(xuse(i).gt.xpkg(iv(i)+1,1)) then
           isrch=1
           i_sign=+1
           imina(i)=min((nx-1),(iv(i)+1))
        else
           iok(i)=.TRUE.
        end if
        !
        !  third iteration
        !
        if(.not.iok(i)) then
           hloci=xpkg(iv(i),3)
           hloc=xpkg(iv(i),2)
           zdelta=(xuse(i)-xpkg(iv(i)+isrch,1))
           if(i_sign*zdelta.le.hloc) then
              inc=i_sign
           else
              inc=zdelta*hloci
              inc=inc+i_sign
           end if
           iv(i)=iv(i)+inc
           iv(i)=max(imina(i),min(imaxa(i),iv(i)))
           if(xuse(i).lt.xpkg(iv(i),1)) then
              isrch=0
              i_sign=-1
              imaxa(i)=max(1,(iv(i)-1))
           else if(xuse(i).gt.xpkg(iv(i)+1,1)) then
              isrch=1
              i_sign=+1
              imina(i)=min((nx-1),(iv(i)+1))
           else
              iok(i)=.TRUE.
           end if
           !
           !  fourth iteration
           !
           if(.not.iok(i)) then
              hloci=xpkg(iv(i),3)
              hloc=xpkg(iv(i),2)
              zdelta=(xuse(i)-xpkg(iv(i)+isrch,1))
              if(i_sign*zdelta.le.hloc) then
                 inc=i_sign
              else
                 inc=zdelta*hloci
                 inc=inc+i_sign
              end if
              iv(i)=iv(i)+inc
              iv(i)=max(imina(i),min(imaxa(i),iv(i)))
              if(xuse(i).lt.xpkg(iv(i),1)) then
                 isrch=0
                 i_sign=-1
                 imaxa(i)=max(1,(iv(i)-1))
              else if(xuse(i).gt.xpkg(iv(i)+1,1)) then
                 isrch=1
                 i_sign=+1
                 imina(i)=min((nx-1),(iv(i)+1))
              else
                 iok(i)=.TRUE.
              end if
              !
              !  fifth iteration
              !
              if(.not.iok(i)) then
                 hloci=xpkg(iv(i),3)
                 hloc=xpkg(iv(i),2)
                 zdelta=(xuse(i)-xpkg(iv(i)+isrch,1))
                 if(i_sign*zdelta.le.hloc) then
                    inc=i_sign
                 else
                    inc=zdelta*hloci
                    inc=inc+i_sign
                 end if
                 iv(i)=iv(i)+inc
                 iv(i)=max(imina(i),min(imaxa(i),iv(i)))
                 if(xuse(i).lt.xpkg(iv(i),1)) then
                    isrch=0
                    i_sign=-1
                    imaxa(i)=max(1,(iv(i)-1))
                 else if(xuse(i).gt.xpkg(iv(i)+1,1)) then
                    isrch=1
                    i_sign=+1
                    imina(i)=min((nx-1),(iv(i)+1))
                 else
                    iok(i)=.TRUE.
                 end if
                 if(.not.iok(i)) iprob=iprob+1
                 !
                 !  end chain of iteration if-then-else blocks
                 !
              end if
           end if
        end if
     end if
     !
     !  end of loop
     !
     iprev=iv(i)
  end do

  go to 500
  !-----------------------------------------------------------
  !  Binary search algorithm
  !
200 continue
  do i=1,ivec
     !
     !  first iteration
     !
     if(init.eq.0) then
        iok(i)=.FALSE.
        imina(i)=1
        imaxa(i)=nx-1
        iv(i)=nx/2
        if(xuse(i).lt.xpkg(iv(i),1)) then
           imaxa(i)=max(1,(iv(i)-1))
        else if(xuse(i).gt.xpkg(iv(i)+1,1)) then
           imina(i)=min((nx-1),(iv(i)+1))
        else
           iok(i)=.TRUE.
        end if
      end if
      !
      !  second iteration
      !
      if(.not.iok(i)) then
        iv(i)=(imina(i)+imaxa(i))/2
        if(xuse(i).lt.xpkg(iv(i),1)) then
          imaxa(i)=max(1,(iv(i)-1))
        else if(xuse(i).gt.xpkg(iv(i)+1,1)) then
          imina(i)=min((nx-1),(iv(i)+1))
        else
          iok(i)=.TRUE.
        end if
        !
        !  third iteration
        !
        if(.not.iok(i)) then
           iv(i)=(imina(i)+imaxa(i))/2
           if(xuse(i).lt.xpkg(iv(i),1)) then
              imaxa(i)=max(1,(iv(i)-1))
           else if(xuse(i).gt.xpkg(iv(i)+1,1)) then
              imina(i)=min((nx-1),(iv(i)+1))
           else
              iok(i)=.TRUE.
           end if
           !
           !  fourth iteration
           !
           if(.not.iok(i)) then
              iv(i)=(imina(i)+imaxa(i))/2
              if(xuse(i).lt.xpkg(iv(i),1)) then
                 imaxa(i)=max(1,(iv(i)-1))
              else if(xuse(i).gt.xpkg(iv(i)+1,1)) then
                 imina(i)=min((nx-1),(iv(i)+1))
              else
                 iok(i)=.TRUE.
              end if
              !
              !  fifth iteration
              !
              if(.not.iok(i)) then
                 iv(i)=(imina(i)+imaxa(i))/2
                 if(xuse(i).lt.xpkg(iv(i),1)) then
                    imaxa(i)=max(1,(iv(i)-1))
                 else if(xuse(i).gt.xpkg(iv(i)+1,1)) then
                    imina(i)=min((nx-1),(iv(i)+1))
                 else
                    iok(i)=.TRUE.
                 end if
                 if(.not.iok(i)) iprob=iprob+1
                 !
                 !  end chain of iteration if-then-else blocks
                 !
              end if
           end if
        end if
     end if
     !
     !  end of loop
     !
  end do

  go to 500
  !-----------------------------------------------------------
  !  Binary search algorithm
  !
250 continue
  do i=1,ivec
     !
     !  first iteration
     !
     if(init.eq.0) then
        if((xpkg(iprev,1).le.xuse(i)).and. &
             (xuse(i).le.xpkg(iprev+1,1))) then
           iok(i)=.TRUE.
           iv(i)=iprev
        else
           iok(i)=.FALSE.
           imina(i)=1
           imaxa(i)=nx-1
           iv(i)=nx/2
           if(xuse(i).lt.xpkg(iv(i),1)) then
              imaxa(i)=max(1,(iv(i)-1))
           else if(xuse(i).gt.xpkg(iv(i)+1,1)) then
              imina(i)=min((nx-1),(iv(i)+1))
           else
              iok(i)=.TRUE.
           end if
        end if
     end if
     !
     !  second iteration
     !
     if(.not.iok(i)) then
        iv(i)=(imina(i)+imaxa(i))/2
        if(xuse(i).lt.xpkg(iv(i),1)) then
           imaxa(i)=max(1,(iv(i)-1))
        else if(xuse(i).gt.xpkg(iv(i)+1,1)) then
           imina(i)=min((nx-1),(iv(i)+1))
        else
           iok(i)=.TRUE.
        end if
        !
        !  third iteration
        !
        if(.not.iok(i)) then
           iv(i)=(imina(i)+imaxa(i))/2
           if(xuse(i).lt.xpkg(iv(i),1)) then
              imaxa(i)=max(1,(iv(i)-1))
           else if(xuse(i).gt.xpkg(iv(i)+1,1)) then
              imina(i)=min((nx-1),(iv(i)+1))
           else
              iok(i)=.TRUE.
           end if
           !
           !  fourth iteration
           !
           if(.not.iok(i)) then
              iv(i)=(imina(i)+imaxa(i))/2
              if(xuse(i).lt.xpkg(iv(i),1)) then
                 imaxa(i)=max(1,(iv(i)-1))
              else if(xuse(i).gt.xpkg(iv(i)+1,1)) then
                 imina(i)=min((nx-1),(iv(i)+1))
              else
                 iok(i)=.TRUE.
              end if
              !
              !  fifth iteration
              !
              if(.not.iok(i)) then
                 iv(i)=(imina(i)+imaxa(i))/2
                 if(xuse(i).lt.xpkg(iv(i),1)) then
                    imaxa(i)=max(1,(iv(i)-1))
                 else if(xuse(i).gt.xpkg(iv(i)+1,1)) then
                    imina(i)=min((nx-1),(iv(i)+1))
                 else
                    iok(i)=.TRUE.
                 end if
                 if(.not.iok(i)) iprob=iprob+1
                 !
                 !  end chain of iteration if-then-else blocks
                 !
              end if
           end if
        end if
     end if
     !
     !  end of loop
     !
     iprev=iv(i)
  end do
  !
  go to 500
  !-----------------------------------------------------------
  !  algorithm:  piecewise linear indexing function lookup & correction
  !
300 continue
  do i=1,ivec
     !
     !  first iteration
     !
     if(init.eq.0) then
        iok(i)=.FALSE.
        !
        !  piecewise linear indexing function on even spaced grid
        !  (not same grid as x axis itself)
        !
        iv(i)=1+havi*(xuse(i)-xpkg(1,1))
        iv(i)=max(1,min(nx-1,iv(i)))
        xfac=(xuse(i)-(xpkg(1,1)+(iv(i)-1)*hav))*havi
        zindx0=xpkg(iv(i),2)
        if(iv(i).lt.nx-1) then
           zdindx=xpkg(iv(i)+1,2)-zindx0
        else
           zdindx=nx-zindx0
        end if
        zindex=zindx0+xfac*zdindx
        iv(i)=zindex
        iv(i)=max(1,min(nx-1,iv(i)))
        if(xuse(i).lt.xpkg(iv(i),1)) then
           i_sign=-1
        else if(xuse(i).gt.xpkg(iv(i)+1,1)) then
           i_sign=+1
        else
           iok(i)=.TRUE.
        end if
     end if
     !
     !  second iteration
     !
     if(.not.iok(i)) then
        iv(i)=iv(i)+i_sign
        if(xuse(i).lt.xpkg(iv(i),1)) then
           i_sign=-1
        else if(xuse(i).gt.xpkg(iv(i)+1,1)) then
           i_sign=+1
        else
           iok(i)=.TRUE.
        end if
        !
        !  third iteration
        !
        if(.not.iok(i)) then
           iv(i)=iv(i)+i_sign
           if(xuse(i).lt.xpkg(iv(i),1)) then
              i_sign=-1
           else if(xuse(i).gt.xpkg(iv(i)+1,1)) then
              i_sign=+1
           else
              iok(i)=.TRUE.
           end if
           !
           !  fourth iteration
           !
           if(.not.iok(i)) then
              iv(i)=iv(i)+i_sign
              if(xuse(i).lt.xpkg(iv(i),1)) then
                 i_sign=-1
              else if(xuse(i).gt.xpkg(iv(i)+1,1)) then
                 i_sign=+1
              else
                 iok(i)=.TRUE.
              end if
              !
              !  fifth iteration
              !
              if(.not.iok(i)) then
                 iv(i)=iv(i)+i_sign
                 if(xuse(i).lt.xpkg(iv(i),1)) then
                    i_sign=-1
                 else if(xuse(i).gt.xpkg(iv(i)+1,1)) then
                    i_sign=+1
                 else
                    iok(i)=.TRUE.
                 end if
                 if(.not.iok(i)) iprob=iprob+1
                 !
                 !  end chain of iteration if-then-else blocks
                 !
              end if
           end if
        end if
     end if
     !
     !  end of loop
     !
  end do

  go to 500
  !-----------------------------------------------------------
  !  algorithm:  piecewise linear indexing function lookup & correction
  !
350 continue
  do i=1,ivec
     !
     !  first iteration
     !
     if(init.eq.0) then
        if((xpkg(iprev,1).le.xuse(i)).and. &
             (xuse(i).le.xpkg(iprev+1,1))) then
           iok(i)=.TRUE.
           iv(i)=iprev
        else
           iok(i)=.FALSE.
           !
           !  piecewise linear indexing function on even spaced grid
           !  (not same grid as x axis itself)
           !
           iv(i)=1+havi*(xuse(i)-xpkg(1,1))
           iv(i)=max(1,min(nx-1,iv(i)))
           xfac=(xuse(i)-(xpkg(1,1)+(iv(i)-1)*hav))*havi
           zindx0=xpkg(iv(i),2)
           if(iv(i).lt.nx-1) then
              zdindx=xpkg(iv(i)+1,2)-zindx0
           else
              zdindx=nx-zindx0
           end if
           zindex=zindx0+xfac*zdindx
           iv(i)=zindex
           iv(i)=max(1,min(nx-1,iv(i)))
           if(xuse(i).lt.xpkg(iv(i),1)) then
              i_sign=-1
           else if(xuse(i).gt.xpkg(iv(i)+1,1)) then
              i_sign=+1
           else
              iok(i)=.TRUE.
           end if
        end if
     end if
     !
     !  second iteration
     !
     if(.not.iok(i)) then
        iv(i)=iv(i)+i_sign
        if(xuse(i).lt.xpkg(iv(i),1)) then
           i_sign=-1
        else if(xuse(i).gt.xpkg(iv(i)+1,1)) then
           i_sign=+1
        else
           iok(i)=.TRUE.
        end if
        !
        !  third iteration
        !
        if(.not.iok(i)) then
           iv(i)=iv(i)+i_sign
           if(xuse(i).lt.xpkg(iv(i),1)) then
              i_sign=-1
           else if(xuse(i).gt.xpkg(iv(i)+1,1)) then
              i_sign=+1
           else
              iok(i)=.TRUE.
           end if
           !
           !  fourth iteration
           !
           if(.not.iok(i)) then
              iv(i)=iv(i)+i_sign
              if(xuse(i).lt.xpkg(iv(i),1)) then
                 i_sign=-1
              else if(xuse(i).gt.xpkg(iv(i)+1,1)) then
                 i_sign=+1
              else
                 iok(i)=.TRUE.
              end if
              !
              !  fifth iteration
              !
              if(.not.iok(i)) then
                 iv(i)=iv(i)+i_sign
                 if(xuse(i).lt.xpkg(iv(i),1)) then
                    i_sign=-1
                 else if(xuse(i).gt.xpkg(iv(i)+1,1)) then
                    i_sign=+1
                 else
                    iok(i)=.TRUE.
                 end if
                 if(.not.iok(i)) iprob=iprob+1
                 !
                 !  end chain of iteration if-then-else blocks
                 !
              end if
           end if
        end if
     end if
     !
     !  end of loop
     !
     iprev=iv(i)
  end do

  go to 500
  !--------------------------------------------------------------------
  !  any "problems" left? if so, re-enter the loop...
  !
500 continue
  if(iprob.gt.0) go to 5
  !
  !  OK -- all zones found; complete output stats
  !
  if(imode.eq.1) then
     dxn=(xuse-xpkg(iv,1))          ! un-normalized
  else if(ialg.ne.3) then
     dxn=(xuse-xpkg(iv,1))*xpkg(iv,3) ! normalized to 1/h in each zone
     hv=xpkg(iv,2)
     hiv=xpkg(iv,3)
  else
     dxn=(xuse-xpkg(iv,1))*xpkg(iv,3) ! normalized to 1/h in each zone
     hv=xpkg(iv+1,1)-xpkg(iv,1)
     hiv=xpkg(iv,3)
  end if

  deallocate(iok)
  if(ialg.lt.3) then
     deallocate(imina)
     deallocate(imaxa)
  end if
  !
  !  all done -- generalized lookup
  !
1000 continue
  deallocate(xuse)

1100 continue
  return
end subroutine xlookup
