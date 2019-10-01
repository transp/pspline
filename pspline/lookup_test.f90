program lookup_test
  use precision_mod, only: fp
  !
  !  timing test:  various lookup algorithms
  !
  !============
  IMPLICIT NONE
  INTEGER imul,istart,i,icur,k,j,i2,ig,ilinx,ier,idum1,idum2, &
       idum,iv,icount,iper,ia1,ia2,ia,is,ivs,ivecs,iloop,ir,iwarn, &
       igen,igenc,igen2,igen3
  !============
  real(fp) :: zfac,zscan,znorm,zval,zval2,zinc,zcur,zprev,zdum,y,x,z, &
       ztime0,ztime1,difmax,ztime
  !============
  integer gsize,g3size,tsize,oddstep,repeat
  parameter (gsize=100)             ! grid
  parameter (g3size=40)             ! grid for 3d object
  parameter (tsize=1000)            ! target vector size
  parameter (oddstep=37)            ! stepsize for scattered vector
  parameter (repeat=500000)         ! repeat count for timing
  !
  !  test grids
  !
  real(fp) :: grids(gsize,3),fspline(4,gsize,3),wk(tsize),ref(tsize)
  real(fp) :: gpkgs(gsize,4,3)
  !
  real(fp) :: gpkg2x(gsize,4),gpkg2y(gsize,4)
  real(fp) :: f2d(4,gsize,gsize)
  real(fp) :: gpkg3x(g3size,4),gpkg3y(g3size,4),gpkg3z(g3size,4)
  real(fp) :: f3d(8,g3size,g3size,g3size)
  !
  real(fp) :: wk2(tsize,2)
  real(fp) :: wk3(tsize,3)
  !
  !  test vectors
  !
  real(fp) :: tvecs(tsize,2)
  !
  real(fp) :: ztest1,ztest2,zdiff
  !
  integer ivecsize(4)
  !
  integer istop(6)
  !
  character(len=20) :: glbls(3)
  character(len=30) :: tlbls(2)
  character(len=20) :: albl(3)
  !
  integer ict(10)
  !
  ivecsize(1) = 1
  ivecsize(2) = 10
  ivecsize(3) = 100
  ivecsize(4) = 1000
  !
  albl(1) = 'pseudo-Newton'
  albl(2) = 'binary search'
  albl(3) = 'linear indexing fcn'
  !
  glbls(1) = 'evenly spaced grid'
  glbls(2) = 'smooth varying grid'
  glbls(3) = 'irregular grid'
  !
  tlbls(1) = 'modestly changing vector'
  tlbls(2) = 'rapidly changing vector'
  !
  ict = 0
  ict(1) = 1
  !
  ztest1 = 1.0_fp
  ztest2 = 1.0000000001_fp
  !
  !---------------------------------
  !
  !  form test vectors
  !
  write(6,*) ' 1d spline / lookup test program run on...'
#ifndef _PGC
  call execute_command_line('date')
  call execute_command_line('uname -a')
#endif

  zdiff=ztest2-ztest1
  write(6,*) ' ...test: ',ztest1,ztest2,' diff=',zdiff
  write(6,*) ' '
  write(6,*) ' total spline evaluations per test loop:  ',repeat
  !
  imul=0
  istart=1
  !
  do i=1,tsize
     zfac=(i-1.0D0)/(tsize-1.0D0)
     zscan=zfac*zfac*(3.0D0-2.0D0*zfac)
     tvecs(i,1)=max(0.0D0,min(100.0D0,(100.0D0*zscan)))
     !
     imul=imul+1
     if(imul.gt.4) imul=1
     istart=istart+imul*oddstep
     if(istart.gt.100) istart=istart-100
     if(istart.gt.100) istart=istart-100
     if(istart.gt.100) istart=istart-100
     tvecs(i,2)=max(0.0D0,min(100.0D0,(istart-0.5D0)))
  enddo
  !
  !  form test grids & packages
  !
  znorm=100.0D0/sqrt(100.0D0)
  do i=1,gsize
     zval=(i-1)*100.0D0/(gsize-1)
     zval2=sqrt(zval)*znorm
     grids(i,1)=zval
     grids(i,2)=zval2
     fspline(1,i,1)=10.0D0+cos(0.1D0*zval)
     fspline(1,i,2)=10.0D0+cos(0.1D0*zval2)
  enddo
  !
  zinc=100.0D0/(gsize-20-1)
  do i=1,6
     istop(i)=((gsize-20)*i)/7
  enddo
  !
  icur=0
  zcur=-zinc
  do i=1,gsize-20
     zprev=zcur
     zcur=zcur+zinc
     k=0
     do j=1,6
        if(i.eq.istop(j)) k=j
     enddo
     if(k.eq.0) then
        icur=icur+1
        grids(icur,3)=zcur
     else if(k.eq.1) then
        do i2=1,2
           icur=icur+1
           grids(icur,3)=zprev+0.5D0*i2*zinc
        enddo
     else if(k.eq.2) then
        do i2=1,4
           icur=icur+1
           grids(icur,3)=zprev+0.25D0*i2*zinc
        enddo
     else if(k.eq.3) then
        do i2=1,8
           icur=icur+1
           grids(icur,3)=zprev+0.125D0*i2*zinc
        enddo
     else if((k.eq.4).or.(k.eq.6)) then
        do i2=1,2
           icur=icur+1
           grids(icur,3)=zprev+0.25D0*i2*zinc
        enddo
        icur=icur+1
        grids(icur,3)=zcur
     else if(k.eq.5) then
        do i2=1,4
           icur=icur+1
           grids(icur,3)=zprev+0.125D0*i2*zinc
        enddo
        icur=icur+1
        grids(icur,3)=zcur
     end if
  enddo
  grids(icur+1,3)=grids(icur,3)+5.0D0
  do i=1,gsize
     fspline(1,i,3)=10.0D0+cos(0.1D0*grids(i,3))
  enddo
  !
  !  compute spline coeffs...
  !
  write(6,*) ' ... setup 1d splines ... '
  do ig=1,3
     call cspline(grids(1,ig),gsize,fspline(1,1,ig), &
          0,zdum,0,zdum,wk,gsize,ilinx,ier)
  enddo
  !
  !  setup grid packages for 2d and 3d tests
  !
  write(6,*) ' ... setup 2d & 3d test grids ... '
  do i=1,gsize
     wk(i)=grids(i,2)
  enddo
  call genxpkg(gsize,wk,gpkg2x,1,1,0,zdum,-3,ier)
  !
  do i=1,gsize
     wk(i)=grids(i,3)
  enddo
  call genxpkg(gsize,wk,gpkg2y,1,1,0,zdum,-3,ier)
  !
  do i=1,g3size
     wk(i)=grids(i,1)
  enddo
  call genxpkg(g3size,wk,gpkg3x,1,1,0,zdum,-3,ier)
  !
  do i=1,g3size
     wk(i)=grids(i,2)
  enddo
  call genxpkg(g3size,wk,gpkg3y,1,1,0,zdum,-3,ier)
  !
  do i=1,g3size
     wk(i)=grids(i,3)
  enddo
  call genxpkg(g3size,wk,gpkg3z,1,1,0,zdum,-3,ier)
  !
  write(6,*) ' ... setup 2d and 3d splines ... '
  do j=1,gsize
     y=gpkg2y(j,1)
     do i=1,gsize
        x=gpkg2x(i,1)
        f2d(1,i,j)=(10.0D0+cos(x))*(10.0D0+sin(2.0D0*y))
     enddo
  enddo
  !
  do k=1,g3size
     z=gpkg3z(k,1)
     do j=1,g3size
        y=gpkg3y(j,1)
        do i=1,g3size
           x=gpkg3x(i,1)
           f3d(1,i,j,k)=(10.0D0+cos(x))*(10.0D0+sin(2.0D0*y))* &
                exp(-0.1D0*z)
        enddo
     enddo
  enddo
  !
  call mkbicub(gpkg2x,gsize,gpkg2y,gsize,f2d,gsize, &
       0,zdum,0,zdum, &
       0,zdum,0,zdum, &
       idum1,idum2,ier)
  !
  call mktricub(gpkg3x,g3size,gpkg3y,g3size,gpkg3z,g3size, &
       f3d,g3size,g3size, &
       0,zdum,0,zdum,g3size, &
       0,zdum,0,zdum,g3size, &
       0,zdum,0,zdum,g3size, &
       idum,idum,idum,ier)
  !
  write(6,*) ' ... setup completed ... '
  !
  !-----------------------------------------------------------------
  !
  do ig=1,3
     write(6,*) '++++++++++++++++++++++++++++++++++++++++++++++++'
     write(6,*) ' ==> grid:  ',glbls(ig)
     write(6,*) ' '
     do iv=1,2
        icount=0
        write(6,*) '---------------------------------------------'
        write(6,*) ' ==> vector:  ',tlbls(iv)
        write(6,*) ' '
        do iper=0,1
           if(iper.eq.0) then
              write(6,*) ' ================='
              write(6,*) ' ...periodic:  NO'
           else
              write(6,*) ' ================='
              write(6,*) ' ...periodic:  YES'
           end if
           if(ig.eq.1) then
              ia1=1
              ia2=1
           else
              ia1=1
              ia2=3
           end if
           do ia=ia1,ia2
              if(ig.eq.1) then
                 write(6,*) ' ...even spacing:  YES'
              else
                 write(6,*) ' ...even spacing:  NO'
                 write(6,*) ' ...algorithm:  ',albl(ia)
              end if
              do is=1,-1,-2
                 write(6,*) ' '
                 if(is.lt.0) then
                    write(6,*) ' ...USE previous result as start ', &
                         'point for next search.'
                 else
                    write(6,*) ' ...DO NOT use previous search ', &
                         'result.'
                 end if
                 write(6,*) ' '
                 !
                 call genxpkg(gsize,grids(1,ig),gpkgs(1,1,ig), &
                      iper,1,0,zdum,is*ia,ier)
                 write(6,*) ' %genpkg(ngrid,grid,gpkg,',iper, &
                      ',1,0,zdum,',is*ia,',ier):  ier=',ier
                 !
                 write(6,*) ' '
                 do ivs=4,1,-1
                    ivecs=ivecsize(ivs)
                    !
                    call cpu_time(ztime0)
                    !
                    do iloop=1,repeat/tsize
                       do ir=1,tsize,ivecs
                          call spvec(ict,ivecs,tvecs(ir,iv), &
                               ivecs,wk(ir),gsize,gpkgs(1,1,ig), &
                               fspline(1,1,ig),iwarn,ier)
                          if(iwarn.ne.0) write(6,*) &
                               ' %xlookup iwarn=',iwarn,'; ir=',ir
                          if(ier.ne.0) write(6,*) &
                               ' %xlookup ier=',iwarn,'; ir=',ir
                       enddo
                    enddo
                    !
                    call cpu_time(ztime1)
                    !
                    icount=icount+1
                    if(icount.eq.1) then
                       difmax=0.0D0
                       do ir=1,tsize
                          !xx                              write(6,*) ir,wk(ir)
                          ref(ir)=wk(ir)
                       enddo
                    else
                       !  check variance
                       difmax=0.0D0
                       do ir=1,tsize
                          difmax=max(difmax,abs(ref(ir)-wk(ir)))
                       enddo
                       difmax=0.1D0*difmax
                       !
                    end if
                    !
                    ztime=ztime1-ztime0
                    write(6,1001) ztime,ivecs,difmax
1001                format(' $$$$$> cpu time (secs) = ',f8.4, &
                         ', vector size = ',i4,', variance =',1pe10.3)
1002                format(' $$$$$> cpu time (secs) = ',f8.4, &
                         ', vector size = ',i4)
                    !
                 enddo
                 !
                 write(6,*) ' '
                 !
              enddo
           enddo
        enddo
        !
     enddo
  enddo
  !
  write(6,*) ' ------------------------------------- '
  write(6,*) ' 2d test (linear fcn lookup algorithm)'
  do ivs=4,1,-1
     ivecs=ivecsize(ivs)
     !
     call cpu_time(ztime0)
     !
     do iloop=1,repeat/tsize
        igen = 1 + iloop - 2*(iloop/2)
        igenc=3-igen
        do ir=1,tsize,ivecs
           call vecbicub(ict,ivecs,tvecs(ir,igen),tvecs(ir,igenc), &
                ivecs,wk2(ir,igen),gsize,gpkg2x,gsize,gpkg2y, &
                f2d,gsize,iwarn,ier)
           if(iwarn.ne.0) write(6,*) &
                ' %xlookup iwarn=',iwarn,'; ir=',ir
           if(ier.ne.0) write(6,*) &
                ' %xlookup ier=',iwarn,'; ir=',ir
        enddo
     enddo
     !
     call cpu_time(ztime1)
     ztime=ztime1-ztime0
     write(6,1002) ztime,ivecs
  enddo
  !
  write(6,*) ' ------------------------------------- '
  write(6,*) ' 3d test (linear fcn lookup algorithm)'
  do ivs=4,1,-1
     ivecs=ivecsize(ivs)
     !
     call cpu_time(ztime0)
     !
     do iloop=1,repeat/tsize
        igen = 1 + iloop - 3*(iloop/3)
        igen2=igen+1
        igen3=igen+2
        if(igen2.gt.3) igen2=igen2-3
        if(igen3.gt.3) igen3=igen3-3
        igen=min(2,igen)            ! (1,2,2) -> (2,2,1) -> (2,1,2)
        igen2=min(2,igen2)
        igen3=min(2,igen3)
        do ir=1,tsize,ivecs
           call vectricub(ict,ivecs, &
                tvecs(ir,igen),tvecs(ir,igen2),tvecs(ir,igen3), &
                ivecs,wk3(ir,igen), &
                g3size,gpkg3x,g3size,gpkg3y,g3size,gpkg3z, &
                f3d,g3size,g3size,iwarn,ier)
           if(iwarn.ne.0) write(6,*) &
                ' %xlookup iwarn=',iwarn,'; ir=',ir
           if(ier.ne.0) write(6,*) &
                ' %xlookup ier=',iwarn,'; ir=',ir
        enddo
     enddo
     !
     call cpu_time(ztime1)
     ztime=ztime1-ztime0
     write(6,1002) ztime,ivecs
  enddo
  !
  stop
end program lookup_test
