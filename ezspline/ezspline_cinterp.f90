subroutine EZspline_cinterp(spline_o, k, p1, p2, p3, f, ier)
  ! list of coordinate triplets
  ! this version of interp that has been inlined for better performance on the
  ! crays (pletzer@pppl.gov) Thu Jun 15 08:52:46 PDT 2000
  ! NOTE:
  ! This routine does not require the -dp switch on Crays since
  ! all constant have been made kind=r8.

  ! DMC Apr 2007 -- only the tricubic SPLINE evaluation has been inlined;
  ! if a different method interpolation is used, the standard f77-style
  ! pspline subroutine is called.

  use precision_mod, only: fp
  use EZspline_obj
  implicit none
  type(EZspline3) spline_o
  integer, intent(in) :: k
  real(fp), intent(in) :: p1(k), p2(k), p3(k)  ! location arrays
  real(fp), intent(out):: f(k)  ! interpolant array
  integer, intent(out) :: ier

  integer :: ifail
  integer, parameter :: ict(10)=(/1,0,0,0,0,0,0,0,0,0/)
  integer:: iwarn=0

  integer ix(k)                  ! zone indices {j}
  REAL(fp) dxn(k)       ! normalized displacements w/in zones
  REAL(fp) hx(k)        ! h(j) vector
  REAL(fp) hxi(k)       ! 1/h(j) vector
  !
  integer iy(k)                  ! zone indices {j}
  REAL(fp) dyn(k)       ! normalized displacements w/in zones
  REAL(fp) hy(k)        ! h(j) vector
  REAL(fp) hyi(k)       ! 1/h(j) vector
  !
  integer iz(k)                  ! zone indices {j}
  REAL(fp) dzn(k)       ! normalized displacements w/in zones
  REAL(fp) hz(k)        ! h(j) vector
  REAL(fp) hzi(k)       ! 1/h(j) vector
  integer v
  !
  REAL(fp) , parameter ::  sixth=0.166666666666666667_fp
  REAL(fp) , parameter ::  one=1.0_fp
  REAL(fp) , parameter ::  three=3.0_fp

  integer j1, j2 , j3
  REAL(fp) z36th,z216th,xp,xpi,xp2,xpi2,cx,cxi,hx2,cxd,cxdi,yp, &
       ypi,yp2,ypi2,cy,cyi,hy2,cyd,cydi,zp,zpi,zp2,zpi2,cz,czi,hz2, &
       czd,czdi,somme

  ier = 0
  ifail=0
  if(spline_o%isReady /= 1) then
    ier = 94
    return
  end if

  if (spline_o%isLinear == 1) then

    call vecpc3(ict, k, p1, p2, p3, k, f, &
         spline_o%n1, spline_o%x1pkg(1,1), &
         spline_o%n2, spline_o%x2pkg(1,1), &
         spline_o%n3, spline_o%x3pkg(1,1), &
         spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
         iwarn,ifail)

  else if (spline_o%isHybrid == 1) then

    call vecintrp3d(ict, k, p1, p2, p3, k, f, &
         spline_o%n1, spline_o%x1pkg(1,1), &
         spline_o%n2, spline_o%x2pkg(1,1), &
         spline_o%n3, spline_o%x3pkg(1,1), &
         spline_o%hspline, &
         spline_o%fspl(1,1,1,1), size(spline_o%fspl,1), &
         size(spline_o%fspl,2), size(spline_o%fspl,3), &
         size(spline_o%fspl,4), iwarn, ifail)

  else if (spline_o%isHermite == 0) then

    call xlookup(k,p1,spline_o%n1,spline_o%x1pkg(1,1),2,ix,dxn,hx,hxi,ifail)
    call xlookup(k,p2,spline_o%n2,spline_o%x2pkg(1,1),2,iy,dyn,hy,hyi,ifail)
    call xlookup(k,p3,spline_o%n3,spline_o%x3pkg(1,1),2,iz,dzn,hz,hzi,ifail)

    z36th=sixth*sixth
    z216th=sixth*sixth*sixth
    !
    !  prepare useful parameters...
    !
    do v=1,k
      j1=ix(v)
      j2=iy(v)
      j3=iz(v)
      !
      !   ...in x direction
      !
      xp=dxn(v)
      xpi=one-xp
      xp2=xp*xp
      xpi2=xpi*xpi
      !
      cx=xp*(xp2-ONE)
      cxi=xpi*(xpi2-ONE)
      hx2=hx(v)*hx(v)
      !
      !   ...and in y direction
      !
      yp=dyn(v)
      ypi=ONE-yp
      yp2=yp*yp
      ypi2=ypi*ypi
      !
      cy=yp*(yp2-ONE)
      cyi=ypi*(ypi2-ONE)
      hy2=hy(v)*hy(v)
      !
      !   ...and in z direction
      !
      zp=dzn(v)
      zpi=ONE-zp
      zp2=zp*zp
      zpi2=zpi*zpi
      !
      cz=zp*(zp2-ONE)
      czi=zpi*(zpi2-ONE)
      hz2=hz(v)*hz(v)
      !
      !  get desired values:
      !
      if(ict(1).eq.1) then
        !
        !  function value:
        !
        somme=(                                                       &
             zpi*(                                                    &
             xpi*(ypi*spline_o%fspl(1,j1,j2,j3)  +yp*spline_o%fspl(1,j1,j2+1,j3))+            &
             xp*(ypi*spline_o%fspl(1,j1+1,j2,j3)+yp*spline_o%fspl(1,j1+1,j2+1,j3)))          &
             +zp*(                                                    &
             xpi*(ypi*spline_o%fspl(1,j1,j2,j3+1)  +yp*spline_o%fspl(1,j1,j2+1,j3+1))+        &
             xp*(ypi*spline_o%fspl(1,j1+1,j2,j3+1)+yp*spline_o%fspl(1,j1+1,j2+1,j3+1))))
        !
        somme=somme+sixth*hx2*(                                       &
             zpi*(                                                    &
             cxi*(ypi*spline_o%fspl(2,j1,j2,j3)  +yp*spline_o%fspl(2,j1,j2+1,j3))+            &
             cx*(ypi*spline_o%fspl(2,j1+1,j2,j3)+yp*spline_o%fspl(2,j1+1,j2+1,j3)))          &
             +zp*(                                                    &
             cxi*(ypi*spline_o%fspl(2,j1,j2,j3+1)  +yp*spline_o%fspl(2,j1,j2+1,j3+1))+        &
             cx*(ypi*spline_o%fspl(2,j1+1,j2,j3+1)+yp*spline_o%fspl(2,j1+1,j2+1,j3+1))))
        !
        somme=somme+sixth*hy2*(                                       &
             zpi*(                                                    &
             xpi*(cyi*spline_o%fspl(3,j1,j2,j3)  +cy*spline_o%fspl(3,j1,j2+1,j3))+            &
             xp*(cyi*spline_o%fspl(3,j1+1,j2,j3)+cy*spline_o%fspl(3,j1+1,j2+1,j3)))          &
             +zp*(                                                    &
             xpi*(cyi*spline_o%fspl(3,j1,j2,j3+1)  +cy*spline_o%fspl(3,j1,j2+1,j3+1))+        &
             xp*(cyi*spline_o%fspl(3,j1+1,j2,j3+1)+cy*spline_o%fspl(3,j1+1,j2+1,j3+1))))
        !
        somme=somme+sixth*hz2*(                                       &
             czi*(                                                    &
             xpi*(ypi*spline_o%fspl(4,j1,j2,j3)  +yp*spline_o%fspl(4,j1,j2+1,j3))+            &
             xp*(ypi*spline_o%fspl(4,j1+1,j2,j3)+yp*spline_o%fspl(4,j1+1,j2+1,j3)))          &
             +cz*(                                                    &
             xpi*(ypi*spline_o%fspl(4,j1,j2,j3+1)  +yp*spline_o%fspl(4,j1,j2+1,j3+1))+        &
             xp*(ypi*spline_o%fspl(4,j1+1,j2,j3+1)+yp*spline_o%fspl(4,j1+1,j2+1,j3+1))))
        !
        somme=somme+z36th*hx2*hy2*(                                   &
             zpi*(                                                    &
             cxi*(cyi*spline_o%fspl(5,j1,j2,j3)  +cy*spline_o%fspl(5,j1,j2+1,j3))+            &
             cx*(cyi*spline_o%fspl(5,j1+1,j2,j3)+cy*spline_o%fspl(5,j1+1,j2+1,j3)))          &
             +zp*(                                                    &
             cxi*(cyi*spline_o%fspl(5,j1,j2,j3+1)  +cy*spline_o%fspl(5,j1,j2+1,j3+1))+        &
             cx*(cyi*spline_o%fspl(5,j1+1,j2,j3+1)+cy*spline_o%fspl(5,j1+1,j2+1,j3+1))))
           !
        somme=somme+z36th*hx2*hz2*(                                   &
             czi*(                                                    &
             cxi*(ypi*spline_o%fspl(6,j1,j2,j3)  +yp*spline_o%fspl(6,j1,j2+1,j3))+            &
             cx*(ypi*spline_o%fspl(6,j1+1,j2,j3)+yp*spline_o%fspl(6,j1+1,j2+1,j3)))          &
             +cz*(                                                    &
             cxi*(ypi*spline_o%fspl(6,j1,j2,j3+1)  +yp*spline_o%fspl(6,j1,j2+1,j3+1))+        &
             cx*(ypi*spline_o%fspl(6,j1+1,j2,j3+1)+yp*spline_o%fspl(6,j1+1,j2+1,j3+1))))
        !
        somme=somme+z36th*hy2*hz2*(                                   &
             czi*(                                                    &
             xpi*(cyi*spline_o%fspl(7,j1,j2,j3)  +cy*spline_o%fspl(7,j1,j2+1,j3))+            &
             xp*(cyi*spline_o%fspl(7,j1+1,j2,j3)+cy*spline_o%fspl(7,j1+1,j2+1,j3)))          &
             +cz*(                                                    &
             xpi*(cyi*spline_o%fspl(7,j1,j2,j3+1)  +cy*spline_o%fspl(7,j1,j2+1,j3+1))+        &
             xp*(cyi*spline_o%fspl(7,j1+1,j2,j3+1)+cy*spline_o%fspl(7,j1+1,j2+1,j3+1))))
        !
        somme=somme+z216th*hx2*hy2*hz2*(                              &
             czi*(                                                    &
             cxi*(cyi*spline_o%fspl(8,j1,j2,j3)  +cy*spline_o%fspl(8,j1,j2+1,j3))+            &
             cx*(cyi*spline_o%fspl(8,j1+1,j2,j3)+cy*spline_o%fspl(8,j1+1,j2+1,j3)))          &
             +cz*(                                                    &
             cxi*(cyi*spline_o%fspl(8,j1,j2,j3+1)  +cy*spline_o%fspl(8,j1,j2+1,j3+1))+        &
             cx*(cyi*spline_o%fspl(8,j1+1,j2,j3+1)+cy*spline_o%fspl(8,j1+1,j2+1,j3+1))))
        !
        f(v)=somme
      else
        write(6,*) 'ERROR ict(1) must be 1 (interpolation)'
      end if
      !
      !
    end do                             ! vector loop


    if(ifail /= 0) ier = 97

  else

    call vecherm3(ict, k, p1, p2, p3, k, f, &
         spline_o%n1, spline_o%x1pkg(1,1), &
         spline_o%n2, spline_o%x2pkg(1,1), &
         spline_o%n3, spline_o%x3pkg(1,1), &
         spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
         iwarn,ifail)

  end if

  if(ifail /= 0) ier = 97

end subroutine EZspline_cinterp
