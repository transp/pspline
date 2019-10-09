
!!!
!!! 1-d
!!!

subroutine EZspline_interp1(spline_o, p1, f, ier)
  use precision_mod, only: fp
  use EZspline_obj
  implicit none
  type(EZspline1) spline_o
  real(fp) p1  ! the location where the interpolation is sought
  real(fp) f   ! the interpolation
  integer, intent(out) :: ier

  integer ifail
  integer, parameter :: ict(3)=(/1, 0, 0/)
  real(fp) :: ansr(1)

  ier = 0
  ifail = 0
  if( .not.EZspline_allocated(spline_o) .or. spline_o%isReady /= 1) then
    ier =94
    return
  end if

  if (spline_o%isLinear == 1) then

    call pc1ev(p1, &
         spline_o%x1(1), spline_o%n1, &
         spline_o%ilin1, &
         spline_o%fspl(1,1), &
         ict, ansr, ifail)

  else if (spline_o%isHermite == 0) then

    call evspline(p1, &
         spline_o%x1(1), spline_o%n1, &
         spline_o%ilin1, &
         spline_o%fspl(1,1), &
         ict, ansr, ifail)

  else

    call herm1ev(p1, &
         spline_o%x1(1), spline_o%n1, &
         spline_o%ilin1,&
         spline_o%fspl(1,1),  &
         ict, ansr, ifail)

  end if

  f=ansr(1)

  if(ifail /= 0) ier = 97

end subroutine EZspline_interp1


subroutine EZspline_interp1_array(spline_o, k, p1, f, ier)
  use precision_mod, only: fp
  use EZspline_obj
  implicit none
  type(EZspline1) spline_o
  integer, intent(in) :: k
  real(fp), intent(in) :: p1(k) ! location arrays
  real(fp), intent(out):: f(k)  ! interpolant array
  integer, intent(out) :: ier

  integer :: ifail
  integer, parameter :: ict(3)=(/1,0,0/)
  integer:: iwarn=0

  ier = 0
  ifail = 0
  if( .not.EZspline_allocated(spline_o) .or. spline_o%isReady /= 1) then
    ier = 94
    return
  end if

  if (spline_o%isLinear == 1) then

    call vecpc1(ict, k, p1, k, f, &
         spline_o%n1,spline_o%x1pkg(1,1), &
         spline_o%fspl(1,1), &
         iwarn, ifail)

  else if (spline_o%isHermite == 0) then

    call vecspline(ict, k, p1, k, f, &
         spline_o%n1,spline_o%x1pkg(1,1), &
         spline_o%fspl(1,1), &
         iwarn, ifail)

  else

    call vecherm1(ict, k, p1, k, f, &
         spline_o%n1, spline_o%x1pkg(1,1), &
         spline_o%fspl(1,1), &
         iwarn,ifail)

  end if

  if(ifail /= 0) ier = 97

end subroutine EZspline_interp1_array

!!!
!!! 2-d
!!!

subroutine EZspline_interp2(spline_o, p1, p2, f, ier)
  use precision_mod, only: fp
  use EZspline_obj
  implicit none
  type(EZspline2) :: spline_o
  real(fp) :: p1, p2  ! the location where the interpolation is sought
  real(fp) :: f       ! the interpolation
  integer, intent(out) :: ier
  integer :: ifail
  integer, parameter :: ict(6)=(/1, 0, 0, 0, 0, 0 /)
  real(fp) :: ansr(1), ansr6(6)

  ier = 0
  ifail = 0
  if( .not.EZspline_allocated(spline_o) .or. spline_o%isReady /= 1) then
    ier =94
    return
  end if

  if (spline_o%isHybrid == 1) then

    call evintrp2d(p1, p2, &
         spline_o%x1(1), spline_o%n1, &
         spline_o%x2(1), spline_o%n2, &
         spline_o%hspline, spline_o%fspl(1,1,1), &
         size(spline_o%fspl,1), size(spline_o%fspl,2), &
         size(spline_o%fspl,3),  &
         ict, ansr6, ifail)
    ansr(1) = ansr6(1)

  else if (spline_o%isLinear == 1) then

    call pc2ev(p1, p2,  &
         spline_o%x1(1), spline_o%n1, &
         spline_o%x2(1), spline_o%n2, &
         spline_o%ilin1, spline_o%ilin2, &
         spline_o%fspl(1,1,1), spline_o%n1, &
         ict, ansr, ifail)

  else if (spline_o%isHermite == 0) then

    call evbicub(p1, p2,  &
         spline_o%x1(1), spline_o%n1, &
         spline_o%x2(1), spline_o%n2, &
         spline_o%ilin1, spline_o%ilin2, &
         spline_o%fspl(1,1,1), spline_o%n1, &
         ict, ansr, ifail)

  else

    call herm2ev(p1, p2,  &
         spline_o%x1(1), spline_o%n1, &
         spline_o%x2(1), spline_o%n2, &
         spline_o%ilin1, spline_o%ilin2, &
         spline_o%fspl(1,1,1), spline_o%n1, &
         ict, ansr, ifail)

  end if

  f=ansr(1)

  if(ifail /= 0) ier = 97

end subroutine EZspline_interp2

subroutine EZspline_interp2_array(spline_o, k1, k2, p1, p2, f, ier)
  use precision_mod, only: fp
  use EZspline_obj
  implicit none
  type(EZspline2) spline_o
  integer, intent(in) :: k1, k2
  real(fp) :: p1(k1), p2(k2) ! location arrays
  real(fp) :: f(k1,k2)  ! interpolated function array
  integer, intent(out) :: ier

  integer :: ifail
  integer:: iwarn=0

  ier = 0
  ifail = 0
  if( .not.EZspline_allocated(spline_o) .or. spline_o%isReady /= 1) then
    ier =94
    return
  end if

  if (spline_o%isHybrid == 1) then

    call gridintrp2d( &
         p1, k1, &
         p2, k2, &
         f, k1, &
         spline_o%n1, spline_o%x1pkg(1,1), &
         spline_o%n2, spline_o%x2pkg(1,1), &
         spline_o%hspline, spline_o%fspl(1,1,1), &
         size(spline_o%fspl,1), size(spline_o%fspl,2), &
         size(spline_o%fspl,3),  &
         iwarn, ifail)

  else if (spline_o%isLinear == 1) then

    call gridpc2( &
         p1, k1, &
         p2, k2, &
         f, k1, &
         spline_o%n1, spline_o%x1pkg(1,1), &
         spline_o%n2, spline_o%x2pkg(1,1), &
         spline_o%fspl(1,1,1), spline_o%n1, &
         iwarn, ifail)

  else if (spline_o%isHermite == 0) then

    call gridbicub( &
         p1, k1, &
         p2, k2, &
         f, k1, &
         spline_o%n1, spline_o%x1pkg(1,1), &
         spline_o%n2, spline_o%x2pkg(1,1), &
         spline_o%fspl(1,1,1), spline_o%n1, &
         iwarn, ifail)

  else

    call gridherm2( &
         p1, k1, &
         p2, k2, &
         f, k1, &
         spline_o%n1, spline_o%x1pkg(1,1), &
         spline_o%n2, spline_o%x2pkg(1,1), &
         spline_o%fspl(1,1,1), spline_o%n1, &
         iwarn, ifail)

  end if

  if(ifail /= 0) ier = 97

end subroutine EZspline_interp2_array

subroutine EZspline_interp2_cloud(spline_o, k, p1, p2, f, ier)
  use precision_mod, only: fp
  ! list of coordinate doublets
  use EZspline_obj
  implicit none
  type(EZspline2) spline_o
  integer, intent(in) :: k
  real(fp), intent(in) :: p1(k), p2(k) ! location arrays
  real(fp), intent(out):: f(k)  ! interpolant array
  integer, intent(out) :: ier
  integer :: ifail
  integer, parameter :: ict(6) = (/1,0,0,0,0,0/)
  integer:: iwarn = 0

  ier = 0
  ifail = 0

  if( .not.EZspline_allocated(spline_o) .or. spline_o%isReady /= 1) then
    ier = 94
    return
  end if


  if (spline_o%isHybrid == 1) then

    call vecintrp2d(ict, k, p1, p2, k, f, &
         spline_o%n1, spline_o%x1pkg(1,1), &
         spline_o%n2, spline_o%x2pkg(1,1), &
         spline_o%hspline, spline_o%fspl(1,1,1), &
         size(spline_o%fspl,1), size(spline_o%fspl,2), &
         size(spline_o%fspl,3),  &
         iwarn, ifail)

  else if (spline_o%isLinear == 1) then

    call vecpc2(ict, k, p1, p2, k, f, &
         spline_o%n1, spline_o%x1pkg(1,1), &
         spline_o%n2, spline_o%x2pkg(1,1), &
         spline_o%fspl(1,1,1), spline_o%n1, &
         iwarn, ifail)

  else if (spline_o%isHermite == 0) then

    call vecbicub(ict, k, p1, p2, k, f, &
         spline_o%n1, spline_o%x1pkg(1,1), &
         spline_o%n2, spline_o%x2pkg(1,1), &
         spline_o%fspl(1,1,1), spline_o%n1, &
         iwarn, ifail)

  else

    call vecherm2(ict, k, p1, p2, k, f, &
         spline_o%n1, spline_o%x1pkg(1,1), &
         spline_o%n2, spline_o%x2pkg(1,1), &
         spline_o%fspl(1,1,1), spline_o%n1, &
         iwarn, ifail)

  end if

  if(ifail /= 0) ier = 97

end subroutine EZspline_interp2_cloud


!!!
!!! 3-d
!!!

subroutine EZspline_interp3(spline_o, p1, p2, p3, f, ier)
  use precision_mod, only: fp
  use EZspline_obj
  implicit none
  type(EZspline3) spline_o
  real(fp) p1, p2, p3 ! the location where the interpolation is sought
  real(fp) f          ! the interpolation

  integer, intent(out) :: ier
  integer ifail
  integer, parameter :: ict(10)=(/1,0,0,0,0,0,0,0,0,0/)
  real(fp) :: ansr(1), ansr10(10)

  ier = 0
  ifail=0
  if( .not.EZspline_allocated(spline_o) .or. spline_o%isReady /= 1) then
    ier =94
    return
  end if

  if (spline_o%isHybrid == 1) then

    call evintrp3d(p1, p2, p3, &
         spline_o%x1(1), spline_o%n1, &
         spline_o%x2(1), spline_o%n2, &
         spline_o%x3(1), spline_o%n3, &
         spline_o%hspline, spline_o%fspl(1,1,1,1), &
         size(spline_o%fspl,1), size(spline_o%fspl,2), &
         size(spline_o%fspl,3), size(spline_o%fspl,4), &
         ict, ansr10, ifail)
    ansr(1) = ansr10(1)

  else if (spline_o%isLinear == 1) then

    call pc3ev(p1, p2, p3, &
         spline_o%x1(1), spline_o%n1, &
         spline_o%x2(1), spline_o%n2, &
         spline_o%x3(1), spline_o%n3, &
         spline_o%ilin1, spline_o%ilin2, spline_o%ilin3, &
         spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
         ict, ansr, ifail)

  else if (spline_o%isHermite == 0) then

    call evtricub(p1, p2, p3, &
         spline_o%x1(1), spline_o%n1, &
         spline_o%x2(1), spline_o%n2, &
         spline_o%x3(1), spline_o%n3, &
         spline_o%ilin1, spline_o%ilin2, spline_o%ilin3, &
         spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
         ict, ansr, ifail)

  else

    call herm3ev(p1, p2, p3, &
         spline_o%x1(1), spline_o%n1, &
         spline_o%x2(1), spline_o%n2, &
         spline_o%x3(1), spline_o%n3, &
         spline_o%ilin1, spline_o%ilin2, spline_o%ilin3, &
         spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
         ict, ansr, ifail)

  end if

  f=ansr(1)

  if(ifail /= 0) ier = 97

end subroutine EZspline_interp3

subroutine EZspline_interp3_array(spline_o, k1, k2, k3, p1, p2, p3, f, ier)
  use precision_mod, only: fp
  use EZspline_obj
  implicit none
  type(EZspline3) spline_o
  integer :: k1, k2, k3
  real(fp) :: p1(k1), p2(k2), p3(k3)  ! location arrays
  real(fp) :: f(k1,k2,k3)  ! interpolant array
  integer, intent(out) :: ier

  integer ifail
  integer:: iwarn=0

  ier = 0
  ifail=0
  if( .not.EZspline_allocated(spline_o) .or. spline_o%isReady /= 1) then
    ier =94
    return
  end if

  if (spline_o%isHybrid == 1) then

    call gridintrp3d( &
         p1, k1, &
         p2, k2, &
         p3, k3, &
         f, k1, k2, &
         spline_o%n1, spline_o%x1pkg(1,1), &
         spline_o%n2, spline_o%x2pkg(1,1), &
         spline_o%n3, spline_o%x3pkg(1,1), &
         spline_o%hspline, spline_o%fspl(1,1,1,1), &
         size(spline_o%fspl,1), size(spline_o%fspl,2), &
         size(spline_o%fspl,3), size(spline_o%fspl,4), &
         iwarn, ifail)

  else if (spline_o%isLinear == 1) then

    call gridpc3( &
         p1, k1, &
         p2, k2, &
         p3, k3, &
         f, k1, k2, &
         spline_o%n1, spline_o%x1pkg(1,1), &
         spline_o%n2, spline_o%x2pkg(1,1), &
         spline_o%n3, spline_o%x3pkg(1,1), &
         spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
         iwarn, ifail)

  else if (spline_o%isHermite == 0) then
    !
    call gridtricub( &
         p1, k1, &
         p2, k2, &
         p3, k3, &
         f, k1, k2, &
         spline_o%n1, spline_o%x1pkg(1,1), &
         spline_o%n2, spline_o%x2pkg(1,1), &
         spline_o%n3, spline_o%x3pkg(1,1), &
         spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
         iwarn,ifail)

  else

    call gridherm3( &
         p1, k1, &
         p2, k2, &
         p3, k3, &
         f, k1, k2, &
         spline_o%n1, spline_o%x1pkg(1,1), &
         spline_o%n2, spline_o%x2pkg(1,1), &
         spline_o%n3, spline_o%x3pkg(1,1), &
         spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
         iwarn, ifail)

  end if

  if(ifail /= 0) ier = 97

end subroutine EZspline_interp3_array

subroutine EZspline_interp3_cloud(spline_o, k, p1, p2, p3, f, ier)
  use precision_mod, only: fp
  ! list of coordinate triplets
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

  ier = 0
  ifail=0
  if( .not.EZspline_allocated(spline_o) .or. spline_o%isReady /= 1) then
    ier = 94
    return
  end if

  if (spline_o%isHybrid == 1) then

    call vecintrp3d(ict, k, p1, p2, p3, k, f, &
         spline_o%n1, spline_o%x1pkg(1,1), &
         spline_o%n2, spline_o%x2pkg(1,1), &
         spline_o%n3, spline_o%x3pkg(1,1), &
         spline_o%hspline, spline_o%fspl(1,1,1,1), &
         size(spline_o%fspl,1), size(spline_o%fspl,2), &
         size(spline_o%fspl,3), size(spline_o%fspl,4), &
         iwarn, ifail)

  else if (spline_o%isLinear == 1) then

    call vecpc3(ict, k, p1, p2, p3, k, f, &
         spline_o%n1, spline_o%x1pkg(1,1), &
         spline_o%n2, spline_o%x2pkg(1,1), &
         spline_o%n3, spline_o%x3pkg(1,1), &
         spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
         iwarn,ifail)

  else if (spline_o%isHermite == 0) then
    !
    call vectricub(ict, k, p1, p2, p3, k, f, &
         spline_o%n1,spline_o%x1pkg(1,1), &
         spline_o%n2,spline_o%x2pkg(1,1), &
         spline_o%n3,spline_o%x3pkg(1,1), &
         spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
         iwarn, ifail)

  else

    call vecherm3(ict, k, p1, p2, p3, k, f, &
         spline_o%n1, spline_o%x1pkg(1,1), &
         spline_o%n2, spline_o%x2pkg(1,1), &
         spline_o%n3, spline_o%x3pkg(1,1), &
         spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
         iwarn,ifail)

  end if

  if(ifail /= 0) ier = 97

end subroutine EZspline_interp3_cloud
