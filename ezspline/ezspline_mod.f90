module EZspline_type
  use precision_mod, only: fp

  real(fp), parameter :: ezspline_twopi = 6.2831853071795865_fp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! EZspline data types
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  type EZspline3
    !
    ! 3-d Spline/Akima Hermite/Piecewise Linear interpolations
    !
    ! Grid
    !
    real(fp), dimension(:), allocatable :: x1, x2, x3
    !
    ! The boundary condition values (for slope and 2nd derivative).
    ! Can be optionally set by the user. Not used for periodic and
    ! not a knot boundary conditions.
    !
    real(fp), dimension(:,:), allocatable :: bcval1min, bcval1max
    real(fp), dimension(:,:), allocatable :: bcval2min, bcval2max
    real(fp), dimension(:,:), allocatable :: bcval3min, bcval3max
    !
    ! Select between spline (0) and Akima spline (1); default=0 (spline)
    !
    integer :: isHermite  ! set after EZspline_init call...
    !
    ! set =0 for Spline, Akima or Hybrid; =1 for piecewise linear: this is set
    ! by EZspline_init, EZhybrid_init, or EZlinear_init; DO NOT SET DIRECTLY:
    !
    integer :: isLinear
    !
    ! set =0 by init routines other than EZhybrid_init which sets it =1:
    integer :: isHybrid
    !
    ! the following is set by EZhybrid_init; other EZ*_init routines clear:
    integer :: hspline(3)  ! interpolation code along each dimension
    !        -1: zonal step fcn; =0: pc linear; =1: Akima Hermite; =2: Spline
    !
    ! Grid sizes (set during EZ*_init call).
    !
    integer :: n1, n2, n3
    !
    ! Grid zone lookup method
    !
    integer :: klookup1,klookup2,klookup3
    !
    ! Type of boundary conditions (set during EZspline_init call) on left
    ! and right hand side. Possible values are:
    !
    ! -1 periodic
    ! 0 not a knot
    ! +1 1st derivative imposed
    ! +2 2nd derivative imposed
    !
    ! For instance, ibctype1 =(/1, 0/) for 1st derivative set on left-hand
    ! and not a knot boundary conditions on right-hand side. The values of
    ! the derivatives a set via  bcval1min. (See above.)
    !
    integer ibctype1(2), ibctype2(2), ibctype3(2)
    !
    ! Grid lengths. DO NOT SET.
    !
    real(fp) :: x1min, x1max, x2min, x2max, x3min, x3max
    !
    ! Compact cubic coefficient arrays. DO NOT SET.
    !
    real(fp), dimension(:,:,:,:), allocatable :: fspl
    !
    ! Control/Other. DO NOT SET.
    !
    integer :: isReady

    integer :: ilin1, ilin2, ilin3
    real(fp), dimension(:,:), allocatable :: x1pkg, x2pkg, x3pkg
    !
    integer :: nguard
  end type EZspline3

  type EZspline2
    !
    ! 2-d Spline/Akima Hermite/Piecewise Linear interpolation
    !
    ! Grid
    !
    real(fp), dimension(:), allocatable :: x1, x2
    !
    ! The boundary condition values (for slope and 2nd derivative).
    ! Can be optionally set by the user. Not used for periodic and
    ! not a knot boundary conditions.
    !
    real(fp), dimension(:), allocatable :: bcval1min, bcval1max
    real(fp), dimension(:), allocatable :: bcval2min, bcval2max
    !
    ! Select between spline (0) and Akima spline (1); default=0 (spline)
    !
    integer :: isHermite  ! set after EZspline_init call...
    !
    ! set =0 for Spline, Akima or Hybrid; =1 for piecewise linear: this is set
    ! by EZspline_init, EZhybrid_init, or EZlinear_init; DO NOT SET DIRECTLY:
    !
    integer :: isLinear
    !
    ! set =0 by init routines other than EZhybrid_init which sets it =1:
    integer :: isHybrid
    !
    ! the following is set by EZhybrid_init; other EZ*_init routines clear:
    integer :: hspline(2)  ! interpolation code along each dimension
    !        -1: zonal step fcn; =0: pc linear; =1: Akima Hermite; =2: Spline
    !
    ! Grid sizes (set during EZ*_init call).
    !
    integer :: n1, n2
    !
    ! Grid zone lookup method
    !
    integer :: klookup1,klookup2
    !
    ! Type of boundary conditions (set during EZspline_init call) on left
    ! and right hand side. Possible values are:
    !
    ! -1 periodic
    ! 0 not a knot
    ! +1 1st derivative imposed
    ! +2 2nd derivative imposed
    !
    ! For instance, ibctype1 =(/1, 0/) for 1st derivative set on left-hand
    ! and not a knot boundary conditions on right-hand side. The values of
    ! the derivatives are set via  bcval1min. (See above)
    !
    integer ibctype1(2), ibctype2(2)
    !
    ! Grid lengths. DO NOT SET.
    !
    real(fp) :: x1min, x1max, x2min, x2max
    !
    ! Compact cubic coefficient arrays. DO NOT SET.
    !
    real(fp), dimension(:,:,:), allocatable :: fspl
    !
    ! Control/Other. DO NOT SET.
    !
    integer :: isReady

    integer :: ilin1, ilin2
    real(fp), dimension(:,:), allocatable :: x1pkg, x2pkg
    !
    integer :: nguard
  end type EZspline2

  type EZspline1
    !
    ! 1-d Spline/Akima Hermite/Piecewise Linear interpolation
    !
    ! Grid
    !
    real(fp), dimension(:), allocatable :: x1
    !
    ! The boundary condition values (for slope and 2nd derivative).
    ! Can be optionally set by the user. Not used for periodic and
    ! not a knot boundary conditions.
    !
    real(fp) :: bcval1min, bcval1max
    !
    ! Select between spline (0) and Akima spline (1); default=0 (spline)
    !
    integer :: isHermite  ! set after EZspline_init call...
    !
    ! set =0 for Spline or Akima; =1 for piecewise linear: this is set
    ! by EZspline_init or EZlinear_init; DO NOT SET DIRECTLY:
    !
    integer :: isLinear
    !
    ! Grid sizes (set during EZ*_init call).
    !
    integer :: n1
    !
    ! Grid zone lookup method
    !
    integer :: klookup1
    !
    ! Type of boundary conditions (set during EZspline_init call) on left
    ! and right hand side. Possible values are:
    !
    ! -1 periodic
    ! 0 not a knot
    ! +1 1st derivative imposed
    ! +2 2nd derivative imposed
    !
    ! For instance, ibctype1 =(/1, 0/) for 1st derivative set on left-hand
    ! and not a knot boundary conditions on right-hand side. The values of
    ! the derivatives are set via  bcval1min. (See above)
    !
    integer ibctype1(2)
    !
    ! Grid lengths. DO NOT SET.
    !
    real(fp) :: x1min, x1max
    !
    ! Compact cubic coefficient arrays. DO NOT SET.
    !
    real(fp), dimension(:,:), allocatable :: fspl
    !
    ! Control/Other. DO NOT SET.
    !
    integer :: isReady

    integer :: ilin1
    real(fp), dimension(:,:), allocatable :: x1pkg
    !
    integer :: nguard
  end type EZspline1

  !=========================================================================
  ! End type
end module EZspline_type

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! EZspline reflective methods (first argument is an EZspline1,2,3 type).
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module EZspline_obj
  use precision_mod, only: fp
  use Ezspline_type
  interface EZspline_preInit
    !
    ! usage:  call EZspline_preinit(spline_o)
    !   ...where spline_o is a 1d, 2d, or 3d spline object.

    module procedure &
         EZspline_preInit1, &
         EZspline_preInit2, &
         EZspline_preInit3

  end interface

  interface EZspline_allocated
    ! logical function returns TRUE if allocated, FALSE otherwise
    !
    ! usage:  
    !   logical :: answer
    !   answer = EZspline_allocated(spline_o)
    !   ...where spline_o is a 1d, 2d, or 3d spline object.

    module procedure &
         EZspline_allocated1, &
         EZspline_allocated2, &
         EZspline_allocated3
  end interface

contains

  subroutine EZspline_preInit1(spline_o)
    use EZspline_type
    type(EZspline1) :: spline_o
    spline_o%nguard=123456789
    spline_o%isReady=0
    spline_o%ibctype1=0 
  end subroutine EZspline_preInit1

  subroutine EZspline_preInit2(spline_o)
    use EZspline_type
    type(EZspline2) :: spline_o
    spline_o%nguard=123456789
    spline_o%isReady=0
    spline_o%ibctype1=0 ; spline_o%ibctype2=0
  end subroutine EZspline_preInit2

  subroutine EZspline_preInit3(spline_o)
    use EZspline_type
    type(EZspline3) spline_o
    spline_o%nguard=123456789
    spline_o%isReady=0
    spline_o%ibctype1=0 ; spline_o%ibctype2=0 ; spline_o%ibctype3=0
  end subroutine EZspline_preInit3

  logical function EZspline_allocated1(spline_o)
    use EZspline_type
    type(EZspline1) spline_o
    EZspline_allocated1 = allocated(spline_o%fspl) &
         & .and. allocated(spline_o%x1) .and. allocated(spline_o%x1pkg) &
         & .and. (spline_o%nguard == 123456789) ! check that ezspline_init has been called
  end function EZspline_allocated1

  logical function EZspline_allocated2(spline_o)
    use EZspline_type
    type(EZspline2) spline_o
    EZspline_allocated2 = allocated(spline_o%fspl) &
         & .and. allocated(spline_o%x1) .and. allocated(spline_o%x1pkg) &
         & .and. allocated(spline_o%x2) .and. allocated(spline_o%x2pkg) &
         & .and. (spline_o%nguard == 123456789) ! check that ezspline_init has been called
  end function EZspline_allocated2

  logical function EZspline_allocated3(spline_o)
    use ezspline_type
    type(EZspline3) spline_o
    EZspline_allocated3 = allocated(spline_o%fspl) &
         & .and. allocated(spline_o%x1) .and. allocated(spline_o%x1pkg) &
         & .and. allocated(spline_o%x2) .and. allocated(spline_o%x2pkg) &
         & .and. allocated(spline_o%x3) .and. allocated(spline_o%x3pkg) &
         & .and. (spline_o%nguard == 123456789) ! check that ezspline_init has been called
  end function EZspline_allocated3

  subroutine ezmake_ict1(i,ict)
    !  (private utility for ezspline derivative2 subroutines)
    !  make ict(1:6) array
    !  for higher derivatives; d[i1+i2]f/dx[i1]dy[i2]
    !  expecting i in range [0:3] (NOT CHECKED)

    implicit NONE
    integer, intent(in) :: i
    integer, intent(out) :: ict(3)

    if(i.eq.0) then
      ict = (/1, 0, 0 /) ! seek f @ (p1)
    else if(i.eq.1) then
      ict = (/0, 1, 0 /) ! df/dx
    else if(i.eq.2) then
      ict = (/0, 0, 1 /) ! d2f/dx2
    else
      ict = (/3, 0, 0 /) ! d3f/dx3
    end if

  end subroutine ezmake_ict1

  subroutine ezmake_ict2(i1,i2,ict)
    !  (private utility for ezspline derivative2 subroutines)
    !  make ict(1:6) array
    !  for higher derivatives; d[i1+i2]f/dx[i1]dy[i2]
    !  expecting i1 & i2 in range [0:3] (NOT CHECKED)

    implicit NONE
    integer, intent(in) :: i1,i2
    integer, intent(out) :: ict(6)

    integer :: imark,isum,iii

    !  this generates the control argument needed by evbicub & similar
    !  routines...
    !----------------------

    isum = i1+i2
    ict(1)=isum

    imark=0

    if(isum.eq.0) then
      ict = (/1, 0, 0, 0, 0, 0 /) ! seek f @ (p1, p2, p3)
    else if(isum.eq.1) then
      if(i1.eq.1) then
        ict = (/0, 1, 0, 0, 0, 0 /) ! df/dx
      else
        ict = (/0, 0, 1, 0, 0, 0 /) ! df/dy
      end if
    else if(isum.eq.2) then
      if(i1.eq.2) then
        ict = (/0, 0, 0, 1, 0, 0 /) ! d2f/dx2
      else if(i2.eq.2) then
        ict = (/0, 0, 0, 0, 1, 0 /) ! d2f/dy2
      else
        ict = (/0, 0, 0, 0, 0, 1 /) ! d2f/dxdy
      end if
    else if(isum.eq.3) then
      if(i1.eq.3) then
        imark=2  ! fxxx
      else if(i1.eq.2) then
        imark=3  ! fxxy
      else if(i1.eq.1) then
        imark=4  ! fxyy
      else
        imark=5  ! fyyy
      end if
    else if(isum.eq.4) then
      if(i1.eq.3) then
        imark=2  ! fxxxy
      else if(i2.eq.3) then
        imark=4  ! fxyyy
      else
        imark=3  ! fxxyy
      end if
    else if(isum.eq.5) then
      if(i1.eq.3) then
        imark=2  ! fxxxyy
      else if(i2.eq.3) then
        imark=3  ! fxxyyy
      end if
    end if

    !  isum=6 --> fxxxyyy

    if(isum.gt.2) then
      do iii=2,6
        if(iii.eq.imark) then
          ict(iii)=1
        else
          ict(iii)=0
        end if
      end do
    end if

  end subroutine ezmake_ict2

  subroutine ezmake_ict3(i1,i2,i3,ict)
    !  (private utility for ezspline derivative3 subroutines)
    !  make ict(1:10) array
    !  for higher derivatives; d[i1+i2+i3]f/dx[i1]dy[i2]dz[i3]
    !  i1 & i2 & i3 in range [0:3] (NOT CHECKED)

    implicit NONE
    integer, intent(in) :: i1,i2,i3
    integer, intent(out) :: ict(10)

    integer :: imark,isum,iii

    !  this generates the control argument needed by evtricub & similar
    !  routines...
    !----------------------

    isum = i1+i2+i3
    if(max(i1,i2,i3).eq.3) then
      isum=-isum
    end if
    ict(1)=isum

    imark=0

    if(isum.eq.0) then
      ict = (/1, 0, 0, 0, 0, 0, 0, 0, 0, 0 /) ! seek f @ (p1, p2, p3)

    else if(isum.eq.1) then
      !  1st derivatives
      if(i1.eq.1) then
        ict = (/0, 1, 0, 0, 0, 0, 0, 0, 0, 0 /)  ! df/dx
      else if(i2.eq.1) then
        ict = (/0, 0, 1, 0, 0, 0, 0, 0, 0, 0 /)  ! df/dy
      else
        ict = (/0, 0, 0, 1, 0, 0, 0, 0, 0, 0 /)  ! df/dz
      end if

    else if(isum.eq.2) then
      !  2nd derivatives -- legacy ordering; x-precedence ordering for all
      !  higher derivatives...

      if(i1.eq.2) then
        ict = (/0, 0, 0, 0, 1, 0, 0, 0, 0, 0 /)  ! d2f/dx2
      else if(i2.eq.2) then
        ict = (/0, 0, 0, 0, 0, 1, 0, 0, 0, 0 /)  ! d2f/dy2
      else if(i3.eq.2) then
        ict = (/0, 0, 0, 0, 0, 0, 1, 0, 0, 0 /)  ! d2f/dz2
      else if(i3.eq.0) then
        ict = (/0, 0, 0, 0, 0, 0, 0, 1, 0, 0 /)  ! d2f/dxdy
      else if(i2.eq.0) then
        ict = (/0, 0, 0, 0, 0, 0, 0, 0, 1, 0 /)  ! d2f/dxdz
      else
        ict = (/0, 0, 0, 0, 0, 0, 0, 0, 0, 1 /)  ! d2f/dydz
      end if

    else if(isum.eq.3) then
      !  3rd derivative, continuous: max(i1,i2,i3)<3
      if(i1.eq.2) then
        if(i2.eq.1) then
          imark=2     ! fxxy
        else
          imark=3     ! fxxz
        end if
      else if(i1.eq.1) then
        if(i2.eq.2) then
          imark=4     ! fxyy
        else if(i2.eq.1) then
          imark=5     ! fxyz
        else
          imark=6     ! fxzz
        end if
      else
        if(i2.eq.2) then
          imark=7     ! fyyz
        else
          imark=8     ! fyzz
        end if
      end if

    else if(isum.eq.-3) then
      !  3rd derivative
      if(i1.eq.3) then
        imark=2        ! fxxx
      else if(i2.eq.3) then
        imark=3        ! fyyy
      else if(i3.eq.3) then
        imark=4        ! fzzz
      end if

    else if(isum.eq.4) then
      !  4th derivative, continuous: max(i1,i2,i3)<3
      if(i1.eq.2) then
        if(i2.eq.2) then
          imark=2     ! fxxyy
        else if(i2.eq.1) then
          imark=3     ! fxxyz
        else
          imark=4     ! fxxzz
        end if
      else if(i1.eq.1) then
        if(i2.eq.2) then
          imark=5     ! fxyyz
        else
          imark=6     ! fxyzz
        end if
      else
        imark=7        ! fyyzz
      end if

    else if(isum.eq.-4) then
      !  4th derivative
      if(i1.eq.3) then
        if(i2.eq.1) then
          imark=2     ! fxxxy
        else
          imark=3     ! fxxxz
        end if
      else if(i1.eq.1) then
        if(i2.eq.3) then
          imark=4     ! fxyyy
        else
          imark=5     ! fxzzz
        end if
      else
        if(i2.eq.3) then
          imark=6     ! fyyyz
        else
          imark=7     ! fyzzz
        end if
      end if

    else if(isum.eq.5) then
      !  5th derivative, continuous: max(i1,i2,i3)<3
      if(i3.eq.1) then
        imark=2     ! fxxyyz
      else if(i2.eq.1) then
        imark=3     ! fxxyzz
      else
        imark=4     ! fxyyzz
      end if

    else if(isum.eq.-5) then
      !  5th derivative
      if(i1.eq.3) then
        if(i2.eq.2) then
          imark=2  ! fxxxyy
        else if(i2.eq.1) then
          imark=3  ! fxxxyz
        else
          imark=4  ! fxxxzz
        end if
      else if(i1.eq.2) then
        if(i2.eq.3) then
          imark=5  ! fxxyyy
        else
          imark=6  ! fxxzzz
        end if
      else if(i1.eq.1) then
        if(i2.eq.3) then
          imark=7  ! fxyyyz
        else
          imark=8  ! fxyzzz
        end if
      else
        if(i2.eq.3) then
          imark=9  ! fyyyzz
        else
          imark=10 ! fyyzzz
        end if
      end if

      !  isum=6 --> fxxyyzz  (i1=i2=i3=2)
    else if(isum.eq.-6) then
      !  6th derivative
      if(i1.eq.3) then
        if(i2.eq.3) then
          imark=2  ! fxxxyyy
        else if(i2.eq.2) then
          imark=3  ! fxxxyyz
        else if(i2.eq.1) then
          imark=4  ! fxxxyzz
        else
          imark=5  ! fxxxzzz
        end if
      else if(i1.eq.2) then
        if(i2.eq.3) then
          imark=6  ! fxxyyyz
        else if(i2.eq.1) then
          imark=7  ! fxxyzzz
        end if
      else if(i1.eq.1) then
        if(i2.eq.3) then
          imark=8  ! fxyyyzz
        else
          imark=9  ! fxyyzzz
        end if
      else
        imark=10    ! fyyyzzz
      end if

      !  isum=7 not possible
    else if(isum.eq.-7) then
      !  7th derivative
      if(i1.eq.3) then
        if(i2.eq.3) then
          imark=2  ! fxxxyyyz
        else if(i2.eq.2) then
          imark=3  ! fxxxyyzz
        else
          imark=4  ! fxxxyzzz
        end if
      else if(i1.eq.2) then
        if(i2.eq.3) then
          imark=5  ! fxxyyyzz
        else
          imark=6  ! fxxyyzzz
        end if
      else
        imark=7     ! fxyyyzzz
      end if

      !  isum=8 not possible
    else if(isum.eq.-8) then
      !  8th derivative
      if(i3.eq.2) then 
        imark=2  ! fxxxyyyzz
      else if(i2.eq.2) then
        imark=3  ! fxxxyyzzz
      else
        imark=4  ! fxxyyyzzz
      end if

      !  isum=9 not possible
      !  isum=-9 --> fxxxyyyzzz

    end if

    if(abs(isum).gt.2) then
      do iii=2,10
        if(iii.eq.imark) then
          ict(iii)=1
        else
          ict(iii)=0
        end if
      end do
    end if

  end subroutine ezmake_ict3

end module EZspline_obj

!=======================================================================

module EZspline
  use precision_mod, only: fp

  interface EZspline_init
    !
    ! Initialize and allocate memory. BCS1,2,3 determine the type of boundary
    ! conditions on either end of the x1,2,3 grids: eg (/0, 0/) for not-a-knot
    ! on the left and (/-1, -1/) periodic. Other BCs such as imposed slope or
    ! second derivative can also be applied by setting '1' or '2' respectively
    ! on either side. For instance (/1, 2/) for 1st derivative on the left and
    ! 2nd derivative on the right. The value of the 1st/2nd derivative must be
    ! set explicitely set through the bcval1,2,3min and bcval1,2,3max arrays.
    !
    subroutine EZspline_init3(spline_o, n1, n2, n3, BCS1, BCS2, BCS3, ier)
      use EZspline_obj
      type(EZspline3) :: spline_o
      integer, intent(in) :: n1, n2, n3
      integer, intent(in) :: BCS1(2), BCS2(2), BCS3(2)
      integer, intent(out) :: ier
    end subroutine EZspline_init3

    subroutine EZspline_init2(spline_o, n1, n2, BCS1, BCS2, ier)
      use EZspline_obj
      type(EZspline2) :: spline_o
      integer, intent(in) :: n1, n2
      integer, intent(in) :: BCS1(2), BCS2(2)
      integer, intent(out) :: ier
    end subroutine EZspline_init2

    subroutine EZspline_init1(spline_o, n1, BCS1, ier)
      use EZspline_obj
      type(EZspline1) :: spline_o
      integer, intent(in) :: n1
      integer, intent(in) :: BCS1(2)
      integer, intent(out) :: ier
    end subroutine EZspline_init1

  end interface


  interface EZlinear_init
    !
    ! Initialize and allocate memory for piecewise LINEAR interpolation
    ! object.  This is C0.  For C2 cubic spline or C1 Akima Hermite spline,
    ! please see EZspline_init.
    !
    ! No boundary conditions are needed for piecewise linear interpolation
    !
    subroutine EZlinear_init3(spline_o, n1, n2, n3, ier)
      use EZspline_obj
      type(EZspline3) spline_o
      integer, intent(in) :: n1, n2, n3
      integer, intent(out) :: ier
    end subroutine EZlinear_init3

    subroutine EZlinear_init2(spline_o, n1, n2, ier)
      use EZspline_obj
      type(EZspline2) spline_o
      integer, intent(in) :: n1, n2
      integer, intent(out) :: ier
    end subroutine EZlinear_init2

    subroutine EZlinear_init1(spline_o, n1, ier)
      use EZspline_obj
      type(EZspline1) spline_o
      integer, intent(in) :: n1
      integer, intent(out) :: ier
    end subroutine EZlinear_init1

  end interface


  interface EZhybrid_init
    !
    ! Initialize and allocate memory for hybrid interpolation object:
    ! dimensionality > 1 only.  Interpolation method is specified separately
    ! for each dimension.  At present, Akima Hermite and Spline interpolation
    ! cannot be mixed.
    !
    ! Boundary condition arguments are optional.  They are appropriate only
    ! for the end points of dimensions for which Hermite or Spline cubic
    ! interpolation is used.
    !
    ! hspline(...) specifies the interpolation method for each dimension,
    ! according to the code: -1 for step function, 0 for piecewise linear,
    ! 1 for Akima Hermite, 2 for cubic Spline.
    !
    subroutine EZhybrid_init3(spline_o, n1, n2, n3, hspline, ier, &
         BCS1, BCS2, BCS3)
      use EZspline_obj
      type(EZspline3) spline_o
      integer, intent(in) :: n1, n2, n3
      integer, intent(in) :: hspline(3)
      integer, intent(out) :: ier
      integer, intent(in), OPTIONAL :: BCS1(2), BCS2(2), BCS3(2)
    end subroutine EZhybrid_init3

    subroutine EZhybrid_init2(spline_o, n1, n2, hspline, ier, &
         BCS1, BCS2)
      use EZspline_obj
      type(EZspline2) spline_o
      integer, intent(in) :: n1, n2
      integer, intent(in) :: hspline(2)
      integer, intent(out) :: ier
      integer, intent(in), OPTIONAL :: BCS1(2), BCS2(2)
    end subroutine EZhybrid_init2

  end interface

  interface EZspline_free
    !
    ! Reset and free the memory. This method must be called to avoid
    ! memory leaks after all interpolations have been computed.
    !
    subroutine EZspline_free3(spline_o, ier)
      use EZspline_obj
      type(EZspline3) spline_o
      integer, intent(out) :: ier
    end subroutine EZspline_free3

    subroutine EZspline_free2(spline_o, ier)
      use EZspline_obj
      type(EZspline2) spline_o
      integer, intent(out) :: ier
    end subroutine EZspline_free2

    subroutine EZspline_free1(spline_o, ier)
      use EZspline_obj
      type(EZspline1) spline_o
      integer, intent(out) :: ier
    end subroutine EZspline_free1

  end interface

  interface EZspline_setup
    !
    ! Compute the cubic spline coefficients. Note: the grid and the
    ! boundary conditions should be properly set prior to this call.
    !
    ! NEW optional argument: exact_dim=TRUE to requre f dimensions to
    ! match higher dimensions of spline_o%fspl exactly; default or FALSE
    ! means f dimensions can match or exceed dimensions of spline_o%fspl.
    !
    ! array arguments are now declared with f90 style dimensioning; the
    ! module interface must be used (if not feasible see ezspline_setupx.f90).
    ! 
    subroutine EZspline_setup3(spline_o, f, ier, exact_dim)
      use EZspline_obj
      type(EZspline3) spline_o
      real(fp), dimension(:,:,:), intent(in) :: f
      integer, intent(out) :: ier
      logical, intent(in), optional :: exact_dim
    end subroutine EZspline_setup3

    subroutine EZspline_setup2(spline_o, f, ier, exact_dim)
      use EZspline_obj
      type(EZspline2) spline_o
      real(fp), dimension(:,:), intent(in) :: f
      integer, intent(out) :: ier
      logical, intent(in), optional :: exact_dim
    end subroutine EZspline_setup2

    subroutine EZspline_setup1(spline_o, f, ier, exact_dim)
      use EZspline_obj
      type(EZspline1) spline_o
      real(fp), dimension(:), intent(in) :: f
      integer, intent(out) :: ier
      logical, intent(in), optional :: exact_dim
    end subroutine EZspline_setup1

  end interface



  interface EZspline_interp
    !
    ! Interpolation at grid point(s) p1, [p2, [p3]]. Result is returned in
    ! f. Interpolation can be sought at a single point (p1, [p2, [p3]] are
    ! scalars), on an unordered list of points (p1, [p2, [p3]] have dimension
    ! k), or on a structured grid (p1, [p2, [p3]] have dimension k1, [k2, [k3]]
    ! respectively).
    !
    subroutine EZspline_interp3(spline_o, p1, p2, p3, f, ier)
      ! single point evaluation
      use EZspline_obj
      type(EZspline3) spline_o
      real(fp) :: p1, p2, p3
      real(fp) f
      integer, intent(out) :: ier
    end subroutine EZspline_interp3

    subroutine EZspline_interp3_array(spline_o, k1, k2, k3, p1, p2, p3, f, ier)
      use EZspline_obj
      type(EZspline3) spline_o
      integer :: k1, k2, k3
      real(fp), intent(in) :: p1(k1), p2(k2), p3(k3)
      real(fp), intent(out):: f(k1,k2,*)
      integer, intent(out) :: ier
    end subroutine EZspline_interp3_array

    subroutine EZspline_interp3_cloud(spline_o, k, p1, p2, p3, f, ier)
      use EZspline_obj
      type(EZspline3) spline_o
      integer, intent(in) :: k
      real(fp), intent(in) :: p1(k), p2(k), p3(k)
      real(fp), intent(out):: f(k)
      integer, intent(out) :: ier
    end subroutine EZspline_interp3_cloud

    subroutine EZspline_interp2(spline_o, p1, p2, f, ier)
      use EZspline_obj
      type(EZspline2) spline_o
      real(fp) :: p1, p2
      real(fp) f
      integer, intent(out) :: ier
    end subroutine EZspline_interp2

    subroutine EZspline_interp2_array(spline_o, k1, k2, p1, p2, f, ier)
      use EZspline_obj
      type(EZspline2) spline_o
      integer :: k1, k2
      real(fp), intent(in) :: p1(k1), p2(k2)
      real(fp), intent(out):: f(k1,*)
      integer, intent(out) :: ier
    end subroutine EZspline_interp2_array

    subroutine EZspline_interp2_cloud(spline_o, k, p1, p2, f, ier)
      use EZspline_obj
      type(EZspline2) spline_o
      integer, intent(in) :: k
      real(fp), intent(in) :: p1(k), p2(k)
      real(fp), intent(out):: f(k)
      integer, intent(out) :: ier
    end subroutine EZspline_interp2_cloud

    subroutine EZspline_interp1(spline_o, p1, f, ier)
      use EZspline_obj
      type(EZspline1) spline_o
      real(fp) :: p1
      real(fp) f
      integer, intent(out) :: ier
    end subroutine EZspline_interp1

    subroutine EZspline_interp1_array(spline_o, k1, p1, f, ier)
      use EZspline_obj
      type(EZspline1) spline_o
      integer :: k1
      real(fp), intent(in) :: p1(k1)
      real(fp), intent(out):: f(k1)
      integer, intent(out) :: ier
    end subroutine EZspline_interp1_array

  end interface

  interface EZspline_derivative
    !
    ! Evaluate the spline/Akima Hermite derivative of order
    ! d^{i1} d^{i2} d^{i3} f / d x1^{i1} d x2^{i2} d x2^{i2}
    ! at p1, [p2, [p3]]. The sum of i1+[i2+[i3]] should be <=2 for spline, or
    ! <=1 for Akima Hermite or Piecewise Linear. 
    ! Result is return in f. The evaluation can
    ! be sought at a single point (p1, [p2, [p3]] are scalars), on
    ! an unordered list of points (p1, [p2, [p3]] have dimension k), or
    ! on a structured grid (p1, [p2, [p3]] have dimension k1, [k2, [k3]]
    ! respectively).
    !
    !
    subroutine EZspline_derivative3(spline_o, i1, i2, i3, p1, p2, p3, f, ier)
      use EZspline_obj
      type(EZspline3) spline_o
      integer, intent(in) :: i1, i2, i3
      real(fp), intent(in) :: p1, p2, p3
      real(fp), intent(out) :: f
      integer, intent(out) :: ier
    end subroutine EZspline_derivative3

    subroutine EZspline_derivative3_array(spline_o, i1, i2, i3, &
         & k1, k2, k3, p1, p2, p3, f, ier)
      use EZspline_obj
      type(EZspline3) spline_o
      integer, intent(in) :: i1, i2, i3, k1, k2, k3
      real(fp), intent(in) :: p1(k1), p2(k2), p3(k3)
      real(fp), intent(out) :: f(k1, k2, *)
      integer, intent(out) :: ier
    end subroutine EZspline_derivative3_array

    subroutine EZspline_derivative3_cloud(spline_o, i1, i2, i3, &
         & k, p1, p2, p3, f, ier)
      use EZspline_obj
      type(EZspline3) spline_o
      integer, intent(in) :: i1, i2, i3, k
      real(fp), intent(in) :: p1(k), p2(k), p3(k)
      real(fp), intent(out) :: f(k)
      integer, intent(out) :: ier
    end subroutine EZspline_derivative3_cloud

    subroutine EZspline_derivative2(spline_o, i1, i2, p1, p2, f, ier)
      use EZspline_obj
      type(EZspline2) spline_o
      integer, intent(in) :: i1, i2
      real(fp), intent(in) :: p1, p2
      real(fp), intent(out) :: f
      integer, intent(out) :: ier
    end subroutine EZspline_derivative2

    subroutine EZspline_derivative2_array(spline_o, i1, i2, k1, k2, p1, p2, f, ier)
      use EZspline_obj
      type(EZspline2) spline_o
      integer, intent(in) :: i1, i2, k1, k2
      real(fp), intent(in) :: p1(k1), p2(k2)
      real(fp), intent(out) :: f(k1,k2)
      integer, intent(out) :: ier
    end subroutine EZspline_derivative2_array

    subroutine EZspline_derivative2_cloud(spline_o, i1, i2, k, p1, p2, f, ier)
      use EZspline_obj
      type(EZspline2) spline_o
      integer, intent(in) :: i1, i2, k
      real(fp), intent(in) :: p1(k), p2(k)
      real(fp), intent(out) :: f(k)
      integer, intent(out) :: ier
    end subroutine EZspline_derivative2_cloud

    subroutine EZspline_derivative1(spline_o, i1, p1, f, ier)
      use EZspline_obj
      type(EZspline1) spline_o
      integer, intent(in) :: i1
      real(fp), intent(in) :: p1
      real(fp), intent(out) :: f
      integer, intent(out) :: ier
    end subroutine EZspline_derivative1

    subroutine EZspline_derivative1_array(spline_o, i1, k1, p1, f, ier)
      use EZspline_obj
      type(EZspline1) spline_o
      integer, intent(in) :: i1, k1
      real(fp), intent(in) :: p1(k1)
      real(fp), intent(out) :: f(k1)
      integer, intent(out) :: ier
    end subroutine EZspline_derivative1_array

  end interface

  interface EZspline_gradient
    !
    ! Return the gradient in df. When the dimensionality is 1 then
    ! df is df/dx. In more than one dimension, df has rank >=1 with
    ! the last index yielding df/dx, df/dy ... etc. Subsequent indices
    ! reflect the node positions when cloud or array evaluation is
    ! sought, as in df(k1, k2, k3, 1) for df/dx at x(k1), y(k2), z(k3).
    !
    subroutine EZspline_gradient3(spline_o, p1, p2, p3, df, ier)
      use EZspline_obj
      type(EZspline3) spline_o
      real(fp), intent(in) :: p1, p2, p3
      real(fp), intent(out) :: df(3)
      integer, intent(out) :: ier
    end subroutine EZspline_gradient3

    subroutine EZspline_gradient3_array(spline_o, k1, k2, k3, &
         & p1, p2, p3, df, ier)
      use EZspline_obj
      type(EZspline3) spline_o
      integer, intent(in) :: k1, k2, k3
      real(fp), intent(in) :: p1(k1), p2(k2), p3(k3)
      real(fp), intent(out) :: df(k1,k2,k3,3)
      integer, intent(out) :: ier
    end subroutine EZspline_gradient3_array

    subroutine EZspline_gradient3_cloud(spline_o, k, &
         & p1, p2, p3, df, ier)
      use EZspline_obj
      type(EZspline3) spline_o
      integer, intent(in) :: k
      real(fp), intent(in) :: p1(k), p2(k), p3(k)
      real(fp), intent(out) :: df(k,3)
      integer, intent(out) :: ier
    end subroutine EZspline_gradient3_cloud

    subroutine EZspline_gradient2(spline_o, p1, p2, df, ier)
      use EZspline_obj
      type(EZspline2) spline_o
      real(fp), intent(in) :: p1, p2
      real(fp), intent(out) :: df(2)
      integer, intent(out) :: ier
    end subroutine EZspline_gradient2

    subroutine EZspline_gradient2_array(spline_o, k1, k2, &
         & p1, p2, df, ier)
      use EZspline_obj
      type(EZspline2) spline_o
      integer, intent(in)  :: k1, k2
      real(fp), intent(in) :: p1(k1), p2(k2)
      real(fp), intent(out) :: df(k1, k2, 2)
      integer, intent(out) :: ier
    end subroutine EZspline_gradient2_array

    subroutine EZspline_gradient2_cloud(spline_o, k, &
         & p1, p2, df, ier)
      use EZspline_obj
      type(EZspline2) spline_o
      integer, intent(in)  :: k
      real(fp), intent(in) :: p1(k), p2(k)
      real(fp), intent(out) :: df(k, 2)
      integer, intent(out) :: ier
    end subroutine EZspline_gradient2_cloud

    subroutine EZspline_gradient1(spline_o, p1, df, ier)
      use EZspline_obj
      type(EZspline1) spline_o
      real(fp), intent(in) :: p1
      real(fp), intent(out) :: df
      integer, intent(out) :: ier
    end subroutine EZspline_gradient1

    subroutine EZspline_gradient1_array(spline_o, k1, p1, df, ier)
      use EZspline_obj
      type(EZspline1) spline_o
      integer, intent(in) :: k1
      real(fp), intent(in) :: p1(k1)
      real(fp), intent(out) :: df(k1)
      integer, intent(out) :: ier
    end subroutine EZspline_gradient1_array

  end interface

  interface EZspline_isInDomain
    !
    ! Return error code if position (p1, [p2, [p3]]) is outside domain.
    ! The evaluation can be sought at a single point (p1, [p2, [p3]]
    ! are scalars), on an unordered list of points (p1, [p2, [p3]] have
    ! dimension k), or on a structured grid (p1, [p2, [p3]] have dimension
    ! k1, [k2, [k3]] respectively).
    !
    subroutine EZspline_isInDomain3(spline_o, p1, p2, p3, ier)
      use EZspline_obj
      type(EZspline3) :: spline_o
      real(fp), intent(in) :: p1, p2, p3
      integer, intent(out) :: ier
    end subroutine EZspline_isInDomain3

    subroutine EZspline_isInDomain3_array(spline_o, k1, k2, k3, p1, p2, p3, ier)
      use EZspline_obj
      type(EZspline3) :: spline_o
      integer, intent(in) :: k1, k2, k3
      real(fp), intent(in) :: p1(k1), p2(k2), p3(k3)
      integer, intent(out) :: ier
    end subroutine EZspline_isInDomain3_array

    subroutine EZspline_isInDomain3_cloud(spline_o, k, p1, p2, p3, ier)
      use EZspline_obj
      type(EZspline3) :: spline_o
      integer, intent(in) :: k
      real(fp), intent(in) :: p1(k), p2(k), p3(k)
      integer, intent(out) :: ier
    end subroutine EZspline_isInDomain3_cloud

    subroutine EZspline_isInDomain2(spline_o, p1, p2, ier)
      use EZspline_obj
      type(EZspline2) :: spline_o
      real(fp), intent(in) :: p1, p2
      integer, intent(out) :: ier
    end subroutine EZspline_isInDomain2

    subroutine EZspline_isInDomain2_array(spline_o, k1, k2, p1, p2, ier)
      use EZspline_obj
      type(EZspline2) :: spline_o
      integer, intent(in) :: k1, k2
      real(fp), intent(in) :: p1(k1), p2(k2)
      integer, intent(out) :: ier
    end subroutine EZspline_isInDomain2_array

    subroutine EZspline_isInDomain2_cloud(spline_o, k, p1, p2, ier)
      use EZspline_obj
      type(EZspline2) :: spline_o
      integer, intent(in) :: k
      real(fp), intent(in) :: p1(k), p2(k)
      integer, intent(out) :: ier
    end subroutine EZspline_isInDomain2_cloud

    subroutine EZspline_isInDomain1(spline_o, p1, ier)
      use EZspline_obj
      type(EZspline1) :: spline_o
      real(fp), intent(in) :: p1
      integer, intent(out) :: ier
    end subroutine EZspline_isInDomain1

    subroutine EZspline_isInDomain1_array(spline_o, k1, p1, ier)
      use EZspline_obj
      type(EZspline1) :: spline_o
      integer, intent(in) :: k1
      real(fp), intent(in) :: p1(k1)
      integer, intent(out) :: ier
    end subroutine EZspline_isInDomain1_array

  end interface

  interface EZspline_isGridRegular
    !
    ! Return error code if grid is not regular (not strictly  increasing).
    ! The evaluation can be sought at a single point (p1, [p2, [p3]]
    ! are scalars), on an unordered list of points (p1, [p2, [p3]] have
    ! dimension k), or on a structured grid (p1, [p2, [p3]] have dimension
    ! k1, [k2, [k3]] respectively).
    !
    subroutine EZspline_isGridRegular3(spline_o, ier)
      use EZspline_obj
      type(EZspline3) :: spline_o
      integer, intent(out) :: ier
    end subroutine EZspline_isGridRegular3

    subroutine EZspline_isGridRegular2(spline_o, ier)
      use EZspline_obj
      type(EZspline2) :: spline_o
      integer, intent(out) :: ier
    end subroutine EZspline_isGridRegular2

    subroutine EZspline_isGridRegular1(spline_o, ier)
      use EZspline_obj
      type(EZspline1) :: spline_o
      integer, intent(out) :: ier
    end subroutine EZspline_isGridRegular1

  end interface

#ifdef _EZCDF
  interface EZspline_save
    !
    ! Save spline/Akima Hermite/Linear object in netcdf file 'filename'. Use
    ! EZspline_load to load spline/Akima Hermite/Linear object from netcdf 
    ! file.
    !
    ! Mod DMC March 2006 -- optionally, by giving each spline a name, 
    ! multiple splines can be saved in a single file.  Also, by specifying
    ! "fullsave=.TRUE." the spline coefficients can be saved as well in the
    ! file, so that they do not have to be recomputed at ezspline_load time.
    !
    ! When creating a single file with multiple spline objects, it is the
    ! user's responsibilitly to make sure that a different name is used for
    ! each spline that is to be saved.  Names can consist of upper or lower
    ! case letters, numerals, and "_", but must not start with a numeral.
    ! Imbedded blanks are not allowed, and the length of the name must be
    ! no more than 20 characters long-- to allow the names of the spline
    ! object elements to be appended.
    !
    subroutine EZspline_save3(spline_o, filename, ier, &
         spl_name,fullsave)
      use EZspline_obj
      use EZcdf
      type(EZspline3) :: spline_o
      character(len=*) :: filename
      integer, intent(out) :: ier

      character(len=*), intent(in),optional :: spl_name
      logical, intent(in), optional :: fullsave
    end subroutine EZspline_save3

    subroutine EZspline_save2(spline_o, filename, ier, &
         spl_name,fullsave)
      use EZspline_obj
      use EZcdf
      type(EZspline2) :: spline_o
      character(len=*) :: filename
      integer, intent(out) :: ier

      character(len=*), intent(in),optional :: spl_name
      logical, intent(in), optional :: fullsave
    end subroutine EZspline_save2

    subroutine EZspline_save1(spline_o, filename, ier, &
         spl_name,fullsave)
      use EZspline_obj
      use EZcdf
      type(EZspline1) :: spline_o
      character(len=*) :: filename
      integer, intent(out) :: ier

      character(len=*), intent(in),optional :: spl_name
      logical, intent(in), optional :: fullsave
    end subroutine EZspline_save1

  end interface
#endif

#ifdef _EZCDF
  interface EZspline_load
    !
    ! Load spline/Akima Hermite object from netcdf file 'filename'. Use
    ! EZspline_save to save spline/Akima Hermite/Linear object in netcdf file.
    !
    ! MOD DMC March 2006-- a single NetCDF file can now contain multiple
    ! *named* spline objects.  If accessing such an object, the name must
    ! be supplied via the optional argument "spl_name".  If this is omitted,
    ! the default is to read the contents of the file which is presumed to
    ! consist of but a single spline object.

    ! Each call still opens and closes the file; if users find this 
    ! inefficient, an improvement to the control interface may be built.

    subroutine EZspline_load3(spline_o, filename, ier, spl_name)
      use EZspline_obj
      use EZcdf
      type(EZspline3) :: spline_o
      character(len=*) :: filename
      integer, intent(out) :: ier

      character(len=*), intent(in),optional :: spl_name
    end subroutine EZspline_load3

    subroutine EZspline_load2(spline_o, filename, ier, spl_name)
      use EZspline_obj
      use EZcdf
      type(EZspline2) :: spline_o
      character(len=*) :: filename
      integer, intent(out) :: ier

      character(len=*), intent(in),optional :: spl_name
    end subroutine EZspline_load2

    subroutine EZspline_load1(spline_o, filename, ier, spl_name)
      use EZspline_obj
      use EZcdf
      type(EZspline1) :: spline_o
      character(len=*) :: filename
      integer, intent(out) :: ier

      character(len=*), intent(in),optional :: spl_name
    end subroutine EZspline_load1

  end interface
#endif

  interface EZspline_modulo
    !
    ! Map argument to (x1,2,3min, x1,2,3max) interval. This is useful to avoid
    ! an out-of-grid error when periodic boundary conditions are applied.
    ! The evaluation can be sought at a single point (p1, [p2, [p3]]
    ! are scalars), on an unordered list of points (p1, [p2, [p3]] have
    ! dimension k), or on a structured grid (p1, [p2, [p3]] have dimension
    ! k1, [k2, [k3]] respectively).
    !
    subroutine EZspline_modulo3(spline_o, p1, p2, p3, ier)
      use EZspline_obj
      type(EZspline3) spline_o
      real(fp) :: p1, p2, p3
      integer, intent(out) :: ier
    end subroutine EZspline_modulo3

    subroutine EZspline_modulo_array3(spline_o, k1, k2, k3, p1, p2, p3, ier)
      use EZspline_obj
      type(EZspline3) spline_o
      integer, intent(in) :: k1, k2, k3
      real(fp) :: p1(k1), p2(k2), p3(k3)
      integer, intent(out) :: ier
    end subroutine EZspline_modulo_array3

    subroutine EZspline_modulo_cloud3(spline_o, k, p1, p2, p3, ier)
      use EZspline_obj
      type(EZspline3) spline_o
      integer, intent(in) :: k
      real(fp) :: p1(k), p2(k), p3(k)
      integer, intent(out) :: ier
    end subroutine EZspline_modulo_cloud3

    subroutine EZspline_modulo2(spline_o, p1, p2, ier)
      use EZspline_obj
      type(EZspline2) spline_o
      real(fp) :: p1, p2
      integer, intent(out) :: ier
    end subroutine EZspline_modulo2

    subroutine EZspline_modulo_array2(spline_o, k1, k2, p1, p2, ier)
      use EZspline_obj
      type(EZspline2) spline_o
      integer, intent(in) :: k1, k2
      real(fp) :: p1(k1), p2(k2)
      integer, intent(out) :: ier
    end subroutine EZspline_modulo_array2

    subroutine EZspline_modulo_cloud2(spline_o, k, p1, p2, ier)
      use EZspline_obj
      type(EZspline2) spline_o
      integer, intent(in) :: k
      real(fp) :: p1(k), p2(k)
      integer, intent(out) :: ier
    end subroutine EZspline_modulo_cloud2

    subroutine EZspline_modulo1(spline_o, p1, ier)
      use EZspline_obj
      type(EZspline1) spline_o
      real(fp) :: p1
      integer, intent(out) :: ier
    end subroutine EZspline_modulo1

    subroutine EZspline_modulo_array1(spline_o, k1, p1, ier)
      use EZspline_obj
      type(EZspline3) spline_o
      integer, intent(in) :: k1
      real(fp) :: p1(k1)
      integer, intent(out) :: ier
    end subroutine EZspline_modulo_array1

  end interface

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! EZspline object-less methods (first argument is NOT an EZspline1,2,3 type).
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef _EZCDF
  interface EZspline_2NetCDF
    !
    ! Save data in netCDF file 'filename'. To save an EZspline1,2,3 object use
    ! EZspline_save method.
    !
    subroutine EZspline_2NetCDF_array3(n1, n2, n3, x1, x2, x3, f, filename, ier)
      use EZspline_obj
      use EZcdf
      implicit none
      integer, intent(in) :: n1, n2, n3
      real(fp), intent(in) :: x1(:), x2(:), x3(:), f(:, :, :)
      character(len=*), intent(in) :: filename
      integer, intent(out) :: ier
    end subroutine EZspline_2NetCDF_array3

    subroutine EZspline_2NetCDF_array2(n1, n2, x1, x2, f, filename, ier)
      use EZspline_obj
      use EZcdf
      implicit none
      integer, intent(in) :: n1, n2
      real(fp), intent(in) ::  x1(:), x2(:), f(:,:)
      character(len=*), intent(in) :: filename
      integer, intent(out) :: ier
    end subroutine EZspline_2NetCDF_array2

    subroutine EZspline_2NetCDF1(n1, x1, f, filename, ier)
      use EZspline_obj
      use EZcdf
      implicit none
      integer, intent(in) :: n1
      real(fp), intent(in) :: x1(:), f(:)
      character(len=*), intent(in) :: filename
      integer, intent(out) :: ier
    end subroutine EZspline_2NetCDF1

    subroutine EZspline_2NetCDF_cloud3(n, x1, x2, x3, f, filename, ier)
      use EZspline_obj
      use EZcdf
      implicit none
      integer, intent(in) :: n
      real(fp), intent(in) :: x1(:), x2(:), x3(:), f(:)
      character(len=*), intent(in) :: filename
      integer, intent(out) :: ier
    end subroutine EZspline_2NetCDF_cloud3

    subroutine EZspline_2NetCDF_cloud2(n, x1, x2, f, filename, ier)
      use EZspline_obj
      use EZcdf
      implicit none
      integer, intent(in) :: n
      real(fp), intent(in) :: x1(:), x2(:), f(:)
      character(len=*), intent(in) :: filename
      integer, intent(out) :: ier
    end subroutine EZspline_2NetCDF_cloud2

  end interface
#endif

contains


  subroutine EZspline_error(ier)
    !
    ! Error handling routine. Maps error ier code to a meaningful message.
    ! Note: does not abort nor stop if ier/=0.
    !
    implicit none
    integer, intent(in) :: ier

    if(ier == 0) return
    write(6,*) '**EZspline** ERROR/WARNING #', ier,' occurred'

    select case(ier)
    case(1)
      write(6,*) '**EZspline** allocation error'
    case(2)
      write(6,*) '**EZspline** wrong BCS1 code'
    case(3)
      write(6,*) '**EZspline** wrong BCS2 code'
    case(4)
      write(6,*) '**EZspline** wrong BCS3 code'
    case(5)
      write(6,*) '**EZspline** Que??'
    case(6)
      write(6,*) '**EZspline** out of interval p1 < min(x1)'
    case(7)
      write(6,*) '**EZspline** out of interval p1 > max(x1)'
    case(8)
      write(6,*) '**EZspline** out of interval p2 < min(x2)'
    case(9)
      write(6,*) '**EZspline** out of interval p2 > max(x2)'
    case(10)
      write(6,*) '**EZspline** out of interval p3 < min(x3)'
    case(11)
      write(6,*) '**EZspline** out of interval p3 > max(x3)'
    case(12)
      write(6,*) '**EZspline** negative derivative order'
    case(13)
      write(6,*) '**EZspline** derivative order too high'
    case(14)
      write(6,*) '**EZspline** x1 grid is not strictly increasing'
    case(15)
      write(6,*) '**EZspline** x2 grid is not strictly increasing'
    case(16)
      write(6,*) '**EZspline** x3 grid is not strictly increasing'
    case(17)
      write(6,*) '**EZspline** could not save spline object in file '
    case(18)
      write(6,*) '**EZspline** memory allocation failure in coefficient setup'

    case(20)
      write(6,*) '**EZspline** attempt to load spline object with wrong rank.'
    case(21)
      write(6,*) '**EZspline** could not load spline object from file '
    case(22)
      write(6,*) '**EZspline** loaded spline object from file but failed at coefficient set-up'
    case(23)
      write(6,*) '**EZspline** failed to free spline object'
    case(24)
      write(6,*) '**EZspline** 2nd order derivative not supported for Akima-Hermite (isHermite=1)'
    case(25)
      write(6,*) '**EZspline** not supported for Akima-Hermite (isHermite=1)'
    case(26)
      write(6,*) '**EZspline** memory allocation error in EZspline_interp'
    case(27)
      write(6,*) '**EZspline** an error ocurred in genxpkg'
    case(28)
      write(6,*) '**EZspline** memory allocation failure in ezspline_interp'
    case(29)
      write(6,*) '**EZspline** memory deallocation failure in ezspline_interp'
    case(30)
      write(6,*) '**EZspline** memory allocation error in EZspline_gradient'
    case(31)
      write(6,*) '**EZspline** memory deallocation error in EZspline_gradient'
    case(32)
      write(6,*) '**EZspline** memory allocation error in EZspline_derivative'
    case(33)
      write(6,*) '**EZspline** memory deallocation error in EZspline_derivative'
    case(34)
      write(6,*) '**EZspline** could not open netCDF file in EZspline_2netcdf'
    case(35)
      write(6,*) '**EZspline** could not write into netCDF file in EZspline_2netcdf'
    case(36)
      write(6,*) '**EZspline** could not read from netCDF file in EZspline_2netcdf'
    case(37)
      write(6,*) '**EZspline** could not close netCDF file in EZspline_2netcdf'
    case(38)
      write(6,*) '**EZspline** could not define variable (cdfDefVar) in EZspline_2netcdf'
    case(39)
      write(6,*) '**EZspline** could not open netCDF file in EZspline_save'
    case(40)
      write(6,*) '**EZspline** could not write into netCDF file in EZspline_save'
    case(41)
      write(6,*) '**EZspline** could not close netCDF file in EZspline_save'
    case(42)
      write(6,*) '**EZspline** could not define variable (cdfDefVar) in EZspline_save'
    case(43)
      write(6,*) '**EZspline** could not open netCDF file in EZspline_load'
    case(44)
      write(6,*) '**EZspline** could not read from netCDF file in EZspline_load'
    case(45)
      write(6,*) '**EZspline** could not close netCDF file in EZspline_load'
    case(46)
      write(6,*) '**EZspline** 2nd order derivative not supported for Piecewise Linear Interpolation (isLinear=1)'
    case(47)
      write(6,*) '**EZspline** not supported for Piecewise Linear Interpolation (isLinear=1)'

    case(50)
      write(6,*) '**EZspline** ezspline_save (optional) spline name is blank.'
    case(51)
      write(6,*) '**EZspline** ezspline_save (optional) spline name too long (max 20 characters).'
    case(52)
      write(6,*) '**EZspline** ezspline_save (optional) spline name contains'
      write(6,*) '             imbedded blanks or other illegal characters.'
    case(53)
      write(6,*) '**EZspline** attempt to write named spline object to NetCDF'
      write(6,*) '             file with change of dimensionality or data type.'

    case(54)
      write(6,*) '**EZspline** hybrid interpolation specification not in range -1:2'
      write(6,*) '             error in EZhybrid_init.'
    case(55)
      write(6,*) '**EZspline** hybrid interpolation cannot mix Hermite and Spline interpolation.'
      write(6,*) '             hspline(i)=1 and hspline(j)=2 in EZhybrid_init.'
    case(56)
      write(6,*) '**EZspline** non-default boundary condition unavailable: zonal or piecewise linear dimension.'
      write(6,*) '             in EZhybrid_init.'

    case(57)
      write(6,*) '**EZspline** dimension of "f" smaller than corresponding "fspl"'
      write(6,*) '             dimension in "spline_o".'
    case(58)
      write(6,*) '**EZspline** dimension of "f" larger than corresponding "fspl"'
      write(6,*) '             dimension in "spline_o".'

    case(90)
      write(6,*) '**EZspline** an error occurred after attempting to evaluate the'
      write(6,*) '             Hermite polynomials'
    case(91)
      write(6,*) '**EZspline** an error occurred after attempting to set up the'
      write(6,*) '             Hermite polynomial coefficients'
    case(92)
      write(6,*) '**EZspline** warning in EZspline_load. Looks like saved object '
      write(6,*) '             was not properly set-up (isReady=0).'
    case(93)
      write(6,*) '**EZspline** warning in EZspline_save. Looks like saved object '
      write(6,*) '             was not properly set-up (isReady=0).'
    case(94)
      write(6,*) '**EZspline** an error occurred in EZspline_interp. Did you forget'
      write(6,*) '             to set up the cubic spline coefficients by calling'
      write(6,*) '             call EZspline_setup(spline_o, f, ier)'
      write(6,*) '             ?'
    case(95)
      write(6,*) '**EZspline** some error occurred in EZspline_gradient'
    case(96)
      write(6,*) '**EZspline** some error occurred in EZspline_derivative'
    case(97)
      write(6,*) '**EZspline** some error occurred in EZspline_interp apparently'
      write(6,*) '             related to a PSPLINE routine. Check if argument is '
      write(6,*) '             outside interpolation domain by calling'
      write(6,*) '             call EZspline_isInDomain(spline_o, [[k1, k2, k3,] .OR. k,] p1, p2, p3, ier ,ier)'
      write(6,*) '             call EZspline_error(ier)'
    case(98)
      write(6,*) '**EZspline** error occurred in EZspline_setup'
      write(6,*) '  if no other explanation-- ezspline_init call never made.'
    case(99)
      write(6,*) '**EZspline** some error occurred in EZspline_init, EZhybrid_init,  or EZlinear_init'
    case(100)
      write(6,*) '**EZSPLINE** EZspline_init, EZhybrid_init,  or EZlinear_init -- object already allocated.'
    case(101)
      write(6,*) '**EZSPLINE** object was never allocated.'
    case default
      write(6,*) '**EZspline** '
    end select

    return
  end subroutine EZspline_error

end module EZspline
