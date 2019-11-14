module czspline3d
  use iso_c_binding
  use EZspline_type, only: EZspline3

  type(EZspline3), private :: spl

  public :: czspline_init3
  public :: czspline_set_axes3
  public :: czspline_set_ishermite3
  public :: czspline_set_bctypes3
  public :: czspline_set_bcvals3
  public :: czspline_setup3
  public :: czspline_interp3
  public :: czspline_interp3_cloud
  public :: czspline_interp3_array
  public :: czspline_free3
  public :: czspline_save3
  public :: czspline_load3
  public :: czspline_isindomain3
  public :: czspline_isindomain3_cloud
  public :: czspline_isindomain3_array
  public :: czspline_gradient3
  public :: czspline_gradient3_cloud
  public :: czspline_gradient3_array
  public :: czspline_derivative3
  public :: czspline_derivative3_cloud
  public :: czspline_derivative3_array
  public :: czspline_isgridregular3

contains

  ! initialization
  subroutine czspline_init3(n1, n2, n3, bcs1, bcs2, bcs3, ier) bind(C)
    use ezspline
    implicit none
    integer(kind=c_int),               intent(in)  :: n1
    integer(kind=c_int),               intent(in)  :: n2
    integer(kind=c_int),               intent(in)  :: n3
    integer(kind=c_int), dimension(*), intent(in)  :: bcs1
    integer(kind=c_int), dimension(*), intent(in)  :: bcs2
    integer(kind=c_int), dimension(*), intent(in)  :: bcs3
    integer(kind=c_int),               intent(out) :: ier
    call ezspline_init3(spl,n1,n2,n3,bcs1,bcs2,bcs3,ier)
    return
  end subroutine czspline_init3

  ! member setters
  subroutine czspline_set_axes3(n1, n2, n3, x1, x2, x3, ier) bind(C)
    use ezspline
    implicit none
    integer(kind=c_int),                intent(in)  :: n1
    integer(kind=c_int),                intent(in)  :: n2
    integer(kind=c_int),                intent(in)  :: n3
    real(kind=c_double), dimension(n1), intent(in)  :: x1
    real(kind=c_double), dimension(n2), intent(in)  :: x2
    real(kind=c_double), dimension(n3), intent(in)  :: x3
    integer(kind=c_int),                intent(out) :: ier
    ier = 0
    spl%x1 = x1
    spl%x2 = x2
    spl%x3 = x3
    return
  end subroutine czspline_set_axes3

  subroutine czspline_set_ishermite3(flag, ier) bind(C)
    use ezspline
    implicit none
    integer(kind=c_int), intent(in)  :: flag
    integer(kind=c_int), intent(out) :: ier
    ier = 0
    spl%isHermite = flag
    return
  end subroutine czspline_set_ishermite3

  ! boundary types
  subroutine czspline_set_bctypes3(bctype1, bctype2, bctype3, ier) bind(C)
    use ezspline
    implicit none
    integer(kind=c_int), dimension(2), intent(in)  :: bctype1
    integer(kind=c_int), dimension(2), intent(in)  :: bctype2
    integer(kind=c_int), dimension(2), intent(in)  :: bctype3
    integer(kind=c_int),               intent(out) :: ier
    ier = 0
    spl%ibctype1 = bctype1
    spl%ibctype2 = bctype2
    spl%ibctype3 = bctype3
    return
  end subroutine czspline_set_bctypes3

  ! boundary conditions bcval1min, bcval1max
  subroutine czspline_set_bcvals3(bcval1, bcval2, bcval3, ier) bind(C)
    use ezspline
    implicit none
    real(kind=c_double), dimension(2), intent(in)  :: bcval1
    real(kind=c_double), dimension(2), intent(in)  :: bcval2
    real(kind=c_double), dimension(2), intent(in)  :: bcval3
    integer(kind=c_int),               intent(out) :: ier
    ier = 0
    spl%bcval1min = bcval1(1)
    spl%bcval1max = bcval1(2)
    spl%bcval2min = bcval2(1)
    spl%bcval2max = bcval2(2)
    spl%bcval3min = bcval3(1)
    spl%bcval3max = bcval3(2)
    return
  end subroutine czspline_set_bcvals3

  ! compute spline coefficients
  subroutine czspline_setup3(n1, n2, n3, f, ier) bind(C)
    use ezspline
    implicit none
    integer(kind=c_int),                      intent(in)  :: n1
    integer(kind=c_int),                      intent(in)  :: n2
    integer(kind=c_int),                      intent(in)  :: n3
    real(kind=c_double), dimension(n1,n2,n3), intent(in)  :: f
    integer(kind=c_int),                      intent(out) :: ier
    call ezspline_setup3(spl,f,ier)
    return
  end subroutine czspline_setup3

  ! point interpolation
  subroutine czspline_interp3(y1, y2, y3, g, ier) bind(C)
    use ezspline
    implicit none
    real(kind=c_double), intent(in)  :: y1
    real(kind=c_double), intent(in)  :: y2
    real(kind=c_double), intent(in)  :: y3
    real(kind=c_double), intent(out) :: g
    integer(kind=c_int), intent(out) :: ier
    call ezspline_interp3(spl,y1,y2,y3,g,ier)
    return
  end subroutine czspline_interp3

  ! cloud interpolation
  subroutine czspline_interp3_cloud(m, y1, y2, y3, g, ier) bind(C)
    use ezspline
    implicit none
    integer(kind=c_int),               intent(in)  :: m
    real(kind=c_double), dimension(m), intent(in)  :: y1
    real(kind=c_double), dimension(m), intent(in)  :: y2
    real(kind=c_double), dimension(m), intent(in)  :: y3
    real(kind=c_double), dimension(m), intent(out) :: g
    integer(kind=c_int),               intent(out) :: ier
    call ezspline_interp3_cloud(spl,m,y1,y2,y3,g,ier)
    return
  end subroutine czspline_interp3_cloud

  ! array interpolation
  subroutine czspline_interp3_array(m1, m2, m3, y1, y2, y3, g, ier) bind(C)
    use ezspline
    implicit none
    integer(kind=c_int),                      intent(in)  :: m1
    integer(kind=c_int),                      intent(in)  :: m2
    integer(kind=c_int),                      intent(in)  :: m3
    real(kind=c_double), dimension(m1),       intent(in)  :: y1
    real(kind=c_double), dimension(m2),       intent(in)  :: y2
    real(kind=c_double), dimension(m3),       intent(in)  :: y3
    real(kind=c_double), dimension(m1,m2,m3), intent(out) :: g
    integer(kind=c_int),                      intent(out) :: ier
    call ezspline_interp3_array(spl,m1,m2,m3,y1,y2,y3,g,ier)
    return
  end subroutine czspline_interp3_array

  ! finalization
  subroutine czspline_free3(ier) bind(C)
    use ezspline
    implicit none
    integer(kind=c_int), intent(out) :: ier
    call ezspline_free3(spl,ier)
    return
  end subroutine czspline_free3

#ifdef _NETCDF
  ! save to file
  subroutine czspline_save3(cname, ier) bind(C)
    use ezspline
    implicit none
    character(kind=c_char,len=1), dimension(*), intent(in) :: cname
    integer(kind=c_int), intent(out) :: ier
    character(len=128) :: fname
    integer :: l
    l = 0
    fname = " "
    do
      if(cname(l+1) == C_NULL_CHAR) exit
      l = l + 1
      fname(l:l) = cname(l)
    end do
    call ezspline_save3(spl,trim(fname),ier)
    return
  end subroutine czspline_save3

  ! load from file
  subroutine czspline_load3(cname, ier) bind(C)
    use ezspline
    character(kind=c_char,len=1), dimension(*), intent(in) :: cname
    integer(kind=c_int), intent(out) :: ier
    character(len=128) :: fname
    integer :: l
    l = 0
    fname = " "
    do
      if(cname(l+1) == C_NULL_CHAR) exit
      l = l + 1
      fname(l:l) = cname(l)
    end do
    call ezspline_load3(spl,trim(fname),ier)
    return
  end subroutine czspline_load3
#endif

  ! point is in domain check
  subroutine czspline_isindomain3(y1, y2, y3, ier) bind(C)
    use ezspline
    implicit none
    real(kind=c_double), intent(in)  :: y1
    real(kind=c_double), intent(in)  :: y2
    real(kind=c_double), intent(in)  :: y3
    integer(kind=c_int), intent(out) :: ier
    call ezspline_isindomain3(spl,y1,y2,y3,ier)
    return
  end subroutine czspline_isindomain3

  ! cloud is in domain check
  subroutine czspline_isindomain3_cloud(m, y1, y2, y3, ier) bind(C)
    use ezspline
    implicit none
    integer(kind=c_int),               intent(in)  :: m
    real(kind=c_double), dimension(m), intent(in)  :: y1
    real(kind=c_double), dimension(m), intent(in)  :: y2
    real(kind=c_double), dimension(m), intent(in)  :: y3
    integer(kind=c_int),               intent(out) :: ier
    call ezspline_isindomain3_cloud(spl,m,y1,y2,y3,ier)
    return
  end subroutine czspline_isindomain3_cloud

  ! array is in domain check
  subroutine czspline_isindomain3_array(m1, m2, m3, y1, y2, y3, ier) bind(C)
    use ezspline
    implicit none
    integer(kind=c_int),                intent(in)  :: m1
    integer(kind=c_int),                intent(in)  :: m2
    integer(kind=c_int),                intent(in)  :: m3
    real(kind=c_double), dimension(m1), intent(in)  :: y1
    real(kind=c_double), dimension(m2), intent(in)  :: y2
    real(kind=c_double), dimension(m3), intent(in)  :: y3
    integer(kind=c_int),                intent(out) :: ier
    call ezspline_isindomain3_array(spl,m1,m2,m3,y1,y2,y3,ier)
    return
  end subroutine czspline_isindomain3_array

  ! gradient at a point
  subroutine czspline_gradient3(y1, y2, y3, df, ier) bind(C)
    use ezspline
    implicit none
    real(kind=c_double),               intent(in)  :: y1
    real(kind=c_double),               intent(in)  :: y2
    real(kind=c_double),               intent(in)  :: y3
    real(kind=c_double), dimension(3), intent(out) :: df
    integer(kind=c_int),               intent(out) :: ier
    call ezspline_gradient3(spl,y1,y2,y3,df,ier)
    return
  end subroutine czspline_gradient3

  ! cloud
  subroutine czspline_gradient3_cloud(m, y1, y2, y3, df, ier) bind(C)
    use ezspline
    implicit none
    integer(kind=c_int),                 intent(in)  :: m
    real(kind=c_double), dimension(m),   intent(in)  :: y1
    real(kind=c_double), dimension(m),   intent(in)  :: y2
    real(kind=c_double), dimension(m),   intent(in)  :: y3
    real(kind=c_double), dimension(m,3), intent(out) :: df
    integer(kind=c_int),                 intent(out) :: ier
    call ezspline_gradient3_cloud(spl,m,y1,y2,y3,df,ier)
    return
  end subroutine czspline_gradient3_cloud

  ! array
  subroutine czspline_gradient3_array(m1, m2, m3, y1, y2, y3, df, ier) bind(C)
    use ezspline
    implicit none
    integer(kind=c_int),                        intent(in)  :: m1
    integer(kind=c_int),                        intent(in)  :: m2
    integer(kind=c_int),                        intent(in)  :: m3
    real(kind=c_double), dimension(m1),         intent(in)  :: y1
    real(kind=c_double), dimension(m2),         intent(in)  :: y2
    real(kind=c_double), dimension(m3),         intent(in)  :: y3
    real(kind=c_double), dimension(m1,m2,m3,3), intent(out) :: df
    integer(kind=c_int),                        intent(out) :: ier
    call ezspline_gradient3_array(spl,m1,m2,m3,y1,y2,y3,df,ier)
    return
  end subroutine czspline_gradient3_array

  ! compute derivatives at a point
  subroutine czspline_derivative3(i1, i2, i3, y1, y2, y3, f, ier) bind(C)
    use ezspline
    implicit none
    integer(kind=c_int), intent(in)  :: i1
    integer(kind=c_int), intent(in)  :: i2
    integer(kind=c_int), intent(in)  :: i3
    real(kind=c_double), intent(in)  :: y1
    real(kind=c_double), intent(in)  :: y2
    real(kind=c_double), intent(in)  :: y3
    real(kind=c_double), intent(out) :: f
    integer(kind=c_int), intent(out) :: ier
    call ezspline_derivative3(spl,i1,i2,i3,y1,y2,y3,f,ier)
    return
  end subroutine czspline_derivative3

  ! cloud
  subroutine czspline_derivative3_cloud(i1, i2, i3, m, y1, y2, y3, f, ier) bind(C)
    use ezspline
    implicit none
    integer(kind=c_int),               intent(in)  :: i1
    integer(kind=c_int),               intent(in)  :: i2
    integer(kind=c_int),               intent(in)  :: i3
    integer(kind=c_int),               intent(in)  :: m
    real(kind=c_double), dimension(m), intent(in)  :: y1
    real(kind=c_double), dimension(m), intent(in)  :: y2
    real(kind=c_double), dimension(m), intent(in)  :: y3
    real(kind=c_double), dimension(m), intent(out) :: f
    integer(kind=c_int),               intent(out) :: ier
    call ezspline_derivative3_cloud(spl,i1,i2,i3,m,y1,y2,y3,f,ier)
    return
  end subroutine czspline_derivative3_cloud

  !array
  subroutine czspline_derivative3_array(i1, i2, i3, m1, m2, m3, y1, y2, y3, f, ier) bind(C)
    use ezspline
    implicit none
    integer(kind=c_int),                      intent(in) :: i1
    integer(kind=c_int),                      intent(in) :: i2
    integer(kind=c_int),                      intent(in) :: i3
    integer(kind=c_int),                      intent(in) :: m1
    integer(kind=c_int),                      intent(in) :: m2
    integer(kind=c_int),                      intent(in) :: m3
    real(kind=c_double), dimension(m1),       intent(in)  :: y1
    real(kind=c_double), dimension(m2),       intent(in)  :: y2
    real(kind=c_double), dimension(m3),       intent(in)  :: y3
    real(kind=c_double), dimension(m1,m2,m3), intent(out) :: f
    integer(kind=c_int),                      intent(out) :: ier
    call ezspline_derivative3_array(spl,i1,i2,i3,m1,m2,m3,y1,y2,y3,f,ier)
    return
  end subroutine czspline_derivative3_array

  ! check to see if grid is strictly increasing
  subroutine czspline_isgridregular3(ier) bind(C)
    use ezspline
    implicit none
    integer(kind=c_int), intent(out) :: ier
    call ezspline_isgridregular3(spl, ier)
    return
  end subroutine czspline_isgridregular3

end module czspline3d
