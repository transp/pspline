module czspline2d
  use iso_c_binding
  use EZspline_type, only: EZspline2

  type(EZspline2), private :: spl

  public :: czspline_init2
  public :: czspline_set_axes2
  public :: czspline_set_ishermite2
  public :: czspline_set_bctypes2
  public :: czspline_set_bcvals2
  public :: czspline_setup2
  public :: czspline_interp2
  public :: czspline_interp2_cloud
  public :: czspline_interp2_array
  public :: czspline_free2
  public :: czspline_save2
  public :: czspline_load2
  public :: czspline_isindomain2
  public :: czspline_isindomain2_cloud
  public :: czspline_isindomain2_array
  public :: czspline_gradient2
  public :: czspline_gradient2_cloud
  public :: czspline_gradient2_array
  public :: czspline_derivative2
  public :: czspline_derivative2_cloud
  public :: czspline_derivative2_array
  public :: czspline_isgridregular2

contains

  ! initialization
  subroutine czspline_init2(n1, n2, bcs1, bcs2, ier) bind(C)
    use ezspline
    implicit none
    integer(kind=c_int),               intent(in)  :: n1
    integer(kind=c_int),               intent(in)  :: n2
    integer(kind=c_int), dimension(*), intent(in)  :: bcs1
    integer(kind=c_int), dimension(*), intent(in)  :: bcs2
    integer(kind=c_int),               intent(out) :: ier
    call ezspline_init2(spl,n1,n2,bcs1,bcs2,ier)
    return
  end subroutine czspline_init2

  ! member setters
  subroutine czspline_set_axes2(n1, n2, x1, x2, ier) bind(C)
    use ezspline
    implicit none
    integer(kind=c_int),                intent(in)  :: n1
    integer(kind=c_int),                intent(in)  :: n2
    real(kind=c_double), dimension(n1), intent(in)  :: x1
    real(kind=c_double), dimension(n2), intent(in)  :: x2
    integer(kind=c_int),                intent(out) :: ier
    ier = 0
    spl%x1 = x1
    spl%x2 = x2
    return
  end subroutine czspline_set_axes2

  subroutine czspline_set_ishermite2(flag, ier) bind(C)
    use ezspline
    implicit none
    integer(kind=c_int), intent(in)  :: flag
    integer(kind=c_int), intent(out) :: ier
    ier = 0
    spl%isHermite = flag
    return
  end subroutine czspline_set_ishermite2

  ! boundary types
  subroutine czspline_set_bctypes2(bctype1, bctype2, ier) bind(C)
    use ezspline
    implicit none
    integer(kind=c_int), dimension(2), intent(in)  :: bctype1
    integer(kind=c_int), dimension(2), intent(in)  :: bctype2
    integer(kind=c_int),               intent(out) :: ier
    ier = 0
    spl%ibctype1 = bctype1
    spl%ibctype2 = bctype2
    return
  end subroutine czspline_set_bctypes2

  ! boundary conditions bcval1min, bcval1max
  subroutine czspline_set_bcvals2(bcval1, bcval2, ier)
    use ezspline
    implicit none
    real(kind=c_double), dimension(2), intent(in)  :: bcval1
    real(kind=c_double), dimension(2), intent(in)  :: bcval2
    integer(kind=c_int),               intent(out) :: ier
    ier = 0
    spl%bcval1min = bcval1(1)
    spl%bcval1max = bcval1(2)
    spl%bcval2min = bcval2(1)
    spl%bcval2max = bcval2(2)
    return
  end subroutine czspline_set_bcvals2

  ! compute spline coefficients
  subroutine czspline_setup2(n1, n2, f, ier) bind(C)
    use ezspline
    implicit none
    integer(kind=c_int),                   intent(in)  :: n1
    integer(kind=c_int),                   intent(in)  :: n2
    real(kind=c_double), dimension(n1,n2), intent(in)  :: f
    integer(kind=c_int),                   intent(out) :: ier
    call ezspline_setup2(spl,f,ier)
    return
  end subroutine czspline_setup2

  ! point interpolation
  subroutine czspline_interp2(y1, y2, g, ier) bind(C)
    use ezspline
    implicit none
    real(kind=c_double), intent(in)  :: y1
    real(kind=c_double), intent(in)  :: y2
    real(kind=c_double), intent(out) :: g
    integer(kind=c_int), intent(out) :: ier
    call ezspline_interp2(spl,y1,y2,g,ier)
    return
  end subroutine czspline_interp2

  ! cloud interpolation
  subroutine czspline_interp2_cloud(m, y1, y2, g, ier) bind(C)
    use ezspline
    implicit none
    integer(kind=c_int),               intent(in)  :: m
    real(kind=c_double), dimension(m), intent(in)  :: y1
    real(kind=c_double), dimension(m), intent(in)  :: y2
    real(kind=c_double), dimension(m), intent(out) :: g
    integer(kind=c_int),               intent(out) :: ier
    call ezspline_interp2_cloud(spl,m,y1,y2,g,ier)
    return
  end subroutine czspline_interp2_cloud

  ! array interpolation
  subroutine czspline_interp2_array(m1, m2, y1, y2, g, ier) bind(C)
    use ezspline
    implicit none
    integer(kind=c_int),                   intent(in)  :: m1
    integer(kind=c_int),                   intent(in)  :: m2
    real(kind=c_double), dimension(m1),    intent(in)  :: y1
    real(kind=c_double), dimension(m2),    intent(in)  :: y2
    real(kind=c_double), dimension(m1,m2), intent(out) :: g
    integer(kind=c_int),                   intent(out) :: ier
    call ezspline_interp2_array(spl,m1,m2,y1,y2,g,ier)
    return
  end subroutine czspline_interp2_array

  ! finalization
  subroutine czspline_free2(ier) bind(C)
    use ezspline
    implicit none
    integer(kind=c_int), intent(out) :: ier
    call ezspline_free2(spl,ier)
    return
  end subroutine czspline_free2

#ifdef _EZCDF
  ! save to file
  subroutine czspline_save2(cname, ier) bind(C)
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
    call ezspline_save2(spl,trim(fname),ier)
    return
  end subroutine czspline_save2

  ! load from file
  subroutine czspline_load2(cname, ier) bind(C)
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
    call ezspline_load2(spl,trim(fname),ier)
    return
  end subroutine czspline_load2
#endif

  ! point is in domain check
  subroutine czspline_isindomain2(y1, y2, ier) bind(C)
    use ezspline
    implicit none
    real(kind=c_double), intent(in)  :: y1
    real(kind=c_double), intent(in)  :: y2
    integer(kind=c_int), intent(out) :: ier
    call ezspline_isindomain2(spl,y1,y2,ier)
    return
  end subroutine czspline_isindomain2

  ! cloud is in domain check
  subroutine czspline_isindomain2_cloud(m, y1, y2, ier) bind(C)
    use ezspline
    implicit none
    integer(kind=c_int),               intent(in)  :: m
    real(kind=c_double), dimension(m), intent(in)  :: y1
    real(kind=c_double), dimension(m), intent(in)  :: y2
    integer(kind=c_int),               intent(out) :: ier
    call ezspline_isindomain2_cloud(spl,m,y1,y2,ier)
    return
  end subroutine czspline_isindomain2_cloud

  ! array is in domain check
  subroutine czspline_isindomain2_array(m1, m2, y1, y2, ier) bind(C)
    use ezspline
    implicit none
    integer(kind=c_int),                intent(in)  :: m1
    integer(kind=c_int),                intent(in)  :: m2
    real(kind=c_double), dimension(m1), intent(in)  :: y1
    real(kind=c_double), dimension(m2), intent(in)  :: y2
    integer(kind=c_int),                intent(out) :: ier
    call ezspline_isindomain2_array(spl,m1,m2,y1,y2,ier)
    return
  end subroutine czspline_isindomain2_array

  ! point
  subroutine czspline_gradient2(y1, y2, df, ier) bind(C)
    use ezspline
    implicit none
    real(kind=c_double),               intent(in)  :: y1
    real(kind=c_double),               intent(in)  :: y2
    real(kind=c_double), dimension(2), intent(out) :: df
    integer(kind=c_int),               intent(out) :: ier
    call ezspline_gradient2(spl,y1,y2,df,ier)
    return
  end subroutine czspline_gradient2

  ! cloud
  subroutine czspline_gradient2_cloud(m, y1, y2, df, ier) bind(C)
    use ezspline
    implicit none
    integer(kind=c_int),                 intent(in)  :: m
    real(kind=c_double), dimension(m),   intent(in)  :: y1
    real(kind=c_double), dimension(m),   intent(in)  :: y2
    real(kind=c_double), dimension(m,2), intent(out) :: df
    integer(kind=c_int),                 intent(out) :: ier
    call ezspline_gradient2_cloud(spl,m,y1,y2,df,ier)
    return
  end subroutine czspline_gradient2_cloud

  ! array
  subroutine czspline_gradient2_array(m1, m2, y1, y2, df, ier) bind(C)
    use ezspline
    implicit none
    integer(kind=c_int),                     intent(in)  :: m1
    integer(kind=c_int),                     intent(in)  :: m2
    real(kind=c_double), dimension(m1),      intent(in)  :: y1
    real(kind=c_double), dimension(m2),      intent(in)  :: y2
    real(kind=c_double), dimension(m1,m2,2), intent(out) :: df
    integer(kind=c_int),                     intent(out) :: ier
    call ezspline_gradient2_array(spl,m1,m2,y1,y2,df,ier)
    return
  end subroutine czspline_gradient2_array

  ! point
  subroutine czspline_derivative2(i1, i2, y1, y2, f, ier) bind(C)
    use ezspline
    implicit none
    integer(kind=c_int), intent(in)  :: i1
    integer(kind=c_int), intent(in)  :: i2
    real(kind=c_double), intent(in)  :: y1
    real(kind=c_double), intent(in)  :: y2
    real(kind=c_double), intent(out) :: f
    integer(kind=c_int), intent(out) :: ier
    call ezspline_derivative2(spl, i1, i2, y1, y2, f, ier)
    return
  end subroutine czspline_derivative2

  ! cloud
  subroutine czspline_derivative2_cloud(i1, i2, m, y1, y2, f, ier) bind(C)
    use ezspline
    implicit none
    integer(kind=c_int),               intent(in)  :: i1
    integer(kind=c_int),               intent(in)  :: i2
    integer(kind=c_int),               intent(in)  :: m
    real(kind=c_double), dimension(m), intent(in)  :: y1
    real(kind=c_double), dimension(m), intent(in)  :: y2
    real(kind=c_double), dimension(m), intent(out) :: f
    integer(kind=c_int),               intent(out) :: ier
    call ezspline_derivative2_cloud(spl,i1,i2,m,y1,y2,f,ier)
    return
  end subroutine czspline_derivative2_cloud

  !array
  subroutine czspline_derivative2_array(i1, i2, m1, m2, y1, y2, f, ier) bind(C)
    use ezspline
    implicit none
    integer(kind=c_int),                   intent(in)  :: i1
    integer(kind=c_int),                   intent(in)  :: i2
    integer(kind=c_int),                   intent(in)  :: m1
    integer(kind=c_int),                   intent(in)  :: m2
    real(kind=c_double), dimension(m1),    intent(in)  :: y1
    real(kind=c_double), dimension(m2),    intent(in)  :: y2
    real(kind=c_double), dimension(m1,m2), intent(out) :: f
    integer(kind=c_int),                   intent(out) :: ier
    call ezspline_derivative2_array(spl,i1,i2,m1,m2,y1,y2,f,ier)
    return
  end subroutine czspline_derivative2_array

  ! check to see if grid is strictly increasing
  subroutine czspline_isgridregular2(ier) bind(C)
    use ezspline
    implicit none
    integer(kind=c_int), intent(out) :: ier
    call ezspline_isgridregular2(spl,ier)
    return
  end subroutine czspline_isgridregular2

end module czspline2d
