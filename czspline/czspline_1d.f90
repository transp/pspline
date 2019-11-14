module czspline1d
  use iso_c_binding
  use EZspline_type, only: EZspline1

  type(EZspline1), private :: spl

  public :: czspline_init1
  public :: czspline_set_axes1
  public :: czspline_set_ishermite1
  public :: czspline_set_bctypes1
  public :: czspline_set_bcvals1
  public :: czspline_setup1
  public :: czspline_interp1
  public :: czspline_interp1_cloud
  public :: czspline_interp1_array
  public :: czspline_free1
  public :: czspline_save1
  public :: czspline_load1
  public :: czspline_isindomain1
  public :: czspline_isindomain1_cloud
  public :: czspline_isindomain1_array
  public :: czspline_gradient1
  public :: czspline_gradient1_cloud
  public :: czspline_gradient1_array
  public :: czspline_derivative1
  public :: czspline_derivative1_cloud
  public :: czspline_derivative1_array
  public :: czspline_isgridregular1

contains

  subroutine czspline_init1(n, bcs, ier) bind(C)
    use ezspline
    implicit none
    integer(kind=c_int),               intent(in)  :: n
    integer(kind=c_int), dimension(*), intent(in)  :: bcs
    integer(kind=c_int),               intent(out) :: ier
    call ezspline_init1(spl,n,bcs,ier)
    return
  end subroutine czspline_init1

  subroutine czspline_set_axes1(n, x, ier) bind(C)
    use ezspline
    implicit none
    integer(kind=c_int),               intent(in)  :: n
    real(kind=c_double), dimension(n), intent(in)  :: x
    integer(kind=c_int),               intent(out) :: ier
    ier = 0
    spl%x1 = x
    return
  end subroutine czspline_set_axes1

  ! isHermite
  subroutine czspline_set_ishermite1(flag, ier) bind(C)
    use ezspline
    implicit none
    integer(kind=c_int), intent(in)  :: flag
    integer(kind=c_int), intent(out) :: ier
    ier = 0
    spl%isHermite = flag
    return
  end subroutine czspline_set_ishermite1

  ! Boundary types
  subroutine czspline_set_bctypes1(bctype1, ier) bind(C)
    use ezspline
    implicit none
    integer(kind=c_int), dimension(2), intent(in)  :: bctype1
    integer(kind=c_int),               intent(out) :: ier
    ier = 0
    spl%ibctype1 = bctype1
    return
  end subroutine czspline_set_bctypes1

  ! Boundary conditions bcval1min, bcval1max
  subroutine czspline_set_bcvals1(bcval1, ier) bind(C)
    use ezspline
    implicit none
    real(kind=c_double), dimension(2), intent(in)  :: bcval1
    integer(kind=c_int),               intent(out) :: ier
    ier = 0
    spl%bcval1min = bcval1(1)
    spl%bcval1max = bcval1(2)
    return
  end subroutine czspline_set_bcvals1

  subroutine czspline_setup1(n, f, ier) bind(C)
    use ezspline
    implicit none
    integer(kind=c_int),               intent(in)  :: n
    real(kind=c_double), dimension(n), intent(in)  :: f
    integer(kind=c_int),               intent(out) :: ier
    call ezspline_setup1(spl,f,ier)
    return
  end subroutine czspline_setup1

  subroutine czspline_interp1(y, g, ier) bind(C)
    use ezspline
    implicit none
    real(kind=c_double), intent(in)  :: y
    real(kind=c_double), intent(out) :: g
    integer(kind=c_int), intent(out) :: ier
    call ezspline_interp1(spl,y,g,ier)
    return
  end subroutine czspline_interp1

  subroutine czspline_interp1_cloud(n, y, g, ier) bind(C)
    use ezspline
    implicit none
    integer(kind=c_int),               intent(in)  :: n
    real(kind=c_double), dimension(n), intent(in)  :: y
    real(kind=c_double), dimension(n), intent(out) :: g
    integer(kind=c_int), intent(out) :: ier
    call ezspline_interp1_array(spl,n,y,g,ier)
    return
  end subroutine czspline_interp1_cloud

  subroutine czspline_interp1_array(n, y, g, ier) bind(C)
    use ezspline
    implicit none
    integer(kind=c_int),               intent(in)  :: n
    real(kind=c_double), dimension(n), intent(in)  :: y
    real(kind=c_double), dimension(n), intent(out) :: g
    integer(kind=c_int), intent(out) :: ier
    call ezspline_interp1_array(spl,n,y,g,ier)
    return
  end subroutine czspline_interp1_array

  subroutine czspline_free1(ier) bind(C)
    use ezspline
    integer(kind=c_int), intent(out) :: ier
    call ezspline_free1(spl,ier)
    return
  end subroutine czspline_free1

#ifdef _NETCDF
  subroutine czspline_save1(cname, ier) bind(C)
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
    call ezspline_save1(spl,trim(fname),ier)
    return
  end subroutine czspline_save1

  subroutine czspline_load1(cname, ier) bind(C)
    use ezspline
    implicit none
    character(kind=c_char,len=1), dimension(*), intent(in)  :: cname
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
    call ezspline_load1(spl,trim(fname),ier)
    return
  end subroutine czspline_load1
#endif

  !point
  subroutine czspline_isindomain1(y, ier) bind(C)
    use ezspline
    implicit none
    real(kind=c_double), intent(in)  :: y
    integer(kind=c_int), intent(out) :: ier
    call ezspline_isindomain(spl, y, ier)
    return
  end subroutine czspline_isindomain1

  !cloud
  subroutine czspline_isindomain1_cloud(m, y, ier) bind(C)
    use ezspline
    implicit none
    integer(kind=c_int),               intent(in)  :: m
    real(kind=c_double), dimension(m), intent(in)  :: y
    integer(kind=c_int),               intent(out) :: ier
    call ezspline_isindomain(spl,m,y,ier)
    return
  end subroutine czspline_isindomain1_cloud

  !array
  subroutine czspline_isindomain1_array(m, y, ier) bind(C)
    use ezspline
    implicit none
    integer(kind=c_int),               intent(in)  :: m
    real(kind=c_double), dimension(m), intent(in)  :: y
    integer(kind=c_int),               intent(out) :: ier
    call ezspline_isindomain(spl,m,y,ier)
    return
  end subroutine czspline_isindomain1_array

  ! point
  subroutine czspline_gradient1(y, df, ier) bind(C)
    use ezspline
    implicit none
    real(kind=c_double), intent(in)  :: y
    real(kind=c_double), intent(out) :: df
    integer(kind=c_int), intent(out) :: ier
    call ezspline_gradient1(spl,y,df,ier)
    return
  end subroutine czspline_gradient1

  ! cloud
  subroutine czspline_gradient1_cloud(m, y, df, ier) bind(C)
    use ezspline
    implicit none
    integer(kind=c_int),               intent(in)  :: m
    real(kind=c_double), dimension(m), intent(in)  :: y
    real(kind=c_double), dimension(m), intent(out) :: df
    integer(kind=c_int),               intent(out) :: ier
    call ezspline_gradient1_array(spl,m,y,df,ier)
    return
  end subroutine czspline_gradient1_cloud

  ! array
  subroutine czspline_gradient1_array(m, y, df, ier) bind(C)
    use ezspline
    implicit none
    integer(kind=c_int),               intent(in)  :: m
    real(kind=c_double), dimension(m), intent(in)  :: y
    real(kind=c_double), dimension(m), intent(out) :: df
    integer(kind=c_int),               intent(out) :: ier
    call ezspline_gradient1_array(spl,m,y,df,ier)
    return
  end subroutine czspline_gradient1_array

  ! point
  subroutine czspline_derivative1(i, y, f, ier) bind(C)
    use ezspline
    implicit none
    integer(kind=c_int), intent(in)  :: i
    real(kind=c_double), intent(in)  :: y
    real(kind=c_double), intent(out) :: f
    integer(kind=c_int), intent(out) :: ier
    call ezspline_derivative1(spl,i,y,f,ier)
    return
  end subroutine czspline_derivative1

  ! cloud
  subroutine czspline_derivative1_cloud(i, m, y, f, ier) bind(C)
    use ezspline
    implicit none
    integer(kind=c_int),               intent(in)  :: i
    integer(kind=c_int),               intent(in)  :: m
    real(kind=c_double), dimension(m), intent(in)  :: y
    real(kind=c_double), dimension(m), intent(out) :: f
    integer(kind=c_int),               intent(out) :: ier
    call ezspline_derivative1_array(spl,i,m,y,f,ier)
    return
  end subroutine czspline_derivative1_cloud

  !array
  subroutine czspline_derivative1_array(i, m, y, f, ier) bind(C)
    use ezspline
    implicit none
    integer(kind=c_int),               intent(in)  :: i
    integer(kind=c_int),               intent(in)  :: m
    real(kind=c_double), dimension(m), intent(in)  :: y
    real(kind=c_double), dimension(m), intent(out) :: f
    integer(kind=c_int),               intent(out) :: ier
    call ezspline_derivative1_array(spl,i,m,y,f,ier)
    return
  end subroutine czspline_derivative1_array

  subroutine czspline_isgridregular1(ier) bind(C)
    use ezspline
    implicit none
    integer(kind=c_int), intent(out) :: ier
    call ezspline_isgridregular(spl,ier)
    return
  end subroutine czspline_isgridregular1

end module czspline1d
