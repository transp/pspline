!---------------------------------------------------------------------------
! This code was developed at Tech-X (www.txcorp.com). It is free for any one
! to use but comes with no warranty whatsoever. Use at your own risk. 
! Thanks for reporting bugs to pletzer@txcorp.com. 
!---------------------------------------------------------------------------
! Boiler plate code to generate pointer to types for 1, 2, and 3 dimensions

module czspline_pointer_types

  use ezspline_obj
  use ezspline

  type czspline1
    type(ezspline1), pointer :: ptr => NULL()
  end type czspline1

  type czspline2
    type(ezspline2), pointer :: ptr => NULL()
  end type czspline2

  type czspline3
    type(ezspline3), pointer :: ptr => NULL()
  end type czspline3

end module czspline_pointer_types
