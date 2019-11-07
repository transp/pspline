subroutine util_bcherm1(fherm,idimx1, &
     jbcxa,jbcxb, &
     zbcxa,zbcxb, &
     x1)
  use psp_precision_mod, only: fp

  !...  insert BCs as needed for Hermite interpolation setup


  implicit none
  integer :: idimx1         ! array dimensions
  real(fp) :: fherm(0:1,idimx1)

  integer :: jbcxa,jbcxb    ! x BC controls
  real(fp) :: zbcxa            ! x(1) BC data
  real(fp) :: zbcxb            ! x(idimx1) BC data

  real(fp) :: x1(idimx1)

  !-------------------------------------------------------------
  real(fp) :: zdx
  !-------------------------------------------------------------

  if((jbcxa.eq.1).or.(jbcxb.eq.1)) then

     if(jbcxa.eq.1) then
        fherm(1,1)=zbcxa
     else
        zdx = x1(2)-x1(1)
        fherm(1,1)=(fherm(0,2)-fherm(0,1))/zdx
     end if

     if(jbcxb.eq.1) then
        fherm(1,idimx1)=zbcxb
     else
        zdx = x1(idimx1)-x1(idimx1-1)
        fherm(1,idimx1)=(fherm(0,idimx1)-fherm(0,idimx1-1))/zdx
     end if

  end if

end subroutine util_bcherm1
