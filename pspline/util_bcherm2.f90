subroutine util_bcherm2(fherm,idimx1,idimx2, &
     jbcx1a,jbcx1b,   jbcx2a,jbcx2b, &
     zbcx1a,zbcx1b,   zbcx2a,zbcx2b, &
     x1, x2)
  use precision_mod, only: fp

  !...  insert BCs as needed for Hermite interpolation setup


  implicit none
  integer :: idimx1,idimx2  ! array dimensions
  real(fp) :: fherm(0:3,idimx1,idimx2)

  integer :: jbcx1a,jbcx1b  ! x1 BC controls
  real(fp) :: zbcx1a(idimx2)    ! x1(1) BC data
  real(fp) :: zbcx1b(idimx2)    ! x1(idimx1) BC data

  integer :: jbcx2a,jbcx2b  ! x2 BC controls
  real(fp) :: zbcx2a(idimx1)    ! x2(1) BC data
  real(fp) :: zbcx2b(idimx1)    ! x2(idimx2) BC data

  real(fp) :: x1(idimx1),x2(idimx2)

  !-------------------------------------------------------------
  integer :: ix
  real(fp) :: zdx
  !-------------------------------------------------------------

  if((jbcx1a.eq.1).or.(jbcx1b.eq.1)) then

     if(jbcx1a.eq.1) then
        fherm(1,1,1:idimx2)=zbcx1a(1:idimx2)
     else
        zdx = x1(2)-x1(1)
        do ix=1,idimx2
           fherm(1,1,ix)=(fherm(0,2,ix)-fherm(0,1,ix))/zdx
        end do
     end if

     if(jbcx1b.eq.1) then
        fherm(1,idimx1,1:idimx2)=zbcx1b(1:idimx2)
     else
        zdx = x1(idimx1)-x1(idimx1-1)
        do ix=1,idimx2
           fherm(1,idimx1,ix)= &
                (fherm(0,idimx1,ix)-fherm(0,idimx1-1,ix))/zdx
        end do
     end if

  end if

  if((jbcx2a.eq.1).or.(jbcx2b.eq.1)) then

     if(jbcx2a.eq.1) then
        fherm(2,1:idimx1,1)=zbcx2a(1:idimx1)
     else
        zdx=x2(2)-x2(1)
        do ix=1,idimx1
           fherm(2,ix,1)=(fherm(0,ix,2)-fherm(0,ix,1))/zdx
        end do
     end if

     if(jbcx2b.eq.1) then
        fherm(2,1:idimx1,idimx2)=zbcx2b(1:idimx1)
     else
        zdx=x2(idimx2)-x2(idimx2-1)
        do ix=1,idimx1
           fherm(2,ix,idimx2)= &
                (fherm(0,ix,idimx2)-fherm(0,ix,idimx2-1))/zdx
        end do
     end if
  end if

end subroutine util_bcherm2
