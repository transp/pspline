!******************** START FILE SPLINE.F90 ; GROUP TRKRLIB ******************
!C
!CCCCCCCCCCCCCCC
!	THE CODES (SPLINE & SEVAL) ARE TAKEN FROM:
!	FORSYTHE,MALCOLM AND MOLER, "COMPUTER METHODS FOR
!	MATHEMATICAL COMPUTATIONS",PRENTICE-HALL, 1977.
!
!	THE CODES (SPLEEN,SPLAAN & SPEVAL) ARE ADAPTATIONS
!	BY R.M. WIELAND FOR SPECIAL CASES ... SEE COMMENTS
!
!  this is nspline -- natural cubic spline
!    second derivative of interpolant set to zero at endpts
!

subroutine nspline(N,X,Y,B,C,D)
  use precision_mod, only: fp
  !============
  implicit none
  integer ip1,ii
  !============
  real(fp) :: zbomb
  !============
  integer N
  real(fp) :: X(N), Y(N), B(N), C(N), D(N)
  !
  !  THE COEFFICIENTS B(I), C(I), AND D(I), I=1,2,...,N ARE COMPUTED
  !  FOR A CUBIC INTERPOLATING SPLINE
  !
  !    S(X) = Y(I) + B(I)*(X-X(I)) + C(I)*(X-X(I))**2 + D(I)*(X-X(I))**3
  !
  !    FOR  X(I) .LE. X .LE. X(I+1)
  !
  !  INPUT..
  !
  !    N = THE NUMBER OF DATA POINTS OR KNOTS (N.GE.2)
  !    X = THE ABSCISSAS OF THE KNOTS IN STRICTLY INCREASING ORDER
  !    Y = THE ORDINATES OF THE KNOTS
  !
  !  OUTPUT..
  !
  !    B, C, D  = ARRAYS OF SPLINE COEFFICIENTS AS DEFINED ABOVE.
  !
  !  USING  P  TO DENOTE DIFFERENTIATION,
  !
  !    Y(I) = S(X(I))
  !    B(I) = SP(X(I))
  !    C(I) = SPP(X(I))/2
  !    D(I) = SPPP(X(I))/6  (DERIVATIVE FROM THE RIGHT)
  !
  !CCCCCCCCCCCCCC
  !  THE ACCOMPANYING FUNCTION SUBPROGRAM  SEVAL  CAN BE USED
  !  TO EVALUATE THE SPLINE.
  !
  !
  integer NM1, NM2, IB, I
  real(fp) :: T
  !
  NM1 = N-1
  IF ( N .LT. 2 ) RETURN
  !
  !
  ! DMC - CHECK ORDINATES
  DO I=1,NM1
    IP1=I+1
    IF(X(IP1).LE.X(I)) then
      write(6,9000)
      write(6,9001) (x(ii),ii=1,N)
      zbomb=sqrt(-1.0_fp/(x(1)*x(1)))
      write(6,*) zbomb
      write(6,*) ' nspline -- code stop.'
      stop 1
9000  FORMAT(/ &
           ' ? UNORDERED ORDINATE ARRAY PASSED TO SPLINE subroutine:')
9001  FORMAT(1X,5(1X,1PE12.5))
    end if
  end do
  !
  IF ( N .LT. 3 ) GO TO 50
  !
  !  SET UP TRIDIAGONAL SYSTEM
  !
  !  B = DIAGONAL, D = OFFDIAGONAL, C = RIGHT HAND SIDE.
  !
  D(1) = X(2) - X(1)
  C(2) = (Y(2) - Y(1))/D(1)
  DO I = 2, NM1
    D(I) = X(I+1) - X(I)
    B(I) = 2._fp*(D(I-1) + D(I))
    C(I+1) = (Y(I+1) - Y(I))/D(I)
    C(I) = C(I+1) - C(I)
  end do
  !
  !  END CONDITIONS:  2nd deriv = 0 at endpts; solve reduced n-2 x n-2
  !  system.
  C(1) = 0._fp
  C(N) = 0._fp
  !
  !  FORWARD ELIMINATION
  !
15 DO I = 3, NM1
    T = D(I-1)/B(I-1)
    B(I) = B(I) - T*D(I-1)
    C(I) = C(I) - T*C(I-1)
  end do
  !
  !  BACK SUBSTITUTION
  !
  NM2=NM1-1
  C(NM1) = C(NM1)/B(NM1)
  DO IB = 2, NM2
    I = N-IB
    C(I) = (C(I) - D(I)*C(I+1))/B(I)
  end do
  !
  !  C(I) IS NOW THE SIGMA(I) OF THE TEXT
  !
  !  COMPUTE POLYNOMIAL COEFFICIENTS
  !
  B(N) = (Y(N) - Y(NM1))/D(NM1) + D(NM1)*(C(NM1) + 2._fp*C(N))
  DO I = 1, NM1
    B(I) = (Y(I+1) - Y(I))/D(I) - D(I)*(C(I+1) + 2._fp*C(I))
    D(I) = (C(I+1) - C(I))/D(I)
    C(I) = 3._fp*C(I)
  end do
  C(N) = 3._fp*C(N)
  D(N) = D(N-1)
  RETURN
  !
50 B(1) = (Y(2)-Y(1))/(X(2)-X(1))
  C(1) = 0._fp
  D(1) = 0._fp
  B(2) = B(1)
  C(2) = 0._fp
  D(2) = 0._fp
  RETURN
end subroutine nspline
