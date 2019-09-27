!  PSPLINE -- modified from SPLINE.F90 -- dmc 3 Oct 1995
!
!	THE CODES (SPLINE & SEVAL) ARE TAKEN FROM:
!	FORSYTHE,MALCOLM AND MOLER, "COMPUTER METHODS FOR
!	MATHEMATICAL COMPUTATIONS",PRENTICE-HALL, 1977.
!
!  PSPLINE -- periodic spline -- adaptation of SPLINE.F90 by D. McCune
!  October 1995:  The function to be splined is taken as having periodic
!  derivatives (it need not be periodic itself), so that
!  the interpolating function satisfies
!            y'(x(1))=y'(x(N))
!           y''(x(1))=y''(x(N))
!  this results in a fully specified system (no addl BC assumptions
!  required); it is not a tridiagonal system but a modified NM1xNM1
!  tridiagonal system with non-zero offdiagonal corner (1,NM1), (NM1,1)
!  elements.  It remains symmetric & diagonally dominant, and is
!  solved by Gaussian elimination w/o pivoting.  NM1=N-1.
!
subroutine pspline (N, X, Y, B, C, D, WK)
  use precision_mod, only: fp
  !============
  implicit none
  real(fp) :: q
  !============
  integer N
  real(fp) :: X(N), Y(N), B(N), C(N), D(N), WK(N)
  !
  !  THE COEFFICIENTS B(I), C(I), AND D(I), I=1,2,...,N ARE COMPUTED
  !  FOR A periodic CUBIC INTERPOLATING SPLINE
  !
  !    S(X) = Y(I) + B(I)*(X-X(I)) + C(I)*(X-X(I))**2 + D(I)*(X-X(I))**3
  !
  !    FOR  X(I) .LE. X .LE. X(I+1)
  !
  !    WK(...) is used as a workspace during the Guassian elimination.
  !    it represents the rightmost column of the array
  !
  !    note B(N),C(N),D(N) give coeffs for a periodic continuation into
  !    the next period, up to X(N)+(X(2)-X(1)) only.
  !
  !    SEVAL can be used to evaluate the spline, but, it is up to the
  !    SEVAL caller to bring the argument into the normal range X(1) to X(N).
  !
  !  INPUT..
  !
  !    N = THE NUMBER OF DATA POINTS OR KNOTS (N.GE.2)
  !        includes a complete period of the data, the last point being
  !        a repeat of the first point.
  !    X = THE ABSCISSAS OF THE KNOTS IN STRICTLY INCREASING ORDER
  !        X(N)-X(1) is the period of the periodic function represented.
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
  !-------------------------------------
  !
  NM1 = N-1
  NM2 = NM1-1
  !
  IF ( N .LT. 4 ) THEN
    write(6,9901)
9901 format(/ &
         ' ?PSPLINE -- at least 4 pts required for periodic spline.')
    return
  end if
  !
  !  SET UP MODIFIED NM1 x NM1 TRIDIAGONAL SYSTEM:
  !  B = DIAGONAL, D = OFFDIAGONAL, C = RIGHT HAND SIDE.
  !  WK(1:NM2) = rightmost column above diagonal
  !  WK(NM1)   = lower left corner element
  !
  D(1) = X(2) - X(1)
  D(NM1) = X(N) - X(NM1)
  B(1) = 2.0_fp*(D(1)+D(NM1))
  C(1) = (Y(N) - Y(NM1))/D(NM1)
  C(2) = (Y(2) - Y(1))/D(1)
  C(1) = C(2) - C(1)
  WK(1) = D(NM1)
  DO I = 2, NM1
    D(I) = X(I+1) - X(I)
    B(I) = 2._fp*(D(I-1) + D(I))
    C(I+1) = (Y(I+1) - Y(I))/D(I)
    C(I) = C(I+1) - C(I)
    WK(I) = 0.0_fp
  end do
  WK(NM2) = D(NM2)
  WK(NM1) = D(NM1)
  !
  !  END CONDITIONS -- implied by periodicity
  !     C(1)=C(N)  B(1)=B(N)  D(1)=D(N) -- no additional assumption needed.
  !
  !  FORWARD ELIMINATION
  !   WK(1)--WK(NM2) represent the rightmost column above the
  !   diagonal; WK(NM1) represents the non-zero lower left corner element
  !   which in each step is moved one column to the right.
  !
15 DO I = 2, NM2
    T = D(I-1)/B(I-1)
    B(I) = B(I) - T*D(I-1)
    C(I) = C(I) - T*C(I-1)
    WK(I) = WK(I) - T*WK(I-1)
    Q = WK(NM1)/B(I-1)
    WK(NM1) = -Q*D(I-1)
    B(NM1) = B(NM1) - Q*WK(I-1)
    C(NM1) = C(NM1) - Q*C(I-1)
  end do
  !
  !  correct the (NM1,NM2) element
  !
  WK(NM1) = WK(NM1) + D(NM2)
  !
  !  complete the forward elimination:  now WK(NM1) and WK(NM2) are
  !  the off diagonal elements of the 2x2 at the lower right hand corner
  !
  T = WK(NM1)/B(NM2)
  B(NM1) = B(NM1) - T*WK(NM2)
  C(NM1) = C(NM1) - T*C(NM2)
  !
  !  BACK SUBSTITUTION
  !
  C(NM1) = C(NM1)/B(NM1)
  C(NM2) = (C(NM2) - WK(NM2)*C(NM1))/B(NM2)
  !
  DO IB = 3, NM1
    I = N-IB
    C(I) = (C(I) - D(I)*C(I+1) - WK(I)*C(NM1))/B(I)
  end do
  !
  !  PERIODICITY:
  !
  C(N)=C(1)
  !
  !  C(I) IS NOW THE SIGMA(I) OF THE TEXT
  !
  !  COMPUTE POLYNOMIAL COEFFICIENTS
  !
  DO I = 1, NM1
    B(I) = (Y(I+1) - Y(I))/D(I) - D(I)*(C(I+1) + 2._fp*C(I))
    D(I) = (C(I+1) - C(I))/D(I)
    C(I) = 3._fp*C(I)
  end do
  ! periodic continuation...
  B(N) = B(1)
  C(N) = C(1)
  D(N) = D(1)
  !
  return
end subroutine pspline
