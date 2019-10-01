!******************** START FILE SPEVAL.F90 ; GROUP TRKRLIB ************
!....................................................
!CCCCCCCCCCCCCC
real(fp) function speval(N, U, X, Y, B, C, D)
  use precision_mod, only: fp
  implicit none
  integer N
  real(fp) ::  U, X(N), Y(N), B(N), C(N), D(N)
  !
  !CCCCCCCCCCCCCC
  !CCCCCCCCCCCCCC
  !  THIS subroutine EVALUATES THE DERIVATIVE OF THE CUBIC SPLINE FUNCTION
  !
  !
  !    WHERE  X(I) .LT. U .LT. X(I+1), USING HORNER'S RULE
  !
  !  IF  U .LT. X(1) THEN  I = 1  IS USED.
  !  IF  U .GE. X(N) THEN  I = N  IS USED.
  !
  !  INPUT..
  !
  !    N = THE NUMBER OF DATA POINTS
  !    U = THE ABSCISSA AT WHICH THE SPLINE IS TO BE EVALUATED
  !    X,Y = THE ARRAYS OF DATA ABSCISSAS AND ORDINATES
  !    B,C,D = ARRAYS OF SPLINE COEFFICIENTS COMPUTED BY SPLINE
  !
  !  IF  U  IS NOT IN THE SAME INTERVAL AS THE PREVIOUS CALL, THEN A
  !  BINARY SEARCH IS PERFORMED TO DETERMINE THE PROPER INTERVAL.
  !
  integer I, J, K
  real(fp) :: DX
  I = 1
  IF ( I .GE. N ) I = 1
  IF ( U .LT. X(I) ) GO TO 10
  IF ( U .LE. X(I+1) ) GO TO 30
  !
  !  BINARY SEARCH
  !
10 I = 1
  J = N+1
20 K = (I+J)/2
  IF ( U .LT. X(K) ) J = K
  IF ( U .GE. X(K) ) I = K
  IF ( J .GT. I+1 ) GO TO 20
  !
  !  EVALUATE SPLINE
  !
30 DX = U - X(I)
  speval = B(I) + DX*(2._fp*C(I) + 3._fp*DX*D(I))
  RETURN
END FUNCTION speval
!******************** END FILE SPEVAL.F90 ; GROUP TRKRLIB **************
