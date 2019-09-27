!******************** START FILE SEVAL2.F90 ; GROUP TRKRLIB ************
!
real(fp) FUNCTION seval2(N, U, X, Y, B, C, D, deriv)
  use precision_mod, only: fp
  !============
  implicit none
  real(fp) :: zslop
  !============
  integer N
  real(fp) ::  U, X(N), Y(N), B(N), C(N), D(N), deriv
  !
  !CCCCCCCCCCCCCC
  !CCCCCCCCCCCCCC
  !  THIS subroutine EVALUATES THE CUBIC SPLINE FUNCTION
  !
  !    SEVAL2 = Y(I) + B(I)*(U-X(I)) + C(I)*(U-X(I))**2 + D(I)*(U-X(I))**3
  !
  !  and the derivative
  !
  !    deriv = B(i) + 2*C(i)*(U-X(i)) + 3*D(i)*(U-X(i))**2
  !
  !  WHERE  X(I) .LT. U .LT. X(I+1), USING HORNER'S RULE
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
  integer I, J, K, ILIN
  real(fp) :: DX
  DATA I/1/
  DATA ILIN/0/
  !
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
  IF(ILIN.EQ.0) THEN
    seval2 = Y(I) + DX*(B(I) + DX*(C(I) + DX*D(I)))
    deriv = B(i) + DX*(2.0_fp*C(i) + DX*3.0_fp*D(i))
  ELSE
    IF(I.EQ.N) THEN
      ZSLOP=(Y(N)-Y(N-1))/(X(N)-X(N-1))
    ELSE
      ZSLOP=(Y(I+1)-Y(I))/(X(I+1)-X(I))
    end if
    seval2=Y(I)+DX*ZSLOP
    deriv=ZSLOP
  end if
  !
  RETURN
END FUNCTION seval2
!******************** END FILE SEVAL2.F90 ; GROUP TRKRLIB **************
