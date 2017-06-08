C     ***CUBIC************************************************08.11.1986
C     Solution of a cubic equation
C     Equations of lesser degree are solved by the appropriate formulas.
C     The solutions are arranged in ascending order.
C     NO WARRANTY, ALWAYS TEST THIS SUBROUTINE AFTER DOWNLOADING
C     ******************************************************************
C     A(0:3)      (i)  vector containing the polynomial coefficients
C     X(1:L)      (o)  results
C     L           (o)  number of valid solutions (beginning with X(1))
C     ==================================================================
      SUBROUTINE CUBIC(A,X,L)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(0:3),X(3),U(3)
      PARAMETER(PI=3.1415926535897932D+0,THIRD=1.D+0/3.D+0)
      INTRINSIC MIN,MAX,ACOS
C
C     define cubic root as statement function
      CBRT(Z)=SIGN(ABS(Z)**THIRD,Z)
C
C     ==== determine the degree of the polynomial ====
C
      IF (A(3).NE.0.D+0) THEN
C
C       cubic problem
	W=A(2)/A(3)*THIRD
	P=(A(1)/A(3)*THIRD-W**2)**3
	Q=-.5D+0*(2.D+0*W**3-(A(1)*W-A(0))/A(3))
	DIS=Q**2+P
	IF (DIS.LT.0.D+0) THEN
C         three real solutions!
C         Confine the argument of ACOS to the interval [-1;1]!
	  PHI=ACOS(MIN(1.D+0,MAX(-1.D+0,Q/SQRT(-P))))
	  P=2.D+0*(-P)**(5.D-1*THIRD)
	  DO 100 I=1,3
  100       U(I)=P*COS((PHI+DBLE(2*I)*PI)*THIRD)-W
	  X(1)=MIN(U(1),U(2),U(3))
	  X(2)=MAX(MIN(U(1),U(2)),MIN(U(1),U(3)),MIN(U(2),U(3)))
	  X(3)=MAX(U(1),U(2),U(3))
	  L=3
	ELSE
C         only one real solution!
	  DIS=SQRT(DIS)
	  X(1)=CBRT(Q+DIS)+CBRT(Q-DIS)-W
	  L=1
	END IF
C
      ELSE IF (A(2).NE.0.D+0) THEN
C
C       quadratic problem
	P=5.D-1*A(1)/A(2)
	DIS=P**2-A(0)/A(2)
	IF (DIS.GE.0.D+0) THEN
C         two real solutions!
	  X(1)=-P-SQRT(DIS)
	  X(2)=-P+SQRT(DIS)
	  L=2
	ELSE
C         no real solution!
	  L=0
	END IF
C
      ELSE IF (A(1).NE.0.D+0) THEN
C
C       linear equation
	X(1)=-A(0)/A(1)
	L=1
C
      ELSE
C       no equation
	L=0
      END IF
C
C     ==== perform one step of a newton iteration in order to minimize
C          round-off errors ====
      DO 110 I=1,L
	X(I)=X(I)-(A(0)+X(I)*(A(1)+X(I)*(A(2)+X(I)*A(3))))
     *  /(A(1)+X(I)*(2.D+0*A(2)+X(I)*3.D+0*A(3)))
  110 CONTINUE
      RETURN
      END
