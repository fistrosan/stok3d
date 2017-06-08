MODULE SPLINES
  !
  ! J M Borrero
  ! July 23, 2013
  ! KIS, Freiburg
  ! Based on Numerical Recipes
  USE CONS_PARAM
CONTAINS
  ! splie2
  ! splin2
  ! spline
  ! splint
  !
  ! -----------------
  ! Subroutine splie2
  ! -----------------
  SUBROUTINE SPLIE2(X1A,X2A,YA,M,N,Y2A)
    IMPLICIT NONE
    INTEGER,       INTENT(IN)    :: M,N
    REAL(DP),      INTENT(IN)    :: X1A(M), X2A(N), YA(M,N)
    REAL(DP),      INTENT(INOUT) :: Y2A(M,N)
    !
    REAL(DP),          PARAMETER :: EPS=1.E30
    INTEGER,           PARAMETER :: NN=100
    REAL(DP)                     :: YTMP(NN),Y2TMP(NN)
    INTEGER                      :: J,K
    !
    DO J=1,M
       DO K=1,N
          YTMP(K)=YA(J,K)
       ENDDO
       CALL SPLINE(X2A,YTMP,N,EPS,EPS,Y2TMP)
       DO  K=1,N
          Y2A(J,K)=Y2TMP(K)
       ENDDO
    ENDDO
    !
    RETURN
    !
  END SUBROUTINE SPLIE2
  ! -----------------
  ! Subroutine splin2
  ! -----------------
  SUBROUTINE SPLIN2(X1A,X2A,YA,Y2A,M,N,X1,X2,Y)
    IMPLICIT NONE
    INTEGER,        INTENT(IN)  :: M,N
    REAL(DP),       INTENT(IN)  :: X1A(M),X2A(N),YA(M,N),Y2A(M,N),X1,X2 
    REAL(DP),       INTENT(OUT) :: Y
    !
    REAL(DP),         PARAMETER :: EPS=1.E30
    INTEGER,          PARAMETER :: NN=100
    REAL(DP)                    :: YTMP(NN),Y2TMP(NN),YYTMP(NN)
    INTEGER                     :: J,K
    !
    DO J=1,M
       DO K=1,N
          YTMP(K)=YA(J,K)
          Y2TMP(K)=Y2A(J,K)
       ENDDO
       CALL SPLINT(X2A,YTMP,Y2TMP,N,X2,YYTMP(J))
    ENDDO
    !
    CALL SPLINE(X1A,YYTMP,M,EPS,EPS,Y2TMP)
    CALL SPLINT(X1A,YYTMP,Y2TMP,M,X1,Y)
    !
  END SUBROUTINE SPLIN2
  ! -----------------
  ! Subroutine spline
  ! -----------------
  SUBROUTINE SPLINE(X,Y,N,YP1,YPN,Y2)
    IMPLICIT NONE
    INTEGER,        INTENT(IN)  :: N
    REAL(DP),       INTENT(IN)  :: X(N),Y(N), YP1, YPN
    REAL(DP),       INTENT(OUT) :: Y2(N)
    !
    INTEGER,          PARAMETER :: NMAX=100
    REAL(DP)                    :: U(NMAX), P, QN, UN, SIG
    INTEGER                     :: I,K
    !
    IF (YP1.GT..99E30) THEN
       Y2(1)=0D0
       U(1)=0D0
    ELSE
       Y2(1)=-0.5D0
       U(1)=(3D0/(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
    ENDIF
    !
    DO I=2,N-1
       SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
       P=SIG*Y2(I-1)+2D0
       Y2(I)=(SIG-1.)/P
       U(I)=(6D0*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1)) &
            /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
    ENDDO
    !
    IF (YPN.GT..99E30) THEN
       QN=0D0
       UN=0D0
    ELSE
       QN=0.5D0
       UN=(3D0/(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
    ENDIF
    !
    Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+1.)
    DO K=N-1,1,-1
       Y2(K)=Y2(K)*Y2(K+1)+U(K)
    ENDDO
    !
  END SUBROUTINE SPLINE
  !
  ! Subroutine splint
  !
  SUBROUTINE SPLINT(XA,YA,Y2A,N,X,Y)
    IMPLICIT NONE
    INTEGER,                 INTENT(IN)  :: N
    REAL(DP),                INTENT(IN)  :: XA(N),YA(N),Y2A(N),X
    REAL(DP),                INTENT(OUT) :: Y
    !
    REAL(DP)                             :: A,B,H
    INTEGER                              :: K,KLO, KHI
    !
    KLO=1
    KHI=N
    DO
       IF (KHI-KLO.LE.1) EXIT
       K=(KHI+KLO)/2
       IF(XA(K).GT.X)THEN
          KHI=K
       ELSE
          KLO=K
       ENDIF
    ENDDO
    H=XA(KHI)-XA(KLO)
    IF (H.EQ.0.) THEN
       PRINT*,'Error during spline interpolation in Subroutine SPLINT. STOP'
       STOP
    ENDIF
    A=(XA(KHI)-X)/H
    B=(X-XA(KLO))/H
    Y=A*YA(KLO)+B*YA(KHI)+((A**3D0-A)*Y2A(KLO)+(B**3D0-B)*Y2A(KHI))*(H**2D0)/6D0
  END SUBROUTINE SPLINT
  !
END MODULE SPLINES
