MODULE RTESOLVER
  !
  ! May 22, 2014
  ! KIS, Freiburg
  !
  USE CONS_PARAM
  USE SIMUL_BOX
  USE ABSORPTION_MATRIX
  !
  IMPLICIT NONE
  !
CONTAINS
  ! analytical_solver
  ! evol_operator
  ! source function
  ! source function_der
  ! set_boundary
  !-----------------------------
  ! Subroutine analytical_solver
  !-----------------------------
  SUBROUTINE ANALYTICAL_SOLVER(I,J,K,L)
    !
    IMPLICIT NONE
    !
    INTEGER,    INTENT(IN)       :: I,J,K,L
    REAL(SP), DIMENSION(4,4)     :: EVOL
    REAL(SP), DIMENSION(4)       :: SF
    REAL(SP), DIMENSION(4,1)     :: SFM, IM, IN
    REAL(SP)                     :: DSF
    ! For derivatives
    INTEGER                      :: M
    REAL(SP), DIMENSION(NZ,4,5)  :: DER
    !REAL(SP), DIMENSION(4,1)     :: DIX
    REAL(SP), DIMENSION(4,4)     :: ABSMAT
    !
    CALL EVOL_OPERATOR(I,J,K,L,EVOL)
    !
    CALL SOURCE_FUNCTION(I,J,K,L,SF)
    SFM(:,1)=SF(:)
    !
    IN(:,1)=STOKES(I,J,K-1,L,:)
    IM=MATMUL(IMAT-EVOL,SFM)+MATMUL(EVOL,IN)
    !
    STOKES(I,J,K,L,:)=IM(:,1)
    !
  END SUBROUTINE ANALYTICAL_SOLVER
  !-----------------------------
  ! Subroutine evol_operator
  !-----------------------------
  SUBROUTINE EVOL_OPERATOR(I,J,K,L,EVOL)
    !
    IMPLICIT NONE
    !
    INTEGER,   INTENT(IN)                     :: I,J,K,L
    REAL(SP), INTENT(OUT), DIMENSION(4,4)     :: EVOL
    REAL(DP)                                  :: ETA_I, ETA_Q, ETA_U, ETA_V
    REAL(DP)                                  :: RHO_Q, RHO_U, RHO_V
    REAL(DP)                                  :: ETA2, RHO2, ETARHO, SIGMA, THETA, RPART, EXPETAI
    REAL(DP)                                  :: LAMBDA1, LAMBDA2, COSHL1, SINHL1, COSL2, SINL2
    REAL(SP), DIMENSION(4,4)                  :: M1, M2, M3, M4
    INTEGER                                   :: M
    LOGICAL                                   :: NAN, CHECKNAN
    !
    M1(:,:)=0.
    M2(:,:)=0.
    M3(:,:)=0.
    M4(:,:)=0.
    SIGMA=0.
    !
    ETA_I=DBLE(ETAI(I,J,K,L))
    ETA_Q=DBLE(ETAQ(I,J,K,L))
    ETA_U=DBLE(ETAU(I,J,K,L))
    ETA_V=DBLE(ETAV(I,J,K,L))
    RHO_Q=DBLE(RHOQ(I,J,K,L))
    RHO_U=DBLE(RHOU(I,J,K,L))
    RHO_V=DBLE(RHOV(I,J,K,L))
    ! Checking for NaN
    NAN=.FALSE.
    NAN=CHECKNAN(ETA_I)
    NAN=CHECKNAN(ETA_Q)
    NAN=CHECKNAN(ETA_U)
    NAN=CHECKNAN(ETA_V)
    NAN=CHECKNAN(RHO_Q)
    NAN=CHECKNAN(RHO_U)
    NAN=CHECKNAN(RHO_V)
    IF (NAN.EQV..TRUE.) THEN
       PRINT*,'NaN detected in EVOL_OPERATOR'
       PRINT*,ETA_I,ETA_Q,ETA_U,ETA_V,RHO_Q,RHO_U,RHO_V
       PRINT*,TEM(I,J,K),RHO(I,J,K),PG(I,J,K),KC(I,J,K,:),KLIN(I,J,K,:)
       PRINT*,I,J,K,L
       STOP
    ENDIF
    !
    ! This part follows Landi Degl'Innocenti & Landi Degl'Innocenti
    ! Solar Physics, 1985, 239-250
    !
    ETA2=ETA_Q**2.+ETA_U**2.+ETA_V**2.
    RHO2=RHO_Q**2.+RHO_U**2.+RHO_V**2.
    ETARHO=ETA_Q*RHO_Q+ETA_U*RHO_U+ETA_V*RHO_V
    RPART=SQRT(((ETA2-RHO2)**2.)/4.+ETARHO**2.)+(ETA2-RHO2)/2.
    LAMBDA1=SQRT(RPART)
    !
    RPART=SQRT(((ETA2-RHO2)**2.)/4.+ETARHO**2.)-(ETA2-RHO2)/2.
    LAMBDA2=SQRT(RPART)
    IF (ETARHO.GT.0) SIGMA=1.
    IF (ETARHO.LT.0) SIGMA=-1.
    THETA=2.*SQRT(((ETA2-RHO2)**2.)/4.+ETARHO**2.)
    ! M1 matrix is the indentiy matrix
    M1=IMAT
    ! M2 matrix
    M2(1,2)=LAMBDA2*ETA_Q-SIGMA*LAMBDA1*RHO_Q
    M2(1,3)=LAMBDA2*ETA_U-SIGMA*LAMBDA1*RHO_U
    M2(1,4)=LAMBDA2*ETA_V-SIGMA*LAMBDA1*RHO_V
    M2(2,1)=LAMBDA2*ETA_Q-SIGMA*LAMBDA1*RHO_Q
    M2(2,3)=SIGMA*LAMBDA1*ETA_V+LAMBDA2*RHO_V
    M2(2,4)=-SIGMA*LAMBDA1*ETA_U-LAMBDA2*RHO_U
    M2(3,1)=LAMBDA2*ETA_U-SIGMA*LAMBDA1*RHO_U
    M2(3,2)=-SIGMA*LAMBDA1*ETA_V-LAMBDA2*RHO_V
    M2(3,4)=SIGMA*LAMBDA1*ETA_Q+LAMBDA2*RHO_Q
    M2(4,1)=LAMBDA2*ETA_V-SIGMA*LAMBDA1*RHO_V
    M2(4,2)=-SIGMA*LAMBDA1*ETA_U+LAMBDA2*RHO_U
    M2(4,3)=-SIGMA*LAMBDA1*ETA_Q-LAMBDA1*RHO_Q
    M2(:,:)=M2(:,:)/THETA
    ! M3 matrix
    M3(1,2)=LAMBDA1*ETA_Q+SIGMA*LAMBDA2*RHO_Q
    M3(1,3)=LAMBDA1*ETA_U+SIGMA*LAMBDA2*RHO_U
    M3(1,4)=LAMBDA1*ETA_V+SIGMA*LAMBDA2*RHO_V
    M3(2,1)=LAMBDA1*ETA_Q+SIGMA*LAMBDA2*RHO_Q
    M3(2,3)=-SIGMA*LAMBDA2*ETA_V+LAMBDA1*RHO_V
    M3(2,4)=SIGMA*LAMBDA2*ETA_U-LAMBDA1*RHO_U
    M3(3,1)=LAMBDA1*ETA_U+SIGMA*LAMBDA2*RHO_U
    M3(3,2)=SIGMA*LAMBDA2*ETA_V-LAMBDA1*RHO_V
    M3(3,4)=-SIGMA*LAMBDA2*ETA_Q+LAMBDA1*RHO_Q
    M3(4,1)=LAMBDA1*ETA_V+SIGMA*LAMBDA2*RHO_V
    M3(4,2)=-SIGMA*LAMBDA2*ETA_U+LAMBDA1*RHO_U
    M3(4,3)=SIGMA*LAMBDA2*ETA_Q-LAMBDA1*RHO_Q
    M3(:,:)=M3(:,:)/THETA
    ! M4 matrix
    M4(1,1)=(ETA2+RHO2)/2.
    M4(1,2)=ETA_V*RHO_U-ETA_U*RHO_V
    M4(1,3)=ETA_Q*RHO_V-ETA_V*RHO_Q
    M4(1,4)=ETA_U*RHO_Q-ETA_Q*RHO_U
    M4(2,1)=ETA_U*RHO_V-ETA_V*RHO_U
    M4(2,2)=ETA_Q**2.+RHO_Q**2.-(ETA2+RHO2)/2.
    M4(2,3)=ETA_Q*ETA_U+RHO_Q*RHO_U
    M4(2,4)=ETA_V*ETA_Q+RHO_V*RHO_Q
    M4(3,1)=ETA_V*RHO_Q-ETA_Q*RHO_V
    M4(3,2)=ETA_Q*ETA_U+RHO_Q*RHO_U
    M4(3,3)=ETA_U**2.+RHO_U**2.-(ETA2+RHO2)/2.
    M4(3,4)=ETA_U*ETA_V+RHO_U*RHO_V
    M4(4,1)=ETA_Q*RHO_U-ETA_U*RHO_Q
    M4(4,2)=ETA_V*ETA_Q+RHO_V*RHO_Q
    M4(4,3)=ETA_U*ETA_V+RHO_U*RHO_V
    M4(4,4)=ETA_V**2.+RHO_V**2.-(ETA2+RHO2)/2.
    M4(:,:)=2.*M4(:,:)/THETA
    ! Multiplicative factors
    COSHL1=COSH(LAMBDA1*1E5*ABS(ZZ(K+1)-ZZ(K)))
    SINHL1=SINH(LAMBDA1*1E5*ABS(ZZ(K+1)-ZZ(K)))
    COSL2=COS(LAMBDA2*1E5*ABS(ZZ(K+1)-ZZ(K)))
    SINL2=SIN(LAMBDA2*1E5*ABS(ZZ(K+1)-ZZ(K)))
    EXPETAI=EXP(-ETA_I*1E5*ABS(ZZ(K+1)-ZZ(K)))
    ! Evolution operator
    EVOL=EXPETAI*(0.5*(COSHL1+COSL2)*M1-SINL2*M2-SINHL1*M3+0.5*(COSHL1-COSL2)*M4)
    !    
  END SUBROUTINE EVOL_OPERATOR
  !-----------------------------------
  ! Subroutine SOURCE_FUNCTION
  !-----------------------------------
  SUBROUTINE SOURCE_FUNCTION(I,J,K,L,SF)
    !
    IMPLICIT NONE
    !
    INTEGER,    INTENT(IN)                  :: I,J,K,L
    REAL(SP), INTENT(OUT), DIMENSION(4)     :: SF
    INTEGER                                 :: M
    REAL(DP)                                :: WAVELENGTH
    !
    DO M=1,NUML
       IF (L.GE.PIXEL_INI(M).AND.M.LE.PIXEL_END(M)) WAVELENGTH=LINE_L0(IND_LINE(M))+WAVE(L)/1E3
    ENDDO
    ! cm !
    WAVELENGTH=WAVELENGTH*1E-8
    ! We assume black body radiation (LTE) and no polarization for the source function
    SF(:)=0.
    SF(1)=2.*HPLA*LIGHT**2./(WAVELENGTH**5.)*(1./(EXP(HPLA*LIGHT/(WAVELENGTH*KBOL*TEM(I,J,K)))-1.))
    !
  END SUBROUTINE SOURCE_FUNCTION
  !-------------------------------
  ! Subroutine SOURCE_FUNCTION_DER
  !-------------------------------
  SUBROUTINE SOURCE_FUNCTION_DER(I,J,K,L,DSF)
    ! Calculates the derivative of the source function with respect to the temperature
    !
    IMPLICIT NONE
    !
    INTEGER,    INTENT(IN)                  :: I,J,K,L
    REAL(SP), INTENT(OUT)                   :: DSF
    INTEGER                                 :: M
    REAL(DP)                                :: WAVELENGTH
    REAL(SP), PARAMETER                     :: FACTOR = 8.5685641E-6 ! h^2*c^3/K
    !
    DO M=1,NUML
       IF (L.GE.PIXEL_INI(M).AND.M.LE.PIXEL_END(M)) WAVELENGTH=LINE_L0(IND_LINE(M))+WAVE(L)/1E3
    ENDDO
    ! cm !
    WAVELENGTH=WAVELENGTH*1E-8
    !
    DSF = (2.*FACTOR/(WAVELENGTH**6.*TEM(I,J,K)**2.))*(EXP(2.*HPLA*LIGHT/(WAVELENGTH*KBOL*TEM(I,J,K))))**2.

!/(EXP(HPLA*LIGHT/(WAVELENGTH*KBOL*TEM(I,J,K)))-1.))*2.
   PRINT*,DSF
    STOP
    !
  END SUBROUTINE SOURCE_FUNCTION_DER
  !-----------------------------------
  ! Subroutine SET_BOUDARY
  !-----------------------------------
  SUBROUTINE SET_BOUNDARY(I,J,K,L)
    !
    IMPLICIT NONE
    !
    INTEGER,    INTENT(IN)                  :: I,J,K,L
    REAL(SP),   DIMENSION(4)                :: SF
    ! At the lower boundary we assume the emergent radiation is the source function
    CALL SOURCE_FUNCTION(I,J,K,L,SF)
    STOKES(I,J,K,L,:)=SF
    !
  END SUBROUTINE SET_BOUNDARY
  !  

END MODULE RTESOLVER
