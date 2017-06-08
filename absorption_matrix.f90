MODULE ABSORPTION_MATRIX
  !
  ! May 20, 2014
  ! KIS, Freiburg
  !
  USE CONS_PARAM
  USE LINES_DATABASE
  USE SIMUL_BOX
  USE MAGSPLIT
  USE ATOM_DATABASE
  !
  REAL(SP),     ALLOCATABLE    :: ETAI(:,:,:,:), ETAQ(:,:,:,:), ETAU(:,:,:,:), ETAV(:,:,:,:)
  REAL(SP),     ALLOCATABLE    :: RHOQ(:,:,:,:), RHOU(:,:,:,:), RHOV(:,:,:,:)
  !
CONTAINS
  ! abs_mat_init
  ! get_abs_mat
  !
  !------------------------
  ! Subroutine ABS_MAT_INIT
  !------------------------
  SUBROUTINE ABS_MAT_INIT
    !
    IMPLICIT NONE
    INTEGER                        :: IERR
    !
    ALLOCATE(ETAI(NX,NY,NZ,NUMW), STAT=IERR)
    IF (IERR.NE.0) THEN
       PRINT*,'Not enought memory to allocate ETAI array. STOP.'
       STOP
    ENDIF
    !
    ALLOCATE(ETAQ(NX,NY,NZ,NUMW), STAT=IERR)
    IF (IERR.NE.0) THEN
       PRINT*,'Not enought memory to allocate ETAQ array. STOP.'
       STOP
    ENDIF
    !
    ALLOCATE(ETAU(NX,NY,NZ,NUMW), STAT=IERR)
    IF (IERR.NE.0) THEN
       PRINT*,'Not enought memory to allocate ETAU array. STOP.'
       STOP
    ENDIF
    !
    ALLOCATE(ETAV(NX,NY,NZ,NUMW), STAT=IERR)
    IF (IERR.NE.0) THEN
       PRINT*,'Not enought memory to allocate ETAV array. STOP.'
       STOP
    ENDIF
    !
    ALLOCATE(RHOQ(NX,NY,NZ,NUMW), STAT=IERR)
    IF (IERR.NE.0) THEN
       PRINT*,'Not enought memory to allocate RHOQ array. STOP.'
       STOP
    ENDIF
    !
    ALLOCATE(RHOU(NX,NY,NZ,NUMW), STAT=IERR)
    IF (IERR.NE.0) THEN
       PRINT*,'Not enought memory to allocate RHOU array. STOP.'
       STOP
    ENDIF
    !
    ALLOCATE(RHOV(NX,NY,NZ,NUMW), STAT=IERR)
    IF (IERR.NE.0) THEN
       PRINT*,'Not enought memory to allocate RHOV array. STOP.'
       STOP
    ENDIF
    !
    ! Set all arrays to zero
    !
    ETAI(:,:,:,:)=0.
    ETAQ(:,:,:,:)=0.
    ETAU(:,:,:,:)=0.
    ETAV(:,:,:,:)=0.
    RHOQ(:,:,:,:)=0.
    RHOU(:,:,:,:)=0.
    RHOV(:,:,:,:)=0.
    !
  END SUBROUTINE ABS_MAT_INIT
  !-----------------------
  ! Subroutine GET_ABS_MAT
  !-----------------------
  SUBROUTINE GET_ABS_MAT (I,J,K,L,DAMP)
    !
    IMPLICIT NONE
    !
    INTEGER,   INTENT(IN)        :: I, J, K, L
    REAL(SP),  INTENT(IN)        :: DAMP
    !
    INTEGER                      :: M, N, COUNT
    REAL(SP)                     :: VDOPPLER, DLDOPPLER, BFIELD
    REAL(SP)                     :: DELTAL_V, DELTAL_B, BBX, BBY, BBZ
    REAL(SP)                     :: COSGAMMA, COSGAMMAFAC, SINGAMMAFAC, COS2PHI, SIN2PHI
    REAL(DP), ALLOCATABLE        :: HH(:), FF(:), UU(:)
    LOGICAL                      :: NAN, CHECKNAN
    !
    ! Allocate arrays for Voigt and Faraday functions
    !
    ALLOCATE(HH(NUMWAVE(L)),FF(NUMWAVE(L)), UU(NUMWAVE(L)))
    !
    ! Doppler width of the spectral line
    ! Note:  VDOPPLER here is different from that of the module collisional_damping
    ! Here, the factor PI^(-1/2) does not appear. It is only a matter of definition
    !
    VDOPPLER = SQRT(2.*KBOL*TEM(I,J,K)/(MATOM(LINE_ZN(IND_LINE(L)))*MAMU)) ! cm/s
    DLDOPPLER = 1E3*(VDOPPLER*LINE_L0(IND_LINE(L))/LIGHT)                  ! mA
    !
    ! Magnetic field strength
    !
    BFIELD=SQRT(BX(I,J,K)**2.+BY(I,J,K)**2.+BZ(I,J,K)**2.)
    !
    ! Cartesian coordinates of the magnetic field
    !
    BBX=BX(I,J,K)
    BBY=BY(I,J,K)
    BBZ=BZ(I,J,K)
    !
    ! Factors in terms of the spherical coordinates of the magnetic field
    !
    COSGAMMA=0.
    COSGAMMAFAC=1.
    SINGAMMAFAC=1.
    COS2PHI=0.
    SIN2PHI=0.
    IF (BFIELD.GT.0) THEN
       COSGAMMA=BBZ/BFIELD                              ! COS(GAMMA)
       COSGAMMAFAC=1+BBZ**2./BFIELD**2.                 ! 1+COS(GAMMA)^2
       SINGAMMAFAC=1-BBZ**2./BFIELD**2.                 ! SIN(GAMMA)^2
    ENDIF
    IF (BBX.NE.0.OR.BBY.NE.0) THEN
       COS2PHI=(BBY**2.-BBX**2.)/(BBX**2.+BBY**2.)      ! COS(2PHI)
       SIN2PHI=2.*BBX*BBY/(BBX**2.+BBY**2.)             ! SIN(2PHI)
    ENDIF
    !
    ! Checking for NaN
    !
    NAN=.FALSE.
    NAN=CHECKNAN(DBLE(COSGAMMA))
    NAN=CHECKNAN(DBLE(COSGAMMAFAC))
    NAN=CHECKNAN(DBLE(SINGAMMAFAC))
    NAN=CHECKNAN(DBLE(COS2PHI))
    NAN=CHECKNAN(DBLE(SIN2PHI))
    IF (NAN.EQV..TRUE.) THEN
       PRINT*,'NaN detected in gamma,phi factors'
       PRINT*,COSGAMMA,COSGAMMAFAC,SINGAMMAFAC,COS2PHI,SIN2PHI
       STOP
    ENDIF
    !
    ! Initialize Zeeman pattern for spectral line L
    !
    CALL ZEEMAN_COMP(L)
    !
    ! Start loop over components of the Zeeman pattern
    !
    DO M=1,NTOT
       COUNT=1
       UU(:)=0D0
       DO N=PIXEL_INI(L),PIXEL_END(L)
          DELTAL_B=LAMBDA_SC(M)*BFIELD                       ! mA
          DELTAL_V=-1E3*LINE_L0(IND_LINE(L))*VZ(I,J,K)/LIGHT ! mA (if Vz in cm/s); Vz < 0 is a blueshift
          UU(COUNT)=DBLE((WAVE(N)+DELTAL_B+DELTAL_V)/DLDOPPLER)
          IF (ABS(UU(COUNT)).GT.1E4) THEN
             PRINT*,UU(COUNT),WAVE(N),DELTAL_B,DELTAL_V,DLDOPPLER
             PRINT*,LINE_L0(IND_LINE(L)),VZ(I,J,K)
             STOP
          ENDIF
          NAN=.FALSE.
          NAN=CHECKNAN(DBLE(UU(COUNT)))
          IF (NAN.EQV..TRUE.) THEN
             PRINT*,'NaN detected in UU'
             PRINT*,UU(COUNT),WAVE(N),DELTAL_B,DELTAL_V,DLDOPPLER
             STOP
          ENDIF
          COUNT=COUNT+1
       ENDDO
       NAN=.FALSE.
       NAN=CHECKNAN(DBLE(STRENGTH_SC(M)))
       IF (NAN.EQV..TRUE.) THEN
          PRINT*,'NaN detected in STRENGTH_SC'
          PRINT*,STRENGTH_SC(M),TYPE_SC(M),LAMBDA_SC(M)
          STOP
       ENDIF
       !
       ! Determine Voigt and Faraday functions
       !
       CALL VOIGT(NUMWAVE(L),DBLE(DAMP),UU,HH,FF)
       !
       ! Add contribution to the different elements of the absorption matrix
       !
       SELECT CASE(TYPE_SC(M))
       CASE(-1) ! SIGMA_RED
          ETAI(I,J,K,PIXEL_INI(L):PIXEL_END(L))=ETAI(I,J,K,PIXEL_INI(L):PIXEL_END(L))+STRENGTH_SC(M)*0.25*HH(:)*COSGAMMAFAC
          ETAV(I,J,K,PIXEL_INI(L):PIXEL_END(L))=ETAV(I,J,K,PIXEL_INI(L):PIXEL_END(L))+STRENGTH_SC(M)*0.5*HH(:)*COSGAMMA
          ETAQ(I,J,K,PIXEL_INI(L):PIXEL_END(L))=ETAQ(I,J,K,PIXEL_INI(L):PIXEL_END(L))-STRENGTH_SC(M)*0.25*HH(:)*SINGAMMAFAC*COS2PHI
          ETAU(I,J,K,PIXEL_INI(L):PIXEL_END(L))=ETAU(I,J,K,PIXEL_INI(L):PIXEL_END(L))-STRENGTH_SC(M)*0.25*HH(:)*SINGAMMAFAC*SIN2PHI
          RHOV(I,J,K,PIXEL_INI(L):PIXEL_END(L))=RHOV(I,J,K,PIXEL_INI(L):PIXEL_END(L))+STRENGTH_SC(M)*0.5*FF(:)*COSGAMMA
          RHOQ(I,J,K,PIXEL_INI(L):PIXEL_END(L))=RHOQ(I,J,K,PIXEL_INI(L):PIXEL_END(L))-STRENGTH_SC(M)*0.25*FF(:)*SINGAMMAFAC*COS2PHI
          RHOU(I,J,K,PIXEL_INI(L):PIXEL_END(L))=RHOU(I,J,K,PIXEL_INI(L):PIXEL_END(L))-STRENGTH_SC(M)*0.25*FF(:)*SINGAMMAFAC*SIN2PHI
       CASE(0) ! PI
          ETAI(I,J,K,PIXEL_INI(L):PIXEL_END(L))=ETAI(I,J,K,PIXEL_INI(L):PIXEL_END(L))+STRENGTH_SC(M)*0.5*HH(:)*SINGAMMAFAC
          ETAQ(I,J,K,PIXEL_INI(L):PIXEL_END(L))=ETAQ(I,J,K,PIXEL_INI(L):PIXEL_END(L))+STRENGTH_SC(M)*0.5*HH(:)*SINGAMMAFAC*COS2PHI
          ETAU(I,J,K,PIXEL_INI(L):PIXEL_END(L))=ETAU(I,J,K,PIXEL_INI(L):PIXEL_END(L))+STRENGTH_SC(M)*0.5*HH(:)*SINGAMMAFAC*SIN2PHI
          RHOQ(I,J,K,PIXEL_INI(L):PIXEL_END(L))=RHOQ(I,J,K,PIXEL_INI(L):PIXEL_END(L))+STRENGTH_SC(M)*0.5*FF(:)*SINGAMMAFAC*COS2PHI
          RHOU(I,J,K,PIXEL_INI(L):PIXEL_END(L))=RHOU(I,J,K,PIXEL_INI(L):PIXEL_END(L))+STRENGTH_SC(M)*0.5*FF(:)*SINGAMMAFAC*SIN2PHI
       CASE(1) ! SIGMA_BLUE
          ETAI(I,J,K,PIXEL_INI(L):PIXEL_END(L))=ETAI(I,J,K,PIXEL_INI(L):PIXEL_END(L))+STRENGTH_SC(M)*0.25*HH(:)*COSGAMMAFAC
          ETAV(I,J,K,PIXEL_INI(L):PIXEL_END(L))=ETAV(I,J,K,PIXEL_INI(L):PIXEL_END(L))-STRENGTH_SC(M)*0.5*HH(:)*COSGAMMA
          ETAQ(I,J,K,PIXEL_INI(L):PIXEL_END(L))=ETAQ(I,J,K,PIXEL_INI(L):PIXEL_END(L))-STRENGTH_SC(M)*0.25*HH(:)*SINGAMMAFAC*COS2PHI
          ETAU(I,J,K,PIXEL_INI(L):PIXEL_END(L))=ETAU(I,J,K,PIXEL_INI(L):PIXEL_END(L))-STRENGTH_SC(M)*0.25*HH(:)*SINGAMMAFAC*SIN2PHI
          RHOV(I,J,K,PIXEL_INI(L):PIXEL_END(L))=RHOV(I,J,K,PIXEL_INI(L):PIXEL_END(L))-STRENGTH_SC(M)*0.5*FF(:)*COSGAMMA
          RHOQ(I,J,K,PIXEL_INI(L):PIXEL_END(L))=RHOQ(I,J,K,PIXEL_INI(L):PIXEL_END(L))-STRENGTH_SC(M)*0.25*FF(:)*SINGAMMAFAC*COS2PHI
          RHOU(I,J,K,PIXEL_INI(L):PIXEL_END(L))=RHOU(I,J,K,PIXEL_INI(L):PIXEL_END(L))-STRENGTH_SC(M)*0.25*FF(:)*SINGAMMAFAC*SIN2PHI
       END SELECT
       HH(:)=0D0
       FF(:)=0D0
    ENDDO
    ! Add line and continuum contribution
    ETAI(I,J,K,PIXEL_INI(L):PIXEL_END(L))=KC(I,J,K,L)+KLIN(I,J,K,L)*ETAI(I,J,K,PIXEL_INI(L):PIXEL_END(L))
    ETAV(I,J,K,PIXEL_INI(L):PIXEL_END(L))=KLIN(I,J,K,L)*ETAV(I,J,K,PIXEL_INI(L):PIXEL_END(L))
    ETAQ(I,J,K,PIXEL_INI(L):PIXEL_END(L))=KLIN(I,J,K,L)*ETAQ(I,J,K,PIXEL_INI(L):PIXEL_END(L))
    ETAU(I,J,K,PIXEL_INI(L):PIXEL_END(L))=KLIN(I,J,K,L)*ETAU(I,J,K,PIXEL_INI(L):PIXEL_END(L))
    RHOV(I,J,K,PIXEL_INI(L):PIXEL_END(L))=KLIN(I,J,K,L)*RHOV(I,J,K,PIXEL_INI(L):PIXEL_END(L))
    RHOQ(I,J,K,PIXEL_INI(L):PIXEL_END(L))=KLIN(I,J,K,L)*RHOQ(I,J,K,PIXEL_INI(L):PIXEL_END(L))
    RHOU(I,J,K,PIXEL_INI(L):PIXEL_END(L))=KLIN(I,J,K,L)*RHOU(I,J,K,PIXEL_INI(L):PIXEL_END(L))
    !
    !PRINT*,MAXVAL(ETAI(I,J,K,PIXEL_INI(L):PIXEL_END(L))),MAXVAL(ETAQ(I,J,K,PIXEL_INI(L):PIXEL_END(L))), &
    !     MAXVAL(ETAU(I,J,K,PIXEL_INI(L):PIXEL_END(L))),MAXVAL(ETAV(I,J,K,PIXEL_INI(L):PIXEL_END(L)))
    !PRINT*,'-------------'
    !PAUSE
    !
    DEALLOCATE(HH, FF, UU)
    DEALLOCATE(TYPE_SC, LAMBDA_SC, STRENGTH_SC)
    !
  END SUBROUTINE GET_ABS_MAT





END MODULE ABSORPTION_MATRIX
