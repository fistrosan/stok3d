PROGRAM TEST
  ! J M Borrero
  ! June 17, 2014
  ! KIS, Freiburg
  USE LOG
  USE TIME_PARAM
  USE CONS_PARAM
  USE LINES_DATABASE
  USE ATOM_DATABASE
  USE SIMUL_BOX
  USE CHEMICAL_EQUILIBRIUM
  USE DAMPING
  USE ABSORPTION_MATRIX
  USE RTESOLVER
  USE BACKGROUND_OPACITY_MODULE
  USE DERIV_TEST
  USE DERIVVAR
  !
  IMPLICIT NONE
  !
  INTEGER                          :: I, J, K, L, M
  REAL(DP)                         :: NHYD, NELEC, PELEC
  REAL(DP)                         :: N, NI, NII, NIII, N_PI
  REAL(DP)                         :: DNELEC, DNI, DNII, DNIII, DN_PI
  REAL(DP)                         :: N2, NI2, NII2, NIII2, DNI2, DNII2, DNIII2
  REAL(SP)                         :: UI, UII, UIII, U_PI
  REAL(SP)                         :: DUI, DUII, DUIII, DU_PI
  REAL(SP)                         :: DAMP, MU
  LOGICAL                          :: NAN, CHECKNAN
  !
  ! Initialize logfile
  CALL LOG_INIT
  ! Initialize clock
  CALL SYSTEM_CLOCK(COUNT=START,COUNT_RATE=RATE,COUNT_MAX=MAX)
  ! Initialize atomic database
  CALL ATOM_INIT
  ! Read atomic line data base
  CALL READ_LINES_DATABASE
  ! Read box data and data of spectral lines to be used
  CALL READ_BOXFILE
  ! Initialize wavelength array
  CALL WAVE_INIT
  ! Initialize simulation box
  CALL BOX_INIT
  ! Initialize absorption matrix
  CALL ABS_MAT_INIT
  ! Initialize null iteration physical parameters
  CALL ZERO_ITER_HYDROSTATIC
  ! Test derivatives (?)
  !CALL DERIV_NE_RHO
  !CALL DERIV_NE_PG
  !CALL DERIV_DAMP
  CALL DERIV_NIR
  STOP
  !-----------------------------------------------------------------------------------
  ! Determine the line and continuum absorption coefficient, and the absorption matrix
  ! at every grid point in the box
  !-----------------------------------------------------------------------------------
  DO I=1,NX
     DO J=1,NY
        DO K=1,NZ
           !-------------------------------------------------------------------------------------
           ! Calculate number of hydrogen atoms per cm^3, electron pressure and molecular weight
           ! assuming temperature, gas pressure and density are known.
           ! If iter=0 then we already did this through hydrostatic equilibrium.
           ! If iter > 0 then we assume T, PG and RHO are given by the MHS equations.
           !-------------------------------------------------------------------------------------
           CALL GET_NHNEMU(TEM(I,J,K),PG(I,J,K),RHO(I,J,K),NHYD,NELEC,MU)
           MW(I,J,K)=MU
           PEL(I,J,K)=NELEC*KBOL*TEM(I,J,K)
           ! Calculate the derivative of the number of hydrogen atoms per cm^3 with respect
           ! to temperature at a) constant density and b) constant gas pressure.
           CALL GET_DNHYDDTEMP(TEM(I,J,K),PG(I,J,K),RHO(I,J,K),MW(I,J,K))
           ! Determening wavelength independent part of the line absorption coefficient.
           DO L=1,NUML
              ! Number of atoms per cm^3 of the specie whose spectral line we are interested in
              N=NHYD*10D0**(ABUND(LINE_ZN(IND_LINE(L)))-12D0)
              !print*,nhyd,n
              ! Determine line absorption coefficient and its derivatives with respect to T
              ! at constant gas pressure and constant density despectively.
              !print*,LINE_ZN(IND_LINE(L))
              CALL GET_KLIN(IND_LINE(L),TEM(I,J,K),NELEC,N,KLIN(I,J,K,L))
              STOP
              ! Determine number of neutrals, single and double ionized atoms and their derivatives with respect
              ! to the temperature at constant density and gas pressure, respectively.
              !CALL ION_CALC_TEMPX(LINE_ZN(IND_LINE(L)),LINE_ION(IND_LINE(L)),TEM(I,J,K),NELEC,N,NI,DNIDT_PG,DNIDT_RHO)
              ! Number of atoms per cm^3 of a given specie with different ionization stages
              !CALL PARTITION_FUNCTION(LINE_ZN(IND_LINE(L)),TEM(I,J,K),UI,UII,UIII,DUI,DUII,DUIII)
              !SELECT CASE (LINE_ION(IND_LINE(L)))
              !   CASE(1)
              !      N_PI = NI
              !      U_PI = UI
              !      DN_PI = DNI
              !      DU_PI = DUI
              !   CASE(2)
              !      N_PI = NII
              !      U_PI = UII
              !      DN_PI= DNII
              !      DU_PI = DUII
              !   CASE(3)
              !      N_PI = NIII
              !      U_PI = UIII
              !      DN_PI = DNIII
              !      DU_PI = DUIII
              !END SELECT
              ! This is the wavelength independet part of the absorption coefficient
              ! It depends on the spectral line however
              KLIN(I,J,K,L)=(SQRT(SPI)*QELE**2./(MELE*LIGHT))*(N_PI/U_PI)* &
                   EXP(-EPLOW(IND_LINE(L))*EVOLT/(KBOL*TEM(I,J,K)))*10.**(LOGGF(IND_LINE(L)))* &
                   (1.-EXP(-HPLA*LIGHT/(1E-8*LINE_L0(IND_LINE(L))*KBOL*TEM(I,J,K))))*LINE_L0(IND_LINE(L))*1E-8* &
                   SQRT(MATOM(LINE_ZN(IND_LINE(L)))*MAMU/(2.*KBOL*TEM(I,J,K)))
              ! Determine damping: radiative + colissional
              NAN=.FALSE.
              CALL GET_DAMPING(LINE_ZN(IND_LINE(L)), LINE_L0(IND_LINE(L))*1E-8, LINE_ION(IND_LINE(L)) &
                  ,EPLOW(IND_LINE(L)), NHYD, TEM(I,J,K), SIGMA(IND_LINE(L)), ALPHA(IND_LINE(L)), DAMP)
              NAN=CHECKNAN(DBLE(DAMP))
              IF (NAN.EQV..TRUE.) THEN
                 PRINT*,'NaN detected in DAMP'
                 PRINT*,DAMP
                 STOP
              ENDIF
              ! Determening the continuum absorption coefficient (1/cm): we need to send partial pressures, not number of atoms per cm^3
              CALL ION_CALC_RHO(1,TEM(I,J,K),NELEC,NHYD,NI,NII,NIII,DNI,DNII,DNIII)
              KC(I,J,K,L)=BACKGROUND_OPACITY(DBLE(TEM(I,J,K)),PELEC,NII*KBOL*TEM(I,J,K),NI*KBOL*TEM(I,J,K),NIII*KBOL*TEM(I,J,K) &
                   ,0D0,0D0,LINE_L0(IND_LINE(L)))
              ! We do not have H2 and H2+ in our chemical equilibrium so we use 0
              ! Get the absorption matrix
              CALL GET_ABS_MAT(I,J,K,L,DAMP)
              ! End loop spectral lines
           ENDDO
           ! End loop Z
        ENDDO
        ! End loop Y
     ENDDO
     ! End loop X
  ENDDO
  !
  ! Storing atmosphere
  ! 
  !OPEN(UNIT=1,FILE='/Users/borrero/Work/newcode/sintesis_test/atmos_nz50001_dz002_loggf.dat',STATUS='UNKNOWN')
  !DO I=1,NZ
  !   WRITE(UNIT=1,FMT='(F8.3,2X,E14.7,2X,E14.7,2X,E14.7,2X,E14.7,2X,E14.7,2X,E14.7,2X,E14.7,2X,E14.7)'),ZZ(I),TEM(1,1,I),PG(1,1,I),BX(1,1,I),BY(1,1,I),BZ(1,1,I) &
  !        ,VZ(1,1,I),RHO(1,1,I)
  !ENDDO
  !CLOSE(UNIT=1)
  !
  ! Solving the radiative transfer equation
  !
  OPEN(UNIT=1,FILE='/Users/borrero/Work/newcode/src_jun14/kk.dat',STATUS='UNKNOWN')
  DO I=1,NX
     DO J=1,NY
        DO K=1,NZ
           ! Loop in Z-axis is from bottom to top
           DO L=1,NUMW
              ! Boundary condition at the bottom
              IF (K.EQ.1) CALL SET_BOUNDARY(I,J,K,L)
              ! Next step
              IF (K.GT.1) CALL ANALYTICAL_SOLVER(I,J,K,L)
              ! End loop wavelengths
           ENDDO
           ! End loop Z
        ENDDO
        ! Write output
        K=NZ
        DO L=1,NUMW
           WRITE(UNIT=1, FMT='(F8.3,2X,E14.7,2X,E14.7,2X,E14.7,2X,E14.7)') WAVE(L),STOKES(I,J,K,L,1)/MAXVAL(STOKES(I,J,K,:,1)) &
                ,STOKES(I,J,K,L,2)/MAXVAL(STOKES(I,J,K,:,1)), STOKES(I,J,K,L,3)/MAXVAL(STOKES(I,J,K,:,1)) &
                ,STOKES(I,J,K,L,4)/MAXVAL(STOKES(I,J,K,:,1))
        ENDDO
        !WRITE(UNIT=1,FMT='(A,2X,I3)') '--------',K
        CLOSE(UNIT=1)
        STOP
        ! End loop Y
     ENDDO
     ! End loop X
  ENDDO
 
  








  STOP
  !
  10 FORMAT(A)
  20 FORMAT (E12.4,2X,E12.4,1X,E12.4,1X,F12.4)
  30 FORMAT (F8.3,2X,E14.7,2X,E14.7,2X,E14.7,2X,E14.7)

END PROGRAM TEST
