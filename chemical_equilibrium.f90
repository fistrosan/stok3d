MODULE CHEMICAL_EQUILIBRIUM
  !
  ! J M Borrero
  ! April 2, 2013
  ! KIS, Freiburg
  !
  USE CONS_PARAM
  USE LINES_DATABASE
  USE ATOM_DATABASE
  USE DERIVVAR
CONTAINS
  ! get_pg
  ! get_rho
  ! get_nhnemu
  ! get_ne_rho
  ! get_ne_temp
  ! ion_calc_temprho
  ! ion_calc_netemp
  ! ion_calc_tempx
  ! get_dnhyddtemp
  ! get_klin
  !----------------------------------------------------------
  !----------------------------------------------------------
  SUBROUTINE GET_PG(TEMP,DENS,NHYD,NELEC,PATOM,PELEC,MOLECW)
    ! This routine obtains the total pressure (electrons+atoms)
    ! from the density and temperature
    IMPLICIT NONE
    !
    INTEGER                     :: I, J, ITER
    REAL(SP),   INTENT(IN)      :: TEMP, DENS
    REAL(DP),   INTENT(OUT)     :: NHYD, NELEC, PELEC, PATOM, MOLECW
    REAL(DP), DIMENSION(NELEM)  :: NT, NTI, NTII, NTIII
    REAL(DP)                    :: NELEC_NEW, NELEC_OLD, NELEC_EST, ERROR
    REAL(DP)                    :: NUME

    MOLECW=0D0
    PATOM=0D0
    PELEC=0D0
    !
    ! Estimation of electron density
    ! 
    IF (NELEC.EQ.0) NELEC = DENS/(150D0*MAMU)
    !
    ERROR=1D0
    ITER=1
    NELEC_NEW=NELEC
    NELEC_OLD=NELEC
    !
    ! Start loop until convergence
    !
    DO WHILE (ERROR.GT.1E-6.AND.ITER.LE.100)
       NELEC_OLD = NELEC_NEW
       ! Estimation of the hydrogen density
       NUME=DENS-NELEC_NEW*MELE
       NHYD=NUME/DENO_SAM
       NELEC_EST=0D0
       ! Initialize arrays
       NT(:)=0D0
       NTI(:)=0D0
       NTII(:)=0D0
       NTIII(:)=0D0
       ! Loop through all elements
       DO J=1,NELEM
          !
          NT(J)=NHYD*10D0**(ABUND(J)-12D0)
          ! Calculate neutrals, single ionized and double ionized atoms of specie J
          CALL ION_CALC_TEMPRHO(J,TEMP,NELEC_NEW,0D0,NT(J),NTI(J),NTII(J),NTIII(J))
          ! In the case of hydrogen (J=1) NTI corresponds to H- (electron capture)
          IF (J.GT.1) NELEC_EST=NELEC_EST+NTII(J)+NTIII(J)
          IF (J.EQ.1) NELEC_EST=NELEC_EST-NTI(J)+NTIII(J)
       ENDDO
       !
       NELEC_NEW=(NELEC_EST+NELEC_OLD)/2D0
       ERROR=ABS((NELEC_NEW-NELEC_OLD)/NELEC_NEW)
       !PRINT*,NELEC_OLD,NELEC_NEW,ERROR
       ITER = ITER + 1
    ENDDO
    !
    ! Calculate pressure due to electrons and due to atoms
    ! Calculate also mean molecuar weight
    !
    NELEC=NELEC_NEW
    PELEC=NELEC*KBOL*TEMP
    !
    DO I=1,NELEM
       PATOM=PATOM+10D0**(ABUND(I)-12D0)
    ENDDO
    PATOM=PATOM*NHYD*KBOL*TEMP
    MOLECW=((DENS*KBOL*TEMP)/MAMU)*(PELEC+PATOM)**(-1D0)
    !
  END SUBROUTINE GET_PG
  !----------------------------------------------------------
  !----------------------------------------------------------
  SUBROUTINE GET_RHO(TEMP,PG,NHYD,NELEC,PATOM,PELEC,DENS,MOLECW)
    ! This routine obtains the density and electron/atom pressure
    ! from the total pressure and temperature
    IMPLICIT NONE
    INTEGER                     :: I, J, ITER
    REAL(SP),   INTENT(IN)      :: TEMP, PG
    REAL(DP),   INTENT(OUT)     :: NHYD, NELEC, PELEC, PATOM, MOLECW, DENS
    REAL(DP), DIMENSION(NELEM)  :: NT, NTI, NTII, NTIII
    REAL(DP)                    :: NELEC_NEW, NELEC_OLD, NELEC_EST, ERROR
    REAL(DP)                    :: NUME
    !
    DENS=0D0
    MOLECW=0D0
    PATOM=0D0
    PELEC=0D0
    !
    ! Estimation of electron density
    ! 
    IF (NELEC.EQ.0) NELEC = PG/(101D0*KBOL*TEMP)
    !
    ERROR=1D0
    ITER=1
    NELEC_NEW=NELEC
    NELEC_OLD=NELEC
    !
    ! Start loop until convergence
    !
    DO WHILE (ERROR.GT.1E-6.AND.ITER.LE.100)
       NELEC_OLD = NELEC_NEW
       ! Estimation of the hydrogen density
       NUME=PG-NELEC_NEW*KBOL*TEMP
       NHYD=NUME/(DENO_SA*KBOL*TEMP)
       NELEC_EST=0D0
       ! Initialize arrays
       NT(:)=0D0
       NTI(:)=0D0
       NTII(:)=0D0
       NTIII(:)=0D0
       ! Loop through all elements
       DO J=1,NELEM
          !
          NT(J)=NHYD*10D0**(ABUND(J)-12D0)
          ! Calculate neutrals, single ionized and double ionized atoms of specie J
          CALL ION_CALC_TEMPRHO(J,TEMP,NELEC_NEW,0D0,NT(J),NTI(J),NTII(J),NTIII(J))
          ! In the case of hydrogen (J=1) NTI corresponds to H- (electron capture)
          IF (J.GT.1) NELEC_EST=NELEC_EST+NTII(J)+NTIII(J)
          IF (J.EQ.1) NELEC_EST=NELEC_EST-NTI(J)+NTIII(J)
       ENDDO
       !
       NELEC_NEW=(NELEC_EST+NELEC_OLD)/2D0
       ERROR=ABS((NELEC_NEW-NELEC_OLD)/NELEC_NEW)
       !PRINT*,NELEC_OLD,NELEC_NEW,ERROR
       ITER = ITER + 1
    ENDDO
    !
    ! Calculate pressure due to electrons and due to atoms
    !
    NELEC=NELEC_NEW
    PELEC=NELEC*KBOL*TEMP
    !
    PATOM=DENO_SA*NHYD*KBOL*TEMP
    !
    ! Calculate density
    !
    DENS=NHYD*DENO_SAM+NELEC*MELE
    MOLECW=((DENS*KBOL*TEMP)/MAMU)*(PELEC+PATOM)**(-1D0)
    !
  END SUBROUTINE GET_RHO
  !----------------------------------------------------------
  !----------------------------------------------------------
  SUBROUTINE GET_NE_RHO(TEMP,RHO,NE, DNE)
    ! Determines the electron pressure and its derivative with respect to T
    ! at constant density: DNELECDTEMP_RHO
    IMPLICIT NONE
    ! Input
    REAL(SP), INTENT(IN)                  :: TEMP, RHO
    ! Output
    REAL(DP), INTENT(OUT)                 :: NE
    REAL(DP), INTENT(OUT), OPTIONAL       :: DNE
    ! Internal
    INTEGER                               :: I, J, ITER
    REAL(DP)                              :: NH, ERROR, DERROR, NELEC, NELEC_NEW, NELEC_OLD
    REAL(DP)                              :: NELEC_EST, DNELEC_NEW, DNELEC_OLD, DNELEC_EST
    REAL(DP), DIMENSION(NELEM)            :: NT, NTI, NTII, NTIII, DNTI, DNTII, DNTIII
!apy 20170605 
    ! Determine NH
    NH=RHO/DENO_SAM
    ! Estimate number of electrons per cm^3
    NELEC=0.1*NH*DENO_SA
    !
    ERROR=1D0
    DERROR=1D0
    ITER=1
    NELEC_NEW=NELEC
    NELEC_OLD=NELEC
    DNELEC_NEW=0D0
    DNELEC_OLD=0D0
    !
    ! Start loop until convergence
    !
    DO WHILE (DERROR.GT.1E-6.AND.ITER.LE.100)
       NELEC_OLD = NELEC_NEW
       DNELEC_OLD = DNELEC_NEW
       NELEC_EST=0D0
       DNELEC_EST=0D0
       ! Initialize arrays
       NT(:)=0D0
       NTI(:)=0D0
       NTII(:)=0D0
       NTIII(:)=0D0
       DNTI(:)=0D0
       DNTII(:)=0D0
       DNTIII(:)=0D0
       ! Loop through all elements
       DO J=1,NELEM
          !
          NT(J)=NH*10D0**(ABUND(J)-12D0)
          ! Calculate neutrals, single ionized and double ionized atoms of specie J
          CALL ION_CALC_TEMPRHO(J,TEMP,NELEC_NEW,DNELEC_NEW,NT(J),NTI(J),NTII(J),NTIII(J),DNTI(J),DNTII(J),DNTIII(J))
          ! In the case of hydrogen (J=1) NTI corresponds to H- (electron capture)
          IF (J.GT.1) THEN
             NELEC_EST=NELEC_EST+NTII(J)+NTIII(J)
             DNELEC_EST=DNELEC_EST+DNTII(J)+DNTIII(J)
          ENDIF
          IF (J.EQ.1) THEN
             NELEC_EST=NELEC_EST-NTI(J)+NTIII(J)
             DNELEC_EST=DNELEC_EST-DNTI(J)+DNTIII(J)
          ENDIF
       ENDDO
       !
       NELEC_NEW=(NELEC_EST+NELEC_OLD)/2D0
       ERROR=ABS((NELEC_NEW-NELEC_OLD)/NELEC_NEW)
       !
       DNELEC_NEW=(DNELEC_EST+DNELEC_OLD)/2D0
       DERROR=ABS((DNELEC_NEW-DNELEC_OLD)/DNELEC_NEW)
       !
       !PRINT*,DNELEC_OLD,DNELEC_NEW,DERROR
       !PRINT*,NELEC_OLD,NELEC_NEW,ERROR
       ITER = ITER + 1
    ENDDO
   
    !
    NE=NELEC_NEW
    DNELECDTEMP_RHO=DNELEC_NEW
    IF (PRESENT(DNE)) DNE=DNELECDTEMP_RHO
    !
  END SUBROUTINE GET_NE_RHO
  !----------------------------------------------------------
  !----------------------------------------------------------
  SUBROUTINE GET_NE_TEMP(TEMP,PGAS,NE,DNE)
    ! Determines the number of electrons per cm^3 and its derivative with respect to the
    ! pressure at constant temperature: DNELECDPG_TEMP
    IMPLICIT NONE
    ! Input
    REAL(SP), INTENT(IN)                :: TEMP, PGAS
    ! Output
    REAL(DP), INTENT(OUT)               :: NE
    REAL(DP), INTENT(OUT), OPTIONAL     :: DNE
    ! Internal
    INTEGER                             :: I, J, K, ITER
    REAL(DP)                            :: NH, ERROR, NELEC, NELEC_NEW, NELEC_OLD
    REAL(DP)                            :: NELEC_EST, DNELEC_NEW, DNELEC_OLD, DNELEC_EST, NUM, DEN
    REAL(DP), DIMENSION(NELEM)          :: NT, NTI, NTII, NTIII, DNTI, DNTII, DNTIII
    REAL(DP), DIMENSION(NELEM)          :: FI, FII, FIII, DFI, DFII, DFIII, DNIDNE
    ! Estimate number of electrons per cm^3
    NELEC=0.1*PGAS/(KBOL*TEMP)
    ! Estimate number of hydrogen atoms per cm^3
    NH=(PGAS-NELEC*KBOL*TEMP)/(KBOL*TEMP*DENO_SA)
    !
    ERROR=1D0
    ITER=1
    NELEC_NEW=NELEC
    NELEC_OLD=NELEC
    !
    ! Start loop until convergence
    !
    DO WHILE (ERROR.GT.1E-6.AND.ITER.LE.100)
       NELEC_OLD = NELEC_NEW
       NELEC_EST=0D0
       ! Initialize arrays
       NT(:)=0D0
       NTI(:)=0D0
       NTII(:)=0D0
       NTIII(:)=0D0
       DNTI(:)=0D0
       DNTII(:)=0D0
       DNTIII(:)=0D0
       FI(:)=0D0
       FII(:)=0D0
       FIII(:)=0D0
       DFI(:)=0D0
       DFII(:)=0D0
       DFIII(:)=0D0
       ! Loop through all elements
       DO J=1,NELEM
          !
          NT(J)=NH*10D0**(ABUND(J)-12D0)
          ! Calculate neutrals, single ionized and double ionized atoms of specie J
          CALL ION_CALC_NETEMP(J,TEMP,NELEC_NEW,FI(J),FII(J),FIII(J),DFI(J),DFII(J),DFIII(J))
          NTI(J)=NT(J)*FI(J)
          NTII(J)=NT(J)*FII(J)
          NTIII(J)=NT(J)*FIII(J)
          ! In the case of hydrogen (J=1) NTI corresponds to H- (electron capture)
          IF (J.GT.1) THEN
             NELEC_EST=NELEC_EST+NTII(J)+NTIII(J)
          ENDIF
          IF (J.EQ.1) THEN
             NELEC_EST=NELEC_EST-NTI(J)+NTIII(J)
          ENDIF
       ENDDO
       !
       NELEC_NEW=(NELEC_EST+NELEC_OLD)/2D0
       ERROR=ABS((NELEC_NEW-NELEC_OLD)/NELEC_NEW)
       ! Update number of hydrogen atoms per cm^3
       NH=(PGAS-NELEC_NEW*KBOL*TEMP)/(KBOL*TEMP*DENO_SA)
       !
       !PRINT*,NELEC_OLD,NELEC_NEW,ERROR
       ITER = ITER + 1
    ENDDO
    !
    NE=NELEC_NEW
    !
    ! After convergence we calculate DNIDNE: derivative of number of atoms
    ! per cm^3 of specie "i" with respect to the number of  electrons
    ! per cm^3 at constant temperature
    DNIDNE(:)=0D0
    DO J=1,NELEM
       NUM=0D0
       DEN=0D0
       DO K=1,NELEM
         IF (K.EQ.1) THEN
            NUM=NUM+NT(K)*(-DFI(K)+DFIII(K))
            DEN=DEN+10D0**(ABUND(K)-ABUND(J))*(-FI(K)+FIII(K))
         ENDIF
         IF (K.GT.1) THEN
            NUM=NUM+NT(K)*(DFII(K)+DFIII(K))
            DEN=DEN+10D0**(ABUND(K)-ABUND(J))*(FII(K)+FIII(K))
         ENDIF
       ENDDO
       DNIDNE(J)=(1D0-NUM)/DEN
    ENDDO
    ! Now derivative of number of electrons per cm^3 with respect
    ! to gas pressure at constant temperature
    DNELECDPG_TEMP=1D0/((SUM(DNIDNE(:))+1D0)*KBOL*TEMP)
    IF (PRESENT(DNE)) DNE=DNELECDPG_TEMP
    !
  END SUBROUTINE GET_NE_TEMP
  !----------------------------------------------------------
  !----------------------------------------------------------
  SUBROUTINE ION_CALC_TEMPRHO(I,TEMP,NE,DNE,N,NI,NII,NIII,DNI,DNII,DNIII)
    !
    IMPLICIT NONE
    ! Input
    INTEGER,  INTENT(IN)            :: I
    REAL(SP), INTENT(IN)            :: TEMP
    REAL(DP), INTENT(IN)            :: NE, N, DNE
    ! Output
    REAL(DP), INTENT(OUT)           :: NI, NII, NIII
    REAL(DP), INTENT(OUT), OPTIONAL :: DNI, DNII, DNIII
    ! Internal
    REAL(SP)                        :: UI, UII, UIII
    REAL(SP)                        :: DUI, DUII, DUIII
    REAL(DP)                        :: LAMELEC, EXPOI, EXPOII, ALPHAI, ALPHAII
    REAL(DP)                        :: DLAMELEC, DEXPOI, DEXPOII, DALPHAI, DALPHAII
    !
    NI=0D0
    NII=0D0
    NIII=0D0
    ! Get partition functions
    CALL PARTITION_FUNCTION(I,TEMP,UI,UII,UIII,DUI,DUII,DUIII)
    ! ALPHAI and ALPHAII factors: these are the factors for the Saha equation
    LAMELEC = CSAHA1*DBLE(TEMP)**(-3D0/2D0)      ! This has units of cm^3
    EXPOI=CSAHA2*XI(I)/TEMP                      ! Dimensionless
    EXPOII=CSAHA2*XII(I)/TEMP                    ! Dimensionless
    ! Check over/under flows
    IF (EXPOI .GT. 0D0) EXPOI=DMIN1(EXPOI,650D0)
    IF (EXPOI .LT. 0D0) EXPOI=DMAX1(EXPOI,-650D0)
    IF (EXPOII .GT. 0D0) EXPOII=DMIN1(EXPOII,650D0)
    IF (EXPOII .LT. 0D0) EXPOII=DMAX1(EXPOII,-650D0)
    !
    EXPOI = EXP(EXPOI)
    EXPOII = EXP(EXPOII)
    !
    ALPHAI=(UI/(2D0*UII))*LAMELEC*EXPOI
    ALPHAII=(UII/(2D0*UIII))*LAMELEC*EXPOII
    !
    NI=N/(1D0+(1D0/(NE*ALPHAI))*(1D0+(1D0/(NE*ALPHAII))))
    NII=N/(NE*ALPHAI+1D0+(1D0/(NE*ALPHAII)))
    NIII=N/(1D0+NE*ALPHAII*(NE*ALPHAI+1D0))
    !
    ! Now derivatives with respect to the temperature at constant density
    !
    IF (PRESENT(DNI).AND.PRESENT(DNII).AND.PRESENT(DNIII)) THEN
       DNI=0D0
       DNII=0D0
       DNIII=0D0
       !
       DLAMELEC = -3D0*LAMELEC/(2D0*TEMP)
       DEXPOI=-EXPOI*CSAHA2*XI(I)/TEMP**2D0
       DEXPOII=-EXPOII*CSAHA2*XII(I)/TEMP**2D0
       !
       !DALPHAI=ALPHAI*(DUI/UI-DUII/UII+DLAMELEC/LAMELEC+DEXPOI/EXPOI)
       !DALPHAII=ALPHAII*(DUII/UII-DUIII/UIII+DLAMELEC/LAMELEC+DEXPOII/EXPOII)
       DALPHAI=ALPHAI*(DUI-DUII+DLAMELEC/LAMELEC+DEXPOI/EXPOI)
       DALPHAII=ALPHAII*(DUII-DUIII+DLAMELEC/LAMELEC+DEXPOII/EXPOII)
       !
       DNI=-N*(1D0+(1D0/(NE*ALPHAI))*(1D0+1D0/(NE*ALPHAII)))**(-2D0)*(1D0/(NE**2D0*ALPHAI))*((-1D0/ALPHAI)* &
            (DNE*ALPHAI+NE*DALPHAI)*(1D0+1D0/(NE*ALPHAII))-(1D0/(NE*ALPHAII**2D0))*(DNE*ALPHAII+NE*DALPHAII))
       !
       DNII=-N*(NE*ALPHAI+1D0+1D0/(NE*ALPHAII))**(-2D0)*(DNE*ALPHAI+NE*DALPHAI-(DNE*ALPHAII+NE*DALPHAII)/ &
            (NE**2D0*ALPHAII**2D0))
       !
       DNIII=-N*(1D0+NE*ALPHAII*(NE*ALPHAI+1D0))**(-2D0)*(2D0*NE*DNE*ALPHAI*ALPHAII+NE**2D0*DALPHAII*ALPHAI+ &
            NE**2D0*ALPHAII*DALPHAI+DNE*ALPHAII+NE*DALPHAII)
       !
    ENDIF
    !----
    !
  END SUBROUTINE ION_CALC_TEMPRHO
  !----------------------------------------------------------
  !----------------------------------------------------------
  SUBROUTINE ION_CALC_NETEMP(I,TEMP,NE,FI,FII,FIII,DFI,DFII,DFIII)
    !
    IMPLICIT NONE
    ! Input
    INTEGER,  INTENT(IN)            :: I
    REAL(SP), INTENT(IN)            :: TEMP
    REAL(DP), INTENT(IN)            :: NE
    ! Output
    REAL(DP), INTENT(OUT)           :: FI, FII, FIII
    REAL(DP), INTENT(OUT), OPTIONAL :: DFI, DFII, DFIII
    ! Internal
    REAL(SP)                        :: UI, UII, UIII
    REAL(SP)                        :: DUI, DUII, DUIII
    REAL(DP)                        :: LAMELEC, EXPOI, EXPOII, ALPHAI, ALPHAII
    !
    FI=0D0
    FII=0D0
    FIII=0D0
    ! Get partition functions
    CALL PARTITION_FUNCTION(I,TEMP,UI,UII,UIII,DUI,DUII,DUIII)
    ! ALPHAI and ALPHAII factors: these are the factors for the Saha equation
    LAMELEC = CSAHA1*DBLE(TEMP)**(-3D0/2D0)      ! This has units of cm^3
    EXPOI=CSAHA2*XI(I)/TEMP                      ! Dimensionless
    EXPOII=CSAHA2*XII(I)/TEMP                    ! Dimensionless
    ! Check over/under flows
    IF (EXPOI .GT. 0D0) EXPOI=DMIN1(EXPOI,650D0)
    IF (EXPOI .LT. 0D0) EXPOI=DMAX1(EXPOI,-650D0)
    IF (EXPOII .GT. 0D0) EXPOII=DMIN1(EXPOII,650D0)
    IF (EXPOII .LT. 0D0) EXPOII=DMAX1(EXPOII,-650D0)
    !
    EXPOI = EXP(EXPOI)
    EXPOII = EXP(EXPOII)
    !
    ALPHAI=(UI/(2D0*UII))*LAMELEC*EXPOI
    ALPHAII=(UII/(2D0*UIII))*LAMELEC*EXPOII
    !
    FI=1D0/(1D0+(1D0/(NE*ALPHAI))*(1D0+(1D0/(NE*ALPHAII))))
    FII=1D0/(NE*ALPHAI+1D0+(1D0/(NE*ALPHAII)))
    FIII=1D0/(1D0+NE*ALPHAII*(NE*ALPHAI+1D0))
    !
    ! Now derivatives of FI, FII and FIII with respect to the electron density at constant temperature
    !
    IF (PRESENT(DFI).AND.PRESENT(DFII).AND.PRESENT(DFIII)) THEN
       DFI=0D0
       DFII=0D0
       DFIII=0D0
       !
       DFI=2D0*NE*ALPHAI*ALPHAII/(NE**2D0*ALPHAI*ALPHAII+NE*ALPHAII+1D0)**2D0
       !
       DFII=ALPHAII*(1D0-NE**2D0*ALPHAI*ALPHAII)/(NE**2D0*ALPHAI*ALPHAII+NE*ALPHAII+1D0)**2D0
       !
       DFIII=-ALPHAII*(2D0*NE*ALPHAI+1D0)/(NE**2D0*ALPHAI*ALPHAII+NE*ALPHAII+1D0)**2D0
       !
    ENDIF
    !----
    !
  END SUBROUTINE ION_CALC_NETEMP
  !----------------------------------------------------------
  !----------------------------------------------------------
  SUBROUTINE ION_CALC_TEMPX(INDEX,TEMP,NE,N,NI,DNIDT_PG,DNIDT_RHO)
    IMPLICIT NONE
    INTEGER,   INTENT(IN)             :: INDEX
    REAL(SP),  INTENT(IN)             :: TEMP
    REAL(DP),  INTENT(IN)             :: NE, N
    REAL(DP),  INTENT(OUT)            :: NI, DNIDT_PG, DNIDT_RHO
    ! Internal
    REAL(SP)                          :: UI, UII, UIII
    REAL(SP)                          :: DUI, DUII, DUIII
    REAL(DP)                          :: LAMELEC, EXPOI, EXPOII, ALPHAI, ALPHAII
    REAL(DP)                          :: DLAMELEC, DEXPOI, DEXPOII, DALPHAI, DALPHAII
    REAL(DP)                          :: DDENODTEMP_PG, DDENODTEMP_RHO
    REAL(DP)                          :: DNUMEDTEMP_PG, DNUMEDTEMP_RHO
    REAL(DP)                          :: DENO, NUME
    INTEGER                           :: ZN, ION
    !
    ZN=LINE_ZN(INDEX)
    ION=LINE_ION(INDEX)
    NI=0D0
    DNIDT_PG=0D0
    DNIDT_RHO=0D0
    !
    ! Get partition functions
    CALL PARTITION_FUNCTION(ZN,TEMP,UI,UII,UIII,DUI,DUII,DUIII)
    ! ALPHAI and ALPHAII factors: these are the factors for the Saha equation
    LAMELEC = CSAHA1*DBLE(TEMP)**(-3D0/2D0)      ! This has units of cm^3
    EXPOI=CSAHA2*XI(ZN)/TEMP                      ! Dimensionless
    EXPOII=CSAHA2*XII(ZN)/TEMP                    ! Dimensionless
    ! Check over/under flows
    IF (EXPOI .GT. 0D0) EXPOI=DMIN1(EXPOI,650D0)
    IF (EXPOI .LT. 0D0) EXPOI=DMAX1(EXPOI,-650D0)
    IF (EXPOII .GT. 0D0) EXPOII=DMIN1(EXPOII,650D0)
    IF (EXPOII .LT. 0D0) EXPOII=DMAX1(EXPOII,-650D0)
    !
    EXPOI = EXP(EXPOI)
    EXPOII = EXP(EXPOII)
    !
    ALPHAI=(UI/(2D0*UII))*LAMELEC*EXPOI
    ALPHAII=(UII/(2D0*UIII))*LAMELEC*EXPOII
    !
    DLAMELEC = -3D0*LAMELEC/(2D0*TEMP)
    DEXPOI=-EXPOI*CSAHA2*XI(ZN)/TEMP**2D0
    DEXPOII=-EXPOII*CSAHA2*XII(ZN)/TEMP**2D0
    !DALPHAI=(ALPHAI/(UI*UII))*(DUI*UII-DUII*UI)+(ALPHAI/LAMELEC)*DLAMELEC+(ALPHAI/EXPOI)*DEXPOI
    !DALPHAII=(ALPHAII/(UII*UIII))*(DUII*UIII-DUIII*UII)+(ALPHAII/LAMELEC)*DLAMELEC+(ALPHAII/EXPOII)*DEXPOII
    DALPHAI=(ALPHAI)*(DUI-DUII)+(ALPHAI/LAMELEC)*DLAMELEC+(ALPHAI/EXPOI)*DEXPOI
    DALPHAII=(ALPHAII)*(DUII-DUIII)+(ALPHAII/LAMELEC)*DLAMELEC+(ALPHAII/EXPOII)*DEXPOII
    !
    DENO=NE**2D0*ALPHAI*ALPHAII+NE*ALPHAII+1D0
    DDENODTEMP_PG=2D0*NE*DNELECDTEMP_PG*ALPHAI*ALPHAII+NE**2D0*DALPHAI*ALPHAII+NE**2D0*ALPHAI*DALPHAII+ &
         DNELECDTEMP_PG*ALPHAII+NE*DALPHAII
    DDENODTEMP_RHO=2D0*NE*DNELECDTEMP_RHO*ALPHAI*ALPHAII+NE**2D0*DALPHAI*ALPHAII+NE**2D0*ALPHAI*DALPHAII+ &
         DNELECDTEMP_RHO*ALPHAII+NE*DALPHAII
    !
    SELECT CASE(ION)
       CASE(1)
          NUME=NE**2D0*ALPHAI*ALPHAII
          DNUMEDTEMP_PG=2D0*NE*DNELECDTEMP_PG*ALPHAI*ALPHAII+NE**2D0*DALPHAI*ALPHAII+ &
               NE**2D0*ALPHAI*DALPHAII
          DNUMEDTEMP_RHO=2D0*NE*DNELECDTEMP_RHO*ALPHAI*ALPHAII+NE**2D0*DALPHAI*ALPHAII+ &
               NE**2D0*ALPHAI*DALPHAII
       CASE(2)
          NUME=NE*ALPHAII
          DNUMEDTEMP_PG=DNELECDTEMP_PG*ALPHAII+NE*DALPHAII
          DNUMEDTEMP_RHO=DNELECDTEMP_RHO*ALPHAII+NE*DALPHAII
       CASE(3)
          NUME=1D0
          DNUMEDTEMP_PG=0D0
          DNUMEDTEMP_RHO=0D0
    END SELECT
    !
    NI=N*NUME/DENO
    DNIDT_PG=(NUME/DENO)*10D0**(ABUND(ZN)-12D0)*DNHYDDTEMP_PG+N*(DNUMEDTEMP_PG*DENO-DDENODTEMP_PG*NUME)/DENO**2D0
    DNIDT_RHO=(NUME/DENO)*10D0**(ABUND(ZN)-12D0)*DNHYDDTEMP_RHO+N*(DNUMEDTEMP_RHO*DENO-DDENODTEMP_RHO*NUME)/DENO**2D0
    !
  END SUBROUTINE ION_CALC_TEMPX
  !----------------------------------------------------------
  !----------------------------------------------------------
  PURE SUBROUTINE GET_NHNEMU(T,P,D,NH,NE,MU)
    !  This routine assumes temperature, gas pressure and density are known
    !  and calculates number of hydrogen atoms per cm^3, number of electrons per cm^3
    !  and molecular weight.
    IMPLICIT NONE
    REAL(SP), INTENT(IN)     :: T, P, D
    REAL(DP), INTENT(OUT)    :: NH, NE
    REAL(SP), INTENT(OUT)    :: MU
    ! Internal
    INTEGER                  :: I
    !
    NH=(DBLE(P)-DBLE(D)*KBOL*DBLE(T)/MELE)/(KBOL*DBLE(T)*(DENO_SA-DENO_SAM/MELE))
    NE=(DBLE(D)-NH*DENO_SAM)/MELE
    MU=D*KBOL*T/(REAL(MHYD)*P)
    !
  END SUBROUTINE GET_NHNEMU
  !----------------------------------------------------------
  !----------------------------------------------------------
  SUBROUTINE GET_DNHYDDTEMP(T,P,D,M,DNHDT_P,DNHDT_R)
             
    IMPLICIT NONE
    REAL(SP), INTENT(IN)              :: T, P, D, M
    REAL(DP), INTENT(OUT), OPTIONAL   :: DNHDT_P,DNHDT_R
    ! Internal
    INTEGER                           :: I
    REAL(DP)                          :: NELEC
    ! This gives the derivative of ne with respect to PG at constant temperature
    CALL GET_NE_TEMP(T, P, NELEC)
    ! Now the derivative of ne with respect to T at constant density
    CALL GET_NE_RHO(T, D, NELEC)
    ! Derivative of ne with repect to T at constant gas pressure
    DNELECDTEMP_PG = DNELECDTEMP_RHO - ((DBLE(D)*KBOL)/(DBLE(M)*MHYD))*DNELECDPG_TEMP
    ! Derivative of NHYD with repect to T at constant gas pressure
    DNHYDDTEMP_PG=(-1D0/DENO_SA)*(DBLE(P)/(KBOL*DBLE(T)**2D0)-DNELECDTEMP_PG)
    ! Derivative of NHYD with repect to T at constant density
    DNHYDDTEMP_RHO= (-MELE/DENO_SAM)*DNELECDTEMP_RHO
    !
    IF (PRESENT(DNHDT_P).AND.PRESENT(DNHDT_R)) THEN
       DNHDT_P=DNHYDDTEMP_PG
       DNHDT_R=DNHYDDTEMP_RHO
    ENDIF
    !
  END SUBROUTINE GET_DNHYDDTEMP
  !----------------------------------------------------------
  !----------------------------------------------------------
  SUBROUTINE GET_KLIN(INDEX,TEMP,NELEC,N,KLIN,DKDT_PG,DKDT_RHO)
    !
    IMPLICIT NONE
    ! Input
    INTEGER, INTENT(IN)                :: INDEX
    REAL(SP), INTENT(IN)               :: TEMP
    REAL(DP), INTENT(IN)               :: NELEC, N
    ! Output
    REAL(SP), INTENT(OUT)              :: KLIN
    REAL(DP), INTENT(OUT),OPTIONAL     :: DKDT_PG, DKDT_RHO
    ! Internal
    REAL(DP)                           :: NIR, DNIRDT_PG, DNIRDT_RHO
    !
    CALL ION_CALC_TEMPX(INDEX,TEMP,NELEC,N,NIR,DNIRDT_PG,DNIRDT_RHO)
    !PRINT*,NIR,DNIRDT_PG,DNIRDT_RHO
    !STOP
    

  END SUBROUTINE GET_KLIN

  !----------------------------------------------------------
  !----------------------------------------------------------
END MODULE CHEMICAL_EQUILIBRIUM
