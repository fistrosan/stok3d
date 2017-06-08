MODULE SIMUL_BOX
  !
  ! J M Borrero
  ! March 22, 2013
  !
  USE LOG
  USE CONS_PARAM
  USE STRINGS
  CHARACTER*100                               :: INPUTFILE
  INTEGER                                     :: NUML, NUMW, NX, NY, NZ
  REAL(DP)                                    :: DX, DY, DZ
  REAL(DP),         ALLOCATABLE               :: DELTA_LAMBDA(:), WAVE_INI(:)
  INTEGER,          ALLOCATABLE               :: IND_LINE(:), NUMWAVE(:), PIXEL_INI(:), PIXEL_END(:)
  !-------------------------------------------------------------------------------------
  REAL(DP),         ALLOCATABLE               :: WAVE(:)
  !-------------------------------------------------------------------------------------
  REAL(SP),         ALLOCATABLE               :: TEM(:,:,:), BX(:,:,:), BY(:,:,:), BZ(:,:,:)
  REAL(SP),         ALLOCATABLE               :: RHO(:,:,:), VX(:,:,:), VY(:,:,:), VZ(:,:,:)
  REAL(SP),         ALLOCATABLE               :: PG(:,:,:), MW(:,:,:), PEL(:,:,:), XX(:), YY(:), ZZ(:)
  REAL(SP),         ALLOCATABLE               :: KLIN(:,:,:,:), KC(:,:,:,:), STOKES(:,:,:,:,:)
  REAL(SP),         ALLOCATABLE               :: DSTOKES(:,:,:,:,:,:)
  REAL(SP),         DIMENSION(4,4)            :: IMAT
  !-------------------------------------------------------------------------------------
  !
CONTAINS
  ! read_boxfile
  ! read_boxdata
  ! read_lines
  ! wave_init
  ! box_init
  ! zero_iter_hydrostatic
  !------------------------
  !------------------------
  SUBROUTINE READ_BOXFILE
    !
    USE LINES_DATABASE
    IMPLICIT NONE
    INTEGER, PARAMETER                         :: REGIONS = 2 !BOX:,LINES:
    CHARACTER*100                              :: LINE
    INTEGER                                    :: IERR, I, J, K, LINES_READ
    INTEGER                                    :: FOUND, FX, IND
    INTEGER                                    :: NUM_LINES_OBS
    INTEGER,    DIMENSION(NUML_DATABASE)       :: DIF
    !
    ! Does input file exist ?
    ! 
    LINES_READ = 0
    CALL GETARG(1,INPUTFILE)
    OPEN(UNIT=1, FILE=TRIM(INPUTFILE), STATUS='OLD', IOSTAT=IERR, FORM='FORMATTED')
    IF (IERR.NE.0) THEN
       PRINT*,'The file containing information about observed spectral lines '//TRIM(INPUTFILE)
       PRINT*,'could not be fonud in the directory of the source code. STOP'
       CLOSE(UNIT=1)
       STOP
    ENDIF
    CLOSE(1)
    FOUND = 0
    !
    DO WHILE (FOUND.LT.REGIONS)
       LINES_READ = LINES_READ+1
       OPEN(UNIT=1, FILE=TRIM(INPUTFILE), STATUS='OLD', IOSTAT=IERR, FORM='FORMATTED')
       ! Read until last read line + 1
       DO I=1,LINES_READ
          READ(UNIT=1, FMT='(A)',IOSTAT=IERR) LINE
       ENDDO
       !
       SELECT CASE (TRIM(LINE))
       CASE('BOX:')
          CLOSE(1)
          CALL READ_BOXDATA(LINES_READ)
          FOUND = FOUND +1
       CASE('LINES:')
          CLOSE(1)
          CALL READ_LINES(LINES_READ,NUM_LINES_OBS)
          FOUND = FOUND +1
       CLOSE(1)
       ENDSELECT
    ENDDO
    ! Check that indexes of spectral lines exist in database
    DO I=1,NUML
       DIF(:)=1E3
       DO J=1,NUML_DATABASE
          DIF(J)=ABS(LINE_NUM(J)-IND_LINE(I))
       ENDDO
       IF (MINVAL(DIF).NE.0) THEN
          PRINT*,'Spectral line index in '//TRIM(INPUTFILE)//' cannot'
          PRINT*,'be found in databse (lines_database.dat). STOP'
          STOP
       ENDIF
    ENDDO
    ! Write in log file
    !CALL LOGW_CHAR(65+LEN_TRIM(TRIM(INPUTFILE)),'File containing observed spectral lines '//TRIM(INPUTFILE)//' has been read correctly.')
    DO I=1,NUML
       CALL LOGW_CHAR_REAL(20,'Spectral line found:',LINE_L0(IND_LINE(I)))
    ENDDO
    !
  END SUBROUTINE READ_BOXFILE
  !--------------------------
  !-------------------------- 
  SUBROUTINE READ_BOXDATA(LINES_READ)
    !
    IMPLICIT NONE
    INTEGER,           PARAMETER     :: NARGS = 2
    INTEGER                          :: LINES_READ, I, IERR, IOS
    INTEGER,       DIMENSION(NARGS)  :: IND_INI,IND_END
    CHARACTER*100                    :: LINE
    ! Open
    OPEN(UNIT=1, FILE=TRIM(INPUTFILE), STATUS='OLD', IOSTAT=IERR, FORM='FORMATTED')
    ! Read until last line
    DO I=1,LINES_READ
       READ(UNIT=1, FMT='(A)',IOSTAT=IERR) LINE
    ENDDO
    ! Read new line: NX
    READ(UNIT=1, FMT='(A)',IOSTAT=IERR) LINE
    CALL SPLIT_LINE(LINE,NARGS,IND_INI,IND_END)
    CALL VALUE_SI(LINE(IND_INI(1):IND_END(1)),NX,IOS)
    CALL VALUE_DR(LINE(IND_INI(2):IND_END(2)),DX,IOS)
    LINES_READ = LINES_READ + 1
    IF (IERR.NE.0) THEN
       PRINT*,'Error while reading '//TRIM(INPUTFILE)//'. STOP.'
       PRINT*,'Error in line #',LINES_READ
       STOP
    ENDIF
    ! Read new line: NY
    READ(UNIT=1, FMT='(A)',IOSTAT=IERR) LINE
    CALL SPLIT_LINE(LINE,NARGS,IND_INI,IND_END)
    CALL VALUE_SI(LINE(IND_INI(1):IND_END(1)),NY,IOS)
    CALL VALUE_DR(LINE(IND_INI(2):IND_END(2)),DY,IOS)
    LINES_READ = LINES_READ + 1
    IF (IERR.NE.0) THEN
       PRINT*,'Error while reading '//TRIM(INPUTFILE)//'. STOP.'
       PRINT*,'Error in line #',LINES_READ
       STOP
    ENDIF
    ! Read new line: NZ
    READ(UNIT=1, FMT='(A)',IOSTAT=IERR) LINE
    CALL SPLIT_LINE(LINE,NARGS,IND_INI,IND_END)
    CALL VALUE_SI(LINE(IND_INI(1):IND_END(1)),NZ,IOS)
    CALL VALUE_DR(LINE(IND_INI(2):IND_END(2)),DZ,IOS)
    LINES_READ = LINES_READ + 1
    IF (IERR.NE.0) THEN
       PRINT*,'Error while reading '//TRIM(INPUTFILE)//'. STOP.'
       PRINT*,'Error in line #',LINES_READ
       STOP
    ENDIF
    !
    CLOSE(1)
    !
  END SUBROUTINE READ_BOXDATA
  !------------------------
  !------------------------
  SUBROUTINE READ_LINES(LINES_READ, NUM_LINES_OBS)
    !
    IMPLICIT NONE
    INTEGER,   PARAMETER                  :: NARGS = 4
    INTEGER                               :: IERR, I, J, IOS, LINES_READ, NUM_LINES_OBS
    CHARACTER*100                         :: LINE
    INTEGER,      DIMENSION(NARGS)        :: IND_INI, IND_END
    ! Open
    OPEN(UNIT=1, FILE=TRIM(INPUTFILE), STATUS='OLD', IOSTAT=IERR, FORM='FORMATTED')
    ! Read until last line
    DO I=1,LINES_READ
       READ(UNIT=1, FMT='(A)',IOSTAT=IERR) LINE
    ENDDO
    !
    ! Determine number of lines inside and allocate arrays
    !
    NUM_LINES_OBS = 0
    DO WHILE (IERR.EQ.0)
       READ(UNIT=1, FMT='(A)',IOSTAT=IERR) LINE
       IF (IERR.NE.0) EXIT
       IF (TRIM(LINE).EQ.'BOX:') EXIT
       !IF (TRIM(LINE).EQ.'OTHER:') EXIT
       NUM_LINES_OBS=NUM_LINES_OBS+1
    ENDDO
    CLOSE(1)
    !
    NUML=NUM_LINES_OBS
    ALLOCATE(IND_LINE(NUML), WAVE_INI(NUML), DELTA_LAMBDA(NUML), NUMWAVE(NUML))
    !
    ! Read data for observed spectral lines
    !
    ! Read until last line
    OPEN(UNIT=1, FILE=TRIM(INPUTFILE), STATUS='OLD', IOSTAT=IERR, FORM='FORMATTED')
    DO I=1,LINES_READ
       READ(UNIT=1, FMT='(A)',IOSTAT=IERR) LINE
    ENDDO
    ! Read info about spectral lines
    DO I=1,NUML
       READ(UNIT=1, FMT='(A)',IOSTAT=IERR) LINE
       CALL SPLIT_LINE(LINE,NARGS,IND_INI,IND_END)
       DO J=1,NARGS
          SELECT CASE (J)
          CASE(1)
             CALL VALUE_SI(LINE(IND_INI(J):IND_END(J)),IND_LINE(I),IOS)
          CASE(2)
             CALL VALUE_SI(LINE(IND_INI(J):IND_END(J)),NUMWAVE(I),IOS)
          CASE(3)
             CALL VALUE_DR(LINE(IND_INI(J):IND_END(J)),WAVE_INI(I),IOS)
          CASE(4)
             CALL VALUE_DR(LINE(IND_INI(J):IND_END(J)),DELTA_LAMBDA(I),IOS)
          END SELECT
       ENDDO
    ENDDO
    CLOSE(1)
    LINES_READ=LINES_READ+NUM_LINES_OBS
  END SUBROUTINE READ_LINES
  !-------------------------------------------
  !-------------------------------------------
  SUBROUTINE WAVE_INIT
    USE LINES_DATABASE
    IMPLICIT NONE
    INTEGER                                    :: I, J
    ! Total number of wavelenghs
    NUMW=SUM(NUMWAVE)
    ! Allocate localizers
    ALLOCATE( PIXEL_INI(NUML), PIXEL_END(NUML) )
    PIXEL_INI(1)=1
    PIXEL_END(1)=NUMWAVE(1)
    DO I=2,NUML
       PIXEL_INI(I)=PIXEL_END(I-1)+1
       PIXEL_END(I)=PIXEL_END(I-1)+NUMWAVE(I)
    ENDDO
    ! Allocate and create wavelength array
    ALLOCATE(WAVE(NUMW))
    DO I=1,NUML
       DO J=PIXEL_INI(I),PIXEL_END(I)
          WAVE(J)=WAVE_INI(I)+DELTA_LAMBDA(I)*(J-PIXEL_INI(I))
       ENDDO
    ENDDO
    ! Message in log file
    CALL LOGW_CHAR(25,'Wavelength array created.')
    !------------------------------------------------
  END SUBROUTINE WAVE_INIT
  !--------------------------------------------------
  !--------------------------------------------------
  SUBROUTINE BOX_INIT
    IMPLICIT NONE
    INTEGER                   :: IERR, I
    ! Geometrical arrays
    ALLOCATE (XX(NX), YY(NY), ZZ(NZ))
    DO I=1,NX
       XX(I)=REAL(I-1)*DX
    ENDDO
    DO I=1,NY
       YY(I)=REAL(I-1)*DY
    ENDDO
    DO I=1,NZ
       ZZ(I)=REAL(I-1)*DZ
    ENDDO
    ! Temperature
    ALLOCATE(TEM(NX,NY,NZ), STAT=IERR)
    IF (IERR.NE.0) THEN
       PRINT*,'Not enought memory to allocate TEM array. STOP.'
       STOP
    ENDIF
    ! Density
    ALLOCATE(RHO(NX,NY,NZ), STAT=IERR)
    IF (IERR.NE.0) THEN
       PRINT*,'Not enought memory to allocate RHO array. STOP.'
       STOP
    ENDIF
    ! Gas pressure
    ALLOCATE(PG(NX,NY,NZ), STAT=IERR)
    IF (IERR.NE.0) THEN
       PRINT*,'Not enought memory to allocate PG array. STOP.'
       STOP
    ENDIF
    ! Electron pressure
    ALLOCATE(PEL(NX,NY,NZ), STAT=IERR)
    IF (IERR.NE.0) THEN
       PRINT*,'Not enought memory to allocate PEL array. STOP.'
       STOP
    ENDIF
    ! Molecular weight
    ALLOCATE(MW(NX,NY,NZ), STAT=IERR)
    IF (IERR.NE.0) THEN
       PRINT*,'Not enought memory to allocate MW array. STOP.'
       STOP
    ENDIF
    ! Opactiy
    ALLOCATE(KC(NX,NY,NZ,NUML), STAT=IERR)
    IF (IERR.NE.0) THEN
       PRINT*,'Not enought memory to allocate KC array. STOP.'
       STOP
    ENDIF
    ! BX
    ALLOCATE(BX(NX,NY,NZ), STAT=IERR)
    IF (IERR.NE.0) THEN
       PRINT*,'Not enought memory to allocate BX array. STOP.'
       STOP
    ENDIF
    ! BY
    ALLOCATE(BY(NX,NY,NZ), STAT=IERR)
    IF (IERR.NE.0) THEN
       PRINT*,'Not enought memory to allocate BY array. STOP.'
       STOP
    ENDIF
    ! BZ
    ALLOCATE(BZ(NX,NY,NZ), STAT=IERR)
    IF (IERR.NE.0) THEN
       PRINT*,'Not enought memory to allocate BZ array. STOP.'
       STOP
    ENDIF
    ! VX
    ALLOCATE(VX(NX,NY,NZ), STAT=IERR)
    IF (IERR.NE.0) THEN
       PRINT*,'Not enought memory to allocate VX array. STOP.'
       STOP
    ENDIF
    ! VY
    ALLOCATE(VY(NX,NY,NZ), STAT=IERR)
    IF (IERR.NE.0) THEN
       PRINT*,'Not enought memory to allocate VY array. STOP.'
       STOP
    ENDIF
    ! VZ
    ALLOCATE(VZ(NX,NY,NZ), STAT=IERR)
    IF (IERR.NE.0) THEN
       PRINT*,'Not enought memory to allocate VZ array. STOP.'
       STOP
    ENDIF
    ! KLIN
    ALLOCATE(KLIN(NX,NY,NZ,NUML),STAT=IERR)
    IF (IERR.NE.0) THEN
       PRINT*,'Not enought memory to allocate KLIN1 array. STOP.'
       STOP
    ENDIF
    ! Zero some arrays
    KLIN(:,:,:,:)=0.
    KC(:,:,:,:)=0.
    VX(:,:,:)=0.
    VY(:,:,:)=0.
    VZ(:,:,:)=0.
    ! STOKES
    ALLOCATE(STOKES(NX,NY,NZ,NUMW,4),STAT=IERR)
     IF (IERR.NE.0) THEN
       PRINT*,'Not enought memory to allocate STOKES array. STOP.'
       STOP
    ENDIF
    ! DSTOKES
    ALLOCATE(DSTOKES(NX,NY,NZ,NUMW,4,5),STAT=IERR)
     IF (IERR.NE.0) THEN
       PRINT*,'Not enought memory to allocate DSTOKES array. STOP.'
       STOP
    ENDIF
    ! IMAT
    IMAT(:,:)=0.
    DO I=1,4
       IMAT(I,I)=1.
    ENDDO
    !
  END SUBROUTINE BOX_INIT
  !------------------------------------------------
  !------------------------------------------------
  SUBROUTINE ZERO_ITER_HYDROSTATIC
    !
    USE CHEMICAL_EQUILIBRIUM
    IMPLICIT NONE
    !
    INTEGER          :: I, J, K
    REAL(DP)         :: NHYD, NELEC, PATOM, PELEC, MOLECW, DENS
    ! Filling temperature
    DO I=1,NZ
       TEM(:,:,I)=(ZZ(I)/(REAL(NZ-1)*DZ))*(3500.-7500.)+7500.
    ENDDO
    ! Filling magnetic field
    BX(:,:,:)=150.
    BY(:,:,:)=150.
    BZ(:,:,:)=150.
    ! Filling vertical component of the velocity
    VZ(:,:,:)=5E4
    ! Filling upper boundary condition for gas pressure
    PG(:,:,NZ)=1.E2
    ! Filling density, electron pressure and mean molecular weight
    DO I=1,NX
       DO J=1,NY
          DO K=NZ,1,-1
             ! Calculate number of hydrogen and electrons per cm^3 using chemical equilibrium
             CALL GET_RHO(TEM(I,J,K),PG(I,J,K),NHYD,NELEC,PATOM,PELEC,DENS,MOLECW)
             MW(I,J,K)=MOLECW
             RHO(I,J,K)=DENS
             PEL(I,J,K)=PELEC
             IF (K.GT.1) THEN
                ! Apply hydrostatic equilibirum to determine gas pressure on grid point below current one
                PG(I,J,K-1)=PG(I,J,K)*EXP(-MHYD*GRAV*MW(I,J,K)*(1./TEM(I,J,K-1)+1./TEM(I,J,K))*1E5*(ZZ(K-1)-ZZ(K))/(2.*KBOL))
                ! Here we consider the mean molecular weight to be constant between K and K-1
                ! It will be re-evaluated at the start of the next K-loop
                !PRINT*,TEM(I,J,K),PG(I,J,K),RHO(I,J,K),PEL(I,J,K)
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    !
  END SUBROUTINE ZERO_ITER_HYDROSTATIC

  !
END MODULE SIMUL_BOX
