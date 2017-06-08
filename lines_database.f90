MODULE LINES_DATABASE
  USE CONS_PARAM
  USE DAMPING
  !
  IMPLICIT NONE
  INTEGER                                       :: NUML_DATABASE
  INTEGER,    ALLOCATABLE                       :: LINE_NUM(:), LINE_ZN(:), LINE_ION(:), OETRANSITION(:,:)
  REAL(DP),   ALLOCATABLE                       :: SL(:), LL(:), JL(:), SU(:), LU(:), JU(:), LINE_L0(:)
  REAL(DP),   ALLOCATABLE                       :: LOGGF(:), EPLOW(:), ALPHA(:), SIGMA(:)
  !
CONTAINS
  ! Subroutine read_lines_database
  ! Subroutine elecconf
  ! Subrourine optical_elec_tran
  !-----------------------------
  !-----------------------------
  SUBROUTINE READ_LINES_DATABASE
    !
    USE LOG
    IMPLICIT NONE
    INTEGER                                     :: I, NL, IERR, FLAG
    CHARACTER*2                                 :: SHELL
    CHARACTER*7                                 :: LOW, UPP
    CHARACTER*100                               :: FIN
    REAL(DP)                                    :: NLOW_EFF, NUPP_EFF
    !
    ! Does lines_database exist ?
    !
    OPEN(UNIT=1, FILE='lines_database.dat', STATUS='OLD', IOSTAT=IERR, FORM='FORMATTED')
    IF (IERR.NE.0) THEN
       PRINT*,'The file containing database about spectral lines (lines_database.dat)'
       PRINT*,'could not be fonud in the directory of the source code. STOP'
       CLOSE(UNIT=1)
       STOP
    ENDIF
    !
    ! Determine number of lines inside
    !
    NL=0
    IERR=0 ! reset ierr
    ALLOCATE(LINE_NUM(1), LINE_L0(1), LINE_ZN(1), LINE_ION(1))
    ALLOCATE(LOGGF(1), EPLOW(1), ALPHA(1), SIGMA(1))
    DO WHILE (IERR.EQ.0)
       READ(UNIT=1, FMT=10, IOSTAT=IERR) LINE_NUM, LINE_ZN, LINE_ION, LINE_L0, LOW, UPP, LOGGF &
            , EPLOW, ALPHA, SIGMA, SHELL
       IF (IERR.EQ.0) NL=NL+1
       IF (IERR.NE.0) EXIT
    ENDDO
    CLOSE(1)
    NUML_DATABASE=NL
    DEALLOCATE(LINE_NUM, LINE_L0, LINE_ZN, LINE_ION, LOGGF, EPLOW, ALPHA, SIGMA)
    !
    ! Read atomic data
    ! 
    ALLOCATE(LINE_NUM(NL), LINE_L0(NL), LINE_ZN(NL), LINE_ION(NL), LOGGF(NL), EPLOW(NL), ALPHA(NL), SIGMA(NL))
    ALLOCATE(SL(NL), SU(NL), LL(NL), LU(NL), JL(NL), JU(NL), OETRANSITION(NL,2))
    !
    ! Zero the elements of created arrays
    !
    LINE_NUM(:)=0
    LINE_L0(:)=0D0
    LINE_ZN(:)=0
    LINE_ION(:)=0
    LOGGF(:)=0D0
    EPLOW(:)=0D0
    ALPHA(:)=0D0
    SIGMA(:)=0D0
    SL(:)=0
    SU(:)=0
    LL(:)=0
    LU(:)=0
    JL(:)=0
    JU(:)=0
    OETRANSITION(:,:)=0
    !
    OPEN(UNIT=1, FILE='lines_database.dat', STATUS='OLD', IOSTAT=IERR, FORM='FORMATTED')
    DO I=1,NL
       !
       FLAG = 0
       IERR = 0
       !
       READ(UNIT=1, FMT=10, IOSTAT=IERR) LINE_NUM(I), LINE_ZN(I), LINE_ION(I), LINE_L0(I), LOW, UPP &
            , LOGGF(I), EPLOW(I), ALPHA(I), SIGMA(I), SHELL
       !
       ! Determine SLJ quantum numbers from electronic configuration
       !
       CALL ELECCONF(LOW,SL(I),LL(I),JL(I))
       CALL ELECCONF(UPP,SU(I),LU(I),JU(I))
       !
       ! Determine the transition type for the optical electron
       !
       CALL OPTICAL_ELEC_TRAN(SHELL,OETRANSITION(I,:),FLAG)
       !----------------------------------------
       ! We determine the collisional parameters
       !----------------------------------------
       ! Case 1: parameters provided in lines_database.dat
       IF (ALPHA(I) .NE.0 .AND. SIGMA(I) .NE. 0) THEN
          CALL LOGW_CHAR_REAL(44,'Collisional parameters provided by user for ',LINE_L0(I))
          CALL LOGW_CHAR_REAL(29,'Collisional cross section is ',SIGMA(I))
          CALL LOGW_CHAR_REAL(25,'Temperature parameter is ',ALPHA(I))
       ENDIF
       ! Case 2: parameters not provided by lines_database.dat
       ! but transition of the optical electron is available
       ! so we attempt an interpolation from Anstee, Barklem and O'Mara tables
       IF ( (ALPHA(I) .EQ. 0D0 .OR. SIGMA(I).EQ.0D0) .AND. (FLAG .EQ.0)) THEN
          ! Determine effective quantum n numbers
          CALL GET_NEFF(LINE_ZN(I),EPLOW(I),LINE_L0(I),LINE_ION(I),NLOW_EFF,NUPP_EFF)
          ! Interpolate from tables
          CALL GET_COLLISIONAL_PARAM(NLOW_EFF,NUPP_EFF,OETRANSITION(I,:),SIGMA(I),ALPHA(I),IERR)
          ! Check for errors:
          IF (IERR .NE. 0) CALL LOGW_CHAR_REAL(41,'Collisional parameters not available for ',LINE_L0(I))
          IF (IERR .EQ. 1) CALL LOGW_CHAR_CHAR(2,SHELL,42,'transition does not have tabulated tables.')
          IF (IERR .EQ. 2) CALL LOGW_CHAR_REAL(49,'Lower effective quantum number is outside tables ',NLOW_EFF)
          IF (IERR .EQ. 3) CALL LOGW_CHAR_REAL(49,'Upper effective quantum number is outside tables ',NUPP_EFF)
          IF (IERR .EQ. 0) THEN
             CALL LOGW_CHAR_REAL(40,'Collisional parameters interpolated for ',LINE_L0(I))
             CALL LOGW_CHAR_REAL(29,'Collisional cross section is ',SIGMA(I))
             CALL LOGW_CHAR_REAL(25,'Temperature parameter is ',ALPHA(I))
          ENDIF
       ENDIF
       ! Case 3: parameters not provided by lines_database.dat
       ! and transition of the optical electron is not available
       IF ( (ALPHA(I) .EQ. 0D0 .OR. SIGMA(I).EQ.0D0) .AND. (FLAG .NE.0)) THEN
          CALL LOGW_CHAR_REAL(41,'Collisional parameters not available for ',LINE_L0(I))
          CALL LOGW_CHAR(50,'For this line we will use Van der Waals broadening')
       ENDIF
    ENDDO
    READ(UNIT=1,FMT=20, IOSTAT=IERR) FIN
    IF (TRIM(FIN).NE.'END') THEN
       PRINT*,'File lines_database.dat could not be properly read. Check format. STOP'
       STOP
    ENDIF
    CLOSE(1)
    CALL LOGW_CHAR(88,'File containing database of spectral lines (lines_database.dat) has been correctly read.')
    !-------------------------------------
10  FORMAT(I4,2X,I3,2X,I1,2X,F9.4,2X,A7,2X,A7,2X,F6.3,2X,F5.3,2X,F5.3,2X,F5.1,2X,A2)
20  FORMAT(A)
  END SUBROUTINE READ_LINES_DATABASE
  !----------------------------
  SUBROUTINE ELECCONF(ELC,SOUT,LOUT,JOUT)
    !
    IMPLICIT NONE
    !
    CHARACTER*7                     :: ELC
    CHARACTER*1                     :: L
    INTEGER                         :: IERR
    REAL(8)                         :: SOUT, LOUT, JOUT
    REAL(4)                         :: S, J
    !
    READ(ELC(1:3),'(F3.1)') S
    READ(ELC(4:4),'(A1)') L
    READ(ELC(5:7),'(F3.1)') J
    JOUT=DBLE(J)
    SOUT=(DBLE(S)-1D0)/2D0
    IERR=1
    !
    !
    SELECT CASE (L)
    CASE('S')
       LOUT=0D0
       IERR=0
    CASE('P')
       LOUT=1D0
       IERR=0
    CASE('D')
       LOUT=2D0
       IERR=0
    CASE('F')
       LOUT=3D0
       IERR=0
    CASE('G')
       LOUT=4D0
       IERR=0
    CASE('H')
       LOUT=5D0
       IERR=0
    CASE('I')
       LOUT=6D0
       IERR=0
    END SELECT
    !
    IF (IERR.NE.0) THEN
       PRINT*,'Unknown electronic configuration '//TRIM(L)//' in lines_database.dat. STOP'
       STOP
    ENDIF
    !
  ENDSUBROUTINE ELECCONF
  !-----------------------------
  !-----------------------------
  SUBROUTINE OPTICAL_ELEC_TRAN(SHELL,TRANSITION,FLAG)
    ! This routine determines the transition type for the optical electron
    ! that is, the last electron in the electronic configuration. In general
    ! it is different from the transition that gives raise to the spectral line
    !
    IMPLICIT NONE
    !
    CHARACTER*2,               INTENT(IN)     :: SHELL
    INTEGER,                   INTENT(OUT)    :: TRANSITION(2)
    INTEGER,                   INTENT(INOUT)  :: FLAG
    CHARACTER*1                               :: ORBITAL
    INTEGER                                   :: I
    !
    DO I=1,2
       ORBITAL=SHELL(I:I)
       FLAG = 1
       SELECT CASE (ORBITAL)
          CASE('s')
             TRANSITION(I)=0
             FLAG=0
          CASE('S')
             TRANSITION(I)=0
             FLAG=0
          CASE('p')
             TRANSITION(I)=1
             FLAG=0
          CASE('P')
             TRANSITION(I)=1
             FLAG=0
          CASE('d')
             TRANSITION(I)=2
             FLAG=0
          CASE('D')
             TRANSITION(I)=2
             FLAG=0
          CASE('f')
             TRANSITION(I)=3
             FLAG=0
          CASE('F')
             TRANSITION(I)=3
             FLAG=0
       END SELECT
    ENDDO
    ! If flag = 1 we will not be able to determine the collisional damping with ABO theory
  END SUBROUTINE OPTICAL_ELEC_TRAN
  !-----------------------------------
  !-----------------------------------
END MODULE LINES_DATABASE
