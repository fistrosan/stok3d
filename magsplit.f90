MODULE MAGSPLIT
  !
  ! J M Borrero
  ! Jan 18, 2012
  !
  USE CONS_PARAM
  USE LINES_DATABASE
  USE SIMUL_BOX
  !
  REAL(SP),       ALLOCATABLE     :: LAMBDA_SC(:)         ! Lambda Split Component
  REAL(SP),       ALLOCATABLE     :: STRENGTH_SC(:)       ! Strength of the Split Component
  INTEGER,        ALLOCATABLE     :: TYPE_SC(:)           ! Type of the Split Component
  INTEGER                         :: NPI, NSR, NSB, NTOT
  !
CONTAINS
  ! zeeman_comp
  !
  !--------------------------------------------
  !--------------------------------------------
  SUBROUTINE ZEEMAN_COMP(I)
    !
    IMPLICIT NONE
    INTEGER                                           :: I
    !
    REAL(SP)                                          :: GLOW, GUPP, GEFF
    REAL(SP)                                          :: SLOW, SUPP, LLOW, LUPP, JLOW, JUPP
    REAL(SP)                                          :: DELTALB
    REAL(SP)                                          :: SUM_PI, SUM_SB, SUM_SR
    INTEGER                                           :: J, K, COUNT, NMJL, NMJU, DELTAM, DELTAJ
    INTEGER,         ALLOCATABLE                      :: ML(:), MU(:)
    ! Number of transitions for each type
    NPI = 2*MIN(JL(IND_LINE(I)),JU(IND_LINE(I)))+1    ! Pi-components
    NSR = JL(IND_LINE(I))+JU(IND_LINE(I))             ! Sigma-red components
    NSB = JL(IND_LINE(I))+JU(IND_LINE(I))             ! Sigma-blue components 
    NTOT = NPI + NSR + NSB
    DELTAJ = JU(IND_LINE(I))-JL(IND_LINE(I))
    ! Selection rule in quantum number J
    IF (ABS(DELTAJ).GT.1) THEN
       PRINT*,'Selection rules imply that DELTAJ=-1,0,+1. However,'
       PRINT*,'in '//TRIM(INPUTFILE)//' the spectral line', LINE_L0(IND_LINE(I))
       PRINT*,'does not verify such rule. Please correct the electronic configuration. STOP.'
       STOP
    ENDIF
    !
    ! We use double precision because this will be employed to determine wavelength shifts
    !
    SLOW=REAL(SL(IND_LINE(I)))
    SUPP=REAL(SU(IND_LINE(I)))
    LLOW=REAL(LL(IND_LINE(I)))
    LUPP=REAL(LU(IND_LINE(I)))
    JLOW=REAL(JL(IND_LINE(I)))
    JUPP=REAL(JU(IND_LINE(I)))
    ! Lande factor lower level
    IF (JLOW.NE.0) THEN
       GLOW=1.5+(SLOW*(SLOW+1.)-LLOW*(LLOW+1.))/(2.*JLOW*(JLOW+1.))
    ELSE
       GLOW=0.
    ENDIF
    ! Lande factor upper level
    IF (JUPP.NE.0) THEN
       GUPP=1.5+(SUPP*(SUPP+1.)-LUPP*(LUPP+1.))/(2.*JUPP*(JUPP+1.))
    ELSE
       GUPP=0.
    ENDIF
    ! Effective Lande factor
    GEFF=0.5*(GLOW+GUPP)+0.25*(GLOW-GUPP)*(JLOW*(JLOW+1.)-JUPP*(JUPP+1.))
    ! Allocate array for wavelength shifts
    ALLOCATE( LAMBDA_SC(NTOT), STRENGTH_SC(NTOT), TYPE_SC(NTOT))
    ! Allocate array for number of M-sublevels that split out of one J level
    NMJL=2*JL(IND_LINE(I))+1
    NMJU=2*JU(IND_LINE(I))+1
    ALLOCATE (ML(NMJL), MU(NMJU))
    ! Create M quantum numbers
    DO J=1,NMJL
       ML(J)=J-1-JLOW
    ENDDO
    DO J=1,NMJU
       MU(J)=J-1-JUPP
    ENDDO
    ! Start Loop to determine allowed transitions and their wavelength shift and strength
    COUNT = 0
    DO J=1,NMJL
       DO K=1,NMJU
          DELTAM=MU(K)-ML(J)
          IF (ABS(DELTAM).LE.1) THEN
             COUNT=COUNT+1
             ! Units are mA/Gauss
             DELTALB=4.66685E-10*(ML(J)*GLOW-MU(K)*GUPP)*LINE_L0(IND_LINE(I))**2.
             TYPE_SC(COUNT)=DELTAM
             LAMBDA_SC(COUNT)=DELTALB
             SELECT CASE (DELTAJ)
                !
             CASE(-1)
                SELECT CASE (DELTAM)
                CASE(-1)
                   STRENGTH_SC(COUNT) = (JLOW+MU(K))*(JUPP+MU(K)+2.)
                CASE(0)
                   STRENGTH_SC(COUNT) = 2.*(JLOW**2.-MU(K)**2.)
                CASE(1)
                   STRENGTH_SC(COUNT) = (JLOW-MU(K))*(JUPP-MU(K)+2.)
                END SELECT
                !
             CASE(0)
                SELECT CASE (DELTAM)
                CASE(-1)
                   STRENGTH_SC(COUNT) = (JUPP-MU(K))*(JUPP+MU(K)+1.)
                CASE(0)
                   STRENGTH_SC(COUNT) = 2.*MU(K)**2.
                CASE(1)
                   STRENGTH_SC(COUNT) = (JUPP+MU(K))*(JUPP-MU(K)+1.)
                END SELECT
             CASE(1)
                SELECT CASE (DELTAM)
                CASE(-1)
                   STRENGTH_SC(COUNT) = (JUPP-MU(K))*(JLOW-MU(K))
                CASE(0)
                   STRENGTH_SC(COUNT) = 2.*(JUPP**2.-MU(K)**2.)
                CASE(1)
                   STRENGTH_SC(COUNT) = (JUPP+MU(K))*(JUPP+MU(K))
                END SELECT
             END SELECT
          ENDIF
       ENDDO
    ENDDO
    !
    DEALLOCATE( ML, MU )
    !
    !
    ! Normalize the strength of the Zeeman pattern individually for PI, SIGMA_BLUE and SIGMA_RED components
    ! 
    SUM_PI=0.
    SUM_SR=0.
    SUM_SB=0.
    DO J=1,NTOT
       SELECT CASE (TYPE_SC(J))
       CASE(-1)
          SUM_SR=SUM_SR+STRENGTH_SC(J)
       CASE(0)
          SUM_PI=SUM_PI+STRENGTH_SC(J)
       CASE(1)
          SUM_SB=SUM_SB+STRENGTH_SC(J)
       END SELECT
    ENDDO
    !
    DO J=1,NTOT
       SELECT CASE (TYPE_SC(J))
       CASE(-1)
          STRENGTH_SC(J)=STRENGTH_SC(J)/SUM_SR
       CASE(0)
          STRENGTH_SC(J)=STRENGTH_SC(J)/SUM_PI
       CASE(1)
          STRENGTH_SC(J)=STRENGTH_SC(J)/SUM_SB
       END SELECT
    ENDDO
    !
  END SUBROUTINE ZEEMAN_COMP
  !--------------------------------------------
  !--------------------------------------------
END MODULE MAGSPLIT
