MODULE DERIV_TEST
  ! J M Borrero
  ! September 6, 2016
  ! KIS, Freiburg
  USE CONS_PARAM
  USE LINES_DATABASE
  USE ATOM_DATABASE
  USE SIMUL_BOX
  USE CHEMICAL_EQUILIBRIUM
  USE DAMPING
CONTAINS
  ! deriv_ne_rho
  ! deriv_nh_pg
  ! deriv_damp
  ! deriv_nir
  SUBROUTINE DERIV_NE_RHO
    ! calculates derivative of ne with respect to T at constant density
    IMPLICIT NONE
    INTEGER                       :: I, J, K
    REAL(SP)                      :: T_1, T_2, D_1, D_2
    REAL(DP)                      :: NELEC_1, DNELEC_1, NELEC_2, DNELEC_2
    !
    PRINT*,'   NUMERICAL                  ANALYTICAL 1                ANALYTICAL 2'
    DO I=1,NX
       DO J=1,NY
          DO K=1,NZ
             ! First atmosphere
             T_1=TEM(I,J,K)
             D_1=RHO(I,J,K)
             CALL GET_NE_RHO(T_1, D_1, NELEC_1, DNELEC_1)
             ! Perturbed atmosphere
             T_2=TEM(I,J,K)*1.003D0
             D_2=D_1
             CALL GET_NE_RHO(T_2, D_2, NELEC_2, DNELEC_2)
             PRINT*,(NELEC_2-NELEC_1)/(T_2-T_1),DNELEC_1,DNELEC_2
          ENDDO
       ENDDO
    ENDDO
    STOP
  END SUBROUTINE DERIV_NE_RHO
  !---------------------------------------------------
  !---------------------------------------------------
  SUBROUTINE DERIV_NE_PG
    ! calculates derivative of ne with respect to T at constant gas pressure
    IMPLICIT NONE
    INTEGER                       :: I, J, K
    REAL(SP)                      :: T_1, T_2, P_1, P_2, D_1, D_2, M_1, M_2
    REAL(SP)                      :: T_3, T_4, P_3, P_4, D_3, D_4, M_3, M_4
    REAL(DP)                      :: NELEC_1, DNELEC_1,  NELEC_2, DNELEC_2 
    REAL(DP)                      :: NELEC_3, DNELEC_3,  NELEC_4, DNELEC_4
    REAL(DP)                      :: NELEC, DNELEC
    !
    DO I=1,NX
       DO J=1,NY
          DO K=1,NZ
             ! First atmosphere
             T_1=TEM(I,J,K)
             P_1=PG(I,J,K)
             D_1=RHO(I,J,K)
             M_1=MW(I,J,K)
             ! This gives the derivative of ne with respect to PG at constant temperature
             CALL GET_NE_TEMP(T_1, P_1, NELEC_1, DNELEC_1)
             ! Perturbed atmosphere
             T_2=TEM(I,J,K)
             P_2=1.003*PG(I,J,K)
             D_2=RHO(I,J,K)
             M_2=MW(I,J,K)
             CALL GET_NE_TEMP(T_2, P_2, NELEC_2, DNELEC_2)
             !PRINT*,'   NUMERICAL                  ANALYTICAL 1                ANALYTICAL 2'
             !PRINT*,(NELEC_2-NELEC_1)/(P_2-P_1),DNELEC_1,DNELEC_2
             ! Now the derivative of ne with respect to T at constant density
             T_3=TEM(I,J,K)
             P_3=PG(I,J,K)
             D_3=RHO(I,J,K)
             M_3=MW(I,J,K)
             CALL GET_NE_RHO(T_3, D_3, NELEC_3, DNELEC_3)
             T_4=1.003*TEM(I,J,K)
             P_4=PG(I,J,K)
             D_4=RHO(I,J,K)
             M_4=MW(I,J,K)
             CALL GET_NE_RHO(T_4, D_4, NELEC_4, DNELEC_4)
             !PRINT*,'   NUMERICAL                  ANALYTICAL 1                ANALYTICAL 2'
             !PRINT*,(NELEC_4-NELEC_3)/(T_4-T_3),DNELEC_3,DNELEC_4
             ! Finally, derivative of ne with repect to T at constant gas pressure
             PRINT*,'   ANALITYCAL                  NUMERICAL'
             PRINT*,DNELEC_3-DNELEC_1*D_1*KBOL/(M_1*MHYD),(NELEC_4-NELEC_3)/(T_4-T_3)- & 
                  (NELEC_2-NELEC_1)/(P_2-P_1)*D_1*KBOL/(M_1*MHYD)
             !
          ENDDO
       ENDDO
    ENDDO
    STOP
  END SUBROUTINE DERIV_NE_PG
  !-----------------------------------------------
  !-----------------------------------------------
  SUBROUTINE DERIV_DAMP
    ! derivative of the damping with respect to temperature
    IMPLICIT NONE
    INTEGER                       :: I, J, K
    REAL(SP)                      :: T_1, T_2, P_1, D_1, M_1
    REAL(SP)                      :: DAM_1, DAM_2
    REAL(DP)                      :: NH_1, NH_2, NE_1, NE_2
    REAL(DP)                      :: DNHDT_P, DNHDT_R
    REAL(DP)                      :: DDAMDT_PG_1, DDAMDT_RHO_1
    REAL(DP)                      :: DDAMDT_PG_2, DDAMDT_RHO_2
    ! Derivative at constant gas pressure
    DO I=1,NX
       DO J=1,NY
          DO K=1,NZ
             T_1=TEM(I,J,K)
             P_1=PG(I,J,K)
             D_1=RHO(I,J,K)
             CALL GET_NHNEMU(T_1,P_1,D_1,NH_1,NE_1,M_1)
             CALL GET_DNHYDDTEMP(T_1,P_1,D_1,M_1,DNHDT_P,DNHDT_R)
             CALL GET_DAMPING(LINE_ZN(IND_LINE(1)),LINE_L0(IND_LINE(1))*1E-8,LINE_ION(IND_LINE(1)) &
                  ,EPLOW(IND_LINE(1)), NH_1, T_1, SIGMA(IND_LINE(1)), ALPHA(IND_LINE(1)),DAM_1 &
                  ,DDAMDT_PG_1, DDAMDT_RHO_1)
             T_2=1.003*TEM(I,J,K)
             NH_2=NH_1+DNHDT_P*(T_2-T_1)
             CALL GET_DAMPING(LINE_ZN(IND_LINE(1)),LINE_L0(IND_LINE(1))*1E-8,LINE_ION(IND_LINE(1)) &
                  ,EPLOW(IND_LINE(1)), NH_2, T_2, SIGMA(IND_LINE(1)), ALPHA(IND_LINE(1)),DAM_2 &
                  ,DDAMDT_PG_2, DDAMDT_RHO_2)
             PRINT*,(DAM_2-DAM_1)/(T_2-T_1),DDAMDT_PG_1,DDAMDT_PG_2
          ENDDO
       ENDDO
    ENDDO
    PRINT*,'----------------------------------------------------------------'
    ! Derivative at constant density
    DO I=1,NX
       DO J=1,NY
          DO K=1,NZ
             T_1=TEM(I,J,K)
             P_1=PG(I,J,K)
             D_1=RHO(I,J,K)
             CALL GET_NHNEMU(T_1,P_1,D_1,NH_1,NE_1,M_1)
             CALL GET_DNHYDDTEMP(T_1,P_1,D_1,M_1,DNHDT_P,DNHDT_R)
             CALL GET_DAMPING(LINE_ZN(IND_LINE(1)),LINE_L0(IND_LINE(1))*1E-8,LINE_ION(IND_LINE(1)) &
                  ,EPLOW(IND_LINE(1)), NH_1, T_1, SIGMA(IND_LINE(1)), ALPHA(IND_LINE(1)),DAM_1 &
                  ,DDAMDT_PG_1, DDAMDT_RHO_1)
             T_2=1.003*TEM(I,J,K)
             NH_2=NH_1+DNHDT_R*(T_2-T_1)
             CALL GET_DAMPING(LINE_ZN(IND_LINE(1)),LINE_L0(IND_LINE(1))*1E-8,LINE_ION(IND_LINE(1)) &
                  ,EPLOW(IND_LINE(1)), NH_2, T_2, SIGMA(IND_LINE(1)), ALPHA(IND_LINE(1)),DAM_2 &
                  ,DDAMDT_PG_2, DDAMDT_RHO_2)
             PRINT*,(DAM_2-DAM_1)/(T_2-T_1),DDAMDT_RHO_1,DDAMDT_RHO_2
          ENDDO
       ENDDO
    ENDDO
    STOP
  END SUBROUTINE DERIV_DAMP
  !-----------------------------------------------
  !-----------------------------------------------
  SUBROUTINE DERIV_NIR
    ! Calculates the derivative of the number of atoms of specie I with ionization stage R with respect to temperature
    IMPLICIT NONE
    INTEGER                                    :: I, J, K
    REAL(SP)                                   :: T_1, P_1, R_1, D_1, M_1, T_2
    REAL(DP)                                   :: NH_1, NE_1, NI_1, NIR_1, NH_2, NE_2, NI_2, NIR_2
    REAL(DP)                                   :: DNIDT_PG_1, DNIDT_RHO_1, DNIDT_PG_2, DNIDT_RHO_2
    REAL(DP)                                   :: DNEDP_T, DNEDT_R, DNEDT_P, DNHDT_P, DNHDT_R
    ! Constant gas pressure
    DO I=1,NX
       DO J=1,NY
          DO K=1,NZ
             T_1=TEM(I,J,K)
             P_1=PG(I,J,K)
             D_1=RHO(I,J,K)
             CALL GET_NHNEMU(T_1,P_1,D_1,NH_1,NE_1,M_1)
             CALL GET_NE_TEMP(T_1,P_1,NE_1,DNEDP_T)
             CALL GET_NE_RHO(T_1,D_1,NE_1,DNEDT_R)
             !
             DNEDT_P=DNEDT_R-((DBLE(D_1)*KBOL)/(DBLE(M_1)*MHYD))*DNEDP_T
             !
             CALL GET_DNHYDDTEMP(T_1,P_1,D_1,M_1,DNHDT_P,DNHDT_R)
             NI_1=NH_1*10D0**(ABUND(LINE_ZN(IND_LINE(1)))-12D0)
             CALL ION_CALC_TEMPX(IND_LINE(1),T_1,NE_1,NI_1,NIR_1,DNIDT_PG_1,DNIDT_RHO_1)
             T_2=1.003*TEM(I,J,K)
             NE_2=NE_1+DNEDT_P*(T_2-T_1)
             NH_2=NH_1+DNHDT_P*(T_2-T_1)
             NI_2=NH_2*10D0**(ABUND(LINE_ZN(IND_LINE(1)))-12D0)
             CALL ION_CALC_TEMPX(IND_LINE(1),T_2,NE_2,NI_2,NIR_2,DNIDT_PG_2,DNIDT_RHO_2)
             PRINT*,(NIR_2-NIR_1)/(T_2-T_1),DNIDT_PG_1,DNIDT_PG_2
          ENDDO
       ENDDO
    ENDDO
    PRINT*,'------------------------------------------------------------'
    ! Constant density
    DO I=1,NX
       DO J=1,NY
          DO K=1,NZ
             T_1=TEM(I,J,K)
             P_1=PG(I,J,K)
             D_1=RHO(I,J,K)
             CALL GET_NHNEMU(T_1,P_1,D_1,NH_1,NE_1,M_1)
             CALL GET_NE_RHO(T_1,D_1,NE_1,DNEDT_R)
             CALL GET_DNHYDDTEMP(T_1,P_1,D_1,M_1,DNHDT_P,DNHDT_R)
             NI_1=NH_1*10D0**(ABUND(LINE_ZN(IND_LINE(1)))-12D0)
             CALL ION_CALC_TEMPX(IND_LINE(1),T_1,NE_1,NI_1,NIR_1,DNIDT_PG_1,DNIDT_RHO_1)
             T_2=1.003*TEM(I,J,K)
             NE_2=NE_1+DNEDT_R*(T_2-T_1)
             NH_2=NH_1+DNHDT_R*(T_2-T_1)
             NI_2=NH_2*10D0**(ABUND(LINE_ZN(IND_LINE(1)))-12D0)
             CALL ION_CALC_TEMPX(IND_LINE(1),T_2,NE_2,NI_2,NIR_2,DNIDT_PG_2,DNIDT_RHO_2)
             PRINT*,(NIR_2-NIR_1)/(T_2-T_1),DNIDT_RHO_1,DNIDT_RHO_2
          ENDDO
       ENDDO
    ENDDO
    !
    STOP
  END SUBROUTINE DERIV_NIR
  !-----------------------------------------------
  !-----------------------------------------------


END MODULE DERIV_TEST
