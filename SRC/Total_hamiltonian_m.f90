MODULE Total_hamiltonian_m
  !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
  USE QDUtil_m
  USE Operator_1D_m
  IMPLICIT NONE

  
  CONTAINS
  
    
  SUBROUTINE MolecCav_Action_Total_Hamiltonian_1D(Result_total_WF, Matter_hamiltonianSystem_WF, &
                                                 & Cavity_hamiltonian_1D, Cavity_position_1D, &
                                                 & Matter_dipolar_moment, System_WF, &
                                                 & lambda_cavity_mode, w_cavity_mode)
    USE QDUtil_m
    USE Operator_1D_m
    IMPLICIT NONE

    real(kind=Rkind), intent(inout)    :: Result_total_WF(:,:)                 ! already allocated !
    real(kind=Rkind), intent(in)       :: Matter_hamiltonianSystem_WF(:,:)     ! size Nb_M*Nb_C. |H_MatterSystem_WF(:,i_C)> = H_Matter|System_WF(:,i_C)>
    TYPE(Operator_1D_t), intent(in)    :: Cavity_hamiltonian_1D
    TYPE(Operator_1D_t), intent(in)    :: Cavity_position_1D
    TYPE(Operator_1D_t), intent(in)    :: Matter_dipolar_moment                ! \hat{\mu}_{M}(R) = Cte.\hat{R} selon hypothèses
    real(kind=Rkind), intent(in)       :: System_WF(:,:)                       ! size Nb_M*Nb_C
    real(kind=Rkind), intent(in)       :: lambda_cavity_mode, w_cavity_mode    ! coupling strenght and eigenpulsation

    real(kind=Rkind), allocatable      :: Cavity_hamiltonian_1DSystem_WF(:)    ! only one line of System_WF
    real(kind=Rkind), allocatable      :: Intermediary(:,:)
    real(kind=Rkind), allocatable      :: Matter_cavity_coupling_hamiltonian_1DSystem_WF(:,:) ! only one line of System_WF
    integer                            :: Nb_M, Nb_C, i_M, i_C
  
    Nb_M = Size(System_WF, 1)
    Nb_C = Size(System_WF, 2)

    Result_total_WF = ZERO

    !-----------H_tot = H_Matter + H_Cavity + H_MatterCavityCoupling-----------
    !----------H_Matter|System_WF> = H_Matter(System_WF) already known---------
    Result_total_WF = Matter_hamiltonianSystem_WF
  
    !-----------H_Cavity|System_WF> = H_Cavity(System_WF) to compute-----------
    ALLOCATE(Cavity_hamiltonian_1DSystem_WF(Nb_C))
    DO i_M = 1, Nb_M
      CALL Action_Operator_1D(Cavity_hamiltonian_1DSystem_WF, &
                                     & Cavity_hamiltonian_1D, &
                                     & System_WF(i_M,:))
      Result_total_WF(i_M,:) = Result_total_WF(i_M,:) + Cavity_hamiltonian_1DSystem_WF(:)
    END DO

    !---------------------H_MatterCavityCoupling|System_WF>--------------------
    !------------------------H_CavityCoupling(System_WF)-----------------------
    ALLOCATE(Intermediary(Nb_M,Nb_C))
    Intermediary = ZERO
    DO i_M = 1, Nb_M
      CALL Action_Operator_1D(Intermediary(i_M,:), &
                                     & Cavity_position_1D, &
                                     & System_WF(i_M,:))
    END DO
    Intermediary = lambda_cavity_mode*w_cavity_mode*Intermediary

    !---------H_MatterCoupling(Intermediary) = Cte.pos_op of the matter--------
    ALLOCATE(Matter_cavity_coupling_hamiltonian_1DSystem_WF(Nb_M,Nb_C))
    CALL MolecCav_Action_matter_dipolar_moment_1D(Matter_cavity_coupling_hamiltonian_1DSystem_WF, &
                                                & Matter_dipolar_moment, &
                                                & Intermediary)                ! the matter_dipolar_moment is assumed to already contains its intensity constant (\frac{d\mu}{dq}) within its matrix /!\
                                            
    !-----------------------------H_tot = summation----------------------------
    Result_total_WF = Result_total_WF + Matter_cavity_coupling_hamiltonian_1DSystem_WF

    DEALLOCATE(Intermediary)

  
  END SUBROUTINE
  

  SUBROUTINE MolecCav_Action_matter_dipolar_moment_1D(Psi_result, Matter_dipolar_moment, Psi_argument) ! /!\ only for the app i.e. for a diatomic molecule with the approximation of a linear dip. mom. with the matter position /!\
    USE QDUtil_m
    USE Operator_1D_m
    IMPLICIT NONE

    real(kind=Rkind),       intent(inout) :: Psi_result(:,:)
    TYPE(Operator_1D_t),    intent(in)    :: Matter_dipolar_moment             ! supposed to contains the constant (\frac{d\mu}{dq}) within its matrix
    real(kind=Rkind),       intent(in)    :: Psi_argument(:,:)

    integer                               :: Nb_C, i_C 

    Nb_C = Size(Psi_result, 2)

    IF (Size(Psi_argument, 1) /= Size(Psi_result, 1) .OR. Size(Psi_argument, 2) /= Nb_C) THEN
      STOP "The size of Psi_argument does not match the size of Psi_result"
    END IF

    DO i_C = 1, Nb_C
      CALL Action_Operator_1D(Psi_result(:,i_C), &                    ! just a matricial product vector by vector to abide by the arguments constraints of the subroutine (i think we could replace by matmul(matter_dip_mom), Psi_argu)
                                     & Matter_dipolar_moment, &
                                     & Psi_argument(:,i_C))
    END DO

  END SUBROUTINE


  SUBROUTINE MolecCav_Construct_H_tot(H_tot, Nb_M, Nb_C, Matter_hamiltonianSystem_WF, &
                                    & Cavity_hamiltonian_1D, Cavity_position_1D, &
                                    & Matter_dipolar_moment, &
                                    & lambda_cavity_mode, w_cavity_mode)
    USE QDUtil_m
    USE Operator_1D_m
    IMPLICIT NONE

    real(kind=Rkind),    intent(inout) :: H_tot(:,:)
    integer,             intent(in)    :: Nb_M, Nb_C
    real(kind=Rkind),    intent(in)    :: Matter_hamiltonianSystem_WF(:,:)     ! size Nb_M*Nb_C. |H_MatterSystem_WF(:,i_C)> = H_Matter|System_WF(:,i_C)>
    TYPE(Operator_1D_t), intent(in)    :: Cavity_hamiltonian_1D
    TYPE(Operator_1D_t), intent(in)    :: Cavity_position_1D
    TYPE(Operator_1D_t), intent(in)    :: Matter_dipolar_moment                ! \hat{\mu}_{M}(R) = Cte.\hat{R} selon hypothèses
    real(kind=Rkind),    intent(in)    :: lambda_cavity_mode, w_cavity_mode    ! coupling strenght and eigenpulsation


    real(kind=Rkind), allocatable      :: Psi_basis(:,:), Psi_result(:,:)
    integer                            :: NB, I, J, j_M, j_C, i_M, i_C

    J = 0
    NB = Size(H_tot, 1)
    H_tot(:,:) = ZERO
    ALLOCATE(Psi_basis(Nb_M, Nb_C))
    ALLOCATE(Psi_result(Nb_M, Nb_C))

    IF (NB /= Size(H_tot, 2) .OR. NB /= Nb_M*Nb_C) THEN
      STOP "wrong allocation of the H_tot matrix"
    END IF

    DO j_M = 1, Nb_M
      DO j_C = 1, Nb_C
        J = J + 1
        I = 0
        Psi_basis = ZERO
        Psi_basis(j_M, j_C) = ONE
        Psi_result = ZERO

        CALL MolecCav_Action_Total_Hamiltonian_1D(Psi_result, Matter_hamiltonianSystem_WF, &
                                                  & Cavity_hamiltonian_1D, Cavity_position_1D, &
                                                  & Matter_dipolar_moment, &
                                                  & Psi_basis, lambda_cavity_mode, w_cavity_mode)
          
        DO i_M = 1, Nb_M
          DO i_C = 1, Nb_C
            I = I + 1
            H_tot(I,J) = Psi_result(i_M, i_C)
          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE


  SUBROUTINE MolecCav_Mapping_WF_2DTO1D(Psi_1D, Psi_2D)
    USE QDUtil_m
    IMPLICIT NONE

    real(kind=Rkind), intent(inout)    :: Psi_1D(:)
    real(kind=Rkind), intent(in)       :: Psi_2D(:,:)

    integer                            :: Nb_M, Nb_C, i_M, i_C, NB, I

    Nb_M   = Size(Psi_2D,1)
    Nb_C   = Size(Psi_2D,2)
    NB     = Size(Psi_1D)
    I      = 0
    Psi_1D = ZERO

    IF (Nb_M*Nb_C /= NB) THEN
      STOP "Wrong size of the matrices"
    END IF

    DO i_M = 1, Nb_M
      DO i_C = 1, Nb_C
        I = I + 1
        Psi_1D(I) = Psi_2D(i_M, i_C)
      END DO
    END DO

  END SUBROUTINE


  SUBROUTINE MolecCav_Average_value_H_tot(Value, H_tot, Psi_argument)   ! /!\ FOR NOW EVERYTHING IS REAL /!\ compute the resulting vector Psi_result(:) from the action of the operator of the cavity mode on the photon state vector Psi_argument(:) written in the Eigenbasis of H_ho
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    IMPLICIT NONE
  
    real(kind=Rkind), intent(inout)    :: Value
    real(kind=Rkind), intent(in)       :: H_tot(:,:)    
    real(kind=Rkind), intent(in)       :: Psi_argument(:)
  
    real(kind=Rkind), allocatable      :: Intermediary(:)
    integer                            :: NB
  
    NB = Size(Psi_argument)
    ALLOCATE(Intermediary(NB))

    Intermediary(:) = MATMUL(H_tot, Psi_argument)
    Value = DOT_PRODUCT(Psi_argument, Intermediary) 

    DEALLOCATE(Intermediary)
  
  END SUBROUTINE
  

END MODULE
  