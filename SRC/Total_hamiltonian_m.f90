MODULE Total_hamiltonian_m
  !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
  USE QDUtil_m
  IMPLICIT NONE

  
  PRIVATE

  PUBLIC Action_total_hamiltonian_1p1D, Construct_total_hamiltonian_1p1D, Mapping_WF_2DTO1D, Average_value_H_tot

  INTERFACE Action_total_hamiltonian_1p1D
    MODULE PROCEDURE MolecCav_Action_total_hamiltonian_1p1D_old, MolecCav_Action_total_hamiltonian_1p1D
  END INTERFACE
  INTERFACE Action_Matter_1p1D
    MODULE PROCEDURE MolecCav_Action_Matter_1p1D
  END INTERFACE
  INTERFACE Action_Cavity_1p1D
    MODULE PROCEDURE MolecCav_Action_Cavity_1p1D
  END INTERFACE
  INTERFACE Action_matter_dipolar_moment_1p1D
    MODULE PROCEDURE MolecCav_Action_matter_dipolar_moment_1p1D
  END INTERFACE
  INTERFACE Construct_total_hamiltonian_1p1D
  MODULE PROCEDURE MolecCav_Construct_total_hamiltonian_1p1D
  END INTERFACE
  INTERFACE Mapping_WF_2DTO1D
  MODULE PROCEDURE MolecCav_Mapping_WF_2DTO1D
  END INTERFACE
  INTERFACE Average_value_H_tot
  MODULE PROCEDURE MolecCav_Average_value_H_tot
  END INTERFACE
    

  CONTAINS
  
    
  SUBROUTINE MolecCav_Action_total_hamiltonian_1p1D(TotH_psi, CavPosition, CavH, Mat_dipolar_moment, MatH, Psi, Debug)
    USE QDUtil_m
    USE Cavity_mode_m
    USE Operator_1D_m
    IMPLICIT NONE

    real(kind=Rkind),    intent(inout) :: TotH_psi(:,:)                        ! already allocated !
    TYPE(Operator_1D_t), intent(in)    :: CavPosition                          ! position operator associated to the cavity mode
    TYPE(Operator_1D_t), intent(in)    :: CavH                                 ! Hamiltonian associated to the cavity mode
    TYPE(Operator_1D_t), intent(in)    :: Mat_dipolar_moment                   ! dipolar moment operator associated to the matter mode. \hat{\mu}_{M}(R) = Cte.\hat{R} selon hypothèses
    TYPE(Operator_1D_t), intent(in)    :: MatH                                 ! Hamiltonian associated to the matter mode
    real(kind=Rkind),    intent(in)    :: Psi(:,:)                             ! size Nb_M*Nb_C
    logical, optional,   intent(in)    :: Debug

    real(kind=Rkind), allocatable      :: Psi_1(:,:), Psi_2(:,:), Psi_3(:,:)   ! cf. below for meaning
    integer                            :: Nb_M, Nb_C
    logical                            :: Debug_local = .TRUE.
  
    ! TotH_psi = [MatH x CavIdentity]_psi + [MatIdentity x CavH]_psi + lambda_C.w_C.[Mat_dipolar_moment x CavPosition]_psi
    !          = [CavIdentity]_psi_1 + [CavH]_psi_2 + lambda_C.w_C.[CavPosition]_psi_3
    !    Psi_1 = [MatH]_psi
    !    Psi_2 = [MatIdentity]_psi
    !    Psi_3 = [Mat_dipolar_moment]_psi

    !------------------------------Initialization------------------------------
    IF (PRESENT(Debug)) Debug_local = Debug
    Nb_M     = Size(Psi, 1)
    Nb_C     = Size(Psi, 2)

    ALLOCATE(Psi_1(Nb_M, Nb_C))
    ALLOCATE(Psi_2(Nb_M, Nb_C))
    ALLOCATE(Psi_3(Nb_M, Nb_C))
    Psi_1    = ZERO
    Psi_2    = ZERO
    Psi_3    = ZERO
    TotH_psi = ZERO

    !----------------------------Checking dimensions---------------------------
    IF (Nb_M /= Size(TotH_psi,1) .OR. Nb_C /= Size(TotH_psi,2)) THEN
      WRITE(out_unit,*) "The dimensions of the operand Psi's matrix does not match the resulting TotH_psi matrix's& 
                       & size. Please check initialization."
      STOP "The dimensions of the operand Psi's matrix does not match the resulting TotH_psi matrix's size. Please&
                       & check initialization."
    END IF 

    IF (Nb_M /= MatH%Nb) THEN
      WRITE(out_unit,*) "The dimensions of the operand Psi's matrix does not match the dimensions of the Operator &
                       &MatH's matrix representation. Please check initialization."
      STOP "The dimensions of the operand Psi's matrix does not match the dimensions of the Operator MatH's matrix&
                       & representation. Please check initialization."
    END IF 

    IF (Nb_M /= Mat_dipolar_moment%Nb) THEN
      WRITE(out_unit,*) "The dimensions of the operand Psi's matrix does not match the dimensions of the Operator &
                       &Mat_dipolar_moment's matrix representation. Please check initialization."
      STOP "The dimensions of the operand Psi's matrix does not match the dimensions of the Operator Mat_dipolar_m&
                       &oment's matrix representation. Please check initialization."
    END IF 

    IF (Nb_C /= CavH%Nb) THEN
      WRITE(out_unit,*) "The dimensions of the operand Psi's matrix does not match the dimensions of the Operator &
                       &CavH's matrix representation. Please check initialization."
      STOP "The dimensions of the operand Psi's matrix does not match the dimensions of the Operator CavH's matrix&
                       & representation. Please check initialization."
    END IF 

    IF (Nb_C /= CavPosition%Nb) THEN
      WRITE(out_unit,*) "The dimensions of the operand Psi's matrix does not match the dimensions of the Operator &
                       &CavPosition's matrix representation. Please check initialization."
      STOP "The dimensions of the operand Psi's matrix does not match the dimensions of the Operator CavPosition's&
                       & matrix representation. Please check initialization."
    END IF 

    !--------------------------------Computation-------------------------------
    CALL Action_Matter_1p1D(Psi_1, Psi_2, Psi_3, Mat_dipolar_moment, MatH, Psi)

    IF (Debug_local .AND. Nb_C <= 10) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "-----------------Computing the action of the total Hamiltonian over the 1p1D system-----------------"
      CALL Write_Mat(Psi,   out_unit, Nb_C, info="Psi")
      WRITE(out_unit,*)
      CALL Write_Mat(Psi_1, out_unit, Nb_C, info="Psi_1")
      WRITE(out_unit,*)
      CALL Write_Mat(Psi_2, out_unit, Nb_C, info="Psi_2")
      WRITE(out_unit,*)
      CALL Write_Mat(Psi_3, out_unit, Nb_C, info="Psi_3")
    ELSE IF (Debug_local .AND. Nb_C > 10) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "-----------------Computing the action of the total Hamiltonian over the 1p1D system-----------------"
      CALL Write_Mat(Psi(1:10,1:10),   out_unit, 10, info="Psi(10:10sliced)")
      WRITE(out_unit,*)
      CALL Write_Mat(Psi_1(1:10,1:10), out_unit, 10, info="Psi_1(10:10sliced)")
      WRITE(out_unit,*)
      CALL Write_Mat(Psi_2(1:10,1:10), out_unit, 10, info="Psi_2(10:10sliced)")
      WRITE(out_unit,*)
      CALL Write_Mat(Psi_3(1:10,1:10), out_unit, 10, info="Psi_3(10:10sliced)")
    END IF

    CALL Action_Cavity_1p1D(TotH_psi, CavPosition, CavH, Psi_1, Psi_2, Psi_3, Debug_local)

  END SUBROUTINE MolecCav_Action_total_hamiltonian_1p1D
  

  SUBROUTINE MolecCav_Action_Matter_1p1D(Psi_1, Psi_2, Psi_3, Mat_dipolar_moment, MatH, Psi)
    USE QDUtil_m
    USE Cavity_mode_m
    USE Operator_1D_m
    IMPLICIT NONE

    real(kind=Rkind),    intent(inout) :: Psi_1(:,:)                           ! cf. below for meaning
    real(kind=Rkind),    intent(inout) :: Psi_2(:,:)                           ! cf. below for meaning
    real(kind=Rkind),    intent(inout) :: Psi_3(:,:)                           ! cf. below for meaning
    TYPE(Operator_1D_t), intent(in)    :: Mat_dipolar_moment                   ! dipolar moment operator associated to the matter mode. \hat{\mu}_{M}(R) = Cte.\hat{R} selon hypothèses
    TYPE(Operator_1D_t), intent(in)    :: MatH                                 ! Hamiltonian associated to the matter mode
    real(kind=Rkind),    intent(in)    :: Psi(:,:)                             ! size Nb_M*Nb_C

    integer                            :: i_C
  
    ! TotH_psi = [MatH x CavIdentity]_psi + [MatIdentity x CavH]_psi + lambda_C.w_C.[Mat_dipolar_moment x CavPosition]_psi
    !          = [CavIdentity]_psi_1 + [CavH]_psi_2 + lambda_C.w_C.[CavPosition]_psi_3
    !    Psi_1 = [MatH]_psi
    !    Psi_2 = [MatIdentity]_psi
    !    Psi_3 = [Mat_dipolar_moment]_psi

    DO i_C = 1, Size(Psi, 2)                                                   ! initialize Matter_hamiltonianSystem_WF by applying the matter hamiltonian to each column of the matrix of the total system WF Psi
      CALL Action_Operator_1D(Psi_1(:,i_C), MatH, Psi(:,i_C))
    END DO
    
    Psi_2(:,:) = Psi(:,:)

    CALL Action_matter_dipolar_moment_1p1D(Psi_3, Mat_dipolar_moment, Psi)

  END SUBROUTINE MolecCav_Action_Matter_1p1D

  
  SUBROUTINE MolecCav_Action_Cavity_1p1D(TotH_psi, CavPosition, CavH, Psi_1, Psi_2, Psi_3, Debug_local)
    USE QDUtil_m
    USE Cavity_mode_m
    USE Operator_1D_m
    IMPLICIT NONE

    real(kind=Rkind),    intent(inout) :: TotH_psi(:,:)                        ! already allocated !
    TYPE(Operator_1D_t), intent(in)    :: CavPosition                          ! position operator associated to the cavity mode
    TYPE(Operator_1D_t), intent(in)    :: CavH                                 ! Hamiltonian associated to the cavity mode
    real(kind=Rkind),    intent(in)    :: Psi_1(:,:)                           ! cf. below for meaning
    real(kind=Rkind),    intent(in)    :: Psi_2(:,:)                           ! cf. below for meaning
    real(kind=Rkind),    intent(in)    :: Psi_3(:,:)                           ! cf. below for meaning
    logical,             intent(in)    :: Debug_local

    real(kind=Rkind), allocatable      :: Intermediary(:,:)
    integer                            :: Nb_C, i_M
  
    ! TotH_psi = [MatH x CavIdentity]_psi + [MatIdentity x CavH]_psi + lambda_C.w_C.[Mat_dipolar_moment x CavPosition]_psi
    !          = [CavIdentity]_psi_1 + [CavH]_psi_2 + lambda_C.w_C.[CavPosition]_psi_3
    !    Psi_1 = [MatH]_psi
    !    Psi_2 = [MatIdentity]_psi
    !    Psi_3 = [Mat_dipolar_moment]_psi

    Nb_C = Size(TotH_psi,2)
    
    TotH_psi(:,:) = Psi_1(:,:)
    IF (Debug_local .AND. Nb_C <= 10) THEN
      WRITE(out_unit,*)
      CALL Write_Mat(TotH_psi, out_unit, Nb_C, info="[MatHxCavIdentity]_psi")
    ELSE IF (Debug_local .AND. Nb_C > 10) THEN
      WRITE(out_unit,*)
      CALL Write_Mat(TotH_psi(1:10,1:10), out_unit, 10, info="[MatHxCavIdentity]_psi(10:10sliced)")
    END IF

    ALLOCATE(Intermediary(Size(TotH_psi,1), Size(TotH_psi,2)))
    Intermediary = ZERO
    DO i_M = 1, Size(TotH_psi,1)
      CALL Action_Operator_1D(Intermediary(i_M,:), CavH, Psi_2(i_M,:))
    END DO
    IF (Debug_local .AND. Nb_C <= 10) THEN
      WRITE(out_unit,*)
      CALL Write_Mat(Intermediary, out_unit, Nb_C, info="[MatIdentityxCavH]_psi")
    ELSE IF (Debug_local .AND. Nb_C > 10) THEN
      WRITE(out_unit,*)
      CALL Write_Mat(Intermediary(1:10,1:10), out_unit, 10, info="[MatIdentityxCavH]_psi(10:10sliced)")
    END IF

    TotH_psi(:,:) = TotH_psi(:,:) + Intermediary(:,:) 
    IF (Debug_local .AND. Nb_C <= 10) THEN
      WRITE(out_unit,*)
      CALL Write_Mat(TotH_psi, out_unit, Nb_C, info="[MatH+CavH]_psi")
    ELSE IF (Debug_local .AND. Nb_C > 10) THEN
      WRITE(out_unit,*)
      CALL Write_Mat(TotH_psi(1:10,1:10), out_unit, 10, info="[MatH+CavH]_psi(10:10sliced)")
    END IF

    Intermediary = ZERO
    DO i_M = 1, Size(TotH_psi,1)
      CALL Action_Operator_1D(Intermediary(i_M,:), CavPosition, Psi_3(i_M,:))
    END DO
    IF (Debug_local .AND. Nb_C <= 10) THEN
      WRITE(out_unit,*)
      CALL Write_Mat(Intermediary, out_unit, Nb_C, info="[Mat_dipolar_momentxCavPosition]_psi")
    ELSE IF (Debug_local .AND. Nb_C > 10) THEN
      WRITE(out_unit,*)
      CALL Write_Mat(Intermediary(1:10,1:10), out_unit, 10, info="[Mat_dipolar_momentxCavPosition]_psi(10:10sliced)")
    END IF

    TotH_psi(:,:) = TotH_psi(:,:) + CavH%lambda * CavH%w * Intermediary(:,:)
    IF (Debug_local .AND. Nb_C <= 10) THEN
      WRITE(out_unit,*)
      CALL Write_Mat(TotH_psi, out_unit, Nb_C, info="TotH_psi")
      WRITE(out_unit,*) "---------------End computing the action of the total Hamiltonian over the 1p1D system---------------"
    ELSE IF (Debug_local .AND. Nb_C > 10) THEN
      WRITE(out_unit,*)
      CALL Write_Mat(TotH_psi(1:10,1:10), out_unit, 10, info="TotH_psi(10:10sliced)")
      WRITE(out_unit,*) "---------------End computing the action of the total Hamiltonian over the 1p1D system---------------"
    END IF

  END SUBROUTINE MolecCav_Action_Cavity_1p1D

  
  SUBROUTINE MolecCav_Action_total_hamiltonian_1p1D_old(TotH_psi, MatH_psi, &
                                                      & CavH, CavPosition,  &
                                                      & Mat_dipolar_moment, Psi, Mode)
    USE QDUtil_m
    USE Cavity_mode_m
    USE Operator_1D_m
    IMPLICIT NONE

    real(kind=Rkind), intent(inout)    :: TotH_psi(:,:)                 ! already allocated !
    real(kind=Rkind), intent(in)       :: MatH_psi(:,:)     ! size Nb_M*Nb_C. |H_MatterPsi(:,i_C)> = H_Matter|Psi(:,i_C)>
    TYPE(Operator_1D_t), intent(in)    :: CavH
    TYPE(Operator_1D_t), intent(in)    :: CavPosition
    TYPE(Operator_1D_t), intent(in)    :: Mat_dipolar_moment                ! \hat{\mu}_{M}(R) = Cte.\hat{R} selon hypothèses
    real(kind=Rkind), intent(in)       :: Psi(:,:)                       ! size Nb_M*Nb_C
    !real(kind=Rkind), intent(in)       :: lambda_cavity_mode, w_cavity_mode    ! coupling strenght and eigenpulsation
    TYPE(Cavity_mode_t), intent(in)    :: Mode

    real(kind=Rkind), allocatable      :: CavHPsi(:)    ! only one line of Psi
    real(kind=Rkind), allocatable      :: Intermediary(:,:)
    real(kind=Rkind), allocatable      :: Matter_cavity_coupling_hamiltonian_1DPsi(:,:) ! only one line of Psi
    integer                            :: Nb_M, Nb_C, i_M, i_C
  
    Nb_M     = Size(Psi, 1)
    Nb_C     = Size(Psi, 2)

    TotH_psi = ZERO

    !-----------H_tot = H_Matter + H_Cavity + H_MatterCavityCoupling-----------
    !----------H_Matter|Psi> = H_Matter(Psi) already known---------
    TotH_psi = MatH_psi
  
    !-----------H_Cavity|Psi> = H_Cavity(Psi) to compute-----------
    ALLOCATE(CavHPsi(Nb_C))
    DO i_M = 1, Nb_M
      CALL Action_Operator_1D(CavHPsi, &
                                     & CavH, &
                                     & Psi(i_M,:))
      TotH_psi(i_M,:) = TotH_psi(i_M,:) + CavHPsi(:)
    END DO

    !---------------------H_MatterCavityCoupling|Psi>--------------------
    !------------------------H_CavityCoupling(Psi)-----------------------
    ALLOCATE(Intermediary(Nb_M,Nb_C))
    Intermediary = ZERO
    DO i_M = 1, Nb_M
      CALL Action_Operator_1D(Intermediary(i_M,:), &
                                     & CavPosition, &
                                     & Psi(i_M,:))
    END DO
    Intermediary = Mode%lambda*Mode%w*Intermediary

    !---------H_MatterCoupling(Intermediary) = Cte.pos_op of the matter--------
    ALLOCATE(Matter_cavity_coupling_hamiltonian_1DPsi(Nb_M,Nb_C))
    CALL MolecCav_Action_matter_dipolar_moment_1p1D(Matter_cavity_coupling_hamiltonian_1DPsi, &
                                                & Mat_dipolar_moment, &
                                                & Intermediary)                ! the matter_dipolar_moment is assumed to already contains its intensity constant (\frac{d\mu}{dq}) within its matrix /!\
                                            
    !-----------------------------H_tot = summation----------------------------
    TotH_psi = TotH_psi + Matter_cavity_coupling_hamiltonian_1DPsi

  END SUBROUTINE MolecCav_Action_total_hamiltonian_1p1D_old
  

  SUBROUTINE MolecCav_Action_matter_dipolar_moment_1p1D(Mat_dipolar_moment_psi, Mat_dipolar_moment, Psi) ! /!\ only for the app i.e. for a diatomic molecule with the approximation of a linear dip. mom. with the matter position /!\
    USE QDUtil_m
    USE Operator_1D_m
    IMPLICIT NONE

    real(kind=Rkind),       intent(inout) :: Mat_dipolar_moment_psi(:,:)
    TYPE(Operator_1D_t),    intent(in)    :: Mat_dipolar_moment             ! supposed to contains the constant (\frac{d\mu}{dq}) within its matrix
    real(kind=Rkind),       intent(in)    :: Psi(:,:)

    integer                               :: Nb_C, i_C 

    Nb_C = Size(Mat_dipolar_moment_psi, 2)

    IF (Size(Psi, 1) /= Size(Mat_dipolar_moment_psi, 1) .OR. Size(Psi, 2) /= Nb_C) THEN
      STOP "The size of Psi does not match the size of Mat_dipolar_moment_psi"
    END IF

    DO i_C = 1, Nb_C
      CALL Action_Operator_1D(Mat_dipolar_moment_psi(:,i_C), &                    ! just a matricial product vector by vector to abide by the arguments constraints of the subroutine (i think we could replace by matmul(matter_dip_mom), Psi_argu)
                                     & Mat_dipolar_moment, &
                                     & Psi(:,i_C))
    END DO

  END SUBROUTINE MolecCav_Action_matter_dipolar_moment_1p1D


  SUBROUTINE MolecCav_Construct_total_hamiltonian_1p1D(H_tot, Nb_M, Nb_C, MatH_psi, &
                                    & CavH, CavPosition, &
                                    & Mat_dipolar_moment, &
                                    & Mode)
    USE QDUtil_m
    USE Cavity_mode_m
    USE Operator_1D_m
    IMPLICIT NONE

    real(kind=Rkind),    intent(inout) :: H_tot(:,:)
    integer,             intent(in)    :: Nb_M, Nb_C
    real(kind=Rkind),    intent(in)    :: MatH_psi(:,:)     ! size Nb_M*Nb_C. |H_MatterPsi(:,i_C)> = H_Matter|Psi(:,i_C)>
    TYPE(Operator_1D_t), intent(in)    :: CavH
    TYPE(Operator_1D_t), intent(in)    :: CavPosition
    TYPE(Operator_1D_t), intent(in)    :: Mat_dipolar_moment                ! \hat{\mu}_{M}(R) = Cte.\hat{R} selon hypothèses
    !real(kind=Rkind),    intent(in)    :: lambda_cavity_mode, w_cavity_mode    ! coupling strenght and eigenpulsation
    TYPE(Cavity_mode_t), intent(in)    :: Mode


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

        CALL Action_total_hamiltonian_1p1D(Psi_result, MatH_psi, &
                                                  & CavH, CavPosition, &
                                                  & Mat_dipolar_moment, &
                                                  & Psi_basis, Mode)
          
        DO i_M = 1, Nb_M
          DO i_C = 1, Nb_C
            I = I + 1
            H_tot(I,J) = Psi_result(i_M, i_C)
          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE MolecCav_Construct_total_hamiltonian_1p1D


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

  END SUBROUTINE MolecCav_Mapping_WF_2DTO1D


  SUBROUTINE MolecCav_Average_value_H_tot(Value, H_tot, Psi)   ! /!\ FOR NOW EVERYTHING IS REAL /!\ compute the resulting vector Psi_result(:) from the action of the operator of the cavity mode on the photon state vector Psi(:) written in the Eigenbasis of H_ho
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    IMPLICIT NONE
  
    real(kind=Rkind), intent(inout)    :: Value
    real(kind=Rkind), intent(in)       :: H_tot(:,:)    
    real(kind=Rkind), intent(in)       :: Psi(:)
    
    Value = DOT_PRODUCT(Psi, MATMUL(H_tot, Psi)) 
  
  END SUBROUTINE MolecCav_Average_value_H_tot
  

END MODULE
  