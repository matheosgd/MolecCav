!==================================================================================================
!==================================================================================================
! This file is part of MolecCav.
!
!==================================================================================================
! MIT License
!
! Copyright (c) 2025 Mathéo Segaud
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.
!==================================================================================================
! README :
! to be written soon
!==================================================================================================
!==================================================================================================
MODULE Total_hamiltonian_m
  !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
  USE QDUtil_m
  IMPLICIT NONE

  
  PRIVATE

  PUBLIC Action_total_hamiltonian_1p1D, Construct_total_hamiltonian_1p1D, Average_value_H_tot, Initialize_transition_matrix, Comp&
        &ute_transition_matrix, Construct_total_hamiltonian_1p1D_R1

  INTERFACE Action_total_hamiltonian_1p1D
    MODULE PROCEDURE MolecCav_Action_total_hamiltonian_1p1D, MolecCav_Action_total_hamiltonian_1p1D_R1
  END INTERFACE
  INTERFACE Action_matter_1p1D
    MODULE PROCEDURE MolecCav_Action_matter_1p1D, MolecCav_Action_matter_1p1D_R1
  END INTERFACE
  INTERFACE Action_cavity_1p1D
    MODULE PROCEDURE MolecCav_Action_cavity_1p1D, MolecCav_Action_cavity_1p1D_R1
  END INTERFACE
  INTERFACE Construct_total_hamiltonian_1p1D
    MODULE PROCEDURE MolecCav_Construct_total_hamiltonian_1p1D
  END INTERFACE
  INTERFACE Construct_total_hamiltonian_1p1D_R1
    MODULE PROCEDURE MolecCav_Construct_total_hamiltonian_1p1D_R1
  END INTERFACE
  INTERFACE Average_value_H_tot
    MODULE PROCEDURE MolecCav_Average_value_H_tot
  END INTERFACE
  INTERFACE Transition_intensity
    MODULE PROCEDURE MolecCav_Transition_intensity
  END INTERFACE
  INTERFACE Initialize_transition_matrix
    MODULE PROCEDURE MolecCav_Initialize_transition_matrix
  END INTERFACE
  INTERFACE Compute_transition_matrix
    MODULE PROCEDURE MolecCav_Compute_transition_matrix
  END INTERFACE
    

  CONTAINS
  
    
  SUBROUTINE MolecCav_Action_total_hamiltonian_1p1D(TotH_psi, CavPosition, CavH, MatDipMomt, MatH, Psi, Debug)
    USE QDUtil_m
    USE Cavity_mode_m
    USE Operator_1D_m
    IMPLICIT NONE

    real(kind=Rkind),    intent(inout) :: TotH_psi(:,:)                        ! already allocated !
    TYPE(Operator_1D_t), intent(in)    :: CavPosition                          ! position operator associated to the cavity mode
    TYPE(Operator_1D_t), intent(in)    :: CavH                                 ! Hamiltonian associated to the cavity mode
    TYPE(Operator_1D_t), intent(in)    :: MatDipMomt                           ! dipolar moment operator associated to the matter mode. \hat{\mu}_{M}(R) = Cte.\hat{R} selon hypothèses
    TYPE(Operator_1D_t), intent(in)    :: MatH                                 ! Hamiltonian associated to the matter mode
    real(kind=Rkind),    intent(in)    :: Psi(:,:)                             ! size Nb_M*Nb_C
    logical, optional,   intent(in)    :: Debug

    real(kind=Rkind), allocatable      :: Psi_1(:,:), Psi_2(:,:), Psi_3(:,:)   ! cf. below for meaning
    integer                            :: Nb_M, Nb_C
    logical                            :: Debug_local = .TRUE.
  
    ! TotH_psi = [MatH x CavIdentity]_psi + [MatIdentity x CavH]_psi + lambda_C.w_C.[MatDipMomt x CavPosition]_psi
    !          = [CavIdentity]_psi_1 + [CavH]_psi_2 + lambda_C.w_C.[CavPosition]_psi_3
    !    Psi_1 = [MatH]_psi
    !    Psi_2 = [MatIdentity]_psi
    !    Psi_3 = [MatDipMomt]_psi

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
      STOP "### The dimensions of the operand Psi's matrix does not match the resulting TotH_psi matrix's size. Please&
                       & check initialization."
    END IF 

    IF (Nb_M /= MatH%Nb) THEN
      WRITE(out_unit,*) "The dimensions of the operand Psi's matrix does not match the dimensions of the Operator &
                       &MatH's matrix representation. Please check initialization."
      STOP "### The dimensions of the operand Psi's matrix does not match the dimensions of the Operator MatH's matrix&
                       & representation. Please check initialization."
    END IF 

    IF (Nb_M /= MatDipMomt%Nb) THEN
      WRITE(out_unit,*) "The dimensions of the operand Psi's matrix does not match the dimensions of the Operator &
                       &MatDipMomt's matrix representation. Please check initialization."
      STOP "### The dimensions of the operand Psi's matrix does not match the dimensions of the Operator Mat_dipolar_m&
                       &oment's matrix representation. Please check initialization."
    END IF 

    IF (Nb_C /= CavH%Nb) THEN
      WRITE(out_unit,*) "The dimensions of the operand Psi's matrix does not match the dimensions of the Operator &
                       &CavH's matrix representation. Please check initialization."
      STOP "### The dimensions of the operand Psi's matrix does not match the dimensions of the Operator CavH's matrix&
                       & representation. Please check initialization."
    END IF 

    IF (Nb_C /= CavPosition%Nb) THEN
      WRITE(out_unit,*) "The dimensions of the operand Psi's matrix does not match the dimensions of the Operator &
                       &CavPosition's matrix representation. Please check initialization."
      STOP "### The dimensions of the operand Psi's matrix does not match the dimensions of the Operator CavPosition's&
                       & matrix representation. Please check initialization."
    END IF 

    !--------------------------------Computation-------------------------------
    CALL Action_matter_1p1D(Psi_1, Psi_2, Psi_3, MatDipMomt, MatH, Psi)

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

    CALL Action_cavity_1p1D(TotH_psi, CavPosition, CavH, Psi_1, Psi_2, Psi_3, Debug_local)

  END SUBROUTINE MolecCav_Action_total_hamiltonian_1p1D
  

  SUBROUTINE MolecCav_Action_matter_1p1D(Psi_1, Psi_2, Psi_3, MatDipMomt, MatH, Psi)
    USE QDUtil_m
    USE Cavity_mode_m
    USE Operator_1D_m
    IMPLICIT NONE

    real(kind=Rkind),    intent(inout) :: Psi_1(:,:)                           ! cf. below for meaning
    real(kind=Rkind),    intent(inout) :: Psi_2(:,:)                           ! cf. below for meaning
    real(kind=Rkind),    intent(inout) :: Psi_3(:,:)                           ! cf. below for meaning
    TYPE(Operator_1D_t), intent(in)    :: MatDipMomt                   ! dipolar moment operator associated to the matter mode. \hat{\mu}_{M}(R) = Cte.\hat{R} selon hypothèses
    TYPE(Operator_1D_t), intent(in)    :: MatH                                 ! Hamiltonian associated to the matter mode
    real(kind=Rkind),    intent(in)    :: Psi(:,:)                             ! size Nb_M*Nb_C

    integer                            :: i_C
  
    ! TotH_psi = [MatH x CavIdentity]_psi + [MatIdentity x CavH]_psi + lambda_C.w_C.[MatDipMomt x CavPosition]_psi
    !          = [CavIdentity]_psi_1 + [CavH]_psi_2 + lambda_C.w_C.[CavPosition]_psi_3
    !    Psi_1 = [MatH]_psi
    !    Psi_2 = [MatIdentity]_psi
    !    Psi_3 = [MatDipMomt]_psi

    DO i_C = 1, Size(Psi, 2)                                                   ! initialize Matter_hamiltonianSystem_WF by applying the matter hamiltonian to each column of the matrix of the total system WF Psi
      CALL Action_Operator_1D(Psi_1(:,i_C), MatH, Psi(:,i_C))
    END DO
    
    Psi_2(:,:) = Psi(:,:)

    DO i_C = 1, Size(Psi, 2)
      CALL Action_Operator_1D(Psi_3(:,i_C), MatDipMomt, Psi(:,i_C))
    END DO

  END SUBROUTINE MolecCav_Action_matter_1p1D

  
  SUBROUTINE MolecCav_Action_cavity_1p1D(TotH_psi, CavPosition, CavH, Psi_1, Psi_2, Psi_3, Debug_local)
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
  
    ! TotH_psi = [MatH x CavIdentity]_psi + [MatIdentity x CavH]_psi + lambda_C.w_C.[MatDipMomt x CavPosition]_psi
    !          = [CavIdentity]_psi_1 + [CavH]_psi_2 + lambda_C.w_C.[CavPosition]_psi_3
    !    Psi_1 = [MatH]_psi
    !    Psi_2 = [MatIdentity]_psi
    !    Psi_3 = [MatDipMomt]_psi

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
      CALL Write_Mat(Intermediary, out_unit, Nb_C, info="[MatDipMomtxCavPosition]_psi")
    ELSE IF (Debug_local .AND. Nb_C > 10) THEN
      WRITE(out_unit,*)
      CALL Write_Mat(Intermediary(1:10,1:10), out_unit, 10, info="[MatDipMomtxCavPosition]_psi(10:10sliced)")
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

  END SUBROUTINE MolecCav_Action_cavity_1p1D


  SUBROUTINE MolecCav_Action_total_hamiltonian_1p1D_R1(TotH_psi, CavPosition, CavH, MatDipMomt, MatH, Psi, Debug)
    USE QDUtil_m
    USE Cavity_mode_m
    USE Operator_1D_m
    IMPLICIT NONE

    real(kind=Rkind),    intent(inout) :: TotH_psi(:)                          ! already allocated !
    TYPE(Operator_1D_t), intent(in)    :: CavPosition                          ! position operator associated to the cavity mode
    TYPE(Operator_1D_t), intent(in)    :: CavH                                 ! Hamiltonian associated to the cavity mode
    TYPE(Operator_1D_t), intent(in)    :: MatDipMomt                           ! dipolar moment operator associated to the matter mode. \hat{\mu}_{M}(R) = Cte.\hat{R} selon hypothèses
    TYPE(Operator_1D_t), intent(in)    :: MatH                                 ! Hamiltonian associated to the matter mode
    real(kind=Rkind),    intent(in)    :: Psi(:)                               ! size Nb_M*Nb_C
    logical, optional,   intent(in)    :: Debug

    real(kind=Rkind), allocatable      :: Psi_1(:), Psi_2(:), Psi_3(:)         ! cf. below for meaning
    integer                            :: Nb_M, Nb_C, NB
    logical                            :: Debug_local = .TRUE.
  
    ! TotH_psi = [MatH x CavIdentity]_psi + [MatIdentity x CavH]_psi + lambda_C.w_C.[MatDipMomt x CavPosition]_psi
    !          = [CavIdentity]_psi_1 + [CavH]_psi_2 + lambda_C.w_C.[CavPosition]_psi_3
    !    Psi_1 = [MatH]_psi
    !    Psi_2 = [MatIdentity]_psi
    !    Psi_3 = [MatDipMomt]_psi

    !------------------------------Initialization------------------------------
    IF (PRESENT(Debug)) Debug_local = Debug
    Nb_M     = MatH%Nb
    Nb_C     = CavH%Nb
    NB       = Size(Psi, dim=1)

    ALLOCATE(Psi_1(NB))
    ALLOCATE(Psi_2(NB))
    ALLOCATE(Psi_3(NB))
    Psi_1    = ZERO
    Psi_2    = ZERO
    Psi_3    = ZERO
    TotH_psi = ZERO

    !----------------------------Checking dimensions---------------------------
    IF (NB /= Size(TotH_psi,dim=1)) THEN
      WRITE(out_unit,*) "The dimension of the operand Psi's matrix does not match the resulting TotH_psi matrix's& 
                       & size. Please check initialization."
      STOP "### The dimension of the operand Psi's matrix does not match the resulting TotH_psi matrix's size. Please&
                       & check initialization."
    END IF 

    IF (Nb_M*Nb_C /= NB) THEN
      WRITE(out_unit,*) "The dimensions of the operand Psi's matrix and resulting TotH_psi's one do not match the dimensions of t&
                        &he two subsystem's Hamiltonian's matrix representation. Please check initialization."
      STOP "### The dimensions of the operand Psi's matrix and resulting TotH_psi's one do not match the dimensions of the two su&
                        &bsystem's Hamiltonian's matrix representation. Please check initialization."
    END IF 

    IF (Nb_M /= MatDipMomt%Nb) THEN
      WRITE(out_unit,*) "The dimensions of the operand Psi's matrix does not match the dimensions of the Operator &
                       &MatDipMomt's matrix representation. Please check initialization."
      STOP "### The dimensions of the operand Psi's matrix does not match the dimensions of the Operator Mat_dipolar_m&
                       &oment's matrix representation. Please check initialization."
    END IF 

    IF (Nb_C /= CavPosition%Nb) THEN
      WRITE(out_unit,*) "The dimensions of the operand Psi's matrix does not match the dimensions of the Operator &
                       &CavPosition's matrix representation. Please check initialization."
      STOP "### The dimensions of the operand Psi's matrix does not match the dimensions of the Operator CavPosition's&
                       & matrix representation. Please check initialization."
    END IF 

    !--------------------------------Computation-------------------------------
    CALL Action_matter_1p1D(Psi_1, Psi_2, Psi_3, MatDipMomt, MatH, Psi, Debug_local)

    IF (Debug_local .AND. Nb_C <= 50) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "-----------------Computing the action of the total Hamiltonian over the 1p1D system-----------------"
      CALL Write_Vec(Psi,   out_unit, 1, info="Psi")
      WRITE(out_unit,*)
      CALL Write_Vec(Psi_1, out_unit, 1, info="Psi_1")
      WRITE(out_unit,*)
      CALL Write_Vec(Psi_2, out_unit, 1, info="Psi_2")
      WRITE(out_unit,*)
      CALL Write_Vec(Psi_3, out_unit, 1, info="Psi_3")
    ELSE IF (Debug_local .AND. Nb_C > 10) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "-----------------Computing the action of the total Hamiltonian over the 1p1D system-----------------"
      CALL Write_Vec(Psi(1:10),   out_unit, 1, info="Psi(10sliced)")
      WRITE(out_unit,*)
      CALL Write_Vec(Psi_1(1:10), out_unit, 1, info="Psi_1(10sliced)")
      WRITE(out_unit,*)
      CALL Write_Vec(Psi_2(1:10), out_unit, 1, info="Psi_2(10sliced)")
      WRITE(out_unit,*)
      CALL Write_Vec(Psi_3(1:10), out_unit, 1, info="Psi_3(10sliced)")
    END IF

    CALL Action_cavity_1p1D(TotH_psi, CavPosition, CavH, Psi_1, Psi_2, Psi_3, Debug_local)

  END SUBROUTINE MolecCav_Action_total_hamiltonian_1p1D_R1
  

  SUBROUTINE MolecCav_Action_matter_1p1D_R1(Psi_1, Psi_2, Psi_3, MatDipMomt, MatH, Psi, Debug_local)
    USE QDUtil_m
    USE Cavity_mode_m
    USE Operator_1D_m
    USE Mapping_m
    IMPLICIT NONE

    real(kind=Rkind),    intent(inout) :: Psi_1(:)                           ! cf. below for meaning
    real(kind=Rkind),    intent(inout) :: Psi_2(:)                           ! cf. below for meaning
    real(kind=Rkind),    intent(inout) :: Psi_3(:)                           ! cf. below for meaning
    TYPE(Operator_1D_t), intent(in)    :: MatDipMomt                         ! dipolar moment operator associated to the matter mode. \hat{\mu}_{M}(R) = Cte.\hat{R} selon hypothèses
    TYPE(Operator_1D_t), intent(in)    :: MatH                               ! Hamiltonian associated to the matter mode
    real(kind=Rkind),    intent(in)    :: Psi(:)                             ! size Nb_M*Nb_C
    logical,             intent(in)    :: Debug_local

    integer                            :: i_C, Nb_C
    real(kind=Rkind), allocatable      :: Psi_R2(:,:), Psi_1_R2(:,:), Psi_3_R2(:,:)
  
    ! TotH_psi = [MatH x CavIdentity]_psi + [MatIdentity x CavH]_psi + lambda_C.w_C.[MatDipMomt x CavPosition]_psi
    !          = [CavIdentity]_psi_1 + [CavH]_psi_2 + lambda_C.w_C.[CavPosition]_psi_3
    !    Psi_1 = [MatH]_psi
    !    Psi_2 = [MatIdentity]_psi
    !    Psi_3 = [MatDipMomt]_psi

    Nb_C = Size(Psi, dim=1)/MatH%Nb

    ALLOCATE(Psi_R2(MatH%Nb, Nb_C))
    CALL Mapping_WF_R1TOR2(Psi_R2, Psi, Debug=Debug_local)
    ALLOCATE(Psi_1_R2(MatH%Nb, Nb_C))
    ALLOCATE(Psi_3_R2(MatH%Nb, Nb_C))
    Psi_1_R2 = ZERO
    Psi_3_R2 = ZERO

    DO i_C = 1, Nb_C                                                   ! initialize Matter_hamiltonianSystem_WF by applying the matter hamiltonian to each column of the matrix of the total system WF Psi
      CALL Action_Operator_1D(Psi_1_R2(:,i_C), MatH, Psi_R2(:,i_C))
    END DO
    
    Psi_2(:) = Psi(:)

    DO i_C = 1, Nb_C
      CALL Action_Operator_1D(Psi_3_R2(:,i_C), MatDipMomt, Psi_R2(:,i_C))
    END DO

    CALL Mapping_WF_R2TOR1(Psi_1, Psi_1_R2)
    CALL Mapping_WF_R2TOR1(Psi_3, Psi_3_R2)

  END SUBROUTINE MolecCav_Action_matter_1p1D_R1

  
  SUBROUTINE MolecCav_Action_cavity_1p1D_R1(TotH_psi, CavPosition, CavH, Psi_1, Psi_2, Psi_3, Debug_local)
    USE QDUtil_m
    USE Cavity_mode_m
    USE Operator_1D_m
    USE Mapping_m
    IMPLICIT NONE

    real(kind=Rkind),    intent(inout) :: TotH_psi(:)                        ! already allocated !
    TYPE(Operator_1D_t), intent(in)    :: CavPosition                          ! position operator associated to the cavity mode
    TYPE(Operator_1D_t), intent(in)    :: CavH                                 ! Hamiltonian associated to the cavity mode
    real(kind=Rkind),    intent(in)    :: Psi_1(:)                           ! cf. below for meaning
    real(kind=Rkind),    intent(in)    :: Psi_2(:)                           ! cf. below for meaning
    real(kind=Rkind),    intent(in)    :: Psi_3(:)                           ! cf. below for meaning
    logical,             intent(in)    :: Debug_local

    real(kind=Rkind), allocatable      :: Intermediary_R1(:)
    integer                            :: Nb_C, Nb_M, i_M
    real(kind=Rkind), allocatable      :: Psi_2_R2(:,:), Psi_3_R2(:,:), Intermediary_R2(:,:)

    ! TotH_psi = [MatH x CavIdentity]_psi + [MatIdentity x CavH]_psi + lambda_C.w_C.[MatDipMomt x CavPosition]_psi
    !          = [CavIdentity]_psi_1 + [CavH]_psi_2 + lambda_C.w_C.[CavPosition]_psi_3
    !    Psi_1 = [MatH]_psi
    !    Psi_2 = [MatIdentity]_psi
    !    Psi_3 = [MatDipMomt]_psi

    Nb_C = CavH%Nb
    Nb_M = Size(Psi_1, dim=1)/Nb_C
    
    TotH_psi(:) = Psi_1(:)
    IF (Debug_local .AND. Nb_C <= 50) THEN
      WRITE(out_unit,*)
      CALL Write_Vec(TotH_psi, out_unit, 1, info="[MatHxCavIdentity]_psi")
    ELSE IF (Debug_local .AND. Nb_C > 10) THEN
      WRITE(out_unit,*)
      CALL Write_Vec(TotH_psi(1:10), out_unit, 1, info="[MatHxCavIdentity]_psi(10sliced)")
    END IF

    ALLOCATE(Intermediary_R2(Nb_M, Nb_C))
    ALLOCATE(Psi_2_R2(Nb_M, Nb_C))
    Intermediary_R2 = ZERO
    Psi_2_R2        = ZERO
    CALL Mapping_WF_R1TOR2(Psi_2_R2, Psi_2, Debug=Debug_local)
    DO i_M = 1, Nb_M
      CALL Action_Operator_1D(Intermediary_R2(i_M,:), CavH, Psi_2_R2(i_M,:))
    END DO
    IF (Debug_local .AND. Nb_C <= 50) THEN
      WRITE(out_unit,*)
      CALL Write_Mat(Intermediary_R2, out_unit, Nb_C, info="[MatIdentityxCavH]_psi")
    ELSE IF (Debug_local .AND. Nb_C > 10) THEN
      WRITE(out_unit,*)
      CALL Write_Mat(Intermediary_R2(1:10,1:10), out_unit, 10, info="[MatIdentityxCavH]_psi(10:10sliced)")
    END IF

    ALLOCATE(Intermediary_R1(Nb_M*Nb_C))
    CALL Mapping_WF_R2TOR1(Intermediary_R1, Intermediary_R2, Debug=Debug_local)
    TotH_psi(:) = TotH_psi(:) + Intermediary_R1(:) 
    IF (Debug_local .AND. Nb_C <= 50) THEN
      WRITE(out_unit,*)
      CALL Write_Vec(TotH_psi, out_unit, 1, info="[MatH+CavH]_psi")
    ELSE IF (Debug_local .AND. Nb_C > 10) THEN
      WRITE(out_unit,*)
      CALL Write_Vec(TotH_psi(1:10), out_unit, 1, info="[MatH+CavH]_psi(10:10sliced)")
    END IF

    Intermediary_R2 = ZERO
    ALLOCATE(Psi_3_R2(Nb_M, Nb_C))
    CALL Mapping_WF_R1TOR2(Psi_3_R2, Psi_3)
    DO i_M = 1, Nb_M
      CALL Action_Operator_1D(Intermediary_R2(i_M,:), CavPosition, Psi_3_R2(i_M,:))
    END DO
    IF (Debug_local .AND. Nb_C <= 50) THEN
      WRITE(out_unit,*)
      CALL Write_Mat(Intermediary_R2, out_unit, Nb_C, info="[MatDipMomtxCavPosition]_psi")
    ELSE IF (Debug_local .AND. Nb_C > 10) THEN
      WRITE(out_unit,*)
      CALL Write_Mat(Intermediary_R2(1:10,1:10), out_unit, 10, info="[MatDipMomtxCavPosition]_psi(10:10sliced)")
    END IF

    CALL Mapping_WF_R2TOR1(Intermediary_R1, Intermediary_R2, Debug=Debug_local)
    TotH_psi(:) = TotH_psi(:) + CavH%lambda * CavH%w * Intermediary_R1(:)
    IF (Debug_local .AND. Nb_C <= 50) THEN
      WRITE(out_unit,*)
      CALL Write_Vec(TotH_psi, out_unit, 1, info="TotH_psi")
      WRITE(out_unit,*) "---------------End computing the action of the total Hamiltonian over the 1p1D system---------------"
    ELSE IF (Debug_local .AND. Nb_C > 10) THEN
      WRITE(out_unit,*)
      CALL Write_Vec(TotH_psi(1:10), out_unit, 1, info="TotH_psi(10:10sliced)")
      WRITE(out_unit,*) "---------------End computing the action of the total Hamiltonian over the 1p1D system---------------"
    END IF

  END SUBROUTINE MolecCav_Action_cavity_1p1D_R1


  SUBROUTINE MolecCav_Construct_total_hamiltonian_1p1D(TotH, CavPosition, CavH, MatDipMomt, MatH, Debug)
    USE QDUtil_m
    USE Cavity_mode_m
    USE Operator_1D_m
    IMPLICIT NONE

    real(kind=Rkind),    intent(inout) :: TotH(:,:)                                                ! already allocated !
    TYPE(Operator_1D_t), intent(in)    :: CavPosition
    TYPE(Operator_1D_t), intent(in)    :: CavH
    TYPE(Operator_1D_t), intent(in)    :: MatDipMomt                                       ! \hat{\mu}_{M}(R) = Cte.\hat{R} according to hypothesis
    TYPE(Operator_1D_t), intent(in)    :: MatH
    logical, optional,   intent(in)    :: Debug
   
    real(kind=Rkind), allocatable      :: Phi(:,:), TotH_phi(:,:)
    integer                            :: Nb_M, Nb_C, i_M, i_C, j_M, j_C, NB, I, J
    logical                            :: Debug_local = .FALSE.

    IF (PRESENT(Debug)) Debug_local = Debug

    IF (ALLOCATED(MatH%Diag_val_R)) THEN
      Nb_M = Size(MatH%Diag_val_R)
    ELSE IF (ALLOCATED(MatH%Dense_val_R)) THEN
      Nb_M = Size(MatH%Dense_val_R, dim=1)
    ELSE 
      WRITE(out_unit,*) "The matter Hamiltonian does not seem to have been initialized (matrices not allocated)."
      STOP "### The matter Hamiltonian does not seem to have been initialized (matrices not allocated)."
    END IF

    IF (ALLOCATED(CavH%Diag_val_R)) THEN
      Nb_C = Size(CavH%Diag_val_R)
    ELSE IF (ALLOCATED(CavH%Dense_val_R)) THEN
      Nb_C = Size(CavH%Dense_val_R, dim=1)
    ELSE 
      WRITE(out_unit,*) "The matter Hamiltonian does not seem to have been initialized (matrices not allocated)."
      STOP "### The matter Hamiltonian does not seem to have been initialized (matrices not allocated)."
    END IF

    NB = Size(TotH, dim=1)
    IF (NB /= Size(TotH, dim=2) .OR. NB /= Nb_M*Nb_C) THEN
      WRITE(out_unit,*) "The TotH matrix seems badly initialized, please check allocation."
      WRITE(out_unit,*) "Size(TotH, dim=1) = "//TO_string(NB)//" = NB"
      WRITE(out_unit,*) "Size(TotH, dim=2) = "//TO_string(Size(TotH, dim=2))
      WRITE(out_unit,*) "Nb_M*Nb_C = "//TO_string(Nb_M)//"*"//TO_string(Nb_C)//" = "//TO_string(Nb_M*Nb_C)
      STOP "### The TotH matrix seems badly initialized, please check allocation."
    END IF

    TotH(:,:) = ZERO
    ALLOCATE(Phi(Nb_M, Nb_C))
    ALLOCATE(TotH_phi(Nb_M, Nb_C))

    J = 0
    DO j_C = 1, Nb_C
      DO j_M = 1, Nb_M
        J = J + 1

        Phi = ZERO
        Phi(j_M, j_C) = ONE

        TotH_phi = ZERO
        CALL Action_total_hamiltonian_1p1D(TotH_phi, CavPosition, CavH, MatDipMomt, MatH, Phi, Debug=Debug_local)
          
        I = 0
        DO i_C = 1, Nb_C
          DO i_M = 1, Nb_M
            I = I + 1
            TotH(I,J) = TotH_phi(i_M, i_C)
          END DO
        END DO
      END DO
    END DO

    IF (Debug_local .AND. NB <= 10) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "------------------------------The total Hamiltonian matrix constructed------------------------------"
      WRITE(out_unit,*) "w_M = "//TO_string(MatH%w)//"; w_C = "//TO_string(CavH%w)//"; lambda_C = "//TO_string(CavH%lambda)
      CALL Write_Mat(TotH, out_unit, Size(TotH, dim=2), info="TotH")
      WRITE(out_unit,*) "---------------------------------End of the total Hamiltonian matrix--------------------------------"
    ELSE IF (Debug_local .AND. NB > 10) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "------------------------------The total Hamiltonian matrix constructed------------------------------"
      WRITE(out_unit,*) "w_M = "//TO_string(MatH%w)//"; w_C = "//TO_string(CavH%w)//"; lambda_C = "//TO_string(CavH%lambda)
      CALL Write_Mat(TotH(1:10,1:10), out_unit, 10, info="TotH(10:10sliced)")
      WRITE(out_unit,*) "---------------------------------End of the total Hamiltonian matrix--------------------------------"
    END IF

  END SUBROUTINE MolecCav_Construct_total_hamiltonian_1p1D


  SUBROUTINE MolecCav_Construct_total_hamiltonian_1p1D_R1(TotH, CavPosition, CavH, MatDipMomt, MatH, Debug)
    USE QDUtil_m
    USE Mapping_m
    USE Cavity_mode_m
    USE Operator_1D_m
    IMPLICIT NONE

    real(kind=Rkind),    intent(inout) :: TotH(:,:)                                                ! already allocated !
    TYPE(Operator_1D_t), intent(in)    :: CavPosition
    TYPE(Operator_1D_t), intent(in)    :: CavH
    TYPE(Operator_1D_t), intent(in)    :: MatDipMomt                                               ! \hat{\mu}_{M}(R) = Cte.\hat{R} according to hypothesis
    TYPE(Operator_1D_t), intent(in)    :: MatH
    logical, optional,   intent(in)    :: Debug
   
    real(kind=Rkind), allocatable      :: Phi_R1(:), TotH_phi_R1(:)
    real(kind=Rkind), allocatable      :: TotH_phi_R2(:,:)
    integer                            :: Nb_M, Nb_C, i_M, i_C, j_M, j_C, NB, I, J
    logical                            :: Debug_local = .FALSE.

    IF (PRESENT(Debug)) Debug_local = Debug

    IF (ALLOCATED(MatH%Diag_val_R)) THEN
      Nb_M = Size(MatH%Diag_val_R)
    ELSE IF (ALLOCATED(MatH%Dense_val_R)) THEN
      Nb_M = Size(MatH%Dense_val_R, dim=1)
    ELSE 
      WRITE(out_unit,*) "The matter Hamiltonian does not seem to have been initialized (matrices not allocated)."
      STOP "### The matter Hamiltonian does not seem to have been initialized (matrices not allocated)."
    END IF

    IF (ALLOCATED(CavH%Diag_val_R)) THEN
      Nb_C = Size(CavH%Diag_val_R)
    ELSE IF (ALLOCATED(CavH%Dense_val_R)) THEN
      Nb_C = Size(CavH%Dense_val_R, dim=1)
    ELSE 
      WRITE(out_unit,*) "The matter Hamiltonian does not seem to have been initialized (matrices not allocated)."
      STOP "### The matter Hamiltonian does not seem to have been initialized (matrices not allocated)."
    END IF

    IF (.NOT. ALLOCATED(MatDipMomt%Band_val_R) .AND. .NOT. ALLOCATED(MatDipMomt%Dense_val_R)) THEN
      WRITE(out_unit,*) "The matter dipole moment operator does not seem to have been initialized (matrices not allocated)."
      STOP "### The matter dipole moment operator does not seem to have been initialized (matrices not allocated)."
    END IF

    IF (.NOT. ALLOCATED(CavPosition%Band_val_R) .AND. .NOT. ALLOCATED(CavPosition%Dense_val_R)) THEN
      WRITE(out_unit,*) "The cavity position operator does not seem to have been initialized (matrices not allocated)."
      STOP "### The cavity position operator does not seem to have been initialized (matrices not allocated)."
    END IF

    NB = Size(TotH, dim=1)
    IF (NB /= Size(TotH, dim=2) .OR. NB /= Nb_M*Nb_C) THEN
      WRITE(out_unit,*) "The TotH matrix seems badly initialized, please check allocation."
      WRITE(out_unit,*) "Size(TotH, dim=1) = "//TO_string(NB)//" = NB"
      WRITE(out_unit,*) "Size(TotH, dim=2) = "//TO_string(Size(TotH, dim=2))
      WRITE(out_unit,*) "Nb_M*Nb_C = "//TO_string(Nb_M)//"*"//TO_string(Nb_C)//" = "//TO_string(Nb_M*Nb_C)
      STOP "### The TotH matrix seems badly initialized, please check allocation."
    END IF

    TotH(:,:) = ZERO
    ALLOCATE(Phi_R1(NB))
    ALLOCATE(TotH_phi_R1(NB))
    ALLOCATE(TotH_phi_R2(Nb_M, Nb_C))

    DO J = 1, NB
      Phi_R1 = ZERO
      Phi_R1(J) = ONE

      TotH_phi_R1 = ZERO
      CALL Action_total_hamiltonian_1p1D(TotH_phi_R1, CavPosition, CavH, MatDipMomt, MatH, Phi_R1, Debug=.FALSE.)
!      CALL Mapping_WF_R1TOR2(TotH_phi_R2, TotH_phi_R1)  
      
!      I = 0
!      DO i_C = 1, Nb_C
!        DO i_M = 1, Nb_M
!          I = I + 1
          TotH(:,J) = TotH_phi_R1(:)
!          TotH(I,J) = TotH_phi_R2(i_M, i_C)
!        END DO
!      END DO
    END DO

    IF (Debug_local .AND. NB <= 20) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "------------------------------The total Hamiltonian matrix constructed------------------------------"
      WRITE(out_unit,*) "w_M = "//TO_string(MatH%w)//"; w_C = "//TO_string(CavH%w)//"; lambda_C = "//TO_string(CavH%lambda)
      CALL Write_Mat(TotH, out_unit, Size(TotH, dim=2), info="TotH")
      WRITE(out_unit,*) "---------------------------------End of the total Hamiltonian matrix--------------------------------"
    ELSE IF (Debug_local .AND. NB > 20) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "------------------------------The total Hamiltonian matrix constructed------------------------------"
      WRITE(out_unit,*) "w_M = "//TO_string(MatH%w)//"; w_C = "//TO_string(CavH%w)//"; lambda_C = "//TO_string(CavH%lambda)
      CALL Write_Mat(TotH(1:20,1:20), out_unit, 20, info="TotH(20:20sliced)")
      WRITE(out_unit,*) "---------------------------------End of the total Hamiltonian matrix--------------------------------"
    END IF

  END SUBROUTINE MolecCav_Construct_total_hamiltonian_1p1D_R1


  SUBROUTINE MolecCav_Average_value_H_tot(Value, H_tot, Psi)   ! /!\ FOR NOW EVERYTHING IS REAL /!\ compute the resulting vector Psi_result(:) from the action of the operator of the cavity mode on the photon state vector Psi(:) written in the Eigenbasis of H_ho
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    IMPLICIT NONE
  
    real(kind=Rkind), intent(inout)    :: Value
    real(kind=Rkind), intent(in)       :: H_tot(:,:)    
    real(kind=Rkind), intent(in)       :: Psi(:)
    
    Value = DOT_PRODUCT(Psi, MATMUL(H_tot, Psi)) 
  
  END SUBROUTINE MolecCav_Average_value_H_tot
  

  SUBROUTINE MolecCav_Transition_intensity(Intensity, InitPsi, MatDipMomt, FinPsi, Rank_size1, Rank_size2)   ! /!\ FOR NOW DESIGNED FOR 1p1D
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    USE Algebra_m
    USE Mapping_m
    USE Operator_1D_m
    IMPLICIT NONE
  
    real(kind=Rkind),    intent(inout) :: Intensity
    real(kind=Rkind),    intent(in)    :: InitPsi(:)                                                                                ! 1p1D so far !
    TYPE(Operator_1D_t), intent(in)    :: MatDipMomt
    real(kind=Rkind),    intent(in)    :: FinPsi(:)

    integer,             intent(in)    :: Rank_size1
    integer,             intent(in)    :: Rank_size2

    real(kind=Rkind), allocatable      :: InitPsi_1p1D(:,:)
    real(kind=Rkind), allocatable      :: FinPsi_1p1D(:,:)
    real(kind=Rkind), allocatable      :: Intermediary(:,:)  
    integer                            :: Nb_M, Nb_C, i_C
    logical, parameter                 :: Debug = .FALSE.

    Nb_M = Rank_size1
    Nb_C = Rank_size2

    IF (Nb_M*Nb_C /= Size(InitPsi, dim=1)) THEN
      WRITE(out_unit,*) "Unconsistent arguments of Transition intensity : InitPsi dimensions do not match the provided Rank_sizes"
      WRITE(out_unit,*) "Rank_size1 = "//TO_string(Nb_M)//";  Rank_size2 = "//TO_string(Nb_C)//"; Size(InitPsi, dim=1) = "//TO_st&
                        &ring(Size(InitPsi, dim=1))
      WRITE(out_unit,*) "Please check arguments"
      STOP "### Unconsistent arguments of Transition_intensity"
    END IF
    IF (Nb_M*Nb_C /= Size(FinPsi, dim=1)) THEN
      WRITE(out_unit,*) "Unconsistent arguments of Transition intensity : InitPsi dimensions do not match the provided Rank_sizes"
      WRITE(out_unit,*) "Rank_size1 = "//TO_string(Nb_M)//";  Rank_size2 = "//TO_string(Nb_C)//"; Size(InitPsi, dim=1) = "//TO_st&
                        &ring(Size(FinPsi, dim=1))
      WRITE(out_unit,*) "Please check arguments"
      STOP "### Unconsistent arguments of Transition_intensity"
    END IF
    IF (Nb_M /= MatDipMomt%Nb) THEN
      WRITE(out_unit,*) "Unconsistent arguments of Transition intensity : MatDipMomt sizes do not match the provided Rank_sizes1"
      WRITE(out_unit,*) "Rank_size1 = "//TO_string(Nb_M)//";  Rank_size2 = "//TO_string(MatDipMomt%Nb)
      WRITE(out_unit,*) "Please check arguments"
      STOP "### Unconsistent arguments of Transition_intensity"
    END IF

    ALLOCATE(InitPsi_1p1D(Nb_M, Nb_C))
    ALLOCATE(FinPsi_1p1D(Nb_M, Nb_C))
    ALLOCATE(Intermediary(Nb_M, Nb_C))

    CALL Mapping_WF_R1TOR2(FinPsi_1p1D,  FinPsi,  Debug=Debug)
    CALL Mapping_WF_R1TOR2(InitPsi_1p1D, InitPsi, Debug=Debug)

    IF (Debug) CALL Write_operator_1D(MatDipMomt)

    DO i_C = 1, Nb_C
      CALL Action_Operator_1D(Intermediary(:,i_C), MatDipMomt, FinPsi_1p1D(:,i_C))
    END DO
    IF (Debug) CALL Write_Mat(Intermediary, out_unit, Size(Intermediary, dim=2), info="\hat{\mu}_{mat}\ket{\Psi_f}")

    CALL Scalar_product(Intensity, InitPsi_1p1D, Intermediary)
    IF (Debug) WRITE(out_unit,*) "\bra{\Psi_i}\hat{\mu}_{mat}\ket{\Psi_f}"//TO_string(Intensity)

    Intensity = ABS(Intensity)**2
    IF (Debug) WRITE(out_unit,*) "|\bra{\Psi_i}\hat{\mu}_{mat}\ket{\Psi_f}|**2"//TO_string(Intensity)
  
  END SUBROUTINE MolecCav_Transition_intensity
  

  SUBROUTINE MolecCav_Initialize_transition_matrix(Intensities, MatDipMomt, REigvec, Energy_threshold, REigval, Nb_states, Debug)   ! /!\ FOR NOW DESIGNED FOR 1p1D
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    USE Operator_1D_m
    IMPLICIT NONE
  
    real(kind=Rkind), allocatable, intent(inout) :: Intensities(:,:)
    TYPE(Operator_1D_t),           intent(in)    :: MatDipMomt
    real(kind=Rkind),              intent(in)    :: REigvec(:,:)
    real(kind=Rkind), optional,    intent(in)    :: REigval(:)
    real(kind=Rkind), optional,    intent(in)    :: Energy_threshold
    integer,          optional,    intent(in)    :: Nb_states
    logical,          optional,    intent(in)    :: Debug

    integer                                      :: N, I
    logical                                      :: Debug_local = .FALSE.

    IF (PRESENT(Debug)) Debug_local = Debug

    IF (MOD(Size(REigvec, dim=1),MatDipMomt%Nb) /= 0) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "Unconsistent arguments of Transition intensity : NB /= Nb_M*Integer"
      WRITE(out_unit,*) "NB = "//TO_string(Size(REigvec, dim=1))//";  Nb_M = "//TO_string(MatDipMomt%Nb)
      WRITE(out_unit,*) "Please check arguments"
      STOP "### Unconsistent arguments of Transition_intensity_matrix"
    END IF

    IF (PRESENT(Energy_threshold) .AND. .NOT. PRESENT(REigval)) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "The list of the total Hamiltonian Eigenenergies (REigval) is expected when the selection criterion is en&
                        &ergy-based (Energy_threshold provided). Please check the arguments."
      STOP "### Missing REigval argument in Transition_intensity_matrix"

    ELSE IF (PRESENT(Energy_threshold)) THEN
      I = 0 ! /!\ REMPLACER PAR COUNT(Reigval(i+1) > Reigval(1) + Energy_threshold) /!\
      DO 
        I = I + 1
        IF (REigval(I+1) - REigval(1) > Energy_threshold) EXIT
      END DO
      N = I

    ELSE IF (PRESENT(Nb_states)) THEN
      N = Nb_states

    ELSE IF ((.NOT. PRESENT(Energy_threshold)) .AND. (.NOT. PRESENT(Nb_states))) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "########################## WARNING ########################## WARNING ########################## WARNING&
                       & #########################"
      WRITE(out_unit,*) "               No criterion provided to select the states with which the transition intensities have to &
                       &be computed"
      WRITE(out_unit,*) "                                              All Eigenstates will thus be considered"
      WRITE(out_unit,*) "########################## WARNING ########################## WARNING ########################## WARNING&
                       & #########################"
      N = Size(REigvec, dim=1)
      WRITE(out_unit,*) N
    END IF 

    IF (Debug_local) WRITE(out_unit,*)
    IF (Debug_local) WRITE(out_unit,*) TO_string(N)//" States will be taken into account to compute the transition intensities."

    ALLOCATE(Intensities(N, N))
    Intensities = ZERO
    
    IF (Debug_local) WRITE(out_unit,*)
    IF (Debug_local) CALL Write_Mat(Intensities, out_unit, Size(Intensities), info="Initialized intensity matrix")
    FLUSH(out_unit)

  END SUBROUTINE MolecCav_Initialize_transition_matrix
  

  SUBROUTINE MolecCav_Compute_transition_matrix(Intensities, MatDipMomt, REigvec, Debug)   ! /!\ FOR NOW DESIGNED FOR 1p1D
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    USE Operator_1D_m
    IMPLICIT NONE
  
    real(kind=Rkind),           intent(inout) :: Intensities(:,:)
    TYPE(Operator_1D_t),        intent(in)    :: MatDipMomt
    real(kind=Rkind),           intent(in)    :: REigvec(:,:)
    logical,          optional, intent(in)    :: Debug

    real(kind=Rkind), allocatable      :: InitPsi_1p1D(:,:)
    real(kind=Rkind), allocatable      :: FinPsi_1p1D(:,:)
    real(kind=Rkind), allocatable      :: Intermediary(:,:)  
    integer                            :: Nb_M, Nb_C, i_C, N, I, J
    logical                            :: Debug_local = .FALSE.

    IF (PRESENT(Debug)) Debug_local = Debug

    IF (MOD(Size(REigvec, dim=1),MatDipMomt%Nb) /= 0) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "Unconsistent arguments of Transition intensity : NB /= Nb_M*Integer"
      WRITE(out_unit,*) "NB = "//TO_string(Size(REigvec, dim=1))//";  Nb_M = "//TO_string(MatDipMomt%Nb)
      WRITE(out_unit,*) "Please check arguments"
      STOP "### Unconsistent arguments of Transition_intensity_matrix"
    END IF

    N    = Size(Intensities, dim=1)
    Nb_M = MatDipMomt%Nb             ! will be ND_indexes soon
    Nb_C = Size(REigvec, dim=1)/Nb_M

    IF (Debug_local) WRITE(out_unit,*)
    DO I = 1, N
      DO J =  1, N
        CALL Transition_intensity(Intensities(I,J), REigvec(:,I), MatDipMomt, REigvec(:,J), Nb_M, Nb_C)
        IF (Debug_local) WRITE(out_unit,*) "Transition \overrightarrow{VP}_"//TO_string(I)//" --> \overrightarrow{VP}_"//TO_strin&
                                           &g(J)//" = "//TO_string(Intensities(I,J))
      END DO
    END DO 
    
    IF (Debug_local) WRITE(out_unit,*)
    IF (Debug_local) CALL Write_Mat(Intensities, out_unit, Size(Intensities), info="Intensities matrix")
    
  END SUBROUTINE MolecCav_Compute_transition_matrix
  

END MODULE
  