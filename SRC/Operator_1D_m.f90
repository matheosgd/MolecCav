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
!==================================================================================================
MODULE Operator_1D_m
  !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
  USE QDUtil_m
  USE Cavity_mode_m
  IMPLICIT NONE


  TYPE, EXTENDS(Cavity_mode_t) :: Operator_1D_t
    character(len=:),    allocatable :: Operator_type                          ! ex : "Hamiltonian", "Position", "Nb_photons", etc
    logical                          :: Dense = .FALSE.                        ! if .TRUE. then the matrix storage will not be optimized and it will be stored as a Dense matrix
    integer                          :: Upper_bandwidth = 0                    ! if type = "Band". Gives the number of additional bands to consider above the diagonal.
    integer                          :: Lower_bandwidth = 0                    ! if type = "Band". Gives the number of additional bands to consider below the diagonal. Ex : Upper_bandwidth=Lower_bandwidth=1 would give a tridiagonal matrix
    real(kind=Rkind),    allocatable :: Dense_val_R(:,:)                       ! if type = "Dense"
    real(kind=Rkind),    allocatable :: Diag_val_R(:)                          ! if type = "Diagonal"
    real(kind=Rkind),    allocatable :: Band_val_R(:,:)                        ! if subtype = "Band". The number of columns will be the number of diagonals to consider
  END TYPE


  PRIVATE

  PUBLIC Operator_1D_t, Construct_Operator_1D, Action_Operator_1D, Average_value_operator_1D, & 
       !& MolecCav_Action_operator_2D, MolecCav_Average_value_operator_2D, &
       & Write_operator_1D

  INTERFACE Construct_Operator_1D
    MODULE PROCEDURE MolecCav_Construct_Operator_1D
  END INTERFACE
  INTERFACE Action_Operator_1D
    MODULE PROCEDURE MolecCav_Action_Operator_1D
  END INTERFACE
  INTERFACE Average_value_operator_1D
    MODULE PROCEDURE MolecCav_Average_value_operator_1D
  END INTERFACE
  INTERFACE Write_operator_1D
    MODULE PROCEDURE MolecCav_Write_operator_1D
  END INTERFACE
    

  CONTAINS


  SUBROUTINE MolecCav_Construct_Operator_1D(Operator, Operator_type, Dense, Mode, Debug)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
    USE QDUtil_m
    IMPLICIT NONE

    TYPE(Operator_1D_t), intent(inout) :: Operator                             ! the object of type Operator_t to be constructed here
    character(len=*),    intent(in)    :: Operator_type                        ! ex : "Hamiltonian", "Position", etc. (len=:) Expects to be allocatable, while (len=*) is dedicated to a procedure argument.
    logical, optional,   intent(in)    :: Dense                                ! if .TRUE. then the matrix storage will not be optimized and it will be stored as a Dense matrix
    TYPE(Cavity_mode_t), intent(in)    :: Mode                                 ! the HO/Cavity mode which the operator is relative to
    logical, optional,   intent(in)    :: Debug

    logical                            :: Debug_local = .FALSE.

    !-----------------------------Debugging options----------------------------
    IF (PRESENT(Debug)) Debug_local = Debug
    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "-------------------Arguments of MolecCav_Construct_Operator_1D------------------"
      WRITE(out_unit,*) "The <<Operator>> argument :"
      CALL Write_operator_1D(Operator)
      WRITE(out_unit,*) "The <<Operator_type>> argument : ", Operator_type
      IF (present(Dense)) WRITE(out_unit,*) "Dense :", Dense
      WRITE(out_unit,*) "The <<Mode>> argument : ", Mode
      WRITE(out_unit,*) "-----------------End arguments of MolecCav_Construct_Operator_1D----------------"
      FLUSH(out_unit)
    END IF
    
    !--------------First steps of the construction of the Operator-------------
    ALLOCATE(character(len=LEN_TRIM(Operator_type)) :: Operator%Operator_type) ! /!\ strings cannot be allocated the exact same way as tables ! /!\
    Operator%Operator_type = TO_lowercase(TRIM(Operator_type))                 ! allocation on assignement (not anymore). Operator_type has the right lengths (no spaces added) thanks to len=* at declaration and it will fit the Op%op_type thanks to len=:, allocatable at declaration of the derived type. 
    Operator%Cavity_mode_t = Mode                                              ! no need to have a variable of type Cavity_mode_t to use the "%" writing !!!

    IF (PRESENT(Dense)) THEN
      Operator%Dense = Dense
    END IF

    !--------------------Construction of the matrix Operator-------------------
    SELECT CASE (Operator%Operator_type)                         ! TO_lowercase avoid case sensitivity issues
      CASE ("hamiltonian")
        CALL MolecCav_Construct_H_cavity_mode(Hamiltonian=Operator)
    
      CASE ("position")
        CALL MolecCav_Construct_x_cavity_mode(Position_Op=Operator)
      
      CASE ("nb_photons")
        CALL MolecCav_Construct_N_cavity_mode(Nb_photon_Op=Operator)

      CASE DEFAULT
        WRITE(out_unit,*) "No Operator type recognized, please verify the input of Construct_Operator_1D subroutine"
        STOP "### No Operator type recognized, please verify the input of Construct_Operator_1D subroutine"

    END SELECT

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "-------------Operator constructed by MolecCav_Construct_Operator_1D-------------"
      CALL Write_operator_1D(Operator)
      WRITE(out_unit,*) "-----------End operator constructed by MolecCav_Construct_Operator_1D-----------"
    END IF

  END SUBROUTINE MolecCav_Construct_Operator_1D


  SUBROUTINE MolecCav_Construct_H_cavity_mode(Hamiltonian)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    IMPLICIT NONE
    
    TYPE(Operator_1D_t), intent(inout) :: Hamiltonian                          ! matrix of the one-dimensional harmonic Hamiltonian associated with HO D

    integer                            :: i                                    ! loop increments /!\ Fortran counts from 1 to Nb !!! /!\

    WRITE(out_unit,*) ''
    WRITE(out_unit,*) '*******************************************************'
    WRITE(out_unit,*) '******* CONSTRUCTING THE HAMILTONIAN OF THE HO ********'

    IF (.NOT. Hamiltonian%Dense) THEN
      !---------------------Initialization to default values-------------------
      ALLOCATE(Hamiltonian%Diag_val_R(Hamiltonian%Nb))
      !------------------------Construction of the matrix----------------------
      DO i = 1, Hamiltonian%Nb                                                             ! /!\ Fortran counts from 1 to Nb !!! /!\
        Hamiltonian%Diag_val_R(i) = Hamiltonian%w*(i - ONE + HALF)                         ! "-1" because the first Fortran vector is the fundamental eigenvector of the HO i.e. the 0^{th} ket 
      END DO
      
    ELSE
      !---------------------Initialization to default values-------------------
      ALLOCATE(Hamiltonian%Dense_val_R(Hamiltonian%Nb, Hamiltonian%Nb))
      Hamiltonian%Dense_val_R = ZERO
      !------------------------Construction of the matrix----------------------
      DO i = 1, Hamiltonian%Nb                                                             ! /!\ Fortran counts from 1 to Nb !!! /!\
        Hamiltonian%Dense_val_R(i,i) = Hamiltonian%w*(i - ONE + HALF)                      ! "-1" because the first Fortran vector is the fundamental eigenvector of the HO i.e. the 0^{th} ket 
      END DO
    END IF
    
    WRITE(out_unit,*) '********* HAMILTONIAN OF THE HO CONSTRUCTED ***********'
    WRITE(out_unit,*) '*******************************************************'

  END SUBROUTINE MolecCav_Construct_H_cavity_mode


  SUBROUTINE MolecCav_Construct_x_cavity_mode(Position_Op)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    IMPLICIT NONE
    
    TYPE(Operator_1D_t), intent(inout) :: Position_Op

    integer                            :: i                                 ! loop increments /!\ Fortran counts from 1 to Nb !!! /!\

    WRITE(out_unit,*) ''
    WRITE(out_unit,*) '*******************************************************'
    WRITE(out_unit,*) '******* CONSTRUCTING THE POSITION OP OF THE HO ********'

    IF ((.NOT. Position_Op%Dense) .AND. Position_Op%Nb > 1) THEN
      !----------Initialization of the characteristics of the operator---------
      Position_Op%Upper_bandwidth   = 1
      Position_Op%Lower_bandwidth   = 1
      !---------------------Initialization to default values-------------------
      ALLOCATE(Position_Op%Band_val_R(Position_Op%Nb,3))                                   ! Nb lines (number of diagonal elements) and 3 columns because 3 bands to consider : the diagonal, and the two bands above and below it
      Position_Op%Band_val_R = ZERO
      !------------------------Construction of the matrix----------------------
      DO i = 1, Position_Op%Nb - 1                                                         ! /!\ Fortran counts from 1 to Nb !!! /!\ Nb-1 not to have Band_val_R(i+1) out of range
        Position_Op%Band_val_R(i,1)   = SQRT(REAL(i,kind=Rkind))
        Position_Op%Band_val_R(i+1,3) = SQRT(REAL(i,kind=Rkind))
      END DO
      Position_Op%Band_val_R = Position_Op%Band_val_R / SQRT(TWO * Position_Op%w * Position_Op%m)
        
    ELSE IF (.NOT. Position_Op%Dense) THEN
            !---------------------Initialization to default values-------------------
      ALLOCATE(Position_Op%Diag_val_R(Position_Op%Nb))
      !------------------------Construction of the matrix----------------------
      DO i = 1, Position_Op%Nb                                                             ! /!\ Fortran counts from 1 to Nb !!! /!\
        Position_Op%Diag_val_R(i) = ZERO                                                   ! the position operator matrix has first value (i.e. only value in the Nb = 0 case) 0 
      END DO

    ELSE 
      !---------------------Initialization to default values-------------------
      ALLOCATE(Position_Op%Dense_val_R(Position_Op%Nb, Position_Op%Nb))
      Position_Op%Dense_val_R = ZERO
      !------------------------Construction of the matrix----------------------
      DO i = 1, Position_Op%Nb - 1                                                         ! /!\ Fortran counts from 1 to Nb !!! /!\
        Position_Op%Dense_val_R(i,i+1) = SQRT(REAL(i,kind=Rkind))
        Position_Op%Dense_val_R(i+1,i) = SQRT(REAL(i,kind=Rkind))
      END DO
      Position_Op%Dense_val_R = Position_Op%Dense_val_R / SQRT(TWO * Position_Op%w * Position_Op%m)
    END IF
    
    WRITE(out_unit,*) '********* POSITION OP OF THE HO CONSTRUCTED ***********'
    WRITE(out_unit,*) '*******************************************************'

  END SUBROUTINE MolecCav_Construct_x_cavity_mode


  SUBROUTINE MolecCav_Construct_N_cavity_mode(Nb_photon_Op)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    IMPLICIT NONE
    
    TYPE(Operator_1D_t), intent(inout) :: Nb_photon_Op

    integer                            :: i                                 ! loop increments /!\ Fortran counts from 1 to Nb !!! /!\

    WRITE(out_unit,*) ''
    WRITE(out_unit,*) '*******************************************************'
    WRITE(out_unit,*) '******* CONSTRUCTING THE NB QUANTA OP OF THE HO *******'

    IF (.NOT. Nb_photon_Op%Dense) THEN
      !---------------------Initialization to default values-------------------
      ALLOCATE(Nb_photon_Op%Diag_val_R(Nb_photon_Op%Nb))
      !------------------------Construction of the matrix----------------------
      DO i = 1, Nb_photon_Op%Nb                                                             ! /!\ Fortran counts from 1 to Nb !!! /!\
        Nb_photon_Op%Diag_val_R(i) = i - 1
      END DO
      
    ELSE
      !---------------------Initialization to default values-------------------
      ALLOCATE(Nb_photon_Op%Dense_val_R(Nb_photon_Op%Nb, Nb_photon_Op%Nb))
      Nb_photon_Op%Dense_val_R = ZERO
      !------------------------Construction of the matrix----------------------
      DO i = 1, Nb_photon_Op%Nb                                                             ! /!\ Fortran counts from 1 to Nb !!! /!\
        Nb_photon_Op%Dense_val_R(i,i) = i - 1
      END DO
    END IF
    
    WRITE(out_unit,*) '********* NB QUANTA OP OF THE HO CONSTRUCTED **********'
    WRITE(out_unit,*) '*******************************************************'

  END SUBROUTINE MolecCav_Construct_N_cavity_mode

  
  SUBROUTINE MolecCav_Action_Operator_1D(Op_psi, Operator, Psi)   ! /!\ FOR NOW EVERYTHING IS REAL /!\ compute the resulting vector Psi_result(:) from the action of the operator of the cavity mode on the photon state vector Psi_argument(:) written in the Eigenbasis of H_ho
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    IMPLICIT NONE

    real(kind=Rkind),    intent(inout) :: Op_psi(:)
    TYPE(Operator_1D_t), intent(in)    :: Operator
    real(kind=Rkind),    intent(in)    :: Psi(:)

    !----------------------------Checking dimensions---------------------------
    IF (Operator%Nb /= Size(Psi)) THEN
      WRITE(out_unit,*) "The dimensions of the Operator's matrix representation does not match the operand Psi's vector&
                       & size. Please check initialization."
      STOP "### The dimensions of the Operator's matrix representation does not match the operand Psi's vector size. Please& 
                       & check initialization."
    END IF 

    IF (Operator%Nb /= Size(Op_psi)) THEN
      WRITE(out_unit,*) "The dimensions of the Operator's matrix representation does not match the resulting Op_psi vector's& 
                       & size. Please check initialization."
      STOP "### The dimensions of the Operator's matrix representation does not match the resulting Op_psi vector's size. Please& 
                       & check initialization."
    END IF 

    !--------------------Selection of the calculation method-------------------
    IF      (ALLOCATED(Operator%Diag_val_R))   THEN
      CALL MolecCav_Action_Diag_Operator_1D(Op_psi=Op_psi, Operator=Operator, Psi=Psi)

    ELSE IF (ALLOCATED(Operator%Band_val_R))   THEN
      CALL MolecCav_Action_Band_Operator_1D(Op_psi=Op_psi, Operator=Operator, Psi=Psi)

    ELSE IF (ALLOCATED(Operator%Dense_val_R)) THEN
      CALL MolecCav_Action_Dense_Operator_1D(Op_psi=Op_psi, Operator=Operator, Psi=Psi)

    ELSE
      WRITE(out_unit,*) "None of this operator's matrices are allocated. Please check its initalization."
      STOP "### None of this operator's matrices are allocated. Please check its initalization."
    
    END IF

  END SUBROUTINE MolecCav_Action_Operator_1D

  
  SUBROUTINE MolecCav_Action_Dense_Operator_1D(Op_psi, Operator, Psi)  ! _R for the case where Psi is a real vector. 
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    IMPLICIT NONE
    
    real(kind=Rkind),    intent(inout) :: Op_psi(:)
    TYPE(Operator_1D_t), intent(in)    :: Operator
    real(kind=Rkind),    intent(in)    :: Psi(:)

    Op_psi(:) = matmul(Operator%Dense_val_R, Psi)

  END SUBROUTINE MolecCav_Action_Dense_Operator_1D

  
  SUBROUTINE MolecCav_Action_Diag_Operator_1D(Op_psi, Operator, Psi) 
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    IMPLICIT NONE
    
    real(kind=Rkind),    intent(inout) :: Op_psi(:)
    TYPE(Operator_1D_t), intent(in)    :: Operator
    real(kind=Rkind),    intent(in)    :: Psi(:)

    Op_psi = Operator%Diag_val_R * Psi

  END SUBROUTINE MolecCav_Action_Diag_Operator_1D

  
  SUBROUTINE MolecCav_Action_Band_Operator_1D(Op_psi, Operator, Psi)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    IMPLICIT NONE
    
    real(kind=Rkind),    intent(inout) :: Op_psi(:)
    TYPE(Operator_1D_t), intent(in)    :: Operator
    real(kind=Rkind),    intent(in)    :: Psi(:)

    integer                            :: i, Nb

    Nb = size(Op_psi)
    IF (Nb /= Operator%Nb) THEN
      WRITE(out_unit,*) "The size of the operator's matrix Band_val_R does not match the size of the &
          & resulting vector Op_psi. Please check their initialization."
      STOP "### The size of the operator's matrix Band_val_R does not match the size of the &
          & resulting vector Op_psi. Please check their initialization."
    ELSE IF (Nb /= Size(Psi)) THEN
      WRITE(out_unit,*) "The size of the resulting vector Op_psi does not match the size of the operand & 
          & vector Psi. Please check their initialization."
      STOP "### The size of the resulting vector Op_psi does not match the size of the operand & 
          & vector Psi. Please check their initialization."
    END IF

    Op_psi     = ZERO
    Op_psi     = Operator%Band_val_R(:,2) * Psi
    Op_psi(1)  = Op_psi(1)  + Operator%Band_val_R(2,3)    * Psi(2)
    Op_psi(Nb) = Op_psi(Nb) + Operator%Band_val_R(Nb-1,1) * Psi(Nb-1)
    DO i = 2, Nb-1
      Op_psi(i) = Op_psi(i) + &
                & Operator%Band_val_R(i-1,1) * Psi(i-1) + &
                & Operator%Band_val_R(i+1,3) * Psi(i+1)
    END DO

  END SUBROUTINE MolecCav_Action_Band_Operator_1D
  

  SUBROUTINE MolecCav_Average_value_operator_1D(Value, Operator, Psi)   ! /!\ FOR NOW EVERYTHING IS REAL /!\
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    USE Algebra_m
    IMPLICIT NONE

    real(kind=Rkind),    intent(inout) :: Value
    TYPE(Operator_1D_t), intent(in)    :: Operator
    real(kind=Rkind),    intent(in)    :: Psi(:)

    real(kind=Rkind), allocatable      :: Intermediary(:)
    integer                            :: Nb

    Nb = Size(Psi)
    ALLOCATE(Intermediary(Nb))

    CALL MolecCav_Action_Operator_1D(Intermediary, Operator, Psi)
    CALL Scalar_product(Value, Psi, Intermediary)

    DEALLOCATE(Intermediary)
    
  END SUBROUTINE MolecCav_Average_value_operator_1D


  SUBROUTINE MolecCav_Write_operator_1D(Operator)
    TYPE(Operator_1D_t), intent(in) :: Operator

    WRITE(out_unit,*) "---------------------------WRITING THE OPERATOR TYPE---------------------------"
    CALL Write_cavity_mode(Operator%Cavity_mode_t)
    IF (ALLOCATED(Operator%Operator_type)) THEN
      WRITE(out_unit,*) "_______________________________________the operator parameters______________________________________"
      WRITE(out_unit,*) "|The operator's nature <<Operator_type>> information do is allocated, and is | ", Operator%Operator_type
      FLUSH(out_unit)
      WRITE(out_unit,*) "|____________________________________________________________________________|______________________"
    ELSE 
      WRITE(out_unit,*) "___________________________the operator parameters___________________________"
      WRITE(out_unit,*) "|The operator's nature <<Operator_type>> information is NOT allocated.       |"
      FLUSH(out_unit)
      WRITE(out_unit,*) "|____________________________________________________________________________|___&
                      &_________________________"
    END IF

    WRITE(out_unit,*) "|Is its matrix supposed to be represented as a dense one ? Operator%Dense    | ", Operator%Dense
    WRITE(out_unit,*) "|____________________________________________________________________________|_____&
                      &_______________________"
    WRITE(out_unit,*) "|Case of a band matrix, (Operator%Upper_bandwidth, Operator%Lower_bandwidth) | (", & 
                     & Operator%Upper_bandwidth, Operator%Lower_bandwidth, ")"
    FLUSH(out_unit)
    WRITE(out_unit,*) "|____________________________________________________________________________|_____&
                      &_______________________"
    WRITE(out_unit,*) "________________________the operator's representations_______________________"
    IF (ALLOCATED(Operator%Diag_val_R)) THEN                                   ! we assume that the code is supposed to be used only allocating one of the matrices of each Operator_t object
      WRITE(out_unit,*) "|The operator's Diagonal matrix representation has been used, and is         |"
      WRITE(out_unit,*) "|____________________________________________________________________________|"
      CALL Write_Vec(Operator%Diag_val_R, out_unit, 1, info="Diag_val_R")
      WRITE(out_unit,*) "|____________________________________________________________________________|"
      FLUSH(out_unit)
    ELSE 
      WRITE(out_unit,*) "|The operator's Diagonal matrix representation is NOT allocated.             |"
      FLUSH(out_unit)
      WRITE(out_unit,*) "|____________________________________________________________________________|"
    END IF

    IF (ALLOCATED(Operator%Band_val_R)) THEN
      WRITE(out_unit,*) "|The operator's Band matrix representation has been used, and is             |"
      WRITE(out_unit,*) "|____________________________________________________________________________|"
      CALL Write_Mat(Operator%Band_val_R, out_unit, 3, info="Band_val_R")
      WRITE(out_unit,*) "|____________________________________________________________________________|"
      FLUSH(out_unit)
    ELSE 
      WRITE(out_unit,*) "|The operator's Band matrix representation is NOT allocated.                 |"
      FLUSH(out_unit)
      WRITE(out_unit,*) "|____________________________________________________________________________|"
    END IF

    IF (ALLOCATED(Operator%Dense_val_R)) THEN
      WRITE(out_unit,*) "|____________________________________________________________________________|"
      WRITE(out_unit,*) "|The operator's Dense matrix representation has been used, and is            |"
      CALL Write_Mat(Operator%Dense_val_R, out_unit, Operator%Nb, info="Dense_val_R")
      WRITE(out_unit,*) "|____________________________________________________________________________|"
      FLUSH(out_unit)
    ELSE 
      WRITE(out_unit,*) "|The operator's Dense matrix representation is NOT allocated.                |"
      FLUSH(out_unit)
      WRITE(out_unit,*) "|____________________________________________________________________________|"
    END IF

    WRITE(out_unit,*) "-----------------------------END WRITE OPERATOR TYPE----------------------------"
    FLUSH(out_unit)

  END SUBROUTINE MolecCav_Write_operator_1D


END MODULE
