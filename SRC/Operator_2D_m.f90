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
MODULE Operator_2D_m
  !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
  USE QDUtil_m
  IMPLICIT NONE

  
  PRIVATE
  
  PUBLIC MolecCav_Action_operator_2D, MolecCav_Average_value_operator_2D
  
  INTERFACE Action_operator_2D
    MODULE PROCEDURE MolecCav_Action_operator_2D
  END INTERFACE
  INTERFACE Average_value_operator_2D
    MODULE PROCEDURE MolecCav_Average_value_operator_2D
  END INTERFACE
      
  
  CONTAINS
  
  
  SUBROUTINE MolecCav_Action_operator_2D(Op_psi, Operator, Psi)   ! /!\ FOR NOW EVERYTHING IS REAL /!\ compute the resulting vector Psi_result(:) from the action of the operator of the cavity mode on the photon state vector Psi_argument(:) written in the Eigenbasis of H_ho
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    USE Cavity_mode_m
    USE Operator_1D_m
    IMPLICIT NONE
  
    real(kind=Rkind),    intent(inout) :: Op_psi(:,:)                      ! already allocated !
    TYPE(Operator_1D_t), intent(in)    :: Operator
    real(kind=Rkind),    intent(in)    :: Psi(:,:)
  
    integer                            :: Nb_M, Nb_C, i_M
 
    Nb_M = Size(Psi, 1)
    Nb_C = Size(Psi, 2)
    
    IF (Nb_M /= Size(Op_psi, 1) .OR. Nb_C /= Size(Op_psi, 2)) THEN
      WRITE(out_unit,*) "The size of the operand Psi's vector does not match the size of the resulting vector Op_psi."
      STOP "### The size of the operand Psi's vector does not match the size of the resulting vector Op_psi."
    END IF

    DO i_M = 1, Nb_M
      CALL Action_Operator_1D(Op_psi(i_M, :), Operator, Psi(i_M, :))
    END DO

  END SUBROUTINE MolecCav_Action_operator_2D


  SUBROUTINE MolecCav_Average_value_operator_2D(Value, Operator, Psi)   ! /!\ FOR NOW EVERYTHING IS REAL /!\ compute the resulting vector Psi_result(:) from the action of the operator of the cavity mode on the photon state vector Psi_argument(:) written in the Eigenbasis of H_ho
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    USE Algebra_m
    USE Cavity_mode_m
    USE Operator_1D_m
    IMPLICIT NONE

    real(kind=Rkind),    intent(inout) :: Value
    TYPE(Operator_1D_t), intent(in)    :: Operator
    real(kind=Rkind),    intent(in)    :: Psi(:,:)

    real(kind=Rkind), allocatable      :: Intermediary(:,:)
    integer                            :: Nb_M, Nb_C

    Nb_M = Size(Psi, 1)
    Nb_C = Size(Psi, 2)
    ALLOCATE(Intermediary(Nb_M, Nb_C))

    CALL Action_operator_2D(Intermediary, Operator, Psi)
    CALL Scalar_product(Value, Intermediary, Psi)

  END SUBROUTINE MolecCav_Average_value_operator_2D
  

END MODULE
