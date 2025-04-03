!==================================================================================================
!==================================================================================================
! This file is part of MolecCav.
!
!==================================================================================================
! MIT License
!
! Copyright (c) 2525 Mathéo Segaud
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
! The module that accounts for the actions of the operators related to a HO over an any statevector
! (described by rank-1 tensors) of this HO.  
! Action_HO1D_operator_R1        : computes the resulting vector Op_psi (Rank-1 tensor) from the a-
! ction of the operator of the 1D HO on the state vector Psi(rank-1 tensor) written in the Eigenba-
! sis of the H of the HO1D. N.B. The wavefunction must be a statefunction OF THIS HO. Calls the ne-
! -xt three procedures.
! Action_dense_HO1D_operator_R1  : the procedure that actually computes the action of Action_HO1D_-
! operator_R1, in the case of an operator represented by its full (analytical) matrix.
! Action_diag_HO1D_operator_R1   : the procedure that actually computes the action of Action_HO1D_-
! operator_R1, in the case of an operator represented by by the vector (rank-1 tensor) of the diag-
! onal elements of its full (analytical) matrix.
! Action_band_HO1D_operator_R1   : the procedure that actually computes the action of Action_HO1D_-
! operator_R1, in the case of an operator represented by matrix (rank-2 tensor) with a column for 
! each of its non-null diagonal.
! Average_value_HO1D_operator_R1 : computes the avarage value of an operator of the 1D HO over a s-
! tatevector (rank-1 tensor), calling first the Action_HO1D_operator_R1 and then projecting the th-
! ereby obtained wavevector on the initial one.
!==================================================================================================
!==================================================================================================
MODULE Actions_HO1D_operators_R1_m
  !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
  USE QDUtil_m
  IMPLICIT NONE


  PRIVATE

  PUBLIC Action_HO1D_operator_R1, Average_value_HO1D_operator_R1

  INTERFACE Action_HO1D_operator_R1
    MODULE PROCEDURE MolecCav_Action_HO1D_operator_R1_real, MolecCav_Action_HO1D_operator_R1_complex
  END INTERFACE
  INTERFACE Action_dense_HO1D_operator_R1
    MODULE PROCEDURE MolecCav_Action_dense_HO1D_operator_R1_real, MolecCav_Action_dense_HO1D_operator_R1_complex
  END INTERFACE
  INTERFACE Action_diag_HO1D_operator_R1
    MODULE PROCEDURE MolecCav_Action_diag_HO1D_operator_R1_real, MolecCav_Action_diag_HO1D_operator_R1_complex
  END INTERFACE
  INTERFACE Action_band_HO1D_operator_R1
    MODULE PROCEDURE MolecCav_Action_band_HO1D_operator_R1_real, MolecCav_Action_band_HO1D_operator_R1_complex
  END INTERFACE
  INTERFACE Average_value_HO1D_operator_R1
    MODULE PROCEDURE MolecCav_Average_value_HO1D_operator_R1_real, MolecCav_Average_value_HO1D_operator_R1_complex
  END INTERFACE
    

  CONTAINS


  SUBROUTINE MolecCav_Action_HO1D_operator_R1_real(Op_psi, Operator, Psi, Verbose, Debug)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    USE HO1D_operator_m
    IMPLICIT NONE

    real(kind=Rkind),      intent(inout) :: Op_psi(:)
    TYPE(HO1D_operator_t), intent(in)    :: Operator
    real(kind=Rkind),      intent(in)    :: Psi(:)
    integer, optional,     intent(in)    :: Verbose                                                                              ! cf. comments in HO1D_parameters_m
    logical, optional,     intent(in)    :: Debug                                                                                ! cf. comments in HO1D_parameters_m

    integer                              :: Nb
    integer                              :: Verbose_local = 25                                                                   ! goes from 25 (= 0 verbose) to 29 (= maximum verbose) at this layer
    logical                              :: Debug_local   = .FALSE.

    !------------------------------------------------------Debugging options-----------------------------------------------------
    IF (PRESENT(Verbose)) Verbose_local = Verbose
    IF (PRESENT(Debug))   Debug_local   = Debug

    IF (Verbose_local > 25) WRITE(out_unit,*) 
    IF (Verbose_local > 25) WRITE(out_unit,*) "---------------------------------------COMPUTING ACTION OF THE HO1D OPERATOR OVER &
                                              &THE R1 WF---------------------------------------"; FLUSH(out_unit)

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "--- Arguments of MolecCav_Action_HO1D_operator_R1_real :"
      WRITE(out_unit,*) "The <<Operator>> argument :"
      CALL Write_HO1D_operator(Operator)
      WRITE(out_unit,*) "The <<Psi>> argument : "
      CALL Write_Vec(Psi, out_unit, 1, info="Psi")
      WRITE(out_unit,*) "The size of its vector : "//TO_string(Size(Psi))
      WRITE(out_unit,*) "--- End arguments of MolecCav_Action_HO1D_operator_R1_real"
      FLUSH(out_unit)
    END IF
    
    !-----------------------------------------------------Checking dimensions----------------------------------------------------
    IF (ALLOCATED(Operator%Dense_val)) Nb = Size(Operator%Dense_val, dim=1)
    IF (ALLOCATED(Operator%Diag_val) ) Nb = Size(Operator%Diag_val,  dim=1)
    IF (ALLOCATED(Operator%Band_val) ) Nb = Size(Operator%Band_val,  dim=1)

    IF (Nb /= Size(Psi)) THEN
      WRITE(out_unit,*) "### The dimensions of the Operator's matrix representation does not match the operand Psi's vector size.&
                       & Please check initialization."
      WRITE(out_unit,*) "    Size(Operator%matrix, dim=1) = "//TO_string(Nb)//"; Size(Psi) = "//TO_string(Size(Psi))
      STOP "### The dimensions of the Operator's matrix representation does not match the operand Psi's vector size. Please check& 
                       & initialization."
    END IF 

    IF (Nb /= Size(Op_psi)) THEN
      WRITE(out_unit,*) "### The dimensions of the Operator's matrix representation does not match the resulting Op_psi vector's & 
                       &size. Please check initialization."
      WRITE(out_unit,*) "    Size(Operator%matrix, dim=1) = "//TO_string(Nb)//"; Size(Op_psi) = "//TO_string(Size(Op_psi))
      STOP "### The dimensions of the Operator's matrix representation does not match the resulting Op_psi vector's size. Please & 
                       &check initialization."
    END IF 

    !---------------------------------------------Selection of the calculation method--------------------------------------------
    IF      (ALLOCATED(Operator%Diag_val))   THEN
      IF (Debug_local) WRITE(out_unit,*) "--- Operator%Diag_val is allocated. The diagonal operator action called."
      CALL Action_diag_HO1D_operator_R1(Op_psi=Op_psi, Operator=Operator, Psi=Psi, Verbose=Verbose_local, Debug=Debug_local)

    ELSE IF (ALLOCATED(Operator%Band_val))   THEN
      IF (Debug_local) WRITE(out_unit,*) "--- Operator%Diag_val is not allocated."
      IF (Debug_local) WRITE(out_unit,*) "--- Operator%Band_val is allocated. The band operator action called."
      CALL Action_band_HO1D_operator_R1(Op_psi=Op_psi, Operator=Operator, Psi=Psi, Verbose=Verbose_local, Debug=Debug_local)

    ELSE IF (ALLOCATED(Operator%Dense_val)) THEN
      IF (Debug_local) WRITE(out_unit,*) "Operator%Diag_val and Operator%Band_val are not allocated."
      IF (Debug_local) WRITE(out_unit,*) "Operator%Dense_val is allocated. The dense operator action called."
      CALL Action_dense_HO1D_operator_R1(Op_psi=Op_psi, Operator=Operator, Psi=Psi, Verbose=Verbose_local, Debug=Debug_local)

    ELSE
      WRITE(out_unit,*) "### None of this operator's matrices are allocated. Please check its initalization."
      STOP "### None of this operator's matrices are allocated. Please check its initalization."
    END IF

    IF (Verbose_local > 26) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "--- Resulting statevector from the action of the HO1D Operator on the Psi statevector operand, computed &
                        &by Action_HO1D_operator_R1 :"
      CALL Write_Vec(Op_psi, out_unit, 1, info="Op_Psi")
      WRITE(out_unit,*) "--- End resulting statevector computed by Action_HO1D_operator_R1"
    END IF
  
    IF (Verbose_local > 25) WRITE(out_unit,*) 
    IF (Verbose_local > 25) WRITE(out_unit,*) "----------------------------------------ACTION OF THE HO1D OPERATOR OVER THE R1 WF&
                                              & COMPUTED---------------------------------------"; FLUSH(out_unit)
  
  END SUBROUTINE MolecCav_Action_HO1D_operator_R1_real

  
  SUBROUTINE MolecCav_Action_dense_HO1D_operator_R1_real(Op_psi, Operator, Psi, Verbose, Debug)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    USE HO1D_operator_m
    IMPLICIT NONE
    
    real(kind=Rkind),      intent(inout) :: Op_psi(:)
    TYPE(HO1D_operator_t), intent(in)    :: Operator
    real(kind=Rkind),      intent(in)    :: Psi(:)
    integer, optional,     intent(in)    :: Verbose                                                                              ! cf. comments in HO1D_parameters_m
    logical, optional,     intent(in)    :: Debug                                                                                ! cf. comments in HO1D_parameters_m

    integer                              :: Verbose_local = 25                                                                   ! goes from 25 (= 0 verbose) to 29 (= maximum verbose) at this layer
    logical                              :: Debug_local   = .FALSE.

    !------------------------------------------------------Debugging options-----------------------------------------------------
    IF (PRESENT(Verbose)) Verbose_local = Verbose
    IF (PRESENT(Debug))   Debug_local   = Debug

    IF (Verbose_local > 27) WRITE(out_unit,*) 
    IF (Verbose_local > 27) WRITE(out_unit,*) "---------------------------------------Computing the HO1D operator using the dense&
                                             & procedure--------------------------------------"; FLUSH(out_unit)
    
    !----------------------------------------------------Computing the action----------------------------------------------------
    Op_psi(:) = matmul(Operator%Dense_val, Psi)

  END SUBROUTINE MolecCav_Action_dense_HO1D_operator_R1_real

  
  SUBROUTINE MolecCav_Action_diag_HO1D_operator_R1_real(Op_psi, Operator, Psi, Verbose, Debug) 
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    USE HO1D_operator_m
    IMPLICIT NONE
    
    real(kind=Rkind),      intent(inout) :: Op_psi(:)
    TYPE(HO1D_operator_t), intent(in)    :: Operator
    real(kind=Rkind),      intent(in)    :: Psi(:)
    integer, optional,     intent(in)    :: Verbose                                                                              ! cf. comments in HO1D_parameters_m
    logical, optional,     intent(in)    :: Debug                                                                                ! cf. comments in HO1D_parameters_m

    integer                              :: Verbose_local = 25                                                                   ! goes from 25 (= 0 verbose) to 29 (= maximum verbose) at this layer
    logical                              :: Debug_local   = .FALSE.

    !------------------------------------------------------Debugging options-----------------------------------------------------
    IF (PRESENT(Verbose)) Verbose_local = Verbose
    IF (PRESENT(Debug))   Debug_local   = Debug

    IF (Verbose_local > 27) WRITE(out_unit,*) 
    IF (Verbose_local > 27) WRITE(out_unit,*) "---------------------------------------Computing the HO1D operator using the diago&
                                             &nal procedure--------------------------------------"; FLUSH(out_unit)
    
    !----------------------------------------------------Computing the action----------------------------------------------------
    Op_psi = Operator%Diag_val * Psi

  END SUBROUTINE MolecCav_Action_diag_HO1D_operator_R1_real

  
  SUBROUTINE MolecCav_Action_band_HO1D_operator_R1_real(Op_psi, Operator, Psi, Verbose, Debug)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    USE HO1D_operator_m
    IMPLICIT NONE
    
    real(kind=Rkind),      intent(inout) :: Op_psi(:)
    TYPE(HO1D_operator_t), intent(in)    :: Operator
    real(kind=Rkind),      intent(in)    :: Psi(:)
    integer, optional,     intent(in)    :: Verbose                                                                              ! cf. comments in HO1D_parameters_m
    logical, optional,     intent(in)    :: Debug                                                                                ! cf. comments in HO1D_parameters_m

    integer                              :: i, Nb
    integer                              :: Verbose_local = 25                                                                   ! goes from 25 (= 0 verbose) to 29 (= maximum verbose) at this layer
    logical                              :: Debug_local   = .FALSE.

    !------------------------------------------------------Debugging options-----------------------------------------------------
    IF (PRESENT(Verbose)) Verbose_local = Verbose
    IF (PRESENT(Debug))   Debug_local   = Debug

    IF (Verbose_local > 27) WRITE(out_unit,*) 
    IF (Verbose_local > 27) WRITE(out_unit,*) "---------------------------------------Computing the HO1D operator using the band &
                                             &procedure--------------------------------------"; FLUSH(out_unit)
    
    !----------------------------------------------------Computing the action----------------------------------------------------
    Nb = size(Op_psi)

    Op_psi     = ZERO
    Op_psi     = Operator%Band_val(:,2) * Psi
    Op_psi(1)  = Op_psi(1)  + Operator%Band_val(2,3)    * Psi(2)
    Op_psi(Nb) = Op_psi(Nb) + Operator%Band_val(Nb-1,1) * Psi(Nb-1)
    DO i = 2, Nb-1
      Op_psi(i) = Op_psi(i) + &
                & Operator%Band_val(i-1,1) * Psi(i-1) + &
                & Operator%Band_val(i+1,3) * Psi(i+1)
    END DO

  END SUBROUTINE MolecCav_Action_band_HO1D_operator_R1_real
  

  SUBROUTINE MolecCav_Average_value_HO1D_operator_R1_real(Value, Operator, Psi, Verbose, Debug)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    USE Algebra_m
    USE HO1D_operator_m
    IMPLICIT NONE

    real(kind=Rkind),      intent(inout) :: Value
    TYPE(HO1D_operator_t), intent(in)    :: Operator
    real(kind=Rkind),      intent(in)    :: Psi(:)
    integer, optional,     intent(in)    :: Verbose                                                                              ! cf. comments in HO1D_parameters_m
    logical, optional,     intent(in)    :: Debug                                                                                ! cf. comments in HO1D_parameters_m

    real(kind=Rkind), allocatable        :: Intermediary(:)
    integer                              :: Nb
    integer                              :: Verbose_local = 25                                                                   ! goes from 25 (= 0 verbose) to 29 (= maximum verbose) at this layer
    logical                              :: Debug_local   = .FALSE.

    !------------------------------------------------------Debugging options-----------------------------------------------------
    IF (PRESENT(Verbose)) Verbose_local = Verbose
    IF (PRESENT(Debug))   Debug_local   = Debug

    IF (Verbose_local > 25) WRITE(out_unit,*) 
    IF (Verbose_local > 25) WRITE(out_unit,*) "---------------------------------------COMPUTING THE AVERAGE VALUE OF THE HO1D OPE&
                                              &RATOR OVER THE R1 WF---------------------------------------"; FLUSH(out_unit)

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "--- Arguments of MolecCav_Average_value_HO1D_operator_R1_real :"
      WRITE(out_unit,*) "The <<Operator>> argument :"
      CALL Write_HO1D_operator(Operator)
      WRITE(out_unit,*) "The <<Psi>> argument : "
      CALL Write_Vec(Psi, out_unit, 1, info="Psi")
      WRITE(out_unit,*) "The size of its vector : "//TO_string(Size(Psi))
      WRITE(out_unit,*) "--- End arguments of MolecCav_Action_HO1D_operator_R1_real"
      FLUSH(out_unit)
    END IF
    
    !-----------------------------------------------------Checking dimensions----------------------------------------------------
    IF (ALLOCATED(Operator%Dense_val)) Nb = Size(Operator%Dense_val, dim=1)
    IF (ALLOCATED(Operator%Diag_val )) Nb = Size(Operator%Diag_val,  dim=1)
    IF (ALLOCATED(Operator%Band_val )) Nb = Size(Operator%Band_val,  dim=1)

    IF (Nb /= Size(Psi)) THEN
      WRITE(out_unit,*) "### The dimensions of the Operator's matrix representation does not match the operand Psi's vector size.&
                       & Please check initialization."
      WRITE(out_unit,*) "    Size(Operator%matrix, dim=1) = "//TO_string(Nb)//"; Size(Psi) = "//TO_string(Size(Psi))
      STOP "### The dimensions of the Operator's matrix representation does not match the operand Psi's vector size. Please check& 
                       & initialization."
    END IF 

    !-------------------------------------------------Computing the average value------------------------------------------------
    ALLOCATE(Intermediary(Nb))

    CALL Action_HO1D_operator_R1(Intermediary, Operator, Psi)
    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "--- Computed Intermediary statevector = \hat{Operator}|Psi> :"
      CALL Write_Vec(Intermediary, out_unit, 1, info="Intermediary")
      WRITE(out_unit,*) "--- End Intermediary statevector"
    END IF 

    CALL Scalar_product(Value, Psi, Intermediary)
    IF (Verbose_local > 26) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "--- Resulting average value of the HO1D Operator on the Psi statevector operand, computed by MolecCav_Av&
                        &erage_value_HO1D_operator_R1_real :"//TO_string(Value)
      WRITE(out_unit,*) "--- End resulting average value computed by MolecCav_Average_value_HO1D_operator_R1_real"
    END IF

    IF (Verbose_local > 25) WRITE(out_unit,*) 
    IF (Verbose_local > 25) WRITE(out_unit,*) "------------------------------------------------------AVERAGE VALUE COMPUTED------&
                                              &------------------------------------------------"; FLUSH(out_unit)

    DEALLOCATE(Intermediary)
    
  END SUBROUTINE MolecCav_Average_value_HO1D_operator_R1_real


  SUBROUTINE MolecCav_Action_HO1D_operator_R1_complex(Op_psi, Operator, Psi, Verbose, Debug)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    USE HO1D_operator_m
    IMPLICIT NONE

    complex(kind=Rkind),   intent(inout) :: Op_psi(:)
    TYPE(HO1D_operator_t), intent(in)    :: Operator
    complex(kind=Rkind),   intent(in)    :: Psi(:)
    integer, optional,     intent(in)    :: Verbose                                                                              ! cf. comments in HO1D_parameters_m
    logical, optional,     intent(in)    :: Debug                                                                                ! cf. comments in HO1D_parameters_m

    integer                              :: Nb
    integer                              :: Verbose_local = 25                                                                   ! goes from 25 (= 0 verbose) to 29 (= maximum verbose) at this layer
    logical                              :: Debug_local   = .FALSE.

    !------------------------------------------------------Debugging options-----------------------------------------------------
    IF (PRESENT(Verbose)) Verbose_local = Verbose
    IF (PRESENT(Debug))   Debug_local   = Debug

    IF (Verbose_local > 25) WRITE(out_unit,*) 
    IF (Verbose_local > 25) WRITE(out_unit,*) "---------------------------------------COMPUTING ACTION OF THE HO1D OPERATOR OVER &
                                              &THE R1 WF---------------------------------------"; FLUSH(out_unit)

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "--- Arguments of MolecCav_Action_HO1D_operator_R1_complex :"
      WRITE(out_unit,*) "The <<Operator>> argument :"
      CALL Write_HO1D_operator(Operator)
      WRITE(out_unit,*) "The <<Psi>> argument : "
      CALL Write_Vec(Psi, out_unit, 1, info="Psi")
      WRITE(out_unit,*) "The size of its vector : "//TO_string(Size(Psi))
      WRITE(out_unit,*) "--- End arguments of MolecCav_Action_HO1D_operator_R1_complex"
      FLUSH(out_unit)
    END IF
    
    !-----------------------------------------------------Checking dimensions----------------------------------------------------
    IF (ALLOCATED(Operator%Dense_val)) Nb = Size(Operator%Dense_val, dim=1)
    IF (ALLOCATED(Operator%Diag_val) ) Nb = Size(Operator%Diag_val,  dim=1)
    IF (ALLOCATED(Operator%Band_val) ) Nb = Size(Operator%Band_val,  dim=1)

    IF (Nb /= Size(Psi)) THEN
      WRITE(out_unit,*) "### The dimensions of the Operator's matrix representation does not match the operand Psi's vector size.&
                       & Please check initialization."
      WRITE(out_unit,*) "    Size(Operator%matrix, dim=1) = "//TO_string(Nb)//"; Size(Psi) = "//TO_string(Size(Psi))
      STOP "### The dimensions of the Operator's matrix representation does not match the operand Psi's vector size. Please check& 
                       & initialization."
    END IF 

    IF (Nb /= Size(Op_psi)) THEN
      WRITE(out_unit,*) "### The dimensions of the Operator's matrix representation does not match the resulting Op_psi vector's & 
                       &size. Please check initialization."
      WRITE(out_unit,*) "    Size(Operator%matrix, dim=1) = "//TO_string(Nb)//"; Size(Op_psi) = "//TO_string(Size(Op_psi))
      STOP "### The dimensions of the Operator's matrix representation does not match the resulting Op_psi vector's size. Please & 
                       &check initialization."
    END IF 

    !---------------------------------------------Selection of the calculation method--------------------------------------------
    IF      (ALLOCATED(Operator%Diag_val))   THEN
      IF (Debug_local) WRITE(out_unit,*) "--- Operator%Diag_val is allocated. The diagonal operator action called."
      CALL Action_diag_HO1D_operator_R1(Op_psi=Op_psi, Operator=Operator, Psi=Psi, Verbose=Verbose_local, Debug=Debug_local)

    ELSE IF (ALLOCATED(Operator%Band_val))   THEN
      IF (Debug_local) WRITE(out_unit,*) "--- Operator%Diag_val is not allocated."
      IF (Debug_local) WRITE(out_unit,*) "--- Operator%Band_val is allocated. The band operator action called."
      CALL Action_band_HO1D_operator_R1(Op_psi=Op_psi, Operator=Operator, Psi=Psi, Verbose=Verbose_local, Debug=Debug_local)

    ELSE IF (ALLOCATED(Operator%Dense_val)) THEN
      IF (Debug_local) WRITE(out_unit,*) "Operator%Diag_val and Operator%Band_val are not allocated."
      IF (Debug_local) WRITE(out_unit,*) "Operator%Dense_val is allocated. The dense operator action called."
      CALL Action_dense_HO1D_operator_R1(Op_psi=Op_psi, Operator=Operator, Psi=Psi, Verbose=Verbose_local, Debug=Debug_local)

    ELSE
      WRITE(out_unit,*) "### None of this operator's matrices are allocated. Please check its initalization."
      STOP "### None of this operator's matrices are allocated. Please check its initalization."
    END IF

    IF (Verbose_local > 26) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "--- Resulting statevector from the action of the HO1D Operator on the Psi statevector operand, computed &
                        &by Action_HO1D_operator_R1 :"
      CALL Write_Vec(Op_psi, out_unit, 1, info="Op_Psi")
      WRITE(out_unit,*) "--- End resulting statevector computed by Action_HO1D_operator_R1"
    END IF
  
    IF (Verbose_local > 25) WRITE(out_unit,*) 
    IF (Verbose_local > 25) WRITE(out_unit,*) "----------------------------------------ACTION OF THE HO1D OPERATOR OVER THE R1 WF&
                                              & COMPUTED---------------------------------------"; FLUSH(out_unit)
  
  END SUBROUTINE MolecCav_Action_HO1D_operator_R1_complex

  
  SUBROUTINE MolecCav_Action_dense_HO1D_operator_R1_complex(Op_psi, Operator, Psi, Verbose, Debug)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    USE HO1D_operator_m
    IMPLICIT NONE
    
    complex(kind=Rkind),   intent(inout) :: Op_psi(:)
    TYPE(HO1D_operator_t), intent(in)    :: Operator
    complex(kind=Rkind),   intent(in)    :: Psi(:)
    integer, optional,     intent(in)    :: Verbose                                                                              ! cf. comments in HO1D_parameters_m
    logical, optional,     intent(in)    :: Debug                                                                                ! cf. comments in HO1D_parameters_m

    integer                              :: Verbose_local = 25                                                                   ! goes from 25 (= 0 verbose) to 29 (= maximum verbose) at this layer
    logical                              :: Debug_local   = .FALSE.

    !------------------------------------------------------Debugging options-----------------------------------------------------
    IF (PRESENT(Verbose)) Verbose_local = Verbose
    IF (PRESENT(Debug))   Debug_local   = Debug

    IF (Verbose_local > 27) WRITE(out_unit,*) 
    IF (Verbose_local > 27) WRITE(out_unit,*) "---------------------------------------Computing the HO1D operator using the dense&
                                             & procedure--------------------------------------"; FLUSH(out_unit)
    
    !----------------------------------------------------Computing the action----------------------------------------------------
    Op_psi(:) = matmul(Operator%Dense_val, Psi)

  END SUBROUTINE MolecCav_Action_dense_HO1D_operator_R1_complex

  
  SUBROUTINE MolecCav_Action_diag_HO1D_operator_R1_complex(Op_psi, Operator, Psi, Verbose, Debug) 
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    USE HO1D_operator_m
    IMPLICIT NONE
    
    complex(kind=Rkind),   intent(inout) :: Op_psi(:)
    TYPE(HO1D_operator_t), intent(in)    :: Operator
    complex(kind=Rkind),   intent(in)    :: Psi(:)
    integer, optional,     intent(in)    :: Verbose                                                                              ! cf. comments in HO1D_parameters_m
    logical, optional,     intent(in)    :: Debug                                                                                ! cf. comments in HO1D_parameters_m

    integer                              :: Verbose_local = 25                                                                   ! goes from 25 (= 0 verbose) to 29 (= maximum verbose) at this layer
    logical                              :: Debug_local   = .FALSE.

    !------------------------------------------------------Debugging options-----------------------------------------------------
    IF (PRESENT(Verbose)) Verbose_local = Verbose
    IF (PRESENT(Debug))   Debug_local   = Debug

    IF (Verbose_local > 27) WRITE(out_unit,*) 
    IF (Verbose_local > 27) WRITE(out_unit,*) "---------------------------------------Computing the HO1D operator using the diago&
                                             &nal procedure--------------------------------------"; FLUSH(out_unit)
    
    !----------------------------------------------------Computing the action----------------------------------------------------
    Op_psi = Operator%Diag_val * Psi

  END SUBROUTINE MolecCav_Action_diag_HO1D_operator_R1_complex

  
  SUBROUTINE MolecCav_Action_band_HO1D_operator_R1_complex(Op_psi, Operator, Psi, Verbose, Debug)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    USE HO1D_operator_m
    IMPLICIT NONE
    
    complex(kind=Rkind),   intent(inout) :: Op_psi(:)
    TYPE(HO1D_operator_t), intent(in)    :: Operator
    complex(kind=Rkind),   intent(in)    :: Psi(:)
    integer, optional,     intent(in)    :: Verbose                                                                              ! cf. comments in HO1D_parameters_m
    logical, optional,     intent(in)    :: Debug                                                                                ! cf. comments in HO1D_parameters_m

    integer                              :: i, Nb
    integer                              :: Verbose_local = 25                                                                   ! goes from 25 (= 0 verbose) to 29 (= maximum verbose) at this layer
    logical                              :: Debug_local   = .FALSE.

    !------------------------------------------------------Debugging options-----------------------------------------------------
    IF (PRESENT(Verbose)) Verbose_local = Verbose
    IF (PRESENT(Debug))   Debug_local   = Debug

    IF (Verbose_local > 27) WRITE(out_unit,*) 
    IF (Verbose_local > 27) WRITE(out_unit,*) "---------------------------------------Computing the HO1D operator using the band &
                                             &procedure--------------------------------------"; FLUSH(out_unit)
    
    !----------------------------------------------------Computing the action----------------------------------------------------
    Nb = size(Op_psi)

    Op_psi     = ZERO
    Op_psi     = Operator%Band_val(:,2) * Psi
    Op_psi(1)  = Op_psi(1)  + Operator%Band_val(2,3)    * Psi(2)
    Op_psi(Nb) = Op_psi(Nb) + Operator%Band_val(Nb-1,1) * Psi(Nb-1)
    DO i = 2, Nb-1
      Op_psi(i) = Op_psi(i) + &
                & Operator%Band_val(i-1,1) * Psi(i-1) + &
                & Operator%Band_val(i+1,3) * Psi(i+1)
    END DO

  END SUBROUTINE MolecCav_Action_band_HO1D_operator_R1_complex
  

  SUBROUTINE MolecCav_Average_value_HO1D_operator_R1_complex(Value, Operator, Psi, Verbose, Debug)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    USE Algebra_m
    USE HO1D_operator_m
    IMPLICIT NONE

    complex(kind=Rkind),   intent(inout) :: Value
    TYPE(HO1D_operator_t), intent(in)    :: Operator
    complex(kind=Rkind),   intent(in)    :: Psi(:)
    integer, optional,     intent(in)    :: Verbose                                                                              ! cf. comments in HO1D_parameters_m
    logical, optional,     intent(in)    :: Debug                                                                                ! cf. comments in HO1D_parameters_m

    complex(kind=Rkind), allocatable     :: Intermediary(:)
    integer                              :: Nb
    integer                              :: Verbose_local = 25                                                                   ! goes from 25 (= 0 verbose) to 29 (= maximum verbose) at this layer
    logical                              :: Debug_local   = .FALSE.

    !------------------------------------------------------Debugging options-----------------------------------------------------
    IF (PRESENT(Verbose)) Verbose_local = Verbose
    IF (PRESENT(Debug))   Debug_local   = Debug

    IF (Verbose_local > 25) WRITE(out_unit,*) 
    IF (Verbose_local > 25) WRITE(out_unit,*) "---------------------------------------COMPUTING THE AVERAGE VALUE OF THE HO1D OPE&
                                              &RATOR OVER THE R1 WF---------------------------------------"; FLUSH(out_unit)

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "--- Arguments of MolecCav_Average_value_HO1D_operator_R1_complex :"
      WRITE(out_unit,*) "The <<Operator>> argument :"
      CALL Write_HO1D_operator(Operator)
      WRITE(out_unit,*) "The <<Psi>> argument : "
      CALL Write_Vec(Psi, out_unit, 1, info="Psi")
      WRITE(out_unit,*) "The size of its vector : "//TO_string(Size(Psi))
      WRITE(out_unit,*) "--- End arguments of MolecCav_Action_HO1D_operator_R1_complex"
      FLUSH(out_unit)
    END IF
    
    !-----------------------------------------------------Checking dimensions----------------------------------------------------
    IF (ALLOCATED(Operator%Dense_val)) Nb = Size(Operator%Dense_val, dim=1)
    IF (ALLOCATED(Operator%Diag_val )) Nb = Size(Operator%Diag_val,  dim=1)
    IF (ALLOCATED(Operator%Band_val )) Nb = Size(Operator%Band_val,  dim=1)

    IF (Nb /= Size(Psi)) THEN
      WRITE(out_unit,*) "### The dimensions of the Operator's matrix representation does not match the operand Psi's vector size.&
                       & Please check initialization."
      WRITE(out_unit,*) "    Size(Operator%matrix, dim=1) = "//TO_string(Nb)//"; Size(Psi) = "//TO_string(Size(Psi))
      STOP "### The dimensions of the Operator's matrix representation does not match the operand Psi's vector size. Please check& 
                       & initialization."
    END IF 

    !-------------------------------------------------Computing the average value------------------------------------------------
    ALLOCATE(Intermediary(Nb))

    CALL Action_HO1D_operator_R1(Intermediary, Operator, Psi)
    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "--- Computed Intermediary statevector = \hat{Operator}|Psi> :"
      CALL Write_Vec(Intermediary, out_unit, 1, info="Intermediary")
      WRITE(out_unit,*) "--- End Intermediary statevector"
    END IF 

    CALL Scalar_product(Value, Psi, Intermediary)
    IF (Verbose_local > 26) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "--- Resulting average value of the HO1D Operator on the Psi statevector operand, computed by MolecCav_Av&
                        &erage_value_HO1D_operator_R1_complex :"//TO_string(Value)
      WRITE(out_unit,*) "--- End resulting average value computed by MolecCav_Average_value_HO1D_operator_R1_complex"
    END IF

    IF (Verbose_local > 25) WRITE(out_unit,*) 
    IF (Verbose_local > 25) WRITE(out_unit,*) "------------------------------------------------------AVERAGE VALUE COMPUTED------&
                                              &------------------------------------------------"; FLUSH(out_unit)

    DEALLOCATE(Intermediary)
    
  END SUBROUTINE MolecCav_Average_value_HO1D_operator_R1_complex


END MODULE
