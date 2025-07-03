!==================================================================================================
!==================================================================================================
!
! This file is part of MolecCav.
!
!==================================================================================================
!
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
!
!==================================================================================================
!
! README :
! The module to initialize the operators related to a HO.
! Initialize_HO1D_operator : constructs the operator using parameters of the HO1D_para object from 
! the so called derived type. Calls the next three procedures. 
! Initialize_H_HO1D        : constructs the Hamiltonian operator using parameters of the HO1D_para 
! object from the so called derived type.  
! Initialize_x_HO1D        : constructs the Position operator using parameters of the HO1D_para ob-
! ject from the so called derived type.
! Initialize_N_HO1D        : constructs the Number of excitation Quanta operator using parameters 
! of the HO1D_para object from the so called derived type.
! Write      : displays values of the type in the output.
! Deallocate_HO1D_operator : deallocates all tables of the type.
! README :
! The module that accounts for the actions of the operators related to a HO over an any statevector
! (described by rank-1 tensors) of this HO.  
! Action_HO1D_operator_R1        : computes the resulting vector Op_psi (Rank-1 tensor) from the a-
! ction of the operator of the 1D HO on the state vector Psi(rank-1 tensor) written in the Eigenba-
! sis of the H of the HO1D. N.B. The wavefunction must be a statefunction OF THIS HO. Calls the ne-
! -xt three procedures.
! Action_dense  : the procedure that actually computes the action of Action_HO1D_-
! operator_R1, in the case of an operator represented by its full (analytical) matrix.
! Action_diag   : the procedure that actually computes the action of Action_HO1D_-
! operator_R1, in the case of an operator represented by by the vector (rank-1 tensor) of the diag-
! onal elements of its full (analytical) matrix.
! Action_band   : the procedure that actually computes the action of Action_HO1D_-
! operator_R1, in the case of an operator represented by matrix (rank-2 tensor) with a column for 
! each of its non-null diagonal.
! Average_value_HO1D_operator_R1 : computes the avarage value of an operator of the 1D HO over a s-
! tatevector (rank-1 tensor), calling first the Action_HO1D_operator_R1 and then projecting the th-
! ereby obtained wavevector on the initial one.
!==================================================================================================
!==================================================================================================
MODULE Elem_op_m
  !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
  USE QDUtil_m
  IMPLICIT NONE


  TYPE                               :: Elem_op_t
    character(len=:),    allocatable :: Operator_type                                                                            ! ex : "Hamiltonian", "Position", "NbQuanta", etc
    logical                          :: Dense = .FALSE.                                                                          ! if .TRUE. then the matrix storage will not be optimized and it will be stored as a Dense matrix. Otherwise, the elements matrix will be stored in a table such as taking advantage of the sparcity of the operator's analytical matrix in the HO Eigenbasis
    logical                          :: Grid  = .FALSE.                                                                          ! if .TRUE. then the matrix storage will not be optimized and it will be stored as a Dense matrix. Otherwise, the elements matrix will be stored in a table such as taking advantage of the sparcity of the operator's analytical matrix in the HO Eigenbasis
    integer                          :: Upper_bandwidth = 0                                                                      ! if type = "Band". Gives the number of additional bands to consider above the diagonal.
    integer                          :: Lower_bandwidth = 0                                                                      ! if type = "Band". Gives the number of additional bands to consider below the diagonal. Ex : Upper_bandwidth=Lower_bandwidth=1 would give a tridiagonal matrix
    real(kind=Rkind),    allocatable :: Dense_val(:,:)                                                                           ! if Dense == .TRUE.
    real(kind=Rkind),    allocatable :: Diag_val(:)                                                                              ! if Dense == .FALSE. .AND. the operator analytical matrix in the HO1D Eigenbasis is diagonal. The diagonal elements are stored in a vector (rank-1 tensor)
    real(kind=Rkind),    allocatable :: Band_val(:,:)                                                                            ! if Dense == .FALSE. .AND. the operator analytical matrix in the HO1D Eigenbasis is band. The number of columns will be the number of diagonals to consider : each considered diagonal is stored in a column
  END TYPE


  PRIVATE

  PUBLIC Elem_op_t, Action, Write, Dealloc


  INTERFACE Action
    MODULE PROCEDURE MolecCav_Action_elem_op_R1_real, MolecCav_Action_elem_op_R1_complex!&
  END INTERFACE
  INTERFACE Action_dense
    MODULE PROCEDURE MolecCav_Action_dense_elem_op_R1_real, MolecCav_Action_dense_elem_op_R1_complex
  END INTERFACE
  INTERFACE Action_diag
    MODULE PROCEDURE MolecCav_Action_diag_elem_op_R1_real,  MolecCav_Action_diag_elem_op_R1_complex
  END INTERFACE
  INTERFACE Action_band
    MODULE PROCEDURE MolecCav_Action_band_elem_op_R1_real,  MolecCav_Action_band_elem_op_R1_complex
  END INTERFACE
  INTERFACE Write
    MODULE PROCEDURE MolecCav_Write_elem_op_R1
  END INTERFACE
  INTERFACE Dealloc
    MODULE PROCEDURE MolecCav_Deallocate_elem_op
  END INTERFACE


  CONTAINS


  SUBROUTINE MolecCav_Action_elem_op_R1_real(Op_psi, Elem_op, Psi, Verbose, Debug)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    IMPLICIT NONE

    real(kind=Rkind),      intent(inout) :: Op_psi(:)
    TYPE(Elem_op_t),       intent(in)    :: Elem_op
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
      WRITE(out_unit,*) "--- Arguments of MolecCav_Action_elem_op_R1_real_real :"
      WRITE(out_unit,*) "The <<Elem_op>> argument :"
      CALL Write(Elem_op)
      WRITE(out_unit,*) "The <<Psi>> argument : "
      CALL Write_Vec(Psi, out_unit, 1, info="Psi")
      WRITE(out_unit,*) "The size of its vector : "//TO_string(Size(Psi))
      WRITE(out_unit,*) "--- End arguments of MolecCav_Action_elem_op_R1_real_real"
      FLUSH(out_unit)
    END IF
    
    !-----------------------------------------------------Checking dimensions----------------------------------------------------
    IF (ALLOCATED(Elem_op%Dense_val)) Nb = Size(Elem_op%Dense_val, dim=1)
    IF (ALLOCATED(Elem_op%Diag_val) ) Nb = Size(Elem_op%Diag_val,  dim=1)
    IF (ALLOCATED(Elem_op%Band_val) ) Nb = Size(Elem_op%Band_val,  dim=1)

    IF (Nb /= Size(Psi)) THEN
      WRITE(out_unit,*) "### The dimensions of the Elem_op's matrix representation does not match the operand Psi's vector size.&
                       & Please check initialization."
      WRITE(out_unit,*) "    Size(Elem_op%matrix, dim=1) = "//TO_string(Nb)//"; Size(Psi) = "//TO_string(Size(Psi))
      STOP "### The dimensions of the Elem_op's matrix representation does not match the operand Psi's vector size. Please check& 
                       & initialization."
    END IF 

    IF (Nb /= Size(Op_psi)) THEN
      WRITE(out_unit,*) "### The dimensions of the Elem_op's matrix representation does not match the resulting Op_psi vector's & 
                       &size. Please check initialization."
      WRITE(out_unit,*) "    Size(Elem_op%matrix, dim=1) = "//TO_string(Nb)//"; Size(Op_psi) = "//TO_string(Size(Op_psi))
      STOP "### The dimensions of the Elem_op's matrix representation does not match the resulting Op_psi vector's size. Please & 
                       &check initialization."
    END IF 

    !---------------------------------------------Selection of the calculation method--------------------------------------------
    IF      (ALLOCATED(Elem_op%Diag_val))   THEN
      IF (Debug_local) WRITE(out_unit,*) "--- Elem_op%Diag_val is allocated. The diagonal operator action called."
      CALL Action_diag(Op_psi=Op_psi, Elem_op=Elem_op, Psi=Psi, Verbose=Verbose_local, Debug=Debug_local)

    ELSE IF (ALLOCATED(Elem_op%Band_val))   THEN
      IF (Debug_local) WRITE(out_unit,*) "--- Elem_op%Diag_val is not allocated."
      IF (Debug_local) WRITE(out_unit,*) "--- Elem_op%Band_val is allocated. The band operator action called."
      CALL Action_band(Op_psi=Op_psi, Elem_op=Elem_op, Psi=Psi, Verbose=Verbose_local, Debug=Debug_local)

    ELSE IF (ALLOCATED(Elem_op%Dense_val)) THEN
      IF (Debug_local) WRITE(out_unit,*) "Elem_op%Diag_val and Elem_op%Band_val are not allocated."
      IF (Debug_local) WRITE(out_unit,*) "Elem_op%Dense_val is allocated. The dense operator action called."
      CALL Action_dense(Op_psi=Op_psi, Elem_op=Elem_op, Psi=Psi, Verbose=Verbose_local, Debug=Debug_local)

    ELSE
      WRITE(out_unit,*) "### None of this operator's matrices are allocated. Please check its initalization."
      STOP "### None of this operator's matrices are allocated. Please check its initalization."
    END IF

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "--- Resulting statevector from the action of the HO1D Elem_op on the Psi statevector operand, computed &
                        &by Action_HO1D_operator_R1 :"
      CALL Write_Vec(Op_psi, out_unit, 1, info="Op_Psi")
      WRITE(out_unit,*) "--- End resulting statevector computed by Action_HO1D_operator_R1"
    END IF
  
    IF (Verbose_local > 25) WRITE(out_unit,*) 
    IF (Verbose_local > 25) WRITE(out_unit,*) "----------------------------------------ACTION OF THE HO1D OPERATOR OVER THE R1 WF&
                                              & COMPUTED---------------------------------------"; FLUSH(out_unit)
  
  END SUBROUTINE MolecCav_Action_elem_op_R1_real

  
  SUBROUTINE MolecCav_Action_dense_elem_op_R1_real(Op_psi, Elem_op, Psi, Verbose, Debug)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    IMPLICIT NONE
    
    real(kind=Rkind),      intent(inout) :: Op_psi(:)
    TYPE(Elem_op_t),       intent(in)    :: Elem_op
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
    Op_psi(:) = matmul(Elem_op%Dense_val, Psi)

  END SUBROUTINE MolecCav_Action_dense_elem_op_R1_real

  
  SUBROUTINE MolecCav_Action_diag_elem_op_R1_real(Op_psi, Elem_op, Psi, Verbose, Debug) 
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    IMPLICIT NONE
    
    real(kind=Rkind),      intent(inout) :: Op_psi(:)
    TYPE(Elem_op_t),       intent(in)    :: Elem_op
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
    Op_psi = Elem_op%Diag_val * Psi

  END SUBROUTINE MolecCav_Action_diag_elem_op_R1_real

  
  SUBROUTINE MolecCav_Action_band_elem_op_R1_real(Op_psi, Elem_op, Psi, Verbose, Debug)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    IMPLICIT NONE
    
    real(kind=Rkind),      intent(inout) :: Op_psi(:)
    TYPE(Elem_op_t),       intent(in)    :: Elem_op
    real(kind=Rkind),      intent(in)    :: Psi(:)
    integer, optional,     intent(in)    :: Verbose                                                                              ! cf. comments in HO1D_parameters_m
    logical, optional,     intent(in)    :: Debug                                                                                ! cf. comments in HO1D_parameters_m

    integer                              :: i, Nb
    integer                              :: Verbose_local                                                                   ! goes from 25 (= 0 verbose) to 29 (= maximum verbose) at this layer
    logical                              :: Debug_local

    !------------------------------------------------------Debugging options-----------------------------------------------------
    IF (PRESENT(Verbose)) THEN; Verbose_local = Verbose
    ELSE; Verbose_local = 20; END IF 
    IF (PRESENT(Debug))   THEN; Debug_local   = Debug
    ELSE; Debug_local = .FALSE.; END IF

    IF (Verbose_local > 27) WRITE(out_unit,*) 
    IF (Verbose_local > 27) WRITE(out_unit,*) "---------------------------------------Computing the HO1D operator using the band &
                                             &procedure--------------------------------------"; FLUSH(out_unit)
    
    !----------------------------------------------------Computing the action----------------------------------------------------
    Nb = size(Op_psi)

    Op_psi     = ZERO
    Op_psi     = Elem_op%Band_val(:,2) * Psi
    Op_psi(1)  = Op_psi(1)  + Elem_op%Band_val(2,3)    * Psi(2)
    Op_psi(Nb) = Op_psi(Nb) + Elem_op%Band_val(Nb-1,1) * Psi(Nb-1)
    DO i = 2, Nb-1
      Op_psi(i) = Op_psi(i) + &
                & Elem_op%Band_val(i-1,1) * Psi(i-1) + &
                & Elem_op%Band_val(i+1,3) * Psi(i+1)
    END DO

  END SUBROUTINE MolecCav_Action_band_elem_op_R1_real


  SUBROUTINE MolecCav_Action_elem_op_R1_complex(Op_psi, Elem_op, Psi, Verbose, Debug)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    IMPLICIT NONE

    complex(kind=Rkind),   intent(inout) :: Op_psi(:)
    TYPE(Elem_op_t),       intent(in)    :: Elem_op
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
      WRITE(out_unit,*) "--- Arguments of MolecCav_Action_elem_op_R1_complex :"
      WRITE(out_unit,*) "The <<Elem_op>> argument :"
      CALL Write(Elem_op)
      WRITE(out_unit,*) "The <<Psi>> argument : "
      CALL Write_Vec(Psi, out_unit, 1, info="Psi")
      WRITE(out_unit,*) "The size of its vector : "//TO_string(Size(Psi))
      WRITE(out_unit,*) "--- End arguments of MolecCav_Action_elem_op_R1_complex"
      FLUSH(out_unit)
    END IF
    
    !-----------------------------------------------------Checking dimensions----------------------------------------------------
    IF (ALLOCATED(Elem_op%Dense_val)) Nb = Size(Elem_op%Dense_val, dim=1)
    IF (ALLOCATED(Elem_op%Diag_val) ) Nb = Size(Elem_op%Diag_val,  dim=1)
    IF (ALLOCATED(Elem_op%Band_val) ) Nb = Size(Elem_op%Band_val,  dim=1)

    IF (Nb /= Size(Psi)) THEN
      WRITE(out_unit,*) "### The dimensions of the Elem_op's matrix representation does not match the operand Psi's vector size.&
                       & Please check initialization."
      WRITE(out_unit,*) "    Size(Elem_op%matrix, dim=1) = "//TO_string(Nb)//"; Size(Psi) = "//TO_string(Size(Psi))
      STOP "### The dimensions of the Elem_op's matrix representation does not match the operand Psi's vector size. Please check& 
                       & initialization."
    END IF 

    IF (Nb /= Size(Op_psi)) THEN
      WRITE(out_unit,*) "### The dimensions of the Elem_op's matrix representation does not match the resulting Op_psi vector's & 
                       &size. Please check initialization."
      WRITE(out_unit,*) "    Size(Elem_op%matrix, dim=1) = "//TO_string(Nb)//"; Size(Op_psi) = "//TO_string(Size(Op_psi))
      STOP "### The dimensions of the Elem_op's matrix representation does not match the resulting Op_psi vector's size. Please & 
                       &check initialization."
    END IF 

    !---------------------------------------------Selection of the calculation method--------------------------------------------
    IF      (ALLOCATED(Elem_op%Diag_val))   THEN
      IF (Debug_local) WRITE(out_unit,*) "--- Elem_op%Diag_val is allocated. The diagonal operator action called."
      CALL Action_diag(Op_psi=Op_psi, Elem_op=Elem_op, Psi=Psi, Verbose=Verbose_local, Debug=Debug_local)

    ELSE IF (ALLOCATED(Elem_op%Band_val))   THEN
      IF (Debug_local) WRITE(out_unit,*) "--- Elem_op%Diag_val is not allocated."
      IF (Debug_local) WRITE(out_unit,*) "--- Elem_op%Band_val is allocated. The band operator action called."
      CALL Action_band(Op_psi=Op_psi, Elem_op=Elem_op, Psi=Psi, Verbose=Verbose_local, Debug=Debug_local)

    ELSE IF (ALLOCATED(Elem_op%Dense_val)) THEN
      IF (Debug_local) WRITE(out_unit,*) "Elem_op%Diag_val and Elem_op%Band_val are not allocated."
      IF (Debug_local) WRITE(out_unit,*) "Elem_op%Dense_val is allocated. The dense operator action called."
      CALL Action_dense(Op_psi=Op_psi, Elem_op=Elem_op, Psi=Psi, Verbose=Verbose_local, Debug=Debug_local)

    ELSE
      WRITE(out_unit,*) "### None of this operator's matrices are allocated. Please check its initalization."
      STOP "### None of this operator's matrices are allocated. Please check its initalization."
    END IF

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "--- Resulting statevector from the action of the HO1D Elem_op on the Psi statevector operand, computed &
                        &by Action_HO1D_operator_R1 :"
      CALL Write_Vec(Op_psi, out_unit, 1, info="Op_Psi")
      WRITE(out_unit,*) "--- End resulting statevector computed by Action_HO1D_operator_R1"
    END IF
  
    IF (Verbose_local > 25) WRITE(out_unit,*) 
    IF (Verbose_local > 25) WRITE(out_unit,*) "----------------------------------------ACTION OF THE HO1D OPERATOR OVER THE R1 WF&
                                              & COMPUTED---------------------------------------"; FLUSH(out_unit)
  
  END SUBROUTINE MolecCav_Action_elem_op_R1_complex

  
  SUBROUTINE MolecCav_Action_dense_elem_op_R1_complex(Op_psi, Elem_op, Psi, Verbose, Debug)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    IMPLICIT NONE
    
    complex(kind=Rkind),   intent(inout) :: Op_psi(:)
    TYPE(Elem_op_t),       intent(in)    :: Elem_op
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
    Op_psi(:) = matmul(Elem_op%Dense_val, Psi)

  END SUBROUTINE MolecCav_Action_dense_elem_op_R1_complex

  
  SUBROUTINE MolecCav_Action_diag_elem_op_R1_complex(Op_psi, Elem_op, Psi, Verbose, Debug) 
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    IMPLICIT NONE
    
    complex(kind=Rkind),   intent(inout) :: Op_psi(:)
    TYPE(Elem_op_t),       intent(in)    :: Elem_op
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
    Op_psi = Elem_op%Diag_val * Psi

  END SUBROUTINE MolecCav_Action_diag_elem_op_R1_complex

  
  SUBROUTINE MolecCav_Action_band_elem_op_R1_complex(Op_psi, Elem_op, Psi, Verbose, Debug)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    IMPLICIT NONE
    
    complex(kind=Rkind),   intent(inout) :: Op_psi(:)
    TYPE(Elem_op_t),       intent(in)    :: Elem_op
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
    Op_psi     = Elem_op%Band_val(:,2) * Psi
    Op_psi(1)  = Op_psi(1)  + Elem_op%Band_val(2,3)    * Psi(2)
    Op_psi(Nb) = Op_psi(Nb) + Elem_op%Band_val(Nb-1,1) * Psi(Nb-1)
    DO i = 2, Nb-1
      Op_psi(i) = Op_psi(i) + &
                & Elem_op%Band_val(i-1,1) * Psi(i-1) + &
                & Elem_op%Band_val(i+1,3) * Psi(i+1)
    END DO

  END SUBROUTINE MolecCav_Action_band_elem_op_R1_complex
  

  SUBROUTINE MolecCav_Write_elem_op_R1(Elem_op)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
    USE QDUtil_m
    IMPLICIT NONE 

    TYPE(Elem_op_t), intent(in) :: Elem_op

    IF (ALLOCATED(Elem_op%Operator_type)) THEN
      WRITE(out_unit,*) "_______________________________________the operator parameters______________________________________"
      WRITE(out_unit,*) "|The operator's nature <<Operator_type>> information do is allocated, and is   | ", Elem_op%Operator_type
      WRITE(out_unit,*) "|______________________________________________________________________________|____________________"
      FLUSH(out_unit)
    ELSE 
      WRITE(out_unit,*) "_____________________________the operator parameters____________________________"
      WRITE(out_unit,*) "|The operator's nature <<Operator_type>> information is NOT allocated.         |"
      WRITE(out_unit,*) "|______________________________________________________________________________|____________________"
      FLUSH(out_unit)
    END IF

    WRITE(out_unit,*) "|Is its matrix supposed to be represented as a dense one ? Elem_op%Dense       | "//TO_string(Elem_op%Dense)
    WRITE(out_unit,*) "|______________________________________________________________________________|____________________"
    
    WRITE(out_unit,*) "|Case of a band matrix, (Elem_op%Upper_bandwidth, Elem_op%Lower_bandwidth)     | ("//TO_string(Elem_op%Up&
                      &per_bandwidth)//","//TO_string(Elem_op%Lower_bandwidth)//")"
    WRITE(out_unit,*) "|______________________________________________________________________________|____________________"
    FLUSH(out_unit)

    WRITE(out_unit,*) "_________________________the operator's representations_____________________________________________"
    WRITE(out_unit,*) "|Is its matrix supposed to be represented on the grid ? Elem_op%Grid           | "//TO_string(Elem_op%Grid)
    WRITE(out_unit,*) "|______________________________________________________________________________|____________________"

    IF (ALLOCATED(Elem_op%Diag_val)) THEN                                                                                     ! we assume that the code is supposed to be used only allocating one of the matrices of each Elem_op_t object
      WRITE(out_unit,*) "|The operator is represented using a matrix of size (Elem_op%Nb) :             | ", SIZE(Elem_op%Diag_val)
      WRITE(out_unit,*) "|______________________________________________________________________________|____________________"
      FLUSH(out_unit)
      
      WRITE(out_unit,*) "|The operator's Diagonal matrix representation has been used, and is           |"
      WRITE(out_unit,*) "|______________________________________________________________________________|"
      CALL Write_Vec(Elem_op%Diag_val, out_unit, 1, info="Diag_val")
      WRITE(out_unit,*) "|______________________________________________________________________________|"
      FLUSH(out_unit)

    ELSE 
      WRITE(out_unit,*) "|The operator's Diagonal matrix representation is NOT allocated.               |"
      WRITE(out_unit,*) "|______________________________________________________________________________|"
      FLUSH(out_unit)
    END IF

    IF (ALLOCATED(Elem_op%Band_val)) THEN
      WRITE(out_unit,*) "|The operator is represented using a matrix of size (Elem_op%Nb) :             | "//TO_string(SIZE(Elem_&
      &op%Band_val, 1))//" * "//TO_string(SIZE(Elem_op%Band_val, 2))
      WRITE(out_unit,*) "|______________________________________________________________________________|_____________________"
      FLUSH(out_unit)

      WRITE(out_unit,*) "|The operator's Band matrix representation has been used, and is               |"
      WRITE(out_unit,*) "|______________________________________________________________________________|"
      CALL Write_Mat(Elem_op%Band_val, out_unit, 3, info="Band_val")
      WRITE(out_unit,*) "|______________________________________________________________________________|"
      FLUSH(out_unit)
  
    ELSE 
      WRITE(out_unit,*) "|The operator's Band matrix representation is NOT allocated.                   |"
      WRITE(out_unit,*) "|______________________________________________________________________________|"
      FLUSH(out_unit)
    END IF

    IF (ALLOCATED(Elem_op%Dense_val)) THEN
      WRITE(out_unit,*) "|The operator is represented using a matrix of size (Elem_op%Nb) :             | "//TO_string(SIZE(Elem_&
      &op%Dense_val, 1))//" * "//TO_string(SIZE(Elem_op%Dense_val, 2))
      WRITE(out_unit,*) "|______________________________________________________________________________|_____________________"
      FLUSH(out_unit)

      WRITE(out_unit,*) "|The operator's Dense matrix representation has been used, and is              |"
      CALL Write_Mat(Elem_op%Dense_val, out_unit, Size(Elem_op%Dense_val, dim=2), info="Dense_val")
      WRITE(out_unit,*) "|______________________________________________________________________________|"
      FLUSH(out_unit)
  
    ELSE 
      WRITE(out_unit,*) "|The operator's Dense matrix representation is NOT allocated.                  |"
      WRITE(out_unit,*) "|__________________________End HO1D operator object____________________________|"
      FLUSH(out_unit)
    END IF
  
  END SUBROUTINE MolecCav_Write_elem_op_R1


  SUBROUTINE MolecCav_Deallocate_elem_op(Elem_op, Verbose, Debug)
    USE QDUtil_m
    IMPLICIT NONE 

    TYPE(Elem_op_t),       intent(inout) :: Elem_op
    integer, optional,     intent(in)    :: Verbose                                                                                 ! cf. comments in HO1D_parameters_m
    logical, optional,     intent(in)    :: Debug                                                                                   ! cf. comments in HO1D_parameters_m

    integer                              :: Verbose_local = 25                                                                      ! goes from 25 (= 0 verbose) to 29 (= maximum verbose) at this layer
    logical                              :: Debug_local   = .FALSE.

    !------------------------------------------------------Debugging options-----------------------------------------------------
    IF (PRESENT(Verbose)) Verbose_local = Verbose
    IF (PRESENT(Debug))   Debug_local   = Debug

    IF (Verbose_local > 27) WRITE(out_unit,*)
    IF (Verbose_local > 27) WRITE(out_unit,*) "-----------------------------------------------Deallocating the HO1D_operator obje&
                                              &ct----------------------------------------------"
    IF (Debug_local) THEN
      WRITE(out_unit,*) "--- The HO1D operator to be deallocated :"
      CALL Write(Elem_op)
      WRITE(out_unit,*) "--- End HO1D operator to be deallocated"
    END IF 

    !-----------------------------Deallocating the HO1D operator object----------------------------
    IF (ALLOCATED(Elem_op%Operator_type)) DEALLOCATE(Elem_op%Operator_type)
    Elem_op%Dense           = .FALSE.
    Elem_op%Grid            = .FALSE.
    Elem_op%Upper_bandwidth = 0
    Elem_op%Lower_bandwidth = 0
    IF (ALLOCATED(Elem_op%Dense_val))     DEALLOCATE(Elem_op%Dense_val)
    IF (ALLOCATED(Elem_op%Diag_val))      DEALLOCATE(Elem_op%Diag_val)
    IF (ALLOCATED(Elem_op%Band_val))      DEALLOCATE(Elem_op%Band_val)

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "--- The HO1D operator object after having been deallocated :"
      CALL Write(Elem_op)
      WRITE(out_unit,*) "--- End dellocated HO1D operator"
    END IF

  END SUBROUTINE MolecCav_Deallocate_elem_op


END MODULE
