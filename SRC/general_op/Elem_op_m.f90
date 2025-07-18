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
! Module designed to manage elementary operators of quantum dynamics. They are considered as opera-
! tors because of the derived type this module takes care of but it actually does not care about p-
! hysics at all. It just manages the computation a tensor's action upon a rank-1 real tensor, choo-
! sing a different method according to the kind of object the "operator" is represented by. This w-
! ay allows to deal with operator's actions in way totally independant of the basis they are expre-
! ssed on, which gives more generality and flexibility. However, because of this choice it only co-
! ntains, in addition to the derived type Elem_op it takes care of, procedures to process Actions, 
! Writing, and Deallocation with this type. Obviouly, from the physical point of view, in all of t-
! he hereiafter action procedure, the wavefunctions are assumed to be expressed in the same basis 
! set used to expressed the elementary operator.
! This module is in the "-5" level. At this level, Verbose is ranging from 21 (degree 0) to 24 (de-
! gree 4).
!
! Elem_op_t : The derived typed used to represent a mono-dimensional operator of quantum mechanics.
! It contains only informations about its nature and about the way the values of its tensor repres-
! entation are actually stored in memory, but nothing about the context of this tensor representat-
! ion. This type is totally unrelated to the basis set choice and just manage values and tables, m-
! ostly. Note that, confusingly enough, there is no procedures within this module to initialize an
! object of this type. This is consistent with the spirit of a module dedicated to the management 
! of operators independantly on the basis set choice done to express their matrices (i.e. managing-
! only their actions assuming their matrices have been provided with some information about their 
! "encoding").  The corresponding initialization procedures are implemented in Quantum_HO1D, which
! takes care of the 1D HO basis. If any other module is implemented to manage another kind of basi-
! s, it should hold inside its own procedures to initialize Elem_op according to its basis.
!   - Operator type : the name of the operator (string). This parameter is not used at all in this
! module since the action depends only on the shape given to the tensor and not on its physical me-
! aning, but is useful to distinguish operators in more external modules. In order to recognised b-
! y the code, please use the following nomenclature :"Hamiltonian", "Position", "NbQuanta", or "Di-
! pMomt". Note that this is NOT case-sensitive.
!   - Dense : a logical parameter that indicated whether the storage of the tensor representation's
! values of the operator has been optimised or not. In other words, if it is set to .TRUE., the st-
! orage will not be optimized and it will be stored as a "Dense" matrix i.e. storing the full/comp-
! lete matrix of the operator on the chosed basis, that is in the Dense_val parameter. This case i-
! mplies that the Dense_val parameter is allocated and should have been initialised. Otherwise, if
! Dense == .FALSE., the matrix' elements are stored in a table (Diag_val or Band_val) such as taki-
! ng advantage of the sparcity of the operator's analytical matrix in the basis. This case implies 
! that the Diag_val or Band_val parameter is allocated and should have been initialised. The defau-
! lt case is the optimised representation.
!   - Grid : a logical parameter that indicate whether the thensor representation is expressed in a 
! discretised spacial representation (Grid == .TRUE.) or in a basis set decomposition (.FALSE.). H-
! owever the grid representation has not been implemented yet and thus is not available.
!   - Upper_bandwidth (resp. Lower_bandwidth) : An integer used only if Dense == .FALSE. and if th-
! e Band_val parameter is used, which corresponds to the situation where the "optimal" way of stor-
! age is used for an operator with non-diagonal band matrix on the used basis. This parameter indi-
! cates the number of additional bands to consider above (resp. below) the diagonal. For example, 
! Upper_bandwidth=Lower_bandwidth=1 would give a tridiagonal matrix. These parameters allows more 
! flexibility to the code if other operators/basis set are implemented, but, for now, only the tri-
! diagonal case is used.
!   - Dense_val : The table table that is expected to be allocated and to have been initialised if 
! Dense == .TRUE. . It is supposed to be the rank-2 tensor representation associated to the operat-
! or on the chosen representation (basis set or space grid), namely its matrix in full.
!   - Diag_val : The table table that is expected to be allocated and to have been initialised if 
! Dense == .FALSE. and if the rank-2 tensor of the operator on the chosen representation is diagon-
! al. In this case, Diag_val is designed to be the rank-1 tensor that holds the diagonal elements 
! of the rank-2 tensor, in a way that Diag_val(i) = Dense_val(i,i).
!   - Band_val : The table table that is expected to be allocated and to have been initialised if 
! Dense == .FALSE. and if the rank-2 tensor of the operator on the chosen representation is band d-
! iagonal. In this case, Band_val is designed to be the collection of rank-1 tensors that holds th-
! e diagonal elements of each non-null diagonal from the rank-2 tensor. If this parameter is used, 
! Upper/Lower_bandwidth should not both be 0. Then, Band_val is supposed to have been initialised 
! in such a way that the N first/leftmost columns (N=Lower_bandwidth) are the bands below the diag-
! onal - from the lowest to the diagonal -, the N+1 columns is the diagonal, and the (SIZE(Band_va-
! l, dim=2)-N-1 = Upper_bandwidth) rightmost columns are the bands above the diagonal - from the d-
! iagonal to the highest. Obviously the non-diagonal bands are shorter than the diagonal. Therefor-
! e the rank-1 tensors that represent the bands below (resp. above) the diagonal are filled with z-
! eros at the bottom (resp. top). Since the case Upper_bandwidth = Lower_bandwidth = 1 reduces its-
! elf to the "diagonal" case, the "band" case may be seen as an extension of the latter, just cons-
! idering more bands i.e. adding more columns around Diag_val.
! As a final comment about the use of the Elem_op derived type, please note that it is not designe-
! d to hold different representations at once. In other words, do NEVER initialise the Dense_val a-
! nd another table parameter in the same Elem_op object. Even if it is numerically possible, it wo-
! uld not be good practice, and would lead to unforseen situations and likely unwanted behaviour. 
!
! Action : Interface for the MolecCav_Action_elem_op_R1_* and procedure. Takes as arguments a tens-
! or of rank 1 with real or complex values (Op_psi) that does not need to have been initialised, a-
! n object of derived type Elem_op_t (Elem_op), and an other tensor of rank 1 with real or complex 
! values (Psi). It computes the action of the allocated tensor of Elem_op upon the Psi one by redi-
! recting to one of the associated procedure according to the arguments. It computes the action by 
! calling the Action_* procedure according to which *_val parameter of Elem_op is allocated, and a-
! ffects the result to Op_psi. From a physical point of view, this is intended in the code to repr-
! esent the computation of the 1D wavefunction Op_psi by taking the action of the 1D quantum mecha-
! nics operator upon the 1D wavefunction Psi.
!
! Write : Interface for the MolecCav_Write_elem_op_R1 procedure. cf. MolecCav_Write_elem_op_R1 
! for details.
!
! Dealloc : Interface for the MolecCav_Deallocate_elem_op procedure. cf. MolecCav_Deallocate_e-
! lem_op for details.
!
! MolecCav_Action_elem_op_R1_real : Takes as arguments a tensor of rank 1 with real values (Op_psi)
! that does not need to have been initialised, an object of derived type Elem_op_t (Elem_op), and 
! an other tensor of rank 1 with real values (Psi). It computes the action of the allocated tensor 
! of Elem_op upon the Psi one by calling the Action_* procedure according to which *_val parameter 
! of Elem_op is allocated, and affects the result to Op_psi. From a physical point of view, this is 
! intended in the code to represent the computation of the 1D wavefunction Op_psi by taking the ac-
! tion of the 1D quantum mechanics operator upon the 1D wavefunction Psi.
!
! MolecCav_Action_elem_op_R1_complex : Same as MolecCav_Action_elem_op_R1_real but for two tensors 
! with complex values.
!
! Action_dense : Interface for the MolecCav_Action_dense_elem_op_R1_* procedures. Takes as argumen-
! ts a tensor of rank 1 with real or complex values (Op_psi) that does not need to have been initi-
! alised, an object of derived type Elem_op_t (Elem_op) whom Dense_val is assumed to be allocated 
! and to have been initialized, and an other tensor of rank 1 with real or complex values (Psi). I-
! t computes the action of the Dense_val parameter of Elem_op upon Psi by by redirecting to one of 
! the associated procedure according to the arguments. It computes the action of the Dense_val par-
! ameter of Elem_op upon Psi by a matrix-vector multiplication, and affects the result to Op_psi.
!
! Action_diag : Interface for the MolecCav_Action_diag_elem_op_R1_* procedures. Takes as arguments 
! a tensor of rank 1 with real values (Op_psi) that does not need to have been initialised, an obj-
! ect of derived type Elem_op_t (Elem_op) whom Diag_val is assumed to be allocated and to have bee-
! n initialized, and an other tensor of rank 1 with real values (Psi). It computes the action of t-
! he Diag_val parameter of Elem_op upon Psi by redirecting to one of the associated procedure acco-
! rding to the arguments. It computes the action by multiplicating two-by-two their respective ele-
! ments with same index, and affects the result to Op_psi.
!
! Action_band : Interface for the MolecCav_Action_band_elem_op_R1_* procedures. Takes as arguments 
! a tensor of rank 1 with real values (Op_psi) that does not need to have been initialised, an obj-
! ect of derived type Elem_op_t (Elem_op) whom Diag_val is assumed to be allocated and to have bee-
! n initialized, and an other tensor of rank 1 with real values (Psi). It computes the action of t-
! he Diag_val parameter of Elem_op upon Psi by redirecting to one of the associated procedure acco-
! rding to the arguments. It computes the action by an element-by-element multiplication that can 
! be demonstrated analytically by writing, and affects the result to Op_psi. However, only the cas-
! e of a Band_val corresponding to a triband matrix is implemented yet and not the general computa-
! tion.
!
! MolecCav_Action_dense_elem_op_R1_real : Takes as arguments a tensor of rank 1 with real values (-
! Op_psi) that does not need to have been initialised, an object of derived type Elem_op_t (Elem_o-
! p) whom Dense_val is assumed to be allocated and to have been initialized, and an other tensor o-
! f rank 1 with real values (Psi). It computes the action of the Dense_val parameter of Elem_op up-
! on Psi by a matrix-vector multiplication, and affects the result to Op_psi.
!
! MolecCav_Action_dense_elem_op_R1_complex : Same as MolecCav_Action_dense_elem_op_R1_real but for
! two tensors with complex values.
!
! MolecCav_Action_diag_elem_op_R1_real : Takes as arguments a tensor of rank 1 with real values (-
! Op_psi) that does not need to have been initialised, an object of derived type Elem_op_t (Elem_o-
! p) whom Diag_val is assumed to be allocated and to have been initialized, and an other tensor of
! rank 1 with real values (Psi). It computes the action of the Diag_val parameter of Elem_op upon 
! Psi by multiplicating two-by-two their respective elements with same index, and affects the resu-
! lt to Op_psi.
!
! MolecCav_Action_diag_elem_op_R1_complex : Same as MolecCav_Action_diag_elem_op_R1_real but for t-
! wo tensors with complex values.
!
! MolecCav_Action_band_elem_op_R1_real : Takes as arguments a tensor of rank 1 with real values (-
! Op_psi) that does not need to have been initialised, an object of derived type Elem_op_t (Elem_o-
! p) whom Band_val is assumed to be allocated and to have been initialized, and an other tensor of
! rank 1 with real values (Psi). It computes the action of the Band_val parameter of Elem_op upon 
! Psi by an element-by-element multiplication that can be demonstrated analytically by writing, and
! affects the result to Op_psi. However, only the case of a Band_val corresponding to a triband ma-
! trix is implemented yet and not the general computation.
!
! MolecCav_Action_band_elem_op_R1_complex : Same as MolecCav_Action_band_elem_op_R1_real but for t-
! wo tensors with complex values.
!
! MolecCav_Write_elem_op_R1 : Takes as argument an object of derived type Elem_op_t (Elem_op) and 
! an optional string messsage (Info), and write in the standard output all of the Elem_op paramete-
! rs, below the Info if provided.
!
! MolecCav_Deallocate_elem_op : Takes as argument an object of derived type Elem_op_t (Elem_op) and
! reset it by deallocating all tables/strings and setting to defaut values the other parameters.
!
!==================================================================================================
!==================================================================================================
MODULE Elem_op_m
  !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
  USE QDUtil_m
  IMPLICIT NONE


  TYPE                               :: Elem_op_t
    character(len=:),    allocatable :: Operator_type                                                   ! ex : "Hamiltonian", "Position", "NbQuanta", etc
    logical                          :: Dense = .FALSE.                                                 ! if .TRUE. then the matrix storage will not be optimized and it will be stored as a Dense matrix. Otherwise, the elements matrix will be stored in a table such as taking advantage of the sparcity of the operator's analytical matrix in the HO Eigenbasis
    logical                          :: Grid  = .FALSE.                                                 ! if .TRUE. then the matrix storage will not be optimized and it will be stored as a Dense matrix. Otherwise, the elements matrix will be stored in a table such as taking advantage of the sparcity of the operator's analytical matrix in the HO Eigenbasis
    integer                          :: Upper_bandwidth = 0                                             ! if type = "Band". Gives the number of additional bands to consider above the diagonal.
    integer                          :: Lower_bandwidth = 0                                             ! if type = "Band". Gives the number of additional bands to consider below the diagonal. Ex : Upper_bandwidth=Lower_bandwidth=1 would give a tridiagonal matrix
    real(kind=Rkind),    allocatable :: Dense_val(:,:)                                                  ! if Dense == .TRUE.
    real(kind=Rkind),    allocatable :: Diag_val(:)                                                     ! if Dense == .FALSE. .AND. the operator analytical matrix in the HO1D Eigenbasis is diagonal. The diagonal elements are stored in a vector (rank-1 tensor)
    real(kind=Rkind),    allocatable :: Band_val(:,:)                                                   ! if Dense == .FALSE. .AND. the operator analytical matrix in the HO1D Eigenbasis is band. The number of columns will be the number of diagonals to consider : each considered diagonal is stored in a column
  END TYPE


  PRIVATE

  PUBLIC Elem_op_t, Action, Write, Dealloc


  INTERFACE Action
    MODULE PROCEDURE MolecCav_Action_elem_op_R1_real, MolecCav_Action_elem_op_R1_complex
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
    integer, optional,     intent(in)    :: Verbose                                                     ! cf. comments in HO1D_parameters_m
    logical, optional,     intent(in)    :: Debug                                                       ! cf. comments in HO1D_parameters_m

    integer                              :: Nb
    integer                              :: Verbose_local
    logical                              :: Debug_local

    !--- Debugging options --------------------------------
    IF (PRESENT(Verbose)) THEN; Verbose_local = Verbose
    ELSE; Verbose_local = 21; END IF 
    IF (PRESENT(Debug)) THEN; Debug_local = Debug
    ELSE; Debug_local = .FALSE.; END IF
    IF (Debug_local) Verbose_local = 24

    IF (Verbose_local > 21) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "o o Arguments of MolecCav_Action_elem_op_R1_real :"
      WRITE(out_unit,*) "The <<Elem_op>> argument :"
      CALL Write(Elem_op)
      WRITE(out_unit,*) "The <<Psi>> argument     : "
      CALL Write_Vec(Psi, out_unit, 1, info="Psi")
      WRITE(out_unit,*) "The size of its vector   : "//TO_string(Size(Psi))
      FLUSH(out_unit)
    END IF
    
    !--- Checking dimensions ------------------------------
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

    !--- Selection of the calculation method --------------
    IF      (ALLOCATED(Elem_op%Diag_val))   THEN
      IF (Debug_local) WRITE(out_unit,*) "--- Elem_op%Diag_val is allocated. The diagonal operator action called."
      CALL Action_diag(Op_psi=Op_psi, Elem_op=Elem_op, Psi=Psi, Verbose=Verbose_local, Debug=Debug_local)

    ELSE IF (ALLOCATED(Elem_op%Band_val))   THEN
      IF (Debug_local) WRITE(out_unit,*) "--- Elem_op%Diag_val is not allocated."
      IF (Debug_local) WRITE(out_unit,*) "--- Elem_op%Band_val is allocated. The band operator action called."
      CALL Action_band(Op_psi=Op_psi, Elem_op=Elem_op, Psi=Psi, Verbose=Verbose_local, Debug=Debug_local)

    ELSE IF (ALLOCATED(Elem_op%Dense_val)) THEN
      IF (Debug_local) WRITE(out_unit,*) "--- Elem_op%Diag_val and Elem_op%Band_val are not allocated."
      IF (Debug_local) WRITE(out_unit,*) "--- Elem_op%Dense_val is allocated. The dense operator action called."
      CALL Action_dense(Op_psi=Op_psi, Elem_op=Elem_op, Psi=Psi, Verbose=Verbose_local, Debug=Debug_local)

    ELSE
      WRITE(out_unit,*) "### None of this operator's matrices are allocated. Please check its initalization."
      STOP "### None of this operator's matrices are allocated. Please check its initalization."
    END IF

    IF (Debug_local) THEN
      WRITE(out_unit,*) "--- Resulting statevector from the action of the HO1D Elem_op on the Psi statevector operand, computed &
                        &by Action_HO1D_operator_R1 :"
      CALL Write_Vec(Op_psi, out_unit, 1, info="Op_Psi")
    END IF
    
  END SUBROUTINE MolecCav_Action_elem_op_R1_real

  
  SUBROUTINE MolecCav_Action_dense_elem_op_R1_real(Op_psi, Elem_op, Psi, Verbose, Debug)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    IMPLICIT NONE
    
    real(kind=Rkind),      intent(inout) :: Op_psi(:)
    TYPE(Elem_op_t),       intent(in)    :: Elem_op
    real(kind=Rkind),      intent(in)    :: Psi(:)
    integer, optional,     intent(in)    :: Verbose                                                     ! cf. comments in HO1D_parameters_m
    logical, optional,     intent(in)    :: Debug                                                       ! cf. comments in HO1D_parameters_m

    integer                              :: Verbose_local
    logical                              :: Debug_local

    !--- Debugging options --------------------------------
    IF (PRESENT(Verbose)) THEN; Verbose_local = Verbose
    ELSE; Verbose_local = 21; END IF 
    IF (PRESENT(Debug)) THEN; Debug_local = Debug
    ELSE; Debug_local = .FALSE.; END IF
    IF (Debug_local) Verbose_local = 24

    IF (Verbose_local > 21) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "o o Arguments of MolecCav_Action_dense_elem_op_R1_real :"
      WRITE(out_unit,*) "The <<Elem_op>> argument :"
      CALL Write(Elem_op)
      WRITE(out_unit,*) "The <<Psi>> argument     : "
      CALL Write_Vec(Psi, out_unit, 1, info="Psi")
      WRITE(out_unit,*) "The size of its vector   : "//TO_string(Size(Psi))
      FLUSH(out_unit)
    END IF
    
    !--- Computing the action -----------------------------
    Op_psi(:) = matmul(Elem_op%Dense_val, Psi)

  END SUBROUTINE MolecCav_Action_dense_elem_op_R1_real

  
  SUBROUTINE MolecCav_Action_diag_elem_op_R1_real(Op_psi, Elem_op, Psi, Verbose, Debug) 
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    IMPLICIT NONE
    
    real(kind=Rkind),      intent(inout) :: Op_psi(:)
    TYPE(Elem_op_t),       intent(in)    :: Elem_op
    real(kind=Rkind),      intent(in)    :: Psi(:)
    integer, optional,     intent(in)    :: Verbose                                                     ! cf. comments in HO1D_parameters_m
    logical, optional,     intent(in)    :: Debug                                                       ! cf. comments in HO1D_parameters_m

    integer                              :: Verbose_local
    logical                              :: Debug_local

    !--- Debugging options --------------------------------
    IF (PRESENT(Verbose)) THEN; Verbose_local = Verbose
    ELSE; Verbose_local = 21; END IF 
    IF (PRESENT(Debug)) THEN; Debug_local = Debug
    ELSE; Debug_local = .FALSE.; END IF
    IF (Debug_local) Verbose_local = 24

    IF (Verbose_local > 21) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "o o Arguments of MolecCav_Action_diag_elem_op_R1_real :"
      WRITE(out_unit,*) "The <<Elem_op>> argument :"
      CALL Write(Elem_op)
      WRITE(out_unit,*) "The <<Psi>> argument     : "
      CALL Write_Vec(Psi, out_unit, 1, info="Psi")
      WRITE(out_unit,*) "The size of its vector   : "//TO_string(Size(Psi))
      FLUSH(out_unit)
    END IF
    
    !--- Computing the action -----------------------------
    Op_psi = Elem_op%Diag_val * Psi

  END SUBROUTINE MolecCav_Action_diag_elem_op_R1_real

  
  SUBROUTINE MolecCav_Action_band_elem_op_R1_real(Op_psi, Elem_op, Psi, Verbose, Debug)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    IMPLICIT NONE
    
    real(kind=Rkind),      intent(inout) :: Op_psi(:)
    TYPE(Elem_op_t),       intent(in)    :: Elem_op
    real(kind=Rkind),      intent(in)    :: Psi(:)
    integer, optional,     intent(in)    :: Verbose                                                     ! cf. comments in HO1D_parameters_m
    logical, optional,     intent(in)    :: Debug                                                       ! cf. comments in HO1D_parameters_m

    integer                              :: i, Nb
    integer                              :: Verbose_local
    logical                              :: Debug_local

    !--- Debugging options --------------------------------
    IF (PRESENT(Verbose)) THEN; Verbose_local = Verbose
    ELSE; Verbose_local = 21; END IF 
    IF (PRESENT(Debug)) THEN; Debug_local = Debug
    ELSE; Debug_local = .FALSE.; END IF
    IF (Debug_local) Verbose_local = 24

    IF (Verbose_local > 21) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "o o Arguments of MolecCav_Action_band_elem_op_R1_real :"
      WRITE(out_unit,*) "The <<Elem_op>> argument :"
      CALL Write(Elem_op)
      WRITE(out_unit,*) "The <<Psi>> argument     : "
      CALL Write_Vec(Psi, out_unit, 1, info="Psi")
      WRITE(out_unit,*) "The size of its vector   : "//TO_string(Size(Psi))
      FLUSH(out_unit)
    END IF

    !--- Computing the action -----------------------------
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
    integer, optional,     intent(in)    :: Verbose                                                     ! cf. comments in HO1D_parameters_m
    logical, optional,     intent(in)    :: Debug                                                       ! cf. comments in HO1D_parameters_m

    integer                              :: Nb
    integer                              :: Verbose_local
    logical                              :: Debug_local

    !--- Debugging options --------------------------------
    IF (PRESENT(Verbose)) THEN; Verbose_local = Verbose
    ELSE; Verbose_local = 21; END IF 
    IF (PRESENT(Debug)) THEN; Debug_local = Debug
    ELSE; Debug_local = .FALSE.; END IF
    IF (Debug_local) Verbose_local = 24

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "o o Arguments of MolecCav_Action_elem_op_R1_complex :"
      WRITE(out_unit,*) "The <<Elem_op>> argument :"
      CALL Write(Elem_op)
      WRITE(out_unit,*) "The <<Psi>> argument     : "
      CALL Write_Vec(Psi, out_unit, 1, info="Psi")
      WRITE(out_unit,*) "The size of its vector   : "//TO_string(Size(Psi))
      FLUSH(out_unit)
    END IF
    
    !--- Checking dimensions ------------------------------
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

    !--- Selection of the calculation method --------------
    IF      (ALLOCATED(Elem_op%Diag_val))   THEN
      IF (Debug_local) WRITE(out_unit,*) "--- Elem_op%Diag_val is allocated. The diagonal operator action called."
      CALL Action_diag(Op_psi=Op_psi, Elem_op=Elem_op, Psi=Psi, Verbose=Verbose_local, Debug=Debug_local)

    ELSE IF (ALLOCATED(Elem_op%Band_val))   THEN
      IF (Debug_local) WRITE(out_unit,*) "--- Elem_op%Diag_val is not allocated."
      IF (Debug_local) WRITE(out_unit,*) "--- Elem_op%Band_val is allocated. The band operator action called."
      CALL Action_band(Op_psi=Op_psi, Elem_op=Elem_op, Psi=Psi, Verbose=Verbose_local, Debug=Debug_local)

    ELSE IF (ALLOCATED(Elem_op%Dense_val)) THEN
      IF (Debug_local) WRITE(out_unit,*) "--- Elem_op%Diag_val and Elem_op%Band_val are not allocated."
      IF (Debug_local) WRITE(out_unit,*) "--- Elem_op%Dense_val is allocated. The dense operator action called."
      CALL Action_dense(Op_psi=Op_psi, Elem_op=Elem_op, Psi=Psi, Verbose=Verbose_local, Debug=Debug_local)

    ELSE
      WRITE(out_unit,*) "### None of this operator's matrices are allocated. Please check its initalization."
      STOP "### None of this operator's matrices are allocated. Please check its initalization."
    END IF

    IF (Debug_local) THEN
      WRITE(out_unit,*) "--- Resulting statevector from the action of the HO1D Elem_op on the Psi statevector operand, computed &
                        &by Action_HO1D_operator_R1 :"
      CALL Write_Vec(Op_psi, out_unit, 1, info="Op_Psi")
    END IF
    
  END SUBROUTINE MolecCav_Action_elem_op_R1_complex

  
  SUBROUTINE MolecCav_Action_dense_elem_op_R1_complex(Op_psi, Elem_op, Psi, Verbose, Debug)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    IMPLICIT NONE
    
    complex(kind=Rkind),   intent(inout) :: Op_psi(:)
    TYPE(Elem_op_t),       intent(in)    :: Elem_op
    complex(kind=Rkind),   intent(in)    :: Psi(:)
    integer, optional,     intent(in)    :: Verbose                                                     ! cf. comments in HO1D_parameters_m
    logical, optional,     intent(in)    :: Debug                                                       ! cf. comments in HO1D_parameters_m

    integer                              :: Verbose_local
    logical                              :: Debug_local

    !--- Debugging options --------------------------------
    IF (PRESENT(Verbose)) THEN; Verbose_local = Verbose
    ELSE; Verbose_local = 21; END IF 
    IF (PRESENT(Debug)) THEN; Debug_local = Debug
    ELSE; Debug_local = .FALSE.; END IF
    IF (Debug_local) Verbose_local = 24

    IF (Verbose_local > 21) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "o o Arguments of MolecCav_Action_dense_elem_op_R1_complex :"
      WRITE(out_unit,*) "The <<Elem_op>> argument :"
      CALL Write(Elem_op)
      WRITE(out_unit,*) "The <<Psi>> argument     : "
      CALL Write_Vec(Psi, out_unit, 1, info="Psi")
      WRITE(out_unit,*) "The size of its vector   : "//TO_string(Size(Psi))
      FLUSH(out_unit)
    END IF
    
    !--- Computing the action -----------------------------
    Op_psi(:) = matmul(Elem_op%Dense_val, Psi)

  END SUBROUTINE MolecCav_Action_dense_elem_op_R1_complex

  
  SUBROUTINE MolecCav_Action_diag_elem_op_R1_complex(Op_psi, Elem_op, Psi, Verbose, Debug) 
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    IMPLICIT NONE
    
    complex(kind=Rkind),   intent(inout) :: Op_psi(:)
    TYPE(Elem_op_t),       intent(in)    :: Elem_op
    complex(kind=Rkind),   intent(in)    :: Psi(:)
    integer, optional,     intent(in)    :: Verbose                                                     ! cf. comments in HO1D_parameters_m
    logical, optional,     intent(in)    :: Debug                                                       ! cf. comments in HO1D_parameters_m

    integer                              :: Verbose_local
    logical                              :: Debug_local

    !--- Debugging options --------------------------------
    IF (PRESENT(Verbose)) THEN; Verbose_local = Verbose
    ELSE; Verbose_local = 21; END IF 
    IF (PRESENT(Debug)) THEN; Debug_local = Debug
    ELSE; Debug_local = .FALSE.; END IF
    IF (Debug_local) Verbose_local = 24

    IF (Verbose_local > 21) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "o o Arguments of MolecCav_Action_diag_elem_op_R1_complex :"
      WRITE(out_unit,*) "The <<Elem_op>> argument :"
      CALL Write(Elem_op)
      WRITE(out_unit,*) "The <<Psi>> argument     : "
      CALL Write_Vec(Psi, out_unit, 1, info="Psi")
      WRITE(out_unit,*) "The size of its vector   : "//TO_string(Size(Psi))
      FLUSH(out_unit)
    END IF
    
    !--- Computing the action -----------------------------
    Op_psi = Elem_op%Diag_val * Psi

  END SUBROUTINE MolecCav_Action_diag_elem_op_R1_complex

  
  SUBROUTINE MolecCav_Action_band_elem_op_R1_complex(Op_psi, Elem_op, Psi, Verbose, Debug)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    IMPLICIT NONE
    
    complex(kind=Rkind),   intent(inout) :: Op_psi(:)
    TYPE(Elem_op_t),       intent(in)    :: Elem_op
    complex(kind=Rkind),   intent(in)    :: Psi(:)
    integer, optional,     intent(in)    :: Verbose                                                     ! cf. comments in HO1D_parameters_m
    logical, optional,     intent(in)    :: Debug                                                       ! cf. comments in HO1D_parameters_m

    integer                              :: i, Nb
    integer                              :: Verbose_local
    logical                              :: Debug_local

    !--- Debugging options --------------------------------
    IF (PRESENT(Verbose)) THEN; Verbose_local = Verbose
    ELSE; Verbose_local = 21; END IF 
    IF (PRESENT(Debug)) THEN; Debug_local = Debug
    ELSE; Debug_local = .FALSE.; END IF
    IF (Debug_local) Verbose_local = 24

    IF (Verbose_local > 21) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "o o Arguments of MolecCav_Action_band_elem_op_R1_complex :"
      WRITE(out_unit,*) "The <<Elem_op>> argument :"
      CALL Write(Elem_op)
      WRITE(out_unit,*) "The <<Psi>> argument     : "
      CALL Write_Vec(Psi, out_unit, 1, info="Psi")
      WRITE(out_unit,*) "The size of its vector   : "//TO_string(Size(Psi))
      FLUSH(out_unit)
    END IF
    
    !--- Computing the action -----------------------------
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
  

  SUBROUTINE MolecCav_Write_elem_op_R1(Elem_op, Info)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
    USE QDUtil_m
    IMPLICIT NONE 

    TYPE(Elem_op_t),            intent(in) :: Elem_op                                                   ! ex : "Hamiltonian", "Position", etc. (len=:) Expects to be allocatable, while (len=*) is dedicated to a procedure argument
    character(len=*), optional, intent(in) :: Info

    IF (PRESENT(Info)) THEN
      WRITE(out_unit,*) "--- Parameters associated to the object of derived type Elem_op_t ("//Info//") :"
    ELSE
      WRITE(out_unit,*) "--- Parameters associated to the object of derived type Elem_op_t :"
    END IF

    IF (ALLOCATED(Elem_op%Operator_type)) THEN
      WRITE(out_unit,*) "Operator_type         = "//Elem_op%Operator_type
    ELSE 
      WRITE(out_unit,*) "Operator_type is NOT allocated"
    END IF 

    WRITE(out_unit,*) "Dense                 = "//TO_string(Elem_op%Dense)
    WRITE(out_unit,*) "Upper,Lower_bandwidth = "//TO_string(Elem_op%Upper_bandwidth)//","//TO_string(Elem_op%Lower_bandwidth)
    WRITE(out_unit,*) "Grid                  = "//TO_string(Elem_op%Grid)

    IF (ALLOCATED(Elem_op%Dense_val)) THEN                                                              ! we assume that the code is supposed to be used only allocating one of the matrices of each Elem_op_t object
      WRITE(out_unit,*) "Dense_val IS allocated, and SIZE(Dense_val, dim=1), SIZE(Dense_val, dim=2) = "//TO_string(SIZE(Elem_op%D&
      &ense_val, dim=1))//","//TO_string(SIZE(Elem_op%Dense_val, dim=2))   
      CALL Write_Mat(Elem_op%Dense_val, out_unit, SIZE(Elem_op%Dense_val, dim=2), info="Dense_val")
    ELSE 
      WRITE(out_unit,*) "Dense_val is NOT allocated"      
    END IF
  
    IF (ALLOCATED(Elem_op%Diag_val)) THEN                                                               ! we assume that the code is supposed to be used only allocating one of the matrices of each Elem_op_t object
      WRITE(out_unit,*) "Diag_val IS allocated, and SIZE(Diag_val, dim=1) = "//TO_string(SIZE(Elem_op%Diag_val, dim=1))      
      CALL Write_Vec(Elem_op%Diag_val, out_unit, 1, info="Diag_val")
    ELSE 
      WRITE(out_unit,*) "Diag_val is NOT allocated"      
    END IF

    IF (ALLOCATED(Elem_op%Band_val)) THEN                                                               ! we assume that the code is supposed to be used only allocating one of the matrices of each Elem_op_t object
      WRITE(out_unit,*) "Band_val IS allocated, and SIZE(Band_val, dim=1), SIZE(Band_val, dim=2) = "//TO_string(SIZE(Elem_op%Band&
      &_val, dim=1))//","//TO_string(SIZE(Elem_op%Band_val, dim=2))   
      CALL Write_Mat(Elem_op%Band_val, out_unit, SIZE(Elem_op%Band_val, dim=2), info="Band_val")
    ELSE 
      WRITE(out_unit,*) "Band_val is NOT allocated"      
    END IF
    FLUSH(out_unit)
    
  END SUBROUTINE MolecCav_Write_elem_op_R1


  SUBROUTINE MolecCav_Deallocate_elem_op(Elem_op, Verbose, Debug)
    USE QDUtil_m
    IMPLICIT NONE 

    TYPE(Elem_op_t),       intent(inout) :: Elem_op
    integer, optional,     intent(in)    :: Verbose                                                     ! cf. comments in HO1D_parameters_m
    logical, optional,     intent(in)    :: Debug                                                       ! cf. comments in HO1D_parameters_m

    integer                              :: Verbose_local
    logical                              :: Debug_local

    !--- Debugging options --------------------------------
    IF (PRESENT(Verbose)) THEN; Verbose_local = Verbose
    ELSE; Verbose_local = 21; END IF 
    IF (PRESENT(Debug)) THEN; Debug_local = Debug
    ELSE; Debug_local = .FALSE.; END IF
    IF (Debug_local) Verbose_local = 24

    IF (Verbose_local > 21) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "o o Arguments of MolecCav_Deallocate_elem_op :"
      WRITE(out_unit,*) "The <<Elem_op>> argument :"
      CALL Write(Elem_op)
      FLUSH(out_unit)
    END IF

    !--- Deallocating the HO1D operator object ------------
    IF (ALLOCATED(Elem_op%Operator_type)) DEALLOCATE(Elem_op%Operator_type)
    Elem_op%Dense           = .FALSE.
    Elem_op%Grid            = .FALSE.
    Elem_op%Upper_bandwidth = 0
    Elem_op%Lower_bandwidth = 0
    IF (ALLOCATED(Elem_op%Dense_val))     DEALLOCATE(Elem_op%Dense_val)
    IF (ALLOCATED(Elem_op%Diag_val))      DEALLOCATE(Elem_op%Diag_val)
    IF (ALLOCATED(Elem_op%Band_val))      DEALLOCATE(Elem_op%Band_val)

    IF (Debug_local) THEN
      WRITE(out_unit,*) "--- The deallocated HO1D operator object :"
      CALL Write(Elem_op)
    END IF

  END SUBROUTINE MolecCav_Deallocate_elem_op


END MODULE
