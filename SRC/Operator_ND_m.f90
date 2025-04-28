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
! The only module related to general HO that the others modules will need to call in a "USE".  
! Initialize_operator_ND : reads the namelist and initialize the type, then constructs the operat-
! or using parameters of the HO1D_para object from the so called derived type.
! Append_operator_ND     : add an HO operator to a already initialized object of Operator_ND_t type.
! Write_operator_ND      : display values of the type in the output
! Deallocate_operator_ND : deallocate all tables of the type
! The module to initialize the HO by reading its parameters from the namelist.  
! Read_HO1D_parameters  : reads the namelist and initialize the type.
! Write_HO1D_parameters : displays values of the type in the output.
!==================================================================================================
!==================================================================================================
MODULE Operator_ND_m
  !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT, real64
  USE QDUtil_m                                                                 ! gives Rkind=real64; out_unit=OUTPUT_UNIT; INPUT_UNIT=in_unit; EYE=i and other numbers; TO_LOWERCASE; TO_UPPERCASE;... We thereby use ZERO instead of 0.0_real64
  USE Cavity_mode_m
  USE Matter_mode_m
  IMPLICIT NONE


  TYPE                 :: Operator_ND_t                                        ! N_mat, N_cav canNOT be part of the derived type because they are not specific of one operator_ND, but parameters of the whole system/calculation. All OpND will have the same. 
  integer, allocatable :: tab_indexes_mat_op(:)                                ! N_mat : nb of DOF of the matter subsystem = nb of vibrational modes
  integer, allocatable :: tab_indexes_cav_op(:)                                ! N_cav : nb of DOF of the cavity subsystem = nb of cavity modes
  END TYPE

  TYPE(Matter_mode_t), allocatable :: tab_mat_op(:)
  TYPE(Cavity_mode_new_t), allocatable :: tab_cav_op(:)

  PRIVATE

  PUBLIC Operator_ND_t, Initialize, Action, Write, Dealloc

  INTERFACE Initialize
    MODULE PROCEDURE MolecCav_Initialize_operator_ND
  END INTERFACE
  INTERFACE Initialize_tabs
    MODULE PROCEDURE MolecCav_Initialize_tab_operator
  END INTERFACE
  INTERFACE Action
    MODULE PROCEDURE MolecCav_Action_operator_ND_R1_real, MolecCav_Action_operator_ND_R1_complex
  END INTERFACE
  INTERFACE Write
    MODULE PROCEDURE MolecCav_Write_operator_ND
  END INTERFACE
  INTERFACE Dealloc
    MODULE PROCEDURE MolecCav_Deallocate_operator_ND
  END INTERFACE
    

  CONTAINS


  SUBROUTINE MolecCav_Initialize_operator_ND(OpND, N_mat, N_cav, Operators_on_mat, Operators_on_cav, nio, Dense, Verbose, Debug) 
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
    USE QDUtil_m
    USE Cavity_mode_m
    USE Matter_mode_m
    IMPLICIT NONE
  
    TYPE(Operator_ND_t), intent(inout) :: OpND
    integer,             intent(in)    :: N_mat
    integer,             intent(in)    :: N_cav
    character(len=*),    intent(in)    :: Operators_on_mat(:) ! syntax : 'hxI', not case sensitive, 'x' <=> \otimes
    character(len=*),    intent(in)    :: Operators_on_cav(:) ! syntax : 'h',   not case sensitive, 'x' <=> \otimes. This exemple means OpND = H_mat_1\otimesI_mat_2\otimesH_cav
    integer,             intent(in)    :: nio
    logical, optional,   intent(in)    :: Dense                                                                        ! cf. comments in HO1D_parameters_m
    integer, optional,   intent(in)    :: Verbose                                                                      ! cf. comments in HO1D_parameters_m
    logical, optional,   intent(in)    :: Debug                                                                        ! cf. comments in HO1D_parameters_m

    logical                            :: Dense_local                                                                  ! goes from 20 (= 0 verbose) to 24 (= maximum verbose) at this layer
    integer                            :: Verbose_local                                                                ! goes from 20 (= 0 verbose) to 24 (= maximum verbose) at this layer
    logical                            :: Debug_local
    
    !------------------------------------------------------Debugging options-----------------------------------------------------
    IF (PRESENT(Verbose)) THEN; Verbose_local = Verbose
    ELSE; Verbose_local = 20; END IF 
    IF (PRESENT(Debug))   THEN; Debug_local   = Debug
    ELSE; Debug_local = .FALSE.; END IF

    IF (Verbose_local > 20) WRITE(out_unit,*) 
    IF (Verbose_local > 20) WRITE(out_unit,*) "-------------------------------------------------INITIALIZING THE OPERATOR_ND OBJE&
                                              &CT-------------------------------------------------"; FLUSH(out_unit)

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "--- Arguments of MolecCav_Initialize_operator_ND :"
      WRITE(out_unit,*) "The <<OpND>> argument :"
      CALL Write(OpND)
      WRITE(out_unit,*) "The <<N_mat>> argument :"//TO_string(N_mat)
      WRITE(out_unit,*) "The <<N_cav>>  argument :"//TO_string(N_cav)
      IF (PRESENT(Dense)) WRITE(out_unit,*) "The <<Dense>> argument : "//TO_string(Dense)
      WRITE(out_unit,*) "--- End arguments of MolecCav_Initialize_operator_ND"
      FLUSH(out_unit)
    END IF
    
    !------------------------------------------Initializing the parameters of the 1D QHO------------------------------------------
    ! ON EN EST LA
    IF (PRESENT(Nb_op)) THEN; OpND%Nb_op  = Nb_op
    ELSE; OpND%Nb_op  = 4; END IF 
    ALLOCATE(OpND%Tab_op(0:OpND%Nb_op-1))
    
    !--------------------------------------Constructing the operators to build-------------------------------------
    IF (PRESENT(Dense)) THEN; Dense_local = Dense
    ELSE; Dense_local = .FALSE.; END IF

    CALL Initialize(OpND%Tab_op(0), "Identity",    Nb=Nb,           Dense=Dense_local, Verbose=Verbose_local, Debug=Debug_local)
    CALL Initialize(OpND%Tab_op(1), "Hamiltonian", Nb=Nb, w=w,      Dense=Dense_local, Verbose=Verbose_local, Debug=Debug_local)
    CALL Initialize(OpND%Tab_op(2), "Position",    Nb=Nb, w=w, m=m, Dense=Dense_local, Verbose=Verbose_local, Debug=Debug_local)
    CALL Initialize(OpND%Tab_op(3), "NbQuanta",    Nb=Nb,           Dense=Dense_local, Verbose=Verbose_local, Debug=Debug_local)

    IF (Verbose_local > 20) WRITE(out_unit,*) 
    IF (Verbose_local > 20) WRITE(out_unit,*) "--------------------------------------------------OPERATOR_ND OBJECT INITIALIZED--&
                                              &-----------------------------------------------"; FLUSH(out_unit)

  END SUBROUTINE MolecCav_Initialize_operator_ND


  SUBROUTINE MolecCav_Initialize_tab_operator(Elem_op, Operator_type, Nb, w, m, Dense, Verbose, Debug)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
    USE QDUtil_m
    USE Cavity_mode_m
    USE Matter_mode_m
    IMPLICIT NONE

    TYPE(Elem_op_t),            intent(inout) :: Elem_op                                                                         ! the object of type Elem_op_t to be constructed here
    character(len=*),           intent(in)    :: Operator_type                                                                   ! ex : "Hamiltonian", "Position", etc. (len=:) Expects to be allocatable, while (len=*) is dedicated to a procedure argument.
    integer,                    intent(in)    :: Nb                                                                              ! the HO/Cavity mode which the operator is relative to
    real(kind=Rkind), optional, intent(in)    :: w                                                                               ! the HO/Cavity mode which the operator is relative to
    real(kind=Rkind), optional, intent(in)    :: m                                                                               ! the HO/Cavity mode which the operator is relative to
    logical,          optional, intent(in)    :: Dense                                                                           ! if .TRUE. then the matrix storage will not be optimized and it will be stored as a Dense matrix
    integer,          optional, intent(in)    :: Verbose                                                                         ! cf. comments in HO1D_parameters_m
    logical,          optional, intent(in)    :: Debug                                                                           ! cf. comments in HO1D_parameters_m

    integer                                   :: Verbose_local                                                              ! goes from 25 (= 0 verbose) to 29 (= maximum verbose) at this layer
    logical                                   :: Debug_local

    !------------------------------------------------------Debugging options-----------------------------------------------------
    IF (PRESENT(Verbose)) THEN; Verbose_local = Verbose
    ELSE; Verbose_local = 20; END IF 
    IF (PRESENT(Debug))   THEN; Debug_local   = Debug
    ELSE; Debug_local = .FALSE.; END IF

    IF (Verbose_local > 25) WRITE(out_unit,*) 
    IF (Verbose_local > 25) WRITE(out_unit,*) "--------------------------------------------------INITIALIZING THE HO1D OPERATOR--&
                                              &------------------------------------------------"; FLUSH(out_unit)

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "--- Arguments of MolecCav_Initialize_HO1D_operator :"
      WRITE(out_unit,*) "The <<Elem_op>> argument :"
      CALL Write(Elem_op)
      WRITE(out_unit,*) "The <<Operator_type>> argument : "//Operator_type
      WRITE(out_unit,*) "The <<Nb>> argument : "//TO_string(Nb)
      IF (PRESENT(w)) WRITE(out_unit,*) "The <<w>> argument : "//TO_string(w)
      IF (PRESENT(m)) WRITE(out_unit,*) "The <<m>> argument : "//TO_string(m)
      IF (PRESENT(Dense)) WRITE(out_unit,*) "The <<Dense>> argument : "//TO_string(Dense)
      WRITE(out_unit,*) "--- End arguments of MolecCav_Construct_Operator_1D"
      FLUSH(out_unit)
    END IF

    IF (.NOT. PRESENT(w) .AND. TRIM(TO_lowercase(Operator_type)) == "hamiltonian") THEN
      WRITE(out_unit,*) "### missing w for the H"
      STOP "### missing w for the H"
    END IF 
    IF ((.NOT. PRESENT(w) .OR. .NOT. PRESENT(m)) .AND. TRIM(TO_lowercase(Operator_type)) == "position") THEN
      WRITE(out_unit,*) "### missing w or m for the x"
      STOP "### missing w or m for the x"
    END IF 
    
    !---------------------------------------First steps of the construction of the Operator--------------------------------------
    ALLOCATE(character(len=LEN_TRIM(Operator_type)) :: Elem_op%Operator_type)                                                   ! /!\ strings cannot be allocated the exact same way as tables ! /!\
    Elem_op%Operator_type = TO_lowercase(TRIM(Operator_type))                                                                   ! allocation on assignement (not anymore : supposed to work but caused dynamic allocation random errors at execution). Elem_op_type has the right lengths (no spaces added) thanks to len=* at declaration and it will fit the Op%op_type thanks to len=:, allocatable at declaration of the derived type. 

    IF (PRESENT(Dense)) Elem_op%Dense = Dense

    !---------------------------------------------Construction of the matrix Operator--------------------------------------------
    IF (Debug_local) THEN
      WRITE(out_unit,*); WRITE(out_unit,*) "--- The Elem_op_t object just before construction of its matrix representation"
      CALL Write(Elem_op)
      WRITE(out_unit,*) "--- End Elem_op_t object (just before construction of its matrix representation)"
    END IF 

    SELECT CASE (Elem_op%Operator_type)                                                                                         ! TO_lowercase avoid case sensitivity issues
      CASE ("identity")
        CALL Initialize_I(Identity=Elem_op,    Nb=Nb,           Verbose=Verbose_local, Debug=Debug_local)
    
      CASE ("hamiltonian")
        CALL Initialize_H(Hamiltonian=Elem_op, Nb=Nb, w=w,      Verbose=Verbose_local, Debug=Debug_local)
    
      CASE ("position")
        CALL Initialize_x(PositionOp=Elem_op,  Nb=Nb, w=w, m=m, Verbose=Verbose_local, Debug=Debug_local)
      
      CASE ("nbquanta")
        CALL Initialize_N(NbQuanta=Elem_op,    Nb=Nb,           Verbose=Verbose_local, Debug=Debug_local)

      CASE DEFAULT
        WRITE(out_unit,*) "### No Operator type recognized, please check the input of Initialize_HO1D_operator subroutine"
        STOP "### No Operator type recognized, please verify the input of Initialize_HO1D_operator subroutine"
    END SELECT

    IF (Verbose_local > 26) THEN
      IF (Verbose_local < 28) WRITE(out_unit,*)
      WRITE(out_unit,*) "--- HO1D operator constructed by MolecCav_Initialize_HO1D_operator :"
      CALL Write(Elem_op)
      WRITE(out_unit,*) "--- End HO1D operator constructed by MolecCav_Initialize_HO1D_operator"
    END IF

    IF (Verbose_local > 25) WRITE(out_unit,*) 
    IF (Verbose_local > 25) WRITE(out_unit,*) "-----------------------------------------------------HO1D OPERATOR INITIALIZED----&
                                              &------------------------------------------------"; FLUSH(out_unit)

  END SUBROUTINE MolecCav_Initialize_tab_operator

  
  SUBROUTINE MolecCav_Action_operator_ND_R1_real(Op_psi, OpND, i_op, Psi, Verbose, Debug)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    USE ND_indexes_m
    USE Cavity_mode_m
    USE Matter_mode_m
    IMPLICIT NONE

    real(kind=Rkind),     intent(inout) :: Op_psi(:)
    TYPE(Operator_ND_t), intent(in)    :: OpND
    integer,              intent(in)    :: i_op
    real(kind=Rkind),     intent(in)    :: Psi(:)
    integer, optional,    intent(in)    :: Verbose                                                                              ! cf. comments in HO1D_parameters_m
    logical, optional,    intent(in)    :: Debug                                                                                ! cf. comments in HO1D_parameters_m

    integer                             :: Nb
    integer                             :: Verbose_local                                                                   ! goes from 25 (= 0 verbose) to 29 (= maximum verbose) at this layer
    logical                             :: Debug_local

    !------------------------------------------------------Debugging options-----------------------------------------------------
    IF (PRESENT(Verbose)) THEN; Verbose_local = Verbose
    ELSE; Verbose_local = 20; END IF 
    IF (PRESENT(Debug))   THEN; Debug_local   = Debug
    ELSE; Debug_local = .FALSE.; END IF

    IF (Verbose_local > 25) WRITE(out_unit,*) 
    IF (Verbose_local > 25) WRITE(out_unit,*) "---------------------------------------COMPUTING ACTION OF THE HO1D OPERATOR OVER &
                                              &THE R1 WF---------------------------------------"; FLUSH(out_unit)

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "--- Arguments of MolecCav_Action_operator_ND :"
      WRITE(out_unit,*) "The <<OpND>> argument :"
      CALL Write(OpND)
      WRITE(out_unit,*) "The <<i_op>> argument :"//TO_string(i_op)
      WRITE(out_unit,*) "The <<Psi>> argument : "
      CALL Write_Vec(Psi, out_unit, 1, info="Psi")
      WRITE(out_unit,*) "The size of its vector : "//TO_string(Size(Psi))
      WRITE(out_unit,*) "--- End arguments of MolecCav_Action_operator_ND"
      FLUSH(out_unit)
    END IF
    
    !-----------------------------------------------------Checking dimensions----------------------------------------------------
    ! ALREADY CHECKED IN THE ACTIONS CODED IN ELEM_OP_M ! (fortunately btw, otherwise the \hat{I}d case should have test above)

    !---------------------------------------------Selection of the calculation method--------------------------------------------
    IF (OpND%Tab_op(i_op)%Operator_type == "identity") THEN ! N.B. we should have tested i_op == 0, easier and consistent with the algorithmic choices made so far
      Op_psi = Psi
    ELSE
      CALL Action(Op_psi=Op_psi, Elem_op=OpND%Tab_op(i_op), Psi=Psi, Verbose=Verbose_local, Debug=Debug_local)
    END IF

    IF (Verbose_local > 26) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "--- Resulting statevector from the action of the HO1D Elem_op on the Psi statevector operand, computed &
                        &by Action_HO1D_operator_R1 :"
      CALL Write_Vec(Op_psi, out_unit, 1, info="Op_Psi")
      WRITE(out_unit,*) "--- End resulting statevector computed by Action_HO1D_operator_R1"
    END IF
  
    IF (Verbose_local > 25) WRITE(out_unit,*) 
    IF (Verbose_local > 25) WRITE(out_unit,*) "----------------------------------------ACTION OF THE HO1D OPERATOR OVER THE R1 WF&
                                              & COMPUTED---------------------------------------"; FLUSH(out_unit)
  
  END SUBROUTINE MolecCav_Action_operator_ND_R1_real

  
  SUBROUTINE MolecCav_Action_operator_ND_R1_complex(Op_psi, OpND, i_op, Psi, Verbose, Debug)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    USE ND_indexes_m
    USE Cavity_mode_m
    USE Matter_mode_m
    IMPLICIT NONE

    complex(kind=Rkind),  intent(inout) :: Op_psi(:)
    TYPE(Operator_ND_t), intent(in)    :: OpND
    integer,              intent(in)    :: i_op
    complex(kind=Rkind),  intent(in)    :: Psi(:)
    integer, optional,    intent(in)    :: Verbose                                                                              ! cf. comments in HO1D_parameters_m
    logical, optional,    intent(in)    :: Debug                                                                                ! cf. comments in HO1D_parameters_m

    integer                             :: Nb
    integer                             :: Verbose_local                                                                   ! goes from 25 (= 0 verbose) to 29 (= maximum verbose) at this layer
    logical                             :: Debug_local

    !------------------------------------------------------Debugging options-----------------------------------------------------
    IF (PRESENT(Verbose)) THEN; Verbose_local = Verbose
    ELSE; Verbose_local = 20; END IF 
    IF (PRESENT(Debug))   THEN; Debug_local   = Debug
    ELSE; Debug_local = .FALSE.; END IF

    IF (Verbose_local > 25) WRITE(out_unit,*) 
    IF (Verbose_local > 25) WRITE(out_unit,*) "---------------------------------------COMPUTING ACTION OF THE HO1D OPERATOR OVER &
                                              &THE R1 WF---------------------------------------"; FLUSH(out_unit)

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "--- Arguments of MolecCav_Action_operator_ND :"
      WRITE(out_unit,*) "The <<OpND>> argument :"
      CALL Write(OpND)
      WRITE(out_unit,*) "The <<i_op>> argument :"//TO_string(i_op)
      WRITE(out_unit,*) "The <<Psi>> argument : "
      CALL Write_Vec(Psi, out_unit, 1, info="Psi")
      WRITE(out_unit,*) "The size of its vector : "//TO_string(Size(Psi))
      WRITE(out_unit,*) "--- End arguments of MolecCav_Action_operator_ND"
      FLUSH(out_unit)
    END IF
    
    !-----------------------------------------------------Checking dimensions----------------------------------------------------
    ! ALREADY CHECKED IN THE ACTIONS CODED IN ELEM_OP_M !

    !---------------------------------------------Selection of the calculation method--------------------------------------------
    IF (OpND%Tab_op(i_op)%Operator_type == "identity") THEN 
      Op_psi = Psi
    ELSE
      CALL Action(Op_psi=Op_psi, Elem_op=OpND%Tab_op(i_op), Psi=Psi, Verbose=Verbose_local, Debug=Debug_local)
    END IF

    IF (Verbose_local > 26) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "--- Resulting statevector from the action of the HO1D Elem_op on the Psi statevector operand, computed &
                        &by Action_HO1D_operator_R1 :"
      CALL Write_Vec(Op_psi, out_unit, 1, info="Op_Psi")
      WRITE(out_unit,*) "--- End resulting statevector computed by Action_HO1D_operator_R1"
    END IF
  
    IF (Verbose_local > 25) WRITE(out_unit,*) 
    IF (Verbose_local > 25) WRITE(out_unit,*) "----------------------------------------ACTION OF THE HO1D OPERATOR OVER THE R1 WF&
                                              & COMPUTED---------------------------------------"; FLUSH(out_unit)
  
  END SUBROUTINE MolecCav_Action_operator_ND_R1_complex

  
  SUBROUTINE MolecCav_Write_operator_ND(OpND)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
    USE QDUtil_m
    USE Cavity_mode_m
    USE Matter_mode_m
    IMPLICIT NONE 
    
    TYPE(Operator_ND_t), intent(in) :: OpND

    integer                          :: i_op

    WRITE(out_unit,*) "____________________________________The parameters of the 1D QHO____________________________________"
    WRITE(out_unit,*) "|Basis set size of the HO (OpND%Nb)                                           | "//TO_string(OpND%Nb)
    WRITE(out_unit,*) "|______________________________________________________________________________|______________________"
    WRITE(out_unit,*) "|Eigenpulsation of the HO (OpND%w)                                            | "//TO_string(OpND%w)
    WRITE(out_unit,*) "|______________________________________________________________________________|______________________"
    WRITE(out_unit,*) "|Mass associated with the HO (OpND%m)                                         | "//TO_string(OpND%m)
    IF (ALLOCATED(OpND%Tab_op)) THEN 
      WRITE(out_unit,*) "|______________________________________________________________________________|______________________"
      WRITE(out_unit,*) "|QHO's Tab_op holds the following nb of operators (Size(OpND%Tab_op))         | ", Size(OpND%Tab_op)
      WRITE(out_unit,*) "|______________________________________________________________________________|______________________"
      WRITE(out_unit,*) "|The "//TO_string(i_op)//"^{th} operator associated with the HO (OpND%Tab_op("//TO_string(i_op)//")) : &
                        &               |"
      WRITE(out_unit,*) "|______________________________________________________________________________|"
      DO i_op = 0, Size(OpND%Tab_op)-1
        CALL Write(OpND%Tab_op(i_op))
      END DO 
    ELSE 
      WRITE(out_unit,*) "|______________________________________________________________________________|______________________"
      WRITE(out_unit,*) "|QHO's Tab_op is not allocated (OpND%Tab_op)                                  | /"
      WRITE(out_unit,*) "|______________________________________________________________________________|______________________"
    END IF
    WRITE(out_unit,*) "|Number of grid points of the QHO DOF (OpND%Nq)                               | "//TO_string(OpND%Nq)
    WRITE(out_unit,*) "|______________________________________________________________________________|______________________"
    WRITE(out_unit,*) "|Equilibrium position of the HO (OpND%Eq_pos)                                 | "//TO_string(OpND%Eq_pos)
    WRITE(out_unit,*) "|______________________________________________________________________________|______________________"
    WRITE(out_unit,*) "|Change in variable coefficient for the DOF grid (OpND%Scale_q)               | "//TO_string(OpND%Scale_q)
    WRITE(out_unit,*) "|________________________________________End HO1D parameters___________________|______________________"
    FLUSH(out_unit)

  END SUBROUTINE MolecCav_Write_operator_ND


  SUBROUTINE MolecCav_Deallocate_operator_ND(OpND, Verbose, Debug)
    USE QDUtil_m
    USE Cavity_mode_m
    USE Matter_mode_m
    IMPLICIT NONE 

    TYPE(Operator_ND_t), intent(inout) :: OpND
    integer, optional,    intent(in)    :: Verbose                                                                                 ! cf. comments in HO1D_parameters_m
    logical, optional,    intent(in)    :: Debug                                                                                   ! cf. comments in HO1D_parameters_m

    integer                             :: i_op
    integer                             :: Verbose_local                                                                      ! goes from 25 (= 0 verbose) to 29 (= maximum verbose) at this layer
    logical                             :: Debug_local

    !------------------------------------------------------Debugging options-----------------------------------------------------
    IF (PRESENT(Verbose)) THEN; Verbose_local = Verbose
    ELSE; Verbose_local = 20; END IF 
    IF (PRESENT(Debug))   THEN; Debug_local   = Debug
    ELSE; Debug_local = .FALSE.; END IF

    IF (Debug_local) THEN
      WRITE(out_unit,*) "--- The OpND to be deallocated :"
      CALL Write(OpND)
      WRITE(out_unit,*) "--- End OpND to be deallocated"
    END IF 

    !-----------------------------Deallocating the HO1D operator object----------------------------
    IF (Verbose_local > 27) WRITE(out_unit,*)
    IF (Verbose_local > 27) WRITE(out_unit,*) "-----------------------------------------------Deallocating the OpND obje&
                                              &ct----------------------------------------------"
  
    OpND%Nb = 0
    OpND%w  = ZERO
    OpND%m  = ZERO
    DO i_op = 0, SIZE(OpND%Tab_op)-1
      CALL Dealloc(OpND%Tab_op(i_op), Verbose=Verbose_local, Debug=Debug_local)
    END DO
    IF (ALLOCATED(OpND%Tab_op)) DEALLOCATE(OpND%Tab_op)
    OpND%Nq      = 0
    OpND%Eq_pos  = -ONE
    OpND%Scale_q = HUGE(ONE)

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "--- The OpND object after having been deallocated :"
      CALL Write(OpND)
      WRITE(out_unit,*) "--- End dellocating OpND"
    END IF

  END SUBROUTINE MolecCav_Deallocate_operator_ND


END MODULE