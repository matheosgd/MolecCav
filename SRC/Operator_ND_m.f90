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


  TYPE                   :: Operator_ND_t                                        ! N_mat, N_cav canNOT be part of the derived type because they are not specific of one operator_ND, but parameters of the whole system/calculation. All OpND will have the same. 
    integer, allocatable :: tab_indexes_mat_op(:)                                ! N_mat : nb of DOF of the matter subsystem = nb of vibrational modes
    integer, allocatable :: tab_indexes_cav_op(:)                                ! N_cav : nb of DOF of the cavity subsystem = nb of cavity modes
  END TYPE

  TYPE(Matter_mode_t), allocatable :: tab_mat_ops(:)
  TYPE(Cavity_mode_new_t), allocatable :: tab_cav_ops(:)


  PRIVATE

  PUBLIC Operator_ND_t, Initialize, Action, Write, Dealloc

  INTERFACE Initialize
    MODULE PROCEDURE MolecCav_Initialize_operator_ND
  END INTERFACE
  INTERFACE Initialize_tabs_ops
    MODULE PROCEDURE MolecCav_Initialize_tabs_operators
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


  SUBROUTINE MolecCav_Initialize_operator_ND(OpND, Mat_operators, Cav_operators, nio, Dense, Verbose, Debug) ! no need for N_mat and N_cav explicitly : they are SIZE(Mat_op and Cav_op)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
    USE QDUtil_m
    USE Cavity_mode_m
    USE Matter_mode_m
    IMPLICIT NONE
  
    TYPE(Operator_ND_t), intent(inout) :: OpND
    character(len=*),    intent(in)    :: Mat_operators ! syntax : '<op_mode_1>, <op_mode_2>, ..., <op_mode_N_mat>', ex : 'hamiltonian, Identity'. Not case sensitive, ' ' <=> \otimes
    character(len=*),    intent(in)    :: Cav_operators ! syntax : '<op_mode_1>, <op_mode_2>, ..., <op_mode_N_cav>', ex : 'hamiltonian'.   Not case sensitive, ' ' <=> \otimes. This exemple means OpND = H_mat1\otimesI_mat2\otimesH_cav
    integer,             intent(in)    :: nio
    logical, optional,   intent(in)    :: Dense                                                                        ! cf. comments in HO1D_parameters_m
    integer, optional,   intent(in)    :: Verbose                                                                      ! cf. comments in HO1D_parameters_m
    logical, optional,   intent(in)    :: Debug                                                                        ! cf. comments in HO1D_parameters_m

    integer                            :: i_mode, i_op
    character(len=25), allocatable     :: Mat_operators_local(:) ! allocatable does not refer to the len=* of each str of the table but to the size of the table i.e. the nb of str within !
    character(len=25), allocatable     :: Cav_operators_local(:) ! /!\ Cannot declare (len=*) because syntaxe dedicated for dummy var (so intent(in)) or parameters of the program. Here an length needs to be provided : just take arbitrary long, and will automatically be completed w. blanks
    integer                            :: N_mat
    integer                            :: N_cav
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
      WRITE(out_unit,*) "The <<Mat_operators>>  argument :"//Mat_operators
      WRITE(out_unit,*) "The <<Cav_operators>>  argument :"//Cav_operators
      IF (PRESENT(Dense)) WRITE(out_unit,*) "The <<Dense>> argument : "//TO_string(Dense)
      WRITE(out_unit,*) "Are the module's <<tab_mat/cav_ops>> allocated ? "//TO_string(ALLOCATED(tab_mat_ops))//TO_string(ALLOCAT&
      &ED(tab_cav_ops))
      WRITE(out_unit,*) "--- End arguments of MolecCav_Initialize_operator_ND"
      FLUSH(out_unit)
    END IF
    
    !------------------------------------------Initializing the procedure------------------------------------------
    IF (PRESENT(Dense)) THEN; Dense_local = Dense
    ELSE; Dense_local = .FALSE.; END IF

    IF (LEN_TRIM(Mat_operators)==0) THEN
      N_mat = 0
    ELSE
      N_mat = 1
      DO i_mode = 1, LEN_TRIM(Mat_operators)
        IF (Mat_operators(i_mode:i_mode)==',') N_mat = N_mat + 1
      END DO 
    END IF 
    IF (LEN_TRIM(Cav_operators)==0) THEN
      N_cav = 0
    ELSE
      N_cav = 1
      DO i_mode = 1, LEN_TRIM(Cav_operators)
        IF (Cav_operators(i_mode:i_mode)==',') N_cav = N_cav + 1
      END DO 
    END IF 
    
    IF (.NOT. ALLOCATED(tab_mat_ops) .OR. .NOT. ALLOCATED(tab_cav_ops)) THEN
      CALL Initialize_tabs_ops(N_mat=N_mat, N_cav=N_cav, Dense=Dense_local, Verbose=Verbose_local, Debug=Debug_local)
    END IF
    
    !--------------------------------------Constructing the OpND = parsing the Mat/Cav_operators strings-------------------------------------
    ALLOCATE(Mat_operators_local(N_mat))
    ALLOCATE(Cav_operators_local(N_cav))
    ALLOCATE(OpND%tab_indexes_mat_op(N_mat))
    ALLOCATE(OpND%tab_indexes_cav_op(N_cav))

    READ(unit=Mat_operators, fmt=*) Mat_operators_local ! /!\ neither trimmed nor lowercased so far /!\
    READ(unit=Cav_operators, fmt=*) Cav_operators_local

    DO i_mode = 1, N_mat
      DO i_op = 0, tab_mat_ops(i_mode)%Nb_op-1
        IF (TO_lowercase(TRIM(Mat_operators_local(i_mode))) == tab_mat_ops(i_mode)%Tab_op(i_op)%Operator_type) THEN
          OpND%tab_indexes_mat_op(i_mode) = i_op
          EXIT 
        ELSE IF (i_op == tab_mat_ops(i_mode)%Nb_op-1) THEN
          WRITE(out_unit,*) "### No operator name recognized. Please check arguments of MolecCav_Initialize_operator_ND"
          STOP "### No operator name recognized in MolecCav_Initialize_operator_ND"
        END IF 
      END DO
    END DO 

    DO i_mode = 1, N_cav
      DO i_op = 0, tab_cav_ops(i_mode)%Nb_op-1
        IF (TO_lowercase(TRIM(Cav_operators_local(i_mode))) == tab_cav_ops(i_mode)%Tab_op(i_op)%Operator_type) THEN
          OpND%tab_indexes_cav_op(i_mode) = i_op
          EXIT 
        ELSE IF (i_op == tab_mat_ops(i_mode)%Nb_op-1) THEN
          WRITE(out_unit,*) "### No operator name recognized. Please check arguments of MolecCav_Initialize_operator_ND"
          STOP "### No operator name recognized in MolecCav_Initialize_operator_ND"
        END IF 
      END DO
    END DO

    IF (Verbose_local > 20) WRITE(out_unit,*) 
    IF (Verbose_local > 20) WRITE(out_unit,*) "--------------------------------------------------OPERATOR_ND OBJECT INITIALIZED--&
    &-----------------------------------------------"; FLUSH(out_unit)

  END SUBROUTINE MolecCav_Initialize_operator_ND


  SUBROUTINE MolecCav_Initialize_tabs_operators(N_mat, N_cav, Dense, Verbose, Debug) ! tab_mat_ops, tab_cav_ops are shared in all the module. Modified in fly in this sub. => no need to pass them in argument here
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
    USE QDUtil_m
    USE Cavity_mode_m
    USE Matter_mode_m
    IMPLICIT NONE

    integer,           intent(in)    :: N_mat                                                                                    ! the HO/Cavity mode which the operator is relative to
    integer,           intent(in)    :: N_cav                                                                                    ! the HO/Cavity mode which the operator is relative to
    logical, optional, intent(in)    :: Dense                                                                                    ! if .TRUE. then the matrix storage will not be optimized and it will be stored as a Dense matrix
    integer, optional, intent(in)    :: Verbose                                                                                  ! cf. comments in HO1D_parameters_m
    logical, optional, intent(in)    :: Debug                                                                                    ! cf. comments in HO1D_parameters_m

    integer                          :: i
    logical                          :: Dense_local
    integer                          :: Verbose_local                                                                            ! goes from 25 (= 0 verbose) to 29 (= maximum verbose) at this layer
    logical                          :: Debug_local

    !------------------------------------------------------Debugging options-----------------------------------------------------
    IF (PRESENT(Verbose)) THEN; Verbose_local = Verbose
    ELSE; Verbose_local = 20; END IF 
    IF (PRESENT(Debug))   THEN; Debug_local   = Debug
    ELSE; Debug_local = .FALSE.; END IF

    IF (Verbose_local > 25) WRITE(out_unit,*) 
    IF (Verbose_local > 25) WRITE(out_unit,*) "--------------------------------------------------INITIALIZING THE OPERATORS' TABL&
    &ES--------------------------------------------------"; FLUSH(out_unit)

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "--- Arguments of MolecCav_Initialize_HO1D_operator :"
      WRITE(out_unit,*) "The module's <<tab_mat_ops>> allocated ? "//TO_string(ALLOCATED(tab_mat_ops))
      WRITE(out_unit,*) "The module's <<tab_cav_ops>> allocated ? "//TO_string(ALLOCATED(tab_cav_ops))
      WRITE(out_unit,*) "The <<N_mat>> argument : "//TO_string(N_mat)
      WRITE(out_unit,*) "The <<N_cav>> argument : "//TO_string(N_cav)
      IF (PRESENT(Dense)) WRITE(out_unit,*) "The <<Dense>> argument : "//TO_string(Dense)
      WRITE(out_unit,*) "--- End arguments of MolecCav_Construct_Operator_1D"
      FLUSH(out_unit)
    END IF
    
    !---------------------------------------First steps of the construction of the tables--------------------------------------
    IF (PRESENT(Dense)) THEN; Dense_local = Dense
    ELSE; Dense_local = .FALSE.; END IF

    ALLOCATE(tab_mat_ops(N_mat))
    ALLOCATE(tab_cav_ops(N_cav))

    !---------------------------------------------Construction of the whole tables--------------------------------------------
    DO i = 1, N_mat
      CALL Initialize(MatMode=tab_mat_ops(i), nio=in_unit, Dense=Dense_local, Verbose=Verbose_local, Debug=Debug_local) ! N.B. => nml shall be constructed as Matmode1\Matmode2\...\Cavmode1\...\CavmodeN_cav\
    END DO

    DO i = 1, N_cav
      CALL Initialize(CavMode=tab_cav_ops(i), nio=in_unit, Dense=Dense_local, Verbose=Verbose_local, Debug=Debug_local)
    END DO 

    IF (Verbose_local > 26) THEN
      IF (Verbose_local < 28) WRITE(out_unit,*)
      WRITE(out_unit,*) "--- tabs_ops constructed by MolecCav_Initialize_tabs_operators :"
      DO i = 1, N_mat
        WRITE(out_unit,*); WRITE(out_unit,*) "--- tabs_mat_ops"//TO_string(i)//" :"
        CALL Write(tab_mat_ops(i))
      END DO 
      DO i = 1, N_cav
        WRITE(out_unit,*); WRITE(out_unit,*) "--- tabs_cav_ops"//TO_string(i)//" :"
        CALL Write(tab_cav_ops(i))
      END DO 
      WRITE(out_unit,*) "--- End tabs_ops constructed by MolecCav_Initialize_tabs_operators"
    END IF

    IF (Verbose_local > 25) WRITE(out_unit,*) 
    IF (Verbose_local > 25) WRITE(out_unit,*) "-----------------------------------------------------OPERATOR'S TABLES INITIALIZED&
    &----------------------------------------------------"; FLUSH(out_unit)

  END SUBROUTINE MolecCav_Initialize_tabs_operators

  
  SUBROUTINE MolecCav_Action_operator_ND_R1_real(Op_psi, OpND, Psi, Verbose, Debug) ! Psi is ND AND R1
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    USE ND_indexes_m
    USE Cavity_mode_m
    USE Matter_mode_m
    IMPLICIT NONE

    real(kind=Rkind),     intent(inout) :: Op_psi(:)
    TYPE(Operator_ND_t),  intent(in)    :: OpND
    real(kind=Rkind),     intent(in)    :: Psi(:)
    integer, optional,    intent(in)    :: Verbose                                                                              ! cf. comments in HO1D_parameters_m
    logical, optional,    intent(in)    :: Debug                                                                                ! cf. comments in HO1D_parameters_m

    integer                             :: N_tot = 1, i_mode
    integer                             :: Verbose_local                                                                   ! goes from 25 (= 0 verbose) to 29 (= maximum verbose) at this layer
    logical                             :: Debug_local

    !------------------------------------------------------Debugging options-----------------------------------------------------
    IF (PRESENT(Verbose)) THEN; Verbose_local = Verbose
    ELSE; Verbose_local = 20; END IF 
    IF (PRESENT(Debug))   THEN; Debug_local   = Debug
    ELSE; Debug_local = .FALSE.; END IF

    IF (Verbose_local > 25) WRITE(out_unit,*) 
    IF (Verbose_local > 25) WRITE(out_unit,*) "---------------------------------------COMPUTING ACTION OF THE OpND OVER &
                                              &THE R1 ND WF---------------------------------------"; FLUSH(out_unit)

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "--- Arguments of MolecCav_Action_operator_ND :"
      WRITE(out_unit,*) "The <<OpND>> argument :"
      CALL Write(OpND)
      WRITE(out_unit,*) "The <<Psi>> argument : "
      CALL Write_Vec(Psi, out_unit, 1, info="Psi")
      WRITE(out_unit,*) "The size of its vector : "//TO_string(Size(Psi))
      WRITE(out_unit,*) "--- End arguments of MolecCav_Action_operator_ND"
      FLUSH(out_unit)
    END IF
    
    !-----------------------------------------------------Checking dimensions----------------------------------------------------
    ! THE DIMENSIONS OF EACH 1D MATMUL WILL BE TESTED IN THE ACTIONS CODED IN ELEM_OP_M ! (fortunately btw, otherwise the \hat{I}d case should have test above)
    IF (SIZE(OpND%tab_indexes_mat_op) /= 0) THEN
      DO i_mode = 1, SIZE(OpND%tab_indexes_mat_op)
        N_tot = N_tot * tab_mat_ops(OpND%tab_indexes_mat_op(i_mode))%Nb
      END DO
    END IF 
    IF (SIZE(OpND%tab_indexes_cav_op) /= 0) THEN
      DO i_mode = 1, SIZE(OpND%tab_indexes_cav_op)
        N_tot = N_tot * tab_cav_ops(OpND%tab_indexes_cav_op(i_mode))%Nb
      END DO
    END IF 

    IF (SIZE(Psi, dim=1) /= N_tot) THEN
      WRITE(out_unit,*) "### The dimension of the wavevector Psi does not match the dimension product of the operators to"
      WRITE(out_unit,*) "   Size(Psi, dim=1) = "//TO_string(Size(Psi, dim=1))//"; PRODUCT(Dims) = "//TO_string(N_tot)
      STOP "### The dimension of Psi does not match the dimension product of the operators."
    END IF
    
    !---------------------------------------------Selection of the calculation method--------------------------------------------
    
    
    
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
    IMPLICIT NONE 
    
    TYPE(Operator_ND_t), intent(in) :: OpND

    IF (ALLOCATED(OpND%tab_indexes_mat_op)) THEN
      WRITE(out_unit,*) "_______________________________________The ND Operator object_______________________________________"
      WRITE(out_unit,*) "|The number of matter modes (SIZE(OpND%tab_indexes_mat_op))                    | "//TO_string(SIZE(OpND%ta&
      &b_indexes_mat_op))
      WRITE(out_unit,*) "|______________________________________________________________________________|______________________"
      FLUSH(out_unit)
      IF (SIZE(OpND%tab_indexes_cav_op) /= 0) THEN
        WRITE(out_unit,*) "|The operators taking action on the matter modes (OpND%tab_indexes_mat_op) :   | "
        CALL Write_Vec(OpND%tab_indexes_mat_op, out_unit, SIZE(OpND%tab_indexes_mat_op), info="tab_indexes_mat_op")
      ELSE
        WRITE(out_unit,*) "|No cavity mode, so no operators to write                                      | "
      END IF
    WRITE(out_unit,*) "|______________________________________________________________________________|______________________"
    ELSE 
      WRITE(out_unit,*) "|_____________________________The ND Operator object___________________________|"
      WRITE(out_unit,*) "|The table of matter operators is NOT allocated (OpND%tab_indexes_mat_op)      |"
      WRITE(out_unit,*) "|______________________________________________________________________________|"
    END IF
    FLUSH(out_unit)
    IF (ALLOCATED(OpND%tab_indexes_mat_op)) THEN
      WRITE(out_unit,*) "|The number of cavity modes (SIZE(OpND%tab_indexes_cav_op))                    | "//TO_string(SIZE(OpND%ta&
      &b_indexes_cav_op))
      WRITE(out_unit,*) "|______________________________________________________________________________|______________________"
      FLUSH(out_unit)
      IF (SIZE(OpND%tab_indexes_cav_op) /= 0) THEN
        WRITE(out_unit,*) "|The operators taking action on the cavity modes (OpND%tab_indexes_cav_op) :   | "
        CALL Write_Vec(OpND%tab_indexes_cav_op, out_unit, SIZE(OpND%tab_indexes_cav_op), info="tab_indexes_cav_op")
      ELSE
        WRITE(out_unit,*) "|No cavity mode, so no operators to write                                      | "
      END IF
      WRITE(out_unit,*) "|_____________________________________End ND Operator object___________________|"
    ELSE 
      WRITE(out_unit,*) "|The table of cavity operators is NOT allocated (OpND%tab_indexes_cav_op)      |"
      WRITE(out_unit,*) "|_____________________________________End ND Operator object___________________|"
    END IF
    FLUSH(out_unit)

  END SUBROUTINE MolecCav_Write_operator_ND


  SUBROUTINE MolecCav_Deallocate_operator_ND(OpND, All, Verbose, Debug)
    USE QDUtil_m
    USE Cavity_mode_m
    USE Matter_mode_m
    IMPLICIT NONE 

    TYPE(Operator_ND_t), intent(inout) :: OpND
    logical, optional,   intent(in)    :: All
    integer, optional,   intent(in)    :: Verbose                                                                                 ! cf. comments in HO1D_parameters_m
    logical, optional,   intent(in)    :: Debug                                                                                   ! cf. comments in HO1D_parameters_m

    integer                            :: i_op
    logical                            :: All_local
    integer                            :: Verbose_local                                                                      ! goes from 25 (= 0 verbose) to 29 (= maximum verbose) at this layer
    logical                            :: Debug_local

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
  
    IF (ALLOCATED(OpND%tab_indexes_mat_op)) DEALLOCATE(OpND%tab_indexes_mat_op)
    IF (ALLOCATED(OpND%tab_indexes_cav_op)) DEALLOCATE(OpND%tab_indexes_cav_op)

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "--- The OpND object after having been deallocated :"
      CALL Write(OpND)
      WRITE(out_unit,*) "--- End dellocating OpND"
    END IF

  END SUBROUTINE MolecCav_Deallocate_operator_ND


END MODULE