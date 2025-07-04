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
  TYPE(Cavity_mode_t), allocatable :: tab_cav_ops(:)


  PRIVATE

  PUBLIC Operator_ND_t, Initialize, Action, Get, Write, Dealloc&
  &, tab_mat_ops

  INTERFACE Initialize
    MODULE PROCEDURE MolecCav_Initialize_operator_ND
  END INTERFACE
  INTERFACE Initialize_tabs_ops
    MODULE PROCEDURE MolecCav_Initialize_tabs_operators
  END INTERFACE
  INTERFACE Action
    MODULE PROCEDURE MolecCav_Action_operator_ND_R1_real, MolecCav_Action_operator_ND_R1_complex
  END INTERFACE
  INTERFACE Get
    MODULE PROCEDURE MolecCav_Get_OpND_parameter_integer, MolecCav_Get_OpND_parameter_real
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
    integer,             intent(in)    :: nio           ! Used only to initialize tab_mat_ops and tab_cav_ops, i.e. ONLY ONCE !!!
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

    !######### Here has only purpose to compute N_mat and N_cav using the OpND declaration. Maybe these two paramaters will passed as arguments later, through a potential Basis_ND_t general object #########
    IF (LEN_TRIM(Mat_operators)==0) THEN
      N_mat = 0
      WRITE(out_unit,*) "########################## WARNING ########################## WARNING ########################## WARNING #&
                      &########################"
      WRITE(out_unit,*) "                          The code is now used without any matter mode to compute : cavity alone "
      WRITE(out_unit,*) "########################## WARNING ########################## WARNING ########################## WARNING #&
                      &########################"
    ELSE
      N_mat = 1
      DO i_mode = 1, LEN_TRIM(Mat_operators)
        IF (Mat_operators(i_mode:i_mode)==',') N_mat = N_mat + 1
      END DO 
    END IF 
    IF (LEN_TRIM(Cav_operators)==0) THEN
      N_cav = 0
      WRITE(out_unit,*) "########################## WARNING ########################## WARNING ########################## WARNING #&
                      &########################"
      WRITE(out_unit,*) "                          The code is now used without any cavity mode to compute : matter alone "
      WRITE(out_unit,*) "########################## WARNING ########################## WARNING ########################## WARNING #&
                      &########################"
    ELSE
      N_cav = 1
      DO i_mode = 1, LEN_TRIM(Cav_operators)
        IF (Cav_operators(i_mode:i_mode)==',') N_cav = N_cav + 1
      END DO 
    END IF 
    !######################################################################################################################################################################################################### 
  
    IF (ALLOCATED(tab_mat_ops) .AND. ALLOCATED(tab_cav_ops)) THEN
      IF (SIZE(tab_mat_ops)/=N_mat) THEN
        WRITE(out_unit,*) "### The number of declared matter modes is not consistent with the one previously declared at last cal&
        &l of Initialize_operator_ND."
        WRITE(out_unit,*) "    Please check arguments to keep consistency in the system's size along the simulation."
        STOP "### The number of declared matter modes is not consistent with the one previously declared."
      ELSE IF (SIZE(tab_cav_ops)/=N_cav) THEN
        WRITE(out_unit,*) "### The number of declared cavity modes is not consistent with the one previously declared at last cal&
        &l of Initialize_operator_ND."
        WRITE(out_unit,*) "    Please check arguments to keep consistency in the system's size along the simulation."
        STOP "### The number of declared cavity modes is not consistent with the one previously declared."
      END IF
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
        ELSE IF (i_op == tab_cav_ops(i_mode)%Nb_op-1) THEN
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
    ! integer,           intent(in)    :: nio ! TO ADD /!\ in the call initialize also and in the call of this sub just above
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
      CALL Initialize(Matmode=tab_mat_ops(i), nio=in_unit, Dense=Dense_local, Verbose=Verbose_local, Debug=Debug_local) ! N.B. => nml shall be constructed as OpND1\OpND2\...\Cavmode1\...\CavmodeN_cav\
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
    USE Cavity_mode_m
    USE Matter_mode_m
    IMPLICIT NONE

    real(kind=Rkind),     intent(inout) :: Op_psi(:)
    TYPE(Operator_ND_t),  intent(in)    :: OpND
    real(kind=Rkind),     intent(in)    :: Psi(:)
    integer, optional,    intent(in)    :: Verbose                                                                              ! cf. comments in HO1D_parameters_m
    logical, optional,    intent(in)    :: Debug                                                                                ! cf. comments in HO1D_parameters_m

    integer                             :: N_mat, N_cav, NB, N1, N2, N3, i_mode, i_1, i_3 ! /!\ N_mat, N_cav are not the basis sets sizes but the respective number of DOF of the matter and cavity subsystems /!\
    integer,          allocatable       :: Ranks_sizes(:)
    real(kind=Rkind), allocatable       :: Cube(:,:,:), Op_cube(:,:,:)
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
      !CALL Write(OpND)
      IF (SIZE(OpND%tab_indexes_mat_op)==0) WRITE(out_unit,*) "tab mat op : $\empty$"; FLUSH(out_unit)
      IF (SIZE(OpND%tab_indexes_mat_op)/=0) WRITE(out_unit,*) "tab mat op : "//TO_string(OpND%tab_indexes_mat_op); FLUSH(out_unit)
      IF (SIZE(OpND%tab_indexes_cav_op)==0) WRITE(out_unit,*) "tab cav op : $\empty$"
      IF (SIZE(OpND%tab_indexes_cav_op)/=0) WRITE(out_unit,*) "tab cav op : "//TO_string(OpND%tab_indexes_cav_op)
      WRITE(out_unit,*) "The <<Psi>> argument : "
      CALL Write_Vec(Psi, out_unit, 1, info="Psi")
      WRITE(out_unit,*) "The size of its vector : "//TO_string(Size(Psi))
      WRITE(out_unit,*) "--- End arguments of MolecCav_Action_operator_ND"
      FLUSH(out_unit)
    END IF
    
    !-----------------------------------------------------Checking dimensions----------------------------------------------------
    ! THE DIMENSIONS OF EACH 1D MATMUL WILL BE TESTED IN THE ACTIONS CODED IN ELEM_OP_M !
    N_mat = SIZE(OpND%tab_indexes_mat_op)
    N_cav = SIZE(OpND%tab_indexes_cav_op)
    ALLOCATE(Ranks_sizes(N_mat + N_cav))

    IF (N_mat /= 0) THEN
      Ranks_sizes(1:N_mat) = tab_mat_ops(1:N_mat)%Nb
    END IF 
    IF (N_cav /= 0) THEN
      Ranks_sizes(1+N_mat:N_cav+N_mat) = tab_cav_ops(1:N_cav)%Nb
    END IF 
    NB = PRODUCT(Ranks_sizes) ! better not to have to calculate again each time needed

    IF (SIZE(Psi, dim=1) /= NB) THEN
      WRITE(out_unit,*) "### The dimension of the wavevector Psi does not match the dimension product of the operators to"
      WRITE(out_unit,*) "   Size(Psi, dim=1) = "//TO_string(Size(Psi, dim=1))//"; PRODUCT(Dims) = "//TO_string(NB)
      STOP "### The dimension of Psi does not match the dimension product of the operators."
    END IF

    !----------------------------Computation using reshape----------------------------------    
      !-----------------------Initialization befor the first loop---------------------------    
!--to compute from first mode to last-!
!    N1 = 1                           !
!    N2 = Ranks_sizes(1)              ! 
!    N3 = NB / N2                     !
!-------------------------------------!
!--to compute from last mode to first-!
    N3 = 1                            !
    N2 = Ranks_sizes(N_mat + N_cav)   !
    N1 = NB / N2                      !
!-------------------------------------!
    ALLOCATE(Cube(   N1, N2, N3))
    ALLOCATE(Op_cube(N1, N2, N3))

    Cube = RESHAPE(Psi, [N1, N2, N3])


!    DO i_mode = 1, N_mat + N_cav     ! from first mode to last
    DO i_mode = N_mat + N_cav, 1, -1 ! from last mode to first
      IF (Debug_local) WRITE(out_unit,*)
      IF (Debug_local) WRITE(out_unit,*) "--- i_mode = "//TO_string(i_mode)
      IF (Debug_local) WRITE(out_unit,*)
      IF (Debug_local) WRITE(out_unit,*) "N1, N2, N3 = "//TO_string(N1)//", "//TO_string(N2)//", "//TO_string(N3)

      !-----------------------Action-----------------------------    
      IF (i_mode <= N_mat) THEN
        DO i_3 = 1, N3
          DO i_1 = 1, N1
            CALL Action(Op_psi=Op_cube(i_1,:,i_3), MatMode=tab_mat_ops(i_mode), i_op=OpND%tab_indexes_mat_op(i_mode), Psi=Cube(i_1,&
            &:,i_3), Verbose=Verbose, Debug=Debug_local)
          END DO 
        END DO 
      ELSE 
        DO i_3 = 1, N3
          DO i_1 = 1, N1
            CALL Action(Op_psi=Op_cube(i_1,:,i_3), CavMode=tab_cav_ops(i_mode-N_mat), i_op=OpND%tab_indexes_cav_op(i_mode-N_mat),&
            & Psi=Cube(i_1,:,i_3), Verbose=Verbose, Debug=Debug_local)
          END DO 
        END DO 
      END IF
!      IF (i_mode == N_mat + N_cav) EXIT ! from first mode to last
      IF (i_mode == 1) EXIT ! from last mode to first
  
      !-----------------------Reinitialization for next loop-----------------------------    
      DEALLOCATE(Cube)
!--to compute from first mode to last--!
!      N1 = N1 * N2                    !
!      N2 = Ranks_sizes(i_mode + 1)    !
!      N3 = N3 / N2                    !
!--------------------------------------!
      N3 = N3 * N2                     !
      N2 = Ranks_sizes(i_mode - 1)     !
      N1 = N1 / N2                     !
!--------------------------------------!
      ALLOCATE(Cube(N1,N2,N3))
      Cube = RESHAPE(Op_cube, [N1, N2, N3])
      DEALLOCATE(Op_cube)
      ALLOCATE(Op_cube(N1,N2,N3))
    END DO 

      !-----------------------Recovering Op_psi (R1)-----------------------------    
    Op_psi = RESHAPE(Op_cube, [NB])
  
    !--------------Conclusion----------------
    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "--- Resulting statevector from the action of the ND Operator on the Psi statevector operand, computed &
                        &by MolecCav_Action_operator_ND_R1_real :"
      CALL Write_Vec(Op_psi, out_unit, 1, info="Op_Psi")
      WRITE(out_unit,*) "--- End resulting statevector computed by MolecCav_Action_operator_ND_R1_real"
    END IF
  
    IF (Verbose_local > 25) WRITE(out_unit,*) 
    IF (Verbose_local > 25) WRITE(out_unit,*) "----------------------------------------ACTION OF THE ND OPERATOR OVER THE R1 WF&
                                              & COMPUTED---------------------------------------"; FLUSH(out_unit)
  
  END SUBROUTINE MolecCav_Action_operator_ND_R1_real

  
  SUBROUTINE MolecCav_Action_operator_ND_R1_complex(Op_psi, OpND, Psi, Verbose, Debug) ! Psi is ND AND R1
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    USE Cavity_mode_m
    USE Matter_mode_m
    IMPLICIT NONE

    complex(kind=Rkind),  intent(inout) :: Op_psi(:)
    TYPE(Operator_ND_t),  intent(in)    :: OpND
    complex(kind=Rkind),  intent(in)    :: Psi(:)
    integer, optional,    intent(in)    :: Verbose                                                                              ! cf. comments in HO1D_parameters_m
    logical, optional,    intent(in)    :: Debug                                                                                ! cf. comments in HO1D_parameters_m

    integer                             :: N_mat, N_cav, NB, N1, N2, N3, i_mode, i_1, i_3 ! /!\ N_mat, N_cav are not the basis sets sizes but the respective number of DOF of the matter and cavity subsystems /!\
    integer,             allocatable    :: Ranks_sizes(:)
    complex(kind=Rkind), allocatable    :: Cube(:,:,:), Op_cube(:,:,:)
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
      !CALL Write(OpND)
      IF (SIZE(OpND%tab_indexes_mat_op)==0) WRITE(out_unit,*) "tab mat op : $\empty$"
      IF (SIZE(OpND%tab_indexes_mat_op)/=0) WRITE(out_unit,*) "tab mat op : "//TO_string(OpND%tab_indexes_mat_op)
      IF (SIZE(OpND%tab_indexes_cav_op)==0) WRITE(out_unit,*) "tab cav op : $\empty$"
      IF (SIZE(OpND%tab_indexes_cav_op)/=0) WRITE(out_unit,*) "tab cav op : "//TO_string(OpND%tab_indexes_cav_op)
      WRITE(out_unit,*) "The <<Psi>> argument : "
      CALL Write_Vec(Psi, out_unit, 1, info="Psi")
      WRITE(out_unit,*) "The size of its vector : "//TO_string(Size(Psi))
      WRITE(out_unit,*) "--- End arguments of MolecCav_Action_operator_ND"
      FLUSH(out_unit)
    END IF
    
    !-----------------------------------------------------Checking dimensions----------------------------------------------------
    ! THE DIMENSIONS OF EACH 1D MATMUL WILL BE TESTED IN THE ACTIONS CODED IN ELEM_OP_M !
    N_mat = SIZE(OpND%tab_indexes_mat_op)
    N_cav = SIZE(OpND%tab_indexes_cav_op)
    ALLOCATE(Ranks_sizes(N_mat + N_cav))

    IF (N_mat /= 0) THEN
      Ranks_sizes(1:N_mat) = tab_mat_ops(1:N_mat)%Nb
    END IF 
    IF (N_cav /= 0) THEN
      Ranks_sizes(1+N_mat:N_cav+N_mat) = tab_cav_ops(1:N_cav)%Nb
    END IF 
    NB = PRODUCT(Ranks_sizes) ! better not to have to calculate again each time needed

    IF (SIZE(Psi, dim=1) /= NB) THEN
      WRITE(out_unit,*) "### The dimension of the wavevector Psi does not match the dimension product of the operators to"
      WRITE(out_unit,*) "   Size(Psi, dim=1) = "//TO_string(Size(Psi, dim=1))//"; PRODUCT(Dims) = "//TO_string(NB)
      STOP "### The dimension of Psi does not match the dimension product of the operators."
    END IF

    !----------------------------Computation using reshape----------------------------------    
      !-----------------------Initialization befor the first loop---------------------------    
    N1 = 1
    N2 = Ranks_sizes(1)
    N3 = NB / N2

    ALLOCATE(Cube(   N1, N2, N3))
    ALLOCATE(Op_cube(N1, N2, N3))

    Cube = RESHAPE(Psi, [N1, N2, N3])


    DO i_mode = 1, N_mat + N_cav
      IF (Debug_local) WRITE(out_unit,*)
      IF (Debug_local) WRITE(out_unit,*) "--- i_mode = "//TO_string(i_mode)
      IF (Debug_local) WRITE(out_unit,*)
      IF (Debug_local) WRITE(out_unit,*) "N1, N2, N3 = "//TO_string(N1)//", "//TO_string(N2)//", "//TO_string(N3)

      !-----------------------Action-----------------------------    
      IF (i_mode <= N_mat) THEN
        DO i_3 = 1, N3
          DO i_1 = 1, N1
            CALL Action(Op_psi=Op_cube(i_1,:,i_3), Matmode=tab_mat_ops(i_mode), i_op=OpND%tab_indexes_mat_op(i_mode), Psi=Cube(i_1,&
            &:,i_3), Verbose=Verbose, Debug=Debug_local)
          END DO 
        END DO 
      ELSE 
        DO i_3 = 1, N3
          DO i_1 = 1, N1
            CALL Action(Op_psi=Op_cube(i_1,:,i_3), CavMode=tab_cav_ops(i_mode-N_mat), i_op=OpND%tab_indexes_cav_op(i_mode-N_mat),&
            & Psi=Cube(i_1,:,i_3), Verbose=Verbose, Debug=Debug_local)
          END DO 
        END DO 
        IF (i_mode == N_mat + N_cav) EXIT
      END IF
  
      !-----------------------Reinitialization for next loop-----------------------------    
      DEALLOCATE(Cube)
      N1 = N1 * N2
      N2 = Ranks_sizes(i_mode + 1)
      N3 = N3 / N2
      ALLOCATE(Cube(N1,N2,N3))
      Cube = RESHAPE(Op_cube, [N1, N2, N3])
      DEALLOCATE(Op_cube)
      ALLOCATE(Op_cube(N1,N2,N3))
    END DO 

      !-----------------------Recovering Op_psi (R1)-----------------------------    
    Op_psi = RESHAPE(Op_cube, [NB])
  
    !--------------Conclusion----------------
    IF (Verbose_local > 0) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "--- Resulting statevector from the action of the ND Operator on the Psi statevector operand, computed &
                        &by MolecCav_Action_operator_ND_R1_complex :"
      CALL Write_Vec(Op_psi, out_unit, 1, info="Op_Psi")
      WRITE(out_unit,*) "--- End resulting statevector computed by MolecCav_Action_operator_ND_R1_complex"
    END IF
  
    IF (Verbose_local > 25) WRITE(out_unit,*) 
    IF (Verbose_local > 25) WRITE(out_unit,*) "----------------------------------------ACTION OF THE HO1D OPERATOR OVER THE R1 WF&
                                              & COMPUTED---------------------------------------"; FLUSH(out_unit)
  
  END SUBROUTINE MolecCav_Action_operator_ND_R1_complex

  
  SUBROUTINE MolecCav_Get_OpND_parameter_integer(Parameter_value, Parameter_name, Subsystem, i_mode)
    USE QDUtil_m
    USE Matter_mode_m
    USE Cavity_mode_m
    IMPLICIT NONE 

    integer,             intent(inout) :: Parameter_value                                                                            ! the current values of the indexes for each dimension
    character(len=*),    intent(in)    :: Parameter_name
    character(len=*),    intent(in)    :: Subsystem
    integer,             intent(in)    :: i_mode ! the number of the mode INSIDE the subsystem, not among the total system modes ! (the number of the 2nd cavity mode is 2, not 2+the number of matter modes)

    IF (.NOT. ALLOCATED(tab_mat_ops) .OR. .NOT. ALLOCATED(tab_cav_ops)) THEN
      WRITE(out_unit,*) "### The tab_<subsystem>_ops are not allocated, and the Get procedure can be used only after the tab_<sub&
      &system>_ops have been initialised. Please CALL Initialize an OpND before using Get."
      STOP "### The Get procedure cannot be used as long as the tab_<subsystem>_ops are not allocated."
    END IF

    IF (TO_lowercase(TRIM(Subsystem)) == "matter") THEN
      CALL Get(Parameter_value, tab_mat_ops(i_mode), Parameter_name)

    ELSE IF (TO_lowercase(TRIM(Subsystem)) == "cavity") THEN
      CALL Get(Parameter_value, tab_cav_ops(i_mode), Parameter_name)
    
    ELSE 
      WRITE(out_unit,*) "### No Subsystem name recognized, please verify the input of Get_OpND_parameter_integer subroutine"
      WRITE(out_unit,*) "    Expected : Matter or Cavity ; Subsystem = "//Subsystem
      STOP "### No Operator type recognized, please verify the input of Get_OpND_parameter_integer subroutine"
    
    END IF

  END SUBROUTINE MolecCav_Get_OpND_parameter_integer


  SUBROUTINE MolecCav_Get_OpND_parameter_real(Parameter_value, Parameter_name, Subsystem, i_mode)
    USE QDUtil_m
    USE Matter_mode_m
    USE Cavity_mode_m
    IMPLICIT NONE 

    real(kind=Rkind),    intent(inout) :: Parameter_value                                                                            ! the current values of the indexes for each dimension
    character(len=*),    intent(in)    :: Parameter_name
    character(len=*),    intent(in)    :: Subsystem
    integer,             intent(in)    :: i_mode

    IF (.NOT. ALLOCATED(tab_mat_ops) .OR. .NOT. ALLOCATED(tab_cav_ops)) THEN
      WRITE(out_unit,*) "### The tab_<subsystem>_ops are not allocated, and the Get procedure can be used only after the tab_<sub&
      &system>_ops have been initialised. Please CALL Initialize an OpND before using Get."
      STOP "### The Get procedure cannot be used as long as the tab_<subsystem>_ops are not allocated."
    END IF
    
    IF (TO_lowercase(TRIM(Subsystem)) == "matter") THEN
      CALL Get(Parameter_value, tab_mat_ops(i_mode), Parameter_name)

    ELSE IF (TO_lowercase(TRIM(Subsystem)) == "cavity") THEN
      CALL Get(Parameter_value, tab_cav_ops(i_mode), Parameter_name)
    
    ELSE 
      WRITE(out_unit,*) "### No Subsystem name recognized, please verify the input of Get_OpND_parameter_integer subroutine"
      WRITE(out_unit,*) "    Expected : Matter or Cavity ; Subsystem = "//Subsystem
      STOP "### No Operator type recognized, please verify the input of Get_OpND_parameter_integer subroutine"
    
    END IF

  END SUBROUTINE MolecCav_Get_OpND_parameter_real


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
      IF (SIZE(OpND%tab_indexes_mat_op) /= 0) THEN
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
    IF (ALLOCATED(OpND%tab_indexes_cav_op)) THEN
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


  SUBROUTINE MolecCav_Deallocate_operator_ND(OpND, Dealloc_all, Verbose, Debug)
    USE QDUtil_m
    USE Cavity_mode_m
    USE Matter_mode_m
    IMPLICIT NONE 

    TYPE(Operator_ND_t), intent(inout) :: OpND
    logical, optional,   intent(in)    :: Dealloc_all !shall be improved later : separate sub to dealloc tab_subsys_ops, not to need to target an any opnd to access it
    integer, optional,   intent(in)    :: Verbose                                                                                 ! cf. comments in HO1D_parameters_m
    logical, optional,   intent(in)    :: Debug                                                                                   ! cf. comments in HO1D_parameters_m

    integer                            :: i_op
    logical                            :: Dealloc_all_local
    integer                            :: Verbose_local                                                                      ! goes from 25 (= 0 verbose) to 29 (= maximum verbose) at this layer
    logical                            :: Debug_local

    !------------------------------------------------------Debugging options-----------------------------------------------------
    IF (PRESENT(Dealloc_all)) THEN; Dealloc_all_local = Dealloc_all
    ELSE; Dealloc_all_local = .FALSE.; END IF 
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

    IF (Dealloc_all_local .AND. ALLOCATED(tab_mat_ops)) THEN
      DO i_op = 1, SIZE(tab_mat_ops)
        CALL Dealloc(tab_mat_ops(i_op), Verbose=Verbose_local, Debug=Debug_local)
      END DO
      DEALLOCATE(tab_mat_ops)
    END IF 
    IF (Dealloc_all_local .AND. ALLOCATED(tab_cav_ops)) THEN
      DO i_op = 1, SIZE(tab_cav_ops)
        CALL Dealloc(tab_cav_ops(i_op), Verbose=Verbose_local, Debug=Debug_local)
      END DO
      DEALLOCATE(tab_cav_ops)
      WRITE(out_unit,*) "### WARNING : the tab_mat_ops and the tab_cav_ops tables have been deallocated. WARNING ###"
    END IF 

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "--- The OpND object after having been deallocated :"
      CALL Write(OpND)
      WRITE(out_unit,*) "--- End dellocating OpND"
    END IF

  END SUBROUTINE MolecCav_Deallocate_operator_ND


END MODULE