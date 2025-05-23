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
! Initialize_quantum_HO1D : reads the namelist and initialize the type, then constructs the operat-
! or using parameters of the HO1D_para object from the so called derived type.
! Append_quantum_HO1D     : add an HO operator to a already initialized object of Matter_mode_t type.
! Write_quantum_HO1D      : display values of the type in the output
! Deallocate_quantum_HO1D : deallocate all tables of the type
! The module to initialize the HO by reading its parameters from the namelist.  
! Read_HO1D_parameters  : reads the namelist and initialize the type.
! Write_HO1D_parameters : displays values of the type in the output.
!==================================================================================================
!==================================================================================================
MODULE Matter_mode_m
  !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT, real64
  USE QDUtil_m                                                                 ! gives Rkind=real64; out_unit=OUTPUT_UNIT; INPUT_UNIT=in_unit; EYE=i and other numbers; TO_LOWERCASE; TO_UPPERCASE;... We thereby use ZERO instead of 0.0_real64
  USE Quantum_HO1D_m
  IMPLICIT NONE


  TYPE, EXTENDS(Quantum_HO1D_t)  :: Matter_mode_t
    real(kind=Rkind)             :: lambda = -ONE
    real(kind=Rkind)             :: CoeffsDipMomt(5) !/!\ if want to make it allocatable and move the allocation in the initialize sub, then mind changing the dealloc sub with a deallocate cmd ! (it is not supposed to change anything else in this file)                                                                         ! size of the Taylor expansion arbitrarily decided /hard-coded here 
  END TYPE


  PRIVATE 
  
  PUBLIC Matter_mode_t, Initialize, Action, Get, Write, Dealloc

  INTERFACE Initialize
    MODULE PROCEDURE MolecCav_Initialize_matter_mode
  END INTERFACE
  INTERFACE Initialize_mu
    MODULE PROCEDURE MolecCav_Initialize_mu_matter_mode
  END INTERFACE
  INTERFACE Action
    MODULE PROCEDURE MolecCav_Action_matter_mode_R1_real, MolecCav_Action_matter_R1_complex
  END INTERFACE
  INTERFACE Get
    MODULE PROCEDURE MolecCav_Get_MatMode_parameter_integer, MolecCav_Get_MatMode_parameter_real
  END INTERFACE
  INTERFACE Write
    MODULE PROCEDURE MolecCav_Write_matter_mode
  END INTERFACE
  INTERFACE Dealloc
    MODULE PROCEDURE MolecCav_Deallocate_matter_mode
  END INTERFACE
    

  CONTAINS


  SUBROUTINE MolecCav_Initialize_matter_mode(MatMode, nio, Dense, Verbose, Debug) ! here init on the 1D QHO basis : don't need Nq etc. will write later Init_grid or Init_other_basis, maybe called by this one, in this case, the "call" will be determined by the optional arg (Nb, Nq etc.)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
    USE QDUtil_m
    USE Quantum_HO1D_m
    IMPLICIT NONE
  
    TYPE(Matter_mode_t), intent(inout)  :: MatMode
    integer,              intent(in)    :: nio
    logical, optional,    intent(in)    :: Dense                                                                         ! cf. comments in HO1D_parameters_m
    integer, optional,    intent(in)    :: Verbose                                                                         ! cf. comments in HO1D_parameters_m
    logical, optional,    intent(in)    :: Debug                                                                           ! cf. comments in HO1D_parameters_m

    real(kind=Rkind)                    :: lambda
    real(kind=Rkind), allocatable       :: CoeffsDipMomt(:)
    real(kind=Rkind)                    :: w, m, Eq_pos, Scale_q
    integer                             :: Nb, Nq
    integer                             :: err_io
    logical                             :: Dense_local                                                              ! goes from 20 (= 0 verbose) to 24 (= maximum verbose) at this layer
    integer                             :: Verbose_local                                                              ! goes from 20 (= 0 verbose) to 24 (= maximum verbose) at this layer
    logical                             :: Debug_local

    NAMELIST /Matter_mode/ lambda, CoeffsDipMomt, Nb, w, m, Nq, Eq_pos, Scale_q

    !------------------------------------------------------Debugging options-----------------------------------------------------
    IF (PRESENT(Verbose)) THEN; Verbose_local = Verbose
    ELSE; Verbose_local = 20; END IF 
    IF (PRESENT(Debug))   THEN; Debug_local   = Debug
    ELSE; Debug_local = .FALSE.; END IF

    IF (Verbose_local > 20) WRITE(out_unit,*) 
    IF (Verbose_local > 20) WRITE(out_unit,*) "-------------------------------------------------INITIALIZING THE MATTER MODE OBJE&
                                              &CT-------------------------------------------------"; FLUSH(out_unit)

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "--- Arguments of MolecCav_Initialize_matter_mode :"
      WRITE(out_unit,*) "The <<MatMode>> argument :"
      CALL Write(MatMode)
      WRITE(out_unit,*) "The <<nio>> argument :"//TO_string(nio)
      IF (PRESENT(Dense)) WRITE(out_unit,*) "The <<Dense>> argument : "//TO_string(Dense)
      WRITE(out_unit,*) "--- End arguments of MolecCav_Initialize_matter_mode"
      FLUSH(out_unit)
    END IF
    
    !----------------------Initialization to default values--------------------
    lambda        = -ONE
    ALLOCATE(CoeffsDipMomt(SIZE(MatMode%CoeffsDipMomt)))
    CoeffsDipMomt = ZERO
    Nb            = 0
    w             = ZERO
    m             = ZERO
    Nq            = 0
    Eq_pos        = -ONE
    Scale_q       = HUGE(ONE)

    !------------------------------------------Initializing the 1D QHO associated to the matter mode------------------------------------------
      !------------------------------Reading of the nml--------------------------
    WRITE(out_unit,*) 
    WRITE(out_unit,*) '********************************************************************************'
    WRITE(out_unit,*) '**************************** READING THE QHO1D PARAMETERS ***************************'
    WRITE(out_unit,*) '********************************************************************************'
    
    READ(nio, nml = Matter_mode, iostat = err_io)                                     ! assign the values read in the nml to the declared list of parameters

    IF (Debug) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "-----------------------The namelist parameters are read as----------------------"
      WRITE(out_unit, nml = Matter_mode)
      WRITE(out_unit,*) "-------------------------End of the namelist parameters-------------------------"
    END IF
    
      !------------------------------Check reading error-------------------------
    IF(err_io < 0) THEN
      WRITE(out_unit,*) ''
      WRITE(out_unit,*) '#######################################################'
      WRITE(out_unit,*) '######## Error in Initialize_matter_mode (err_io<0) #########'
      WRITE(out_unit,*) '#######################################################'
      WRITE(out_unit,*) '################ err_io = ', err_io, '################'
      STOP '################# Check basis data ################'

    ELSE IF( err_io > 0) THEN
      WRITE(out_unit,*) ''
      WRITE(out_unit,*) '#######################################################'
      WRITE(out_unit,*) '######## Error in Initialize_matter_mode (err_io>0) #########'
      WRITE(out_unit,*) '#######################################################'
      WRITE(out_unit,*) '################ err_io = ', err_io, '################'
      STOP '################# Check basis data ################'

    END IF

    IF (Nb == 0) THEN
      WRITE(out_unit,*) "The number of basis vector associated to any HO CANNOT be 0 (what are are you going to study if there is&
                       & no system ???). Please check the data file '.nml'"
      STOP "### The number of basis vector associated to any HO CANNOT be 0 (what are are you going to study if there is&
                       & no system ???). Please check the data file '.nml'"
    END IF
    
      !---------------Construction of the Quantum_HO1D_t type object-----------
    CALL Initialize(MatMode%Quantum_HO1D_t, Nb, w, m, Nb_op=5, Dense=Dense_local, Verbose=Verbose_local, Debug=Debug_local)
    !### TEMPORARY ### (until we implement the grid)
    MatMode%Nq      = Nq ! so far they are not affected in the Initialize_quantum_HO1D, which takes care of representations upon the basis, later these 3 lines will be in a sub in QHO1D_m
    MatMode%Eq_pos  = Eq_pos
    MatMode%Scale_q = Scale_q
    !#################

    WRITE(out_unit,*) 
    WRITE(out_unit,*) '********************************************************************************'
    WRITE(out_unit,*) '************************** QHO1D CONSTRUCTED *************************'
    WRITE(out_unit,*) '********************************************************************************'

    IF (Debug) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "--------------Matter mode constructed by MolecCav_Initialize_matter_mode--------------"
      CALL Write(MatMode)
      WRITE(out_unit,*) "------------End Matter mode constructed by MolecCav_Initialize_matter_mode------------"
    END IF

    !------------------------------------------Completing the matter mode------------------------------------------
    MatMode%lambda        = lambda
    MatMode%CoeffsDipMomt = CoeffsDipMomt
    
    CALL Initialize_mu(MatMode, Verbose=Verbose_local, Debug=Debug_local) ! cannot be non dense ! All MatMode is passed and not anly Mat%Tab(4) because its other parameters are needed for the construction

    IF (Verbose_local > 20) WRITE(out_unit,*) 
    IF (Verbose_local > 20) WRITE(out_unit,*) "--------------------------------------------------QUANTUM HO1D OBJECT INITIALIZED-&
                                              &------------------------------------------------"; FLUSH(out_unit)

  END SUBROUTINE MolecCav_Initialize_matter_mode


  SUBROUTINE MolecCav_Initialize_mu_matter_mode(MatMode, Verbose, Debug)                         ! a matchup of the Initialize_elem_op and the Initialize_x from QHO1D_m. Stays here rather than the latter module for the sake of the consistency (a HO do not have a dipole moment by itself)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
    USE QDUtil_m
    USE Quantum_HO1D_m
    IMPLICIT NONE

    TYPE(Matter_mode_t),        intent(inout) :: MatMode                                                                         ! the object of type Elem_op_t to be constructed here
    integer,          optional, intent(in)    :: Verbose                                                                         ! cf. comments in HO1D_parameters_m
    logical,          optional, intent(in)    :: Debug                                                                           ! cf. comments in HO1D_parameters_m

    TYPE(Quantum_HO1D_t)                      :: QHO1D_dense
    integer                                   :: i
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
      WRITE(out_unit,*) "--- Arguments of MolecCav_Initialize_mu_matter_mode :"
      WRITE(out_unit,*) "The <<MatMode>> argument :"
      CALL Write(MatMode)
      WRITE(out_unit,*) "--- End arguments of MolecCav_Initialize_mu_matter_mode"
      FLUSH(out_unit)
    END IF
    
    !---------------------------------------First steps of the construction of the Operator--------------------------------------
    ALLOCATE(character(len=LEN_TRIM("DipMomt")) :: MatMode%Tab_op(4)%Operator_type)                                                   ! /!\ strings cannot be allocated the exact same way as tables ! /!\
    MatMode%Tab_op(4)%Operator_type = TO_lowercase(TRIM("DipMomt")) ! yes it is useless here to lowercase/trim instead of just writing "dipmomt", but can be convenient if want to change the name and do "research/replace". Allocation on assignement (not anymore : supposed to work but caused dynamic allocation random errors at execution). Elem_op_type has the right lengths (no spaces added) thanks to len=* at declaration and it will fit the Op%op_type thanks to len=:, allocatable at declaration of the derived type. 
    MatMode%Tab_op(4)%Dense         = .TRUE.

!---------------------------------------------Construction of the matrix Operator--------------------------------------------
    IF (Debug_local) THEN
      WRITE(out_unit,*); WRITE(out_unit,*) "--- The MatMode object just before construction of its DipMomt matrix's representation"
      CALL Write(MatMode)
      WRITE(out_unit,*) "--- End Elem_op_t object (just before construction of its DipMomt's matrix representation)"
    END IF 

    CALL Initialize(QHO1D_dense, MatMode%Nb, MatMode%w, MatMode%m, Dense=.TRUE., Verbose=Verbose_local, Debug=Debug_local)
    ALLOCATE(MatMode%Tab_op(4)%Dense_val(MatMode%Nb, MatMode%Nb))

    MatMode%Tab_op(4)%Dense_val = MatMode%CoeffsDipMomt(1) * QHO1D_dense%Tab_op(0)%Dense_val
    DO i = 1, SIZE(MatMode%CoeffsDipMomt) - 1
      QHO1D_dense%Tab_op(0)%Dense_val = MATMUL(QHO1D_dense%Tab_op(2)%Dense_val, QHO1D_dense%Tab_op(0)%Dense_val)
      MatMode%Tab_op(4)%Dense_val = MatMode%Tab_op(4)%Dense_val + MatMode%CoeffsDipMomt(i+1) * QHO1D_dense%Tab_op(0)%Dense_val
    END DO
    
    IF (Verbose_local > 26) THEN
      IF (Verbose_local < 28) WRITE(out_unit,*)
      WRITE(out_unit,*) "--- HO1D operator constructed by MolecCav_Initialize_HO1D_operator :"
      CALL Write(MatMode)
      WRITE(out_unit,*) "--- End HO1D operator constructed by MolecCav_Initialize_HO1D_operator"
    END IF

    IF (Verbose_local > 25) WRITE(out_unit,*) 
    IF (Verbose_local > 25) WRITE(out_unit,*) "-----------------------------------------------------HO1D OPERATOR INITIALIZED----&
                                              &------------------------------------------------"; FLUSH(out_unit)

  END SUBROUTINE MolecCav_Initialize_mu_matter_mode


  SUBROUTINE MolecCav_Action_matter_mode_R1_real(Op_psi, MatMode, i_op, Psi, Verbose, Debug)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    USE Quantum_HO1D_m
    IMPLICIT NONE

    real(kind=Rkind),    intent(inout) :: Op_psi(:)
    TYPE(Matter_mode_t), intent(in)    :: MatMode
    integer,             intent(in)    :: i_op
    real(kind=Rkind),    intent(in)    :: Psi(:)
    integer, optional,   intent(in)    :: Verbose                                                                              ! cf. comments in HO1D_parameters_m
    logical, optional,   intent(in)    :: Debug                                                                                ! cf. comments in HO1D_parameters_m

    integer                             :: Nb
    integer                             :: Verbose_local                                                                   ! goes from 25 (= 0 verbose) to 29 (= maximum verbose) at this layer
    logical                             :: Debug_local

    !------------------------------------------------------Debugging options-----------------------------------------------------
    IF (PRESENT(Verbose)) THEN; Verbose_local = Verbose
    ELSE; Verbose_local = 20; END IF 
    IF (PRESENT(Debug))   THEN; Debug_local   = Debug
    ELSE; Debug_local = .FALSE.; END IF

    IF (Verbose_local > 25) WRITE(out_unit,*) 
    IF (Verbose_local > 25) WRITE(out_unit,*) "------------------------COMPUTING ACTION OF THE "//TO_string(i_op)//"^{th} OPERATO&
    &R OF THE  MATTER MODE OVER THE R1 WF------------------------"; FLUSH(out_unit)

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "--- Arguments of MolecCav_Action_matter_mode_R1_real :"
      WRITE(out_unit,*) "The <<MatMode>> argument :"
      CALL Write(MatMode)
      WRITE(out_unit,*) "The <<i_op>> argument :"//TO_string(i_op)
      WRITE(out_unit,*) "The <<Psi>> argument : "
      CALL Write_Vec(Psi, out_unit, 1, info="Psi")
      WRITE(out_unit,*) "The size of its vector : "//TO_string(Size(Psi))
      WRITE(out_unit,*) "--- End arguments of MolecCav_Action_matter_mode_R1_real"
      FLUSH(out_unit)
    END IF
    
    !-----------------------------------------------------Checking dimensions----------------------------------------------------
    ! ALREADY CHECKED IN THE ACTIONS CODED IN ELEM_OP_M ! 

    !---------------------------------------------Selection of the calculation method--------------------------------------------
    CALL Action(Op_psi=Op_psi, QHO1D=MatMode%Quantum_HO1D_t, i_op=i_op, Psi=Psi, Verbose=Verbose_local, Debug=Debug_local)

    IF (Verbose_local > 26) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "--- Resulting statevector from the action of the HO1D Elem_op on the Psi statevector operand, computed &
                        &by Action_HO1D_operator_R1 :"
      CALL Write_Vec(Op_psi, out_unit, 1, info="Op_Psi")
      WRITE(out_unit,*) "--- End resulting statevector computed by Action_HO1D_operator_R1"
    END IF
  
    IF (Verbose_local > 25) WRITE(out_unit,*) 
    IF (Verbose_local > 25) WRITE(out_unit,*) "-------------------------ACTION OF THE "//TO_string(i_op)//"^{th} OPERATO&
    &R OF THE  MATTER MODE OVER THE R1 WF COMPUTED------------------------"; FLUSH(out_unit)
  
  END SUBROUTINE MolecCav_Action_matter_mode_R1_real


  SUBROUTINE MolecCav_Action_matter_R1_complex(Op_psi, MatMode, i_op, Psi, Verbose, Debug)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    USE Quantum_HO1D_m
    IMPLICIT NONE

    complex(kind=Rkind), intent(inout) :: Op_psi(:)
    TYPE(Matter_mode_t), intent(in)    :: MatMode
    integer,             intent(in)    :: i_op
    complex(kind=Rkind), intent(in)    :: Psi(:)
    integer, optional,   intent(in)    :: Verbose                                                                              ! cf. comments in HO1D_parameters_m
    logical, optional,   intent(in)    :: Debug                                                                                ! cf. comments in HO1D_parameters_m

    integer                             :: Nb
    integer                             :: Verbose_local                                                                   ! goes from 25 (= 0 verbose) to 29 (= maximum verbose) at this layer
    logical                             :: Debug_local

    !------------------------------------------------------Debugging options-----------------------------------------------------
    IF (PRESENT(Verbose)) THEN; Verbose_local = Verbose
    ELSE; Verbose_local = 20; END IF 
    IF (PRESENT(Debug))   THEN; Debug_local   = Debug
    ELSE; Debug_local = .FALSE.; END IF

    IF (Verbose_local > 25) WRITE(out_unit,*) 
    IF (Verbose_local > 25) WRITE(out_unit,*) "------------------------COMPUTING ACTION OF THE "//TO_string(i_op)//"^{th} OPERATO&
    &R OF THE  MATTER MODE OVER THE R1 WF------------------------"; FLUSH(out_unit)

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "--- Arguments of MolecCav_Action_matter_R1_complex :"
      WRITE(out_unit,*) "The <<MatMode>> argument :"
      CALL Write(MatMode)
      WRITE(out_unit,*) "The <<i_op>> argument :"//TO_string(i_op)
      WRITE(out_unit,*) "The <<Psi>> argument : "
      CALL Write_Vec(Psi, out_unit, 1, info="Psi")
      WRITE(out_unit,*) "The size of its vector : "//TO_string(Size(Psi))
      WRITE(out_unit,*) "--- End arguments of MolecCav_Action_matter_R1_complex"
      FLUSH(out_unit)
    END IF
    
    !-----------------------------------------------------Checking dimensions----------------------------------------------------
    ! ALREADY CHECKED IN THE ACTIONS CODED IN ELEM_OP_M !

    !---------------------------------------------Selection of the calculation method--------------------------------------------
    CALL Action(Op_psi=Op_psi, QHO1D=MatMode%Quantum_HO1D_t, i_op=i_op, Psi=Psi, Verbose=Verbose_local, Debug=Debug_local)

    IF (Verbose_local > 26) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "--- Resulting statevector from the action of the HO1D Elem_op on the Psi statevector operand, computed &
                        &by Action_HO1D_operator_R1 :"
      CALL Write_Vec(Op_psi, out_unit, 1, info="Op_Psi")
      WRITE(out_unit,*) "--- End resulting statevector computed by Action_HO1D_operator_R1"
    END IF
  
    IF (Verbose_local > 25) WRITE(out_unit,*) 
    IF (Verbose_local > 25) WRITE(out_unit,*) "-------------------------ACTION OF THE "//TO_string(i_op)//"^{th} OPERATO&
    &R OF THE  MATTER MODE OVER THE R1 WF COMPUTED------------------------"; FLUSH(out_unit)

  END SUBROUTINE MolecCav_Action_matter_R1_complex

  
  SUBROUTINE MolecCav_Get_MatMode_parameter_integer(Parameter_value, MatMode, Parameter_name)
    USE QDUtil_m
    IMPLICIT NONE 

    integer,             intent(inout) :: Parameter_value                                                                            ! the current values of the indexes for each dimension
    TYPE(Matter_mode_t), intent(in)    :: MatMode
    character(len=*),    intent(in)    :: Parameter_name

    CALL Get(Parameter_value, MatMode%Quantum_HO1D_t, Parameter_name)

  END SUBROUTINE MolecCav_Get_MatMode_parameter_integer


  SUBROUTINE MolecCav_Get_MatMode_parameter_real(Parameter_value, MatMode, Parameter_name)
    USE QDUtil_m
    IMPLICIT NONE 

    real(kind=Rkind),     intent(inout) :: Parameter_value                                                                            ! the current values of the indexes for each dimension
    TYPE(Matter_mode_t), intent(in)    :: MatMode
    character(len=*),     intent(in)    :: Parameter_name

    IF (TO_lowercase(TRIM(Parameter_name)) == "lambda") THEN
      Parameter_value = MatMode%lambda
    ELSE 
      CALL Get(Parameter_value, MatMode%Quantum_HO1D_t, Parameter_name)
    END IF

  END SUBROUTINE MolecCav_Get_MatMode_parameter_real


  SUBROUTINE MolecCav_Write_matter_mode(MatMode)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
    USE QDUtil_m
    IMPLICIT NONE 
    
    TYPE(Matter_mode_t), intent(in) :: MatMode

    integer                          :: i_op

    WRITE(out_unit,*) "___________________________________The parameters of the Matter mode___________"
    WRITE(out_unit,*) "|The associated 1D Quantum HO (MatMode%Quantum_HO1D_t) :                       |"
    CALL Write(MatMode%Quantum_HO1D_t)
    WRITE(out_unit,*) "_______________________________The peculiar parameters to the Matter mode_____________________________"
    WRITE(out_unit,*) "|The strength parameter of its coupling with a cavity mode (MatMode%lambda)    | "//TO_string(MatMode%lambda)
    WRITE(out_unit,*) "|______________________________________________________________________________|______________________"
    WRITE(out_unit,*) "|The n-derivatives of the matter dipole moment with respect to this mode's DOF | "
    WRITE(out_unit,*) "| (i.e. the coefficients of its Taylor expansion) (MatMode%CoeffsDipMomt)      | " 
    CALL Write_Vec(MatMode%CoeffsDipMomt, out_unit, SIZE(MatMode%CoeffsDipMomt), info="CoeffsDipMomt")
    WRITE(out_unit,*) "|_____________________________________End Matter mode parameters_______________|"
    FLUSH(out_unit)

  END SUBROUTINE MolecCav_Write_matter_mode


  SUBROUTINE MolecCav_Deallocate_matter_mode(MatMode, Verbose, Debug)
    USE QDUtil_m
    USE Quantum_HO1D_m
    IMPLICIT NONE 

    TYPE(Matter_mode_t), intent(inout) :: MatMode
    integer, optional,   intent(in)    :: Verbose                                                                                 ! cf. comments in HO1D_parameters_m
    logical, optional,   intent(in)    :: Debug                                                                                   ! cf. comments in HO1D_parameters_m

    integer                            :: i_op
    integer                            :: Verbose_local                                                                      ! goes from 25 (= 0 verbose) to 29 (= maximum verbose) at this layer
    logical                            :: Debug_local

    !------------------------------------------------------Debugging options-----------------------------------------------------
    IF (PRESENT(Verbose)) THEN; Verbose_local = Verbose
    ELSE; Verbose_local = 20; END IF 
    IF (PRESENT(Debug))   THEN; Debug_local   = Debug
    ELSE; Debug_local   = .FALSE.; END IF

    IF (Debug_local) THEN
      WRITE(out_unit,*) "--- The Matter mode to be deallocated :"
      CALL Write(MatMode)
      WRITE(out_unit,*) "--- End Matter mode to be deallocated"
    END IF 

    !-----------------------------Deallocating the HO1D operator object----------------------------
    IF (Verbose_local > 27) WRITE(out_unit,*)
    IF (Verbose_local > 27) WRITE(out_unit,*) "-----------------------------------------------Deallocating the matter mode object&
    &----------------------------------------------"
  
    MatMode%lambda        = -ONE
    MatMode%CoeffsDipMomt  = ZERO
    CALL Dealloc(MatMode%Quantum_HO1D_t, Verbose=Verbose_local, Debug=Debug_local)

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "--- The Matter mode object after having been deallocated :"
      CALL Write(MatMode)
      WRITE(out_unit,*) "--- End dellocating Matter mode"
    END IF

  END SUBROUTINE MolecCav_Deallocate_matter_mode


END MODULE