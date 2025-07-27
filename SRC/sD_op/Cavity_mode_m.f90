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
! Initialize_quantum_HO1D : reads the namelist and initialize the type, then constructs the operat-
! or using parameters of the HO1D_para object from the so called derived type.
! Append_quantum_HO1D     : add an HO operator to a already initialized object of Cavity_mode_t type.
! Write_quantum_HO1D      : display values of the type in the output
! Deallocate_quantum_HO1D : deallocate all tables of the type
! The module to initialize the HO by reading its parameters from the namelist.  
! Read_HO1D_parameters  : reads the namelist and initialize the type.
! Write_HO1D_parameters : displays values of the type in the output.
!==================================================================================================
!==================================================================================================
MODULE Cavity_mode_m
  !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT, real64
  USE QDUtil_m                                                                 ! gives Rkind=real64; out_unit=OUTPUT_UNIT; INPUT_UNIT=in_unit; EYE=i and other numbers; TO_LOWERCASE; TO_UPPERCASE;... We thereby use ZERO instead of 0.0_real64
  USE Quantum_HO1D_m
  IMPLICIT NONE


  TYPE, EXTENDS(Quantum_HO1D_t)  :: Cavity_mode_t
    real(kind=Rkind)             :: lambda        = -ONE
  END TYPE


  PRIVATE 
  
  PUBLIC Cavity_mode_t, Initialize, Action, Get, Write, Dealloc

  INTERFACE Initialize
    MODULE PROCEDURE MolecCav_Initialize_Cavity_mode
  END INTERFACE
  INTERFACE Action
    MODULE PROCEDURE MolecCav_Action_Cavity_mode_R1_real, MolecCav_Action_Cavity_R1_complex
  END INTERFACE
  INTERFACE Get
    MODULE PROCEDURE MolecCav_Get_CavMode_parameter_integer, MolecCav_Get_CavMOde_parameter_real
  END INTERFACE
  INTERFACE Write
    MODULE PROCEDURE MolecCav_Write_Cavity_mode
  END INTERFACE
  INTERFACE Dealloc
    MODULE PROCEDURE MolecCav_Deallocate_Cavity_mode
  END INTERFACE
    

  CONTAINS


  SUBROUTINE MolecCav_Initialize_Cavity_mode(CavMode, nio, Dense, Verbose, Debug) ! here init on the 1D QHO basis : don't need Nq etc. will write later Init_grid or Init_other_basis, maybe called by this one, in this case, the "call" will be determined by the optional arg (Nb, Nq etc.)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
    USE QDUtil_m
    USE Quantum_HO1D_m
    IMPLICIT NONE
  
    TYPE(Cavity_mode_t), intent(inout)  :: CavMode
    integer,              intent(in)    :: nio
    logical, optional,    intent(in)    :: Dense                                                                              ! cf. comments in HO1D_parameters_m
    integer, optional,    intent(in)    :: Verbose                                                                         ! cf. comments in HO1D_parameters_m
    logical, optional,    intent(in)    :: Debug                                                                           ! cf. comments in HO1D_parameters_m

    real(kind=Rkind)                    :: lambda
    real(kind=Rkind)                    :: w, m, Eq_pos, Scale_q
    integer                             :: Nb, Nq
    integer                             :: err_io
    logical                             :: Dense_local                                                              ! goes from 16 (= 0 verbose) to 24 (= maximum verbose) at this layer
    integer                             :: Verbose_local                                                              ! goes from 16 (= 0 verbose) to 24 (= maximum verbose) at this layer
    logical                             :: Debug_local

    NAMELIST /Cavity_mode/ lambda, Nb, w, m, Nq, Eq_pos, Scale_q

    !------------------------------------------------------Debugging options-----------------------------------------------------
    IF (PRESENT(Verbose)) THEN; Verbose_local = Verbose
    ELSE; Verbose_local = 16; END IF 
    IF (PRESENT(Debug))   THEN; Debug_local   = Debug
    ELSE; Debug_local = .FALSE.; END IF
    IF (Debug_local) Verbose_local = 24

    IF (Verbose_local > 16) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "o o o Arguments of MolecCav_Initialize_Cavity_mode :"
      WRITE(out_unit,*) "The <<CavMode>> argument :"
      CALL Write(CavMode)
      WRITE(out_unit,*) "The <<nio>> argument     : "//TO_string(nio)
      IF (PRESENT(Dense)) WRITE(out_unit,*) "The <<Dense>> argument   : "//TO_string(Dense)
      FLUSH(out_unit)
    END IF
    
    !----------------------Initialization to default values--------------------
    lambda        = -ONE
    Nb            = 0
    w             = ZERO
    m             = ZERO
    Nq            = 0
    Eq_pos        = -ONE
    Scale_q       = HUGE(ONE)

    !------------------------------------------Initializing the 1D QHO associated to the Cavity mode------------------------------------------    
    READ(nio, nml = Cavity_mode, iostat = err_io)                                     ! assign the values read in the nml to the declared list of parameters

    IF (Debug) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "--- The namelist parameters read :"
      WRITE(out_unit, nml = Cavity_mode)
    END IF
    
      !------------------------------Check reading error-------------------------
    IF(err_io < 0) THEN
      WRITE(out_unit,*) ''
      WRITE(out_unit,*) '#######################################################'
      WRITE(out_unit,*) '######## Error in Initialize_Cavity_mode (err_io<0) #########'
      WRITE(out_unit,*) '#######################################################'
      WRITE(out_unit,*) '################ err_io = ', err_io, '################'
      STOP '################# Check basis data ################'

    ELSE IF( err_io > 0) THEN
      WRITE(out_unit,*) ''
      WRITE(out_unit,*) '#######################################################'
      WRITE(out_unit,*) '######## Error in Initialize_Cavity_mode (err_io>0) #########'
      WRITE(out_unit,*) '#######################################################'
      WRITE(out_unit,*) '################ err_io = ', err_io, '################'
      STOP '################# Check basis data ################'

    END IF

    IF (Nb == 0) THEN
      WRITE(out_unit,*) "### The number of basis vector associated to any HO CANNOT be 0 (what are are you going to study if ther&
      &e is no system ???). Please check the data file '.nml'"
      STOP "### The number of basis vector associated to any HO CANNOT be 0 (what are are you going to study if there is no syste&
      &m ???). Please check the data file '.nml'"
    END IF
    
      !---------------Construction of the Quantum_HO1D_t type object-----------
    IF (PRESENT(Dense)) THEN; Dense_local = Dense
    ELSE; Dense_local = .FALSE.; END IF

    CALL Initialize(CavMode%Quantum_HO1D_t, Nb, w, m, Nb_op=4, Dense=Dense_local, Verbose=Verbose_local, Debug=Debug_local)
    !### TEMPORARY ###
    CavMode%Nq      = Nq ! so far they are not affected in the Initialize_quantum_HO1D, which takes care of representations upon the basis, later these 3 lines will be in a sub in QHO1D_m
    CavMode%Eq_pos  = Eq_pos
    CavMode%Scale_q = Scale_q
    !#################

    IF (Debug) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "--- Cavity mode constructed by MolecCav_Initialize_Cavity_mode"
      CALL Write(CavMode)
    END IF

    !------------------------------------------Completing the Cavity mode------------------------------------------
    CavMode%lambda        = lambda
    
  END SUBROUTINE MolecCav_Initialize_Cavity_mode


  SUBROUTINE MolecCav_Action_Cavity_mode_R1_real(Op_psi, CavMode, i_op, Psi, Verbose, Debug)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    USE Quantum_HO1D_m
    IMPLICIT NONE

    real(kind=Rkind),    intent(inout) :: Op_psi(:)
    TYPE(Cavity_mode_t), intent(in)    :: CavMode
    integer,             intent(in)    :: i_op
    real(kind=Rkind),    intent(in)    :: Psi(:)
    integer, optional,   intent(in)    :: Verbose                                                                              ! cf. comments in HO1D_parameters_m
    logical, optional,   intent(in)    :: Debug                                                                                ! cf. comments in HO1D_parameters_m

    integer                            :: Nb
    integer                            :: Verbose_local                                                                   ! goes from 25 (= 0 verbose) to 29 (= maximum verbose) at this layer
    logical                            :: Debug_local

    !------------------------------------------------------Debugging options-----------------------------------------------------
    IF (PRESENT(Verbose)) THEN; Verbose_local = Verbose
    ELSE; Verbose_local = 16; END IF 
    IF (PRESENT(Debug))   THEN; Debug_local   = Debug
    ELSE; Debug_local = .FALSE.; END IF
    IF (Debug_local) Verbose_local = 24

    IF (Verbose_local > 16) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "o o Arguments of MolecCav_Action_Cavity_mode_R1_real :"
      WRITE(out_unit,*) "The <<CavMode>> argument :"
      CALL Write(CavMode)
      WRITE(out_unit,*) "The <<i_op>> argument    : "//TO_string(i_op)
      WRITE(out_unit,*) "The <<Psi>> argument     :"
      CALL Write_Vec(Psi, out_unit, 1, info="Psi")
      WRITE(out_unit,*) "The size of its vector   : "//TO_string(Size(Psi))
      FLUSH(out_unit)
    END IF
    
    !-----------------------------------------------------Checking dimensions----------------------------------------------------
    ! ALREADY CHECKED IN THE ACTIONS CODED IN ELEM_OP_M ! 

    !---------------------------------------------Selection of the calculation method--------------------------------------------
    CALL Action(Op_psi=Op_psi, QHO1D=CavMode%Quantum_HO1D_t, i_op=i_op, Psi=Psi, Verbose=Verbose_local, Debug=Debug_local)

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "--- Resulting statevector from the action of the HO1D Elem_op on the Psi statevector operand, computed &
                        &by Action_HO1D_operator_R1 :"
      CALL Write_Vec(Op_psi, out_unit, 1, info="Op_Psi")
      WRITE(out_unit,*) "--- End resulting statevector computed by Action_HO1D_operator_R1"
    END IF
    
  END SUBROUTINE MolecCav_Action_Cavity_mode_R1_real


  SUBROUTINE MolecCav_Action_Cavity_R1_complex(Op_psi, CavMode, i_op, Psi, Verbose, Debug)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    USE Quantum_HO1D_m
    IMPLICIT NONE

    complex(kind=Rkind), intent(inout) :: Op_psi(:)
    TYPE(Cavity_mode_t), intent(in)    :: CavMode
    integer,             intent(in)    :: i_op
    complex(kind=Rkind), intent(in)    :: Psi(:)
    integer, optional,   intent(in)    :: Verbose                                                                              ! cf. comments in HO1D_parameters_m
    logical, optional,   intent(in)    :: Debug                                                                                ! cf. comments in HO1D_parameters_m

    integer                             :: Nb
    integer                             :: Verbose_local                                                                   ! goes from 25 (= 0 verbose) to 29 (= maximum verbose) at this layer
    logical                             :: Debug_local

    !------------------------------------------------------Debugging options-----------------------------------------------------
    IF (PRESENT(Verbose)) THEN; Verbose_local = Verbose
    ELSE; Verbose_local = 16; END IF 
    IF (PRESENT(Debug))   THEN; Debug_local   = Debug
    ELSE; Debug_local = .FALSE.; END IF
    IF (Debug_local) Verbose_local = 24

    IF (Verbose_local > 16) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "o o Arguments of MolecCav_Action_Cavity_R1_complex :"
      WRITE(out_unit,*) "The <<CavMode>> argument :"
      CALL Write(CavMode)
      WRITE(out_unit,*) "The <<i_op>> argument    : "//TO_string(i_op)
      WRITE(out_unit,*) "The <<Psi>> argument     :"
      CALL Write_Vec(Psi, out_unit, 1, info="Psi")
      WRITE(out_unit,*) "The size of its vector   :"//TO_string(Size(Psi))
      FLUSH(out_unit)
    END IF
    
    !-----------------------------------------------------Checking dimensions----------------------------------------------------
    ! ALREADY CHECKED IN THE ACTIONS CODED IN ELEM_OP_M !

    !---------------------------------------------Selection of the calculation method--------------------------------------------
    CALL Action(Op_psi=Op_psi, QHO1D=CavMode%Quantum_HO1D_t, i_op=i_op, Psi=Psi, Verbose=Verbose_local, Debug=Debug_local)

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "--- Resulting statevector from the action of the HO1D Elem_op on the Psi statevector operand, computed &
                        &by Action_HO1D_operator_R1 :"
      CALL Write_Vec(Op_psi, out_unit, 1, info="Op_Psi")
      WRITE(out_unit,*) "--- End resulting statevector computed by Action_HO1D_operator_R1"
    END IF
  
  END SUBROUTINE MolecCav_Action_Cavity_R1_complex

  
  SUBROUTINE MolecCav_Get_CavMode_parameter_integer(Parameter_value, CavMode, Parameter_name)
    USE QDUtil_m
    IMPLICIT NONE 

    integer,              intent(inout) :: Parameter_value                                                                            ! the current values of the indexes for each dimension
    TYPE(Cavity_mode_t), intent(in)    :: CavMode
    character(len=*),     intent(in)    :: Parameter_name

    CALL Get(Parameter_value, CavMode%Quantum_HO1D_t, Parameter_name)

  END SUBROUTINE MolecCav_Get_CavMode_parameter_integer


  SUBROUTINE MolecCav_Get_CavMode_parameter_real(Parameter_value, CavMode, Parameter_name)
    USE QDUtil_m
    IMPLICIT NONE 

    real(kind=Rkind),     intent(inout) :: Parameter_value                                                                            ! the current values of the indexes for each dimension
    TYPE(Cavity_mode_t), intent(in)    :: CavMode
    character(len=*),     intent(in)    :: Parameter_name

    IF (TO_lowercase(TRIM(Parameter_name)) == "lambda") THEN
      Parameter_value = CavMode%lambda
    ELSE 
      CALL Get(Parameter_value, CavMode%Quantum_HO1D_t, Parameter_name)
    END IF

  END SUBROUTINE MolecCav_Get_CavMode_parameter_real


  SUBROUTINE MolecCav_Write_Cavity_mode(CavMode)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
    USE QDUtil_m
    IMPLICIT NONE 
    
    TYPE(Cavity_mode_t), intent(in) :: CavMode

    integer                          :: i_op

    WRITE(out_unit,*) "___________________________________The parameters of the Cavity mode___________"
    WRITE(out_unit,*) "|The associated 1D Quantum HO (CavMode%Quantum_HO1D_t) :                       |"
    CALL Write(CavMode%Quantum_HO1D_t)
    WRITE(out_unit,*) "_______________________________The peculiar parameters to the Cavity mode_____________________________"
    WRITE(out_unit,*) "|The strength parameter of its coupling with a cavity mode (CavMode%lambda)    | "//TO_string(CavMode%lambda)
    WRITE(out_unit,*) "|_____________________________________End Cavity mode parameters_______________|"
    FLUSH(out_unit)

  END SUBROUTINE MolecCav_Write_Cavity_mode


  SUBROUTINE MolecCav_Deallocate_Cavity_mode(CavMode, Verbose, Debug)
    USE QDUtil_m
    USE Quantum_HO1D_m
    IMPLICIT NONE 

    TYPE(Cavity_mode_t), intent(inout) :: CavMode
    integer, optional,   intent(in)    :: Verbose                                                                                 ! cf. comments in HO1D_parameters_m
    logical, optional,   intent(in)    :: Debug                                                                                   ! cf. comments in HO1D_parameters_m

    integer                            :: i_op
    integer                            :: Verbose_local                                                                      ! goes from 25 (= 0 verbose) to 29 (= maximum verbose) at this layer
    logical                            :: Debug_local

    !------------------------------------------------------Debugging options-----------------------------------------------------
    IF (PRESENT(Verbose)) THEN; Verbose_local = Verbose
    ELSE; Verbose_local = 16; END IF 
    IF (PRESENT(Debug))   THEN; Debug_local   = Debug
    ELSE; Debug_local   = .FALSE.; END IF
    IF (Debug_local) Verbose_local = 24

    IF (Verbose_local > 16) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "o o Arguments of MolecCav_Deallocate_Cavity_mode :"
      WRITE(out_unit,*) "The <<CavMode>> argument :"
      CALL Write(CavMode)
      FLUSH(out_unit)
    END IF

    !-----------------------------Deallocating the HO1D operator object----------------------------  
    CavMode%lambda        = -ONE

  END SUBROUTINE MolecCav_Deallocate_Cavity_mode


END MODULE