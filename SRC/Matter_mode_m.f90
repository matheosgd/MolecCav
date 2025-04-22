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
  USE Elem_op_m
  USE Quantum_HO1D_m
  IMPLICIT NONE


  TYPE, EXTENDS(Quantum_HO1D_t)  :: Matter_mode_t
    real(kind=Rkind)             :: lambda       = -ONE
    real(kind=Rkind)             :: CoeffDipMomt = -ONE
  END TYPE


  PUBLIC Matter_mode_t, Initialize, Action, Write, Dealloc

  INTERFACE Initialize
    MODULE PROCEDURE MolecCav_Initialize_matter_mode, MolecCav_Initialize_matter_mode_QHO1D
  END INTERFACE
  INTERFACE Action
    MODULE PROCEDURE MolecCav_Action_matter_mode_R1_real, MolecCav_Action_matter_R1_complex
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
    USE Elem_op_m
    IMPLICIT NONE
  
    TYPE(Matter_mode_t), intent(inout)  :: MatMode
    integer,              intent(in)    :: nio
    logical, optional,    intent(in)    :: Dense                                                                         ! cf. comments in HO1D_parameters_m
    integer, optional,    intent(in)    :: Verbose                                                                         ! cf. comments in HO1D_parameters_m
    logical, optional,    intent(in)    :: Debug                                                                           ! cf. comments in HO1D_parameters_m

    real(kind=Rkind)                    :: lambda, CoeffDipMomt, w, m, Eq_pos, Scale_q
    integer                             :: Nb, Nq
    logical                             :: Dense_local                                                              ! goes from 20 (= 0 verbose) to 24 (= maximum verbose) at this layer
    integer                             :: Verbose_local                                                              ! goes from 20 (= 0 verbose) to 24 (= maximum verbose) at this layer
    logical                             :: Debug_local

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
    
    !------------------------------------------Initializing the parameters of the 1D QHO------------------------------------------
    CALL MolecCav_Initialize_matter_mode_QHO1D(MatMode%Quantum_HO1D_t, lambda, CoeffDipMomt, nio, Dense_local, Verbose=Verbose_lo&
    &cal, Debug=Debug_local)
    
    IF (Verbose_local > 20) WRITE(out_unit,*) 
    IF (Verbose_local > 20) WRITE(out_unit,*) "--------------------------------------------------QUANTUM HO1D OBJECT INITIALIZED-&
                                              &------------------------------------------------"; FLUSH(out_unit)

  END SUBROUTINE MolecCav_Initialize_matter_mode


  SUBROUTINE MolecCav_Initialize_matter_mode_QHO1D(QHO1D, lambda, CoeffDipMomt, nio, Dense, Verbose, Debug)                              ! nio is the label of the file from which the values have to be drawn.
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE Quantum_HO1D_m
    IMPLICIT NONE
    
    TYPE(Quantum_HO1D_t), intent(inout) :: QHO1D   
    real(kind=Rkind),     intent(inout) :: lambda   
    real(kind=Rkind),     intent(inout) :: CoeffDipMomt   
    integer,              intent(in)    :: nio
    logical, optional,    intent(in)    :: Dense                                                                         ! cf. comments in HO1D_parameters_m
    integer, optional,    intent(in)    :: Verbose                                                                         ! cf. comments in HO1D_parameters_m
    logical, optional,    intent(in)    :: Debug                                                                           ! cf. comments in HO1D_parameters_m

    real(kind=Rkind)                    :: w, m, Eq_pos, Scale_q
    integer                             :: Nb, Nq
    integer                             :: err_io
    logical                             :: Dense_local                                                              ! goes from 20 (= 0 verbose) to 24 (= maximum verbose) at this layer
    integer                             :: Verbose_local                                                              ! goes from 20 (= 0 verbose) to 24 (= maximum verbose) at this layer
    logical                             :: Debug_local

    NAMELIST /Matter_mode/ lambda, CoeffDipMomt, Nb, w, m, Nq, Eq_pos, Scale_q

    !------------------------------------------------------Debugging options-----------------------------------------------------
    IF (PRESENT(Verbose)) THEN; Verbose_local = Verbose
    ELSE; Verbose_local = 20; END IF 
    IF (PRESENT(Debug))   THEN; Debug_local   = Debug
    ELSE; Debug_local = .FALSE.; END IF

    IF (Verbose_local > 20) WRITE(out_unit,*) 
    IF (Verbose_local > 20) WRITE(out_unit,*) "-------------------------------------------------INITIALIZING THE HO1D OF THE MATT&
                                              &ER MODE OBJECT-------------------------------------------------"; FLUSH(out_unit)

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "--- Arguments of MolecCav_Initialize_matter_mode_QHO1D :"
      WRITE(out_unit,*) "The <<QHO1D>> argument :"
      CALL Write(QHO1D)
      WRITE(out_unit,*) "The <<nio>> argument :"//TO_string(nio)
      IF (PRESENT(Dense))   WRITE(out_unit,*) "The <<Dense>> argument : "//TO_string(Dense)
      IF (PRESENT(Verbose)) WRITE(out_unit,*) "The <<Verbose>> argument : "//TO_string(Verbose)
      IF (PRESENT(Debug))   WRITE(out_unit,*) "The <<Debug>> argument : "//TO_string(Debug)
      WRITE(out_unit,*) "--- End arguments of MolecCav_Initialize_matter_mode_QHO1D"
      FLUSH(out_unit)
    END IF
    
    !----------------------Initialization to default values--------------------
    lambda       = -ONE
    CoeffDipMomt = -ONE
    Nb           = 0
    w            = ZERO
    m            = ZERO
    Nq           = 0
    Eq_pos       = -ONE
    Scale_q      = ZERO

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
      WRITE(out_unit,*) '######## Error in Read_cavity_mode (err_io<0) #########'
      WRITE(out_unit,*) '#######################################################'
      WRITE(out_unit,*) '################ err_io = ', err_io, '################'
      STOP '################# Check basis data ################'

    ELSE IF( err_io > 0) THEN
      WRITE(out_unit,*) ''
      WRITE(out_unit,*) '#######################################################'
      WRITE(out_unit,*) '######## Error in Read_cavity_mode (err_io>0) #########'
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
    QHO1D%Nb      = Nb
    QHO1D%w       = w
    QHO1D%m       = m
    QHO1D%Nq      = Nq
    QHO1D%Eq_pos  = Eq_pos
    QHO1D%Scale_q = Scale_q

    WRITE(out_unit,*) 
    WRITE(out_unit,*) '********************************************************************************'
    WRITE(out_unit,*) '************************** QHO1D CONSTRUCTED *************************'
    WRITE(out_unit,*) '********************************************************************************'

    IF (Debug) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "--------------Cavity mode constructed by MolecCav_Read_cavity_mode--------------"
      CALL Write(QHO1D)
      WRITE(out_unit,*) "------------End Cavity mode constructed by MolecCav_Read_cavity_mode------------"
    END IF

  END SUBROUTINE MolecCav_Initialize_matter_mode_QHO1D


  SUBROUTINE MolecCav_Action_matter_mode_R1_real(Op_psi, QHO1D, i_op, Psi, Verbose, Debug)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    USE Elem_op_m
    IMPLICIT NONE

    real(kind=Rkind),     intent(inout) :: Op_psi(:)
    TYPE(Matter_mode_t), intent(in)    :: QHO1D
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
      WRITE(out_unit,*) "--- Arguments of MolecCav_Action_quantum_HO1D :"
      WRITE(out_unit,*) "The <<QHO1D>> argument :"
      CALL Write(QHO1D)
      WRITE(out_unit,*) "The <<i_op>> argument :"//TO_string(i_op)
      WRITE(out_unit,*) "The <<Psi>> argument : "
      CALL Write_Vec(Psi, out_unit, 1, info="Psi")
      WRITE(out_unit,*) "The size of its vector : "//TO_string(Size(Psi))
      WRITE(out_unit,*) "--- End arguments of MolecCav_Action_quantum_HO1D"
      FLUSH(out_unit)
    END IF
    
    !-----------------------------------------------------Checking dimensions----------------------------------------------------
    ! ALREADY CHECKED IN THE ACTIONS CODED IN ELEM_OP_M ! (fortunately btw, otherwise the \hat{I}d case should have test above)

    !---------------------------------------------Selection of the calculation method--------------------------------------------
    IF (QHO1D%Tab_op(i_op)%Operator_type == "identity") THEN ! N.B. we should have tested i_op == 0, easier and consistent with the algorithmic choices made so far
      Op_psi = Psi
    ELSE
      CALL Action(Op_psi=Op_psi, Elem_op=QHO1D%Tab_op(i_op), Psi=Psi, Verbose=Verbose_local, Debug=Debug_local)
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
  
  END SUBROUTINE MolecCav_Action_matter_mode_R1_real


  SUBROUTINE MolecCav_Action_matter_R1_complex(Op_psi, QHO1D, i_op, Psi, Verbose, Debug)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    USE Elem_op_m
    IMPLICIT NONE

    complex(kind=Rkind),  intent(inout) :: Op_psi(:)
    TYPE(Matter_mode_t), intent(in)    :: QHO1D
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
      WRITE(out_unit,*) "--- Arguments of MolecCav_Action_quantum_HO1D :"
      WRITE(out_unit,*) "The <<QHO1D>> argument :"
      CALL Write(QHO1D)
      WRITE(out_unit,*) "The <<i_op>> argument :"//TO_string(i_op)
      WRITE(out_unit,*) "The <<Psi>> argument : "
      CALL Write_Vec(Psi, out_unit, 1, info="Psi")
      WRITE(out_unit,*) "The size of its vector : "//TO_string(Size(Psi))
      WRITE(out_unit,*) "--- End arguments of MolecCav_Action_quantum_HO1D"
      FLUSH(out_unit)
    END IF
    
    !-----------------------------------------------------Checking dimensions----------------------------------------------------
    ! ALREADY CHECKED IN THE ACTIONS CODED IN ELEM_OP_M !

    !---------------------------------------------Selection of the calculation method--------------------------------------------
    IF (QHO1D%Tab_op(i_op)%Operator_type == "identity") THEN 
      Op_psi = Psi
    ELSE
      CALL Action(Op_psi=Op_psi, Elem_op=QHO1D%Tab_op(i_op), Psi=Psi, Verbose=Verbose_local, Debug=Debug_local)
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
  
  END SUBROUTINE MolecCav_Action_matter_R1_complex

  
  SUBROUTINE MolecCav_Write_matter_mode(QHO1D)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
    USE QDUtil_m
    USE Elem_op_m
    IMPLICIT NONE 
    
    TYPE(Matter_mode_t), intent(in) :: QHO1D

    integer                          :: i_op

    WRITE(out_unit,*) "____________________________________The parameters of the 1D QHO____________________________________"
    WRITE(out_unit,*) "|Basis set size of the HO (QHO1D%Nb)                                           | "//TO_string(QHO1D%Nb)
    WRITE(out_unit,*) "|______________________________________________________________________________|______________________"
    WRITE(out_unit,*) "|Eigenpulsation of the HO (QHO1D%w)                                            | "//TO_string(QHO1D%w)
    WRITE(out_unit,*) "|______________________________________________________________________________|______________________"
    WRITE(out_unit,*) "|Mass associated with the HO (QHO1D%m)                                         | "//TO_string(QHO1D%m)
    IF (ALLOCATED(QHO1D%Tab_op)) THEN 
      WRITE(out_unit,*) "|______________________________________________________________________________|______________________"
      WRITE(out_unit,*) "|The "//TO_string(i_op)//"^{th} operator associated with the HO (QHO1D%Tab_op("//TO_string(i_op)//")) : &
                        &               |"
      WRITE(out_unit,*) "|______________________________________________________________________________|"
      DO i_op = 0, Size(QHO1D%Tab_op)-1
        CALL Write(QHO1D%Tab_op(i_op))
      END DO 
    END IF
    WRITE(out_unit,*) "|______________________________________________________________________________|______________________"
    WRITE(out_unit,*) "|Number of grid points of the QHO DOF (QHO1D%Nq)                               | "//TO_string(QHO1D%Nq)
    WRITE(out_unit,*) "|______________________________________________________________________________|______________________"
    WRITE(out_unit,*) "|Equilibrium position of the HO (QHO1D%Eq_pos)                                 | "//TO_string(QHO1D%Eq_pos)
    WRITE(out_unit,*) "|______________________________________________________________________________|______________________"
    WRITE(out_unit,*) "|Change in variable coefficient for the DOF grid (QHO1D%Scale_q)               | "//TO_string(QHO1D%Scale_q)
    WRITE(out_unit,*) "|________________________________________End HO1D parameters___________________|______________________"
    FLUSH(out_unit)

  END SUBROUTINE MolecCav_Write_matter_mode


  SUBROUTINE MolecCav_Deallocate_matter_mode(QHO1D, Verbose, Debug)
    USE QDUtil_m
    USE Elem_op_m
    IMPLICIT NONE 

    TYPE(Matter_mode_t), intent(inout) :: QHO1D
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
      WRITE(out_unit,*) "--- The QHO1D to be deallocated :"
      CALL Write(QHO1D)
      WRITE(out_unit,*) "--- End QHO1D to be deallocated"
    END IF 

    !-----------------------------Deallocating the HO1D operator object----------------------------
    IF (Verbose_local > 27) WRITE(out_unit,*)
    IF (Verbose_local > 27) WRITE(out_unit,*) "-----------------------------------------------Deallocating the QHO1D obje&
                                              &ct----------------------------------------------"
  
    QHO1D%Nb = 0
    QHO1D%w  = ZERO
    QHO1D%m  = ZERO
    DO i_op = 0, SIZE(QHO1D%Tab_op)-1
      IF (QHO1D%Tab_op(i_op)%Dense         .AND. ALLOCATED(QHO1D%Tab_op(i_op)%Dense_val)) DEALLOCATE(QHO1D%Tab_op(i_op)%Dense_val)
      IF ((.NOT. QHO1D%Tab_op(i_op)%Dense) .AND. ALLOCATED(QHO1D%Tab_op(i_op)%Diag_val))  DEALLOCATE(QHO1D%Tab_op(i_op)%Diag_val)
      IF ((.NOT. QHO1D%Tab_op(i_op)%Dense) .AND. ALLOCATED(QHO1D%Tab_op(i_op)%Band_val))  DEALLOCATE(QHO1D%Tab_op(i_op)%Band_val)
    END DO
    IF (ALLOCATED(QHO1D%Tab_op)) DEALLOCATE(QHO1D%Tab_op)
    QHO1D%Nq      = 0
    QHO1D%Eq_pos  = -ONE
    QHO1D%Scale_q = ZERO

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "--- The QHO1D object after having been deallocated :"
      CALL Write(QHO1D)
      WRITE(out_unit,*) "--- End dellocating QHO1D"
    END IF

  END SUBROUTINE MolecCav_Deallocate_matter_mode


END MODULE