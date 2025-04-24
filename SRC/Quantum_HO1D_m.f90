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
! Append_quantum_HO1D     : add an HO operator to a already initialized object of Quantum_HO1D_t type.
! Write_quantum_HO1D      : display values of the type in the output
! Deallocate_quantum_HO1D : deallocate all tables of the type
! The module to initialize the HO by reading its parameters from the namelist.  
! Read_HO1D_parameters  : reads the namelist and initialize the type.
! Write_HO1D_parameters : displays values of the type in the output.
!==================================================================================================
!==================================================================================================
MODULE Quantum_HO1D_m
  !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT, real64
  USE QDUtil_m                                                                 ! gives Rkind=real64; out_unit=OUTPUT_UNIT; INPUT_UNIT=in_unit; EYE=i and other numbers; TO_LOWERCASE; TO_UPPERCASE;... We thereby use ZERO instead of 0.0_real64
  USE Elem_op_m
  IMPLICIT NONE


  TYPE :: Cavity_mode_t ! OLD                                                  ! MC = MolecCav NB: everything is initialized at values that are not supposed to make it possible of the cavity mode lecture/creation have successfully been executed
    integer          :: D      = 0                                             ! label of the HO/mode/dimension/associated basis set
    integer          :: Nb     = 0                                             ! number of basis vectors associated with the HO D
    real(kind=Rkind) :: w      = ZERO                                          ! eigenpulsation associated with the HO D
    real(kind=Rkind) :: m      = ZERO                                          ! mass associated with the HO D
    real(kind=Rkind) :: lambda = -ONE                                          ! strength parameter of the coupling between the mode D and the molecule
    real(kind=Rkind) :: eq_pos = -ONE                                          ! equilibrium position of the HO
  END TYPE

  TYPE                           :: Quantum_HO1D_t
    integer                      :: Nb      = 0
    real(kind=Rkind)             :: w       = ZERO
    real(kind=Rkind)             :: m       = ZERO
    integer                      :: Nb_op   = 0
    TYPE(Elem_op_t), allocatable :: Tab_op(:)                                  ! 0 : \hat{Id} ; 1 : \hat{H} ; 2 : \hat{x} ; 3 : \hat{N} ; 4 : \hat{\mu_{mat}} ; 5 : \hat{who knows ?}
    integer                      :: Nq      = 0
    real(kind=Rkind)             :: Eq_pos  = -ONE
    real(kind=Rkind)             :: Scale_q = HUGE(ONE)
  END TYPE


  PUBLIC Cavity_mode_t, Read_cavity_mode, Write_cavity_mode,& ! OLD 
       & Quantum_HO1D_t, Initialize, Action, Write, Dealloc

  INTERFACE Initialize
    MODULE PROCEDURE MolecCav_Initialize_quantum_HO1D, MolecCav_Initialize_QHO1D_Elem_op
  END INTERFACE
  INTERFACE Initialize_I
    MODULE PROCEDURE MolecCav_Initialize_I_QHO1D
  END INTERFACE
  INTERFACE Initialize_H
    MODULE PROCEDURE MolecCav_Initialize_H_QHO1D
  END INTERFACE
  INTERFACE Initialize_x
    MODULE PROCEDURE MolecCav_Initialize_x_QHO1D
  END INTERFACE
  INTERFACE Initialize_N
    MODULE PROCEDURE MolecCav_Initialize_N_QHO1D
  END INTERFACE
  INTERFACE Action
    MODULE PROCEDURE MolecCav_Action_quantum_HO1D_R1_real, MolecCav_Action_quantum_HO1D_R1_complex
  END INTERFACE
  INTERFACE Write
    MODULE PROCEDURE MolecCav_Write_quantum_HO1D
  END INTERFACE
  INTERFACE Dealloc
    MODULE PROCEDURE MolecCav_Deallocate_quantum_HO1D
  END INTERFACE
    

  INTERFACE Read_cavity_mode ! OLD
    MODULE PROCEDURE MolecCav_Read_cavity_mode
  END INTERFACE
  INTERFACE Write_cavity_mode ! OLD
    MODULE PROCEDURE MolecCav_Write_cavity_mode_old
  END INTERFACE


  CONTAINS


  SUBROUTINE MolecCav_Initialize_quantum_HO1D(QHO1D, Nb, w, m, Nb_op, Dense, Verbose, Debug) ! here init on the 1D QHO basis : don't need Nq etc. will write later Init_grid or Init_other_basis, maybe called by this one, in this case, the "call" will be determined by the optional arg (Nb, Nq etc.)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
    USE QDUtil_m
    USE Elem_op_m
    IMPLICIT NONE
  
    TYPE(Quantum_HO1D_t), intent(inout) :: QHO1D
    integer,              intent(in)    :: Nb
    real(kind=Rkind),     intent(in)    :: w 
    real(kind=Rkind),     intent(in)    :: m 
    integer, optional,    intent(in)    :: Nb_op
    logical, optional,    intent(in)    :: Dense                                                                         ! cf. comments in HO1D_parameters_m
    integer, optional,    intent(in)    :: Verbose                                                                         ! cf. comments in HO1D_parameters_m
    logical, optional,    intent(in)    :: Debug                                                                           ! cf. comments in HO1D_parameters_m

    logical                             :: Dense_local                                                              ! goes from 20 (= 0 verbose) to 24 (= maximum verbose) at this layer
    integer                             :: Verbose_local                                                              ! goes from 20 (= 0 verbose) to 24 (= maximum verbose) at this layer
    logical                             :: Debug_local
    
    !------------------------------------------------------Debugging options-----------------------------------------------------
    IF (PRESENT(Verbose)) THEN; Verbose_local = Verbose
    ELSE; Verbose_local = 20; END IF 
    IF (PRESENT(Debug))   THEN; Debug_local   = Debug
    ELSE; Debug_local = .FALSE.; END IF

    IF (Verbose_local > 20) WRITE(out_unit,*) 
    IF (Verbose_local > 20) WRITE(out_unit,*) "-------------------------------------------------INITIALIZING THE QUANTUM HO1D OBJ&
                                              &ECT-------------------------------------------------"; FLUSH(out_unit)

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "--- Arguments of MolecCav_Initialize_quantum_HO1D :"
      WRITE(out_unit,*) "The <<QHO1D>> argument :"
      CALL Write(QHO1D)
      WRITE(out_unit,*) "The <<Nb>> argument :"//TO_string(Nb)
      WRITE(out_unit,*) "The <<w>>  argument :"//TO_string(w)
      WRITE(out_unit,*) "The <<m>>  argument :"//TO_string(m)
      IF (PRESENT(Dense)) WRITE(out_unit,*) "The <<Dense>> argument : "//TO_string(Dense)
      WRITE(out_unit,*) "--- End arguments of MolecCav_Initialize_quantum_HO1D"
      FLUSH(out_unit)
    END IF
    
    !------------------------------------------Initializing the parameters of the 1D QHO------------------------------------------
    QHO1D%Nb = Nb
    QHO1D%w  = w
    QHO1D%m  = m
    IF (PRESENT(Nb_op)) THEN; QHO1D%Nb_op  = Nb_op
    ELSE; QHO1D%Nb_op  = 4; END IF 
    ALLOCATE(QHO1D%Tab_op(0:QHO1D%Nb_op-1))
    
    !--------------------------------------Constructing the operators to build-------------------------------------
    IF (PRESENT(Dense)) THEN; Dense_local = Dense
    ELSE; Dense_local = .FALSE.; END IF

    CALL Initialize(QHO1D%Tab_op(0), "Identity",    Nb=Nb,           Dense=Dense_local, Verbose=Verbose_local, Debug=Debug_local)
    CALL Initialize(QHO1D%Tab_op(1), "Hamiltonian", Nb=Nb, w=w,      Dense=Dense_local, Verbose=Verbose_local, Debug=Debug_local)
    CALL Initialize(QHO1D%Tab_op(2), "Position",    Nb=Nb, w=w, m=m, Dense=Dense_local, Verbose=Verbose_local, Debug=Debug_local)
    CALL Initialize(QHO1D%Tab_op(3), "NbQuanta",    Nb=Nb,           Dense=Dense_local, Verbose=Verbose_local, Debug=Debug_local)

    IF (Verbose_local > 20) WRITE(out_unit,*) 
    IF (Verbose_local > 20) WRITE(out_unit,*) "--------------------------------------------------QUANTUM HO1D OBJECT INITIALIZED-&
                                              &------------------------------------------------"; FLUSH(out_unit)

  END SUBROUTINE MolecCav_Initialize_quantum_HO1D


  SUBROUTINE MolecCav_Initialize_QHO1D_Elem_op(Elem_op, Operator_type, Nb, w, m, Dense, Verbose, Debug)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
    USE QDUtil_m
    USE Elem_op_m
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

  END SUBROUTINE MolecCav_Initialize_QHO1D_Elem_op


  SUBROUTINE MolecCav_Initialize_I_QHO1D(Identity, Nb, Verbose, Debug)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    USE Elem_op_m
    IMPLICIT NONE
    
    TYPE(Elem_op_t),   intent(inout) :: Identity                                                                        ! matrix of the one-dimensional harmonic Hamiltonian associated with HO D
    integer,           intent(in)    :: Nb                                                                            ! cf. comments in HO1D_parameters_m
    integer, optional, intent(in)    :: Verbose                                                                            ! cf. comments in HO1D_parameters_m
    logical, optional, intent(in)    :: Debug                                                                              ! cf. comments in HO1D_parameters_m

    integer                          :: i
    integer                          :: Verbose_local                                                                 ! goes from 25 (= 0 verbose) to 29 (= maximum verbose) at this layer
    logical                          :: Debug_local

    !------------------------------------------------------Debugging options-----------------------------------------------------
    IF (PRESENT(Verbose)) THEN; Verbose_local = Verbose
    ELSE; Verbose_local = 20; END IF 
    IF (PRESENT(Debug))   THEN; Debug_local   = Debug
    ELSE; Debug_local = .FALSE.; END IF

    IF (Verbose_local > 27) WRITE(out_unit,*) 
    IF (Verbose_local > 27) WRITE(out_unit,*) "----------------------------------Constructing the matrix representation of the 1D&
                                              & HO Identity---------------------------------"
 
    !---------------------------------------------Construction of the matrix Operator--------------------------------------------
    IF (.NOT. Identity%Dense) THEN
      IF (Verbose_local > 28) WRITE(out_unit,*) "--- The Dense parameter of the Identity is .FALSE., so the 1D HO Identity'&
                                               &s matrix representation will be a rank-0 tensor of ONE"
      !---------------------------------------------Initialization to default values---------------------------------------------
      ALLOCATE(Identity%Diag_val(1))

      !------------------------------------------------Construction of the matrix------------------------------------------------
      Identity%Diag_val = ONE                                                                 ! "-1" because the first Fortran vector is the fundamental eigenvector of the HO i.e. the 0^{th} ket 

      IF (Debug_local) WRITE(out_unit,*)
      IF (Debug_local) CALL Write_Vec(Identity%Diag_val, out_unit, 1, info="QHO1DIdentity")

    ELSE
      IF (Verbose_local > 28) WRITE(out_unit,*) "--- The Dense parameter of the Identity is .TRUE., so the full 1D HO Identity&
                                                &'s matrix will be constructed (in Eigenbasis) for the representation, as if t&
                                                &he analytical matrix was a dense one"
      !---------------------------------------------Initialization to default values---------------------------------------------
      ALLOCATE(Identity%Dense_val(Nb, Nb))
      Identity%Dense_val = ZERO

      !------------------------------------------------Construction of the matrix------------------------------------------------
      DO i = 1, Nb                                                                                                     ! /!\ Fortran counts from 1 to Nb !!! /!\
        Identity%Dense_val(i,i) = ONE                                                              ! "-1" because the first Fortran vector is the fundamental eigenvector of the HO i.e. the 0^{th} ket 
      END DO

      IF (Debug_local) WRITE(out_unit,*)
      IF (Debug_local) CALL Write_Mat(Identity%Dense_val, out_unit, Size(Identity%Dense_val), info="QHO1DIdentity")
    END IF
      
  END SUBROUTINE MolecCav_Initialize_I_QHO1D


  SUBROUTINE MolecCav_Initialize_H_QHO1D(Hamiltonian, Nb, w, Verbose, Debug)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    USE Elem_op_m
    IMPLICIT NONE
    
    TYPE(Elem_op_t),   intent(inout) :: Hamiltonian                                                                        ! matrix of the one-dimensional harmonic Hamiltonian associated with HO D
    integer,           intent(in)    :: Nb                                                                            ! cf. comments in HO1D_parameters_m
    real(kind=Rkind),  intent(in)    :: w                                                                            ! cf. comments in HO1D_parameters_m
    integer, optional, intent(in)    :: Verbose                                                                            ! cf. comments in HO1D_parameters_m
    logical, optional, intent(in)    :: Debug                                                                              ! cf. comments in HO1D_parameters_m

    integer                          :: i                                                                                  ! loop increments /!\ Fortran counts from 1 to Nb !!! /!\
    integer                          :: Verbose_local                                                                 ! goes from 25 (= 0 verbose) to 29 (= maximum verbose) at this layer
    logical                          :: Debug_local

    !------------------------------------------------------Debugging options-----------------------------------------------------
    IF (PRESENT(Verbose)) THEN; Verbose_local = Verbose
    ELSE; Verbose_local = 20; END IF 
    IF (PRESENT(Debug))   THEN; Debug_local   = Debug
    ELSE; Debug_local = .FALSE.; END IF

    IF (Verbose_local > 27) WRITE(out_unit,*) 
    IF (Verbose_local > 27) WRITE(out_unit,*) "----------------------------------Constructing the matrix representation of the 1D&
                                              & HO Hamiltonian---------------------------------"
 
    !---------------------------------------------Construction of the matrix Operator--------------------------------------------
    IF (.NOT. Hamiltonian%Dense) THEN
      IF (Verbose_local > 28) WRITE(out_unit,*) "--- The Dense parameter of the Hamiltonian is .FALSE., so the 1D HO Hamiltonian'&
                                               &s matrix representation will be a rank-1 tensor of the diagonal elementsof its an&
                                               &alytical matrix (in Eigenbasis)"
      !---------------------------------------------Initialization to default values---------------------------------------------
      ALLOCATE(Hamiltonian%Diag_val(Nb))

      !------------------------------------------------Construction of the matrix------------------------------------------------
      DO i = 1, Nb                                                                                                     ! /!\ Fortran counts from 1 to Nb !!! /!\
        Hamiltonian%Diag_val(i) = w*(i - ONE + HALF)                                                                 ! "-1" because the first Fortran vector is the fundamental eigenvector of the HO i.e. the 0^{th} ket 
      END DO

      IF (Debug_local) WRITE(out_unit,*)
      IF (Debug_local) CALL Write_Vec(Hamiltonian%Diag_val, out_unit, 1, info="HO1DHamiltonian")

    ELSE
      IF (Verbose_local > 28) WRITE(out_unit,*) "--- The Dense parameter of the Hamiltonian is .TRUE., so the full 1D HO Hamilton&
                                                &ian's matrix will be constructed (in Eigenbasis) for the representation, as if t&
                                                &he analytical matrix was a dense one"
      !---------------------------------------------Initialization to default values---------------------------------------------
      ALLOCATE(Hamiltonian%Dense_val(Nb, Nb))
      Hamiltonian%Dense_val = ZERO

      !------------------------------------------------Construction of the matrix------------------------------------------------
      DO i = 1, Nb                                                                                                     ! /!\ Fortran counts from 1 to Nb !!! /!\
        Hamiltonian%Dense_val(i,i) = w*(i - ONE + HALF)                                                              ! "-1" because the first Fortran vector is the fundamental eigenvector of the HO i.e. the 0^{th} ket 
      END DO

      IF (Debug_local) WRITE(out_unit,*)
      IF (Debug_local) CALL Write_Mat(Hamiltonian%Dense_val, out_unit, Size(Hamiltonian%Dense_val), info="HO1DHamiltonian")
    END IF
      
  END SUBROUTINE MolecCav_Initialize_H_QHO1D


  SUBROUTINE MolecCav_Initialize_x_QHO1D(PositionOp, Nb, w, m, Verbose, Debug)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT, OUTPUT_UNIT, real64 
    USE QDUtil_m
    USE Elem_op_m
   IMPLICIT NONE
    
    TYPE(Elem_op_t),   intent(inout) :: PositionOp
    integer,           intent(in)    :: Nb                                                                            ! cf. comments in HO1D_parameters_m
    real(kind=Rkind),  intent(in)    :: w                                                                            ! cf. comments in HO1D_parameters_m
    real(kind=Rkind),  intent(in)    :: m                                                                            ! cf. comments in HO1D_parameters_m
    integer, optional, intent(in)    :: Verbose                                                                            ! cf. comments in HO1D_parameters_m
    logical, optional, intent(in)    :: Debug                                                                              ! cf. comments in HO1D_parameters_m

    integer                          :: i                                                                                  ! loop increments /!\ Fortran counts from 1 to Nb !!! /!\
    integer                          :: Verbose_local                                                                 ! goes from 25 (= 0 verbose) to 29 (= maximum verbose) at this layer
    logical                          :: Debug_local

    !------------------------------------------------------Debugging options-----------------------------------------------------
    IF (PRESENT(Verbose)) THEN; Verbose_local = Verbose
    ELSE; Verbose_local = 20; END IF 
    IF (PRESENT(Debug))   THEN; Debug_local   = Debug
    ELSE; Debug_local = .FALSE.; END IF

    IF (Verbose_local > 27) WRITE(out_unit,*) 
    IF (Verbose_local > 27) WRITE(out_unit,*) "-------------------------------Constructing the matrix representation of the 1D HO&
                                             & Position operator------------------------------"

    !---------------------------------------------Construction of the matrix Operator--------------------------------------------
    IF ((.NOT. PositionOp%Dense) .AND. Nb > 1) THEN
      IF (Verbose_local > 28) WRITE(out_unit,*) "--- The Dense parameter of the Position operator is .FALSE., so the 1D HO Positi&
                                                &on operator's matrix representation will be a rank-2 tensor of the tridiagonal e&
                                                &lements of its analytical matrix (in Eigenbasis)"
      !-----------------------------------Initialization of the characteristics of the operator----------------------------------
      PositionOp%Upper_bandwidth   = 1
      PositionOp%Lower_bandwidth   = 1

      !---------------------------------------------Initialization to default values---------------------------------------------
      ALLOCATE(PositionOp%Band_val(Nb,3))                                                                            ! Nb lines (number of diagonal elements) and 3 columns because 3 bands to consider : the diagonal, and the two bands above and below it
      PositionOp%Band_val = ZERO

      !------------------------------------------------Construction of the matrix------------------------------------------------
      DO i = 1, Nb - 1                                                                                                 ! /!\ Fortran counts from 1 to Nb !!! /!\ Nb-1 not to have Band_val(i+1) out of range
        PositionOp%Band_val(i,1)   = SQRT(REAL(i,kind=Rkind))
        PositionOp%Band_val(i+1,3) = SQRT(REAL(i,kind=Rkind))
      END DO
      PositionOp%Band_val = PositionOp%Band_val / SQRT(TWO * w * m)
    
      IF (Debug_local) WRITE(out_unit,*)
      IF (Debug_local) CALL Write_Mat(PositionOp%Band_val, out_unit, 3, info="HO1DPositionOp")

    ELSE IF (.NOT. PositionOp%Dense) THEN
      IF (Verbose_local > 28) WRITE(out_unit,*) "--- The Dense parameter of the Position operator is .FALSE. BUT the basis set si&
                                                &ze is only 1, so the 1D HO Position operator's matrix representation will use th&
                                                &e Diag_val rank-1 tensor to store the only element of the analytical matrix (i&
                                                &n Eigenbasis)"
      !---------------------------------------------Initialization to default values---------------------------------------------
      ALLOCATE(PositionOp%Diag_val(Nb))

      !------------------------------------------------Construction of the matrix------------------------------------------------
      DO i = 1, Nb                                                                                                     ! /!\ Fortran counts from 1 to Nb !!! /!\
        PositionOp%Diag_val(i) = ZERO                                                                                          ! the position operator matrix has first value (i.e. only value in the Nb = 0 case) 0 
      END DO

      IF (Debug_local) WRITE(out_unit,*)
      IF (Debug_local) CALL Write_Vec(PositionOp%Diag_val, out_unit, 1, info="HO1DPositionOp")

    ELSE 
      IF (Verbose_local > 28) WRITE(out_unit,*) "--- The Dense parameter of the Position operator is .TRUE., so the full 1D HO Po&
                                                &sition operator's matrix will be constructed (in Eigenbasis) for the representat&
                                                &ion, as if the analytical matrix was a dense one"
      !---------------------------------------------Initialization to default values---------------------------------------------
      ALLOCATE(PositionOp%Dense_val(Nb, Nb))
      PositionOp%Dense_val = ZERO
      
      !------------------------------------------------Construction of the matrix------------------------------------------------
      DO i = 1, Nb - 1                                                                                                 ! /!\ Fortran counts from 1 to Nb !!! /!\
        PositionOp%Dense_val(i,i+1) = SQRT(REAL(i,kind=Rkind))
        PositionOp%Dense_val(i+1,i) = SQRT(REAL(i,kind=Rkind))
      END DO
      PositionOp%Dense_val = PositionOp%Dense_val / SQRT(TWO * w * m)
    
      IF (Debug_local) WRITE(out_unit,*)
      IF (Debug_local) CALL Write_Mat(PositionOp%Dense_val, out_unit, Size(PositionOp%Dense_val), info="HO1DPositionOp")
    END IF
      
  END SUBROUTINE MolecCav_Initialize_x_QHO1D


  SUBROUTINE MolecCav_Initialize_N_QHO1D(NbQuanta, Nb, Verbose, Debug)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    USE Elem_op_m
    IMPLICIT NONE
    
    TYPE(Elem_op_t),   intent(inout) :: NbQuanta
    integer,           intent(in)    :: Nb                                                                            ! cf. comments in HO1D_parameters_m
    integer, optional, intent(in)    :: Verbose                                                                            ! cf. comments in HO1D_parameters_m
    logical, optional, intent(in)    :: Debug                                                                              ! cf. comments in HO1D_parameters_m

    integer                          :: i                                                                                  ! loop increments /!\ Fortran counts from 1 to Nb !!! /!\
    integer                          :: Verbose_local                                                                 ! goes from 25 (= 0 verbose) to 29 (= maximum verbose) at this layer
    logical                          :: Debug_local

    !------------------------------------------------------Debugging options-----------------------------------------------------
    IF (PRESENT(Verbose)) THEN; Verbose_local = Verbose
    ELSE; Verbose_local = 20; END IF 
    IF (PRESENT(Debug))   THEN; Debug_local   = Debug
    ELSE; Debug_local = .FALSE.; END IF

    IF (Verbose_local > 27) WRITE(out_unit,*) 
    IF (Verbose_local > 27) WRITE(out_unit,*) "---------------------Constructing the matrix representation of the 1D HO Number of&
                                             & excitation Quanta operator---------------------"
  
    !---------------------------------------------Construction of the matrix Operator--------------------------------------------
    IF (.NOT. NbQuanta%Dense) THEN
      IF (Verbose_local > 28) WRITE(out_unit,*) "--- The Dense parameter of the NbQuanta operator is .FALSE., so the 1D HO NbQuan&
                                                &ta's matrix representation will be a rank-1 tensor of the diagonal elements of i&
                                                &ts analytical matrix (in Eigenbasis)"
      !---------------------------------------------Initialization to default values---------------------------------------------
      ALLOCATE(NbQuanta%Diag_val(Nb))

      !------------------------------------------------Construction of the matrix------------------------------------------------
      DO i = 1, Nb                                                                                                     ! /!\ Fortran counts from 1 to Nb !!! /!\
        NbQuanta%Diag_val(i) = i - 1
      END DO
  
      IF (Debug_local) WRITE(out_unit,*)
      IF (Debug_local) CALL Write_Vec(NbQuanta%Diag_val, out_unit, 1, info="HO1DNbQuanta")

    ELSE
      IF (Verbose_local > 28) WRITE(out_unit,*) "--- The Dense parameter of the NbQuanta is .TRUE., so the full 1D HO NbQuanta's &
                                                &matrix will be constructed (in Eigenbasis) for the representation, as if the ana&
                                                &lytical matrix was a dense one"
      !---------------------------------------------Initialization to default values---------------------------------------------
      ALLOCATE(NbQuanta%Dense_val(Nb, Nb))
      NbQuanta%Dense_val = ZERO

      !------------------------------------------------Construction of the matrix------------------------------------------------
      DO i = 1, Nb                                                                            ! /!\ Fortran counts from 1 to Nb !!! /!\
        NbQuanta%Dense_val(i,i) = i - 1
      END DO

      IF (Debug_local) WRITE(out_unit,*)
      IF (Debug_local) CALL Write_Mat(NbQuanta%Dense_val, out_unit, Size(NbQuanta%Dense_val), info="HO1DNbQuanta")
    END IF
      
  END SUBROUTINE MolecCav_Initialize_N_QHO1D

  
  SUBROUTINE MolecCav_Action_quantum_HO1D_R1_real(Op_psi, QHO1D, i_op, Psi, Verbose, Debug)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    USE Elem_op_m
    IMPLICIT NONE

    real(kind=Rkind),     intent(inout) :: Op_psi(:)
    TYPE(Quantum_HO1D_t), intent(in)    :: QHO1D
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
  
  END SUBROUTINE MolecCav_Action_quantum_HO1D_R1_real

  
  SUBROUTINE MolecCav_Action_quantum_HO1D_R1_complex(Op_psi, QHO1D, i_op, Psi, Verbose, Debug)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    USE Elem_op_m
    IMPLICIT NONE

    complex(kind=Rkind),  intent(inout) :: Op_psi(:)
    TYPE(Quantum_HO1D_t), intent(in)    :: QHO1D
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
  
  END SUBROUTINE MolecCav_Action_quantum_HO1D_R1_complex

  
  SUBROUTINE MolecCav_Write_quantum_HO1D(QHO1D)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
    USE QDUtil_m
    USE Elem_op_m
    IMPLICIT NONE 
    
    TYPE(Quantum_HO1D_t), intent(in) :: QHO1D

    integer                          :: i_op

    WRITE(out_unit,*) "____________________________________The parameters of the 1D QHO____________________________________"
    WRITE(out_unit,*) "|Basis set size of the HO (QHO1D%Nb)                                           | "//TO_string(QHO1D%Nb)
    WRITE(out_unit,*) "|______________________________________________________________________________|______________________"
    WRITE(out_unit,*) "|Eigenpulsation of the HO (QHO1D%w)                                            | "//TO_string(QHO1D%w)
    WRITE(out_unit,*) "|______________________________________________________________________________|______________________"
    WRITE(out_unit,*) "|Mass associated with the HO (QHO1D%m)                                         | "//TO_string(QHO1D%m)
    IF (ALLOCATED(QHO1D%Tab_op)) THEN 
      WRITE(out_unit,*) "|______________________________________________________________________________|______________________"
      WRITE(out_unit,*) "|QHO's Tab_op holds the following nb of operators (Size(QHO1D%Tab_op))         | ", Size(QHO1D%Tab_op)
      WRITE(out_unit,*) "|______________________________________________________________________________|______________________"
      WRITE(out_unit,*) "|The "//TO_string(i_op)//"^{th} operator associated with the HO (QHO1D%Tab_op("//TO_string(i_op)//")) : &
                        &               |"
      WRITE(out_unit,*) "|______________________________________________________________________________|"
      DO i_op = 0, Size(QHO1D%Tab_op)-1
        CALL Write(QHO1D%Tab_op(i_op))
      END DO 
    ELSE 
      WRITE(out_unit,*) "|______________________________________________________________________________|______________________"
      WRITE(out_unit,*) "|QHO's Tab_op is not allocated (QHO1D%Tab_op)                                  | /"
      WRITE(out_unit,*) "|______________________________________________________________________________|______________________"
    END IF
    WRITE(out_unit,*) "|Number of grid points of the QHO DOF (QHO1D%Nq)                               | "//TO_string(QHO1D%Nq)
    WRITE(out_unit,*) "|______________________________________________________________________________|______________________"
    WRITE(out_unit,*) "|Equilibrium position of the HO (QHO1D%Eq_pos)                                 | "//TO_string(QHO1D%Eq_pos)
    WRITE(out_unit,*) "|______________________________________________________________________________|______________________"
    WRITE(out_unit,*) "|Change in variable coefficient for the DOF grid (QHO1D%Scale_q)               | "//TO_string(QHO1D%Scale_q)
    WRITE(out_unit,*) "|________________________________________End HO1D parameters___________________|______________________"
    FLUSH(out_unit)

  END SUBROUTINE MolecCav_Write_quantum_HO1D


  SUBROUTINE MolecCav_Deallocate_quantum_HO1D(QHO1D, Verbose, Debug)
    USE QDUtil_m
    USE Elem_op_m
    IMPLICIT NONE 

    TYPE(Quantum_HO1D_t), intent(inout) :: QHO1D
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
      CALL Dealloc(QHO1D%Tab_op(i_op), Verbose=Verbose_local, Debug=Debug_local)
    END DO
    IF (ALLOCATED(QHO1D%Tab_op)) DEALLOCATE(QHO1D%Tab_op)
    QHO1D%Nq      = 0
    QHO1D%Eq_pos  = -ONE
    QHO1D%Scale_q = HUGE(ONE)

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "--- The QHO1D object after having been deallocated :"
      CALL Write(QHO1D)
      WRITE(out_unit,*) "--- End dellocating QHO1D"
    END IF

  END SUBROUTINE MolecCav_Deallocate_quantum_HO1D


  !############################################################################################
  !############################################## OLD #########################################
  !############################################################################################
  SUBROUTINE MolecCav_Read_cavity_mode(Mode, nio)                              ! nio is the label of the file from which the values have to be drawn.
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    IMPLICIT NONE
    
    TYPE(Cavity_mode_t),    intent(inout) :: Mode   
    integer,                intent(in)    :: nio

    integer                               :: D, Nb, err_io                     ! label of the basis/HO/mode/dimension, its number of basis vectors, and an error control variable
    real(kind=Rkind)                      :: w, m, lambda, eq_pos              ! eigenpulsation, mass, molecule-coupling strength, and equilibrium position associated with this HO
    logical, parameter                    :: Debug = .FALSE.

    NAMELIST /HO_1/ D, Nb, w, m, lambda, eq_pos                                ! declare the nml HO_1 and specify the parameter's list to be found within

    !----------------------Initialization to default values--------------------
    D      = 0
    Nb     = 1
    w      = ZERO
    m      = ZERO
    lambda = -ONE
    eq_pos = -ONE
 
    !------------------------------Reading of the nml--------------------------
    WRITE(out_unit,*) 
    WRITE(out_unit,*) '********************************************************************************'
    WRITE(out_unit,*) '**************************** READING BASIS OF THE HO ***************************'
    WRITE(out_unit,*) '********************************************************************************'
    
    READ(nio, nml = HO_1, iostat = err_io)                                     ! assign the values read in the nml to the declared list of parameters

    IF (Debug) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "-----------------------The namelist parameters are read as----------------------"
      WRITE(out_unit, nml = HO_1)
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
    
    !---------------Construction of the Cavity_mode_t type object-----------
    Mode%D      = D
    Mode%Nb     = Nb
    Mode%w      = w
    Mode%m      = m
    Mode%lambda = lambda
    Mode%eq_pos = eq_pos

    WRITE(out_unit,*) 
    WRITE(out_unit,*) '********************************************************************************'
    WRITE(out_unit,*) '************************** BASIS OF THE HO CONSTRUCTED *************************'
    WRITE(out_unit,*) '********************************************************************************'

    IF (Debug) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "--------------Cavity mode constructed by MolecCav_Read_cavity_mode--------------"
      CALL Write_cavity_mode(Mode)
      WRITE(out_unit,*) "------------End Cavity mode constructed by MolecCav_Read_cavity_mode------------"
    END IF

  END SUBROUTINE MolecCav_Read_cavity_mode

  
  SUBROUTINE MolecCav_Write_cavity_mode_old(Mode)
    TYPE(Cavity_mode_t), intent(in) :: Mode

    WRITE(out_unit,*) "____________________________________The associated HO cavity mode___________________________________"
    WRITE(out_unit,*) "|Index of the cavity mode Mode%D                                             | ", Mode%D
    WRITE(out_unit,*) "|____________________________________________________________________________|______________________"
    WRITE(out_unit,*) "|Basis set size of the cavity mode Mode%Nb                                   | ", Mode%Nb
    WRITE(out_unit,*) "|____________________________________________________________________________|______________________"
    WRITE(out_unit,*) "|Eigenpulsation of the cavity mode Mode%w                                    | ", Mode%w
    WRITE(out_unit,*) "|____________________________________________________________________________|______________________"
    WRITE(out_unit,*) "|Mass of the associated HO Mode%m                                            | ", Mode%m
    WRITE(out_unit,*) "|____________________________________________________________________________|______________________"
    WRITE(out_unit,*) "|Coupling strength between the matter and this cavity mode Mode%lambda       | ", Mode%lambda
    WRITE(out_unit,*) "|____________________________________________________________________________|______________________"
    WRITE(out_unit,*) "|Equilibrium position of this cavity mode's HO Mode%eq_pos                   | ", Mode%eq_pos
    WRITE(out_unit,*) "|____________________________________________________________________________|______________________"
    FLUSH(out_unit)

  END SUBROUTINE MolecCav_Write_cavity_mode_old


END MODULE