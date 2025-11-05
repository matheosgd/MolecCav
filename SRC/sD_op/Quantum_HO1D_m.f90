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
! Module designed to manage GENERAL 1D Quantum Harmonic Oscillators (QHO1D), independantly on any 
! other parts of the code, except for Elem_op_m which is used to manage each operator associated to
! the de scribed QHO1D and their action on a wavefunction representing a state on this QHO1D. ONE 
! QHO1D is represented by an object of derived type Quantum_HO1D_t In other words, this module may 
! be seen as the implementation of one of the possible basis to construct the Elem_op_t objects on.
! If one want to implement another representation (basis, grid) for the elementary operators, it s-
! hould be done in a module built in an analogous way. It is aimed to use as general methods as po-
! ssible in order to be useful not only in the MolecCav program but any other one that would need 
! to describe a 1D HO within a few modifications only. It contains a set of initialisation subrout-
! ines for the QHO1D and the associated operators - but one just has to call Initialize(QHO1D, Nb, 
! w, m) once and all others will be automatically chained - subroutines to compute the Action of a
! QHO1D operator, to Write and Deallocate a QHO1D object, and to Get some of its parameter to pass
! it to less internal module. 
! This module is in the "-4" level. At this level, Verbose is ranging from 16 (degree 0) to 20 (de-
! gree 4).
!
! Quantum_HO1D_t : The derived type used to represent a Quantum Harmonic Oscillator of 1 Dimension 
! (QHO1D). More rigorously, it represents the basis set associated with a QHO1D, through the physi-
! cal parameters of the potential, the list gathering all the quantum mechanics operators associat-
! ed to this HO expanded as rank-2 tensors unpon the HO's eigenvectors basis set, and some spatial
! paramaters that are unused for now but left available for later grid implementation.
!   - Nb : The number (integer) of vectors that one want to build the HO's basis set with. It will
! thus decide the sizes of the Elem_op_t operator's matrices associated with this QHO1D, as well as
! the sizes associated to the wavefunction matrix representation on which the operators are expect-
! ed to take action. 
!   - w : The real number that represents the eigenpulsation associated to the QHO1D potential.
!   - m : The real number that represents the mass associated to the QHO1D potential.
!   - Nb_op : The number (integer) of operators (Elem_op_t type objects) that are expected to be a-
! ssociated to this QHO1D and thus constructed for it in this module. In other words it is the len-
! gth the parameter Tab_op is expected to be allocated to. 
!   - Tab_op : The table in which the QHO1D's operators are actually constructed and stored. It is 
! a rank-1 table of length Nb_op and of elements' type Elem_op. Here is an crucial key of the code 
! framework : in aim to be able to access easily to the operator an order had been defined in whic-
! h they are ALWAYS put in Tab_op, such as the same index will ALWAYS correspond to the same Opera-
! tors's type EVRYWHERE in the code. We give here the correspondance of the nomenclature : "0" : \-
! hat{Identity} ; "1" : \hat{Hamiltonian} ; "2" : \hat{Position} ; "3" : \hat{Number of quanta exc-
! itation} ; 4 : \hat{Dipole moment} (only if the QHO1D is associated with a matter mode of the sy-
! stem) ; 5 : \hat{To be implemented ?}.
!   - Nq : The number (integer) of grid points associated to the discretisation of the spatial coo-
! rdinate of this QHO1D (<=> \hat{Position}).
!   - Eq_pos : The real number corresponding to the equilibrium spatial position of this QHO1D's p-
! otential.
!   - Scale_q : A real number used for the change of coordinate of the space grid representation (-
! the scale factor).
! 
! Initialize : Interface for the MolecCav_Initialize_quantum_HO1D and the MolecCav_Initialize_QHO1-
! D_Elem_op and procedures. Takes as arguments
!
! Action : Interface for the MolecCav_Action_quantum_HO1D_R1_* procedures. Takes as arguments
!
! Get : Interface for the MolecCav_Get_QHO1D_parameter_* procedures. Takes as arguments
!
! Write : Interface for the MolecCav_Write_quantum_HO1D procedure. cf. MolecCav_Write_quantum_HO1D 
! for details.
!
! Dealloc : Interface for the MolecCav_Deallocate_quantum_HO1D procedure. cf. MolecCav_Deallocate_-
! quantum_HO1D for details.
!
! MolecCav_Initialize_quantum_HO1D : 
!
! MolecCav_Initialize_QHO1D_Elem_op :
!
! Initialize_I :
!
! Initialize_H :
!
! Initialize_x :
!
! Initialize_N :
!
! MolecCav_Action_quantum_HO1D_R1_real :
!
! MolecCav_Action_quantum_HO1D_R1_complex :
!
! MolecCav_Get_QHO1D_parameter_integer :
!
! MolecCav_Get_QHO1D_parameter_real :
!
! MolecCav_Write_quantum_HO1D : Takes as argument an object of derived type Quantum_HO1D_t (QHO1D),
! an optional string messsage (Info), and an optional logical (More), and write in the standard ou-
! tput all of the QHO1D parameters, below the Info if provided. Yet, if More is not provided (or set to .F.), the elements inside the Tab_op parameters will not be developed because it would be most of the time redundant. 
!
! MolecCav_Deallocate_quantum_HO1D : Takes as argument an object of derived type Quantum_HO1D_t (Q-
! HO1D) and reset it by deallocating all tables/strings and setting to defaut values the other par-
! ameters.
!
!==================================================================================================
!==================================================================================================
MODULE Quantum_HO1D_m
  !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT, real64
  USE QDUtil_m                                                                                          ! gives Rkind=real64; out_unit=OUTPUT_UNIT; INPUT_UNIT=in_unit; EYE=i and other numbers; TO_LOWERCASE; TO_UPPERCASE;... We thereby use ZERO instead of 0.0_real64
  USE Elem_op_m
  IMPLICIT NONE


  TYPE                           :: Quantum_HO1D_t
    integer                      :: Nb      = 0
    real(kind=Rkind)             :: w       = ZERO
    real(kind=Rkind)             :: m       = ZERO
    integer                      :: Nb_op   = 0
    TYPE(Elem_op_t), allocatable :: Tab_op(:)                                                           ! 0 : \hat{Id} ; 1 : \hat{H} ; 2 : \hat{x} ; 3 : \hat{N} ; 4 : \hat{\mu_{mat}} ; 5 : \hat{who knows ?}
    integer                      :: Nq      = 0
    real(kind=Rkind)             :: Eq_pos  = -ONE
    real(kind=Rkind)             :: Scale_q = HUGE(ONE)
  END TYPE


  PRIVATE 

  PUBLIC Quantum_HO1D_t, Initialize, Action, Get, Write, Dealloc

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
  INTERFACE Get
    MODULE PROCEDURE MolecCav_Get_QHO1D_parameter_integer, MolecCav_Get_QHO1D_parameter_real
  END INTERFACE
  INTERFACE Write
    MODULE PROCEDURE MolecCav_Write_quantum_HO1D
  END INTERFACE
  INTERFACE Dealloc
    MODULE PROCEDURE MolecCav_Deallocate_quantum_HO1D
  END INTERFACE
    

  CONTAINS


  SUBROUTINE MolecCav_Initialize_quantum_HO1D(QHO1D, Nb, w, m, Nb_op, Dense, Verbose, Debug)            ! here init on the 1D QHO basis : don't need Nq etc. will write later Init_grid or Init_other_basis, maybe called by this one, in this case, the "call" will be determined by the optional arg (Nb, Nq etc.)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
    USE QDUtil_m
    USE Elem_op_m
    IMPLICIT NONE
  
    TYPE(Quantum_HO1D_t), intent(inout) :: QHO1D
    integer,              intent(in)    :: Nb
    real(kind=Rkind),     intent(in)    :: w 
    real(kind=Rkind),     intent(in)    :: m 
    integer, optional,    intent(in)    :: Nb_op
    logical, optional,    intent(in)    :: Dense                                                        ! cf. comments in HO1D_parameters_m
    integer, optional,    intent(in)    :: Verbose                                                      ! cf. comments in HO1D_parameters_m
    logical, optional,    intent(in)    :: Debug                                                        ! cf. comments in HO1D_parameters_m

    logical                             :: Dense_local
    integer                             :: Verbose_local                                                ! goes from 20 (= 0 verbose) to 24 (= maximum verbose) at this layer
    logical                             :: Debug_local

    !--- Debugging options --------------------------------
    IF (PRESENT(Verbose)) THEN; Verbose_local = Verbose
    ELSE; Verbose_local = 16; END IF 
    IF (PRESENT(Debug)) THEN; Debug_local = Debug
    ELSE; Debug_local = .FALSE.; END IF
    IF (Debug_local) Verbose_local = 24

    IF (Verbose_local > 16) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "o o o Arguments of MolecCav_Initialize_quantum_HO1D :"
      WRITE(out_unit,*) "The <<QHO1D>> argument :"
      CALL Write(QHO1D, More=Debug_local)
      WRITE(out_unit,*) "The <<Nb>> argument    : "//TO_string(Nb)
      WRITE(out_unit,*) "The <<w>>  argument    : "//TO_string(w)
      WRITE(out_unit,*) "The <<m>>  argument    : "//TO_string(m)
      IF (PRESENT(Dense)) WRITE(out_unit,*) "The <<Dense>> argument : "//TO_string(Dense)
      FLUSH(out_unit)
    END IF
    
    !--- Initializing QHO1D parameters --------------------
    QHO1D%Nb = Nb
    QHO1D%w  = w
    QHO1D%m  = m
    IF (PRESENT(Nb_op)) THEN; QHO1D%Nb_op  = Nb_op
    ELSE; QHO1D%Nb_op  = 4; END IF 
    ALLOCATE(QHO1D%Tab_op(0:QHO1D%Nb_op-1))
    
    !--- Constructing the operators to build --------------
    IF (PRESENT(Dense)) THEN; Dense_local = Dense
    ELSE; Dense_local = .FALSE.; END IF

    IF (Verbose_local > 18) WRITE(out_unit,*) "--- Initialising the Elem_op objects associated to the QHO1D's operators..."
    CALL Initialize(QHO1D%Tab_op(0), "Identity",    Nb=Nb,           Dense=Dense_local, Verbose=Verbose_local, Debug=Debug_local)
    CALL Initialize(QHO1D%Tab_op(1), "Hamiltonian", Nb=Nb, w=w,      Dense=Dense_local, Verbose=Verbose_local, Debug=Debug_local)
    CALL Initialize(QHO1D%Tab_op(2), "Position",    Nb=Nb, w=w, m=m, Dense=Dense_local, Verbose=Verbose_local, Debug=Debug_local)
    CALL Initialize(QHO1D%Tab_op(3), "NbQuanta",    Nb=Nb,           Dense=Dense_local, Verbose=Verbose_local, Debug=Debug_local)
    IF (Verbose_local > 17) WRITE(out_unit,*) "    ...back to MolecCav_Initialize_quantum_HO1D"

    !--- Conclusion ---------------------------------------
    IF (Verbose_local > 17) THEN
      WRITE(out_unit,*) "--- Initialised QHO1D object :"
      CALL Write(QHO1D, More=.TRUE.)
    END IF

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

    integer                                   :: Verbose_local
    logical                                   :: Debug_local

    !--- Debugging options --------------------------------
    IF (PRESENT(Verbose)) THEN; Verbose_local = Verbose
    ELSE; Verbose_local = 16; END IF 
    IF (PRESENT(Debug)) THEN; Debug_local = Debug
    ELSE; Debug_local = .FALSE.; END IF
    IF (Debug_local) Verbose_local = 24

    IF (Verbose_local > 16) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "o o o Arguments of MolecCav_Initialize_HO1D_operator :"
      WRITE(out_unit,*) "The <<Elem_op>> argument       :"
      CALL Write(Elem_op)
      WRITE(out_unit,*) "The <<Operator_type>> argument : "//Operator_type
      WRITE(out_unit,*) "The <<Nb>> argument            : "//TO_string(Nb)
      IF (PRESENT(w)) WRITE(out_unit,*) "The <<w>> argument             : "//TO_string(w)
      IF (PRESENT(m)) WRITE(out_unit,*) "The <<m>> argument             : "//TO_string(m)
      IF (PRESENT(Dense)) WRITE(out_unit,*) "The <<Dense>> argument         : "//TO_string(Dense)
      FLUSH(out_unit)
    END IF

    IF (.NOT. PRESENT(w) .AND. TRIM(TO_lowercase(Operator_type)) == "hamiltonian") THEN
      WRITE(out_unit,*) "### Missing w for the H at MolecCav_Initialize_QHO1D_Elem_op"
      WRITE(out_unit,*) "    Please check arguments"
      STOP "### Missing w for the H"
    END IF 
    IF ((.NOT. PRESENT(w) .OR. .NOT. PRESENT(m)) .AND. TRIM(TO_lowercase(Operator_type)) == "position") THEN
      WRITE(out_unit,*) "### Missing w or m for the x at MolecCav_Initialize_QHO1D_Elem_op"
      WRITE(out_unit,*) "    Please check arguments"
      STOP "### Missing w or m for the x"
    END IF 
    
    !--- First steps of the construction of the Operator --
    ALLOCATE(character(len=LEN_TRIM(Operator_type)) :: Elem_op%Operator_type)                                                   ! /!\ strings cannot be allocated the exact same way as tables ! /!\
    Elem_op%Operator_type = TO_lowercase(TRIM(Operator_type))                                                                   ! allocation on assignement (not anymore : supposed to work but caused dynamic allocation random errors at execution). Elem_op_type has the right lengths (no spatials added) thanks to len=* at declaration and it will fit the Op%op_type thanks to len=:, allocatable at declaration of the derived type. 

    IF (PRESENT(Dense)) Elem_op%Dense = Dense

    !--- Construction of the matrix Operator --------------
    IF (Debug_local) THEN
      WRITE(out_unit,*); WRITE(out_unit,*) "--- The Elem_op_t object before constructing its matrix representation :"
      CALL Write(Elem_op)
    END IF 

    SELECT CASE (Elem_op%Operator_type)                                                                                         ! TO_lowercase avoid case sensitivity issues
      CASE ("identity")
        IF (Verbose_local > 18) WRITE(out_unit,*) "--- Initialising the Identity Elem_op inside QHO1D%Tab_op..."
        CALL Initialize_I(Identity=Elem_op,    Nb=Nb,           Verbose=Verbose_local, Debug=Debug_local)
        IF (Verbose_local > 18) WRITE(out_unit,*) "    ...back to MolecCav_Initialize_HO1D_operator"
    
      CASE ("hamiltonian")
        IF (Verbose_local > 18) WRITE(out_unit,*) "--- Initialising the Hamiltonian Elem_op inside QHO1D%Tab_op..."
        CALL Initialize_H(Hamiltonian=Elem_op, Nb=Nb, w=w,      Verbose=Verbose_local, Debug=Debug_local)
        IF (Verbose_local > 18) WRITE(out_unit,*) "    ...back to MolecCav_Initialize_HO1D_operator"
    
      CASE ("position")
        IF (Verbose_local > 18) WRITE(out_unit,*) "--- Initialising the Position Elem_op inside QHO1D%Tab_op..."
        CALL Initialize_x(Position=Elem_op,  Nb=Nb, w=w, m=m, Verbose=Verbose_local, Debug=Debug_local)
        IF (Verbose_local > 18) WRITE(out_unit,*) "    ...back to MolecCav_Initialize_HO1D_operator"
      
      CASE ("nbquanta")
        IF (Verbose_local > 18) WRITE(out_unit,*) "--- Initialising the Number of excitation quanta Elem_op inside QHO1D%Tab_op..."
        CALL Initialize_N(NbQuanta=Elem_op,    Nb=Nb,           Verbose=Verbose_local, Debug=Debug_local)
        IF (Verbose_local > 18) WRITE(out_unit,*) "    ...back to MolecCav_Initialize_HO1D_operator"

      CASE DEFAULT
        WRITE(out_unit,*) "### No Operator type recognized, please check the input of Initialize_HO1D_operator subroutine"
        STOP "### No Operator type recognized, please verify the input of Initialize_HO1D_operator subroutine"
    END SELECT

    IF (Verbose_local > 17) THEN
      WRITE(out_unit,*) "--- The initialised Elem_op of the QHO1D :"
      CALL Write(Elem_op)
    END IF

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
    integer                          :: Verbose_local
    logical                          :: Debug_local

    !--- Debugging options --------------------------------
    IF (PRESENT(Verbose)) THEN; Verbose_local = Verbose
    ELSE; Verbose_local = 16; END IF 
    IF (PRESENT(Debug)) THEN; Debug_local = Debug
    ELSE; Debug_local = .FALSE.; END IF
    IF (Debug_local) Verbose_local = 24

    IF (Verbose_local > 16) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "o o o Arguments of MolecCav_Initialize_I_QHO1D :"
      WRITE(out_unit,*) "The <<Identity>> argument      :"
      CALL Write(Identity, Info="QHO1D%Tab(0)")
      WRITE(out_unit,*) "The <<Nb>> argument            : "//TO_string(Nb)

    END IF 

    !--- Construction of the operator's matrix ------------
    IF (.NOT. Identity%Dense) THEN
      IF (Verbose_local > 18) WRITE(out_unit,*) "--- The Dense parameter is passed as .FALSE., the QHO1D Identity's matrix repres&
      &entation will be a rank-0 tensor of value ONE"
      !--- Initialization to default values ---------------
      ALLOCATE(Identity%Diag_val(1))

      !--- Construction of the matrix ---------------------
      Identity%Diag_val = ONE                                                                           ! "-1" because the first Fortran vector is the fundamental eigenvector of the HO i.e. the 0^{th} ket 

      IF (Debug_local) CALL Write_Vec(Identity%Diag_val, out_unit, 1, info="QHO1DIdentity")

    ELSE
      IF (Verbose_local > 18) WRITE(out_unit,*) "--- The Dense parameter is passed as .TRUE., the full QHO1D Identity's matrix wi&
      &ll be constructed for the representation, regardless of the sparce character of the analytical matrix."
      !---------------------------------------------Initialization to default values---------------------------------------------
      ALLOCATE(Identity%Dense_val(Nb, Nb))
      Identity%Dense_val = ZERO

      !------------------------------------------------Construction of the matrix------------------------------------------------
      DO i = 1, Nb                                                                                                     ! /!\ Fortran counts from 1 to Nb !!! /!\
        Identity%Dense_val(i,i) = ONE                                                              ! "-1" because the first Fortran vector is the fundamental eigenvector of the HO i.e. the 0^{th} ket 
      END DO

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
    integer                          :: Verbose_local
    logical                          :: Debug_local

    !--- Debugging options --------------------------------
    IF (PRESENT(Verbose)) THEN; Verbose_local = Verbose
    ELSE; Verbose_local = 16; END IF 
    IF (PRESENT(Debug)) THEN; Debug_local = Debug
    ELSE; Debug_local = .FALSE.; END IF
    IF (Debug_local) Verbose_local = 24

    IF (Verbose_local > 16) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "o o o Arguments of MolecCav_Initialize_H_QHO1D :"
      WRITE(out_unit,*) "The <<Hamiltonian>> argument   :"
      CALL Write(Hamiltonian, Info="Hamiltonian")
      WRITE(out_unit,*) "The <<Nb>> argument            : "//TO_string(Nb)

    END IF 

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

      IF (Debug_local) CALL Write_Mat(Hamiltonian%Dense_val, out_unit, Size(Hamiltonian%Dense_val), info="HO1DHamiltonian")
    END IF
      
  END SUBROUTINE MolecCav_Initialize_H_QHO1D


  SUBROUTINE MolecCav_Initialize_x_QHO1D(Position, Nb, w, m, Verbose, Debug)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT, OUTPUT_UNIT, real64 
    USE QDUtil_m
    USE Elem_op_m
   IMPLICIT NONE
    
    TYPE(Elem_op_t),   intent(inout) :: Position
    integer,           intent(in)    :: Nb                                                                            ! cf. comments in HO1D_parameters_m
    real(kind=Rkind),  intent(in)    :: w                                                                            ! cf. comments in HO1D_parameters_m
    real(kind=Rkind),  intent(in)    :: m                                                                            ! cf. comments in HO1D_parameters_m
    integer, optional, intent(in)    :: Verbose                                                                            ! cf. comments in HO1D_parameters_m
    logical, optional, intent(in)    :: Debug                                                                              ! cf. comments in HO1D_parameters_m

    integer                          :: i                                                                                  ! loop increments /!\ Fortran counts from 1 to Nb !!! /!\
    integer                          :: Verbose_local
    logical                          :: Debug_local

    !--- Debugging options --------------------------------
    IF (PRESENT(Verbose)) THEN; Verbose_local = Verbose
    ELSE; Verbose_local = 16; END IF 
    IF (PRESENT(Debug)) THEN; Debug_local = Debug
    ELSE; Debug_local = .FALSE.; END IF
    IF (Debug_local) Verbose_local = 24

    IF (Verbose_local > 16) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "o o o Arguments of MolecCav_Initialize_x_QHO1D :"
      WRITE(out_unit,*) "The <<Position>> argument      :"
      CALL Write(Position, Info="Position")
      WRITE(out_unit,*) "The <<Nb>> argument            : "//TO_string(Nb)

    END IF 
    
    !---------------------------------------------Construction of the matrix Operator--------------------------------------------
    IF ((.NOT. Position%Dense) .AND. Nb > 1) THEN
      IF (Verbose_local > 28) WRITE(out_unit,*) "--- The Dense parameter of the Position operator is .FALSE., so the 1D HO Positi&
                                                &on operator's matrix representation will be a rank-2 tensor of the tridiagonal e&
                                                &lements of its analytical matrix (in Eigenbasis)"
      !-----------------------------------Initialization of the characteristics of the operator----------------------------------
      Position%Upper_bandwidth   = 1
      Position%Lower_bandwidth   = 1

      !---------------------------------------------Initialization to default values---------------------------------------------
      ALLOCATE(Position%Band_val(Nb,3))                                                                            ! Nb lines (number of diagonal elements) and 3 columns because 3 bands to consider : the diagonal, and the two bands above and below it
      Position%Band_val = ZERO

      !------------------------------------------------Construction of the matrix------------------------------------------------
      DO i = 1, Nb - 1                                                                                                 ! /!\ Fortran counts from 1 to Nb !!! /!\ Nb-1 not to have Band_val(i+1) out of range
        Position%Band_val(i,1)   = SQRT(REAL(i,kind=Rkind))
        Position%Band_val(i+1,3) = SQRT(REAL(i,kind=Rkind))
      END DO
      Position%Band_val = Position%Band_val / SQRT(TWO * w * m)
    
      IF (Debug_local) CALL Write_Mat(Position%Band_val, out_unit, 3, info="HO1DPosition")

    ELSE IF (.NOT. Position%Dense) THEN
      IF (Verbose_local > 28) WRITE(out_unit,*) "--- The Dense parameter of the Position operator is .FALSE. BUT the basis set si&
                                                &ze is only 1, so the 1D HO Position operator's matrix representation will use th&
                                                &e Diag_val rank-1 tensor to store the only element of the analytical matrix (i&
                                                &n Eigenbasis)"
      !---------------------------------------------Initialization to default values---------------------------------------------
      ALLOCATE(Position%Diag_val(Nb))

      !------------------------------------------------Construction of the matrix------------------------------------------------
      DO i = 1, Nb                                                                                                     ! /!\ Fortran counts from 1 to Nb !!! /!\
        Position%Diag_val(i) = ZERO                                                                                          ! the position operator matrix has first value (i.e. only value in the Nb = 0 case) 0 
      END DO

      IF (Debug_local) CALL Write_Vec(Position%Diag_val, out_unit, 1, info="HO1DPosition")

    ELSE 
      IF (Verbose_local > 28) WRITE(out_unit,*) "--- The Dense parameter of the Position operator is .TRUE., so the full 1D HO Po&
                                                &sition operator's matrix will be constructed (in Eigenbasis) for the representat&
                                                &ion, as if the analytical matrix was a dense one"
      !---------------------------------------------Initialization to default values---------------------------------------------
      ALLOCATE(Position%Dense_val(Nb, Nb))
      Position%Dense_val = ZERO
      
      !------------------------------------------------Construction of the matrix------------------------------------------------
      DO i = 1, Nb - 1                                                                                                 ! /!\ Fortran counts from 1 to Nb !!! /!\
        Position%Dense_val(i,i+1) = SQRT(REAL(i,kind=Rkind))
        Position%Dense_val(i+1,i) = SQRT(REAL(i,kind=Rkind))
      END DO
      Position%Dense_val = Position%Dense_val / SQRT(TWO * w * m)
    
      IF (Debug_local) CALL Write_Mat(Position%Dense_val, out_unit, Size(Position%Dense_val), info="HO1DPosition")
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
    integer                          :: Verbose_local
    logical                          :: Debug_local

    !--- Debugging options --------------------------------
    IF (PRESENT(Verbose)) THEN; Verbose_local = Verbose
    ELSE; Verbose_local = 16; END IF 
    IF (PRESENT(Debug)) THEN; Debug_local = Debug
    ELSE; Debug_local = .FALSE.; END IF
    IF (Debug_local) Verbose_local = 24

    IF (Verbose_local > 16) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "o o o Arguments of MolecCav_Initialize_N_QHO1D :"
      WRITE(out_unit,*) "The <<NbQuanta>> argument      :"
      CALL Write(NbQuanta, Info="NbQuanta")
      WRITE(out_unit,*) "The <<Nb>> argument            : "//TO_string(Nb)
    END IF 

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
    integer                             :: Verbose_local
    logical                             :: Debug_local

    !--- Debugging options --------------------------------
    IF (PRESENT(Verbose)) THEN; Verbose_local = Verbose
    ELSE; Verbose_local = 16; END IF 
    IF (PRESENT(Debug)) THEN; Debug_local = Debug
    ELSE; Debug_local = .FALSE.; END IF
    IF (Debug_local) Verbose_local = 24

    IF (Verbose_local > 16) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "o o o Arguments of MolecCav_Action_quantum_HO1D :"
      WRITE(out_unit,*) "The <<QHO1D>> argument :"
      CALL Write(QHO1D, More=Debug_local)
      WRITE(out_unit,*) "The <<i_op>> argument  : "//TO_string(i_op)
      WRITE(out_unit,*) "The <<Psi>> argument   : "
      CALL Write_Vec(Psi, out_unit, 1, info="Psi")
      WRITE(out_unit,*) "The size of its vector : "//TO_string(Size(Psi))
      FLUSH(out_unit)
    END IF
    
    !-----------------------------------------------------Checking dimensions----------------------------------------------------
    ! ALREADY CHECKED IN THE ACTIONS CODED IN ELEM_OP_M ! (fortunately btw, otherwise the \hat{I}d case should have test above)

    !---------------------------------------------Selection of the calculation method--------------------------------------------
    IF (QHO1D%Tab_op(i_op)%Operator_type == "identity") THEN ! N.B. we should have tested i_op == 0, easier and consistent with the algorithmic choices made so far
      Op_psi = Psi
    ELSE
      CALL Action(Op_psi=Op_psi, Elem_op=QHO1D%Tab_op(i_op), Psi=Psi, Verbose=Verbose_local, Debug=Debug_local)
      WRITE(out_unit,*) "    ...back to MolecCav_Action_quantum_HO1D"
    END IF

    IF (Debug_local) THEN
      WRITE(out_unit,*) "--- Resulting statevector from the action of the HO1D Elem_op on the Psi statevector operand, computed &
                        &by Action_HO1D_operator_R1 :"
      CALL Write_Vec(Op_psi, out_unit, 1, info="Op_Psi")
    END IF
    
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
    integer                             :: Verbose_local
    logical                             :: Debug_local

    !--- Debugging options --------------------------------
    IF (PRESENT(Verbose)) THEN; Verbose_local = Verbose
    ELSE; Verbose_local = 16; END IF 
    IF (PRESENT(Debug)) THEN; Debug_local = Debug
    ELSE; Debug_local = .FALSE.; END IF
    IF (Debug_local) Verbose_local = 24

    IF (Verbose_local > 16) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "o o o Arguments of MolecCav_Action_quantum_HO1D :"
      WRITE(out_unit,*) "The <<QHO1D>> argument  :"
      CALL Write(QHO1D, More=Debug_local)
      WRITE(out_unit,*) "The <<i_op>> argument   : "//TO_string(i_op)
      WRITE(out_unit,*) "The <<Psi>> argument    :"
      CALL Write_Vec(Psi, out_unit, 1, info="Psi")
      WRITE(out_unit,*) "The size of its vector  : "//TO_string(Size(Psi))
      FLUSH(out_unit)
    END IF
    
    !-----------------------------------------------------Checking dimensions----------------------------------------------------
    ! ALREADY CHECKED IN THE ACTIONS CODED IN ELEM_OP_M !

    !---------------------------------------------Selection of the calculation method--------------------------------------------
    IF (QHO1D%Tab_op(i_op)%Operator_type == "identity") THEN 
      Op_psi = Psi
    ELSE
      CALL Action(Op_psi=Op_psi, Elem_op=QHO1D%Tab_op(i_op), Psi=Psi, Verbose=Verbose_local, Debug=Debug_local)
      WRITE(out_unit,*) "    ...back to MolecCav_Action_quantum_HO1D_R1_complex"
    END IF

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "--- Resulting statevector from the action of the HO1D Elem_op on the Psi statevector operand, computed &
                        &by Action_HO1D_operator_R1 :"
      CALL Write_Vec(Op_psi, out_unit, 1, info="Op_Psi")
    END IF
    
  END SUBROUTINE MolecCav_Action_quantum_HO1D_R1_complex

  
  SUBROUTINE MolecCav_Get_QHO1D_parameter_integer(Parameter_value, QHO1D, Parameter_name)
    USE QDUtil_m
    IMPLICIT NONE 

    integer,              intent(inout) :: Parameter_value                                                                            ! the current values of the indexes for each dimension
    TYPE(Quantum_HO1D_t), intent(in)    :: QHO1D
    character(len=*),     intent(in)    :: Parameter_name

    SELECT CASE (TO_lowercase(TRIM(Parameter_name)))                         ! TO_lowercase avoid case sensitivity issues
      CASE ("nb")
        Parameter_value = QHO1D%Nb
  
      CASE ("nb_op")
        Parameter_value = QHO1D%Nb_op
    
      CASE ("nq")
        Parameter_value = QHO1D%Nq

      CASE DEFAULT
        WRITE(out_unit,*) "No Parameter name recognized, please verify the input of Get_QHO1D_parameter_integer subroutine"
        STOP "### No Operator type recognized, please verify the input of Get_QHO1D_parameter_integer subroutine"

    END SELECT

  END SUBROUTINE MolecCav_Get_QHO1D_parameter_integer


  SUBROUTINE MolecCav_Get_QHO1D_parameter_real(Parameter_value, QHO1D, Parameter_name)
    USE QDUtil_m
    IMPLICIT NONE 

    real(kind=Rkind),     intent(inout) :: Parameter_value                                                                            ! the current values of the indexes for each dimension
    TYPE(Quantum_HO1D_t), intent(in)    :: QHO1D
    character(len=*),     intent(in)    :: Parameter_name

    SELECT CASE (TO_lowercase(TRIM(Parameter_name)))                         ! TO_lowercase avoid case sensitivity issues
      CASE ("w")
        Parameter_value = QHO1D%w
  
      CASE ("m")
        Parameter_value = QHO1D%m
    
      CASE ("eq_pos")
        Parameter_value = QHO1D%Eq_pos

      CASE ("scale_q")
        Parameter_value = QHO1D%Scale_q

      CASE DEFAULT
        WRITE(out_unit,*) "No Parameter name recognized, please verify the input of Get_QHO1D_parameter_integer subroutine"
        STOP "### No Operator type recognized, please verify the input of Get_QHO1D_parameter_integer subroutine"

    END SELECT

  END SUBROUTINE MolecCav_Get_QHO1D_parameter_real


  SUBROUTINE MolecCav_Write_quantum_HO1D(QHO1D, Info, More)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
    USE QDUtil_m
    USE Elem_op_m
    IMPLICIT NONE 
    
    TYPE(Quantum_HO1D_t),       intent(in) :: QHO1D
    character(len=*), optional, intent(in) :: Info
    logical,          optional, intent(in) :: More

    integer                                :: i_op
    logical                                :: More_local

    IF (PRESENT(More)) THEN; More_local = More
    ELSE; More_local = .FALSE.; END IF

    IF (PRESENT(Info)) THEN
      WRITE(out_unit,*) "--- Parameters associated to the object of derived type Quantum_HO1D_t ("//Info//") :"
    ELSE
      WRITE(out_unit,*) "--- Parameters associated to the object of derived type Quantum_HO1D_t :"
    END IF

    WRITE(out_unit,*) "Nb                  = "//TO_string(QHO1D%Nb)
    WRITE(out_unit,*) "w                   = "//TO_string(QHO1D%w)
    WRITE(out_unit,*) "m                   = "//TO_string(QHO1D%m)
    WRITE(out_unit,*) "Nb_op               = "//TO_string(QHO1D%Nb_op)

    IF (ALLOCATED(QHO1D%Tab_op)) THEN 
      WRITE(out_unit,*) "SIZE(Tab_op, dim=1) = "//TO_string(SIZE(QHO1D%Tab_op, dim=1))
      IF (More_local) THEN 
        WRITE(out_unit,*) "--- Writing Tab_op..."
        DO i_op = 0, SIZE(QHO1D%Tab_op)-1
          CALL Write(QHO1D%Tab_op(i_op), Info="Tab_op("//TO_string(i_op)//")")
        END DO 
        WRITE(out_unit,*) "    ...back to MolecCav_Write_quantum_HO1D"
      END IF 
    ELSE 
      WRITE(out_unit,*) "Tab_op is NOT allocated"
    END IF

    WRITE(out_unit,*) "Nq                  = "//TO_string(QHO1D%Nq)
    WRITE(out_unit,*) "Eq_pos              = "//TO_string(QHO1D%Eq_pos)
    WRITE(out_unit,*) "Scale_q             = "//TO_string(QHO1D%Scale_q)
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
    integer                             :: Verbose_local
    logical                             :: Debug_local

    !--- Debugging options --------------------------------
    IF (PRESENT(Verbose)) THEN; Verbose_local = Verbose
    ELSE; Verbose_local = 16; END IF 
    IF (PRESENT(Debug)) THEN; Debug_local = Debug
    ELSE; Debug_local = .FALSE.; END IF
    IF (Debug_local) Verbose_local = 24

    IF (Verbose_local > 16) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "o o o Arguments of MolecCav_Deallocate_quantum_HO1D :"
      CALL Write(QHO1D, More=Debug_local)
    END IF 

    !-----------------------------Deallocating the HO1D operator object----------------------------  
    QHO1D%Nb = 0
    QHO1D%w  = ZERO
    QHO1D%m  = ZERO
    DO i_op = 0, SIZE(QHO1D%Tab_op)-1
      CALL Dealloc(QHO1D%Tab_op(i_op), Verbose=Verbose_local, Debug=Debug_local)
      WRITE(out_unit,*) "    ...back to MolecCav_Deallocate_quantum_HO1D"
    END DO
    IF (ALLOCATED(QHO1D%Tab_op)) DEALLOCATE(QHO1D%Tab_op)
    QHO1D%Nq      = 0
    QHO1D%Eq_pos  = -ONE
    QHO1D%Scale_q = HUGE(ONE)

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "--- The deallocated QHO1D object :"
      CALL Write(QHO1D, More=Debug_local)
    END IF

  END SUBROUTINE MolecCav_Deallocate_quantum_HO1D


END MODULE