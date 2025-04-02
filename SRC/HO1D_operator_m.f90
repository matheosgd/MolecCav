!==================================================================================================
!==================================================================================================
! This file is part of MolecCav.
!
!==================================================================================================
! MIT License
!
! Copyright (c) 2525 Mathéo Segaud
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
! The module to initialize the operators related to a HO.
! Initialize_HO1D_operator : constructs the operator using parameters of the HO1D_para object from 
! the so called derived type. Calls the next three procedures. 
! Initialize_H_HO1D        : constructs the Hamiltonian operator using parameters of the HO1D_para 
! object from the so called derived type.  
! Initialize_x_HO1D        : constructs the Position operator using parameters of the HO1D_para ob-
! ject from the so called derived type.
! Initialize_N_HO1D        : constructs the Number of excitation Quanta operator using parameters 
! of the HO1D_para object from the so called derived type.
! Write_HO1D_operator      : displays values of the type in the output.
! Deallocate_HO1D_operator : deallocates all tables of the type.
!==================================================================================================
!==================================================================================================
MODULE HO1D_operator_m
  !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
  USE QDUtil_m
  IMPLICIT NONE


  TYPE                               :: HO1D_operator_t
    character(len=:),    allocatable :: Operator_type                                                                            ! ex : "Hamiltonian", "Position", "NbQuanta", etc
    logical                          :: Dense = .FALSE.                                                                          ! if .TRUE. then the matrix storage will not be optimized and it will be stored as a Dense matrix. Otherwise, the elements matrix will be stored in a table such as taking advantage of the sparcity of the operator's analytical matrix in the HO Eigenbasis
    integer                          :: Upper_bandwidth = 0                                                                      ! if type = "Band". Gives the number of additional bands to consider above the diagonal.
    integer                          :: Lower_bandwidth = 0                                                                      ! if type = "Band". Gives the number of additional bands to consider below the diagonal. Ex : Upper_bandwidth=Lower_bandwidth=1 would give a tridiagonal matrix
    real(kind=Rkind),    allocatable :: Dense_val_R(:,:)                                                                         ! if Dense == .TRUE.
    real(kind=Rkind),    allocatable :: Diag_val_R(:)                                                                            ! if Dense == .FALSE. .AND. the operator analytical matrix in the HO1D Eigenbasis is diagonal. The diagonal elements are stored in a vector (rank-1 tensor)
    real(kind=Rkind),    allocatable :: Band_val_R(:,:)                                                                          ! if Dense == .FALSE. .AND. the operator analytical matrix in the HO1D Eigenbasis is band. The number of columns will be the number of diagonals to consider : each considered diagonal is stored in a column
  END TYPE
  
  
  PRIVATE

  PUBLIC HO1D_operator_t, Initialize_HO1D_operator, Write_HO1D_operator, Deallocate_HO1D_operator
 
  INTERFACE Initialize_HO1D_operator
    MODULE PROCEDURE MolecCav_Initialize_HO1D_operator
  END INTERFACE
  INTERFACE Initialize_H_HO1D
    MODULE PROCEDURE MolecCav_Initialize_H_HO1D
  END INTERFACE
  INTERFACE Initialize_x_HO1D
    MODULE PROCEDURE MolecCav_Initialize_x_HO1D
  END INTERFACE
  INTERFACE Initialize_N_HO1D
    MODULE PROCEDURE MolecCav_Initialize_N_HO1D
  END INTERFACE
  INTERFACE Write_HO1D_operator
    MODULE PROCEDURE MolecCav_Write_HO1D_operator
  END INTERFACE
  INTERFACE Deallocate_HO1D_operator
    MODULE PROCEDURE MolecCav_Deallocate_HO1D_operator
  END INTERFACE
    

  CONTAINS


  SUBROUTINE MolecCav_Initialize_HO1D_operator(Operator, Operator_type, HO1D_para, Dense, Verbose, Debug)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
    USE QDUtil_m
    USE HO1D_parameters_m
    IMPLICIT NONE

    TYPE(HO1D_operator_t),   intent(inout) :: Operator                                                                           ! the object of type Operator_t to be constructed here
    character(len=*),        intent(in)    :: Operator_type                                                                      ! ex : "Hamiltonian", "Position", etc. (len=:) Expects to be allocatable, while (len=*) is dedicated to a procedure argument.
    TYPE(HO1D_parameters_t), intent(in)    :: HO1D_para                                                                          ! the HO/Cavity mode which the operator is relative to
    logical, optional,       intent(in)    :: Dense                                                                              ! if .TRUE. then the matrix storage will not be optimized and it will be stored as a Dense matrix
    integer, optional,       intent(in)    :: Verbose                                                                            ! cf. comments in HO1D_parameters_m
    logical, optional,       intent(in)    :: Debug                                                                              ! cf. comments in HO1D_parameters_m

    integer                                :: Verbose_local = 25                                                                 ! goes from 25 (= 0 verbose) to 29 (= maximum verbose) at this layer
    logical                                :: Debug_local   = .FALSE.

    !------------------------------------------------------Debugging options-----------------------------------------------------
    IF (PRESENT(Verbose)) Verbose_local = Verbose
    IF (PRESENT(Debug))   Debug_local   = Debug

    IF (Verbose_local > 25) WRITE(out_unit,*) 
    IF (Verbose_local > 25) WRITE(out_unit,*) "--------------------------------------------------INITIALIZING THE HO1D OPERATOR--&
                                              &------------------------------------------------"; FLUSH(out_unit)

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "--- Arguments of MolecCav_Initialize_HO1D_operator :"
      WRITE(out_unit,*) "The <<Operator>> argument :"
      CALL Write_HO1D_operator(Operator)
      WRITE(out_unit,*) "The <<Operator_type>> argument : "//Operator_type
      WRITE(out_unit,*) "The <<HO1D_para>> argument :"
      CALL Write_HO1D_parameters(HO1D_para)
      IF (PRESENT(Dense)) WRITE(out_unit,*) "The <<Dense>> argument : "//TO_string(Dense)
      WRITE(out_unit,*) "--- End arguments of MolecCav_Construct_Operator_1D"
      FLUSH(out_unit)
    END IF
    
    !---------------------------------------First steps of the construction of the Operator--------------------------------------
    ALLOCATE(character(len=LEN_TRIM(Operator_type)) :: Operator%Operator_type)                                                   ! /!\ strings cannot be allocated the exact same way as tables ! /!\
    Operator%Operator_type = TO_lowercase(TRIM(Operator_type))                                                                   ! allocation on assignement (not anymore : supposed to work but caused dynamic allocation random errors at execution). Operator_type has the right lengths (no spaces added) thanks to len=* at declaration and it will fit the Op%op_type thanks to len=:, allocatable at declaration of the derived type. 

    IF (PRESENT(Dense)) Operator%Dense = Dense

    !---------------------------------------------Construction of the matrix Operator--------------------------------------------
    IF (Debug_local) THEN
      WRITE(out_unit,*); WRITE(out_unit,*) "--- The HO1D_operator_t object just before construction of its matrix representation"
      CALL Write_HO1D_operator(Operator)
      WRITE(out_unit,*) "--- End HO1D_operator_t object (just before construction of its matrix representation)"
    END IF 

    SELECT CASE (Operator%Operator_type)                                                                                         ! TO_lowercase avoid case sensitivity issues
      CASE ("hamiltonian")
        CALL MolecCav_Initialize_H_HO1D(Hamiltonian=Operator, HO1D_para=HO1D_para, Verbose=Verbose_local, Debug=Debug_local)
    
      CASE ("position")
        CALL MolecCav_Initialize_x_HO1D(PositionOp=Operator,  HO1D_para=HO1D_para, Verbose=Verbose_local, Debug=Debug_local)
      
      CASE ("nbquanta")
        CALL MolecCav_Initialize_N_HO1D(NbQuanta=Operator,    HO1D_para=HO1D_para, Verbose=Verbose_local, Debug=Debug_local)

      CASE DEFAULT
        WRITE(out_unit,*) "### No Operator type recognized, please check the input of Initialize_HO1D_operator subroutine"
        STOP "### No Operator type recognized, please verify the input of Initialize_HO1D_operator subroutine"
    END SELECT

    IF (Verbose_local > 26) THEN
      IF (Verbose_local < 28) WRITE(out_unit,*)
      WRITE(out_unit,*) "--- HO1D operator constructed by MolecCav_Initialize_HO1D_operator :"
      CALL Write_HO1D_operator(Operator)
      WRITE(out_unit,*) "--- End HO1D operator constructed by MolecCav_Initialize_HO1D_operator"
    END IF

    IF (Verbose_local > 25) WRITE(out_unit,*) 
    IF (Verbose_local > 25) WRITE(out_unit,*) "-----------------------------------------------------HO1D OPERATOR INITIALIZED----&
                                              &------------------------------------------------"; FLUSH(out_unit)

  END SUBROUTINE MolecCav_Initialize_HO1D_operator


  SUBROUTINE MolecCav_Initialize_H_HO1D(Hamiltonian, HO1D_para, Verbose, Debug)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    USE HO1D_parameters_m
    IMPLICIT NONE
    
    TYPE(HO1D_operator_t),   intent(inout) :: Hamiltonian                                                                        ! matrix of the one-dimensional harmonic Hamiltonian associated with HO D
    TYPE(HO1D_parameters_t), intent(in)    :: HO1D_para                                                                          ! the HO/Cavity mode which the operator is relative to
    integer, optional,       intent(in)    :: Verbose                                                                            ! cf. comments in HO1D_parameters_m
    logical, optional,       intent(in)    :: Debug                                                                              ! cf. comments in HO1D_parameters_m

    integer                                :: i                                                                                  ! loop increments /!\ Fortran counts from 1 to Nb !!! /!\
    integer                                :: Verbose_local = 25                                                                 ! goes from 25 (= 0 verbose) to 29 (= maximum verbose) at this layer
    logical                                :: Debug_local   = .FALSE.

    !------------------------------------------------------Debugging options-----------------------------------------------------
    IF (PRESENT(Verbose)) Verbose_local = Verbose
    IF (PRESENT(Debug))   Debug_local   = Debug

    IF (Verbose_local > 27) WRITE(out_unit,*) 
    IF (Verbose_local > 27) WRITE(out_unit,*) "----------------------------------Constructing the matrix representation of the 1D&
                                              & HO Hamiltonian---------------------------------"
 
    !---------------------------------------------Construction of the matrix Operator--------------------------------------------
    IF (.NOT. Hamiltonian%Dense) THEN
      IF (Verbose_local > 28) WRITE(out_unit,*) "--- The Dense parameter of the Hamiltonian is .FALSE., so the 1D HO Hamiltonian'&
                                               &s matrix representation will be a rank-1 tensor of the diagonal elementsof its an&
                                               &alytical matrix (in Eigenbasis)"
      !---------------------------------------------Initialization to default values---------------------------------------------
      ALLOCATE(Hamiltonian%Diag_val_R(HO1D_para%Nb))

      !------------------------------------------------Construction of the matrix------------------------------------------------
      DO i = 1, HO1D_para%Nb                                                                                                     ! /!\ Fortran counts from 1 to Nb !!! /!\
        Hamiltonian%Diag_val_R(i) = HO1D_para%w*(i - ONE + HALF)                                                                 ! "-1" because the first Fortran vector is the fundamental eigenvector of the HO i.e. the 0^{th} ket 
      END DO

      IF (Debug_local) WRITE(out_unit,*)
      IF (Debug_local) CALL Write_Vec(Hamiltonian%Diag_val_R, out_unit, 1, info="HO1DHamiltonian")

    ELSE
      IF (Verbose_local > 28) WRITE(out_unit,*) "--- The Dense parameter of the Hamiltonian is .TRUE., so the full 1D HO Hamilton&
                                                &ian's matrix will be constructed (in Eigenbasis) for the representation, as if t&
                                                &he analytical matrix was a dense one"
      !---------------------------------------------Initialization to default values---------------------------------------------
      ALLOCATE(Hamiltonian%Dense_val_R(HO1D_para%Nb, HO1D_para%Nb))
      Hamiltonian%Dense_val_R = ZERO

      !------------------------------------------------Construction of the matrix------------------------------------------------
      DO i = 1, HO1D_para%Nb                                                                                                     ! /!\ Fortran counts from 1 to Nb !!! /!\
        Hamiltonian%Dense_val_R(i,i) = HO1D_para%w*(i - ONE + HALF)                                                              ! "-1" because the first Fortran vector is the fundamental eigenvector of the HO i.e. the 0^{th} ket 
      END DO

      IF (Debug_local) WRITE(out_unit,*)
      IF (Debug_local) CALL Write_Mat(Hamiltonian%Dense_val_R, out_unit, Size(Hamiltonian%Dense_val_R), info="HO1DHamiltonian")
    END IF
      
  END SUBROUTINE MolecCav_Initialize_H_HO1D


  SUBROUTINE MolecCav_Initialize_x_HO1D(PositionOp, HO1D_para, Verbose, Debug)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT, OUTPUT_UNIT, real64 
    USE QDUtil_m
     USE HO1D_parameters_m
   IMPLICIT NONE
    
    TYPE(HO1D_operator_t),   intent(inout) :: PositionOp
    TYPE(HO1D_parameters_t), intent(in)    :: HO1D_para                                                                          ! the HO/Cavity mode which the operator is relative to
    integer, optional,       intent(in)    :: Verbose                                                                            ! cf. comments in HO1D_parameters_m
    logical, optional,       intent(in)    :: Debug                                                                              ! cf. comments in HO1D_parameters_m

    integer                                :: i                                                                                  ! loop increments /!\ Fortran counts from 1 to Nb !!! /!\
    integer                                :: Verbose_local = 25                                                                 ! goes from 25 (= 0 verbose) to 29 (= maximum verbose) at this layer
    logical                                :: Debug_local   = .FALSE.

    !------------------------------------------------------Debugging options-----------------------------------------------------
    IF (PRESENT(Verbose)) Verbose_local = Verbose
    IF (PRESENT(Debug))   Debug_local   = Debug

    IF (Verbose_local > 27) WRITE(out_unit,*) 
    IF (Verbose_local > 27) WRITE(out_unit,*) "-------------------------------Constructing the matrix representation of the 1D HO&
                                             & Position operator------------------------------"

    !---------------------------------------------Construction of the matrix Operator--------------------------------------------
    IF ((.NOT. PositionOp%Dense) .AND. HO1D_para%Nb > 1) THEN
      IF (Verbose_local > 28) WRITE(out_unit,*) "--- The Dense parameter of the Position operator is .FALSE., so the 1D HO Positi&
                                                &on operator's matrix representation will be a rank-2 tensor of the tridiagonal e&
                                                &lements of its analytical matrix (in Eigenbasis)"
      !-----------------------------------Initialization of the characteristics of the operator----------------------------------
      PositionOp%Upper_bandwidth   = 1
      PositionOp%Lower_bandwidth   = 1

      !---------------------------------------------Initialization to default values---------------------------------------------
      ALLOCATE(PositionOp%Band_val_R(HO1D_para%Nb,3))                                                                            ! Nb lines (number of diagonal elements) and 3 columns because 3 bands to consider : the diagonal, and the two bands above and below it
      PositionOp%Band_val_R = ZERO

      !------------------------------------------------Construction of the matrix------------------------------------------------
      DO i = 1, HO1D_para%Nb - 1                                                                                                 ! /!\ Fortran counts from 1 to Nb !!! /!\ Nb-1 not to have Band_val_R(i+1) out of range
        PositionOp%Band_val_R(i,1)   = SQRT(REAL(i,kind=Rkind))
        PositionOp%Band_val_R(i+1,3) = SQRT(REAL(i,kind=Rkind))
      END DO
      PositionOp%Band_val_R = PositionOp%Band_val_R / SQRT(TWO * HO1D_para%w * HO1D_para%m)
    
      IF (Debug_local) WRITE(out_unit,*)
      IF (Debug_local) CALL Write_Mat(PositionOp%Band_val_R, out_unit, 3, info="HO1DPositionOp")

    ELSE IF (.NOT. PositionOp%Dense) THEN
      IF (Verbose_local > 28) WRITE(out_unit,*) "--- The Dense parameter of the Position operator is .FALSE. BUT the basis set si&
                                                &ze is only 1, so the 1D HO Position operator's matrix representation will use th&
                                                &e Diag_val_R rank-1 tensor to store the only element of the analytical matrix (i&
                                                &n Eigenbasis)"
      !---------------------------------------------Initialization to default values---------------------------------------------
      ALLOCATE(PositionOp%Diag_val_R(HO1D_para%Nb))

      !------------------------------------------------Construction of the matrix------------------------------------------------
      DO i = 1, HO1D_para%Nb                                                                                                     ! /!\ Fortran counts from 1 to Nb !!! /!\
        PositionOp%Diag_val_R(i) = ZERO                                                                                          ! the position operator matrix has first value (i.e. only value in the Nb = 0 case) 0 
      END DO

      IF (Debug_local) WRITE(out_unit,*)
      IF (Debug_local) CALL Write_Vec(PositionOp%Diag_val_R, out_unit, 1, info="HO1DPositionOp")

    ELSE 
      IF (Verbose_local > 28) WRITE(out_unit,*) "--- The Dense parameter of the Position operator is .TRUE., so the full 1D HO Po&
                                                &sition operator's matrix will be constructed (in Eigenbasis) for the representat&
                                                &ion, as if the analytical matrix was a dense one"
      !---------------------------------------------Initialization to default values---------------------------------------------
      ALLOCATE(PositionOp%Dense_val_R(HO1D_para%Nb, HO1D_para%Nb))
      PositionOp%Dense_val_R = ZERO
      
      !------------------------------------------------Construction of the matrix------------------------------------------------
      DO i = 1, HO1D_para%Nb - 1                                                                                                 ! /!\ Fortran counts from 1 to Nb !!! /!\
        PositionOp%Dense_val_R(i,i+1) = SQRT(REAL(i,kind=Rkind))
        PositionOp%Dense_val_R(i+1,i) = SQRT(REAL(i,kind=Rkind))
      END DO
      PositionOp%Dense_val_R = PositionOp%Dense_val_R / SQRT(TWO * HO1D_para%w * HO1D_para%m)
    
      IF (Debug_local) WRITE(out_unit,*)
      IF (Debug_local) CALL Write_Mat(PositionOp%Dense_val_R, out_unit, Size(PositionOp%Dense_val_R), info="HO1DPositionOp")
    END IF
      
  END SUBROUTINE MolecCav_Initialize_x_HO1D


  SUBROUTINE MolecCav_Initialize_N_HO1D(NbQuanta, HO1D_para, Verbose, Debug)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    USE HO1D_parameters_m
    IMPLICIT NONE
    
    TYPE(HO1D_operator_t),   intent(inout) :: NbQuanta
    TYPE(HO1D_parameters_t), intent(in)    :: HO1D_para                                                                          ! the HO/Cavity mode which the operator is relative to
    integer, optional,       intent(in)    :: Verbose                                                                            ! cf. comments in HO1D_parameters_m
    logical, optional,       intent(in)    :: Debug                                                                              ! cf. comments in HO1D_parameters_m

    integer                                :: i                                                                                  ! loop increments /!\ Fortran counts from 1 to Nb !!! /!\
    integer                                :: Verbose_local = 25                                                                 ! goes from 25 (= 0 verbose) to 29 (= maximum verbose) at this layer
    logical                                :: Debug_local   = .FALSE.

    !------------------------------------------------------Debugging options-----------------------------------------------------
    IF (PRESENT(Verbose)) Verbose_local = Verbose
    IF (PRESENT(Debug))   Debug_local   = Debug

    IF (Verbose_local > 27) WRITE(out_unit,*) 
    IF (Verbose_local > 27) WRITE(out_unit,*) "---------------------Constructing the matrix representation of the 1D HO Number of&
                                             & excitation Quanta operator---------------------"
  
    !---------------------------------------------Construction of the matrix Operator--------------------------------------------
    IF (.NOT. NbQuanta%Dense) THEN
      IF (Verbose_local > 28) WRITE(out_unit,*) "--- The Dense parameter of the NbQuanta operator is .FALSE., so the 1D HO NbQuan&
                                                &ta's matrix representation will be a rank-1 tensor of the diagonal elements of i&
                                                &ts analytical matrix (in Eigenbasis)"
      !---------------------------------------------Initialization to default values---------------------------------------------
      ALLOCATE(NbQuanta%Diag_val_R(HO1D_para%Nb))

      !------------------------------------------------Construction of the matrix------------------------------------------------
      DO i = 1, HO1D_para%Nb                                                                                                     ! /!\ Fortran counts from 1 to Nb !!! /!\
        NbQuanta%Diag_val_R(i) = i - 1
      END DO
  
      IF (Debug_local) WRITE(out_unit,*)
      IF (Debug_local) CALL Write_Vec(NbQuanta%Diag_val_R, out_unit, 1, info="HO1DNbQuanta")

    ELSE
      IF (Verbose_local > 28) WRITE(out_unit,*) "--- The Dense parameter of the NbQuanta is .TRUE., so the full 1D HO NbQuanta's &
                                                &matrix will be constructed (in Eigenbasis) for the representation, as if the ana&
                                                &lytical matrix was a dense one"
      !---------------------------------------------Initialization to default values---------------------------------------------
      ALLOCATE(NbQuanta%Dense_val_R(HO1D_para%Nb, HO1D_para%Nb))
      NbQuanta%Dense_val_R = ZERO

      !------------------------------------------------Construction of the matrix------------------------------------------------
      DO i = 1, HO1D_para%Nb                                                                            ! /!\ Fortran counts from 1 to Nb !!! /!\
        NbQuanta%Dense_val_R(i,i) = i - 1
      END DO

      IF (Debug_local) WRITE(out_unit,*)
      IF (Debug_local) CALL Write_Mat(NbQuanta%Dense_val_R, out_unit, Size(NbQuanta%Dense_val_R), info="HO1DNbQuanta")
    END IF
      
  END SUBROUTINE MolecCav_Initialize_N_HO1D

  
  FUNCTION MolecCav_Allocated_HO1D_operator(HO1D_operator) RESULT(Alloc)
    USE QDUtil_m
    IMPLICIT NONE 

    integer, allocatable           :: List_indexes(:)                                                                            ! the current values of the indexes for each dimension
    TYPE(ND_indexes_t), intent(in) :: ND_indexes

    integer                        :: Zero_index

    List_indexes = ND_indexes%Starting_indexes                                                                                   ! dynamic allocation : allows to call the function for any allocatable object, already allocated or not, and even for a tabular not declared allocatable it seems

    Zero_index = 1 + ND_indexes%Begin_right*(ND_indexes%N_dim-1)
    List_indexes(Zero_index) = List_indexes(Zero_index) - 1                                                                      ! to prepare the initial point I = 0 (state before first loop)

  END FUNCTION MolecCav_Allocated_HO1D_operator


  SUBROUTINE MolecCav_Write_HO1D_operator(Operator)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
    USE QDUtil_m
    USE HO1D_parameters_m
    IMPLICIT NONE 

    TYPE(HO1D_operator_t), intent(in) :: Operator

    IF (ALLOCATED(Operator%Operator_type)) THEN
      WRITE(out_unit,*) "_______________________________________the operator parameters______________________________________"
      WRITE(out_unit,*) "|The operator's nature <<Operator_type>> information do is allocated, and is   | ", Operator%Operator_type
      WRITE(out_unit,*) "|______________________________________________________________________________|____________________"
      FLUSH(out_unit)
    ELSE 
      WRITE(out_unit,*) "_____________________________the operator parameters____________________________"
      WRITE(out_unit,*) "|The operator's nature <<Operator_type>> information is NOT allocated.         |"
      WRITE(out_unit,*) "|______________________________________________________________________________|____________________"
      FLUSH(out_unit)
    END IF

    WRITE(out_unit,*) "|Is its matrix supposed to be represented as a dense one ? Operator%Dense      | "//TO_string(Operator%Dense)
    WRITE(out_unit,*) "|______________________________________________________________________________|____________________"
    
    WRITE(out_unit,*) "|Case of a band matrix, (Operator%Upper_bandwidth, Operator%Lower_bandwidth)   | ("//TO_string(Operator%Up&
                      &per_bandwidth)//TO_string(Operator%Lower_bandwidth)//")"
    WRITE(out_unit,*) "|______________________________________________________________________________|____________________"
    FLUSH(out_unit)

    WRITE(out_unit,*) "_________________________the operator's representations_________________________"
    IF (ALLOCATED(Operator%Diag_val_R)) THEN                                                                                     ! we assume that the code is supposed to be used only allocating one of the matrices of each Operator_t object
      WRITE(out_unit,*) "|The operator is represented using the H Eigenbasis with a basis size of       | ", Operator%Operator_type
      WRITE(out_unit,*) "|______________________________________________________________________________|____________________"
      FLUSH(out_unit)
      
      WRITE(out_unit,*) "|The operator's Diagonal matrix representation has been used, and is           |"
      WRITE(out_unit,*) "|______________________________________________________________________________|"
      CALL Write_Vec(Operator%Diag_val_R, out_unit, 1, info="Diag_val_R")
      WRITE(out_unit,*) "|______________________________________________________________________________|"
      FLUSH(out_unit)

    ELSE 
      WRITE(out_unit,*) "|The operator's Diagonal matrix representation is NOT allocated.               |"
      WRITE(out_unit,*) "|______________________________________________________________________________|"
      FLUSH(out_unit)
    END IF

    IF (ALLOCATED(Operator%Band_val_R)) THEN
      WRITE(out_unit,*) "|The operator is represented using the H Eigenbasis with a basis size of      | ", Operator%Operator_type
      WRITE(out_unit,*) "|_____________________________________________________________________________|_____________________"
      FLUSH(out_unit)

      WRITE(out_unit,*) "|The operator's Band matrix representation has been used, and is               |"
      WRITE(out_unit,*) "|______________________________________________________________________________|"
      CALL Write_Mat(Operator%Band_val_R, out_unit, 3, info="Band_val_R")
      WRITE(out_unit,*) "|______________________________________________________________________________|"
      FLUSH(out_unit)
  
    ELSE 
      WRITE(out_unit,*) "|The operator's Band matrix representation is NOT allocated.                   |"
      WRITE(out_unit,*) "|______________________________________________________________________________|"
      FLUSH(out_unit)
    END IF

    IF (ALLOCATED(Operator%Dense_val_R)) THEN
      WRITE(out_unit,*) "|The operator is represented using the H Eigenbasis with a basis size of      | ", Operator%Operator_type
      WRITE(out_unit,*) "|_____________________________________________________________________________|_____________________"
      FLUSH(out_unit)

      WRITE(out_unit,*) "|______________________________________________________________________________|"
      WRITE(out_unit,*) "|The operator's Dense matrix representation has been used, and is              |"
      CALL Write_Mat(Operator%Dense_val_R, out_unit, Size(Operator%Dense_val_R, dim=2), info="Dense_val_R")
      WRITE(out_unit,*) "|______________________________________________________________________________|"
      FLUSH(out_unit)
  
    ELSE 
      WRITE(out_unit,*) "|The operator's Dense matrix representation is NOT allocated.                  |"
      WRITE(out_unit,*) "|__________________________End HO1D operator object____________________________|"
      FLUSH(out_unit)
    END IF
  
  END SUBROUTINE MolecCav_Write_HO1D_operator

  SUBROUTINE MolecCav_Deallocate_HO1D_operator(Operator, Verbose, Debug)
    USE QDUtil_m
    IMPLICIT NONE 

    TYPE(HO1D_operator_t), intent(inout) :: Operator
    integer, optional,     intent(in)    :: Verbose                                                                                 ! cf. comments in HO1D_parameters_m
    logical, optional,     intent(in)    :: Debug                                                                                   ! cf. comments in HO1D_parameters_m

    integer                              :: Verbose_local = 25                                                                      ! goes from 25 (= 0 verbose) to 29 (= maximum verbose) at this layer
    logical                              :: Debug_local   = .FALSE.

    !------------------------------------------------------Debugging options-----------------------------------------------------
    IF (PRESENT(Verbose)) Verbose_local = Verbose
    IF (PRESENT(Debug))   Debug_local   = Debug

    !-----------------------------Deallocating the HO1D operator object----------------------------
    IF (Verbose_local > 27) WRITE(out_unit,*)
    IF (Verbose_local > 27) WRITE(out_unit,*) "-----------------------------------------------Deallocating the HO1D_operator obje&
                                              &ct----------------------------------------------"
  
    IF (Debug_local) THEN
      WRITE(out_unit,*) "--- The HO1D operator to be deallocated :"
      CALL Write_HO1D_operator(Operator)
      WRITE(out_unit,*) "--- End HO1D operator to be deallocated"
    END IF 

    IF (ALLOCATED(Operator%Operator_type)) DEALLOCATE(Operator%Operator_type)
    IF (ALLOCATED(Operator%Dense_val_R))   DEALLOCATE(Operator%Dense_val_R)
    IF (ALLOCATED(Operator%Diag_val_R))    DEALLOCATE(Operator%Diag_val_R)
    IF (ALLOCATED(Operator%Band_val_R))    DEALLOCATE(Operator%Band_val_R)

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "--- The HO1D operator object after having been deallocated :"
      CALL Write_HO1D_operator(Operator)
      WRITE(out_unit,*) "--- End dellocated HO1D operator"
    END IF

  END SUBROUTINE MolecCav_Deallocate_HO1D_operator


END MODULE
  