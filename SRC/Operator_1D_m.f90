MODULE Operator_1D_m
  !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
  USE QDUtil_m
  USE Cavity_mode_m
  IMPLICIT NONE


  TYPE, EXTENDS(Cavity_mode_t) :: Operator_1D_t
    character(len=:),    allocatable :: Operator_type                          ! ex : "Hamiltonian", "Position", "Nb_photons", etc
    !character(len=:),    allocatable :: matrix_shape_type                      ! ex : "Dense", "Diagonal", "Band", etc
    logical                          :: Dense = .FALSE.                        ! if .TRUE. then the matrix storage will not be optimized and it will be stored as a Dense matrix
    integer                          :: Upper_bandwidth = 0                    ! if type = "Band". Gives the number of additional bands to consider above the diagonal.
    integer                          :: Lower_bandwidth = 0                    ! if type = "Band". Gives the number of additional bands to consider below the diagonal. Ex : Upper_bandwidth=Lower_bandwidth=1 would give a tridiagonal matrix
    real(kind=Rkind),    allocatable :: Dense_val_R(:,:)                       ! if type = "Dense"
    real(kind=Rkind),    allocatable :: Diag_val_R(:)                          ! if type = "Diagonal"
    real(kind=Rkind),    allocatable :: Band_val_R(:,:)                        ! if subtype = "Band". The number of columns will be the number of diagonals to consider
  END TYPE


  PRIVATE

  PUBLIC Operator_1D_t, Construct_Operator, Action_Operator_1D, Average_value_operator_1D, & 
       & MolecCav_Action_cavity_operator_2D, MolecCav_Average_value_cavity_operator_2D

  INTERFACE Construct_Operator
    MODULE PROCEDURE MolecCav_Construct_Operator
  END INTERFACE
  INTERFACE Action_Operator_1D
    MODULE PROCEDURE MolecCav_Action_Operator_1D
  END INTERFACE
  INTERFACE Average_value_operator_1D
  MODULE PROCEDURE MolecCav_Average_value_operator_1D
  END INTERFACE
    

  CONTAINS


  SUBROUTINE MolecCav_Construct_Operator(Operator, Operator_type, Dense, Mode)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
    USE QDUtil_m
    IMPLICIT NONE

    TYPE(Operator_1D_t), intent(inout) :: Operator                             ! the object of type Operator_t to be constructed here
    character(len=*),    intent(in)    :: Operator_type                        ! ex : "Hamiltonian", "Position", etc
    logical, optional,   intent(in)    :: Dense                                ! if .TRUE. then the matrix storage will not be optimized and it will be stored as a Dense matrix
    TYPE(Cavity_mode_t), intent(in)    :: Mode                                 ! the HO/Cavity mode which the operator is relative to

    !--------------first steps of the construction of the Operator-------------
    Operator%Operator_type = Operator_type                                     ! allocation on assignement. Operator_type has the right lengths (no spaces added) thanks to len=* at declaration and it will fit the Op%op_type thanks to len=:, allocatable at declaration of the derived type. 
    Operator%D             = Mode%D
    Operator%Nb            = Mode%Nb
    Operator%w             = Mode%w
    Operator%m             = Mode%m
    Operator%lambda        = Mode%lambda

    IF (PRESENT(Dense)) THEN
      Operator%Dense = Dense
    END IF

    !--------------------construction of the matrix Operator-------------------
    SELECT CASE (TO_lowercase(Operator%Operator_type))                         ! TO_lowercase avoid case sensitivity issues
      CASE ("hamiltonian")
        CALL MolecCav_Construct_H_cavity_mode(Hamiltonian=Operator)
    
      CASE ("position")
        CALL MolecCav_Construct_x_cavity_mode(Position_Op=Operator)
      
      CASE ("nb_photons")
        CALL MolecCav_Construct_N_cavity_mode(Nb_photon_Op=Operator)

      CASE DEFAULT
        STOP "No Operator type recognized, please verify the input of Construct_Operator subroutine"

    END SELECT

  END SUBROUTINE


  SUBROUTINE MolecCav_Construct_H_cavity_mode(Hamiltonian)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    IMPLICIT NONE
    
    TYPE(Operator_1D_t), intent(inout) :: Hamiltonian                          ! matrix of the one-dimensional harmonic Hamiltonian associated with HO D

    integer                            :: i                                    ! loop increments /!\ Fortran counts from 1 to Nb !!! /!\

    WRITE(out_unit,*) ''
    WRITE(out_unit,*) '*******************************************************'
    WRITE(out_unit,*) '******* CONSTRUCTING THE HAMILTONIAN OF THE HO ********'

    IF (.NOT. Hamiltonian%Dense) THEN
      !---------------------Initialization to default values-------------------
      ALLOCATE(Hamiltonian%Diag_val_R(Hamiltonian%Nb))
      Hamiltonian%Diag_val_R = ZERO
      !------------------------Construction of the matrix----------------------
      DO i = 1, Hamiltonian%Nb                                                             ! /!\ Fortran counts from 1 to Nb !!! /!\
        Hamiltonian%Diag_val_R(i) = Hamiltonian%w*(i - ONE + HALF)                         ! "-1" because the first Fortran vector is the fundamental eigenvector of the HO i.e. the 0^{th} ket 
      END DO
      
    ELSE
      !---------------------Initialization to default values-------------------
      ALLOCATE(Hamiltonian%Dense_val_R(Hamiltonian%Nb, Hamiltonian%Nb))
      Hamiltonian%Dense_val_R = ZERO
      !------------------------Construction of the matrix----------------------
      DO i = 1, Hamiltonian%Nb                                                             ! /!\ Fortran counts from 1 to Nb !!! /!\
        Hamiltonian%Dense_val_R(i,i) = Hamiltonian%w*(i - ONE + HALF)                      ! "-1" because the first Fortran vector is the fundamental eigenvector of the HO i.e. the 0^{th} ket 
      END DO
    END IF
    
    WRITE(out_unit,*) '********* HAMILTONIAN OF THE HO CONSTRUCTED ***********'
    WRITE(out_unit,*) '*******************************************************'

  END SUBROUTINE


  SUBROUTINE MolecCav_Construct_x_cavity_mode(Position_Op)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    IMPLICIT NONE
    
    TYPE(Operator_1D_t), intent(inout) :: Position_Op

    integer                            :: i                                 ! loop increments /!\ Fortran counts from 1 to Nb !!! /!\

    WRITE(out_unit,*) ''
    WRITE(out_unit,*) '*******************************************************'
    WRITE(out_unit,*) '******* CONSTRUCTING THE POSITION OP OF THE HO ********'

    IF (.NOT. Position_Op%Dense) THEN
      !----------Initialization of the characteristics of the operator---------
      Position_Op%Upper_bandwidth   = 1
      Position_Op%Lower_bandwidth   = 1
      !---------------------Initialization to default values-------------------
      ALLOCATE(Position_Op%Band_val_R(Position_Op%Nb,3))                                   ! Nb lines (number of diagonal elements) and 3 columns because 3 bands to consider : the diagonal, and the two bands above and below it
      Position_Op%Band_val_R = ZERO
      !------------------------Construction of the matrix----------------------
      DO i = 1, Position_Op%Nb - 1                                                         ! /!\ Fortran counts from 1 to Nb !!! /!\ Nb-1 not to have Band_val_R(i+1) out of range
        Position_Op%Band_val_R(i,1)   = SQRT(REAL(i,kind=Rkind))
        Position_Op%Band_val_R(i+1,3) = SQRT(REAL(i,kind=Rkind))
      END DO
      Position_Op%Band_val_R = Position_Op%Band_val_R / SQRT(TWO * Position_Op%w * Position_Op%m)
        
    ELSE
      !---------------------Initialization to default values-------------------
      ALLOCATE(Position_Op%Dense_val_R(Position_Op%Nb, Position_Op%Nb))
      Position_Op%Dense_val_R = ZERO
      !------------------------Construction of the matrix----------------------
      DO i = 1, Position_Op%Nb - 1                                                         ! /!\ Fortran counts from 1 to Nb !!! /!\
        Position_Op%Dense_val_R(i,i+1) = SQRT(REAL(i,kind=Rkind))
        Position_Op%Dense_val_R(i+1,i) = SQRT(REAL(i,kind=Rkind))
      END DO
      Position_Op%Dense_val_R = Position_Op%Dense_val_R / SQRT(TWO * Position_Op%w * Position_Op%m)
    END IF
    
    WRITE(out_unit,*) '********* POSITION OP OF THE HO CONSTRUCTED ***********'
    WRITE(out_unit,*) '*******************************************************'

  END SUBROUTINE


  SUBROUTINE MolecCav_Construct_N_cavity_mode(Nb_photon_Op)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    IMPLICIT NONE
    
    TYPE(Operator_1D_t), intent(inout) :: Nb_photon_Op

    integer                            :: i                                 ! loop increments /!\ Fortran counts from 1 to Nb !!! /!\

    WRITE(out_unit,*) ''
    WRITE(out_unit,*) '*******************************************************'
    WRITE(out_unit,*) '******* CONSTRUCTING THE NB QUANTA OP OF THE HO *******'

    IF (.NOT. Nb_photon_Op%Dense) THEN
      !---------------------Initialization to default values-------------------
      ALLOCATE(Nb_photon_Op%Diag_val_R(Nb_photon_Op%Nb))
      Nb_photon_Op%Diag_val_R = ZERO
      !------------------------Construction of the matrix----------------------
      DO i = 1, Nb_photon_Op%Nb                                                             ! /!\ Fortran counts from 1 to Nb !!! /!\
        Nb_photon_Op%Diag_val_R(i) = i - 1
      END DO
      
    ELSE
      !---------------------Initialization to default values-------------------
      ALLOCATE(Nb_photon_Op%Dense_val_R(Nb_photon_Op%Nb, Nb_photon_Op%Nb))
      Nb_photon_Op%Dense_val_R = ZERO
      !------------------------Construction of the matrix----------------------
      DO i = 1, Nb_photon_Op%Nb                                                             ! /!\ Fortran counts from 1 to Nb !!! /!\
        Nb_photon_Op%Dense_val_R(i,i) = i - 1
      END DO
    END IF
    
    WRITE(out_unit,*) '********* NB QUANTA OP OF THE HO CONSTRUCTED **********'
    WRITE(out_unit,*) '*******************************************************'

  END SUBROUTINE

  
  SUBROUTINE MolecCav_Action_Operator_1D(Op_psi, Operator, Psi)   ! /!\ FOR NOW EVERYTHING IS REAL /!\ compute the resulting vector Psi_result(:) from the action of the operator of the cavity mode on the photon state vector Psi_argument(:) written in the Eigenbasis of H_ho
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    IMPLICIT NONE

    real(kind=Rkind),    intent(inout) :: Op_psi(:)
    TYPE(Operator_1D_t), intent(in)    :: Operator
    real(kind=Rkind),    intent(in)    :: Psi(:)

    !--------------------Selection of the calculation method-------------------
    IF      (ALLOCATED(Operator%Diag_val_R))   THEN
      CALL MolecCav_Action_Diag_Operator_1D(Op_psi=Op_psi, Operator=Operator, Psi=Psi)

    ELSE IF (ALLOCATED(Operator%Band_val_R))   THEN
      CALL MolecCav_Action_Band_Operator_1D(Op_psi=Op_psi, Operator=Operator, Psi=Psi)

    ELSE IF (ALLOCATED(Operator%Dense_val_R)) THEN
      CALL MolecCav_Action_Dense_Operator_1D(Op_psi=Op_psi, Operator=Operator, Psi=Psi)

    ELSE
      STOP "None of this operator's matrices are allocated. Please check its initalization."
    
    END IF

  END SUBROUTINE

  
  SUBROUTINE MolecCav_Action_Dense_Operator_1D(Op_psi, Operator, Psi)  ! _R for the case where Psi is a real vector. 
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    IMPLICIT NONE
    
    real(kind=Rkind),    intent(inout) :: Op_psi(:)
    TYPE(Operator_1D_t), intent(in)    :: Operator
    real(kind=Rkind),    intent(in)    :: Psi(:)

    Op_psi = ZERO
    Op_psi(:) = matmul(Operator%Dense_val_R, Psi)

  END SUBROUTINE

  
  SUBROUTINE MolecCav_Action_Diag_Operator_1D(Op_psi, Operator, Psi) 
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    IMPLICIT NONE
    
    real(kind=Rkind),    intent(inout) :: Op_psi(:)
    TYPE(Operator_1D_t), intent(in)    :: Operator
    real(kind=Rkind),    intent(in)    :: Psi(:)

    Op_psi = ZERO
    Op_psi = Operator%Diag_val_R * Psi

  END SUBROUTINE

  
  SUBROUTINE MolecCav_Action_Band_Operator_1D(Op_psi, Operator, Psi)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    IMPLICIT NONE
    
    real(kind=Rkind),    intent(inout) :: Op_psi(:)
    TYPE(Operator_1D_t), intent(in)    :: Operator
    real(kind=Rkind),    intent(in)    :: Psi(:)

    integer                            :: i, Nb

    Nb = size(Op_psi)
    IF (Nb /= Operator%Nb) THEN
      STOP "The size of the operator's matrix Band_val_R does not match the size of the &
          & resulting vector Op_psi. Please check their initialization."
    ELSE IF (Nb /= Size(Psi)) THEN
      STOP "The size of the resulting vector Op_psi does not match the size of the operand & 
          & vector Psi. Please check their initialization."
    END IF

    Op_psi     = ZERO
    Op_psi     = Operator%Band_val_R(:,2) * Psi
    Op_psi(1)  = Op_psi(1)  + Operator%Band_val_R(2,3)    * Psi(2)
    Op_psi(Nb) = Op_psi(Nb) + Operator%Band_val_R(Nb-1,1) * Psi(Nb-1)
    DO i = 2, Nb-1
      Op_psi(i) = Op_psi(i) + &
                & Operator%Band_val_R(i-1,1) * Psi(i-1) + &
                & Operator%Band_val_R(i+1,3) * Psi(i+1)
    END DO

  END SUBROUTINE
  

  SUBROUTINE MolecCav_Average_value_operator_1D(Value, Operator, Psi)   ! /!\ FOR NOW EVERYTHING IS REAL /!\ compute the resulting vector Psi_result(:) from the action of the operator of the cavity mode on the photon state vector Psi_argument(:) written in the Eigenbasis of H_ho
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    IMPLICIT NONE

    real(kind=Rkind),    intent(inout) :: Value
    TYPE(Operator_1D_t), intent(in)    :: Operator
    real(kind=Rkind),    intent(in)    :: Psi(:)

    real(kind=Rkind), allocatable      :: Intermediary(:)
    integer                            :: Nb

    Nb = Size(Psi)
    ALLOCATE(Intermediary(Nb))

    CALL MolecCav_Action_Operator_1D(Intermediary, Operator, Psi)
    Value = DOT_PRODUCT(Psi, Intermediary) 

    DEALLOCATE(Intermediary)
    
  END SUBROUTINE


  SUBROUTINE MolecCav_Action_cavity_operator_2D(Op_psi, Operator, Psi)   ! /!\ FOR NOW EVERYTHING IS REAL /!\ compute the resulting vector Psi_result(:) from the action of the operator of the cavity mode on the photon state vector Psi_argument(:) written in the Eigenbasis of H_ho
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    IMPLICIT NONE

    real(kind=Rkind),    intent(inout) :: Op_psi(:,:)                      ! already allocated !
    TYPE(Operator_1D_t), intent(in)    :: Operator
    real(kind=Rkind),    intent(in)    :: Psi(:,:)

    integer                            :: Nb_M, Nb_C, i_M

    Nb_M = Size(Psi, 1)
    Nb_C = Size(Psi, 2)
    
    IF (Nb_M /= Size(Op_psi, 1) .OR. Nb_C /= Size(Op_psi, 2)) THEN
      STOP "The size of the operand Psi's vector does not match the size of the resulting vector Op_psi."
    END IF

    DO i_M = 1, Nb_M
      CALL MolecCav_Action_Operator_1D(Op_psi(i_M, :), Operator, Psi(i_M, :))
    END DO

  END SUBROUTINE


  SUBROUTINE MolecCav_Average_value_cavity_operator_2D(Value, Operator, Psi)   ! /!\ FOR NOW EVERYTHING IS REAL /!\ compute the resulting vector Psi_result(:) from the action of the operator of the cavity mode on the photon state vector Psi_argument(:) written in the Eigenbasis of H_ho
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    USE Algebra_m
    IMPLICIT NONE

    real(kind=Rkind),    intent(inout) :: Value
    TYPE(Operator_1D_t), intent(in)    :: Operator
    real(kind=Rkind),    intent(in)    :: Psi(:,:)

    real(kind=Rkind), allocatable      :: Intermediary(:,:)
    integer                            :: Nb_M, Nb_C

    Nb_M = Size(Psi, 1)
    Nb_C = Size(Psi, 2)
    ALLOCATE(Intermediary(Nb_M, Nb_C))

    CALL MolecCav_Action_cavity_operator_2D(Intermediary, Operator, Psi)
    CALL Scalar_product(Value, Intermediary, Psi)

  END SUBROUTINE
  

END MODULE
