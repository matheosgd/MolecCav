MODULE MC_operator_1D_m
  !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
  USE QDUtil_m
  IMPLICIT NONE


  TYPE :: MC_operator_1D_t
    character(len=:),    allocatable :: operator_type                          ! ex : "Hamiltonian", "Position", "Nb_photons", etc
    character(len=:),    allocatable :: scalar_space                           ! ex : "Real" or "Complex"
    character(len=:),    allocatable :: matrix_shape_type                      ! ex : "Dense", "Diagonal", "Band", etc
    integer                          :: Upper_bandwidth = 0                    ! if type = "Band". Gives the number of additional bands to consider above the diagonal.
    integer                          :: Lower_bandwidth = 0                    ! if type = "Band". Gives the number of additional bands to consider below the diagonal. Ex : Upper_bandwidth=Lower_bandwidth=1 would give a tridiagonal matrix
    real(kind=Rkind),    allocatable :: Dense_val_R(:,:)                       ! if type = "Dense"
    complex(kind=Rkind), allocatable :: Dense_val_C(:,:)                       ! if type = "Dense"
    real(kind=Rkind),    allocatable :: Diag_val_R(:)                          ! if type = "Diagonal"
    complex(kind=Rkind), allocatable :: Diag_val_C(:)                          ! if type = "Diagonal"
    real(kind=Rkind),    allocatable :: Band_val_R(:,:)                        ! if subtype = "Band". The number of columns will be the number of diagonals to consider
    complex(kind=Rkind), allocatable :: Band_val_C(:,:)                        ! if subtype = "Band"
  END TYPE


  CONTAINS


  SUBROUTINE MolecCav_Construct_Operator(Operator, operator_type, scalar_space, matrix_shape_type, Nb, w, m)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
    USE QDUtil_m
    IMPLICIT NONE

    TYPE(MC_operator_1D_t), intent(inout) :: Operator                          ! the object of type Operator_t to be constructed here
    character(len=*),       intent(in)    :: operator_type                     ! ex : "Hamiltonian", "Position", etc
    character(len=*),       intent(in)    :: scalar_space                      ! ex : "Real" or "Complex"
    character(len=*),       intent(in)    :: matrix_shape_type                 ! the assignement to use the analytical shape of the operator's matrix, or the dense shape
    integer,                intent(in)    :: Nb                                ! number of basis vectors associated with the HO of whom Operator is the operator
    real(kind=Rkind),       intent(in)    :: w, m                              ! eigenpulsation and mass associated with this HO

    !--------------first steps of the construction of the Operator-------------
    Operator%operator_type = operator_type                                     ! allocation on assignement. operator_type has the right lengths (no spaces added) thanks to len=* at declaration and it will fit the Op%op_type thanks to len=:, allocatable at declaration of the derived type. 
    Operator%scalar_space  = scalar_space                                      ! allocation on assignement too.

    !--------------------construction of the matrix Operator-------------------
    SELECT CASE (TO_lowercase(Operator%operator_type))                         ! TO_lowercase avoid case sensitivity issues
      CASE ("hamiltonian")
        CALL MolecCav_Construct_H_cavity_mode(Hamiltonian=Operator, matrix_shape_type=matrix_shape_type, Nb=Nb, w=w)
      
      CASE ("position")
        CALL MolecCav_Construct_x_cavity_mode(Position_Op=Operator, matrix_shape_type=matrix_shape_type, Nb=Nb, w=w, m=m)
      
      CASE ("nb_photons")
        CALL MolecCav_Construct_N_cavity_mode(Nb_photon_Op=Operator, matrix_shape_type=matrix_shape_type, Nb=Nb)

      CASE DEFAULT
        STOP "No Operator type recognized, please verify the input of Construct_Operator subroutine"

    END SELECT

  END SUBROUTINE


  SUBROUTINE MolecCav_Construct_H_cavity_mode(Hamiltonian, matrix_shape_type, Nb, w)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    IMPLICIT NONE
    
    TYPE(MC_operator_1D_t), intent(inout) :: Hamiltonian                       ! matrix of the one-dimensional harmonic Hamiltonian associated with HO D
    character(len=*),       intent(in)    :: matrix_shape_type                 ! the assignement to use the analytical shape of the operator's matrix, or the dense shape
    integer,                intent(in)    :: Nb                                ! number of basis vectors associated with the HO of whom Operator is the operator
    real(kind=Rkind),       intent(in)    :: w                                 ! eigenpulsation associated with this HO

    integer                               :: i                                 ! loop increments /!\ Fortran counts from 1 to Nb !!! /!\

    WRITE(out_unit,*) ''
    WRITE(out_unit,*) '*******************************************************'
    WRITE(out_unit,*) '******* CONSTRUCTING THE HAMILTONIAN OF THE HO ********'

    IF (TO_lowercase(matrix_shape_type) == "opt") THEN
      !----------Initialization of the characteristics of the operator---------
      Hamiltonian%matrix_shape_type = "diagonal"
      !---------------------Initialization to default values-------------------
      ALLOCATE(Hamiltonian%Diag_val_R(Nb))
      Hamiltonian%Diag_val_R = ZERO
      !------------------------Construction of the matrix----------------------
      DO i = 1, Nb                                                             ! /!\ Fortran counts from 1 to Nb !!! /!\
        Hamiltonian%Diag_val_R(i) = w*(i - ONE + HALF)                         ! "-1" because the first Fortran vector is the fundamental eigenvector of the HO i.e. the 0^{th} ket 
      END DO
      
    ELSE IF (TO_lowercase(matrix_shape_type) == "non_opt") THEN
      !----------Initialization of the characteristics of the operator---------
      Hamiltonian%matrix_shape_type = "dense"
      !---------------------Initialization to default values-------------------
      ALLOCATE(Hamiltonian%Dense_val_R(Nb, Nb))
      Hamiltonian%Dense_val_R = ZERO
      !------------------------Construction of the matrix----------------------
      DO i = 1, Nb                                                             ! /!\ Fortran counts from 1 to Nb !!! /!\
        Hamiltonian%Dense_val_R(i,i) = w*(i - ONE + HALF)                      ! "-1" because the first Fortran vector is the fundamental eigenvector of the HO i.e. the 0^{th} ket 
      END DO
    END IF
    
    WRITE(out_unit,*) '********* HAMILTONIAN OF THE HO CONSTRUCTED ***********'
    WRITE(out_unit,*) '*******************************************************'

  END SUBROUTINE


  SUBROUTINE MolecCav_Construct_x_cavity_mode(Position_Op, matrix_shape_type, Nb, w, m)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    IMPLICIT NONE
    
    TYPE(MC_operator_1D_t), intent(inout) :: Position_Op
    character(len=*),       intent(in)    :: matrix_shape_type                 ! the assignement to use the analytical shape of the operator's matrix, or the dense shape
    integer,                intent(in)    :: Nb                                ! number of basis vectors associated with the HO of whom Operator is the operator
    real(kind=Rkind),       intent(in)    :: w, m                              ! eigenpulsation and mass associated with this HO (m=1 for a cavity HO)

    integer                               :: i                                 ! loop increments /!\ Fortran counts from 1 to Nb !!! /!\

    WRITE(out_unit,*) ''
    WRITE(out_unit,*) '*******************************************************'
    WRITE(out_unit,*) '******* CONSTRUCTING THE POSITION OP OF THE HO ********'

    IF (TO_lowercase(matrix_shape_type) == "opt") THEN
      !----------Initialization of the characteristics of the operator---------
      Position_Op%matrix_shape_type = "band"
      Position_Op%Upper_bandwidth   = 1
      Position_Op%Lower_bandwidth   = 1
      !---------------------Initialization to default values-------------------
      ALLOCATE(Position_Op%Band_val_R(Nb,3))                                   ! Nb lines (number of diagonal elements) and 3 columns because 3 bands to consider : the diagonal, and the two bands above and below it
      Position_Op%Band_val_R = ZERO
      !------------------------Construction of the matrix----------------------
      DO i = 1, Nb - 1                                                         ! /!\ Fortran counts from 1 to Nb !!! /!\ Nb-1 not to have Band_val_R(i+1) out of range
        Position_Op%Band_val_R(i,1)   = SQRT(REAL(i,kind=Rkind))
        Position_Op%Band_val_R(i+1,3) = SQRT(REAL(i,kind=Rkind))
      END DO
      Position_Op%Band_val_R = Position_Op%Band_val_R / SQRT(TWO * w * m)
        
    ELSE IF (TO_lowercase(matrix_shape_type) == "non_opt") THEN
      !----------Initialization of the characteristics of the operator---------
      Position_Op%matrix_shape_type = "dense"
      !---------------------Initialization to default values-------------------
      ALLOCATE(Position_Op%Dense_val_R(Nb, Nb))
      Position_Op%Dense_val_R = ZERO
      !------------------------Construction of the matrix----------------------
      DO i = 1, Nb - 1                                                         ! /!\ Fortran counts from 1 to Nb !!! /!\
        Position_Op%Dense_val_R(i,i+1) = SQRT(REAL(i,kind=Rkind))
        Position_Op%Dense_val_R(i+1,i) = SQRT(REAL(i,kind=Rkind))
      END DO
      Position_Op%Dense_val_R = Position_Op%Dense_val_R / SQRT(TWO * w * m)
    END IF
    
    WRITE(out_unit,*) '********* POSITION OP OF THE HO CONSTRUCTED ***********'
    WRITE(out_unit,*) '*******************************************************'

  END SUBROUTINE


  SUBROUTINE MolecCav_Construct_N_cavity_mode(Nb_photon_Op, matrix_shape_type, Nb)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    IMPLICIT NONE
    
    TYPE(MC_operator_1D_t), intent(inout) :: Nb_photon_Op
    character(len=*),       intent(in)    :: matrix_shape_type                 ! the assignement to use the analytical shape of the operator's matrix, or the dense shape
    integer,                intent(in)    :: Nb                                ! number of basis vectors associated with the HO of whom Operator is the operator

    integer                               :: i                                 ! loop increments /!\ Fortran counts from 1 to Nb !!! /!\

    WRITE(out_unit,*) ''
    WRITE(out_unit,*) '*******************************************************'
    WRITE(out_unit,*) '******* CONSTRUCTING THE NB QUANTA OP OF THE HO *******'

    IF (TO_lowercase(matrix_shape_type) == "opt") THEN
      !----------Initialization of the characteristics of the operator---------
      Nb_photon_Op%matrix_shape_type = "diagonal"
      !---------------------Initialization to default values-------------------
      ALLOCATE(Nb_photon_Op%Diag_val_R(Nb))
      Nb_photon_Op%Diag_val_R = ZERO
      !------------------------Construction of the matrix----------------------
      DO i = 1, Nb                                                             ! /!\ Fortran counts from 1 to Nb !!! /!\
        Nb_photon_Op%Diag_val_R(i) = i - 1
      END DO
      
    ELSE IF (TO_lowercase(matrix_shape_type) == "non_opt") THEN
      !----------Initialization of the characteristics of the operator---------
      Nb_photon_Op%matrix_shape_type = "dense"
      !---------------------Initialization to default values-------------------
      ALLOCATE(Nb_photon_Op%Dense_val_R(Nb, Nb))
      Nb_photon_Op%Dense_val_R = ZERO
      !------------------------Construction of the matrix----------------------
      DO i = 1, Nb                                                             ! /!\ Fortran counts from 1 to Nb !!! /!\
        Nb_photon_Op%Dense_val_R(i,i) = i - 1
      END DO
    END IF
    
    WRITE(out_unit,*) '********* NB QUANTA OP OF THE HO CONSTRUCTED **********'
    WRITE(out_unit,*) '*******************************************************'

  END SUBROUTINE

  
  SUBROUTINE MolecCav_Action_Operator_1D(Psi_result, Operator, Psi_argument)   ! /!\ FOR NOW EVERYTHING IS REAL /!\ compute the resulting vector Psi_result(:) from the action of the operator of the cavity mode on the photon state vector Psi_argument(:) written in the Eigenbasis of H_ho
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    IMPLICIT NONE

    real(kind=Rkind), intent(inout)    :: Psi_result(:)
    TYPE(MC_operator_1D_t), intent(in) :: Operator
    real(kind=Rkind), intent(in)       :: Psi_argument(:)

    !--------------------Selection of the calculation method-------------------
    SELECT CASE (TO_lowercase(Operator%matrix_shape_type))
      CASE ("dense")
        CALL MolecCav_Action_Dense_Operator_1D(Psi_result=Psi_result, Operator=Operator, Psi_argument=Psi_argument)
      
      CASE ("diagonal")
        CALL MolecCav_Action_Diag_Operator_1D(Psi_result=Psi_result, Operator=Operator, Psi_argument=Psi_argument)
      
      CASE ("band")
        CALL MolecCav_Action_Band_Operator_1D(Psi_result=Psi_result, Operator=Operator, Psi_argument=Psi_argument)

      CASE DEFAULT
        STOP "No Operator type recognized, please verify the input of Construct_Operator subroutine"

    END SELECT

  END SUBROUTINE

  
  SUBROUTINE MolecCav_Action_Dense_Operator_1D(Psi_result, Operator, Psi_argument)  ! _R for the case where Psi is a real vector. 
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    IMPLICIT NONE
    
    real(kind=Rkind), intent(inout)    :: Psi_result(:)
    TYPE(MC_operator_1D_t), intent(in) :: Operator
    real(kind=Rkind), intent(in)       :: Psi_argument(:)

    Psi_result = ZERO
    Psi_result(:) = matmul(Operator%Dense_val_R, Psi_argument)

  END SUBROUTINE

  
  SUBROUTINE MolecCav_Action_Diag_Operator_1D(Psi_result, Operator, Psi_argument) 
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    IMPLICIT NONE
    
    real(kind=Rkind), intent(inout)    :: Psi_result(:)
    TYPE(MC_operator_1D_t), intent(in) :: Operator
    real(kind=Rkind), intent(in)       :: Psi_argument(:)

    Psi_result = ZERO
    Psi_result = Operator%Diag_val_R * Psi_argument

  END SUBROUTINE

  
  SUBROUTINE MolecCav_Action_Band_Operator_1D(Psi_result, Operator, Psi_argument)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    IMPLICIT NONE
    
    real(kind=Rkind), intent(inout)    :: Psi_result(:)
    TYPE(MC_operator_1D_t), intent(in) :: Operator
    real(kind=Rkind), intent(in)       :: Psi_argument(:)

    integer                            :: i, Nb

    Nb = size(Psi_result)

    Psi_result     = ZERO
    Psi_result     = Operator%Band_val_R(:,2) * Psi_argument
    Psi_result(1)  = Psi_result(1)  + Operator%Band_val_R(2,3)    * Psi_argument(2)
    Psi_result(Nb) = Psi_result(Nb) + Operator%Band_val_R(Nb-1,1) * Psi_argument(Nb-1)
    DO i = 2, Nb-1
      Psi_result(i) = Psi_result(i) + &
                    & Operator%Band_val_R(i-1,1) * Psi_argument(i-1) + &
                    & Operator%Band_val_R(i+1,3) * Psi_argument(i+1)
    END DO

  END SUBROUTINE
  

  SUBROUTINE MolecCav_Average_value_operator_1D(Value, Operator, Psi_argument)   ! /!\ FOR NOW EVERYTHING IS REAL /!\ compute the resulting vector Psi_result(:) from the action of the operator of the cavity mode on the photon state vector Psi_argument(:) written in the Eigenbasis of H_ho
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    IMPLICIT NONE

    real(kind=Rkind), intent(inout)    :: Value
    TYPE(MC_operator_1D_t), intent(in) :: Operator
    real(kind=Rkind), intent(in)       :: Psi_argument(:)

    real(kind=Rkind), allocatable      :: Intermediary(:)
    integer                            :: Nb

    Nb = Size(Psi_argument)
    ALLOCATE(Intermediary(Nb))

    CALL MolecCav_Action_Operator_1D(Intermediary, Operator, Psi_argument)
    Value = DOT_PRODUCT(Psi_argument, Intermediary) 

    DEALLOCATE(Intermediary)
    

  END SUBROUTINE


  SUBROUTINE MolecCav_Action_cavity_operator_2D(Psi_result, Operator, Psi_argument)   ! /!\ FOR NOW EVERYTHING IS REAL /!\ compute the resulting vector Psi_result(:) from the action of the operator of the cavity mode on the photon state vector Psi_argument(:) written in the Eigenbasis of H_ho
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    IMPLICIT NONE

    real(kind=Rkind), intent(inout)    :: Psi_result(:,:)                      ! already allocated !
    TYPE(MC_operator_1D_t), intent(in) :: Operator
    real(kind=Rkind), intent(in)       :: Psi_argument(:,:)

    integer                            :: Nb_M, Nb_C, i_M

    Nb_M = Size(Psi_argument, 1)
    Nb_C = Size(Psi_argument, 2)
    
    IF (Nb_M /= Size(Psi_result, 1) .OR. Nb_C /= Size(Psi_result, 2)) THEN
      STOP "The size of the Psi_argument does not match the size of the Psi_result."
    END IF

    DO i_M = 1, Nb_M
      CALL MolecCav_Action_Operator_1D(Psi_result(i_M, :), Operator, Psi_argument(i_M, :))
    END DO

  END SUBROUTINE


  SUBROUTINE MolecCav_Average_value_cavity_operator_2D(Value, Operator, Psi_argument)   ! /!\ FOR NOW EVERYTHING IS REAL /!\ compute the resulting vector Psi_result(:) from the action of the operator of the cavity mode on the photon state vector Psi_argument(:) written in the Eigenbasis of H_ho
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    USE MC_algebra_m
    IMPLICIT NONE

    real(kind=Rkind), intent(inout)    :: Value
    TYPE(MC_operator_1D_t), intent(in) :: Operator
    real(kind=Rkind), intent(in)       :: Psi_argument(:,:)

    real(kind=Rkind), allocatable      :: Intermediary(:,:)
    integer                            :: Nb_M, Nb_C

    Nb_M = Size(Psi_argument, 1)
    Nb_C = Size(Psi_argument, 2)
    ALLOCATE(Intermediary(Nb_M, Nb_C))

    CALL MolecCav_Action_cavity_operator_2D(Intermediary, Operator, Psi_argument)
    CALL MolecCav_scalar_product_2D(Value, Intermediary, Psi_argument)

  END SUBROUTINE
  

END MODULE
