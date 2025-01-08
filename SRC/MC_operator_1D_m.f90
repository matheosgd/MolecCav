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
    SELECT CASE (Operator%operator_type)
      CASE ("Hamiltonian")
        CALL MolecCav_Construct_H_cavity_mode(Hamiltonian=Operator, matrix_shape_type=matrix_shape_type, Nb=Nb, w=w)
      
      CASE ("Position")
        CALL MolecCav_Construct_x_cavity_mode(Position_Op=Operator, matrix_shape_type=matrix_shape_type, Nb=Nb, w=w, m=m)
      
      CASE ("Nb_photons")
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

    IF (matrix_shape_type == "Opt") THEN
      !----------Initialization of the characteristics of the operator---------
      Hamiltonian%matrix_shape_type = "Diagonal"
      !---------------------Initialization to default values-------------------
      ALLOCATE(Hamiltonian%Diag_val_R(Nb))
      Hamiltonian%Diag_val_R = ZERO
      !------------------------Construction of the matrix----------------------
      DO i = 1, Nb                                                             ! /!\ Fortran counts from 1 to Nb !!! /!\
        Hamiltonian%Diag_val_R(i) = w*(i - ONE + HALF)                         ! "-1" because the first Fortran vector is the fundamental eigenvector of the HO i.e. the 0^{th} ket 
      END DO
      
    ELSE IF (matrix_shape_type == "Non_opt") THEN
      !----------Initialization of the characteristics of the operator---------
      Hamiltonian%matrix_shape_type = "Dense"
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

    IF (matrix_shape_type == "Opt") THEN
      !----------Initialization of the characteristics of the operator---------
      Position_Op%matrix_shape_type = "Band"
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
        
    ELSE IF (matrix_shape_type == "Non_opt") THEN
      !----------Initialization of the characteristics of the operator---------
      Position_Op%matrix_shape_type = "Dense"
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

    IF (matrix_shape_type == "Opt") THEN
      !----------Initialization of the characteristics of the operator---------
      Nb_photon_Op%matrix_shape_type = "Diagonal"
      !---------------------Initialization to default values-------------------
      ALLOCATE(Nb_photon_Op%Diag_val_R(Nb))
      Nb_photon_Op%Diag_val_R = ZERO
      !------------------------Construction of the matrix----------------------
      DO i = 1, Nb                                                             ! /!\ Fortran counts from 1 to Nb !!! /!\
        Nb_photon_Op%Diag_val_R(i) = i - 1
      END DO
      
    ELSE IF (matrix_shape_type == "Non_opt") THEN
      !----------Initialization of the characteristics of the operator---------
      Nb_photon_Op%matrix_shape_type = "Dense"
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
    SELECT CASE (Operator%matrix_shape_type)
      CASE ("Dense")
        CALL MolecCav_Action_Dense_Operator_1D(Psi_result=Psi_result, Operator=Operator, Psi_argument=Psi_argument)
      
      CASE ("Diagonal")
        CALL MolecCav_Action_Diag_Operator_1D(Psi_result=Psi_result, Operator=Operator, Psi_argument=Psi_argument)
      
      CASE ("Band")
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

  SUBROUTINE MolecCav_Action_Total_Hamiltonian_1D(Result_total_WF, Matter_hamiltonianSystem_WF, &
                                                 & Cavity_hamiltonian_1D, Cavity_position_1D, &
                                                 & Matter_dipolar_moment, System_WF, &
                                                 & lambda_cavity_mode, w_cavity_mode)
    USE QDUtil_m
    IMPLICIT NONE

    real(kind=Rkind), intent(inout)    :: Result_total_WF(:,:)                 ! already allocated !
    real(kind=Rkind), intent(in)       :: Matter_hamiltonianSystem_WF(:,:)     ! size Nb_M*Nb_C. |H_MatterSystem_WF(:,i_C)> = H_Matter|System_WF(:,i_C)>
    TYPE(MC_operator_1D_t), intent(in) :: Cavity_hamiltonian_1D
    TYPE(MC_operator_1D_t), intent(in) :: Cavity_position_1D
    TYPE(MC_operator_1D_t), intent(in) :: Matter_dipolar_moment
    real(kind=Rkind), intent(in)       :: System_WF(:,:)                       ! size Nb_M*Nb_C
    real(kind=Rkind), intent(in)       :: lambda_cavity_mode, w_cavity_mode    ! coupling strenght and eigenpulsation

    real(kind=Rkind), allocatable      :: Cavity_hamiltonian_1DSystem_WF(:)    ! only one line of System_WF
    real(kind=Rkind), allocatable      :: Intermediary(:,:)
    real(kind=Rkind), allocatable      :: Matter_cavity_coupling_hamiltonian_1DSystem_WF(:,:) ! only one line of System_WF
    integer                            :: Nb_M, Nb_C, i_M, i_C

    Nb_M = Size(System_WF, 1)
    Nb_C = Size(System_WF, 2)

    Result_total_WF = ZERO

    !-----------H_tot = H_Matter + H_Cavity + H_MatterCavityCoupling-----------
    !----------H_Matter|System_WF> = H_Matter(System_WF) already known---------
    Result_total_WF = Matter_hamiltonianSystem_WF

    !-----------H_Cavity|System_WF> = H_MatterSystem_WF already known----------
    ALLOCATE(Cavity_hamiltonian_1DSystem_WF(Nb_C))
    DO i_M = 1, Nb_M
      CALL MolecCav_Action_Operator_1D(Cavity_hamiltonian_1DSystem_WF, &
                                     & Cavity_hamiltonian_1D, &
                                     & System_WF(i_M,:))
      Result_total_WF(i_M,:) = Result_total_WF(i_M,:) + Cavity_hamiltonian_1DSystem_WF(:)
    END DO

    !---------------------H_MatterCavityCoupling|System_WF>--------------------
    !------------------------H_CavityCoupling(System_WF)-----------------------
    ALLOCATE(Intermediary(Nb_M,Nb_C))
    DO i_M = 1, Nb_M
      CALL MolecCav_Action_Operator_1D(Intermediary(i_M,:), &
                                     & Cavity_position_1D, &
                                     & System_WF(i_M,:))
    END DO
    Intermediary = lambda_cavity_mode*w_cavity_mode*Intermediary

    !----------------------H_MatterCoupling(Intermediary)----------------------
    ALLOCATE(Matter_cavity_coupling_hamiltonian_1DSystem_WF(Nb_M,Nb_C))
    DO i_C = 1, Nb_C
      CALL MolecCav_Action_Operator_1D(Matter_cavity_coupling_hamiltonian_1DSystem_WF(:,i_C), &
                                     & Matter_dipolar_moment, &
                                     & Intermediary(:,i_C))
    END DO

    !-----------------------------H_tot = summation----------------------------
    Result_total_WF = Result_total_WF + Matter_cavity_coupling_hamiltonian_1DSystem_WF
  
  END SUBROUTINE


END MODULE
