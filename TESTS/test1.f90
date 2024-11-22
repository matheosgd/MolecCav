PROGRAM test1
  USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
  IMPLICIT NONE


  TYPE :: MC_cavity_mode_t                                                     ! MC = MolecCav
    integer                        :: D      = 1                               ! label of the HO/mode/dimension/associated basis set
    integer                        :: Nb     = 1                               ! number of basis vectors associated with the HO D
    real(kind=real64)              :: w      = 1.0_real64                      ! eigenpulsation associated with the HO D
    real(kind=real64)              :: m      = 1.0_real64                      ! mass associated with the HO D
    real(kind=real64)              :: lambda = 1.0_real64                      ! Strength parameter of the coupling between the mode D and the molecule
  END TYPE


  TYPE :: MC_operator_1D_t
    character(len=:),     allocatable :: operator_type                         ! ex : "Hamiltonian", "Position", "Nb_photons", etc
    character(len=:),     allocatable :: scalar_space                          ! ex : "Real" or "Complex"
    !### NB : The following types concern mainly the two-dimension matrices ###
    character(len=:),     allocatable :: matrix_shape_type                     ! ex : "Dense", "Diagonal", "Band", etc
    integer                           :: Upper_bandwidth = 0                   ! if type = "Band". Gives the number of additional bands to consider above the diagonal.
    integer                           :: Lower_bandwidth = 0                   ! if type = "Band". Gives the number of additional bands to consider below the diagonal. Ex : Upper_bandwidth=Lower_bandwidth=1 would give a tridiagonal matrix
    real(kind=real64),    allocatable :: Dense_val_R(:,:)                      ! if type = "Dense"
    complex(kind=real64), allocatable :: Dense_val_C(:,:)                      ! if type = "Dense"
    real(kind=real64),    allocatable :: Diag_val_R(:)                         ! if type = "Diagonal"
    complex(kind=real64), allocatable :: Diag_val_C(:)                         ! if type = "Diagonal"
    real(kind=real64),    allocatable :: Band_val_R(:,:)                       ! if subtype = "Band". The number of columns will be the number of diagonals to consider
    complex(kind=real64), allocatable :: Band_val_C(:,:)                       ! if subtype = "Band"
  END TYPE


!--------------------------------First cavity mode-----------------------------
  TYPE(MC_cavity_mode_t)  :: Cavity_mode_1
  TYPE(MC_operator_1D_t)  :: H_ho_cavity_mode_1                                ! matrix of the one-dimensional harmonic Hamiltonian associated with HO D
  !TYPE(MC_operator_1D_t)  :: x_ho_cavity_1
  !TYPE(MC_operator_1D_t)  :: N_ho_cavity_1


  !--------------------------------Just for test purposes----------------------
  integer                         :: j
  TYPE(MC_operator_1D_t)          :: H_ho_cavity_mode_1_diag
  TYPE(MC_operator_1D_t)          :: H_ho_cavity_mode_1_dense
  !real(kind=real64), allocatable :: Psi_Cavity_mode_1(:)                                ! a one dimension vector describing the photon state of the first mode. Vector on the HO_1 Eigenbasis. It's a superposition of the Eigenvector states
  !real(kind=real64), allocatable :: Psi1_result(:), Psi2_result(:), Psi3_result(:)
  !----------------------------------------------------------------------------


  CALL MolecCav_Read_cavity_mode(Mode=Cavity_mode_1, nio=INPUT_UNIT)

  CALL MolecCav_Construct_Operator(Operator=H_ho_cavity_mode_1, &
                                 & operator_type="Hamiltonian", &
                                 & scalar_space="Real", &
                                 & matrix_shape_type="Opt", &                  ! opt => get analytical shape. non_opt => get dense shape
                                 & Nb=Cavity_mode_1%Nb, &
                                 & w=Cavity_mode_1%w, &
                                 & m=Cavity_mode_1%m)

  !--------------------------------Just for test purposes----------------------
  WRITE(OUTPUT_UNIT,*) "H_{cavity mode 1}"
  DO j = 1, Cavity_mode_1%Nb
    WRITE(OUTPUT_UNIT,*) H_ho_cavity_mode_1%Diag_val_R(j)
  END DO

  CALL MolecCav_Construct_Operator(Operator=H_ho_cavity_mode_1_diag, &              ! to compare with non_opt
                                 & operator_type="Hamiltonian", &
                                 & scalar_space="Real", &
                                 & matrix_shape_type="Opt", &
                                 & Nb=Cavity_mode_1%Nb, &
                                 & w=Cavity_mode_1%w, &
                                 & m=Cavity_mode_1%m)

  WRITE(OUTPUT_UNIT,*) "H_{cavity mode 1 Optimized}"
  DO j = 1, Cavity_mode_1%Nb
    WRITE(OUTPUT_UNIT,*) H_ho_cavity_mode_1_diag%Diag_val_R(j)
  END DO

  CALL MolecCav_Construct_Operator(Operator=H_ho_cavity_mode_1_dense, &
                                 & operator_type="Hamiltonian", &
                                 & scalar_space="Real", &
                                 & matrix_shape_type="Non_opt", &                  ! opt => get analytical shape. non_opt => get dense shape
                                 & Nb=Cavity_mode_1%Nb, &
                                 & w=Cavity_mode_1%w, &
                                 & m=Cavity_mode_1%m)

  WRITE(OUTPUT_UNIT,*) "H_{cavity mode 1 Non Optimized}"
  DO j = 1, Cavity_mode_1%Nb
    WRITE(OUTPUT_UNIT,*) H_ho_cavity_mode_1_dense%Dense_val_R(j,:)
  END DO
  !----------------------------------------------------------------------------


  !CALL MolecCav_Construct_Operator(Operator=x_ho_cavity_1, operator_type="Position", &
  !                               & Nb_dimensions_array=2, scalar_space="Real", Mode=Cavity_mode_1)

  !WRITE(OUTPUT_UNIT,*) "x_{cavity mode 1}"
  !DO j = 1, Cavity_mode_1%Nb
  !  WRITE(OUTPUT_UNIT,*) x_ho_cavity_1%Diag_val_R(j)
  !END DO 

  !CALL MolecCav_Construct_Operator(Operator=N_ho_cavity_1, operator_type="Nb_photons", &
  !                               & Nb_dimensions_array=2, scalar_space="Real", Mode=Cavity_mode_1)

  !WRITE(OUTPUT_UNIT,*) "N_{cavity_mode_1}"
  !DO j = 1, Cavity_mode_1%Nb
  !  WRITE(OUTPUT_UNIT,*) N_ho_cavity_1%Diag_val_R(j)
  !END DO 

  !ALLOCATE(Psi_Cavity_mode_1(Cavity_mode_1%Nb))
  !DO j = Cavity_mode_1%Nb, 1, -1
  !  Psi_Cavity_mode_1(Cavity_mode_1%Nb - j + 1) = j
  !END DO
  !WRITE(OUTPUT_UNIT,*) "Psi_{cavity mode 1}"
  !DO j = 1, Cavity_mode_1%Nb
  !  WRITE(OUTPUT_UNIT,*) Psi_Cavity_mode_1(j)
  !END DO

  !ALLOCATE(Psi1_result(Cavity_mode_1%Nb))
  !CALL MolecCav_Action_H_cavity_mode_R(Psi_H_Cavity_R=Psi1_result, Mode=Cavity_mode_1, Psi_R=Psi_Cavity_mode_1)
  !WRITE(OUTPUT_UNIT,*) "Psi_{Hc}"
  !DO j = 1, Cavity_mode_1%Nb
  !  WRITE(OUTPUT_UNIT,*) Psi1_result(j)
  !END DO 


  !ALLOCATE(Psi2_result(Cavity_mode_1%Nb))
  !CALL MolecCav_Action_x_cavity_mode_R(Psi2_result, x_ho_cavity_1, Psi_Cavity_mode_1)
  !WRITE(OUTPUT_UNIT,*) "Psi_{xc}"
  !DO j = 1, Cavity_mode_1%Nb
  !  WRITE(OUTPUT_UNIT,*) Psi2_result(j)
  !END DO 

  !ALLOCATE(Psi3_result(Cavity_mode_1%Nb))
  !CALL MolecCav_Action_N_cavity_mode_R(Psi3_result, N_ho_cavity_1, Psi_Cavity_mode_1)
  !WRITE(OUTPUT_UNIT,*) "Psi_{Nc}"
  !DO j = 1, Cavity_mode_1%Nb
  !  WRITE(OUTPUT_UNIT,*) Psi3_result(j)
  !END DO 


  WRITE(OUTPUT_UNIT,*) ''
  WRITE(OUTPUT_UNIT,*) 'Test 1 checked !'



  
  CONTAINS




  SUBROUTINE MolecCav_Read_cavity_mode(Mode, nio)                              ! nio is the label of the file from which the values have to be drawn.
    USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    IMPLICIT NONE
    
    TYPE(MC_cavity_mode_t), intent(inout) :: Mode   
    integer,                intent(in)    :: nio

    integer                               :: D, Nb, err_io                     ! label of the basis/HO/mode/dimension, its number of basis vectors, and an error control variable
    real(kind=real64)                     :: w, m, lambda                      ! eigenpulsation, mass, and molecule-coupling strength associated with this HO

    NAMELIST /HO_1/ D, Nb, w, m, lambda                                        ! declare the nml HO_1 and specify the parameter's list to be found within

    !----------------------Initialization to default values---------------------
    D      = 1
    Nb     = 1
    w      = 1.0
    m      = 1.0
    lambda = 1.0
 
    !------------------------------Reading of the nml---------------------------
    WRITE(OUTPUT_UNIT,*) ''
    WRITE(OUTPUT_UNIT,*) '*******************************************************'
    WRITE(OUTPUT_UNIT,*) '************* READING BASIS OF THE HO *****************'
    WRITE(OUTPUT_UNIT,*) '*******************************************************'
    
    READ(nio, nml = HO_1, iostat = err_io)                                     ! assign the values read in the nml to the declared list of parameters

    WRITE(OUTPUT_UNIT, nml = HO_1)                                             ! just to have it in the output file
    
    !------------------------------Check reading error-------------------------
    IF(err_io < 0) then
      WRITE(OUTPUT_UNIT,*) ''
      WRITE(OUTPUT_UNIT,*) '#######################################################'
      WRITE(OUTPUT_UNIT,*) '######## Error in Read_cavity_mode (err_io<0) #########'
      WRITE(OUTPUT_UNIT,*) '#######################################################'
      WRITE(OUTPUT_UNIT,*) '################ err_io = ', err_io, '################'
      STOP '################# Check basis data ################'

    ELSE IF( err_io > 0) then
      WRITE(OUTPUT_UNIT,*) ''
      WRITE(OUTPUT_UNIT,*) '#######################################################'
      WRITE(OUTPUT_UNIT,*) '######## Error in Read_cavity_mode (err_io>0) #########'
      WRITE(OUTPUT_UNIT,*) '#######################################################'
      WRITE(OUTPUT_UNIT,*) '################ err_io = ', err_io, '################'
      STOP '################# Check basis data ################'

    END IF

    !---------------Construction of the MC_cavity_mode_t type object-----------
    Mode%D      = D
    Mode%Nb     = Nb
    Mode%w      = w
    Mode%m      = m
    Mode%lambda = lambda

    WRITE(OUTPUT_UNIT,*) ''
    WRITE(OUTPUT_UNIT,*) '*******************************************************'
    WRITE(OUTPUT_UNIT,*) '*********** BASIS OF THE HO CONSTRUCTED ***************'
    WRITE(OUTPUT_UNIT,*) '*******************************************************'

  END SUBROUTINE
 

  SUBROUTINE MolecCav_Construct_Operator(Operator, operator_type, scalar_space, matrix_shape_type, Nb, w, m)
    USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    IMPLICIT NONE

    TYPE(MC_operator_1D_t), intent(inout) :: Operator                          ! the object of type Operator_t to be constructed here
    character(len=*),       intent(in)    :: operator_type                     ! ex : "Hamiltonian", "Position", etc
    character(len=*),       intent(in)    :: scalar_space                      ! ex : "Real" or "Complex"
    character(len=*),       intent(in)    :: matrix_shape_type                 ! the assignement to use the analytical shape of the operator's matrix, or the dense shape
    integer,                intent(in)    :: Nb                                ! number of basis vectors associated with the HO of whom Operator is the operator
    real(kind=real64),      intent(in)    :: w, m                              ! eigenpulsation and mass associated with this HO

    !--------------first steps of the construction of the Operator-------------
    Operator%operator_type = operator_type                                     ! allocation on assignement. operator type has the right lengths (no spaces added) thanks to len=* at declaration and it will fit the Op%op_type thanks to len=:, allocatable at declaration of the derived type. 
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
    USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    IMPLICIT NONE
    
    TYPE(MC_operator_1D_t), intent(inout) :: Hamiltonian                       ! matrix of the one-dimensional harmonic Hamiltonian associated with HO D
    character(len=*),       intent(in)    :: matrix_shape_type                 ! the assignement to use the analytical shape of the operator's matrix, or the dense shape
    integer,                intent(in)    :: Nb                                ! number of basis vectors associated with the HO of whom Operator is the operator
    real(kind=real64),      intent(in)    :: w                                 ! eigenpulsation associated with this HO

    integer                               :: i                                 ! loop increments /!\ Fortran counts from 1 to Nb !!! /!\

    WRITE(OUTPUT_UNIT,*) ''
    WRITE(OUTPUT_UNIT,*) '*******************************************************'
    WRITE(OUTPUT_UNIT,*) '******* CONSTRUCTING THE HAMILTONIAN OF THE HO ********'

    IF (matrix_shape_type == "Opt") THEN
      !----------Initialization of the characteristics of the operator---------
      Hamiltonian%matrix_shape_type = "Diagonal"
      !---------------------Initialization to default values-------------------
      ALLOCATE(Hamiltonian%Diag_val_R(Nb))
      Hamiltonian%Diag_val_R = 0.0_real64
      !------------------------Construction of the matrix----------------------
      DO i = 1, Nb                                                             ! /!\ Fortran counts from 1 to Nb !!! /!\
        Hamiltonian%Diag_val_R(i) = w*(i - 1 + 0.5_real64)                     ! "-1" because the first Fortran vector is the fundamental eigenvector of the HO i.e. the 0^{th} ket 
      END DO
      
    ELSE IF (matrix_shape_type == "Non_opt") THEN
      !----------Initialization of the characteristics of the operator---------
      Hamiltonian%matrix_shape_type = "Dense"
      !---------------------Initialization to default values-------------------
      ALLOCATE(Hamiltonian%Dense_val_R(Nb, Nb))
      Hamiltonian%Dense_val_R = 0.0_real64
      !------------------------Construction of the matrix----------------------
      DO i = 1, Nb                                                             ! /!\ Fortran counts from 1 to Nb !!! /!\
        Hamiltonian%Dense_val_R(i,i) = w*(i - 1 + 0.5_real64)                     ! "-1" because the first Fortran vector is the fundamental eigenvector of the HO i.e. the 0^{th} ket 
      END DO
    END IF
    
    WRITE(OUTPUT_UNIT,*) '********* HAMILTONIAN OF THE HO CONSTRUCTED ***********'
    WRITE(OUTPUT_UNIT,*) '*******************************************************'

  END SUBROUTINE


  SUBROUTINE MolecCav_Construct_x_cavity_mode(Position_Op, matrix_shape_type, Nb, w, m)
    USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    IMPLICIT NONE
    
    TYPE(MC_operator_1D_t), intent(inout) :: Position_Op
    character(len=*),       intent(in)    :: matrix_shape_type                 ! the assignement to use the analytical shape of the operator's matrix, or the dense shape
    integer,                intent(in)    :: Nb                                ! number of basis vectors associated with the HO of whom Operator is the operator
    real(kind=real64),      intent(in)    :: w, m                              ! eigenpulsation and mass associated with this HO

    WRITE(OUTPUT_UNIT,*) '*******************************************************'

  END SUBROUTINE


  SUBROUTINE MolecCav_Construct_N_cavity_mode(Nb_photon_Op, matrix_shape_type, Nb)
    USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    IMPLICIT NONE
    
    TYPE(MC_operator_1D_t), intent(inout) :: Nb_photon_Op
    character(len=*),       intent(in)    :: matrix_shape_type                 ! the assignement to use the analytical shape of the operator's matrix, or the dense shape
    integer,                intent(in)    :: Nb                                ! number of basis vectors associated with the HO of whom Operator is the operator

    WRITE(OUTPUT_UNIT,*) '*******************************************************'

  END SUBROUTINE


END PROGRAM
