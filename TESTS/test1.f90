PROGRAM test1
  USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
  IMPLICIT NONE

  TYPE :: MC_cavity_mode_t                                                     ! MC = MolecCav
    integer                        :: D                                        ! label of the HO/mode/dimension/associated basis set
    integer                        :: Nb                                       ! number of basis vectors associated with the HO D
    real(kind=real64)              :: w                                        ! eigenpulsation associated with the HO D
    real(kind=real64)              :: m                                        ! mass associated with the HO D
    real(kind=real64)              :: lambda                                   ! Strength parameter of the coupling between the mode D and the molecule
    real(kind=real64), allocatable :: H_ho(:,:)                                ! matrix of the one-dimensional harmonic Hamiltonian associated with HO D
  END TYPE

  TYPE :: MC_operator_t
    character(len=20)                 :: operator_type                         ! ex : "Hamiltonian", "Position", etc
    integer                           :: Nb_dimensions_array                   ! /!\ This is the array dimension and NOT the number of cavity modes !!! Consequently, this integer is not supposed to be lower than 2 for any matrix operator /!\
    integer,              allocatable :: Dimensions(:)                         ! the shape of the matrix. Ex : for a 2D_array, Dimensions = (Cavity_mode%Nb-1, Cavity_mode%Nb-1)
    character(len=10)                 :: Scalar_space                          ! ex : "Real" or "Complex"
    !### NB : The following types concern mainly the two-dimension matrices ###
    character(len=20)                 :: matrix_shape_type                     ! ex : "Dense", "Diagonal", "Band", etc
    character(len=20)                 :: matrix_shape_subtype                  ! if type = "Band". Ex : "Tridiagonal", "Up_Triangular", "Low_Triangular", "Other", etc
    integer                           :: Upper_bandwidth                       ! if subtype = "Other". Ex : Upper_bandwidth=Lower_bandwidth=1 would give a tridiagonal matrix
    integer                           :: Lower_bandwidth                       ! if subtype = "Other"
    real(kind=real64),    allocatable :: Dense_val_R(:,:)                      ! if type = "Dense"
    complex(kind=real64), allocatable :: Dense_val_C(:,:)                      ! if type = "Dense"
    real(kind=real64),    allocatable :: Diag_val_R(:)                         ! if type = "Diagonal"
    complex(kind=real64), allocatable :: Diag_val_C(:)                         ! if type = "Diagonal"
    real(kind=real64),    allocatable :: Tridiag_val_R(:,:,:)                  ! if subtype = "Tridiagonal"
    complex(kind=real64), allocatable :: Tridiag_val_C(:,:,:)                  ! if subtype = "Tridiagonal"
    real(kind=real64),    allocatable :: Other_val_R(:)                        ! if subtype = "Other"
    complex(kind=real64), allocatable :: Other_val_C(:)                        ! if subtype = "Other"
  END TYPE
!--------------------------------First cavity mode-----------------------------

  TYPE(MC_cavity_mode_t)         :: Cavity_mode_1
  TYPE(MC_operator_t)            :: H_ho_cavity_mode_1
  real(kind=real64), allocatable :: x_ho_cavity_1(:,:), N_ho_cavity_1(:,:)
  real(kind=real64), allocatable :: Psi_Cavity_mode_1(:)                                      ! for tests purposes. A one dimension vector describing the photon state of the first mode. Vector on the HO_1 Eigenbasis. It's a superposition of the Eigenvector states
  real(kind=real64), allocatable :: Psi1_result(:), Psi2_result(:), Psi3_result(:)            ! just to store the result of some tests
  integer                        :: j

  CALL MolecCav_Read_cavity_mode(Mode=Cavity_mode_1, nio=INPUT_UNIT)

  CALL MolecCav_Construct_H_cavity_mode(Operator=H_ho_cavity_mode_1, Mode=Cavity_mode_1)

  WRITE(OUTPUT_UNIT,*) "H_{cavity mode 1}"
  DO j = 1, Cavity_mode_1%Nb
    WRITE(OUTPUT_UNIT,*) H_ho_cavity_mode_1(j,:)
  END DO

  ALLOCATE(Psi_Cavity_mode_1(Cavity_mode_1%Nb))
  DO j = Cavity_mode_1%Nb, 1, -1
    Psi_Cavity_mode_1(Cavity_mode_1%Nb - j + 1) = j
  END DO
  WRITE(OUTPUT_UNIT,*) "Psi_{cavity mode 1}"
  DO j = 1, Cavity_mode_1%Nb
    WRITE(OUTPUT_UNIT,*) Psi_Cavity_mode_1(j)
  END DO

  ALLOCATE(Psi1_result(Cavity_mode_1%Nb))
  CALL MolecCav_Action_H_cavity_mode_R(Psi_H_Cavity_R=Psi1_result, , Mode=Cavity_mode_1, Psi_R=Psi_Cavity_mode_1)
  WRITE(OUTPUT_UNIT,*) "Psi_{Hc}"
  DO j = 1, Cavity_mode_1%Nb
    WRITE(OUTPUT_UNIT,*) Psi1_result(j)
  END DO 

  CALL MolecCav_Construct_x_cavity_mode(x_ho_cavity_1, Cavity_mode_1)
  WRITE(OUTPUT_UNIT,*) "x_{cavity mode 1}"
  DO j = 1, Cavity_mode_1%Nb
    WRITE(OUTPUT_UNIT,*) x_ho_cavity_1(j,:)
  END DO 

  ALLOCATE(Psi2_result(Cavity_mode_1%Nb))
  CALL MolecCav_Action_x_cavity_mode_R(Psi2_result, x_ho_cavity_1, Psi_Cavity_mode_1)
  WRITE(OUTPUT_UNIT,*) "Psi_{xc}"
  DO j = 1, Cavity_mode_1%Nb
    WRITE(OUTPUT_UNIT,*) Psi2_result(j)
  END DO 

  CALL MolecCav_Construct_N_cavity_mode(N_ho_cavity_1, Cavity_mode_1)
  WRITE(OUTPUT_UNIT,*) "N_{cavity_mode_1}"
  DO j = 1, Cavity_mode_1%Nb
    WRITE(OUTPUT_UNIT,*) N_ho_cavity_1(j,:)
  END DO 

  ALLOCATE(Psi3_result(Cavity_mode_1%Nb))
  CALL MolecCav_Action_N_cavity_mode_R(Psi3_result, N_ho_cavity_1, Psi_Cavity_mode_1)
  WRITE(OUTPUT_UNIT,*) "Psi_{Nc}"
  DO j = 1, Cavity_mode_1%Nb
    WRITE(OUTPUT_UNIT,*) Psi3_result(j)
  END DO 

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
 

  SUBROUTINE MolecCav_Construct_H_cavity_mode(Operator, Mode)
    USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    IMPLICIT NONE
    
    TYPE(MC_operator_t), intent(inout) :: Operator
    TYPE(MC_cavity_mode_t), intent(in) :: Mode

    integer                               :: i                                 ! loop increments /!\ Fortran counts from 1 to Nb !!! /!\

    !----------------------Initialization to default values--------------------
    WRITE(OUTPUT_UNIT,*) ''
    WRITE(OUTPUT_UNIT,*) '*******************************************************'
    WRITE(OUTPUT_UNIT,*) '******* CONSTRUCTING THE HAMILTONIAN OF THE HO ********'
    WRITE(OUTPUT_UNIT,*) '*******************************************************'

    Operator%operator_type = "Hamiltonian"
    Operator%Nb_dimensions_array = 2
    ALLOCATE(Operator%Dimensions(Operator%Nb_dimensions_array))
    Operator%Dimensions(:) = Mode%Nb
    Operator%Scalar_space = "Real"
    Operator%matrix_shape_type = "Diagonal"
    ALLOCATE(Operator%Diag_val_R(Mode%Nb))

    Operator%Diag_val_R = 0.0_real64
    
    !WRITE(OUTPUT_UNIT,*) "Initialized matrix :"
    !DO i = 1, Mode%Nb
    !  WRITE(OUTPUT_UNIT,*) Mode%H_ho(i,:)
    !END DO 

    !-------------------------Construction of the matrix-----------------------
    DO i = 1, Mode%Nb                                                          ! /!\ Fortran counts from 1 to Nb !!! /!\
      Operator%Diag_val_R(i) = Mode%w*(i - 1 + 0.5_real64)                     ! "-1" because the first Fortran vector is the fundamental eigenvector of the HO i.e. the 0^{th} ket 
    END DO

    !WRITE(OUTPUT_UNIT,*) "Constructed matrix :"
    !DO i = 1, Mode%Nb
    !  WRITE(OUTPUT_UNIT,*) Mode%H_ho(i,:)
    !END DO

    WRITE(OUTPUT_UNIT,*) ''
    WRITE(OUTPUT_UNIT,*) '*******************************************************'
    WRITE(OUTPUT_UNIT,*) '********* HAMILTONIAN OF THE HO CONSTRUCTED ***********'
    WRITE(OUTPUT_UNIT,*) '*******************************************************'

  END SUBROUTINE


  SUBROUTINE MolecCav_Construct_x_cavity_mode(x_ho, Mode)
    USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    IMPLICIT NONE
    
    real(kind=real64), allocatable, intent(inout) :: x_ho(:,:)
    TYPE(MC_cavity_mode_t), intent(in)            :: Mode

    integer                                       :: i                         ! loop increments /!\ Fortran counts from 1 to Nb !!! /!\

    !----------------------Initialization to default values--------------------
    WRITE(OUTPUT_UNIT,*) ''
    WRITE(OUTPUT_UNIT,*) '*******************************************************'
    WRITE(OUTPUT_UNIT,*) '******* CONSTRUCTING THE POSITION OP OF THE HO ********'
    WRITE(OUTPUT_UNIT,*) '*******************************************************'

    ALLOCATE(x_ho(Mode%Nb,Mode%Nb))
    x_ho(:,:) = 0.0_real64

    !WRITE(OUTPUT_UNIT,*) "Initialized matrix :"
    !DO i = 1, Mode%Nb
    !  WRITE(OUTPUT_UNIT,*) Mode%H_ho(i,:)
    !END DO 

    !-------------------------Construction of the matrix-----------------------
    DO i = 1, Mode%Nb - 1                                                      ! /!\ Fortran counts from 1 to Nb !!! /!\
      x_ho(i,i+1) = SQRT(REAL(i,kind=real64))
      x_ho(i+1,i) = SQRT(REAL(i,kind=real64))
    END DO
    x_ho(:,:) = x_ho(:,:) / SQRT(2.0_real64 * Mode%w)

    !WRITE(OUTPUT_UNIT,*) "Constructed matrix :"
    !DO i = 1, Mode%Nb
    !  WRITE(OUTPUT_UNIT,*) Mode%x_ho(i,:)
    !END DO

    WRITE(OUTPUT_UNIT,*) ''
    WRITE(OUTPUT_UNIT,*) '*******************************************************'
    WRITE(OUTPUT_UNIT,*) '********* POSITION OP OF THE HO CONSTRUCTED ***********'
    WRITE(OUTPUT_UNIT,*) '*******************************************************'

  END SUBROUTINE


  SUBROUTINE MolecCav_Construct_N_cavity_mode(N_ho, Mode)
    USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    IMPLICIT NONE
    
    real(kind=real64), allocatable, intent(inout) :: N_ho(:,:)
    TYPE(MC_cavity_mode_t), intent(in)            :: Mode

    integer                                       :: i                         ! loop increments /!\ Fortran counts from 1 to Nb !!! /!\

    !----------------------Initialization to default values--------------------
    WRITE(OUTPUT_UNIT,*) ''
    WRITE(OUTPUT_UNIT,*) '*******************************************************'
    WRITE(OUTPUT_UNIT,*) '******* CONSTRUCTING THE NB QUANTA OP OF THE HO *******'
    WRITE(OUTPUT_UNIT,*) '*******************************************************'

    ALLOCATE(N_ho(Mode%Nb,Mode%Nb))
    N_ho(:,:) = 0.0_real64

    !WRITE(OUTPUT_UNIT,*) "Initialized matrix :"
    !DO i = 1, Mode%Nb
    !  WRITE(OUTPUT_UNIT,*) N_ho(i,:)
    !END DO 

    !-------------------------Construction of the matrix-----------------------
    DO i = 1, Mode%Nb                                                          ! /!\ Fortran counts from 1 to Nb !!! /!\
      N_ho(i,i) = i - 1
    END DO

    !WRITE(OUTPUT_UNIT,*) "Constructed matrix :"
    !DO i = 1, Mode%Nb
    !  WRITE(OUTPUT_UNIT,*) N_ho(i,:)
    !END DO

    WRITE(OUTPUT_UNIT,*) ''
    WRITE(OUTPUT_UNIT,*) '*******************************************************'
    WRITE(OUTPUT_UNIT,*) '********* NB QUANTA OP OF THE HO CONSTRUCTED **********'
    WRITE(OUTPUT_UNIT,*) '*******************************************************'

  END SUBROUTINE


  SUBROUTINE MolecCav_Action_H_cavity_mode_R(Psi_H_Cavity_R, Mode, Psi_R)      ! _R for the case where Psi is a real vector. Compute the resulting vector Psi_H_Cavity_R(:) from the action of the cavity Hamiltonian of the mode Mode%H_ho(:,:) on the photon state vector Psi(:) associated to the mode of the Hamiltonian (written in the Eigenbasis of H_ho)
    USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    IMPLICIT NONE
    
    real(kind=real64), intent(inout)                :: Psi_H_Cavity_R(:)
    TYPE(MC_cavity_mode_t), intent(in)              :: Mode
    real(kind=real64), intent(in)                   :: Psi_R(:)

    Psi_H_Cavity_R(:) = matmul(Mode%H_ho, Psi_R)

  END SUBROUTINE


  SUBROUTINE MolecCav_Action_x_cavity_mode_R(Psi_x_Cavity_R, x_ho, Psi_R)      ! _R for the case where Psi is a real vector. Compute the resulting vector Psi_x_Cavity_R(:) from the action of the 1D-position operator of the mode x_ho(:,:) on the photon state vector Psi(:) associated to the mode of the Hamiltonian (written in the Eigenbasis of H_ho)
    USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    IMPLICIT NONE
    
    real(kind=real64), intent(inout)                :: Psi_x_Cavity_R(:)
    real(kind=real64), intent(in)                   :: Psi_R(:), x_ho(:,:)
        
    Psi_x_Cavity_R(:) = matmul(x_ho, Psi_R)
    
  END SUBROUTINE
  

  SUBROUTINE MolecCav_Action_N_cavity_mode_R(Psi_N_Cavity_R, N_ho, Psi_R)      ! _R for the case where Psi is a real vector. Compute the resulting vector Psi_H_Cavity_R(:) from the action of the cavity Hamiltonian of the mode Mode%H_ho(:,:) on the photon state vector Psi(:) associated to the mode of the Hamiltonian (written in the Eigenbasis of H_ho)
    USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    IMPLICIT NONE
    
    real(kind=real64), intent(inout)                :: Psi_N_Cavity_R(:)
    real(kind=real64), intent(in)                   :: Psi_R(:), N_ho(:,:)

    Psi_N_Cavity_R(:) = matmul(N_ho, Psi_R)

  END SUBROUTINE

END PROGRAM
