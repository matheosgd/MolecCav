PROGRAM test1
  USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
  IMPLICIT NONE

  TYPE :: MC_cavity_mode_t                           ! MC = MolecCav
    integer                        :: D              ! label of the OHL/mode/dimension/associated basis set
    integer                        :: Nb_D           ! number of basis vectors associated with the OHL D
    real(kind=real64)              :: w_D            ! eigenpulsation associated with the OHL D
    real(kind=real64)              :: m_D            ! mass associated with the OHL D
    real(kind=real64), allocatable :: H_1Dohl_D(:,:) ! matrix of the one-dimensional harmonic Hamiltonian associated with OHL D
  END TYPE


!------------------------------First cavity mode---------------------------

  TYPE(MC_cavity_mode_t) :: Cavity_mode_1

  CALL MolecCav_Read_cavity_mode(Mode=Cavity_mode_1, nio=INPUT_UNIT)

  CALL MolecCav_Construct_H_cavity_mode(Mode=Cavity_mode_1)

  WRITE(OUTPUT_UNIT,*) ''
  WRITE(OUTPUT_UNIT,*) 'Test 1 checked !'



  CONTAINS



  SUBROUTINE MolecCav_Read_cavity_mode(Mode, nio)                           ! nio is the label of the file from which the values have to be drawn.
    USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    IMPLICIT NONE
    
    TYPE(MC_cavity_mode_t), intent(inout) :: Mode   
    integer,                intent(in)    :: nio

    integer                               :: D, Nb_D, err_io                ! label of the basis/OHL/mode/dimension, its number of basis vectors, and an error control variable
    real(kind=real64)                     :: w_D, m_D                       ! eigenpulsation and mass associated with this OHL

    NAMELIST /PARAMETERS_OHL_1/ D, Nb_D, w_D, m_D                           ! declare the nml PARAMETERS_OHL_1 and specify the parameter's list to be found within

!----------------------Initialization to default values---------------------
    D    = 1
    Nb_D = 10
    w_D  = 1.0
    m_D  = 1.0
 
!------------------------------Reading of the nml---------------------------
    WRITE(OUTPUT_UNIT,*) ''
    WRITE(OUTPUT_UNIT,*) '*******************************************************'
    WRITE(OUTPUT_UNIT,*) '************* READING BASIS OF THE OHL ****************'
    WRITE(OUTPUT_UNIT,*) '*******************************************************'
    
    READ(nio, nml = PARAMETERS_OHL_1, iostat = err_io)                      ! assign the values read in the nml to the declared list of parameters

    WRITE(OUTPUT_UNIT, nml = PARAMETERS_OHL_1)                              ! just to have it in the output file
    
!------------------------------Check reading error--------------------------
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

!---------------Construction of the MC_cavity_mode_t type object------------
    Mode%D    = D
    Mode%Nb_D = Nb_D
    Mode%w_D  = w_D
    Mode%m_D  = m_D

    WRITE(OUTPUT_UNIT,*) ''
    WRITE(OUTPUT_UNIT,*) '*******************************************************'
    WRITE(OUTPUT_UNIT,*) '*********** BASIS OF THE OHL CONSTRUCTED **************'
    WRITE(OUTPUT_UNIT,*) '*******************************************************'

  END SUBROUTINE
 
  SUBROUTINE MolecCav_Construct_H_cavity_mode(Mode)
    USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    IMPLICIT NONE
    
    TYPE(MC_cavity_mode_t), intent(inout) :: Mode

    integer                               :: i, j                           ! loop increments /!\ Fortran counts from 1 to Nb !!! /!\

!----------------------Initialization to default values---------------------
    WRITE(OUTPUT_UNIT,*) ''
    WRITE(OUTPUT_UNIT,*) '*******************************************************'
    WRITE(OUTPUT_UNIT,*) '******* CONSTRUCTING THE HAMILTONIAN OF THE OHL *******'
    WRITE(OUTPUT_UNIT,*) '*******************************************************'

    ALLOCATE(Mode%H_1Dohl_D(Mode%Nb_D,Mode%Nb_D))
    Mode%H_1Dohl_D(:,:) = 0
    
    !WRITE(OUTPUT_UNIT,*) "Initialized matrix :"
    !DO i = 1, Mode%Nb_D
    !  WRITE(OUTPUT_UNIT,*) Mode%H_1Dohl_D(i,:)
    !END DO 

!-------------------------Construction of the matrix------------------------
    DO i = 1, Mode%Nb_D                                                     ! /!\ Fortran counts from 1 to Nb !!! /!\
      Mode%H_1Dohl_D(i,i) = Mode%w_D*(i - 1 + 0.5_real64)                   ! -1 because the first Fortran vector is the fundamental eigenvector of the OHL i.e. the 0^{th} ket 
    END DO

    !WRITE(OUTPUT_UNIT,*) "Constructed matrix :"
    !DO i = 1, Mode%Nb_D
    !  WRITE(OUTPUT_UNIT,*) Mode%H_1Dohl_D(i,:)
    !END DO

    WRITE(OUTPUT_UNIT,*) ''
    WRITE(OUTPUT_UNIT,*) '*******************************************************'
    WRITE(OUTPUT_UNIT,*) '********* HAMILTONIAN OF THE OHL CONSTRUCTED **********'
    WRITE(OUTPUT_UNIT,*) '*******************************************************'

  END SUBROUTINE

END PROGRAM
