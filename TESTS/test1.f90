PROGRAM test1
  USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
  IMPLICIT NONE

  TYPE :: MC_harmo_basis_t    ! MC = MolecCav
    integer           :: D    ! label of the basis/OHL/mode/dimension
    integer           :: Nb_D ! number of basis vectors associated with the OHL D
    real(kind=real64) :: w_D  ! eigenpulsation associated with the OHL D
    real(kind=real64) :: m_D  ! mass associated with the OHL D
  END TYPE

  TYPE(MC_harmo_basis_t) :: Basis_mode_1

  CALL MolecCav_read_cav_mode_Basis(Basis=Basis_mode_1, nio=INPUT_UNIT)

  WRITE(OUTPUT_UNIT,*) ''
  WRITE(OUTPUT_UNIT,*) 'Test 1 checked !'



  CONTAINS



  SUBROUTINE MolecCav_read_cav_mode_Basis(Basis, nio)                                         ! nio le fichier dans lequel prendre les valeurs
    USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    IMPLICIT NONE
    
    TYPE(MC_harmo_basis_t), intent(inout) :: Basis   
    integer,                intent(in)    :: nio

    integer                               :: D, Nb_D, err_io                ! label of the basis/OHL/mode/dimension, its number of basis vectors, and an error control variable
    real(kind=real64)                     :: w_D, m_D                       ! eigenpulsation and mass associated with this OHL

    NAMELIST /BASIS_OHL_1/ D, Nb_D, w_D, m_D                                ! declare the nml BASIS_OHL_1 and specify the parameter's list to be found within

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

    READ(nio, nml = BASIS_OHL_1, iostat = err_io)                           ! assign the values read in the nml to the declared list of parameters

    WRITE(OUTPUT_UNIT, nml = BASIS_OHL_1)                                   ! just to have it in the output file
    
!------------------------------Check reading error--------------------------
    IF(err_io < 0) then
      WRITE(OUTPUT_UNIT,*) ''
      WRITE(OUTPUT_UNIT,*) '#######################################################'
      WRITE(OUTPUT_UNIT,*) '################# Error in read_Basis #################'
      WRITE(OUTPUT_UNIT,*) '#######################################################'
      STOP '################### Check basis data ##################'

    ELSE IF( err_io > 0) then
      WRITE(OUTPUT_UNIT,*) ''
      WRITE(OUTPUT_UNIT,*) '#######################################################'
      WRITE(OUTPUT_UNIT,*) '################# Error in read_Basis #################'
      WRITE(OUTPUT_UNIT,*) '#######################################################'
      STOP '################### Check basis data ##################'

    END IF

!---------------Construction of the MC_harmo_basis_t type object------------
    Basis%D    = D
    Basis%Nb_D = Nb_D
    Basis%w_D  = w_D
    Basis%m_D  = m_D

    WRITE(OUTPUT_UNIT,*) ''
    WRITE(OUTPUT_UNIT,*) '*******************************************************'
    WRITE(OUTPUT_UNIT,*) '*********** BASIS OF THE OHL CONSTRUCTED **************'
    WRITE(OUTPUT_UNIT,*) '*******************************************************'

  END SUBROUTINE
 
END PROGRAM
