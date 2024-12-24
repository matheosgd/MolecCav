MODULE MC_cavity_mode_m
  !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
  USE QDUtil_m                                                                 ! gives Rkind=real64; out_unit=OUTPUT_UNIT; INPUT_UNIT=in_unit; EYE=i and other numbers; TO_LOWERCASE; TO_UPPERCASE;... We thereby use ZERO instead of 0.0_real64
  IMPLICIT NONE


  TYPE :: MC_cavity_mode_t                                                     ! MC = MolecCav NB: everything is initialized at values that are not supposed to make it possible of the cavity mode lecture/creation have successfully been executed
    integer          :: D      = 0                                             ! label of the HO/mode/dimension/associated basis set
    integer          :: Nb     = 1                                             ! number of basis vectors associated with the HO D
    real(kind=Rkind) :: w      = ZERO                                          ! eigenpulsation associated with the HO D
    real(kind=Rkind) :: m      = ZERO                                          ! mass associated with the HO D
    real(kind=Rkind) :: lambda = ZERO                                          ! strength parameter of the coupling between the mode D and the molecule
  END TYPE


  CONTAINS


  SUBROUTINE MolecCav_Read_cavity_mode(Mode, nio)                              ! nio is the label of the file from which the values have to be drawn.
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    IMPLICIT NONE
    
    TYPE(MC_cavity_mode_t), intent(inout) :: Mode   
    integer,                intent(in)    :: nio

    integer                               :: D, Nb, err_io                     ! label of the basis/HO/mode/dimension, its number of basis vectors, and an error control variable
    real(kind=Rkind)                      :: w, m, lambda                      ! eigenpulsation, mass, and molecule-coupling strength associated with this HO

    NAMELIST /HO_1/ D, Nb, w, m, lambda                                        ! declare the nml HO_1 and specify the parameter's list to be found within

    !----------------------Initialization to default values---------------------
    D      = 1
    Nb     = 1
    w      = ONE
    m      = ONE
    lambda = ONE
 
    !------------------------------Reading of the nml---------------------------
    WRITE(out_unit,*) ''
    WRITE(out_unit,*) '*******************************************************'
    WRITE(out_unit,*) '************* READING BASIS OF THE HO *****************'
    WRITE(out_unit,*) '*******************************************************'
    
    READ(nio, nml = HO_1, iostat = err_io)                                     ! assign the values read in the nml to the declared list of parameters

    WRITE(out_unit, nml = HO_1)                                             ! just to have it in the output file
    
    !------------------------------Check reading error-------------------------
    IF(err_io < 0) then
      WRITE(out_unit,*) ''
      WRITE(out_unit,*) '#######################################################'
      WRITE(out_unit,*) '######## Error in Read_cavity_mode (err_io<0) #########'
      WRITE(out_unit,*) '#######################################################'
      WRITE(out_unit,*) '################ err_io = ', err_io, '################'
      STOP '################# Check basis data ################'

    ELSE IF( err_io > 0) then
      WRITE(out_unit,*) ''
      WRITE(out_unit,*) '#######################################################'
      WRITE(out_unit,*) '######## Error in Read_cavity_mode (err_io>0) #########'
      WRITE(out_unit,*) '#######################################################'
      WRITE(out_unit,*) '################ err_io = ', err_io, '################'
      STOP '################# Check basis data ################'

    END IF

    !---------------Construction of the MC_cavity_mode_t type object-----------
    Mode%D      = D
    Mode%Nb     = Nb
    Mode%w      = w
    Mode%m      = m
    Mode%lambda = lambda

    WRITE(out_unit,*) ''
    WRITE(out_unit,*) '*******************************************************'
    WRITE(out_unit,*) '*********** BASIS OF THE HO CONSTRUCTED ***************'
    WRITE(out_unit,*) '*******************************************************'

  END SUBROUTINE

  
END MODULE