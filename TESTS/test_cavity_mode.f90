PROGRAM test_cavity_mode
  USE QDUtil_m
  USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
  USE MC_cavity_mode_m
  IMPLICIT NONE


!--------------------------------First cavity mode-----------------------------
  TYPE(MC_cavity_mode_t)  :: Cavity_mode_1
  integer                 :: error


  error = 0
  CALL MolecCav_Read_cavity_mode(Mode=Cavity_mode_1, nio=INPUT_UNIT)


  IF (Cavity_mode_1%D == 0) THEN
    WRITE(OUTPUT_UNIT,*) ''
    WRITE(OUTPUT_UNIT,*) 'Mode%D failed to initialize'
    error = 1
  END IF
  IF (Cavity_mode_1%Nb == 1) THEN
    WRITE(OUTPUT_UNIT,*) ''
    WRITE(OUTPUT_UNIT,*) 'Mode%Nb failed to initialize'
    error = 1
  END IF
  IF (Cavity_mode_1%w == 0) THEN
    WRITE(OUTPUT_UNIT,*) ''
    WRITE(OUTPUT_UNIT,*) 'Mode%w failed to initialize'
    error = 1
  END IF
  IF (Cavity_mode_1%m == 0) THEN
    WRITE(OUTPUT_UNIT,*) ''
    WRITE(OUTPUT_UNIT,*) 'Mode%m failed to initialize'
    error = 1
  END IF
  IF (Cavity_mode_1%lambda == 0) THEN
    WRITE(OUTPUT_UNIT,*) ''
    WRITE(OUTPUT_UNIT,*) 'Mode%lambda failed to initialize'
    error = 1
  END IF

  IF (error == 0) THEN
    WRITE(OUTPUT_UNIT,*) ''
    WRITE(OUTPUT_UNIT,*) 'Test 1 checked ! The Cavity mode initializes successfully !'
  ELSE IF (error == 1) THEN
    WRITE(OUTPUT_UNIT,*) ''
    WRITE(OUTPUT_UNIT,*) 'Test 1 failed ! The Cavity mode did not initialize successfully...'
  ELSE
    WRITE(OUTPUT_UNIT,*) ''
    WRITE(OUTPUT_UNIT,*) 'Test 1 stopped ! The file did not execute as it should have !'
  END IF
  
  
END PROGRAM
