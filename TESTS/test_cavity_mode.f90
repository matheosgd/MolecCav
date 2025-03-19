!==================================================================================================
!==================================================================================================
! This file is part of MolecCav.
!
!==================================================================================================
! MIT License
!
! Copyright (c) 2025 Mathéo Segaud
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
!==================================================================================================
PROGRAM test_cavity_mode
  !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
  USE QDUtil_m
  USE QDUtil_Test_m
  USE Cavity_mode_m
  IMPLICIT NONE


  TYPE(Cavity_mode_t) :: Cavity_mode_1
  logical             :: Debug = .FALSE.
  logical             :: error = .FALSE.
  TYPE(test_t)        :: test_cavity_init


  !-----------------------------Test initialization----------------------------
  CALL Initialize_Test(test_cavity_init, test_name="OUT/test_file_cav_mode")


  !-------------------------Cavity mode initialization-------------------------
  CALL MolecCav_Read_cavity_mode(Mode=Cavity_mode_1, nio=in_unit)
  IF (Debug) THEN
    WRITE(out_unit,*)
    WRITE(out_unit,*) "--------------Cavity mode constructed by MolecCav_Read_cavity_mode--------------"
    CALL Write_cavity_mode(Cavity_mode_1)
    WRITE(out_unit,*) "------------End Cavity mode constructed by MolecCav_Read_cavity_mode------------"
  END IF


  !-----------------------------------The tests--------------------------------
  error = (Cavity_mode_1%D == 0)
  CALL Logical_Test(test_cavity_init, error, test2=.FALSE., info="Cavity_mode_1%D = 0")
  IF (error .AND. Debug) WRITE(out_unit,*) 'Mode%D failed to initialize'

  error = (Cavity_mode_1%Nb == 0)
  CALL Logical_Test(test_cavity_init, error, test2=.FALSE., info="Cavity_mode_1%Nb == 0")
  IF (error .AND. Debug) WRITE(out_unit,*) 'Mode%Nb failed to initialize'

  error = (Cavity_mode_1%w == 0)
  CALL Logical_Test(test_cavity_init, error, test2=.FALSE., info="Cavity_mode_1%w == 0")
  IF (error .AND. Debug) WRITE(out_unit,*) 'Mode%w failed to initialize'

  error = (Cavity_mode_1%m == 0)
  CALL Logical_Test(test_cavity_init, error, test2=.FALSE., info="Cavity_mode_1%m == 0")
  IF (error .AND. Debug) WRITE(out_unit,*) 'Mode%m failed to initialize'

  error = (Cavity_mode_1%lambda < 0)
  CALL Logical_Test(test_cavity_init, error, test2=.FALSE., info="Cavity_mode_1%lambda == 0")
  IF (error .AND. Debug) WRITE(out_unit,*) 'Mode%lambda failed to initialize'

  error = (Cavity_mode_1%eq_pos < 0)
  CALL Logical_Test(test_cavity_init, error, test2=.FALSE., info="Cavity_mode_1%eq_pos == 0")
  IF (error .AND. Debug) WRITE(out_unit,*) 'Mode%eq_pos failed to initialize'

  CALL Finalize_Test(test_cavity_init)

  !IF (error == 0) THEN
  !  WRITE(out_unit,*) ''
  !  WRITE(out_unit,*) 'Test 1 checked ! The Cavity mode initializes successfully !'
  !ELSE IF (error == 1) THEN
  !  WRITE(out_unit,*) ''
  !  WRITE(out_unit,*) 'Test 1 failed ! The Cavity mode did not initialize successfully...'
  !ELSE
  !  WRITE(out_unit,*) ''
  !  WRITE(out_unit,*) 'Test 1 stopped ! The file did not execute as it should have !'
  !END IF
  
  
END PROGRAM
