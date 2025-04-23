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
PROGRAM Quantum_HO1D
  !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
  USE QDUtil_m
  USE Tests_m
  USE Matter_mode_m
  IMPLICIT NONE


  integer             :: Verbose = 0
  logical             :: Debug   = .TRUE.

  TYPE(Matter_mode_t) :: MatMode
  logical             :: Dense   = .FALSE.

  TYPE(test_t)        :: test_matmode
  logical             :: error_matmode = .FALSE.

  integer             :: i, i_op


  !-----------------------------Test initialization----------------------------
  CALL Initialize_Test(test_matmode, test_name="OUT/test_file_mat_mode")


  !-------------------------Matter mode initialization-------------------------
  CALL Initialize(MatMode, nio=in_unit, Dense=Dense, Verbose=Verbose, Debug=Debug)
  IF (Debug) THEN
    WRITE(out_unit,*)
    WRITE(out_unit,*) "--------------Cavity mode constructed by MolecCav_Read_cavity_mode--------------"
    CALL Write(MatMode)
    WRITE(out_unit,*) "------------End Cavity mode constructed by MolecCav_Read_cavity_mode------------"
  END IF


  !-----------------------------------The tests--------------------------------

  CALL Finalize_Test(test_matmode)

  
END PROGRAM
