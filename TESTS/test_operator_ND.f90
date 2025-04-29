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
PROGRAM test_operator_ND
  !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
  USE QDUtil_m
  USE Tests_m
  USE Algebra_m
  !USE Operator_ND_m
  IMPLICIT NONE


  integer             :: Verbose = 40
  logical             :: Debug   = .TRUE.

  character(len=39)   :: test
  character(len=20)   :: Ops(3)

  TYPE(test_t)        :: test_opnd
  logical             :: error_opnd = .FALSE.

  integer             :: i, i_op


  !-----------------------------Test initialization----------------------------
  CALL Initialize_Test(test_opnd, test_name="OUT/test_file_opnd")

  test="hamiltonian ,   position ,Identity     "
  WRITE(out_unit,*) test 
  WRITE(out_unit,*) test(4:4)
  WRITE(out_unit,*) (test=='h')
  READ(unit=test, fmt=*) Ops

  WRITE(out_unit,*) Ops 
  WRITE(out_unit,*) Ops(1) 
  WRITE(out_unit,*) Ops(2)
  WRITE(out_unit,*) Ops(3)
  WRITE(out_unit,*) LEN(Ops)
  WRITE(out_unit,*) LEN_trim(Ops)
  WRITE(out_unit,*) LEN(Ops(1))
  WRITE(out_unit,*) LEN_trim(Ops(1))
  WRITE(out_unit,*) LEN(Ops(2))
  WRITE(out_unit,*) LEN_trim(Ops(2))
  WRITE(out_unit,*) LEN(Ops(3))
  WRITE(out_unit,*) LEN_trim(Ops(3))

  !-----------------------------------The tests--------------------------------

  CALL Finalize_Test(test_opnd)

  
END PROGRAM
