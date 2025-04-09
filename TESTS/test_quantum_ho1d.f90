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
PROGRAM test_quantum_ho1d
  !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
  USE QDUtil_m
  USE Tests_m
  USE Cavity_mode_old_m
  USE Elem_op_m
  USE Quantum_HO1D_m
  IMPLICIT NONE


  integer                       :: Verbose = 0
  logical                       :: Debug   = .FALSE.

  TYPE(Quantum_HO1D_t)          :: QHO1D_opt_1_6_1
  TYPE(Quantum_HO1D_t)          :: QHO1D_opt_14_6_1
  TYPE(Quantum_HO1D_t)          :: QHO1D_opt_1_17_1
  TYPE(Quantum_HO1D_t)          :: QHO1D_opt_1_6_7
  TYPE(Quantum_HO1D_t)          :: QHO1D_dense_1_6_1
  TYPE(Quantum_HO1D_t)          :: QHO1D_dense_14_6_1
  TYPE(Quantum_HO1D_t)          :: QHO1D_dense_1_17_1
  TYPE(Quantum_HO1D_t)          :: QHO1D_dense_1_6_7

  real(kind=Rkind), allocatable :: H_HO1D_diag_ana_14_17(:)
  real(kind=Rkind), allocatable :: H_HO1D_dense_ana_1_6(:,:)
  real(kind=Rkind), allocatable :: H_HO1D_dense_ana_14_17(:,:)
  real(kind=Rkind), allocatable :: x_HO1D_band_ana_1_17_1(:,:)
  real(kind=Rkind), allocatable :: x_HO1D_dense_ana_1_6_1(:,:)
  real(kind=Rkind), allocatable :: x_HO1D_dense_ana_14_17_7(:,:)
  real(kind=Rkind), allocatable :: N_HO1D_dense_ana_17(:,:)

  TYPE(test_t)                  :: test_construct
  logical                       :: error_construct = .FALSE.

  integer                       :: i


  !-----------------------------Test initialization----------------------------
  CALL Initialize_Test(test_construct, test_name="OUT/test_file_qntm_ho1d")

  !---------------------------Construct Quantum HO1D to test---------------------------
  CALL Initialize(QHO1D_opt_1_6_1,    6,  ONE,    ONE,   Verbose=Verbose, Debug=Debug)
  CALL Initialize(QHO1D_opt_14_6_1,   6,  14*ONE, ONE,   Verbose=Verbose, Debug=Debug)
  CALL Initialize(QHO1D_opt_1_17_1,   17, ONE,    ONE,   Verbose=Verbose, Debug=Debug)
  CALL Initialize(QHO1D_opt_1_6_7,    6,  ONE,  SEVEN,   Verbose=Verbose, Debug=Debug)
  CALL Initialize(QHO1D_dense_1_6_1,  6,  ONE,    ONE,   Verbose=Verbose, Debug=Debug)
  CALL Initialize(QHO1D_dense_14_6_1, 6,  14*ONE, ONE,   Verbose=Verbose, Debug=Debug)
  CALL Initialize(QHO1D_dense_1_17_1, 17, ONE,    ONE,   Verbose=Verbose, Debug=Debug)
  CALL Initialize(QHO1D_dense_1_6_7,  6,  ONE,  SEVEN,   Verbose=Verbose, Debug=Debug)


  !-----------------------------------sum up-----------------------------------
  CALL Finalize_Test(test_construct)
  

END PROGRAM

