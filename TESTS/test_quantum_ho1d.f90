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

  real(kind=Rkind), allocatable :: I_HO1D_dense_ana_17(:)
  real(kind=Rkind), allocatable :: H_HO1D_diag_ana_14_17(:)
  real(kind=Rkind), allocatable :: H_HO1D_dense_ana_14_17(:,:)
  real(kind=Rkind), allocatable :: x_HO1D_band_ana_14_17_7(:,:)
  real(kind=Rkind), allocatable :: x_HO1D_dense_ana_14_17_7(:,:)
  real(kind=Rkind), allocatable :: N_HO1D_diag_ana_17(:)
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


  !-------------------------Construct reference matricies-----------------------
    !-------------------------H matricies-----------------------
  ALLOCATE(H_HO1D_diag_ana_14_17(17))
  ALLOCATE(H_HO1D_dense_ana_14_17(17,17))

  H_HO1D_diag_ana_14_17  = ZERO
  H_HO1D_dense_ana_14_17 = ZERO

  DO i = 1, 17
    H_HO1D_diag_ana_14_17(i)    = 14*(i - ONE + HALF)
    H_HO1D_dense_ana_14_17(i,i) = 14*(i - ONE + HALF)
  END DO

  IF (Debug) THEN
    WRITE(out_unit,*) 
    WRITE(out_unit,*) "H_{Reference 1 (diagonal, w = 14, Nb = 17)}"
    CALL Write_Vec(H_HO1D_diag_ana_14_17,  out_unit, 1, info="H_{Reference 1}")
  
    WRITE(out_unit,*) 
    WRITE(out_unit,*) "H_{Reference 2 (dense, w = 14, Nb = 17)}"
    CALL Write_Mat(H_HO1D_dense_ana_14_17, out_unit, Size(H_HO1D_dense_ana_14_17, dim=2), info="H_{Reference 2}")
  END IF

  !-----------------------x matricies----------------------
  ALLOCATE(x_HO1D_band_ana_14_17_7(17,3))
  ALLOCATE(x_HO1D_dense_ana_14_17_7(17,17))

  x_HO1D_band_ana_14_17_7  = ZERO
  x_HO1D_dense_ana_14_17_7 = ZERO

  DO i = 1, 16                                                                 ! /!\ Fortran counts from 1 to Nb !!! /!\ Nb-1 not to have Band_val(i+1) out of range
    x_HO1D_band_ana_14_17_7(i,1)    = SQRT(REAL(i,kind=Rkind))
    x_HO1D_band_ana_14_17_7(i+1,3)  = SQRT(REAL(i,kind=Rkind))
    x_HO1D_dense_ana_14_17_7(i,i+1) = SQRT(REAL(i,kind=Rkind))
    x_HO1D_dense_ana_14_17_7(i+1,i) = SQRT(REAL(i,kind=Rkind))
  END DO

  x_HO1D_band_ana_14_17_7  = x_HO1D_band_ana_14_17_7  / SQRT(TWO * 14.0_Rkind * 7.0_Rkind)
  x_HO1D_dense_ana_14_17_7 = x_HO1D_dense_ana_14_17_7 / SQRT(TWO * 14.0_Rkind * 7.0_Rkind)

  IF (Debug) THEN
    WRITE(out_unit,*) 
    WRITE(out_unit,*) "x_{Reference 1 (band, w = 1.0, Nb = 17, m = 1.0)}"
    CALL Write_Mat(x_HO1D_band_ana_14_17_7,  out_unit, 3, info="x_{Reference 1}")

    WRITE(out_unit,*) 
    WRITE(out_unit,*) "x_{Reference 2 (dense, w = 14.0, Nb = 17, m = 7.0)}"
    CALL Write_Mat(x_HO1D_dense_ana_14_17_7, out_unit, Size(H_HO1D_dense_ana_14_17, dim=2), info="x_{Reference 2}")
  END IF

  !-----------------------N matricies----------------------
  ALLOCATE(N_HO1D_diag_ana_17(17))
  ALLOCATE(N_HO1D_dense_ana_17(17,17))
  N_HO1D_diag_ana_17  = ZERO
  N_HO1D_dense_ana_17 = ZERO

  DO i = 1, 17                                                             ! /!\ Fortran counts from 1 to Nb !!! /!\
    N_HO1D_diag_ana_17(i) = i - 1
    N_HO1D_dense_ana_17(i,i) = i - 1
  END DO

  IF (Debug) THEN
    WRITE(out_unit,*) 
    WRITE(out_unit,*) "N_{Reference 1 (diagonal, Nb = 17)}"
    CALL Write_Mat(N_HO1D_dense_ana_17, out_unit, Size(N_HO1D_dense_ana_17, dim=2), info="N_{Reference 1}")

    WRITE(out_unit,*) 
    WRITE(out_unit,*) "N_{Reference 2 (dense, Nb = 17)}"
    CALL Write_Mat(N_HO1D_dense_ana_17, out_unit, Size(N_HO1D_dense_ana_17, dim=2), info="N_{Reference 2}")
  END IF
  FLUSH(out_unit)


  !--------------------------------Comparisons H-------------------------------
  !QHO1D_opt_1_17_1%Tab_op(1)%Diag_val(1) = 0
  !WRITE(out_unit,*) 
  !WRITE(out_unit,*) "H_{Spurious 1 (diagonal, w = 14, Nb = 17)}"
  !DO i = 1, 17
  !  WRITE(out_unit,*) H_HO1D_diag_ana_14_17(i)
  !END DO

  CALL Equal_tensor(error_construct, H_HO1D_diag_ana_14_17, 14*QHO1D_opt_1_17_1%Tab_op(1)%Diag_val)
  CALL Logical_Test(test_construct, test1=error_construct, test2=.FALSE., info="QHO1D_opt_1_17_1%H%Diag_val well initialized ?")
  IF (error_construct .AND. Debug) THEN
    CALL Write_Vec(H_HO1D_diag_ana_14_17, out_unit, Size(H_HO1D_diag_ana_14_17), info="H_HO1D_diag_ana_14_17")
    CALL Write_Vec(QHO1D_opt_1_17_1%Tab_op(1)%Diag_val, out_unit, Size(QHO1D_opt_1_17_1%Tab_op(1)%Diag_val, dim=1), in&
                  &fo="QHO1D_opt_1_17_1%Tab_op(1)%Diag_val")
  END IF

  !QHO1D_opt_14_6_1%Tab_op(1)%Diag_val(1,2) = 1
  !WRITE(out_unit,*) 
  !WRITE(out_unit,*) "H_{Spurious 2 (diagonal, w = 1, Nb = 6)}"
  !DO i = 1, 17
  !  WRITE(out_unit,*) QHO1D_opt_14_6_1%Tab_op(1)%Diag_val(i,:)
  !END DO

  CALL Equal_tensor(error_construct, H_HO1D_diag_ana_14_17(1:6), QHO1D_opt_14_6_1%Tab_op(1)%Diag_val)
  CALL Logical_Test(test_construct, test1=error_construct, test2=.FALSE., info="QHO1D_opt_14_6_1%H%Diag_val well initialized ?")
  IF (error_construct .AND. Debug) THEN
    CALL Write_Vec(H_HO1D_diag_ana_14_17, out_unit, Size(H_HO1D_diag_ana_14_17, dim=1), info="H_HO1D_diag_ana_14_17")
    CALL Write_Vec(QHO1D_opt_14_6_1%Tab_op(1)%Diag_val, out_unit, Size(QHO1D_opt_14_6_1%Tab_op(1)%Diag_val, dim=1), in&
                  &fo="QHO1D_opt_14_6_1%Tab_op(1)%Diag_val")
  END IF

  !QHO1D_opt_1_6_1%Tab_op(1)%Diag_val(1,2) = 1
  !WRITE(out_unit,*) 
  !WRITE(out_unit,*) "H_{Spurious 2 (diagonal, w = 1, Nb = 6)}"
  !DO i = 1, 17
  !  WRITE(out_unit,*) QHO1D_opt_1_6_1%Tab_op(1)%Diag_val(i,:)
  !END DO

  CALL Equal_tensor(error_construct, H_HO1D_diag_ana_14_17(1:6), 14*QHO1D_opt_1_6_1%Tab_op(1)%Diag_val)
  CALL Logical_Test(test_construct, test1=error_construct, test2=.FALSE., info="QHO1D_opt_1_6_1%H%Diag_val well initialized ?")
  IF (error_construct .AND. Debug) THEN
    CALL Write_Vec(H_HO1D_diag_ana_14_17, out_unit, Size(H_HO1D_diag_ana_14_17, dim=1), info="H_HO1D_diag_ana_14_17")
    CALL Write_Vec(14*QHO1D_opt_1_6_1%Tab_op(1)%Diag_val, out_unit, Size(QHO1D_opt_1_6_1%Tab_op(1)%Diag_val, dim=1), in&
                  &fo="14*QHO1D_opt_1_6_1%Tab_op(1)%Diag_val")
  END IF

  !QHO1D_dense_1_17_1%Tab_op(1)%Dense_val(1) = 0
  !WRITE(out_unit,*) 
  !WRITE(out_unit,*) "H_{Spurious 1 (diagonal, w = 14, Nb = 17)}"
  !DO i = 1, 17
  !  WRITE(out_unit,*) H_HO1D_dense_ana_14_17(i)
  !END DO

  WRITE(out_unit,*) "H_HO1D_dense_ana_14_17 : "//TO_string(SIZE(H_HO1D_dense_ana_14_17))
  WRITE(out_unit,*) "QHO1D_dense_1_17_1%Tab_op(1)%Dense_val : "//TO_string(SIZE(QHO1D_dense_1_17_1%Tab_op(1)%Dense_val))
  STOP "STOPSTOPSTOP"
  CALL Equal_tensor(error_construct, H_HO1D_dense_ana_14_17, 14*QHO1D_dense_1_17_1%Tab_op(1)%Dense_val)
  CALL Logical_Test(test_construct, test1=error_construct, test2=.FALSE., info="QHO1D_dense_1_17_1%H%Dense_val well initialized ?")
  IF (error_construct .AND. Debug) THEN
    CALL Write_Mat(H_HO1D_dense_ana_14_17, out_unit, Size(H_HO1D_dense_ana_14_17), info="H_HO1D_dense_ana_14_17")
    CALL Write_Mat(QHO1D_dense_1_17_1%Tab_op(1)%Dense_val, out_unit, Size(QHO1D_dense_1_17_1%Tab_op(1)%Dense_val), in&
                  &fo="QHO1D_dense_1_17_1%Tab_op(1)%Dense_val")
  END IF
  
  !QHO1D_dense_14_6_1%Tab_op(1)%Dense_val(1,2) = 1
  !WRITE(out_unit,*) 
  !WRITE(out_unit,*) "H_{Spurious 2 (diagonal, w = 1, Nb = 6)}"
  !DO i = 1, 17
  !  WRITE(out_unit,*) QHO1D_dense_14_6_1%Tab_op(1)%Dense_val(i,:)
  !END DO

  CALL Equal_tensor(error_construct, H_HO1D_dense_ana_14_17(1:6,1:6), QHO1D_dense_14_6_1%Tab_op(1)%Dense_val)
  CALL Logical_Test(test_construct, test1=error_construct, test2=.FALSE., info="QHO1D_dense_14_6_1%H%Dense_val well initialized ?")
  IF (error_construct .AND. Debug) THEN
    CALL Write_Mat(H_HO1D_dense_ana_14_17, out_unit, Size(H_HO1D_dense_ana_14_17, dim=2), info="H_HO1D_dense_ana_14_17")
    CALL Write_Mat(QHO1D_dense_14_6_1%Tab_op(1)%Dense_val, out_unit, Size(QHO1D_dense_14_6_1%Tab_op(1)%Dense_val, dim=2), in&
                  &fo="QHO1D_dense_14_6_1%Tab_op(1)%Dense_val")
  END IF

  !QHO1D_dense_1_6_1%Tab_op(1)%Dense_val(1,2) = 1
  !WRITE(out_unit,*) 
  !WRITE(out_unit,*) "H_{Spurious 2 (diagonal, w = 1, Nb = 6)}"
  !DO i = 1, 17
  !  WRITE(out_unit,*) QHO1D_dense_1_6_1%Tab_op(1)%Dense_val(i,:)
  !END DO

  CALL Equal_tensor(error_construct, H_HO1D_dense_ana_14_17(1:6,1:6), 14*QHO1D_dense_1_6_1%Tab_op(1)%Dense_val)
  CALL Logical_Test(test_construct, test1=error_construct, test2=.FALSE., info="QHO1D_dense_1_6_1%H%Dense_val well initialized ?")
  IF (error_construct .AND. Debug) THEN
    CALL Write_Mat(H_HO1D_dense_ana_14_17, out_unit, Size(H_HO1D_dense_ana_14_17, dim=2), info="H_HO1D_dense_ana_14_17")
    CALL Write_Mat(14*QHO1D_dense_1_6_1%Tab_op(1)%Dense_val, out_unit, Size(QHO1D_dense_1_6_1%Tab_op(1)%Dense_val, dim=2), in&
                  &fo="14*QHO1D_dense_1_6_1%Tab_op(1)%Dense_val")
  END IF


  !-----------------------------------sum up-----------------------------------
  CALL Finalize_Test(test_construct)
  

END PROGRAM

