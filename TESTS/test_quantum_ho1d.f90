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
  USE Algebra_m
  USE Elem_op_m
  USE Quantum_HO1D_m
  IMPLICIT NONE


  integer                       :: Verbose = 0
  logical                       :: Debug   = .TRUE.

  TYPE(Quantum_HO1D_t)          :: QHO1D_non_alloc

  TYPE(Quantum_HO1D_t)          :: QHO1D_opt_1_6_1
  TYPE(Quantum_HO1D_t)          :: QHO1D_opt_14_6_1
  TYPE(Quantum_HO1D_t)          :: QHO1D_opt_1_17_1
  TYPE(Quantum_HO1D_t)          :: QHO1D_opt_1_6_7
  TYPE(Quantum_HO1D_t)          :: QHO1D_dense_1_6_1
  TYPE(Quantum_HO1D_t)          :: QHO1D_dense_14_6_1
  TYPE(Quantum_HO1D_t)          :: QHO1D_dense_1_17_1
  TYPE(Quantum_HO1D_t)          :: QHO1D_dense_1_6_7

  real(kind=Rkind), allocatable :: I_HO1D_diag_ana_17(:)
  real(kind=Rkind), allocatable :: I_HO1D_dense_ana_17(:,:)
  real(kind=Rkind), allocatable :: H_HO1D_diag_ana_14_17(:)
  real(kind=Rkind), allocatable :: H_HO1D_dense_ana_14_17(:,:)
  real(kind=Rkind), allocatable :: x_HO1D_band_ana_14_17_7(:,:)
  real(kind=Rkind), allocatable :: x_HO1D_dense_ana_14_17_7(:,:)
  real(kind=Rkind), allocatable :: N_HO1D_diag_ana_17(:)
  real(kind=Rkind), allocatable :: N_HO1D_dense_ana_17(:,:)

  real(kind=Rkind)              :: Psi_real(6)
  real(kind=Rkind)              :: Op_psi_real_qho1d(6)
  real(kind=Rkind)              :: Op_psi_real_elem_op(6)
  complex(kind=Rkind)           :: Psi_complex(6)
  complex(kind=Rkind)           :: Op_psi_complex_qho1d(6)
  complex(kind=Rkind)           :: Op_psi_complex_elem_op(6)
  
  integer                      :: Param_integer
  real(kind=Rkind)             :: Param_real

  TYPE(test_t)                  :: test_qho1d
  logical                       :: error_qho1d = .FALSE.

  integer                       :: i, i_op


  !-----------------------------Test initialization----------------------------
  CALL Initialize_Test(test_qho1d, test_name="OUT/test_file_qntm_ho1d")


  !---------------------------Construct Quantum HO1D to test---------------------------
  CALL Initialize(QHO1D_opt_1_6_1,    6,  ONE,    ONE,                 Verbose=Verbose, Debug=Debug)
  CALL Initialize(QHO1D_opt_14_6_1,   6,  14*ONE, ONE,                 Verbose=Verbose, Debug=Debug)
  CALL Initialize(QHO1D_opt_1_17_1,   17, ONE,    ONE,                 Verbose=Verbose, Debug=Debug)
  CALL Initialize(QHO1D_opt_1_6_7,    6,  ONE,  SEVEN,                 Verbose=Verbose, Debug=Debug)
  CALL Initialize(QHO1D_dense_1_6_1,  6,  ONE,    ONE,   Dense=.TRUE., Verbose=Verbose, Debug=Debug)
  CALL Initialize(QHO1D_dense_14_6_1, 6,  14*ONE, ONE,   Dense=.TRUE., Verbose=Verbose, Debug=Debug)
  CALL Initialize(QHO1D_dense_1_17_1, 17, ONE,    ONE,   Dense=.TRUE., Verbose=Verbose, Debug=Debug)
  CALL Initialize(QHO1D_dense_1_6_7,  6,  ONE,  SEVEN,   Dense=.TRUE., Verbose=Verbose, Debug=Debug)

    !-----------------------Construct Wavepackets to take action on----------------------
  Psi_real    = [ONE, SQRT(TWO), PI, SQRT(TWO), ONE, HALF]
  CALL Normalize(Psi_real)
  Psi_complex = [ONE, SQRT(TWO), PI, SQRT(TWO), ONE, HALF]
  Psi_complex(1) = Psi_complex(2) + Psi_complex(1)*EYE
  Psi_complex(2) = Psi_complex(2)*EYE
  Psi_complex(3) = Psi_complex(1) + Psi_complex(3)*EYE
  Psi_complex(5) = Psi_complex(5)*EYE
  CALL Normalize(Psi_complex)
  
  IF (Debug) THEN
    WRITE(out_unit,*)
    CALL Write_Vec(Psi_real,    out_unit, 1, info="Psi_real")
    WRITE(out_unit,*)
    CALL Write_Vec(Psi_complex, out_unit, 1, info="Psi_complex")
  END IF

  !-------------------------Construct reference matricies-----------------------
    !-------------------------I matricies-----------------------
  ALLOCATE(I_HO1D_diag_ana_17(1))
  ALLOCATE(I_HO1D_dense_ana_17(17,17))

  I_HO1D_diag_ana_17  = ONE
  I_HO1D_dense_ana_17 = ZERO

  DO i = 1, 17
    I_HO1D_dense_ana_17(i,i) = ONE
  END DO

  IF (Debug) THEN
    WRITE(out_unit,*) 
    WRITE(out_unit,*) "I_{Reference 1 (diagonal, Nb = 17)}"
    CALL Write_Vec(I_HO1D_diag_ana_17,  out_unit, 1, info="I_{Reference 1}")
  
    WRITE(out_unit,*) 
    WRITE(out_unit,*) "I_{Reference 2 (dense, Nb = 17)}"
    CALL Write_Mat(I_HO1D_dense_ana_17, out_unit, Size(I_HO1D_dense_ana_17, dim=2), info="I_{Reference 2}")
  END IF

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

  DO i = 1, 16                                                                 !/!\ Fortran counts from 1 to Nb !!!/!\ Nb-1 not to have Band_val(i+1) out of range
    x_HO1D_band_ana_14_17_7(i,1)    = SQRT(REAL(i,kind=Rkind))
    x_HO1D_band_ana_14_17_7(i+1,3)  = SQRT(REAL(i,kind=Rkind))
    x_HO1D_dense_ana_14_17_7(i,i+1) = SQRT(REAL(i,kind=Rkind))
    x_HO1D_dense_ana_14_17_7(i+1,i) = SQRT(REAL(i,kind=Rkind))
  END DO

  x_HO1D_band_ana_14_17_7  = x_HO1D_band_ana_14_17_7 /SQRT(TWO * 14.0_Rkind * 7.0_Rkind)
  x_HO1D_dense_ana_14_17_7 = x_HO1D_dense_ana_14_17_7/SQRT(TWO * 14.0_Rkind * 7.0_Rkind)

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

  DO i = 1, 17                                                             !/!\ Fortran counts from 1 to Nb !!!/!\
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


  !-------------------------------Comparisons I-------------------------------
  !N_HO1D_dense_17%Dense_val(1,1) = 5

  CALL Equal_tensor(error_qho1d, I_HO1D_diag_ana_17, QHO1D_opt_1_17_1%Tab_op(0)%Diag_val)
  CALL Logical_Test(test_qho1d, test1=error_qho1d, test2=.FALSE., info="QHO1D_opt_1_17_1%Tab_op(0)%Diag_val ?")
  IF (error_qho1d .AND. Debug) THEN
    CALL Write_Vec(N_HO1D_diag_ana_17, out_unit, Size(N_HO1D_diag_ana_17, dim=1), info="I_HO1D_diag_ana_17")
    CALL Write_Vec(QHO1D_opt_1_17_1%Tab_op(3)%Diag_val, out_unit, Size(QHO1D_opt_1_17_1%Tab_op(3)%Diag_val, dim=&
                  &1), info="QHO1D_opt_1_17_1%Tab_op(3)%Diag_val")
  END IF

  !N_HO1D_dense_17%Dense_val(1,1) = 5

  CALL Equal_tensor(error_qho1d, I_HO1D_dense_ana_17, QHO1D_dense_1_17_1%Tab_op(0)%Dense_val)
  CALL Logical_Test(test_qho1d, test1=error_qho1d, test2=.FALSE., info="I_HO1D_dense_17%De&
                   &nse_val well initialized ?")
  IF (error_qho1d .AND. Debug) THEN
    CALL Write_Vec(N_HO1D_diag_ana_17, out_unit, Size(N_HO1D_diag_ana_17, dim=1), info="I_HO1D_diag_ana_17")
    CALL Write_Vec(QHO1D_opt_1_6_1%Tab_op(3)%Diag_val, out_unit, Size(QHO1D_opt_1_6_1%Tab_op(3)%Diag_val, dim=&
                  &1), info="QHO1D_opt_1_17_1%Tab_op(1)%Diag_val")
  END IF

  !N_HO1D_diag_17%Diag_val(1) = 5
  !N_HO1D_dense_6%Dense_val(1,1) = 5
  !N_HO1D_diag_6%Diag_val(1) = 5
  
  CALL Equal_tensor(error_qho1d, I_HO1D_dense_ana_17(1:6,1:6), QHO1D_dense_1_6_1%Tab_op(0)%Dense_val)
  CALL Logical_Test(test_qho1d, test1=error_qho1d, test2=.FALSE., info="QHO1D_opt_1_17_1%I%Diag_val well initialized ?")
  IF (error_qho1d .AND. Debug) THEN
    CALL Write_Mat(N_HO1D_dense_ana_17(1:6,1:6), out_unit, Size(N_HO1D_dense_ana_17(1:6,1:6), d&
                  &im=2), info="I_HO1D_dense_ana_17(1:6,1:6)")
    CALL Write_Mat(QHO1D_dense_1_17_1%Tab_op(3)%Dense_val, out_unit, Size(QHO1D_dense_1_17_1%Tab_op(3)%Dense_val, dim=2),&
                 & info="_densQHO1D_opt_1_17_1%Tab_op(1)%Diag_val")
  END IF

  !N_HO1D_diag_17%Diag_val(1) = 5
  !N_HO1D_dense_6%Dense_val(1,1) = 5
  !N_HO1D_diag_6%Diag_val(1) = 5
  
  CALL Equal_tensor(error_qho1d, N_HO1D_dense_ana_17(1:6,1:6), QHO1D_dense_1_6_1%Tab_op(3)%Dense_val)
  CALL Logical_Test(test_qho1d, test1=error_qho1d, test2=.FALSE., info="N_HO1D_dense_6%Dense_val well initialized ?")
  IF (error_qho1d .AND. Debug) THEN
    CALL Write_Mat(N_HO1D_dense_ana_17(1:6,1:6), out_unit, Size(N_HO1D_dense_ana_17(1:6,1:6), d&
                  &im=2), info="N_HO1D_dense_ana_17(1:6,1:6)")
    CALL Write_Mat(QHO1D_dense_1_6_1%Tab_op(3)%Dense_val, out_unit, Size(QHO1D_dense_1_6_1%Tab_op(3)%Dense_val, dim=2),&
                 & info="QHO1D_dense_1_6_1%Tab_op(3)%Dense_val")
  END IF


  !--------------------------------Comparisons H-------------------------------
  !QHO1D_opt_1_17_1%Tab_op(1)%Diag_val(1) = 0
  !WRITE(out_unit,*) 
  !WRITE(out_unit,*) "H_{Spurious 1 (diagonal, w = 14, Nb = 17)}"
  !DO i = 1, 17
  !  WRITE(out_unit,*) H_HO1D_diag_ana_14_17(i)
  !END DO

  CALL Equal_tensor(error_qho1d, H_HO1D_diag_ana_14_17, 14*QHO1D_opt_1_17_1%Tab_op(1)%Diag_val)
  CALL Logical_Test(test_qho1d, test1=error_qho1d, test2=.FALSE., info="QHO1D_opt_1_17_1%H%Diag_val well initialized ?")
  IF (error_qho1d .AND. Debug) THEN
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

  CALL Equal_tensor(error_qho1d, H_HO1D_diag_ana_14_17(1:6), QHO1D_opt_14_6_1%Tab_op(1)%Diag_val)
  CALL Logical_Test(test_qho1d, test1=error_qho1d, test2=.FALSE., info="QHO1D_opt_14_6_1%H%Diag_val well initialized ?")
  IF (error_qho1d .AND. Debug) THEN
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

  CALL Equal_tensor(error_qho1d, H_HO1D_diag_ana_14_17(1:6), 14*QHO1D_opt_1_6_1%Tab_op(1)%Diag_val)
  CALL Logical_Test(test_qho1d, test1=error_qho1d, test2=.FALSE., info="QHO1D_opt_1_6_1%H%Diag_val well initialized ?")
  IF (error_qho1d .AND. Debug) THEN
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
  CALL Equal_tensor(error_qho1d, H_HO1D_dense_ana_14_17, 14*QHO1D_dense_1_17_1%Tab_op(1)%Dense_val)
  CALL Logical_Test(test_qho1d, test1=error_qho1d, test2=.FALSE., info="QHO1D_dense_1_17_1%H%Dense_val well initialized ?")
  IF (error_qho1d .AND. Debug) THEN
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

  CALL Equal_tensor(error_qho1d, H_HO1D_dense_ana_14_17(1:6,1:6), QHO1D_dense_14_6_1%Tab_op(1)%Dense_val)
  CALL Logical_Test(test_qho1d, test1=error_qho1d, test2=.FALSE., info="QHO1D_dense_14_6_1%H%Dense_val well initialized ?")
  IF (error_qho1d .AND. Debug) THEN
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

  CALL Equal_tensor(error_qho1d, H_HO1D_dense_ana_14_17(1:6,1:6), 14*QHO1D_dense_1_6_1%Tab_op(1)%Dense_val)
  CALL Logical_Test(test_qho1d, test1=error_qho1d, test2=.FALSE., info="QHO1D_dense_1_6_1%H%Dense_val well initialized ?")
  IF (error_qho1d .AND. Debug) THEN
    CALL Write_Mat(H_HO1D_dense_ana_14_17, out_unit, Size(H_HO1D_dense_ana_14_17, dim=2), info="H_HO1D_dense_ana_14_17")
    CALL Write_Mat(14*QHO1D_dense_1_6_1%Tab_op(1)%Dense_val, out_unit, Size(QHO1D_dense_1_6_1%Tab_op(1)%Dense_val, dim=2), in&
                  &fo="14*QHO1D_dense_1_6_1%Tab_op(1)%Dense_val")
  END IF


  !-------------------------------Comparisons x-------------------------------
  !QHO1D_opt_1_17_1%Tab_op(2)%Band_val(6,1) = 15
  !WRITE(out_unit,*) "x_{Spurious 4 (Band, w = 1, Nb = 6, m = 1)}"
  !DO i = 1, 6
  !  WRITE(out_unit,*) QHO1D_opt_1_17_1%Tab_op(2)%Band_val(i,:)
  !END DO

  CALL Equal_tensor(error_qho1d, x_HO1D_band_ana_14_17_7, QHO1D_opt_1_17_1%Tab_op(2)%Band_val/SQRT(14*SEVEN))
  CALL Logical_Test(test_qho1d, test1=error_qho1d, test2=.FALSE., info="x_HO1D_band_1_17_1&
                   &%Band_val well initialized ?")
  IF (error_qho1d .AND. Debug) THEN
    CALL Write_Mat(x_HO1D_band_ana_14_17_7, out_unit, Size(x_HO1D_band_ana_14_17_7, dim=2), info="x_HO1D_band_ana_14_17_7")
    CALL Write_Mat(QHO1D_opt_1_17_1%Tab_op(2)%Band_val/SQRT(14*SEVEN), out_unit, Size(QHO1D_opt_1_17_1%Tab_op(2)%Band_val, dim=&
                 &2), info="QHO1D_opt_1_17_1%Tab_op(2)%Band_val/SQRT(14*SEVEN)")
  END IF

  !QHO1D_opt_1_6_1%Tab_op(2)%Band_val(1,1) = 5
  !WRITE(out_unit,*)
  !WRITE(out_unit,*) "x_{Spurious 2 (dense, w = 1, Nb = 6, m = 1)}"
  !DO i = 1, 6
  !  WRITE(out_unit,*) QHO1D_opt_1_6_1%Tab_op(2)%Band_val(i,:)
  !END DO

  CALL Equal_tensor(error_qho1d, x_HO1D_band_ana_14_17_7(1:5,:), QHO1D_opt_1_6_1%Tab_op(2)%Band_val(1:5,:)/SQRT(14*SEVEN)) !/!\ have to truncate at 5 because on the 6-vectors basis set, the tridiagonal matrix have only 5 values on the lower/upper diagonals
  CALL Logical_Test(test_qho1d, test1=error_qho1d, test2=.FALSE., info="QHO1D_opt_1_17_1%Tab_op(2)%Band_val/SQRT(14*SEV&
                   &EN) well initialized ?")
  IF (error_qho1d .AND. Debug) THEN
    CALL Write_Mat(x_HO1D_band_ana_14_17_7, out_unit, Size(x_HO1D_band_ana_14_17_7, dim=2), info="x_HO1D_band_ana_14_17_7")
    CALL Write_Mat(QHO1D_opt_1_6_1%Tab_op(2)%Band_val/SQRT(14*SEVEN), out_unit, Size(QHO1D_opt_1_6_1%Tab_op(2)%Band_val/SQRT(&
                  &14*SEVEN), dim=2), info="QHO1D_opt_1_6_1%Tab_op(2)%Band_val/SQRT(14*SEVEN)")
  END IF

  !QHO1D_opt_1_6_1%Tab_op(2)%Band_val/SQRT(14*SEVEN)(1,1) = 5
  !WRITE(out_unit,*)
  !WRITE(out_unit,*) "x_{Spurious 3 (dense, w = 14, Nb = 17, m = 7)}"
  !DO i = 1, 17
  !  WRITE(out_unit,*) QHO1D_opt_1_6_1%Tab_op(2)%Band_val/SQRT(14*SEVEN)(i,:)
  !END DO

  CALL Equal_tensor(error_qho1d, x_HO1D_band_ana_14_17_7(1:5,:), QHO1D_opt_14_6_1%Tab_op(2)%Band_val(1:5,:)/SQRT(SEVEN))
  CALL Logical_Test(test_qho1d, test1=error_qho1d, test2=.FALSE., info="QHO1D_opt_14_6_1%x(2)%Band_val well initialized ?")
  IF (error_qho1d .AND. Debug) THEN
    CALL Write_Mat(x_HO1D_band_ana_14_17_7, out_unit, Size(x_HO1D_band_ana_14_17_7, dim=2), info="x_HO1D_band_ana_14_17_7")
    CALL Write_Mat(QHO1D_opt_14_6_1%Tab_op(2)%Band_val/SQRT(SEVEN), out_unit, Size(QHO1D_opt_14_6_1%Tab_op(2)%Band_val/SQRT(S&
                  &EVEN), dim=2), info="QHO1D_opt_14_6_1%Tab_op(2)%Band_val/SQRT(14*SEVEN)")
  END IF

  !QHO1D_opt_1_6_7%Tab_op(2)%Band_val/SQRT(14*ONE)(1,1) = 5
  !WRITE(out_unit,*)
  !WRITE(out_unit,*) "x_{Spurious 3 (dense, w = 14, Nb = 17, m = 7)}"
  !DO i = 1, 17
  !  WRITE(out_unit,*) QHO1D_opt_1_6_7%Tab_op(2)%Band_val/SQRT(14*ONE)(i,:)
  !END DO

  CALL Equal_tensor(error_qho1d, x_HO1D_band_ana_14_17_7(1:5,:), QHO1D_opt_1_6_7%Tab_op(2)%Band_val(1:5,:)/SQRT(14*ONE))
  CALL Logical_Test(test_qho1d, test1=error_qho1d, test2=.FALSE., info="x_HO1D_dense_14_17&
                  &_7%Dense_val well initialized ?")
  IF (error_qho1d .AND. Debug) THEN
    CALL Write_Mat(x_HO1D_band_ana_14_17_7, out_unit, Size(x_HO1D_band_ana_14_17_7, dim=2), info="x_HO1D_band_ana_14_17_7")
    CALL Write_Mat(QHO1D_opt_1_6_7%Tab_op(2)%Band_val/SQRT(14*ONE), out_unit, Size(QHO1D_opt_1_6_7%Tab_op(2)%Band_val/SQRT(14&
                  &*ONE), dim=2), info="QHO1D_opt_1_6_7%Tab_op(2)%Band_val/SQRT(14*ONE)")
  END IF


  !QHO1D_dense_1_17_1%Tab_op(2)%Dense_val(6,1) = 15
  !WRITE(out_unit,*) "x_{Spurious 4 (Band, w = 1, Nb = 6, m = 1)}"
  !DO i = 1, 6
  !  WRITE(out_unit,*) QHO1D_dense_1_17_1%Tab_op(2)%Dense_val(i,:)
  !END DO

  CALL Equal_tensor(error_qho1d, x_HO1D_dense_ana_14_17_7, QHO1D_dense_1_17_1%Tab_op(2)%Dense_val/SQRT(14*SEVEN))
  CALL Logical_Test(test_qho1d, test1=error_qho1d, test2=.FALSE., info="x_HO1D_band_1_17_1&
                   &%Dense_val well initialized ?")
  IF (error_qho1d .AND. Debug) THEN
    CALL Write_Mat(x_HO1D_dense_ana_14_17_7, out_unit, Size(x_HO1D_dense_ana_14_17_7, dim=2), info="x_HO1D_dense_ana_14_17_7")
    CALL Write_Mat(QHO1D_dense_1_17_1%Tab_op(2)%Dense_val/SQRT(14*SEVEN), out_unit, Size(QHO1D_dense_1_17_1%Tab_op(2)%Dense_val&
                 &, dim=2), info="QHO1D_dense_1_17_1%Tab_op(2)%Dense_val/SQRT(14*SEVEN)")
  END IF

  !QHO1D_dense_1_6_1%Tab_op(2)%Dense_val(1,1) = 5
  !WRITE(out_unit,*)
  !WRITE(out_unit,*) "x_{Spurious 2 (dense, w = 1, Nb = 6, m = 1)}"
  !DO i = 1, 6
  !  WRITE(out_unit,*) QHO1D_dense_1_6_1%Tab_op(2)%Dense_val(i,:)
  !END DO

  CALL Equal_tensor(error_qho1d,x_HO1D_dense_ana_14_17_7(1:5,1:5),QHO1D_dense_1_6_1%Tab_op(2)%Dense_val(1:5,1:5)/SQRT(14*SEVEN)) !/!\ have to truncate at 5 because on the 6-vectors basis set, the tridiagonal matrix have only 5 values on the lower/upper diagonals
  CALL Logical_Test(test_qho1d, test1=error_qho1d, test2=.FALSE., info="QHO1D_dense_1_17_1%Tab_op(2)%Dense_val/SQRT(14*&
                   &SEVEN) well initialized ?")
  IF (error_qho1d .AND. Debug) THEN
    CALL Write_Mat(x_HO1D_dense_ana_14_17_7, out_unit, Size(x_HO1D_dense_ana_14_17_7, dim=2), info="x_HO1D_dense_ana_14_17_7")
    CALL Write_Mat(QHO1D_dense_1_6_1%Tab_op(2)%Dense_val/SQRT(14*SEVEN), out_unit, Size(QHO1D_dense_1_6_1%Tab_op(2)%Dense_val/&
                 & SQRT(14*SEVEN), dim=2), info="QHO1D_dense_1_6_1%Tab_op(2)%Dense_val/SQRT(14*SEVEN)")
  END IF

  !QHO1D_dense_1_6_1%Tab_op(2)%Dense_val/SQRT(14*SEVEN)(1,1) = 5
  !WRITE(out_unit,*)
  !WRITE(out_unit,*) "x_{Spurious 3 (dense, w = 14, Nb = 17, m = 7)}"
  !DO i = 1, 17
  !  WRITE(out_unit,*) QHO1D_dense_1_6_1%Tab_op(2)%Dense_val/SQRT(14*SEVEN)(i,:)
  !END DO

  CALL Equal_tensor(error_qho1d, x_HO1D_dense_ana_14_17_7(1:5,1:5), QHO1D_dense_14_6_1%Tab_op(2)%Dense_val(1:5,1:5)/SQRT(SEVEN))
  CALL Logical_Test(test_qho1d, test1=error_qho1d, test2=.FALSE., info="QHO1D_dense_14_6_1%x%Dense_val well initialized ?")
  IF (error_qho1d .AND. Debug) THEN
    CALL Write_Mat(x_HO1D_dense_ana_14_17_7, out_unit, Size(x_HO1D_dense_ana_14_17_7, dim=2), info="x_HO1D_dense_ana_14_17_7")
    CALL Write_Mat(QHO1D_dense_14_6_1%Tab_op(2)%Dense_val/SQRT(SEVEN), out_unit, Size(QHO1D_dense_14_6_1%Tab_op(2)%Dense_val/&
                  &SQRT(SEVEN), dim=2), info="QHO1D_dense_14_6_1%Tab_op(2)%Dense_val/SQRT(14*SEVEN)")
  END IF

  !QHO1D_dense_1_6_7%Tab_op(2)%Dense_val/SQRT(14*ONE)(1,1) = 5
  !WRITE(out_unit,*)
  !WRITE(out_unit,*) "x_{Spurious 3 (dense, w = 14, Nb = 17, m = 7)}"
  !DO i = 1, 17
  !  WRITE(out_unit,*) QHO1D_dense_1_6_7%Tab_op(2)%Dense_val/SQRT(14*ONE)(i,:)
  !END DO

  CALL Equal_tensor(error_qho1d, x_HO1D_dense_ana_14_17_7(1:5,1:5), QHO1D_dense_1_6_7%Tab_op(2)%Dense_val(1:5,1:5)/SQRT(14*ONE))
  CALL Logical_Test(test_qho1d, test1=error_qho1d, test2=.FALSE., info="x_HO1D_dense_14_17&
                  &_7%Dense_val well initialized ?")
  IF (error_qho1d .AND. Debug) THEN
    CALL Write_Mat(x_HO1D_dense_ana_14_17_7, out_unit, Size(x_HO1D_dense_ana_14_17_7, dim=2), info="x_HO1D_dense_ana_14_17_7")
    CALL Write_Mat(QHO1D_dense_1_6_7%Tab_op(2)%Dense_val/SQRT(14*ONE), out_unit, Size(QHO1D_dense_1_6_7%Tab_op(2)%Dense_val/S&
                  &QRT(14*ONE), dim=2), info="QHO1D_dense_1_6_7%Tab_op(2)%Dense_val/SQRT(14*ONE)")
  END IF


  !-------------------------------Comparisons N-------------------------------
  !N_HO1D_dense_17%Dense_val(1,1) = 5

  CALL Equal_tensor(error_qho1d, N_HO1D_diag_ana_17, QHO1D_opt_1_17_1%Tab_op(3)%Diag_val)
  CALL Logical_Test(test_qho1d, test1=error_qho1d, test2=.FALSE., info="N_HO1D_dense_17%De&
                   &nse_val well initialized ?")
  IF (error_qho1d .AND. Debug) THEN
    CALL Write_Vec(N_HO1D_diag_ana_17, out_unit, Size(N_HO1D_diag_ana_17, dim=1), info="N_HO1D_diag_ana_17")
    CALL Write_Vec(QHO1D_opt_1_17_1%Tab_op(3)%Diag_val, out_unit, Size(QHO1D_opt_1_17_1%Tab_op(3)%Diag_val, dim=&
                  &1), info="QHO1D_opt_1_17_1%Tab_op(3)%Diag_val")
  END IF

  !N_HO1D_dense_17%Dense_val(1,1) = 5

  CALL Equal_tensor(error_qho1d, N_HO1D_diag_ana_17(1:6), QHO1D_opt_1_6_1%Tab_op(3)%Diag_val)
  CALL Logical_Test(test_qho1d, test1=error_qho1d, test2=.FALSE., info="N_HO1D_dense_17%De&
                   &nse_val well initialized ?")
  IF (error_qho1d .AND. Debug) THEN
    CALL Write_Vec(N_HO1D_diag_ana_17, out_unit, Size(N_HO1D_diag_ana_17, dim=1), info="N_HO1D_diag_ana_17")
    CALL Write_Vec(QHO1D_opt_1_6_1%Tab_op(3)%Diag_val, out_unit, Size(QHO1D_opt_1_6_1%Tab_op(3)%Diag_val, dim=&
                  &1), info="QHO1D_opt_1_6_1%Tab_op(3)%Diag_val")
  END IF

  !N_HO1D_diag_17%Diag_val(1) = 5
  !N_HO1D_dense_6%Dense_val(1,1) = 5
  !N_HO1D_diag_6%Diag_val(1) = 5
  
  CALL Equal_tensor(error_qho1d, N_HO1D_dense_ana_17, QHO1D_dense_1_17_1%Tab_op(3)%Dense_val)
  CALL Logical_Test(test_qho1d, test1=error_qho1d, test2=.FALSE., info="N_HO1D_dense_6%Dense_val well initialized ?")
  IF (error_qho1d .AND. Debug) THEN
    CALL Write_Mat(N_HO1D_dense_ana_17(1:6,1:6), out_unit, Size(N_HO1D_dense_ana_17(1:6,1:6), d&
                  &im=2), info="N_HO1D_dense_ana_17(1:6,1:6)")
    CALL Write_Mat(QHO1D_dense_1_17_1%Tab_op(3)%Dense_val, out_unit, Size(QHO1D_dense_1_17_1%Tab_op(3)%Dense_val, dim=2),&
                 & info="_dense_1_17_1%Tab_op(3)%Dense_val")
  END IF

  !N_HO1D_diag_17%Diag_val(1) = 5
  !N_HO1D_dense_6%Dense_val(1,1) = 5
  !N_HO1D_diag_6%Diag_val(1) = 5
  
  CALL Equal_tensor(error_qho1d, N_HO1D_dense_ana_17(1:6,1:6), QHO1D_dense_1_6_1%Tab_op(3)%Dense_val)
  CALL Logical_Test(test_qho1d, test1=error_qho1d, test2=.FALSE., info="N_HO1D_dense_6%Dense_val well initialized ?")
  IF (error_qho1d .AND. Debug) THEN
    CALL Write_Mat(N_HO1D_dense_ana_17(1:6,1:6), out_unit, Size(N_HO1D_dense_ana_17(1:6,1:6), d&
                  &im=2), info="N_HO1D_dense_ana_17(1:6,1:6)")
    CALL Write_Mat(QHO1D_dense_1_6_1%Tab_op(3)%Dense_val, out_unit, Size(QHO1D_dense_1_6_1%Tab_op(3)%Dense_val, dim=2),&
                 & info="QHO1D_dense_1_6_1%Tab_op(3)%Dense_val")
  END IF


  !-------------------------Action of a Quantum_HO1D object-----------------------
  CALL Action(Op_psi_real_qho1d,     QHO1D_opt_14_6_1, 0,           Psi_real, Verbose=Verbose, Debug=Debug)
  CALL Equal_tensor(error_qho1d, Op_psi_real_qho1d, Psi_real)
  CALL Logical_Test(test_qho1d,   test1=error_qho1d, test2=.FALSE., info="QHO1D"//TO_string(0)//"^{th} operator action on&
    & a real rank-one tensor")
  DO i_op = 1, SIZE(QHO1D_dense_14_6_1%Tab_op)-1
    CALL Action(Op_psi_real_qho1d,   QHO1D_opt_14_6_1, i_op,        Psi_real, Verbose=Verbose, Debug=Debug)
    CALL Action(Op_psi_real_elem_op, QHO1D_opt_14_6_1%Tab_op(i_op), Psi_real, Verbose=Verbose, Debug=Debug)
    CALL Equal_tensor(error_qho1d, Op_psi_real_qho1d, Op_psi_real_elem_op)
    CALL Logical_Test(test_qho1d, test1=error_qho1d, test2=.FALSE., info="QHO1D"//TO_string(i_op)//"^{th} operator action&
    & on a real rank-one tensor")
  END DO

  CALL Action(Op_psi_complex_qho1d,     QHO1D_opt_14_6_1, 0,           Psi_complex, Verbose=Verbose, Debug=Debug)
  CALL Equal_tensor(error_qho1d, Op_psi_complex_qho1d, Psi_complex)
  CALL Logical_Test(test_qho1d,   test1=error_qho1d, test2=.FALSE., info="QHO1D"//TO_string(0)//"^{th} operator action on&
    & a complex rank-one tensor")
  DO i_op = 1, SIZE(QHO1D_dense_14_6_1%Tab_op)-1
    CALL Action(Op_psi_complex_qho1d,   QHO1D_opt_14_6_1, i_op,        Psi_complex, Verbose=Verbose, Debug=Debug)
    CALL Action(Op_psi_complex_elem_op, QHO1D_opt_14_6_1%Tab_op(i_op), Psi_complex, Verbose=Verbose, Debug=Debug)
    CALL Equal_tensor(error_qho1d, Op_psi_complex_qho1d, Op_psi_complex_elem_op)
    CALL Logical_Test(test_qho1d, test1=error_qho1d, test2=.FALSE., info="QHO1D"//TO_string(i_op)//"^{th} operator action&
    & on a complex rank-one tensor")
  END DO


  !-------------------------Get a parameter from a Quantum_HO1D object-----------------------
  CALL Get(Param_integer,      QHO1D_opt_14_6_1, "Nb")
  CALL Logical_Test(test_qho1d, (Param_integer/=QHO1D_opt_14_6_1%Nb),test2=.FALSE.,info="Get Nb ?")
  CALL Get(Param_real,       QHO1D_opt_14_6_1, "w")
  CALL Logical_Test(test_qho1d, (Param_real/=QHO1D_opt_14_6_1%w),test2=.FALSE.,info="Get w ?")
  CALL Get(Param_real,       QHO1D_opt_14_6_1, "m")
  CALL Logical_Test(test_qho1d, (Param_real/=QHO1D_opt_14_6_1%m),test2=.FALSE.,info="Get m ?")
  CALL Get(Param_integer,   QHO1D_opt_14_6_1, "Nb_op")
  CALL Logical_Test(test_qho1d, (Param_integer/=QHO1D_opt_14_6_1%Nb_op),test2=.FALSE.,info="Get Nb_op ?")
  CALL Get(Param_integer,      QHO1D_opt_14_6_1, "Nq")
  CALL Logical_Test(test_qho1d, (Param_integer/=QHO1D_opt_14_6_1%Nq),test2=.FALSE.,info="Get Nq ?")
  CALL Get(Param_real,  QHO1D_opt_14_6_1, "Eq_pos")
  CALL Logical_Test(test_qho1d, (Param_real/=QHO1D_opt_14_6_1%Eq_pos),test2=.FALSE.,info="Get Eq_pos ?")
  CALL Get(Param_real, QHO1D_opt_14_6_1, "Scale_q")
  CALL Logical_Test(test_qho1d, (Param_real/=QHO1D_opt_14_6_1%Scale_q),test2=.FALSE.,info="Get Scale_q ?")


  !-------------------------Deallocate Quantum_HO1D object-----------------------
  CALL Dealloc(QHO1D_opt_14_6_1, Verbose=Verbose, Debug=Debug)
  CALL Logical_Test(test_qho1d, (QHO1D_opt_14_6_1%Nb/=QHO1D_non_alloc%Nb),test2=.FALSE.,info="QHO1D%Nq well deallocated ?")
  CALL Equal_tensor(error_qho1d, QHO1D_opt_14_6_1%w, QHO1D_non_alloc%w)
  CALL Logical_Test(test_qho1d, test1=error_qho1d,                    test2=.FALSE.,info="QHO1D%w well deallocated ?")
  CALL Equal_tensor(error_qho1d, QHO1D_opt_14_6_1%m, QHO1D_non_alloc%m)
  CALL Logical_Test(test_qho1d, test1=error_qho1d,                    test2=.FALSE.,info="QHO1D%m well deallocated ?")
  CALL Logical_Test(test_qho1d, test1=ALLOCATED(QHO1D_opt_14_6_1%Tab_op), test2=.FALSE.,info="Tab_op well deallocated ?")
  CALL Logical_Test(test_qho1d, (QHO1D_opt_14_6_1%Nq/=QHO1D_non_alloc%Nq),test2=.FALSE.,info="QHO1D%Nq well deallocated ?")
  CALL Equal_tensor(error_qho1d, QHO1D_opt_14_6_1%Eq_pos, QHO1D_non_alloc%Eq_pos)
  CALL Logical_Test(test_qho1d, test1=error_qho1d,                    test2=.FALSE.,info="QHO1D%Eq_pos well deallocated ?")
  CALL Equal_tensor(error_qho1d, QHO1D_opt_14_6_1%Scale_q, QHO1D_non_alloc%Scale_q)
  CALL Logical_Test(test_qho1d, test1=error_qho1d,                    test2=.FALSE.,info="QHO1D%Scale_q well deallocated ?")

  CALL Dealloc(QHO1D_dense_1_17_1, Verbose=Verbose, Debug=Debug)
  CALL Logical_Test(test_qho1d, (QHO1D_dense_1_17_1%Nb/=QHO1D_non_alloc%Nb),test2=.FALSE.,info="QHO1D%Nq well deallocated ?")
  CALL Equal_tensor(error_qho1d, QHO1D_dense_1_17_1%w, QHO1D_non_alloc%w)
  CALL Logical_Test(test_qho1d, test1=error_qho1d,                      test2=.FALSE.,info="QHO1D%w well deallocated ?")
  CALL Equal_tensor(error_qho1d, QHO1D_dense_1_17_1%m, QHO1D_non_alloc%m)
  CALL Logical_Test(test_qho1d, test1=error_qho1d,                      test2=.FALSE.,info="QHO1D%m well deallocated ?")
  CALL Logical_Test(test_qho1d, test1=ALLOCATED(QHO1D_dense_1_17_1%Tab_op), test2=.FALSE.,info="Tab_op well deallocated ?")
  CALL Logical_Test(test_qho1d, (QHO1D_dense_1_17_1%Nq/=QHO1D_non_alloc%Nq),test2=.FALSE.,info="QHO1D%Nq well deallocated ?")
  CALL Equal_tensor(error_qho1d, QHO1D_dense_1_17_1%Eq_pos, QHO1D_non_alloc%Eq_pos)
  CALL Logical_Test(test_qho1d, test1=error_qho1d,                      test2=.FALSE.,info="QHO1D%Eq_pos  deallocated ?")
  CALL Equal_tensor(error_qho1d, QHO1D_dense_1_17_1%Scale_q, QHO1D_non_alloc%Scale_q)
  CALL Logical_Test(test_qho1d, test1=error_qho1d,                      test2=.FALSE.,info="QHO1D%Scale_q deallocated ?")


  !-----------------------------------sum up-----------------------------------
  CALL Finalize_Test(test_qho1d)
  

END PROGRAM

