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
PROGRAM test_construct_op_1D
  !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
  USE QDUtil_m
  USE QDUtil_Test_m
  USE Cavity_mode_m
  USE Operator_1D_m
  IMPLICIT NONE


  logical                       :: Debug = .FALSE.

  TYPE(Cavity_mode_t)           :: Mode 

  TYPE(Operator_1D_t)           :: H_ho_1D_diag_1_6                            ! matrix of the one-dimensional harmonic Hamiltonian associated with HO D. The first number = eigenpulsation of the HO, the second one = the number of basis functions to expand the matrix in.
  TYPE(Operator_1D_t)           :: H_ho_1D_diag_14_6                           ! 14 = w (eigenpulsation); 6 = Nb (number of Ho basis vectors)
  TYPE(Operator_1D_t)           :: H_ho_1D_diag_1_17
  TYPE(Operator_1D_t)           :: H_ho_1D_diag_14_17
  TYPE(Operator_1D_t)           :: H_ho_1D_dense_1_6
  TYPE(Operator_1D_t)           :: H_ho_1D_dense_14_17
  real(kind=Rkind), allocatable :: H_ho_1D_diag_theo_14_17(:)
  real(kind=Rkind), allocatable :: H_ho_1D_dense_theo_1_6(:,:)
  real(kind=Rkind), allocatable :: H_ho_1D_dense_theo_14_17(:,:)

  TYPE(Operator_1D_t)           :: x_ho_1D_band_1_6_1
  TYPE(Operator_1D_t)           :: x_ho_1D_band_14_6_1                         ! 14 = w (eigenpulsation); 6 = Nb (number of Ho basis vectors); 1 = m (mass associated to the HO)
  TYPE(Operator_1D_t)           :: x_ho_1D_band_1_17_1
  TYPE(Operator_1D_t)           :: x_ho_1D_band_1_6_7
  TYPE(Operator_1D_t)           :: x_ho_1D_band_14_17_7
  TYPE(Operator_1D_t)           :: x_ho_1D_dense_1_6_1
  TYPE(Operator_1D_t)           :: x_ho_1D_dense_14_17_7
  real(kind=Rkind), allocatable :: x_ho_1D_band_theo_1_17_1(:,:)
  real(kind=Rkind), allocatable :: x_ho_1D_dense_theo_1_6_1(:,:)
  real(kind=Rkind), allocatable :: x_ho_1D_dense_theo_14_17_7(:,:)

  TYPE(Operator_1D_t)           :: N_ho_1D_diag_6
  TYPE(Operator_1D_t)           :: N_ho_1D_diag_17
  TYPE(Operator_1D_t)           :: N_ho_1D_dense_6
  TYPE(Operator_1D_t)           :: N_ho_1D_dense_17
  real(kind=Rkind), allocatable :: N_ho_1D_dense_theo_17(:,:)

  TYPE(test_t)                  :: test_construct
  logical                       :: error_construct = .FALSE.

  integer                       :: i


  !-----------------------------Test initialization----------------------------
  CALL Initialize_Test(test_construct, test_name="OUT/test_file_cnstrct_op_1D")
  IF (Debug) THEN
    WRITE(out_unit,*)
    WRITE(out_unit,*) "--------------Cavity mode constructed by MolecCav_Read_cavity_mode--------------"
    CALL Write_cavity_mode(Mode)
    WRITE(out_unit,*) "------------End Cavity mode constructed by MolecCav_Read_cavity_mode------------"
  END IF

  
  !----------------------------Construct Cavity mode---------------------------
  CALL Read_cavity_mode(Mode, nio=in_unit)


  !------------------------Construct H matricies to test-----------------------
  CALL Construct_Operator_1D(H_ho_1D_diag_1_6,    "Hamiltonian", Mode=Mode, Debug=Debug)

  Mode%w  = 14*ONE
  CALL Construct_Operator_1D(H_ho_1D_diag_14_6,   "Hamiltonian", Mode=Mode, Debug=Debug)

  Mode%w  = ONE
  Mode%Nb = 17
  CALL Construct_Operator_1D(H_ho_1D_diag_1_17,   "Hamiltonian", Mode=Mode, Debug=Debug)

  Mode%w  = 14*ONE
  CALL Construct_Operator_1D(H_ho_1D_diag_14_17,  "Hamiltonian", Mode=Mode, Debug=Debug)

  Mode%w  = ONE
  Mode%Nb = 6
  CALL Construct_Operator_1D(H_ho_1D_dense_1_6,   "Hamiltonian", Dense=.TRUE., Mode=Mode, Debug=Debug)
                                 
  Mode%w  = 14*ONE
  Mode%Nb = 17
  CALL Construct_Operator_1D(H_ho_1D_dense_14_17, "Hamiltonian", Dense=.TRUE., Mode=Mode, Debug=Debug)


  !------------------------Construct x matricies to test-----------------------
  Mode%w  = ONE
  Mode%Nb = 6
  Mode%m  = ONE
  CALL Construct_Operator_1D(x_ho_1D_band_1_6_1,    "Position", Mode=Mode, Debug=Debug)

  Mode%w  = 14*ONE
  CALL Construct_Operator_1D(x_ho_1D_band_14_6_1,   "Position", Mode=Mode, Debug=Debug)

  Mode%w  = ONE
  Mode%Nb = 17
  CALL Construct_Operator_1D(x_ho_1D_band_1_17_1,   "Position", Mode=Mode, Debug=Debug)

  Mode%Nb = 6
  Mode%m  = SEVEN
  CALL Construct_Operator_1D(x_ho_1D_band_1_6_7,    "Position", Mode=Mode, Debug=Debug)

  Mode%w  = 14*ONE
  Mode%Nb = 17
  Mode%m  = SEVEN
  CALL Construct_Operator_1D(x_ho_1D_band_14_17_7,  "Position", Mode=Mode, Debug=Debug)

  Mode%w  = ONE
  Mode%Nb = 6
  Mode%m  = ONE
  CALL Construct_Operator_1D(x_ho_1D_dense_1_6_1,   "Position", Dense=.TRUE., Mode=Mode, Debug=Debug)

  Mode%w  = 14*ONE
  Mode%Nb = 17
  Mode%m  = SEVEN
  CALL Construct_Operator_1D(x_ho_1D_dense_14_17_7, "Position", Dense=.TRUE., Mode=Mode, Debug=Debug)


  !------------------------Construct N matricies to test-----------------------
  Mode%Nb = 6
  CALL Construct_Operator_1D(N_ho_1D_diag_6,   "Nb_photons", Mode=Mode, Debug=Debug)

  Mode%Nb = 17
  CALL Construct_Operator_1D(N_ho_1D_diag_17,  "Nb_photons", Mode=Mode, Debug=Debug)

  Mode%Nb = 6
  CALL Construct_Operator_1D(N_ho_1D_dense_6,  "Nb_photons", Dense=.TRUE., Mode=Mode, Debug=Debug)

  Mode%Nb = 17
  CALL Construct_Operator_1D(N_ho_1D_dense_17, "Nb_photons", Dense=.TRUE., Mode=Mode, Debug=Debug)


  !-----------------------Construct reference H matricies----------------------
  ALLOCATE(H_ho_1D_diag_theo_14_17(17))
  ALLOCATE(H_ho_1D_dense_theo_1_6(6,6))
  ALLOCATE(H_ho_1D_dense_theo_14_17(17,17))

  H_ho_1D_diag_theo_14_17  = ZERO
  H_ho_1D_dense_theo_1_6   = ZERO
  H_ho_1D_dense_theo_14_17 = ZERO

  DO i = 1, 17
    H_ho_1D_diag_theo_14_17(i)    = 14*(i - ONE + HALF)
    H_ho_1D_dense_theo_14_17(i,i) = 14*(i - ONE + HALF)
    IF (i <= 6) THEN
      H_ho_1D_dense_theo_1_6(i,i) = (i - ONE + HALF)
    END IF
  END DO

  IF (Debug) THEN
    WRITE(out_unit,*) 
    WRITE(out_unit,*) "H_{Reference 1 (diagonal, w = 14, Nb = 17)}"
    CALL Write_Vec(H_ho_1D_diag_theo_14_17,  out_unit, 1, info="H_{Reference 1}")

    WRITE(out_unit,*) 
    WRITE(out_unit,*) "H_{Reference 2 (dense, w = 1, Nb = 6)}"
    CALL Write_Mat(H_ho_1D_dense_theo_1_6,   out_unit, Size(H_ho_1D_dense_theo_1_6,   dim=2), info="H_{Reference 2}")
  
    WRITE(out_unit,*) 
    WRITE(out_unit,*) "H_{Reference 3 (dense, w = 14, Nb = 17)}"
    CALL Write_Mat(H_ho_1D_dense_theo_14_17, out_unit, Size(H_ho_1D_dense_theo_14_17, dim=2), info="H_{Reference 3}")
  
  END IF


  !-----------------------Construct reference x matricies----------------------
  ALLOCATE(x_ho_1D_band_theo_1_17_1(17,3))
  ALLOCATE(x_ho_1D_dense_theo_1_6_1(6,6))
  ALLOCATE(x_ho_1D_dense_theo_14_17_7(17,17))

  x_ho_1D_band_theo_1_17_1   = ZERO
  x_ho_1D_dense_theo_1_6_1   = ZERO
  x_ho_1D_dense_theo_14_17_7 = ZERO

  DO i = 1, 16                                                                 ! /!\ Fortran counts from 1 to Nb !!! /!\ Nb-1 not to have Band_val_R(i+1) out of range
    x_ho_1D_band_theo_1_17_1(i,1)     = SQRT(REAL(i,kind=Rkind))
    x_ho_1D_band_theo_1_17_1(i+1,3)   = SQRT(REAL(i,kind=Rkind))
    x_ho_1D_dense_theo_14_17_7(i,i+1) = SQRT(REAL(i,kind=Rkind))
    x_ho_1D_dense_theo_14_17_7(i+1,i) = SQRT(REAL(i,kind=Rkind))
    IF (i <= 5) THEN
      x_ho_1D_dense_theo_1_6_1(i,i+1) = SQRT(REAL(i,kind=Rkind))
      x_ho_1D_dense_theo_1_6_1(i+1,i) = SQRT(REAL(i,kind=Rkind))
    END IF
  END DO

  x_ho_1D_band_theo_1_17_1   = x_ho_1D_band_theo_1_17_1   / SQRT(TWO * 1.0_Rkind  * 1.0_Rkind)
  x_ho_1D_dense_theo_1_6_1   = x_ho_1D_dense_theo_1_6_1   / SQRT(TWO * 1.0_Rkind  * 1.0_Rkind)
  x_ho_1D_dense_theo_14_17_7 = x_ho_1D_dense_theo_14_17_7 / SQRT(TWO * 14.0_Rkind * 7.0_Rkind)

  IF (Debug) THEN
    WRITE(out_unit,*) 
    WRITE(out_unit,*) "x_{Reference 1 (band, w = 1.0, Nb = 17, m = 1.0)}"
    CALL Write_Mat(x_ho_1D_band_theo_1_17_1,   out_unit, 3, info="x_{Reference 1}")

    WRITE(out_unit,*) 
    WRITE(out_unit,*) "x_{Reference 2 (dense, w = 1.0, Nb = 6, m = 1.0)}"
    CALL Write_Mat(x_ho_1D_dense_theo_1_6_1,   out_unit, Size(H_ho_1D_dense_theo_1_6,   dim=2), info="x_{Reference 2}")
  
    WRITE(out_unit,*) 
    WRITE(out_unit,*) "x_{Reference 3 (dense, w = 14.0, Nb = 17, m = 7.0)}"
    CALL Write_Mat(x_ho_1D_dense_theo_14_17_7, out_unit, Size(H_ho_1D_dense_theo_14_17, dim=2), info="x_{Reference 3}")
  
  END IF


  !-----------------------Construct reference N matricies----------------------
  ALLOCATE(N_ho_1D_dense_theo_17(17,17))
  N_ho_1D_dense_theo_17 = ZERO

  DO i = 1, 17                                                             ! /!\ Fortran counts from 1 to Nb !!! /!\
    N_ho_1D_dense_theo_17(i,i) = i - 1
  END DO

  IF (Debug) THEN
    WRITE(out_unit,*) 
    WRITE(out_unit,*) "N_{Reference (diagonal, Nb = 17)}"
    CALL Write_Mat(N_ho_1D_dense_theo_17, out_unit, Size(N_ho_1D_dense_theo_17, dim=2), info="N_{Reference}")
  END IF
  FLUSH(out_unit)

  !--------------------------------Comparisons H-------------------------------
  !H_ho_1D_diag_14_17%Diag_val_R(1) = 0
  !WRITE(out_unit,*) 
  !WRITE(out_unit,*) "H_{Spurious 1 (diagonal, w = 14, Nb = 17)}"
  !DO i = 1, 17
  !  WRITE(out_unit,*) H_ho_1D_diag_theo_14_17(i)
  !END DO

  CALL Equal_R_R_vector(error_construct, H_ho_1D_diag_theo_14_17, H_ho_1D_diag_14_17%Diag_val_R)
  CALL Logical_Test(test_construct, test1=error_construct, test2=.FALSE., info="H_ho_1D_diag_14_17%Diag_val_R well initialized ?")
  IF (error_construct .AND. Debug) THEN
    CALL Write_Vec(H_ho_1D_diag_theo_14_17, out_unit, Size(H_ho_1D_diag_theo_14_17), info="H_ho_1D_diag_theo_14_17")
    CALL Write_Vec(H_ho_1D_diag_14_17%Diag_val_R, out_unit, Size(H_ho_1D_diag_14_17%Diag_val_R), in&
                  &fo="H_ho_1D_diag_14_17%Diag_val_R")
  END IF
  
  !H_ho_1D_dense_1_6%Dense_val_R(1,2) = 1
  !WRITE(out_unit,*) 
  !WRITE(out_unit,*) "H_{Spurious 2 (diagonal, w = 1, Nb = 6)}"
  !DO i = 1, 17
  !  WRITE(out_unit,*) H_ho_1D_dense_1_6%Dense_val_R(i,:)
  !END DO

  CALL Equal_R_R_matrix(error_construct, H_ho_1D_dense_theo_1_6, H_ho_1D_dense_1_6%Dense_val_R)
  CALL Logical_Test(test_construct, test1=error_construct, test2=.FALSE., info="H_ho_1D_dense_1_6%Dense_val_R well initialized ?")
  IF (error_construct .AND. Debug) THEN
    CALL Write_Mat(H_ho_1D_dense_theo_1_6, out_unit, Size(H_ho_1D_dense_theo_1_6, dim=2), info="H_ho_1D_dense_theo_1_6")
    CALL Write_Mat(H_ho_1D_dense_1_6%Dense_val_R, out_unit, Size(H_ho_1D_dense_1_6%Dense_val_R, dim=2), in&
                  &fo="H_ho_1D_dense_1_6%Dense_val_R")
  END IF

  !H_ho_1D_dense_14_17%Dense_val_R(1,2) = 1
  !WRITE(out_unit,*) 
  !WRITE(out_unit,*) "H_{Spurious 3 (diagonal, w = 14, Nb = 17)}"
  !DO i = 1, 17
  !  WRITE(out_unit,*) H_ho_1D_dense_14_17%Dense_val_R(i,:)
  !END DO

  CALL Equal_R_R_matrix(error_construct, H_ho_1D_dense_theo_14_17, H_ho_1D_dense_14_17%Dense_val_R)
  CALL Logical_Test(test_construct, test1=error_construct, test2=.FALSE., info="H_ho_1D_dense_14_17&
                   &%Dense_val_R well initialized ?")
  IF (error_construct .AND. Debug) THEN
    CALL Write_Mat(H_ho_1D_dense_theo_14_17, out_unit, Size(H_ho_1D_dense_theo_14_17, dim=2), info="H_ho_1D_dense_theo_14_17")
    CALL Write_Mat(H_ho_1D_dense_14_17%Dense_val_R, out_unit, Size(H_ho_1D_dense_14_17%Dense_val_R, dim=2), in&
                  &fo="H_ho_1D_dense_14_17%Dense_val_R")
  END IF

  !H_ho_1D_diag_1_6%Diag_val_R(1) = 1
  !WRITE(out_unit,*) 
  !WRITE(out_unit,*) "H_{Spurious 4 (diagonal, w = 1, Nb = 6)}"
  !DO i = 1, 17
  !  WRITE(out_unit,*) H_ho_1D_diag_1_6%Dense_val_R(i,:)
  !END DO

  DO i = 1, 6
    CALL Equal_R_R_scalar(error_construct, H_ho_1D_dense_1_6%Dense_val_R(i,i), H_ho_1D_diag_1_6%Diag_val_R(i))
    CALL Logical_Test(test_construct, test1=error_construct, test2=.FALSE., info="H_ho_1D_diag_1_6%&
                     &Diag_val_R("//TO_string(i)//") well initialized ?")
    IF (error_construct .AND. Debug) WRITE(out_unit,*) "H_ho_1D_dense_1_6%Dense_val_R("//TO_string(&
                                                       &i)//","//TO_string(i)//"), H_ho_1D_diag_1_6&
                                                       &%Diag_val_R("//TO_string(i)//")"
  END DO

  CALL Equal_R_R_vector(error_construct, H_ho_1D_diag_1_17%Diag_val_R(1:6), H_ho_1D_diag_1_6%Diag_val_R)
  CALL Logical_Test(test_construct, test1=error_construct, test2=.FALSE., info="H_ho_1D_diag_1_17%D&
                     &iag_val_R well initialized ?")
  IF (error_construct .AND. Debug) THEN
    CALL Write_Vec(H_ho_1D_diag_1_17%Diag_val_R, out_unit, Size(H_ho_1D_diag_1_17%Diag_val_R), info="H_ho_1D_diag_1_17%Diag_val_R")
    CALL Write_Vec(H_ho_1D_diag_1_6%Diag_val_R, out_unit, Size(H_ho_1D_diag_1_6%Diag_val_R), in&
                  &fo="H_ho_1D_diag_1_6%Diag_val_R")
  END IF

  !H_ho_1D_diag_14_6%Diag_val_R(1) = 1
  !WRITE(out_unit,*) 
  !WRITE(out_unit,*) "H_{Spurious 5 (diagonal, w = 14, Nb = 6)}"

  CALL Equal_R_R_vector(error_construct, H_ho_1D_diag_14_6%Diag_val_R, H_ho_1D_diag_14_17%Diag_val_R(1:6))
  CALL Logical_Test(test_construct, test1=error_construct, test2=.FALSE., info="H_ho_1D_diag_14_6%D&
                   &iag_val_R well initialized ?")
  IF (error_construct .AND. Debug) THEN
    CALL Write_Vec(H_ho_1D_diag_14_6%Diag_val_R, out_unit, Size(H_ho_1D_diag_14_6%Diag_val_R), info="H_ho_1D_diag_14_6%Diag_val_R")
    CALL Write_Vec(H_ho_1D_diag_14_17%Diag_val_R, out_unit, Size(H_ho_1D_diag_14_17%Diag_val_R), in&
                  &fo="H_ho_1D_diag_14_17%Diag_val_R")
  END IF


  !-------------------------------Comparisons x-------------------------------
  CALL Equal_R_R_matrix(error_construct, x_ho_1D_band_theo_1_17_1, x_ho_1D_band_1_17_1%Band_val_R)
  CALL Logical_Test(test_construct, test1=error_construct, test2=.FALSE., info="x_ho_1D_band_1_17_1&
                   &%Band_val_R well initialized ?")
  IF (error_construct .AND. Debug) THEN
    CALL Write_Mat(x_ho_1D_band_theo_1_17_1, out_unit, Size(x_ho_1D_band_theo_1_17_1, dim=2), info="x_ho_1D_band_theo_1_17_1")
    CALL Write_Mat(x_ho_1D_band_theo_1_17_1, out_unit, Size(x_ho_1D_band_theo_1_17_1, dim=2), in&
                  &fo="x_ho_1D_band_theo_1_17_1")
  END IF

  !x_ho_1D_dense_1_6_1%Dense_val_R(1,1) = 5
  !WRITE(out_unit,*)
  !WRITE(out_unit,*) "x_{Spurious 2 (dense, w = 1, Nb = 6, m = 1)}"
  !DO i = 1, 6
  !  WRITE(out_unit,*) x_ho_1D_dense_1_6_1%Dense_val_R(i,:)
  !END DO

  CALL Equal_R_R_matrix(error_construct, x_ho_1D_dense_theo_1_6_1, x_ho_1D_dense_1_6_1%Dense_val_R)
  CALL Logical_Test(test_construct, test1=error_construct, test2=.FALSE., info="x_ho_1D_dense_1_6_1&
                   &%Dense_val_R well initialized ?")
  IF (error_construct .AND. Debug) THEN
    CALL Write_Mat(x_ho_1D_dense_theo_1_6_1, out_unit, Size(x_ho_1D_dense_theo_1_6_1, dim=2), info="x_ho_1D_dense_theo_1_6_1")
    CALL Write_Mat(x_ho_1D_dense_1_6_1%Dense_val_R, out_unit, Size(x_ho_1D_dense_1_6_1%Dense_val_R, dim=2), in&
                  &fo="x_ho_1D_dense_1_6_1%Dense_val_R")
  END IF

  !x_ho_1D_dense_14_17_7%Dense_val_R(1,1) = 5
  !WRITE(out_unit,*)
  !WRITE(out_unit,*) "x_{Spurious 3 (dense, w = 14, Nb = 17, m = 7)}"
  !DO i = 1, 17
  !  WRITE(out_unit,*) x_ho_1D_dense_14_17_7%Dense_val_R(i,:)
  !END DO

  CALL Equal_R_R_matrix(error_construct, x_ho_1D_dense_theo_14_17_7, x_ho_1D_dense_14_17_7%Dense_val_R)
  CALL Logical_Test(test_construct, test1=error_construct, test2=.FALSE., info="x_ho_1D_dense_14_17&
                  &_7%Dense_val_R well initialized ?")
  IF (error_construct .AND. Debug) THEN
    CALL Write_Mat(x_ho_1D_dense_theo_14_17_7, out_unit, Size(x_ho_1D_dense_theo_14_17_7, dim=2), info="x_ho_1D_dense_theo_14_17_7")
    CALL Write_Mat(x_ho_1D_dense_14_17_7%Dense_val_R, out_unit, Size(x_ho_1D_dense_14_17_7%Dense_val_R, dim=2), in&
                  &fo="x_ho_1D_dense_14_17_7%Dense_val_R")
  END IF

  !x_ho_1D_band_1_6_1%Band_val_R(6,1) = 15
  !WRITE(out_unit,*) "x_{Spurious 4 (Band, w = 1, Nb = 6, m = 1)}"
  !DO i = 1, 6
  !  WRITE(out_unit,*) x_ho_1D_band_1_6_1%Band_val_R(i,:)
  !END DO

  CALL Equal_R_R_vector(error_construct, x_ho_1D_band_1_17_1%Band_val_R(1:5,1), x_ho_1D_band_1_6_1%Band_val_R(1:5,1))
  CALL Logical_Test(test_construct, test1=error_construct, test2=.FALSE., info="x_ho_1D_band_1_6_1%&
                   &Band_val_R(:,1) well initialized ?")
  IF (error_construct .AND. Debug) THEN
    CALL Write_Vec(x_ho_1D_band_1_17_1%Band_val_R(1:5,1), out_unit, Size(x_ho_1D_band_1_17_1%Band_v&
                  &al_R(1:5,1)), info="x_ho_1D_band_1_17_1%Band_val_R(1:5,1)")
    CALL Write_Vec(x_ho_1D_band_1_6_1%Band_val_R(1:5,1), out_unit, Size(x_ho_1D_band_1_6_1%Band_val&
                  &_R(1:5,1)), info="x_ho_1D_band_1_6_1%Band_val_R(1:5,1)")
  END IF
  CALL Equal_R_R_vector(error_construct, x_ho_1D_band_1_17_1%Band_val_R(1:6,2), x_ho_1D_band_1_6_1%Band_val_R(:,2))
  CALL Logical_Test(test_construct, test1=error_construct, test2=.FALSE., info="x_ho_1D_band_1_6_1%&
                   &Band_val_R(:,2) well initialized ?")
  IF (error_construct .AND. Debug) THEN
    CALL Write_Vec(x_ho_1D_band_1_17_1%Band_val_R(1:6,2), out_unit, Size(x_ho_1D_band_1_17_1%Band_v&
                  &al_R(1:6,2)), info="x_ho_1D_band_1_17_1%Band_val_R(1:6,2)")
    CALL Write_Vec(x_ho_1D_band_1_6_1%Band_val_R(1:6,2), out_unit, Size(x_ho_1D_band_1_6_1%Band_val&
                  &_R(1:6,2)), info="x_ho_1D_band_1_6_1%Band_val_R(1:6,2)")
  END IF
  CALL Equal_R_R_vector(error_construct, x_ho_1D_band_1_17_1%Band_val_R(1:6,3), x_ho_1D_band_1_6_1%Band_val_R(:,3))
  CALL Logical_Test(test_construct, test1=error_construct, test2=.FALSE., info="x_ho_1D_band_1_6_1%&
                  &Band_val_R(:,3) well initialized ?")
  IF (error_construct .AND. Debug) THEN
    CALL Write_Vec(x_ho_1D_band_1_17_1%Band_val_R(1:6,3), out_unit, Size(x_ho_1D_band_1_17_1%Band_v&
                  &al_R(1:6,3)), info="x_ho_1D_band_1_17_1%Band_val_R(1:6,3)")
    CALL Write_Vec(x_ho_1D_band_1_6_1%Band_val_R(1:6,3), out_unit, Size(x_ho_1D_band_1_6_1%Band_val&
                  &_R(1:6,3)), info="x_ho_1D_band_1_6_1%Band_val_R(1:6,3)")
  END IF
  CALL Equal_R_R_scalar(error_construct, ZERO, x_ho_1D_band_1_6_1%Band_val_R(6,1))
  CALL Logical_Test(test_construct, test1=error_construct, test2=.FALSE., info="x_ho_1D_band_1_6_1%&
                  &Band_val_R(6,1) well initialized ?")
  IF (error_construct .AND. Debug) WRITE(out_unit,*) "x_ho_1D_band_1_6_1%Band_val_R(6,1)", x_ho_1D_&
                  &band_1_6_1%Band_val_R(6,1), "/= 0"

  !x_ho_1D_band_14_17_7%Band_val_R(1,3) = 5
  DO i = 1, 16
    CALL Equal_R_R_scalar(error_construct, x_ho_1D_dense_theo_14_17_7(i+1,i), x_ho_1D_band_14_17_7%Band_val_R(i,1))
    CALL Logical_Test(test_construct, test1=error_construct, test2=.FALSE., info="x_ho_1D_band_14_1&
                     &7_7%Band_val_R("//TO_string(i)//",1) well initialized ?")
    IF (error_construct .AND. Debug) WRITE(out_unit,*) x_ho_1D_band_14_17_7%Band_val_R(i,1), x_ho_1D_dense_theo_14_17_7(i+1,i)

    CALL Equal_R_R_scalar(error_construct, x_ho_1D_dense_theo_14_17_7(i,i+1), x_ho_1D_band_14_17_7%Band_val_R(i+1,3))
    CALL Logical_Test(test_construct, test1=error_construct, test2=.FALSE., info="x_ho_1D_band_14_1&
                     &7_7%Band_val_R("//TO_string(i+1)//",3) well initialized ?")
    IF (error_construct .AND. Debug) WRITE(out_unit,*) x_ho_1D_band_14_17_7%Band_val_R(i+1,3), x_ho_1D_dense_theo_14_17_7(i,i+1)

  END DO

  CALL Equal_R_R_scalar(error_construct, ZERO, x_ho_1D_band_14_17_7%Band_val_R(17,1))
  CALL Logical_Test(test_construct, test1=error_construct, test2=.FALSE., info="x_ho_1D_band_14_17_&
                   &7%Band_val_R(17,1) well initialized ?")
  IF (error_construct .AND. Debug) WRITE(out_unit,*) x_ho_1D_band_14_17_7%Band_val_R(17,1), "/= 0"
  CALL Equal_R_R_scalar(error_construct, ZERO, x_ho_1D_band_14_17_7%Band_val_R(1,3))
  CALL Logical_Test(test_construct, test1=error_construct, test2=.FALSE., info="x_ho_1D_band_14_17_&
                   &7%Band_val_R(1,3) well initialized ?")
  IF (error_construct .AND. Debug) WRITE(out_unit,*) x_ho_1D_band_14_17_7%Band_val_R(1,3), "/= 0"

  !x_ho_1D_band_14_6_1%Band_val_R(6,1) = 5
  !WRITE(out_unit,*) "x_{Spurious 4 (Band, w = 1, Nb = 6, m = 1)}"
  !DO i = 1, 6
  !  WRITE(out_unit,*) x_ho_1D_band_1_6_1%Band_val_R(i,:)
  !END DO

  CALL Equal_R_R_matrix(error_construct, x_ho_1D_band_1_6_1%Band_val_R/(SQRT(14.0_Rkind)), x_ho_1D_band_14_6_1%Band_val_R)
  CALL Logical_Test(test_construct, test1=error_construct, test2=.FALSE., info="x_ho_1D_band_14_6_1&
                   &%Band_val_R well initialized ?")
  IF (error_construct .AND. Debug) THEN
    CALL Write_Mat(x_ho_1D_band_1_6_1%Band_val_R/(SQRT(14.0_Rkind)), out_unit, Size(x_ho_1D_band_1_&
                  &6_1%Band_val_R/(SQRT(14.0_Rkind)), dim=2), info="x_ho_1D_band_1_6_1%Band_val_R/(SQRT(14.0_Rkind))")
    CALL Write_Mat(x_ho_1D_band_14_6_1%Band_val_R, out_unit, Size(x_ho_1D_band_14_6_1%Band_val_R, d&
                  &im=2), info="x_ho_1D_band_14_6_1%Band_val_R")
  END IF

  !x_ho_1D_band_1_6_7%Band_val_R(6,1) = 5
  !WRITE(out_unit,*) "x_{Spurious 4 (Band, w = 1, Nb = 6, m = 1)}"
  !DO i = 1, 6
  !  WRITE(out_unit,*) x_ho_1D_band_1_6_1%Band_val_R(i,:)
  !END DO

  CALL Equal_R_R_matrix(error_construct, x_ho_1D_band_1_6_1%Band_val_R/(SQRT(SEVEN)), x_ho_1D_band_1_6_7%Band_val_R)
  CALL Logical_Test(test_construct, test1=error_construct, test2=.FALSE., info="x_ho_1D_band_1_6_7%&
                   &Band_val_R well initialized ?")
  IF (error_construct .AND. Debug) THEN
    CALL Write_Mat(x_ho_1D_band_1_6_1%Band_val_R/(SQRT(SEVEN)), out_unit, Size(x_ho_1D_band_1_6_1%B&
                  &and_val_R/(SQRT(SEVEN)), dim=2), info="x_ho_1D_band_1_6_1%Band_val_R/(SQRT(SEVEN))")
    CALL Write_Mat(x_ho_1D_band_1_6_7%Band_val_R, out_unit, Size(x_ho_1D_band_1_6_7%Band_val_R, dim&
                  &=2), info="x_ho_1D_band_1_6_7%Band_val_R")
  END IF


  !-------------------------------Comparisons N-------------------------------
  !N_ho_1D_dense_17%Dense_val_R(1,1) = 5

  CALL Equal_R_R_matrix(error_construct, N_ho_1D_dense_theo_17, N_ho_1D_dense_17%Dense_val_R)
  CALL Logical_Test(test_construct, test1=error_construct, test2=.FALSE., info="N_ho_1D_dense_17%De&
                   &nse_val_R well initialized ?")
  IF (error_construct .AND. Debug) THEN
    CALL Write_Mat(N_ho_1D_dense_theo_17, out_unit, Size(N_ho_1D_dense_theo_17, dim=2), info="N_ho_1D_dense_theo_17")
    CALL Write_Mat(N_ho_1D_dense_17%Dense_val_R, out_unit, Size(N_ho_1D_dense_17%Dense_val_R, dim=&
                  &2), info="N_ho_1D_dense_17%Dense_val_R")
  END IF

  !N_ho_1D_diag_17%Diag_val_R(1) = 5
  !N_ho_1D_dense_6%Dense_val_R(1,1) = 5
  !N_ho_1D_diag_6%Diag_val_R(1) = 5
  
  DO i = 1, 6
    CALL Equal_R_R_scalar(error_construct, N_ho_1D_dense_theo_17(i,i), N_ho_1D_diag_17%Diag_val_R(i))
    CALL Logical_Test(test_construct, test1=error_construct, test2=.FALSE., info="N_ho_1D_diag_17%D&
                     &iag_val_R("//TO_string(i)//") well initialized ?")
    IF (error_construct .AND. Debug) WRITE(out_unit,*) N_ho_1D_dense_theo_17(i,i), N_ho_1D_diag_17%Diag_val_R(i)
  END DO

  CALL Equal_R_R_matrix(error_construct, N_ho_1D_dense_theo_17(1:6,1:6), N_ho_1D_dense_6%Dense_val_R)
  CALL Logical_Test(test_construct, test1=error_construct, test2=.FALSE., info="N_ho_1D_dense_6%Dense_val_R well initialized ?")
  IF (error_construct .AND. Debug) THEN
    CALL Write_Mat(N_ho_1D_dense_theo_17(1:6,1:6), out_unit, Size(N_ho_1D_dense_theo_17(1:6,1:6), d&
                  &im=2), info="N_ho_1D_dense_theo_17(1:6,1:6)")
    CALL Write_Mat(N_ho_1D_dense_6%Dense_val_R, out_unit, Size(N_ho_1D_dense_6%Dense_val_R, dim=2),&
                 & info="N_ho_1D_dense_6%Dense_val_R")
  END IF

  DO i = 1, 6
    CALL Equal_R_R_scalar(error_construct, N_ho_1D_dense_theo_17(i,i), N_ho_1D_diag_6%Diag_val_R(i))
    CALL Logical_Test(test_construct, test1=error_construct, test2=.FALSE., info="N_ho_1D_diag_6%Di&
                     &ag_val_R("//TO_string(i)//") well initialized ?")
    IF (error_construct .AND. Debug) WRITE(out_unit,*) N_ho_1D_dense_theo_17(i,i), N_ho_1D_diag_6%Diag_val_R(i)
  END DO


  !-----------------------------------sum up-----------------------------------
  CALL Finalize_Test(test_construct)
  

  CONTAINS


  SUBROUTINE Equal_R_R_scalar(error, Rl_1, Rl_2)
    USE QDUtil_m
    IMPLICIT NONE 

    logical,          intent(inout) :: error
    real(kind=Rkind), intent(in)    :: Rl_1
    real(kind=Rkind), intent(in)    :: Rl_2
    
    real(kind=Rkind), parameter     :: Threshold   = 1E-10_Rkind
    logical, parameter              :: Debug_local = .FALSE.

    IF (ABS(Rl_1 - Rl_2) > Threshold) THEN
      error = .TRUE.
      IF (Debug_local) WRITE(out_unit,*) "The two numbers are not close enough to be considered equ&
                                         &al : R_1 =", Rl_1, "R_2 =", Rl_2, "|R_1-R_2| = ", ABS(Rl_1 - Rl_2)
    ELSE 
      error = .FALSE.
      IF (Debug_local) WRITE(out_unit,*) "The two numbers are close enough to be considered equal :&
                                         & R_1 =", Rl_1, "R_2 =", Rl_2, "|R_1-R_2| = ", ABS(Rl_1 - Rl_2)
    END IF

  END SUBROUTINE Equal_R_R_scalar
  

  SUBROUTINE Equal_R_R_vector(error, Rl_1, Rl_2)
    USE QDUtil_m
    IMPLICIT NONE 

    logical,          intent(inout) :: error
    real(kind=Rkind), intent(in)    :: Rl_1(:)
    real(kind=Rkind), intent(in)    :: Rl_2(:)
    
    real(kind=Rkind), parameter     :: Threshold   = 1E-10_Rkind
    logical, parameter              :: Debug_local = .FALSE.
    integer                         :: Nb_1

    Nb_1 = Size(Rl_1)
    IF (Nb_1 /= Size(Rl_2)) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "The two vectors must have same dimensions to compare them. Please, check initialization."
      WRITE(out_unit,*) "Size(Rl_1) = "//TO_string(Size(Rl_1))//"; Size(Rl_2) = "//TO_string(Size(Rl_2))
      STOP "### The two vectors must have same dimensions to compare them. Please, check initialization."
    END IF 

    IF (ANY(ABS(Rl_1 - Rl_2) > Threshold)) THEN
      error = .TRUE.
      IF (Debug_local) THEN
        WRITE(out_unit,*) "The two vectors are not close enough to be considered equal :"
        CALL Write_Vec(Rl_1, out_unit, Nb_1, info="R_1(:)")
        CALL Write_Vec(Rl_2, out_unit, Nb_1, info="R_2(:)")
        CALL Write_Vec(ABS(Rl_1 - Rl_2), out_unit, Nb_1, info="|R_1-R_2| = ")
      END IF 

    ELSE 
      error = .FALSE.
      IF (Debug_local) THEN
        WRITE(out_unit,*) "The two vectors are close enough to be considered equal :"
        CALL Write_Vec(Rl_1, out_unit, Nb_1, info="R_1(:)")
        CALL Write_Vec(Rl_2, out_unit, Nb_1, info="R_2(:)")
        CALL Write_Vec(ABS(Rl_1 - Rl_2), out_unit, Nb_1, info="|R_1-R_2| = ")
      END IF
    END IF 

  END SUBROUTINE Equal_R_R_vector
  

  SUBROUTINE Equal_R_R_matrix(error, Rl_1, Rl_2)
    USE QDUtil_m
    IMPLICIT NONE 

    logical,          intent(inout) :: error
    real(kind=Rkind), intent(in)    :: Rl_1(:,:)
    real(kind=Rkind), intent(in)    :: Rl_2(:,:)
    
    real(kind=Rkind), parameter     :: Threshold   = 1E-10_Rkind
    logical, parameter              :: Debug_local = .FALSE.
    integer                         :: Nb_1, Nb_2

    Nb_1 = Size(Rl_1, dim=1)
    Nb_2 = Size(Rl_2, dim=2)
    IF (Nb_1 /= Size(Rl_2, dim=1) .OR. Nb_2 /= Size(Rl_2, dim=2)) THEN
      WRITE(out_unit,*) "The two matrices must have same dimensions to compare them. Please, check initialization."
      STOP "### The two matrices must have same dimensions to compare them. Please, check initialization."
    END IF 

    IF (ANY(ABS(Rl_1 - Rl_2) > Threshold)) THEN
      error = .TRUE.
      IF (Debug_local) THEN
        WRITE(out_unit,*) "The two matrices are not close enough to be considered equal :"
        CALL Write_Mat(Rl_1, out_unit, Nb_2, info="R_1(:,:)")
        CALL Write_Mat(Rl_2, out_unit, Nb_2, info="R_2(:,:)")
        CALL Write_Mat(ABS(Rl_1 - Rl_2), out_unit, Nb_2, info="|R_1-R_2| = ")
      END IF 

    ELSE 
      error = .FALSE.
      IF (Debug_local) THEN
        WRITE(out_unit,*) "The two matrices are close enough to be considered equal :"
        CALL Write_Mat(Rl_1, out_unit, Nb_2, info="R_1(:,:)")
        CALL Write_Mat(Rl_2, out_unit, Nb_2, info="R_2(:,:)")
        CALL Write_Mat(ABS(Rl_1 - Rl_2), out_unit, Nb_2, info="|R_1-R_2| = ")
      END IF
    END IF 

  END SUBROUTINE Equal_R_R_matrix


END PROGRAM

