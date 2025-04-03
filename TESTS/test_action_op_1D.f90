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
PROGRAM test_action_op_1D
  !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
  USE QDUtil_m
  USE QDUtil_Test_m
  USE Algebra_m
  USE Cavity_mode_m
  USE Elem_op_m
  IMPLICIT NONE


  logical                       :: Debug = .FALSE.

  TYPE(Cavity_mode_t)           :: Mode                                                            ! The well construction of the Operator's matrices is assumed checked by the dedicated test, so hard-coded references matrices are not needed here

  TYPE(Operator_1D_t)           :: H_diag_half                                                     ! Nomenclature : H_<shape>_<eigenpulsation>
  TYPE(Operator_1D_t)           :: H_dense_half
  TYPE(Operator_1D_t)           :: H_diag_3
  TYPE(Operator_1D_t)           :: H_dense_3

  TYPE(Operator_1D_t)           :: x_band_half_1                                                   ! Nomenclature : x_<shape>_<eigenpulsation>_<mass>
  TYPE(Operator_1D_t)           :: x_dense_half_1
  TYPE(Operator_1D_t)           :: x_band_3_1
  TYPE(Operator_1D_t)           :: x_dense_3_1
  TYPE(Operator_1D_t)           :: x_band_3_4
  TYPE(Operator_1D_t)           :: x_dense_3_4

  TYPE(Operator_1D_t)           :: N_diag                                                          ! Nomenclature : N_<shape>
  TYPE(Operator_1D_t)           :: N_dense

  real(kind=Rkind)              :: b_0(3) = [ONE, ZERO, ZERO]                                      ! three vectors of the HO basis set |0>, |1>, |2> 
  real(kind=Rkind)              :: b_1(3) = [ZERO, ONE, ZERO]
  real(kind=Rkind)              :: b_2(3) = [ZERO, ZERO, ONE]
  real(kind=Rkind)              :: Coeff_0 = ONE
  real(kind=Rkind)              :: Coeff_1 = HALF
  real(kind=Rkind)              :: Coeff_2 = THREE
  real(kind=Rkind)              :: Psi_1D(3)                                                       ! an any vector representing the excitation state/wavefunction of the HO = a linear combination of the basis functions /!\ Not normalized yet !
  real(kind=Rkind)              :: Norm                                                            ! SQRT(Coeff_0**2 + Coeff_1**2 + COeff_2**2)
  real(kind=Rkind)              :: Op_psi(3)                                                       ! the resulting vector from the action of a 1D operator upon Psi_1D
  real(kind=Rkind)              :: Op_psi_ana(3)                                                   ! the analytical action = hard-coded reference for comparison

  TYPE(test_t)                  :: test_action
  logical                       :: error_action = .FALSE.

  !IF (Debug) THEN
  !  WRITE(out_unit,*)
  !  WRITE(out_unit,*) "--------------Cavity mode constructed by MolecCav_Read_cavity_mode--------------"
  !  CALL Write_cavity_mode(Mode)
  !  WRITE(out_unit,*) "------------End Cavity mode constructed by MolecCav_Read_cavity_mode------------"
  !END IF


  !-----------------------------Test initialization----------------------------
  CALL Initialize_Test(test_action, test_name="OUT/test_file_actn_op_1D")
  

  !-------------------------Wavefunction initialization------------------------
  Psi_1D(:) = [Coeff_0, Coeff_1, Coeff_2]
  CALL Norm_of(Norm, Psi_1D)
  CALL Normalize(Psi_1D)
  IF (Debug) THEN
    WRITE(out_unit,*)
    WRITE(out_unit,*) "----------------The cavity wavefunction has been initialized as :---------------"
    CALL Write_Vec(Psi_1D, out_unit, Size(Psi_1D), info="Psi_1D")
    WRITE(out_unit,*) "-----------------------------End cavity wavefunction----------------------------"
  END IF


  !-------------------------Cavity mode initialization-------------------------
  CALL Read_cavity_mode(Mode, nio=in_unit)
  Mode%Nb = 3
  IF (Debug) THEN
    WRITE(out_unit,*)
    WRITE(out_unit,*) "--------------Cavity mode constructed by MolecCav_Read_cavity_mode--------------"
    CALL Write_cavity_mode(Mode)
    WRITE(out_unit,*) "------------End Cavity mode constructed by MolecCav_Read_cavity_mode------------"
  END IF

  
  !-------------------------H matricies initialization-------------------------
  Mode%w  = HALF
  CALL Construct_Operator_1D(H_diag_half,  "Hamiltonian", Mode=Mode, Debug=Debug)
  CALL Construct_Operator_1D(H_dense_half, "Hamiltonian", Dense=.TRUE., Mode=Mode, Debug=Debug)

  Mode%w  = THREE
  CALL Construct_Operator_1D(H_diag_3,  "Hamiltonian", Mode=Mode, Debug=Debug)
  CALL Construct_Operator_1D(H_dense_3, "Hamiltonian", Dense=.TRUE., Mode=Mode, Debug=Debug)


  !----------------------------Testing the H actions---------------------------
    !---------------------------------H_{w=1/2}--------------------------------
  CALL Action(Op_psi, H_diag_half, b_0)
  Op_psi_ana(:) = [H_diag_half%Diag_val(1), ZERO, ZERO]
  CALL Equal_R_R_vector(error_action, Op_psi, Op_psi_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="H_{diag, w=1/2}|0>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi, out_unit, Size(Op_psi), info="H_{diag, w=1/2}|0>")
    CALL Write_Vec(Op_psi_ana, out_unit, Size(Op_psi_ana), info="H_{diag, w=1/2}|0>(Analitical)")
  END IF

  CALL Action(Op_psi, H_diag_half, b_1)
  Op_psi_ana(:) = [ZERO, H_diag_half%Diag_val(2), ZERO]
  CALL Equal_R_R_vector(error_action, Op_psi, Op_psi_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="H_{diag, w=1/2}|1>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi, out_unit, Size(Op_psi), info="H_{diag, w=1/2}|1>")
    CALL Write_Vec(Op_psi_ana, out_unit, Size(Op_psi_ana), info="H_{diag, w=1/2}|1>(Analitical)")
  END IF

  CALL Action(Op_psi, H_diag_half, b_2)
  Op_psi_ana(:) = [ZERO, ZERO, H_diag_half%Diag_val(3)]
  CALL Equal_R_R_vector(error_action, Op_psi, Op_psi_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="H_{diag, w=1/2}|2>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi, out_unit, Size(Op_psi), info="H_{diag, w=1/2}|2>")
    CALL Write_Vec(Op_psi_ana, out_unit, Size(Op_psi_ana), info="H_{diag, w=1/2}|2>(Analitical)")
  END IF

  CALL Action(Op_psi, H_dense_half, b_0)
  Op_psi_ana(:) = [H_dense_half%Dense_val(1,1), ZERO, ZERO]
  CALL Equal_R_R_vector(error_action, Op_psi, Op_psi_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="H_{dense, w=1/2}|0>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi, out_unit, Size(Op_psi), info="H_{dense, w=1/2}|0>")
    CALL Write_Vec(Op_psi_ana, out_unit, Size(Op_psi_ana), info="H_{dense, w=1/2}|0>(Analitical)")
  END IF

  CALL Action(Op_psi, H_dense_half, b_1)
  Op_psi_ana(:) = [ZERO, H_dense_half%Dense_val(2,2), ZERO]
  CALL Equal_R_R_vector(error_action, Op_psi, Op_psi_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="H_{dense, w=1/2}|1>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi, out_unit, Size(Op_psi), info="H_{dense, w=1/2}|1>")
    CALL Write_Vec(Op_psi_ana, out_unit, Size(Op_psi_ana), info="H_{dense, w=1/2}|1>(Analitical)")
  END IF

  CALL Action(Op_psi, H_dense_half, b_2)
  Op_psi_ana(:) = [ZERO, ZERO, H_dense_half%Dense_val(3,3)]
  CALL Equal_R_R_vector(error_action, Op_psi, Op_psi_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="H_{dense, w=1/2}|2>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi, out_unit, Size(Op_psi), info="H_{dense, w=1/2}|2>")
    CALL Write_Vec(Op_psi_ana, out_unit, Size(Op_psi_ana), info="H_{dense, w=1/2}|2>(Analitical)")
  END IF

  CALL Action(Op_psi, H_diag_half, Psi_1D)
  Op_psi_ana(:) = [Coeff_0*H_diag_half%Diag_val(1), Coeff_1*H_diag_half%Diag_val(2), Coeff_2*H_diag_half%Diag_val(3)] / Norm
  CALL Equal_R_R_vector(error_action, Op_psi, Op_psi_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="H_{diag, w=1/2}|Psi_1D>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi*Norm, out_unit, Size(Op_psi), info="H_{diag, w=1/2}|Psi_1D>*Norm")
    CALL Write_Vec(Op_psi_ana*Norm, out_unit, Size(Op_psi_ana), info="H_{diag, w=1/2}|Psi_1D>*Norm(Analitical)")
  END IF

  CALL Action(Op_psi, H_dense_half, Psi_1D)
  Op_psi_ana(:) = [Coeff_0*H_dense_half%Dense_val(1,1), Coeff_1*H_dense_half%Dense_val(2,2), Co&
                  &eff_2*H_dense_half%Dense_val(3,3)] / Norm
  CALL Equal_R_R_vector(error_action, Op_psi, Op_psi_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="H_{dense, w=1/2}|Psi_1D>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi*Norm, out_unit, Size(Op_psi), info="H_{dense, w=1/2}|Psi_1D>*Norm")
    CALL Write_Vec(Op_psi_ana*Norm, out_unit, Size(Op_psi_ana), info="H_{dense, w=1/2}|Psi_1D>*Norm(Analitical)")
  END IF

    !----------------------------------H_{w=3}---------------------------------
  CALL Action(Op_psi, H_diag_3, b_0)
  Op_psi_ana(:) = [H_diag_3%Diag_val(1), ZERO, ZERO]
  CALL Equal_R_R_vector(error_action, Op_psi, Op_psi_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="H_{diag, w=3}|0>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi, out_unit, Size(Op_psi), info="H_{diag, w=3}|0>")
    CALL Write_Vec(Op_psi_ana, out_unit, Size(Op_psi_ana), info="H_{diag, w=3}|0>(Analitical)")
  END IF

  CALL Action(Op_psi, H_diag_3, b_1)
  Op_psi_ana(:) = [ZERO, H_diag_3%Diag_val(2), ZERO]
  CALL Equal_R_R_vector(error_action, Op_psi, Op_psi_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="H_{diag, w=3}|1>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi, out_unit, Size(Op_psi), info="H_{diag, w=3}|1>")
    CALL Write_Vec(Op_psi_ana, out_unit, Size(Op_psi_ana), info="H_{diag, w=3}|1>(Analitical)")
  END IF

  CALL Action(Op_psi, H_diag_3, b_2)
  Op_psi_ana(:) = [ZERO, ZERO, H_diag_3%Diag_val(3)]
  CALL Equal_R_R_vector(error_action, Op_psi, Op_psi_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="H_{diag, w=3}|2>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi, out_unit, Size(Op_psi), info="H_{diag, w=3}|2>")
    CALL Write_Vec(Op_psi_ana, out_unit, Size(Op_psi_ana), info="H_{diag, w=3}|2>(Analitical)")
  END IF

  CALL Action(Op_psi, H_dense_3, b_0)
  Op_psi_ana(:) = [H_dense_3%Dense_val(1,1), ZERO, ZERO]
  CALL Equal_R_R_vector(error_action, Op_psi, Op_psi_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="H_{dense, w=3}|0>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi, out_unit, Size(Op_psi), info="H_{dense, w=3}|0>")
    CALL Write_Vec(Op_psi_ana, out_unit, Size(Op_psi_ana), info="H_{dense, w=3}|0>(Analitical)")
  END IF

  CALL Action(Op_psi, H_dense_3, b_1)
  Op_psi_ana(:) = [ZERO, H_dense_3%Dense_val(2,2), ZERO]
  CALL Equal_R_R_vector(error_action, Op_psi, Op_psi_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="H_{dense, w=3}|1>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi, out_unit, Size(Op_psi), info="H_{dense, w=3}|1>")
    CALL Write_Vec(Op_psi_ana, out_unit, Size(Op_psi_ana), info="H_{dense, w=3}|1>(Analitical)")
  END IF

  CALL Action(Op_psi, H_dense_3, b_2)
  Op_psi_ana(:) = [ZERO, ZERO, H_dense_3%Dense_val(3,3)]
  CALL Equal_R_R_vector(error_action, Op_psi, Op_psi_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="H_{dense, w=3}|2>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi, out_unit, Size(Op_psi), info="H_{dense, w=3}|2>")
    CALL Write_Vec(Op_psi_ana, out_unit, Size(Op_psi_ana), info="H_{dense, w=3}|2>(Analitical)")
  END IF

  CALL Action(Op_psi, H_diag_3, Psi_1D)
  Op_psi_ana(:) = [Coeff_0*H_diag_3%Diag_val(1), Coeff_1*H_diag_3%Diag_val(2), Coeff_2*H_diag_3%Diag_val(3)] / Norm
  CALL Equal_R_R_vector(error_action, Op_psi, Op_psi_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="H_{diag, w=3}|Psi_1D>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi*Norm, out_unit, Size(Op_psi), info="H_{diag, w=3}|Psi_1D>*Norm")
    CALL Write_Vec(Op_psi_ana*Norm, out_unit, Size(Op_psi_ana), info="H_{diag, w=3}|Psi_1D>*Norm(Analitical)")
  END IF

  CALL Action(Op_psi, H_dense_3, Psi_1D)
  Op_psi_ana(:) = [Coeff_0*H_dense_3%Dense_val(1,1), Coeff_1*H_dense_3%Dense_val(2,2), Coeff_2*&
                  &H_dense_3%Dense_val(3,3)] / Norm
  CALL Equal_R_R_vector(error_action, Op_psi, Op_psi_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="H_{dense, w=3}|Psi_1D>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi*Norm, out_unit, Size(Op_psi), info="H_{dense, w=3}|Psi_1D>*Norm")
    CALL Write_Vec(Op_psi_ana*Norm, out_unit, Size(Op_psi_ana), info="H_{dense, w=3}|Psi_1D>*Norm(Analitical)")
  END IF

  !-------------------------x matricies initialization-------------------------
  Mode%w  = HALF
  Mode%m = ONE
  CALL Construct_Operator_1D(x_band_half_1,  "Position", Mode=Mode, Debug=Debug)
  CALL Construct_Operator_1D(x_dense_half_1, "Position", Dense=.TRUE., Mode=Mode, Debug=Debug)

  Mode%w  = THREE
  CALL Construct_Operator_1D(x_band_3_1,  "Position", Mode=Mode, Debug=Debug)
  CALL Construct_Operator_1D(x_dense_3_1, "Position", Dense=.TRUE., Mode=Mode, Debug=Debug)

  Mode%m = FOUR
  CALL Construct_Operator_1D(x_band_3_4,  "Position", Mode=Mode, Debug=Debug)
  CALL Construct_Operator_1D(x_dense_3_4, "Position", Dense=.TRUE., Mode=Mode, Debug=Debug)


  !----------------------------Testing the x actions---------------------------
    !------------------------------x_{w=1/2, m=1}------------------------------
  CALL Action(Op_psi, x_band_half_1, b_0)
  Op_psi_ana(:) = [ZERO, x_band_half_1%Band_val(1,1), ZERO]
  CALL Equal_R_R_vector(error_action, Op_psi, Op_psi_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="x_{band, w=1/2, m=1}|0>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi, out_unit, Size(Op_psi), info="x_{band, w=1/2, m=1}|0>")
    CALL Write_Vec(Op_psi_ana, out_unit, Size(Op_psi_ana), info="x_{band, w=1/2, m=1}|0>(Analitical)")
  END IF

  CALL Action(Op_psi, x_band_half_1, b_1)
  Op_psi_ana(:) = [x_band_half_1%Band_val(2,3), ZERO, x_band_half_1%Band_val(2,1)]
  CALL Equal_R_R_vector(error_action, Op_psi, Op_psi_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="x_{band, w=1/2, m=1}|1>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi, out_unit, Size(Op_psi), info="x_{band, w=1/2, m=1}|1>")
    CALL Write_Vec(Op_psi_ana, out_unit, Size(Op_psi_ana), info="x_{band, w=1/2, m=1}|1>(Analitical)")
  END IF

  CALL Action(Op_psi, x_band_half_1, b_2)
  Op_psi_ana(:) = [ZERO, x_band_half_1%Band_val(3,3), ZERO]
  CALL Equal_R_R_vector(error_action, Op_psi, Op_psi_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="x_{band, w=1/2, m=1}|2>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi, out_unit, Size(Op_psi), info="x_{band, w=1/2, m=1}|2>")
    CALL Write_Vec(Op_psi_ana, out_unit, Size(Op_psi_ana), info="x_{band, w=1/2, m=1}|2>(Analitical)")
  END IF

  CALL Action(Op_psi, x_dense_half_1, b_0)
  Op_psi_ana(:) = [ZERO, x_dense_half_1%Dense_val(2,1), ZERO]
  CALL Equal_R_R_vector(error_action, Op_psi, Op_psi_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="x_{dense, w=1/2, m=1}|0>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi, out_unit, Size(Op_psi), info="x_{dense, w=1/2, m=1}|0>")
    CALL Write_Vec(Op_psi_ana, out_unit, Size(Op_psi_ana), info="x_{dense, w=1/2, m=1}|0>(Analitical)")
  END IF

  CALL Action(Op_psi, x_dense_half_1, b_1)
  Op_psi_ana(:) = [x_dense_half_1%Dense_val(1,2), ZERO, x_dense_half_1%Dense_val(3,2)]
  CALL Equal_R_R_vector(error_action, Op_psi, Op_psi_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="x_{dense, w=1/2, m=1}|1>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi, out_unit, Size(Op_psi), info="x_{dense, w=1/2, m=1}|1>")
    CALL Write_Vec(Op_psi_ana, out_unit, Size(Op_psi_ana), info="x_{dense, w=1/2, m=1}|1>(Analitical)")
  END IF

  CALL Action(Op_psi, x_dense_half_1, b_2)
  Op_psi_ana(:) = [ZERO, x_dense_half_1%Dense_val(2,3), ZERO]
  CALL Equal_R_R_vector(error_action, Op_psi, Op_psi_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="x_{dense, w=1/2, m=1}|2>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi, out_unit, Size(Op_psi), info="x_{dense, w=1/2, m=1}|2>")
    CALL Write_Vec(Op_psi_ana, out_unit, Size(Op_psi_ana), info="x_{dense, w=1/2, m=1}|2>(Analitical)")
  END IF

  CALL Action(Op_psi, x_band_half_1, Psi_1D)
  Op_psi_ana(:) = [Coeff_1*x_band_half_1%Band_val(2,3), Coeff_0*x_band_half_1%Band_val(1,1) + C&
                  &oeff_2*x_band_half_1%Band_val(3,3), Coeff_1*x_band_half_1%Band_val(2,1)] / Norm
  CALL Equal_R_R_vector(error_action, Op_psi, Op_psi_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="x_{band, w=1/2, m=1}|Psi_1D>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi*Norm, out_unit, Size(Op_psi), info="x_{band, w=1/2, m=1}|Psi_1D>*Norm")
    CALL Write_Vec(Op_psi_ana*Norm, out_unit, Size(Op_psi_ana), info="x_{band, w=1/2, m=1}|Psi_1D>*Norm(Analitical)")
  END IF

  CALL Action(Op_psi, x_dense_half_1, Psi_1D)
  Op_psi_ana(:) = [Coeff_1*x_dense_half_1%Dense_val(1,2), Coeff_0*x_dense_half_1%Dense_val(2,1)&
                 & + Coeff_2*x_dense_half_1%Dense_val(2,3), Coeff_1*x_dense_half_1%Dense_val(3,2)] / Norm
  CALL Equal_R_R_vector(error_action, Op_psi, Op_psi_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="x_{dense, w=1/2, m=1}|Psi_1D>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi*Norm, out_unit, Size(Op_psi), info="x_{dense, w=1/2, m=1}|Psi_1D>*Norm")
    CALL Write_Vec(Op_psi_ana*Norm, out_unit, Size(Op_psi_ana), info="x_{dense, w=1/2, m=1}|Psi_1D>*Norm(Analitical)")
  END IF

    !-------------------------------x_{w=3, m=1}-------------------------------
  CALL Action(Op_psi, x_band_3_1, b_0)
  Op_psi_ana(:) = [ZERO, x_band_3_1%Band_val(1,1), ZERO]
  CALL Equal_R_R_vector(error_action, Op_psi, Op_psi_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="x_{band, w=3, m=1}|0>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi, out_unit, Size(Op_psi), info="x_{band, w=3, m=1}|0>")
    CALL Write_Vec(Op_psi_ana, out_unit, Size(Op_psi_ana), info="x_{band, w=3, m=1}|0>(Analitical)")
  END IF

  CALL Action(Op_psi, x_band_3_1, b_1)
  Op_psi_ana(:) = [x_band_3_1%Band_val(2,3), ZERO, x_band_3_1%Band_val(2,1)]
  CALL Equal_R_R_vector(error_action, Op_psi, Op_psi_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="x_{band, w=3, m=1}|1>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi, out_unit, Size(Op_psi), info="x_{band, w=3, m=1}|1>")
    CALL Write_Vec(Op_psi_ana, out_unit, Size(Op_psi_ana), info="x_{band, w=3, m=1}|1>(Analitical)")
  END IF

  CALL Action(Op_psi, x_band_3_1, b_2)
  Op_psi_ana(:) = [ZERO, x_band_3_1%Band_val(3,3), ZERO]
  CALL Equal_R_R_vector(error_action, Op_psi, Op_psi_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="x_{band, w=3, m=1}|2>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi, out_unit, Size(Op_psi), info="x_{band, w=3, m=1}|2>")
    CALL Write_Vec(Op_psi_ana, out_unit, Size(Op_psi_ana), info="x_{band, w=3, m=1}|2>(Analitical)")
  END IF

  CALL Action(Op_psi, x_dense_3_1, b_0)
  Op_psi_ana(:) = [ZERO, x_dense_3_1%Dense_val(2,1), ZERO]
  CALL Equal_R_R_vector(error_action, Op_psi, Op_psi_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="x_{dense, w=3, m=1}|0>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi, out_unit, Size(Op_psi), info="x_{dense, w=3, m=1}|0>")
    CALL Write_Vec(Op_psi_ana, out_unit, Size(Op_psi_ana), info="x_{dense, w=3, m=1}|0>(Analitical)")
  END IF

  CALL Action(Op_psi, x_dense_3_1, b_1)
  Op_psi_ana(:) = [x_dense_3_1%Dense_val(1,2), ZERO, x_dense_3_1%Dense_val(3,2)]
  CALL Equal_R_R_vector(error_action, Op_psi, Op_psi_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="x_{dense, w=3, m=1}|1>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi, out_unit, Size(Op_psi), info="x_{dense, w=3, m=1}|1>")
    CALL Write_Vec(Op_psi_ana, out_unit, Size(Op_psi_ana), info="x_{dense, w=3, m=1}|1>(Analitical)")
  END IF

  CALL Action(Op_psi, x_dense_3_1, b_2)
  Op_psi_ana(:) = [ZERO, x_dense_3_1%Dense_val(2,3), ZERO]
  CALL Equal_R_R_vector(error_action, Op_psi, Op_psi_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="x_{dense, w=3, m=1}|2>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi, out_unit, Size(Op_psi), info="x_{dense, w=3, m=1}|2>")
    CALL Write_Vec(Op_psi_ana, out_unit, Size(Op_psi_ana), info="x_{dense, w=3, m=1}|2>(Analitical)")
  END IF

  CALL Action(Op_psi, x_band_3_1, Psi_1D)
  Op_psi_ana(:) = [Coeff_1*x_band_3_1%Band_val(2,3), Coeff_0*x_band_3_1%Band_val(1,1) + Coeff_2&
                  &*x_band_3_1%Band_val(3,3), Coeff_1*x_band_3_1%Band_val(2,1)] / Norm
  CALL Equal_R_R_vector(error_action, Op_psi, Op_psi_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="x_{band, w=3, m=1}|Psi_1D>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi*Norm, out_unit, Size(Op_psi), info="x_{band, w=3, m=1}|Psi_1D>*Norm")
    CALL Write_Vec(Op_psi_ana*Norm, out_unit, Size(Op_psi_ana), info="x_{band, w=3, m=1}|Psi_1D>*Norm(Analitical)")
  END IF

  CALL Action(Op_psi, x_dense_3_1, Psi_1D)
  Op_psi_ana(:) = [Coeff_1*x_dense_3_1%Dense_val(1,2), Coeff_0*x_dense_3_1%Dense_val(2,1) + Coe&
                  &ff_2*x_dense_3_1%Dense_val(2,3), Coeff_1*x_dense_3_1%Dense_val(3,2)] / Norm
  CALL Equal_R_R_vector(error_action, Op_psi, Op_psi_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="x_{dense, w=3, m=1}|Psi_1D>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi*Norm, out_unit, Size(Op_psi), info="x_{dense, w=3, m=1}|Psi_1D>*Norm")
    CALL Write_Vec(Op_psi_ana*Norm, out_unit, Size(Op_psi_ana), info="x_{dense, w=3, m=1}|Psi_1D>*Norm(Analitical)")
  END IF

    !-------------------------------x_{w=3, m=4}-------------------------------
  CALL Action(Op_psi, x_band_3_4, b_0)
  Op_psi_ana(:) = [ZERO, x_band_3_4%Band_val(1,1), ZERO]
  CALL Equal_R_R_vector(error_action, Op_psi, Op_psi_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="x_{band, w=3, m=4}|0>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi, out_unit, Size(Op_psi), info="x_{band, w=3, m=4}|0>")
    CALL Write_Vec(Op_psi_ana, out_unit, Size(Op_psi_ana), info="x_{band, w=3, m=4}|0>(Analitical)")
  END IF

  CALL Action(Op_psi, x_band_3_4, b_1)
  Op_psi_ana(:) = [x_band_3_4%Band_val(2,3), ZERO, x_band_3_4%Band_val(2,1)]
  CALL Equal_R_R_vector(error_action, Op_psi, Op_psi_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="x_{band, w=3, m=4}|1>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi, out_unit, Size(Op_psi), info="x_{band, w=3, m=4}|1>")
    CALL Write_Vec(Op_psi_ana, out_unit, Size(Op_psi_ana), info="x_{band, w=3, m=4}|1>(Analitical)")
  END IF

  CALL Action(Op_psi, x_band_3_4, b_2)
  Op_psi_ana(:) = [ZERO, x_band_3_4%Band_val(3,3), ZERO]
  CALL Equal_R_R_vector(error_action, Op_psi, Op_psi_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="x_{band, w=3, m=4}|2>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi, out_unit, Size(Op_psi), info="x_{band, w=3, m=4}|2>")
    CALL Write_Vec(Op_psi_ana, out_unit, Size(Op_psi_ana), info="x_{band, w=3, m=4}|2>(Analitical)")
  END IF

  CALL Action(Op_psi, x_dense_3_4, b_0)
  Op_psi_ana(:) = [ZERO, x_dense_3_4%Dense_val(2,1), ZERO]
  CALL Equal_R_R_vector(error_action, Op_psi, Op_psi_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="x_{dense, w=3, m=4}|0>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi, out_unit, Size(Op_psi), info="x_{dense, w=3, m=4}|0>")
    CALL Write_Vec(Op_psi_ana, out_unit, Size(Op_psi_ana), info="x_{dense, w=3, m=4}|0>(Analitical)")
  END IF

  CALL Action(Op_psi, x_dense_3_4, b_1)
  Op_psi_ana(:) = [x_dense_3_4%Dense_val(1,2), ZERO, x_dense_3_4%Dense_val(3,2)]
  CALL Equal_R_R_vector(error_action, Op_psi, Op_psi_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="x_{dense, w=3, m=4}|1>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi, out_unit, Size(Op_psi), info="x_{dense, w=3, m=4}|1>")
    CALL Write_Vec(Op_psi_ana, out_unit, Size(Op_psi_ana), info="x_{dense, w=3, m=4}|1>(Analitical)")
  END IF

  CALL Action(Op_psi, x_dense_3_4, b_2)
  Op_psi_ana(:) = [ZERO, x_dense_3_4%Dense_val(2,3), ZERO]
  CALL Equal_R_R_vector(error_action, Op_psi, Op_psi_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="x_{dense, w=3, m=4}|2>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi, out_unit, Size(Op_psi), info="x_{dense, w=3, m=4}|2>")
    CALL Write_Vec(Op_psi_ana, out_unit, Size(Op_psi_ana), info="x_{dense, w=3, m=4}|2>(Analitical)")
  END IF

  CALL Action(Op_psi, x_band_3_4, Psi_1D)
  Op_psi_ana(:) = [Coeff_1*x_band_3_4%Band_val(2,3), Coeff_0*x_band_3_4%Band_val(1,1) + Coeff_2&
                  &*x_band_3_4%Band_val(3,3), Coeff_1*x_band_3_4%Band_val(2,1)] / Norm
  CALL Equal_R_R_vector(error_action, Op_psi, Op_psi_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="x_{band, w=3, m=4}|Psi_1D>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi*Norm, out_unit, Size(Op_psi), info="x_{band, w=3, m=4}|Psi_1D>*Norm")
    CALL Write_Vec(Op_psi_ana*Norm, out_unit, Size(Op_psi_ana), info="x_{band, w=3, m=4}|Psi_1D>*Norm(Analitical)")
  END IF

  CALL Action(Op_psi, x_dense_3_4, Psi_1D)
  Op_psi_ana(:) = [Coeff_1*x_dense_3_4%Dense_val(1,2), Coeff_0*x_dense_3_4%Dense_val(2,1) + Coe&
                  &ff_2*x_dense_3_4%Dense_val(2,3), Coeff_1*x_dense_3_4%Dense_val(3,2)] / Norm
  CALL Equal_R_R_vector(error_action, Op_psi, Op_psi_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="x_{dense, w=3, m=4}|Psi_1D>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi*Norm, out_unit, Size(Op_psi), info="x_{dense, w=3, m=4}|Psi_1D>*Norm")
    CALL Write_Vec(Op_psi_ana*Norm, out_unit, Size(Op_psi_ana), info="x_{dense, w=3, m=4}|Psi_1D>*Norm(Analitical)")
  END IF


  !-------------------------N matricies initialization-------------------------
  CALL Construct_Operator_1D(N_diag,  "Nb_photons", Mode=Mode, Debug=Debug)
  CALL Construct_Operator_1D(N_dense, "Nb_photons", Dense=.TRUE., Mode=Mode, Debug=Debug)


  !----------------------------Testing the N actions---------------------------
  CALL Action(Op_psi, N_diag, b_0)
  Op_psi_ana(:) = ZERO
  CALL Equal_R_R_vector(error_action, Op_psi, Op_psi_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="N_{diag}|0>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi, out_unit, Size(Op_psi), info="N_{diag}|0>")
    CALL Write_Vec(Op_psi_ana, out_unit, Size(Op_psi_ana), info="N_{diag}|0>(Analitical)")
  END IF

  CALL Action(Op_psi, N_diag, b_1)
  Op_psi_ana(:) = [ZERO, N_diag%Diag_val(2), ZERO]
  CALL Equal_R_R_vector(error_action, Op_psi, Op_psi_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="N_{diag}|1>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi, out_unit, Size(Op_psi), info="N_{diag}|1>")
    CALL Write_Vec(Op_psi_ana, out_unit, Size(Op_psi_ana), info="N_{diag}|1>(Analitical)")
  END IF

  CALL Action(Op_psi, N_diag, b_2)
  Op_psi_ana(:) = [ZERO, ZERO, N_diag%Diag_val(3)]
  CALL Equal_R_R_vector(error_action, Op_psi, Op_psi_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="N_{diag}|2>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi, out_unit, Size(Op_psi), info="N_{diag}|2>")
    CALL Write_Vec(Op_psi_ana, out_unit, Size(Op_psi_ana), info="N_{diag}|2>(Analitical)")
  END IF

  CALL Action(Op_psi, N_dense, b_0)
  Op_psi_ana(:) = ZERO
  CALL Equal_R_R_vector(error_action, Op_psi, Op_psi_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="N_{dense}|0>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi, out_unit, Size(Op_psi), info="N_{dense}|0>")
    CALL Write_Vec(Op_psi_ana, out_unit, Size(Op_psi_ana), info="N_{dense}|0>(Analitical)")
  END IF

  CALL Action(Op_psi, N_dense, b_1)
  Op_psi_ana(:) = [ZERO, N_dense%Dense_val(2,2), ZERO]
  CALL Equal_R_R_vector(error_action, Op_psi, Op_psi_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="N_{dense}|1>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi, out_unit, Size(Op_psi), info="N_{dense}|1>")
    CALL Write_Vec(Op_psi_ana, out_unit, Size(Op_psi_ana), info="N_{dense}|1>(Analitical)")
  END IF

  CALL Action(Op_psi, N_dense, b_2)
  Op_psi_ana(:) = [ZERO, ZERO, N_dense%Dense_val(3,3)]
  CALL Equal_R_R_vector(error_action, Op_psi, Op_psi_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="N_{dense}|2>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi, out_unit, Size(Op_psi), info="N_{dense}|2>")
    CALL Write_Vec(Op_psi_ana, out_unit, Size(Op_psi_ana), info="N_{dense}|2>(Analitical)")
  END IF

  CALL Action(Op_psi, N_diag, Psi_1D)
  Op_psi_ana(:) = [ZERO, Coeff_1*N_diag%Diag_val(2), Coeff_2*N_diag%Diag_val(3)] / Norm
  CALL Equal_R_R_vector(error_action, Op_psi, Op_psi_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="N_{diag}|Psi_1D>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi*Norm, out_unit, Size(Op_psi), info="N_{diag}|Psi_1D>*Norm")
    CALL Write_Vec(Op_psi_ana*Norm, out_unit, Size(Op_psi_ana), info="N_{diag}|Psi_1D>*Norm(Analitical)")
  END IF

  CALL Action(Op_psi, N_dense, Psi_1D)
  Op_psi_ana(:) = [ZERO, Coeff_1*N_dense%Dense_val(2,2), Coeff_2*N_dense%Dense_val(3,3)] / Norm
  CALL Equal_R_R_vector(error_action, Op_psi, Op_psi_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="N_{dense}|Psi_1D>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi*Norm, out_unit, Size(Op_psi), info="N_{dense}|Psi_1D>*Norm")
    CALL Write_Vec(Op_psi_ana*Norm, out_unit, Size(Op_psi_ana), info="N_{dense}|Psi_1D>*Norm(Analitical)")
  END IF

  
  CALL Finalize_Test(test_action)
  
  
  CONTAINS


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
