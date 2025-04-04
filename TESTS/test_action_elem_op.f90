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
PROGRAM test_action_elem_op
  !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
  USE QDUtil_m
  USE Cavity_mode_m
  USE Tests_m
  USE Algebra_m
  USE Elem_op_m
  IMPLICIT NONE

  logical                       :: Debug = .FALSE.

  TYPE(Cavity_mode_t)           :: Mode ! OLD                                                                                    ! The well construction of the Operator's matrices is assumed checked by the dedicated test, so hard-coded references matrices are not needed here

  TYPE(Operator_1D_t)           :: H_diag_half! OLD                                                                                   ! Nomenclature : H_<shape>_<eigenpulsation>
  TYPE(Operator_1D_t)           :: H_dense_half! OLD
  TYPE(Operator_1D_t)           :: H_diag_3! OLD
  TYPE(Operator_1D_t)           :: H_dense_3! OLD

  TYPE(Operator_1D_t)           :: x_band_half_1 ! OLD                                                                                ! Nomenclature : x_<shape>_<eigenpulsation>_<mass>
  TYPE(Operator_1D_t)           :: x_dense_half_1! OLD
  TYPE(Operator_1D_t)           :: x_band_3_1! OLD
  TYPE(Operator_1D_t)           :: x_dense_3_1! OLD
  TYPE(Operator_1D_t)           :: x_band_3_4! OLD
  TYPE(Operator_1D_t)           :: x_dense_3_4! OLD

  TYPE(Operator_1D_t)           :: N_diag          ! OLD                                                                              ! Nomenclature : N_<shape>
  TYPE(Operator_1D_t)           :: N_dense! OLD

  TYPE(Elem_op_t)               :: EO1_diag_1
  TYPE(Elem_op_t)               :: EO1_dense_1
  TYPE(Elem_op_t)               :: EO2_diag_tenth
  TYPE(Elem_op_t)               :: EO2_dense_tenth
  TYPE(Elem_op_t)               :: EO3_band_1
  TYPE(Elem_op_t)               :: EO3_dense_1
  TYPE(Elem_op_t)               :: EO4_band_pi
  TYPE(Elem_op_t)               :: EO4_dense_pi

  real(kind=Rkind)              :: b_0(3)  = [ONE, ZERO, ZERO]                                                                   ! three vectors of the HO basis set |0>, |1>, |2> 
  real(kind=Rkind)              :: b_1(3)  = [ZERO, ONE, ZERO]
  real(kind=Rkind)              :: b_2(3)  = [ZERO, ZERO, ONE]
  real(kind=Rkind)              :: Coeff_0_real = ONE
  real(kind=Rkind)              :: Coeff_1_real = HALF
  real(kind=Rkind)              :: Coeff_2_real = PI
  real(kind=Rkind)              :: Psi_1D_R1_real(3)                                                                             ! an any vector representing the excitation state/wavefunction = a linear combination of the basis functions from the canonical basis set over \mathbb{R} /!\ Not normalized yet !
  real(kind=Rkind)              :: Op_psi_real(3)                                                                                ! the resulting vector from the action of a 1D operator upon Psi_1D_R1_real
  real(kind=Rkind)              :: Op_psi_real_ana(3)                                                                            ! the analytical action = hard-coded reference for comparison

  real(kind=Rkind)              :: Norm                                                                                          ! SQRT(Coeff_0_real**2 + Coeff_1_real**2 + Coeff_2_real**2) or using MOD(Coeff_i_complex)**2 instead for complex Psi

  complex(kind=Rkind)           :: Coeff_0_complex = ONE*EYE
  complex(kind=Rkind)           :: Coeff_1_complex = HALF
  complex(kind=Rkind)           :: Coeff_2_complex = PI*EYE
  complex(kind=Rkind)           :: Psi_1D_R1_complex(3)                                                                          ! an any vector representing the excitation state/wavefunction = a linear combination of the basis functions from the canonical basis set over \mathbb{R} BUT with complexes expansion coefficients to make complexe WF. N.B. it is the same vector as the real one but some coeff are in the complex space. /!\ Not normalized yet !
  complex(kind=Rkind)           :: Op_psi_complex(3)                                                                             ! the resulting vector from the action of a 1D operator upon Psi_1D_R1_complex
  complex(kind=Rkind)           :: Op_psi_complex_ana(3)                                                                         ! the analytical action = hard-coded reference for comparison

  TYPE(test_t)                  :: test_action
  logical                       :: error_action = .FALSE.


  !-----------------------------Test initialization----------------------------
  CALL Initialize_Test(test_action, test_name="OUT/test_file_actn_op_1D")
  

  !-------------------------Wavefunction initialization (real)------------------------
  Psi_1D_R1_real(:) = [Coeff_0_real, Coeff_1_real, Coeff_2_real]
  CALL Norm_of(Norm, Psi_1D_R1_real)
  CALL Normalize(Psi_1D_R1_real)
  IF (Debug) THEN
    WRITE(out_unit,*)
    WRITE(out_unit,*) "----------------The cavity wavefunction has been initialized as :---------------"
    CALL Write_Vec(Psi_1D_R1_real, out_unit, Size(Psi_1D_R1_real), info="Psi_1D_R1_real")
    WRITE(out_unit,*) "-----------------------------End cavity wavefunction----------------------------"
  END IF


  !-------------------------Wavefunction initialization (complex)------------------------
  Psi_1D_R1_complex(:) = [Coeff_0_complex, Coeff_1_complex, Coeff_2_complex]                                                     ! uses the same basis set as the real WF, but the expansion coefficients are complexes
  CALL Norm_of(Norm, Psi_1D_R1_complex)
  CALL Normalize(Psi_1D_R1_complex)
  IF (Debug) THEN
    WRITE(out_unit,*)
    WRITE(out_unit,*) "----------------The cavity wavefunction has been initialized as :---------------"
    CALL Write_Vec(Psi_1D_R1_complex, out_unit, Size(Psi_1D_R1_complex), info="Psi_1D_R1_complex")
    WRITE(out_unit,*) "-----------------------------End cavity wavefunction----------------------------"
  END IF


  !-------------------------Cavity mode initialization OLD-------------------------
  CALL Read_cavity_mode(Mode, nio=in_unit)
  Mode%Nb = 3
  IF (Debug) THEN
    WRITE(out_unit,*)
    WRITE(out_unit,*) "--------------Cavity mode constructed by MolecCav_Read_cavity_mode--------------"
    CALL Write_cavity_mode(Mode)
    WRITE(out_unit,*) "------------End Cavity mode constructed by MolecCav_Read_cavity_mode------------"
  END IF

  
  !-------------------------H matricies initialization OLD-------------------------
  Mode%w  = HALF
  CALL Construct_Operator_1D(H_diag_half,  "Hamiltonian", Mode=Mode, Debug=Debug) !OLD
  CALL Construct_Operator_1D(H_dense_half, "Hamiltonian", Dense=.TRUE., Mode=Mode, Debug=Debug) !OLD

  Mode%w  = THREE
  CALL Construct_Operator_1D(H_diag_3,  "Hamiltonian", Mode=Mode, Debug=Debug)
  CALL Construct_Operator_1D(H_dense_3, "Hamiltonian", Dense=.TRUE., Mode=Mode, Debug=Debug)


  !-------------------------Diagonal matricies initialization-------------------------
  CALL MolecCav_Construct_Elem_op(EO1_diag_1,      Coeff=ONE,      Operator_type="Hamiltonian", Dense=.FALSE., Debug_opt=Debug)
  CALL MolecCav_Construct_Elem_op(EO1_dense_1,     Coeff=ONE,      Operator_type="Hamiltonian", Dense=.TRUE.,  Debug_opt=Debug)
  CALL MolecCav_Construct_Elem_op(EO2_diag_tenth,  Coeff=ONETENTH, Operator_type="Hamiltonian", Dense=.FALSE., Debug_opt=Debug)
  CALL MolecCav_Construct_Elem_op(EO2_dense_tenth, Coeff=ONETENTH, Operator_type="Hamiltonian", Dense=.TRUE.,  Debug_opt=Debug)
  CALL MolecCav_Construct_Elem_op(EO3_band_1,      Coeff=ONE,      Operator_type="Position",    Dense=.FALSE., Debug_opt=Debug)
  CALL MolecCav_Construct_Elem_op(EO3_dense_1,     Coeff=ONE,      Operator_type="Position",    Dense=.TRUE.,  Debug_opt=Debug)
  CALL MolecCav_Construct_Elem_op(EO4_band_pi,     Coeff=PI,       Operator_type="Position",    Dense=.FALSE., Debug_opt=Debug)
  CALL MolecCav_Construct_Elem_op(EO4_dense_pi,    Coeff=PI,       Operator_type="Position",    Dense=.TRUE.,  Debug_opt=Debug)
  
  !----------------------------Testing the actions of diagonal matrices (real)---------------------------
    !---------------------------------Coeff = 1--------------------------------
  CALL MolecCav_Construct_op_psi_ana_real(Op_psi_real_ana, EO1_diag_1, ONE, b_0, Debug) ! same analytical result as with E01_dense_1

  CALL Action(Op_psi_real, EO1_diag_1, b_0, Debug=Debug)
  CALL Equal_tensor(error_action, Op_psi_real, Op_psi_real_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="\hat{O}_{diag, 1}|0>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_real, out_unit, Size(Op_psi_real), info="\hat{O}_{diag, 1}|0>")
    CALL Write_Vec(Op_psi_real_ana, out_unit, Size(Op_psi_real_ana), info="\hat{O}_{diag, 1}|0>(Analitical)")
  END IF

  CALL Action(Op_psi_real, EO1_dense_1, b_0, Debug=Debug)
  CALL Equal_tensor(error_action, Op_psi_real, Op_psi_real_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="\hat{O}_{dense, 1}|0>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_real, out_unit, Size(Op_psi_real), info="\hat{O}_{dense, 1}|0>")
    CALL Write_Vec(Op_psi_real_ana, out_unit, Size(Op_psi_real_ana), info="\hat{O}_{dense, 1}|0>(Analitical)")
  END IF

  Op_psi_complex_ana = Op_psi_real_ana*EYE                                                               ! make as if the b_0 was a vector of the canonical basis set on \mathbb{C} i.e. i*b_0
  CALL Action(Op_psi_complex, EO1_diag_1, EYE*b_0, Debug=Debug)
  ! we are here and the contains/var declaration/WF initialization is finished. the test works so far

  CALL MolecCav_Construct_op_psi_ana_real(Op_psi_real_ana,    EO1_diag_1, ONE, b_1,               Debug) ! same analytical result as with E01_dense_1
  Op_psi_complex_ana = Op_psi_real_ana*EYE                                                               ! make as if the b_1 was a vector of the canonical basis set on \mathbb{C} i.e. i*b_1
  CALL MolecCav_Construct_op_psi_ana_real(Op_psi_real_ana,    EO1_diag_1, ONE, b_1,               Debug) ! same analytical result as with E01_dense_1
  Op_psi_complex_ana = Op_psi_real_ana*EYE                                                               ! make as if the b_2 was a vector of the canonical basis set on \mathbb{C} i.e. i*b_2
  CALL MolecCav_Construct_op_psi_ana_real(Op_psi_real_ana,    EO1_diag_1, ONE, Psi_1D_R1_real,    Debug) ! same analytical result as with E01_dense_1
  CALL MolecCav_Construct_op_psi_ana_complex(Op_psi_complex_ana, EO1_diag_1, ONE, Psi_1D_R1_complex, Debug) ! no easy product here because not all the coefficients are imaginary for the LC

  CALL Action(Op_psi_real, EO1_diag_1, b_0)
  Op_psi_real_ana(:) = [H_diag_half%Diag_val(1), ZERO, ZERO]

  CALL Action(Op_psi_real, H_diag_half, b_1)
  Op_psi_real_ana(:) = [ZERO, H_diag_half%Diag_val(2), ZERO]
  CALL Equal_tensor(error_action, Op_psi_real, Op_psi_real_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="H_{diag, w=1/2}|1>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_real, out_unit, Size(Op_psi_real), info="H_{diag, w=1/2}|1>")
    CALL Write_Vec(Op_psi_real_ana, out_unit, Size(Op_psi_real_ana), info="H_{diag, w=1/2}|1>(Analitical)")
  END IF

  CALL Action(Op_psi_real, H_diag_half, b_2)
  Op_psi_real_ana(:) = [ZERO, ZERO, H_diag_half%Diag_val(3)]
  CALL Equal_tensor(error_action, Op_psi_real, Op_psi_real_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="H_{diag, w=1/2}|2>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_real, out_unit, Size(Op_psi_real), info="H_{diag, w=1/2}|2>")
    CALL Write_Vec(Op_psi_real_ana, out_unit, Size(Op_psi_real_ana), info="H_{diag, w=1/2}|2>(Analitical)")
  END IF

  CALL Action(Op_psi_real, H_dense_half, b_0)
  Op_psi_real_ana(:) = [H_dense_half%Dense_val(1,1), ZERO, ZERO]
  CALL Equal_tensor(error_action, Op_psi_real, Op_psi_real_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="H_{dense, w=1/2}|0>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_real, out_unit, Size(Op_psi_real), info="H_{dense, w=1/2}|0>")
    CALL Write_Vec(Op_psi_real_ana, out_unit, Size(Op_psi_real_ana), info="H_{dense, w=1/2}|0>(Analitical)")
  END IF

  CALL Action(Op_psi_real, H_dense_half, b_1)
  Op_psi_real_ana(:) = [ZERO, H_dense_half%Dense_val(2,2), ZERO]
  CALL Equal_tensor(error_action, Op_psi_real, Op_psi_real_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="H_{dense, w=1/2}|1>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_real, out_unit, Size(Op_psi_real), info="H_{dense, w=1/2}|1>")
    CALL Write_Vec(Op_psi_real_ana, out_unit, Size(Op_psi_real_ana), info="H_{dense, w=1/2}|1>(Analitical)")
  END IF

  CALL Action(Op_psi_real, H_dense_half, b_2)
  Op_psi_real_ana(:) = [ZERO, ZERO, H_dense_half%Dense_val(3,3)]
  CALL Equal_tensor(error_action, Op_psi_real, Op_psi_real_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="H_{dense, w=1/2}|2>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_real, out_unit, Size(Op_psi_real), info="H_{dense, w=1/2}|2>")
    CALL Write_Vec(Op_psi_real_ana, out_unit, Size(Op_psi_real_ana), info="H_{dense, w=1/2}|2>(Analitical)")
  END IF

  CALL Action(Op_psi_real, H_diag_half, Psi_1D_R1_real)
  Op_psi_real_ana(:) = [Coeff_0_real*H_diag_half%Diag_val(1), Coeff_1_real*H_diag_half%Diag_val(2), Coeff_2_real*H_diag_half%Diag&
  &_val(3)] / Norm
  CALL Equal_tensor(error_action, Op_psi_real, Op_psi_real_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="H_{diag, w=1/2}|Psi_1D_R1_real>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_real*Norm, out_unit, Size(Op_psi_real), info="H_{diag, w=1/2}|Psi_1D_R1_real>*Norm")
    CALL Write_Vec(Op_psi_real_ana*Norm, out_unit, Size(Op_psi_real_ana), info="H_{diag, w=1/2}|Psi_1D_R1_real>*Norm(Analitical)")
  END IF

  CALL Action(Op_psi_real, H_dense_half, Psi_1D_R1_real)
  Op_psi_real_ana(:) = [Coeff_0_real*H_dense_half%Dense_val(1,1), Coeff_1_real*H_dense_half%Dense_val(2,2), Co&
                  &eff_2_real*H_dense_half%Dense_val(3,3)] / Norm
  CALL Equal_tensor(error_action, Op_psi_real, Op_psi_real_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="H_{dense, w=1/2}|Psi_1D_R1_real>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_real*Norm, out_unit, Size(Op_psi_real), info="H_{dense, w=1/2}|Psi_1D_R1_real>*Norm")
    CALL Write_Vec(Op_psi_real_ana*Norm, out_unit, Size(Op_psi_real_ana), info="H_{dense, w=1/2}|Psi_1D_R1_real>*Norm(Analitical)")
  END IF

    !----------------------------------H_{w=3}---------------------------------
  CALL Action(Op_psi_real, H_diag_3, b_0)
  Op_psi_real_ana(:) = [H_diag_3%Diag_val(1), ZERO, ZERO]
  CALL Equal_tensor(error_action, Op_psi_real, Op_psi_real_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="H_{diag, w=3}|0>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_real, out_unit, Size(Op_psi_real), info="H_{diag, w=3}|0>")
    CALL Write_Vec(Op_psi_real_ana, out_unit, Size(Op_psi_real_ana), info="H_{diag, w=3}|0>(Analitical)")
  END IF

  CALL Action(Op_psi_real, H_diag_3, b_1)
  Op_psi_real_ana(:) = [ZERO, H_diag_3%Diag_val(2), ZERO]
  CALL Equal_tensor(error_action, Op_psi_real, Op_psi_real_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="H_{diag, w=3}|1>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_real, out_unit, Size(Op_psi_real), info="H_{diag, w=3}|1>")
    CALL Write_Vec(Op_psi_real_ana, out_unit, Size(Op_psi_real_ana), info="H_{diag, w=3}|1>(Analitical)")
  END IF

  CALL Action(Op_psi_real, H_diag_3, b_2)
  Op_psi_real_ana(:) = [ZERO, ZERO, H_diag_3%Diag_val(3)]
  CALL Equal_tensor(error_action, Op_psi_real, Op_psi_real_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="H_{diag, w=3}|2>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_real, out_unit, Size(Op_psi_real), info="H_{diag, w=3}|2>")
    CALL Write_Vec(Op_psi_real_ana, out_unit, Size(Op_psi_real_ana), info="H_{diag, w=3}|2>(Analitical)")
  END IF

  CALL Action(Op_psi_real, H_dense_3, b_0)
  Op_psi_real_ana(:) = [H_dense_3%Dense_val(1,1), ZERO, ZERO]
  CALL Equal_tensor(error_action, Op_psi_real, Op_psi_real_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="H_{dense, w=3}|0>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_real, out_unit, Size(Op_psi_real), info="H_{dense, w=3}|0>")
    CALL Write_Vec(Op_psi_real_ana, out_unit, Size(Op_psi_real_ana), info="H_{dense, w=3}|0>(Analitical)")
  END IF

  CALL Action(Op_psi_real, H_dense_3, b_1)
  Op_psi_real_ana(:) = [ZERO, H_dense_3%Dense_val(2,2), ZERO]
  CALL Equal_tensor(error_action, Op_psi_real, Op_psi_real_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="H_{dense, w=3}|1>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_real, out_unit, Size(Op_psi_real), info="H_{dense, w=3}|1>")
    CALL Write_Vec(Op_psi_real_ana, out_unit, Size(Op_psi_real_ana), info="H_{dense, w=3}|1>(Analitical)")
  END IF

  CALL Action(Op_psi_real, H_dense_3, b_2)
  Op_psi_real_ana(:) = [ZERO, ZERO, H_dense_3%Dense_val(3,3)]
  CALL Equal_tensor(error_action, Op_psi_real, Op_psi_real_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="H_{dense, w=3}|2>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_real, out_unit, Size(Op_psi_real), info="H_{dense, w=3}|2>")
    CALL Write_Vec(Op_psi_real_ana, out_unit, Size(Op_psi_real_ana), info="H_{dense, w=3}|2>(Analitical)")
  END IF

  CALL Action(Op_psi_real, H_diag_3, Psi_1D_R1_real)
  Op_psi_real_ana(:) = [Coeff_0_real*H_diag_3%Diag_val(1), Coeff_1_real*H_diag_3%Diag_val(2), Coeff_2_real*H_diag_3%Diag_val(3)] &
  &/ Norm
  CALL Equal_tensor(error_action, Op_psi_real, Op_psi_real_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="H_{diag, w=3}|Psi_1D_R1_real>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_real*Norm, out_unit, Size(Op_psi_real), info="H_{diag, w=3}|Psi_1D_R1_real>*Norm")
    CALL Write_Vec(Op_psi_real_ana*Norm, out_unit, Size(Op_psi_real_ana), info="H_{diag, w=3}|Psi_1D_R1_real>*Norm(Analitical)")
  END IF

  CALL Action(Op_psi_real, H_dense_3, Psi_1D_R1_real)
  Op_psi_real_ana(:) = [Coeff_0_real*H_dense_3%Dense_val(1,1), Coeff_1_real*H_dense_3%Dense_val(2,2), Coeff_2_real*&
                  &H_dense_3%Dense_val(3,3)] / Norm
  CALL Equal_tensor(error_action, Op_psi_real, Op_psi_real_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="H_{dense, w=3}|Psi_1D_R1_real>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_real*Norm, out_unit, Size(Op_psi_real), info="H_{dense, w=3}|Psi_1D_R1_real>*Norm")
    CALL Write_Vec(Op_psi_real_ana*Norm, out_unit, Size(Op_psi_real_ana), info="H_{dense, w=3}|Psi_1D_R1_real>*Norm(Analitical)")
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
  CALL Action(Op_psi_real, x_band_half_1, b_0)
  Op_psi_real_ana(:) = [ZERO, x_band_half_1%Band_val(1,1), ZERO]
  CALL Equal_tensor(error_action, Op_psi_real, Op_psi_real_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="x_{band, w=1/2, m=1}|0>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_real, out_unit, Size(Op_psi_real), info="x_{band, w=1/2, m=1}|0>")
    CALL Write_Vec(Op_psi_real_ana, out_unit, Size(Op_psi_real_ana), info="x_{band, w=1/2, m=1}|0>(Analitical)")
  END IF

  CALL Action(Op_psi_real, x_band_half_1, b_1)
  Op_psi_real_ana(:) = [x_band_half_1%Band_val(2,3), ZERO, x_band_half_1%Band_val(2,1)]
  CALL Equal_tensor(error_action, Op_psi_real, Op_psi_real_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="x_{band, w=1/2, m=1}|1>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_real, out_unit, Size(Op_psi_real), info="x_{band, w=1/2, m=1}|1>")
    CALL Write_Vec(Op_psi_real_ana, out_unit, Size(Op_psi_real_ana), info="x_{band, w=1/2, m=1}|1>(Analitical)")
  END IF

  CALL Action(Op_psi_real, x_band_half_1, b_2)
  Op_psi_real_ana(:) = [ZERO, x_band_half_1%Band_val(3,3), ZERO]
  CALL Equal_tensor(error_action, Op_psi_real, Op_psi_real_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="x_{band, w=1/2, m=1}|2>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_real, out_unit, Size(Op_psi_real), info="x_{band, w=1/2, m=1}|2>")
    CALL Write_Vec(Op_psi_real_ana, out_unit, Size(Op_psi_real_ana), info="x_{band, w=1/2, m=1}|2>(Analitical)")
  END IF

  CALL Action(Op_psi_real, x_dense_half_1, b_0)
  Op_psi_real_ana(:) = [ZERO, x_dense_half_1%Dense_val(2,1), ZERO]
  CALL Equal_tensor(error_action, Op_psi_real, Op_psi_real_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="x_{dense, w=1/2, m=1}|0>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_real, out_unit, Size(Op_psi_real), info="x_{dense, w=1/2, m=1}|0>")
    CALL Write_Vec(Op_psi_real_ana, out_unit, Size(Op_psi_real_ana), info="x_{dense, w=1/2, m=1}|0>(Analitical)")
  END IF

  CALL Action(Op_psi_real, x_dense_half_1, b_1)
  Op_psi_real_ana(:) = [x_dense_half_1%Dense_val(1,2), ZERO, x_dense_half_1%Dense_val(3,2)]
  CALL Equal_tensor(error_action, Op_psi_real, Op_psi_real_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="x_{dense, w=1/2, m=1}|1>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_real, out_unit, Size(Op_psi_real), info="x_{dense, w=1/2, m=1}|1>")
    CALL Write_Vec(Op_psi_real_ana, out_unit, Size(Op_psi_real_ana), info="x_{dense, w=1/2, m=1}|1>(Analitical)")
  END IF

  CALL Action(Op_psi_real, x_dense_half_1, b_2)
  Op_psi_real_ana(:) = [ZERO, x_dense_half_1%Dense_val(2,3), ZERO]
  CALL Equal_tensor(error_action, Op_psi_real, Op_psi_real_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="x_{dense, w=1/2, m=1}|2>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_real, out_unit, Size(Op_psi_real), info="x_{dense, w=1/2, m=1}|2>")
    CALL Write_Vec(Op_psi_real_ana, out_unit, Size(Op_psi_real_ana), info="x_{dense, w=1/2, m=1}|2>(Analitical)")
  END IF

  CALL Action(Op_psi_real, x_band_half_1, Psi_1D_R1_real)
  Op_psi_real_ana(:) = [Coeff_1_real*x_band_half_1%Band_val(2,3), Coeff_0_real*x_band_half_1%Band_val(1,1) + C&
                  &oeff_2_real*x_band_half_1%Band_val(3,3), Coeff_1_real*x_band_half_1%Band_val(2,1)] / Norm
  CALL Equal_tensor(error_action, Op_psi_real, Op_psi_real_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="x_{band, w=1/2, m=1}|Psi_1D_R1_real>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_real*Norm, out_unit, Size(Op_psi_real), info="x_{band, w=1/2, m=1}|Psi_1D_R1_real>*Norm")
    CALL Write_Vec(Op_psi_real_ana*Norm, out_unit, Size(Op_psi_real_ana), info="x_{band, w=1/2, m=1}|Psi_1D_R1_real>*Norm(Analiti&
    &cal)")
  END IF

  CALL Action(Op_psi_real, x_dense_half_1, Psi_1D_R1_real)
  Op_psi_real_ana(:) = [Coeff_1_real*x_dense_half_1%Dense_val(1,2), Coeff_0_real*x_dense_half_1%Dense_val(2,1)&
                 & + Coeff_2_real*x_dense_half_1%Dense_val(2,3), Coeff_1_real*x_dense_half_1%Dense_val(3,2)] / Norm
  CALL Equal_tensor(error_action, Op_psi_real, Op_psi_real_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="x_{dense, w=1/2, m=1}|Psi_1D_R1_real>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_real*Norm, out_unit, Size(Op_psi_real), info="x_{dense, w=1/2, m=1}|Psi_1D_R1_real>*Norm")
    CALL Write_Vec(Op_psi_real_ana*Norm, out_unit, Size(Op_psi_real_ana), info="x_{dense, w=1/2, m=1}|Psi_1D_R1_real>*Norm(Analit&
    &ical)")
  END IF

    !-------------------------------x_{w=3, m=1}-------------------------------
  CALL Action(Op_psi_real, x_band_3_1, b_0)
  Op_psi_real_ana(:) = [ZERO, x_band_3_1%Band_val(1,1), ZERO]
  CALL Equal_tensor(error_action, Op_psi_real, Op_psi_real_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="x_{band, w=3, m=1}|0>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_real, out_unit, Size(Op_psi_real), info="x_{band, w=3, m=1}|0>")
    CALL Write_Vec(Op_psi_real_ana, out_unit, Size(Op_psi_real_ana), info="x_{band, w=3, m=1}|0>(Analitical)")
  END IF

  CALL Action(Op_psi_real, x_band_3_1, b_1)
  Op_psi_real_ana(:) = [x_band_3_1%Band_val(2,3), ZERO, x_band_3_1%Band_val(2,1)]
  CALL Equal_tensor(error_action, Op_psi_real, Op_psi_real_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="x_{band, w=3, m=1}|1>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_real, out_unit, Size(Op_psi_real), info="x_{band, w=3, m=1}|1>")
    CALL Write_Vec(Op_psi_real_ana, out_unit, Size(Op_psi_real_ana), info="x_{band, w=3, m=1}|1>(Analitical)")
  END IF

  CALL Action(Op_psi_real, x_band_3_1, b_2)
  Op_psi_real_ana(:) = [ZERO, x_band_3_1%Band_val(3,3), ZERO]
  CALL Equal_tensor(error_action, Op_psi_real, Op_psi_real_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="x_{band, w=3, m=1}|2>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_real, out_unit, Size(Op_psi_real), info="x_{band, w=3, m=1}|2>")
    CALL Write_Vec(Op_psi_real_ana, out_unit, Size(Op_psi_real_ana), info="x_{band, w=3, m=1}|2>(Analitical)")
  END IF

  CALL Action(Op_psi_real, x_dense_3_1, b_0)
  Op_psi_real_ana(:) = [ZERO, x_dense_3_1%Dense_val(2,1), ZERO]
  CALL Equal_tensor(error_action, Op_psi_real, Op_psi_real_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="x_{dense, w=3, m=1}|0>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_real, out_unit, Size(Op_psi_real), info="x_{dense, w=3, m=1}|0>")
    CALL Write_Vec(Op_psi_real_ana, out_unit, Size(Op_psi_real_ana), info="x_{dense, w=3, m=1}|0>(Analitical)")
  END IF

  CALL Action(Op_psi_real, x_dense_3_1, b_1)
  Op_psi_real_ana(:) = [x_dense_3_1%Dense_val(1,2), ZERO, x_dense_3_1%Dense_val(3,2)]
  CALL Equal_tensor(error_action, Op_psi_real, Op_psi_real_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="x_{dense, w=3, m=1}|1>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_real, out_unit, Size(Op_psi_real), info="x_{dense, w=3, m=1}|1>")
    CALL Write_Vec(Op_psi_real_ana, out_unit, Size(Op_psi_real_ana), info="x_{dense, w=3, m=1}|1>(Analitical)")
  END IF

  CALL Action(Op_psi_real, x_dense_3_1, b_2)
  Op_psi_real_ana(:) = [ZERO, x_dense_3_1%Dense_val(2,3), ZERO]
  CALL Equal_tensor(error_action, Op_psi_real, Op_psi_real_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="x_{dense, w=3, m=1}|2>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_real, out_unit, Size(Op_psi_real), info="x_{dense, w=3, m=1}|2>")
    CALL Write_Vec(Op_psi_real_ana, out_unit, Size(Op_psi_real_ana), info="x_{dense, w=3, m=1}|2>(Analitical)")
  END IF

  CALL Action(Op_psi_real, x_band_3_1, Psi_1D_R1_real)
  Op_psi_real_ana(:) = [Coeff_1_real*x_band_3_1%Band_val(2,3), Coeff_0_real*x_band_3_1%Band_val(1,1) + Coeff_2_real&
                  &*x_band_3_1%Band_val(3,3), Coeff_1_real*x_band_3_1%Band_val(2,1)] / Norm
  CALL Equal_tensor(error_action, Op_psi_real, Op_psi_real_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="x_{band, w=3, m=1}|Psi_1D_R1_real>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_real*Norm, out_unit, Size(Op_psi_real), info="x_{band, w=3, m=1}|Psi_1D_R1_real>*Norm")
    CALL Write_Vec(Op_psi_real_ana*Norm, out_unit, Size(Op_psi_real_ana), info="x_{band, w=3, m=1}|Psi_1D_R1_real>*Norm(Analitica&
    &l)")
  END IF

  CALL Action(Op_psi_real, x_dense_3_1, Psi_1D_R1_real)
  Op_psi_real_ana(:) = [Coeff_1_real*x_dense_3_1%Dense_val(1,2), Coeff_0_real*x_dense_3_1%Dense_val(2,1) + Coe&
                  &ff_2_real*x_dense_3_1%Dense_val(2,3), Coeff_1_real*x_dense_3_1%Dense_val(3,2)] / Norm
  CALL Equal_tensor(error_action, Op_psi_real, Op_psi_real_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="x_{dense, w=3, m=1}|Psi_1D_R1_real>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_real*Norm, out_unit, Size(Op_psi_real), info="x_{dense, w=3, m=1}|Psi_1D_R1_real>*Norm")
    CALL Write_Vec(Op_psi_real_ana*Norm, out_unit, Size(Op_psi_real_ana), info="x_{dense, w=3, m=1}|Psi_1D_R1_real>*Norm(Analitic&
    &al)")
  END IF

    !-------------------------------x_{w=3, m=4}-------------------------------
  CALL Action(Op_psi_real, x_band_3_4, b_0)
  Op_psi_real_ana(:) = [ZERO, x_band_3_4%Band_val(1,1), ZERO]
  CALL Equal_tensor(error_action, Op_psi_real, Op_psi_real_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="x_{band, w=3, m=4}|0>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_real, out_unit, Size(Op_psi_real), info="x_{band, w=3, m=4}|0>")
    CALL Write_Vec(Op_psi_real_ana, out_unit, Size(Op_psi_real_ana), info="x_{band, w=3, m=4}|0>(Analitical)")
  END IF

  CALL Action(Op_psi_real, x_band_3_4, b_1)
  Op_psi_real_ana(:) = [x_band_3_4%Band_val(2,3), ZERO, x_band_3_4%Band_val(2,1)]
  CALL Equal_tensor(error_action, Op_psi_real, Op_psi_real_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="x_{band, w=3, m=4}|1>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_real, out_unit, Size(Op_psi_real), info="x_{band, w=3, m=4}|1>")
    CALL Write_Vec(Op_psi_real_ana, out_unit, Size(Op_psi_real_ana), info="x_{band, w=3, m=4}|1>(Analitical)")
  END IF

  CALL Action(Op_psi_real, x_band_3_4, b_2)
  Op_psi_real_ana(:) = [ZERO, x_band_3_4%Band_val(3,3), ZERO]
  CALL Equal_tensor(error_action, Op_psi_real, Op_psi_real_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="x_{band, w=3, m=4}|2>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_real, out_unit, Size(Op_psi_real), info="x_{band, w=3, m=4}|2>")
    CALL Write_Vec(Op_psi_real_ana, out_unit, Size(Op_psi_real_ana), info="x_{band, w=3, m=4}|2>(Analitical)")
  END IF

  CALL Action(Op_psi_real, x_dense_3_4, b_0)
  Op_psi_real_ana(:) = [ZERO, x_dense_3_4%Dense_val(2,1), ZERO]
  CALL Equal_tensor(error_action, Op_psi_real, Op_psi_real_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="x_{dense, w=3, m=4}|0>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_real, out_unit, Size(Op_psi_real), info="x_{dense, w=3, m=4}|0>")
    CALL Write_Vec(Op_psi_real_ana, out_unit, Size(Op_psi_real_ana), info="x_{dense, w=3, m=4}|0>(Analitical)")
  END IF

  CALL Action(Op_psi_real, x_dense_3_4, b_1)
  Op_psi_real_ana(:) = [x_dense_3_4%Dense_val(1,2), ZERO, x_dense_3_4%Dense_val(3,2)]
  CALL Equal_tensor(error_action, Op_psi_real, Op_psi_real_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="x_{dense, w=3, m=4}|1>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_real, out_unit, Size(Op_psi_real), info="x_{dense, w=3, m=4}|1>")
    CALL Write_Vec(Op_psi_real_ana, out_unit, Size(Op_psi_real_ana), info="x_{dense, w=3, m=4}|1>(Analitical)")
  END IF

  CALL Action(Op_psi_real, x_dense_3_4, b_2)
  Op_psi_real_ana(:) = [ZERO, x_dense_3_4%Dense_val(2,3), ZERO]
  CALL Equal_tensor(error_action, Op_psi_real, Op_psi_real_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="x_{dense, w=3, m=4}|2>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_real, out_unit, Size(Op_psi_real), info="x_{dense, w=3, m=4}|2>")
    CALL Write_Vec(Op_psi_real_ana, out_unit, Size(Op_psi_real_ana), info="x_{dense, w=3, m=4}|2>(Analitical)")
  END IF

  CALL Action(Op_psi_real, x_band_3_4, Psi_1D_R1_real)
  Op_psi_real_ana(:) = [Coeff_1_real*x_band_3_4%Band_val(2,3), Coeff_0_real*x_band_3_4%Band_val(1,1) + Coeff_2_real&
                  &*x_band_3_4%Band_val(3,3), Coeff_1_real*x_band_3_4%Band_val(2,1)] / Norm
  CALL Equal_tensor(error_action, Op_psi_real, Op_psi_real_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="x_{band, w=3, m=4}|Psi_1D_R1_real>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_real*Norm, out_unit, Size(Op_psi_real), info="x_{band, w=3, m=4}|Psi_1D_R1_real>*Norm")
    CALL Write_Vec(Op_psi_real_ana*Norm, out_unit, Size(Op_psi_real_ana), info="x_{band, w=3, m=4}|Psi_1D_R1_real>*Norm(Analitica&
    &l)")
  END IF

  CALL Action(Op_psi_real, x_dense_3_4, Psi_1D_R1_real)
  Op_psi_real_ana(:) = [Coeff_1_real*x_dense_3_4%Dense_val(1,2), Coeff_0_real*x_dense_3_4%Dense_val(2,1) + Coe&
                  &ff_2_real*x_dense_3_4%Dense_val(2,3), Coeff_1_real*x_dense_3_4%Dense_val(3,2)] / Norm
  CALL Equal_tensor(error_action, Op_psi_real, Op_psi_real_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="x_{dense, w=3, m=4}|Psi_1D_R1_real>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_real*Norm, out_unit, Size(Op_psi_real), info="x_{dense, w=3, m=4}|Psi_1D_R1_real>*Norm")
    CALL Write_Vec(Op_psi_real_ana*Norm, out_unit, Size(Op_psi_real_ana), info="x_{dense, w=3, m=4}|Psi_1D_R1_real>*Norm(Analitic&
    &al)")
  END IF


  !-------------------------N matricies initialization-------------------------
  CALL Construct_Operator_1D(N_diag,  "Nb_photons", Mode=Mode, Debug=Debug)
  CALL Construct_Operator_1D(N_dense, "Nb_photons", Dense=.TRUE., Mode=Mode, Debug=Debug)


  !----------------------------Testing the N actions---------------------------
  CALL Action(Op_psi_real, N_diag, b_0)
  Op_psi_real_ana(:) = ZERO
  CALL Equal_tensor(error_action, Op_psi_real, Op_psi_real_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="N_{diag}|0>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_real, out_unit, Size(Op_psi_real), info="N_{diag}|0>")
    CALL Write_Vec(Op_psi_real_ana, out_unit, Size(Op_psi_real_ana), info="N_{diag}|0>(Analitical)")
  END IF

  CALL Action(Op_psi_real, N_diag, b_1)
  Op_psi_real_ana(:) = [ZERO, N_diag%Diag_val(2), ZERO]
  CALL Equal_tensor(error_action, Op_psi_real, Op_psi_real_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="N_{diag}|1>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_real, out_unit, Size(Op_psi_real), info="N_{diag}|1>")
    CALL Write_Vec(Op_psi_real_ana, out_unit, Size(Op_psi_real_ana), info="N_{diag}|1>(Analitical)")
  END IF

  CALL Action(Op_psi_real, N_diag, b_2)
  Op_psi_real_ana(:) = [ZERO, ZERO, N_diag%Diag_val(3)]
  CALL Equal_tensor(error_action, Op_psi_real, Op_psi_real_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="N_{diag}|2>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_real, out_unit, Size(Op_psi_real), info="N_{diag}|2>")
    CALL Write_Vec(Op_psi_real_ana, out_unit, Size(Op_psi_real_ana), info="N_{diag}|2>(Analitical)")
  END IF

  CALL Action(Op_psi_real, N_dense, b_0)
  Op_psi_real_ana(:) = ZERO
  CALL Equal_tensor(error_action, Op_psi_real, Op_psi_real_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="N_{dense}|0>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_real, out_unit, Size(Op_psi_real), info="N_{dense}|0>")
    CALL Write_Vec(Op_psi_real_ana, out_unit, Size(Op_psi_real_ana), info="N_{dense}|0>(Analitical)")
  END IF

  CALL Action(Op_psi_real, N_dense, b_1)
  Op_psi_real_ana(:) = [ZERO, N_dense%Dense_val(2,2), ZERO]
  CALL Equal_tensor(error_action, Op_psi_real, Op_psi_real_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="N_{dense}|1>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_real, out_unit, Size(Op_psi_real), info="N_{dense}|1>")
    CALL Write_Vec(Op_psi_real_ana, out_unit, Size(Op_psi_real_ana), info="N_{dense}|1>(Analitical)")
  END IF

  CALL Action(Op_psi_real, N_dense, b_2)
  Op_psi_real_ana(:) = [ZERO, ZERO, N_dense%Dense_val(3,3)]
  CALL Equal_tensor(error_action, Op_psi_real, Op_psi_real_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="N_{dense}|2>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_real, out_unit, Size(Op_psi_real), info="N_{dense}|2>")
    CALL Write_Vec(Op_psi_real_ana, out_unit, Size(Op_psi_real_ana), info="N_{dense}|2>(Analitical)")
  END IF

  CALL Action(Op_psi_real, N_diag, Psi_1D_R1_real)
  Op_psi_real_ana(:) = [ZERO, Coeff_1_real*N_diag%Diag_val(2), Coeff_2_real*N_diag%Diag_val(3)] / Norm
  CALL Equal_tensor(error_action, Op_psi_real, Op_psi_real_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="N_{diag}|Psi_1D_R1_real>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_real*Norm, out_unit, Size(Op_psi_real), info="N_{diag}|Psi_1D_R1_real>*Norm")
    CALL Write_Vec(Op_psi_real_ana*Norm, out_unit, Size(Op_psi_real_ana), info="N_{diag}|Psi_1D_R1_real>*Norm(Analitical)")
  END IF

  CALL Action(Op_psi_real, N_dense, Psi_1D_R1_real)
  Op_psi_real_ana(:) = [ZERO, Coeff_1_real*N_dense%Dense_val(2,2), Coeff_2_real*N_dense%Dense_val(3,3)] / Norm
  CALL Equal_tensor(error_action, Op_psi_real, Op_psi_real_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="N_{dense}|Psi_1D_R1_real>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_real*Norm, out_unit, Size(Op_psi_real), info="N_{dense}|Psi_1D_R1_real>*Norm")
    CALL Write_Vec(Op_psi_real_ana*Norm, out_unit, Size(Op_psi_real_ana), info="N_{dense}|Psi_1D_R1_real>*Norm(Analitical)")
  END IF

    !-------------------------Wavefunction initialization (complex)------------------------
  Psi_1D_R1_complex(:) = [Coeff_0_complex, Coeff_1_complex, Coeff_2_complex]
  CALL Norm_of(Norm, Psi_1D_R1_complex)
  CALL Normalize(Psi_1D_R1_complex)
  IF (Debug) THEN
    WRITE(out_unit,*)
    WRITE(out_unit,*) "----------------The cavity wavefunction has been initialized as :---------------"
    CALL Write_Vec(Psi_1D_R1_complex, out_unit, Size(Psi_1D_R1_complex), info="Psi_1D_R1_complex")
    WRITE(out_unit,*) "-----------------------------End cavity wavefunction----------------------------"
  END IF

  
  CALL Finalize_Test(test_action)
  
  
  CONTAINS


  SUBROUTINE MolecCav_Construct_Elem_op(Operator, Coeff, Operator_type, Dense, Debug_opt)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
    USE QDUtil_m
    USE Elem_op_m
    IMPLICIT NONE

    TYPE(Elem_op_t),   intent(inout) :: Operator                                                                                 ! the object of type Elem_op_t to be constructed here
    real(kind=Rkind),  intent(in)    :: Coeff
    character(len=*),  intent(in)    :: Operator_type                                                                            ! ex : "Hamiltonian", "Position", etc. (len=:) Expects to be allocatable, while (len=*) is dedicated to a procedure argument.
    logical, optional, intent(in)    :: Dense                                                                                    ! if .TRUE. then the matrix storage will not be optimized and it will be stored as a Dense matrix
    logical, optional, intent(in)    :: Debug_opt

    logical                          :: Debug_local = .FALSE.

    !-----------------------------Debugging options----------------------------
    IF (PRESENT(Debug_opt)) Debug_local = Debug_opt
    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "-------------------Arguments of MolecCav_Construct_Elem_op------------------"
      WRITE(out_unit,*) "The <<Operator>> argument :"
      CALL Write(Operator)
      WRITE(out_unit,*) "The <<Operator_type>> argument : "//Operator_type
      WRITE(out_unit,*) "The <<Coeff>> argument :"//TO_string(Coeff)
      IF (PRESENT(Dense)) WRITE(out_unit,*) "Dense :", Dense
      WRITE(out_unit,*) "-----------------End Arguments of MolecCav_Construct_Elem_op----------------"
    END IF
    
    !--------------First steps of the construction of the Operator-------------
    ALLOCATE(character(len=LEN_TRIM(Operator_type)) :: Operator%Operator_type)                                                   ! /!\ strings cannot be allocated the exact same way as tables ! /!\
    Operator%Operator_type = TO_lowercase(TRIM(Operator_type))                                                                   ! allocation on assignement (not anymore). Operator_type has the right lengths (no spaces added) thanks to len=* at declaration and it will fit the Op%op_type thanks to len=:, allocatable at declaration of the derived type. 

    IF (PRESENT(Dense)) THEN
      Operator%Dense = Dense
    END IF

    !--------------------Construction of the matrix Operator-------------------
    SELECT CASE (Operator%Operator_type)                                                                                         ! we actually do not care about the operator type here but only about the shape of the matrix representation, that is why only H and x are possible to use to have one diag and one band
      CASE ("hamiltonian")
        CALL MolecCav_Construct_diag_elem_op(Operator=Operator, Coeff=Coeff)                                                     ! /!\ contrary to the true procedures, diag and band means here the shape of the actual full analytical matrix and not the representation used to store the operator's matrices as in the module Quantum_HO1D_m /!\ 
    
      CASE ("position")
        CALL MolecCav_Construct_band_elem_op(Operator=Operator, Coeff=Coeff)                                                     ! /!\ contrary to the true procedures, diag and band means here the shape of the actual full analytical matrix and not the representation used to store the operator's matrices as in the module Quantum_HO1D_m /!\
      
      CASE DEFAULT
        WRITE(out_unit,*) "No Operator type recognized, please verify the input of Construct_Elem_op subroutine"
        STOP "### No Operator type recognized, please verify the input of Construct_Elem_op subroutine"
    END SELECT

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "-------------Operator constructed by MolecCav_Construct_Elem_op-------------"
      CALL Write(Operator)
      WRITE(out_unit,*) "-----------End operator constructed by MolecCav_Construct_Elem_op-----------"
    END IF

  END SUBROUTINE MolecCav_Construct_Elem_op


  SUBROUTINE MolecCav_Construct_diag_elem_op(Operator, Coeff)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    USE Elem_op_m
    IMPLICIT NONE
    
    TYPE(Elem_op_t),  intent(inout) :: Operator
    real(kind=Rkind), intent(in)    :: Coeff
 
    integer                         :: i                                                                                         ! loop increments /!\ Fortran counts from 1 to Nb !!! /!\

    WRITE(out_unit,*)
    WRITE(out_unit,*) '*******************************************************'
    WRITE(out_unit,*) '********** CONSTRUCTING THE DIAGONAL MATRIX ***********'

    IF (.NOT. Operator%Dense) THEN
      !---------------------Initialization to default values-------------------
      ALLOCATE(Operator%Diag_val(3))
      !------------------------Construction of the matrix----------------------
      DO i = 1, 3                                                                                                                ! /!\ Fortran counts from 1 to Nb !!! /!\
        Operator%Diag_val(i) = Coeff*(i - ONE + HALF)                                                                            ! "-1" because the first Fortran vector is the fundamental eigenvector of the HO i.e. the 0^{th} ket 
      END DO
      
    ELSE
      !---------------------Initialization to default values-------------------
      ALLOCATE(Operator%Dense_val(3, 3))
      Operator%Dense_val = ZERO
      !------------------------Construction of the matrix----------------------
      DO i = 1, 3                                                                                                                ! /!\ Fortran counts from 1 to Nb !!! /!\
        Operator%Dense_val(i,i) = Coeff*(i - ONE + HALF)                                                                         ! "-1" because the first Fortran vector is the fundamental eigenvector of the HO i.e. the 0^{th} ket 
      END DO
    END IF
    
    WRITE(out_unit,*) '************* DIAGONAL MATRIX CONSTRUCTED *************'
    WRITE(out_unit,*) '*******************************************************'

  END SUBROUTINE MolecCav_Construct_diag_elem_op


  SUBROUTINE MolecCav_Construct_band_elem_op(Operator, Coeff)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    USE Elem_op_m
    IMPLICIT NONE
    
    TYPE(Elem_op_t),  intent(inout) :: Operator
    real(kind=Rkind), intent(in)    :: Coeff

    integer                         :: i                                                                                         ! loop increments /!\ Fortran counts from 1 to Nb !!! /!\

    WRITE(out_unit,*) ''
    WRITE(out_unit,*) '*******************************************************'
    WRITE(out_unit,*) '************ CONSTRUCTING THE BAND MATRIX *************'

    IF ((.NOT. Operator%Dense) .AND. 3 > 1) THEN
      !----------Initialization of the characteristics of the operator---------
      Operator%Upper_bandwidth   = 1
      Operator%Lower_bandwidth   = 1
      !---------------------Initialization to default values-------------------
      ALLOCATE(Operator%Band_val(3,3))                                                                                           ! Nb lines (number of diagonal elements) and 3 columns because 3 bands to consider : the diagonal, and the two bands above and below it
      Operator%Band_val = ZERO
      !------------------------Construction of the matrix----------------------
      DO i = 1, 3 - 1                                                                                                            ! /!\ Fortran counts from 1 to Nb !!! /!\ Nb-1 not to have Band_val(i+1) out of range
        Operator%Band_val(i,1)   = SQRT(REAL(i,kind=Rkind))
        Operator%Band_val(i+1,3) = SQRT(REAL(i,kind=Rkind))
      END DO
      Operator%Band_val = Operator%Band_val / Coeff
        
    ELSE IF (.NOT. Operator%Dense) THEN
            !---------------------Initialization to default values-------------------
      ALLOCATE(Operator%Diag_val(3))
      !------------------------Construction of the matrix----------------------
      DO i = 1, 3                                                                                                                ! /!\ Fortran counts from 1 to Nb !!! /!\
        Operator%Diag_val(i) = ZERO                                                                                              ! the position operator matrix has first value (i.e. only value in the Nb = 0 case) 0 
      END DO

    ELSE 
      !---------------------Initialization to default values-------------------
      ALLOCATE(Operator%Dense_val(3, 3))
      Operator%Dense_val = ZERO
      !------------------------Construction of the matrix----------------------
      DO i = 1, 3 - 1                                                                                                            ! /!\ Fortran counts from 1 to Nb !!! /!\
        Operator%Dense_val(i,i+1) = SQRT(REAL(i,kind=Rkind))
        Operator%Dense_val(i+1,i) = SQRT(REAL(i,kind=Rkind))
      END DO
      Operator%Dense_val = Operator%Dense_val / Coeff
    END IF
    
    WRITE(out_unit,*) '************** BAND MATRIX CONSTRUCTED ****************'
    WRITE(out_unit,*) '*******************************************************'

  END SUBROUTINE MolecCav_Construct_band_elem_op


  SUBROUTINE MolecCav_Construct_op_psi_ana_real(Op_psi_ana, Operator, Op_coeff, Psi, Debug_opt)                                  ! /!\ Can be used only in this module /!\ (needs the op to have been constructed using the above procedures) The Op_coeff is the argument of the Construct procedures of this module.
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
    USE QDUtil_m
    USE Elem_op_m
    IMPLICIT NONE

    real(kind=Rkind),  intent(inout) :: Op_psi_ana(:)
    TYPE(Elem_op_t),   intent(in)    :: Operator                                                                                 ! the object of type Elem_op_t to be constructed here
    real(kind=Rkind),  intent(in)    :: Op_coeff
    real(kind=Rkind),  intent(in)    :: Psi(:)
    logical, optional, intent(in)    :: Debug_opt

    logical                          :: Debug_local = .FALSE.

    !-----------------------------Debugging options----------------------------
    IF (PRESENT(Debug_opt)) Debug_local = Debug_opt
    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "-------------------Arguments of MolecCav_Construct_Elem_op------------------"
      WRITE(out_unit,*) "The <<Operator>> argument :"
      CALL Write(Operator)
      WRITE(out_unit,*) "The <<Op_coeff>> argument : "//TO_string(Op_coeff)
      WRITE(out_unit,*) "The <<Psi>> argument :"
      CALL Write_Vec(Psi, out_unit, 1, info="Psi")
      WRITE(out_unit,*) "-----------------End Arguments of MolecCav_Construct_Elem_op----------------"
    END IF
    
    !--------------------Construction of the Op_psi_ana matrix-------------------
    SELECT CASE (Operator%Operator_type)                                                                                         ! We actually do not care about the operator type here but only about the shape of the matrix representation, that is why only H and x are possible to use to have one diag and one band. We do not use the allocated(Operator%<shape>_val) criterium because it would not allows to know the shape of the matrices constructed as dense ones
      CASE ("hamiltonian")
        CALL MolecCav_Construct_diag_op_psi_ana_real(Op_psi_ana, Op_coeff, Psi, Debug_local)
    
      CASE ("position")
        CALL MolecCav_Construct_band_op_psi_ana_real(Op_psi_ana, Op_coeff, Psi, Debug_local)
      
      CASE DEFAULT
        WRITE(out_unit,*) "No Operator type recognized, please verify the input of MolecCav_Construct_op_psi_ana_real subroutine"
        STOP "### No Operator type recognized, please verify the input of MolecCav_Construct_op_psi_ana_real subroutine"
    END SELECT

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "-------------Operator constructed by MolecCav_Construct_Elem_op-------------"
      CALL Write(Operator)
      WRITE(out_unit,*) "-----------End operator constructed by MolecCav_Construct_Elem_op-----------"
    END IF

  END SUBROUTINE MolecCav_Construct_op_psi_ana_real


  SUBROUTINE MolecCav_Construct_diag_op_psi_ana_real(Op_psi_ana, Op_coeff, Psi, Debug_opt)                                  ! /!\ Can be used only in this module /!\ (needs the op to have been constructed using the above procedures) The Op_coeff is the argument of the Construct procedures of this module.
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
    USE QDUtil_m
    USE Elem_op_m
    IMPLICIT NONE

    real(kind=Rkind),  intent(inout) :: Op_psi_ana(:)
    real(kind=Rkind),  intent(in)    :: Op_coeff
    real(kind=Rkind),  intent(in)    :: Psi(:)
    logical, optional, intent(in)    :: Debug_opt

    logical                          :: Debug_local = .FALSE.

    !-----------------------------Debugging options----------------------------
    IF (PRESENT(Debug_opt)) Debug_local = Debug_opt
    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "-------------------Arguments of MolecCav_Construct_Elem_op------------------"
      WRITE(out_unit,*) "The <<Op_coeff>> argument : "//TO_string(Op_coeff)
      WRITE(out_unit,*) "The <<Psi>> argument :"
      CALL Write_Vec(Psi, out_unit, 1, info="Psi")
      WRITE(out_unit,*) "-----------------End Arguments of MolecCav_Construct_Elem_op----------------"
    END IF
    
    !--------------------Construction of the Op_psi_ana matrix-------------------
    Op_psi_ana = Op_coeff * [Psi(1), 3*Psi(2), 5*Psi(3)] / TWO 

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "-------------Op_psi_ana constructed by MolecCav_Construct_Elem_op-------------"
      CALL Write_Vec(Op_psi_ana, out_unit, 1, info="Op_psi")
      WRITE(out_unit,*) "-----------End Op_psi_ana constructed by MolecCav_Construct_Elem_op-----------"
    END IF

  END SUBROUTINE MolecCav_Construct_diag_op_psi_ana_real


  SUBROUTINE MolecCav_Construct_band_op_psi_ana_real(Op_psi_ana, Op_coeff, Psi, Debug_opt)                                  ! /!\ Can be used only in this module /!\ (needs the op to have been constructed using the above procedures) The Op_coeff is the argument of the Construct procedures of this module.
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
    USE QDUtil_m
    USE Elem_op_m
    IMPLICIT NONE

    real(kind=Rkind),  intent(inout) :: Op_psi_ana(:)
    real(kind=Rkind),  intent(in)    :: Op_coeff
    real(kind=Rkind),  intent(in)    :: Psi(:)
    logical, optional, intent(in)    :: Debug_opt

    logical                          :: Debug_local = .FALSE.

    !-----------------------------Debugging options----------------------------
    IF (PRESENT(Debug_opt)) Debug_local = Debug_opt
    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "-------------------Arguments of MolecCav_Construct_Elem_op------------------"
      WRITE(out_unit,*) "The <<Op_coeff>> argument : "//TO_string(Op_coeff)
      WRITE(out_unit,*) "The <<Psi>> argument :"
      CALL Write_Vec(Psi, out_unit, 1, info="Psi")
      WRITE(out_unit,*) "-----------------End Arguments of MolecCav_Construct_Elem_op----------------"
    END IF
    
    !--------------------Construction of the Op_psi_ana matrix-------------------
    Op_psi_ana = [Psi(2), Psi(1) + SQRT(TWO)*Psi(3), SQRT(TWO)*Psi(2)] / Op_coeff

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "-------------Operator constructed by MolecCav_Construct_Elem_op-------------"
      CALL Write_Vec(Op_psi_ana, out_unit, 1, info="Op_psi")
      WRITE(out_unit,*) "-----------End operator constructed by MolecCav_Construct_Elem_op-----------"
    END IF

  END SUBROUTINE MolecCav_Construct_band_op_psi_ana_real


  SUBROUTINE MolecCav_Construct_op_psi_ana_complex(Op_psi_ana, Operator, Op_coeff, Psi, Debug_opt)                                  ! /!\ Can be used only in this module /!\ (needs the op to have been constructed using the above procedures) The Op_coeff is the argument of the Construct procedures of this module.
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
    USE QDUtil_m
    USE Elem_op_m
    IMPLICIT NONE

    complex(kind=Rkind), intent(inout) :: Op_psi_ana(:)
    TYPE(Elem_op_t),     intent(in)    :: Operator                                                                                 ! the object of type Elem_op_t to be constructed here
    real(kind=Rkind),    intent(in)    :: Op_coeff
    complex(kind=Rkind), intent(in)    :: Psi(:)
    logical, optional,   intent(in)    :: Debug_opt

    logical                            :: Debug_local = .FALSE.

    !-----------------------------Debugging options----------------------------
    IF (PRESENT(Debug_opt)) Debug_local = Debug_opt
    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "-------------------Arguments of MolecCav_Construct_Elem_op------------------"
      WRITE(out_unit,*) "The <<Operator>> argument :"
      CALL Write(Operator)
      WRITE(out_unit,*) "The <<Op_coeff>> argument : "//TO_string(Op_coeff)
      WRITE(out_unit,*) "The <<Psi>> argument :"
      CALL Write_Vec(Psi, out_unit, 1, info="Psi")
      WRITE(out_unit,*) "-----------------End Arguments of MolecCav_Construct_Elem_op----------------"
    END IF
    
    !--------------------Construction of the Op_psi_ana matrix-------------------
    SELECT CASE (Operator%Operator_type)                                                                                         ! We actually do not care about the operator type here but only about the shape of the matrix representation, that is why only H and x are possible to use to have one diag and one band. We do not use the allocated(Operator%<shape>_val) criterium because it would not allows to know the shape of the matrices constructed as dense ones
      CASE ("hamiltonian")
        CALL MolecCav_Construct_diag_op_psi_ana_complex(Op_psi_ana, Op_coeff, Psi, Debug_local)
    
      CASE ("position")
        CALL MolecCav_Construct_band_op_psi_ana_complex(Op_psi_ana, Op_coeff, Psi, Debug_local)
      
      CASE DEFAULT
        WRITE(out_unit,*) "No Operator type recognized, please verify the input of MolecCav_Construct_op_psi_ana_complex subroutine"
        STOP "### No Operator type recognized, please verify the input of MolecCav_Construct_op_psi_ana_complex subroutine"
    END SELECT

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "-------------Operator constructed by MolecCav_Construct_Elem_op-------------"
      CALL Write(Operator)
      WRITE(out_unit,*) "-----------End operator constructed by MolecCav_Construct_Elem_op-----------"
    END IF

  END SUBROUTINE MolecCav_Construct_op_psi_ana_complex


  SUBROUTINE MolecCav_Construct_diag_op_psi_ana_complex(Op_psi_ana, Op_coeff, Psi, Debug_opt)                                  ! /!\ Can be used only in this module /!\ (needs the op to have been constructed using the above procedures) The Op_coeff is the argument of the Construct procedures of this module.
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
    USE QDUtil_m
    USE Elem_op_m
    IMPLICIT NONE

    complex(kind=Rkind), intent(inout) :: Op_psi_ana(:)
    real(kind=Rkind),    intent(in)    :: Op_coeff
    complex(kind=Rkind), intent(in)    :: Psi(:)
    logical, optional,   intent(in)    :: Debug_opt

    logical                            :: Debug_local = .FALSE.

    !-----------------------------Debugging options----------------------------
    IF (PRESENT(Debug_opt)) Debug_local = Debug_opt
    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "-------------------Arguments of MolecCav_Construct_Elem_op------------------"
      WRITE(out_unit,*) "The <<Op_coeff>> argument : "//TO_string(Op_coeff)
      WRITE(out_unit,*) "The <<Psi>> argument :"
      CALL Write_Vec(Psi, out_unit, 1, info="Psi")
      WRITE(out_unit,*) "-----------------End Arguments of MolecCav_Construct_Elem_op----------------"
    END IF
    
    !--------------------Construction of the Op_psi_ana matrix-------------------
    Op_psi_ana = Op_coeff * [Psi(1), 3*Psi(2), 5*Psi(3)] / TWO 

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "-------------Op_psi_ana constructed by MolecCav_Construct_Elem_op-------------"
      CALL Write_Vec(Op_psi_ana, out_unit, 1, info="Op_psi")
      WRITE(out_unit,*) "-----------End Op_psi_ana constructed by MolecCav_Construct_Elem_op-----------"
    END IF

  END SUBROUTINE MolecCav_Construct_diag_op_psi_ana_complex


  SUBROUTINE MolecCav_Construct_band_op_psi_ana_complex(Op_psi_ana, Op_coeff, Psi, Debug_opt)                                  ! /!\ Can be used only in this module /!\ (needs the op to have been constructed using the above procedures) The Op_coeff is the argument of the Construct procedures of this module.
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
    USE QDUtil_m
    USE Elem_op_m
    IMPLICIT NONE

    complex(kind=Rkind), intent(inout) :: Op_psi_ana(:)
    real(kind=Rkind),    intent(in)    :: Op_coeff
    complex(kind=Rkind), intent(in)    :: Psi(:)
    logical, optional,   intent(in)    :: Debug_opt

    logical                            :: Debug_local = .FALSE.

    !-----------------------------Debugging options----------------------------
    IF (PRESENT(Debug_opt)) Debug_local = Debug_opt
    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "-------------------Arguments of MolecCav_Construct_Elem_op------------------"
      WRITE(out_unit,*) "The <<Op_coeff>> argument : "//TO_string(Op_coeff)
      WRITE(out_unit,*) "The <<Psi>> argument :"
      CALL Write_Vec(Psi, out_unit, 1, info="Psi")
      WRITE(out_unit,*) "-----------------End Arguments of MolecCav_Construct_Elem_op----------------"
    END IF
    
    !--------------------Construction of the Op_psi_ana matrix-------------------
    Op_psi_ana = [Psi(2), Psi(1) + SQRT(TWO)*Psi(3), SQRT(TWO)*Psi(2)] / Op_coeff

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "-------------Operator constructed by MolecCav_Construct_Elem_op-------------"
      CALL Write_Vec(Op_psi_ana, out_unit, 1, info="Op_psi")
      WRITE(out_unit,*) "-----------End operator constructed by MolecCav_Construct_Elem_op-----------"
    END IF

  END SUBROUTINE MolecCav_Construct_band_op_psi_ana_complex


END PROGRAM
