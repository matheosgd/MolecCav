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
PROGRAM test_elem_op
  !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
  USE QDUtil_m
  USE Tests_m
  USE Algebra_m
  USE Elem_op_m
  IMPLICIT NONE

  integer                       :: Verbose = 0
  logical                       :: Debug   = .FALSE.

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

  complex(kind=Rkind)           :: Coeff_0_complex = ONE*EYE + SQRT(TWO)
  complex(kind=Rkind)           :: Coeff_1_complex = HALF
  complex(kind=Rkind)           :: Coeff_2_complex = PI*EYE
  complex(kind=Rkind)           :: Psi_1D_R1_complex(3)                                                                          ! an any vector representing the excitation state/wavefunction = a linear combination of the basis functions from the canonical basis set over \mathbb{R} BUT with complexes expansion coefficients to make complexe WF. /!\ Not normalized yet !
  complex(kind=Rkind)           :: Op_psi_complex(3)                                                                             ! the resulting vector from the action of a 1D operator upon Psi_1D_R1_complex
  complex(kind=Rkind)           :: Op_psi_complex_ana(3)                                                                         ! the analytical action = hard-coded reference for comparison

  TYPE(test_t)                  :: test_action
  logical                       :: error_action = .FALSE.


  !-----------------------------Test initialization----------------------------
  CALL Initialize_Test(test_action, test_name="OUT/test_file_elem_op")
  

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


  !-------------------------Operators matricies initialization-------------------------
  CALL MolecCav_Construct_Elem_op(EO1_diag_1,      Coeff=ONE,      Operator_type="Hamiltonian", Dense=.FALSE., Debug_opt=Debug)
  CALL MolecCav_Construct_Elem_op(EO1_dense_1,     Coeff=ONE,      Operator_type="Hamiltonian", Dense=.TRUE.,  Debug_opt=Debug)
  CALL MolecCav_Construct_Elem_op(EO2_diag_tenth,  Coeff=ONETENTH, Operator_type="Hamiltonian", Dense=.FALSE., Debug_opt=Debug)
  CALL MolecCav_Construct_Elem_op(EO2_dense_tenth, Coeff=ONETENTH, Operator_type="Hamiltonian", Dense=.TRUE.,  Debug_opt=Debug)
  CALL MolecCav_Construct_Elem_op(EO3_band_1,      Coeff=ONE,      Operator_type="Position",    Dense=.FALSE., Debug_opt=Debug)
  CALL MolecCav_Construct_Elem_op(EO3_dense_1,     Coeff=ONE,      Operator_type="Position",    Dense=.TRUE.,  Debug_opt=Debug)
  CALL MolecCav_Construct_Elem_op(EO4_band_pi,     Coeff=PI,       Operator_type="Position",    Dense=.FALSE., Debug_opt=Debug)
  CALL MolecCav_Construct_Elem_op(EO4_dense_pi,    Coeff=PI,       Operator_type="Position",    Dense=.TRUE.,  Debug_opt=Debug)
  
  
  !----------------------------Testing the actions---------------------------
    !---------------------------------(Diag/Dense) Coeff = 1 (EO1)--------------------------------
      !---------------------------------first basis vector (i*)b_0--------------------------------
  CALL MolecCav_Construct_op_psi_ana_real(Op_psi_real_ana, EO1_diag_1, ONE, b_0, Debug) ! same analytical result as with E01_dense_1
  CALL Action(Op_psi_real, EO1_diag_1, b_0, Verbose=Verbose, Debug=Debug)
  CALL Equal_tensor(error_action, Op_psi_real, Op_psi_real_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="\hat{O}_{diag, 1}|0>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_real, out_unit, Size(Op_psi_real), info="\hat{O}_{diag, 1}|0>")
    CALL Write_Vec(Op_psi_real_ana, out_unit, Size(Op_psi_real_ana), info="\hat{O}_{diag, 1}|0>(Analitical)")
  END IF
  CALL Action(Op_psi_real, EO1_dense_1, b_0, Verbose=Verbose, Debug=Debug)
  CALL Equal_tensor(error_action, Op_psi_real, Op_psi_real_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="\hat{O}_{dense, 1}|0>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_real, out_unit, Size(Op_psi_real), info="\hat{O}_{dense, 1}|0>")
    CALL Write_Vec(Op_psi_real_ana, out_unit, Size(Op_psi_real_ana), info="\hat{O}_{dense, 1}|0>(Analitical)")
  END IF

  Op_psi_complex_ana = Op_psi_real_ana*EYE                                                               ! make as if the b_0 was a vector of the canonical basis set on \mathbb{C} i.e. i*b_0
  CALL Action(Op_psi_complex, EO1_diag_1, EYE*b_0, Verbose=Verbose, Debug=Debug)
  CALL Equal_tensor(error_action, Op_psi_complex, Op_psi_complex_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="\hat{O}_{diag, 1}|0>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_complex, out_unit, Size(Op_psi_complex), info="\hat{O}_{diag, 1}|0>")
    CALL Write_Vec(Op_psi_complex_ana, out_unit, Size(Op_psi_complex_ana), info="\hat{O}_{diag, 1}|0>(Analitical)")
  END IF
  CALL Action(Op_psi_complex, EO1_dense_1, EYE*b_0, Verbose=Verbose, Debug=Debug)
  CALL Equal_tensor(error_action, Op_psi_complex, Op_psi_complex_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="\hat{O}_{dense, 1}|0>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_complex, out_unit, Size(Op_psi_complex), info="\hat{O}_{dense, 1}|0>")
    CALL Write_Vec(Op_psi_complex_ana, out_unit, Size(Op_psi_complex_ana), info="\hat{O}_{dense, 1}|0>(Analitical)")
  END IF
      !---------------------------------second basis vector (i*)b_1--------------------------------
  CALL MolecCav_Construct_op_psi_ana_real(Op_psi_real_ana, EO1_diag_1, ONE, b_1, Debug) ! same analytical result as with E01_dense_1
  CALL Action(Op_psi_real, EO1_diag_1, b_1, Verbose=Verbose, Debug=Debug)
  CALL Equal_tensor(error_action, Op_psi_real, Op_psi_real_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="\hat{O}_{diag, 1}|1>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_real, out_unit, Size(Op_psi_real), info="\hat{O}_{diag, 1}|1>")
    CALL Write_Vec(Op_psi_real_ana, out_unit, Size(Op_psi_real_ana), info="\hat{O}_{diag, 1}|1>(Analitical)")
  END IF
  CALL Action(Op_psi_real, EO1_dense_1, b_1, Verbose=Verbose, Debug=Debug)
  CALL Equal_tensor(error_action, Op_psi_real, Op_psi_real_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="\hat{O}_{dense, 1}|1>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_real, out_unit, Size(Op_psi_real), info="\hat{O}_{dense, 1}|1>")
    CALL Write_Vec(Op_psi_real_ana, out_unit, Size(Op_psi_real_ana), info="\hat{O}_{dense, 1}|1>(Analitical)")
  END IF

  Op_psi_complex_ana = Op_psi_real_ana*EYE                                                               ! make as if the b_1 was a vector of the canonical basis set on \mathbb{C} i.e. i*b_1
  CALL Action(Op_psi_complex, EO1_diag_1, EYE*b_1, Verbose=Verbose, Debug=Debug)
  CALL Equal_tensor(error_action, Op_psi_complex, Op_psi_complex_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="\hat{O}_{diag, 1}|1>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_complex, out_unit, Size(Op_psi_complex), info="\hat{O}_{diag, 1}|1>")
    CALL Write_Vec(Op_psi_complex_ana, out_unit, Size(Op_psi_complex_ana), info="\hat{O}_{diag, 1}|1>(Analitical)")
  END IF
  CALL Action(Op_psi_complex, EO1_dense_1, EYE*b_1, Verbose=Verbose, Debug=Debug)
  CALL Equal_tensor(error_action, Op_psi_complex, Op_psi_complex_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="\hat{O}_{dense, 1}|1>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_complex, out_unit, Size(Op_psi_complex), info="\hat{O}_{dense, 1}|1>")
    CALL Write_Vec(Op_psi_complex_ana, out_unit, Size(Op_psi_complex_ana), info="\hat{O}_{dense, 1}|1>(Analitical)")
  END IF
      !---------------------------------third basis vector (i*)b_2--------------------------------
  CALL MolecCav_Construct_op_psi_ana_real(Op_psi_real_ana, EO1_diag_1, ONE, b_2, Debug) ! same analytical result as with E01_dense_1
  CALL Action(Op_psi_real, EO1_diag_1, b_2, Verbose=Verbose, Debug=Debug)
  CALL Equal_tensor(error_action, Op_psi_real, Op_psi_real_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="\hat{O}_{diag, 1}|2>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_real, out_unit, Size(Op_psi_real), info="\hat{O}_{diag, 1}|2>")
    CALL Write_Vec(Op_psi_real_ana, out_unit, Size(Op_psi_real_ana), info="\hat{O}_{diag, 1}|2>(Analitical)")
  END IF
  CALL Action(Op_psi_real, EO1_dense_1, b_2, Verbose=Verbose, Debug=Debug)
  CALL Equal_tensor(error_action, Op_psi_real, Op_psi_real_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="\hat{O}_{dense, 1}|2>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_real, out_unit, Size(Op_psi_real), info="\hat{O}_{dense, 1}|2>")
    CALL Write_Vec(Op_psi_real_ana, out_unit, Size(Op_psi_real_ana), info="\hat{O}_{dense, 1}|2>(Analitical)")
  END IF

  Op_psi_complex_ana = Op_psi_real_ana*EYE                                                               ! make as if the b_2 was a vector of the canonical basis set on \mathbb{C} i.e. i*b_2
  CALL Action(Op_psi_complex, EO1_diag_1, EYE*b_2, Verbose=Verbose, Debug=Debug)
  CALL Equal_tensor(error_action, Op_psi_complex, Op_psi_complex_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="\hat{O}_{diag, 1}|2>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_complex, out_unit, Size(Op_psi_complex), info="\hat{O}_{diag, 1}|2>")
    CALL Write_Vec(Op_psi_complex_ana, out_unit, Size(Op_psi_complex_ana), info="\hat{O}_{diag, 1}|2>(Analitical)")
  END IF
  CALL Action(Op_psi_complex, EO1_dense_1, EYE*b_2, Verbose=Verbose, Debug=Debug)
  CALL Equal_tensor(error_action, Op_psi_complex, Op_psi_complex_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="\hat{O}_{dense, 1}|2>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_complex, out_unit, Size(Op_psi_complex), info="\hat{O}_{dense, 1}|2>")
    CALL Write_Vec(Op_psi_complex_ana, out_unit, Size(Op_psi_complex_ana), info="\hat{O}_{dense, 1}|2>(Analitical)")
  END IF      
      !---------------------------------WF LC of the basis vectors--------------------------------
  CALL MolecCav_Construct_op_psi_ana_real(Op_psi_real_ana, EO1_diag_1, ONE, Psi_1D_R1_real, Debug) ! same analytical result as with E01_dense_1
  CALL Action(Op_psi_real, EO1_diag_1, Psi_1D_R1_real, Verbose=Verbose, Debug=Debug)
  CALL Equal_tensor(error_action, Op_psi_real, Op_psi_real_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="\hat{O}_{diag, 1}|\Psi_{real}>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_real, out_unit, Size(Op_psi_real), info="\hat{O}_{diag, 1}|\Psi_{real}>")
    CALL Write_Vec(Op_psi_real_ana, out_unit, Size(Op_psi_real_ana), info="\hat{O}_{diag, 1}|\Psi_{real}>(Analitical)")
  END IF
  CALL Action(Op_psi_real, EO1_dense_1, Psi_1D_R1_real, Verbose=Verbose, Debug=Debug)
  CALL Equal_tensor(error_action, Op_psi_real, Op_psi_real_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="\hat{O}_{dense, 1}|\Psi_{real}>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_real, out_unit, Size(Op_psi_real), info="\hat{O}_{dense, 1}|\Psi_{real}>")
    CALL Write_Vec(Op_psi_real_ana, out_unit, Size(Op_psi_real_ana), info="\hat{O}_{dense, 1}|\Psi_{real}>(Analitical)")
  END IF

  CALL MolecCav_Construct_op_psi_ana_complex(Op_psi_complex_ana, EO1_diag_1, ONE, Psi_1D_R1_complex, Debug) ! same analytical result as with E01_dense_1
  CALL Action(Op_psi_complex, EO1_diag_1, Psi_1D_R1_complex, Verbose=Verbose, Debug=Debug)
  CALL Equal_tensor(error_action, Op_psi_complex, Op_psi_complex_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="\hat{O}_{diag, 1}|\Psi_{complex}>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_complex, out_unit, Size(Op_psi_complex), info="\hat{O}_{diag, 1}|\Psi_{complex}>")
    CALL Write_Vec(Op_psi_complex_ana, out_unit, Size(Op_psi_complex_ana), info="\hat{O}_{diag, 1}|\Psi_{complex}>(Analitical)")
  END IF
  CALL Action(Op_psi_complex, EO1_dense_1, Psi_1D_R1_complex, Verbose=Verbose, Debug=Debug)
  CALL Equal_tensor(error_action, Op_psi_complex, Op_psi_complex_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="\hat{O}_{dense, 1}|\Psi_{complex}>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_complex, out_unit, Size(Op_psi_complex), info="\hat{O}_{dense, 1}|\Psi_{complex}>")
    CALL Write_Vec(Op_psi_complex_ana, out_unit, Size(Op_psi_complex_ana), info="\hat{O}_{dense, 1}|\Psi_{complex}>(Analitical)")
  END IF
    !---------------------------------(Diag/Dense) Coeff = 1/10 (EO2)--------------------------------
      !---------------------------------first basis vector (i*)b_0--------------------------------
  CALL MolecCav_Construct_op_psi_ana_real(Op_psi_real_ana, EO2_diag_tenth, ONETENTH, b_0, Debug) ! same analytical result as with E01_dense_1
  CALL Action(Op_psi_real, EO2_diag_tenth, b_0, Verbose=Verbose, Debug=Debug)
  CALL Equal_tensor(error_action, Op_psi_real, Op_psi_real_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="\hat{O}_{diag, 1/10}|0>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_real, out_unit, Size(Op_psi_real), info="\hat{O}_{diag, 1/10}|0>")
    CALL Write_Vec(Op_psi_real_ana, out_unit, Size(Op_psi_real_ana), info="\hat{O}_{diag, 1/10}|0>(Analitical)")
  END IF
  CALL Action(Op_psi_real, EO2_dense_tenth, b_0, Verbose=Verbose, Debug=Debug)
  CALL Equal_tensor(error_action, Op_psi_real, Op_psi_real_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="\hat{O}_{dense, 1/10}|0>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_real, out_unit, Size(Op_psi_real), info="\hat{O}_{dense, 1/10}|0>")
    CALL Write_Vec(Op_psi_real_ana, out_unit, Size(Op_psi_real_ana), info="\hat{O}_{dense, 1/10}|0>(Analitical)")
  END IF

  Op_psi_complex_ana = Op_psi_real_ana*EYE                                                               ! make as if the b_0 was a vector of the canonical basis set on \mathbb{C} i.e. i*b_0
  CALL Action(Op_psi_complex, EO2_diag_tenth, EYE*b_0, Verbose=Verbose, Debug=Debug)
  CALL Equal_tensor(error_action, Op_psi_complex, Op_psi_complex_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="\hat{O}_{diag, 1/10}|0>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_complex, out_unit, Size(Op_psi_complex), info="\hat{O}_{diag, 1/10}|0>")
    CALL Write_Vec(Op_psi_complex_ana, out_unit, Size(Op_psi_complex_ana), info="\hat{O}_{diag, 1/10}|0>(Analitical)")
  END IF
  CALL Action(Op_psi_complex, EO2_dense_tenth, EYE*b_0, Verbose=Verbose, Debug=Debug)
  CALL Equal_tensor(error_action, Op_psi_complex, Op_psi_complex_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="\hat{O}_{dense, 1/10}|0>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_complex, out_unit, Size(Op_psi_complex), info="\hat{O}_{dense, 1/10}|0>")
    CALL Write_Vec(Op_psi_complex_ana, out_unit, Size(Op_psi_complex_ana), info="\hat{O}_{dense, 1/10}|0>(Analitical)")
  END IF
      !---------------------------------second basis vector (i*)b_1--------------------------------
  CALL MolecCav_Construct_op_psi_ana_real(Op_psi_real_ana, EO2_diag_tenth, ONETENTH, b_1, Debug) ! same analytical result as with E01_dense_1
  CALL Action(Op_psi_real, EO2_diag_tenth, b_1, Verbose=Verbose, Debug=Debug)
  CALL Equal_tensor(error_action, Op_psi_real, Op_psi_real_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="\hat{O}_{diag, 1/10}|1>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_real, out_unit, Size(Op_psi_real), info="\hat{O}_{diag, 1/10}|1>")
    CALL Write_Vec(Op_psi_real_ana, out_unit, Size(Op_psi_real_ana), info="\hat{O}_{diag, 1/10}|1>(Analitical)")
  END IF
  CALL Action(Op_psi_real, EO2_dense_tenth, b_1, Verbose=Verbose, Debug=Debug)
  CALL Equal_tensor(error_action, Op_psi_real, Op_psi_real_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="\hat{O}_{dense, 1/10}|1>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_real, out_unit, Size(Op_psi_real), info="\hat{O}_{dense, 1/10}|1>")
    CALL Write_Vec(Op_psi_real_ana, out_unit, Size(Op_psi_real_ana), info="\hat{O}_{dense, 1/10}|1>(Analitical)")
  END IF

  Op_psi_complex_ana = Op_psi_real_ana*EYE                                                               ! make as if the b_1 was a vector of the canonical basis set on \mathbb{C} i.e. i*b_1
  CALL Action(Op_psi_complex, EO2_diag_tenth, EYE*b_1, Verbose=Verbose, Debug=Debug)
  CALL Equal_tensor(error_action, Op_psi_complex, Op_psi_complex_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="\hat{O}_{diag, 1/10}|1>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_complex, out_unit, Size(Op_psi_complex), info="\hat{O}_{diag, 1/10}|1>")
    CALL Write_Vec(Op_psi_complex_ana, out_unit, Size(Op_psi_complex_ana), info="\hat{O}_{diag, 1/10}|1>(Analitical)")
  END IF
  CALL Action(Op_psi_complex, EO2_dense_tenth, EYE*b_1, Verbose=Verbose, Debug=Debug)
  CALL Equal_tensor(error_action, Op_psi_complex, Op_psi_complex_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="\hat{O}_{dense, 1/10}|1>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_complex, out_unit, Size(Op_psi_complex), info="\hat{O}_{dense, 1/10}|1>")
    CALL Write_Vec(Op_psi_complex_ana, out_unit, Size(Op_psi_complex_ana), info="\hat{O}_{dense, 1/10}|1>(Analitical)")
  END IF
      !---------------------------------third basis vector (i*)b_2--------------------------------
  CALL MolecCav_Construct_op_psi_ana_real(Op_psi_real_ana, EO2_diag_tenth, ONETENTH, b_2, Debug) ! same analytical result as with E01_dense_1
  CALL Action(Op_psi_real, EO2_diag_tenth, b_2, Verbose=Verbose, Debug=Debug)
  CALL Equal_tensor(error_action, Op_psi_real, Op_psi_real_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="\hat{O}_{diag, 1/10}|2>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_real, out_unit, Size(Op_psi_real), info="\hat{O}_{diag, 1/10}|2>")
    CALL Write_Vec(Op_psi_real_ana, out_unit, Size(Op_psi_real_ana), info="\hat{O}_{diag, 1/10}|2>(Analitical)")
  END IF
  CALL Action(Op_psi_real, EO2_dense_tenth, b_2, Verbose=Verbose, Debug=Debug)
  CALL Equal_tensor(error_action, Op_psi_real, Op_psi_real_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="\hat{O}_{dense, 1/10}|2>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_real, out_unit, Size(Op_psi_real), info="\hat{O}_{dense, 1/10}|2>")
    CALL Write_Vec(Op_psi_real_ana, out_unit, Size(Op_psi_real_ana), info="\hat{O}_{dense, 1/10}|2>(Analitical)")
  END IF

  Op_psi_complex_ana = Op_psi_real_ana*EYE                                                               ! make as if the b_2 was a vector of the canonical basis set on \mathbb{C} i.e. i*b_2
  CALL Action(Op_psi_complex, EO2_diag_tenth, EYE*b_2, Verbose=Verbose, Debug=Debug)
  CALL Equal_tensor(error_action, Op_psi_complex, Op_psi_complex_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="\hat{O}_{diag, 1/10}|2>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_complex, out_unit, Size(Op_psi_complex), info="\hat{O}_{diag, 1/10}|2>")
    CALL Write_Vec(Op_psi_complex_ana, out_unit, Size(Op_psi_complex_ana), info="\hat{O}_{diag, 1/10}|2>(Analitical)")
  END IF
  CALL Action(Op_psi_complex, EO2_dense_tenth, EYE*b_2, Verbose=Verbose, Debug=Debug)
  CALL Equal_tensor(error_action, Op_psi_complex, Op_psi_complex_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="\hat{O}_{dense, 1/10}|2>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_complex, out_unit, Size(Op_psi_complex), info="\hat{O}_{dense, 1/10}|2>")
    CALL Write_Vec(Op_psi_complex_ana, out_unit, Size(Op_psi_complex_ana), info="\hat{O}_{dense, 1/10}|2>(Analitical)")
  END IF      
      !---------------------------------WF LC of the basis vectors--------------------------------
  CALL MolecCav_Construct_op_psi_ana_real(Op_psi_real_ana, EO2_diag_tenth, ONETENTH, Psi_1D_R1_real, Debug) ! same analytical result as with E01_dense_1
  CALL Action(Op_psi_real, EO2_diag_tenth, Psi_1D_R1_real, Verbose=Verbose, Debug=Debug)
  CALL Equal_tensor(error_action, Op_psi_real, Op_psi_real_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="\hat{O}_{diag, 1/10}|\Psi_{real}>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_real, out_unit, Size(Op_psi_real), info="\hat{O}_{diag, 1/10}|\Psi_{real}>")
    CALL Write_Vec(Op_psi_real_ana, out_unit, Size(Op_psi_real_ana), info="\hat{O}_{diag, 1/10}|\Psi_{real}>(Analitical)")
  END IF
  CALL Action(Op_psi_real, EO2_dense_tenth, Psi_1D_R1_real, Verbose=Verbose, Debug=Debug)
  CALL Equal_tensor(error_action, Op_psi_real, Op_psi_real_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="\hat{O}_{dense, 1/10}|\Psi_{real}>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_real, out_unit, Size(Op_psi_real), info="\hat{O}_{dense, 1/10}|\Psi_{real}>")
    CALL Write_Vec(Op_psi_real_ana, out_unit, Size(Op_psi_real_ana), info="\hat{O}_{dense, 1/10}|\Psi_{real}>(Analitical)")
  END IF

  CALL MolecCav_Construct_op_psi_ana_complex(Op_psi_complex_ana, EO2_diag_tenth, ONETENTH, Psi_1D_R1_complex, Debug) ! same analytical result as with E01_dense_1
  CALL Action(Op_psi_complex, EO2_diag_tenth, Psi_1D_R1_complex, Verbose=Verbose, Debug=Debug)
  CALL Equal_tensor(error_action, Op_psi_complex, Op_psi_complex_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="\hat{O}_{diag, 1/10}|\Psi_{complex}>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_complex, out_unit, Size(Op_psi_complex), info="\hat{O}_{diag, 1/10}|\Psi_{complex}>")
    CALL Write_Vec(Op_psi_complex_ana, out_unit, Size(Op_psi_complex_ana), info="\hat{O}_{diag, 1/10}|\Psi_{complex}>(Analitical)")
  END IF
  CALL Action(Op_psi_complex, EO2_dense_tenth, Psi_1D_R1_complex, Verbose=Verbose, Debug=Debug)
  CALL Equal_tensor(error_action, Op_psi_complex, Op_psi_complex_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="\hat{O}_{dense, 1/10}|\Psi_{complex}>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_complex, out_unit, Size(Op_psi_complex), info="\hat{O}_{dense, 1/10}|\Psi_{complex}>")
    CALL Write_Vec(Op_psi_complex_ana, out_unit, Size(Op_psi_complex_ana), info="\hat{O}_{dense, 1/10}|\Psi_{complex}>(Analitical)")
  END IF
    !--------------------------------(Band/Dense) Coeff = 1 (EO3)--------------------------------
      !---------------------------------first basis vector (i*)b_0--------------------------------
  CALL MolecCav_Construct_op_psi_ana_real(Op_psi_real_ana, EO3_band_1, ONE, b_0, Debug) ! same analytical result as with E01_dense_1
  CALL Action(Op_psi_real, EO3_band_1, b_0, Verbose=Verbose, Debug=Debug)
  CALL Equal_tensor(error_action, Op_psi_real, Op_psi_real_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="\hat{O}_{band, 1}|0>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_real, out_unit, Size(Op_psi_real), info="\hat{O}_{band, 1}|0>")
    CALL Write_Vec(Op_psi_real_ana, out_unit, Size(Op_psi_real_ana), info="\hat{O}_{band, 1}|0>(Analitical)")
  END IF
  CALL Action(Op_psi_real, EO3_dense_1, b_0, Verbose=Verbose, Debug=Debug)
  CALL Equal_tensor(error_action, Op_psi_real, Op_psi_real_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="\hat{O}_{dense, 1}|0>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_real, out_unit, Size(Op_psi_real), info="\hat{O}_{dense, 1}|0>")
    CALL Write_Vec(Op_psi_real_ana, out_unit, Size(Op_psi_real_ana), info="\hat{O}_{dense, 1}|0>(Analitical)")
  END IF

  Op_psi_complex_ana = Op_psi_real_ana*EYE                                                               ! make as if the b_0 was a vector of the canonical basis set on \mathbb{C} i.e. i*b_0
  CALL Action(Op_psi_complex, EO3_band_1, EYE*b_0, Verbose=Verbose, Debug=Debug)
  CALL Equal_tensor(error_action, Op_psi_complex, Op_psi_complex_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="\hat{O}_{band, 1}|0>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_complex, out_unit, Size(Op_psi_complex), info="\hat{O}_{band, 1}|0>")
    CALL Write_Vec(Op_psi_complex_ana, out_unit, Size(Op_psi_complex_ana), info="\hat{O}_{band, 1}|0>(Analitical)")
  END IF
  CALL Action(Op_psi_complex, EO3_dense_1, EYE*b_0, Verbose=Verbose, Debug=Debug)
  CALL Equal_tensor(error_action, Op_psi_complex, Op_psi_complex_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="\hat{O}_{dense, 1}|0>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_complex, out_unit, Size(Op_psi_complex), info="\hat{O}_{dense, 1}|0>")
    CALL Write_Vec(Op_psi_complex_ana, out_unit, Size(Op_psi_complex_ana), info="\hat{O}_{dense, 1}|0>(Analitical)")
  END IF
      !---------------------------------second basis vector (i*)b_1--------------------------------
  CALL MolecCav_Construct_op_psi_ana_real(Op_psi_real_ana, EO3_band_1, ONE, b_1, Debug) ! same analytical result as with E01_dense_1
  CALL Action(Op_psi_real, EO3_band_1, b_1, Verbose=Verbose, Debug=Debug)
  CALL Equal_tensor(error_action, Op_psi_real, Op_psi_real_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="\hat{O}_{band, 1}|1>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_real, out_unit, Size(Op_psi_real), info="\hat{O}_{band, 1}|1>")
    CALL Write_Vec(Op_psi_real_ana, out_unit, Size(Op_psi_real_ana), info="\hat{O}_{band, 1}|1>(Analitical)")
  END IF
  CALL Action(Op_psi_real, EO3_dense_1, b_1, Verbose=Verbose, Debug=Debug)
  CALL Equal_tensor(error_action, Op_psi_real, Op_psi_real_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="\hat{O}_{dense, 1}|1>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_real, out_unit, Size(Op_psi_real), info="\hat{O}_{dense, 1}|1>")
    CALL Write_Vec(Op_psi_real_ana, out_unit, Size(Op_psi_real_ana), info="\hat{O}_{dense, 1}|1>(Analitical)")
  END IF

  Op_psi_complex_ana = Op_psi_real_ana*EYE                                                               ! make as if the b_1 was a vector of the canonical basis set on \mathbb{C} i.e. i*b_1
  CALL Action(Op_psi_complex, EO3_band_1, EYE*b_1, Verbose=Verbose, Debug=Debug)
  CALL Equal_tensor(error_action, Op_psi_complex, Op_psi_complex_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="\hat{O}_{band, 1}|1>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_complex, out_unit, Size(Op_psi_complex), info="\hat{O}_{band, 1}|1>")
    CALL Write_Vec(Op_psi_complex_ana, out_unit, Size(Op_psi_complex_ana), info="\hat{O}_{band, 1}|1>(Analitical)")
  END IF
  CALL Action(Op_psi_complex, EO3_dense_1, EYE*b_1, Verbose=Verbose, Debug=Debug)
  CALL Equal_tensor(error_action, Op_psi_complex, Op_psi_complex_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="\hat{O}_{dense, 1}|1>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_complex, out_unit, Size(Op_psi_complex), info="\hat{O}_{dense, 1}|1>")
    CALL Write_Vec(Op_psi_complex_ana, out_unit, Size(Op_psi_complex_ana), info="\hat{O}_{dense, 1}|1>(Analitical)")
  END IF
      !---------------------------------third basis vector (i*)b_2--------------------------------
  CALL MolecCav_Construct_op_psi_ana_real(Op_psi_real_ana, EO3_band_1, ONE, b_2, Debug) ! same analytical result as with E01_dense_1
  CALL Action(Op_psi_real, EO3_band_1, b_2, Verbose=Verbose, Debug=Debug)
  CALL Equal_tensor(error_action, Op_psi_real, Op_psi_real_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="\hat{O}_{band, 1}|2>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_real, out_unit, Size(Op_psi_real), info="\hat{O}_{band, 1}|2>")
    CALL Write_Vec(Op_psi_real_ana, out_unit, Size(Op_psi_real_ana), info="\hat{O}_{band, 1}|2>(Analitical)")
  END IF
  CALL Action(Op_psi_real, EO3_dense_1, b_2, Verbose=Verbose, Debug=Debug)
  CALL Equal_tensor(error_action, Op_psi_real, Op_psi_real_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="\hat{O}_{dense, 1}|2>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_real, out_unit, Size(Op_psi_real), info="\hat{O}_{dense, 1}|2>")
    CALL Write_Vec(Op_psi_real_ana, out_unit, Size(Op_psi_real_ana), info="\hat{O}_{dense, 1}|2>(Analitical)")
  END IF

  Op_psi_complex_ana = Op_psi_real_ana*EYE                                                               ! make as if the b_2 was a vector of the canonical basis set on \mathbb{C} i.e. i*b_2
  CALL Action(Op_psi_complex, EO3_band_1, EYE*b_2, Verbose=Verbose, Debug=Debug)
  CALL Equal_tensor(error_action, Op_psi_complex, Op_psi_complex_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="\hat{O}_{band, 1}|2>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_complex, out_unit, Size(Op_psi_complex), info="\hat{O}_{band, 1}|2>")
    CALL Write_Vec(Op_psi_complex_ana, out_unit, Size(Op_psi_complex_ana), info="\hat{O}_{band, 1}|2>(Analitical)")
  END IF
  CALL Action(Op_psi_complex, EO3_dense_1, EYE*b_2, Verbose=Verbose, Debug=Debug)
  CALL Equal_tensor(error_action, Op_psi_complex, Op_psi_complex_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="\hat{O}_{dense, 1}|2>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_complex, out_unit, Size(Op_psi_complex), info="\hat{O}_{dense, 1}|2>")
    CALL Write_Vec(Op_psi_complex_ana, out_unit, Size(Op_psi_complex_ana), info="\hat{O}_{dense, 1}|2>(Analitical)")
  END IF      
      !---------------------------------WF LC of the basis vectors--------------------------------
  CALL MolecCav_Construct_op_psi_ana_real(Op_psi_real_ana, EO3_band_1, ONE, Psi_1D_R1_real, Debug) ! same analytical result as with E01_dense_1
  CALL Action(Op_psi_real, EO3_band_1, Psi_1D_R1_real, Verbose=Verbose, Debug=Debug)
  CALL Equal_tensor(error_action, Op_psi_real, Op_psi_real_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="\hat{O}_{band, 1}|\Psi_{real}>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_real, out_unit, Size(Op_psi_real), info="\hat{O}_{band, 1}|\Psi_{real}>")
    CALL Write_Vec(Op_psi_real_ana, out_unit, Size(Op_psi_real_ana), info="\hat{O}_{band, 1}|\Psi_{real}>(Analitical)")
  END IF
  CALL Action(Op_psi_real, EO3_dense_1, Psi_1D_R1_real, Verbose=Verbose, Debug=Debug)
  CALL Equal_tensor(error_action, Op_psi_real, Op_psi_real_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="\hat{O}_{dense, 1}|\Psi_{real}>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_real, out_unit, Size(Op_psi_real), info="\hat{O}_{dense, 1}|\Psi_{real}>")
    CALL Write_Vec(Op_psi_real_ana, out_unit, Size(Op_psi_real_ana), info="\hat{O}_{dense, 1}|\Psi_{real}>(Analitical)")
  END IF

  CALL MolecCav_Construct_op_psi_ana_complex(Op_psi_complex_ana, EO3_band_1, ONE, Psi_1D_R1_complex, Debug) ! same analytical result as with E01_dense_1
  CALL Action(Op_psi_complex, EO3_band_1, Psi_1D_R1_complex, Verbose=Verbose, Debug=Debug)
  CALL Equal_tensor(error_action, Op_psi_complex, Op_psi_complex_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="\hat{O}_{band, 1}|\Psi_{complex}>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_complex, out_unit, Size(Op_psi_complex), info="\hat{O}_{band, 1}|\Psi_{complex}>")
    CALL Write_Vec(Op_psi_complex_ana, out_unit, Size(Op_psi_complex_ana), info="\hat{O}_{band, 1}|\Psi_{complex}>(Analitical)")
  END IF
  CALL Action(Op_psi_complex, EO3_dense_1, Psi_1D_R1_complex, Verbose=Verbose, Debug=Debug)
  CALL Equal_tensor(error_action, Op_psi_complex, Op_psi_complex_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="\hat{O}_{dense, 1}|\Psi_{complex}>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_complex, out_unit, Size(Op_psi_complex), info="\hat{O}_{dense, 1}|\Psi_{complex}>")
    CALL Write_Vec(Op_psi_complex_ana, out_unit, Size(Op_psi_complex_ana), info="\hat{O}_{dense, 1}|\Psi_{complex}>(Analitical)")
  END IF
    !---------------------------------(Band/Dense) Coeff = PI (EO4)--------------------------------
      !---------------------------------first basis vector (i*)b_0--------------------------------
  CALL MolecCav_Construct_op_psi_ana_real(Op_psi_real_ana, EO4_band_pi, PI, b_0, Debug) ! same analytical result as with E01_dense_1
  CALL Action(Op_psi_real, EO4_band_pi, b_0, Verbose=Verbose, Debug=Debug)
  CALL Equal_tensor(error_action, Op_psi_real, Op_psi_real_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="\hat{O}_{band, \pi}|0>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_real, out_unit, Size(Op_psi_real), info="\hat{O}_{band, \pi}|0>")
    CALL Write_Vec(Op_psi_real_ana, out_unit, Size(Op_psi_real_ana), info="\hat{O}_{band, \pi}|0>(Analitical)")
  END IF
  CALL Action(Op_psi_real, EO4_dense_pi, b_0, Verbose=Verbose, Debug=Debug)
  CALL Equal_tensor(error_action, Op_psi_real, Op_psi_real_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="\hat{O}_{dense, \pi}|0>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_real, out_unit, Size(Op_psi_real), info="\hat{O}_{dense, \pi}|0>")
    CALL Write_Vec(Op_psi_real_ana, out_unit, Size(Op_psi_real_ana), info="\hat{O}_{dense, \pi}|0>(Analitical)")
  END IF

  Op_psi_complex_ana = Op_psi_real_ana*EYE                                                               ! make as if the b_0 was a vector of the canonical basis set on \mathbb{C} i.e. i*b_0
  CALL Action(Op_psi_complex, EO4_band_pi, EYE*b_0, Verbose=Verbose, Debug=Debug)
  CALL Equal_tensor(error_action, Op_psi_complex, Op_psi_complex_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="\hat{O}_{band, \pi}|0>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_complex, out_unit, Size(Op_psi_complex), info="\hat{O}_{band, \pi}|0>")
    CALL Write_Vec(Op_psi_complex_ana, out_unit, Size(Op_psi_complex_ana), info="\hat{O}_{band, \pi}|0>(Analitical)")
  END IF
  CALL Action(Op_psi_complex, EO4_dense_pi, EYE*b_0, Verbose=Verbose, Debug=Debug)
  CALL Equal_tensor(error_action, Op_psi_complex, Op_psi_complex_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="\hat{O}_{dense, \pi}|0>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_complex, out_unit, Size(Op_psi_complex), info="\hat{O}_{dense, \pi}|0>")
    CALL Write_Vec(Op_psi_complex_ana, out_unit, Size(Op_psi_complex_ana), info="\hat{O}_{dense, \pi}|0>(Analitical)")
  END IF
      !---------------------------------second basis vector (i*)b_1--------------------------------
  CALL MolecCav_Construct_op_psi_ana_real(Op_psi_real_ana, EO4_band_pi, PI, b_1, Debug) ! same analytical result as with E01_dense_1
  CALL Action(Op_psi_real, EO4_band_pi, b_1, Verbose=Verbose, Debug=Debug)
  CALL Equal_tensor(error_action, Op_psi_real, Op_psi_real_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="\hat{O}_{band, \pi}|1>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_real, out_unit, Size(Op_psi_real), info="\hat{O}_{band, \pi}|1>")
    CALL Write_Vec(Op_psi_real_ana, out_unit, Size(Op_psi_real_ana), info="\hat{O}_{band, \pi}|1>(Analitical)")
  END IF
  CALL Action(Op_psi_real, EO4_dense_pi, b_1, Verbose=Verbose, Debug=Debug)
  CALL Equal_tensor(error_action, Op_psi_real, Op_psi_real_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="\hat{O}_{dense, \pi}|1>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_real, out_unit, Size(Op_psi_real), info="\hat{O}_{dense, \pi}|1>")
    CALL Write_Vec(Op_psi_real_ana, out_unit, Size(Op_psi_real_ana), info="\hat{O}_{dense, \pi}|1>(Analitical)")
  END IF

  Op_psi_complex_ana = Op_psi_real_ana*EYE                                                               ! make as if the b_1 was a vector of the canonical basis set on \mathbb{C} i.e. i*b_1
  CALL Action(Op_psi_complex, EO4_band_pi, EYE*b_1, Verbose=Verbose, Debug=Debug)
  CALL Equal_tensor(error_action, Op_psi_complex, Op_psi_complex_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="\hat{O}_{band, \pi}|1>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_complex, out_unit, Size(Op_psi_complex), info="\hat{O}_{band, \pi}|1>")
    CALL Write_Vec(Op_psi_complex_ana, out_unit, Size(Op_psi_complex_ana), info="\hat{O}_{band, \pi}|1>(Analitical)")
  END IF
  CALL Action(Op_psi_complex, EO4_dense_pi, EYE*b_1, Verbose=Verbose, Debug=Debug)
  CALL Equal_tensor(error_action, Op_psi_complex, Op_psi_complex_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="\hat{O}_{dense, \pi}|1>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_complex, out_unit, Size(Op_psi_complex), info="\hat{O}_{dense, \pi}|1>")
    CALL Write_Vec(Op_psi_complex_ana, out_unit, Size(Op_psi_complex_ana), info="\hat{O}_{dense, \pi}|1>(Analitical)")
  END IF
      !---------------------------------third basis vector (i*)b_2--------------------------------
  CALL MolecCav_Construct_op_psi_ana_real(Op_psi_real_ana, EO4_band_pi, PI, b_2, Debug) ! same analytical result as with E01_dense_1
  CALL Action(Op_psi_real, EO4_band_pi, b_2, Verbose=Verbose, Debug=Debug)
  CALL Equal_tensor(error_action, Op_psi_real, Op_psi_real_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="\hat{O}_{band, \pi}|2>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_real, out_unit, Size(Op_psi_real), info="\hat{O}_{band, \pi}|2>")
    CALL Write_Vec(Op_psi_real_ana, out_unit, Size(Op_psi_real_ana), info="\hat{O}_{band, \pi}|2>(Analitical)")
  END IF
  CALL Action(Op_psi_real, EO4_dense_pi, b_2, Verbose=Verbose, Debug=Debug)
  CALL Equal_tensor(error_action, Op_psi_real, Op_psi_real_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="\hat{O}_{dense, \pi}|2>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_real, out_unit, Size(Op_psi_real), info="\hat{O}_{dense, \pi}|2>")
    CALL Write_Vec(Op_psi_real_ana, out_unit, Size(Op_psi_real_ana), info="\hat{O}_{dense, \pi}|2>(Analitical)")
  END IF

  Op_psi_complex_ana = Op_psi_real_ana*EYE                                                               ! make as if the b_2 was a vector of the canonical basis set on \mathbb{C} i.e. i*b_2
  CALL Action(Op_psi_complex, EO4_band_pi, EYE*b_2, Verbose=Verbose, Debug=Debug)
  CALL Equal_tensor(error_action, Op_psi_complex, Op_psi_complex_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="\hat{O}_{band, \pi}|2>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_complex, out_unit, Size(Op_psi_complex), info="\hat{O}_{band, \pi}|2>")
    CALL Write_Vec(Op_psi_complex_ana, out_unit, Size(Op_psi_complex_ana), info="\hat{O}_{band, \pi}|2>(Analitical)")
  END IF
  CALL Action(Op_psi_complex, EO4_dense_pi, EYE*b_2, Verbose=Verbose, Debug=Debug)
  CALL Equal_tensor(error_action, Op_psi_complex, Op_psi_complex_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="\hat{O}_{dense, \pi}|2>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_complex, out_unit, Size(Op_psi_complex), info="\hat{O}_{dense, \pi}|2>")
    CALL Write_Vec(Op_psi_complex_ana, out_unit, Size(Op_psi_complex_ana), info="\hat{O}_{dense, \pi}|2>(Analitical)")
  END IF      
      !---------------------------------WF LC of the basis vectors--------------------------------
  CALL MolecCav_Construct_op_psi_ana_real(Op_psi_real_ana, EO4_band_pi, PI, Psi_1D_R1_real, Debug) ! same analytical result as with E01_dense_1
  CALL Action(Op_psi_real, EO4_band_pi, Psi_1D_R1_real, Verbose=Verbose, Debug=Debug)
  CALL Equal_tensor(error_action, Op_psi_real, Op_psi_real_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="\hat{O}_{band, \pi}|\Psi_{real}>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_real, out_unit, Size(Op_psi_real), info="\hat{O}_{band, \pi}|\Psi_{real}>")
    CALL Write_Vec(Op_psi_real_ana, out_unit, Size(Op_psi_real_ana), info="\hat{O}_{band, \pi}|\Psi_{real}>(Analitical)")
  END IF
  CALL Action(Op_psi_real, EO4_dense_pi, Psi_1D_R1_real, Verbose=Verbose, Debug=Debug)
  CALL Equal_tensor(error_action, Op_psi_real, Op_psi_real_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="\hat{O}_{dense, \pi}|\Psi_{real}>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_real, out_unit, Size(Op_psi_real), info="\hat{O}_{dense, \pi}|\Psi_{real}>")
    CALL Write_Vec(Op_psi_real_ana, out_unit, Size(Op_psi_real_ana), info="\hat{O}_{dense, \pi}|\Psi_{real}>(Analitical)")
  END IF

  CALL MolecCav_Construct_op_psi_ana_complex(Op_psi_complex_ana, EO4_band_pi, PI, Psi_1D_R1_complex, Debug) ! same analytical result as with E01_dense_1
  CALL Action(Op_psi_complex, EO4_band_pi, Psi_1D_R1_complex, Verbose=Verbose, Debug=Debug)
  CALL Equal_tensor(error_action, Op_psi_complex, Op_psi_complex_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="\hat{O}_{band, \pi}|\Psi_{complex}>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_complex, out_unit, Size(Op_psi_complex), info="\hat{O}_{band, \pi}|\Psi_{complex}>")
    CALL Write_Vec(Op_psi_complex_ana, out_unit, Size(Op_psi_complex_ana), info="\hat{O}_{band, \pi}|\Psi_{complex}>(Analitical)")
  END IF
  CALL Action(Op_psi_complex, EO4_dense_pi, Psi_1D_R1_complex, Verbose=Verbose, Debug=Debug)
  CALL Equal_tensor(error_action, Op_psi_complex, Op_psi_complex_ana)
  CALL Logical_Test(test_action, error_action, test2=.FALSE., info="\hat{O}_{dense, \pi}|\Psi_{complex}>")
  IF (error_action .AND. Debug) THEN
    CALL Write_Vec(Op_psi_complex, out_unit, Size(Op_psi_complex), info="\hat{O}_{dense, \pi}|\Psi_{complex}>")
    CALL Write_Vec(Op_psi_complex_ana, out_unit, Size(Op_psi_complex_ana), info="\hat{O}_{dense, \pi}|\Psi_{complex}>(Analitical)")
  END IF

  !----------------------------Testing the Deallocation--------------------------- ! maybe should use a non initialized Elem_op so as not to have to change the test if the elem_op type is modified
      !----------------------------Diagonal guy---------------------------
  CALL Dealloc(EO1_diag_1, Verbose=Verbose, Debug=Debug)
  CALL Logical_Test(test_action, EO1_diag_1%Dense, test2=.FALSE., info="Dense")
  CALL Logical_Test(test_action, EO1_diag_1%Grid, test2=.FALSE., info="Grid")
  CALL Logical_Test(test_action, EO1_diag_1%Upper_bandwidth /= 0, test2=.FALSE., info="\hat{O}_{dense, 1}|0>")
  CALL Logical_Test(test_action, EO1_diag_1%Lower_bandwidth /= 0, test2=.FALSE., info="\hat{O}_{dense, 1}|0>")
  CALL Logical_Test(test_action, ALLOCATED(EO1_diag_1%Operator_type),  test2=.FALSE., info="\hat{O}_{dense, 1}|0>")
  CALL Logical_Test(test_action, ALLOCATED(EO1_diag_1%Diag_val),  test2=.FALSE., info="\hat{O}_{dense, 1}|0>")
  CALL Logical_Test(test_action, ALLOCATED(EO1_diag_1%Band_val),  test2=.FALSE., info="\hat{O}_{dense, 1}|0>")
  CALL Logical_Test(test_action, ALLOCATED(EO1_diag_1%Dense_val), test2=.FALSE., info="\hat{O}_{dense, 1}|0>")
      !----------------------------Dense guy---------------------------
  CALL Dealloc(EO1_dense_1, Verbose=Verbose, Debug=Debug)
  CALL Logical_Test(test_action, EO1_dense_1%Dense, test2=.FALSE., info="Dense")
  CALL Logical_Test(test_action, EO1_dense_1%Grid, test2=.FALSE., info="Grid")
  CALL Logical_Test(test_action, EO1_dense_1%Upper_bandwidth /= 0, test2=.FALSE., info="\hat{O}_{dense, 1}|0>")
  CALL Logical_Test(test_action, EO1_dense_1%Lower_bandwidth /= 0, test2=.FALSE., info="\hat{O}_{dense, 1}|0>")
  CALL Logical_Test(test_action, ALLOCATED(EO1_dense_1%Operator_type),  test2=.FALSE., info="\hat{O}_{dense, 1}|0>")
  CALL Logical_Test(test_action, ALLOCATED(EO1_dense_1%Diag_val),  test2=.FALSE., info="\hat{O}_{dense, 1}|0>")
  CALL Logical_Test(test_action, ALLOCATED(EO1_dense_1%Band_val),  test2=.FALSE., info="\hat{O}_{dense, 1}|0>")
  CALL Logical_Test(test_action, ALLOCATED(EO1_dense_1%Dense_val), test2=.FALSE., info="\hat{O}_{dense, 1}|0>")
      !----------------------------Band guy---------------------------
  CALL Dealloc(EO3_band_1, Verbose=Verbose, Debug=Debug)
  CALL Logical_Test(test_action, EO3_band_1%Dense, test2=.FALSE., info="Dense")
  CALL Logical_Test(test_action, EO3_band_1%Grid, test2=.FALSE., info="Grid")
  CALL Logical_Test(test_action, EO3_band_1%Upper_bandwidth /= 0, test2=.FALSE., info="\hat{O}_{dense, 1}|0>")
  CALL Logical_Test(test_action, EO3_band_1%Lower_bandwidth /= 0, test2=.FALSE., info="\hat{O}_{dense, 1}|0>")
  CALL Logical_Test(test_action, ALLOCATED(EO3_band_1%Operator_type),  test2=.FALSE., info="\hat{O}_{dense, 1}|0>")
  CALL Logical_Test(test_action, ALLOCATED(EO3_band_1%Diag_val),  test2=.FALSE., info="\hat{O}_{dense, 1}|0>")
  CALL Logical_Test(test_action, ALLOCATED(EO3_band_1%Band_val),  test2=.FALSE., info="\hat{O}_{dense, 1}|0>")
  CALL Logical_Test(test_action, ALLOCATED(EO3_band_1%Dense_val), test2=.FALSE., info="\hat{O}_{dense, 1}|0>")


  CALL Finalize_Test(test_action)
  
  
  CONTAINS


  SUBROUTINE MolecCav_Construct_Elem_op(Operator, Coeff, Operator_type, Dense, Debug_opt)                                        ! needed in this test file because Elem_op_m does not contains any Initialization procedure to stay independant on any basis representation choice (and cares only about the maths/matrices/actions of the operators) 
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
