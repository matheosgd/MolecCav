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
PROGRAM test_operator_ND_1p1D
  !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
  USE QDUtil_m
  USE Tests_m
  USE Algebra_m
  USE Operator_ND_m
  IMPLICIT NONE


  integer             :: Verbose = 40
  logical             :: Debug   = .TRUE.

  logical             :: Dense   = .FALSE.
  TYPE(Operator_ND_t) :: MatHxCavI
  TYPE(Operator_ND_t) :: MatIxCavH
  TYPE(Operator_ND_t) :: DipMomtxCavPos
  real(kind=Rkind)    :: Matw
  real(kind=Rkind)    :: Matm
  real(kind=Rkind)    :: Cavw
  real(kind=Rkind)    :: Cavlambda
  real(kind=Rkind)    :: CoeffDipMomt

  real(kind=Rkind)    :: b_0(6)                                                                ! six vectors of the HO basis set |00>, |10>, |20>, |01>, |11>, |21> 
  real(kind=Rkind)    :: b_1(6)
  real(kind=Rkind)    :: b_2(6)
  real(kind=Rkind)    :: b_3(6)
  real(kind=Rkind)    :: b_4(6)
  real(kind=Rkind)    :: b_5(6)

  real(kind=Rkind)    :: Coeff_0 = SQRT(TWO)
  real(kind=Rkind)    :: Coeff_1 = HALF
  real(kind=Rkind)    :: Coeff_2 = PI
  real(kind=Rkind)    :: Coeff_3 = ONETENTH
  real(kind=Rkind)    :: Coeff_4 = TWELVE
  real(kind=Rkind)    :: Coeff_5 = ONE
  real(kind=Rkind)    :: Coeffs(0:5)                                                             ! /!\ the indexes are here renamed to match the indexes of the basis vectors and coefficients ! the elements starts from 0 to 5 and not from 1 to 6 !!! /!\

  real(kind=Rkind)    :: Psi_R1_real(6)                                                                             ! an any vector representing the excitation state/wavefunction = a linear combination of the basis functions from the canonical basis set over \mathbb{R} /!\ Not normalized yet !
  real(kind=Rkind)    :: Op_psi_R1_real(6)                                                                                ! the resulting vector from the action of a 1D operator upon Psi_ND_R1_real
  real(kind=Rkind)    :: Ana_op_psi_R1_real(6)                                                                                ! the resulting vector from the action of a 1D operator upon Psi_ND_R1_real
  complex(kind=Rkind) :: Psi_R1_complex(6)                                                                          ! an any vector representing the excitation state/wavefunction = a linear combination of the basis functions from the canonical basis set over \mathbb{R} BUT with complexes expansion coefficients to make complexe WF. /!\ Not normalized yet !
  complex(kind=Rkind) :: Op_psi_R1_complex(6)                                                                             ! the resulting vector from the action of a 1D operator upon Psi_ND_R1_complex
  complex(kind=Rkind) :: Ana_op_psi_R1_complex(6)                                                                             ! the resulting vector from the action of a 1D operator upon Psi_ND_R1_complex

  TYPE(test_t)        :: test_opnd
  logical             :: error_opnd = .FALSE.

  integer             :: i, i_op


  !-----------------------------Test initialization----------------------------
  CALL Initialize_Test(test_opnd, test_name="OUT/test_file_opnd")


  !----------------------------------Wavefunction initialization---------------------------------
  b_5 = ZERO

  b_0 = b_5
  b_0(1) = ONE
  b_1 = b_5
  b_1(2) = ONE
  b_2 = b_5
  b_2(3) = ONE
  b_3 = b_5
  b_3(4) = ONE
  b_4 = b_5
  b_4(5) = ONE

  b_5(6) = ONE
  
  Coeffs = [Coeff_0, Coeff_1, Coeff_2, Coeff_3, Coeff_4, Coeff_5]                                ! /!\ the indexes are here renamed to match the indexes of the basis vectors and coefficients ! the elements starts from 0 to 5 and not from 1 to 6 !!! /!\

  Psi_R1_real    = Coeffs(0)*b_0 + Coeffs(1)*b_1 + Coeffs(2)*b_2 + Coeffs(3)*b_3 + Coeffs(4)*b_4 + Coeffs(5)*b_5
  Psi_R1_complex = Coeffs(0)*b_0 + Coeffs(1)*b_1 + Coeffs(2)*b_2 + Coeffs(3)*b_3 + Coeffs(4)*b_4 + Coeffs(5)*b_5
  CALL Normalize(Psi_R1_real)
  CALL Normalize(Psi_R1_complex)

  IF (Debug) THEN
    WRITE(out_unit,*)
    WRITE(out_unit,*); WRITE(out_unit,*) "---------------Basis vectors of the [Matter x Cavity] 1p1D system---------------"
    CALL Write_Vec(b_0, out_unit, Size(b_0, dim=1), info="b_0"); WRITE(out_unit,*)
    CALL Write_Vec(b_1, out_unit, Size(b_0, dim=1), info="b_1"); WRITE(out_unit,*)
    CALL Write_Vec(b_2, out_unit, Size(b_0, dim=1), info="b_2"); WRITE(out_unit,*)
    CALL Write_Vec(b_3, out_unit, Size(b_0, dim=1), info="b_3"); WRITE(out_unit,*)
    CALL Write_Vec(b_4, out_unit, Size(b_0, dim=1), info="b_4"); WRITE(out_unit,*)
    CALL Write_Vec(b_5, out_unit, Size(b_0, dim=1), info="b_5"); WRITE(out_unit,*)
    WRITE(out_unit,*); WRITE(out_unit,*) "-------------End basis vectors of the [Matter x Cavity] 1p1D system-------------"
    WRITE(out_unit,*); WRITE(out_unit,*) "----------------The any linear combinations of the basis functions---------------"
    CALL Write_Vec(Psi_R1_real,    out_unit, Size(Psi_R1_real,    dim=1), info="Psi_R1_real(normalized)")
    CALL Write_Vec(Psi_R1_complex, out_unit, Size(Psi_R1_complex, dim=1), info="Psi_R1_complex(normalized)")
    WRITE(out_unit,*); WRITE(out_unit,*) "------------------------End of the any linear combination-----------------------"
  END IF

  !-------------------------Operator_ND object initialization-------------------------
  CALL Initialize(MatHxCavI, "Hamiltonian", " identity ", in_unit, Dense=Dense, Verbose=Verbose, Debug=Debug)
  IF (Debug) THEN
    WRITE(out_unit,*)
    WRITE(out_unit,*) "--------------Operator_ND object constructed by MolecCav_Initialize_operator_ND--------------"
    CALL Write(MatHxCavI)
    WRITE(out_unit,*) "------------End Operator_ND object constructed by MolecCav_Initialize_operator_ND------------"
  END IF

  CALL Initialize(MatIxCavH, "Identity ", " hamiltonian", in_unit, Dense=Dense, Verbose=Verbose, Debug=Debug)
  IF (Debug) THEN
    WRITE(out_unit,*)
    WRITE(out_unit,*) "--------------Operator_ND object constructed by MolecCav_Initialize_operator_ND--------------"
    CALL Write(MatIxCavH)
    WRITE(out_unit,*) "------------End Operator_ND object constructed by MolecCav_Initialize_operator_ND------------"
  END IF

  CALL Initialize(DipMomtxCavPos, "DipMomt", "Position ", in_unit, Dense=Dense, Verbose=Verbose, Debug=Debug)
  IF (Debug) THEN
    WRITE(out_unit,*)
    WRITE(out_unit,*) "--------------Operator_ND object constructed by MolecCav_Initialize_operator_ND--------------"
    CALL Write(DipMomtxCavPos)
    WRITE(out_unit,*) "------------End Operator_ND object constructed by MolecCav_Initialize_operator_ND------------"
  END IF

  !----------------------------Testing the initialization---------------------------
  CALL Logical_Test(test_opnd, ANY(MatHxCavI%tab_indexes_mat_op/=[1]), test2=.FALSE., info="MatHxCavI%tab_mat_op")
  CALL Logical_Test(test_opnd, ANY(MatHxCavI%tab_indexes_cav_op/=[0]), test2=.FALSE., info="MatHxCavI%tab_cav_op")

  CALL Logical_Test(test_opnd, ANY(MatIxCavH%tab_indexes_mat_op/=[0]), test2=.FALSE., info="MatIxCavH%tab_mat_op")
  CALL Logical_Test(test_opnd, ANY(MatIxCavH%tab_indexes_cav_op/=[1]), test2=.FALSE., info="MatIxCavH%tab_cav_op")

  CALL Logical_Test(test_opnd, ANY(DipMomtxCavPos%tab_indexes_mat_op/=[4]), test2=.FALSE., info="DipMomtxCavPos%tab_mat_op")
  CALL Logical_Test(test_opnd, ANY(DipMomtxCavPos%tab_indexes_cav_op/=[2]), test2=.FALSE., info="DipMomtxCavPos%tab_cav_op")


  !----------------------------Testing the actions---------------------------
  CALL Action(Op_psi_R1_real, MatHxCavI,      Psi_R1_real, Debug=Debug)
  CALL Action(Op_psi_R1_real, MatIxCavH,      Op_psi_R1_real, Debug=Debug)
  CALL Action(Op_psi_R1_real, DipMomtxCavPos, Op_psi_R1_real, Debug=Debug)
  CALL Get(Matw, "w", "Matter", 1)
  CALL Get(Matm, "m", "Matter", 1)
  CALL Get(Cavw_loc, "w", "Cavity", 1)
  CALL Get(Cavlambda_loc, "lambda", "Cavity", 1)
  CoeffDipMomt_loc = ONE ! assumed linear here
  CALL Construct_ana_TotH_psi(Ana_op_psi_R1_real, Matw, Matm, Cavw_loc, Cavlambda_loc, CoeffDipMomt_loc, Psi_R1_real, Debug_opt=Debug)
  CALL Equal_R_R_matrix(error_opnd, Ana_op_psi_R1_real, Psi_R1_real)
  CALL Logical_Test(test_opnd, error_opnd, test2=.FALSE., info="TotH_Psi_R1_real")

  CALL Action(Op_psi_R1_complex, MatHxCavI,      Psi_R1_complex, Debug=Debug)
  CALL Action(Op_psi_R1_complex, MatIxCavH,      Op_psi_R1_complex, Debug=Debug)
  CALL Action(Op_psi_R1_complex, DipMomtxCavPos, Op_psi_R1_complex, Debug=Debug)
  CALL Construct_ana_TotH_psi(Ana_op_psi_R1_complex, Matw, Matm, Cavw_loc, Cavlambda_loc, CoeffDipMomt_loc, Psi_R1_complex, Debug_opt=Debug)
  CALL Equal_R_R_matrix(error_opnd, Ana_op_psi_R1_complex, Psi_R1_complex)
  CALL Logical_Test(test_opnd, error_opnd, test2=.FALSE., info="TotH_Psi_R1_complex")


  !----------------------------Testing the writing---------------------------
  CALL Write(MatHxCavI)

  
  !----------------------------Testing the deallocation---------------------------
  CALL Dealloc(MatHxCavI, Verbose=Verbose, Debug=Debug)
  CALL Logical_Test(test_opnd, ALLOCATED(MatHxCavI%tab_indexes_mat_op), test2=.FALSE., info="OpND%tab_mat deallocated ?")
  CALL Logical_Test(test_opnd, ALLOCATED(MatHxCavI%tab_indexes_cav_op), test2=.FALSE., info="OpND%tab_cav deallocated ?")
  
  CALL Finalize_Test(test_opnd)


  CONTAINS


  SUBROUTINE Construct_ana_TotH_psi_real(Ana_TotH_psi_R1_real, Matw_loc, Matm_loc, Cavw_loc, Cavlambda_loc, CoeffDipMomt_loc, Psi&
    &, Debug_opt)
    USE QDUtil_m
    USE Algebra_m
    IMPLICIT NONE 

    real(kind=Rkind),    intent(inout) :: Ana_TotH_psi_R1_real(3,2)
    real(kind=Rkind),    intent(in)    :: Matw_loc
    real(kind=Rkind),    intent(in)    :: Matm_loc
    real(kind=Rkind),    intent(in)    :: Cavw_loc
    real(kind=Rkind),    intent(in)    :: Cavlambda_loc
    real(kind=Rkind),    intent(in)    :: CoeffDipMomt_loc
    real(kind=Rkind),    intent(in)    :: Psi(3,2)
    logical, optional,   intent(in)    :: Debug_opt

    logical                            :: Debug_local = .FALSE.
    real(kind=Rkind)                   :: Norm_local

    IF (PRESENT(Debug_opt)) Debug_local = Debug_opt

    CALL Norm_of(Norm_local, Psi)

    Ana_TotH_psi_R1_real(1,1) = (  Matw_loc +   Cavm)*Psi(1,1) + Cavlambda_loc*CoeffDipMomt_loc*Psi(2,2)*SQRT(Cavm/(Matw_loc*Matm&
    &_loc))
    Ana_TotH_psi_R1_real(2,1) = (3*Matw_loc +   Cavm)*Psi(2,1) + Cavlambda_loc*CoeffDipMomt_loc*SQRT(Cavm/(Matw_loc*Matm_loc))*(P&
    &si(1,2) + SQRT(TWO)*Psi(3,2))
    Ana_TotH_psi_R1_real(3,1) = (5*Matw_loc +   Cavm)*Psi(3,1) + Cavlambda_loc*CoeffDipMomt_loc*SQRT(TWO)*Psi(2,2)*SQRT(Cavm/(Mat&
    &w_loc*Matm))
    Ana_TotH_psi_R1_real(1,2) = (  Matw_loc + 3*Cavm)*Psi(1,2) + Cavlambda_loc*CoeffDipMomt_loc*Psi(2,1)*SQRT(Cavm/(Matw_loc*Matm&
    &_loc))
    Ana_TotH_psi_R1_real(2,2) = (3*Matw_loc + 3*Cavm)*Psi(2,2) + Cavlambda_loc*CoeffDipMomt_loc*SQRT(Cavm/(Matw_loc*Matm_loc))*(P&
    &si(1,1) + SQRT(TWO)*Psi(3,1))
    Ana_TotH_psi_R1_real(3,2) = (5*Matw_loc + 3*Cavm)*Psi(3,2) + Cavlambda_loc*CoeffDipMomt_loc*SQRT(TWO)*Psi(2,1)*SQRT(Cavm/(Mat&
    &w_loc*Matm))
    
    Ana_TotH_psi_R1_real = Ana_TotH_psi_R1_real / (2*Norm_local)

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "Norm of the operand Psi = ", TO_string(Norm_local)
      WRITE(out_unit,*) "Resulting matrix of the TotH action over the operand :"
      CALL Write_Mat(Ana_TotH_psi_R1_real, out_unit, Size(Ana_TotH_psi_R1_real), info="Ana_TotH_psi_R1_real")
    END IF

  END SUBROUTINE


  SUBROUTINE Construct_ana_TotH_psi_complex(Ana_TotH_psi_R1_complex, Matw_loc, Matm_loc, Cavw_loc, Cavlambda_loc, CoeffDipMomt_loc, Psi&
    &, Debug_opt)
    USE QDUtil_m
    USE Algebra_m
    IMPLICIT NONE 

    complex(kind=Rkind),    intent(inout) :: Ana_TotH_psi_R1_complex(3,2)
    real(kind=Rkind),    intent(in)    :: Matw_loc
    real(kind=Rkind),    intent(in)    :: Matm_loc
    real(kind=Rkind),    intent(in)    :: Cavw_loc
    real(kind=Rkind),    intent(in)    :: Cavlambda_loc
    real(kind=Rkind),    intent(in)    :: CoeffDipMomt_loc
    real(kind=Rkind),    intent(in)    :: Psi(3,2)
    logical, optional,   intent(in)    :: Debug_opt

    logical                            :: Debug_local = .FALSE.
    real(kind=Rkind)                   :: Norm_local

    IF (PRESENT(Debug_opt)) Debug_local = Debug_opt

    CALL Norm_of(Norm_local, Psi)

    Ana_TotH_psi_R1_complex(1,1) = (  Matw_loc +   Cavm)*Psi(1,1) + Cavlambda_loc*CoeffDipMomt_loc*Psi(2,2)*SQRT(Cavm/(Matw_loc*M&
    &atm_loc))
    Ana_TotH_psi_R1_complex(2,1) = (3*Matw_loc +   Cavm)*Psi(2,1) + Cavlambda_loc*CoeffDipMomt_loc*SQRT(Cavm/(Matw_loc*Matm_loc))&
    &*(Psi(1,2) + SQRT(TWO)*Psi(3,2))
    Ana_TotH_psi_R1_complex(3,1) = (5*Matw_loc +   Cavm)*Psi(3,1) + Cavlambda_loc*CoeffDipMomt_loc*SQRT(TWO)*Psi(2,2)*SQRT(Cavm/(&
    &Matw_loc*Matm))
    Ana_TotH_psi_R1_complex(1,2) = (  Matw_loc + 3*Cavm)*Psi(1,2) + Cavlambda_loc*CoeffDipMomt_loc*Psi(2,1)*SQRT(Cavm/(Matw_loc*M&
    &atm_loc))
    Ana_TotH_psi_R1_complex(2,2) = (3*Matw_loc + 3*Cavm)*Psi(2,2) + Cavlambda_loc*CoeffDipMomt_loc*SQRT(Cavm/(Matw_loc*Matm_loc))&
    &*(Psi(1,1) + SQRT(TWO)*Psi(3,1))
    Ana_TotH_psi_R1_complex(3,2) = (5*Matw_loc + 3*Cavm)*Psi(3,2) + Cavlambda_loc*CoeffDipMomt_loc*SQRT(TWO)*Psi(2,1)*SQRT(Cavm/(&
    &Matw_loc*Matm))
    
    Ana_TotH_psi_R1_complex = Ana_TotH_psi_R1_complex / (2*Norm_local)

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "Norm of the operand Psi = ", TO_string(Norm_local)
      WRITE(out_unit,*) "Resulting matrix of the TotH action over the operand :"
      CALL Write_Mat(Ana_TotH_psi_R1_complex, out_unit, Size(Ana_TotH_psi_R1_complex), info="Ana_TotH_psi_R1_complex")
    END IF

  END SUBROUTINE


END PROGRAM
