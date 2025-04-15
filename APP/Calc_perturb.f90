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
! README:
! to be written soon
!==================================================================================================
!==================================================================================================
PROGRAM Calc_perturb
  USE QDUtil_m
  USE Algebra_m
  USE ND_indexes_m
  USE Mapping_m
  USE Cavity_mode_m
  USE Operator_1D_m
  USE Operator_2D_m
  USE Total_hamiltonian_m
  USE Psi_analysis_m
  IMPLICIT NONE


  logical, parameter            :: Debug = .TRUE.
  integer, parameter            :: Verbose = 0

  real(kind=Rkind)              :: w_mat, DT, m_mat, lambda, Cte

  real(kind=Rkind)              :: H_tot(3,3)

  real(kind=Rkind)              :: H_0_RC(3,3)
  real(kind=Rkind)              :: W_detu(3,3)
  real(kind=Rkind)              :: H_0_hRnC(3,3)
  real(kind=Rkind)              :: W_cplg(3,3)
  real(kind=Rkind)              :: H_0_G(3,3)
  real(kind=Rkind)              :: W_G(3,3)

  real(kind=Rkind)              :: REigval(3), Buffer(3)
  real(kind=Rkind)              :: REigvec(3,3), Corrections(3,3)

  real(kind=Rkind)              :: Energy


  WRITE(out_unit,*) "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
  WRITE(out_unit,*) "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx A(0) >> DT xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
  WRITE(out_unit,*) "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"

  w_mat  = 0.0058_Rkind
  DT     = 0.00002_Rkind
  m_mat  = 1744.60504565_Rkind
  lambda = 0.008_Rkind
  Cte    = ONE
  WRITE(out_unit,*) "--- w_mat   = "//TO_string(w_mat)
  WRITE(out_unit,*) "--- DT      = "//TO_string(DT)
  WRITE(out_unit,*) "--- m_mat   = "//TO_string(m_mat)
  WRITE(out_unit,*) "--- lambda  = "//TO_string(lambda)
  WRITE(out_unit,*) "--- Cte     = "//TO_string(Cte)
  WRITE(out_unit,*) "--- A(0)    = "//TO_string(  Couplings(lambda, Cte_loc=Cte, w_mat_loc=w_mat, m_mat_loc=m_mat, DT_loc=ZERO))
  WRITE(out_unit,*) "--- 2*A(0)  = "//TO_string(2*Couplings(lambda, Cte_loc=Cte, w_mat_loc=w_mat, m_mat_loc=m_mat, DT_loc=ZERO))
  WRITE(out_unit,*) "--- A(DT)   = "//TO_string(  Couplings(lambda, Cte_loc=Cte, w_mat_loc=w_mat, m_mat_loc=m_mat, DT_loc=DT))
  WRITE(out_unit,*) "--- 2*A(DT) = "//TO_string(2*Couplings(lambda, Cte_loc=Cte, w_mat_loc=w_mat, m_mat_loc=m_mat, DT_loc=DT))

  WRITE(out_unit,*)
  WRITE(out_unit,*) "-------------------------------------------------------------------------------"
  WRITE(out_unit,*) "-------------- _                                     _-------------------------"
  WRITE(out_unit,*) "--------------| w_mat+DT/2        0             0     |------------------------"
  WRITE(out_unit,*) "------H^tot = |     0      2*w_mat+3DT/2      A(DT)   |------------------------"
  WRITE(out_unit,*) "--------------|     0            A(DT)   2*w_mat+DT/2 |------------------------"
  WRITE(out_unit,*) "--------------|_                                     _|------------------------"
  WRITE(out_unit,*) "-------------------------------------------------------------------------------"

  H_tot      = ZERO
  H_tot(1,1) =   w_mat +   DT/2
  H_tot(2,2) = 2*w_mat + 3*DT/2
  H_tot(3,3) = 2*w_mat +   DT/2
  H_tot(2,3) = Couplings(lambda_loc=lambda, Cte_loc=Cte, w_mat_loc=w_mat, m_mat_loc=m_mat, DT_loc=DT)
  H_tot(3,2) = Couplings(lambda_loc=lambda, Cte_loc=Cte, w_mat_loc=w_mat, m_mat_loc=m_mat, DT_loc=DT)
  CALL Write_Mat(H_tot, out_unit, 3, info="H_tot")

  CALL diagonalization(H_tot, REigVal, REigVec)
  WRITE(out_unit,*); CALL Write_Vec(REigval, out_unit, 3, info="Energy levels tot")
  WRITE(out_unit,*) " Energy gap E_2 - E_1 = "//TO_string(REigval(3) - REigval(2)) 
  WRITE(out_unit,*); CALL Write_Mat(REigvec, out_unit, 3, info="\Psi_tot")

  WRITE(out_unit,*)
  WRITE(out_unit,*) "------------------------------------------------------------------------------------------------------------"
  WRITE(out_unit,*) "----------- _                             _            _                               _ -------------------"
  WRITE(out_unit,*) "-----------| w_mat       0             0   |          | DT/2        0            0      |-------------------"
  WRITE(out_unit,*) "---H^(0) = |   0      2*w_mat        A(0)  | W_detu = |   0       3DT/2     A(DT)-A(0)  |-------------------"
  WRITE(out_unit,*) "-----------|   0        A(0)       2*w_mat |          |   0     A(DT)-A(0)     DT/2     |-------------------"
  WRITE(out_unit,*) "-----------|_                             _|          |_                               _|-------------------"
  WRITE(out_unit,*) "------------------------------------------------------------------------------------------------------------"

  H_0_RC      = ZERO
  H_0_RC(1,1) = w_mat
  H_0_RC(2,2) = 2*w_mat
  H_0_RC(3,3) = 2*w_mat
  H_0_RC(2,3) = Couplings(lambda_loc=lambda, Cte_loc=Cte, w_mat_loc=w_mat, m_mat_loc=m_mat, DT_loc=ZERO)
  H_0_RC(3,2) = Couplings(lambda_loc=lambda, Cte_loc=Cte, w_mat_loc=w_mat, m_mat_loc=m_mat, DT_loc=ZERO)
  CALL Write_Mat(H_0_RC, out_unit, 3, info="H_0_RC")

  CALL diagonalization(H_0_RC, REigVal, REigVec)
  WRITE(out_unit,*); CALL Write_Vec(REigval, out_unit, 3, info="Energy levels (0)")
  WRITE(out_unit,*) " Energy gap (0) : |E_2 - E_1| = "//TO_string(REigval(3) - REigval(2)) 
  WRITE(out_unit,*); CALL Write_Mat(REigvec, out_unit, 3, info="\Psi^(0)")

  W_detu      = ZERO
  W_detu(1,1) =   DT / 2
  W_detu(2,2) = 3*DT / 2
  W_detu(3,3) =   DT / 2
  W_detu(2,3) = Couplings(lambda_loc=lambda, Cte_loc=Cte, w_mat_loc=w_mat, m_mat_loc=m_mat, DT_loc=DT) - &
              & Couplings(lambda_loc=lambda, Cte_loc=Cte, w_mat_loc=w_mat, m_mat_loc=m_mat, DT_loc=ZERO)
  W_detu(3,2) = Couplings(lambda_loc=lambda, Cte_loc=Cte, w_mat_loc=w_mat, m_mat_loc=m_mat, DT_loc=DT) - &
              & Couplings(lambda_loc=lambda, Cte_loc=Cte, w_mat_loc=w_mat, m_mat_loc=m_mat, DT_loc=ZERO)
  WRITE(out_unit,*); CALL Write_Mat(W_detu, out_unit, 3, info="W_detu")

  Corrections = matmul(TRANSPOSE(REigvec), matmul(W_detu, REigvec))
  WRITE(out_unit,*); CALL Write_Mat(Corrections, out_unit, 3, info="<Psi^(0)|W_detu|Psi^(0)>")
  CALL Write_Vec(REigval + [Corrections(1,1), Corrections(2,2), Corrections(3,3)], out_unit, 3, info="Energy levels (1) :")
  WRITE(out_unit,*) " Energy gap (1) : |E_2 - E_1| = "//TO_string(REigval(3)+Corrections(3,3) - REigval(2)-Corrections(2,2)) 


  WRITE(out_unit,*); WRITE(out_unit,*)
  WRITE(out_unit,*) "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
  WRITE(out_unit,*) "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx DT >> A(0) xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
  WRITE(out_unit,*) "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"

  w_mat  = 0.0058_Rkind
  DT     = 0.0002_Rkind
  m_mat  = 1744.60504565_Rkind
  lambda = 0.004_Rkind
  Cte    = ONE
  WRITE(out_unit,*) "--- w_mat   = "//TO_string(w_mat)
  WRITE(out_unit,*) "--- DT      = "//TO_string(DT)
  WRITE(out_unit,*) "--- m_mat   = "//TO_string(m_mat)
  WRITE(out_unit,*) "--- lambda  = "//TO_string(lambda)
  WRITE(out_unit,*) "--- Cte     = "//TO_string(Cte)
  WRITE(out_unit,*) "--- A(0)    = "//TO_string(  Couplings(lambda, Cte_loc=Cte, w_mat_loc=w_mat, m_mat_loc=m_mat, DT_loc=ZERO))
  WRITE(out_unit,*) "--- 2*A(0)  = "//TO_string(2*Couplings(lambda, Cte_loc=Cte, w_mat_loc=w_mat, m_mat_loc=m_mat, DT_loc=ZERO))
  WRITE(out_unit,*) "--- A(DT)   = "//TO_string(  Couplings(lambda, Cte_loc=Cte, w_mat_loc=w_mat, m_mat_loc=m_mat, DT_loc=DT))
  WRITE(out_unit,*) "--- 2*A(DT) = "//TO_string(2*Couplings(lambda, Cte_loc=Cte, w_mat_loc=w_mat, m_mat_loc=m_mat, DT_loc=DT))

  WRITE(out_unit,*)
  WRITE(out_unit,*) "-------------------------------------------------------------------------------"
  WRITE(out_unit,*) "-------------- _                                     _-------------------------"
  WRITE(out_unit,*) "--------------| w_mat+DT/2        0             0     |------------------------"
  WRITE(out_unit,*) "------H^tot = |     0      2*w_mat+3DT/2      A(DT)   |------------------------"
  WRITE(out_unit,*) "--------------|     0            A(DT)   2*w_mat+DT/2 |------------------------"
  WRITE(out_unit,*) "--------------|_                                     _|------------------------"
  WRITE(out_unit,*) "-------------------------------------------------------------------------------"

  H_tot      = ZERO
  H_tot(1,1) =   w_mat +   DT/2
  H_tot(2,2) = 2*w_mat + 3*DT/2
  H_tot(3,3) = 2*w_mat +   DT/2
  H_tot(2,3) = Couplings(lambda_loc=lambda, Cte_loc=Cte, w_mat_loc=w_mat, m_mat_loc=m_mat, DT_loc=DT)
  H_tot(3,2) = Couplings(lambda_loc=lambda, Cte_loc=Cte, w_mat_loc=w_mat, m_mat_loc=m_mat, DT_loc=DT)
  CALL Write_Mat(H_tot, out_unit, 3, info="H_tot")

  CALL diagonalization(H_tot, REigVal, REigVec)
  WRITE(out_unit,*); CALL Write_Vec(REigval, out_unit, 3, info="Energy levels tot")
  WRITE(out_unit,*) " Energy gap E_2 - E_1 = "//TO_string(REigval(3) - REigval(2)) 
  WRITE(out_unit,*); CALL Write_Mat(REigvec, out_unit, 3, info="\Psi_tot")

  WRITE(out_unit,*)
  WRITE(out_unit,*) "---------------------------------------------------------------------------------------------------------"
  WRITE(out_unit,*) "----------- _                                     _            _                         _ --------------"
  WRITE(out_unit,*) "-----------| w_mat+DT/2       0             0      |          | 0        0            0   |--------------"
  WRITE(out_unit,*) "---H^(0) = |   0        2*w_mat+3DT/2       0      | W_cplg = | 0        0          A(DT) |--------------"
  WRITE(out_unit,*) "-----------|   0              0       2*w_mat+DT/2 |          | 0      A(DT)          0   |--------------"
  WRITE(out_unit,*) "-----------|_                                     _|          |_                         _|--------------"
  WRITE(out_unit,*) "---------------------------------------------------------------------------------------------------------"

  H_0_hRnC      = ZERO
  H_0_hRnC(1,1) = w_mat   +   DT/2
  H_0_hRnC(2,2) = 2*w_mat + 3*DT/2
  H_0_hRnC(3,3) = 2*w_mat +   DT/2
  CALL Write_Mat(H_0_hRnC, out_unit, 3, info="H_0_hRnC")

  CALL diagonalization(H_0_hRnC, REigVal, REigVec)
  WRITE(out_unit,*); CALL Write_Vec(REigval, out_unit, 3, info="Energy levels (0)")
  WRITE(out_unit,*) " Energy gap (0) : |E_2 - E_1| = "//TO_string(REigval(3) - REigval(2)) 
  WRITE(out_unit,*); CALL Write_Mat(REigvec, out_unit, 3, info="\Psi^(0)")

  W_cplg      = ZERO
  W_cplg(2,3) = Couplings(lambda_loc=lambda, Cte_loc=Cte, w_mat_loc=w_mat, m_mat_loc=m_mat, DT_loc=DT)
  W_cplg(3,2) = Couplings(lambda_loc=lambda, Cte_loc=Cte, w_mat_loc=w_mat, m_mat_loc=m_mat, DT_loc=DT)
  WRITE(out_unit,*); CALL Write_Mat(W_cplg, out_unit, 3, info="W_cplg")

  Corrections = matmul(TRANSPOSE(REigvec), matmul(W_cplg, REigvec))
  WRITE(out_unit,*); CALL Write_Mat(Corrections, out_unit, 3, info="<Psi^(0)|W_cplg|Psi^(0)>")
  CALL Write_Vec(REigval + [Corrections(1,1), Corrections(2,2), Corrections(3,3)], out_unit, 3, info="Energy levels (1) :")
  WRITE(out_unit,*) " Energy gap (1) : |E_2 - E_1| = "//TO_string(REigval(3)+Corrections(3,3) - REigval(2)-Corrections(2,2)) 
  WRITE(out_unit,*) " Pas de correction à l'ordre 1 si la dégénérescence est déjà levée ? " 


  WRITE(out_unit,*); WRITE(out_unit,*)
  WRITE(out_unit,*) "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
  WRITE(out_unit,*) "xxxxxxxxxxxxxxxxxxxxxxxxx CAS GENERAL : 2 PERTURBATIONS xxxxxxxxxxxxxxxxxxxxxxxx"
  WRITE(out_unit,*) "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"

  w_mat  = 0.0058_Rkind
  DT     = 0.0002_Rkind
  m_mat  = 1744.60504565_Rkind
  lambda = 0.004_Rkind
  Cte    = ONE
  WRITE(out_unit,*) "--- w_mat   = "//TO_string(w_mat)
  WRITE(out_unit,*) "--- DT      = "//TO_string(DT)
  WRITE(out_unit,*) "--- m_mat   = "//TO_string(m_mat)
  WRITE(out_unit,*) "--- lambda  = "//TO_string(lambda)
  WRITE(out_unit,*) "--- Cte     = "//TO_string(Cte)
  WRITE(out_unit,*) "--- A(0)    = "//TO_string(  Couplings(lambda, Cte_loc=Cte, w_mat_loc=w_mat, m_mat_loc=m_mat, DT_loc=ZERO))
  WRITE(out_unit,*) "--- 2*A(0)  = "//TO_string(2*Couplings(lambda, Cte_loc=Cte, w_mat_loc=w_mat, m_mat_loc=m_mat, DT_loc=ZERO))
  WRITE(out_unit,*) "--- A(DT)   = "//TO_string(  Couplings(lambda, Cte_loc=Cte, w_mat_loc=w_mat, m_mat_loc=m_mat, DT_loc=DT))
  WRITE(out_unit,*) "--- 2*A(DT) = "//TO_string(2*Couplings(lambda, Cte_loc=Cte, w_mat_loc=w_mat, m_mat_loc=m_mat, DT_loc=DT))

  WRITE(out_unit,*)
  WRITE(out_unit,*) "-------------------------------------------------------------------------------"
  WRITE(out_unit,*) "-------------- _                                     _-------------------------"
  WRITE(out_unit,*) "--------------| w_mat+DT/2        0             0     |------------------------"
  WRITE(out_unit,*) "------H^tot = |     0      2*w_mat+3DT/2      A(DT)   |------------------------"
  WRITE(out_unit,*) "--------------|     0            A(DT)   2*w_mat+DT/2 |------------------------"
  WRITE(out_unit,*) "--------------|_                                     _|------------------------"
  WRITE(out_unit,*) "-------------------------------------------------------------------------------"

  H_tot      = ZERO
  H_tot(1,1) =   w_mat +   DT/2
  H_tot(2,2) = 2*w_mat + 3*DT/2
  H_tot(3,3) = 2*w_mat +   DT/2
  H_tot(2,3) = Couplings(lambda_loc=lambda, Cte_loc=Cte, w_mat_loc=w_mat, m_mat_loc=m_mat, DT_loc=DT)
  H_tot(3,2) = Couplings(lambda_loc=lambda, Cte_loc=Cte, w_mat_loc=w_mat, m_mat_loc=m_mat, DT_loc=DT)
  CALL Write_Mat(H_tot, out_unit, 3, info="H_tot")

  CALL diagonalization(H_tot, REigVal, REigVec)
  WRITE(out_unit,*); CALL Write_Vec(REigval, out_unit, 3, info="Energy levels tot")
  WRITE(out_unit,*) " Energy gap E_2 - E_1 = "//TO_string(REigval(3) - REigval(2)) 
  WRITE(out_unit,*); CALL Write_Mat(REigvec, out_unit, 3, info="\Psi_tot")

  WRITE(out_unit,*)
  WRITE(out_unit,*) "---------------------------------------------------------------------------------------------------------"
  WRITE(out_unit,*) "----------- _                                _         _                            _ -------------------"
  WRITE(out_unit,*) "-----------| w_mat       0             0      |       | DT/2        0            0   |-------------------"
  WRITE(out_unit,*) "---H^(0) = |   0        2*w_mat        0      | W_G = | 0         3DT/2        A(DT) |-------------------"
  WRITE(out_unit,*) "-----------|   0         0            2*w_mat |       | 0         A(DT)        DT/2  |-------------------"
  WRITE(out_unit,*) "-----------|_                                _|       |_                            _|-------------------"
  WRITE(out_unit,*) "---------------------------------------------------------------------------------------------------------"

  H_0_G      = ZERO
  H_0_G(1,1) = w_mat
  H_0_G(2,2) = 2*w_mat
  H_0_G(3,3) = 2*w_mat
  CALL Write_Mat(H_0_G, out_unit, 3, info="H_0_G")

  CALL diagonalization(H_0_G, REigVal, REigVec)
  WRITE(out_unit,*); CALL Write_Vec(REigval, out_unit, 3, info="Energy levels (0)")
  WRITE(out_unit,*) " Energy gap (0) : |E_2 - E_1| = "//TO_string(REigval(3) - REigval(2)) 
  WRITE(out_unit,*); CALL Write_Mat(REigvec, out_unit, 3, info="\Psi^(0)")
  Buffer = REigval

  W_G      = ZERO
  W_G(1,1) =   DT / 2
  W_G(2,2) = 3*DT / 2
  W_G(3,3) =   DT / 2
  W_G(2,3) = Couplings(lambda_loc=lambda, Cte_loc=Cte, w_mat_loc=w_mat, m_mat_loc=m_mat, DT_loc=DT)
  W_G(3,2) = Couplings(lambda_loc=lambda, Cte_loc=Cte, w_mat_loc=w_mat, m_mat_loc=m_mat, DT_loc=DT)
  WRITE(out_unit,*); CALL Write_Mat(W_G, out_unit, 3, info="W_G")

  WRITE(out_unit,*)
  CALL diagonalization(W_G, REigVal, REigVec)
  WRITE(out_unit,*); CALL Write_Vec(REigval, out_unit, 3, info="Energy corrections (1)")
  CALL Write_Mat(REigvec, out_unit, 3, info="REigvec W_G")

  WRITE(out_unit,*)
  CALL Write_Vec(REigval + Buffer, out_unit, 3, info="Energy levels (1)")
  WRITE(out_unit,*) " Energy gap (1) : |E_2 - E_1| = "//TO_string(REigval(3)+Buffer(3) - REigval(2)-Buffer(2)) 


  CONTAINS


  FUNCTION Couplings(lambda_loc, Cte_loc, w_mat_loc, m_mat_loc, DT_loc) RESULT (A_loc)
    USE QDUtil_m
    IMPLICIT NONE

    real(kind=Rkind), intent(in)  :: lambda_loc, Cte_loc, w_mat_loc, m_mat_loc, DT_loc  

    real(kind=Rkind)              :: A_loc

    A_loc = lambda_loc * Cte_loc * SQRT( (w_mat_loc + DT_loc) / (w_mat_loc * m_mat_loc) ) / 2

  END FUNCTION Couplings


END PROGRAM
