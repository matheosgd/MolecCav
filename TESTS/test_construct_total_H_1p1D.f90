!==================================================================================================
!==================================================================================================
!This file is part of MolecCav.
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
PROGRAM test_construct_total_H_1p1D
  !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
  USE QDUtil_m
  USE QDUtil_Test_m
  USE Algebra_m
  USE Cavity_mode_m
  USE Operator_1D_m
  USE Total_hamiltonian_m
  IMPLICIT NONE
  
  
  logical             :: Debug = .FALSE.
  
  TYPE(Cavity_mode_t) :: Matter_DOF                                                                ! DOF = the only Degree Of Freedom of the matter part of the system consider so far
  TYPE(Cavity_mode_t) :: Cavity_mode                                                               ! The well construction of the Operator's matrices is assumed checked by the dedicated test, so hard-coded references matrices are not needed here

  TYPE(Operator_1D_t) :: CavPosition                                                               ! position operator of the cavity mode
  TYPE(Operator_1D_t) :: CavH                                                                      ! Hamiltonian operator of the cavity mode
  TYPE(Operator_1D_t) :: Mat_dipolar_moment                                                        ! dipolar moment operator of the matter subsystem
  real(kind=Rkind)    :: Coeff_dipole_moment                                                       ! variation coefficient of the dipole moment with the matter DOF
  TYPE(Operator_1D_t) :: MatH                                                                      ! Hamiltonian operator of the matter subsystem
    
  real(kind=Rkind)    :: TotH(6,6)                                                                 ! size : NBxNB; NB = Nb_M*Nb_C (tensor product of the basis sets)
  real(kind=Rkind)    :: TotH_ana(6,6)                                                             ! the analytical matrix = hard-coded reference for comparison
  
  TYPE(test_t)        :: test_cnstrct_tot_H
  logical             :: error_cnstrct_tot_H = .FALSE.
    
  
  !--------------------------------------Test initialization-------------------------------------
  CALL Initialize_Test(test_cnstrct_tot_H, test_name="OUT/test_file_cnstrct_tot_H_1p1D")
  
  
  !-------------------------------------System initialization------------------------------------
  CALL Read_cavity_mode(Matter_DOF, nio=in_unit)
  Matter_DOF%Nb = 3
  IF (Debug) THEN
    WRITE(out_unit,*)
    WRITE(out_unit,*); WRITE(out_unit,*) "--------------Cavity mode constructed by MolecCav_Read_cavity_mode--------------"
    CALL Write_cavity_mode(Matter_DOF)
    WRITE(out_unit,*) "------------End Cavity mode constructed by MolecCav_Read_cavity_mode------------"
  END IF
  
  CALL Read_cavity_mode(Cavity_mode, nio=in_unit)
  Cavity_mode%Nb = 2
  IF (Debug) THEN
    WRITE(out_unit,*)
    WRITE(out_unit,*); WRITE(out_unit,*) "--------------Cavity mode constructed by MolecCav_Read_cavity_mode--------------"
    CALL Write_cavity_mode(Cavity_mode)
    WRITE(out_unit,*) "------------End Cavity mode constructed by MolecCav_Read_cavity_mode------------"
  END IF
  

  !-----------------TotH_{1p1D_(w_M=1)_(m_M=1)_(w_C=1)_(Coeff_\mu=1)_(lambda=0)}-----------------
    !---------------------------------1D operators initialization--------------------------------
  WRITE(out_unit,*)
  WRITE(out_unit,*) "-----------------TotH_{1p1D_(w_M=1)_(m_M=1)_(w_C=1)_(Coeff_\mu=1)_(lambda=0)}-----------------"
  Matter_DOF%w        = ONE
  Matter_DOF%m        = ONE
  Cavity_mode%w       = ONE
  Coeff_dipole_moment = ONE
  Cavity_mode%lambda  = ZERO

  CALL Construct_Operator_1D(Mat_dipolar_moment, "Position",    Mode=Matter_DOF,  Debug=Debug)     ! hypothesis : \mu(q) = \frac{\partial\mu}{\partial q}*\q ; \frac{\partial\mu}{\partial q} = Coeff_dip_momt
  Mat_dipolar_moment%Band_val_R = Coeff_dipole_moment*Mat_dipolar_moment%Band_val_R                ! => \hat{\mu} = Coeff_dip_momt*\hat{q}
  CALL Construct_Operator_1D(MatH,               "Hamiltonian", Mode=Matter_DOF,  Debug=Debug)
  CALL Construct_Operator_1D(CavPosition,        "Position",    Mode=Cavity_mode, Debug=Debug)
  CALL Construct_Operator_1D(CavH,               "Hamiltonian", Mode=Cavity_mode, Debug=Debug)
  
    !----------------------------------Testing the H construction----------------------------------
  CALL Construct_total_hamiltonian_1p1D(TotH, CavPosition, CavH, Mat_dipolar_moment, MatH, Debug=Debug)
  CALL Construct_ref_total_H_matrix_1p1D(TotH_ana, Matter_DOF, Cavity_mode, Coeff_dipole_moment, Debug_opt=Debug)
  CALL Equal_R_R_matrix(error_cnstrct_tot_H, TotH_ana, TotH)
  CALL Logical_Test(test_cnstrct_tot_H, error_cnstrct_tot_H, test2=.FALSE., info="TotH test")

  
  !----------------TotH_{1p1D_(w_M=1)_(m_M=1)_(w_C=1)_(Coeff_\mu=1)_(lambda=1/2)}----------------
    !---------------------------------1D operators initialization--------------------------------
  WRITE(out_unit,*)
  WRITE(out_unit,*) "-------------------TotH_{1p1D_(w_M=1)_(m_M=1)_(w_C=1)_(Coeff_\mu=1)_(lambda=1/2)}-------------------"
  CavH%lambda  = HALF; Cavity_mode%lambda = HALF                                                   ! /!\/!\/!\ The procedure use CavH%lambda and CavH%w as lambda and w_C, and not Mode%<> /!\ (but when we change w_C we rebuild CavH anyway). We still change the Mode one for the ref matrix building
  
    !------------------------------------Testing the H actions-----------------------------------
  CALL Construct_total_hamiltonian_1p1D(TotH, CavPosition, CavH, Mat_dipolar_moment, MatH, Debug=Debug)
  CALL Construct_ref_total_H_matrix_1p1D(TotH_ana, Matter_DOF, Cavity_mode, Coeff_dipole_moment, Debug_opt=Debug)
  CALL Equal_R_R_matrix(error_cnstrct_tot_H, TotH_ana, TotH)
  CALL Logical_Test(test_cnstrct_tot_H, error_cnstrct_tot_H, test2=.FALSE., info="TotH test")


  !------------TotH_{1p1D_(w_M=3*\sqrt(2))_(m_M=1)_(w_C=1)_(Coeff_\mu=1)_(lambda=1/2)}-----------
    !---------------------------------1D operators initialization--------------------------------
  WRITE(out_unit,*)
  WRITE(out_unit,*) "---------------TotH_{1p1D_(w_M=3*\sqrt(2))_(m_M=1)_(w_C=1)_(Coeff_\mu=1)_(lambda=1/2)}--------------"
  Matter_DOF%w    = THREE*SQRT(TWO)

  DEALLOCATE(Mat_dipolar_moment%Band_val_R); DEALLOCATE(Mat_dipolar_moment%Operator_type)
  DEALLOCATE(MatH%Diag_val_R);               DEALLOCATE(MatH%Operator_type)

  CALL Construct_Operator_1D(Mat_dipolar_moment, "Position",    Mode=Matter_DOF,  Debug=Debug)     ! hypothesis : \mu(q) = \frac{\partial\mu}{\partial q}*\q ; \frac{\partial\mu}{\partial q} = Coeff_dip_momt
  Mat_dipolar_moment%Band_val_R = Coeff_dipole_moment*Mat_dipolar_moment%Band_val_R                ! => \hat{\mu} = Coeff_dip_momt*\hat{q}
  CALL Construct_Operator_1D(MatH,               "Hamiltonian", Mode=Matter_DOF,  Debug=Debug)
  
    !------------------------------------Testing the H actions-----------------------------------
  CALL Construct_total_hamiltonian_1p1D(TotH, CavPosition, CavH, Mat_dipolar_moment, MatH, Debug=Debug)
  CALL Construct_ref_total_H_matrix_1p1D(TotH_ana, Matter_DOF, Cavity_mode, Coeff_dipole_moment, Debug_opt=Debug)
  CALL Equal_R_R_matrix(error_cnstrct_tot_H, TotH_ana, TotH)
  CALL Logical_Test(test_cnstrct_tot_H, error_cnstrct_tot_H, test2=.FALSE., info="TotH test")


  !----------TotH_{1p1D_(w_M=3*\sqrt(2))_(m_M=m(HF))_(w_C=1)_(Coeff_\mu=1)_(lambda=1/2)}---------
    !---------------------------------1D operators initialization--------------------------------
  WRITE(out_unit,*)
  WRITE(out_unit,*) "-------------TotH_{1p1D_(w_M=3*\sqrt(2))_(m_M=m(HF))_(w_C=1)_(Coeff_\mu=1)_(lambda=1/2)}------------"
  Matter_DOF%m    = 1744.60504565_Rkind

  DEALLOCATE(Mat_dipolar_moment%Band_val_R); DEALLOCATE(Mat_dipolar_moment%Operator_type)
  DEALLOCATE(MatH%Diag_val_R);               DEALLOCATE(MatH%Operator_type)

  CALL Construct_Operator_1D(Mat_dipolar_moment, "Position",    Mode=Matter_DOF,  Debug=Debug)     ! hypothesis : \mu(q) = \frac{\partial\mu}{\partial q}*\q ; \frac{\partial\mu}{\partial q} = Coeff_dip_momt
  Mat_dipolar_moment%Band_val_R = Coeff_dipole_moment*Mat_dipolar_moment%Band_val_R                ! => \hat{\mu} = Coeff_dip_momt*\hat{q}
  CALL Construct_Operator_1D(MatH,               "Hamiltonian", Mode=Matter_DOF,  Debug=Debug)
  
    !------------------------------------Testing the H actions-----------------------------------
  CALL Construct_total_hamiltonian_1p1D(TotH, CavPosition, CavH, Mat_dipolar_moment, MatH, Debug=Debug)
  CALL Construct_ref_total_H_matrix_1p1D(TotH_ana, Matter_DOF, Cavity_mode, Coeff_dipole_moment, Debug_opt=Debug)
  CALL Equal_R_R_matrix(error_cnstrct_tot_H, TotH_ana, TotH)
  CALL Logical_Test(test_cnstrct_tot_H, error_cnstrct_tot_H, test2=.FALSE., info="TotH test")
    

  !------TotH_{1p1D_(w_M=3*\sqrt(2))_(m_M=m(HF))_(w_C=\sqrt(2))_(Coeff_\mu=1)_(lambda=1/2)}------
    !---------------------------------1D operators initialization--------------------------------
  WRITE(out_unit,*)
  WRITE(out_unit,*) "---------TotH_{1p1D_(w_M=3*\sqrt(2))_(m_M=m(HF))_(w_C=\sqrt(2))_(Coeff_\mu=1)_(lambda=1/2)}---------"
  Cavity_mode%w     = SQRT(TWO)

  DEALLOCATE(CavPosition%Band_val_R); DEALLOCATE(CavPosition%Operator_type)
  DEALLOCATE(CavH%Diag_val_R);        DEALLOCATE(CavH%Operator_type)

  CALL Construct_Operator_1D(CavPosition, "Position",    Mode=Cavity_mode, Debug=Debug)
  CALL Construct_Operator_1D(CavH,        "Hamiltonian", Mode=Cavity_mode, Debug=Debug)
  
    !------------------------------------Testing the H actions-----------------------------------
  CALL Construct_total_hamiltonian_1p1D(TotH, CavPosition, CavH, Mat_dipolar_moment, MatH, Debug=Debug)
  CALL Construct_ref_total_H_matrix_1p1D(TotH_ana, Matter_DOF, Cavity_mode, Coeff_dipole_moment, Debug_opt=Debug)
  CALL Equal_R_R_matrix(error_cnstrct_tot_H, TotH_ana, TotH)
  CALL Logical_Test(test_cnstrct_tot_H, error_cnstrct_tot_H, test2=.FALSE., info="TotH test")
  

  !-----TotH_{1p1D_(w_M=3*\sqrt(2))_(m_M=m(HF))_(w_C=\sqrt(2))_(Coeff_\mu=\pi)_(lambda=1/2)}-----
    !---------------------------------1D operators initialization--------------------------------
  WRITE(out_unit,*)
  WRITE(out_unit,*) "--------TotH_{1p1D_(w_M=3*\sqrt(2))_(m_M=m(HF))_(w_C=\sqrt(2))_(Coeff_\mu=\pi)_(lambda=1/2)}--------"
  Mat_dipolar_moment%Band_val_R = Mat_dipolar_moment%Band_val_R / Coeff_dipole_moment
  
  Coeff_dipole_moment           = pi

  Mat_dipolar_moment%Band_val_R = Coeff_dipole_moment*Mat_dipolar_moment%Band_val_R                ! => \hat{\mu} = Coeff_dip_momt*\hat{q}
  
    !------------------------------------Testing the H actions-----------------------------------
  CALL Construct_total_hamiltonian_1p1D(TotH, CavPosition, CavH, Mat_dipolar_moment, MatH, Debug=Debug)
  CALL Construct_ref_total_H_matrix_1p1D(TotH_ana, Matter_DOF, Cavity_mode, Coeff_dipole_moment, Debug_opt=Debug)
  CALL Equal_R_R_matrix(error_cnstrct_tot_H, TotH_ana, TotH)
  CALL Logical_Test(test_cnstrct_tot_H, error_cnstrct_tot_H, test2=.FALSE., info="TotH test")
  

  CALL Finalize_Test(test_cnstrct_tot_H)
  
  
  CONTAINS
  
  
  SUBROUTINE Construct_ref_total_H_matrix_1p1D(TotH_ana_loc, Matter_DOF_loc, Cavity_mode_loc, Coeff_dip_mo_loc, Debug_opt)
    USE QDUtil_m
    USE Algebra_m
    IMPLICIT NONE 
  
    real(kind=Rkind),    intent(inout) :: TotH_ana_loc(6,6)
    TYPE(Cavity_mode_t), intent(in)    :: Matter_DOF_loc
    TYPE(Cavity_mode_t), intent(in)    :: Cavity_mode_loc
    real(kind=Rkind),    intent(in)    :: Coeff_dip_mo_loc
    logical, optional,   intent(in)    :: Debug_opt

    logical                            :: Debug_local = .FALSE.
    logical, parameter                 :: Debug_loc_submatrices = .FALSE.
    real(kind=Rkind)                   :: Norm_local
    real(kind=Rkind)                   :: MatterHO(6,6), CavityHO(6,6), Couplings(6,6)
    integer                            :: i

    IF (PRESENT(Debug_opt)) Debug_local = Debug_opt

    MatterHO = ZERO; CavityHO = ZERO; Couplings = ZERO

    MatterHO(1,1) = ONE
    MatterHO(2,2) = THREE
    MatterHO(3,3) = FIVE
    MatterHO(4,4) = ONE
    MatterHO(5,5) = THREE
    MatterHO(6,6) = FIVE

    CavityHO(1,1) = ONE
    CavityHO(2,2) = ONE
    CavityHO(3,3) = ONE
    CavityHO(4,4) = THREE
    CavityHO(5,5) = THREE
    CavityHO(6,6) = THREE

    Couplings(5,1) = ONE
    Couplings(4,2) = ONE
    Couplings(2,4) = ONE
    Couplings(1,5) = ONE

    Couplings(6,2) = SQRT(TWO)
    Couplings(5,3) = SQRT(TWO)
    Couplings(3,5) = SQRT(TWO)
    Couplings(2,6) = SQRT(TWO)

    IF (Debug_loc_submatrices) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "Intermediary Sub matrices for the TotH :"
      CALL Write_Mat(MatterHO, out_unit, Size(TotH_ana_loc), info="MatterHO")
      FLUSH(out_unit)
      WRITE(out_unit,*) "Intermediary Sub matrices for the TotH :"
      CALL Write_Mat(CavityHO, out_unit, Size(TotH_ana_loc), info="CavityHO")
      FLUSH(out_unit)
      WRITE(out_unit,*) "Intermediary Sub matrices for the TotH :"
      CALL Write_Mat(Couplings, out_unit, Size(TotH_ana_loc), info="Couplings")
      FLUSH(out_unit)
    END IF
  
    TotH_ana_loc = ( Matter_DOF_loc%w*MatterHO + Cavity_mode_loc%w*CavityHO + Cavity_mode_loc%lambd&
                    &a*Coeff_dip_mo_loc*SQRT(Cavity_mode_loc%w / (Matter_DOF_loc%w*Matter_DOF_loc%m&
                    &))*Couplings ) / 2

    IF (Debug_local) THEN
    WRITE(out_unit,*)
    WRITE(out_unit,*) "Constructed matrix of the TotH over the coupled basis [matter x cavity] (tensor producted) :"
    CALL Write_Mat(TotH_ana_loc, out_unit, Size(TotH_ana_loc), info="TotH_ana_loc")
    END IF

  END SUBROUTINE
  
  
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
    STOP "The two matrices must have same dimensions to compare them. Please, check initialization."
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
  