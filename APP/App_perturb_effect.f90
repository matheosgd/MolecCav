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
PROGRAM App_perturbations
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

  !--------------------------------------Diatomic molecule in a harmonic electonic potential-------------------------------------
  TYPE(Cavity_mode_t)           :: Molecule_1
  TYPE(Operator_1D_t)           :: Mol1H                                                                                         ! matrix of the one-dimensional harmonic Hamiltonian associated with HO D
  TYPE(Operator_1D_t)           :: Mol1DipMomt
  real(kind=Rkind)              :: CteMol1DipMomt = ONE                                                                          ! the intensity of the variation of the dipole moment with a variation of the matter DOF
  
  !-----------------------------------------------Cavity mode----------------------------------------------
  TYPE(Cavity_mode_t)           :: Cavity_mode
  TYPE(Operator_1D_t)           :: CavH                                                                                           ! matrix of the one-dimensional harmonic Hamiltonian associated with HO D
  TYPE(Operator_1D_t)           :: CavPosition
  
  real(kind=Rkind)              :: DT, A_0, A_DT, Ratio
  
  !-------------------------------------------------------Total Hamiltonian------------------------------------------------------
  real(kind=Rkind), allocatable :: TotH(:,:)

  !--------------------------------------------------Results - system properties-------------------------------------------------
  real(kind=Rkind), allocatable :: REigval(:)
  real(kind=Rkind), allocatable :: REigvec(:,:)

  !-----------------------------------------------------------Utilities----------------------------------------------------------
  integer                       :: I, Nb_M, Nb_C, NB


  !-----------------------------------------------------SYSTEM INITIALIZATION----------------------------------------------------
    !-------------------------------------Diatomic molecule in a harmonic electonic potential------------------------------------
  WRITE(out_unit,*) "-------------------------------------------------------SYSTEM INITIALIZATION--------------------------------&
                    &----------------------"
  WRITE(out_unit,*) "  ---------------------------------------Diatomic molecule in a harmonic electonic potential----------------&
                    &----------------------"
  CALL MolecCav_Read_cavity_mode(Mode=Molecule_1, nio=in_unit)

  WRITE(out_unit,*) "Molecular Hamiltonian"
  CALL Construct_Operator_1D(Operator=Mol1H,        operator_type="Hamiltonian", Mode=Molecule_1, Debug=.FALSE.)
  WRITE(out_unit,*) "Molecular Dipole moment"
  CALL Construct_Operator_1D(Operator=Mol1DipMomt,  operator_type="Position",    Mode=Molecule_1, Debug=.FALSE.)    ! initialized as a position operator because of approximation over its expression (cf. readme.md or manual)

  IF (ALLOCATED(Mol1DipMomt%Diag_val_R))  Mol1DipMomt%Diag_val_R  = Mol1DipMomt%Diag_val_R *CteMol1DipMomt                       ! /!\ so that the matrix already contains the intensity constant of the dipolar moment with the position of the matter (cf. manual for formulas)
  IF (ALLOCATED(Mol1DipMomt%Dense_val_R)) Mol1DipMomt%Dense_val_R = Mol1DipMomt%Dense_val_R*CteMol1DipMomt                       ! /!\ so that the matrix already contains the intensity constant of the dipolar moment with the position of the matter (cf. manual for formulas)
  IF (ALLOCATED(Mol1DipMomt%Band_val_R))  Mol1DipMomt%Band_val_R  = Mol1DipMomt%Band_val_R *CteMol1DipMomt                       ! /!\ so that the matrix already contains the intensity constant of the dipolar moment with the position of the matter (cf. manual for formulas)
    
  IF (.FALSE. .AND. ALLOCATED(Mol1DipMomt%Diag_val_R )) CALL Write_Vec(Mol1DipMomt%Diag_val_R,  out_unit, 3, info="Mol1DipMomt")
  IF (.FALSE. .AND. ALLOCATED(Mol1DipMomt%Band_val_R )) CALL Write_Mat(Mol1DipMomt%Band_val_R,  out_unit, 3, info="Mol1DipMomt")
  IF (.FALSE. .AND. ALLOCATED(Mol1DipMomt%Dense_val_R)) CALL Write_Mat(Mol1DipMomt%Dense_val_R, out_unit, 3, info="Mol1DipMomt")
  FLUSH(out_unit)

  !-----------------------------------------------Cavity mode----------------------------------------------
  WRITE(out_unit,*); WRITE(out_unit,*) "  -----------------------------------------------Zeroth cavity mode Reso/uncoupled-------&
  &----------------------------------------"

  CALL MolecCav_Read_cavity_mode(Mode=Cavity_mode, nio=in_unit)

  WRITE(out_unit,*) "Cavity mode Hamiltonian"
  CALL Construct_Operator_1D(Operator=CavH,        operator_type="Hamiltonian", Mode=Cavity_mode, Debug=.FALSE.)
  WRITE(out_unit,*) "Cavity mode Position"
  CALL Construct_Operator_1D(Operator=CavPosition, operator_type="Position",    Mode=Cavity_mode, Debug=.FALSE.)
  FLUSH(out_unit)

  Nb_M = Molecule_1%Nb
  Nb_C = Cavity_mode%Nb
  NB   = Molecule_1%Nb * Cavity_mode%Nb  

  DT = Cavity_mode%w-Molecule_1%w
  A_0  = Couplings(Cavity_mode%lambda, CteMol1DipMomt, Molecule_1%w, Molecule_1%m, ZERO)
  A_DT = Couplings(Cavity_mode%lambda, CteMol1DipMomt, Molecule_1%w, Molecule_1%m, DT)
  Ratio = DT/A_0

  !-------------------------------Construction of the Total Hamiltonian matrix------------------------------
  WRITE(out_unit,*); WRITE(out_unit,*)
  WRITE(out_unit,*) "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
  WRITE(out_unit,*) "xxxxxxxxxxxxxxxxxxxxxx DT/A(DT) = "//TO_string(Ratio)
  WRITE(out_unit,*) "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
  WRITE(out_unit,*); WRITE(out_unit,*) "  ----------------------------------Construction of the Total Hamiltonian matrix---------&
  &-------------------------"
  WRITE(out_unit,*) "--- w_mat   = "//TO_string(Molecule_1%w)
  WRITE(out_unit,*) "--- DT      = "//TO_string(DT)
  WRITE(out_unit,*) "--- m_mat   = "//TO_string(Molecule_1%m)
  WRITE(out_unit,*) "--- lambda  = "//TO_string(Cavity_mode%lambda)
  WRITE(out_unit,*) "--- Cte     = "//TO_string(CteMol1DipMomt)
  WRITE(out_unit,*) "--- A(0)    = "//TO_string(A_0)
  WRITE(out_unit,*) "--- 2xA(0)  = "//TO_string(2*A_0)
  WRITE(out_unit,*) "--- A(DT)   = "//TO_string(A_DT)//" (1st order energy correction of the 1st excited level)"
  WRITE(out_unit,*) "--- 2xA(DT) = "//TO_string(2*A_DT)

  ALLOCATE(TotH(NB, NB))
  CALL Construct_total_hamiltonian_1p1D_R1(TotH, CavPosition, CavH, Mol1DipMomt, Mol1H, Debug=.FALSE.)
  WRITE(out_unit,*); CALL Write_Mat(TotH, out_unit, NB, info="TotH")

    !-------------------------------------------------Computation of Eigenstates-------------------------------------------------
  WRITE(out_unit,*); WRITE(out_unit,*) "  ------------------------------------------------Computation of Eigenstates--------------&
  &----------------------------------"
  ALLOCATE(REigval(NB))
  ALLOCATE(REigvec(NB,NB))

  CALL diagonalization(TotH, REigval, REigvec)
  CALL Write_Vec(REigval, out_unit, 10, info="Energy levels TotH")
  WRITE(out_unit,*) " Energy level of the first excited state = "//TO_string(REigval(2))
  WRITE(out_unit,*) " Energy level of the secnd excited state = "//TO_string(REigval(3))
  WRITE(out_unit,*) " Energy gap : |E_2 - E_1| = "//TO_string(REigval(3) - REigval(2)) 
  WRITE(out_unit,*); CALL Write_Mat(REigvec, out_unit, SIZE(REigvec), info="\Psi_tot")
  

  CONTAINS


  FUNCTION Couplings(lambda_loc, Cte_loc, w_mat_loc, m_mat_loc, DT_loc) RESULT (A_loc)
    USE QDUtil_m
    IMPLICIT NONE

    real(kind=Rkind), intent(in)  :: lambda_loc, Cte_loc, w_mat_loc, m_mat_loc, DT_loc  

    real(kind=Rkind)              :: A_loc

    A_loc = lambda_loc * Cte_loc * SQRT( (w_mat_loc + DT_loc) / (w_mat_loc * m_mat_loc) ) / 2

  END FUNCTION Couplings


END PROGRAM
