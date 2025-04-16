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
  
  !-----------------------------------------------Zeroth cavity mode Reso/uncoupled----------------------------------------------
  TYPE(Cavity_mode_t)           :: Cavity_mode_RuC
  TYPE(Operator_1D_t)           :: CavH_RuC                                                                                         ! matrix of the one-dimensional harmonic Hamiltonian associated with HO D
  TYPE(Operator_1D_t)           :: CavPosition_RuC
  
  !--------------------------------------------------First cavity mode A(0) >> DT------------------------------------------------
  TYPE(Cavity_mode_t)           :: Cavity_mode_RC
  TYPE(Operator_1D_t)           :: CavH_RC                                                                                         ! matrix of the one-dimensional harmonic Hamiltonian associated with HO D
  TYPE(Operator_1D_t)           :: CavPosition_RC
  
  !--------------------------------------------------Secnd cavity mode DT >> A(0)------------------------------------------------
  TYPE(Cavity_mode_t)           :: Cavity_mode_oRuC
  TYPE(Operator_1D_t)           :: CavH_oRuC                                                                                         ! matrix of the one-dimensional harmonic Hamiltonian associated with HO D
  TYPE(Operator_1D_t)           :: CavPosition_oRuC

  real(kind=Rkind)              :: DT
  
  !-------------------------------------------------------Total Hamiltonian------------------------------------------------------
  real(kind=Rkind), allocatable :: TotH_RuC(:,:)
  real(kind=Rkind), allocatable :: TotH_RC(:,:)
  real(kind=Rkind), allocatable :: TotH_oRuC(:,:)

  !--------------------------------------------------Results - system properties-------------------------------------------------
  real(kind=Rkind), allocatable :: REigval(:)
  real(kind=Rkind), allocatable :: REigvec(:,:)
  real(kind=Rkind), allocatable :: Intensities(:,:)

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

  !-----------------------------------------------Zeroth cavity mode Reso/uncoupled----------------------------------------------
  WRITE(out_unit,*); WRITE(out_unit,*) "  -----------------------------------------------Zeroth cavity mode Reso/uncoupled-------&
  &----------------------------------------"

  CALL MolecCav_Read_cavity_mode(Mode=Cavity_mode_RuC, nio=in_unit)

  WRITE(out_unit,*) "Cavity mode Hamiltonian"
  CALL Construct_Operator_1D(Operator=CavH_RuC,        operator_type="Hamiltonian", Mode=Cavity_mode_RuC, Debug=.FALSE.)
  WRITE(out_unit,*) "Cavity mode Position"
  CALL Construct_Operator_1D(Operator=CavPosition_RuC, operator_type="Position",    Mode=Cavity_mode_RuC, Debug=.FALSE.)
  FLUSH(out_unit)

  !--------------------------------------------------First cavity mode A(0) >> DT------------------------------------------------
  WRITE(out_unit,*); WRITE(out_unit,*) "  --------------------------------------------------First cavity mode A(0) >> DT------&
  &------------------------------------------"
  
  CALL MolecCav_Read_cavity_mode(Mode=Cavity_mode_RC, nio=in_unit)

  WRITE(out_unit,*) "Cavity mode Hamiltonian"
  CALL Construct_Operator_1D(Operator=CavH_RC,        operator_type="Hamiltonian", Mode=Cavity_mode_RC, Debug=.FALSE.)
  WRITE(out_unit,*) "Cavity mode Position"
  CALL Construct_Operator_1D(Operator=CavPosition_RC, operator_type="Position",    Mode=Cavity_mode_RC, Debug=.FALSE.)
  FLUSH(out_unit)

  !--------------------------------------------------Secnd cavity mode DT >> A(0)------------------------------------------------
  WRITE(out_unit,*); WRITE(out_unit,*) "  --------------------------------------------------Secnd cavity mode DT >> A(0)----&
  &--------------------------------------------"
  
  CALL MolecCav_Read_cavity_mode(Mode=Cavity_mode_oRuC, nio=in_unit)

  WRITE(out_unit,*) "Cavity mode Hamiltonian"
  CALL Construct_Operator_1D(Operator=CavH_oRuC,           operator_type="Hamiltonian", Mode=Cavity_mode_oRuC, Debug=.FALSE.)
  WRITE(out_unit,*) "Cavity mode Position"
  CALL Construct_Operator_1D(Operator=CavPosition_oRuC, operator_type="Position",    Mode=Cavity_mode_oRuC, Debug=.FALSE.)
  FLUSH(out_unit)

  Nb_M = Molecule_1%Nb
  Nb_C = Cavity_mode_RC%Nb
  NB   = Molecule_1%Nb * Cavity_mode_RC%Nb  
                   
  !-------------------------------Construction of the Total Hamiltonian matrix with lambda = DT = 0------------------------------
  WRITE(out_unit,*); WRITE(out_unit,*)
  WRITE(out_unit,*) "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
  WRITE(out_unit,*) "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx lambda = DT = 0 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
  WRITE(out_unit,*) "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
  WRITE(out_unit,*); WRITE(out_unit,*) "  ----------------------------------Construction of the Total Hamiltonian matrix with lam&
  &bda = DT = 0----------------------------------"
  WRITE(out_unit,*) "--- w_mat   = "//TO_string(Molecule_1%w)
  DT = Cavity_mode_RC%w-Molecule_1%w
  WRITE(out_unit,*) "--- DT      = "//TO_string(DT)
  WRITE(out_unit,*) "--- m_mat   = "//TO_string(Molecule_1%m)
  WRITE(out_unit,*) "--- lambda  = "//TO_string(Cavity_mode_RuC%lambda)
  WRITE(out_unit,*) "--- Cte     = "//TO_string(CteMol1DipMomt)
  WRITE(out_unit,*) "--- A(0)    = "//TO_string(  Couplings(Cavity_mode_RuC%lambda, CteMol1DipMomt, Molecule_1%w, Molecule_1%m, ZE&
  &RO))
  WRITE(out_unit,*) "--- 2*A(0)  = "//TO_string(2*Couplings(Cavity_mode_RuC%lambda, CteMol1DipMomt, Molecule_1%w, Molecule_1%m, ZE&
  &RO))
  WRITE(out_unit,*) "--- A(DT)   = "//TO_string(  Couplings(Cavity_mode_RuC%lambda, CteMol1DipMomt, Molecule_1%w, Molecule_1%m, DT&
  &))//" (1st order energy correction of the 1st excited level)"
  WRITE(out_unit,*) "--- 2*A(DT) = "//TO_string(2*Couplings(Cavity_mode_RuC%lambda, CteMol1DipMomt, Molecule_1%w, Molecule_1%m, DT))

  ALLOCATE(TotH_RuC(NB, NB))
  CALL Construct_total_hamiltonian_1p1D_R1(TotH_RuC, CavPosition_RuC, CavH_RuC, Mol1DipMomt, Mol1H, Debug=.FALSE.)
  WRITE(out_unit,*); CALL Write_Mat(TotH_RuC, out_unit, NB, info="TotH_RuC")

    !-------------------------------------------------Computation of Eigenstates-------------------------------------------------
  WRITE(out_unit,*); WRITE(out_unit,*) "  ------------------------------------------------Computation of Eigenstates--------------&
  &----------------------------------"
  ALLOCATE(REigval(NB))
  ALLOCATE(REigvec(NB,NB))

  CALL diagonalization(TotH_RuC, REigval, REigvec)
  CALL Write_Vec(REigval, out_unit, 10, info="Energy levels TotH_RuC")
  WRITE(out_unit,*) " Energy gap : |E_2 - E_1| (RuC) = "//TO_string(REigval(3) - REigval(2)) 
  WRITE(out_unit,*); CALL Write_Mat(REigvec, out_unit, SIZE(REigvec), info="\Psi_tot")
  
    !--------------------------------------------Computation of transition intensities-------------------------------------------
!  WRITE(out_unit,*); WRITE(out_unit,*) "--------------------------------------------Computation of transition intensities------&
!                                       &-------------------------------------"
!  CALL Initialize_transition_matrix(Intensities, Mol1DipMomt, REigvec, Nb_states=10, Debug=.FALSE.)
!  CALL Compute_transition_matrix(Intensities,    Mol1DipMomt, REigvec, Debug=.FALSE.)
!
!  DO I = 1, 3
!    WRITE(out_unit,*) "Transition energy GSto"//TO_string(I)//" = "//TO_string((REigval(I+1)-REigval(1)))
!  END DO
!  
  DEALLOCATE(TotH_RuC); DEALLOCATE(REigval); DEALLOCATE(REigvec)!; DEALLOCATE(Intensities)


  !---------------------------------Construction of the Total Hamiltonian matrix with A(0) >> DT---------------------------------
  WRITE(out_unit,*); WRITE(out_unit,*)
  WRITE(out_unit,*) "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
  WRITE(out_unit,*) "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx A(0) >> DT xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
  WRITE(out_unit,*) "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
  WRITE(out_unit,*); WRITE(out_unit,*) "  ----------------------------------Construction of the Total Hamiltonian matrix with A(0)&
  & >> DT----------------------------------"
  WRITE(out_unit,*) "--- w_mat   = "//TO_string(Molecule_1%w)
  DT = Cavity_mode_RC%w-Molecule_1%w
  WRITE(out_unit,*) "--- DT      = "//TO_string(DT)
  WRITE(out_unit,*) "--- m_mat   = "//TO_string(Molecule_1%m)
  WRITE(out_unit,*) "--- lambda  = "//TO_string(Cavity_mode_RC%lambda)
  WRITE(out_unit,*) "--- Cte     = "//TO_string(CteMol1DipMomt)
  WRITE(out_unit,*) "--- A(0)    = "//TO_string(  Couplings(Cavity_mode_RC%lambda, CteMol1DipMomt, Molecule_1%w, Molecule_1%m, ZE&
  &RO))
  WRITE(out_unit,*) "--- 2*A(0)  = "//TO_string(2*Couplings(Cavity_mode_RC%lambda, CteMol1DipMomt, Molecule_1%w, Molecule_1%m, ZE&
  &RO))
  WRITE(out_unit,*) "--- A(DT)   = "//TO_string(  Couplings(Cavity_mode_RC%lambda, CteMol1DipMomt, Molecule_1%w, Molecule_1%m, DT&
  &))//" (1st order energy correction of the 1st excited level)"
  WRITE(out_unit,*) "--- 2*A(DT) = "//TO_string(2*Couplings(Cavity_mode_RC%lambda, CteMol1DipMomt, Molecule_1%w, Molecule_1%m, DT))

  ALLOCATE(TotH_RC(NB, NB))
  CALL Construct_total_hamiltonian_1p1D_R1(TotH_RC, CavPosition_RC, CavH_RC, Mol1DipMomt, Mol1H, Debug=.FALSE.)
  WRITE(out_unit,*); CALL Write_Mat(TotH_RC, out_unit, NB, info="TotH_RC")

    !-------------------------------------------------Computation of Eigenstates-------------------------------------------------
  WRITE(out_unit,*); WRITE(out_unit,*) "  ------------------------------------------------Computation of Eigenstates--------------&
  &----------------------------------"
  ALLOCATE(REigval(NB))
  ALLOCATE(REigvec(NB,NB))

  CALL diagonalization(TotH_RC, REigval, REigvec)
  CALL Write_Vec(REigval, out_unit, 10, info="Energy levels TotH_RC")
  WRITE(out_unit,*) " Energy gap : |E_2 - E_1| (RC) = "//TO_string(REigval(3) - REigval(2)) 
  WRITE(out_unit,*); CALL Write_Mat(REigvec, out_unit, 3, info="\Psi_tot")
  
    !--------------------------------------------Computation of transition intensities-------------------------------------------
!  WRITE(out_unit,*); WRITE(out_unit,*) "--------------------------------------------Computation of transition intensities------&
!                                       &-------------------------------------"
!  CALL Initialize_transition_matrix(Intensities, Mol1DipMomt, REigvec, Nb_states=10, Debug=.FALSE.)
!  CALL Compute_transition_matrix(Intensities,    Mol1DipMomt, REigvec, Debug=.FALSE.)
!
!  DO I = 1, 3
!    WRITE(out_unit,*) "Transition energy GSto"//TO_string(I)//" = "//TO_string((REigval(I+1)-REigval(1)))
!  END DO
!  
  DEALLOCATE(TotH_RC); DEALLOCATE(REigval); DEALLOCATE(REigvec)!; DEALLOCATE(Intensities)


  !---------------------------------Construction of the Total Hamiltonian matrix with A(0) >> DT---------------------------------
  WRITE(out_unit,*); WRITE(out_unit,*)
  WRITE(out_unit,*) "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
  WRITE(out_unit,*) "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx DT >> A(0) xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
  WRITE(out_unit,*) "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
  WRITE(out_unit,*); WRITE(out_unit,*) "  ----------------------------------Construction of the Total Hamiltonian matrix with A(0)&
  & >> DT----------------------------------"
  WRITE(out_unit,*) "--- w_mat  = "//TO_string(Molecule_1%w)
  DT = Cavity_mode_oRuC%w-Molecule_1%w
  WRITE(out_unit,*) "--- DT     = "//TO_string(DT)
  WRITE(out_unit,*) "--- m_mat  = "//TO_string(Molecule_1%m)
  WRITE(out_unit,*) "--- lambda = "//TO_string(Cavity_mode_oRuC%lambda)
  WRITE(out_unit,*) "--- Cte    = "//TO_string(CteMol1DipMomt)
  WRITE(out_unit,*) "--- A(0)   = "//TO_string(  Couplings(Cavity_mode_oRuC%lambda, CteMol1DipMomt, Molecule_1%w, Molecule_1%m, Z&
  &ERO))
  WRITE(out_unit,*) "--- 2*A(0) = "//TO_string(2*Couplings(Cavity_mode_oRuC%lambda, CteMol1DipMomt, Molecule_1%w, Molecule_1%m, Z&
  &ERO))
  WRITE(out_unit,*) "--- A(DT)  = "//TO_string(  Couplings(Cavity_mode_oRuC%lambda, CteMol1DipMomt, Molecule_1%w, Molecule_1%m, D&
  &T))//" (1st order energy correction of the 1st excited level)"
  WRITE(out_unit,*) "--- 2*A(DT) = "//TO_string(2*Couplings(Cavity_mode_oRuC%lambda, CteMol1DipMomt, Molecule_1%w, Molecule_1%m, &
  &DT))

  ALLOCATE(TotH_oRuC(NB, NB))
  CALL Construct_total_hamiltonian_1p1D_R1(TotH_oRuC, CavPosition_oRuC, CavH_oRuC, Mol1DipMomt, Mol1H, Debug=.FALSE.)
  WRITE(out_unit,*); CALL Write_Mat(TotH_oRuC, out_unit, NB, info="TotH_oRuC")

    !-------------------------------------------------Computation of Eigenstates-------------------------------------------------
  WRITE(out_unit,*); WRITE(out_unit,*) "  ------------------------------------------------Computation of Eigenstates-------------&
  &-----------------------------------"
  ALLOCATE(REigval(NB))
  ALLOCATE(REigvec(NB,NB))

  CALL diagonalization(TotH_oRuC, REigval, REigvec)
  CALL Write_Vec(REigval, out_unit, 10, info="Energy levels TotH_oRuC")
  WRITE(out_unit,*) " Energy gap : |E_2 - E_1| (oRuC) = "//TO_string(REigval(3) - REigval(2)) 
  WRITE(out_unit,*); CALL Write_Mat(REigvec, out_unit, 3, info="\Psi_tot")
  
    !--------------------------------------------Computation of transition intensities-------------------------------------------
!  WRITE(out_unit,*); WRITE(out_unit,*) "--------------------------------------------Computation of transition intensities------&
!                                       &-------------------------------------"
!  CALL Initialize_transition_matrix(Intensities, Mol1DipMomt, REigvec, Nb_states=10, Debug=.FALSE.)
!  CALL Compute_transition_matrix(Intensities,    Mol1DipMomt, REigvec, Debug=.FALSE.)
!
!  DO I = 1, 3
!    WRITE(out_unit,*) "Transition energy GSto"//TO_string(I)//" = "//TO_string((REigval(I+1)-REigval(1)))
!  END DO
  
  DEALLOCATE(TotH_oRuC); DEALLOCATE(REigval); DEALLOCATE(REigvec)!; DEALLOCATE(Intensities)


  CONTAINS


  FUNCTION Couplings(lambda_loc, Cte_loc, w_mat_loc, m_mat_loc, DT_loc) RESULT (A_loc)
    USE QDUtil_m
    IMPLICIT NONE

    real(kind=Rkind), intent(in)  :: lambda_loc, Cte_loc, w_mat_loc, m_mat_loc, DT_loc  

    real(kind=Rkind)              :: A_loc

    A_loc = lambda_loc * Cte_loc * SQRT( (w_mat_loc + DT_loc) / (w_mat_loc * m_mat_loc) ) / 2

  END FUNCTION Couplings


END PROGRAM
