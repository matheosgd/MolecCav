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
PROGRAM App_MolecCav
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


  logical, parameter            :: Debug = .FALSE.
  integer, parameter            :: Verbose = 0

  !--------------------------------------Diatomic molecule in a harmonic electonic potential-------------------------------------
  TYPE(Cavity_mode_t)           :: Molecule_1

  TYPE(Operator_1D_t)           :: Mol1H_opt                                                                                     ! matrix of the one-dimensional harmonic Hamiltonian associated with HO D
  TYPE(Operator_1D_t)           :: Mol1DipMomt_opt

  TYPE(Operator_1D_t)           :: Mol1H_dense                                                                                   ! matrix of the one-dimensional harmonic Hamiltonian associated with HO D
  TYPE(Operator_1D_t)           :: Mol1DipMomt_dense

  real(kind=Rkind)              :: CteMol1DipMomt = ONE                                                                          ! the intensity of the variation of the dipole moment with a variation of the matter DOF
  
  !-------------------------------------------------------First cavity mode------------------------------------------------------
  TYPE(Cavity_mode_t)           :: Cavity_mode_1

  TYPE(Operator_1D_t)           :: Cav1H_opt                                                                                     ! matrix of the one-dimensional harmonic Hamiltonian associated with HO D
  TYPE(Operator_1D_t)           :: Cav1Position_opt

  TYPE(Operator_1D_t)           :: Cav1H_dense                                                                                   ! matrix of the one-dimensional harmonic Hamiltonian associated with HO D
  TYPE(Operator_1D_t)           :: Cav1Position_dense

  !-------------------------------------------------------Total Hamiltonian------------------------------------------------------
  real(kind=Rkind), allocatable :: TotH_opt(:,:)
  real(kind=Rkind), allocatable :: TotH_dense(:,:)

  !--------------------------------------------------Results - system properties-------------------------------------------------
  real(kind=Rkind), allocatable :: REigval(:)
  real(kind=Rkind), allocatable :: REigvec(:,:)

  !-----------------------------------------------------------Utilities----------------------------------------------------------
  integer                       :: i, Nb_M, Nb_C, NB, N, j, k


  !--------------------------------------------------SYSTEM INITIALIZATION [OPT]-------------------------------------------------
    !-------------------------------------Diatomic molecule in a harmonic electonic potential------------------------------------
  WRITE(out_unit,*) "----------------------------------------------------SYSTEM INITIALIZATION [OPT]-----------------------------&
  &----------------------"
  WRITE(out_unit,*) "  ---------------------------------------Diatomic molecule in a harmonic electonic potential----------------&
  &----------------------"
  CALL time_perso("Beginning of time (opt)")

  CALL MolecCav_Read_cavity_mode(Mode=Molecule_1, nio=in_unit)

  WRITE(out_unit,*) "Molecular Hamiltonian   :"
  CALL Construct_Operator_1D(Operator=Mol1H_opt,        operator_type="Hamiltonian",                Mode=Molecule_1, Debug=Debug)
  WRITE(out_unit,*) "Molecular Dipole moment :"
  CALL Construct_Operator_1D(Operator=Mol1DipMomt_opt,  operator_type="Position",    Dense=.FALSE., Mode=Molecule_1, Debug=Debug)    ! initialized as a position operator because of approximation over its expression (cf. readme.md or manual)

  IF (ALLOCATED(Mol1DipMomt_opt%Diag_val_R))  Mol1DipMomt_opt%Diag_val_R  = Mol1DipMomt_opt%Diag_val_R  *CteMol1DipMomt                       ! /!\ so that the matrix already contains the intensity constant of the dipolar moment with the position of the matter (cf. manual for formulas)
  IF (ALLOCATED(Mol1DipMomt_opt%Dense_val_R)) Mol1DipMomt_opt%Dense_val_R = Mol1DipMomt_opt%Dense_val_R *CteMol1DipMomt                       ! /!\ so that the matrix already contains the intensity constant of the dipolar moment with the position of the matter (cf. manual for formulas)
  IF (ALLOCATED(Mol1DipMomt_opt%Band_val_R))  Mol1DipMomt_opt%Band_val_R  = Mol1DipMomt_opt%Band_val_R  *CteMol1DipMomt                       ! /!\ so that the matrix already contains the intensity constant of the dipolar moment with the position of the matter (cf. manual for formulas)
    
  IF (Debug .AND. ALLOCATED(Mol1DipMomt_opt%Diag_val_R )) CALL Write_Vec(Mol1DipMomt_opt%Diag_val_R,  out_unit, 3, info="Mol1DipM&
  &omt_opt")
  IF (Debug .AND. ALLOCATED(Mol1DipMomt_opt%Band_val_R )) CALL Write_Mat(Mol1DipMomt_opt%Band_val_R,  out_unit, 3, info="Mol1DipM&
  &omt_opt")
  IF (Debug .AND. ALLOCATED(Mol1DipMomt_opt%Dense_val_R)) CALL Write_Mat(Mol1DipMomt_opt%Dense_val_R, out_unit, 3, info="Mol1DipM&
  &omt_opt")
  FLUSH(out_unit)

    !------------------------------------------------------First cavity mode-----------------------------------------------------
  WRITE(out_unit,*); WRITE(out_unit,*) "  --------------------------------------------------------First cavity mode--------------&
                                       &-----------------------------------------"
  CALL MolecCav_Read_cavity_mode(Mode=Cavity_mode_1, nio=in_unit)

  WRITE(out_unit,*) "Cavity mode Hamiltonian :"
  CALL Construct_Operator_1D(Operator=Cav1H_opt,        operator_type="Hamiltonian",               Mode=Cavity_mode_1, Debug=Debug)
  WRITE(out_unit,*) "Cavity mode Position    :"
  CALL Construct_Operator_1D(Operator=Cav1Position_opt, operator_type="Position",   Dense=.FALSE., Mode=Cavity_mode_1, Debug=Debug)
  FLUSH(out_unit)

  Nb_M = Molecule_1%Nb
  Nb_C = Cavity_mode_1%Nb
  NB   = Molecule_1%Nb * Cavity_mode_1%Nb

  !------------------------------------------------------SYSTEM INITIALIZED------------------------------------------------------
  WRITE(out_unit,*) "------------------------------------------------------SYSTEM INITIALIZED------------------------------------&
  &------------------"

  
  !--------------------------------Construction of the Total Hamiltonian matrix with CM-couplings--------------------------------
  WRITE(out_unit,*); WRITE(out_unit,*) "----------------------------------Construction of the Total Hamiltonian matrix with CM-co&
                                       &uplings----------------------------------"

  ALLOCATE(TotH_opt(NB, NB))

  CALL Construct_total_hamiltonian_1p1D_R1(TotH_opt, Cav1Position_opt, Cav1H_opt, Mol1DipMomt_opt, Mol1H_opt, Debug=.FALSE.)
  IF (Debug) THEN
    WRITE(out_unit,*); WRITE(out_unit,*) "Total Hamiltonian 1p1D (lambda = "//TO_string(Cavity_mode_1%lambda)//" ; w_M = "//TO_st&
    &ring(Molecule_1%lambda)//" ; w_C = "//TO_string(Cavity_mode_1%w)//")"
    CALL Write_Mat(TotH_opt, out_unit, NB, info="TotH_opt")
  END IF

    !-------------------------------------------------Computation of Eigenstates-------------------------------------------------
  ALLOCATE(REigval(NB))
  ALLOCATE(REigvec(NB,NB))

  CALL diagonalization(TotH_opt, REigval, REigvec)
  WRITE(out_unit,*); WRITE(out_unit,*) 'EIGENVALUES'
  IF (Debug) CALL WRITE_Vec(REigval, out_unit, 10, info = 'VP_TotH_opt[Ha]')
  
  IF (Debug) THEN
    WRITE(out_unit,*); WRITE(out_unit,*) 'EIGENVECTORS'
    CALL WRITE_Mat(REigvec, out_unit, 6, info = 'Eigenvectors')
  END IF 

  CALL time_perso("End of computation (opt)")
  DEALLOCATE(Mol1H_opt%Diag_val_R, Mol1DipMomt_opt%Band_val_R, Cav1H_opt%Diag_val_R, Cav1Position_opt%Band_val_R, TotH_opt, REigv&
  &al, REigvec)

  !-------------------------------------------------SYSTEM INITIALIZATION [DENSE]------------------------------------------------
    !-------------------------------------Diatomic molecule in a harmonic electonic potential------------------------------------
  WRITE(out_unit,*) "------------------------------------------------------SYSTEM INITIALIZATION [DENSE]-------------------------&
  &-------------------------"
  WRITE(out_unit,*) "  ---------------------------------------Diatomic molecule in a harmonic electonic potential----------------&
  &----------------------"
!  CALL time_perso("Beginning of time (dense)") ! no need

  CALL MolecCav_Read_cavity_mode(Mode=Molecule_1, nio=in_unit)

  WRITE(out_unit,*) "Molecular Hamiltonian   :"
  CALL Construct_Operator_1D(Operator=Mol1H_dense,        operator_type="Hamiltonian", Dense=.TRUE., Mode=Molecule_1, Debug=Debug)
  WRITE(out_unit,*) "Molecular Dipole moment :"
  CALL Construct_Operator_1D(Operator=Mol1DipMomt_dense,  operator_type="Position",    Dense=.TRUE., Mode=Molecule_1, Debug=Debug)    ! initialized as a position operator because of approximation over its expression (cf. readme.md or manual)

  IF (ALLOCATED(Mol1DipMomt_dense%Diag_val_R))  Mol1DipMomt_dense%Diag_val_R  = Mol1DipMomt_dense%Diag_val_R  *CteMol1DipMomt         ! /!\ so that the matrix already contains the intensity constant of the dipolar moment with the position of the matter (cf. manual for formulas)
  IF (ALLOCATED(Mol1DipMomt_dense%Dense_val_R)) Mol1DipMomt_dense%Dense_val_R = Mol1DipMomt_dense%Dense_val_R *CteMol1DipMomt         ! /!\ so that the matrix already contains the intensity constant of the dipolar moment with the position of the matter (cf. manual for formulas)
  IF (ALLOCATED(Mol1DipMomt_dense%Band_val_R))  Mol1DipMomt_dense%Band_val_R  = Mol1DipMomt_dense%Band_val_R  *CteMol1DipMomt         ! /!\ so that the matrix already contains the intensity constant of the dipolar moment with the position of the matter (cf. manual for formulas)
    
  IF (Debug .AND. ALLOCATED(Mol1DipMomt_dense%Diag_val_R )) CALL Write_Vec(Mol1DipMomt_dense%Diag_val_R,  out_unit, 3, info="Mol1&
  &DipMomt_dense")
  IF (Debug .AND. ALLOCATED(Mol1DipMomt_dense%Band_val_R )) CALL Write_Mat(Mol1DipMomt_dense%Band_val_R,  out_unit, 3, info="Mol1&
  &DipMomt_dense")
  IF (Debug .AND. ALLOCATED(Mol1DipMomt_dense%Dense_val_R)) CALL Write_Mat(Mol1DipMomt_dense%Dense_val_R, out_unit, 3, info="Mol1&
  &DipMomt_dense")
  FLUSH(out_unit)

    !------------------------------------------------------First cavity mode-----------------------------------------------------
  WRITE(out_unit,*); WRITE(out_unit,*) "  --------------------------------------------------------First cavity mode--------------&
                                       &-----------------------------------------"
  CALL MolecCav_Read_cavity_mode(Mode=Cavity_mode_1, nio=in_unit)

  WRITE(out_unit,*) "Cavity mode Hamiltonian :"
  CALL Construct_Operator_1D(Operator=Cav1H_dense,        operator_type="Hamiltonian",Dense=.TRUE., Mode=Cavity_mode_1, Debug=Debug)
  WRITE(out_unit,*) "Cavity mode Position    :"
  CALL Construct_Operator_1D(Operator=Cav1Position_dense, operator_type="Position",   Dense=.TRUE., Mode=Cavity_mode_1, Debug=Debug)
  FLUSH(out_unit)

  Nb_M = Molecule_1%Nb
  Nb_C = Cavity_mode_1%Nb
  NB   = Molecule_1%Nb * Cavity_mode_1%Nb

  !------------------------------------------------------SYSTEM INITIALIZED------------------------------------------------------
  WRITE(out_unit,*) "------------------------------------------------------SYSTEM INITIALIZED------------------------------------&
  &------------------"

  
  !--------------------------------Construction of the Total Hamiltonian matrix with CM-couplings--------------------------------
  WRITE(out_unit,*); WRITE(out_unit,*) "----------------------------------Construction of the Total Hamiltonian matrix with CM-co&
                                       &uplings----------------------------------"

  ALLOCATE(TotH_dense(NB, NB))

  CALL Construct_total_hamiltonian_1p1D_R1(TotH_dense, Cav1Position_dense, Cav1H_dense, Mol1DipMomt_dense, Mol1H_dense)
  IF (Debug) THEN
    WRITE(out_unit,*); WRITE(out_unit,*) "Total Hamiltonian 1p1D (lambda = "//TO_string(Cavity_mode_1%lambda)//" ; w_M = "//TO_st&
    &ring(Molecule_1%lambda)//" ; w_C = "//TO_string(Cavity_mode_1%w)//")"
    CALL Write_Mat(TotH_dense, out_unit, NB, info="TotH_dense")
  END IF

    !-------------------------------------------------Computation of Eigenstates-------------------------------------------------
  ALLOCATE(REigval(NB))
  ALLOCATE(REigvec(NB,NB))

  CALL diagonalization(TotH_dense, REigval, REigvec)
  WRITE(out_unit,*); WRITE(out_unit,*) 'EIGENVALUES'
  IF (Debug) CALL WRITE_Vec(REigval, out_unit, 10, info = 'VP_TotH_dense[Ha]')
  
  IF (Debug) THEN
    WRITE(out_unit,*); WRITE(out_unit,*) 'EIGENVECTORS'
    CALL WRITE_Mat(REigvec, out_unit, 6, info = 'Eigenvectors')
  END IF 

  CALL time_perso("End of computation (dense)")


END PROGRAM