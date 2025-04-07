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


  logical, parameter            :: Debug = .TRUE.
  integer, parameter            :: Verbose = 0

  !--------------------------------------Diatomic molecule in a harmonic electonic potential-------------------------------------
  TYPE(Cavity_mode_t)           :: Molecule_1
  TYPE(Operator_1D_t)           :: Mol1H                                                                                         ! matrix of the one-dimensional harmonic Hamiltonian associated with HO D
  TYPE(Operator_1D_t)           :: Mol1Position
  TYPE(Operator_1D_t)           :: Mol1N                                                                                         ! here N does not count the number of photons but of excitation quanta in the vibrational state
  TYPE(Operator_1D_t)           :: Mol1DipMomt
  real(kind=Rkind)              :: CteMol1DipMomt = ONE                                                                          ! the intensity of the variation of the dipole moment with a variation of the matter DOF
  
  !-------------------------------------------------------First cavity mode------------------------------------------------------
  TYPE(Cavity_mode_t)           :: Cavity_mode_1
  TYPE(Operator_1D_t)           :: Cav1H                                                                                         ! matrix of the one-dimensional harmonic Hamiltonian associated with HO D
  TYPE(Operator_1D_t)           :: Cav1Position
  TYPE(Operator_1D_t)           :: Cav1N

  TYPE(Cavity_mode_t)           :: Cavity_mode_1_uncoupled
  TYPE(Operator_1D_t)           :: Cav1H_uncloupled                                                                              ! /!\ just need to change this one because Action TotH takes \lambda from the CavH

  !---------------------------------------------------------Wavefunctions--------------------------------------------------------
  real(kind=Rkind), allocatable :: Psi_1p1D_R2(:,:)                                                                              ! the total system (matter-cavity) wavefunction. Size Nb_M*Nb_C. |Psi_1p1D_R2> = |Molecule_WF>.TENSOR.|Cavity_WF> 
  real(kind=Rkind), allocatable :: Psi_1p1D_R1(:)
  real(kind=Rkind), allocatable :: CavPsi(:)

  !-------------------------------------------------------Total Hamiltonian------------------------------------------------------
  real(kind=Rkind), allocatable :: TotH_uncoupled(:,:)
  real(kind=Rkind)              :: MWH_uncoupled(2,2)                                                                            ! the mass-weighted Hessian matrix of the total, 1p1D, coupled system [matter x cavity] in harmonic approximation of the matter                     
  real(kind=Rkind), allocatable :: TotH(:,:)
  real(kind=Rkind)              :: MWH(2,2)                                                                                      ! the mass-weighted Hessian matrix of the total, 1p1D, coupled system [matter x cavity] in harmonic approximation of the matter                     

  !--------------------------------------------------Results - system properties-------------------------------------------------
  real(kind=Rkind), allocatable :: REigval_uncoupled(:)
  real(kind=Rkind), allocatable :: REigvec_uncoupled(:,:)
  real(kind=Rkind)              :: Normal_modes_uncoupled(2)                                                                     ! VP of the MWH
  real(kind=Rkind)              :: Normal_coordinates_uncoupled(2,2)                                                             ! \overrightarrow{VP} pf the MWH
  real(kind=Rkind), allocatable :: REigval(:)
  real(kind=Rkind), allocatable :: REigvec(:,:)
  real(kind=Rkind)              :: Normal_modes(2)                                                                               ! VP of the MWH
  real(kind=Rkind)              :: Normal_coordinates(2,2)                                                                       ! \overrightarrow{VP} pf the MWH

  !------------------------------------------------Results - temporary computation-----------------------------------------------
  real(kind=Rkind), allocatable :: Result_psi_1p1D_R1(:)
  real(kind=Rkind), allocatable :: Result_psi_1p1D_R2(:,:)
  real(kind=Rkind)              :: Average
  real(kind=Rkind), allocatable :: Intensities(:,:)
  real(kind=Rkind), allocatable :: Mol1Weights(:)
  real(kind=Rkind), allocatable :: Cav1Weights(:)
  real(kind=Rkind), allocatable :: TotH_bis(:,:)
  real(kind=Rkind), allocatable :: REigval_bis(:)
  real(kind=Rkind), allocatable :: REigvec_bis(:,:)
  real(kind=Rkind)              :: A(4,4)
  real(kind=Rkind)              :: B(16)

  !-----------------------------------------------------------Utilities----------------------------------------------------------
  integer                       :: i, Nb_M, Nb_C, NB, N, j, k
  TYPE(ND_indexes_t)            :: ND_indexes
  integer                       :: Ranks_sizes(2), List_indexes(2)
  logical                       :: Continue_loop = .TRUE.


  !-----------------------------------------------------SYSTEM INITIALIZATION----------------------------------------------------
    !-------------------------------------Diatomic molecule in a harmonic electonic potential------------------------------------
  WRITE(out_unit,*) "-------------------------------------------------------SYSTEM INITIALIZATION--------------------------------&
                    &----------------------"
  WRITE(out_unit,*) "  ---------------------------------------Diatomic molecule in a harmonic electonic potential----------------&
                    &----------------------"
  CALL MolecCav_Read_cavity_mode(Mode=Molecule_1, nio=in_unit)

  WRITE(out_unit,*) "Molecular Hamiltonian   :"
  CALL Construct_Operator_1D(Operator=Mol1H,        operator_type="Hamiltonian",                Mode=Molecule_1, Debug=Debug)
  WRITE(out_unit,*) "Molecular Position      :"
  CALL Construct_Operator_1D(Operator=Mol1Position, operator_type="Position",    Dense=.FALSE., Mode=Molecule_1, Debug=Debug)
  WRITE(out_unit,*) "Molecular Nb_photon (vibrationnal excitation quanta)"
  CALL Construct_Operator_1D(Operator=Mol1N,        operator_type="Nb_photons",                 Mode=Molecule_1, Debug=Debug)
  WRITE(out_unit,*) "Molecular Dipole moment :"
  CALL Construct_Operator_1D(Operator=Mol1DipMomt,  operator_type="Position",    Dense=.FALSE., Mode=Molecule_1, Debug=Debug)    ! initialized as a position operator because of approximation over its expression (cf. readme.md or manual)

  IF (ALLOCATED(Mol1DipMomt%Dense_val_R)) Mol1DipMomt%Dense_val_R = Mol1DipMomt%Dense_val_R*CteMol1DipMomt                       ! /!\ so that the matrix already contains the intensity constant of the dipolar moment with the position of the matter (cf. manual for formulas)
  IF (ALLOCATED(Mol1DipMomt%Band_val_R))  Mol1DipMomt%Band_val_R  = Mol1DipMomt%Band_val_R *CteMol1DipMomt                       ! /!\ so that the matrix already contains the intensity constant of the dipolar moment with the position of the matter (cf. manual for formulas)
    
  IF (Debug .AND. ALLOCATED(Mol1DipMomt%Diag_val_R )) CALL Write_Vec(Mol1DipMomt%Diag_val_R, out_unit, 3, info="Mol1DipMomt")
  IF (Debug .AND. ALLOCATED(Mol1DipMomt%Band_val_R )) CALL Write_Mat(Mol1DipMomt%Band_val_R, out_unit, 3, info="Mol1DipMomt")
  IF (Debug .AND. ALLOCATED(Mol1DipMomt%Dense_val_R)) CALL Write_Mat(Mol1DipMomt%Band_val_R, out_unit, 3, info="Mol1DipMomt")
  FLUSH(out_unit)

    !------------------------------------------------------First cavity mode-----------------------------------------------------
  WRITE(out_unit,*); WRITE(out_unit,*) "  --------------------------------------------------------First cavity mode--------------&
                                       &-----------------------------------------"
  CALL MolecCav_Read_cavity_mode(Mode=Cavity_mode_1, nio=in_unit)

  WRITE(out_unit,*) "Cavity mode Hamiltonian :"
  CALL Construct_Operator_1D(Operator=Cav1H,        operator_type="Hamiltonian",               Mode=Cavity_mode_1, Debug=Debug)
  WRITE(out_unit,*) "Cavity mode Position    :"
  CALL Construct_Operator_1D(Operator=Cav1Position, operator_type="Position",   Dense=.FALSE., Mode=Cavity_mode_1, Debug=.TRUE.)
  WRITE(out_unit,*) "Cavity mode Nb_photons  :"
  CALL Construct_Operator_1D(Operator=Cav1N,        operator_type="Nb_photons",                Mode=Cavity_mode_1, Debug=Debug)
  FLUSH(out_unit)

  Nb_M = Molecule_1%Nb
  Nb_C = Cavity_mode_1%Nb
  NB   = Molecule_1%Nb * Cavity_mode_1%Nb

  !--------------------------------Construction of the Total Hamiltonian matrix with CM-couplings--------------------------------
  WRITE(out_unit,*); WRITE(out_unit,*) "----------------------------------Construction of the Total Hamiltonian matrix with CM-co&
                                       &uplings----------------------------------"
  CALL time_perso("Beginning of time")

  ALLOCATE(TotH(NB, NB))
  CALL Construct_total_hamiltonian_1p1D_R1(TotH, Cav1Position, Cav1H, Mol1DipMomt, Mol1H, Debug=.FALSE.)
  CALL time_perso("TotH constructed")

  IF (Debug .AND. NB <= 20) THEN
    WRITE(out_unit,*); WRITE(out_unit,*) "Total Hamiltonian 1p1D (lambda /= 0, w_C /= w_M)"
    CALL Write_Mat(TotH, out_unit, NB, info="TotH")
  ELSE IF (Debug .AND. NB > 20) THEN
    WRITE(out_unit,*); WRITE(out_unit,*) "Total Hamiltonian 1p1D (lambda /= 0, w_C /= w_M) (50:50 slicing)"
    CALL Write_Mat(TotH(1:20,1:20), out_unit, 20, info="TotH (sliced)")
  END IF

  IF (Verbose > 0 ) WRITE(out_unit,*)
  IF (Verbose > 0 ) CALL Write_Mat(TotH(1:NB, 1:NB), out_unit, Size(TotH), info="TotH(FULL)")

    !-------------------------------------------------Computation of Eigenstates-------------------------------------------------
  ALLOCATE(REigval(NB))
  ALLOCATE(REigvec(NB,NB))
  CALL time_perso("Beginning diagonalization")

  CALL diagonalization(TotH, REigval, REigvec)
  CALL time_perso("TotH diagonalized")

  WRITE(out_unit,*); WRITE(out_unit,*) 'EIGENVALUES'
  CALL WRITE_Vec(REigval, out_unit, 10, info = 'VP_TotH[Ha]')

  IF (Debug) THEN
    WRITE(out_unit,*); WRITE(out_unit,*) 'EIGENVECTORS'
    CALL WRITE_Mat(REigvec, out_unit, 6, info = 'Eigenvectors')
  END IF 

    !--------------------------------------------Computation of transition intensities-------------------------------------------
  CALL time_perso("Beginning to compute the transition intensities matrix")
  
  CALL Initialize_transition_matrix(Intensities, Mol1DipMomt, REigvec, Nb_states=10, Debug=.TRUE.)
  CALL Compute_transition_matrix(Intensities,    Mol1DipMomt, REigvec, Debug=.TRUE.)
  CALL time_perso("Transition matrix computed")

  DO I = 1, 4
    WRITE(out_unit,*) "Transition energy GSto"//TO_string(I)//" = "//TO_string((REigval(I+1)-REigval(1)))
  END DO
  
  DEALLOCATE(TotH); DEALLOCATE(REigval); DEALLOCATE(REigvec)

END PROGRAM