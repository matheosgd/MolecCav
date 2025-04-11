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
  TYPE(Operator_1D_t)           :: Mol1Position
  TYPE(Operator_1D_t)           :: Mol1N                                                                                         ! here N does not count the number of photons but of excitation quanta in the vibrational state
  TYPE(Operator_1D_t)           :: Mol1DipMomt
  real(kind=Rkind)              :: CteMol1DipMomt = ONE                                                                          ! the intensity of the variation of the dipole moment with a variation of the matter DOF
  
  !-------------------------------------------------------First cavity mode------------------------------------------------------
  TYPE(Cavity_mode_t)           :: Cavity_mode_1
  TYPE(Operator_1D_t)           :: Cav1H                                                                                         ! matrix of the one-dimensional harmonic Hamiltonian associated with HO D
  TYPE(Operator_1D_t)           :: Cav1Position
  TYPE(Operator_1D_t)           :: Cav1N
  
  !---------------------------------------------------------Wavefunctions--------------------------------------------------------
  real(kind=Rkind), allocatable :: Psi_1p1D_R2(:,:)                                                                              ! the total system (matter-cavity) wavefunction. Size Nb_M*Nb_C. |Psi_1p1D_R2> = |Molecule_WF>.TENSOR.|Cavity_WF> 
  real(kind=Rkind), allocatable :: Psi_1p1D_R1(:)
  real(kind=Rkind), allocatable :: CavPsi(:)

  !-------------------------------------------------------Total Hamiltonian------------------------------------------------------
  real(kind=Rkind), allocatable :: TotH(:,:)
  real(kind=Rkind)              :: MWH(2,2)                                                                                      ! the mass-weighted Hessian matrix of the total, 1p1D, coupled system [matter x cavity] in harmonic approximation of the matter                     

  !--------------------------------------------------Results - system properties-------------------------------------------------
  real(kind=Rkind), allocatable :: REigval(:)
  real(kind=Rkind), allocatable :: REigvec(:,:)
  real(kind=Rkind)              :: Normal_modes(2)                                                                               ! VP of the MWH
  real(kind=Rkind)              :: Normal_coordinates(2,2)                                                                       ! \overrightarrow{VP} pf the MWH
  real(kind=Rkind), allocatable :: Intensities(:,:)

  !-----------------------------------------------------------Utilities----------------------------------------------------------
  integer                       :: i, Nb_M, Nb_C, NB, N, j, k


  !-----------------------------------------------------SYSTEM INITIALIZATION----------------------------------------------------
    !-------------------------------------Diatomic molecule in a harmonic electonic potential------------------------------------
  WRITE(out_unit,*) "-------------------------------------------------------SYSTEM INITIALIZATION--------------------------------&
                    &----------------------"
  WRITE(out_unit,*) "  ---------------------------------------Diatomic molecule in a harmonic electonic potential----------------&
                    &----------------------"
  CALL MolecCav_Read_cavity_mode(Mode=Molecule_1, nio=in_unit)

  WRITE(out_unit,*) "Molecular Hamiltonian"
  CALL Construct_Operator_1D(Operator=Mol1H,        operator_type="Hamiltonian", Mode=Molecule_1, Debug=.FALSE.)
  WRITE(out_unit,*) "Molecular Position"
  CALL Construct_Operator_1D(Operator=Mol1Position, operator_type="Position",    Mode=Molecule_1, Debug=.FALSE.)
  WRITE(out_unit,*) "Molecular Nb_photon (vibrationnal excitation quanta)"
  CALL Construct_Operator_1D(Operator=Mol1N,        operator_type="Nb_photons",  Mode=Molecule_1, Debug=.FALSE.)
  WRITE(out_unit,*) "Molecular Dipole moment"
  CALL Construct_Operator_1D(Operator=Mol1DipMomt,  operator_type="Position",    Mode=Molecule_1, Debug=.FALSE.)    ! initialized as a position operator because of approximation over its expression (cf. readme.md or manual)

  IF (ALLOCATED(Mol1DipMomt%Dense_val_R)) Mol1DipMomt%Dense_val_R = Mol1DipMomt%Dense_val_R*CteMol1DipMomt                       ! /!\ so that the matrix already contains the intensity constant of the dipolar moment with the position of the matter (cf. manual for formulas)
  IF (ALLOCATED(Mol1DipMomt%Band_val_R))  Mol1DipMomt%Band_val_R  = Mol1DipMomt%Band_val_R *CteMol1DipMomt                       ! /!\ so that the matrix already contains the intensity constant of the dipolar moment with the position of the matter (cf. manual for formulas)
    
  IF (.FALSE. .AND. ALLOCATED(Mol1DipMomt%Diag_val_R )) CALL Write_Vec(Mol1DipMomt%Diag_val_R, out_unit, 3, info="Mol1DipMomt")
  IF (.FALSE. .AND. ALLOCATED(Mol1DipMomt%Band_val_R )) CALL Write_Mat(Mol1DipMomt%Band_val_R, out_unit, 3, info="Mol1DipMomt")
  IF (.FALSE. .AND. ALLOCATED(Mol1DipMomt%Dense_val_R)) CALL Write_Mat(Mol1DipMomt%Band_val_R, out_unit, 3, info="Mol1DipMomt")
  FLUSH(out_unit)

    !------------------------------------------------------First cavity mode-----------------------------------------------------
  WRITE(out_unit,*); WRITE(out_unit,*) "  --------------------------------------------------------First cavity mode--------------&
                                       &-----------------------------------------"
  CALL MolecCav_Read_cavity_mode(Mode=Cavity_mode_1, nio=in_unit)

  WRITE(out_unit,*) "Cavity mode Hamiltonian"
  CALL Construct_Operator_1D(Operator=Cav1H,        operator_type="Hamiltonian", Mode=Cavity_mode_1, Debug=.FALSE.)
  WRITE(out_unit,*) "Cavity mode Position"
  CALL Construct_Operator_1D(Operator=Cav1Position, operator_type="Position",    Mode=Cavity_mode_1, Debug=.FALSE.)
  WRITE(out_unit,*) "Cavity mode Nb_photons"
  CALL Construct_Operator_1D(Operator=Cav1N,        operator_type="Nb_photons",  Mode=Cavity_mode_1, Debug=.FALSE.)
  FLUSH(out_unit)

  Nb_M = Molecule_1%Nb
  Nb_C = Cavity_mode_1%Nb
  NB   = Molecule_1%Nb * Cavity_mode_1%Nb  
                   
  !--------------------------------Construction of the Total Hamiltonian matrix with CM-couplings--------------------------------
  WRITE(out_unit,*); WRITE(out_unit,*) "----------------------------------Construction of the Total Hamiltonian matrix with CM-co&
                                       &uplings----------------------------------"

  ALLOCATE(TotH(NB, NB))
  CALL Construct_total_hamiltonian_1p1D_R1(TotH, Cav1Position, Cav1H, Mol1DipMomt, Mol1H, Debug=.FALSE.)

  IF (Debug) THEN
    WRITE(out_unit,*); WRITE(out_unit,*) "Total Hamiltonian 1p1D (lambda = "//TO_string(Cav1H%lambda)//" ; w_C = "//TO_string(Ca&
                                        &v1H%w)//" ; w_M = "//TO_string(Mol1H%w)//")"
    CALL Write_Mat(TotH, out_unit, NB, info="TotH")
  END IF

  WRITE(out_unit,*); WRITE(out_unit,*) "A = "//TO_string(Cav1H%lambda * CteMol1DipMomt / (2*SQRT(Molecule_1%m)))//" (1st order &
                                      &energy correction of the 1st excited level)"
  !WRITE(out_unit,*); WRITE(out_unit,*) "A detu = "//TO_string(Cav1H%lambda * CteMol1DipMomt * SQRT(Cavity_mode_1%w) / (2*SQRT(M&
   !                                  &olecule_1%w * Molecule_1%m)))//" (1st order energy correction of the 1st excited level)"
  WRITE(out_unit,*); WRITE(out_unit,*) "2A = "//TO_string(2 * Cav1H%lambda * CteMol1DipMomt / (2*SQRT(Molecule_1%m)))
  !WRITE(out_unit,*); WRITE(out_unit,*) "2A detu = "//TO_string(2 * Cav1H%lambda * CteMol1DipMomt * SQRT(Cavity_mode_1%w) / (2*SQR&
    !                                 &T(Molecule_1%w * Molecule_1%m)))

    !-------------------------------------------------Computation of Eigenstates-------------------------------------------------
  WRITE(out_unit,*); WRITE(out_unit,*) "-------------------------------------------------Computation of Eigenstates--------------&
                                       &-----------------------------------"
  ALLOCATE(REigval(NB))
  ALLOCATE(REigvec(NB,NB))
  CALL diagonalization(TotH, REigval, REigvec)

  WRITE(out_unit,*); WRITE(out_unit,*) 'EIGENVALUES'
  CALL WRITE_Vec(REigval, out_unit, 10, info = 'VP_TotH[Ha]')

  WRITE(out_unit,*); WRITE(out_unit,*) "Energy gap = "//TO_string(REigval(3) - REigval(2))

  IF (Debug) THEN
    WRITE(out_unit,*); WRITE(out_unit,*) 'EIGENVECTORS'
    CALL WRITE_Mat(REigvec, out_unit, 6, info = 'Eigenvectors')
  END IF 
  
    !--------------------------------------------Computation of transition intensities-------------------------------------------
  WRITE(out_unit,*); WRITE(out_unit,*) "--------------------------------------------Computation of transition intensities------&
                                       &-------------------------------------"
  CALL Initialize_transition_matrix(Intensities, Mol1DipMomt, REigvec, Nb_states=10, Debug=.FALSE.)
  CALL Compute_transition_matrix(Intensities,    Mol1DipMomt, REigvec, Debug=.FALSE.)

  DO I = 1, 3
    WRITE(out_unit,*) "Transition energy GSto"//TO_string(I)//" = "//TO_string((REigval(I+1)-REigval(1)))
  END DO
  
  DEALLOCATE(TotH); DEALLOCATE(REigval); DEALLOCATE(REigvec)

    !-----------------------------Construction of the 1p1D total system Mass-weighted Hessian matrix-----------------------------
  WRITE(out_unit,*); WRITE(out_unit,*) "--------------------------------Construction of the 1p1D total system Mass-weighted Hessi&
                                       &an matrix--------------------------------"
  MWH(1,1) = Cavity_mode_1%w**2
  MWH(2,2) = Molecule_1%w**2
  MWH(1,2) = Cavity_mode_1%lambda*Cavity_mode_1%w*CteMol1DipMomt / SQRT(Molecule_1%m)
  MWH(2,1) = MWH(1,2)

  IF (Debug) THEN
    CALL Write_Mat(MWH, out_unit, Size(MWH, dim=2), info="MWH")
  END IF

  CALL diagonalization(MWH, Normal_modes, Normal_coordinates)
  CALL Write_Vec(Normal_modes, out_unit, Size(Normal_modes), info="Normal modes")
  CALL Write_Mat(Normal_coordinates, out_unit, Size(Normal_coordinates, dim=2), info="Normal coordniates")

  WRITE(out_unit,*)
  DO i = 1, Size(Normal_modes)
    IF (Normal_modes(i) >= 0) THEN
      WRITE(out_unit,*) TO_string(i)//"^{th} Normal coordinate has positive squared frequency, lead&
                       &ing to w_"//TO_string(i)//" = "//TO_string(SQRT(Normal_modes(i)))
    ELSE
      WRITE(out_unit,*) TO_string(i)//"^{th} Normal coordinate has NEGATIVE squared frequency, lead&
      &ing to w_"//TO_string(i)//" = "//TO_string(EYE*SQRT(-Normal_modes(i)))
    END IF
  END DO
  WRITE(out_unit,*) "Expected ZPE by half-sum of the total system eigenpulsations (from MWH): "//TO&
                    &_string( ( SQRT(Normal_modes(1))+SQRT(Normal_modes(2)) )/2 )


END PROGRAM
