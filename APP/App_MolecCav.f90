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
  real(kind=Rkind), allocatable :: Psi_1p1D_R2(:,:)                                                                                 ! the total system (matter-cavity) wavefunction. Size Nb_M*Nb_C. |Psi_1p1D_R2> = |Molecule_WF>.TENSOR.|Cavity_WF> 
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

  !-----------------------------------------------------------Utilities----------------------------------------------------------
  integer                       :: i, Nb_M, Nb_C, NB, N


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

    !-----------------------------------------------------Second cavity mode-----------------------------------------------------
  Cavity_mode_1_uncoupled        = Cavity_mode_1

  Cavity_mode_1_uncoupled%lambda = ZERO

  CALL Construct_Operator_1D(Cav1H_uncloupled,       "Hamiltonian", Mode=Cavity_mode_1_uncoupled, Debug=.FALSE.)

    !-----------------------------------------------------Total Wavefunctions----------------------------------------------------
  WRITE(out_unit,*); WRITE(out_unit,*) "  -------------------------------------------------------Total Wavefunctions-------------&
                                       &-----------------------------------------"
  ALLOCATE(Psi_1p1D_R2(Molecule_1%Nb,Cavity_mode_1%Nb))
  Psi_1p1D_R2(:,:) = ZERO
  DO i = 1, MIN(Molecule_1%Nb, Cavity_mode_1%Nb)                                                                                 ! initialize Systel_WF arbitrarily
    Psi_1p1D_R2(i,i) = i
  END DO

  IF (Debug .AND. Nb_C <= 10) THEN
    WRITE(out_unit,*) "Not normalized Psi_1p1D_R2 total wavefunction"
    CALL Write_Mat(Psi_1p1D_R2, out_unit, Cavity_mode_1%Nb, info="NN Psi_1p1D_R2")
  ELSE IF (Debug .AND. Nb_C > 10) THEN
    WRITE(out_unit,*); WRITE(out_unit,*) "Not normalized Psi_1p1D_R2 total wavefunction (10:10 slicing)"
    CALL Write_Mat(Psi_1p1D_R2(1:MIN(Nb_M,10),1:10), out_unit, 10, info="NN Psi_1p1D_R2 (sliced)")
  END IF 

  ALLOCATE(Psi_1p1D_R1(Molecule_1%Nb * Cavity_mode_1%Nb))
  CALL Mapping_WF_2DTO1D(Psi_1p1D_R1, Psi_1p1D_R2)

  IF (Debug .AND. NB <= 10) THEN
    WRITE(out_unit,*); WRITE(out_unit,*) "Not normalized Psi_1p1D_R1"
    CALL Write_Vec(Psi_1p1D_R1, out_unit, 1, info="NN Psi_1p1D_R1")
  ELSE IF (Debug .AND. NB > 10) THEN
    WRITE(out_unit,*); WRITE(out_unit,*) "Not normalized Psi_1p1D_R1 total wavefunction (10:10 slicing)"
    CALL Write_Vec(Psi_1p1D_R1(1:10), out_unit, 1, info="NN Psi_1p1D_R1 (sliced)")
  END IF 

  CALL Normalize(Psi_1p1D_R2)
  IF (Nb_C <= 10) THEN
    WRITE(out_unit,*); WRITE(out_unit,*) "Normalized Psi_1p1D_R2 wavefunction"
    CALL Write_Mat(Psi_1p1D_R2, out_unit, Cavity_mode_1%Nb, info="Psi_1p1D_R2")
  ELSE 
    WRITE(out_unit,*); WRITE(out_unit,*) "Normalized Psi_1p1D_R2 wavefunction (10:10 slicing)"
    CALL Write_Mat(Psi_1p1D_R2(1:MIN(Nb_M,10),1:10), out_unit, Cavity_mode_1%Nb, info="Psi_1p1D_R2 (sliced)")
  END IF 

  CALL Normalize(Psi_1p1D_R1)
  IF (NB <= 10) THEN
    WRITE(out_unit,*); WRITE(out_unit,*) "Normalized Psi_1p1D_R1"
    CALL Write_Vec(Psi_1p1D_R1, out_unit, 1, info="Psi_1p1D_R1")
  ELSE
    WRITE(out_unit,*); WRITE(out_unit,*) "Normalized Psi_1p1D_R1 (10:10 slicing)"
    CALL Write_Vec(Psi_1p1D_R1(1:10), out_unit, 1, info="Psi_1p1D_R1 (sliced)")
  END IF

  WRITE(out_unit,*) "--------------------------------------------------FINISHED SYSTEM INITIALIZATION----------------------------&
                    &----------------------"


  !----------------------------------------Computation of the photon number of Psi_1p1D_R2---------------------------------------
  WRITE(out_unit,*); WRITE(out_unit,*) "-------------------------------------------Computation of the photon number of Psi_1p1D_R&
                                       &2-------------------------------------------"

  ALLOCATE(Result_psi_1p1D_R2(Molecule_1%Nb, Cavity_mode_1%Nb))

  CALL MolecCav_Action_operator_2D(Result_psi_1p1D_R2, Cav1N, Psi_1p1D_R2)                                                             ! /!\ do not trust so much this subroutine elsewhere than here (has not been refined yet)
  IF (Debug .AND. Nb_C <= 10) THEN
    WRITE(out_unit,*) "Action of Nb photon operator over Psi_1p1D_R2 (for the first method) :"
    CALL Write_Mat(Result_psi_1p1D_R2, out_unit, Cavity_mode_1%Nb, info="N_psi_1p1D")
  ELSE IF (Debug .AND. Nb_C > 10) THEN
    WRITE(out_unit,*) "Action of Nb photon operator over Psi_1p1D_R2 (for the first method) (1:10,1:10 slicing) :"
    CALL Write_Mat(Result_psi_1p1D_R2(1:MIN(Nb_M,10),1:10), out_unit, Cavity_mode_1%Nb, info="N_psi_1p1D (sliced)")
  END IF 

  CALL Scalar_product(Average, Result_psi_1p1D_R2, Psi_1p1D_R2)
  WRITE(out_unit,*) "The average nb of photons of the normalised Psi_1p1D_R2 is (first method) : ", Average

  CALL MolecCav_Average_value_operator_2D(Average, Cav1N, Psi_1p1D_R2)
  WRITE(out_unit,*) "The average nb of photons of the normalised Psi_1p1D_R2 is (second method) : ", Average

  DEALLOCATE(Result_psi_1p1D_R2)

  !------------------------------------------------An experiment on the Nb_photons-----------------------------------------------
  WRITE(out_unit,*); WRITE(out_unit,*) "--------------------------------------------------An experiment on the Nb_photons--------&
                                       &-----------------------------------------"
  IF (Nb_C == 1) THEN
    WRITE(out_unit,*) "########################## WARNING ########################## WARNING ########################## WARNING #&
                      &########################"
    WRITE(out_unit,*) "                          This test is NOT likely to return the expected outcome when Nb_C = 1 is chosen "
    WRITE(out_unit,*) "########################## WARNING ########################## WARNING ########################## WARNING #&
                      &########################"
  END IF 

  ALLOCATE(CavPsi(Cavity_mode_1%Nb))
  CavPsi    = ZERO
  CavPsi(1) = ONE
  IF (Nb_C > 1) CavPsi(2) = ONE                                                                                                  ! \ket{Psi} = \ket{0} + \ket{1}
  IF (Debug .AND. Size(CavPsi) <= 10) THEN
    WRITE(out_unit,*) 'Not normalized cavity mode WF (1D)'
    CALL Write_Vec(CavPsi, out_unit, 1, info="NN CavPsi")
  ELSE IF (Debug .AND. Size(CavPsi) > 10) THEN
    WRITE(out_unit,*) 'Not normalized cavity mode WF (1D) (1:10 slicing)'
    CALL Write_Vec(CavPsi(1:10), out_unit, 1, info="NN CavPsi (sliced)")
  END IF

  CALL Normalize(CavPsi)
  IF (Debug .AND. Size(CavPsi) <= 10) THEN
    WRITE(out_unit,*) 'Normalized cavity mode WF (1D) (N.B. \frac{1}{\sqrt{2}} = 0.7071067811865475)'
    CALL Write_Vec(CavPsi, out_unit, 1, info="CavPsi")
  ELSE IF (Debug .AND. Size(CavPsi) > 10) THEN
    WRITE(out_unit,*) 'Normalized cavity mode WF (1D) (1:10 slicing) (N.B. \frac{1}{\sqrt{2}} = 0.7071067811865475)'
    CALL Write_Vec(CavPsi(1:10), out_unit, 1, info="CavPsi (sliced)")
  END IF

  ALLOCATE(Result_psi_1p1D_R1(Cavity_mode_1%Nb))
  CALL Action_Operator_1D(Result_psi_1p1D_R1, Cav1N, CavPsi)
  IF (Debug .AND. Size(CavPsi) <= 10) THEN
    WRITE(out_unit,*) 'Action of Nb_photons over the cavity mode WF (1D)'
    CALL Write_Vec(Result_psi_1p1D_R1, out_unit, 1, info="N_CavPsi")
  ELSE IF (Debug .AND. Size(CavPsi) > 10) THEN
    WRITE(out_unit,*) 'Action of Nb_photons over the cavity mode WF (1D) (1:10 slicing)'
    CALL Write_Vec(Result_psi_1p1D_R1(1:10), out_unit, 1, info="N_CavPsi (sliced)")
  END IF

  CALL Average_value_operator_1D(Average, Cav1N, CavPsi)
  WRITE(out_unit,*) 'The avegeraged number of photons of that CavPsi is analytically expected to be &
                    & 0.5 and is = ', Average
  
  DEALLOCATE(Result_psi_1p1D_R1); DEALLOCATE(CavPsi)

  !--------------------------------Construction of a Total Hamiltonian matrix without CM-couplings-------------------------------
  WRITE(out_unit,*); WRITE(out_unit,*) "----------------------------------Construction of a Total Hamiltonian matrix without CM-c&
                                       &ouplings---------------------------------"
  ALLOCATE(TotH_uncoupled(NB, NB))
  
  CALL Construct_total_hamiltonian_1p1D(TotH_uncoupled, Cav1Position, Cav1H_uncloupled, Mol1DipMomt, Mol1H, Debug=.FALSE.)

  IF (Debug .AND. NB <= 10) THEN
    WRITE(out_unit,*); WRITE(out_unit,*) "Total Hamiltonian 1p1D (lambda = 0, w_C /= w_M)"
    CALL Write_Mat(TotH_uncoupled, out_unit, NB, info="TotH")
  ELSE IF (Debug .AND. NB > 10) THEN
    WRITE(out_unit,*); WRITE(out_unit,*) "Total Hamiltonian 1p1D (lambda = 0, w_C /= w_M) (10:10 slicing)"
    CALL Write_Mat(TotH_uncoupled(1:10,1:10), out_unit, 10, info="TotH (sliced)")
  END IF

    !----------------------------------------------------Computation Eigenstates---------------------------------------------------
  ALLOCATE(REigval_uncoupled(NB))
  ALLOCATE(REigvec_uncoupled(NB,NB))
  CALL diagonalization(TotH_uncoupled, REigval_uncoupled, REigvec_uncoupled)

  WRITE(out_unit,*); WRITE(out_unit,*) 'EIGENVALUES'
  IF (NB <= 10) CALL WRITE_Vec(REigval_uncoupled,       out_unit, 10, info = 'VP_TotH[Ha]')
  IF (NB > 10)  CALL WRITE_Vec(REigval_uncoupled(1:10), out_unit, 10, info = 'Ten_first_VP_TotH[Ha]')

  IF (Debug .AND. NB <= 10) THEN
    WRITE(out_unit,*); WRITE(out_unit,*) 'EIGENVECTORS'
    CALL WRITE_Mat(REigvec_uncoupled, out_unit, 6, info = 'Eigenvectors')
  ELSE IF (Debug .AND. NB > 10) THEN
    WRITE(out_unit,*); WRITE(out_unit,*) 'EIGENVECTORS'
    CALL WRITE_Mat(REigvec_uncoupled(1:10,1:10), out_unit, 6, info = 'Ten first Eigenvectors (1:10 slicing)')
  END IF 

  DEALLOCATE(TotH_uncoupled); DEALLOCATE(REigval_uncoupled); DEALLOCATE(REigvec_uncoupled)

    !---------------------------Construction of the 1p1D uncoupled system Mass-weighted Hessian matrix---------------------------
  WRITE(out_unit,*); WRITE(out_unit,*) "------------------------------Construction of the 1p1D uncoupled system Mass-weighted Hes&
                                       &sian matrix------------------------------"
  MWH_uncoupled(1,1) = Cavity_mode_1_uncoupled%w**2
  MWH_uncoupled(2,2) = Molecule_1%w**2
  MWH_uncoupled(1,2) = Cavity_mode_1_uncoupled%lambda*Cavity_mode_1%w*CteMol1DipMomt / SQRT(Molecule_1%m)
  MWH_uncoupled(2,1) = MWH_uncoupled(1,2)

  IF (Debug) THEN
    CALL Write_Mat(MWH_uncoupled, out_unit, Size(MWH, dim=2), info="MWH")
  END IF

  CALL diagonalization(MWH_uncoupled, Normal_modes_uncoupled, Normal_coordinates_uncoupled)
  CALL Write_Vec(Normal_modes_uncoupled,       out_unit, Size(Normal_modes_uncoupled),              info="Normal modes")
  CALL Write_Mat(Normal_coordinates_uncoupled, out_unit, Size(Normal_coordinates_uncoupled, dim=2), info="Normal coordniates")

  WRITE(out_unit,*)
  DO i = 1, Size(Normal_modes_uncoupled)
    IF (Normal_modes_uncoupled(i) >= 0) THEN
      WRITE(out_unit,*) TO_string(i)//"^{th} Normal coordinate has positive squared frequency, lead&
                       &ing to w_"//TO_string(i)//" = "//TO_string(SQRT(Normal_modes_uncoupled(i)))
    ELSE
      WRITE(out_unit,*) TO_string(i)//"^{th} Normal coordinate has NEGATIVE squared frequency, lead&
      &ing to w_"//TO_string(i)//" = "//TO_string(EYE*SQRT(-Normal_modes_uncoupled(i)))
    END IF
  END DO
  WRITE(out_unit,*) "Expected ZPE by half-sum of the total system eigenpulsations : "//TO_string( (&
                   & SQRT(Normal_modes_uncoupled(1))+SQRT(Normal_modes_uncoupled(2)) )/2 )

                   
  !--------------------------------Construction of the Total Hamiltonian matrix with CM-couplings--------------------------------
  WRITE(out_unit,*); WRITE(out_unit,*) "----------------------------------Construction of the Total Hamiltonian matrix with CM-co&
                                       &uplings----------------------------------"
  CALL time_perso("Beginning of time")

  ALLOCATE(TotH(NB, NB))
  CALL Construct_total_hamiltonian_1p1D(TotH, Cav1Position, Cav1H, Mol1DipMomt, Mol1H, Debug=.FALSE.)
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

  ALLOCATE(TotH_bis(NB, NB))
  CALL Construct_total_hamiltonian_1p1D_R1(TotH_bis, Cav1Position, Cav1H, Mol1DipMomt, Mol1H, Debug=.FALSE.)
    !-----------------------------------------------Computation of some observables----------------------------------------------
!  CALL Average_value_TotH(Average, TotH, Psi_1p1D_R1)
!  WRITE(out_unit,*) "Average E_tot = ", Average, "Ha"

    !-------------------------------------------------Computation of Eigenstates-------------------------------------------------
  ALLOCATE(REigval(NB))
  ALLOCATE(REigvec(NB,NB))
  CALL time_perso("Beginning diagonalization")
  CALL diagonalization(TotH, REigval, REigvec)
  CALL time_perso("TotH diagonalized")

  WRITE(out_unit,*); WRITE(out_unit,*) 'EIGENVALUES'
  IF (NB <= 20) CALL WRITE_Vec(REigval, out_unit, 10, info = 'VP_TotH[Ha]')
  IF (NB > 20)  CALL WRITE_Vec(REigval(1:20), out_unit, 20, info = 'Twenty_first_VP_TotH[Ha]')

  IF (Debug .AND. NB <= 20) THEN
    WRITE(out_unit,*); WRITE(out_unit,*) 'EIGENVECTORS'
    CALL WRITE_Mat(REigvec, out_unit, 6, info = 'Eigenvectors')
  ELSE IF (Debug .AND. NB > 10) THEN
    WRITE(out_unit,*); WRITE(out_unit,*) 'EIGENVECTORS'
    CALL WRITE_Mat(REigvec(1:20,1:20), out_unit, 6, info = 'Twenty first Eigenvectors (1:10 slicing)')
  END IF 

  ALLOCATE(REigval_bis(NB))
  ALLOCATE(REigvec_bis(NB,NB))
  CALL diagonalization(TotH_bis, REigval_bis, REigvec_bis)
  WRITE(out_unit,*); WRITE(out_unit,*) 'EIGENVALUES'
  IF (NB <= 20) CALL WRITE_Vec(REigval_bis, out_unit, 10, info = 'VP_TotH[Ha]')
  IF (NB > 20)  CALL WRITE_Vec(REigval_bis(1:20), out_unit, 20, info = 'Twenty_first_VP_TotH[Ha]')

    !--------------------------------------------Computation of transition intensities-------------------------------------------
  CALL time_perso("Beginning to compute the transition intensities matrix")
  CALL Initialize_transition_matrix(Intensities, Mol1DipMomt, REigvec, Nb_states=10, Debug=.TRUE.)
  CALL Compute_transition_matrix(Intensities,    Mol1DipMomt, REigvec, Debug=.TRUE.)
  CALL time_perso("Transition matrix computed")

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

                    
  !-------------------------------------------An experiment on the Psi_1p1D_R2 analysis------------------------------------------
  ALLOCATE(Mol1Weights(Nb_M))
  ALLOCATE(Cav1Weights(Nb_C))

  CALL Reduced_density_psi_R(Mol1Weights, Cav1Weights, Psi_1p1D_R2, Debug=.TRUE.)
  WRITE(out_unit,*); CALL Write_Vec(Mol1Weights, out_unit, 1, info="Mol1Weights")
  WRITE(out_unit,*); CALL Write_Vec(Cav1Weights, out_unit, 1, info="Cav1Weights")
  CALL Reduced_density_psi_R(Mol1Weights, Cav1Weights, Psi_1p1D_R1, Debug=.TRUE.)
  WRITE(out_unit,*); CALL Write_Vec(Mol1Weights, out_unit, 1, info="Mol1Weights")
  WRITE(out_unit,*); CALL Write_Vec(Cav1Weights, out_unit, 1, info="Cav1Weights")
  DEALLOCATE(Mol1Weights); DEALLOCATE(Cav1Weights)

  !-----------------------------------------------An experiment on the TotH action-----------------------------------------------
  ALLOCATE(Result_psi_1p1D_R1(NB))
  ALLOCATE(Result_psi_1p1D_R2(Nb_M, Nb_C))
!  CALL Write_Mat(Psi_1p1D_R2, out_unit, 10, info="Psi_R2")
!  CALL Write_Vec(Psi_1p1D_R1, out_unit, 1,  info="Psi_R1")
  CALL Action_total_hamiltonian_1p1D(Result_psi_1p1D_R2, Cav1Position, Cav1H, Mol1DipMomt, Mol1H, Psi_1p1D_R2, Debug=.FALSE.)
  CALL Action_total_hamiltonian_1p1D(Result_psi_1p1D_R1, Cav1Position, Cav1H, Mol1DipMomt, Mol1H, Psi_1p1D_R1, Debug=.FALSE.)
!  CALL Write_Mat(Result_psi_1p1D_R2, out_unit, 10, info="Res_R2")
!  CALL Write_Vec(Result_psi_1p1D_R1, out_unit, 1,  info="Res_R1")
  DEALLOCATE(Result_psi_1p1D_R2); DEALLOCATE(Result_psi_1p1D_R1)


END PROGRAM