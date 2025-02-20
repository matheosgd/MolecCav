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
PROGRAM App_MolecCav
  USE QDUtil_m
  USE Algebra_m
  USE Cavity_mode_m
  USE Operator_1D_m
  USE Operator_2D_m
  USE Total_hamiltonian_m
  IMPLICIT NONE


  !-------------Diatomic molecule in a harmonic electonic potential------------
  TYPE(Cavity_mode_t)           :: Molecule_1
  TYPE(Operator_1D_t)           :: Mol1H                                       ! matrix of the one-dimensional harmonic Hamiltonian associated with HO D
  TYPE(Operator_1D_t)           :: Mol1Position
  TYPE(Operator_1D_t)           :: Mol1N                                       ! here N does not count the number of photons but of excitation quanta in the vibrational state
  TYPE(Operator_1D_t)           :: Mol1_dipolar_moment
  real(kind=Rkind)              :: Cte_dipole_moment = ONE                     ! the intensity of the variation of the dipole moment with a variation of the matter DOF
  
  !------------------------------First cavity mode-----------------------------
  TYPE(Cavity_mode_t)           :: Cavity_mode_1
  TYPE(Operator_1D_t)           :: Cav1H                                       ! matrix of the one-dimensional harmonic Hamiltonian associated with HO D
  TYPE(Operator_1D_t)           :: Cav1Position
  TYPE(Operator_1D_t)           :: Cav1N

  !--------------------------------Wavefunctions-------------------------------
  real(kind=Rkind), allocatable :: Psi_1p1D(:,:)                               ! The total system (matter-cavity) wavefunction. Size Nb_M*Nb_C. |Psi_1p1D> = |Molecule_WF>.TENSOR.|Cavity_WF> 
  real(kind=Rkind), allocatable :: Psi_1D_mapped(:)
  real(kind=Rkind), allocatable :: CavPsi(:)

  !------------------------------Total Hamiltonian-----------------------------
  real(kind=Rkind), allocatable :: TotH(:,:)                               

  !-----------------------------------Results----------------------------------
  real(kind=Rkind)              :: Average
  real(kind=Rkind), allocatable :: Result_psi_1D(:)
  real(kind=Rkind), allocatable :: Result_psi_1p1D(:,:)
  real(kind=Rkind), allocatable :: REigval(:)
  real(kind=Rkind), allocatable :: REigvec(:,:)

  !----------------------------------Utilities---------------------------------
  integer                       :: i, NB
  logical, parameter            :: Debug = .TRUE.


!-----------------------------SYSTEM INITIALIZATION----------------------------
  !-------------Diatomic molecule in a harmonic electonic potential------------
  WRITE(out_unit,*) "-----------------------------SYSTEM INITIALIZATION----------------------------"
  WRITE(out_unit,*) "  -------------Diatomic molecule in a harmonic electonic potential------------"
  CALL MolecCav_Read_cavity_mode(Mode=Molecule_1, nio=in_unit)

  WRITE(out_unit,*) "Molecular Hamiltonian :"
  CALL Construct_Operator_1D(Operator=Mol1H, &
                           & operator_type="Hamiltonian", &
                           & Mode=Molecule_1, &
                           & Debug=.TRUE.)
  
  WRITE(out_unit,*) "Molecular Position :"
  CALL Construct_Operator_1D(Operator=Mol1Position, &
                           & operator_type="Position", &
!                          & Dense=.TRUE., &
                           & Mode=Molecule_1, &
                           & Debug=Debug)

  WRITE(out_unit,*) "Molecular Nb_photon (vibrationnal excitation quanta)"
  CALL Construct_Operator_1D(Operator=Mol1N, &
                           & operator_type="Nb_photons", &
                           & Mode=Molecule_1, &
                           & Debug=Debug)
  
  WRITE(out_unit,*) "Molecular Dipole moment :"
  CALL Construct_Operator_1D(Operator=Mol1_dipolar_moment, &
                           & operator_type="Position", &                       ! initialized as a position operator because of approximation over its expression (cf. readme.md or manual)
  !                        & Dense=.TRUE., &                                   !/!\ if initialize as dense MUST change two lines below : Mol1_dipolar_moment%Dense_val_R /!\
                           & Mode=Molecule_1, &
                           & Debug=.FALSE.)

  Mol1_dipolar_moment%Band_val_R = Mol1_dipolar_moment%Band_val_R*Cte_dipole_moment ! /!\ so that the matrix already contains the intensity constant of the dipolar moment with the position of the matter (cf. manual for formulas)
    
  IF (Debug) CALL Write_Mat(Mol1_dipolar_moment%Band_val_R, out_unit, 3, info="Mol1_dipolar_moment")
  FLUSH(out_unit)

    !-----------------------------First cavity mode----------------------------
  WRITE(out_unit,*); WRITE(out_unit,*) "  -----------------------------First cavity mode----------------------------"
  CALL MolecCav_Read_cavity_mode(Mode=Cavity_mode_1, nio=in_unit)

  WRITE(out_unit,*) "Cavity mode Hamiltonian :"
  CALL Construct_Operator_1D(Operator=Cav1H, &
                           & operator_type="Hamiltonian", &
                           & Mode=Cavity_mode_1, &
                           & Debug=.TRUE.)

  WRITE(out_unit,*) "Cavity mode Position :"
  CALL Construct_Operator_1D(Operator=Cav1Position, &
                           & operator_type="Position", &
!                           & Dense=.TRUE., &
                           & Mode=Cavity_mode_1, &
                           & Debug=Debug)

  WRITE(out_unit,*) "Cavity mode Nb_photons :"
  CALL Construct_Operator_1D(Operator=Cav1N, &
                           & operator_type="Nb_photons", &
                           & Mode=Cavity_mode_1, &
                           & Debug=Debug)
  FLUSH(out_unit)

  NB = Molecule_1%Nb * Cavity_mode_1%Nb

    !----------------------------Total Wavefunctions---------------------------
  WRITE(out_unit,*); WRITE(out_unit,*) "  ----------------------------Total Wavefunctions---------------------------"
  ALLOCATE(Psi_1p1D(Molecule_1%Nb,Cavity_mode_1%Nb))
  Psi_1p1D(:,:) = ZERO
  DO i = 1, MIN(Molecule_1%Nb, Cavity_mode_1%Nb)                               ! initialize Systel_WF arbitrarily
    Psi_1p1D(i,i) = i
  END DO

  IF (Debug .AND. NB <= 10) THEN
    WRITE(out_unit,*) "Not normalized Psi_1p1D total wavefunction"
    CALL Write_Mat(Psi_1p1D, out_unit, Cavity_mode_1%Nb, info="NN Psi_1p1D")
  ELSE IF (Debug .AND. NB > 10) THEN
    WRITE(out_unit,*); WRITE(out_unit,*) "Not normalized Psi_1p1D total wavefunction (10:10 slicing)"
    CALL Write_Mat(Psi_1p1D(1:10,1:10), out_unit, 10, info="NN Psi_1p1D (sliced)")
  END IF 

  ALLOCATE(Psi_1D_mapped(Molecule_1%Nb * Cavity_mode_1%Nb))
  CALL Mapping_WF_2DTO1D(Psi_1D_mapped, Psi_1p1D)

  IF (Debug .AND. NB <= 10) THEN
    WRITE(out_unit,*); WRITE(out_unit,*) "Not normalized Psi_1D_mapped"
    CALL Write_Vec(Psi_1D_mapped, out_unit, 1, info="NN Psi_1D_mapped")
  ELSE IF (Debug .AND. NB > 10) THEN
    WRITE(out_unit,*); WRITE(out_unit,*) "Not normalized Psi_1D_mapped total wavefunction (10:10 slicing)"
    CALL Write_Vec(Psi_1D_mapped(1:10), out_unit, 1, info="NN Psi_1D_mapped (sliced)")
  END IF 

  CALL Normalize(Psi_1p1D)
  IF (NB < 10) THEN
    WRITE(out_unit,*); WRITE(out_unit,*) "Normalized Psi_1p1D wavefunction"
    CALL Write_Mat(Psi_1p1D, out_unit, Cavity_mode_1%Nb, info="Psi_1p1D")
  ELSE 
    WRITE(out_unit,*); WRITE(out_unit,*) "Normalized Psi_1p1D wavefunction (10:10 slicing)"
    CALL Write_Mat(Psi_1p1D(1:10,1:10), out_unit, Cavity_mode_1%Nb, info="Psi_1p1D (sliced)")
  END IF 

  CALL Normalize(Psi_1D_mapped)
  IF (NB < 10) THEN
    WRITE(out_unit,*); WRITE(out_unit,*) "Normalized Psi_1D_mapped"
    CALL Write_Vec(Psi_1D_mapped, out_unit, 1, info="Psi_1D_mapped")
  ELSE
    WRITE(out_unit,*); WRITE(out_unit,*) "Normalized Psi_1D_mapped (10:10 slicing)"
    CALL Write_Vec(Psi_1D_mapped(1:10), out_unit, 1, info="Psi_1D_mapped (sliced)")
  END IF

  WRITE(out_unit,*) "------------------------FINISHED SYSTEM INITIALIZATION------------------------"


  !----------------Computation of the photon number of Psi_1p1D----------------
  WRITE(out_unit,*); WRITE(out_unit,*) "----------------Computation of the photon number of Psi_1p1D---------------"

  ALLOCATE(Result_psi_1p1D(Molecule_1%Nb, Cavity_mode_1%Nb))

  CALL MolecCav_Action_operator_2D(Result_psi_1p1D, Cav1N, Psi_1p1D)           ! /!\ do not trust so much this subroutine elsewhere than here (has not been refined yet)
  IF (Debug .AND. NB <= 10) THEN
    WRITE(out_unit,*) "Action of Nb photon operator over Psi_1p1D (for the first method) :"
    CALL Write_Mat(Result_psi_1p1D, out_unit, Cavity_mode_1%Nb, info="N_psi_1p1D")
  ELSE IF (Debug .AND. NB > 10) THEN
    WRITE(out_unit,*) "Action of Nb photon operator over Psi_1p1D (for the first method) (1:10,1:10 slicing) :"
    CALL Write_Mat(Result_psi_1p1D(1:10,1:10), out_unit, Cavity_mode_1%Nb, info="N_psi_1p1D (sliced)")
  END IF 

  CALL Scalar_product(Average, Result_psi_1p1D, Psi_1p1D)
  WRITE(out_unit,*) "The average nb of photons of the normalised Psi_1p1D is (first method) : ", Average

  CALL MolecCav_Average_value_operator_2D(Average, Cav1N, Psi_1p1D)
  WRITE(out_unit,*) "The average nb of photons of the normalised Psi_1p1D is (second method) : ", Average

  DEALLOCATE(Result_psi_1p1D)

  !-----------------------An experiment on the Nb_photons----------------------
  WRITE(out_unit,*); WRITE(out_unit,*) "-----------------------An experiment on the Nb_photons----------------------"

  ALLOCATE(CavPsi(Cavity_mode_1%Nb))
  CavPsi    = ZERO
  CavPsi(1) = ONE
  CavPsi(2) = ONE                                                              ! \ket{Psi} = \ket{0} + \ket{1}
  IF (Debug) WRITE(out_unit,*) 'Not normalized cavity mode WF (1D)'
  IF (Debug) CALL Write_Vec(CavPsi, out_unit, 1, info="NN CavPsi")
  CALL Normalize(CavPsi)
  WRITE(out_unit,*) 'Normalized cavity mode WF (1D) (N.B. \frac{1}{\sqrt{2}} = 0.7071067811865475)'
  CALL Write_Vec(CavPsi, out_unit, 1, info="CavPsi")

  ALLOCATE(Result_psi_1D(Cavity_mode_1%Nb))
  CALL Action_Operator_1D(Result_psi_1D, Cav1N, CavPsi)
  IF (Debug .AND. Size(CavPsi) <= 10) THEN
    WRITE(out_unit,*) 'Action of Nb_photons over the cavity mode WF (1D)'
    CALL Write_Vec(Result_psi_1D, out_unit, 1, info="N_CavPsi")
  ELSE IF (Debug .AND. Size(CavPsi) > 10) THEN
    WRITE(out_unit,*) 'Action of Nb_photons over the cavity mode WF (1D) (1:10 slicing)'
    CALL Write_Vec(Result_psi_1D(1:10), out_unit, 1, info="N_CavPsi (sliced)")
  END IF

  CALL Average_value_operator_1D(Average, Cav1N, CavPsi)
  WRITE(out_unit,*) 'The avegeraged number of photons of that CavPsi is analytically expected to be &
                    & 0.5 and is = ', Average
  
  DEALLOCATE(Result_psi_1D); DEALLOCATE(CavPsi)

  !----------------An experiment on the total Hamiltonian action---------------
  WRITE(out_unit,*); WRITE(out_unit,*) "----------------An experiment on the total Hamiltonian action---------------"

  IF (Debug .AND. NB <= 10) THEN
    WRITE(out_unit,*); WRITE(out_unit,*) "Psi_1p1D"
    CALL Write_Mat(Psi_1p1D, out_unit, Cavity_mode_1%Nb)
  ELSE IF (Debug .AND. NB > 10) THEN
    WRITE(out_unit,*); WRITE(out_unit,*) "Psi_1p1D (1:10,1:10 slicing)"
    CALL Write_Mat(Psi_1p1D(1:10,1:10), out_unit, Cavity_mode_1%Nb, info="Psi_1p1D (sliced)")
  END IF
  FLUSH(out_unit) 

  ALLOCATE(Result_psi_1p1D(Molecule_1%Nb, Cavity_mode_1%Nb))
  Result_psi_1p1D(:,:) = ZERO

  Cavity_mode_1%lambda = ZERO; Cav1H%lambda = ZERO                                                 ! /!\ the action tot H procedure uses Cav1H%lambda as coupling strenght parameter
  CALL Action_total_hamiltonian_1p1D(Result_psi_1p1D, Cav1Position, Cav1H, Mol1_dipolar_moment, Mol1H, Psi_1p1D, Debug=debug)
  Cavity_mode_1%lambda = ONE;  Cav1H%lambda = ONE

  IF (NB <= 10) THEN
    WRITE(out_unit,*); WRITE(out_unit,*) "Action of the uncoupled total Hamiltonian over Psi_1p1D"
    CALL Write_Mat(Result_psi_1p1D, out_unit, Cavity_mode_1%Nb, info="TotH_psi_1p1D")
  ELSE IF (NB > 10) THEN
    WRITE(out_unit,*); WRITE(out_unit,*) "Action of the uncoupled total Hamiltonian over Psi_1p1D (1:10,1:10 slicing)"
    CALL Write_Mat(Result_psi_1p1D(1:10,1:10), out_unit, Cavity_mode_1%Nb, info="TotH_psi_1p1D (sliced)")
  END IF
  FLUSH(out_unit) 

  DEALLOCATE(Result_psi_1p1D)

  !-------Construction of a Total Hamiltonian matrix without CM-couplings------
  WRITE(out_unit,*); WRITE(out_unit,*) "-------Construction of a Total Hamiltonian matrix without CM-couplings------"
  ALLOCATE(TotH(NB, NB))
  
  Cavity_mode_1%lambda = ZERO; Cav1H%lambda = ZERO                                                 ! /!\ the action tot H procedure uses Cav1H%lambda as coupling strenght parameter
  CALL Construct_total_hamiltonian_1p1D(TotH, Cav1Position, Cav1H, Mol1_dipolar_moment, Mol1H, Debug=.FALSE.)
  Cavity_mode_1%lambda = ONE;  Cav1H%lambda = ONE

  IF (Debug .AND. NB <= 10) THEN
    WRITE(out_unit,*); WRITE(out_unit,*) "Total Hamiltonian 1p1D (lambda = 0, w_C /= w_M)"
    CALL Write_Mat(TotH, out_unit, NB, info="TotH")
  ELSE 
    WRITE(out_unit,*); WRITE(out_unit,*) "Total Hamiltonian 1p1D (lambda = 0, w_C /= w_M) (10:10 slicing)"
    CALL Write_Mat(TotH(1:10,1:10), out_unit, 10, info="TotH (sliced)")
  END IF

    !---------------------------Computation Eigenstates--------------------------
  ALLOCATE(REigval(NB))
  ALLOCATE(REigvec(NB,NB))
  CALL diagonalization(TotH, REigval, Reigvec)

  WRITE(out_unit,*); WRITE(out_unit,*) 'EIGENVALUES'
  IF (NB <= 10) CALL WRITE_Vec(Reigval, out_unit, 10, info = 'VP[Ha]')
  IF (NB > 10)  CALL WRITE_Vec(Reigval(1:10), out_unit, 10, info = 'Ten first VP[Ha]')

  IF (Debug .AND. NB <= 10) THEN
    WRITE(out_unit,*); WRITE(out_unit,*) 'EIGENVECTORS'
    CALL WRITE_Mat(Reigvec, out_unit, 6, info = 'Eigenvectors')
  ELSE IF (Debug .AND. NB > 10) THEN
    WRITE(out_unit,*); WRITE(out_unit,*) 'EIGENVECTORS'
    CALL WRITE_Mat(Reigvec(1:10,1:10), out_unit, 6, info = 'Ten first Eigenvectors (1:10 slicing)')
  END IF 

  DEALLOCATE(TotH); DEALLOCATE(REigval); DEALLOCATE(REigvec)

  !-------Construction of the Total Hamiltonian matrix with CM-couplings-------
  WRITE(out_unit,*); WRITE(out_unit,*) "-------Construction of the Total Hamiltonian matrix with CM-couplings-------"

  ALLOCATE(TotH(NB, NB))
  CALL Construct_total_hamiltonian_1p1D(TotH, Cav1Position, Cav1H, Mol1_dipolar_moment, Mol1H, Debug=.FALSE.)

  IF (Debug .AND. NB <= 10) THEN
    WRITE(out_unit,*); WRITE(out_unit,*) "Total Hamiltonian 1p1D (lambda /= 0, w_C /= w_M)"
    CALL Write_Mat(TotH, out_unit, NB, info="TotH")
  ELSE IF (Debug .AND. NB > 10) THEN
    WRITE(out_unit,*); WRITE(out_unit,*) "Total Hamiltonian 1p1D (lambda /= 0, w_C /= w_M) (10:10 slicing)"
    CALL Write_Mat(TotH(1:10,1:10), out_unit, 10, info="TotH (sliced)")
  END IF

  CALL Write_Mat(TotH, out_unit, Size(TotH), info="TotH")
    !----------------------Computation of some observables---------------------

!  CALL Average_value_TotH(Average, TotH, Psi_1D_mapped)
!  WRITE(out_unit,*) "Average E_tot = ", Average, "Ha"

    !------------------------Computation of Eigenstates------------------------
  ALLOCATE(REigval(NB))
  ALLOCATE(REigvec(NB,NB))
  CALL diagonalization(TotH, REigval, Reigvec)

  WRITE(out_unit,*); WRITE(out_unit,*) 'EIGENVALUES'
  IF (NB <= 10) CALL WRITE_Vec(Reigval, out_unit, 10, info = 'VP[Ha]')
  IF (NB > 10)  CALL WRITE_Vec(Reigval(1:10), out_unit, 10, info = 'Ten first VP[Ha]')

  IF (Debug .AND. NB <= 10) THEN
    WRITE(out_unit,*); WRITE(out_unit,*) 'EIGENVECTORS'
    CALL WRITE_Mat(Reigvec, out_unit, 6, info = 'Eigenvectors')
  ELSE IF (Debug .AND. NB > 10) THEN
    WRITE(out_unit,*); WRITE(out_unit,*) 'EIGENVECTORS'
    CALL WRITE_Mat(Reigvec(1:10,1:10), out_unit, 6, info = 'Ten first Eigenvectors (1:10 slicing)')
  END IF 

  DEALLOCATE(TotH); DEALLOCATE(REigval); DEALLOCATE(REigvec)


END PROGRAM