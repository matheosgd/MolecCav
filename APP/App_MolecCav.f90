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
  TYPE(Operator_1D_t)           :: H_ho_molecule_1                             ! matrix of the one-dimensional harmonic Hamiltonian associated with HO D
  TYPE(Operator_1D_t)           :: x_ho_molecule_1
  TYPE(Operator_1D_t)           :: N_ho_molecule_1                             ! here N does not count the number of photons but of excitation quanta in the vibrational state
  TYPE(Operator_1D_t)           :: Matter_dipolar_moment
  real(kind=Rkind)              :: Cte_dipole_moment                           ! the intensity of the variation of the dipole moment with a variation of the matter DOF
  
  !------------------------------First cavity mode-----------------------------
  TYPE(Cavity_mode_t)           :: Cavity_mode_1
  TYPE(Operator_1D_t)           :: H_ho_cavity_mode_1                          ! matrix of the one-dimensional harmonic Hamiltonian associated with HO D
  TYPE(Operator_1D_t)           :: x_ho_cavity_mode_1
  TYPE(Operator_1D_t)           :: N_ho_cavity_mode_1

  !--------------------------------Wavefunctions-------------------------------
  real(kind=Rkind), allocatable :: System_WF(:,:)                              ! The total system (matter-cavity) wavefunction. Size Nb_M*Nb_C. |System_WF> = |Molecule_WF>.TENSOR.|Cavity_WF> 
  real(kind=Rkind), allocatable :: System_WF_mapped(:)
  real(kind=Rkind), allocatable :: Matter_hamiltonianSystem_WF(:,:)            ! The wavefunction resulting from the action of the matter hamiltonian on System_WF. Size Nb_M*Nb_C. |H_MatterSystem_WF(:,i_C)> = H_Matter|System_WF(:,i_C)>
  real(kind=Rkind), allocatable :: H_tot(:,:)                               

  !------------------------------Some expermients------------------------------

  !-----------------------------------Results----------------------------------
  real(kind=Rkind)              :: Average
  real(kind=Rkind), allocatable :: Psi_cavity(:), Psi_result_bis(:)
  real(kind=Rkind), allocatable :: REigval(:)
  real(kind=Rkind), allocatable :: REigvec(:,:)

  !----------------------------------Utilities---------------------------------
  integer                       :: i, NB

  !------------------------------Tests in progress-----------------------------
  real(kind=Rkind)              :: Norm_sys
  real(kind=Rkind), allocatable :: Intermediary(:,:), Intermediary2(:,:)
  real(kind=Rkind), allocatable :: Result_total_WF(:,:)
  real(kind=Rkind), allocatable :: Psi_result(:,:), TotH_psi(:,:)
  !real(kind=Rkind), allocatable :: Psi_mapped_result(:)
  !real(kind=Rkind), allocatable :: Cavity_hamiltonian_1DSystem_WF(:,:)
  real(kind=Rkind), allocatable :: Psi_test1(:,:)
  integer                       :: j
  

  
!-----------------------------SYSTEM INITIALIZATION----------------------------
  !-------------Diatomic molecule in a harmonic electonic potential------------
  WRITE(out_unit,*) "-----------------------------SYSTEM INITIALIZATION----------------------------"
  
  CALL MolecCav_Read_cavity_mode(Mode=Molecule_1, nio=in_unit)

  CALL Construct_Operator_1D(Operator=H_ho_molecule_1, &
                                 & operator_type="Hamiltonian", &
                                 & Mode=Molecule_1)!, &
                                 !& Debug=.TRUE.)

  WRITE(out_unit,*) "Molecular Hamiltonian :"
  CALL Write_Operator_1D(H_ho_molecule_1)
  
  CALL Construct_Operator_1D(Operator=x_ho_molecule_1, &
                                 & operator_type="Position", &
                                 & Dense=.TRUE., &
                                 & Mode=Molecule_1)

  WRITE(out_unit,*) "Molecular Position :"
  CALL Write_Mat(x_ho_molecule_1%Dense_val_R, out_unit, Molecule_1%Nb)

  CALL Construct_Operator_1D(Operator=N_ho_molecule_1, &
                                 & operator_type="Nb_photons", &
                                 & Mode=Molecule_1)
  
  WRITE(out_unit,*) "Molecular Nb_photon (vibrationnal excitation quanta)"
  CALL Write_Vec(N_ho_molecule_1%Diag_val_R, out_unit, 1)

  CALL Construct_Operator_1D(Operator=Matter_dipolar_moment, &
                                 & operator_type="Position", &                 ! initialized as a position operator because of approximation over its expression (cf. readme.md or manual)
  !                               & Dense=.TRUE., &                            !/!\ if initialize as dense MUST change two lines below : Matter_dipolar_moment%Dense_val_R /!\
                                 & Mode=Molecule_1)

  Cte_dipole_moment = ONE
  Matter_dipolar_moment%Band_val_R = Matter_dipolar_moment%Band_val_R*Cte_dipole_moment ! /!\ so that the matrix already contains the intensity constant of the dipolar moment with the position of the matter (cf. manual for formulas)
    
  WRITE(out_unit,*) "Molecular Dipole moment :"
  CALL Write_Mat(Matter_dipolar_moment%Band_val_R, out_unit, 3)
  !CALL Write_Mat(Matter_dipolar_moment%Dense_val_R, out_unit, Molecule_1%Nb)

  !------------------------------First cavity mode-----------------------------
  CALL MolecCav_Read_cavity_mode(Mode=Cavity_mode_1, nio=in_unit)

  CALL Construct_Operator_1D(Operator=H_ho_cavity_mode_1, &
                                 & operator_type="Hamiltonian", &
                                 & Mode=Cavity_mode_1)

  WRITE(out_unit,*) "Cavity mode Hamiltonian :"
  CALL Write_Vec(H_ho_cavity_mode_1%Diag_val_R, out_unit, 1)

  CALL Construct_Operator_1D(Operator=x_ho_cavity_mode_1, &
                                 & operator_type="Position", &
                                 & Dense=.TRUE., &
                                 & Mode=Cavity_mode_1)

  WRITE(out_unit,*) "Cavity mode Position :"
  CALL Write_Mat(x_ho_cavity_mode_1%Dense_val_R, out_unit, Molecule_1%Nb)

  CALL Construct_Operator_1D(Operator=N_ho_cavity_mode_1, &
                                 & operator_type="Nb_photons", &
                                 & Mode=Cavity_mode_1)


  WRITE(out_unit,*) "Cavity mode Nb_photons :"
  CALL Write_Vec(N_ho_cavity_mode_1%Diag_val_R, out_unit, 1)
  FLUSH(out_unit)
  
  !-----------------------------Total Wavefunctions----------------------------
  WRITE(out_unit,*) "-----------------------------Total Wavefunctions----------------------------"

  ALLOCATE(System_WF(Molecule_1%Nb,Cavity_mode_1%Nb))
  System_WF(:,:) = ZERO
  DO i = 1, MIN(Molecule_1%Nb, Cavity_mode_1%Nb)                               ! initialize Systel_WF arbitrarily
    System_WF(i,i) = i
  END DO

  !WRITE(out_unit,*) "System_WF"
  !DO i = 1, Molecule_1%Nb
  !  WRITE(out_unit,*) System_WF(i,:)
  !END DO
  !FLUSH(out_unit)                                                             ! without the flush the whole matrix is not written but stops at the 10th line ? Not anymore now the afterward error is fixed
  WRITE(out_unit,*) "Not normalized System_WF total wavefunction"
  CALL Write_Mat(System_WF, out_unit, Cavity_mode_1%Nb)
  
  ALLOCATE(System_WF_mapped(Molecule_1%Nb * Cavity_mode_1%Nb))
  CALL Mapping_WF_2DTO1D(System_WF_mapped, System_WF)
  WRITE(out_unit,*) "Not normalized System_WF_mapped"
  CALL Write_Vec(System_WF_mapped, out_unit, 1)

  CALL Normalize(System_WF)
  WRITE(out_unit,*) "Normalized System_WF wavefunction"
  CALL Write_Mat(System_WF, out_unit, Cavity_mode_1%Nb)

  CALL Normalize(System_WF_mapped)
  WRITE(out_unit,*) "Normalized System_WF_mapped"
  CALL Write_Vec(System_WF_mapped, out_unit, 1)

  !----------------------Action of the Matter Hamiltonian----------------------
  WRITE(out_unit,*) "----------------------Action of the Matter Hamiltonian----------------------"

  ALLOCATE(Matter_hamiltonianSystem_WF(Molecule_1%Nb, Cavity_mode_1%Nb))
  Matter_hamiltonianSystem_WF = ZERO

  DO i = 1, Cavity_mode_1%Nb                                                   ! initialize Matter_hamiltonianSystem_WF by applying the matter hamiltonian to each column of the matrix of the total system WF but creates an sysmalloc : assertion failed if we initialize like that
  CALL Action_Operator_1D(Matter_hamiltonianSystem_WF(:,i), &
                                     & H_ho_molecule_1, &
                                     & System_WF(:,i))
  END DO

  WRITE(out_unit,*) "Action of the matter hamiltonian over System_WF"
  CALL Write_Mat(Matter_hamiltonianSystem_WF, out_unit, Cavity_mode_1%Nb)

  !----------------Computation of the photon number of System_WF---------------
  WRITE(out_unit,*) "----------------Computation of the photon number of System_WF---------------"

  ALLOCATE(Psi_result(Molecule_1%Nb, Cavity_mode_1%Nb))

  CALL MolecCav_Action_operator_2D(Psi_result, N_ho_cavity_mode_1, System_WF)
  WRITE(out_unit,*) "Action of Nb photon operator over System_WF"
  CALL Write_Mat(Psi_result, out_unit, Cavity_mode_1%Nb)

  CALL Scalar_product(Average, Psi_result, System_WF)
  WRITE(out_unit,*) "The average nb of photons of the normalised System_WF is (first method) : ", Average

  CALL MolecCav_Average_value_operator_2D(Average, N_ho_cavity_mode_1, System_WF)
  WRITE(out_unit,*) "The average nb of photons of the normalised System_WF is (second method) : ", Average

  DEALLOCATE(Psi_result)

  !-----------------------An experiment on the Nb_photons----------------------
  WRITE(out_unit,*) "-----------------------An experiment on the Nb_photons----------------------"

  ALLOCATE(Psi_cavity(Cavity_mode_1%Nb))
  Psi_cavity = ZERO
  Psi_cavity(1) = ONE
  Psi_cavity(2) = ONE                                                          ! \ket{Psi} = \ket{0} + \ket{1}
  WRITE(out_unit,*) 'Not normalized Psi_cavity'
  CALL Write_Vec(Psi_cavity, out_unit, 1)
  CALL Normalize(Psi_cavity)
  WRITE(out_unit,*) 'Normalized Psi_cavity (N.B. \frac{1}{\sqrt{2}} = 0.7071067811865475)'
  CALL Write_Vec(Psi_cavity, out_unit, 1)

  ALLOCATE(Psi_result_bis(Cavity_mode_1%Nb))
  CALL Action_Operator_1D(Psi_result_bis, N_ho_cavity_mode_1, Psi_cavity)
  WRITE(out_unit,*) 'Action of Nb_photons over Psi_cavity'
  CALL Write_Vec(Psi_result_bis, out_unit, 1)

  CALL Average_value_operator_1D(Average, N_ho_cavity_mode_1, Psi_cavity)
  WRITE(out_unit,*) 'The avegeraged number of photons of that Psi_cavity is analytically expected to be &
                    & 0.5 and is = ', Average
  
  DEALLOCATE(Psi_result_bis)

  !----------------An experiment on the total Hamiltonian action---------------
  WRITE(out_unit,*) "----------------------Action of the total Hamiltonian----------------------"

  WRITE(out_unit,*) "System_WF"
  CALL Write_Mat(System_WF, out_unit, Cavity_mode_1%Nb)
  FLUSH(out_unit) 

  ALLOCATE(Psi_result(Molecule_1%Nb, Cavity_mode_1%Nb))
  Psi_result(:,:) = ZERO

  Cavity_mode_1%lambda = ZERO
  CALL Action_total_hamiltonian_1p1D(Psi_result, Matter_hamiltonianSystem_WF, &
                                          & H_ho_cavity_mode_1, x_ho_cavity_mode_1, &
                                          & Matter_dipolar_moment, System_WF, &
                                          & Cavity_mode_1)
  Cavity_mode_1%lambda = ONE

  WRITE(out_unit,*) "Action of the uncoupled total Hamiltonian over System_WF"
  CALL Write_Mat(Psi_result, out_unit, Cavity_mode_1%Nb)
  FLUSH(out_unit) 

  !----------------Construction of the Total Hamiltonian matrix----------------
  WRITE(out_unit,*) "----------------Construction of the Total Hamiltonian matrix----------------"

  NB = Molecule_1%Nb * Cavity_mode_1%Nb
  ALLOCATE(H_tot(NB, NB))

  CALL Construct_H_tot(H_tot, Molecule_1%Nb, Cavity_mode_1%Nb, Matter_hamiltonianSystem_WF, &
                              & H_ho_cavity_mode_1, x_ho_cavity_mode_1, Matter_dipolar_moment, &
                              & Cavity_mode_1)  

  WRITE(out_unit,*) "Total Hamiltonian (mapped, lambda /= 0, w_C /= w_M) (10:10 slicing)"
  !CALL Write_Mat(H_tot, out_unit, NB)
  CALL Write_Mat(H_tot(1:10,1:10), out_unit, 10)

  ! JUST SOME TESTS FOR THE SUBROUTINES/CONSTRUCTION OF HTOT WITH DIFFERENT DIPOLAR MOMENT VALUES :
  !Matter_dipolar_moment%Band_val_R = Matter_dipolar_moment%Band_val_R/Cte_dipole_moment ! 
  !Cte_dipole_moment = ONE
  !Matter_dipolar_moment%Band_val_R = Matter_dipolar_moment%Band_val_R*Cte_dipole_moment ! /!\ so that the matrix already contains the intensity constant of the dipolar moment with the position of the matter (cf. manual for formulas)
  !ALLOCATE(H_tot2(NB, NB))
  !CALL Construct_H_tot(H_tot2, Molecule_1%Nb, Cavity_mode_1%Nb, Matter_hamiltonianSystem_WF, &
  !                            & H_ho_cavity_mode_1, x_ho_cavity_mode_1, Matter_dipolar_moment, &
  !                            & Cte_dipole_moment, Cavity_mode_1%lambda, Cavity_mode_1%w)

  !DO i = 1, NB
  !  DO j =1, NB
  !    IF (ABS(H_tot2(i,j) - H_tot(i,j)) > 1E-08_Rkind) THEN
  !      WRITE(out_unit,*) "i,j,H_tot2(i,j),H_tot(i,j) = ", i, j, H_tot2(i,j), H_tot(i,j)
  !    END IF
  !  END DO
  !END DO

  !-----------------------Computation of some observables----------------------

  CALL Average_value_H_tot(Average, H_tot, System_WF_mapped)
  WRITE(out_unit,*) "Average E_tot = ", Average, "Ha"

  !---------------------------Computation Eigenstates--------------------------

  ALLOCATE(REigval(NB))
  ALLOCATE(REigvec(NB,NB))
  CALL diagonalization(H_tot, REigval, Reigvec)

  WRITE(out_unit,*) 'EIGENVALUES'
  !CALL WRITE_Vec(Reigval, out_unit, 6, info = 'VP[Ha]')
  CALL WRITE_Vec(Reigval(1:10), out_unit, 10, info = 'VP[Ha]')

  !WRITE(out_unit,*) 'EIGENVECTORS'
  !CALL WRITE_Mat(Reigvec, out_unit, 6, info = 'Eigenvectors')

  !-------Construction of a Total Hamiltonian matrix without CM-couplings------
  WRITE(out_unit,*) "-------Construction of a Total Hamiltonian matrix without CM-couplings------"

  DEALLOCATE(H_tot)
  ALLOCATE(H_tot(NB, NB))

  Cavity_mode_1%lambda = ZERO
  CALL Construct_H_tot(H_tot, Molecule_1%Nb, Cavity_mode_1%Nb, Matter_hamiltonianSystem_WF, &
                              & H_ho_cavity_mode_1, x_ho_cavity_mode_1, Matter_dipolar_moment, &
                              & Cavity_mode_1)  

  WRITE(out_unit,*) "Total Hamiltonian (mapped, lambda = 0, w_C /= w_M) (10:10 slicing)"
  !CALL Write_Mat(H_tot, out_unit, NB)
  CALL Write_Mat(H_tot(1:10,1:10), out_unit, 10)

  !---------------------------Computation Eigenstates--------------------------
  DEALLOCATE(REigval)
  DEALLOCATE(REigvec)
  ALLOCATE(REigval(NB))
  ALLOCATE(REigvec(NB,NB))
  CALL diagonalization(H_tot, REigval, Reigvec)
  !Missing : procedure to write diagonal H_tot : no need, assume it worked and diag. elements = Eigenvalues
  !WRITE(out_unit,*) 'Total Hamiltonian (mapped, lambda = 0, diagonalized)'
  !CALL Write_Mat(H_tot, out_unit, 42)

  WRITE(out_unit,*) 'EIGENVALUES'
  !CALL Write_Vec(Reigval, out_unit, 6, info = 'VP[Ha]')
  CALL Write_Vec(Reigval(1:10), out_unit, 10, info = 'VP[Ha]')

  !WRITE(out_unit,*) 'EIGENVECTORS'
  !CALL Write_Mat(Reigvec, out_unit, 6, info = 'Eigenvectors')

  !--------Construction of a Total Hamiltonian matrix with CM-couplings--------
  WRITE(out_unit,*) "-------Construction of a Total Hamiltonian matrix without CM-couplings------"

  DEALLOCATE(H_tot)
  ALLOCATE(H_tot(NB, NB))

  Cavity_mode_1%lambda = TWO
  Cavity_mode_1%w = Molecule_1%w
  CALL Construct_H_tot(H_tot, Molecule_1%Nb, Cavity_mode_1%Nb, Matter_hamiltonianSystem_WF, &
                              & H_ho_cavity_mode_1, x_ho_cavity_mode_1, Matter_dipolar_moment, &
                              & Cavity_mode_1)  

  WRITE(out_unit,*) "Total Hamiltonian (mapped, lambda = 1.0, w_C = w_M) (10:10 slicing)"
  !CALL Write_Mat(H_tot, out_unit, NB)
  CALL Write_Mat(H_tot(1:10,1:10), out_unit, 10)

  !---------------------------Computation Eigenstates--------------------------
  DEALLOCATE(REigval)
  DEALLOCATE(REigvec)
  ALLOCATE(REigval(NB))
  ALLOCATE(REigvec(NB,NB))
  CALL diagonalization(H_tot, REigval, Reigvec)

  WRITE(out_unit,*) 'EIGENVALUES (10 firsts)'
  !CALL Write_Vec(Reigval, out_unit, 6, info = 'VP[Ha]')
  CALL Write_Vec(Reigval(1:10), out_unit, 10, info = 'VP[Ha]')
  
  !WRITE(out_unit,*) 'EIGENVECTORS'
  !CALL Write_Mat(Reigvec, out_unit, 6, info = 'Eigenvectors')


  !---------------------------tests--------------------------
  Norm_sys = ZERO
  ALLOCATE(Psi_test1(2,3))
  DO i = 1, 2
    DO j = 1, 3
    Psi_test1(i,j) = i+j
    END DO 
  END DO

  CALL Scalar_product(Norm_sys, Psi_test1, Psi_test1)

  ALLOCATE(TotH_psi(Molecule_1%Nb, Cavity_mode_1%Nb))
  CALL Action_total_hamiltonian_1p1D(TotH_psi, x_ho_cavity_mode_1, H_ho_cavity_mode_1, &
  &Matter_dipolar_moment, H_ho_molecule_1, System_WF, Debug=.TRUE.)

END PROGRAM