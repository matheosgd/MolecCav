PROGRAM App_MolecCav
  USE QDUtil_m
  USE MC_cavity_mode_m
  USE MC_operator_1D_m
  USE MC_total_hamiltonian_m
  USE MC_algebra_m
  IMPLICIT NONE


  !-------------Diatomic molecule in a harmonic electonic potential------------
  TYPE(MC_cavity_mode_t)        :: Molecule_1
  TYPE(MC_operator_1D_t)        :: H_ho_molecule_1                             ! matrix of the one-dimensional harmonic Hamiltonian associated with HO D
  TYPE(MC_operator_1D_t)        :: x_ho_molecule_1
  TYPE(MC_operator_1D_t)        :: N_ho_molecule_1                             ! here N does not count the number of photons but of excitation quanta in the vibrational state
  TYPE(MC_operator_1D_t)        :: Matter_dipolar_moment
  real(kind=Rkind)              :: Cte_dipole_moment                           ! the intensity of the variation of the dipole moment with a variation of the matter DOF
  
  !------------------------------First cavity mode-----------------------------
  TYPE(MC_cavity_mode_t)        :: Cavity_mode_1
  TYPE(MC_operator_1D_t)        :: H_ho_cavity_mode_1                          ! matrix of the one-dimensional harmonic Hamiltonian associated with HO D
  TYPE(MC_operator_1D_t)        :: x_ho_cavity_mode_1
  TYPE(MC_operator_1D_t)        :: N_ho_cavity_mode_1

  !--------------------------------Wavefunctions-------------------------------
  real(kind=Rkind), allocatable :: System_WF(:,:)                              ! The total system (matter-cavity) wavefunction. Size Nb_M*Nb_C. |System_WF> = |Molecule_WF>.TENSOR.|Cavity_WF> 
  real(kind=Rkind), allocatable :: Matter_hamiltonianSystem_WF(:,:)            ! The wavefunction resulting from the action of the matter hamiltonian on System_WF. Size Nb_M*Nb_C. |H_MatterSystem_WF(:,i_C)> = H_Matter|System_WF(:,i_C)>
  !real(kind=Rkind), allocatable :: Result_total_WF(:,:)
  real(kind=Rkind), allocatable :: H_tot(:,:)                               

  !----------------------------------Utilities---------------------------------
  integer                       :: i, NB

  !------------------------------Tests in progress-----------------------------
  real(kind=Rkind)              :: Norm_sys, Projection, Average
  real(kind=Rkind), allocatable :: System_WF_mapped(:)
  real(kind=Rkind), allocatable :: Psi_result(:,:)
  !real(kind=Rkind), allocatable :: H_tot2(:,:)                               
  !integer                       :: j
  real(kind=Rkind), allocatable :: Psi_mapped_result(:)

  
!-----------------------------SYSTEM INITIALIZATION----------------------------
  !-------------Diatomic molecule in a harmonic electonic potential------------
  CALL MolecCav_Read_cavity_mode(Mode=Molecule_1, nio=in_unit)

  CALL MolecCav_Construct_Operator(Operator=H_ho_molecule_1, &
                                 & operator_type="Hamiltonian", &
                                 & scalar_space="Real", &
                                 & matrix_shape_type="Opt", &                  ! opt => get analytical shape. non_opt => get dense shape
                                 & Nb=Molecule_1%Nb, &
                                 & w=Molecule_1%w, &
                                 & m=Molecule_1%m)

  CALL MolecCav_Construct_Operator(Operator=x_ho_molecule_1, &
                                 & operator_type="Position", &
                                 & scalar_space="Real", &
                                 & matrix_shape_type="Opt", &                  ! opt => get analytical shape. non_opt => get dense shape
                                 & Nb=Molecule_1%Nb, &
                                 & w=Molecule_1%w, &
                                 & m=Molecule_1%m)

  CALL MolecCav_Construct_Operator(Operator=N_ho_molecule_1, &
                                 & operator_type="Nb_photons", &
                                 & scalar_space="Real", &
                                 & matrix_shape_type="Opt", &                  ! opt => get analytical shape. non_opt => get dense shape
                                 & Nb=Molecule_1%Nb, &
                                 & w=Molecule_1%w, &
                                 & m=Molecule_1%m)
  
  CALL MolecCav_Construct_Operator(Operator=Matter_dipolar_moment, &
                                 & operator_type="Position", &                 ! initialized as a position operator because of approximation over its expression (cf. readme.md or manual)
                                 & scalar_space="Real", &
                                 & matrix_shape_type="Opt", &                  ! opt => get analytical shape. non_opt => get dense shape
                                 & Nb=Molecule_1%Nb, &
                                 & w=Molecule_1%w, &
                                 & m=Molecule_1%m)

  Cte_dipole_moment = TWELVE
  Matter_dipolar_moment%Band_val_R = Matter_dipolar_moment%Band_val_R*Cte_dipole_moment ! /!\ so that the matrix already contains the intensity constant of the dipolar moment with the position of the matter (cf. manual for formulas)
                                
  !------------------------------First cavity mode-----------------------------
  CALL MolecCav_Read_cavity_mode(Mode=Cavity_mode_1, nio=in_unit)

  CALL MolecCav_Construct_Operator(Operator=H_ho_cavity_mode_1, &
                                 & operator_type="Hamiltonian", &
                                 & scalar_space="Real", &
                                 & matrix_shape_type="Opt", &                  ! opt => get analytical shape. non_opt => get dense shape
                                 & Nb=Cavity_mode_1%Nb, &
                                 & w=Cavity_mode_1%w, &
                                 & m=Cavity_mode_1%m)

  CALL MolecCav_Construct_Operator(Operator=x_ho_cavity_mode_1, &
                                 & operator_type="Position", &
                                 & scalar_space="Real", &
                                 & matrix_shape_type="Opt", &                  ! opt => get analytical shape. non_opt => get dense shape
                                 & Nb=Cavity_mode_1%Nb, &
                                 & w=Cavity_mode_1%w, &
                                 & m=Cavity_mode_1%m)

  CALL MolecCav_Construct_Operator(Operator=N_ho_cavity_mode_1, &
                                 & operator_type="Nb_photons", &
                                 & scalar_space="Real", &
                                 & matrix_shape_type="Opt", &                  ! opt => get analytical shape. non_opt => get dense shape
                                 & Nb=Cavity_mode_1%Nb, &
                                 & w=Cavity_mode_1%w, &
                                 & m=Cavity_mode_1%m)


!---------------------------------Wavefunctions--------------------------------
!------------------Part of the code using the MolecCav library-----------------
  ALLOCATE(System_WF(Molecule_1%Nb,Cavity_mode_1%Nb))
  ALLOCATE(Matter_hamiltonianSystem_WF(Molecule_1%Nb,Cavity_mode_1%Nb))
  
  System_WF = ZERO
  Matter_hamiltonianSystem_WF = ZERO

  DO i = 1, MIN(Molecule_1%Nb, Cavity_mode_1%Nb)                               ! initialize Systel_WF arbitrarily
    System_WF(i,i) = i
  END DO

  WRITE(out_unit,*) "System_WF"
  DO i = 1, Molecule_1%Nb
    WRITE(out_unit,*) System_WF(i,:)
  END DO
  !FLUSH(out_unit)                                                             ! without the flush the whole matrix is not written but stops at the 10th line ? Not anymore now the afterward error is fixed

  ALLOCATE(System_WF_mapped(Molecule_1%Nb * Cavity_mode_1%Nb))
  CALL MolecCav_Mapping_WF_2DTO1D(System_WF_mapped, System_WF)
  WRITE(out_unit,*) "System_WF_mapped"
  DO i = 1, Size(System_WF_mapped)
    WRITE(out_unit,*) System_WF_mapped(i)
  END DO

  CALL MolecCav_Normalize_WF_2D(System_WF)
  WRITE(out_unit,*) "Normalized System_WF"
  DO i = 1, Molecule_1%Nb
    WRITE(out_unit,*) System_WF(i,:)
  END DO

  CALL MolecCav_Normalize_WF_1D(System_WF_mapped)
  WRITE(out_unit,*) "Normalized System_WF_mapped"
  DO i = 1, Size(System_WF_mapped)
    WRITE(out_unit,*) System_WF_mapped(i)
  END DO

  CALL MolecCav_Norm_WF_1D(Norm_sys, System_WF_mapped)
  WRITE(out_unit, *) "the norm of system_wf_mapped is ", Norm_sys, "but those of system_wf cannot be &
                    & computed here because of <<multiple definition of moleccav_norm_wf_2d>> ??"

  DO i = 1, Cavity_mode_1%Nb                                                   ! initialize Matter_hamiltonianSystem_WF by applying the matter hamiltonian to each column of the matrix of the total system WF but creates an sysmalloc : assertion failed if we initialize like that
    CALL MolecCav_Action_Operator_1D(Matter_hamiltonianSystem_WF(:,i), &
                                     & H_ho_molecule_1, &
                                     & System_WF(:,i))
  END DO

  ALLOCATE(Psi_result(Molecule_1%Nb, Cavity_mode_1%Nb))
  CALL MolecCav_Action_cavity_operator_2D(Psi_result, N_ho_cavity_mode_1, System_WF)
  WRITE(out_unit,*) "Result Nb photon op"
  DO i = 1, Size(Psi_result, 1)
    WRITE(out_unit,*) Psi_result(i,:)
  END DO
  CALL MolecCav_scalar_product_2D(Average, Psi_result, System_WF)
  WRITE(out_unit, *) "The average nb of photons of the normalised System_WF is (first method) : ", Average
  CALL MolecCav_Average_value_cavity_operator_2D(Average, N_ho_cavity_mode_1, System_WF)
  WRITE(out_unit, *) "The average nb of photons of the normalised System_WF is (second method) : ", Average

  WRITE(out_unit,*) "Matter_hamiltonianSystem_WF"
  DO i = 1, Molecule_1%Nb
    WRITE(out_unit,*) Matter_hamiltonianSystem_WF(i,:)
  END DO

!-------------------------Part of the MolecCav library-------------------------
  !ALLOCATE(Result_total_WF(Molecule_1%Nb,Cavity_mode_1%Nb))
  !Result_total_WF = ZERO

  !CALL MolecCav_Action_Total_Hamiltonian_1D(Result_total_WF, Matter_hamiltonianSystem_WF, &
  !                                        & H_ho_cavity_mode_1, x_ho_cavity_mode_1, &
  !                                        & Matter_dipolar_moment, Cte_dipole_moment, System_WF, &
  !                                        & Cavity_mode_1%lambda, Cavity_mode_1%w)
  
  !WRITE(out_unit,*) "Result_total_WF"
  !DO i = 1, Molecule_1%Nb
  !  WRITE(out_unit,*) Result_total_WF(i,:)
  !END DO

  NB = Molecule_1%Nb * Cavity_mode_1%Nb
  ALLOCATE(H_tot(NB, NB))

  CALL MolecCav_Construct_H_tot(H_tot, Molecule_1%Nb, Cavity_mode_1%Nb, Matter_hamiltonianSystem_WF, &
                              & H_ho_cavity_mode_1, x_ho_cavity_mode_1, Matter_dipolar_moment, &
                              & Cavity_mode_1%lambda, Cavity_mode_1%w)  

  WRITE(out_unit,*) "H_tot"
  DO i = 1, NB
    WRITE(out_unit,*) H_tot(i,:)
  END DO

  ! JUST SOME TESTS FOR THE SUBROUTINES :
  !Matter_dipolar_moment%Band_val_R = Matter_dipolar_moment%Band_val_R/Cte_dipole_moment ! /!\ so that the matrix already contains the intensity constant of the dipolar moment with the position of the matter (cf. manual for formulas)
  !Cte_dipole_moment = ONE
  !Matter_dipolar_moment%Band_val_R = Matter_dipolar_moment%Band_val_R*Cte_dipole_moment ! /!\ so that the matrix already contains the intensity constant of the dipolar moment with the position of the matter (cf. manual for formulas)
  !ALLOCATE(H_tot2(NB, NB))
  !CALL MolecCav_Construct_H_tot(H_tot2, Molecule_1%Nb, Cavity_mode_1%Nb, Matter_hamiltonianSystem_WF, &
  !                            & H_ho_cavity_mode_1, x_ho_cavity_mode_1, Matter_dipolar_moment, &
  !                            & Cte_dipole_moment, Cavity_mode_1%lambda, Cavity_mode_1%w)

  !DO i = 1, NB
  !  DO j =1, NB
  !    IF (ABS(H_tot2(i,j) - H_tot(i,j)) > 1E-08_Rkind) THEN
  !      WRITE(out_unit,*) "i,j,H_tot2(i,j),H_tot(i,j) = ", i, j, H_tot2(i,j), H_tot(i,j)
  !    END IF
  !  END DO
  !END DO

  CALL MolecCav_Average_value_H_tot(Average, H_tot, System_WF_mapped)
  WRITE(out_unit, *) "Average E_tot = ", Average, "Ha"


END PROGRAM