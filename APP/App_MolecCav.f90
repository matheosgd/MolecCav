PROGRAM App_MolecCav
  USE QDUtil_m
  USE MC_cavity_mode_m
  USE MC_operator_1D_m
  IMPLICIT NONE


!--------------Diatomic molecule in a harmonic electonic potential-------------
  TYPE(MC_cavity_mode_t)        :: Molecule_1
  TYPE(MC_operator_1D_t)        :: H_ho_molecule_1                             ! matrix of the one-dimensional harmonic Hamiltonian associated with HO D
  TYPE(MC_operator_1D_t)        :: x_ho_molecule_1
  TYPE(MC_operator_1D_t)        :: N_ho_molecule_1                             ! here N does not count the number of photons but of excitation quanta in the vibrational state
  TYPE(MC_operator_1D_t)        :: Matter_dipolar_moment

!-------------------------------First cavity mode------------------------------
  TYPE(MC_cavity_mode_t)        :: Cavity_mode_1
  TYPE(MC_operator_1D_t)        :: H_ho_cavity_mode_1                          ! matrix of the one-dimensional harmonic Hamiltonian associated with HO D
  TYPE(MC_operator_1D_t)        :: x_ho_cavity_mode_1
  TYPE(MC_operator_1D_t)        :: N_ho_cavity_mode_1

!---------------------------------Wavefunctions--------------------------------
  real(kind=Rkind), allocatable :: System_WF(:,:)                              ! size Nb_M*Nb_C. |System_WF> = |Molecule_WF>.TENSOR.|Cavity_WF> 
  real(kind=Rkind), allocatable :: Matter_hamiltonianSystem_WF(:,:)            ! size Nb_M*Nb_C. |H_MatterSystem_WF(:,i_C)> = H_Matter|System_WF(:,i_C)>
  real(kind=Rkind), allocatable :: Result_total_WF(:,:)                        ! already allocated !

!-----------------------------------Utilities----------------------------------
  integer                       :: i


!--------------Diatomic molecule in a harmonic electonic potential-------------
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
                                 & operator_type="Nb_photons", &               ! oui, ce n'est pas un op nb de photons
                                 & scalar_space="Real", &
                                 & matrix_shape_type="Opt", &                  ! opt => get analytical shape. non_opt => get dense shape
                                 & Nb=Molecule_1%Nb, &
                                 & w=Molecule_1%w, &
                                 & m=Molecule_1%m)
                                
!--------------------------------First cavity mode-----------------------------
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

  DO i = 1, MIN(Molecule_1%Nb, Cavity_mode_1%Nb)
    System_WF(i,i) = i
  END DO

  WRITE(out_unit,*) "System_WF"
  DO i = 1, Molecule_1%Nb
    WRITE(out_unit,*) System_WF(i,:)
  END DO

  DO i = 1, MIN(Molecule_1%Nb, Cavity_mode_1%Nb)
    Matter_hamiltonianSystem_WF(i,i) = (i+ONETENTH)*ELEVEN
  END DO

  WRITE(out_unit,*) "Matter_hamiltonianSystem_WF"
  DO i = 1, Molecule_1%Nb
    WRITE(out_unit,*) Matter_hamiltonianSystem_WF(i,:)
  END DO

!-------------------------Part of the MolecCav library-------------------------
  ALLOCATE(Result_total_WF(Molecule_1%Nb,Cavity_mode_1%Nb))
  Result_total_WF = ZERO

  CALL MolecCav_Action_Total_Hamiltonian_1D(Result_total_WF, Matter_hamiltonianSystem_WF, &
                                          & H_ho_cavity_mode_1, x_ho_cavity_mode_1, &
                                          & Matter_dipolar_moment, System_WF, &
                                          & Cavity_mode_1%lambda, Cavity_mode_1%w)
  
  WRITE(out_unit,*) "Result_total_WF"
  DO i = 1, Molecule_1%Nb
    WRITE(out_unit,*) Result_total_WF(i,:)
  END DO


  WRITE(out_unit,*) "youhou"
END PROGRAM