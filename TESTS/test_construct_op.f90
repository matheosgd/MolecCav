PROGRAM test_construct_op
  !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
  USE QDUtil_m
  USE MC_operator_1D_m
  IMPLICIT NONE


!--------------------------------First cavity mode-----------------------------
  TYPE(MC_operator_1D_t)        :: H_ho_1D_diag_1_6                            ! matrix of the one-dimensional harmonic Hamiltonian associated with HO D. The first number = eigenpulsation of the HO, the second one = the number of basis functions to expand the matrix in.
  TYPE(MC_operator_1D_t)        :: H_ho_1D_diag_14_6                           ! 14 = w (eigenpulsation); 6 = Nb (number of Ho basis vectors)
  TYPE(MC_operator_1D_t)        :: H_ho_1D_diag_1_17
  TYPE(MC_operator_1D_t)        :: H_ho_1D_diag_14_17
  TYPE(MC_operator_1D_t)        :: H_ho_1D_dense_1_6
  TYPE(MC_operator_1D_t)        :: H_ho_1D_dense_14_17
  real(kind=Rkind), allocatable :: H_ho_1D_diag_theo_14_17(:)
  real(kind=Rkind), allocatable :: H_ho_1D_dense_theo_1_6(:,:)
  real(kind=Rkind), allocatable :: H_ho_1D_dense_theo_14_17(:,:)

  TYPE(MC_operator_1D_t)        :: x_ho_1D_band_1_6_1
  TYPE(MC_operator_1D_t)        :: x_ho_1D_band_14_6_1                         ! 14 = w (eigenpulsation); 6 = Nb (number of Ho basis vectors); 1 = m (mass associated to the HO)
  TYPE(MC_operator_1D_t)        :: x_ho_1D_band_1_17_1
  TYPE(MC_operator_1D_t)        :: x_ho_1D_band_1_6_7
  TYPE(MC_operator_1D_t)        :: x_ho_1D_band_14_17_7
  TYPE(MC_operator_1D_t)        :: x_ho_1D_dense_1_6_1
  TYPE(MC_operator_1D_t)        :: x_ho_1D_dense_14_17_7
  real(kind=Rkind), allocatable :: x_ho_1D_band_theo_1_17_1(:,:)
  real(kind=Rkind), allocatable :: x_ho_1D_dense_theo_1_6_1(:,:)
  real(kind=Rkind), allocatable :: x_ho_1D_dense_theo_14_17_7(:,:)

  TYPE(MC_operator_1D_t)        :: N_ho_1D_diag_6
  TYPE(MC_operator_1D_t)        :: N_ho_1D_diag_17
  TYPE(MC_operator_1D_t)        :: N_ho_1D_dense_6
  TYPE(MC_operator_1D_t)        :: N_ho_1D_dense_17
  real(kind=Rkind), allocatable :: N_ho_1D_dense_theo_17(:,:)

  integer                       :: i, error_H, error_x, error_N
  real(kind=Rkind)              :: threshold_H, threshold_x, threshold_N       ! the accuracy thresholds for the comparisons between matrices elements of the Hamiltonian, position and number of photon operators


  !-----------------------construct HO matricies to test-----------------------
  CALL MolecCav_Construct_Operator(Operator=H_ho_1D_diag_1_6, &
                                 & operator_type="Hamiltonian", &
                                 & scalar_space="Real", &
                                 & matrix_shape_type="Opt", &                  ! opt => get analytical shape. non_opt => get dense shape
                                 & Nb=6, &
                                 & w=ONE, &
                                 & m=ONE)

  CALL MolecCav_Construct_Operator(Operator=H_ho_1D_diag_14_6, &
                                 & operator_type="Hamiltonian", &
                                 & scalar_space="Real", &
                                 & matrix_shape_type="Opt", &                  ! opt => get analytical shape. non_opt => get dense shape
                                 & Nb=6, &
                                 & w=14.0_Rkind, &
                                 & m=ONE)

  CALL MolecCav_Construct_Operator(Operator=H_ho_1D_diag_1_17, &
                                 & operator_type="Hamiltonian", &
                                 & scalar_space="Real", &
                                 & matrix_shape_type="Opt", &                  ! opt => get analytical shape. non_opt => get dense shape
                                 & Nb=17, &
                                 & w=ONE, &
                                 & m=ONE)

  CALL MolecCav_Construct_Operator(Operator=H_ho_1D_diag_14_17, &
                                 & operator_type="Hamiltonian", &
                                 & scalar_space="Real", &
                                 & matrix_shape_type="Opt", &                  ! opt => get analytical shape. non_opt => get dense shape
                                 & Nb=17, &
                                 & w=14.0_Rkind, &
                                 & m=ONE)

  CALL MolecCav_Construct_Operator(Operator=H_ho_1D_dense_1_6, &
                                 & operator_type="Hamiltonian", &
                                 & scalar_space="Real", &
                                 & matrix_shape_type="Non_opt", &              ! opt => get analytical shape. non_opt => get dense shape
                                 & Nb=6, &
                                 & w=ONE, &
                                 & m=ONE)
                                 
  CALL MolecCav_Construct_Operator(Operator=H_ho_1D_dense_14_17, &
                                 & operator_type="Hamiltonian", &
                                 & scalar_space="Real", &
                                 & matrix_shape_type="Non_opt", &              ! opt => get analytical shape. non_opt => get dense shape
                                 & Nb=17, &
                                 & w=14.0_Rkind, &
                                 & m=ONE)


  !------------------------construct x matricies to test-----------------------
  CALL MolecCav_Construct_Operator(Operator=x_ho_1D_band_1_6_1, &
                                 & operator_type="Position", &
                                 & scalar_space="Real", &
                                 & matrix_shape_type="Opt", &                  ! opt => get analytical shape. non_opt => get dense shape
                                 & Nb=6, &
                                 & w=ONE, &
                                 & m=ONE)

  CALL MolecCav_Construct_Operator(Operator=x_ho_1D_band_14_6_1, &
                                 & operator_type="Position", &
                                 & scalar_space="Real", &
                                 & matrix_shape_type="Opt", &                  ! opt => get analytical shape. non_opt => get dense shape
                                 & Nb=6, &
                                 & w=14.0_Rkind, &
                                 & m=ONE)

  CALL MolecCav_Construct_Operator(Operator=x_ho_1D_band_1_17_1, &
                                 & operator_type="Position", &
                                 & scalar_space="Real", &
                                 & matrix_shape_type="Opt", &                  ! opt => get analytical shape. non_opt => get dense shape
                                 & Nb=17, &
                                 & w=ONE, &
                                 & m=ONE)

  CALL MolecCav_Construct_Operator(Operator=x_ho_1D_band_1_6_7, &
                                 & operator_type="Position", &
                                 & scalar_space="Real", &
                                 & matrix_shape_type="Opt", &                  ! opt => get analytical shape. non_opt => get dense shape
                                 & Nb=6, &
                                 & w=ONE, &
                                 & m=SEVEN)

  CALL MolecCav_Construct_Operator(Operator=x_ho_1D_band_14_17_7, &
                                 & operator_type="Position", &
                                 & scalar_space="Real", &
                                 & matrix_shape_type="Opt", &                  ! opt => get analytical shape. non_opt => get dense shape
                                 & Nb=17, &
                                 & w=14.0_Rkind, &
                                 & m=SEVEN)

  CALL MolecCav_Construct_Operator(Operator=x_ho_1D_dense_1_6_1, &
                                 & operator_type="Position", &
                                 & scalar_space="Real", &
                                 & matrix_shape_type="Non_opt", &              ! opt => get analytical shape. non_opt => get dense shape
                                 & Nb=6, &
                                 & w=ONE, &
                                 & m=ONE)

  CALL MolecCav_Construct_Operator(Operator=x_ho_1D_dense_14_17_7, &
                                 & operator_type="Position", &
                                 & scalar_space="Real", &
                                 & matrix_shape_type="Non_opt", &              ! opt => get analytical shape. non_opt => get dense shape
                                 & Nb=17, &
                                 & w=14.0_Rkind, &
                                 & m=SEVEN)

  !------------------------construct N matricies to test-----------------------
  CALL MolecCav_Construct_Operator(Operator=N_ho_1D_diag_6, &
                                 & operator_type="Nb_photons", &
                                 & scalar_space="Real", &
                                 & matrix_shape_type="Opt", &                  ! opt => get analytical shape. non_opt => get dense shape
                                 & Nb=6, &
                                 & w=ONE, &
                                 & m=ONE)

  CALL MolecCav_Construct_Operator(Operator=N_ho_1D_diag_17, &
                                 & operator_type="Nb_photons", &
                                 & scalar_space="Real", &
                                 & matrix_shape_type="Opt", &                  ! opt => get analytical shape. non_opt => get dense shape
                                 & Nb=17, &
                                 & w=ONE, &
                                 & m=ONE)

  CALL MolecCav_Construct_Operator(Operator=N_ho_1D_dense_6, &
                                 & operator_type="Nb_photons", &
                                 & scalar_space="Real", &
                                 & matrix_shape_type="Non_opt", &              ! opt => get analytical shape. non_opt => get dense shape
                                 & Nb=6, &
                                 & w=ONE, &
                                 & m=ONE)

  CALL MolecCav_Construct_Operator(Operator=N_ho_1D_dense_17, &
                                 & operator_type="Nb_photons", &
                                 & scalar_space="Real", &
                                 & matrix_shape_type="Non_opt", &              ! opt => get analytical shape. non_opt => get dense shape
                                 & Nb=17, &
                                 & w=ONE, &
                                 & m=ONE)


  !----------------------construct reference HO matricies----------------------
  ALLOCATE(H_ho_1D_diag_theo_14_17(17))
  ALLOCATE(H_ho_1D_dense_theo_1_6(6,6))
  ALLOCATE(H_ho_1D_dense_theo_14_17(17,17))

  H_ho_1D_diag_theo_14_17  = ZERO
  H_ho_1D_dense_theo_1_6   = ZERO
  H_ho_1D_dense_theo_14_17 = ZERO

  DO i = 1, 17
    H_ho_1D_diag_theo_14_17(i)    = 14*(i - ONE + HALF)
    H_ho_1D_dense_theo_14_17(i,i) = 14*(i - ONE + HALF)
    IF (i <= 6) THEN
      H_ho_1D_dense_theo_1_6(i,i) = (i - ONE + HALF)
    END IF
  END DO

  WRITE(out_unit,*) 
  WRITE(out_unit,*) "H_{Reference 1 (diagonal, w = 14, Nb = 17)}"
  DO i = 1, 17
    WRITE(out_unit,*) H_ho_1D_diag_theo_14_17(i)
  END DO

  WRITE(out_unit,*) 
  WRITE(out_unit,*) "H_{Reference 2 (dense, w = 1, Nb = 6)}"
  DO i = 1, 6
    WRITE(out_unit,*) H_ho_1D_dense_theo_1_6(i,:)
  END DO

  WRITE(out_unit,*) 
  WRITE(out_unit,*) "H_{Reference 3 (dense, w = 14, Nb = 17)}"
  DO i = 1, 17
    WRITE(out_unit,*) H_ho_1D_dense_theo_14_17(i,:)
  END DO


  !----------------------construct reference x matricies---------------------
  ALLOCATE(x_ho_1D_band_theo_1_17_1(17,3))
  ALLOCATE(x_ho_1D_dense_theo_1_6_1(6,6))
  ALLOCATE(x_ho_1D_dense_theo_14_17_7(17,17))

  x_ho_1D_band_theo_1_17_1   = ZERO
  x_ho_1D_dense_theo_1_6_1   = ZERO
  x_ho_1D_dense_theo_14_17_7 = ZERO

  DO i = 1, 16                                                                 ! /!\ Fortran counts from 1 to Nb !!! /!\ Nb-1 not to have Band_val_R(i+1) out of range
    x_ho_1D_band_theo_1_17_1(i,1)     = SQRT(REAL(i,kind=Rkind))
    x_ho_1D_band_theo_1_17_1(i+1,3)   = SQRT(REAL(i,kind=Rkind))
    x_ho_1D_dense_theo_14_17_7(i,i+1) = SQRT(REAL(i,kind=Rkind))
    x_ho_1D_dense_theo_14_17_7(i+1,i) = SQRT(REAL(i,kind=Rkind))
    IF (i <= 5) THEN
      x_ho_1D_dense_theo_1_6_1(i,i+1) = SQRT(REAL(i,kind=Rkind))
      x_ho_1D_dense_theo_1_6_1(i+1,i) = SQRT(REAL(i,kind=Rkind))
    END IF
  END DO

  x_ho_1D_band_theo_1_17_1   = x_ho_1D_band_theo_1_17_1   / SQRT(TWO * 1.0_Rkind  * 1.0_Rkind)
  x_ho_1D_dense_theo_1_6_1   = x_ho_1D_dense_theo_1_6_1   / SQRT(TWO * 1.0_Rkind  * 1.0_Rkind)
  x_ho_1D_dense_theo_14_17_7 = x_ho_1D_dense_theo_14_17_7 / SQRT(TWO * 14.0_Rkind * 7.0_Rkind)

  WRITE(out_unit,*) 
  WRITE(out_unit,*) "x_{Reference 1 (band, w = 1.0, Nb = 17, m = 1.0)}"
  DO i = 1, 17
    WRITE(out_unit,*) x_ho_1D_band_theo_1_17_1(i,:)
  END DO

  WRITE(out_unit,*) 
  WRITE(out_unit,*) "x_{Reference 2 (dense, w = 1.0, Nb = 6, m = 1.0)}"
  DO i = 1, 6
    WRITE(out_unit,*) x_ho_1D_dense_theo_1_6_1(i,:)
  END DO


  WRITE(out_unit,*) 
  WRITE(out_unit,*) "x_{Reference 3 (dense, w = 14.0, Nb = 17, m = 7.0)}"
  DO i = 1, 17
    WRITE(out_unit,*) x_ho_1D_dense_theo_14_17_7(i,:)
  END DO

  !----------------------construct reference N matricies---------------------
  ALLOCATE(N_ho_1D_dense_theo_17(17,17))
  N_ho_1D_dense_theo_17 = ZERO

  DO i = 1, 17                                                             ! /!\ Fortran counts from 1 to Nb !!! /!\
    N_ho_1D_dense_theo_17(i,i) = i - 1
  END DO

  WRITE(out_unit,*) 
  WRITE(out_unit,*) "N_{Reference 1 (dense, Nb = 17)}"
  DO i = 1, 17
    WRITE(out_unit,*) N_ho_1D_dense_theo_17(i,:)
  END DO


  !-------------------------------comparisons HO-------------------------------
  error_H = 0
  threshold_H = 1E-08_Rkind

  !H_ho_1D_diag_14_17%Diag_val_R(1) = 0
  !WRITE(out_unit,*) 
  !WRITE(out_unit,*) "H_{Spurious 1 (diagonal, w = 14, Nb = 17)}"
  !DO i = 1, 17
  !  WRITE(out_unit,*) H_ho_1D_diag_theo_14_17(i)
  !END DO

  IF (ANY(ABS(H_ho_1D_diag_theo_14_17 - H_ho_1D_diag_14_17%Diag_val_R) > threshold_H)) THEN
    WRITE(out_unit,*) 
    WRITE(out_unit,*) 'H_ho_1D_diag_14_17 failed to initialize'
    error_H = 1
  END IF

  !H_ho_1D_dense_1_6%Dense_val_R(1,2) = 1
  !WRITE(out_unit,*) 
  !WRITE(out_unit,*) "H_{Spurious 2 (diagonal, w = 1, Nb = 6)}"
  !DO i = 1, 17
  !  WRITE(out_unit,*) H_ho_1D_dense_1_6%Dense_val_R(i,:)
  !END DO

  IF (ANY(ABS(H_ho_1D_dense_theo_1_6 - H_ho_1D_dense_1_6%Dense_val_R) > threshold_H)) THEN
    WRITE(out_unit,*) 
    WRITE(out_unit,*) 'H_ho_1D_dense_1_6 failed to initialize'
    error_H = 1
  END IF

  !H_ho_1D_dense_14_17%Dense_val_R(1,2) = 1
  !WRITE(out_unit,*) 
  !WRITE(out_unit,*) "H_{Spurious 3 (diagonal, w = 14, Nb = 17)}"
  !DO i = 1, 17
  !  WRITE(out_unit,*) H_ho_1D_dense_14_17%Dense_val_R(i,:)
  !END DO

  IF (ANY(ABS(H_ho_1D_dense_theo_14_17 - H_ho_1D_dense_14_17%Dense_val_R) > threshold_H)) THEN
    WRITE(out_unit,*) 
    WRITE(out_unit,*) 'H_ho_1D_dense_14_17 failed to initialize'
    error_H = 1
  END IF

  !H_ho_1D_diag_1_6%Diag_val_R(1) = 1
  !WRITE(out_unit,*) 
  !WRITE(out_unit,*) "H_{Spurious 4 (diagonal, w = 1, Nb = 6)}"
  !DO i = 1, 17
  !  WRITE(out_unit,*) H_ho_1D_diag_1_6%Dense_val_R(i,:)
  !END DO

  DO i = 1, 6
    IF (ABS(H_ho_1D_diag_1_6%Diag_val_R(i) - H_ho_1D_dense_1_6%Dense_val_R(i,i)) > threshold_H) THEN
      WRITE(out_unit,*) 
      WRITE(out_unit,*) "|", H_ho_1D_diag_1_6%Diag_val_R(i), "-", H_ho_1D_dense_1_6%Dense_val_R(i,i), "| >", threshold_H
      WRITE(out_unit,*) 'H_ho_1D_diag_1_6 failed to initialize'
      error_H = 1
    END IF

    IF (ABS(H_ho_1D_diag_1_17%Diag_val_R(i) - H_ho_1D_diag_1_6%Diag_val_R(i)) > threshold_H) THEN
      WRITE(out_unit,*) 
      WRITE(out_unit,*) "|", H_ho_1D_diag_1_17%Diag_val_R(i), "-", H_ho_1D_diag_1_6%Diag_val_R(i), "| >", threshold_H
      WRITE(out_unit,*) 'H_ho_1D_diag_1_17 failed to initialize'
      error_H = 1
    END IF

    !H_ho_1D_diag_14_6%Diag_val_R(1) = 1
    !WRITE(out_unit,*) 
    !WRITE(out_unit,*) "H_{Spurious 5 (diagonal, w = 14, Nb = 6)}"
    !DO i = 1, 17
    !  WRITE(out_unit,*) H_ho_1D_dense_14_6%Dense_val_R(i,:)
    !END DO

    IF (ABS(H_ho_1D_diag_14_6%Diag_val_R(i) - H_ho_1D_diag_14_17%Diag_val_R(i)) > threshold_H) THEN
      WRITE(out_unit,*) 
      WRITE(out_unit,*) "|", H_ho_1D_diag_14_6%Diag_val_R(i), "-", H_ho_1D_diag_14_17%Diag_val_R(i), "| >", threshold_H
      WRITE(out_unit,*) 'H_ho_1D_diag_14_6 failed to initialize'
      error_H = 1
    END IF
  END DO


  !-------------------------------comparisons x-------------------------------
  error_x = 0
  threshold_x = 1E-08_Rkind

!  TYPE(MC_operator_1D_t)        :: x_ho_1D_band_14_6_1                         ! 14 = w (eigenpulsation); 6 = Nb (number of Ho basis vectors); 1 = m (mass associated to the HO)
!  TYPE(MC_operator_1D_t)        :: x_ho_1D_band_1_6_7

  !x_ho_1D_band_1_17_1%Band_val_R(1,1) = 5

  IF (ANY(ABS(x_ho_1D_band_1_17_1%Band_val_R - x_ho_1D_band_theo_1_17_1) > threshold_x)) THEN
    WRITE(out_unit,*) 
    WRITE(out_unit,*) 'x_ho_1D_band_1_17_1 failed to initialize'
    error_x = 1
  END IF

  !x_ho_1D_dense_1_6_1%Dense_val_R(1,1) = 5
  !WRITE(out_unit,*)
  !WRITE(out_unit,*) "x_{Spurious 2 (dense, w = 1, Nb = 6, m = 1)}"
  !DO i = 1, 6
  !  WRITE(out_unit,*) x_ho_1D_dense_1_6_1%Dense_val_R(i,:)
  !END DO

  IF (ANY(ABS(x_ho_1D_dense_1_6_1%Dense_val_R - x_ho_1D_dense_theo_1_6_1) > threshold_x)) THEN
    WRITE(out_unit,*) 
    WRITE(out_unit,*) 'x_ho_1D_dense_1_6_1 failed to initialize'
    error_x = 1
  END IF

  !x_ho_1D_dense_14_17_7%Dense_val_R(1,1) = 5
  !WRITE(out_unit,*)
  !WRITE(out_unit,*) "x_{Spurious 3 (dense, w = 14, Nb = 17, m = 7)}"
  !DO i = 1, 17
  !  WRITE(out_unit,*) x_ho_1D_dense_14_17_7%Dense_val_R(i,:)
  !END DO

  IF (ANY(ABS(x_ho_1D_dense_14_17_7%Dense_val_R - x_ho_1D_dense_theo_14_17_7) > threshold_x)) THEN
    WRITE(out_unit,*) 
    WRITE(out_unit,*) 'x_ho_1D_dense_14_17_7 failed to initialize'
    error_x = 1
  END IF

  !x_ho_1D_band_1_6_1%Band_val_R(6,1) = 5
  !WRITE(out_unit,*) "x_{Spurious 4 (Band, w = 1, Nb = 6, m = 1)}"
  !DO i = 1, 6
  !  WRITE(out_unit,*) x_ho_1D_band_1_6_1%Band_val_R(i,:)
  !END DO

  DO i = 1, 16                                                                  ! stops at 16 because "x_ho_1D_dense_theo_14_17_7(i+1,i)"
    IF (i <= 5) THEN
      IF (ANY(ABS(x_ho_1D_band_1_6_1%Band_val_R(i,:) - x_ho_1D_band_1_17_1%Band_val_R(i,:)) > threshold_x)) THEN ! stops at 5 because x_ho_1D_band_1_6_1%Band_val_R(6,1) is necessary 0 due to the end of the matrix. x_ho_1D_band_1_6_1%Band_val_R(6,3) Is tested at the end
        WRITE(out_unit,*) 
        WRITE(out_unit,*) "|", x_ho_1D_band_1_6_1%Band_val_R(i,:), "-", x_ho_1D_band_1_17_1%Band_val_R(i,:), "| >", threshold_x
        WRITE(out_unit,*) 'x_ho_1D_band_1_6_1 failed to initialize'
        error_x = 1
      END IF
    END IF

    !x_ho_1D_band_14_17_7%Band_val_R(1,3) = 5

    IF (ABS(x_ho_1D_band_14_17_7%Band_val_R(i,1) - x_ho_1D_dense_theo_14_17_7(i+1,i)) > threshold_x) THEN
      WRITE(out_unit,*) 
      WRITE(out_unit,*) "|", x_ho_1D_band_14_17_7%Band_val_R(i,1), "-", x_ho_1D_dense_theo_14_17_7(i+1,i), "| >", threshold_x
      WRITE(out_unit,*) 'x_ho_1D_band_14_17_7 failed to initialize'
      error_x = 1
    ELSE IF (ABS(x_ho_1D_band_14_17_7%Band_val_R(i+1,3) - x_ho_1D_dense_theo_14_17_7(i,i+1)) > threshold_x) THEN
      WRITE(out_unit,*) 
      WRITE(out_unit,*) "|", x_ho_1D_band_14_17_7%Band_val_R(i,1), "-", x_ho_1D_dense_theo_14_17_7(i+1,i), "| >", threshold_x
      WRITE(out_unit,*) 'x_ho_1D_band_14_17_7 failed to initialize'
      error_x = 1
    END IF
  END DO

  IF (ABS(x_ho_1D_band_1_6_1%Band_val_R(6,3) - x_ho_1D_band_1_17_1%Band_val_R(6,3)) > threshold_x) THEN
    WRITE(out_unit,*) 
    WRITE(out_unit,*) "|", x_ho_1D_band_1_6_1%Band_val_R(6,3), "-", x_ho_1D_band_1_17_1%Band_val_R(6,3), "| >", threshold_x
    WRITE(out_unit,*) 'x_ho_1D_band_1_6_1 failed to initialize'
    error_x = 1
  ELSE IF (ABS(x_ho_1D_band_1_6_1%Band_val_R(6,1)) > threshold_x) THEN
    WRITE(out_unit,*) 
    WRITE(out_unit,*) "|", x_ho_1D_band_1_6_1%Band_val_R(6,1), "| >", threshold_x
    WRITE(out_unit,*) 'x_ho_1D_band_1_6_1 failed to initialize'
    error_x = 1  
  END IF

  IF (ABS(x_ho_1D_band_14_17_7%Band_val_R(17,1)) > threshold_x) THEN
    WRITE(out_unit,*) 
    WRITE(out_unit,*) "|", x_ho_1D_band_14_17_7%Band_val_R(17,1), "| >", threshold_x
    WRITE(out_unit,*) 'x_ho_1D_band_14_17_7 failed to initialize'
    error_x = 1
  ELSE IF (ABS(x_ho_1D_band_14_17_7%Band_val_R(1,3)) > threshold_x) THEN
    WRITE(out_unit,*) 
    WRITE(out_unit,*) "|", x_ho_1D_band_14_17_7%Band_val_R(1,3), "| >", threshold_x
    WRITE(out_unit,*) 'x_ho_1D_band_14_17_7 failed to initialize'
    error_x = 1
  END IF

  !x_ho_1D_band_14_6_1%Band_val_R(6,1) = 5
  !WRITE(out_unit,*) "x_{Spurious 4 (Band, w = 1, Nb = 6, m = 1)}"
  !DO i = 1, 6
  !  WRITE(out_unit,*) x_ho_1D_band_1_6_1%Band_val_R(i,:)
  !END DO

  IF (ANY(ABS(x_ho_1D_band_14_6_1%Band_val_R - (x_ho_1D_band_1_6_1%Band_val_R)/(SQRT(14.0_Rkind))) > threshold_x)) THEN
    WRITE(out_unit,*) 
    WRITE(out_unit,*) 'x_ho_1D_band_14_6_1 failed to initialize'
    error_x = 1
  END IF

  !x_ho_1D_band_1_6_7%Band_val_R(6,1) = 5
  !WRITE(out_unit,*) "x_{Spurious 4 (Band, w = 1, Nb = 6, m = 1)}"
  !DO i = 1, 6
  !  WRITE(out_unit,*) x_ho_1D_band_1_6_1%Band_val_R(i,:)
  !END DO

  IF (ANY(ABS(x_ho_1D_band_1_6_7%Band_val_R - (x_ho_1D_band_1_6_1%Band_val_R)/(SQRT(SEVEN))) > threshold_x)) THEN
    WRITE(out_unit,*) 
    WRITE(out_unit,*) 'x_ho_1D_band_1_6_7 failed to initialize'
    error_x = 1
  END IF

  !-------------------------------comparisons N-------------------------------
  error_N = 0
  threshold_N = 1E-08_Rkind
  !TYPE(MC_operator_1D_t)        :: N_ho_1D_diag_6

  !N_ho_1D_dense_17%Dense_val_R(1,1) = 5

  IF (ANY(ABS(N_ho_1D_dense_17%Dense_val_R - N_ho_1D_dense_theo_17) > threshold_x)) THEN
    WRITE(out_unit,*) 
    WRITE(out_unit,*) 'N_ho_1D_dense_17 failed to initialize'
    error_N = 1
  END IF

  !N_ho_1D_diag_17%Diag_val_R(1) = 5
  !N_ho_1D_dense_6%Dense_val_R(1,1) = 5
  !N_ho_1D_diag_6%Diag_val_R(1) = 5

  DO i = 1, 17
    IF (ABS(N_ho_1D_diag_17%Diag_val_R(i) - N_ho_1D_dense_theo_17(i,i)) > threshold_x) THEN
      WRITE(out_unit,*) 
      WRITE(out_unit,*) "|", N_ho_1D_diag_17%Diag_val_R(i), "-", N_ho_1D_dense_theo_17(i,i), "| >", threshold_x
      WRITE(out_unit,*) 'N_ho_1D_diag_17 failed to initialize'
      error_N = 1
    END IF

    IF (i <= 6) THEN
      IF (ABS(N_ho_1D_dense_6%Dense_val_R(i,i) - N_ho_1D_dense_theo_17(i,i)) > threshold_x) THEN
        WRITE(out_unit,*) 
        WRITE(out_unit,*) "|", N_ho_1D_dense_6%Dense_val_R(i,i), "-", N_ho_1D_dense_theo_17(i,i), "| >", threshold_x
        WRITE(out_unit,*) 'N_ho_1D_dense_6 failed to initialize'
        error_N = 1
      END IF

      IF (ABS(N_ho_1D_diag_6%Diag_val_R(i) - N_ho_1D_dense_theo_17(i,i)) > threshold_x) THEN
        WRITE(out_unit,*) 
        WRITE(out_unit,*) "|", N_ho_1D_diag_6%Diag_val_R(i), "-", N_ho_1D_dense_theo_17(i,i), "| >", threshold_x
        WRITE(out_unit,*) 'N_ho_1D_diag_6 failed to initialize'
        error_N = 1
      END IF
    END IF
  END DO


  !-----------------------------------sum up-----------------------------------
  !error_x = 1
  IF (error_H == 0) THEN
    WRITE(out_unit,*) 
    WRITE(out_unit,*) 'The 1D Hamiltonian initializes successfully !'
  ELSE IF (error_H == 1) THEN
    WRITE(out_unit,*) 
    WRITE(out_unit,*) 'The 1D Hamiltonian did not initialize successfully...'
  ELSE
    WRITE(out_unit,*) 
    WRITE(out_unit,*) 'The file did not execute as it should have !'
  END IF

  IF (error_x == 0) THEN
    WRITE(out_unit,*) 'The 1D position operator initializes successfully !'
  ELSE IF (error_x == 1) THEN
    WRITE(out_unit,*) 'The 1D position operator did not initialize successfully...'
  ELSE
    WRITE(out_unit,*) 'The file did not execute as it should have !'
  END IF

  IF (error_N == 0) THEN
    WRITE(out_unit,*) 'The 1D number of photons operator initializes successfully !'
  ELSE IF (error_N == 1) THEN
    WRITE(out_unit,*) 'The 1D number of photons operator did not initialize successfully...'
  ELSE
    WRITE(out_unit,*) 'The file did not execute as it should have !'
  END IF

  IF (error_H + error_x + error_N == 0) THEN
    WRITE(out_unit,*) 
    WRITE(out_unit,*) 'Test 2 checked ! The 1D Hamiltonian, the 1D position operator and &
                                       &the number of photons operator initialize successfully !'
  ELSE IF (error_H + error_x + error_N >= 1 .AND. error_H + error_x <= 3) THEN
    WRITE(out_unit,*) 
    WRITE(out_unit,*) 'Test 2 failed ! At least one operator did not initialize successfully, cf. the log file for more details...'
  ELSE
    WRITE(out_unit,*) 
    WRITE(out_unit,*) 'Test 2 stopped ! The file did not execute as it should have !'
  END IF


END PROGRAM

