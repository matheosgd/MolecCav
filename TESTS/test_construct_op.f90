PROGRAM test_construct_op
  USE QDUtil_m
  USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
  USE MC_operator_1D_m
  IMPLICIT NONE

!--------------------------------First cavity mode-----------------------------
  TYPE(MC_operator_1D_t)         :: H_ho_1D_diag_1_6                           ! matrix of the one-dimensional harmonic Hamiltonian associated with HO D. The first number = eigenpulsation of the HO, the second one = the number of basis functions to expand the matrix in.
  TYPE(MC_operator_1D_t)         :: H_ho_1D_diag_14_6
  TYPE(MC_operator_1D_t)         :: H_ho_1D_diag_1_17
  TYPE(MC_operator_1D_t)         :: H_ho_1D_diag_14_17
  TYPE(MC_operator_1D_t)         :: H_ho_1D_dense_1_6
  TYPE(MC_operator_1D_t)         :: H_ho_1D_dense_14_17
  real(kind=real64), allocatable :: H_ho_1D_diag_theo_14_17(:)
  real(kind=real64), allocatable :: H_ho_1D_dense_theo_1_6(:,:)
  real(kind=real64), allocatable :: H_ho_1D_dense_theo_14_17(:,:)
  TYPE(MC_operator_1D_t)         :: x_ho_cavity_mode_1
  TYPE(MC_operator_1D_t)         :: N_ho_cavity_mode_1

  integer                        :: i, error_H


  !-------------------------construct matricies to test------------------------
  CALL MolecCav_Construct_Operator(Operator=H_ho_1D_diag_1_6, &
                                 & operator_type="Hamiltonian", &
                                 & scalar_space="Real", &
                                 & matrix_shape_type="Opt", &                  ! opt => get analytical shape. non_opt => get dense shape
                                 & Nb=6, &
                                 & w=1.0_real64, &
                                 & m=1.0_real64)

  CALL MolecCav_Construct_Operator(Operator=H_ho_1D_diag_14_6, &
                                 & operator_type="Hamiltonian", &
                                 & scalar_space="Real", &
                                 & matrix_shape_type="Opt", &                  ! opt => get analytical shape. non_opt => get dense shape
                                 & Nb=6, &
                                 & w=14.0_real64, &
                                 & m=1.0_real64)

  CALL MolecCav_Construct_Operator(Operator=H_ho_1D_diag_1_17, &
                                 & operator_type="Hamiltonian", &
                                 & scalar_space="Real", &
                                 & matrix_shape_type="Opt", &                  ! opt => get analytical shape. non_opt => get dense shape
                                 & Nb=17, &
                                 & w=1.0_real64, &
                                 & m=1.0_real64)

  CALL MolecCav_Construct_Operator(Operator=H_ho_1D_diag_14_17, &
                                 & operator_type="Hamiltonian", &
                                 & scalar_space="Real", &
                                 & matrix_shape_type="Opt", &                  ! opt => get analytical shape. non_opt => get dense shape
                                 & Nb=17, &
                                 & w=14.0_real64, &
                                 & m=1.0_real64)

  CALL MolecCav_Construct_Operator(Operator=H_ho_1D_dense_1_6, &
                                 & operator_type="Hamiltonian", &
                                 & scalar_space="Real", &
                                 & matrix_shape_type="Non_opt", &                  ! opt => get analytical shape. non_opt => get dense shape
                                 & Nb=6, &
                                 & w=1.0_real64, &
                                 & m=1.0_real64)
                                 
  CALL MolecCav_Construct_Operator(Operator=H_ho_1D_dense_14_17, &
                                 & operator_type="Hamiltonian", &
                                 & scalar_space="Real", &
                                 & matrix_shape_type="Non_opt", &                  ! opt => get analytical shape. non_opt => get dense shape
                                 & Nb=17, &
                                 & w=14.0_real64, &
                                 & m=1.0_real64)


  !------------------------construct reference matricies-----------------------
  ALLOCATE(H_ho_1D_diag_theo_14_17(17))
  ALLOCATE(H_ho_1D_dense_theo_1_6(6,6))
  ALLOCATE(H_ho_1D_dense_theo_14_17(17,17))

  H_ho_1D_diag_theo_14_17  = 0.0_real64
  H_ho_1D_dense_theo_1_6   = 0.0_real64
  H_ho_1D_dense_theo_14_17 = 0.0_real64

  DO i = 1, 17
    H_ho_1D_diag_theo_14_17(i)    = 14*(i - 1 + 0.5_real64)
    H_ho_1D_dense_theo_14_17(i,i) = 14*(i - 1 + 0.5_real64)
    IF (i <= 6) THEN
      H_ho_1D_dense_theo_1_6(i,i) = (i - 1 + 0.5_real64)
    END IF
  END DO

  WRITE(OUTPUT_UNIT,*) ''
  WRITE(OUTPUT_UNIT,*) "H_{Reference 1 (diagonal, w = 14, Nb = 17)}"
  DO i = 1, 17
    WRITE(OUTPUT_UNIT,*) H_ho_1D_diag_theo_14_17(i)
  END DO

  WRITE(OUTPUT_UNIT,*) ''
  WRITE(OUTPUT_UNIT,*) "H_{Reference 2 (dense, w = 1, Nb = 6)}"
  DO i = 1, 6
    WRITE(OUTPUT_UNIT,*) H_ho_1D_dense_theo_1_6(i,:)
  END DO

  WRITE(OUTPUT_UNIT,*) ''
  WRITE(OUTPUT_UNIT,*) "H_{Reference 3 (dense, w = 14, Nb = 17)}"
  DO i = 1, 17
    WRITE(OUTPUT_UNIT,*) H_ho_1D_dense_theo_14_17(i,:)
  END DO


  !---------------------------------comparisons---------------------------------
  IF (ALL(H_ho_1D_diag_theo_14_17 /= H_ho_1D_diag_14_17%Diag_val_R)) THEN
    WRITE(OUTPUT_UNIT,*) ''
    WRITE(OUTPUT_UNIT,*) 'H_ho_1D_diag_14_17 failed to initialize'
    error_H = 1
  END IF

  IF (ALL(H_ho_1D_dense_theo_1_6 /= H_ho_1D_dense_1_6%Dense_val_R)) THEN
    WRITE(OUTPUT_UNIT,*) ''
    WRITE(OUTPUT_UNIT,*) 'H_ho_1D_dense_1_6 failed to initialize'
    error_H = 1
  END IF

  IF (ALL(H_ho_1D_dense_theo_14_17 /= H_ho_1D_dense_14_17%Dense_val_R)) THEN
    WRITE(OUTPUT_UNIT,*) ''
    WRITE(OUTPUT_UNIT,*) 'H_ho_1D_dense_14_17 failed to initialize'
    error_H = 1
  END IF

  DO i = 1, 6
    IF (H_ho_1D_diag_1_6%Diag_val_R(i) /= H_ho_1D_dense_1_6%Dense_val_R(i,i)) THEN
      WRITE(OUTPUT_UNIT,*) ''
      WRITE(OUTPUT_UNIT,*) H_ho_1D_diag_1_6%Diag_val_R(i), "/=", H_ho_1D_dense_1_6%Dense_val_R(i,i)
      WRITE(OUTPUT_UNIT,*) 'H_ho_1D_diag_1_6 failed to initialize'
      error_H = 1
    END IF

    IF (H_ho_1D_diag_1_17%Diag_val_R(i) /= H_ho_1D_diag_1_6%Diag_val_R(i)) THEN
      WRITE(OUTPUT_UNIT,*) ''
      WRITE(OUTPUT_UNIT,*) H_ho_1D_diag_1_17%Diag_val_R(i), "/=", H_ho_1D_diag_1_6%Diag_val_R(i)
      WRITE(OUTPUT_UNIT,*) 'H_ho_1D_diag_1_17 failed to initialize'
      error_H = 1
    END IF

    IF (H_ho_1D_diag_14_6%Diag_val_R(i) /= H_ho_1D_diag_14_17%Diag_val_R(i)) THEN
      WRITE(OUTPUT_UNIT,*) ''
      WRITE(OUTPUT_UNIT,*) H_ho_1D_diag_14_6%Diag_val_R(i), "/=", H_ho_1D_diag_14_17%Diag_val_R(i)
      WRITE(OUTPUT_UNIT,*) 'H_ho_1D_diag_14_6 failed to initialize'
      error_H = 1
    END IF
  END DO

  !---------------------------------sum up---------------------------------
  IF (error_H == 0) THEN
    WRITE(OUTPUT_UNIT,*) ''
    WRITE(OUTPUT_UNIT,*) 'Test 2 checked ! The 1D Hamiltonian initializes successfully !'
  ELSE IF (error_H == 1) THEN
    WRITE(OUTPUT_UNIT,*) ''
    WRITE(OUTPUT_UNIT,*) 'Test 2 failed ! The 1D Hamiltonian did not initialize successfully...'
  ELSE
    WRITE(OUTPUT_UNIT,*) ''
    WRITE(OUTPUT_UNIT,*) 'Test 2 stopped ! The file did not execute as it should have !'
  END IF


END PROGRAM
