MODULE MC_algebra_m
  !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
  USE QDUtil_m
  IMPLICIT NONE

  
  CONTAINS
  
    
  SUBROUTINE MolecCav_Normalize_2D_real(Psi)
    USE QDUtil_m
    IMPLICIT NONE

    real(kind=Rkind), intent(inout)    :: Psi(:,:)                             ! already allocated

    real(kind=Rkind)                   :: Norm
    logical, parameter                 :: debug = .FALSE.

    Norm = ZERO
    CALL MolecCav_Norm_2D_real(Norm, Psi)
    IF (debug) WRITE(out_unit, *) "Computed norm of 2D WF = ", Norm

    IF (Norm < 1E-08_Rkind) THEN
      STOP "Attempt to normalize matrix of norm ZERO"
    ELSE 
      Psi(:,:) = Psi(:,:) / Norm
    END IF 

    IF (debug) THEN
      CALL MolecCav_Norm_2D_real(Norm, Psi)
      WRITE(out_unit, *) "Computed new norm of 2D WF = ", Norm
    END IF 


  END SUBROUTINE

  
  SUBROUTINE MolecCav_Normalize_1D_real(Psi)
    USE QDUtil_m
    IMPLICIT NONE

    real(kind=Rkind), intent(inout)    :: Psi(:)                               ! already allocated

    real(kind=Rkind)                   :: Norm
    logical, parameter                 :: debug = .FALSE.

    Norm = ZERO
    CALL MolecCav_Norm_1D_real(Norm, Psi)
    IF (debug) WRITE(out_unit, *) "Computed norm of 1D WF = ", Norm

    IF (Norm < 1E-08_Rkind) THEN
      STOP "Attempt to normalize matrix of norm ZERO"
    ELSE 
      Psi(:) = Psi(:) / Norm
    END IF 

    IF (debug) THEN
      CALL MolecCav_Norm_1D_real(Norm, Psi)
      WRITE(out_unit, *) "Computed new norm of 1D WF = ", Norm
    END IF 


  END SUBROUTINE


  SUBROUTINE MolecCav_Norm_2D_real(Norm, Psi)
    USE QDUtil_m
    IMPLICIT NONE
  
    real(kind=Rkind), intent(inout) :: Norm
    real(kind=Rkind), intent(in)    :: Psi(:,:)                                ! already allocated

    real(kind=Rkind), allocatable   :: A(:,:)
    integer                         :: dim, i
    
    CALL MolecCav_scalar_product_2D_real(Norm, Psi, Psi)
    Norm = SQRT(Norm)
  
  END SUBROUTINE


  SUBROUTINE MolecCav_Norm_1D_real(Norm, Psi)
    USE QDUtil_m
    IMPLICIT NONE
  
    real(kind=Rkind), intent(inout) :: Norm
    real(kind=Rkind), intent(in)    :: Psi(:)                                  ! already allocated
    
    Norm = SQRT(DOT_PRODUCT(Psi, Psi))
  
  END SUBROUTINE


  SUBROUTINE MolecCav_scalar_product_2D_real(scalar_product, Psi_1, Psi_2)
    USE QDUtil_m
    IMPLICIT NONE
  
    real(kind=Rkind),         intent(inout) :: scalar_product
    real(kind=Rkind), target, intent(in)    :: Psi_1(:,:)                                ! already allocated + "target" means that it can be pointed by a pointer
    real(kind=Rkind), target, intent(in)    :: Psi_2(:,:)                                ! already allocated + "target" means that it can be pointed by a pointer

    real(kind=Rkind), pointer               :: V1(:), V2(:)                              ! <=> real(kind=Rkind), dimension(:), pointer :: V1, V2
    real(kind=Rkind), target, allocatable   :: A(:), B(:)     ! obligé à cause "Rank remapping target must be rank 1 or simply contiguous" error : à corriger
    integer                                 :: dim
    logical, parameter                      :: debug = .TRUE.
    
    dim = Size(Psi_1, dim=2)
    IF (dim /= Size(Psi_2, dim=2) .OR. Size(Psi_2, dim=1) /= Size(Psi_2, dim=1)) THEN
      STOP "The matrices are expected to have same dimensions for the scalar product."
    END IF

    scalar_product = ZERO

    ALLOCATE(A(dim*Size(Psi_1, dim=1)))
    ALLOCATE(B(dim*Size(Psi_1, dim=1)))
    A(:) = reshape(Psi_1, shape=[dim*Size(Psi_1, dim=1)])
    B(:) = reshape(Psi_2, shape=[dim*Size(Psi_1, dim=1)])
    
    V1(1:Size(Psi_1)) => A(:) ! Psi_1(:,:) ideally
    V2(1:Size(Psi_2)) => B(:) ! Psi_2(:,:) //
    scalar_product = DOT_PRODUCT(V1, V2)

    NULLIFY(V1)
    NULLIFY(V2)

    IF (debug) WRITE(out_unit, *) "Computed scalar product : < Psi_1 |  Psi_2 >  =", scalar_product, "supposed to get 79 the 10/02/2025"

  END SUBROUTINE


END MODULE
  