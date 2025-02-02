MODULE MC_algebra_m
  !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
  USE QDUtil_m
  IMPLICIT NONE

  
  CONTAINS
  
    
  SUBROUTINE MolecCav_Normalize_WF_2D(Psi)
    USE QDUtil_m
    IMPLICIT NONE

    real(kind=Rkind), intent(inout)    :: Psi(:,:)                             ! already allocated

    real(kind=Rkind)                   :: Norm

    Norm = ZERO
    CALL MolecCav_Norm_WF_2D(Norm, Psi)
    WRITE(out_unit, *) "Computed norm of 2D WF = ", Norm

    IF (Norm < 1E-08_Rkind) THEN
      STOP "Attempt to normalize matrix of norm ZERO"
    ELSE 
      Psi(:,:) = Psi(:,:) / Norm
    END IF 

    CALL MolecCav_Norm_WF_2D(Norm, Psi)
    WRITE(out_unit, *) "Computed new norm of 2D WF = ", Norm

  END SUBROUTINE

  
  SUBROUTINE MolecCav_Normalize_WF_1D(Psi)
    USE QDUtil_m
    IMPLICIT NONE

    real(kind=Rkind), intent(inout)    :: Psi(:)                               ! already allocated

    real(kind=Rkind)                   :: Norm

    Norm = ZERO
    CALL MolecCav_Norm_WF_1D(Norm, Psi)
    WRITE(out_unit, *) "Computed norm of 1D WF = ", Norm

    IF (Norm < 1E-08_Rkind) THEN
      STOP "Attempt to normalize matrix of norm ZERO"
    ELSE 
      Psi(:) = Psi(:) / Norm
    END IF 

    CALL MolecCav_Norm_WF_1D(Norm, Psi)
    WRITE(out_unit, *) "Computed new norm of 1D WF = ", Norm

  END SUBROUTINE


  SUBROUTINE MolecCav_Norm_WF_2D(Norm, Psi)
    USE QDUtil_m
    IMPLICIT NONE
  
    real(kind=Rkind), intent(inout) :: Norm
    real(kind=Rkind), intent(in)    :: Psi(:,:)                                ! already allocated

    real(kind=Rkind), allocatable   :: A(:,:)
    integer                         :: dim, i
    
    CALL MolecCav_scalar_product_2D(Norm, Psi, Psi)
    Norm = SQRT(Norm)
  
  END SUBROUTINE


  SUBROUTINE MolecCav_Norm_WF_1D(Norm, Psi)
    USE QDUtil_m
    IMPLICIT NONE
  
    real(kind=Rkind), intent(inout) :: Norm
    real(kind=Rkind), intent(in)    :: Psi(:)                                  ! already allocated
    
    Norm = SQRT(DOT_PRODUCT(Psi, Psi))
  
  END SUBROUTINE


  SUBROUTINE MolecCav_scalar_product_2D(Projection, Psi_1, Psi_2)
    USE QDUtil_m
    IMPLICIT NONE
  
    real(kind=Rkind), intent(inout) :: Projection
    real(kind=Rkind), intent(in)    :: Psi_1(:,:)                                ! already allocated
    real(kind=Rkind), intent(in)    :: Psi_2(:,:)                                ! already allocated

    real(kind=Rkind), allocatable   :: A(:,:)
    integer                         :: dim, i
    
    dim = Size(Psi_1,2)
    IF (dim /= Size(Psi_2,2) .OR. Size(Psi_2,1) /= Size(Psi_2,1)) THEN
      STOP "The matrices are expected to have same dimensions for the projection."
    END IF

    ALLOCATE(A(dim,dim))
    Projection = ZERO

    A = MATMUL(TRANSPOSE(Psi_1), Psi_2)
    DO i=1, dim
      Projection = Projection + A(i,i)                                                     ! compute the trace
    END DO
  
  END SUBROUTINE


END MODULE
  