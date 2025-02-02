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

  
  SUBROUTINE MolecCav_Norm_WF_2D(Norm, Psi)
    USE QDUtil_m
    IMPLICIT NONE
  
    real(kind=Rkind), intent(inout) :: Norm
    real(kind=Rkind), intent(in)    :: Psi(:,:)                                ! already allocated

    real(kind=Rkind), allocatable   :: A(:,:)
    integer                         :: dim, i
    
    dim = Size(Psi,2)
    ALLOCATE(A(dim,dim))
    Norm = ZERO

    A = MATMUL(TRANSPOSE(Psi), Psi)
    DO i=1, dim
      Norm = Norm + A(i,i)                                                     ! compute the trace
    END DO
    Norm = SQRT(Norm)
  
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


  SUBROUTINE MolecCav_Norm_WF_1D(Norm, Psi)
    USE QDUtil_m
    IMPLICIT NONE
  
    real(kind=Rkind), intent(inout) :: Norm
    real(kind=Rkind), intent(in)    :: Psi(:)                                  ! already allocated
    
    Norm = SQRT(DOT_PRODUCT(Psi, Psi))
  
  END SUBROUTINE


END MODULE
  