MODULE Algebra_m
  !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
  USE QDUtil_m
  IMPLICIT NONE

  
  PRIVATE

  PUBLIC Normalize, Scalar_product

  INTERFACE Normalize
    MODULE PROCEDURE MolecCav_Normalize_2D_real, MolecCav_Normalize_2D_complex, & 
                   & MolecCav_Normalize_1D_real, MolecCav_Normalize_1D_complex
  END INTERFACE
  INTERFACE Scalar_product
    MODULE PROCEDURE MolecCav_Scalar_product_2D_real, MolecCav_Scalar_product_2D_complex, & 
                   & MolecCav_Scalar_product_1D_real, MolecCav_Scalar_product_1D_complex
  END INTERFACE

    
  CONTAINS


  SUBROUTINE MolecCav_Normalize_2D_real(Psi)
    USE QDUtil_m
    IMPLICIT NONE

    real(kind=Rkind), intent(inout) :: Psi(:,:)                             ! already allocated

    real(kind=Rkind)                :: Norm
    logical, parameter              :: debug = .FALSE.

    Norm = ZERO
    CALL MolecCav_Norm_2D_real(Norm, Psi)
    IF (debug) WRITE(out_unit,*) "Computed norm of 2D WF = ", Norm

    IF (Norm < 1E-08_Rkind) THEN
      STOP "Attempt to normalize matrix of norm ZERO"
    ELSE 
      Psi(:,:) = Psi(:,:) / Norm
    END IF 

    IF (debug) THEN
      CALL MolecCav_Norm_2D_real(Norm, Psi)
      WRITE(out_unit,*) "Computed new norm of 2D WF = ", Norm
    END IF 


  END SUBROUTINE

  
  SUBROUTINE MolecCav_Normalize_2D_complex(Psi)
    USE QDUtil_m
    IMPLICIT NONE

    complex(kind=Rkind), intent(inout) :: Psi(:,:)                             ! already allocated

    real(kind=Rkind)                   :: Norm
    logical, parameter                 :: debug = .FALSE.

    Norm = ZERO
    CALL MolecCav_Norm_2D_complex(Norm, Psi)
    IF (debug) WRITE(out_unit,*) "Computed norm of 2D WF = ", Norm

    IF (Norm < 1E-08_Rkind) THEN
      STOP "Attempt to normalize matrix of norm ZERO"
    ELSE 
      Psi(:,:) = Psi(:,:) / Norm
    END IF 

    IF (debug) THEN
      CALL MolecCav_Norm_2D_complex(Norm, Psi)
      WRITE(out_unit,*) "Computed new norm of 2D WF = ", Norm
    END IF 


  END SUBROUTINE

  
  SUBROUTINE MolecCav_Normalize_1D_real(Psi)
    USE QDUtil_m
    IMPLICIT NONE

    real(kind=Rkind), intent(inout) :: Psi(:)                               ! already allocated

    real(kind=Rkind)                :: Norm
    logical, parameter              :: debug = .FALSE.

    Norm = ZERO
    CALL MolecCav_Norm_1D_real(Norm, Psi)
    IF (debug) WRITE(out_unit,*) "Computed norm of 1D WF = ", Norm

    IF (Norm < 1E-08_Rkind) THEN
      STOP "Attempt to normalize matrix of norm ZERO"
    ELSE 
      Psi(:) = Psi(:) / Norm
    END IF 

    IF (debug) THEN
      CALL MolecCav_Norm_1D_real(Norm, Psi)
      WRITE(out_unit,*) "Computed new norm of 1D WF = ", Norm
    END IF 


  END SUBROUTINE


  SUBROUTINE MolecCav_Normalize_1D_complex(Psi)
    USE QDUtil_m
    IMPLICIT NONE

    complex(kind=Rkind), intent(inout) :: Psi(:)                               ! already allocated

    real(kind=Rkind)                   :: Norm
    logical, parameter                 :: debug = .FALSE.

    Norm = ZERO
    CALL MolecCav_Norm_1D_complex(Norm, Psi)
    IF (debug) WRITE(out_unit,*) "Computed norm of 1D WF = ", Norm

    IF (Norm < 1E-08_Rkind) THEN
      STOP "Attempt to normalize matrix of norm ZERO"
    ELSE 
      Psi(:) = Psi(:) / Norm
    END IF 

    IF (debug) THEN
      CALL MolecCav_Norm_1D_complex(Norm, Psi)
      WRITE(out_unit,*) "Computed new norm of 1D WF = ", Norm
    END IF 


  END SUBROUTINE


  SUBROUTINE MolecCav_Norm_2D_real(Norm, Psi)
    USE QDUtil_m
    IMPLICIT NONE
  
    real(kind=Rkind), intent(inout) :: Norm
    real(kind=Rkind), intent(in)    :: Psi(:,:)                                ! already allocated
    
    CALL MolecCav_Scalar_product_2D_real(Norm, Psi, Psi)
    Norm = SQRT(Norm)
  
  END SUBROUTINE


  SUBROUTINE MolecCav_Norm_2D_complex(Norm, Psi)
    USE QDUtil_m
    IMPLICIT NONE
  
    real(kind=Rkind),    intent(inout) :: Norm
    complex(kind=Rkind), intent(in)    :: Psi(:,:)                                ! already allocated

    CALL MolecCav_Scalar_product_2D_complex(Norm, Psi, Psi)
    Norm = SQRT(Norm)
  
  END SUBROUTINE


  SUBROUTINE MolecCav_Norm_1D_real(Norm, Psi)
    USE QDUtil_m
    IMPLICIT NONE
  
    real(kind=Rkind), intent(inout) :: Norm
    real(kind=Rkind), intent(in)    :: Psi(:)                                  ! already allocated
    
    Norm = SQRT(DOT_PRODUCT(Psi, Psi))
  
  END SUBROUTINE


  SUBROUTINE MolecCav_Norm_1D_complex(Norm, Psi)
    USE QDUtil_m
    IMPLICIT NONE
  
    real(kind=Rkind),    intent(inout) :: Norm
    complex(kind=Rkind), intent(in)    :: Psi(:)                                  ! already allocated
    
    Norm = SQRT(DOT_PRODUCT(Psi, Psi))
  
  END SUBROUTINE


  SUBROUTINE MolecCav_Scalar_product_2D_real(Sca_pdt, Psi_1, Psi_2)
    USE QDUtil_m
    IMPLICIT NONE
  
    real(kind=Rkind),         intent(inout) :: Sca_pdt
    real(kind=Rkind), target, intent(in)    :: Psi_1(:,:)                                ! already allocated + "target" means that it can be pointed by a pointer
    real(kind=Rkind), target, intent(in)    :: Psi_2(:,:)                                ! already allocated + "target" means that it can be pointed by a pointer

    real(kind=Rkind), pointer               :: V1(:), V2(:)                              ! <=> real(kind=Rkind), Dimension(:), pointer :: V1, V2
    real(kind=Rkind), target, allocatable   :: A(:), B(:)     ! obligé à cause "Rank remapping target must be rank 1 or simply contiguous" error : à corriger
    integer                                 :: Dim
    logical, parameter                      :: debug = .FALSE.
    
    Dim = Size(Psi_1, Dim=2)
    IF (Dim /= Size(Psi_2, Dim=2) .OR. Size(Psi_2, Dim=1) /= Size(Psi_2, Dim=1)) THEN
      STOP "The matrices are expected to have same Dimensions for the scalar product."
    END IF

    Sca_pdt = ZERO

    ALLOCATE(A(Dim*Size(Psi_1, Dim=1)))
    ALLOCATE(B(Dim*Size(Psi_1, Dim=1)))
    A(:) = reshape(Psi_1, shape=[Dim*Size(Psi_1, Dim=1)])
    B(:) = reshape(Psi_2, shape=[Dim*Size(Psi_1, Dim=1)])
    
    V1(1:Size(Psi_1)) => A(:) ! Psi_1(:,:) ideally
    V2(1:Size(Psi_2)) => B(:) ! Psi_2(:,:) //
    Sca_pdt = DOT_PRODUCT(V1, V2)

    NULLIFY(V1)
    NULLIFY(V2)

    IF (debug) WRITE(out_unit,*) "Computed scalar product : < Psi_1 |  Psi_2 >  =", Sca_pdt, &
                               & "supposed to get 79 the 10/02/2025"

  END SUBROUTINE


  SUBROUTINE MolecCav_Scalar_product_2D_complex(Sca_pdt, Psi_1, Psi_2)
    USE QDUtil_m
    IMPLICIT NONE
  
    real(kind=Rkind),            intent(inout) :: Sca_pdt
    complex(kind=Rkind), target, intent(in)    :: Psi_1(:,:)                                ! already allocated + "target" means that it can be pointed by a pointer
    complex(kind=Rkind), target, intent(in)    :: Psi_2(:,:)                                ! already allocated + "target" means that it can be pointed by a pointer

    complex(kind=Rkind), pointer               :: V1(:), V2(:)                              ! <=> real(kind=Rkind), Dimension(:), pointer :: V1, V2
    complex(kind=Rkind), target, allocatable   :: A(:), B(:)     ! obligé à cause "Rank remapping target must be rank 1 or simply contiguous" error : à corriger
    integer                                    :: Dim
    logical, parameter                         :: debug = .FALSE.
    
    Dim = Size(Psi_1, Dim=2)
    IF (Dim /= Size(Psi_2, Dim=2) .OR. Size(Psi_2, Dim=1) /= Size(Psi_2, Dim=1)) THEN
      STOP "The matrices are expected to have same Dimensions for the scalar product."
    END IF

    Sca_pdt = ZERO

    ALLOCATE(A(Dim*Size(Psi_1, Dim=1)))
    ALLOCATE(B(Dim*Size(Psi_1, Dim=1)))
    A(:) = reshape(Psi_1, shape=[Dim*Size(Psi_1, Dim=1)])
    B(:) = reshape(Psi_2, shape=[Dim*Size(Psi_1, Dim=1)])
    
    V1(1:Size(Psi_1)) => A(:) ! Psi_1(:,:) ideally
    V2(1:Size(Psi_2)) => B(:) ! Psi_2(:,:) //
    Sca_pdt = DOT_PRODUCT(V1, V2)

    NULLIFY(V1)
    NULLIFY(V2)

    IF (debug) WRITE(out_unit,*) "Computed scalar product : < Psi_1 |  Psi_2 >  =", Sca_pdt, &
                               & "supposed to get 79 the 10/02/2025"

  END SUBROUTINE


  SUBROUTINE MolecCav_Scalar_product_1D_real(Sca_pdt, Psi_1, Psi_2)
    USE QDUtil_m
    IMPLICIT NONE
  
    real(kind=Rkind), intent(inout) :: Sca_pdt
    real(kind=Rkind), intent(in)    :: Psi_1(:)                                ! already allocated + "target" means that it can be pointed by a pointer
    real(kind=Rkind), intent(in)    :: Psi_2(:)                                ! already allocated + "target" means that it can be pointed by a pointer

    integer                         :: Dim
    logical, parameter              :: debug = .FALSE.
    
    Dim = Size(Psi_1)
    IF (Dim /= Size(Psi_2)) THEN
      STOP "The matrices are expected to have same Dimensions for the scalar product."
    END IF

    Sca_pdt = ZERO

    Sca_pdt = DOT_PRODUCT(Psi_1, Psi_2)

    IF (debug) WRITE(out_unit,*) "Computed scalar product : < Psi_1 |  Psi_2 >  =", Sca_pdt, &
                               & "supposed to get 79 the 10/02/2025"

  END SUBROUTINE


  SUBROUTINE MolecCav_Scalar_product_1D_complex(Sca_pdt, Psi_1, Psi_2)
    USE QDUtil_m
    IMPLICIT NONE
  
    real(kind=Rkind),    intent(inout) :: Sca_pdt
    complex(kind=Rkind), intent(in)    :: Psi_1(:)                                ! already allocated + "target" means that it can be pointed by a pointer
    complex(kind=Rkind), intent(in)    :: Psi_2(:)                                ! already allocated + "target" means that it can be pointed by a pointer

    integer                            :: Dim
    logical, parameter                 :: debug = .FALSE.
    
    Dim = Size(Psi_1)
    IF (Dim /= Size(Psi_2)) THEN
      STOP "The matrices are expected to have same Dimensions for the scalar product."
    END IF

    Sca_pdt = ZERO

    Sca_pdt = DOT_PRODUCT(Psi_1, Psi_2)

    IF (debug) WRITE(out_unit,*) "Computed scalar product : < Psi_1 |  Psi_2 >  =", Sca_pdt, &
                               & "supposed to get 79 the 10/02/2025"

  END SUBROUTINE


END MODULE
  