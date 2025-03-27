!==================================================================================================
!==================================================================================================
! This file is part of MolecCav.
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
! README :
! to be written soon
!==================================================================================================
!==================================================================================================
MODULE Algebra_m
  !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
  USE QDUtil_m
  IMPLICIT NONE

  
  PRIVATE

  PUBLIC Normalize, Norm_of, Scalar_product

  INTERFACE Normalize
    MODULE PROCEDURE MolecCav_Normalize_R2_real, MolecCav_Normalize_R2_complex, & 
                   & MolecCav_Normalize_R1_real, MolecCav_Normalize_R1_complex
  END INTERFACE
  INTERFACE Norm_of
    MODULE PROCEDURE MolecCav_Norm_R2_real, MolecCav_Norm_R2_complex, & 
                   & MolecCav_Norm_R1_real, MolecCav_Norm_R1_complex
  END INTERFACE
  INTERFACE Scalar_product
    MODULE PROCEDURE MolecCav_Scalar_product_R2_real, MolecCav_Scalar_product_R2_complex, & 
                   & MolecCav_Scalar_product_R1_real, MolecCav_Scalar_product_R1_complex
  END INTERFACE

  
  CONTAINS


  SUBROUTINE MolecCav_Normalize_R2_real(Psi)
    USE QDUtil_m
    IMPLICIT NONE

    real(kind=Rkind), intent(inout) :: Psi(:,:)                                                                                  ! already allocated

    real(kind=Rkind)                :: Norm, Threshold
    logical, parameter              :: Debug = .FALSE.

    Threshold = 1E-08_Rkind
    CALL MolecCav_Norm_R2_real(Norm, Psi)
    IF (Debug) WRITE(out_unit,*) "Computed norm of R2 WF = "//TO_string(Norm)

    IF (Norm < Threshold) THEN
      WRITE(out_unit,*) "Attempt to normalize matrix of norm ZERO"
      STOP "### Attempt to normalize matrix of norm ZERO"
    ELSE 
      Psi(:,:) = Psi(:,:) / Norm
    END IF 

    IF (Debug) THEN
      CALL MolecCav_Norm_R2_real(Norm, Psi)
      WRITE(out_unit,*) "Computed new norm of R2 WF = "//TO_string(Norm)
    END IF 

  END SUBROUTINE MolecCav_Normalize_R2_real

  
  SUBROUTINE MolecCav_Normalize_R2_complex(Psi)
    USE QDUtil_m
    IMPLICIT NONE

    complex(kind=Rkind), intent(inout) :: Psi(:,:)                                                                               ! already allocated

    real(kind=Rkind)                   :: Norm, Threshold
    logical, parameter                 :: Debug = .FALSE.

    Threshold = 1E-08_Rkind
    CALL MolecCav_Norm_R2_complex(Norm, Psi)
    IF (Debug) WRITE(out_unit,*) "Computed norm of R2 WF = "//TO_string(Norm)

    IF (Norm < Threshold) THEN
      WRITE(out_unit,*) "Attempt to normalize matrix of norm ZERO"
      STOP "### Attempt to normalize matrix of norm ZERO"
    ELSE 
      Psi(:,:) = Psi(:,:) / Norm
    END IF 

    IF (Debug) THEN
      CALL MolecCav_Norm_R2_complex(Norm, Psi)
      WRITE(out_unit,*) "Computed new norm of R2 WF = "//TO_string(Norm)
    END IF 

  END SUBROUTINE MolecCav_Normalize_R2_complex

  
  SUBROUTINE MolecCav_Normalize_R1_real(Psi)
    USE QDUtil_m
    IMPLICIT NONE

    real(kind=Rkind), intent(inout) :: Psi(:)                                                                                    ! already allocated

    real(kind=Rkind)                :: Norm, Threshold
    logical, parameter              :: Debug = .FALSE.

    Threshold = 1E-08_Rkind
    CALL MolecCav_Norm_R1_real(Norm, Psi)
    IF (Debug) WRITE(out_unit,*) "Computed norm of R1 WF = "//TO_string(Norm)

    IF (Norm < Threshold) THEN
      WRITE(out_unit,*) "Attempt to normalize matrix of norm ZERO"
      STOP "### Attempt to normalize matrix of norm ZERO"
    ELSE 
      Psi(:) = Psi(:) / Norm
    END IF 

    IF (Debug) THEN
      CALL MolecCav_Norm_R1_real(Norm, Psi)
      WRITE(out_unit,*) "Computed new norm of R1 WF = "//TO_string(Norm)
    END IF 

  END SUBROUTINE MolecCav_Normalize_R1_real


  SUBROUTINE MolecCav_Normalize_R1_complex(Psi)
    USE QDUtil_m
    IMPLICIT NONE

    complex(kind=Rkind), intent(inout) :: Psi(:)                                                                                 ! already allocated

    real(kind=Rkind)                   :: Norm, Threshold
    logical, parameter                 :: Debug = .FALSE.

    Threshold = 1E-08_Rkind
    CALL MolecCav_Norm_R1_complex(Norm, Psi)
    IF (Debug) WRITE(out_unit,*) "Computed norm of R1 WF = "//TO_string(Norm)

    IF (Norm < Threshold) THEN
      WRITE(out_unit,*) "Attempt to normalize matrix of norm ZERO"
      STOP "### Attempt to normalize matrix of norm ZERO"
    ELSE 
      Psi(:) = Psi(:) / Norm
    END IF 

    IF (Debug) THEN
      CALL MolecCav_Norm_R1_complex(Norm, Psi)
      WRITE(out_unit,*) "Computed new norm of R1 WF = "//TO_string(Norm)
    END IF 

  END SUBROUTINE MolecCav_Normalize_R1_complex


  SUBROUTINE MolecCav_Norm_R2_real(Norm, Psi)
    USE QDUtil_m
    IMPLICIT NONE
  
    real(kind=Rkind), intent(inout) :: Norm
    real(kind=Rkind), intent(in)    :: Psi(:,:)                                                                                  ! already allocated
    
    CALL MolecCav_Scalar_product_R2_real(Norm, Psi, Psi)
    Norm = SQRT(Norm)
  
  END SUBROUTINE MolecCav_Norm_R2_real


  SUBROUTINE MolecCav_Norm_R2_complex(Norm, Psi)
    USE QDUtil_m
    IMPLICIT NONE
  
    real(kind=Rkind),    intent(inout) :: Norm
    complex(kind=Rkind), intent(in)    :: Psi(:,:)                                                                               ! already allocated

    complex(kind=Rkind)                :: Sca_pdt

    CALL MolecCav_Scalar_product_R2_complex(Sca_pdt, Psi, Psi)
    Norm = SQRT(REAL(Sca_pdt))
  
  END SUBROUTINE MolecCav_Norm_R2_complex


  SUBROUTINE MolecCav_Norm_R1_real(Norm, Psi)
    USE QDUtil_m
    IMPLICIT NONE
  
    real(kind=Rkind), intent(inout) :: Norm
    real(kind=Rkind), intent(in)    :: Psi(:)                                                                                    ! already allocated
    
    CALL Scalar_product(Norm, Psi, Psi)
    Norm = SQRT(Norm)
  
  END SUBROUTINE MolecCav_Norm_R1_real


  SUBROUTINE MolecCav_Norm_R1_complex(Norm, Psi)
    USE QDUtil_m
    IMPLICIT NONE
  
    real(kind=Rkind),    intent(inout) :: Norm
    complex(kind=Rkind), intent(in)    :: Psi(:)                                                                                 ! already allocated

    complex(kind=Rkind)                :: Sca_pdt
    
    CALL Scalar_product(Sca_pdt, Psi, Psi)
    Norm = SQRT(REAL(Sca_pdt, kind=Rkind))
  
  END SUBROUTINE MolecCav_Norm_R1_complex


  SUBROUTINE MolecCav_Scalar_product_R2_real(Sca_pdt, Psi_1, Psi_2)
    USE QDUtil_m
    IMPLICIT NONE
  
    real(kind=Rkind), intent(inout) :: Sca_pdt
    real(kind=Rkind), intent(in)    :: Psi_1(:,:)                                                                        ! already allocated
    real(kind=Rkind), intent(in)    :: Psi_2(:,:)                                                                        ! already allocated

    integer                         :: Dim, i_2
    logical, parameter              :: Debug = .FALSE.
    
    Dim = Size(Psi_1, Dim=2)
    IF (Dim /= Size(Psi_2, Dim=2) .OR. Size(Psi_2, Dim=1) /= Size(Psi_2, Dim=1)) THEN
      WRITE(out_unit,*) "The matrices are expected to have same Dimensions for the scalar product."
      STOP "### The matrices are expected to have same Dimensions for the scalar product."
    END IF
    
    Sca_pdt = ZERO
    
    DO i_2 = 1, Dim
      Sca_pdt = Sca_pdt + DOT_PRODUCT(Psi_1(:,i_2), Psi_2(:,i_2))
    END DO

    IF (Debug) WRITE(out_unit,*) "Computed scalar product : < Psi_1 |  Psi_2 >  ="//TO_string(Sca_pdt)

  END SUBROUTINE MolecCav_Scalar_product_R2_real


  SUBROUTINE MolecCav_Scalar_product_R2_complex(Sca_pdt, Psi_1, Psi_2)
    USE QDUtil_m
    IMPLICIT NONE
  
    complex(kind=Rkind), intent(inout) :: Sca_pdt
    complex(kind=Rkind), intent(in)    :: Psi_1(:,:)                                                                     ! already allocated
    complex(kind=Rkind), intent(in)    :: Psi_2(:,:)                                                                     ! already allocated

    integer                            :: Dim, i_2
    logical, parameter                 :: Debug = .FALSE.
    
    Dim = Size(Psi_1, Dim=2)
    IF (Dim /= Size(Psi_2, Dim=2) .OR. Size(Psi_2, Dim=1) /= Size(Psi_2, Dim=1)) THEN
      WRITE(out_unit,*) "The matrices are expected to have same Dimensions for the scalar product."
      STOP "### The matrices are expected to have same Dimensions for the scalar product."
    END IF

    Sca_pdt = ZERO
    
    DO i_2 = 1, Dim
      Sca_pdt = Sca_pdt + DOT_PRODUCT(Psi_1(:,i_2), Psi_2(:,i_2))
    END DO

    IF (Debug) WRITE(out_unit,*) "Computed scalar product : < Psi_1 |  Psi_2 >  ="//TO_string(Sca_pdt)

  END SUBROUTINE MolecCav_Scalar_product_R2_complex


  SUBROUTINE MolecCav_Scalar_product_R1_real(Sca_pdt, Psi_1, Psi_2)
    USE QDUtil_m
    IMPLICIT NONE
  
    real(kind=Rkind), intent(inout) :: Sca_pdt
    real(kind=Rkind), intent(in)    :: Psi_1(:)                                                                                  ! already allocated
    real(kind=Rkind), intent(in)    :: Psi_2(:)                                                                                  ! already allocated

    integer                         :: Dim
    logical, parameter              :: Debug = .FALSE.
    
    Dim = Size(Psi_1)
    IF (Dim /= Size(Psi_2)) THEN
      WRITE(out_unit,*) "The matrices are expected to have same Dimensions for the scalar product."
      STOP "### The matrices are expected to have same Dimensions for the scalar product."
    END IF

    Sca_pdt = DOT_PRODUCT(Psi_1, Psi_2)

    IF (Debug) WRITE(out_unit,*) "Computed scalar product : < Psi_1 |  Psi_2 >  ="//TO_string(Sca_pdt)

  END SUBROUTINE MolecCav_Scalar_product_R1_real


  SUBROUTINE MolecCav_Scalar_product_R1_complex(Sca_pdt, Psi_1, Psi_2)
    USE QDUtil_m
    IMPLICIT NONE
  
    complex(kind=Rkind), intent(inout) :: Sca_pdt
    complex(kind=Rkind), intent(in)    :: Psi_1(:)                                                                               ! already allocated
    complex(kind=Rkind), intent(in)    :: Psi_2(:)                                                                               ! already allocated

    integer                            :: Dim
    logical, parameter                 :: Debug = .FALSE.
    
    Dim = Size(Psi_1)
    IF (Dim /= Size(Psi_2)) THEN
      WRITE(out_unit,*) "The matrices are expected to have same Dimensions for the scalar product."
      STOP "### The matrices are expected to have same Dimensions for the scalar product."
    END IF

    Sca_pdt = DOT_PRODUCT(Psi_1, Psi_2)

    IF (Debug) WRITE(out_unit,*) "Computed scalar product : < Psi_1 |  Psi_2 >  ="//TO_string(Sca_pdt)

  END SUBROUTINE MolecCav_Scalar_product_R1_complex


END MODULE
  