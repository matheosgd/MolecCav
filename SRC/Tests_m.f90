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
MODULE Tests_m
  !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
  USE QDUtil_m
  USE QDUtil_Test_m
  IMPLICIT NONE

  
!  PRIVATE

!  PUBLIC Equal_tensor ! because we want to leave QDUtil_Test_m available juste USEing this module, and idk how to leave public otherwise

  INTERFACE Equal_tensor
    MODULE PROCEDURE MolecCav_Equal_R_R_tensor_R0, MolecCav_Equal_C_C_tensor_R0, &
                    &MolecCav_Equal_R_R_tensor_R1, MolecCav_Equal_C_C_tensor_R1, &
                    &MolecCav_Equal_R_R_tensor_R2, MolecCav_Equal_C_C_tensor_R2
  END INTERFACE

  
  CONTAINS


  SUBROUTINE MolecCav_Equal_R_R_tensor_R0(error, Rl_1, Rl_2)
    USE QDUtil_m
    IMPLICIT NONE 

    logical,          intent(inout) :: error
    real(kind=Rkind), intent(in)    :: Rl_1
    real(kind=Rkind), intent(in)    :: Rl_2
    
    real(kind=Rkind), parameter     :: Threshold   = 1E-10_Rkind
    logical, parameter              :: Debug_local = .FALSE.

    IF (ABS(Rl_1 - Rl_2) > Threshold) THEN
      error = .TRUE.
      IF (Debug_local) WRITE(out_unit,*) "The two numbers are not close enough to be considered equ&
                                         &al : R_1 =", Rl_1, "R_2 =", Rl_2, "|R_1-R_2| = ", ABS(Rl_1 - Rl_2)
    ELSE 
      error = .FALSE.
      IF (Debug_local) WRITE(out_unit,*) "The two numbers are close enough to be considered equal :&
                                         & R_1 =", Rl_1, "R_2 =", Rl_2, "|R_1-R_2| = ", ABS(Rl_1 - Rl_2)
    END IF

  END SUBROUTINE MolecCav_Equal_R_R_tensor_R0
  

  SUBROUTINE MolecCav_Equal_C_C_tensor_R0(error, Cmplx_1, Cmplx_2)
    USE QDUtil_m
    IMPLICIT NONE 
  
    logical,             intent(inout) :: error
    complex(kind=Rkind), intent(in)    :: Cmplx_1
    complex(kind=Rkind), intent(in)    :: Cmplx_2
    
    complex(kind=Rkind)                :: Difference
    real(kind=Rkind), parameter        :: Threshold   = 1E-10_Rkind
    logical,          parameter        :: Debug_local = .FALSE.
  
    Difference = Cmplx_1 - Cmplx_2
  
    IF (ABS(Difference%Re) > Threshold) THEN
      error = .TRUE.
      IF (Debug_local) WRITE(out_unit,*) "The real part of the two numbers are not close enough to &
                                         &be considered equal : Re_1 =", Cmplx_1%Re, "Re_2 =", Cmpl&
                                         &x_2%Re, "|Re_1-Re_2| = ", ABS(Difference%Re)
    ELSE
      error = .FALSE.
      IF (Debug_local) WRITE(out_unit,*) "The real part of the two numbers are close enough to be c&
                                         &onsidered equal : Re_1 =", Cmplx_1%Re, "Re_2 =", Cmplx_2%&
                                         &Re, "|Re_1-Re_2| = ", ABS(Difference%Re)
    END IF
  
    IF (ABS(AIMAG(Difference)) > Threshold) THEN
      error = .TRUE.
      IF (Debug_local) WRITE(out_unit,*) "The cmplx part of the two numbers are not close enough to&
                                         & be considered equal : Im_1 =", AIMAG(Cmplx_1), "Im_2 =",&
                                         & AIMAG(Cmplx_2), "|Im_1-Re_2| = ", ABS(AIMAG(Difference))
    ELSE
      error = .FALSE.
      IF (Debug_local) WRITE(out_unit,*) "The cmplx part of the two numbers are close enough to be &
                                         &considered equal : Im_1 =", AIMAG(Cmplx_1), "Im_2 =", AIM&
                                         &AG(Cmplx_2), "|Im_1-Im_2| = ", ABS(AIMAG(Difference))
    END IF

  END SUBROUTINE MolecCav_Equal_C_C_tensor_R0
    
  
  SUBROUTINE MolecCav_Equal_R_R_tensor_R1(error, Rl_1, Rl_2)
    USE QDUtil_m
    IMPLICIT NONE 

    logical,          intent(inout) :: error
    real(kind=Rkind), intent(in)    :: Rl_1(:)
    real(kind=Rkind), intent(in)    :: Rl_2(:)
    
    real(kind=Rkind), parameter     :: Threshold   = 1E-10_Rkind
    logical, parameter              :: Debug_local = .FALSE.
    integer                         :: Nb_1

    Nb_1 = Size(Rl_1)
    IF (Nb_1 /= Size(Rl_2)) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "The two vectors must have same dimensions to compare them. Please, check initialization."
      WRITE(out_unit,*) "Size(Rl_1) = "//TO_string(Size(Rl_1))//"; Size(Rl_2) = "//TO_string(Size(Rl_2))
      STOP "### The two vectors must have same dimensions to compare them. Please, check initialization."
    END IF 

    IF (ANY(ABS(Rl_1 - Rl_2) > Threshold)) THEN
      error = .TRUE.
      IF (Debug_local) THEN
        WRITE(out_unit,*) "The two vectors are not close enough to be considered equal :"
        CALL Write_Vec(Rl_1, out_unit, Nb_1, info="R_1(:)")
        CALL Write_Vec(Rl_2, out_unit, Nb_1, info="R_2(:)")
        CALL Write_Vec(ABS(Rl_1 - Rl_2), out_unit, Nb_1, info="|R_1-R_2| = ")
      END IF 

    ELSE 
      error = .FALSE.
      IF (Debug_local) THEN
        WRITE(out_unit,*) "The two vectors are close enough to be considered equal :"
        CALL Write_Vec(Rl_1, out_unit, Nb_1, info="R_1(:)")
        CALL Write_Vec(Rl_2, out_unit, Nb_1, info="R_2(:)")
        CALL Write_Vec(ABS(Rl_1 - Rl_2), out_unit, Nb_1, info="|R_1-R_2| = ")
      END IF
    END IF 

  END SUBROUTINE MolecCav_Equal_R_R_tensor_R1
  

  SUBROUTINE MolecCav_Equal_C_C_tensor_R1(error, Cmplx_1, Cmplx_2)
    USE QDUtil_m
    IMPLICIT NONE 
  
    logical,             intent(inout) :: error
    complex(kind=Rkind), intent(in)    :: Cmplx_1(:)
    complex(kind=Rkind), intent(in)    :: Cmplx_2(:)
    
    complex(kind=Rkind), allocatable   :: Difference(:)
    real(kind=Rkind), parameter        :: Threshold   = 1E-10_Rkind
    logical,          parameter        :: Debug_local = .FALSE.
    integer                            :: Nb_1

    Nb_1 = Size(Cmplx_1)
    IF (Nb_1 /= Size(Cmplx_2)) THEN
      WRITE(out_unit,*) "The two vectors must have same dimensions to compare them. Please, check initialization."
      STOP "### The two vectors must have same dimensions to compare them. Please, check initialization."
    END IF 
  
    ALLOCATE(Difference(Nb_1))
    Difference(:) = Cmplx_1(:) - Cmplx_2(:)
  
    IF (ANY(ABS(Difference(:)%Re) > Threshold)) THEN
      error = .TRUE.
      IF (Debug_local) THEN
        WRITE(out_unit,*) "The real part of the two vectors are not close enough to be considered equal :"
        CALL Write_Vec(REAL(Cmplx_1, kind=Rkind), out_unit, Nb_1, info="Re(Cmplx_1(:))")
        CALL Write_Vec(REAL(Cmplx_2, kind=Rkind), out_unit, Nb_1, info="Re(Cmplx_2(:))")
        CALL Write_Vec(ABS(REAL(Difference, kind=Rkind)), out_unit, Nb_1, info="|Re(Cmplx_1-Cmplx_2)| = ")
      END IF 

    ELSE
      error = .FALSE.
      IF (Debug_local) THEN
        WRITE(out_unit,*) "The real part of the two vectors are close enough to be considered equal :"
        CALL Write_Vec(REAL(Cmplx_1, kind=Rkind), out_unit, Nb_1, info="Re(Cmplx_1(:))")
        CALL Write_Vec(REAL(Cmplx_2, kind=Rkind), out_unit, Nb_1, info="Re(Cmplx_2(:))")
        CALL Write_Vec(ABS(REAL(Difference, kind=Rkind)), out_unit, Nb_1, info="|Re(Cmplx_1-Cmplx_2)| = ")
      END IF 
    END IF
  
    IF (ANY(ABS(AIMAG(Difference(:))) > Threshold)) THEN
      error = .TRUE.
      IF (Debug_local) THEN
        WRITE(out_unit,*) "The imaginary part of the two vectors are not close enough to be considered equal :"
        CALL Write_Vec(AIMAG(Cmplx_1), out_unit, Nb_1, info="Cmplx_1(:)")
        CALL Write_Vec(AIMAG(Cmplx_2), out_unit, Nb_1, info="Cmplx_2(:)")
        CALL Write_Vec(ABS(AIMAG(Difference)), out_unit, Nb_1, info="|Im(Cmplx_1-Cmplx_2)| = ")
      END IF 

    ELSE
      error = .FALSE.
      IF (Debug_local) THEN
        WRITE(out_unit,*) "The imaginary part of the two vectors are close enough to be considered equal :"
        CALL Write_Vec(AIMAG(Cmplx_1), out_unit, Nb_1, info="Im(Cmplx_1(:))")
        CALL Write_Vec(AIMAG(Cmplx_2), out_unit, Nb_1, info="Im(Cmplx_2(:))")
        CALL Write_Vec(ABS(AIMAG(Difference)), out_unit, Nb_1, info="|Im(Cmplx_1-Cmplx_2)| = ")
      END IF 
    END IF

  END SUBROUTINE MolecCav_Equal_C_C_tensor_R1
    
  
  SUBROUTINE MolecCav_Equal_R_R_tensor_R2(error, Rl_1, Rl_2)
    USE QDUtil_m
    IMPLICIT NONE 

    logical,          intent(inout) :: error
    real(kind=Rkind), intent(in)    :: Rl_1(:,:)
    real(kind=Rkind), intent(in)    :: Rl_2(:,:)
    
    real(kind=Rkind), parameter     :: Threshold   = 1E-10_Rkind
    logical, parameter              :: Debug_local = .FALSE.
    integer                         :: Nb_1, Nb_2

    Nb_1 = Size(Rl_1, dim=1)
    Nb_2 = Size(Rl_2, dim=2)
    IF (Nb_1 /= Size(Rl_2, dim=1) .OR. Nb_2 /= Size(Rl_2, dim=2)) THEN
      WRITE(out_unit,*) "The two matrices must have same dimensions to compare them. Please, check initialization."
      STOP "### The two matrices must have same dimensions to compare them. Please, check initialization."
    END IF 

    IF (ANY(ABS(Rl_1 - Rl_2) > Threshold)) THEN
      error = .TRUE.
      IF (Debug_local) THEN
        WRITE(out_unit,*) "The two matrices are not close enough to be considered equal :"
        CALL Write_Mat(Rl_1, out_unit, Nb_2, info="R_1(:,:)")
        CALL Write_Mat(Rl_2, out_unit, Nb_2, info="R_2(:,:)")
        CALL Write_Mat(ABS(Rl_1 - Rl_2), out_unit, Nb_2, info="|R_1-R_2| = ")
      END IF 

    ELSE 
      error = .FALSE.
      IF (Debug_local) THEN
        WRITE(out_unit,*) "The two matrices are close enough to be considered equal :"
        CALL Write_Mat(Rl_1, out_unit, Nb_2, info="R_1(:,:)")
        CALL Write_Mat(Rl_2, out_unit, Nb_2, info="R_2(:,:)")
        CALL Write_Mat(ABS(Rl_1 - Rl_2), out_unit, Nb_2, info="|R_1-R_2| = ")
      END IF
    END IF 

  END SUBROUTINE MolecCav_Equal_R_R_tensor_R2
  

  SUBROUTINE MolecCav_Equal_C_C_tensor_R2(error, Cmplx_1, Cmplx_2)
    USE QDUtil_m
    IMPLICIT NONE 
  
    logical,             intent(inout) :: error
    complex(kind=Rkind), intent(in)    :: Cmplx_1(:,:)
    complex(kind=Rkind), intent(in)    :: Cmplx_2(:,:)
    
    complex(kind=Rkind), allocatable   :: Difference(:,:)
    real(kind=Rkind), parameter        :: Threshold   = 1E-10_Rkind
    logical,          parameter        :: Debug_local = .FALSE.
    integer                            :: Nb_1, Nb_2

    Nb_1 = Size(Cmplx_1, dim=1)
    Nb_2 = Size(Cmplx_2, dim=2)
    IF (Nb_1 /= Size(Cmplx_2, dim=1) .OR. Nb_2 /= Size(Cmplx_2, dim=2)) THEN
      WRITE(out_unit,*) "The two matrices must have same dimensions to compare them. Please, check initialization."
      STOP "### The two matrices must have same dimensions to compare them. Please, check initialization."
    END IF 
  
    ALLOCATE(Difference(Nb_1, Nb_2))
    Difference(:,:) = Cmplx_1(:,:) - Cmplx_2(:,:)
  
    IF (ANY(ABS(Difference(:,:)%Re) > Threshold)) THEN
      error = .TRUE.
      IF (Debug_local) THEN
        WRITE(out_unit,*) "The real part of the two matrices are not close enough to be considered equal :"
        CALL Write_Mat(REAL(Cmplx_1, kind=Rkind), out_unit, Nb_2, info="Re(Cmplx_1(:,:))")
        CALL Write_Mat(REAL(Cmplx_2, kind=Rkind), out_unit, Nb_2, info="Re(Cmplx_2(:,:))")
        CALL Write_Mat(ABS(REAL(Difference, kind=Rkind)), out_unit, Nb_2, info="|Re(Cmplx_1-Cmplx_2)| = ")
      END IF 

    ELSE
      error = .FALSE.
      IF (Debug_local) THEN
        WRITE(out_unit,*) "The real part of the two matrices are close enough to be considered equal :"
        CALL Write_Mat(REAL(Cmplx_1, kind=Rkind), out_unit, Nb_2, info="Re(Cmplx_1(:,:))")
        CALL Write_Mat(REAL(Cmplx_2, kind=Rkind), out_unit, Nb_2, info="Re(Cmplx_2(:,:))")
        CALL Write_Mat(ABS(REAL(Difference, kind=Rkind)), out_unit, Nb_2, info="|Re(Cmplx_1-Cmplx_2)| = ")
      END IF 
    END IF
  
    IF (ANY(ABS(AIMAG(Difference(:,:))) > Threshold)) THEN
      error = .TRUE.
      IF (Debug_local) THEN
        WRITE(out_unit,*) "The imaginary part of the two matrices are not close enough to be considered equal :"
        CALL Write_Mat(AIMAG(Cmplx_1), out_unit, Nb_2, info="Cmplx_1(:,:)")
        CALL Write_Mat(AIMAG(Cmplx_2), out_unit, Nb_2, info="Cmplx_2(:,:)")
        CALL Write_Mat(ABS(AIMAG(Difference)), out_unit, Nb_2, info="|Im(Cmplx_1-Cmplx_2)| = ")
      END IF 

    ELSE
      error = .FALSE.
      IF (Debug_local) THEN
        WRITE(out_unit,*) "The imaginary part of the two matrices are close enough to be considered equal :"
        CALL Write_Mat(AIMAG(Cmplx_1), out_unit, Nb_2, info="Im(Cmplx_1(:,:))")
        CALL Write_Mat(AIMAG(Cmplx_2), out_unit, Nb_2, info="Im(Cmplx_2(:,:))")
        CALL Write_Mat(ABS(AIMAG(Difference)), out_unit, Nb_2, info="|Im(Cmplx_1-Cmplx_2)| = ")
      END IF 
    END IF

  END SUBROUTINE MolecCav_Equal_C_C_tensor_R2
    
  
END MODULE