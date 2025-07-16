!==================================================================================================
!==================================================================================================
!
! This file is part of MolecCav.
!
!==================================================================================================
!
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
!
!==================================================================================================
!
! README :
! Module of utilities, designed to be used only in the tests_*.f90 PROGRAMs. It contains procedures
! to test equality between tensors, and, most importantly, a USE of the QDUtil_Test_m module. This 
! makes it possible to use all the Logical_test-related procedures from the QDUtil library.
! 
! Equal_tensor : Interface for the MolecCav_Equal_*_*_tensor_R* procedures. It takes as arguments a
! logical (Error) that does not need to have been initialized, and two tensors of rank 0, 1, or 2, 
! with real values (Real_1 and Real_2) or complex ones (Cmplx_1 and Cmplx_2). It checks whether th-
! ey can be considered equal by comparing their absolute difference (or the absolute real and imag-
! inary part of it) to a Threshold, and affects Error to .T. if they are different, and .F. otherw-
! ise. 
! 
! MolecCav_Equal_R_R_tensor_R0 : Takes as agruments a logical (Error) that does not need to have b-
! een initialized, and two tensors of rank 0 with real values (Real_1 and Real_2). It checks wheth-
! er they can be considered equal by comparing their absolute difference to a Threshold, and affec-
! ts Error to .T. if they are different, and .F. otherwise. 
! 
! MolecCav_Equal_C_C_tensor_R0 : Same as MolecCav_Equal_R_R_tensor_R0, but for two tensors with co-
! mplex values (Cmplx_1 and Cmplx_2). It checks whether their respective real and imaginary parts 
! can be considered equal by comparing to a Threshold the absolute values of the real and imaginary
! parts of the difference between the two complexes, and affects Error to .T. if any of them are d-
! ifferent, and .F. otherwise. 
!
! MolecCav_Equal_R_R_tensor_R1 : Same as MolecCav_Equal_R_R_tensor_R0, but for two tensors of rank 
! 1. It checks whether all values of same index from the two tensors can be considered equal, and 
! affects Error to .T. if any of them are different, and .F. otherwise. 
! 
! MolecCav_Equal_C_C_tensor_R1 : Same as MolecCav_Equal_R_R_tensor_R0, but for two tensors of rank 
! 1 with complex values (Cmplx_1 and Cmplx_2). It checks whether the real and imaginary parts of a-
! ll values of same index from the two tensors can be considered equal by comparing to a Threshold
! the absolute values of the real and imaginary parts of the difference between the two complexes, 
! and affects Error to .T. if any of them are different, and .F. otherwise.
!
! MolecCav_Equal_R_R_tensor_R2 : Same as MolecCav_Equal_R_R_tensor_R0, but for two tensors of rank 
! 2. It checks whether all values of same index from the two tensors can be considered equal, and 
! affects Error to .T. if any of them are different, and .F. otherwise.
!
! MolecCav_Equal_C_C_tensor_R2 : Same as MolecCav_Equal_R_R_tensor_R0, but for two tensors of rank 
! 2 with complex values (Cmplx_1 and Cmplx_2). It checks whether the real and imaginary parts of a-
! ll values of same index from the two tensors can be considered equal by comparing to a Threshold 
! the absolute values of the real and imaginary parts of the difference between the two complexes, 
! and affects Error to .T. if any of them are different, and .F. otherwise.
!
!==================================================================================================
!==================================================================================================
MODULE Tests_m
  !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
  USE QDUtil_m
  USE QDUtil_Test_m
  IMPLICIT NONE


  !========================================================
  !
  ! In this module everything is left public because we wa-
  ! nt to leave QDUtil_Test_m available juste USEing this 
  ! module, and i don't know how to leave it public otherw-
  ! ise.
  !
  !========================================================
  INTERFACE Equal_tensor
    MODULE PROCEDURE MolecCav_Equal_R_R_tensor_R0, MolecCav_Equal_C_C_tensor_R0, &
                    &MolecCav_Equal_R_R_tensor_R1, MolecCav_Equal_C_C_tensor_R1, &
                    &MolecCav_Equal_R_R_tensor_R2, MolecCav_Equal_C_C_tensor_R2
  END INTERFACE

  
  CONTAINS


  SUBROUTINE MolecCav_Equal_R_R_tensor_R0(Error, Real_1, Real_2, Debug)
    USE QDUtil_m
    IMPLICIT NONE 

    logical,           intent(inout) :: Error                                                           ! does not need to have been initialised
    real(kind=Rkind),  intent(in)    :: Real_1
    real(kind=Rkind),  intent(in)    :: Real_2
    logical, optional, intent(in)    :: Debug
    
    real(kind=Rkind), parameter      :: Threshold   = 1E-10_Rkind
    logical                          :: Debug_local

    !--- Debugging options --------------------------------
    IF (PRESENT(Debug)) THEN; Debug_local = Debug
    ELSE; Debug_local = .FALSE.; END IF

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "o Arguments of MolecCav_Equal_R_R_tensor_R0 :"
      WRITE(out_unit,*) "The <<Error>> argument               :"//TO_string(Error)
      WRITE(out_unit,*) "The <<Real_1>> argument              :"//TO_string(Real_1)
      WRITE(out_unit,*) "The <<Real_2>> argument              :"//TO_string(Real_2)
      WRITE(out_unit,*) "The <<Threshold>> internal parameter :"//TO_string(Threshold)
      FLUSH(out_unit)
    END IF
      
    !--- The test -----------------------------------------
    IF (ABS(Real_1 - Real_2) > Threshold) THEN
      Error = .TRUE.
      IF (Debug_local) WRITE(out_unit,*) "--- The two numbers are not close enough to be considered equal : Re_1 =", Real_1, "Re_2 &
      &=", Real_2, "|Re_1-Re_2| = ", ABS(Real_1 - Real_2)
    ELSE 
      Error = .FALSE.
      IF (Debug_local) WRITE(out_unit,*) "--- The two numbers are close enough to be considered equal : Re_1 =", Real_1, "Re_2 =", &
      &Real_2, "|Re_1-Re_2| = ", ABS(Real_1 - Real_2)
    END IF

    IF (Debug_local) WRITE(out_unit,*) "--- The final value of Error = "//TO_string(Error)

  END SUBROUTINE MolecCav_Equal_R_R_tensor_R0
  

  SUBROUTINE MolecCav_Equal_C_C_tensor_R0(Error, Cmplx_1, Cmplx_2, Debug)
    USE QDUtil_m
    IMPLICIT NONE 
  
    logical,             intent(inout) :: Error
    complex(kind=Rkind), intent(in)    :: Cmplx_1
    complex(kind=Rkind), intent(in)    :: Cmplx_2
    logical, optional,   intent(in)    :: Debug
    
    complex(kind=Rkind)                :: Difference
    real(kind=Rkind), parameter        :: Threshold   = 1E-10_Rkind
    logical                            :: Debug_local

    !--- Debugging options --------------------------------
    IF (PRESENT(Debug)) THEN; Debug_local = Debug
    ELSE; Debug_local = .FALSE.; END IF

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "o Arguments of MolecCav_Equal_C_C_tensor_R0 :"
      WRITE(out_unit,*) "The <<Error>> argument               :"//TO_string(Error)
      WRITE(out_unit,*) "The <<Cmplx_1>> argument             :"//TO_string(Cmplx_1)
      WRITE(out_unit,*) "The <<Cmplx_2>> argument             :"//TO_string(Cmplx_2)
      WRITE(out_unit,*) "The <<Threshold>> internal parameter :"//TO_string(Threshold)
      FLUSH(out_unit)
    END IF
      
    !--- The test -----------------------------------------
    Difference = Cmplx_1 - Cmplx_2

    IF (ABS(Difference%Re) > Threshold) THEN
      Error = .TRUE.
      IF (Debug_local) WRITE(out_unit,*) "--- The real part of the two numbers are not close enough to be considered equal : Re_1&
      & =", Cmplx_1%Re, "Re_2 =", Cmplx_2%Re, "|Re_1-Re_2| = ", ABS(Difference%Re)
    ELSE
      Error = .FALSE.
      IF (Debug_local) WRITE(out_unit,*) "--- The real part of the two numbers are close enough to be considered equal : Re_1 =",&
      & Cmplx_1%Re, "Re_2 =", Cmplx_2%Re, "|Re_1-Re_2| = ", ABS(Difference%Re)
    END IF
  
    IF (ABS(AIMAG(Difference)) > Threshold) THEN
      Error = .TRUE.
      IF (Debug_local) WRITE(out_unit,*) "--- The imaginary part of the two numbers are not close enough to be considered equal :&
      & Im_1 =", AIMAG(Cmplx_1), "Im_2 =", AIMAG(Cmplx_2), "|Im_1-Im_2| = ", ABS(AIMAG(Difference))
    ELSE
      Error = .FALSE.
      IF (Debug_local) WRITE(out_unit,*) "--- The cmplx part of the two numbers are close enough to be considered equal : Im_1 ="&
      &, AIMAG(Cmplx_1), "Im_2 =", AIMAG(Cmplx_2), "|Im_1-Im_2| = ", ABS(AIMAG(Difference))
    END IF

    IF (Debug_local) WRITE(out_unit,*) "--- The final value of Error = "//TO_string(Error)

  END SUBROUTINE MolecCav_Equal_C_C_tensor_R0
    
  
  SUBROUTINE MolecCav_Equal_R_R_tensor_R1(Error, Real_1, Real_2, Debug)
    USE QDUtil_m
    IMPLICIT NONE 

    logical,           intent(inout) :: Error
    real(kind=Rkind),  intent(in)    :: Real_1(:)
    real(kind=Rkind),  intent(in)    :: Real_2(:)
    logical, optional, intent(in)    :: Debug
    
    real(kind=Rkind), parameter      :: Threshold   = 1E-10_Rkind
    integer                          :: Nb_1
    logical                          :: Debug_local

    !--- Debugging options --------------------------------
    IF (PRESENT(Debug)) THEN; Debug_local = Debug
    ELSE; Debug_local = .FALSE.; END IF

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "o Arguments of MolecCav_Equal_R_R_tensor_R1 :"
      WRITE(out_unit,*) "The <<Error>> argument               :"//TO_string(Error)
      WRITE(out_unit,*) "The <<Real_1>> argument              :"
      CALL Write_Vec(Real_1, out_unit, 1, info="Real_1")
      WRITE(out_unit,*) "The <<Real_2>> argument              :"
      CALL Write_Vec(Real_2, out_unit, 1, info="Real_2")
      WRITE(out_unit,*) "The <<Threshold>> internal parameter :"//TO_string(Threshold)
      FLUSH(out_unit)
    END IF
      
    IF (SIZE(Real_1, dim=1) /= SIZE(Real_2, dim=1)) THEN
      WRITE(out_unit,*) "### The two tensors must have same dimensions to compare them. Please, check initialization."
      WRITE(out_unit,*) "    SIZE(Real_1) = "//TO_string(SIZE(Real_1))//"; SIZE(Real_2) = "//TO_string(SIZE(Real_2))
      STOP "### The two tensors must have same dimensions to compare them. Please, check initialization."
    END IF 

    !--- The test -----------------------------------------
    Nb_1 = SIZE(Real_1, dim=1)

    IF (ANY(ABS(Real_1 - Real_2) > Threshold)) THEN
      Error = .TRUE.
      IF (Debug_local) THEN
        WRITE(out_unit,*) "--- The two vectors are not close enough to be considered equal :"
        CALL Write_Vec(Real_1, out_unit, Nb_1, info="Re_1(:)")
        CALL Write_Vec(Real_2, out_unit, Nb_1, info="Re_2(:)")
        CALL Write_Vec(ABS(Real_1 - Real_2), out_unit, Nb_1, info="|Re_1-Re_2| = ")
      END IF 

    ELSE 
      Error = .FALSE.
      IF (Debug_local) THEN
        WRITE(out_unit,*) "The two vectors are close enough to be considered equal :"
        CALL Write_Vec(Real_1, out_unit, Nb_1, info="Re_1(:)")
        CALL Write_Vec(Real_2, out_unit, Nb_1, info="Re_2(:)")
        CALL Write_Vec(ABS(Real_1 - Real_2), out_unit, Nb_1, info="|Re_1-Re_2| = ")
      END IF
    END IF 

    IF (Debug_local) WRITE(out_unit,*) "--- The final value of Error = "//TO_string(Error)

  END SUBROUTINE MolecCav_Equal_R_R_tensor_R1
  

  SUBROUTINE MolecCav_Equal_C_C_tensor_R1(Error, Cmplx_1, Cmplx_2, Debug)
    USE QDUtil_m
    IMPLICIT NONE 
  
    logical,             intent(inout) :: Error
    complex(kind=Rkind), intent(in)    :: Cmplx_1(:)
    complex(kind=Rkind), intent(in)    :: Cmplx_2(:)
    logical, optional,   intent(in)    :: Debug
    
    complex(kind=Rkind), allocatable   :: Difference(:)
    real(kind=Rkind), parameter        :: Threshold   = 1E-10_Rkind
    integer                            :: Nb_1
    logical                            :: Debug_local

    !--- Debugging options --------------------------------
    IF (PRESENT(Debug)) THEN; Debug_local = Debug
    ELSE; Debug_local = .FALSE.; END IF

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "o Arguments of MolecCav_Equal_C_C_tensor_R1 :"
      WRITE(out_unit,*) "The <<Error>> argument                :"//TO_string(Error)
      WRITE(out_unit,*) "The <<Cmplx_1>> argument              :"
      CALL Write_Vec(Cmplx_1, out_unit, 1, info="Cmplx_1")
      WRITE(out_unit,*) "The <<Cmplx_2>> argument              :"
      CALL Write_Vec(Cmplx_2, out_unit, 1, info="Cmplx_2")
      WRITE(out_unit,*) "The <<Threshold>> internal parameter  :"//TO_string(Threshold)
      FLUSH(out_unit)
    END IF
      
    IF (SIZE(Cmplx_1, dim=1) /= SIZE(Cmplx_2, dim=1)) THEN
      WRITE(out_unit,*) "### The two tensors must have same dimensions to compare them. Please, check initialization."
      WRITE(out_unit,*) "    SIZE(Cmplx_1) = "//TO_string(SIZE(Cmplx_1))//"; SIZE(Cmplx_2) = "//TO_string(SIZE(Cmplx_2))
      STOP "### The two tensors must have same dimensions to compare them. Please, check initialization."
    END IF 

    !--- The test -----------------------------------------
    Nb_1 = SIZE(Cmplx_1, dim=1)
    ALLOCATE(Difference(Nb_1))
    Difference(:) = Cmplx_1(:) - Cmplx_2(:)
  
    IF (ANY(ABS(Difference(:)%Re) > Threshold)) THEN
      Error = .TRUE.
      IF (Debug_local) THEN
        WRITE(out_unit,*) "--- The real part of the two vectors are not close enough to be considered equal :"
        CALL Write_Vec(REAL(Cmplx_1, kind=Rkind), out_unit, Nb_1, info="Re(Cmplx_1(:))")
        CALL Write_Vec(REAL(Cmplx_2, kind=Rkind), out_unit, Nb_1, info="Re(Cmplx_2(:))")
        CALL Write_Vec(ABS(REAL(Difference, kind=Rkind)), out_unit, Nb_1, info="|Re(Cmplx_1-Cmplx_2)| = ")
      END IF 

    ELSE
      Error = .FALSE.
      IF (Debug_local) THEN
        WRITE(out_unit,*) "--- The real part of the two vectors are close enough to be considered equal :"
        CALL Write_Vec(REAL(Cmplx_1, kind=Rkind), out_unit, Nb_1, info="Re(Cmplx_1(:))")
        CALL Write_Vec(REAL(Cmplx_2, kind=Rkind), out_unit, Nb_1, info="Re(Cmplx_2(:))")
        CALL Write_Vec(ABS(REAL(Difference, kind=Rkind)), out_unit, Nb_1, info="|Re(Cmplx_1-Cmplx_2)| = ")
      END IF 
    END IF
  
    IF (ANY(ABS(AIMAG(Difference(:))) > Threshold)) THEN
      Error = .TRUE.
      IF (Debug_local) THEN
        WRITE(out_unit,*) "--- The imaginary part of the two vectors are not close enough to be considered equal :"
        CALL Write_Vec(AIMAG(Cmplx_1), out_unit, Nb_1, info="Cmplx_1(:)")
        CALL Write_Vec(AIMAG(Cmplx_2), out_unit, Nb_1, info="Cmplx_2(:)")
        CALL Write_Vec(ABS(AIMAG(Difference)), out_unit, Nb_1, info="|Im(Cmplx_1-Cmplx_2)| = ")
      END IF 

    ELSE
      Error = .FALSE.
      IF (Debug_local) THEN
        WRITE(out_unit,*) "--- The imaginary part of the two vectors are close enough to be considered equal :"
        CALL Write_Vec(AIMAG(Cmplx_1), out_unit, Nb_1, info="Im(Cmplx_1(:))")
        CALL Write_Vec(AIMAG(Cmplx_2), out_unit, Nb_1, info="Im(Cmplx_2(:))")
        CALL Write_Vec(ABS(AIMAG(Difference)), out_unit, Nb_1, info="|Im(Cmplx_1-Cmplx_2)| = ")
      END IF 
    END IF

  END SUBROUTINE MolecCav_Equal_C_C_tensor_R1
    
  
  SUBROUTINE MolecCav_Equal_R_R_tensor_R2(Error, Real_1, Real_2, Debug)
    USE QDUtil_m
    IMPLICIT NONE 

    logical,           intent(inout) :: Error
    real(kind=Rkind),  intent(in)    :: Real_1(:,:)
    real(kind=Rkind),  intent(in)    :: Real_2(:,:)
    logical, optional, intent(in)    :: Debug

    real(kind=Rkind), parameter     :: Threshold   = 1E-10_Rkind
    integer                         :: Nb_1, Nb_2
    logical                         :: Debug_local

    !--- Debugging options --------------------------------
    IF (PRESENT(Debug)) THEN; Debug_local = Debug
    ELSE; Debug_local = .FALSE.; END IF

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "o Arguments of MolecCav_Equal_R_R_tensor_R2 :"
      WRITE(out_unit,*) "The <<Error>> argument               :"//TO_string(Error)
      WRITE(out_unit,*) "The <<Real_1>> argument              :"
      CALL Write_Mat(Real_1, out_unit, SIZE(Real_1, dim=2), info="Real_1")
      WRITE(out_unit,*) "The <<Real_2>> argument              :"
      CALL Write_Mat(Real_2, out_unit, SIZE(Real_2, dim=2), info="Real_2")
      WRITE(out_unit,*) "The <<Threshold>> internal parameter :"//TO_string(Threshold)
      FLUSH(out_unit)
    END IF
      
    IF (SIZE(Real_1, dim=1) /= SIZE(Real_2, dim=1) .OR. SIZE(Real_1, dim=2) /= SIZE(Real_2, dim=2)) THEN
      WRITE(out_unit,*) "### The two tensors must have same sizes to compare them. Please, check initialization."
      WRITE(out_unit,*) "    SIZE(Real_1, dim=1) = "//TO_string(SIZE(Real_1, dim=1))//"; SIZE(Real_1, dim=2) = "//TO_string(SIZE(&
      &Real_1, dim=2))
      WRITE(out_unit,*) "    SIZE(Real_2, dim=1) = "//TO_string(SIZE(Real_2, dim=1))//"; SIZE(Real_2, dim=2) = "//TO_string(SIZE(&
      &Real_2, dim=2))
      STOP "### The two tensors must have same dimensions to compare them. Please, check initialization."
    END IF 

    !--- The test -----------------------------------------
    Nb_1 = Size(Real_1, dim=1)
    Nb_2 = Size(Real_2, dim=2)

    IF (ANY(ABS(Real_1 - Real_2) > Threshold)) THEN
      Error = .TRUE.
      IF (Debug_local) THEN
        WRITE(out_unit,*) "--- The two matrices are not close enough to be considered equal :"
        CALL Write_Mat(Real_1, out_unit, Nb_2, info="Re_1(:,:)")
        CALL Write_Mat(Real_2, out_unit, Nb_2, info="Re_2(:,:)")
        CALL Write_Mat(ABS(Real_1 - Real_2), out_unit, Nb_2, info="|Re_1-Re_2| = ")
      END IF 

    ELSE 
      Error = .FALSE.
      IF (Debug_local) THEN
        WRITE(out_unit,*) "--- The two matrices are close enough to be considered equal :"
        CALL Write_Mat(Real_1, out_unit, Nb_2, info="Re_1(:,:)")
        CALL Write_Mat(Real_2, out_unit, Nb_2, info="Re_2(:,:)")
        CALL Write_Mat(ABS(Real_1 - Real_2), out_unit, Nb_2, info="|Re_1-Re_2| = ")
      END IF
    END IF 

    IF (Debug_local) WRITE(out_unit,*) "--- The final value of Error = "//TO_string(Error)

  END SUBROUTINE MolecCav_Equal_R_R_tensor_R2
  

  SUBROUTINE MolecCav_Equal_C_C_tensor_R2(Error, Cmplx_1, Cmplx_2, Debug)
    USE QDUtil_m
    IMPLICIT NONE 
  
    logical,             intent(inout) :: Error
    complex(kind=Rkind), intent(in)    :: Cmplx_1(:,:)
    complex(kind=Rkind), intent(in)    :: Cmplx_2(:,:)
    logical, optional,   intent(in)    :: Debug
    
    complex(kind=Rkind), allocatable   :: Difference(:,:)
    real(kind=Rkind), parameter        :: Threshold   = 1E-10_Rkind
    integer                            :: Nb_1, Nb_2
    logical                            :: Debug_local

    !--- Debugging options --------------------------------
    IF (PRESENT(Debug)) THEN; Debug_local = Debug
    ELSE; Debug_local = .FALSE.; END IF

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "o Arguments of MolecCav_Equal_C_C_tensor_R2 :"
      WRITE(out_unit,*) "The <<Error>> argument               :"//TO_string(Error)
      WRITE(out_unit,*) "The <<Cmplx_1>> argument             :"
      CALL Write_Mat(Cmplx_1, out_unit, SIZE(Cmplx_1, dim=2), info="Cmplx_1")
      WRITE(out_unit,*) "The <<Cmplx_2>> argument             :"
      CALL Write_Mat(Cmplx_2, out_unit, SIZE(Cmplx_2, dim=2), info="Cmplx_2")
      WRITE(out_unit,*) "The <<Threshold>> internal parameter :"//TO_string(Threshold)
      FLUSH(out_unit)
    END IF
      
    IF (SIZE(Cmplx_1, dim=1) /= SIZE(Cmplx_2, dim=1) .OR. SIZE(Cmplx_1, dim=2) /= SIZE(Cmplx_2, dim=2)) THEN
      WRITE(out_unit,*) "### The two tensors must have same sizes to compare them. Please, check initialization."
      WRITE(out_unit,*) "    SIZE(Cmplx_1, dim=1) = "//TO_string(SIZE(Cmplx_1, dim=1))//"; SIZE(Cmplx_1, dim=2) = "//TO_string(SI&
      &ZE(Cmplx_1, dim=2))
      WRITE(out_unit,*) "    SIZE(Cmplx_2, dim=1) = "//TO_string(SIZE(Cmplx_2, dim=1))//"; SIZE(Cmplx_2, dim=2) = "//TO_string(SI&
      &ZE(Cmplx_2, dim=2))
      STOP "### The two tensors must have same dimensions to compare them. Please, check initialization."
    END IF 

    !--- The test -----------------------------------------
    Nb_1 = Size(Cmplx_1, dim=1)
    Nb_2 = Size(Cmplx_2, dim=2)
    ALLOCATE(Difference(Nb_1, Nb_2))
    Difference(:,:) = Cmplx_1(:,:) - Cmplx_2(:,:)
  
    IF (ANY(ABS(Difference(:,:)%Re) > Threshold)) THEN
      Error = .TRUE.
      IF (Debug_local) THEN
        WRITE(out_unit,*) "--- The real part of the two matrices are not close enough to be considered equal :"
        CALL Write_Mat(REAL(Cmplx_1, kind=Rkind), out_unit, Nb_2, info="Re(Cmplx_1(:,:))")
        CALL Write_Mat(REAL(Cmplx_2, kind=Rkind), out_unit, Nb_2, info="Re(Cmplx_2(:,:))")
        CALL Write_Mat(ABS(REAL(Difference, kind=Rkind)), out_unit, Nb_2, info="|Re(Cmplx_1-Cmplx_2)| = ")
      END IF 

    ELSE
      Error = .FALSE.
      IF (Debug_local) THEN
        WRITE(out_unit,*) "--- The real part of the two matrices are close enough to be considered equal :"
        CALL Write_Mat(REAL(Cmplx_1, kind=Rkind), out_unit, Nb_2, info="Re(Cmplx_1(:,:))")
        CALL Write_Mat(REAL(Cmplx_2, kind=Rkind), out_unit, Nb_2, info="Re(Cmplx_2(:,:))")
        CALL Write_Mat(ABS(REAL(Difference, kind=Rkind)), out_unit, Nb_2, info="|Re(Cmplx_1-Cmplx_2)| = ")
      END IF 
    END IF
  
    IF (ANY(ABS(AIMAG(Difference(:,:))) > Threshold)) THEN
      Error = .TRUE.
      IF (Debug_local) THEN
        WRITE(out_unit,*) "--- The imaginary part of the two matrices are not close enough to be considered equal :"
        CALL Write_Mat(AIMAG(Cmplx_1), out_unit, Nb_2, info="Cmplx_1(:,:)")
        CALL Write_Mat(AIMAG(Cmplx_2), out_unit, Nb_2, info="Cmplx_2(:,:)")
        CALL Write_Mat(ABS(AIMAG(Difference)), out_unit, Nb_2, info="|Im(Cmplx_1-Cmplx_2)| = ")
      END IF 

    ELSE
      Error = .FALSE.
      IF (Debug_local) THEN
        WRITE(out_unit,*) "--- The imaginary part of the two matrices are close enough to be considered equal :"
        CALL Write_Mat(AIMAG(Cmplx_1), out_unit, Nb_2, info="Im(Cmplx_1(:,:))")
        CALL Write_Mat(AIMAG(Cmplx_2), out_unit, Nb_2, info="Im(Cmplx_2(:,:))")
        CALL Write_Mat(ABS(AIMAG(Difference)), out_unit, Nb_2, info="|Im(Cmplx_1-Cmplx_2)| = ")
      END IF 
    END IF

    IF (Debug_local) WRITE(out_unit,*) "--- The final value of Error = "//TO_string(Error)

  END SUBROUTINE MolecCav_Equal_C_C_tensor_R2
    
  
END MODULE