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
! Module of mathematics utilities. It contains some procedures for algebraic tools, mostly about m-
! atrix management and a Gram-Schmidt procedure. 
! This module is in the "-6" level.
! 
! Gram_schmidt : Interface for the MolecCav_Gram_schmidt procedure. cf. MolecCav_Gram_schmidt for 
! details. HAVE NOT BEEN TESTED IN THIS PROGRAM.
!
! Normalize : Interface for the MolecCav_Normalize_R*_* procedures. It takes as arguments a tensor 
! of rank 1 or 2 with real or complex values (Psi) and normalises it by redirecting to one of the 
! associated procedure according to the arguments. It computes its norm calling the Norm_of subrou-
! tine, and divide Psi by it.
!
! Norm_of : Interface for the MolecCav_Norm_R*_* procedures. Takes as arguments a real (Norm) that 
! does not need to have been initialized, and a tensor of rank 1 or 2 with real or complex values 
! (Psi). It computes the norm of Psi by redirecting to one of the associated procedure according t-
! o the arguments. It computes the scalar product by calling the Scalar_product subroutine and aff-
! ects sqrt of the result to Norm.
!
! Scalar_product : Interface for the MolecCav_Scalar_product_R*_* procedures. Takes as arguments a
! real or complex (ScaP) that does not need to have been initialized, and two tensors of rank 1 or 
! 2 with real or complex values (Psi_1 and Psi_2). It computes the scalar product between Psi_1 an-
! d Psi_2 by redirecting to one of the associated procedure according to the arguments. It compute-
! s the scalar product by summing the results from the intrinsic DOT_PRODUCT for each column of sa-
! me index of the tensors, and affects the result to ScaP.
!
! MolecCav_Gram_schmidt : Takes as arguments a rank 2 tensor with real values that does not need to
! have been initialized (OrthoBasis) and an other one that need to have actual values in (NonOrtho-
! Basis). It applies the Gram-Schmidt procedure to NonOrthoBasis to orthonormalise it, going throu-
! gh an intermediary rank 1 tensor to manage the column vectors individually, and affect the resul-
! t to OrthoBasis. HAVE NOT BEEN TESTED IN THIS PROGRAM.
!
! MolecCav_Normalize_R2_real : Takes as arguments a tensor of rank 2 with real values (Psi) and no-
! rmalises it. It computes its norm calling the Norm_of subroutine, and divide Psi by it. 
!
! MolecCav_Normalize_R2_complex : Same as MolecCav_Normalize_R2_real, but for a tensor with complex
! values.
!
! MolecCav_Normalize_R1_real : Same as MolecCav_Normalize_R2_real, but for a tensor of rank 1.
!
! MolecCav_Normalize_R1_complex : Same as MolecCav_Normalize_R2_real but for a tensor of rank 1 wi-
! th complex values.
!
! MolecCav_Norm_R2_real : Takes as arguments a real (Norm) that does not need to have been initial-
! ized, and a tensor of rank 2 with real values (Psi). It computes the norm of Psi by calling Scal-
! ar_product and affects the result to Norm.
!
! MolecCav_Norm_R2_complex : Same as MolecCav_Norm_R2_real, but for a tensor with complex values.
!
! MolecCav_Norm_R1_real : Same as MolecCav_Norm_R2_real, but for a tensor of rank 1.
!
! MolecCav_Norm_R1_complex : Same as MolecCav_Norm_R2_real but for a tensor of rank 1 with complex 
! values.
!
! MolecCav_Scalar_product_R2_real : Takes as arguments a real (ScaP) that does not need to have be-
! en initialized, and two tensors of rank 2 with real values (Psi_1 and Psi_2). It computes the sc-
! alar product between Psi_1 and Psi_2 by summing the results from the intrinsic DOT_PRODUCT for e-
! ach column of same index of the tensors, and affects the result to ScaP.
!
! MolecCav_Scalar_product_R2_complex : Same as MolecCav_Scalar_product_R2_real, but for three argu-
! ments with complex values.
!
! MolecCav_Scalar_product_R1_real : Same as MolecCav_Scalar_product_R2_real, but for two tensors of
! rank 1.
!
! MolecCav_Scalar_product_R1_complex : Same as MolecCav_Scalar_product_R2_real but for two tensors 
! of rank 1, and three arguments with complex values.
!
!==================================================================================================
!==================================================================================================
MODULE Algebra_m
  !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
  USE QDUtil_m
  IMPLICIT NONE

  
  PRIVATE

  PUBLIC Gram_schmidt, Normalize, Norm_of, Scalar_product

  INTERFACE Gram_schmidt
    MODULE PROCEDURE MolecCav_Gram_schmidt
  END INTERFACE
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


  SUBROUTINE MolecCav_Gram_schmidt(OrthoBasis, NonOrthoBasis, Debug)
    USE QDUtil_m
    IMPLICIT NONE
    
    real(kind=Rkind),  intent(inout) :: OrthoBasis(:,:)                                                 ! has to have been allocated BEFORE calling MolecCav_Gram_schmidt
    real(kind=Rkind),  intent(in)    :: NonOrthoBasis(:,:)                                              ! has to have been allocated BEFORE calling MolecCav_Gram_schmidt
    logical, optional, intent(in)    :: Debug

    real(kind=Rkind), allocatable    :: Intermediary(:)
    integer                          :: Nb_1, Nb_2, i, j
    logical                          :: Debug_local

    !--- Debugging options --------------------------------
    IF (PRESENT(Debug)) THEN; Debug_local = Debug
    ELSE; Debug_local = .FALSE.; END IF

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "o Arguments of MolecCav_Gram_schmidt :"
      WRITE(out_unit,*) "The <<OrthoBasis>> argument    :" 
      CALL Write_Mat(OrthoBasis, out_unit, SIZE(OrthoBasis, dim=2), info="OrthoBasis")
      WRITE(out_unit,*) "The <<NonOrthoBasis>> argument :" 
      CALL Write_Mat(NonOrthoBasis, out_unit, SIZE(NonOrthoBasis, dim=2), info="NonOrthoBasis")
      FLUSH(out_unit)
    END IF

    !--- Initialization -----------------------------------
    Nb_1 = SIZE(NonOrthoBasis, dim=1)
    Nb_2 = SIZE(NonOrthoBasis, dim=2)    
    ALLOCATE(Intermediary(Nb_1))
    OrthoBasis(:,:) = ZERO
    Intermediary(:) = ZERO
  
    !--- Orthonormalisation -------------------------------
    OrthoBasis(:,1) = NonOrthoBasis(:,1) / SQRT(DOT_PRODUCT(NonOrthoBasis(:,1),NonOrthoBasis(:,1)))
    
    DO i = 2, Nb_2
      Intermediary(:) = NonOrthoBasis(:, i)
      DO j = 1, i-1
        Intermediary(:) = Intermediary(:) - DOT_PRODUCT(NonOrthoBasis(:,i), OrthoBasis(:,j)) * OrthoBasis(:,j) ! removes projection of the new vector upon every vectors of the ortho basis
      END DO
      OrthoBasis(:, i) = Intermediary(:) / SQRT(DOT_PRODUCT(Intermediary(:),Intermediary(:)))
      IF (Debug_local) CALL Write_Mat(OrthoBasis, out_unit, SIZE(OrthoBasis, dim=1), info="Orthobasis(:,"//TO_string(i)//") :")
    END DO
    
    !--- Finishing initialisation -------------------------
    IF (Debug_local) THEN
      WRITE(out_unit,*) "--- The orthonormalised basis :"
      CALL Write_Mat(OrthoBasis, out_unit, SIZE(OrthoBasis, dim=2), info="OrthoBasis")
    END IF 


  END SUBROUTINE MolecCav_Gram_schmidt
  
  
  SUBROUTINE MolecCav_Normalize_R2_real(Psi, Debug)
    USE QDUtil_m
    IMPLICIT NONE

    real(kind=Rkind),  intent(inout) :: Psi(:,:)                                                        ! already allocated
    logical, optional, intent(in)    :: Debug

    real(kind=Rkind)                 :: Norm
    real(kind=Rkind), parameter      :: Threshold = 1E-10_Rkind
    logical                          :: Debug_local

    !--- Debugging options --------------------------------
    IF (PRESENT(Debug)) THEN; Debug_local = Debug
    ELSE; Debug_local = .FALSE.; END IF

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "o Arguments of MolecCav_Normalize_R2_real :"
      WRITE(out_unit,*) "The <<Psi>> argument                 :"
      CALL Write_Mat(Psi, out_unit, SIZE(Psi, dim=2), info="Psi")
      WRITE(out_unit,*) "The <<Threshold>> internal parameter :"//TO_string(Threshold)
      FLUSH(out_unit)
    END IF

    IF (Debug_local) WRITE(out_unit,*) "--- Computing norm of Psi..."
    CALL Norm_of(Norm, Psi, Debug_local)                                                                ! the printing of the norm according to Debug is within the call of Norm_of
    IF (Debug_local) WRITE(out_unit,*) "    ...back to MolecCav_Normalize_R2_real"

    IF (Norm < Threshold) THEN
      WRITE(out_unit,*) "### Attempt to normalize matrix of norm ZERO at MolecCav_Normalize_R2_real."
      STOP "### Attempt to normalize matrix of norm ZERO at MolecCav_Normalize_R2_real."
      
    !--- Normalisation ------------------------------------
    ELSE 
      Psi(:,:) = Psi(:,:) / Norm
    END IF 

    !--- Conclusion ---------------------------------------
    IF (Debug_local) THEN
      WRITE(out_unit,*) "--- The normalized Psi :"
      CALL Write_Mat(Psi, out_unit, SIZE(Psi, dim=2), info="Normalised Psi")
      WRITE(out_unit,*) "--- Computing the new norm of Psi..."
      CALL Norm_of(Norm, Psi, Debug_local)                                                              ! the printing of the norm according to Debug is within the call of Norm_of
      WRITE(out_unit,*) "    ...back to MolecCav_Normalize_R2_real"
    END IF 

  END SUBROUTINE MolecCav_Normalize_R2_real

  
  SUBROUTINE MolecCav_Normalize_R2_complex(Psi, Debug)
    USE QDUtil_m
    IMPLICIT NONE

    complex(kind=Rkind), intent(inout) :: Psi(:,:)                                                      ! already allocated
    logical, optional,   intent(in)    :: Debug

    real(kind=Rkind)                   :: Norm
    real(kind=Rkind), parameter        :: Threshold = 1E-10_Rkind
    logical                            :: Debug_local

    !--- Debugging options --------------------------------
    IF (PRESENT(Debug)) THEN; Debug_local = Debug
    ELSE; Debug_local = .FALSE.; END IF

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "o Arguments of MolecCav_Normalize_R2_complex :"
      WRITE(out_unit,*) "The <<Psi>> argument                 :"
      CALL Write_Mat(Psi, out_unit, SIZE(Psi, dim=2), info="Psi")
      WRITE(out_unit,*) "The <<Threshold>> internal parameter :"//TO_string(Threshold)
      FLUSH(out_unit)
    END IF

    IF (Debug_local) WRITE(out_unit,*) "--- Computing norm of Psi..."
    CALL Norm_of(Norm, Psi, Debug_local)                                                                ! the printing of the norm according to Debug is within the call of Norm_of
    IF (Debug_local) WRITE(out_unit,*) "    ...back to MolecCav_Normalize_R2_complex"

    IF (Norm < Threshold) THEN
      WRITE(out_unit,*) "### Attempt to normalize matrix of norm ZERO at MolecCav_Normalize_R2_complex."
      STOP "### Attempt to normalize matrix of norm ZERO at MolecCav_Normalize_R2_complex."
      
    !--- Normalisation ------------------------------------
    ELSE 
      Psi(:,:) = Psi(:,:) / Norm
    END IF 

    !--- Conclusion ---------------------------------------
    IF (Debug_local) THEN
      WRITE(out_unit,*) "--- The normalized Psi :"
      CALL Write_Mat(Psi, out_unit, SIZE(Psi, dim=2), info="Normalised Psi")
      WRITE(out_unit,*) "--- Computing the new norm of Psi..."
      CALL Norm_of(Norm, Psi, Debug_local)                                                              ! the printing of the norm according to Debug is within the call of Norm_of
      WRITE(out_unit,*) "    ...back to MolecCav_Normalize_R2_complex"
    END IF 

  END SUBROUTINE MolecCav_Normalize_R2_complex

  
  SUBROUTINE MolecCav_Normalize_R1_real(Psi, Debug)
    USE QDUtil_m
    IMPLICIT NONE

    real(kind=Rkind),  intent(inout) :: Psi(:)                                                          ! already allocated
    logical, optional, intent(in)    :: Debug

    real(kind=Rkind)                 :: Norm
    real(kind=Rkind), parameter      :: Threshold = 1E-10_Rkind
    logical                          :: Debug_local

    !--- Debugging options --------------------------------
    IF (PRESENT(Debug)) THEN; Debug_local = Debug
    ELSE; Debug_local = .FALSE.; END IF

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "o Arguments of MolecCav_Normalize_R1_real :"
      WRITE(out_unit,*) "The <<Psi>> argument                 :"
      CALL Write_Vec(Psi, out_unit, 1, info="Psi")
      WRITE(out_unit,*) "The <<Threshold>> internal parameter :"//TO_string(Threshold)
      FLUSH(out_unit)
    END IF

    IF (Debug_local) WRITE(out_unit,*) "--- Computing norm of Psi..."
    CALL Norm_of(Norm, Psi, Debug_local)                                                                ! the printing of the norm according to Debug is within the call of Norm_of
    IF (Debug_local) WRITE(out_unit,*) "    ...back to MolecCav_Normalize_R1_real"

    IF (Norm < Threshold) THEN
      WRITE(out_unit,*) "### Attempt to normalize matrix of norm ZERO at MolecCav_Normalize_R1_real."
      STOP "### Attempt to normalize matrix of norm ZERO at MolecCav_Normalize_R1_real."
      
    !--- Normalisation ------------------------------------
    ELSE 
      Psi(:) = Psi(:) / Norm
    END IF 

    !--- Conclusion ---------------------------------------
    IF (Debug_local) THEN
      WRITE(out_unit,*) "--- The normalized Psi :"
      CALL Write_Vec(Psi, out_unit, 1, info="Normalised Psi")
      WRITE(out_unit,*) "--- Computing the new norm of Psi..."
      CALL Norm_of(Norm, Psi, Debug_local)                                                              ! the printing of the norm according to Debug is within the call of Norm_of
      WRITE(out_unit,*) "    ...back to MolecCav_Normalize_R1_real"
    END IF 

  END SUBROUTINE MolecCav_Normalize_R1_real


  SUBROUTINE MolecCav_Normalize_R1_complex(Psi, Debug)
    USE QDUtil_m
    IMPLICIT NONE

    complex(kind=Rkind), intent(inout) :: Psi(:)                                                        ! already allocated
    logical, optional, intent(in)      :: Debug

    real(kind=Rkind)                   :: Norm
    real(kind=Rkind), parameter        :: Threshold = 1E-10_Rkind
    logical                            :: Debug_local

    !--- Debugging options --------------------------------
    IF (PRESENT(Debug)) THEN; Debug_local = Debug
    ELSE; Debug_local = .FALSE.; END IF

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "o Arguments of MolecCav_Normalize_R1_complex :"
      WRITE(out_unit,*) "The <<Psi>> argument                 :"
      CALL Write_Vec(Psi, out_unit, 1, info="Psi")
      WRITE(out_unit,*) "The <<Threshold>> internal parameter :"//TO_string(Threshold)
      FLUSH(out_unit)
    END IF

    IF (Debug_local) WRITE(out_unit,*) "--- Computing norm of Psi..."
    CALL Norm_of(Norm, Psi, Debug_local)                                                                ! the printing of the norm according to Debug is within the call of Norm_of
    IF (Debug_local) WRITE(out_unit,*) "    ...back to MolecCav_Normalize_R1_complex"

    IF (Norm < Threshold) THEN
      WRITE(out_unit,*) "### Attempt to normalize matrix of norm ZERO at MolecCav_Normalize_R1_complex."
      STOP "### Attempt to normalize matrix of norm ZERO at MolecCav_Normalize_R1_complex."
      
    !--- Normalisation ------------------------------------
    ELSE 
      Psi(:) = Psi(:) / Norm
    END IF 

    !--- Conclusion ---------------------------------------
    IF (Debug_local) THEN
      WRITE(out_unit,*) "--- The normalized Psi :"
      CALL Write_Vec(Psi, out_unit, 1, info="Normalised Psi")
      WRITE(out_unit,*) "--- Computing the new norm of Psi..."
      CALL Norm_of(Norm, Psi, Debug_local)                                                              ! the printing of the norm according to Debug is within the call of Norm_of
      WRITE(out_unit,*) "    ...back to MolecCav_Normalize_R1_complex"
    END IF 

  END SUBROUTINE MolecCav_Normalize_R1_complex


  SUBROUTINE MolecCav_Norm_R2_real(Norm, Psi, Debug)
    USE QDUtil_m
    IMPLICIT NONE
  
    real(kind=Rkind),  intent(inout) :: Norm
    real(kind=Rkind),  intent(in)    :: Psi(:,:)                                                        ! already allocated
    logical, optional, intent(in)    :: Debug

    logical                          :: Debug_local

    !--- Debugging options --------------------------------
    IF (PRESENT(Debug)) THEN; Debug_local = Debug
    ELSE; Debug_local = .FALSE.; END IF

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "o Arguments of MolecCav_Norm_R2_real :"
      WRITE(out_unit,*) "The <<Norm>> argument :"//TO_string(Norm)
      WRITE(out_unit,*) "The <<Psi>> argument  :"
      CALL Write_Mat(Psi, out_unit, SIZE(Psi, dim=2), info="Psi")
      FLUSH(out_unit)
    END IF

    !--- Norm computation ---------------------------------
    IF (Debug_local) WRITE(out_unit,*) "--- Computing squared modulus of Psi..."
    CALL Scalar_product(Norm, Psi, Psi, Debug_local)
    IF (Debug_local) WRITE(out_unit,*) "    ...back to MolecCav_Norm_R2_real"

    Norm = SQRT(Norm)
  
    IF (Debug_local) WRITE(out_unit,*) "--- The computed norm of Psi :"//TO_string(Norm)

  END SUBROUTINE MolecCav_Norm_R2_real


  SUBROUTINE MolecCav_Norm_R2_complex(Norm, Psi, Debug)
    USE QDUtil_m
    IMPLICIT NONE
  
    real(kind=Rkind),    intent(inout) :: Norm
    complex(kind=Rkind), intent(in)    :: Psi(:,:)                                                      ! already allocated
    logical, optional,   intent(in)    :: Debug

    complex(kind=Rkind)                :: ScaP
    logical                            :: Debug_local

    !--- Debugging options --------------------------------
    IF (PRESENT(Debug)) THEN; Debug_local = Debug
    ELSE; Debug_local = .FALSE.; END IF

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "o Arguments of MolecCav_Norm_R2_complex :"
      WRITE(out_unit,*) "The <<Norm>> argument :"//TO_string(Norm)
      WRITE(out_unit,*) "The <<Psi>> argument  :"
      CALL Write_Mat(Psi, out_unit, SIZE(Psi, dim=2), info="Psi")
      FLUSH(out_unit)
    END IF

    !--- Norm computation ---------------------------------
    IF (Debug_local) WRITE(out_unit,*) "--- Computing squared modulus of Psi..."
    CALL Scalar_product(ScaP, Psi, Psi, Debug_local)
    IF (Debug_local) WRITE(out_unit,*) "    ...back to MolecCav_Norm_R2_complex"

    Norm = SQRT(REAL(ScaP, kind=Rkind))
  
    !--- Conclusion ---------------------------------------
    IF (Debug_local) WRITE(out_unit,*) "--- The computed norm of Psi :"//TO_string(Norm)

  END SUBROUTINE MolecCav_Norm_R2_complex


  SUBROUTINE MolecCav_Norm_R1_real(Norm, Psi, Debug)
    USE QDUtil_m
    IMPLICIT NONE
  
    real(kind=Rkind),  intent(inout) :: Norm
    real(kind=Rkind),  intent(in)    :: Psi(:)                                                          ! already allocated
    logical, optional, intent(in)    :: Debug

    logical                          :: Debug_local

    !--- Debugging options --------------------------------
    IF (PRESENT(Debug)) THEN; Debug_local = Debug
    ELSE; Debug_local = .FALSE.; END IF

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "o Arguments of MolecCav_Norm_R1_real :"
      WRITE(out_unit,*) "The <<Norm>> argument :"//TO_string(Norm)
      WRITE(out_unit,*) "The <<Psi>> argument  :"
      CALL Write_Vec(Psi, out_unit, 1, info="Psi")
      FLUSH(out_unit)
    END IF

    !--- Norm computation ---------------------------------
    IF (Debug_local) WRITE(out_unit,*) "--- Computing squared modulus of Psi..."
    CALL Scalar_product(Norm, Psi, Psi, Debug_local)
    IF (Debug_local) WRITE(out_unit,*) "    ...back to MolecCav_Norm_R1_real"

    Norm = SQRT(Norm)
  
    !--- Conclusion ---------------------------------------
    IF (Debug_local) WRITE(out_unit,*) "--- The computed norm of Psi :"//TO_string(Norm)

  END SUBROUTINE MolecCav_Norm_R1_real


  SUBROUTINE MolecCav_Norm_R1_complex(Norm, Psi, Debug)
    USE QDUtil_m
    IMPLICIT NONE
  
    real(kind=Rkind),    intent(inout) :: Norm
    complex(kind=Rkind), intent(in)    :: Psi(:)                                                        ! already allocated
    logical, optional,   intent(in)    :: Debug

    complex(kind=Rkind)                :: ScaP
    logical                            :: Debug_local

    !--- Debugging options --------------------------------
    IF (PRESENT(Debug)) THEN; Debug_local = Debug
    ELSE; Debug_local = .FALSE.; END IF

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "o Arguments of MolecCav_Norm_R1_complex :"
      WRITE(out_unit,*) "The <<Norm>> argument :"//TO_string(Norm)
      WRITE(out_unit,*) "The <<Psi>> argument  :"
      CALL Write_Vec(Psi, out_unit, 1, info="Psi")
      FLUSH(out_unit)
    END IF

    !--- Norm computation ---------------------------------
    IF (Debug_local) WRITE(out_unit,*) "--- Computing squared modulus of Psi..."
    CALL Scalar_product(ScaP, Psi, Psi, Debug_local)
    IF (Debug_local) WRITE(out_unit,*) "    ...back to MolecCav_Norm_R1_complex"

    Norm = SQRT(REAL(ScaP, kind=Rkind))
  
    !--- Conclusion ---------------------------------------
    IF (Debug_local) WRITE(out_unit,*) "--- The computed norm of Psi :"//TO_string(Norm)

  END SUBROUTINE MolecCav_Norm_R1_complex


  SUBROUTINE MolecCav_Scalar_product_R2_real(ScaP, Psi_1, Psi_2, Debug)
    USE QDUtil_m
    IMPLICIT NONE
  
    real(kind=Rkind),  intent(inout) :: ScaP
    real(kind=Rkind),  intent(in)    :: Psi_1(:,:)                                                      ! already allocated
    real(kind=Rkind),  intent(in)    :: Psi_2(:,:)                                                      ! already allocated
    logical, optional, intent(in)    :: Debug

    integer                          :: Dim, i_2
    logical                          :: Debug_local

    !--- Debugging options --------------------------------
    IF (PRESENT(Debug)) THEN; Debug_local = Debug
    ELSE; Debug_local = .FALSE.; END IF

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "o Arguments of MolecCav_Scalar_product_R2_real :"
      WRITE(out_unit,*) "The <<ScaP>> argument  :"//TO_string(ScaP)
      WRITE(out_unit,*) "The <<Psi_1>> argument :"
      CALL Write_Mat(Psi_1, out_unit, SIZE(Psi_1, dim=2), info="Psi_1")
      WRITE(out_unit,*) "The <<Psi_2>> argument :"
      CALL Write_Mat(Psi_2, out_unit, SIZE(Psi_2, dim=2), info="Psi_2")
      FLUSH(out_unit)
    END IF

    IF (SIZE(Psi_1, Dim=2) /= SIZE(Psi_2, Dim=2) .OR. SIZE(Psi_1, Dim=1) /= SIZE(Psi_2, Dim=1)) THEN
      WRITE(out_unit,*) "### The matrices are expected to have same sizes at MolecCav_Scalar_product_R2_real."
      WRITE(out_unit,*) "    SIZE(Psi_1, dim=1) = "//TO_string(SIZE(Psi_1, dim=1))//"; SIZE(Psi_1, dim=2) = "//TO_string(SIZE(Psi&
      &_1, dim=2))
      WRITE(out_unit,*) "    SIZE(Psi_2, dim=1) = "//TO_string(SIZE(Psi_2, dim=1))//"; SIZE(Psi_2, dim=2) = "//TO_string(SIZE(Psi&
      &_2, dim=2))
      STOP "### The matrices are expected to have same sizes MolecCav_Scalar_product_R2_real."
    END IF

    !--- Scalar product computation -----------------------
    Dim = SIZE(Psi_1, Dim=2)

    ScaP = ZERO
    DO i_2 = 1, Dim
      ScaP = ScaP + DOT_PRODUCT(Psi_1(:,i_2), Psi_2(:,i_2))
    END DO

    !--- Conclusion ---------------------------------------
    IF (Debug_local) WRITE(out_unit,*) "--- The computed scalar product < Psi_1 | Psi_2 >  ="//TO_string(ScaP)

  END SUBROUTINE MolecCav_Scalar_product_R2_real


  SUBROUTINE MolecCav_Scalar_product_R2_complex(ScaP, Psi_1, Psi_2, Debug)
    USE QDUtil_m
    IMPLICIT NONE
  
    complex(kind=Rkind),  intent(inout) :: ScaP
    complex(kind=Rkind),  intent(in)    :: Psi_1(:,:)                                                   ! already allocated
    complex(kind=Rkind),  intent(in)    :: Psi_2(:,:)                                                   ! already allocated
    logical, optional, intent(in)       :: Debug

    integer                             :: Dim, i_2
    logical                             :: Debug_local

    !--- Debugging options --------------------------------
    IF (PRESENT(Debug)) THEN; Debug_local = Debug
    ELSE; Debug_local = .FALSE.; END IF

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "o Arguments of MolecCav_Scalar_product_R2_complex :"
      WRITE(out_unit,*) "The <<ScaP>> argument  :"//TO_string(ScaP)
      WRITE(out_unit,*) "The <<Psi_1>> argument :"
      CALL Write_Mat(Psi_1, out_unit, SIZE(Psi_1, dim=2), info="Psi_1")
      WRITE(out_unit,*) "The <<Psi_2>> argument :"
      CALL Write_Mat(Psi_2, out_unit, SIZE(Psi_2, dim=2), info="Psi_2")
      FLUSH(out_unit)
    END IF

    IF (SIZE(Psi_1, Dim=2) /= SIZE(Psi_2, Dim=2) .OR. SIZE(Psi_1, Dim=1) /= SIZE(Psi_2, Dim=1)) THEN
      WRITE(out_unit,*) "### The matrices are expected to have same sizes at MolecCav_Scalar_product_R2_complex."
      WRITE(out_unit,*) "    SIZE(Psi_1, dim=1) = "//TO_string(SIZE(Psi_1, dim=1))//"; SIZE(Psi_1, dim=2) = "//TO_string(SIZE(Psi&
      &_1, dim=2))
      WRITE(out_unit,*) "    SIZE(Psi_2, dim=1) = "//TO_string(SIZE(Psi_2, dim=1))//"; SIZE(Psi_2, dim=2) = "//TO_string(SIZE(Psi&
      &_2, dim=2))
      STOP "### The matrices are expected to have same sizes MolecCav_Scalar_product_R2_complex."
    END IF

    !--- Scalar product computation -----------------------
    Dim = SIZE(Psi_1, Dim=2)

    ScaP = ZERO
    DO i_2 = 1, Dim
      ScaP = ScaP + DOT_PRODUCT(Psi_1(:,i_2), Psi_2(:,i_2))
    END DO

    !--- Conclusion ---------------------------------------
    IF (Debug_local) WRITE(out_unit,*) "--- The computed scalar product < Psi_1 | Psi_2 >  ="//TO_string(ScaP)

  END SUBROUTINE MolecCav_Scalar_product_R2_complex


  SUBROUTINE MolecCav_Scalar_product_R1_real(ScaP, Psi_1, Psi_2, Debug)
    USE QDUtil_m
    IMPLICIT NONE
  
    real(kind=Rkind),  intent(inout) :: ScaP
    real(kind=Rkind),  intent(in)    :: Psi_1(:)                                                        ! already allocated
    real(kind=Rkind),  intent(in)    :: Psi_2(:)                                                        ! already allocated
    logical, optional, intent(in)    :: Debug

    logical                          :: Debug_local

    !--- Debugging options --------------------------------
    IF (PRESENT(Debug)) THEN; Debug_local = Debug
    ELSE; Debug_local = .FALSE.; END IF

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "o Arguments of MolecCav_Scalar_product_R1_real :"
      WRITE(out_unit,*) "The <<ScaP>> argument  :"//TO_string(ScaP)
      WRITE(out_unit,*) "The <<Psi_1>> argument :"
      CALL Write_Vec(Psi_1, out_unit, 1, info="Psi_1")
      WRITE(out_unit,*) "The <<Psi_2>> argument :"
      CALL Write_Vec(Psi_2, out_unit, 1, info="Psi_2")
      FLUSH(out_unit)
    END IF

    IF (SIZE(Psi_1, Dim=1) /= SIZE(Psi_2, Dim=1)) THEN
      WRITE(out_unit,*) "### The matrices are expected to have same sizes at MolecCav_Scalar_product_R1_real."
      WRITE(out_unit,*) "    SIZE(Psi_1, dim=1) = "//TO_string(SIZE(Psi_1, dim=1))//"; SIZE(Psi_2, dim=1) = "//TO_string(SIZE(Psi&
      &_2, dim=1))
      STOP "### The matrices are expected to have same sizes MolecCav_Scalar_product_R1_real."
    END IF

    !--- Scalar product computation -----------------------
    ScaP = DOT_PRODUCT(Psi_1, Psi_2)

    !--- Conclusion ---------------------------------------
    IF (Debug_local) WRITE(out_unit,*) "--- The computed scalar product < Psi_1 | Psi_2 >  ="//TO_string(ScaP)

  END SUBROUTINE MolecCav_Scalar_product_R1_real


  SUBROUTINE MolecCav_Scalar_product_R1_complex(ScaP, Psi_1, Psi_2, Debug)
    USE QDUtil_m
    IMPLICIT NONE
  
    complex(kind=Rkind),  intent(inout) :: ScaP
    complex(kind=Rkind),  intent(in)    :: Psi_1(:)                                                     ! already allocated
    complex(kind=Rkind),  intent(in)    :: Psi_2(:)                                                     ! already allocated
    logical, optional, intent(in)       :: Debug

    logical                             :: Debug_local

    !--- Debugging options --------------------------------
    IF (PRESENT(Debug)) THEN; Debug_local = Debug
    ELSE; Debug_local = .FALSE.; END IF

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "o Arguments of MolecCav_Scalar_product_R1_complex :"
      WRITE(out_unit,*) "The <<ScaP>> argument  :"//TO_string(ScaP)
      WRITE(out_unit,*) "The <<Psi_1>> argument :"
      CALL Write_Vec(Psi_1, out_unit, 1, info="Psi_1")
      WRITE(out_unit,*) "The <<Psi_2>> argument :"
      CALL Write_Vec(Psi_2, out_unit, 1, info="Psi_2")
      FLUSH(out_unit)
    END IF

    IF (SIZE(Psi_1, Dim=1) /= SIZE(Psi_2, Dim=1)) THEN
      WRITE(out_unit,*) "### The matrices are expected to have same sizes at MolecCav_Scalar_product_R1_complex."
      WRITE(out_unit,*) "    SIZE(Psi_1, dim=1) = "//TO_string(SIZE(Psi_1, dim=1))//"; SIZE(Psi_2, dim=1) = "//TO_string(SIZE(Psi&
      &_2, dim=1))
      STOP "### The matrices are expected to have same sizes MolecCav_Scalar_product_R1_complex."
    END IF

    !--- Scalar product computation -----------------------
    ScaP = DOT_PRODUCT(Psi_1, Psi_2)

    !--- Conclusion ---------------------------------------
    IF (Debug_local) WRITE(out_unit,*) "--- The computed scalar product < Psi_1 | Psi_2 >  ="//TO_string(ScaP)

  END SUBROUTINE MolecCav_Scalar_product_R1_complex


END MODULE
  