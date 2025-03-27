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
!==================================================================================================
PROGRAM test_algebra
  !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
  USE QDUtil_m
  USE QDUtil_Test_m
  USE Algebra_m
  IMPLICIT NONE
  
  logical, parameter               :: Debug = .FALSE.

  integer                          :: Dim_1 = 5
  integer                          :: Dim_2 = 6

  real(kind=Rkind),    allocatable :: Basis_R_1D(:)
  complex(kind=Rkind), allocatable :: Basis_C_1D(:)
  real(kind=Rkind),    allocatable :: Basis_R_2D(:,:)
  complex(kind=Rkind), allocatable :: Basis_C_2D(:,:)

  integer                          :: Coeff_1 = 1
  integer                          :: Coeff_2 = 1
  integer                          :: Coeff_3 = 4                                                                                ! for the linear combinations
  real(kind=Rkind),    allocatable :: Any_R_1D(:)
  complex(kind=Rkind), allocatable :: Any_C_1D(:)
  real(kind=Rkind),    allocatable :: Any_R_2D(:,:)
  complex(kind=Rkind), allocatable :: Any_C_2D(:,:)
  real(kind=Rkind)                 :: Normalization_coeff

  real(kind=Rkind)                 :: Sca_pdt_R
  complex(kind=Rkind)              :: Sca_pdt_C
  real(kind=Rkind)                 :: Norm

  real(kind=Rkind),    allocatable :: R_vec(:)                                                               ! Buffer
  real(kind=Rkind),    allocatable :: R_mat(:,:)                                                             ! Buffer
  complex(kind=Rkind), allocatable :: C_vec(:)                                                               ! Buffer
  complex(kind=Rkind), allocatable :: C_mat(:,:)                                                             ! Buffer

  TYPE(test_t)                     :: test_sca_pdt
  TYPE(test_t)                     :: test_norm
  TYPE(test_t)                     :: test_normalization  
  logical                          :: error_sca_pdt       = .FALSE.
  logical                          :: error_norm          = .FALSE.
  logical                          :: error_normalization = .FALSE.


  !----------------------------Tests initialization----------------------------
  CALL Initialize_Test(test_sca_pdt,      test_name='OUT/test_file_sca_pdt')
  CALL Initialize_Test(test_norm,         test_name='OUT/test_file_norm')
  CALL Initialize_Test(test_normalization,test_name='OUT/test_file_normalztn')


  !---------------------------Matrices initialization--------------------------
  ALLOCATE(Basis_R_1D(Dim_1))
  ALLOCATE(Basis_C_1D(Dim_1))
  ALLOCATE(Basis_R_2D(Dim_1, Dim_2))
  ALLOCATE(Basis_C_2D(Dim_1, Dim_2))
  Basis_R_1D = ZERO
  Basis_C_1D = ZERO
  Basis_R_2D = ZERO
  Basis_C_2D = ZERO

  ALLOCATE(Any_R_1D(Dim_1))
  ALLOCATE(Any_C_1D(Dim_1))
  ALLOCATE(Any_R_2D(Dim_1, Dim_2))
  ALLOCATE(Any_C_2D(Dim_1, Dim_2))

  ALLOCATE(R_vec(Dim_1))
  ALLOCATE(C_vec(Dim_1))
  ALLOCATE(R_mat(Dim_1, Dim_2))
  ALLOCATE(C_mat(Dim_1, Dim_2))

  Basis_R_1D(1)   = ONE                                                        ! first basis vector. Supposed to have norm 1
  Basis_C_1D(1)   = EYE                                                        ! first basis vector. Supposed to have norm 1
  Basis_R_2D(1,1) = ONE                                                        ! first basis vector. Supposed to have norm 1
  Basis_C_2D(1,1) = EYE                                                        ! first basis vector. Supposed to have norm 1

  Any_R_1D(:)           = Coeff_1*Basis_R_1D(:)
  Any_R_1D(Dim_1/2)     = Coeff_2*ONE
  Any_R_1D(Dim_1/2 + 1) = Coeff_3*ONE                                          ! an arbitrary superposition of 3 basis vectors. Supposed to have norm 18
  Any_C_1D(:)           = Coeff_1*Basis_C_1D(:)
  Any_C_1D(Dim_1/2)     = Coeff_2*ONE
  Any_C_1D(Dim_1/2 + 1) = Coeff_3*EYE                                          ! an arbitrary superposition of 3 basis vectors. Supposed to have norm 18

  Any_R_2D(:,:)                     = Coeff_1*Basis_R_2D(:,:)
  Any_R_2D(1, Dim_2/2)              = Coeff_2*ONE
  Any_R_2D(Dim_1/2 + 1, Dim_2/2 +1) = Coeff_3*ONE                              ! an arbitrary superposition of 3 basis vectors. Supposed to have norm 18
  Any_C_2D(:,:)                     = Coeff_1*Basis_C_2D(:,:)
  Any_C_2D(1, Dim_2/2)              = Coeff_2*ONE
  Any_C_2D(Dim_1/2 + 1, Dim_2/2 +1) = Coeff_3*EYE                              ! an arbitrary superposition of 3 basis vectors. Supposed to have norm 18

  Normalization_coeff = SQRT(REAL((Coeff_1**2 + Coeff_2**2 + Coeff_3**2), kind=Rkind)) ! the theoretical/analytical normalization coefficient

  IF (Debug) THEN
    WRITE(out_unit,*) "---------------------------Matrices initialization--------------------------"
    CALL Write_Vec(Basis_R_1D, out_unit, 1, info="Basis_R_1D")
    WRITE(out_unit,*)
    CALL Write_Vec(Basis_C_1D, out_unit, 1, info="Basis_C_1D")
    WRITE(out_unit,*)
    CALL Write_Mat(Basis_R_2D, out_unit, Dim_2, info="Basis_R_2D")
    WRITE(out_unit,*)
    CALL Write_Mat(Basis_C_2D, out_unit, Dim_2, info="Basis_C_2D")
    WRITE(out_unit,*)
    CALL Write_Vec(Any_R_1D, out_unit, 1, info="Any_R_1D")
    WRITE(out_unit,*)
    CALL Write_Vec(Any_C_1D, out_unit, 1, info="Any_C_1D")
    WRITE(out_unit,*)
    CALL Write_Mat(Any_R_2D, out_unit, Dim_2, info="Any_R_2D")
    WRITE(out_unit,*)
    CALL Write_Mat(Any_C_2D, out_unit, Dim_2, info="Any_C_2D")
    WRITE(out_unit,*)
    WRITE(out_unit,*) "Their theoretical/analytical normalization coefficient :", Normalization_coeff
  END IF

  
  !-------------------Tests of the Scalar_product procedures-------------------
  WRITE(out_unit,*)
  WRITE(out_unit,*) "------------------------The Scalar_product procedures-----------------------"
  
  IF (Debug) WRITE(out_unit,*) "Expecting", 1
  CALL Scalar_product(Sca_pdt_R, Basis_R_1D, Basis_R_1D)
  IF (Debug) WRITE(out_unit,*) "Computed scalar product : < Basis_R_1D | Basis_R_1D >  =", Sca_pdt_R
  CALL Equal_R_R_scalar(error_sca_pdt, Sca_pdt_R, ONE)
  CALL Logical_Test(test_sca_pdt, test1=error_sca_pdt, test2=.FALSE., info="< Basis_R_1D | Basis_R_1D > = 1")

  CALL Scalar_product(Sca_pdt_C, Basis_C_1D, Basis_C_1D)
  IF (Debug) WRITE(out_unit,*)
  IF (Debug) WRITE(out_unit,*) "Computed scalar product : < Basis_C_1D |  Basis_C_1D >  =", Sca_pdt_C
  CALL Equal_C_C_scalar(error_sca_pdt, Sca_pdt_C, CMPLX(ONE, kind=Rkind))
  CALL Logical_Test(test_sca_pdt, test1=error_sca_pdt, test2=.FALSE., info="< Basis_C_1D | Basis_C_1D > = 1")

  CALL Scalar_product(Sca_pdt_R, Basis_R_2D, Basis_R_2D)
  IF (Debug) WRITE(out_unit,*)
  IF (Debug) WRITE(out_unit,*) "Computed scalar product : < Basis_R_2D |  Basis_R_2D >  =", Sca_pdt_R
  CALL Equal_R_R_scalar(error_sca_pdt, Sca_pdt_R, ONE)
  CALL Logical_Test(test_sca_pdt, test1=error_sca_pdt, test2=.FALSE., info="< Basis_R_2D | Basis_R_2D > = 1")

  CALL Scalar_product(Sca_pdt_C, Basis_C_2D, Basis_C_2D)
  IF (Debug) WRITE(out_unit,*)
  IF (Debug) WRITE(out_unit,*) "Computed scalar product : < Basis_C_2D |  Basis_C_2D >  =", Sca_pdt_C
  CALL Equal_C_C_scalar(error_sca_pdt, Sca_pdt_C, CMPLX(ONE, kind=Rkind))
  CALL Logical_Test(test_sca_pdt, test1=error_sca_pdt, test2=.FALSE., info="< Basis_C_2D | Basis_C_2D > = 1")


  IF (Debug) WRITE(out_unit,*)
  IF (Debug) WRITE(out_unit,*) "Expecting", Coeff_1**2 + Coeff_2**2 + Coeff_3**2
  CALL Scalar_product(Sca_pdt_R, Any_R_1D, Any_R_1D)
  IF (Debug) WRITE(out_unit,*) "Computed scalar product : < Any_R_1D | Any_R_1D >  =", Sca_pdt_R
  CALL Equal_R_R_scalar(error_sca_pdt, Sca_pdt_R, 18.0_Rkind)
  CALL Logical_Test(test_sca_pdt, test1=error_sca_pdt, test2=.FALSE., info="< Any_R_1D | Any_R_1D > = 18")

  CALL Scalar_product(Sca_pdt_C, Any_C_1D, Any_C_1D)
  IF (Debug) WRITE(out_unit,*)
  IF (Debug) WRITE(out_unit,*) "Computed scalar product : < Any_C_1D | Any_C_1D >  =", Sca_pdt_C
  CALL Equal_C_C_scalar(error_sca_pdt, Sca_pdt_C, CMPLX(18.0_Rkind, kind=Rkind))
  CALL Logical_Test(test_sca_pdt, test1=error_sca_pdt, test2=.FALSE., info="< Any_C_1D | Any_C_1D > = 18")

  CALL Scalar_product(Sca_pdt_R, Any_R_2D, Any_R_2D)
  IF (Debug) WRITE(out_unit,*)
  IF (Debug) WRITE(out_unit,*) "Computed scalar product : < Any_R_2D | Any_R_2D >  =", Sca_pdt_R
  CALL Equal_R_R_scalar(error_sca_pdt, Sca_pdt_R, 18.0_Rkind)
  CALL Logical_Test(test_sca_pdt, test1=error_sca_pdt, test2=.FALSE., info=" < Any_R_2D | Any_R_2D > = 18")

  CALL Scalar_product(Sca_pdt_C, Any_C_2D, Any_C_2D)
  IF (Debug) WRITE(out_unit,*)
  IF (Debug) WRITE(out_unit,*) "Computed scalar product : < Any_C_2D | Any_C_2D >  =", Sca_pdt_C
  CALL Equal_C_C_scalar(error_sca_pdt, Sca_pdt_C, CMPLX(18.0_Rkind, kind=Rkind))
  CALL Logical_Test(test_sca_pdt, test1=error_sca_pdt, test2=.FALSE., info="< Any_C_2D | Any_C_2D > = 18")


  IF (Debug) WRITE(out_unit,*)
  IF (Debug) WRITE(out_unit,*) "Expecting", 1, "or", "   i"
  CALL Scalar_product(Sca_pdt_R, Basis_R_1D, Any_R_1D)
  IF (Debug) WRITE(out_unit,*) "Computed scalar product : < Basis_R_1D | Any_R_1D >  =", Sca_pdt_R
  CALL Equal_R_R_scalar(error_sca_pdt, Sca_pdt_R, ONE)
  CALL Logical_Test(test_sca_pdt, test1=error_sca_pdt, test2=.FALSE., info="< Basis_R_1D | Any_R_1D > = 1")

  CALL Scalar_product(Sca_pdt_C, Basis_C_1D, Any_C_1D)
  IF (Debug) WRITE(out_unit,*)
  IF (Debug) WRITE(out_unit,*) "Computed scalar product : < Basis_C_1D | Any_C_1D >  =", Sca_pdt_C
  CALL Equal_C_C_scalar(error_sca_pdt, Sca_pdt_C, CMPLX(ONE, kind=Rkind))
  CALL Logical_Test(test_sca_pdt, test1=error_sca_pdt, test2=.FALSE., info="< Basis_C_1D | Any_C_1D > = 1")

  CALL Scalar_product(Sca_pdt_C, CMPLX(Basis_R_1D(:), kind=Rkind), Any_C_1D)
  IF (Debug) WRITE(out_unit,*)
  IF (Debug) WRITE(out_unit,*) "Computed scalar product : < Basis_R_1D | Any_C_1D >  =", Sca_pdt_C
  CALL Equal_C_C_scalar(error_sca_pdt, Sca_pdt_C, CMPLX(ZERO, ONE, kind=Rkind))
  CALL Logical_Test(test_sca_pdt, test1=error_sca_pdt, test2=.FALSE., info="< Basis_R_1D | Any_C_1D > = i")

  CALL Scalar_product(Sca_pdt_R, Basis_R_2D, Any_R_2D)
  IF (Debug) WRITE(out_unit,*)
  IF (Debug) WRITE(out_unit,*) "Computed scalar product : < Basis_R_2D | Any_R_2D >  =", Sca_pdt_R
  CALL Equal_R_R_scalar(error_sca_pdt, Sca_pdt_R, ONE)
  CALL Logical_Test(test_sca_pdt, test1=error_sca_pdt, test2=.FALSE., info="< Basis_R_2D | Any_R_2D > = 1")

  CALL Scalar_product(Sca_pdt_C, Basis_C_2D, Any_C_2D)
  IF (Debug) WRITE(out_unit,*)
  IF (Debug) WRITE(out_unit,*) "Computed scalar product : < Basis_C_2D | Any_C_2D >  =", Sca_pdt_C
  CALL Equal_C_C_scalar(error_sca_pdt, Sca_pdt_C, CMPLX(ONE, kind=Rkind))
  CALL Logical_Test(test_sca_pdt, test1=error_sca_pdt, test2=.FALSE., info="< Basis_C_2D | Any_C_2D > = 1")

  CALL Scalar_product(Sca_pdt_C, CMPLX(Basis_R_2D(:,:), kind=Rkind), Any_C_2D)
  IF (Debug) WRITE(out_unit,*)
  IF (Debug) WRITE(out_unit,*) "Computed scalar product : < Basis_R_2D | Any_C_2D >  =", Sca_pdt_C
  CALL Equal_C_C_scalar(error_sca_pdt, Sca_pdt_C, CMPLX(ZERO, ONE, kind=Rkind))
  CALL Logical_Test(test_sca_pdt, test1=error_sca_pdt, test2=.FALSE., info="< Basis_R_2D | Any_C_2D > = i")


  !------------------------Tests of the Norm computation-----------------------
  WRITE(out_unit,*)
  WRITE(out_unit,*) "----------------------------The Norm computation----------------------------"

  CALL Norm_of(Norm, Basis_R_1D)
  IF (Debug) WRITE(out_unit,*) "Computed norm : ||Basis_R_1D|| =", Norm
  CALL Equal_R_R_scalar(error_norm, Norm, ONE)
  CALL Logical_Test(test_norm, test1=error_norm, test2=.FALSE., info="||Basis_R_1D|| = 1")

  CALL Norm_of(Norm, Basis_C_1D)
  IF (Debug) WRITE(out_unit,*)
  IF (Debug) WRITE(out_unit,*) "Computed norm : ||Basis_C_1D|| =", Norm
  CALL Equal_R_R_scalar(error_norm, Norm, ONE)
  CALL Logical_Test(test_norm, test1=error_norm, test2=.FALSE., info="||Basis_C_1D|| = 1")

  CALL Norm_of(Norm, Basis_R_2D)
  IF (Debug) WRITE(out_unit,*)
  IF (Debug) WRITE(out_unit,*) "Computed norm : ||Basis_R_2D|| =", Norm
  CALL Equal_R_R_scalar(error_norm, Norm, ONE)
  CALL Logical_Test(test_norm, test1=error_norm, test2=.FALSE., info="||Basis_R_2D|| = 1")

  CALL Norm_of(Norm, Basis_C_2D)
  IF (Debug) WRITE(out_unit,*)
  IF (Debug) WRITE(out_unit,*) "Computed norm : ||Basis_C_2D|| =", Norm
  CALL Equal_R_R_scalar(error_norm, Norm, ONE)
  CALL Logical_Test(test_norm, test1=error_norm, test2=.FALSE., info="||Basis_C_2D|| = 1")

  CALL Norm_of(Norm, Any_R_1D)
  IF (Debug) WRITE(out_unit,*)
  IF (Debug) WRITE(out_unit,*) "Computed norm : ||Any_R_1D|| =", Norm
  CALL Equal_R_R_scalar(error_norm, Norm, Normalization_coeff)
  CALL Logical_Test(test_norm, test1=error_norm, test2=.FALSE., info="||Any_R_1D|| = 18")

  CALL Norm_of(Norm, Any_C_1D)
  IF (Debug) WRITE(out_unit,*)
  IF (Debug) WRITE(out_unit,*) "Computed norm : ||Any_C_1D|| =", Norm
  CALL Equal_R_R_scalar(error_norm, Norm, Normalization_coeff)
  CALL Logical_Test(test_norm, test1=error_norm, test2=.FALSE., info="||Any_C_1D|| = 18")

  CALL Norm_of(Norm, Any_R_2D)
  IF (Debug) WRITE(out_unit,*)
  IF (Debug) WRITE(out_unit,*) "Computed norm : ||Any_R_2D|| =", Norm
  CALL Equal_R_R_scalar(error_norm, Norm, Normalization_coeff)
  CALL Logical_Test(test_norm, test1=error_norm, test2=.FALSE., info="||Any_R_2D|| = 18")

  CALL Norm_of(Norm, Any_C_2D)
  IF (Debug) WRITE(out_unit,*)
  IF (Debug) WRITE(out_unit,*) "Computed norm : ||Any_C_2D|| =", Norm
  CALL Equal_R_R_scalar(error_norm, Norm, Normalization_coeff)
  CALL Logical_Test(test_norm, test1=error_norm, test2=.FALSE., info="||Any_C_2D|| = 18")


  !--------------------Tests of the Normalization procedures-------------------
  WRITE(out_unit,*)
  WRITE(out_unit,*) "------------------------The Normalization procedures------------------------"

  R_vec(:) = Basis_R_1D(:)
  CALL Normalize(Basis_R_1D)
  CALL Norm_of(Norm, Basis_R_1D)
  IF (Debug) CALL Write_Vec(Basis_R_1D, out_unit, Dim_1, info="Normalized_Basis_R_1D")
  CALL Equal_R_R_scalar(error_normalization, Norm, ONE)
  CALL Logical_Test(test_normalization, test1=error_normalization, test2=.FALSE., info="||Normalized_Basis_R_1D|| = 1")
  CALL Equal_R_R_vector(error_normalization, Basis_R_1D, R_vec)
  CALL Logical_Test(test_normalization, test1=error_normalization, test2=.FALSE., info="Normalized_Basis_R_1D normalized ?")

  C_vec(:) = Basis_C_1D(:)
  CALL Normalize(Basis_C_1D)
  CALL Norm_of(Norm, Basis_C_1D)
  IF (Debug) WRITE(out_unit,*)
  IF (Debug) CALL Write_Vec(Basis_C_1D, out_unit, Dim_1, info="Normalized_Basis_C_1D")
  CALL Equal_R_R_scalar(error_normalization, Norm, ONE)
  CALL Logical_Test(test_normalization, test1=error_normalization, test2=.FALSE., info="||Normalized_Basis_C_1D|| = 1")
  CALL Equal_C_C_vector(error_normalization, Basis_C_1D, C_vec)
  CALL Logical_Test(test_normalization, test1=error_normalization, test2=.FALSE., info="Normalized_Basis_C_1D normalized ?")

  R_mat(:,:) = Basis_R_2D(:,:)
  CALL Normalize(Basis_R_2D)
  CALL Norm_of(Norm, Basis_R_2D)
  IF (Debug) WRITE(out_unit,*)
  IF (Debug) CALL Write_Mat(Basis_R_2D, out_unit, Dim_2, info="Normalized_Basis_R_2D")
  CALL Equal_R_R_scalar(error_normalization, Norm, ONE)
  CALL Logical_Test(test_normalization, test1=error_normalization, test2=.FALSE., info="||Normalized_Basis_R_2D|| = 1")
  CALL Equal_R_R_matrix(error_normalization, Basis_R_2D, R_mat)
  CALL Logical_Test(test_normalization, test1=error_normalization, test2=.FALSE., info="Normalized_Basis_R_2D normalized ?")

  C_mat(:,:) = Basis_C_2D(:,:)
  CALL Normalize(Basis_C_2D)
  CALL Norm_of(Norm, Basis_C_2D)
  IF (Debug) WRITE(out_unit,*)
  IF (Debug) CALL Write_Mat(Basis_C_2D, out_unit, Dim_2, info="Normalized_Basis_C_2D")
  CALL Equal_R_R_scalar(error_normalization, Norm, ONE)
  CALL Logical_Test(test_normalization, test1=error_normalization, test2=.FALSE., info="||Normalized_Basis_C_2D|| = 1")
  CALL Equal_C_C_matrix(error_normalization, Basis_C_2D, C_mat)
  CALL Logical_Test(test_normalization, test1=error_normalization, test2=.FALSE., info="Normalized_Basis_C_2D normalized ?")

  R_vec(:) = Any_R_1D(:)
  CALL Normalize(Any_R_1D)
  CALL Norm_of(Norm, Any_R_1D)
  IF (Debug) WRITE(out_unit,*)
  IF (Debug) CALL Write_Vec(Any_R_1D, out_unit, Dim_1, info="Normalized_Any_R_1D")
  CALL Equal_R_R_scalar(error_normalization, Norm, ONE)
  CALL Logical_Test(test_normalization, test1=error_normalization, test2=.FALSE., info="||Normalized_Any_R_1D|| = 1")
  CALL Equal_R_R_vector(error_normalization, Normalization_coeff*Any_R_1D, R_vec)
  CALL Logical_Test(test_normalization, test1=error_normalization, test2=.FALSE., info="Normalized_Any_R_1D normalized ?")

  C_vec(:) = Any_C_1D(:)
  CALL Normalize(Any_C_1D)
  CALL Norm_of(Norm, Any_C_1D)
  IF (Debug) WRITE(out_unit,*)
  IF (Debug) CALL Write_Vec(Any_C_1D, out_unit, Dim_1, info="Normalized_Any_C_1D")
  CALL Equal_R_R_scalar(error_normalization, Norm, ONE)
  CALL Logical_Test(test_normalization, test1=error_normalization, test2=.FALSE., info="||Normalized_Any_C_1D|| = 1")
  CALL Equal_C_C_vector(error_normalization, Normalization_coeff*Any_C_1D, C_vec)
  CALL Logical_Test(test_normalization, test1=error_normalization, test2=.FALSE., info="Normalized_Any_C_1D normalized ?")

  R_mat(:,:) = Any_R_2D(:,:)
  CALL Normalize(Any_R_2D)
  CALL Norm_of(Norm, Any_R_2D)
  IF (Debug) WRITE(out_unit,*)
  IF (Debug) CALL Write_Mat(Any_R_2D, out_unit, Dim_2, info="Normalized_Any_R_2D")
  CALL Equal_R_R_scalar(error_normalization, Norm, ONE)
  CALL Logical_Test(test_normalization, test1=error_normalization, test2=.FALSE., info="||Normalized_Any_R_2D|| = 1")
  CALL Equal_R_R_matrix(error_normalization, Normalization_coeff*Any_R_2D, R_mat)
  CALL Logical_Test(test_normalization, test1=error_normalization, test2=.FALSE., info="Normalized_Any_R_2D normalized ?")

  C_mat(:,:) = Any_C_2D(:,:)
  CALL Normalize(Any_C_2D)
  !Any_C_2D(4,4) = FOUR
  CALL Norm_of(Norm, Any_C_2D)
  IF (Debug) WRITE(out_unit,*)
  IF (Debug) CALL Write_Mat(Any_C_2D, out_unit, Dim_2, info="Normalized_Any_C_2D")
  CALL Equal_R_R_scalar(error_normalization, Norm, ONE)
  CALL Logical_Test(test_normalization, test1=error_normalization, test2=.FALSE., info="||Normalized_Any_C_2D|| = 1")
  CALL Equal_C_C_matrix(error_normalization, Normalization_coeff*Any_C_2D, C_mat)
  CALL Logical_Test(test_normalization, test1=error_normalization, test2=.FALSE., info="Normalized_Any_C_2D normalized ?")


  !-----------------------------------Sum up-----------------------------------
  WRITE(out_unit,*)
  WRITE(out_unit,*) "-----------------------------------Sum up-----------------------------------"

  CALL Finalize_Test(test_sca_pdt)
  CALL Finalize_Test(test_norm)
  CALL Finalize_Test(test_normalization)

  !IF (test_sca_pdt%nb_Err > 0 .OR. test_norm%nb_Err > 0 .OR. test_normalization%nb_Err > 0) THEN
  !  WRITE(out_unit,*) 
  !  WRITE(out_unit,*) "Test 1 failed ! At least one procedure among the Scalar_product, Norm_of and&
  !                   & Normalization ones does NOT work properly ! Please refer to the test_algebra&
  !                   &.log output file for more information."
  !ELSE
  !  WRITE(out_unit,*) 
  !  WRITE(out_unit,*) "Test 1 checked ! The Scalar_product, Norm_of and Normalization procedures do work properly !"
  !
  !END IF


  CONTAINS


  SUBROUTINE Equal_R_R_scalar(error, Rl_1, Rl_2)
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

  END SUBROUTINE Equal_R_R_scalar
  

  SUBROUTINE Equal_C_C_scalar(error, Cmplx_1, Cmplx_2)
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

  END SUBROUTINE Equal_C_C_scalar
    
  
  SUBROUTINE Equal_R_R_vector(error, Rl_1, Rl_2)
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

  END SUBROUTINE Equal_R_R_vector
  

  SUBROUTINE Equal_C_C_vector(error, Cmplx_1, Cmplx_2)
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

  END SUBROUTINE Equal_C_C_vector
    
  
  SUBROUTINE Equal_R_R_matrix(error, Rl_1, Rl_2)
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

  END SUBROUTINE Equal_R_R_matrix
  

  SUBROUTINE Equal_C_C_matrix(error, Cmplx_1, Cmplx_2)
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

  END SUBROUTINE Equal_C_C_matrix
    
  
END PROGRAM