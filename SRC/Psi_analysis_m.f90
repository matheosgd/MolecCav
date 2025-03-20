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
MODULE Psi_analysis_m
  !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
  USE QDUtil_m
  IMPLICIT NONE

  
  PRIVATE

  PUBLIC Reduced_density_psi

  INTERFACE Reduced_density_psi
    MODULE PROCEDURE MolecCav_Reduced_density_psi_R2_R, MolecCav_Reduced_density_psi_R1_R
  END INTERFACE
    

  CONTAINS
  
    
  SUBROUTINE MolecCav_Reduced_density_psi_R2_R(weight_dim_1, weight_dim_2, Psi_R2, Debug)
    USE QDUtil_m
    USE Mapping_m
    IMPLICIT NONE

    real(kind=Rkind),  intent(inout) :: weight_dim_1(:)                        ! already allocated !
    real(kind=Rkind),  intent(inout) :: weight_dim_2(:)                        ! already allocated !
    real(kind=Rkind),  intent(in)    :: Psi_R2(:,:)                            ! already allocated !
    logical, optional, intent(in)    :: Debug

    logical                          :: Debug_local = .FALSE.
    integer                          :: i

    IF (PRESENT(Debug)) Debug_local = Debug

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "----------Computing the wavefunction Psi_R2's reduced density weights---------"
      WRITE(out_unit,*) "The squared Psi_R2 matrix to be summed :"
      CALL Write_Mat(Psi_R2**2, out_unit, Size(Psi_R2), info="Psi_R2**2")
      WRITE(out_unit,*)
    END IF

    weight_dim_1 = ZERO
    weight_dim_2 = ZERO

    DO i = 1, MAX(Size(Psi_R2, dim=1), Size(Psi_R2, dim=2))
      IF (i<=Size(Psi_R2, dim=1)) THEN
        weight_dim_1(i) = SUM((Psi_R2(i,:))**2) ! sum over i_2
        IF (Debug_local) WRITE(out_unit,*) "weight_dim_1("//TO_string(i)//") = "//TO_string(weight_dim_1(i))
      END IF
      IF (i<=Size(Psi_R2, dim=2)) THEN
        weight_dim_2(i) = SUM((Psi_R2(:,i))**2) ! sum over i_1
        IF (Debug_local) WRITE(out_unit,*) "weight_dim_2("//TO_string(i)//") = "//TO_string(weight_dim_2(i))
      END IF
    END DO

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      CALL Write_Vec(weight_dim_1, out_unit, Size(Psi_R2, dim=1), info="weight_dim_1")
      CALL Write_Vec(weight_dim_2, out_unit, Size(Psi_R2, dim=2), info="weight_dim_2")
      WRITE(out_unit,*) "--------End computing the wavefunction Psi_R2's reduced density weights-------"
    END IF

    END SUBROUTINE MolecCav_Reduced_density_psi_R2_R


    SUBROUTINE MolecCav_Reduced_density_psi_R1_R(weight_dim_1, weight_dim_2, Psi_R1, Debug)
      USE QDUtil_m
      USE Mapping_m
      IMPLICIT NONE
  
      real(kind=Rkind),  intent(inout) :: weight_dim_1(:)                        ! already allocated !
      real(kind=Rkind),  intent(inout) :: weight_dim_2(:)                        ! already allocated !
      real(kind=Rkind),  intent(in)    :: Psi_R1(:)                            ! already allocated !
      logical, optional, intent(in)    :: Debug
  
      logical                          :: Debug_local = .FALSE.
      integer                          :: i
      real(kind=Rkind), allocatable    :: Psi_R2(:,:)
  
      IF (PRESENT(Debug)) Debug_local = Debug
    
      IF (Debug_local) THEN
        WRITE(out_unit,*)
        WRITE(out_unit,*) "----------Computing the wavefunction Psi_R1's reduced density weights---------"
!        WRITE(out_unit,*) "The squared Psi_R1 matrix to be summed :"
!        CALL Write_Vec(Psi_R1**2, out_unit, Size(Psi_R1), info="Psi_R1**2")
        WRITE(out_unit,*)
      END IF
    
      ALLOCATE(Psi_R2(Size(weight_dim_1, dim=1), Size(weight_dim_2, dim=1)))
      CALL Mapping_WF_1DTO2D(Psi_R2, Psi_R1, Debug=.FALSE.)

      CALL Reduced_density_psi(weight_dim_1, weight_dim_2, Psi_R2, Debug=Debug_local)

      IF (Debug_local) THEN
!        WRITE(out_unit,*)
!        CALL Write_Vec(weight_dim_1, out_unit, Size(weight_dim_1, dim=1), info="weight_dim_1")
!        CALL Write_Vec(weight_dim_2, out_unit, Size(weight_dim_2, dim=1), info="weight_dim_2")
        WRITE(out_unit,*) "--------End computing the wavefunction Psi_R1's reduced density weights-------"
      END IF
  
      END SUBROUTINE MolecCav_Reduced_density_psi_R1_R
  
  
  END MODULE
