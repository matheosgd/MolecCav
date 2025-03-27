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
MODULE Psi_analysis_m
  !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
  USE QDUtil_m
  IMPLICIT NONE

  
  PRIVATE

  PUBLIC Reduced_density_psi

  INTERFACE Reduced_density_psi
    MODULE PROCEDURE MolecCav_Reduced_density_psi_1p1D_R2_real, MolecCav_Reduced_density_psi_1p1D_R1_real
  END INTERFACE
    

  CONTAINS
  
    
  SUBROUTINE MolecCav_Reduced_density_psi_1p1D_R2_real(weight_dim_1, weight_dim_2, Psi_R2, Debug)
    USE QDUtil_m
    USE Mapping_m
    IMPLICIT NONE

    real(kind=Rkind),  intent(inout) :: weight_dim_1(:)                                                                          ! already allocated !
    real(kind=Rkind),  intent(inout) :: weight_dim_2(:)                                                                          ! already allocated !
    real(kind=Rkind),  intent(in)    :: Psi_R2(:,:)                                                                              ! already allocated !
    logical, optional, intent(in)    :: Debug

    logical                          :: Debug_local = .FALSE.
    integer                          :: i

    IF (PRESENT(Debug)) Debug_local = Debug

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "------------------------------------Computing the wavefunction Psi_R2's reduced density weights---------&
                        &--------------------------"
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
      WRITE(out_unit,*) "----------------------------------End computing the wavefunction Psi_R2's reduced density weights-------&
                        &--------------------------"
    END IF

  END SUBROUTINE MolecCav_Reduced_density_psi_1p1D_R2_real


  SUBROUTINE MolecCav_Reduced_density_psi_1p1D_R1_real(weight_dim_1, weight_dim_2, Psi_R1, Debug)
    USE QDUtil_m
    USE Mapping_m
    IMPLICIT NONE

    real(kind=Rkind),  intent(inout) :: weight_dim_1(:)                                                                          ! already allocated !
    real(kind=Rkind),  intent(inout) :: weight_dim_2(:)                                                                          ! already allocated !
    real(kind=Rkind),  intent(in)    :: Psi_R1(:)                                                                                ! already allocated !
    logical, optional, intent(in)    :: Debug

    logical                          :: Debug_local = .FALSE.
    integer                          :: i
    real(kind=Rkind), allocatable    :: Psi_R2(:,:)

    IF (PRESENT(Debug)) Debug_local = Debug
  
    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "------------------------------------Computing the wavefunction Psi_R1's reduced density weights---------&
                        &--------------------------"
!        WRITE(out_unit,*) "The squared Psi_R1 matrix to be summed :"
!        CALL Write_Vec(Psi_R1**2, out_unit, Size(Psi_R1), info="Psi_R1**2")
      WRITE(out_unit,*)
    END IF
  
    ALLOCATE(Psi_R2(Size(weight_dim_1, dim=1), Size(weight_dim_2, dim=1)))
    CALL Mapping_WF_R1TOR2(Psi_R2, Psi_R1, Debug=.TRUE.)

    CALL Reduced_density_psi(weight_dim_1, weight_dim_2, Psi_R2, Debug=Debug_local)

    IF (Debug_local) THEN
!        WRITE(out_unit,*)
!        CALL Write_Vec(weight_dim_1, out_unit, Size(weight_dim_1, dim=1), info="weight_dim_1")
!        CALL Write_Vec(weight_dim_2, out_unit, Size(weight_dim_2, dim=1), info="weight_dim_2")
      WRITE(out_unit,*) "----------------------------------End computing the wavefunction Psi_R1's reduced density weights-------&
                        &--------------------------"
    END IF
  
  END SUBROUTINE MolecCav_Reduced_density_psi_1p1D_R1_real
  
  SUBROUTINE MolecCav_Reduced_density_psi_1p1D_R1_real_new(weight_dim_1, weight_dim_2, Psi_R1, Ranks_sizes, Debug)
    USE QDUtil_m
    USE ND_indexes_m
    IMPLICIT NONE

    real(kind=Rkind),  intent(inout) :: weight_dim_1(:)                                                                          ! already allocated !
    real(kind=Rkind),  intent(inout) :: weight_dim_2(:)                                                                          ! already allocated !
    real(kind=Rkind),  intent(in)    :: Psi_R1(:)                                                                                ! already allocated !
    integer,           intent(in)    :: Ranks_sizes(2)                                                                               ! the basis set size of each dimension of the system. Already allocated !
    logical, optional, intent(in)    :: Debug

    logical                          :: Debug_local = .FALSE.
    integer                          :: i
    real(kind=Rkind), allocatable    :: Psi_R2(:,:)

    integer                          :: Ranks_sizes_local(3)
    TYPE(ND_indexes_t)               :: ND_indexes
    integer                          :: List_indexes(3)                                                                          ! actually only two indexes are needed (2D here) but we left the "Cube" framework to prepare for the MpKD case
    logical                          :: Continue_loop

    IF (PRESENT(Debug)) Debug_local = Debug
  
    IF (Size(weight_dim_1, dim=1) /= Ranks_sizes(1)) THEN
      WRITE(out_unit,*) "The dimension of the weight vector for the first system's dimension and the one declared in Ranks_sizes &
                        &do not match. Please check initialization."
      WRITE(out_unit,*) "Size(weight_dim_1) = "//TO_string(Size(weight_dim_1))//"; Ranks_sizes(1) = "//TO_string(Ranks_sizes(1))
      STOP "### The dimension of the weight vector 1 and the one declared in Ranks_sizes do not match."
    END IF

    IF (Size(weight_dim_2, dim=1) /= Ranks_sizes(2)) THEN
      WRITE(out_unit,*) "The dimension of the weight vector for the second system's dimension and the one declared in Ranks_sizes&
                        & do not match. Please check initialization."
      WRITE(out_unit,*) "Size(weight_dim_2) = "//TO_string(Size(weight_dim_2))//"; Ranks_sizes(2) = "//TO_string(Ranks_sizes(2))
      STOP "### The dimension of the weight vector 2 and the one declared in Ranks_sizes do not match."
    END IF

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "------------------------------------Computing the wavefunction Psi_R1's reduced density weights---------&
                        &--------------------------"
      WRITE(out_unit,*) "The squared Psi_R1 matrix to be summed :"
      CALL Write_Vec(Psi_R1**2, out_unit, Size(Psi_R1), info="Psi_R1**2")
      WRITE(out_unit,*)
    END IF
  
    Ranks_sizes_local = [1, Ranks_sizes(:)] ! maybe strange writing be careful
    CALL Initialize_ND_indexes(ND_indexes, Ranks_sizes, Starting_indexes=[1,1,1], Debug=Debug_local)
    List_indexes = Initialize_List_indexes(ND_indexes)
    IF (Debug_local) CALL Write_Vec(List_indexes, out_unit, 3, info="List_indexes")


    IF (Debug_local) THEN
!        WRITE(out_unit,*)
!        CALL Write_Vec(weight_dim_1, out_unit, Size(weight_dim_1, dim=1), info="weight_dim_1")
!        CALL Write_Vec(weight_dim_2, out_unit, Size(weight_dim_2, dim=1), info="weight_dim_2")
      WRITE(out_unit,*) "----------------------------------End computing the wavefunction Psi_R1's reduced density weights-------&
                        &--------------------------"
    END IF

  END SUBROUTINE MolecCav_Reduced_density_psi_1p1D_R1_real_new


END MODULE
