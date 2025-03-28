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
    USE ND_indexes_m
    IMPLICIT NONE

    real(kind=Rkind),  intent(inout) :: weight_dim_1(:)                                                                          ! already allocated !
    real(kind=Rkind),  intent(inout) :: weight_dim_2(:)                                                                          ! already allocated !
    real(kind=Rkind),  intent(in)    :: Psi_R1(:)                                                                                ! already allocated !
!    integer,           intent(in)    :: Ranks_sizes(2)                                                                               ! the basis set size of each dimension of the system. Already allocated !
    logical, optional, intent(in)    :: Debug

    integer                          :: Ranks_sizes_local(3)
    TYPE(ND_indexes_t)               :: ND_indexes
    integer                          :: List_indexes(3)                                                                          ! actually only two indexes are needed (2D here) but we left the "Cube" framework to prepare for the MpKD case
    logical                          :: Continue_loop
    real(kind=Rkind), allocatable    :: Cube(:,:,:)
    integer                          :: I = 0, i_1, i_2, i_3
    logical                          :: Debug_local = .FALSE.

    IF (PRESENT(Debug)) Debug_local = Debug
  
!ASK BEFORE DELETE
!    IF (Size(weight_dim_1, dim=1) /= Ranks_sizes(1)) THEN
!      WRITE(out_unit,*) "The dimension of the weight vector for the first system's dimension and the one declared in Ranks_sizes &
!                        &do not match. Please check initialization."
!      WRITE(out_unit,*) "Size(weight_dim_1) = "//TO_string(Size(weight_dim_1))//"; Ranks_sizes(1) = "//TO_string(Ranks_sizes(1))
!      STOP "### The dimension of the weight vector 1 and the one declared in Ranks_sizes do not match."
!    END IF

!    IF (Size(weight_dim_2, dim=1) /= Ranks_sizes(2)) THEN
!      WRITE(out_unit,*) "The dimension of the weight vector for the second system's dimension and the one declared in Ranks_sizes&
!                        & do not match. Please check initialization."
!      WRITE(out_unit,*) "Size(weight_dim_2) = "//TO_string(Size(weight_dim_2))//"; Ranks_sizes(2) = "//TO_string(Ranks_sizes(2))
!      STOP "### The dimension of the weight vector 2 and the one declared in Ranks_sizes do not match."
!    END IF

    Ranks_sizes_local = [1, Size(weight_dim_1), Size(weight_dim_2)]
    IF (Size(Psi_R1, dim=1) /= PRODUCT(Ranks_sizes_local)) THEN
      WRITE(out_unit,*) "### The dimension of the wavevector to analyse Psi_R1 does not match the declared dimensions of the subs&
                        &ystems in weights_dim_1,2. Please check initialization."
      WRITE(out_unit,*) "   Size(Psi_R1, dim=1) = "//TO_string(Size(Psi_R1, dim=1))//"; PRODUCT(Ranks_sizes) = "//TO_string(PRODU&
                        &CT(Ranks_sizes_local))
      CALL Write_Vec(Ranks_sizes_local, out_unit, 3, info="   [1, Size(weight_dim_1), Size(weight_dim_2)] = ")
      STOP "### The dimension of Psi_R1 does not match the values of Ranks_sizes_local."
    END IF

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "------------------------------------Computing the wavefunction Psi_R1's reduced density weights---------&
                        &--------------------------"
      WRITE(out_unit,*) "The squared Psi_R1 matrix to be summed /!\ Rank 1, just written in line /!\ :"
      CALL Write_Vec(Psi_R1**2, out_unit, Size(Psi_R1), info="Psi_R1**2")
      WRITE(out_unit,*)
    END IF
  
    CALL Initialize_ND_indexes(ND_indexes, Ranks_sizes_local, Starting_indexes=[1,1,1], Debug=Debug_local)
    List_indexes = Initialize_List_indexes(ND_indexes)
    IF (Debug_local) CALL Write_Vec(List_indexes, out_unit, 3, info="List_indexes")

    weight_dim_1 = ZERO
    weight_dim_2 = ZERO
!--------------------------------------------------------------------------------------------------------------------------------!
! A first straightforward and easy-reading version of the procedure, improved afterwards.                                        !
!--------------------------------------------------------------------------------------------------------------------------------!
!    ALLOCATE(Cube(ND_indexes%Ranks_sizes(1), ND_indexes%Ranks_sizes(2), ND_indexes%Ranks_sizes(3)))
!
!    DO
!      IF (Debug_local) THEN
!        WRITE(out_unit,*)
!        WRITE(out_unit,*) "After looping "//TO_string(I)//" : I = "//TO_string(I)//"; NB = "//TO_string(ND_indexes%NB)
!        CALL Write_Vec(List_indexes,           out_unit, Size(List_indexes),           info="List_indexes_"//TO_string(I)//" =")
!        CALL Write_Vec(ND_indexes%Ranks_sizes, out_unit, Size(ND_indexes%Ranks_sizes), info="Ranks_sizes    =")
!      END IF
!
!      CALL Increment_indexes(Continue_loop, List_indexes, ND_indexes, Debug=.FALSE.)
!      I = I + 1
!      IF (Debug_local) WRITE(out_unit,*) "--> Continue_loop : "//TO_string(Continue_loop)
!      IF (.NOT. Continue_loop) EXIT 
!      IF (I>ND_indexes%NB) THEN
!        WRITE(out_unit,*) "The looping should not continue when I is greater that NB !"
!        STOP              "The looping should not continue when I is greater that NB !"
!      END IF
!
!      IF (Debug_local) WRITE(out_unit,*)
!      IF (Debug_local) WRITE(out_unit,*) "Looping "//TO_string(I)//" -> I = "//TO_string(I)//"..."
!      Cube(List_indexes(1), List_indexes(2), List_indexes(3)) = Psi_R1(I)
!    END DO
!
!    IF (Debug_local) CALL Write_Mat(Cube(1,:,:), out_unit, ND_indexes%Ranks_sizes(3), info="Cube(1,:,:)")
!
!    I = 0
!    DO i_2 = 1, ND_indexes%Ranks_sizes(2)
!      DO i_1 = 1, ND_indexes%Ranks_sizes(1)
!        I = I + 1
!        weight_dim_1(I) = SUM(Cube(i_1,i_2,:)**2)
!      END DO 
!    END DO 
!    I = 0
!    DO i_3 = 1, ND_indexes%Ranks_sizes(3)
!      DO i_1 = 1, ND_indexes%Ranks_sizes(1)
!        I = I + 1
!        weight_dim_2(I) = SUM(Cube(i_1,:,i_3)**2)
!      END DO 
!    END DO 
!
!--------------------------------------------------------------------------------------------------------------------------------!
! The following version does exactly the same thing as the above commented one but in a slightly less easy-reading way. It does  !
! not generates the Cube tensor and directly sort and sum the elements of Psi_R1, taking advantage of the correspondances create-!
! d by Increment_indexes between I and List_indexes.                                                                             !
!--------------------------------------------------------------------------------------------------------------------------------!
    DO
      IF (Debug_local) THEN
        WRITE(out_unit,*)
        WRITE(out_unit,*) "After looping "//TO_string(I)//" : I = "//TO_string(I)//"; NB = "//TO_string(ND_indexes%NB)
        CALL Write_Vec(List_indexes,           out_unit, Size(List_indexes),           info="List_indexes_"//TO_string(I)//" =")
        CALL Write_Vec(ND_indexes%Ranks_sizes, out_unit, Size(ND_indexes%Ranks_sizes), info="Ranks_sizes    =")
      END IF

      CALL Increment_indexes(Continue_loop, List_indexes, ND_indexes, Debug=.FALSE.)
      I = I + 1
      IF (Debug_local) WRITE(out_unit,*) "--> Continue_loop : "//TO_string(Continue_loop)
      IF (.NOT. Continue_loop) EXIT 
      IF (I>ND_indexes%NB) THEN
        WRITE(out_unit,*) "The looping should not continue when I is greater that NB !"
        STOP              "The looping should not continue when I is greater that NB !"
      END IF

      IF (Debug_local) WRITE(out_unit,*)
      IF (Debug_local) WRITE(out_unit,*) "Looping "//TO_string(I)//" -> I = "//TO_string(I)//"..."
      weight_dim_1(List_indexes(2)) = weight_dim_1(List_indexes(2)) + Psi_R1(I)**2
      weight_dim_2(List_indexes(3)) = weight_dim_2(List_indexes(3)) + Psi_R1(I)**2
    END DO

    IF (Debug_local) THEN
        WRITE(out_unit,*)
        CALL Write_Vec(weight_dim_1, out_unit, Size(weight_dim_1, dim=1), info="weight_dim_1")
        CALL Write_Vec(weight_dim_2, out_unit, Size(weight_dim_2, dim=1), info="weight_dim_2")
      WRITE(out_unit,*) "----------------------------------End computing the wavefunction Psi_R1's reduced density weights-------&
                        &--------------------------"
    END IF

  END SUBROUTINE MolecCav_Reduced_density_psi_1p1D_R1_real


END MODULE
