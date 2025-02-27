!==================================================================================================
!==================================================================================================
!This file is part of MolecCav.
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
MODULE Mapping_m
  !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
  USE QDUtil_m
  USE ND_indexes_m
  IMPLICIT NONE
    

  PRIVATE

  PUBLIC Mapping_WF_2DTO1D, Mapping_WF_1DTO2D

  INTERFACE Mapping_WF_2DTO1D
    MODULE PROCEDURE MolecCav_Mapping_WF_2DTO1D
  END INTERFACE
  INTERFACE Mapping_WF_1DTO2D
    MODULE PROCEDURE MolecCav_Mapping_WF_1DTO2D
  END INTERFACE

  
  CONTAINS


  SUBROUTINE MolecCav_Mapping_WF_2DTO1D(Psi_1D, Psi_2D, Debug)
    USE QDUtil_m
    IMPLICIT NONE

    real(kind=Rkind),  intent(inout) :: Psi_1D(:)
    real(kind=Rkind),  intent(in)    :: Psi_2D(:,:)
    logical, optional, intent(in)    :: Debug

    logical                          :: Debug_local = .FALSE.
    integer                          :: Nb_1, Nb_2, i_1, i_2, NB, I

    IF (PRESENT(Debug)) Debug_local = Debug

    Nb_1   = Size(Psi_2D, dim=1)
    Nb_2   = Size(Psi_2D, dim=2)
    NB     = Size(Psi_1D)
    I      = 0
    Psi_1D = ZERO

    IF (Nb_1*Nb_2 /= NB) THEN
      STOP "### Wrong size of the matrices"
    END IF

    DO i_2 = 1, Nb_2
      DO i_1 = 1, Nb_1
        I = I + 1
        Psi_1D(I) = Psi_2D(i_1, i_2)
      END DO
    END DO

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "The non-tensor producted WF Psi_2D"
      CALL Write_Mat(Psi_2D, out_unit, Nb_2, info="Psi_2D")
      WRITE(out_unit,*) "                     ||"
      WRITE(out_unit,*) "                     \/"
      WRITE(out_unit,*) "The tensor producted (mapped) WF Psi_1D"
      CALL Write_Vec(Psi_1D, out_unit, 1, info="Psi_1D")
    END IF 
 
  END SUBROUTINE MolecCav_Mapping_WF_2DTO1D
  
  
  SUBROUTINE MolecCav_Mapping_WF_1DTO2D(Psi_2D, Psi_1D, Starting_indexes, Debug)
    USE QDUtil_m
    IMPLICIT NONE

    real(kind=Rkind),   intent(inout) :: Psi_2D(:,:)
    real(kind=Rkind),   intent(in)    :: Psi_1D(:)
    integer, optional,  intent(in)    :: Starting_indexes(:)
    logical, optional,  intent(in)    :: Debug

    integer, allocatable              :: Ranks_sizes(:), List_indexes(:)
    TYPE(ND_indexes_t)                :: ND_indexes
    integer                           :: NB, I
    logical                           :: Debug_local = .FALSE., Continue_loop

    IF (PRESENT(Debug)) Debug_local = Debug

    ALLOCATE(Ranks_sizes(RANK(Psi_2D)))
    Ranks_sizes = [(Size(Psi_2D, dim=i), i=1,RANK(Psi_2D))]

    IF (Size(Psi_1D) /= PRODUCT(Ranks_sizes)) THEN
      WRITE(out_unit,*) "Wrong size of the matrices : Size(Psi_2D, dim=1)*...*Size(Psi_2D, dim=rank(psi_2d)) = "//TO_string(Size(&
                        &Psi_2D, dim=1))//"*"//TO_string(Size(Psi_2D, dim=2))//" /= Size(Psi_1D) = "//TO_string(Size(Psi_1D))
      STOP "### Wrong size of the matrices. cf. output file for more information."
    END IF

    IF (PRESENT(Starting_indexes)) THEN
      CALL Initialize_ND_indexes(ND_indexes, Ranks_sizes, Starting_indexes=Starting_indexes, Debug=Debug)
    ELSE 
      CALL Initialize_ND_indexes(ND_indexes, Ranks_sizes, Debug=Debug)
    END IF
    CALL Initialize_List_indexes(List_indexes, ND_indexes)
    I      = 0
    Psi_2D = ZERO

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "Ranks_sizes has been initialized as : "
      CALL Write_Vec(Ranks_sizes, out_unit, Size(Ranks_sizes), info="Ranks_sizes")
      WRITE(out_unit,*) "NB recovered as : "//TO_string(ND_indexes%NB)//"; and the ND_indexes object has been initialized as : "
      CALL Write_ND_indexes(ND_indexes)
    END IF

!    List_indexes = List_indexes + 1 - ND_indexes%Starting_indexes ! List_indexes may start from 0 but not Psi_2D : if start from 0 -> +1, if star from 1 -> +0
    IF (Debug_local) CALL Write_Vec(List_indexes, out_unit, ND_indexes%N_dim, info="List_indexes")
    
    Continue_loop = .TRUE.
    DO WHILE (Continue_loop)
      IF (Debug_local) THEN
        WRITE(out_unit,*)
        WRITE(out_unit,*) "After looping "//TO_string(I)//" : I = "//TO_string(I)//"; NB = "//TO_string(ND_indexes%NB)
        CALL Write_Vec(List_indexes,           out_unit, Size(List_indexes),           info="List_indexes")
        CALL Write_Vec(ND_indexes%Ranks_sizes, out_unit, Size(ND_indexes%Ranks_sizes), info="Ranks_sizes")
      END IF

      I = I + 1; IF (Debug_local) WRITE(out_unit,*) "Looping "//TO_string(I)//" -> I = "//TO_string(I)
      IF (I>ND_indexes%NB) THEN
        WRITE(out_unit,*) "I should not be greater that NB !"
        STOP              "I should not be greater that NB !"
      END IF

      Psi_2D(List_indexes(1), List_indexes(2)) = Psi_1D(I)
      IF (Debug_local) CALL Write_Mat(Psi_2D, out_unit, Ranks_sizes(2), info="Psi_2D")

      CALL Increment_indexes(Continue_loop, List_indexes, ND_indexes, Debug=Debug_local)
      IF (Debug_local) WRITE(out_unit,*) "Continue_loop : "//TO_string(Continue_loop)
    END DO

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "The tensor producted (mapped) WF Psi_1D"
      CALL Write_Vec(Psi_1D, out_unit, 1, info="Psi_1D")
      WRITE(out_unit,*) "                     ||"
      WRITE(out_unit,*) "                     \/"
      WRITE(out_unit,*) "The non-tensor producted WF Psi_2D"
      CALL Write_Mat(Psi_2D, out_unit, Ranks_sizes(2), info="Psi_2D")
    END IF 

    CALL Deallocate_ND_indexes(ND_indexes)

  END SUBROUTINE MolecCav_Mapping_WF_1DTO2D
 

END MODULE
  