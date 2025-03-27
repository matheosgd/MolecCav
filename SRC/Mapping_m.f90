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
MODULE Mapping_m
  !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
  USE QDUtil_m
  USE ND_indexes_m
  IMPLICIT NONE
    

  PRIVATE

  PUBLIC Mapping_WF_R2TOR1, Mapping_WF_R1TOR2

  INTERFACE Mapping_WF_R2TOR1
    MODULE PROCEDURE MolecCav_Mapping_WF_R2TOR1
  END INTERFACE
  INTERFACE Mapping_WF_R1TOR2
    MODULE PROCEDURE MolecCav_Mapping_WF_R1TOR2
  END INTERFACE

  
  CONTAINS


  SUBROUTINE MolecCav_Mapping_WF_R2TOR1(Psi_R1, Psi_R2, Debug)
    USE QDUtil_m
    IMPLICIT NONE

    real(kind=Rkind),  intent(inout) :: Psi_R1(:)
    real(kind=Rkind),  intent(in)    :: Psi_R2(:,:)
    logical, optional, intent(in)    :: Debug

    logical                          :: Debug_local = .FALSE.
    integer                          :: Nb_1, Nb_2, i_1, i_2, NB, I

    IF (PRESENT(Debug)) Debug_local = Debug

    Nb_1   = Size(Psi_R2, dim=1)
    Nb_2   = Size(Psi_R2, dim=2)
    NB     = Size(Psi_R1)
    I      = 0
    Psi_R1 = ZERO

    IF (Nb_1*Nb_2 /= NB) THEN
      STOP "### Wrong size of the matrices"
    END IF

    DO i_2 = 1, Nb_2
      DO i_1 = 1, Nb_1
        I = I + 1
        Psi_R1(I) = Psi_R2(i_1, i_2)
      END DO
    END DO

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "The non-tensor producted WF Psi_R2"
      CALL Write_Mat(Psi_R2, out_unit, Nb_2, info="Psi_R2")
      WRITE(out_unit,*) "                     ||"
      WRITE(out_unit,*) "                     \/"
      WRITE(out_unit,*) "The tensor producted (mapped) WF Psi_R1"
      CALL Write_Vec(Psi_R1, out_unit, 1, info="Psi_R1")
    END IF 
 
  END SUBROUTINE MolecCav_Mapping_WF_R2TOR1
  
  
  SUBROUTINE MolecCav_Mapping_WF_R1TOR2(Psi_R2, Psi_R1, Starting_indexes, Debug)
    USE QDUtil_m
    IMPLICIT NONE

    real(kind=Rkind),   intent(inout) :: Psi_R2(:,:)
    real(kind=Rkind),   intent(in)    :: Psi_R1(:)
    integer, optional,  intent(in)    :: Starting_indexes(:)
    logical, optional,  intent(in)    :: Debug

    integer, allocatable              :: Ranks_sizes(:), List_indexes(:)
    TYPE(ND_indexes_t)                :: ND_indexes
    integer                           :: NB, I
    logical                           :: Debug_local = .FALSE., Continue_loop

    IF (PRESENT(Debug)) Debug_local = Debug

    ALLOCATE(Ranks_sizes(RANK(Psi_R2)))
    Ranks_sizes = [(Size(Psi_R2, dim=i), i=1,RANK(Psi_R2))]

    IF (Size(Psi_R1) /= PRODUCT(Ranks_sizes)) THEN
      WRITE(out_unit,*) "Wrong size of the matrices : Size(Psi_R2, dim=1)*...*Size(Psi_R2, dim=rank(psi_2d)) = "//TO_string(Size(&
                        &Psi_R2, dim=1))//"*"//TO_string(Size(Psi_R2, dim=2))//" /= Size(Psi_R1) = "//TO_string(Size(Psi_R1))
      STOP "### Wrong size of the matrices. cf. output file for more information."
    END IF

    IF (PRESENT(Starting_indexes)) THEN
      CALL Initialize_ND_indexes(ND_indexes, Ranks_sizes, Starting_indexes=Starting_indexes, Debug=Debug_local)
    ELSE 
      CALL Initialize_ND_indexes(ND_indexes, Ranks_sizes, Debug=Debug_local)
    END IF
    List_indexes = Initialize_List_indexes(ND_indexes)
    I      = 0
    Psi_R2 = ZERO

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "Ranks_sizes has been initialized as : "
      CALL Write_Vec(Ranks_sizes, out_unit, Size(Ranks_sizes), info="Ranks_sizes")
      WRITE(out_unit,*) "NB recovered as : "//TO_string(ND_indexes%NB)//"; and the ND_indexes object has been initialized as : "
      CALL Write_ND_indexes(ND_indexes)
    END IF

    IF (Debug_local) CALL Write_Vec(List_indexes, out_unit, ND_indexes%N_dim, info="List_indexes")
    
    Continue_loop = .TRUE.
    DO WHILE (Continue_loop)
      IF (Debug_local) THEN
        WRITE(out_unit,*)
        WRITE(out_unit,*) "After looping "//TO_string(I)//" : I = "//TO_string(I)//"; NB = "//TO_string(ND_indexes%NB)
        CALL Write_Vec(List_indexes,           out_unit, Size(List_indexes),           info="List_indexes")
        CALL Write_Vec(ND_indexes%Ranks_sizes, out_unit, Size(ND_indexes%Ranks_sizes), info="Ranks_sizes")
      END IF

      CALL Increment_indexes(Continue_loop, List_indexes, ND_indexes, Debug=Debug_local)
      IF (Debug_local) WRITE(out_unit,*) "--> Continue_loop : "//TO_string(Continue_loop)

      IF (Continue_loop) THEN
        I = I + 1
        IF (Debug_local) WRITE(out_unit,*)
        IF (Debug_local) WRITE(out_unit,*) "Looping "//TO_string(I)//" -> I = "//TO_string(I)//" :"
      END IF 

      IF (I>ND_indexes%NB) THEN
        WRITE(out_unit,*) "I should not be greater that NB !"
        STOP              "I should not be greater that NB !"
      END IF
      Psi_R2(List_indexes(1), List_indexes(2)) = Psi_R1(I)
      IF (Debug_local) CALL Write_Mat(Psi_R2, out_unit, Ranks_sizes(2), info="Psi_R2")
    END DO

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "The tensor producted (mapped) WF Psi_R1"
      CALL Write_Vec(Psi_R1, out_unit, 1, info="Psi_R1")
      WRITE(out_unit,*) "                     ||"
      WRITE(out_unit,*) "                     \/"
      WRITE(out_unit,*) "The non-tensor producted WF Psi_R2"
      CALL Write_Mat(Psi_R2, out_unit, Ranks_sizes(2), info="Psi_R2")
    END IF 

    CALL Deallocate_ND_indexes(ND_indexes)

  END SUBROUTINE MolecCav_Mapping_WF_R1TOR2
 

END MODULE
  