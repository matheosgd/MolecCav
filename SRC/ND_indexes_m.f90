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
MODULE ND_indexes_m
  !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
  USE QDUtil_m
  IMPLICIT NONE


  TYPE :: ND_indexes_t
    integer              :: N                   = 0                                                ! number of dimensions of the tensor
    integer, allocatable :: Starting_indexes(:)                                                    ! the labelling of the basis vectors may starts from 0 or 1
    integer, allocatable :: Dim_sizes(:)                                                           ! the basis set size for each dimension of the tensor /!\ NOT THE INDEX OF THE LAST BASIS VECTOR !!! /!\
    integer, allocatable :: Indexes(:)                                                             ! the current values of the indexes for each dimension
  END TYPE
  

  PRIVATE

  PUBLIC ND_indexes_t, Initialize_ND_indexes, Mapping_WF_2DTO1D

  INTERFACE Initialize_ND_indexes
    MODULE PROCEDURE MolecCav_Initialize_ND_indexes
  END INTERFACE
  INTERFACE Mapping_WF_2DTO1D
    MODULE PROCEDURE MolecCav_Mapping_WF_2DTO1D
  END INTERFACE
  INTERFACE Write_ND_indexes
    MODULE PROCEDURE MolecCav_Write_ND_indexes
  END INTERFACE
    

  CONTAINS


  SUBROUTINE MolecCav_Initialize_ND_indexes(ND_indexes, Dim_sizes, Starting_indexes, Debug)
    USE QDUtil_m
    IMPLICIT NONE 

    TYPE(ND_indexes_t), intent(inout) :: ND_indexes
    integer,            intent(in)    :: Dim_sizes(:)                                              ! the basis set size for each dimension of the tensor
    integer, optional,  intent(in)    :: Starting_indexes(:)                                       ! the labelling of the basis vectors may starts from 0 or 1 each
    logical, optional,  intent(in)    :: Debug

    logical                           :: Debug_local = .FALSE.
    
    IF (PRESENT(Debug)) Debug_local = Debug

    IF (Debug_local) THEN
      WRITE(out_unit,*) 
      WRITE(out_unit,*) '********************************************************************************'
      WRITE(out_unit,*) '********************** INITIALIZING ND_INDEXES PARAMETERS **********************'
    END IF 

    IF (PRESENT(Starting_indexes) .AND. Size(Starting_indexes) /= Size(Dim_sizes)) THEN
      WRITE(out_unit,*) "The sizes of the optional argument Starting_index ("//TO_string(Size(Starting_indexes))//") and of the D&
                        &im_sizes argument ("//TO_string(Size(Dim_sizes))//") do not match. Please check the arguments."
      STOP "The sizes of the Starting_index and the Dim_sizes arguments do not match. Please check output file for more informati&
           &on."
    ELSE IF (PRESENT(Starting_indexes) .AND. Size(Starting_indexes) == Size(Dim_sizes)) THEN
      ALLOCATE(ND_indexes%Starting_indexes(Size(Dim_sizes)))
      ND_indexes%Starting_indexes = Starting_indexes
    ELSE 
      ALLOCATE(ND_indexes%Starting_indexes(Size(Dim_sizes)))      
      ND_indexes%Starting_indexes = 1                                                              ! as a default arbitrary choice, the labelling of the basis vectors starts as the Fortran convention i.e. from 1
    END IF

    ND_indexes%N = Size(Dim_sizes)
    ALLOCATE(ND_indexes%Dim_sizes(ND_indexes%N))
    ALLOCATE(ND_indexes%Indexes(  ND_indexes%N))

    ND_indexes%Dim_sizes = Dim_sizes
    ND_indexes%Indexes = Starting_indexes

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "------------ND_indexes constructed by MolecCav_Initialize_ND_indexes------------"
      CALL Write_ND_indexes(ND_indexes)
      WRITE(out_unit,*) "----------End ND_indexes constructed by MolecCav_Initialize_ND_indexes----------"
      WRITE(out_unit,*) 
      WRITE(out_unit,*) '*********************** ND_INDEXES PARAMETERS INITIALIZED **********************'
      WRITE(out_unit,*) '********************************************************************************'      
    END IF
      
  END SUBROUTINE MolecCav_Initialize_ND_indexes


  SUBROUTINE MolecCav_Mapping_WF_2DTO1D(Psi_1D, Psi_2D, Debug)
    USE QDUtil_m
    IMPLICIT NONE

    real(kind=Rkind),  intent(inout) :: Psi_1D(:)
    real(kind=Rkind),  intent(in)    :: Psi_2D(:,:)
    logical, optional, intent(in)    :: Debug

    logical                          :: Debug_local = .FALSE.
    integer                          :: Nb_M, Nb_C, i_M, i_C, NB, I

    IF (PRESENT(Debug)) Debug_local = Debug

    Nb_M   = Size(Psi_2D, dim=1)
    Nb_C   = Size(Psi_2D, dim=2)
    NB     = Size(Psi_1D)
    I      = 0
    Psi_1D = ZERO

    IF (Nb_M*Nb_C /= NB) THEN
      STOP "Wrong size of the matrices"
    END IF

    DO i_C = 1, Nb_C
      DO i_M = 1, Nb_M
        I = I + 1
        Psi_1D(I) = Psi_2D(i_M, i_C)
      END DO
    END DO

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "The non-tensor producted WF Psi_2D"
      CALL Write_Mat(Psi_2D, out_unit, Nb_C, info="Psi_2D")
      WRITE(out_unit,*) "The tensor producted (mapped) WF Psi_1D"
      CALL Write_Vec(Psi_1D, out_unit, 1, info="Psi_1D")
    END IF 

  END SUBROUTINE MolecCav_Mapping_WF_2DTO1D


  SUBROUTINE MolecCav_Mapping_WF_1DTO2D(Psi_2D, Psi_1D, Dim_sizes, Debug)
    USE QDUtil_m
    IMPLICIT NONE

    real(kind=Rkind),  intent(inout) :: Psi_2D(:,:)
    real(kind=Rkind),  intent(in)    :: Psi_1D(:)
    integer,           intent(in)    :: Dim_sizes(:)
    logical, optional, intent(in)    :: Debug

    logical                          :: Debug_local = .FALSE.
    integer, allocatable             :: Indexes(:)
    integer                          :: Nb_M, Nb_C, i_M, i_C, NB, I

    IF (PRESENT(Debug)) Debug_local = Debug

    Nb_M   = Size(Psi_2D, dim=1)
    Nb_C   = Size(Psi_2D, dim=2)
    NB     = Size(Psi_1D)
    I      = 0
    Psi_2D = ZERO

    IF (Nb_M*Nb_C /= NB) THEN
      WRITE(out_unit,*) "Wrong size of the matrices : Size(Psi_2D, dim=1)*Size(Psi_2D, dim=1) = "//&
                        &TO_string(Size(Psi_2D, dim=1))//"*"//TO_string(Size(Psi_2D, dim=1))//" /= &
                        &Size(Psi_1D) = "//TO_string(Size(Psi_1D))
      STOP "Wrong size of the matrices. cf. output file for more information."
    END IF
    IF (Size(Dim_sizes) /= NB) THEN
      WRITE(out_unit,*) "Wrong size of the matrices : Size(Psi_2D, dim=1)*Size(Psi_2D, dim=1) = "//&
      &TO_string(Size(Psi_2D, dim=1))//"*"//TO_string(Size(Psi_2D, dim=1))//" /= &
      &Size(Psi_1D) = "//TO_string(Size(Psi_1D))
      STOP "Wrong size of the matrices. cf. output file for more information."
    END IF

    DO i_C = 1, Nb_C
      DO i_M = 1, Nb_M
        I = I + 1
!        Psi_1D(I) = Psi_2D(i_M, i_C)
      END DO
    END DO

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "The non-tensor producted WF Psi_2D"
      CALL Write_Mat(Psi_2D, out_unit, Nb_C, info="Psi_2D")
      WRITE(out_unit,*) "The tensor producted (mapped) WF Psi_1D"
      CALL Write_Vec(Psi_1D, out_unit, 1, info="Psi_1D")
    END IF 

  END SUBROUTINE MolecCav_Mapping_WF_1DTO2D


  SUBROUTINE MolecCav_Write_ND_indexes(ND_indexes)
    TYPE(ND_indexes_t), intent(in) :: ND_indexes

    WRITE(out_unit,*) "___________________________________The associated ND_indexes object_________________________________"
    WRITE(out_unit,*) "|Number of dimension of the tensor to map (ND_indexes%N)                     | ", ND_indexes%N
    WRITE(out_unit,*) "|____________________________________________________________________________|______________________"
    FLUSH(out_unit)
    IF (ALLOCATED(ND_indexes%Starting_indexes)) THEN
      WRITE(out_unit,*) "|The labelling of the basis vectors information do is allocated, and ...     | "
      WRITE(out_unit,*) "|... this labelling shall start from (ND_indexes%Starting_indexes)           | "
      WRITE(out_unit,*) "|____________________________________________________________________________| "
      CALL Write_Vec(ND_indexes%Starting_indexes, out_unit, Size(ND_indexes%Starting_indexes), info="ND_indexes%Starting_indexes")
      FLUSH(out_unit)
      WRITE(out_unit,*) "|____________________________________________________________________________|"
    ELSE 
      WRITE(out_unit,*) "|The labelling of the basis vectors information is NOT allocated.            | "
      WRITE(out_unit,*) "|____________________________________________________________________________| "
    END IF
    FLUSH(out_unit)
    IF (ALLOCATED(ND_indexes%Dim_sizes)) THEN
      WRITE(out_unit,*) "|The list of the basis set size of each of these dimensions do is allocated..| "
      WRITE(out_unit,*) "|... and is initialized as (ND_indexes%Dim_sizes)                            | "
      WRITE(out_unit,*) "|____________________________________________________________________________| "
      CALL Write_Vec(ND_indexes%Dim_sizes, out_unit, Size(ND_indexes%Dim_sizes), info="ND_indexes%Dim_sizes")
      FLUSH(out_unit)
      WRITE(out_unit,*) "|____________________________________________________________________________|"
    ELSE 
      WRITE(out_unit,*) "|The list of the basis set size of each of these dimensions is NOT allocated.| "
      WRITE(out_unit,*) "|____________________________________________________________________________| "
    END IF
    FLUSH(out_unit)
    IF (ALLOCATED(ND_indexes%Indexes)) THEN
      WRITE(out_unit,*) "|The current values list of the indexes for each dimension do is allocated...| "
      WRITE(out_unit,*) "|... and is initialized as (ND_indexes%Indexes)                              | "
      WRITE(out_unit,*) "|____________________________________________________________________________| "
      CALL Write_Vec(ND_indexes%Indexes, out_unit, Size(ND_indexes%Indexes), info="ND_indexes%Indexes")
      FLUSH(out_unit)
      WRITE(out_unit,*) "|____________________________________________________________________________|"
    ELSE 
      WRITE(out_unit,*) "|The current values list of the indexes for each dimension is NOT allocated. | "
      WRITE(out_unit,*) "|____________________________________________________________________________| "
    END IF
    FLUSH(out_unit)
  

  END SUBROUTINE MolecCav_Write_ND_indexes


END MODULE
