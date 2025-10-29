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
! to be written soon. Anyway, this module is (for now) not used in the MolecCav so don't break your
! head too much trying to understand it rather than another one. Stay healthy. And keep hydrated.
! This module is in the "-6" level.
!
!==================================================================================================
!==================================================================================================
MODULE ND_indexes_m
  !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
  USE QDUtil_m
  IMPLICIT NONE


  TYPE :: ND_indexes_t
    integer              :: N_dim               = 0                                                                              ! number of dimensions of the tensor
    integer, allocatable :: Starting_indexes(:)                                                                                  ! the labelling of the basis vectors may starts from 0 or 1
    integer              :: Begin_right         = 0 ! 0=FALSE; 1=TRUE                                                            ! one can either iterate/increment from the first index (leftmost) => Begin_right=FALSE=0 or from the last one (rightmost) => Befin_right=TRUE=1
    integer, allocatable :: Ranks_sizes(:)                                                                                       ! the basis set size for each dimension of the tensor /!\ NOT THE INDEX OF THE LAST BASIS VECTOR !!! /!\
    integer              :: NB                  = 0                                                                              ! dimension of the tensor producted vector = product(Ranks_sizes)
  END TYPE
  

  PRIVATE

  PUBLIC ND_indexes_t, Initialize, Initialize_List_indexes, Increment, Deallocate, Write

  INTERFACE Initialize
    MODULE PROCEDURE MolecCav_Initialize_ND_indexes
  END INTERFACE
  INTERFACE Initialize_List_indexes
    MODULE PROCEDURE MolecCav_Initialize_List_indexes
  END INTERFACE
  INTERFACE Increment
    MODULE PROCEDURE MolecCav_Increment_indexes
  END INTERFACE
  INTERFACE Deallocate
    MODULE PROCEDURE MolecCav_Deallocate_ND_indexes
  END INTERFACE
  INTERFACE Write
    MODULE PROCEDURE MolecCav_Write_ND_indexes
  END INTERFACE
    

  CONTAINS

  
  SUBROUTINE MolecCav_Initialize_ND_indexes(ND_indexes, Ranks_sizes, Starting_indexes, Begin_right, Debug)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
    USE QDUtil_m
    IMPLICIT NONE 

    TYPE(ND_indexes_t), intent(inout) :: ND_indexes
    integer,            intent(in)    :: Ranks_sizes(:)                                                                          ! the basis set size for each dimension of the tensor /!\ i.e. the number of elements in each dimension AND NOT the upper bound of the indexes, which also depends on the starting indexes/labelling ! /!\
    integer, optional,  intent(in)    :: Starting_indexes(:)                                                                     ! the labelling of the basis vectors may starts from 0 or 1 each
    integer, optional,  intent(in)    :: Begin_right
    logical, optional,  intent(in)    :: Debug

    logical                           :: Debug_local
    integer                           :: NB
    
    !--- Debugging options --------------------------------
    IF (PRESENT(Debug)) THEN; Debug_local = Debug
    ELSE; Debug_local = .FALSE.; END IF

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "o Arguments of MolecCav_Initialize_ND_indexes :"
      WRITE(out_unit,*) "The <<ND_indexes>> argument       :" 
      CALL Write(ND_indexes)
      WRITE(out_unit,*) "The <<Ranks_sizes>> argument      :" 
      CALL Write_Vec(Ranks_sizes, out_unit, SIZE(Ranks_sizes, dim=1), info="Ranks_sizes")
      WRITE(out_unit,*) "The <<Starting_indexes>> argument :" 
      CALL Write_Vec(Starting_indexes, out_unit, SIZE(Starting_indexes, dim=1), info="Starting_indexes")
      WRITE(out_unit,*) "The <<Begin_right>> argument      : "//TO_string(Begin_right)
      FLUSH(out_unit)
    END IF

    !--- Checking dimensions and Allocating ND_indexes%Starting_Indexes
    IF (PRESENT(Starting_indexes) .AND. Size(Starting_indexes) /= Size(Ranks_sizes)) THEN
      WRITE(out_unit,*) "### The sizes of the optional argument Starting_index and of the Ranks_sizes argument do not match."
      WRITE(out_unit,*) "    SIZE(Starting_indexes) = "//TO_string(SIZE(Starting_indexes))//"; SIZE(Ranks_sizes) = "//TO_string(S&
      &IZE(Ranks_sizes))
      WRITE(out_unit,*) "    Please check the arguments."
      STOP "### The sizes of the Starting_index and the Ranks_sizes arguments do not match. cf. output file for more information."
    ELSE IF (PRESENT(Starting_indexes) .AND. Size(Starting_indexes) == Size(Ranks_sizes)) THEN
      ALLOCATE(ND_indexes%Starting_indexes(Size(Ranks_sizes)))
      ND_indexes%Starting_indexes = Starting_indexes
    ELSE 
      ALLOCATE(ND_indexes%Starting_indexes(Size(Ranks_sizes)))      
      ND_indexes%Starting_indexes = 1                                                                                            ! as a default arbitrary choice, the labelling of the basis vectors starts as the Fortran convention i.e. from 1
    END IF

    !--- Finishing initialisation -------------------------
    IF (PRESENT(Begin_right)) ND_indexes%Begin_right = Begin_right
    ND_indexes%N_dim       = SIZE(Ranks_sizes)
    ALLOCATE(ND_indexes%Ranks_sizes(ND_indexes%N_dim))
    ND_indexes%Ranks_sizes = Ranks_sizes
    ND_indexes%NB          = PRODUCT(ND_indexes%Ranks_sizes)

    !--- Conclusion ---------------------------------------
    IF (Debug_local) THEN
      WRITE(out_unit,*) "--- The initialised ND_indexes :"
      CALL Write(ND_indexes)
    END IF 
          
  END SUBROUTINE MolecCav_Initialize_ND_indexes


  FUNCTION MolecCav_Initialize_List_indexes(ND_indexes) RESULT(List_indexes)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
    USE QDUtil_m
    IMPLICIT NONE 

    integer, allocatable           :: List_indexes(:)                                                                            ! the current values of the indexes for each dimension
    TYPE(ND_indexes_t), intent(in) :: ND_indexes

    integer                        :: Zero_index

    List_indexes = ND_indexes%Starting_indexes                                                                                   ! dynamic allocation : allows to call the function for any allocatable object, already allocated or not, and even for a tabular not declared allocatable it seems

    Zero_index = 1 + ND_indexes%Begin_right*(ND_indexes%N_dim-1)
    List_indexes(Zero_index) = List_indexes(Zero_index) - 1                                                                      ! to prepare the initial point I = 0 (state before first loop)

  END FUNCTION MolecCav_Initialize_List_indexes


  SUBROUTINE MolecCav_Increment_indexes(Continue_loop, List_indexes, ND_indexes, Debug)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
    USE QDUtil_m
    IMPLICIT NONE 

    logical,            intent(inout) :: Continue_loop
    integer,            intent(inout) :: List_indexes(:)
    TYPE(ND_indexes_t), intent(inout) :: ND_indexes
    logical, optional,  intent(in)    :: Debug 

    logical                           :: Debug_local
    integer                           :: N, i

    !--- Debugging options --------------------------------
    IF (PRESENT(Debug)) THEN; Debug_local = Debug
    ELSE; Debug_local = .FALSE.; END IF

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "o Arguments of MolecCav_Increment_indexes :"
      WRITE(out_unit,*) "The <<Continue_loop>> argument :"//TO_string(Continue_loop)
      WRITE(out_unit,*) "The <<List_indexes>> argument  :" 
      CALL Write_Vec(List_indexes, out_unit, SIZE(List_indexes, dim=1), info="List_indexes")
      WRITE(out_unit,*) "The <<ND_indexes>> argument    :" 
      CALL Write(ND_indexes)
      FLUSH(out_unit)
    END IF

    !--- Initialization -----------------------------------
    N = ND_indexes%N_dim
    Continue_loop = .TRUE.
    i = 1 + ND_indexes%Begin_right*(N-1)                                                                                         ! we chose first i from 0 to N using the index N-i in the loop and the final test i==N (so looping on dim=N then N-1 ... up to 1), but we are forced to loop in the other sens since the construct TotH has to loop on i_M first then i_C, and so does mappping 2DTO1D

    !--- Incrementation -----------------------------------
    DO WHILE (Continue_loop)
      List_indexes(i) = List_indexes(i) + 1
      IF (Debug_local) CALL Write_Vec(List_indexes, out_unit, N, info="-> incremented List_indexes")

      IF (List_indexes(i) > (ND_indexes%Ranks_sizes(i) + ND_indexes%Starting_indexes(i) - 1)) THEN                               ! i<N <=> N-i>0 (when overcome the last dimension, i is still incremented => N-i goes to 0 at the last iteration)
        List_indexes(i) = ND_indexes%Starting_indexes(i)                                                                         ! upperbound of the rank index i = Ranks_sizes(i)(-0)   if the rank index started at  1(=Starting_indexes(i)) (does count the number of elements);
        i = i + 1 - 2*ND_indexes%Begin_right                                                                                     !                                = Ranks_sizes(i) -1    if       ...      started at  0(=Starting_indexes(i)) (there is an offset of 1 towards small values);
        IF (ALL(List_indexes == ND_indexes%Starting_indexes) .AND. i==(N+1 - (N+1)*ND_indexes%Begin_right)) THEN                 !                                = Ranks_sizes(i) -2    if       ...      started at -1(=Starting_indexes(i))
          IF (Debug_local) CALL Write_Vec(List_indexes, out_unit, N, info="-> incremented List_indexes")                         !                                = Ranks_sizes(i)+[k-1] if       ...      started at  k(=Starting_indexes(i)), k\in\mathbb{N}
          List_indexes = ND_indexes%Ranks_sizes + ND_indexes%Starting_indexes - 1
          Continue_loop = .FALSE.
          IF (Debug_local) CALL Write_Vec(List_indexes, out_unit, N, info="-> final List_indexes")
        END IF
      ELSE 
        EXIT
      END IF
    END DO

    !--- Conclusion ---------------------------------------
    IF (Debug_local) THEN
      WRITE(out_unit,*) "The final Continue_loop      : "//TO_string(Continue_loop)
      WRITE(out_unit,*) "The incremented List_indexes :" 
      CALL Write_Vec(List_indexes, out_unit, SIZE(List_indexes, dim=1), info="List_indexes")
      WRITE(out_unit,*) "The final ND_indexes         :" 
      CALL Write(ND_indexes)
    END IF 

  END SUBROUTINE MolecCav_Increment_indexes


  SUBROUTINE MolecCav_Write_ND_indexes(ND_indexes, Info)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
    USE QDUtil_m
    IMPLICIT NONE

    TYPE(ND_indexes_t), intent(in)         :: ND_indexes
    character(len=*), optional, intent(in) :: Info

    IF (PRESENT(Info)) THEN
      WRITE(out_unit,*) "--- Parameters associated to the object of derived type ND_indexes ("//Info//") :"
    ELSE
      WRITE(out_unit,*) "--- Parameters associated to the object of derived type ND_indexes :"
    END IF

    WRITE(out_unit,*) "N_dim = "//TO_string(ND_indexes%N_dim)

    IF (ALLOCATED(ND_indexes%Starting_indexes)) THEN
      WRITE(out_unit,*) "Starting_indexes IS allocated, and SIZE(Starting_indexes, dim=1) = "//TO_string(SIZE(ND_indexes%Starting&
      &_indexes))
      CALL Write_Vec(ND_indexes%Starting_indexes, out_unit, Size(ND_indexes%Starting_indexes), info="Starting_indexes")
    ELSE 
      WRITE(out_unit,*) "Starting_indexes is NOT allocated"
    END IF

    WRITE(out_unit,*) "Begin_right = "//TO_string(ND_indexes%Begin_right)

    IF (ALLOCATED(ND_indexes%Ranks_sizes)) THEN
      WRITE(out_unit,*) "Ranks_sizes IS allocated, and SIZE(Ranks_sizes, dim=1) = "//TO_string(SIZE(ND_indexes%Ranks_sizes))
      CALL Write_Vec(ND_indexes%Ranks_sizes, out_unit, Size(ND_indexes%Ranks_sizes), info="Ranks_sizes")
    ELSE 
      WRITE(out_unit,*) "Ranks_sizes is NOT allocated"
    END IF

    WRITE(out_unit,*) "NB = "//TO_string(ND_indexes%NB)
    FLUSH(out_unit)
  
  END SUBROUTINE MolecCav_Write_ND_indexes


  SUBROUTINE MolecCav_Deallocate_ND_indexes(ND_indexes, Debug)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
    USE QDUtil_m
    IMPLICIT NONE 

    TYPE(ND_indexes_t), intent(inout) :: ND_indexes
    logical, optional,     intent(in) :: Debug                                                       ! cf. comments in HO1D_parameters_m

    logical                           :: Debug_local

    !--- Debugging options --------------------------------
    IF (PRESENT(Debug)) THEN; Debug_local = Debug
    ELSE; Debug_local = .FALSE.; END IF

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "o Arguments of MolecCav_Deallocate_ND_indexes :"
      WRITE(out_unit,*) "The <<ND_indexes>> argument    :" 
      CALL Write(ND_indexes)
      FLUSH(out_unit)
    END IF

    !--- Deallocation -------------------------------------
    ND_indexes%N_dim       = 0
    IF (ALLOCATED(ND_indexes%Ranks_sizes))      DEALLOCATE(ND_indexes%Ranks_sizes)
    ND_indexes%Begin_right = 0
    IF (ALLOCATED(ND_indexes%Starting_indexes)) DEALLOCATE(ND_indexes%Starting_indexes)
    ND_indexes%NB          = 0

    !--- Conclusion ---------------------------------------
    IF (Debug_local) THEN
      CALL Write(ND_indexes)
    END IF 

  END SUBROUTINE MolecCav_Deallocate_ND_indexes


END MODULE
