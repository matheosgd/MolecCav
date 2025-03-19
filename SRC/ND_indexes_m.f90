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
MODULE ND_indexes_m
  !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
  USE QDUtil_m
  IMPLICIT NONE


  TYPE :: ND_indexes_t
    integer              :: N_dim               = 0                                                ! number of dimensions of the tensor
    integer, allocatable :: Starting_indexes(:)                                                    ! the labelling of the basis vectors may starts from 0 or 1
    integer              :: Begin_right         = 0 ! 0=FALSE; 1=TRUE                              ! one can either iterate/increment from the first index (leftmost) => Begin_right=FALSE=0 or from the last one (rightmost) => Befin_right=TRUE=1
    integer, allocatable :: Ranks_sizes(:)                                                         ! the basis set size for each dimension of the tensor /!\ NOT THE INDEX OF THE LAST BASIS VECTOR !!! /!\
    integer              :: NB                  = 0                                                ! dimension of the tensor producted vector = product(Ranks_sizes)
  END TYPE
  

  PRIVATE

  PUBLIC ND_indexes_t, Initialize_ND_indexes, Increment_indexes, Deallocate_ND_indexes, Write_ND_indexes, Initialize_List_indexes

  INTERFACE Initialize_ND_indexes
    MODULE PROCEDURE MolecCav_Initialize_ND_indexes
  END INTERFACE
  INTERFACE Initialize_List_indexes_old
    MODULE PROCEDURE MolecCav_Initialize_List_indexes_old
  END INTERFACE
  INTERFACE Initialize_List_indexes
    MODULE PROCEDURE MolecCav_Initialize_List_indexes
  END INTERFACE
  INTERFACE Increment_indexes
    MODULE PROCEDURE MolecCav_Increment_indexes
  END INTERFACE
  INTERFACE Deallocate_ND_indexes
    MODULE PROCEDURE MolecCav_Deallocate_ND_indexes
  END INTERFACE
  INTERFACE Write_ND_indexes
    MODULE PROCEDURE MolecCav_Write_ND_indexes
  END INTERFACE
    

  CONTAINS

  
  SUBROUTINE MolecCav_Initialize_ND_indexes(ND_indexes, Ranks_sizes, Starting_indexes, Begin_right, Debug)
    USE QDUtil_m
    IMPLICIT NONE 

    TYPE(ND_indexes_t), intent(inout) :: ND_indexes
    integer,            intent(in)    :: Ranks_sizes(:)                                            ! the basis set size for each dimension of the tensor
    integer, optional,  intent(in)    :: Starting_indexes(:)                                       ! the labelling of the basis vectors may starts from 0 or 1 each
    integer, optional,  intent(in)    :: Begin_right
    logical, optional,  intent(in)    :: Debug

    logical                           :: Debug_local = .FALSE.
    integer                           :: NB
    
    IF (PRESENT(Debug)) Debug_local = Debug

    IF (Debug_local) THEN
      WRITE(out_unit,*) 
      WRITE(out_unit,*) '********************************************************************************'
      WRITE(out_unit,*) '********************** INITIALIZING ND_INDEXES PARAMETERS **********************'
    END IF 

    IF (PRESENT(Starting_indexes) .AND. Size(Starting_indexes) /= Size(Ranks_sizes)) THEN
      WRITE(out_unit,*) "The sizes of the optional argument Starting_index ("//TO_string(Size(Starting_indexes))//") and of the D&
                        &im_sizes argument ("//TO_string(Size(Ranks_sizes))//") do not match. Please check the arguments."
      STOP "### The sizes of the Starting_index and the Ranks_sizes arguments do not match. cf. output file for more information."
    ELSE IF (PRESENT(Starting_indexes) .AND. Size(Starting_indexes) == Size(Ranks_sizes)) THEN
      ALLOCATE(ND_indexes%Starting_indexes(Size(Ranks_sizes)))
      ND_indexes%Starting_indexes = Starting_indexes
    ELSE 
      ALLOCATE(ND_indexes%Starting_indexes(Size(Ranks_sizes)))      
      ND_indexes%Starting_indexes = 1                                                              ! as a default arbitrary choice, the labelling of the basis vectors starts as the Fortran convention i.e. from 1
    END IF

    IF (PRESENT(Begin_right)) ND_indexes%Begin_right = Begin_right
    ND_indexes%N_dim       = Size(Ranks_sizes)
    ALLOCATE(ND_indexes%Ranks_sizes(ND_indexes%N_dim))
    ND_indexes%Ranks_sizes = Ranks_sizes
    ND_indexes%NB          = PRODUCT(ND_indexes%Ranks_sizes)

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


  FUNCTION MolecCav_Initialize_List_indexes(ND_indexes) RESULT(List_indexes)
    USE QDUtil_m
    IMPLICIT NONE 

    integer, allocatable           :: List_indexes(:)                                              ! the current values of the indexes for each dimension
    TYPE(ND_indexes_t), intent(in) :: ND_indexes

    integer                        :: Zero_index

    List_indexes = ND_indexes%Starting_indexes                                                     ! dynamic allocattion : allows to call the function for any allocatable object, already allocated or not, and even for a tabular not declared allocatable it seems

    Zero_index = 1 + ND_indexes%Begin_right*(ND_indexes%N_dim-1)
    List_indexes(Zero_index) = List_indexes(Zero_index) - 1                                        ! to prepare the initial point I = 0 (state before first loop)

  END FUNCTION MolecCav_Initialize_List_indexes


  SUBROUTINE MolecCav_Initialize_List_indexes_old(List_indexes, ND_indexes)
    USE QDUtil_m
    IMPLICIT NONE 

    integer, allocatable, intent(inout) :: List_indexes(:)                                         ! the current values of the indexes for each dimension
    TYPE(ND_indexes_t),   intent(in)    :: ND_indexes

    logical, parameter                  :: Debug_local = .FALSE.
    
    IF (Debug_local) THEN
      WRITE(out_unit,*) 
      WRITE(out_unit,*) '********************************************************************************'
      WRITE(out_unit,*) '*************************** INITIALIZING LIST_INDEXES **************************'
    END IF 
    
    ALLOCATE(List_indexes(ND_indexes%N_dim))

!    IF (Size(List_indexes) /= Size(ND_indexes%Starting_indexes)) THEN                             ! meaningless since Starting_indexes is built of size N_dim 
!      WRITE(out_unit,*) "The size of the indexes list ("//TO_string(Size(List_indexes))//") and the size of the starting values o&
!                        &f these indexes ("//TO_string(Size(ND_indexes%Starting_indexes))//") does not match"
!      STOP "### The size of the List_indexes and of the Starting_values does not match. cf. output file for more information."
!    ELSE
    List_indexes = ND_indexes%Starting_indexes
!    END IF

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "----------List_indexes constructed by MolecCav_Initialize_List_indexes----------"
      CALL Write_Vec(List_indexes, out_unit, Size(List_indexes), info="List_indexes")
      WRITE(out_unit,*) "--------End List_indexes constructed by MolecCav_Initialize_List_indexes--------"
      WRITE(out_unit,*) 
      WRITE(out_unit,*) '**************************** LIST_INDEXES INITIALIZED **************************'
      WRITE(out_unit,*) '********************************************************************************'      
    END IF
      
  END SUBROUTINE MolecCav_Initialize_List_indexes_old


  SUBROUTINE MolecCav_Increment_indexes(Continue_loop, List_indexes, ND_indexes, Debug)
    USE QDUtil_m
    IMPLICIT NONE 

    logical,            intent(inout) :: Continue_loop
    integer,            intent(inout) :: List_indexes(:)
    TYPE(ND_indexes_t), intent(inout) :: ND_indexes
    logical, optional,  intent(in)    :: Debug 

    logical                           :: Debug_local = .FALSE.
    integer                           :: N, i

    !/!\ WE DO NOT CONSIDER THE EFFECT OF STARTING INDEX FOR NOW /!\
    IF (PRESENT(Debug)) Debug_local = Debug

    N = ND_indexes%N_dim
    Continue_loop = .TRUE.
    i = 1 + ND_indexes%Begin_right*(N-1)                                                           ! we chosed first i from 0 to N using the index N-i in the loop and the final test i==N (so looping on dim=N then N-1 ... up to 1), but we are forced to loop in the other sens since the construct TotH has to loop on i_M first then i_C, and so does mappping 2DTO1D

    IF (Debug_local) WRITE(out_unit,*)
    DO WHILE (Continue_loop)
      List_indexes(i) = List_indexes(i) + 1
      IF (Debug_local) CALL Write_Vec(List_indexes, out_unit, N, info="-> incremented List_indexes")

      IF (List_indexes(i) > ND_indexes%Ranks_sizes(i)) THEN                                        ! i<N <=> N-i>0 (when overcome the last dimension, i is still incremented => N-i goes to 0 at the last iteration)
        List_indexes(i) = ND_indexes%Starting_indexes(i)
        i = i + 1 - 2*ND_indexes%Begin_right
        IF (ALL(List_indexes == ND_indexes%Starting_indexes) .AND. i==(N+1 - (N+1)*ND_indexes%Begin_right)) THEN
          IF (Debug_local) CALL Write_Vec(List_indexes, out_unit, N, info="-> incremented List_indexes")
          List_indexes = ND_indexes%Ranks_sizes
          Continue_loop = .FALSE.
          IF (Debug_local) CALL Write_Vec(List_indexes, out_unit, N, info="-> final List_indexes")
        END IF
      ELSE 
        EXIT
      END IF
    END DO

  END SUBROUTINE MolecCav_Increment_indexes


  SUBROUTINE MolecCav_Deallocate_ND_indexes(ND_indexes)
    USE QDUtil_m
    IMPLICIT NONE 

    TYPE(ND_indexes_t) :: ND_indexes

    logical, parameter :: Debug_local = .FALSE.

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "------------------Deallocating the following ND_indexes object------------------"
      CALL Write_ND_indexes(ND_indexes)
    END IF

    IF (ALLOCATED(ND_indexes%Ranks_sizes))        DEALLOCATE(ND_indexes%Ranks_sizes)
    IF (ALLOCATED(ND_indexes%Starting_indexes)) DEALLOCATE(ND_indexes%Starting_indexes)

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "------------------------The deallocated ND_indexes object-----------------------"
      CALL Write_ND_indexes(ND_indexes)
      WRITE(out_unit,*) "----------------End Deallocating the previous ND_indexes object---------------"
    END IF

  END SUBROUTINE MolecCav_Deallocate_ND_indexes


  SUBROUTINE MolecCav_Write_ND_indexes(ND_indexes)
    USE QDUtil_m
    IMPLICIT NONE

    TYPE(ND_indexes_t), intent(in) :: ND_indexes

    WRITE(out_unit,*) "___________________________________The associated ND_indexes object_________________________________"
    WRITE(out_unit,*) "|Number of dimension of the tensor to map (ND_indexes%N_dim)                 | ", ND_indexes%N_dim
    WRITE(out_unit,*) "|____________________________________________________________________________|_____________________|"
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
    IF (ALLOCATED(ND_indexes%Ranks_sizes)) THEN
      WRITE(out_unit,*) "|The list of the basis set size of each of these dimensions do is allocated..| "
      WRITE(out_unit,*) "|... and is initialized as (ND_indexes%Ranks_sizes)                          | "
      WRITE(out_unit,*) "|____________________________________________________________________________| "
      CALL Write_Vec(ND_indexes%Ranks_sizes, out_unit, Size(ND_indexes%Ranks_sizes), info="ND_indexes%Ranks_sizes")
      FLUSH(out_unit)
      WRITE(out_unit,*) "|____________________________________________________________________________|______________________"
    ELSE 
      WRITE(out_unit,*) "|The list of the basis set size of each of these dimensions is NOT allocated.| "
      WRITE(out_unit,*) "|____________________________________________________________________________|______________________"
    END IF
    FLUSH(out_unit)
    WRITE(out_unit,*) "|The dimension of the tensor producted vector (ND_indexes%NB)                | ", ND_indexes%NB
    WRITE(out_unit,*) "|____________________________________________________________________________|_____________________|"
    FLUSH(out_unit)
  
  END SUBROUTINE MolecCav_Write_ND_indexes


END MODULE
