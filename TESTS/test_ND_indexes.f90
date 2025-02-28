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
PROGRAM test_ND_indexes
  !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
  USE QDUtil_m
  USE QDUtil_Test_m
  USE ND_indexes_m
  IMPLICIT NONE
   

  logical, parameter            :: Debug = .TRUE.
  

  !-----------------------------The indexes objects----------------------------
  integer, allocatable          :: Starting_indexes(:)
  integer                       :: Ranks_sizes(3) = [3,5,4]
  TYPE(ND_indexes_t)            :: ND_indexes

  integer, allocatable          :: List_indexes(:)
  
  !-----------------------------------Results----------------------------------
  integer, allocatable          :: Table_indexes(:,:)                          ! indexes history through incrementation
  integer, allocatable          :: Table_indexes_ref(:,:)

  !----------------------------------Utilities---------------------------------
  integer                       :: N_dim, NB, I, i_1, i_2, i_3
  logical                       :: Continue_loop = .TRUE.
  
  TYPE(test_t)                  :: test_ND_ind
  logical                       :: error_ND_ind = .FALSE.
    
  
  !---------------------------------------Test initialization--------------------------------------
  CALL Initialize_Test(test_ND_ind, test_name="OUT/test_file_ND_ind")
  

  !--------------------------------------System initialization-------------------------------------
  N_dim = Size(Ranks_sizes)
  ALLOCATE(Starting_indexes(N_dim))
  Starting_indexes = 1
  IF (Debug) WRITE(out_unit,*) "Starting_indexes", Starting_indexes

  CALL Initialize_ND_indexes(ND_indexes, Ranks_sizes, Starting_indexes=Starting_indexes, Debug=Debug)

  List_indexes = Initialize_List_indexes(ND_indexes)
  IF (Debug) WRITE(out_unit,*)
  IF (Debug) WRITE(out_unit,*) "List_indexes", List_indexes

  ALLOCATE(Table_indexes_ref(ND_indexes%N_dim, ND_indexes%NB))
  I = 0
  DO i_3 = 1, Ranks_sizes(3)
    DO i_2 = 1, Ranks_sizes(2)
      DO i_1 = 1, Ranks_sizes(1)
        I = I + 1
        Table_indexes_ref(:,I) = [i_1, i_2, i_3]
      END DO 
    END DO 
  END DO 
  Table_indexes_ref(1,2) = 10


  !----------------------------------------Tests ND_indexes----------------------------------------
    !---------------------------------The ND_indexes initialization--------------------------------
  CALL Equal_I_I_scalar(error_ND_ind, N_dim, ND_indexes%N_dim)
  CALL Logical_Test(test_ND_ind, error_ND_ind, test2=.FALSE., info="initialization ND_indexes%N_dim")
  IF (error_ND_ind .AND. Debug) THEN
    WRITE(out_unit,*)
    WRITE(out_unit,*) "N_dim = "//TO_string(N_dim)//"; ND_indexes%N_dim = "//TO_string(ND_indexes%N_dim)
  END IF

  IF (.NOT. ALLOCATED(ND_indexes%Starting_indexes)) THEN
    WRITE(out_unit,*) "ND_indexes%Starting_indexes failed to initialize : not allocated. Please check initialization."
    STOP "### ND_indexes%Starting_indexes failed to initialize : not allocated. cf. output file."
  END IF
  CALL Equal_I_I_vector(error_ND_ind, Starting_indexes, ND_indexes%Starting_indexes)
  CALL Logical_Test(test_ND_ind, error_ND_ind, test2=.FALSE., info="initialization ND_indexes%Starting_indexes")
  IF (error_ND_ind .AND. Debug) THEN
    WRITE(out_unit,*)
    CALL Write_Vec(Starting_indexes, out_unit, Size(Starting_indexes), info="Starting_indexes")
    WRITE(out_unit,*)
    CALL Write_Vec(ND_indexes%Starting_indexes, out_unit, Size(Starting_indexes), info="ND_indexes%Starting_indexes")
  END IF

  IF (.NOT. ALLOCATED(ND_indexes%Ranks_sizes)) THEN
    WRITE(out_unit,*) "ND_indexes%Ranks_sizes failed to initialize : not allocated. Please check initialization."
    STOP "### ND_indexes%Starting_indexes failed to initialize : not allocated. cf. output file."
  END IF
  CALL Equal_I_I_vector(error_ND_ind, Ranks_sizes, ND_indexes%Ranks_sizes)
  CALL Logical_Test(test_ND_ind, error_ND_ind, test2=.FALSE., info="initialization ND_indexes%Ranks_sizes")
  IF (error_ND_ind .AND. Debug) THEN
    WRITE(out_unit,*)
    CALL Write_Vec(Ranks_sizes, out_unit, Size(Ranks_sizes), info="Ranks_sizes")
    WRITE(out_unit,*)
    CALL Write_Vec(ND_indexes%Ranks_sizes, out_unit, Size(Ranks_sizes), info="ND_indexes%Ranks_sizes")
  END IF

  CALL Equal_I_I_scalar(error_ND_ind, PRODUCT(Ranks_sizes), ND_indexes%NB)
  CALL Logical_Test(test_ND_ind, error_ND_ind, test2=.FALSE., info="initialization ND_indexes%NB")
  IF (error_ND_ind .AND. Debug) THEN
    WRITE(out_unit,*)
    WRITE(out_unit,*) "N_dim = "//TO_string(PRODUCT(Ranks_sizes))//"; ND_indexes%N_dim = "//TO_string(ND_indexes%NB)
  END IF

  CALL Equal_I_I_vector(error_ND_ind, List_indexes, [(1, i = 1, N_dim)])
  CALL Logical_Test(test_ND_ind, error_ND_ind, test2=.FALSE., info="List_indexes well initialized ?")
  IF (error_ND_ind .AND. Debug) THEN
    WRITE(out_unit,*)
    WRITE(out_unit,*) "List_indexes", List_indexes
  END IF    


    !---------------------------------The ND_indexes incrementation--------------------------------
  ALLOCATE(Table_indexes(ND_indexes%N_dim, ND_indexes%NB))

  I = 0
  
  Continue_loop = .TRUE.
  DO WHILE (Continue_loop)
    IF (Debug) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "After looping "//TO_string(I)//" : I = "//TO_string(I)//"; NB = "//TO_string(ND_indexes%NB)
      CALL Write_Vec(List_indexes,           out_unit, Size(List_indexes),           info="List_indexes")
      CALL Write_Vec(ND_indexes%Ranks_sizes, out_unit, Size(ND_indexes%Ranks_sizes), info="Ranks_sizes")
    END IF

    I = I + 1      
    IF (Debug) WRITE(out_unit,*)
    IF (Debug) WRITE(out_unit,*) "Looping "//TO_string(I)//" -> I = "//TO_string(I)//" :"
    IF (I>ND_indexes%NB) THEN
      WRITE(out_unit,*) "I should not be greater that NB !"
      STOP              "I should not be greater that NB !"
    END IF

    Table_indexes(:,I) = List_indexes
      
    CALL Increment_indexes(Continue_loop, List_indexes, ND_indexes, Debug=Debug)
    IF (Debug) WRITE(out_unit,*) " -> Continue_loop : "//TO_string(Continue_loop)
  END DO

  CALL Equal_I_I_matrix(error_ND_ind, Table_indexes_ref, Table_indexes)
  CALL Logical_Test(test_ND_ind, error_ND_ind, test2=.FALSE., info="indexes history through incrementation")
  IF (error_ND_ind .AND. Debug) THEN
    WRITE(out_unit,*)
    CALL Write_Mat(Table_indexes, out_unit, 10, info="Table_indexes")
    WRITE(out_unit,*)
    CALL Write_Mat(Table_indexes_ref, out_unit, 10, info="Table_indexes_ref")
  END IF


    !----------------------------------The ND_indexes deallocation---------------------------------
  CALL Deallocate_ND_indexes(ND_indexes)
  IF (Debug) CALL Write_ND_indexes(ND_indexes)

  CALL Logical_Test(test_ND_ind, ALLOCATED(ND_indexes%Starting_indexes),test2=.FALSE.,info="ND_ind%Starting_indexes deallocated ?")
  CALL Logical_Test(test_ND_ind, ALLOCATED(ND_indexes%Ranks_sizes),     test2=.FALSE.,info="ND_ind%Ranks_sizes      deallocated ?")


  CALL Finalize_Test(test_ND_ind)
  
    
  CONTAINS
  
  
  SUBROUTINE Equal_I_I_scalar(error, Int_1, Int_2)
    USE QDUtil_m
    IMPLICIT NONE 

    logical, intent(inout) :: error
    integer, intent(in)    :: Int_1
    integer, intent(in)    :: Int_2
    
    logical, parameter     :: Debug_local = .FALSE.

    IF (Int_1 /= Int_2) THEN
      error = .TRUE.
      IF (Debug_local) WRITE(out_unit,*) "The two integers are not equal : Int_1 =", Int_1, "Int_2 =", Int_2
    ELSE 
      error = .FALSE.
      IF (Debug_local) WRITE(out_unit,*) "The two numbers are close enough to be considered equal : Int_1 =", Int_1, "Int_2 =", I&
                                         &nt_2
    END IF

  END SUBROUTINE Equal_I_I_scalar
  

  SUBROUTINE Equal_I_I_vector(error, Int_1, Int_2)
    USE QDUtil_m
    IMPLICIT NONE 

    logical, intent(inout) :: error
    integer, intent(in)    :: Int_1(:)
    integer, intent(in)    :: Int_2(:)
    
    logical, parameter     :: Debug_local = .FALSE.
    integer                :: Nb_1

    Nb_1 = Size(Int_1)
    IF (Nb_1 /= Size(Int_2)) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "The two vectors must have same dimensions to compare them. Please, check initialization."
      WRITE(out_unit,*) "Size(Int_1) = "//TO_string(Size(Int_1))//"; Size(Int_2) = "//TO_string(Size(Int_2))
      STOP "### The two vectors must have same dimensions to compare them. Please, check initialization."
    END IF 

    IF (ANY(Int_1 /= Int_2)) THEN
      error = .TRUE.
      IF (Debug_local) THEN
        WRITE(out_unit,*) "The two vectors are not equal :"
        CALL Write_Vec(Int_1, out_unit, Nb_1, info="Int_1(:)")
        CALL Write_Vec(Int_2, out_unit, Nb_1, info="Int_2(:)")
      END IF 

    ELSE 
      error = .FALSE.
      IF (Debug_local) THEN
        WRITE(out_unit,*) "The two vectors are equal :"
        CALL Write_Vec(Int_1, out_unit, Nb_1, info="Int_1(:)")
        CALL Write_Vec(Int_2, out_unit, Nb_1, info="Int_2(:)")
      END IF
    END IF 

  END SUBROUTINE Equal_I_I_vector
  

  SUBROUTINE Equal_I_I_matrix(error, Int_1, Int_2)
    USE QDUtil_m
    IMPLICIT NONE 
  
    logical, intent(inout) :: error
    integer, intent(in)    :: Int_1(:,:)
    integer, intent(in)    :: Int_2(:,:)
    
    logical, parameter              :: Debug_local = .FALSE.
    integer                         :: Nb_1_local, Nb_2_local

    Nb_1_local = Size(Int_1, dim=1)
    Nb_2_local = Size(Int_2, dim=2)
    IF (Nb_1_local /= Size(Int_1, dim=1) .OR. Nb_2_local /= Size(Int_2, dim=2)) THEN
      WRITE(out_unit,*) "The two matrices must have same dimensions to compare them. Please, check initialization."
      STOP "### The two matrices must have same dimensions to compare them. Please, check initialization."
    END IF 
  
    IF (ANY(Int_1 /= Int_2)) THEN
      error = .TRUE.
      IF (Debug_local) THEN
        WRITE(out_unit,*) "The two matrices are not equal :"
        CALL Write_Mat(Int_1, out_unit, Nb_2_local, info="Int_1(:,:)")
        CALL Write_Mat(Int_2, out_unit, Nb_2_local, info="Int_2(:,:)")
      END IF 

    ELSE 
      error = .FALSE.
      IF (Debug_local) THEN
        WRITE(out_unit,*) "The two matrices are considered equal :"
        CALL Write_Mat(Int_1, out_unit, Nb_2_local, info="Int_1(:,:)")
        CALL Write_Mat(Int_2, out_unit, Nb_2_local, info="Int_2(:,:)")
      END IF
    END IF 

  END SUBROUTINE Equal_I_I_matrix


END PROGRAM
  