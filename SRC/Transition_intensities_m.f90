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
MODULE Transition_spectrum_m
  !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
  USE QDUtil_m
  USE Algebra_m
  USE Sum_of_products_m
  IMPLICIT NONE

  
  TYPE                            :: Transition_spectrum_t                  ! N_mat, N_cav canNOT be part of the derived type because they are not specific of one operator_ND, but parameters of the whole system/calculation. All OpND will have the same. (therefore only one namelist per mode is needed, and not one per mode and oOND)
    integer                       :: N_trstns = 0                           ! the number of terms in the sum of products operator
    real(kind=Rkind), allocatable :: tab_ints(:)                            ! the list of the products terms in the sum, which one being an OpND
    real(kind=Rkind), allocatable :: tab_energies(:)
    ! real(kind=Rkind), allocatable :: tab_ints_HighTemp(:,:)                            ! the list of the products terms in the sum, which one being an OpND
    ! real(kind=Rkind), allocatable :: tab_energies_HighTemp(:,:)
  END TYPE


  PRIVATE

  PUBLIC Transition_spectrum_t!, Initialize, Alloc

  INTERFACE Initialize
    MODULE PROCEDURE MolecCav_Initialize_transition_matrix_0K
  END INTERFACE
  INTERFACE Alloc
    MODULE PROCEDURE MolecCav_Allocate_transition_matrix_0K
  END INTERFACE


  CONTAINS


  SUBROUTINE MolecCav_Initialize_transition_matrix_0K(TranSpec, REigvec, DipMomt, E_threshold, REigval, Nb_states, Verbose, Debug)   ! /!\ FOR NOW DESIGNED FOR 1p1D
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    USE Algebra_m
    USE Sum_of_products_m
    IMPLICIT NONE
  
    TYPE(Transition_spectrum_t), intent(inout) :: TranSpec
    real(kind=Rkind),            intent(in)    :: REigvec(:,:)
    TYPE(Sum_of_products_t),     intent(in)    :: DipMomt
    real(kind=Rkind), optional,  intent(in)    :: E_threshold
    real(kind=Rkind), optional,  intent(in)    :: REigval(:)
    integer,          optional,  intent(in)    :: Nb_states
    integer, optional,           intent(in)    :: Verbose
    logical, optional,           intent(in)    :: Debug

    real(kind=Rkind), allocatable      :: InitPsi_1p1D(:,:)
    real(kind=Rkind), allocatable      :: FinPsi_1p1D(:,:)
    real(kind=Rkind), allocatable      :: Intermediary(:,:)  
    integer                            :: Nb_M, Nb_C, i_C, N, I, J
    integer                            :: Verbose_local
    logical                            :: Debug_local

    !###################### WE ARE HEREEEEEEEEEEEEEEEEEEEEEEEEEEEEE ######################
    IF (PRESENT(Debug)) Debug_local = Debug

    IF (MOD(Size(REigvec, dim=1),MatDipMomt%Nb) /= 0) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "Unconsistent arguments of Transition intensity : NB /= Nb_M*Integer"
      WRITE(out_unit,*) "NB = "//TO_string(Size(REigvec, dim=1))//";  Nb_M = "//TO_string(MatDipMomt%Nb)
      WRITE(out_unit,*) "Please check arguments"
      STOP "### Unconsistent arguments of Transition_intensity_matrix"
    END IF

    N    = Size(Intensities, dim=1)
    Nb_M = MatDipMomt%Nb             ! will be ND_indexes soon
    Nb_C = Size(REigvec, dim=1)/Nb_M

    IF (Debug_local) WRITE(out_unit,*)
    DO I = 1, N
      DO J =  1, N
        CALL Transition_intensity(Intensities(I,J), REigvec(:,I), MatDipMomt, REigvec(:,J), Nb_M, Nb_C)
        IF (Debug_local) WRITE(out_unit,*) "Transition \overrightarrow{VP}_"//TO_string(I)//" --> \overrightarrow{VP}_"//TO_strin&
                                           &g(J)//" = "//TO_string(Intensities(I,J))
      END DO
    END DO 
    
    IF (Debug_local) WRITE(out_unit,*)
    IF (Debug_local) CALL Write_Mat(Intensities, out_unit, Size(Intensities), info="Intensities matrix")
    
  END SUBROUTINE MolecCav_Initialize_transition_matrix_0K
  

  SUBROUTINE MolecCav_Allocate_transition_matrix_0K(TranSpec, Energy_threshold, REigval, Nb_states, Verbose, Debug)   ! /!\ FOR NOW DESIGNED FOR 1p1D
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    IMPLICIT NONE
  
    TYPE(Transition_spectrum_t), intent(inout) :: TranSpec
    real(kind=Rkind), optional,  intent(in)    :: Energy_threshold
    real(kind=Rkind), optional,  intent(in)    :: REigval(:)
    integer,          optional,  intent(in)    :: Nb_states
    integer,          optional,  intent(in)    :: Verbose
    logical,          optional,  intent(in)    :: Debug

    integer                                    :: Verbose_local
    logical                                    :: Debug_local

    !------------------------------------------------------Debugging options-----------------------------------------------------
    IF (PRESENT(Verbose)) THEN; Verbose_local = Verbose
    ELSE; Verbose_local = 20; END IF 
    IF (PRESENT(Debug))   THEN; Debug_local   = Debug
    ELSE; Debug_local = .FALSE.; END IF

    IF (Debug_local) WRITE(out_unit,*) 
    IF (Debug_local) WRITE(out_unit,*) "-------------------------------------------------INITIALIZING THE TRANSITION SPECTRUM OBJ&
                                              &ECT-------------------------------------------------"; FLUSH(out_unit)

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "--- Arguments of MolecCav_Allocate_transition_matrix_0K :"
      WRITE(out_unit,*) "The <<TranSpec>> argument :"
    !   CALL Write(TranSpec)
      IF (PRESENT(Energy_threshold)) WRITE(out_unit,*) "The <<Energy_threshold>> argument : "//TO_string(Energy_threshold)
      IF (PRESENT(REigval)) WRITE(out_unit,*) "The <<REigval>> argument : "
      IF (PRESENT(REigval)) CALL Write_Vec(REigval, out_unit, 1, info="REigval")
      IF (PRESENT(Nb_states)) WRITE(out_unit,*) "The <<Nb_states>> argument : "//TO_string(Nb_states)
      WRITE(out_unit,*) "--- End arguments of MolecCav_Allocate_transition_matrix_0K"
      FLUSH(out_unit)
    END IF
    
    IF (PRESENT(Energy_threshold) .AND. .NOT. PRESENT(REigval)) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "### The list of the total Hamiltonian Eigenenergies (REigval) is expected when the selection criterion i&
                        &s energy-based (Energy_threshold provided). Please check the arguments."
      STOP "### Missing REigval argument in Transition_intensity_matrix"

    !-------------------------------------------------Counting the number of states to take into account-------------------------------------------
    ELSE IF (PRESENT(Energy_threshold)) THEN
      TranSpec%N_trstns = COUNT(Reigval > Reigval(1) + Energy_threshold)
    
    ELSE
        TranSpec%N_trstns = 0
    END IF 
    
    IF (PRESENT(Nb_states)) THEN
      IF (TranSpec%N_trstns == 0 .OR. Nb_states < TranSpec%N_trstns) TranSpec%N_trstns = Nb_states
    END IF

    IF ((.NOT. PRESENT(Energy_threshold)) .AND. (.NOT. PRESENT(Nb_states))) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "########################## WARNING ########################## WARNING ########################## WARNING&
                       & #########################"
      WRITE(out_unit,*) "               No criterion provided to select the states with which the transition intensities have to &
                       &be computed"
      WRITE(out_unit,*) "                                              All Eigenstates will thus be considered"
      WRITE(out_unit,*) "########################## WARNING ########################## WARNING ########################## WARNING&
                       & #########################"
      TranSpec%N_trstns = Size(REigval, dim=1)
    END IF 

    IF (Debug_local) WRITE(out_unit,*)
    IF (Debug_local) WRITE(out_unit,*) TO_string(TranSpec%N_trstns)//" States will be taken into account to compute the transition intensities."

    ALLOCATE(TranSpec%tab_ints(    TranSpec%N_trstns))
    ALLOCATE(TranSpec%tab_energies(TranSpec%N_trstns))
    TranSpec%tab_ints     = ZERO
    TranSpec%tab_energies = ZERO
    
    IF (Debug_local) WRITE(out_unit,*)
    IF (Debug_local) CALL Write_Mat(Intensities, out_unit, Size(Intensities), info="Initialized intensity matrix")
    FLUSH(out_unit)

  END SUBROUTINE MolecCav_Allocate_transition_matrix_0K
  

END MODULE