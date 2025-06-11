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

  PUBLIC Transition_spectrum_t, Initialize, Alloc, Write, Dealloc

  INTERFACE Initialize
    MODULE PROCEDURE MolecCav_Initialize_transition_matrix_0K
  END INTERFACE
  INTERFACE Alloc
    MODULE PROCEDURE MolecCav_Allocate_transition_matrix_0K
  END INTERFACE
  INTERFACE Write
    MODULE PROCEDURE MolecCav_Write_transition_matrix_0K
  END INTERFACE
  INTERFACE Dealloc
    MODULE PROCEDURE MolecCav_Deallocate_transition_matrix_0K
  END INTERFACE


  CONTAINS


  SUBROUTINE MolecCav_Initialize_transition_matrix_0K(TranSpec, REigvec, DipMomt, REigval, E_threshold, Nb_states, Verbose, Debug)   ! /!\ FOR NOW DESIGNED FOR 1p1D
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    USE Algebra_m
    USE Sum_of_products_m
    IMPLICIT NONE
  
    TYPE(Transition_spectrum_t), intent(inout) :: TranSpec
    real(kind=Rkind),            intent(in)    :: REigvec(:,:)
    TYPE(Sum_of_products_t),     intent(in)    :: DipMomt
    real(kind=Rkind),            intent(in)    :: REigval(:)
    real(kind=Rkind), optional,  intent(in)    :: E_threshold ! for the state selection
    integer,          optional,  intent(in)    :: Nb_states   ! for the state selection
    integer, optional,           intent(in)    :: Verbose
    logical, optional,           intent(in)    :: Debug

    real(kind=Rkind), allocatable              :: DipMomt_GS(:)  ! action of the dipole moment upon the ground state
    integer                                    :: I
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
      CALL Write(TranSpec)
      WRITE(out_unit,*) "The <<REigvec>> argument :" 
      CALL Write_Mat(REigvec, out_unit, SIZE(REigvec, dim=2), info="REigvec")
      WRITE(out_unit,*) "The <<DipMomt>> argument :" 
      CALL Write(DipMomt)
      WRITE(out_unit,*) "The <<REigval>> argument :" 
      CALL Write_vec(REigval, out_unit, SIZE(REigval, dim=1), info="REigval")
      IF (PRESENT(E_threshold)) WRITE(out_unit,*) "The <<E_threshold>> argument : "//TO_string(E_threshold)
      IF (PRESENT(Nb_states)) WRITE(out_unit,*) "The <<Nb_states>> argument : "//TO_string(Nb_states)
      WRITE(out_unit,*) "--- End arguments of MolecCav_Allocate_transition_matrix_0K"
      FLUSH(out_unit)
    END IF
    
    !-------------------------------------------------Allocating the transition spectrum-------------------------------------------
    IF (PRESENT(E_threshold) .AND. PRESENT(Nb_states)) THEN
      CALL Alloc(TranSpec, E_threshold, REigval=REigval, Nb_states=Nb_states, Verbose=Verbose_local, Debug=Debug_local)
    ELSE IF (PRESENT(E_threshold)) THEN
      CALL Alloc(TranSpec, Energy_threshold=E_threshold, REigval=REigval, Verbose=Verbose_local, Debug=Debug_local)
    ELSE IF (PRESENT(Nb_states)) THEN
      CALL Alloc(TranSpec, Nb_states=Nb_states, Verbose=Verbose_local, Debug=Debug_local)
    ELSE
      CALL Alloc(TranSpec, Verbose=Verbose_local, Debug=Debug_local)
    END IF

    !-------------------------------------------------Computing the transition spectrum-------------------------------------------
    ALLOCATE(DipMomt_GS(SIZE(REigvec, dim=1)))
    CALL Action(DipMomt_GS, DipMomt, REigvec(:,1), Verbose=Verbose_local, Debug=.FALSE.)

    DO I = 2, TranSpec%N_trstns + 1 ! "+1" because if there is N_trstns to compute, then the first one is 2 <- 1 (G.S.) and the last one N_trstns+1 <- 1
      IF (Debug_local) WRITE(out_unit,*) "Transition 1 --> \overrightarrow{VP}_"//TO_string(I); FLUSH(out_unit)
      IF (Debug_local) CALL Write_Vec(TranSpec%tab_energies, out_unit, 10, info="energies"); FLUSH(out_unit)
      IF (Debug_local) CALL Write_Vec(REigval, out_unit, 10, info="eig val"); FLUSH(out_unit)
      
      TranSpec%tab_energies(I-1) = REigval(I) - REigval(1)  
      
      CALL Scalar_product(TranSpec%tab_ints(I-1), REigvec(:,I), DipMomt_GS)
      TranSpec%tab_ints(I-1) = ABS(TranSpec%tab_ints(I-1))**2 !* TranSpec%tab_energies(I-1) * 4.35974E-18 / (6.62607015E-34*2.99792458E10)
      
      IF (Debug_local) WRITE(out_unit,*) "Transition \overrightarrow{VP}_"//TO_string(1)//" --> \overrightarrow{VP}_"//TO_strin&
      &g(I)//" = "//TO_string(TranSpec%tab_ints(I-1))
    END DO 
    
    IF (Debug_local) WRITE(out_unit,*)
    IF (Debug_local) CALL Write_Vec(TranSpec%tab_energies, out_unit, Size(TranSpec%tab_energies), info="Transition energies")
    IF (Debug_local) CALL Write_Vec(TranSpec%tab_ints, out_unit, Size(TranSpec%tab_ints), info="Transition Intensities")
    
  END SUBROUTINE MolecCav_Initialize_transition_matrix_0K
  

  SUBROUTINE MolecCav_Allocate_transition_matrix_0K(TranSpec, Energy_threshold, REigval, Nb_states, Verbose, Debug)
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
    IF (Debug_local) WRITE(out_unit,*) "-------------------------------------------------ALLOCATING THE TRANSITION SPECTRUM OBJ&
                                              &ECT-------------------------------------------------"; FLUSH(out_unit)

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "--- Arguments of MolecCav_Allocate_transition_matrix_0K :"
      WRITE(out_unit,*) "The <<TranSpec>> argument :"
      CALL Write(TranSpec)
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
      STOP "### Missing REigval argument in MolecCav_Allocate_transition_matrix_0K"

    !-------------------------------------------------Counting the number of states to take into account-------------------------------------------
    ELSE IF (PRESENT(Energy_threshold)) THEN
      TranSpec%N_trstns = COUNT(Reigval < Reigval(1) + Energy_threshold) - 1 ! "-1" because G.S. <- G.S. is not a transition
    
    ELSE
        TranSpec%N_trstns = 0
    END IF 
    
    IF (PRESENT(Nb_states)) THEN
      IF (TranSpec%N_trstns == 0 .OR. Nb_states < TranSpec%N_trstns) TranSpec%N_trstns = Nb_states - 1 ! "-1" because G.S. <- G.S. is not a transition
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
      TranSpec%N_trstns = Size(REigval, dim=1) - 1 ! "-1" because G.S. <- G.S. is not a transition
    END IF 

    IF (Debug_local) WRITE(out_unit,*)
    IF (Debug_local) WRITE(out_unit,*) TO_string(TranSpec%N_trstns + 1)//" States will be taken into account to compute the trans&
    &ition intensities."

    ALLOCATE(TranSpec%tab_ints(    TranSpec%N_trstns))
    ALLOCATE(TranSpec%tab_energies(TranSpec%N_trstns))
    TranSpec%tab_ints     = ZERO
    TranSpec%tab_energies = ZERO
    
    IF (Debug_local) WRITE(out_unit,*)
    IF (Debug_local) CALL Write(TranSpec)
    FLUSH(out_unit)

    IF (Debug_local) WRITE(out_unit,*) 
    IF (Debug_local) WRITE(out_unit,*) "-------------------------------------------------TRANSITION SPECTRUM OBJECT ALLOCATED----&
                                              &---------------------------------------------"; FLUSH(out_unit)

  END SUBROUTINE MolecCav_Allocate_transition_matrix_0K
  

  SUBROUTINE MolecCav_Write_transition_matrix_0K(TranSpec)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
    USE QDUtil_m
    IMPLICIT NONE 
    
    TYPE(Transition_spectrum_t), intent(in) :: TranSpec


    WRITE(out_unit,*)
    WRITE(out_unit,*) "_________________________________The Transition_spectrum object__________________________________"
    WRITE(out_unit,*) "|The number of transitions considered (TranSpec%N_trstns)                      | "//TO_string(TranSpec%N_t&
    &rstns)
    WRITE(out_unit,*) "|______________________________________________________________________________|______________________"
    FLUSH(out_unit)
    IF (ALLOCATED(TranSpec%tab_energies)) THEN
      WRITE(out_unit,*) "|The energies associated to the transitions (TranSpec%tab_energies) :          |"
      CALL Write_Vec(TranSpec%tab_energies, out_unit, SIZE(TranSpec%tab_energies), info="Transition_Energies")
    ELSE 
      WRITE(out_unit,*) "|The transition energies list is not allocated (TranSpec%tab_energies)         |"
    END IF
    WRITE(out_unit,*) "|______________________________________________________________________________|"
    FLUSH(out_unit)
    IF (ALLOCATED(TranSpec%tab_ints)) THEN
      WRITE(out_unit,*) "|The associated list of intensities (TranSpec%tab_ints) :                      |"
      CALL Write_Vec(TranSpec%tab_ints, out_unit, SIZE(TranSpec%tab_ints), info="TranSpec%tab_ints")
    ELSE 
      WRITE(out_unit,*) "|The associated list of intensities is not allocated (TranSpec%tab_ints)       |"
    END IF
    FLUSH(out_unit)
    WRITE(out_unit,*) "|_____________________________________End ND Operator object___________________|"
    FLUSH(out_unit)

  END SUBROUTINE MolecCav_Write_transition_matrix_0K


  SUBROUTINE MolecCav_Deallocate_transition_matrix_0K(TranSpec, Verbose, Debug)
    USE QDUtil_m
    IMPLICIT NONE 

    TYPE(Transition_spectrum_t), intent(inout) :: TranSpec
    integer, optional,       intent(in)        :: Verbose                                                                                 ! cf. comments in HO1D_parameters_m
    logical, optional,       intent(in)        :: Debug                                                                                   ! cf. comments in HO1D_parameters_m

    integer                                    :: Verbose_local                                                                      ! goes from 25 (= 0 verbose) to 29 (= maximum verbose) at this layer
    logical                                    :: Debug_local

    !------------------------------------------------------Debugging options-----------------------------------------------------
    IF (PRESENT(Verbose)) THEN; Verbose_local = Verbose
    ELSE; Verbose_local = 20; END IF 
    IF (PRESENT(Debug))   THEN; Debug_local   = Debug
    ELSE; Debug_local = .FALSE.; END IF

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "--- The TranSpec to be deallocated :"
      CALL Write(TranSpec)
      WRITE(out_unit,*) "--- End TranSpec to be deallocated"
    END IF 

    !-----------------------------Deallocating the transition spectrum object----------------------------
    IF (Debug_local) WRITE(out_unit,*)
    IF (Debug_local) WRITE(out_unit,*) "-----------------------------------------------Deallocating the Transition_spectrum objec&
    &t----------------------------------------------"
  
    TranSpec%N_trstns = 0
    IF (ALLOCATED(TranSpec%tab_energies)) DEALLOCATE(TranSpec%tab_energies) 
    IF (ALLOCATED(TranSpec%tab_ints)    ) DEALLOCATE(TranSpec%tab_ints) 

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "--- The TranSpec object after having been deallocated :"
      CALL Write(TranSpec)
      WRITE(out_unit,*) "--- End dellocating TranSpec"
    END IF

  END SUBROUTINE MolecCav_Deallocate_transition_matrix_0K


END MODULE