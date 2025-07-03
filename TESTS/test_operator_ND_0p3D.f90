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
PROGRAM test_operator_ND_0p3D
  !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real124
  USE QDUtil_m
  USE Tests_m
  USE Algebra_m
  USE Operator_ND_m
  IMPLICIT NONE


  integer             :: Verbose = 50
  logical             :: Debug   = .TRUE.

  logical             :: Dense   = .FALSE.
  TYPE(Operator_ND_t) :: HxIxI
  TYPE(Operator_ND_t) :: IxHxI
  TYPE(Operator_ND_t) :: IxIxH
  real(kind=Rkind)    :: Cav1w
  real(kind=Rkind)    :: Cav2w
  real(kind=Rkind)    :: Cav3w
  integer             :: Cav1Nb
  integer             :: Cav2Nb
  integer             :: Cav3Nb
  integer             :: NB

  real(kind=Rkind), allocatable ::    Phi(:)                                                                             ! an any vector representing the excitation state/wavefunction = a linear combination of the basis functions from the canonical basis set over \mathbb{R} /!\ Not normalized yet !
  real(kind=Rkind), allocatable :: Op_phi(:)                                                                             ! an any vector representing the excitation state/wavefunction = a linear combination of the basis functions from the canonical basis set over \mathbb{R} /!\ Not normalized yet !
  real(kind=Rkind), allocatable :: TotH(:,:)                                                                             ! an any vector representing the excitation state/wavefunction = a linear combination of the basis functions from the canonical basis set over \mathbb{R} /!\ Not normalized yet !
  real(kind=Rkind), allocatable :: Eigenenergies(:)                                                                             ! an any vector representing the excitation state/wavefunction = a linear combination of the basis functions from the canonical basis set over \mathbb{R} /!\ Not normalized yet !
  real(kind=Rkind), allocatable :: Eigenstates(:,:)                                                                             ! an any vector representing the excitation state/wavefunction = a linear combination of the basis functions from the canonical basis set over \mathbb{R} /!\ Not normalized yet !

  real(kind=Rkind)              :: Modes_freq(3)
  real(kind=Rkind)              :: N_modes(3)

  TYPE(test_t)        :: test_opnd
  logical             :: error_opnd = .FALSE.

  integer             :: J, i_1, i_2, i_3, min_index


  !-----------------------------Test initialization----------------------------
  CALL Initialize_Test(test_opnd, test_name="OUT/tests/test_file_opnd_0p3D")


  !-------------------------Operator_ND object initialization-------------------------
  CALL Initialize(HxIxI,"" , "Hamiltonian, identity, identity ", in_unit, Dense=Dense, Verbose=Verbose, Debug=.FALSE.)
  IF (Debug) THEN
    WRITE(out_unit,*)
    WRITE(out_unit,*) "--------------Operator_ND object constructed by MolecCav_Initialize_operator_ND--------------"
    CALL Write(HxIxI)
    WRITE(out_unit,*) "------------End Operator_ND object constructed by MolecCav_Initialize_operator_ND------------"
  END IF

  CALL Initialize(IxHxI,"" , "identity, Hamiltonian, identity", in_unit, Dense=Dense, Verbose=Verbose, Debug=.FALSE.)
  IF (Debug) THEN
    WRITE(out_unit,*)
    WRITE(out_unit,*) "--------------Operator_ND object constructed by MolecCav_Initialize_operator_ND--------------"
    CALL Write(IxHxI)
    WRITE(out_unit,*) "------------End Operator_ND object constructed by MolecCav_Initialize_operator_ND------------"
  END IF

  CALL Initialize(IxIxH,"" , "identity, identity, Hamiltonian", in_unit, Dense=Dense, Verbose=Verbose, Debug=.FALSE.)
  IF (Debug) THEN
    WRITE(out_unit,*)
    WRITE(out_unit,*) "--------------Operator_ND object constructed by MolecCav_Initialize_operator_ND--------------"
    CALL Write(IxIxH)
    WRITE(out_unit,*) "------------End Operator_ND object constructed by MolecCav_Initialize_operator_ND------------"
  END IF


  !----------------------------Testing the initialization---------------------------
  CALL Logical_Test(test_opnd, SIZE(HxIxI%tab_indexes_mat_op)/=0,       test2=.FALSE., info="HxIxI%tab_mat_op")
  CALL Logical_Test(test_opnd, ANY( HxIxI%tab_indexes_cav_op/=[1,0,0]), test2=.FALSE., info="HxIxI%tab_cav_op")

  CALL Logical_Test(test_opnd, SIZE(IxHxI%tab_indexes_mat_op)/=0,       test2=.FALSE., info="IxHxI%tab_mat_op")
  CALL Logical_Test(test_opnd, ANY( IxHxI%tab_indexes_cav_op/=[0,1,0]), test2=.FALSE., info="IxHxI%tab_cav_op")

  CALL Logical_Test(test_opnd, SIZE(IxIxH%tab_indexes_mat_op)/=0,       test2=.FALSE., info="IxIxH%tab_mat_op")
  CALL Logical_Test(test_opnd, ANY( IxIxH%tab_indexes_cav_op/=[0,0,1]), test2=.FALSE., info="IxIxH%tab_cav_op")


  !----------------------------Testing the actions---------------------------
  CALL Get(Cav1w,  "w",  "Cavity", 1)
  CALL Get(Cav2w,  "w",  "Cavity", 2)
  CALL Get(Cav3w,  "w",  "Cavity", 3)
  CALL Get(Cav1Nb, "Nb", "Cavity", 1)
  CALL Get(Cav2Nb, "Nb", "Cavity", 2)
  CALL Get(Cav3Nb, "Nb", "Cavity", 3)
  NB = Cav1Nb * Cav2Nb * Cav3Nb

  IF (Debug) THEN
    WRITE(out_unit,*)
    WRITE(out_unit,*) "--- System parameters"
    WRITE(out_unit,*) "Cav1w  = "//TO_string(Cav1w)
    WRITE(out_unit,*) "Cav2w  = "//TO_string(Cav2w)
    WRITE(out_unit,*) "Cav3w  = "//TO_string(Cav3w)
    WRITE(out_unit,*) "Cav1Nb = "//TO_string(Cav1Nb)
    WRITE(out_unit,*) "Cav2Nb = "//TO_string(Cav2Nb)
    WRITE(out_unit,*) "Cav3Nb = "//TO_string(Cav3Nb)
    WRITE(out_unit,*) "NB     = "//TO_string(NB)
  END IF 

  ALLOCATE(Phi(NB))
  ALLOCATE(Op_phi(NB))
  ALLOCATE(TotH(NB,NB))
  ALLOCATE(Eigenenergies(NB))
  ALLOCATE(Eigenstates(NB,NB))
  TotH = ZERO
  DO J = 1, NB
    Phi = ZERO
    Phi(J) = ONE

    Op_phi = ZERO
    CALL Action(Op_phi, HxIxI, Phi, Verbose=Verbose, Debug=.FALSE.)
    TotH(:,J) = TotH(:,J) + Op_phi

    Op_phi = ZERO
    CALL Action(Op_phi, IxHxI, Phi, Verbose=Verbose, Debug=.FALSE.)
    TotH(:,J) = TotH(:,J) + Op_phi
    
    Op_phi = ZERO
    CALL Action(Op_phi, IxIxH, Phi, Verbose=Verbose, Debug=.FALSE.)
    TotH(:,J) = TotH(:,J) + Op_phi
  END DO

  IF (Debug) WRITE(out_unit,*)
  IF (Debug) CALL Write_Mat(TotH, out_unit, NB, info="TotH")

  CALL diagonalization(TotH, Eigenenergies, Eigenstates)
  IF (Debug) WRITE(out_unit,*)
  IF (Debug) CALL Write_Vec(Eigenenergies(1:4), out_unit, 1, info="EigenEnergies(1:4)")
  
  CALL Compute_Normal_modes(N_modes, Eigenenergies(1:4), Debug_opt=.FALSE.)
  
  Modes_freq = [Cav1w, Cav2w, Cav3w]

  DO J = 1, 3
    min_index = MINLOC(Modes_freq, dim=1)
      ! WRITE(out_unit,*) "MININDEX = "//TO_string(min_index)
    CALL Equal_tensor(error_opnd, Modes_freq(min_index), N_modes(J))
    CALL Logical_Test(test_opnd, error_opnd, test2=.FALSE., info="Normal mode "//TO_string(J))
    IF (Debug .OR. error_opnd) THEN
      WRITE(out_unit,*) "J, N_modes(J), Modes_freq(min_index) = "//TO_string(J)//", "//TO_string(N_modes(J))//&
      &", "//TO_string(Modes_freq(min_index))
    END IF
    Modes_freq(min_index) = HUGE(1)
  END DO


  !----------------------------Testing the writing---------------------------
  CALL Write(IxHxI)

  
  !----------------------------Testing the deallocation---------------------------
  CALL Dealloc(HxIxI, Verbose=Verbose, Debug=Debug)
  CALL Logical_Test(test_opnd, ALLOCATED(HxIxI%tab_indexes_mat_op), test2=.FALSE., info="HxIxI%tab_mat deallocated")
  CALL Logical_Test(test_opnd, ALLOCATED(HxIxI%tab_indexes_cav_op), test2=.FALSE., info="HxIxI%tab_cav deallocated")
  

  CALL Finalize_Test(test_opnd)


  CONTAINS


  SUBROUTINE Compute_Normal_modes(N_modes_loc, Eigenenergies_loc, Debug_opt)
    USE QDUtil_m
    IMPLICIT NONE 

    real(kind=Rkind),  intent(inout) :: N_modes_loc(3)
    real(kind=Rkind),  intent(in)    :: Eigenenergies_loc(:)
    logical, optional, intent(in)    :: Debug_opt

    integer                          :: i
    logical                          :: Debug_local = .FALSE.


    IF (PRESENT(Debug_opt)) THEN; Debug_local = Debug_opt
    ELSE; Debug_local = .FALSE.; END IF 


    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "--- Arguments given to Compute_Ana_Eigenenergie"
      WRITE(out_unit,*) "E_0, E_1, E_2, E_3 = "//TO_string(Eigenenergies_loc(1))//", "//TO_string(Eigenenergies_loc(2))//", "//TO&
      &_string(Eigenenergies_loc(3))//", "//TO_string(Eigenenergies_loc(4))
    END IF 

    DO i = 1, SIZE(N_modes_loc)
      N_modes_loc(i) = Eigenenergies_loc(i+1) - Eigenenergies_loc(1)
    END DO

    IF (Debug_local) WRITE(out_unit,*)
    IF (Debug_local) WRITE(out_unit,*) "Normal modes :"
    IF (Debug_local) CALL Write_Vec(N_modes_loc, out_unit, 1, info="Normal modes")

  END SUBROUTINE Compute_Normal_modes


END PROGRAM
