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
PROGRAM test_sum_of_products_1p1D
  !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
  USE QDUtil_m
  USE Algebra_m
  USE Sum_of_products_m
  USE Transition_spectrum_m
  IMPLICIT NONE


  integer                       :: Verbose = 40
  logical                       :: Debug   = .FALSE.

  logical                       :: Dense   = .FALSE.
  TYPE(Sum_of_products_t)       :: TotH
  TYPE(Sum_of_products_t)       :: DipMomt

  integer                       :: Nb_1, Nb_2, NB, J
  real(kind=Rkind), allocatable :: Phi(:)
  real(kind=Rkind), allocatable :: TotH_matrix(:,:)
  real(kind=Rkind), allocatable :: REigvec(:,:)
  real(kind=Rkind), allocatable :: REigval(:)

  real(kind=Rkind)              :: E_threshold
  integer                       :: N_states
  TYPE(Transition_spectrum_t)   :: TranSpec


  !-------------------------Sum_of_products operators initialization-------------------------
  CALL Initialize_totH(TotH, in_unit, Dense=Dense, Verbose=Verbose, Debug=.TRUE.)
  IF (Debug) THEN
    WRITE(out_unit,*)
    WRITE(out_unit,*) "--------------TotH object constructed by MolecCav_Initialize_total_hamiltonian--------------"
    CALL Write(TotH)
    WRITE(out_unit,*) "------------End TotH object constructed by MolecCav_Initialize_total_hamiltonian------------"
  END IF

  CALL Initialize_dipmomt(DipMomt, in_unit, Dense=Dense, Verbose=Verbose, Debug=.TRUE.)
  IF (Debug) THEN
    WRITE(out_unit,*)
    WRITE(out_unit,*) "--------------DipMomt object constructed by MolecCav_Initialize_dipmomt--------------"
    CALL Write(DipMomt)
    WRITE(out_unit,*) "------------End DipMomt object constructed by MolecCav_Initialize_dipmomt------------"
  END IF

  !-------------------------System initialization-------------------------
  CALL Get(Nb_1, "Nb", "Matter", 1)
  CALL Get(Nb_2, "Nb", "Cavity", 1)
  NB = Nb_1 * Nb_2
  ALLOCATE(TotH_matrix(NB, NB))
  ALLOCATE(Phi(NB))
  TotH_matrix = ZERO

  DO J = 1, NB
    Phi = ZERO
    Phi(J) = ONE
    CALL Action(TotH_matrix(:,J), TotH, Phi, Verbose=Verbose, Debug=.FALSE.)
  END DO

  WRITE(out_unit,*)
  WRITE(out_unit,*) "*** RESULTING MATRIX OF TOTH"
  CALL Write_Mat(TotH_matrix, out_unit, NB, info="TotH_matrix")

  ALLOCATE(REigval(NB))
  ALLOCATE(REigvec(NB,NB))

  CALL diagonalization(TotH_matrix, REigval, REigvec)
  WRITE(out_unit,*)
  CALL Write_Vec(REigval, out_unit, 1, info="EigenEnergies(1:12)")
  
  !----------------------------Testing the computation of the spectra---------------------------
  CALL Initialize(TranSpec, REigvec, DipMomt, E_threshold=E_threshold, REigval=REigval, Nb_states=N_states, Verbose=Verbose, Deb&
  &ug=Debug)

  !----------------------------Testing the writing---------------------------
  CALL Write(TotH)

  
  !----------------------------Testing the deallocation---------------------------
  CALL Dealloc(TotH, Verbose=Verbose, Debug=Debug)
  CALL Write(TotH)


END PROGRAM
