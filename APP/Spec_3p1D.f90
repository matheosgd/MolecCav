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
PROGRAM Spec_3p1D
  !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
  USE QDUtil_m
  USE Algebra_m
  USE Sum_of_products_m
  USE Transition_spectrum_m
  IMPLICIT NONE


  integer                       :: Verbose = 40
  logical                       :: Debug   = .FALSE.
  integer                       :: nioint, niospec

  logical                       :: Dense   = .FALSE.
  TYPE(Sum_of_products_t)       :: TotH
  TYPE(Sum_of_products_t)       :: DipMomt

  real(kind=Rkind)              :: Matw, Cavw, Matlambda, Cavlambda, lambda
  integer                       :: Nb_1, Nb_2, Nb_3, Nb_4, NB, J
  real(kind=Rkind), allocatable :: Phi(:)
  real(kind=Rkind), allocatable :: TotH_matrix(:,:)
  real(kind=Rkind), allocatable :: REigvec(:,:)
  real(kind=Rkind), allocatable :: REigval(:)

  real(kind=Rkind)              :: E_threshold
  integer                       :: N_states
  TYPE(Transition_spectrum_t)   :: TranSpec

  real(kind=Rkind)              :: Gamma, Start_plot, Stop_plot, Step_plot, Conversion, Energy, Intensity
  integer                       :: I


  !-------------------------Sum_of_products operators initialization-------------------------
  CALL Initialize_totH(TotH, in_unit, Dense=Dense, Verbose=Verbose, Debug=Debug)
  IF (Debug) THEN
    WRITE(out_unit,*)
    WRITE(out_unit,*) "--------------TotH object constructed by MolecCav_Initialize_total_hamiltonian--------------"
    CALL Write(TotH)
    WRITE(out_unit,*) "------------End TotH object constructed by MolecCav_Initialize_total_hamiltonian------------"
  END IF

  CALL Initialize_dipmomt(DipMomt, in_unit, Dense=Dense, Verbose=Verbose, Debug=Debug)
  IF (Debug) THEN
    WRITE(out_unit,*)
    WRITE(out_unit,*) "--------------DipMomt object constructed by MolecCav_Initialize_dipmomt--------------"
    CALL Write(DipMomt)
    WRITE(out_unit,*) "------------End DipMomt object constructed by MolecCav_Initialize_dipmomt------------"
  END IF

  !-------------------------System initialization-------------------------
  CALL Get(Matw, "w", "Matter", 1) ! all matter are the same so far
  CALL Get(Cavw, "w", "Cavity", 1)
  CALL Get(Matlambda, "lambda", "Matter", 1)
  CALL Get(Cavlambda, "lambda", "Cavity", 1)
  lambda = Matlambda * Cavlambda ! used only for naming the files and as indication (not for the caculations)
  CALL Get(Nb_1, "Nb", "Matter", 1)
  CALL Get(Nb_2, "Nb", "Matter", 2)
  CALL Get(Nb_3, "Nb", "Matter", 3)
  CALL Get(Nb_4, "Nb", "Cavity", 1)
  NB = Nb_1 * Nb_2 * Nb_3 * Nb_4

  WRITE(out_unit,*)
  WRITE(out_unit,*) "--- System parameters :"
  WRITE(out_unit,*) "Matw      = "//TO_string(Matw)
  WRITE(out_unit,*) "Cavw      = "//TO_string(Cavw)
  WRITE(out_unit,*) "Matlambda = "//TO_string(Matlambda)
  WRITE(out_unit,*) "Cavlambda = "//TO_string(Cavlambda)
  WRITE(out_unit,*) "lambda    = "//TO_string(lambda)
  WRITE(out_unit,*) "Mat1Nb    = "//TO_string(Nb_1)
  WRITE(out_unit,*) "Mat2Nb    = "//TO_string(Nb_2)
  WRITE(out_unit,*) "Mat3Nb    = "//TO_string(Nb_3)
  WRITE(out_unit,*) "CavNb     = "//TO_string(Nb_4)

  ALLOCATE(TotH_matrix(NB, NB))
  ALLOCATE(Phi(NB))
  TotH_matrix = ZERO

  DO J = 1, NB
    Phi = ZERO
    Phi(J) = ONE
    CALL Action(TotH_matrix(:,J), TotH, Phi, Verbose=Verbose, Debug=Debug)
  END DO

  WRITE(out_unit,*)
  WRITE(out_unit,*) "*** RESULTING MATRIX OF TOTH"
  !CALL Write_Mat(TotH_matrix, out_unit, NB, info="TotH_matrix")

  ALLOCATE(REigval(NB))
  ALLOCATE(REigvec(NB,NB))

  CALL diagonalization(TotH_matrix, REigval, REigvec)
  WRITE(out_unit,*)
  !CALL Write_Vec(REigval, out_unit, 1, info="EigenEnergies")
  

  !----------------------------Computing of the spectra---------------------------
  ! CALL Initialize(TranSpec, REigvec, DipMomt, E_threshold=E_threshold, REigval=REigval, Nb_states=N_states, Verbose=Verbose, Deb&
  ! &ug=Debug)
  ! CALL Initialize(TranSpec, REigvec, DipMomt, E_threshold=E_threshold, REigval=REigval, Nb_states=N_states, Verbose=Verbose, Deb&
  ! &ug=Debug)
  CALL Initialize(TranSpec, REigvec, DipMomt, REigval, Nb_states=10, Verbose=Verbose, Debug=Debug)

  WRITE(out_unit,*) "Spectrum information"
  CALL Write_Vec(TranSpec%tab_energies, out_unit, SIZE(TranSpec%tab_energies), info="Transition Energies")
  CALL Write_Vec(TranSpec%tab_ints,     out_unit, SIZE(TranSpec%tab_ints),     info="Transition Intensities")


  OPEN(NEWUNIT = nioint, FILE = 'OUT/TransInts_3p1D_wmat'//TO_string(REAL(Matw,kind=RkS))//'_wcav'//TO_string(REAL(Cavw,kind=RkS)&
  &)//'_lamb'//TO_string(REAL(lambda,kind=RkS))//'.txt',  FORM = 'formatted', ACTION = 'write', POSITION = 'rewind')

  OPEN(NEWUNIT = niospec, FILE = 'OUT/Spectrum_3p1D_wmat'//TO_string(REAL(Matw,kind=RkS))//'_wcav'//TO_string(REAL(Cavw,kind=RkS)&
  &)//'_lamb'//TO_string(REAL(lambda,kind=RkS))//'.txt',  FORM = 'formatted', ACTION = 'write', POSITION = 'rewind')

  !----------------------------computing of the spectra---------------------------
  WRITE(nioint, *) "Transition Energy -------- Transition intensity"
  DO J = 1, TranSpec%N_trstns
    IF (TranSpec%tab_ints(J)>1E-10) WRITE(nioint, *) TranSpec%tab_energies(J), TranSpec%tab_ints(J)
    IF (TranSpec%tab_ints(J)<1E-10) WRITE(nioint, *) TranSpec%tab_energies(J), 0
  END DO


  Conversion = 219474.6                       ! 1 Ha = Conversion.cm-1
  Gamma      = 1.0                            ! => 10cm-1
  Start_plot = 5.8154906356422710E-003 - 1E-4 ! in Ha 
  Stop_plot  = 5.9170696435614355E-003 + 1E-4 ! in Ha
  Step_plot  = (Stop_plot - Start_plot) / 900 ! divide by desired number of points

  WRITE(niospec, *) "Energy -------- Transition intensity"
  DO I = 1, 900
    Energy    = (Start_plot + I*Step_plot)*Conversion
    Intensity = TranSpec%tab_ints(1)*Lorentzian(Energy, TranSpec%tab_energies(1)*Conversion, Gamma) &
    &         + TranSpec%tab_ints(4)*Lorentzian(Energy, TranSpec%tab_energies(4)*Conversion, Gamma)
    WRITE(niospec, *) Energy, Intensity
  END DO


END PROGRAM
