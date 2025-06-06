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
PROGRAM test_sum_of_products_2p1D
  !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
  USE QDUtil_m
  USE Sum_of_products_m
  IMPLICIT NONE


  integer                          :: Verbose = 40
  logical                          :: Debug   = .TRUE.

  logical                          :: Dense   = .FALSE.
  TYPE(Sum_of_products_t)          :: TotH
  TYPE(Sum_of_products_t)          :: DipMomt

  integer                          :: Mat1Nb
  integer                          :: Mat2Nb
  integer                          :: CavNb
  integer                          :: NB

  real(kind=Rkind),    allocatable ::    Phi(:)                                                                             ! an any vector representing the excitation state/wavefunction = a linear combination of the basis functions from the canonical basis set over \mathbb{R} /!\ Not normalized yet !
  complex(kind=Rkind), allocatable ::    Phi_complex(:)                                                                             ! an any vector representing the excitation state/wavefunction = a linear combination of the basis functions from the canonical basis set over \mathbb{R} /!\ Not normalized yet !
  complex(kind=Rkind), allocatable :: Op_phi_complex(:)                                                                             ! an any vector representing the excitation state/wavefunction = a linear combination of the basis functions from the canonical basis set over \mathbb{R} /!\ Not normalized yet !

  real(kind=Rkind), allocatable    :: TotH_matrix(:,:)
  real(kind=Rkind), allocatable    :: Eigenenergies(:)                                                                             ! an any vector representing the excitation state/wavefunction = a linear combination of the basis functions from the canonical basis set over \mathbb{R} /!\ Not normalized yet !
  real(kind=Rkind), allocatable    :: Eigenstates(:,:)                                                                             ! an any vector representing the excitation state/wavefunction = a linear combination of the basis functions from the canonical basis set over \mathbb{R} /!\ Not normalized yet !

  integer                          :: J, i_1, i_2, i_3, min_index


  !-------------------------Sum_of_products operators initialization-------------------------
  CALL Initialize_totH(TotH, in_unit, Dense=Dense, Verbose=Verbose, Debug=Debug)
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

  
  !----------------------------------Wavefunction initialization---------------------------------
  CALL Get(Mat1Nb, "Nb", "Matter", 1)
  CALL Get(Mat2Nb, "Nb", "Matter", 2)
  CALL Get(CavNb,  "Nb", "Cavity", 1)
  NB = Mat1Nb * Mat2Nb * CavNb

  ALLOCATE(Phi(NB))
  ALLOCATE(TotH_matrix(NB,NB))
  ALLOCATE(Phi_complex(NB))
  ALLOCATE(Op_phi_complex(NB))

  Phi_complex = ZERO
  DO J = 1, SIZE(Phi_complex)
    Phi_complex(J) = J*HALF + J*EYE  
  END DO


  !----------------------------Testing the actions---------------------------
  TotH_matrix = ZERO

  DO J = 1, NB
    Phi = ZERO
    Phi(J) = ONE
    CALL Action(TotH_matrix(:,J), TotH, Phi, Verbose=Verbose, Debug=.FALSE.)
  END DO

  WRITE(out_unit,*)
  WRITE(out_unit,*) "*** RESULTING MATRIX OF TOTH"
  CALL Write_Mat(TotH_matrix, out_unit, NB, info="TotH_matrix")

  ALLOCATE(Eigenenergies(NB))
  ALLOCATE(Eigenstates(NB,NB))

  CALL diagonalization(TotH_matrix, Eigenenergies, Eigenstates)
  WRITE(out_unit,*)
  CALL Write_Vec(Eigenenergies(1:12), out_unit, 1, info="EigenEnergies(1:12)")

  CALL Action(Op_phi_complex, TotH, Phi_complex, Verbose=Verbose, Debug=Debug)
  WRITE(out_unit,*)
  WRITE(out_unit,*) "*** RESULTING WF VECTOR FROM ACTION OF TOTH ON PHI COMPLEX"
  CALL Write_Vec(Op_phi_complex, out_unit, 1, info="TotH_Phi_complex")


  !----------------------------Testing the writing---------------------------
  CALL Write(TotH)

  
  !----------------------------Testing the deallocation---------------------------
  CALL Dealloc(TotH, Verbose=Verbose, Debug=Debug)


END PROGRAM
