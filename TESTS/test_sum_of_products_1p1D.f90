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
  IMPLICIT NONE


  integer                 :: Verbose = 40
  logical                 :: Debug   = .FALSE.

  logical                 :: Dense   = .FALSE.
  TYPE(Sum_of_products_t) :: TotH
  TYPE(Sum_of_products_t) :: DipMomt

  real(kind=Rkind)        ::    Psi_R1_real(6)                                                                             ! an any vector representing the excitation state/wavefunction = a linear combination of the basis functions from the canonical basis set over \mathbb{R} /!\ Not normalized yet !
  real(kind=Rkind)        :: Op_psi_R1_real(6)                                                                                ! the resulting vector from the action of a 1D operator upon Psi_ND_R1_real
  complex(kind=Rkind)     ::    Psi_R1_complex(6)                                                                          ! an any vector representing the excitation state/wavefunction = a linear combination of the basis functions from the canonical basis set over \mathbb{R} BUT with complexes expansion coefficients to make complexe WF. /!\ Not normalized yet !
  complex(kind=Rkind)     :: Op_psi_R1_complex(6)                                                                             ! the resulting vector from the action of a 1D operator upon Psi_ND_R1_complex


!   !----------------------------------Wavefunction initialization---------------------------------
  Psi_R1_real = [SQRT(TWO), HALF, PI, ONETENTH, TWELVE, ONE]                                ! /!\ the indexes are here renamed to match the indexes of the basis vectors and coefficients ! the elements starts from 0 to 5 and not from 1 to 6 !!! /!\
  CALL Normalize(Psi_R1_real)
  Psi_R1_complex = [SQRT(TWO)*EYE, COMPLEX(HALF, Rkind), COMPLEX(PI, Rkind), COMPLEX(ONETENTH, Rkind)+HALF*EYE, TW&
  &ELVE*EYE, COMPLEX(ONE, Rkind)]
  CALL Normalize(Psi_R1_complex)


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

  
  !----------------------------Testing the actions---------------------------
  CALL Action(Op_psi_R1_real, TotH, Psi_R1_real, Verbose=Verbose, Debug=Debug) ! at the end we have the resulting action of the total 1p1D hamiltonian upon Psi_R1_real  
  WRITE(out_unit,*)
  WRITE(out_unit,*) "*** RESULTING WF VECTOR FROM ACTION OF TOTH ON PSI_R1_REAL"
  CALL Write_Vec(Op_psi_R1_real, out_unit, 1, info="Op_psi_R1_real")

  CALL Action(Op_psi_R1_complex, TotH, Psi_R1_complex, Verbose=Verbose, Debug=Debug)
  WRITE(out_unit,*)
  WRITE(out_unit,*) "*** RESULTING WF VECTOR FROM ACTION OF TOTH ON PSI_R1_COMPLEX"
  CALL Write_Vec(Op_psi_R1_complex, out_unit, 1, info="Op_psi_R1_complex")


  !----------------------------Testing the writing---------------------------
  CALL Write(TotH)

  
  !----------------------------Testing the deallocation---------------------------
  CALL Dealloc(TotH, Verbose=Verbose, Debug=Debug)
  CALL Write(TotH)


END PROGRAM
