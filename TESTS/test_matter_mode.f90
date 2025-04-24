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
PROGRAM Quantum_HO1D
  !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
  USE QDUtil_m
  USE Tests_m
  USE Algebra_m
  USE Matter_mode_m
  IMPLICIT NONE


  integer             :: Verbose = 40
  logical             :: Debug   = .TRUE.

  TYPE(Matter_mode_t) :: MatMode
  logical             :: Dense   = .FALSE.

  real(kind=Rkind)    :: Psi_1D_R1_real(3)                                                                             ! an any vector representing the excitation state/wavefunction = a linear combination of the basis functions from the canonical basis set over \mathbb{R} /!\ Not normalized yet !
  real(kind=Rkind)    :: Coeff_0_real = ONE
  real(kind=Rkind)    :: Coeff_1_real = HALF
  real(kind=Rkind)    :: Coeff_2_real = PI
  real(kind=Rkind)    :: Op_psi_real(3)                                                                                ! the resulting vector from the action of a 1D operator upon Psi_1D_R1_real
  complex(kind=Rkind) :: Psi_1D_R1_complex(3)                                                                          ! an any vector representing the excitation state/wavefunction = a linear combination of the basis functions from the canonical basis set over \mathbb{R} BUT with complexes expansion coefficients to make complexe WF. /!\ Not normalized yet !
  complex(kind=Rkind) :: Coeff_0_complex = ONE*EYE + SQRT(TWO)
  complex(kind=Rkind) :: Coeff_1_complex = HALF
  complex(kind=Rkind) :: Coeff_2_complex = PI*EYE
  complex(kind=Rkind) :: Op_psi_complex(3)                                                                             ! the resulting vector from the action of a 1D operator upon Psi_1D_R1_complex

  TYPE(test_t)        :: test_matmode
  logical             :: error_matmode = .FALSE.

  integer             :: i, i_op


  !-----------------------------Test initialization----------------------------
  CALL Initialize_Test(test_matmode, test_name="OUT/test_file_mat_mode")


  !-------------------------Matter mode initialization-------------------------
  CALL Initialize(MatMode, nio=in_unit, Dense=Dense, Verbose=Verbose, Debug=Debug)
  IF (Debug) THEN
    WRITE(out_unit,*)
    WRITE(out_unit,*) "--------------Matter mode constructed by MolecCav_Initialize_matter_mode--------------"
    CALL Write(MatMode)
    WRITE(out_unit,*) "------------End Matter mode constructed by MolecCav_Initialize_matter_mode------------"
  END IF


  !-------------------------Wavefunction initialization (real)------------------------
  Psi_1D_R1_real(:) = [Coeff_0_real, Coeff_1_real, Coeff_2_real]
  CALL Normalize(Psi_1D_R1_real)
  IF (Debug) THEN
    WRITE(out_unit,*)
    WRITE(out_unit,*) "----------------The matter wavefunction has been initialized as :---------------"
    CALL Write_Vec(Psi_1D_R1_real, out_unit, Size(Psi_1D_R1_real), info="Psi_1D_R1_real")
    WRITE(out_unit,*) "-----------------------------End matter wavefunction----------------------------"
  END IF


  !-------------------------Wavefunction initialization (complex)------------------------
  Psi_1D_R1_complex(:) = [Coeff_0_complex, Coeff_1_complex, Coeff_2_complex]                                                     ! uses the same basis set as the real WF, but the expansion coefficients are complexes
  CALL Normalize(Psi_1D_R1_complex)
  IF (Debug) THEN
    WRITE(out_unit,*)
    WRITE(out_unit,*) "----------------The matter wavefunction has been initialized as :---------------"
    CALL Write_Vec(Psi_1D_R1_complex, out_unit, Size(Psi_1D_R1_complex), info="Psi_1D_R1_complex")
    WRITE(out_unit,*) "-----------------------------End matter wavefunction----------------------------"
  END IF


  !----------------------------Testing the actions---------------------------
  DO i = 0, SIZE(MatMode%Tab_op)-1
    CALL Action(Op_psi_real,    MatMode, i, Psi_1D_R1_real,    Verbose=Verbose, Debug=Debug)
    CALL Action(Op_psi_complex, MatMode, i, Psi_1D_R1_complex, Verbose=Verbose, Debug=Debug)    
  END DO


  !----------------------------Testing the actions---------------------------
  CALL Dealloc(MatMode, Verbose=Verbose, Debug=Debug)

  
  !-----------------------------------The tests--------------------------------

  CALL Finalize_Test(test_matmode)

  
END PROGRAM
