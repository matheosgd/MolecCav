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
PROGRAM test_cavity_mode
  !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
  USE QDUtil_m
  USE Tests_m
  USE Algebra_m
  USE Cavity_mode_m
  IMPLICIT NONE


  integer             :: Verbose = 40
  logical             :: Debug   = .TRUE.

  TYPE(Cavity_mode_new_t) :: CavMode
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

  TYPE(test_t)        :: test_cavmode
  logical             :: error_cavmode = .FALSE.

  integer             :: i, i_op

  !############ OLD ####################
  TYPE(Cavity_mode_t) :: Cavity_mode_1
  logical             :: error = .FALSE.


  !-----------------------------Test initialization----------------------------
  CALL Initialize_Test(test_cavmode, test_name="OUT/test_file_cav_mode")


  !-------------------------Cavity mode initialization-------------------------
  CALL Initialize(CavMode, nio=in_unit, Dense=Dense, Verbose=Verbose, Debug=Debug)
  IF (Debug) THEN
    WRITE(out_unit,*)
    WRITE(out_unit,*) "--------------Cavity mode constructed by MolecCav_Initialize_matter_mode--------------"
    CALL Write(CavMode)
    WRITE(out_unit,*) "------------End Cavity mode constructed by MolecCav_Initialize_matter_mode------------"
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
  DO i = 0, SIZE(CavMode%Tab_op)-1
    CALL Action(Op_psi_real,    CavMode, i, Psi_1D_R1_real,    Verbose=Verbose, Debug=Debug)
    CALL Action(Op_psi_complex, CavMode, i, Psi_1D_R1_complex, Verbose=Verbose, Debug=Debug)    
  END DO


  !----------------------------Testing the actions---------------------------
  CALL Dealloc(CavMode, Verbose=Verbose, Debug=Debug)


  ! ################################### OLD #######################################
  !-------------------------Cavity mode initialization-------------------------
  CALL MolecCav_Read_cavity_mode(Mode=Cavity_mode_1, nio=in_unit)
  IF (Debug) THEN
    WRITE(out_unit,*)
    WRITE(out_unit,*) "--------------Cavity mode constructed by MolecCav_Read_cavity_mode--------------"
    CALL Write_cavity_mode(Cavity_mode_1)
    WRITE(out_unit,*) "------------End Cavity mode constructed by MolecCav_Read_cavity_mode------------"
  END IF


  !-----------------------------------The tests--------------------------------
  error = (Cavity_mode_1%D == 0)
  CALL Logical_Test(test_cavmode, error, test2=.FALSE., info="Cavity_mode_1%D = 0")
  IF (error .AND. Debug) WRITE(out_unit,*) 'Mode%D failed to initialize'

  error = (Cavity_mode_1%Nb == 0)
  CALL Logical_Test(test_cavmode, error, test2=.FALSE., info="Cavity_mode_1%Nb == 0")
  IF (error .AND. Debug) WRITE(out_unit,*) 'Mode%Nb failed to initialize'

  error = (Cavity_mode_1%w == 0)
  CALL Logical_Test(test_cavmode, error, test2=.FALSE., info="Cavity_mode_1%w == 0")
  IF (error .AND. Debug) WRITE(out_unit,*) 'Mode%w failed to initialize'

  error = (Cavity_mode_1%m == 0)
  CALL Logical_Test(test_cavmode, error, test2=.FALSE., info="Cavity_mode_1%m == 0")
  IF (error .AND. Debug) WRITE(out_unit,*) 'Mode%m failed to initialize'

  error = (Cavity_mode_1%lambda < 0)
  CALL Logical_Test(test_cavmode, error, test2=.FALSE., info="Cavity_mode_1%lambda == 0")
  IF (error .AND. Debug) WRITE(out_unit,*) 'Mode%lambda failed to initialize'

  error = (Cavity_mode_1%eq_pos < 0)
  CALL Logical_Test(test_cavmode, error, test2=.FALSE., info="Cavity_mode_1%eq_pos == 0")
  IF (error .AND. Debug) WRITE(out_unit,*) 'Mode%eq_pos failed to initialize'

  CALL Finalize_Test(test_cavmode)

  
END PROGRAM
