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
PROGRAM test_mapping
  !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
  USE QDUtil_m
  USE QDUtil_Test_m
  USE Algebra_m
  USE Mapping_m
  USE Cavity_mode_m
  USE Operator_1D_m
  USE Total_hamiltonian_m
  IMPLICIT NONE


  logical, parameter            :: Debug = .TRUE.
  logical, parameter            :: Ease_visual_check = .FALSE.                                     ! takes (3,2) basis set instead of the data file one

  TYPE(Cavity_mode_t)           :: DOF_1                                                           ! DOF = the only Degree Of Freedom of the matter part of the system consider so far
  TYPE(Cavity_mode_t)           :: DOF_2                                                           ! The well construction of the Operator's matrices is assumed checked by the dedicated test, so hard-coded references matrices are not needed here

  real(kind=Rkind), allocatable :: b_00_2D(:,:)                                                    ! six vectors of the HO basis set |00>, |10>, |20>, |01>, |11>, |21> 
  real(kind=Rkind), allocatable :: b_01_2D(:,:)
  real(kind=Rkind), allocatable :: b_10_2D(:,:)
  real(kind=Rkind), allocatable :: b_11_2D(:,:)
  real(kind=Rkind), allocatable :: b_20_2D(:,:)
  real(kind=Rkind), allocatable :: b_21_2D(:,:)
  real(kind=Rkind)              :: Coeff_00 = ONE
  real(kind=Rkind)              :: Coeff_01 = HALF
  real(kind=Rkind)              :: Coeff_10 = THREE
  real(kind=Rkind)              :: Coeff_11 = TEN
  real(kind=Rkind)              :: Coeff_20 = SEVEN
  real(kind=Rkind)              :: Coeff_21 = ONETENTH
  real(kind=Rkind)              :: Coeffs(0:5)                                                     ! /!\ the indexes are here renamed to match the indexes of the basis vectors and coefficients ! the elements starts from 0 to 5 and not from 1 to 6 !!! /!\
  real(kind=Rkind), allocatable :: Psi_2D(:,:)                                                     ! an any vector representing the excitation state/wavefunction of the HO = a linear combination of the basis functions /!\ Not normalized yet !

  real(kind=Rkind), allocatable :: b_00_2D_bis(:,:)                                                    ! six vectors of the HO basis set |00>, |10>, |20>, |01>, |11>, |21> 
  real(kind=Rkind), allocatable :: b_01_2D_bis(:,:)
  real(kind=Rkind), allocatable :: b_10_2D_bis(:,:)
  real(kind=Rkind), allocatable :: b_11_2D_bis(:,:)
  real(kind=Rkind), allocatable :: b_20_2D_bis(:,:)
  real(kind=Rkind), allocatable :: b_21_2D_bis(:,:)
  real(kind=Rkind), allocatable :: Psi_2D_bis(:,:)                                                     ! an any vector representing the excitation state/wavefunction of the HO = a linear combination of the basis functions /!\ Not normalized yet !

  real(kind=Rkind), allocatable :: b_00_1D(:)                                                      ! six vectors of the HO basis set |00>, |10>, |20>, |01>, |11>, |21> 
  real(kind=Rkind), allocatable :: b_01_1D(:)
  real(kind=Rkind), allocatable :: b_10_1D(:)
  real(kind=Rkind), allocatable :: b_11_1D(:)
  real(kind=Rkind), allocatable :: b_20_1D(:)
  real(kind=Rkind), allocatable :: b_21_1D(:)
  real(kind=Rkind), allocatable :: Psi_1D(:)

  ! for a particular test we will assume DOF_1 = Matter_DOF, DOF_2 = Cavity_mode, to be able to compute the action of a total Hamiltonian 1p1D
  TYPE(Operator_1D_t)           :: CavPosition                                                     ! position operator of the cavity mode
  TYPE(Operator_1D_t)           :: CavH                                                            ! Hamiltonian operator of the cavity mode
  TYPE(Operator_1D_t)           :: Mat_dipolar_moment                                              ! dipolar moment operator of the matter subsystem
  real(kind=Rkind)              :: Coeff_dipole_moment = ONE                                       ! variation coefficient of the dipole moment with the matter DOF
  TYPE(Operator_1D_t)           :: MatH                                                            ! Hamiltonian operator of the matter subsystem  
  real(kind=Rkind), allocatable :: TotH(:,:)                                                       ! size : NBxNB; NB = Nb_M*Nb_C (tensor product of the basis sets)
  real(kind=Rkind), allocatable :: TotH_psi_1p1D                                                   ! the resulting vector from the action of a 1D operator upon Psi_1D

  integer                       :: Nb_1, Nb_2, NB

  TYPE(test_t)                  :: test_mapp
  logical                       :: error_mapp = .FALSE.
    
  
  !---------------------------------------Test initialization--------------------------------------
  CALL Initialize_Test(test_mapp, test_name="OUT/test_file_mapp")
  
  
  !--------------------------------------System initialization-------------------------------------
    !-------------------------------------Modes initialization-------------------------------------
  CALL Read_cavity_mode(DOF_1, nio=in_unit)
  IF (Debug) THEN
    WRITE(out_unit,*)
    WRITE(out_unit,*); WRITE(out_unit,*) "--------------Cavity mode constructed by MolecCav_Read_cavity_mode--------------"
    CALL Write_cavity_mode(DOF_1)
    WRITE(out_unit,*) "------------End Cavity mode constructed by MolecCav_Read_cavity_mode------------"
  END IF
  IF (Ease_visual_check) DOF_1%Nb = 3

  CALL Read_cavity_mode(DOF_2, nio=in_unit)
  IF (Debug) THEN
    WRITE(out_unit,*)
    WRITE(out_unit,*); WRITE(out_unit,*) "--------------Cavity mode constructed by MolecCav_Read_cavity_mode--------------"
    CALL Write_cavity_mode(DOF_2)
    WRITE(out_unit,*) "------------End Cavity mode constructed by MolecCav_Read_cavity_mode------------"
  END IF
  IF (Ease_visual_check) DOF_2%Nb = 2

  Nb_1 = DOF_1%Nb
  Nb_2 = DOF_2%Nb
  NB   = Nb_1*Nb_2

    !----------------------------------Wavefunction initialization---------------------------------
  ALLOCATE(b_00_2D(Nb_1, Nb_2))
  ALLOCATE(b_01_2D(Nb_1, Nb_2))
  ALLOCATE(b_10_2D(Nb_1, Nb_2))
  ALLOCATE(b_11_2D(Nb_1, Nb_2))
  ALLOCATE(b_20_2D(Nb_1, Nb_2))
  ALLOCATE(b_21_2D(Nb_1, Nb_2))
  ALLOCATE(Psi_2D(Nb_1,  Nb_2))

  b_00_2D      = ZERO
  b_00_2D(1,1) = ONE
  b_01_2D      = ZERO
  b_01_2D(1,2) = ONE
  b_10_2D      = ZERO
  b_10_2D(2,1) = ONE
  b_11_2D      = ZERO
  b_11_2D(2,2) = ONE
  b_20_2D      = ZERO
  b_20_2D(3,1) = ONE
  b_21_2D      = ZERO
  b_21_2D(3,2) = ONE
  
  Coeffs = [Coeff_00, Coeff_01, Coeff_10, Coeff_11, Coeff_20, Coeff_21]                            ! /!\ the indexes are here renamed to match the indexes of the basis vectors and coefficients ! the elements starts from 0 to 5 and not from 1 to 6 !!! /!\

  Psi_2D = Coeffs(0)*b_00_2D + Coeffs(1)*b_01_2D + Coeffs(2)*b_10_2D + Coeffs(3)*b_11_2D + Coeffs(4)*b_20_2D + Coeffs(5)*b_21_2D
  
  IF (Debug) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*); WRITE(out_unit,*) "---------------Basis vectors of the [DOF_1 x DOF_2] 2D system---------------"
      CALL Write_Mat(b_00_2D, out_unit, Size(b_00_2D, dim=2), info="b_00_2D"); WRITE(out_unit,*)
      CALL Write_Mat(b_01_2D, out_unit, Size(b_00_2D, dim=2), info="b_01_2D"); WRITE(out_unit,*)
      CALL Write_Mat(b_10_2D, out_unit, Size(b_00_2D, dim=2), info="b_10_2D"); WRITE(out_unit,*)
      CALL Write_Mat(b_11_2D, out_unit, Size(b_00_2D, dim=2), info="b_11_2D"); WRITE(out_unit,*)
      CALL Write_Mat(b_20_2D, out_unit, Size(b_00_2D, dim=2), info="b_20_2D"); WRITE(out_unit,*)
      CALL Write_Mat(b_21_2D, out_unit, Size(b_00_2D, dim=2), info="b_21_2D"); WRITE(out_unit,*)
      WRITE(out_unit,*); WRITE(out_unit,*) "-------------End basis vectors of the [DOF_1 x DOF_2] 2D system-------------"
      WRITE(out_unit,*); WRITE(out_unit,*) "----------------The any linear combination of the basis functions---------------"
      CALL Write_Mat(Psi_2D, out_unit, Size(Psi_2D, dim=2), info="Psi_2D(NOT normalized)")
      WRITE(out_unit,*); WRITE(out_unit,*) "------------------------End of the any linear combination-----------------------"
  END IF

  CALL Normalize(Psi_2D)
  IF (Debug) THEN
    WRITE(out_unit,*)
    WRITE(out_unit,*); WRITE(out_unit,*) "----------------The any linear combination of the basis functions---------------"
    CALL Write_Mat(Psi_2D, out_unit, Size(Psi_2D, dim=2), info="Psi_2D(normalized)")
    WRITE(out_unit,*); WRITE(out_unit,*) "------------------------End of the any linear combination-----------------------"
  END IF

    !-----------------------------------operators initialization-----------------------------------
  CALL Construct_Operator_1D(Mat_dipolar_moment, "Position",    Mode=DOF_1,  Debug=Debug)          ! hypothesis : \mu(q) = \frac{\partial\mu}{\partial q}*\q ; \frac{\partial\mu}{\partial q} = Coeff_dip_momt
  Mat_dipolar_moment%Band_val_R = Coeff_dipole_moment*Mat_dipolar_moment%Band_val_R                ! => \hat{\mu} = Coeff_dip_momt*\hat{q}
  CALL Construct_Operator_1D(MatH,               "Hamiltonian", Mode=DOF_1,  Debug=Debug)
  CALL Construct_Operator_1D(CavPosition,        "Position",    Mode=DOF_2, Debug=Debug)
  CALL Construct_Operator_1D(CavH,               "Hamiltonian", Mode=DOF_2, Debug=Debug)
  ALLOCATE(TotH(NB, NB))
  CALL Construct_total_hamiltonian_1p1D(TotH, CavPosition, CavH, Mat_dipolar_moment, MatH, Debug=.FALSE.)


  !------------------------------------------Tests mapping-----------------------------------------
  WRITE(out_unit,*)
  WRITE(out_unit,*) "------------------------------------------------b_00------------------------------------------------"
  ALLOCATE(b_00_1D(NB))
  CALL Mapping_WF_2DTO1D(b_00_1D, b_00_2D, Debug=Debug)
  ALLOCATE(b_00_2D_bis(Nb_1, Nb_2))
  CALL Mapping_WF_1DTO2D(b_00_2D_bis, b_00_1D, Debug=Debug)
  CALL Equal_R_R_matrix(error_mapp, b_00_2D_bis, b_00_2D)
  CALL Logical_Test(test_mapp, error_mapp, test2=.FALSE., info="b_00_2D recovered ?")

  WRITE(out_unit,*)
  WRITE(out_unit,*) "------------------------------------------------b_01------------------------------------------------"
  ALLOCATE(b_01_1D(NB))
  CALL Mapping_WF_2DTO1D(b_01_1D, b_01_2D, Debug=Debug)
  ALLOCATE(b_01_2D_bis(Nb_1, Nb_2))
  CALL Mapping_WF_1DTO2D(b_01_2D_bis, b_01_1D, Debug=Debug)
  CALL Equal_R_R_matrix(error_mapp, b_01_2D_bis, b_01_2D)
  CALL Logical_Test(test_mapp, error_mapp, test2=.FALSE., info="b_01_2D recovered ?")

  WRITE(out_unit,*)
  WRITE(out_unit,*) "------------------------------------------------b_10------------------------------------------------"
  ALLOCATE(b_10_1D(NB))
  CALL Mapping_WF_2DTO1D(b_10_1D, b_10_2D, Debug=Debug)
  ALLOCATE(b_10_2D_bis(Nb_1, Nb_2))
  CALL Mapping_WF_1DTO2D(b_10_2D_bis, b_10_1D, Debug=Debug)
  CALL Equal_R_R_matrix(error_mapp, b_10_2D_bis, b_10_2D)
  CALL Logical_Test(test_mapp, error_mapp, test2=.FALSE., info="b_10_2D recovered ?")

  WRITE(out_unit,*)
  WRITE(out_unit,*) "------------------------------------------------b_11------------------------------------------------"
  ALLOCATE(b_11_1D(NB))
  CALL Mapping_WF_2DTO1D(b_11_1D, b_11_2D, Debug=Debug)
  ALLOCATE(b_11_2D_bis(Nb_1, Nb_2))
  CALL Mapping_WF_1DTO2D(b_11_2D_bis, b_11_1D, Debug=Debug)
  CALL Equal_R_R_matrix(error_mapp, b_11_2D_bis, b_11_2D)
  CALL Logical_Test(test_mapp, error_mapp, test2=.FALSE., info="b_11_2D recovered ?")

  WRITE(out_unit,*)
  WRITE(out_unit,*) "------------------------------------------------b_20------------------------------------------------"
  ALLOCATE(b_20_1D(NB))
  CALL Mapping_WF_2DTO1D(b_20_1D, b_20_2D, Debug=Debug)
  ALLOCATE(b_20_2D_bis(Nb_1, Nb_2))
  CALL Mapping_WF_1DTO2D(b_20_2D_bis, b_20_1D, Debug=Debug)
  CALL Equal_R_R_matrix(error_mapp, b_20_2D_bis, b_20_2D)
  CALL Logical_Test(test_mapp, error_mapp, test2=.FALSE., info="b_20_2D recovered ?")

  WRITE(out_unit,*)
  WRITE(out_unit,*) "------------------------------------------------b_21------------------------------------------------"
  ALLOCATE(b_21_1D(NB))
  CALL Mapping_WF_2DTO1D(b_21_1D, b_21_2D, Debug=Debug)
  ALLOCATE(b_21_2D_bis(Nb_1, Nb_2))
  CALL Mapping_WF_1DTO2D(b_21_2D_bis, b_21_1D, Debug=Debug)
  CALL Equal_R_R_matrix(error_mapp, b_21_2D_bis, b_21_2D)
  CALL Logical_Test(test_mapp, error_mapp, test2=.FALSE., info="b_21_2D recovered ?")

  WRITE(out_unit,*)
  WRITE(out_unit,*) "-------------------------------------------------Psi------------------------------------------------"
  ALLOCATE(Psi_1D(NB))
  CALL Mapping_WF_2DTO1D(Psi_1D, Psi_2D, Debug=Debug)
  ALLOCATE(Psi_2D_bis(Nb_1, Nb_2 ))
  CALL Mapping_WF_1DTO2D(Psi_2D_bis, Psi_1D, Debug=Debug)
  CALL Equal_R_R_matrix(error_mapp, Psi_2D_bis, Psi_2D)
  CALL Logical_Test(test_mapp, error_mapp, test2=.FALSE., info="Psi_2D recovered ?")

  CALL Finalize_Test(test_mapp)

  
  CONTAINS


  SUBROUTINE Equal_R_R_matrix(error, Rl_1, Rl_2)
    USE QDUtil_m
    IMPLICIT NONE 

    logical,          intent(inout) :: error
    real(kind=Rkind), intent(in)    :: Rl_1(:,:)
    real(kind=Rkind), intent(in)    :: Rl_2(:,:)
    
    real(kind=Rkind), parameter     :: Threshold   = 1E-10_Rkind
    logical, parameter              :: Debug_local = .FALSE.
    integer                         :: Nb_1_local, Nb_2_local

    Nb_1_local = Size(Rl_1, dim=1)
    Nb_2_local = Size(Rl_2, dim=2)
    IF (Nb_1_local /= Size(Rl_2, dim=1) .OR. Nb_2_local /= Size(Rl_2, dim=2)) THEN
      WRITE(out_unit,*) "The two matrices must have same dimensions to compare them. Please, check initialization."
      STOP "### The two matrices must have same dimensions to compare them. Please, check initialization."
    END IF 

    IF (ANY(ABS(Rl_1 - Rl_2) > Threshold)) THEN
      error = .TRUE.
      IF (Debug_local) THEN
        WRITE(out_unit,*) "The two matrices are not close enough to be considered equal :"
        CALL Write_Mat(Rl_1, out_unit, Nb_2_local, info="R_1(:,:)")
        CALL Write_Mat(Rl_2, out_unit, Nb_2_local, info="R_2(:,:)")
        CALL Write_Mat(ABS(Rl_1 - Rl_2), out_unit, Nb_2_local, info="|R_1-R_2| = ")
      END IF 

    ELSE 
      error = .FALSE.
      IF (Debug_local) THEN
        WRITE(out_unit,*) "The two matrices are close enough to be considered equal :"
        CALL Write_Mat(Rl_1, out_unit, Nb_2_local, info="R_1(:,:)")
        CALL Write_Mat(Rl_2, out_unit, Nb_2_local, info="R_2(:,:)")
        CALL Write_Mat(ABS(Rl_1 - Rl_2), out_unit, Nb_2_local, info="|R_1-R_2| = ")
      END IF
    END IF 

  END SUBROUTINE Equal_R_R_matrix


END PROGRAM
