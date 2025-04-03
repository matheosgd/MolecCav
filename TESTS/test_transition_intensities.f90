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
PROGRAM test_transition_intensities
  !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
  USE QDUtil_m
  USE QDUtil_Test_m
  USE Algebra_m
  USE Cavity_mode_m
  USE Elem_op_m
  USE Total_hamiltonian_m
  IMPLICIT NONE
  

  logical, parameter            :: Debug =     .TRUE.
  logical, parameter            :: View_more = .FALSE.

  !--------------------------------------Diatomic molecule in a HARMONIC electonic potential-------------------------------------
  TYPE(Cavity_mode_t)           :: Molecule_1
  TYPE(Operator_1D_t)           :: Mol1H                                                                                         ! matrix of the one-dimensional harmonic Hamiltonian associated with HO D
  TYPE(Operator_1D_t)           :: Mol1DipMomt
  real(kind=Rkind)              :: CteDipMomt = ONE                                                                              ! the intensity of the variation of the dipole moment with a variation of the matter DOF
  
  !-------------------------------------------------------First cavity mode------------------------------------------------------
  TYPE(Cavity_mode_t)           :: Cavity_mode_1
  TYPE(Operator_1D_t)           :: Cav1H                                                                                         ! matrix of the one-dimensional harmonic Hamiltonian associated with HO D
  TYPE(Operator_1D_t)           :: Cav1Position
  
  !-------------------------------------------------------Total Hamiltonian------------------------------------------------------
  real(kind=Rkind), allocatable :: TotH(:,:)

  !------------------------------------------------------------Results-----------------------------------------------------------
  real(kind=Rkind), allocatable :: REigval(:)
  real(kind=Rkind), allocatable :: REigvec(:,:)
  real(kind=Rkind), allocatable :: Intensities1(:,:)
  real(kind=Rkind), allocatable :: Intensities2(:,:)

  !-----------------------------------------------------------Utilities----------------------------------------------------------
  integer                       :: Nb_M, Nb_C, NB, I, N_min

  TYPE(test_t)                  :: test_trnstn_int
  logical                       :: error_trnstn_int = .FALSE.


  !------------------------------------------------------Test initialization-----------------------------------------------------
  CALL Initialize_Test(test_trnstn_int, test_name="OUT/test_file_trnstn_int")
  

  !-----------------------------------------------------SYSTEM INITIALIZATION----------------------------------------------------
    !-------------------------------------Diatomic molecule in a harmonic electonic potential------------------------------------
  WRITE(out_unit,*) "-------------------------------------------------------SYSTEM INITIALIZATION--------------------------------&
                    &----------------------"
  WRITE(out_unit,*) "  ---------------------------------------Diatomic molecule in a harmonic electonic potential----------------&
                    &----------------------"
  CALL MolecCav_Read_cavity_mode(Mode=Molecule_1, nio=in_unit)
  Molecule_1%Nb = 10
  Molecule_1%w  = 0.0058665_Rkind
  Molecule_1%m  = 1744.60504565_Rkind

  WRITE(out_unit,*) "Molecular Hamiltonian :"
  CALL Construct_Operator_1D(Operator=Mol1H, &
                           & operator_type="Hamiltonian", &
                           & Mode=Molecule_1, &
                           & Debug=.TRUE.)
  
  WRITE(out_unit,*) "Molecular Dipole moment :"
  CALL Construct_Operator_1D(Operator=Mol1DipMomt, &
                           & operator_type="Position", &                                                                         ! initialized as a position operator because of approximation over its expression (cf. readme.md or manual)
                           & Mode=Molecule_1, &
                           & Debug=.FALSE.)

  Mol1DipMomt%Band_val = Mol1DipMomt%Band_val*CteDipMomt                                                              ! /!\ so that the matrix already contains the intensity constant of the dipolar moment with the position of the matter (cf. manual for formulas)
    
  IF (Debug) CALL Write_Mat(Mol1DipMomt%Band_val, out_unit, 3, info="Mol1DipMomt")
  FLUSH(out_unit)

    !------------------------------------------------------First cavity mode-----------------------------------------------------
  WRITE(out_unit,*); WRITE(out_unit,*) "  --------------------------------------------------------First cavity mode--------------&
                                     &-----------------------------------------"
  CALL MolecCav_Read_cavity_mode(Mode=Cavity_mode_1, nio=in_unit)
  Cavity_mode_1%Nb     = 11
  Cavity_mode_1%w      = 0.0068665_Rkind
  Cavity_mode_1%m      = ONE
  Cavity_mode_1%lambda = ONETENTH

  WRITE(out_unit,*) "Cavity mode Hamiltonian :"
  CALL Construct_Operator_1D(Operator=Cav1H, &
                           & operator_type="Hamiltonian", &
                           & Mode=Cavity_mode_1, &
                           & Debug=.TRUE.)

  WRITE(out_unit,*) "Cavity mode Position :"
  CALL Construct_Operator_1D(Operator=Cav1Position, &
                           & operator_type="Position", &
                           & Mode=Cavity_mode_1, &
                           & Debug=Debug)

  FLUSH(out_unit)

  Nb_M = Molecule_1%Nb
  Nb_C = Cavity_mode_1%Nb
  NB   = Molecule_1%Nb * Cavity_mode_1%Nb
      
    !-------------------------------Construction of the Total Hamiltonian matrix with CM-couplings-------------------------------
  WRITE(out_unit,*); WRITE(out_unit,*) "  ---------------------------------Construction of the Total Hamiltonian matrix with CM-c&
                                      &ouplings---------------------------------"

  ALLOCATE(TotH(NB, NB))
  CALL Construct_total_hamiltonian_1p1D(TotH, Cav1Position, Cav1H, Mol1DipMomt, Mol1H, Debug=.FALSE.)

  IF (Debug .AND. NB <= 10) THEN
    WRITE(out_unit,*); WRITE(out_unit,*) "Total Hamiltonian 1p1D (lambda /= 0, w_C /= w_M)"
    CALL Write_Mat(TotH, out_unit, NB, info="TotH")
  ELSE IF (Debug .AND. NB > 10) THEN
    WRITE(out_unit,*); WRITE(out_unit,*) "Total Hamiltonian 1p1D (lambda /= 0, w_C /= w_M) (10:10 slicing)"
    CALL Write_Mat(TotH(1:10,1:10), out_unit, 10, info="TotH (sliced)")
  END IF

  IF (View_more) WRITE(out_unit,*)
  IF (View_more) CALL Write_Mat(TotH(1:NB, 1:NB), out_unit, Size(TotH), info="TotH(FULL)")

    !-------------------------------------------------Computation of Eigenstates-------------------------------------------------
  ALLOCATE(REigval(NB))
  ALLOCATE(REigvec(NB,NB))
  CALL diagonalization(TotH, REigval, REigvec)

  WRITE(out_unit,*); WRITE(out_unit,*) 'EIGENVALUES'
  IF (NB <= 10) CALL WRITE_Vec(REigval, out_unit, 10, info = 'VP[Ha]')
  IF (NB > 10)  CALL WRITE_Vec(REigval(1:10), out_unit, 10, info = 'Ten_first_VP[Ha]')

  IF (Debug .AND. NB <= 10) THEN
    WRITE(out_unit,*); WRITE(out_unit,*) 'EIGENVECTORS'
    CALL WRITE_Mat(REigvec, out_unit, 6, info = 'Eigenvectors')
  ELSE IF (Debug .AND. NB > 10) THEN
    WRITE(out_unit,*); WRITE(out_unit,*) 'EIGENVECTORS'
    CALL WRITE_Mat(REigvec(1:10,1:10), out_unit, 6, info = 'Ten first Eigenvectors (1:10 slicing)')
  END IF 

  IF (View_more) WRITE(out_unit,*)
  IF (View_more) CALL Write_Mat(REigvec(1:NB, 1:NB), out_unit, Size(REigvec), info="REigvec(FULL)")

  WRITE(out_unit,*) "--------------------------------------------------FINISHED SYSTEM INITIALIZATION----------------------------&
                    &----------------------"


  !-------------------------------------------Computing the Transition intensity matrix------------------------------------------
  WRITE(out_unit,*); WRITE(out_unit,*) "---------------------------------------------Computing the transition intensity matrix---&
                                       &-----------------------------------------"
  WRITE(out_unit,*)

  CALL Initialize_transition_matrix(Intensities1, Mol1DipMomt, REigvec, Nb_states=NB/12, Debug=Debug)
  CALL Compute_transition_matrix(Intensities1, Mol1DipMomt, REigvec, Debug=Debug)

  CALL Initialize_transition_matrix(Intensities2, Mol1DipMomt, REigvec, Energy_threshold=0.020_Rkind, REigval=REigval, Debug=Debug)
  CALL Compute_transition_matrix(Intensities2, Mol1DipMomt, REigvec, Debug=Debug)


  !--------------------------------------------Testing the Transition intensity matrix-------------------------------------------
  N_min = MIN(Size(Intensities1, dim=1), Size(Intensities2, dim=1))
  CALL Equal_R_R_matrix(error_trnstn_int, Intensities1(1:N_min,1:N_min), Intensities2(1:N_min,1:N_min))
  CALL Logical_Test(test_trnstn_int, error_trnstn_int, test2=.FALSE., info="Both selection methods give same result ?")
  CALL Finalize_Test(test_trnstn_int)
  
  
  CONTAINS


  SUBROUTINE Equal_R_R_matrix(error, Rl_1, Rl_2)
    USE QDUtil_m
    IMPLICIT NONE 

    logical,          intent(inout) :: error
    real(kind=Rkind), intent(in)    :: Rl_1(:,:)
    real(kind=Rkind), intent(in)    :: Rl_2(:,:)
    
    real(kind=Rkind), parameter     :: Threshold   = 1E-10_Rkind
    logical, parameter              :: Debug_local = .FALSE.
    integer                         :: Nb_1, Nb_2

    Nb_1 = Size(Rl_1, dim=1)
    Nb_2 = Size(Rl_2, dim=2)
    IF (Nb_1 /= Size(Rl_2, dim=1) .OR. Nb_2 /= Size(Rl_2, dim=2)) THEN
      WRITE(out_unit,*) "The two matrices must have same dimensions to compare them. Please, check initialization."
      STOP "### (Equal_R_R_matrix) The two matrices must have same dimensions to compare them. Please, check initialization."
    END IF 

    IF (ANY(ABS(Rl_1 - Rl_2) > Threshold)) THEN
      error = .TRUE.
      IF (Debug_local) THEN
        WRITE(out_unit,*) "The two matrices are not close enough to be considered equal :"
        CALL Write_Mat(Rl_1, out_unit, Nb_2, info="R_1(:,:)")
        CALL Write_Mat(Rl_2, out_unit, Nb_2, info="R_2(:,:)")
        CALL Write_Mat(ABS(Rl_1 - Rl_2), out_unit, Nb_2, info="|R_1-R_2| = ")
      END IF 

    ELSE 
      error = .FALSE.
      IF (Debug_local) THEN
        WRITE(out_unit,*) "The two matrices are close enough to be considered equal :"
        CALL Write_Mat(Rl_1, out_unit, Nb_2, info="R_1(:,:)")
        CALL Write_Mat(Rl_2, out_unit, Nb_2, info="R_2(:,:)")
        CALL Write_Mat(ABS(Rl_1 - Rl_2), out_unit, Nb_2, info="|R_1-R_2| = ")
      END IF
    END IF 

  END SUBROUTINE Equal_R_R_matrix


END PROGRAM