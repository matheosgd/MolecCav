!==================================================================================================
!==================================================================================================
!This file is part of MolecCav.
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
PROGRAM test_normal_modes_1p1D
  !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
  USE QDUtil_m
  USE QDUtil_Test_m
  USE Algebra_m
  USE Cavity_mode_m
  USE Operator_1D_m
  USE Total_hamiltonian_m
  IMPLICIT NONE
  

  logical, parameter            :: Debug =     .TRUE.
  logical, parameter            :: View_more = .FALSE.

  !-------------Diatomic molecule in a HARMONIC electonic potential------------
  TYPE(Cavity_mode_t)           :: Molecule_1
  TYPE(Operator_1D_t)           :: Mol1H                                       ! matrix of the one-dimensional harmonic Hamiltonian associated with HO D
  TYPE(Operator_1D_t)           :: Mol1Position
  TYPE(Operator_1D_t)           :: Mol1N                                       ! here N does not count the number of photons but of excitation quanta in the vibrational state
  TYPE(Operator_1D_t)           :: Mol1_dipolar_moment
  real(kind=Rkind)              :: Cte_dipole_moment = ONE                     ! the intensity of the variation of the dipole moment with a variation of the matter DOF
  
  !------------------------------First cavity mode-----------------------------
  TYPE(Cavity_mode_t)           :: Cavity_mode_1
  TYPE(Operator_1D_t)           :: Cav1H                                       ! matrix of the one-dimensional harmonic Hamiltonian associated with HO D
  TYPE(Operator_1D_t)           :: Cav1Position
  TYPE(Operator_1D_t)           :: Cav1N

  TYPE(Cavity_mode_t)           :: Cavity_mode_1_uncoupled
  TYPE(Operator_1D_t)           :: Cav1H_uncloupled                            ! matrix of the one-dimensional harmonic Hamiltonian associated with HO D
  TYPE(Operator_1D_t)           :: Cav1Position_uncoupled

  !------------------------------Total Hamiltonian-----------------------------
  real(kind=Rkind), allocatable :: TotH_uncoupled(:,:)
  real(kind=Rkind), allocatable :: TotH(:,:)
  real(kind=Rkind)              :: MWH_uncoupled(2,2)                          ! the mass-weighted Hessian matrix of the total, 1p1D, coupled system [matter x cavity] in harmonic approximation of the matter                     
  real(kind=Rkind)              :: MWH(2,2)                                    ! the mass-weighted Hessian matrix of the total, 1p1D, coupled system [matter x cavity] in harmonic approximation of the matter                     

  !-----------------------------------Results----------------------------------
  real(kind=Rkind), allocatable :: REigval_uncoupled(:)
  real(kind=Rkind), allocatable :: REigvec_uncoupled(:,:)
  real(kind=Rkind)              :: Normal_modes_uncoupled(2)                   ! VP of the MWH
  real(kind=Rkind)              :: Normal_coordinates_uncoupled(2,2)           ! \overrightarrow{VP} pf the MWH
  real(kind=Rkind), allocatable :: REigval(:)
  real(kind=Rkind), allocatable :: REigvec(:,:)
  real(kind=Rkind)              :: Normal_modes(2)                             ! VP of the MWH
  real(kind=Rkind)              :: Normal_coordinates(2,2)                     ! \overrightarrow{VP} pf the MWH

  !----------------------------------Utilities---------------------------------
  integer                       :: Nb_M, Nb_C, NB, i
  real(kind=Rkind)              :: TotE, w_1, w_2

  TYPE(test_t)                  :: test_nrml_mdes
  logical                       :: error_nrml_mdes = .FALSE.

  !---------------------------------------Test initialization--------------------------------------
  CALL Initialize_Test(test_nrml_mdes, test_name="OUT/test_file_nrml_mdes_1p1D")

  !----------------------------SYSTEM INITIALIZATION---------------------------
    !------------Diatomic molecule in a harmonic electonic potential-----------
  WRITE(out_unit,*) "-----------------------------SYSTEM INITIALIZATION----------------------------"
  WRITE(out_unit,*) "  -------------Diatomic molecule in a harmonic electonic potential------------"
  CALL MolecCav_Read_cavity_mode(Mode=Molecule_1, nio=in_unit)
  Molecule_1%Nb = 10
  Molecule_1%w  = 0.15_Rkind
  Molecule_1%m  = 1744.60504565_Rkind

  WRITE(out_unit,*) "Molecular Hamiltonian :"
  CALL Construct_Operator_1D(Operator=Mol1H, &
                           & operator_type="Hamiltonian", &
                           & Mode=Molecule_1, &
                           & Debug=.TRUE.)
  
  WRITE(out_unit,*) "Molecular Position :"
  CALL Construct_Operator_1D(Operator=Mol1Position, &
                           & operator_type="Position", &
                           & Mode=Molecule_1, &
                           & Debug=Debug)

  WRITE(out_unit,*) "Molecular Nb_photon (vibrationnal excitation quanta)"
  CALL Construct_Operator_1D(Operator=Mol1N, &
                           & operator_type="Nb_photons", &
                           & Mode=Molecule_1, &
                           & Debug=Debug)
  
  WRITE(out_unit,*) "Molecular Dipole moment :"
  CALL Construct_Operator_1D(Operator=Mol1_dipolar_moment, &
                           & operator_type="Position", &                       ! initialized as a position operator because of approximation over its expression (cf. readme.md or manual)
                           & Mode=Molecule_1, &
                           & Debug=.FALSE.)

  Mol1_dipolar_moment%Band_val_R = Mol1_dipolar_moment%Band_val_R*Cte_dipole_moment ! /!\ so that the matrix already contains the intensity constant of the dipolar moment with the position of the matter (cf. manual for formulas)
    
  IF (Debug .AND. ALLOCATED(Mol1_dipolar_moment%Diag_val_R)) CALL Write_Vec(Mol1_dipolar_moment%Diag_val_R, out_unit, 3, info="Mo&
                                                            &l1_dipolar_moment")
  IF (Debug .AND. (ALLOCATED(Mol1_dipolar_moment%Band_val_R) .OR. ALLOCATED(Mol1_dipolar_moment%Dense_val_R))) CALL Write_Mat(Mol&
                                                            &1_dipolar_moment%Band_val_R, out_unit, 3, info="Mol1_dipolar_moment")
  FLUSH(out_unit)

    !-----------------------------First cavity mode----------------------------
  WRITE(out_unit,*); WRITE(out_unit,*) "  -----------------------------First cavity mode----------------------------"
  CALL MolecCav_Read_cavity_mode(Mode=Cavity_mode_1, nio=in_unit)
  Cavity_mode_1%Nb     = 11
  Cavity_mode_1%w      = 0.17_Rkind
  Cavity_mode_1%m      = ONE
  Cavity_mode_1%lambda = ONE


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

  WRITE(out_unit,*) "Cavity mode Nb_photons :"
  CALL Construct_Operator_1D(Operator=Cav1N, &
                           & operator_type="Nb_photons", &
                           & Mode=Cavity_mode_1, &
                           & Debug=Debug)
  FLUSH(out_unit)

  Nb_M = Molecule_1%Nb
  Nb_C = Cavity_mode_1%Nb
  NB   = Molecule_1%Nb * Cavity_mode_1%Nb

    !---------------------------uncoupled cavity mode--------------------------
  WRITE(out_unit,*); WRITE(out_unit,*) "  -----------------------------First cavity mode----------------------------"
  Cavity_mode_1_uncoupled = Cavity_mode_1
  Cavity_mode_1_uncoupled%lambda = ZERO
  CALL Construct_Operator_1D(Cav1H_uncloupled,       "Hamiltonian", Mode=Cavity_mode_1_uncoupled, Debug=.FALSE.)
  CALL Construct_Operator_1D(Cav1Position_uncoupled, "Position",    Mode=Cavity_mode_1_uncoupled, Debug=.FALSE.)

  WRITE(out_unit,*) "------------------------FINISHED SYSTEM INITIALIZATION------------------------"


    !-------Construction of a Total Hamiltonian matrix without CM-couplings------
  WRITE(out_unit,*); WRITE(out_unit,*) "-------Construction of a Total Hamiltonian matrix without CM-couplings------"
  ALLOCATE(TotH_uncoupled(NB, NB))
  
  CALL Construct_total_hamiltonian_1p1D(TotH_uncoupled, Cav1Position_uncoupled, Cav1H_uncloupled, Mol1_dipolar_moment, Mol1H, Deb&
                                       &ug=.FALSE.)

  IF (Debug .AND. NB <= 10) THEN
    WRITE(out_unit,*); WRITE(out_unit,*) "Uncoupled total Hamiltonian 1p1D (lambda = 0, w_C /= w_M)"
    CALL Write_Mat(TotH_uncoupled,            out_unit, NB, info="TotH_uncoupled")
  ELSE 
    WRITE(out_unit,*); WRITE(out_unit,*) "Uncoupled total Hamiltonian 1p1D (lambda = 0, w_C /= w_M) (10:10 slicing)"
    CALL Write_Mat(TotH_uncoupled(1:10,1:10), out_unit, 10, info="TotH_uncoupled (sliced)")
  END IF

    !---------------------------Computation Eigenstates--------------------------
  ALLOCATE(REigval_uncoupled(NB))
  ALLOCATE(REigvec_uncoupled(NB,NB))
  CALL diagonalization(TotH_uncoupled, REigval_uncoupled, REigvec_uncoupled)

  WRITE(out_unit,*); WRITE(out_unit,*) 'EIGENVALUES (\lambda = 0)'
  IF (NB <= 10) CALL WRITE_Vec(REigval_uncoupled,       out_unit, 10, info = 'VP[Ha]')
  IF (NB > 10)  CALL WRITE_Vec(REigval_uncoupled(1:10), out_unit, 10, info = 'Ten_first_VP[Ha]')

  IF (Debug .AND. NB <= 10) THEN
    WRITE(out_unit,*); WRITE(out_unit,*) 'EIGENVECTORS (\lambda = 0)'
    CALL WRITE_Mat(REigvec_uncoupled,            out_unit, 6, info = 'Eigenvectors')
  ELSE IF (Debug .AND. NB > 10) THEN
    WRITE(out_unit,*); WRITE(out_unit,*) 'EIGENVECTORS (\lambda = 0)'
    CALL WRITE_Mat(REigvec_uncoupled(1:10,1:10), out_unit, 6, info = 'Ten first Eigenvectors (1:10 slicing)')
  END IF 

    !--Construction of the 1p1D uncoupled system Mass-weighted Hessian matrix--
  CALL Compute_normal_modes(Normal_modes_uncoupled, Normal_coordinates_uncoupled, MWH_uncoupled, Molecule_1, Cavity_mode_1_uncoup&
                           &led, Cte_dipole_moment)


  !-------Construction of the Total Hamiltonian matrix with CM-couplings-------
  WRITE(out_unit,*); WRITE(out_unit,*) "-------Construction of the Total Hamiltonian matrix with CM-couplings-------"

  ALLOCATE(TotH(NB, NB))
  CALL Construct_total_hamiltonian_1p1D(TotH, Cav1Position, Cav1H, Mol1_dipolar_moment, Mol1H, Debug=.FALSE.)

  IF (Debug .AND. NB <= 10) THEN
    WRITE(out_unit,*); WRITE(out_unit,*) "Total Hamiltonian 1p1D (lambda /= 0, w_C /= w_M)"
    CALL Write_Mat(TotH, out_unit, NB, info="TotH")
  ELSE IF (Debug .AND. NB > 10) THEN
    WRITE(out_unit,*); WRITE(out_unit,*) "Total Hamiltonian 1p1D (lambda /= 0, w_C /= w_M) (10:10 slicing)"
    CALL Write_Mat(TotH(1:10,1:10), out_unit, 10, info="TotH (sliced)")
  END IF

  IF (View_more) WRITE(out_unit,*)
  IF (View_more) CALL Write_Mat(TotH(1:NB, 1:NB), out_unit, Size(TotH), info="TotH(1:150)")

    !------------------------Computation of Eigenstates------------------------
  ALLOCATE(REigval(NB))
  ALLOCATE(REigvec(NB,NB))
  CALL diagonalization(TotH, REigval, Reigvec)

  WRITE(out_unit,*); WRITE(out_unit,*) 'EIGENVALUES'
  IF (NB <= 10) CALL WRITE_Vec(Reigval, out_unit, 10, info = 'VP[Ha]')
  IF (NB > 10)  CALL WRITE_Vec(Reigval(1:10), out_unit, 10, info = 'Ten_first_VP[Ha]')

  IF (Debug .AND. NB <= 10) THEN
    WRITE(out_unit,*); WRITE(out_unit,*) 'EIGENVECTORS'
    CALL WRITE_Mat(Reigvec, out_unit, 6, info = 'Eigenvectors')
  ELSE IF (Debug .AND. NB > 10) THEN
    WRITE(out_unit,*); WRITE(out_unit,*) 'EIGENVECTORS'
    CALL WRITE_Mat(Reigvec(1:10,1:10), out_unit, 6, info = 'Ten first Eigenvectors (1:10 slicing)')
  END IF 

    !-----Construction of the 1p1D total system Mass-weighted Hessian matrix-----
  CALL Compute_normal_modes(Normal_modes, Normal_coordinates, MWH, Molecule_1, Cavity_mode_1, Cte_dipole_moment)

  
  !----------------------------Testing the eigenstates---------------------------
    ! Are the eigenpulsations of the total system's normal modes (w_1, w_2) i.e. the sqrt of the MWH's Eigenvalues equal to the eigenpulsations of the two separated subsystems, when lambda = 0 ?
  CALL Equal_R_R_vector(error_nrml_mdes, SQRT(Normal_modes_uncoupled), [MIN(Molecule_1%w, Cavity_mode_1%w), MAX(Molecule_1%w, Cav&
                                                                       &ity_mode_1%w)])
  CALL Logical_Test(test_nrml_mdes, error_nrml_mdes, test2=.FALSE., info="Uncloupled system: NM' MWH =?= eigenpulsations")

    ! Are the energies of the eigenstates of the total system i.e. the Eigenvalues of TotH(:,:) matching the Energies computed with the analytical expression from the eigenpulsations of the total system's normal modes (w_1, w_2) ?
    ! TotE(i_1, i_2) = \barh.w_1(i_1 + 1/2) + \barh.w_2(i_2 + 1/2)
      ! Without coupling :
  TotE = Total_Eigenenergy_1p1D(0,0,SQRT(Normal_modes_uncoupled(1)),SQRT(Normal_modes_uncoupled(2)))
  CALL Equal_R_R_scalar(error_nrml_mdes, TotE, REigval_uncoupled(1))
  CALL Logical_Test(test_nrml_mdes, error_nrml_mdes, test2=.FALSE., info="Uncoupled system TotE")
  IF (error_nrml_mdes .AND. Debug) THEN 
    WRITE(out_unit,*) "Total_Eigenenergy_1p1D(0,0) = "//TO_string(TotE)//"; E_{00} = "//TO_string(REigval_uncoupled(1))
  END IF 

  TotE = Total_Eigenenergy_1p1D(1,0,SQRT(Normal_modes_uncoupled(1)),SQRT(Normal_modes_uncoupled(2)))
  CALL Equal_R_R_scalar(error_nrml_mdes, TotE, REigval_uncoupled(2))
  CALL Logical_Test(test_nrml_mdes, error_nrml_mdes, test2=.FALSE., info="Uncoupled system TotE")
  IF (error_nrml_mdes .AND. Debug) THEN 
    WRITE(out_unit,*) "Total_Eigenenergy_1p1D(1,0) = "//TO_string(TotE)//"; E_{10} = "//TO_string(REigval_uncoupled(2))
  END IF 

  TotE = Total_Eigenenergy_1p1D(0,1,SQRT(Normal_modes_uncoupled(1)),SQRT(Normal_modes_uncoupled(2)))
  CALL Equal_R_R_scalar(error_nrml_mdes, TotE, REigval_uncoupled(3))
  CALL Logical_Test(test_nrml_mdes, error_nrml_mdes, test2=.FALSE., info="Uncoupled system TotE")
  IF (error_nrml_mdes .AND. Debug) THEN 
    WRITE(out_unit,*) "Total_Eigenenergy_1p1D(0,1) = "//TO_string(TotE)//"; E_{01} = "//TO_string(REigval_uncoupled(3))
  END IF 

      ! With coupling :
  TotE = Total_Eigenenergy_1p1D(0,0,SQRT(Normal_modes(1)),SQRT(Normal_modes(2)))
  CALL Equal_R_R_scalar(error_nrml_mdes, TotE, REigval(1))
  CALL Logical_Test(test_nrml_mdes, error_nrml_mdes, test2=.FALSE., info="Coupled system TotE")
  IF (error_nrml_mdes .AND. Debug) THEN 
    WRITE(out_unit,*) "Total_Eigenenergy_1p1D(0,0) = "//TO_string(TotE)//"; E_{00} = "//TO_string(REigval(1))
  END IF 

  TotE = Total_Eigenenergy_1p1D(1,0,SQRT(Normal_modes(1)),SQRT(Normal_modes(2)))
  CALL Equal_R_R_scalar(error_nrml_mdes, TotE, REigval(2))
  CALL Logical_Test(test_nrml_mdes, error_nrml_mdes, test2=.FALSE., info="Coupled system TotE")
  IF (error_nrml_mdes .AND. Debug) THEN 
    WRITE(out_unit,*) "Total_Eigenenergy_1p1D(1,0) = "//TO_string(TotE)//"; E_{10} = "//TO_string(REigval(2))
  END IF 

  TotE = Total_Eigenenergy_1p1D(0,1,SQRT(Normal_modes(1)),SQRT(Normal_modes(2)))
  CALL Equal_R_R_scalar(error_nrml_mdes, TotE, REigval(3))
  CALL Logical_Test(test_nrml_mdes, error_nrml_mdes, test2=.FALSE., info="Coupled system TotE")
  IF (error_nrml_mdes .AND. Debug) THEN 
    WRITE(out_unit,*) "Total_Eigenenergy_1p1D(0,1) = "//TO_string(TotE)//"; E_{01} = "//TO_string(REigval(3))
  END IF 
    
    ! Are the eigenpulsations of the total system's normal modes (w_1, w_2) computed analytically from the Eigenvalues of TotH(:,:) matching the eigenpulsation computed from the sqrt of the MWH's Eigenvalues:
    ! w_1 = TotE(1,0) - TotE(0,0) = E_1 - E_0; w_2 = TotE(0,1) - TotE(0,0) = E_2 - E_0 
  w_1 = REigval(2) - REigval(1)
  CALL Equal_R_R_scalar(error_nrml_mdes, w_1, SQRT(Normal_modes(1)))
  CALL Logical_Test(test_nrml_mdes, error_nrml_mdes, test2=.FALSE., info="E_1 - E_0 =?= w_1")
  IF (error_nrml_mdes .AND. Debug) THEN 
    WRITE(out_unit,*) "E_1 - E_0 = "//TO_string(w_1)//"; w_1 = "//TO_string(SQRT(Normal_modes(1)))
  END IF 

  w_2 = REigval(3) - REigval(1)
  CALL Equal_R_R_scalar(error_nrml_mdes, w_2, SQRT(Normal_modes(2)))
  CALL Logical_Test(test_nrml_mdes, error_nrml_mdes, test2=.FALSE., info="E_2 - E_0 =?= w_2")
  IF (error_nrml_mdes .AND. Debug) THEN 
    WRITE(out_unit,*) "E_2 - E_0 = "//TO_string(w_1)//"; w_2 = "//TO_string(SQRT(Normal_modes(2)))
  END IF 


  CALL Finalize_Test(test_nrml_mdes)


  CONTAINS


  SUBROUTINE Compute_normal_modes(Nrml_modes, Nrml_coordinates, MWH_local, DOF_1, DOF_2, Coeff_dipole_moment)
    USE QDUtil_m
    USE Cavity_mode_m
    IMPLICIT NONE 

    real(kind=Rkind), intent(inout)              :: Nrml_modes(2)                             ! VP of the MWH
    real(kind=Rkind), intent(inout)              :: Nrml_coordinates(2,2)                     ! \overrightarrow{VP} of the MWH
    real(kind=Rkind), intent(inout)              :: MWH_local(2,2)                            ! the Mass-Weighted Hessian matrix
    TYPE(Cavity_mode_t), intent(in)              :: DOF_1                                     ! the molecule in this case
    TYPE(Cavity_mode_t), intent(in)              :: DOF_2                                     ! the cavity mode in this case
    real(kind=Rkind)                             :: Coeff_dipole_moment
    logical, parameter                           :: Debug_local = .TRUE.

    WRITE(out_unit,*); WRITE(out_unit,*) "-----Construction of the 1p1D total system Mass-weighted Hessian matrix-----"
    MWH_local(1,1) = DOF_2%w**2
    MWH_local(2,2) = DOF_1%w**2
    MWH_local(1,2) = DOF_2%lambda*DOF_2%w*Coeff_dipole_moment / SQRT(DOF_1%m)
    MWH_local(2,1) = MWH_local(1,2)

    IF (Debug_local) THEN
      CALL Write_Mat(MWH_local, out_unit, Size(MWH_local, dim=2), info="MWH")
    END IF

    CALL diagonalization(MWH_local, Nrml_modes, Nrml_coordinates)
    CALL Write_Vec(Nrml_modes,       out_unit, Size(Nrml_modes),              info="Normal modes")
    CALL Write_Mat(Nrml_coordinates, out_unit, Size(Nrml_coordinates, dim=2), info="Normal coordniates")

    WRITE(out_unit,*)
    DO i = 1, Size(Nrml_modes)
      IF (Nrml_modes(i) >= 0) THEN
        WRITE(out_unit,*) TO_string(i)//"^{th} Normal coordinate has positive squared frequency, lead&
                         &ing to w_"//TO_string(i)//" = "//TO_string(SQRT(Nrml_modes(i)))
      ELSE
        WRITE(out_unit,*) TO_string(i)//"^{th} Normal coordinate has NEGATIVE squared frequency, lead&
        &ing to w_"//TO_string(i)//" = "//TO_string(EYE*SQRT(-Nrml_modes(i)))
      END IF
    END DO
    WRITE(out_unit,*) "Expected ZPE by half-sum of the total system eigenpulsations (from MWH): "//TO&
                    &_string( ( SQRT(Nrml_modes(1))+SQRT(Nrml_modes(2)) )/2 )

  END SUBROUTINE Compute_normal_modes


  FUNCTION Total_Eigenenergy_1p1D(i_1, i_2, w_1_local, w_2_local) RESULT(TotE_local)
    USE QDUtil_m
    USE Cavity_mode_m
    IMPLICIT NONE
 
    real(kind=Rkind)    :: TotE_local
    integer,          intent(in) :: i_1
    integer,          intent(in) :: i_2
    real(kind=Rkind), intent(in) :: w_1_local
    real(kind=Rkind), intent(in) :: w_2_local

    TotE = w_1_local*(REAL(i_1, kind=Rkind) + HALF) + w_2_local*(REAL(i_2, kind=Rkind) +  HALF)

  END FUNCTION Total_Eigenenergy_1p1D


  SUBROUTINE Equal_R_R_scalar(error, Rl_1, Rl_2)
    USE QDUtil_m
    IMPLICIT NONE 

    logical,          intent(inout) :: error
    real(kind=Rkind), intent(in)    :: Rl_1
    real(kind=Rkind), intent(in)    :: Rl_2
    
    real(kind=Rkind), parameter     :: Threshold   = 1E-10_Rkind
    logical, parameter              :: Debug_local = .FALSE.

    IF (ABS(Rl_1 - Rl_2) > Threshold) THEN
      error = .TRUE.
      IF (Debug_local) WRITE(out_unit,*) "The two numbers are not close enough to be considered equ&
                                         &al : R_1 =", Rl_1, "R_2 =", Rl_2, "|R_1-R_2| = ", ABS(Rl_1 - Rl_2)
    ELSE 
      error = .FALSE.
      IF (Debug_local) WRITE(out_unit,*) "The two numbers are close enough to be considered equal :&
                                         & R_1 =", Rl_1, "R_2 =", Rl_2, "|R_1-R_2| = ", ABS(Rl_1 - Rl_2)
    END IF

  END SUBROUTINE Equal_R_R_scalar
  

  SUBROUTINE Equal_R_R_vector(error, Rl_1, Rl_2)
    USE QDUtil_m
    IMPLICIT NONE 

    logical,          intent(inout) :: error
    real(kind=Rkind), intent(in)    :: Rl_1(:)
    real(kind=Rkind), intent(in)    :: Rl_2(:)
    
    real(kind=Rkind), parameter     :: Threshold   = 1E-10_Rkind
    logical, parameter              :: Debug_local = .FALSE.
    integer                         :: Nb_1

    Nb_1 = Size(Rl_1)
    IF (Nb_1 /= Size(Rl_2)) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "The two vectors must have same dimensions to compare them. Please, check initialization."
      WRITE(out_unit,*) "Size(Rl_1) = "//TO_string(Size(Rl_1))//"; Size(Rl_2) = "//TO_string(Size(Rl_2))
      STOP "### The two vectors must have same dimensions to compare them. Please, check initialization."
    END IF 

    IF (ANY(ABS(Rl_1 - Rl_2) > Threshold)) THEN
      error = .TRUE.
      IF (Debug_local) THEN
        WRITE(out_unit,*) "The two vectors are not close enough to be considered equal :"
        CALL Write_Vec(Rl_1, out_unit, Nb_1, info="R_1(:)")
        CALL Write_Vec(Rl_2, out_unit, Nb_1, info="R_2(:)")
        CALL Write_Vec(ABS(Rl_1 - Rl_2), out_unit, Nb_1, info="|R_1-R_2| = ")
      END IF 

    ELSE 
      error = .FALSE.
      IF (Debug_local) THEN
        WRITE(out_unit,*) "The two vectors are close enough to be considered equal :"
        CALL Write_Vec(Rl_1, out_unit, Nb_1, info="R_1(:)")
        CALL Write_Vec(Rl_2, out_unit, Nb_1, info="R_2(:)")
        CALL Write_Vec(ABS(Rl_1 - Rl_2), out_unit, Nb_1, info="|R_1-R_2| = ")
      END IF
    END IF 

  END SUBROUTINE Equal_R_R_vector


END PROGRAM