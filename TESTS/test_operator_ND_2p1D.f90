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
PROGRAM test_operator_ND_2p1D
  !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real124
  USE QDUtil_m
  USE Tests_m
  USE Algebra_m
  USE Operator_ND_m
  USE Total_hamiltonian_m
  IMPLICIT NONE


  integer             :: Verbose = 50
  logical             :: Debug   = .TRUE.

  logical             :: Dense   = .FALSE.
  TYPE(Operator_ND_t) :: HxIxI
  TYPE(Operator_ND_t) :: IxHxI
  TYPE(Operator_ND_t) :: IxIxH
  TYPE(Operator_ND_t) :: DipMomtxIxPos
  TYPE(Operator_ND_t) :: IxDipMomtxPos
  real(kind=Rkind)    :: Mat1w
  real(kind=Rkind)    :: Mat2w
  real(kind=Rkind)    :: Mat1m
  real(kind=Rkind)    :: Mat2m
  real(kind=Rkind)    :: Cavw
  real(kind=Rkind)    :: Mat1lambda
  real(kind=Rkind)    :: Mat2lambda
  real(kind=Rkind)    :: Cavlambda
  real(kind=Rkind)    :: CoeffDipMomt1
  real(kind=Rkind)    :: CoeffDipMomt2
  integer             :: Mat1Nb
  integer             :: Mat2Nb
  integer             :: CavNb
  integer             :: NB

  real(kind=Rkind), allocatable ::    Phi(:)                                                                             ! an any vector representing the excitation state/wavefunction = a linear combination of the basis functions from the canonical basis set over \mathbb{R} /!\ Not normalized yet !
  real(kind=Rkind), allocatable :: Op_phi(:)                                                                             ! an any vector representing the excitation state/wavefunction = a linear combination of the basis functions from the canonical basis set over \mathbb{R} /!\ Not normalized yet !
  real(kind=Rkind), allocatable :: TotH(:,:)                                                                             ! an any vector representing the excitation state/wavefunction = a linear combination of the basis functions from the canonical basis set over \mathbb{R} /!\ Not normalized yet !
  real(kind=Rkind), allocatable :: Eigenenergies(:)                                                                             ! an any vector representing the excitation state/wavefunction = a linear combination of the basis functions from the canonical basis set over \mathbb{R} /!\ Not normalized yet !
  real(kind=Rkind), allocatable :: Eigenstates(:,:)                                                                             ! an any vector representing the excitation state/wavefunction = a linear combination of the basis functions from the canonical basis set over \mathbb{R} /!\ Not normalized yet !

  real(kind=Rkind)              :: Ana_Normal_modes(0:2)
  real(kind=Rkind), allocatable :: Ana_Eigenenergies(:)

  TYPE(test_t)        :: test_opnd
  logical             :: error_opnd = .FALSE.

  integer             :: J, i_1, i_2, i_3, min_index


  !-----------------------------Test initialization----------------------------
  CALL Initialize_Test(test_opnd, test_name="OUT/test_file_opnd_2p1D")


  !-------------------------Operator_ND object initialization-------------------------
  CALL Initialize(HxIxI, "Hamiltonian, identity ", "identity", in_unit, Dense=Dense, Verbose=Verbose, Debug=.FALSE.)
  IF (Debug) THEN
    WRITE(out_unit,*)
    WRITE(out_unit,*) "--------------Operator_ND object constructed by MolecCav_Initialize_operator_ND--------------"
    CALL Write(HxIxI)
    WRITE(out_unit,*) "------------End Operator_ND object constructed by MolecCav_Initialize_operator_ND------------"
  END IF

  CALL Initialize(IxHxI, "identity, Hamiltonian ", "identity", in_unit, Dense=Dense, Verbose=Verbose, Debug=.FALSE.)
  IF (Debug) THEN
    WRITE(out_unit,*)
    WRITE(out_unit,*) "--------------Operator_ND object constructed by MolecCav_Initialize_operator_ND--------------"
    CALL Write(IxHxI)
    WRITE(out_unit,*) "------------End Operator_ND object constructed by MolecCav_Initialize_operator_ND------------"
  END IF

  CALL Initialize(IxIxH, "identity, identity ", "Hamiltonian", in_unit, Dense=Dense, Verbose=Verbose, Debug=.FALSE.)
  IF (Debug) THEN
    WRITE(out_unit,*)
    WRITE(out_unit,*) "--------------Operator_ND object constructed by MolecCav_Initialize_operator_ND--------------"
    CALL Write(IxIxH)
    WRITE(out_unit,*) "------------End Operator_ND object constructed by MolecCav_Initialize_operator_ND------------"
  END IF

  CALL Initialize(DipMomtxIxPos, "DipMomt, Identity", "position", in_unit, Dense=Dense, Verbose=Verbose, Debug=.FALSE.)
  IF (Debug) THEN
    WRITE(out_unit,*)
    WRITE(out_unit,*) "--------------Operator_ND object constructed by MolecCav_Initialize_operator_ND--------------"
    CALL Write(DipMomtxIxPos)
    WRITE(out_unit,*) "------------End Operator_ND object constructed by MolecCav_Initialize_operator_ND------------"
  END IF

  CALL Initialize(IxDipMomtxPos, "identity, DipMomt", "Position ", in_unit, Dense=Dense, Verbose=Verbose, Debug=.FALSE.)
  IF (Debug) THEN
    WRITE(out_unit,*)
    WRITE(out_unit,*) "--------------Operator_ND object constructed by MolecCav_Initialize_operator_ND--------------"
    CALL Write(IxDipMomtxPos)
    WRITE(out_unit,*) "------------End Operator_ND object constructed by MolecCav_Initialize_operator_ND------------"
  END IF


  !----------------------------Testing the initialization---------------------------
  CALL Logical_Test(test_opnd, ANY(HxIxI%tab_indexes_mat_op/=[1,0]), test2=.FALSE., info="HxIxI%tab_mat_op")
  CALL Logical_Test(test_opnd, ANY(HxIxI%tab_indexes_cav_op/=[0]  ), test2=.FALSE., info="HxIxI%tab_cav_op")

  CALL Logical_Test(test_opnd, ANY(IxHxI%tab_indexes_mat_op/=[0,1]), test2=.FALSE., info="IxHxI%tab_mat_op")
  CALL Logical_Test(test_opnd, ANY(IxHxI%tab_indexes_cav_op/=[0]  ), test2=.FALSE., info="IxHxI%tab_cav_op")

  CALL Logical_Test(test_opnd, ANY(IxIxH%tab_indexes_mat_op/=[0,0]), test2=.FALSE., info="IxIxH%tab_mat_op")
  CALL Logical_Test(test_opnd, ANY(IxIxH%tab_indexes_cav_op/=[1]  ), test2=.FALSE., info="IxIxH%tab_cav_op")

  CALL Logical_Test(test_opnd, ANY(DipMomtxIxPos%tab_indexes_mat_op/=[4,0]), test2=.FALSE., info="DipMomtxIxPos%tab_mat_op")
  CALL Logical_Test(test_opnd, ANY(DipMomtxIxPos%tab_indexes_cav_op/=[2]  ), test2=.FALSE., info="DipMomtxIxPos%tab_cav_op")

  CALL Logical_Test(test_opnd, ANY(IxDipMomtxPos%tab_indexes_mat_op/=[0,4]), test2=.FALSE., info="IxDipMomtxPos%tab_mat_op")
  CALL Logical_Test(test_opnd, ANY(IxDipMomtxPos%tab_indexes_cav_op/=[2]  ), test2=.FALSE., info="IxDipMomtxPos%tab_cav_op")


  !----------------------------Testing the actions---------------------------
  CALL Get(Mat1w, "w", "Matter", 1)
  CALL Get(Mat2w, "w", "Matter", 2)
  CALL Get(Cavw,  "w", "Cavity", 1)
  CALL Get(Mat1m, "m", "Matter", 1)
  CALL Get(Mat2m, "m", "Matter", 2)
  CALL Get(Mat1lambda, "lambda", "Matter", 1)
  CALL Get(Mat2lambda, "lambda", "Matter", 2)
  CALL Get(Cavlambda,  "lambda", "Cavity", 1)
  CoeffDipMomt1 = ONE ! assumed linear here
  CoeffDipMomt2 = ONE ! assumed linear here
  CALL Get(Mat1Nb, "Nb", "Matter", 1)
  CALL Get(Mat2Nb, "Nb", "Matter", 2)
  CALL Get(CavNb,  "Nb", "Cavity", 1)
  NB = Mat1Nb * Mat2Nb * CavNb

  IF (Debug) THEN
    WRITE(out_unit,*)
    WRITE(out_unit,*) "--- System paramaters"
    WRITE(out_unit,*) "Mat1w         = "//TO_string(Mat1w)
    WRITE(out_unit,*) "Mat2w         = "//TO_string(Mat2w)
    WRITE(out_unit,*) "Mat1m         = "//TO_string(Mat1m)
    WRITE(out_unit,*) "Mat2m         = "//TO_string(Mat2m)
    WRITE(out_unit,*) "Cavw          = "//TO_string(Cavw)
    WRITE(out_unit,*) "Mat1lambda    = "//TO_string(Mat1lambda)
    WRITE(out_unit,*) "Mat2lambda    = "//TO_string(Mat2lambda)
    WRITE(out_unit,*) "Cavlambda     = "//TO_string(Cavlambda)
    WRITE(out_unit,*) "CoeffDipMomt1 = "//TO_string(CoeffDipMomt1)
    WRITE(out_unit,*) "CoeffDipMomt2 = "//TO_string(CoeffDipMomt2)
    WRITE(out_unit,*) "Mat1Nb        = "//TO_string(Mat1Nb)
    WRITE(out_unit,*) "Mat2Nb        = "//TO_string(Mat2Nb)
    WRITE(out_unit,*) "CavNb         = "//TO_string(CavNb)
    WRITE(out_unit,*) "NB            = "//TO_string(NB)
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
    
    Op_phi = ZERO
    CALL Action(Op_phi, DipMomtxIxPos, Phi, Verbose=Verbose, Debug=.FALSE.)
    TotH(:,J) = TotH(:,J) + Op_phi * Cavw * Mat1lambda * Cavlambda * CoeffDipMomt1
    
    Op_phi = ZERO
    CALL Action(Op_phi, IxDipMomtxPos, Phi, Verbose=Verbose, Debug=.FALSE.)
    TotH(:,J) = TotH(:,J) + Op_phi * Cavw * Mat2lambda * Cavlambda * CoeffDipMomt2
  END DO
  IF (Debug) WRITE(out_unit,*)
  IF (Debug) CALL Write_Mat(TotH, out_unit, NB, info="TotH")

  CALL diagonalization(TotH, Eigenenergies, Eigenstates)
  IF (Debug) WRITE(out_unit,*)
  IF (Debug) CALL Write_Vec(Eigenenergies(1:12), out_unit, 1, info="EigenEnergies(1:12)")
  
  CALL Construct_Ana_Normal_modes(Ana_Normal_modes, Mat1w, Mat2w, Mat1m, Mat2m, Cavw, Mat1lambda, Mat2lambda, Cavlambda, CoeffDip&
  &Momt1, CoeffDipMomt2, Debug_opt=.FALSE.)
  IF (Debug) WRITE(out_unit,*)
  IF (Debug) CALL Write_Vec(Ana_Normal_modes, out_unit, 1, info="Ana_Normal_modes")

!----------------------------------------------------------------------------------!
!  /!\            /!\            /!\            /!\            /!\            /!\  !
!   Here we will use the basis vectors in the same order as used in the Action     !
!   procedure of the code :                                                        !
!   \bigl\{\ket{000}, \ket{100}, \ket{200}, \ket{010}, \ket{110}, \ket{210},       !
!         {\ket{001}, \ket{101}, \ket{201}, \ket{011}, \ket{111}, \ket{211}\bigl\}.!
!   EVEN THOUGH the basis set is not the one used on the code /!!!!!\ It is the    !
!   normal coordinates and not the HO basis anymore ! /!\                          !
!  /!\            /!\            /!\            /!\            /!\            /!\  !
!----------------------------------------------------------------------------------!

  ALLOCATE(Ana_Eigenenergies(0:NB-1))
  J = 0
  DO i_3 = 0, CavNb-1
    DO i_2 = 0, Mat2Nb-1
      DO i_1 = 0, Mat1Nb-1
        CALL Compute_Ana_Eigenenergie(Ana_Eigenenergies(J), i_1, i_2, i_3, Ana_Normal_modes(0), Ana_Normal_modes(1)&
        &, Ana_Normal_modes(2), Debug_opt=.FALSE.)
        J = J + 1
      END DO
    END DO
  END DO
  IF (Debug) WRITE(out_unit,*)
  IF (Debug) CALL Write_Vec(Ana_Eigenenergies, out_unit, 1, info="Ana_Eigenenergies")

  DO J = 1, 12
    min_index = MINLOC(Ana_Eigenenergies, dim=1) - 1 ! "-1" because MINLOC does not take into account that the indexes were renamed (?)
    CALL Equal_tensor(error_opnd, EigenEnergies(J), Ana_Eigenenergies(min_index))
    CALL Logical_Test(test_opnd, error_opnd, test2=.FALSE., info="State "//TO_string(J-1))
    IF (Debug .OR. error_opnd) THEN
      WRITE(out_unit,*) "J, EigenEnergies(J), Ana_Eigenenergies(min_index) = "//TO_string(J)//", "//TO_string(EigenEnergies(J))//&
      &", "//TO_string(Ana_Eigenenergies(min_index))
    END IF
    Ana_Eigenenergies(min_index) = HUGE(1)
  END DO


  !----------------------------Testing the writing---------------------------
  CALL Write(IxDipMomtxPos)

  
  !----------------------------Testing the deallocation---------------------------
  CALL Dealloc(IxDipMomtxPos, Verbose=Verbose, Debug=Debug)
  CALL Logical_Test(test_opnd, ALLOCATED(IxDipMomtxPos%tab_indexes_mat_op), test2=.FALSE., info="IxDipMomtxPos%tab_mat deallocated")
  CALL Logical_Test(test_opnd, ALLOCATED(IxDipMomtxPos%tab_indexes_cav_op), test2=.FALSE., info="IxDipMomtxPos%tab_cav deallocated")
  
  CALL Finalize_Test(test_opnd)


  CONTAINS


  SUBROUTINE Construct_Ana_Normal_modes(N_modes, Mat1w_loc, Mat2w_loc, Mat1m_loc, Mat2m_loc, Cavw_loc, Mat1lambda_loc, &
    &Mat2lambda_loc, Cavlambda_loc, CoeffDipMomt1_loc, CoeffDipMomt2_loc, Debug_opt)
    USE QDUtil_m
    IMPLICIT NONE 

    real(kind=Rkind),    intent(inout) :: N_modes(3)
    real(kind=Rkind),    intent(in)    :: Mat1w_loc
    real(kind=Rkind),    intent(in)    :: Mat2w_loc
    real(kind=Rkind),    intent(in)    :: Mat1m_loc
    real(kind=Rkind),    intent(in)    :: Mat2m_loc
    real(kind=Rkind),    intent(in)    :: Cavw_loc
    real(kind=Rkind),    intent(in)    :: Mat1lambda_loc
    real(kind=Rkind),    intent(in)    :: Mat2lambda_loc
    real(kind=Rkind),    intent(in)    :: Cavlambda_loc
    real(kind=Rkind),    intent(in)    :: CoeffDipMomt1_loc
    real(kind=Rkind),    intent(in)    :: CoeffDipMomt2_loc
    logical, optional,   intent(in)    :: Debug_opt

    real(kind=Rkind)                   :: L1, L2 ! analytical coupling terms
    real(kind=Rkind)                   :: MWH(3,3)
    real(kind=Rkind)                   :: N_coos(3,3)
    logical                            :: Debug_local = .FALSE.


    IF (PRESENT(Debug_opt)) THEN; Debug_local = Debug_opt
    ELSE; Debug_local = .FALSE.; END IF 


    L1 = Cavlambda_loc * Mat1lambda_loc * CoeffDipMomt1_loc * Cavw_loc / SQRT(Mat1m_loc)
    L2 = Cavlambda_loc * Mat2lambda_loc * CoeffDipMomt2_loc * Cavw_loc / SQRT(Mat2m_loc)

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "--- System paramaters given to Construct_Ana_Normal_modes"
      WRITE(out_unit,*) "Mat1w_loc         = "//TO_string(Mat1w_loc)
      WRITE(out_unit,*) "Mat2w_loc         = "//TO_string(Mat2w_loc)
      WRITE(out_unit,*) "Mat1m_loc         = "//TO_string(Mat1m_loc)
      WRITE(out_unit,*) "Mat2m_loc         = "//TO_string(Mat2m_loc)
      WRITE(out_unit,*) "Cavw_loc          = "//TO_string(Cavw_loc)
      WRITE(out_unit,*) "Mat1lambda_loc    = "//TO_string(Mat1lambda_loc)
      WRITE(out_unit,*) "Mat2lambda_loc    = "//TO_string(Mat2lambda_loc)
      WRITE(out_unit,*) "Cavlambda_loc     = "//TO_string(Cavlambda_loc)
      WRITE(out_unit,*) "CoeffDipMomt1_loc = "//TO_string(CoeffDipMomt1_loc)
      WRITE(out_unit,*) "CoeffDipMomt2_loc = "//TO_string(CoeffDipMomt2_loc)
      WRITE(out_unit,*) "L1                = "//TO_string(L1)
      WRITE(out_unit,*) "L2                = "//TO_string(L2)
    END IF 

    MWH      = ZERO
    MWH(1,1) = Mat1w**2
    MWH(2,2) = Mat2w**2
    MWH(3,3) = Cavw **2
    MWH(1,3) = L1
    MWH(3,1) = L1
    MWH(2,3) = L2
    MWH(3,2) = L2

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "Analytical MWH :"
      CALL Write_Mat(MWH, out_unit, Size(MWH), info="MWH")
    END IF

    CALL diagonalization(MWH, N_modes, N_coos)
    N_modes = SQRT(N_modes)

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "Normal modes resulting from MWH diagonalization :"
      CALL Write_Vec(N_modes, out_unit, 1, info="Normal modes")
    END IF

  END SUBROUTINE Construct_Ana_Normal_modes


  SUBROUTINE Compute_Ana_Eigenenergie(E, n1, n2, n3, w1, w2, w3, Debug_opt)
    USE QDUtil_m
    IMPLICIT NONE 

    real(kind=Rkind),  intent(inout) :: E
    integer,           intent(in)    :: n1
    integer,           intent(in)    :: n2
    integer,           intent(in)    :: n3
    real(kind=Rkind),  intent(in)    :: w1
    real(kind=Rkind),  intent(in)    :: w2
    real(kind=Rkind),  intent(in)    :: w3
    logical, optional, intent(in)    :: Debug_opt

    logical                          :: Debug_local = .FALSE.


    IF (PRESENT(Debug_opt)) THEN; Debug_local = Debug_opt
    ELSE; Debug_local = .FALSE.; END IF 


    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "--- Arguments given to Compute_Ana_Eigenenergie"
      WRITE(out_unit,*) "n1, n2, n3 = "//TO_string(n1)//", "//TO_string(n2)//", "//TO_string(n3)
      WRITE(out_unit,*) "w1, w2, w3 = "//TO_string(w1)//", "//TO_string(w2)//", "//TO_string(w3)
    END IF 

    E = w1*(REAL(n1, Rkind) + HALF) + w2*(REAL(n2, Rkind) + HALF) + w3*(REAL(n3, Rkind) + HALF)
    IF (Debug_local) WRITE(out_unit,*)
    IF (Debug_local) WRITE(out_unit,*) "E_{|"//TO_string(n1)//TO_string(n2)//TO_string(n3)//">} = "//TO_string(E)

  END SUBROUTINE Compute_Ana_Eigenenergie


END PROGRAM
