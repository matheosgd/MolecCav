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
PROGRAM test_operator_ND_3p0D
  !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real124
  USE QDUtil_m
  USE Tests_m
  USE Algebra_m
  USE Operator_ND_m
  IMPLICIT NONE


  integer             :: Verbose = 50
  logical             :: Debug   = .TRUE.

  logical             :: Dense   = .FALSE.
  TYPE(Operator_ND_t) :: HxIxI
  TYPE(Operator_ND_t) :: IxHxI
  TYPE(Operator_ND_t) :: IxIxH
  TYPE(Operator_ND_t) :: DipMomtxPosxI
  TYPE(Operator_ND_t) :: DipMomtxIxPos
  real(kind=Rkind)    :: Matw
  real(kind=Rkind)    :: Matm
  real(kind=Rkind)    :: Cav1w
  real(kind=Rkind)    :: Cav2w
  real(kind=Rkind)    :: Matlambda
  real(kind=Rkind)    :: Cav1lambda
  real(kind=Rkind)    :: Cav2lambda
  real(kind=Rkind)    :: CoeffDipMomt
  integer             :: MatNb
  integer             :: Cav1Nb
  integer             :: Cav2Nb
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
  CALL Initialize_Test(test_opnd, test_name="OUT/test_file_opnd_1p2D")


  !-------------------------Operator_ND object initialization-------------------------
  CALL Initialize(HxIxI, "Hamiltonian ", "identity, identity", in_unit, Dense=Dense, Verbose=Verbose, Debug=.FALSE.)
  IF (Debug) THEN
    WRITE(out_unit,*)
    WRITE(out_unit,*) "--------------Operator_ND object constructed by MolecCav_Initialize_operator_ND--------------"
    CALL Write(HxIxI)
    WRITE(out_unit,*) "------------End Operator_ND object constructed by MolecCav_Initialize_operator_ND------------"
  END IF

  CALL Initialize(IxHxI, "identity ", "Hamiltonian, identity", in_unit, Dense=Dense, Verbose=Verbose, Debug=.FALSE.)
  IF (Debug) THEN
    WRITE(out_unit,*)
    WRITE(out_unit,*) "--------------Operator_ND object constructed by MolecCav_Initialize_operator_ND--------------"
    CALL Write(IxHxI)
    WRITE(out_unit,*) "------------End Operator_ND object constructed by MolecCav_Initialize_operator_ND------------"
  END IF

  CALL Initialize(IxIxH, "identity ", "identity, Hamiltonian", in_unit, Dense=Dense, Verbose=Verbose, Debug=.FALSE.)
  IF (Debug) THEN
    WRITE(out_unit,*)
    WRITE(out_unit,*) "--------------Operator_ND object constructed by MolecCav_Initialize_operator_ND--------------"
    CALL Write(IxIxH)
    WRITE(out_unit,*) "------------End Operator_ND object constructed by MolecCav_Initialize_operator_ND------------"
  END IF

  CALL Initialize(DipMomtxPosxI, "DipMomt", "position, identity", in_unit, Dense=Dense, Verbose=Verbose, Debug=.FALSE.)
  IF (Debug) THEN
    WRITE(out_unit,*)
    WRITE(out_unit,*) "--------------Operator_ND object constructed by MolecCav_Initialize_operator_ND--------------"
    CALL Write(DipMomtxPosxI)
    WRITE(out_unit,*) "------------End Operator_ND object constructed by MolecCav_Initialize_operator_ND------------"
  END IF

  CALL Initialize(DipMomtxIxPos, "DipMomt", "identity, Position ", in_unit, Dense=Dense, Verbose=Verbose, Debug=.FALSE.)
  IF (Debug) THEN
    WRITE(out_unit,*)
    WRITE(out_unit,*) "--------------Operator_ND object constructed by MolecCav_Initialize_operator_ND--------------"
    CALL Write(DipMomtxIxPos)
    WRITE(out_unit,*) "------------End Operator_ND object constructed by MolecCav_Initialize_operator_ND------------"
  END IF


  !----------------------------Testing the initialization---------------------------
  CALL Logical_Test(test_opnd, ANY(HxIxI%tab_indexes_mat_op/=[1]  ), test2=.FALSE., info="HxIxI%tab_mat_op")
  CALL Logical_Test(test_opnd, ANY(HxIxI%tab_indexes_cav_op/=[0,0]), test2=.FALSE., info="HxIxI%tab_cav_op")

  CALL Logical_Test(test_opnd, ANY(IxHxI%tab_indexes_mat_op/=[0]  ), test2=.FALSE., info="IxHxI%tab_mat_op")
  CALL Logical_Test(test_opnd, ANY(IxHxI%tab_indexes_cav_op/=[1,0]), test2=.FALSE., info="IxHxI%tab_cav_op")

  CALL Logical_Test(test_opnd, ANY(IxIxH%tab_indexes_mat_op/=[0]  ), test2=.FALSE., info="IxIxH%tab_mat_op")
  CALL Logical_Test(test_opnd, ANY(IxIxH%tab_indexes_cav_op/=[0,1]), test2=.FALSE., info="IxIxH%tab_cav_op")

  CALL Logical_Test(test_opnd, ANY(DipMomtxPosxI%tab_indexes_mat_op/=[4]  ), test2=.FALSE., info="DipMomtxPosxI%tab_mat_op")
  CALL Logical_Test(test_opnd, ANY(DipMomtxPosxI%tab_indexes_cav_op/=[2,0]), test2=.FALSE., info="DipMomtxPosxI%tab_cav_op")

  CALL Logical_Test(test_opnd, ANY(DipMomtxIxPos%tab_indexes_mat_op/=[4]  ), test2=.FALSE., info="DipMomtxIxPos%tab_mat_op")
  CALL Logical_Test(test_opnd, ANY(DipMomtxIxPos%tab_indexes_cav_op/=[0,2]), test2=.FALSE., info="DipMomtxIxPos%tab_cav_op")


  !----------------------------Testing the actions---------------------------
  CALL Get(Matw,  "w", "Matter", 1)
  CALL Get(Cav1w, "w", "Cavity", 1)
  CALL Get(Cav2w, "w", "Cavity", 2)
  CALL Get(Matm,  "m", "Matter", 1)
  CALL Get(Matlambda,  "lambda", "Matter", 1)
  CALL Get(Cav1lambda, "lambda", "Cavity", 1)
  CALL Get(Cav2lambda, "lambda", "Cavity", 2)
  CoeffDipMomt = ONE ! assumed linear here
  CALL Get(MatNb,  "Nb", "Matter", 1)
  CALL Get(Cav1Nb, "Nb", "Cavity", 1)
  CALL Get(Cav2Nb, "Nb", "Cavity", 2)
  NB = MatNb * Cav1Nb * Cav2Nb

  IF (Debug) THEN
    WRITE(out_unit,*)
    WRITE(out_unit,*) "--- System paramaters"
    WRITE(out_unit,*) "Matw         = "//TO_string(Matw)
    WRITE(out_unit,*) "Matm         = "//TO_string(Matm)
    WRITE(out_unit,*) "Cav1w        = "//TO_string(Cav1w)
    WRITE(out_unit,*) "Cav2w        = "//TO_string(Cav2w)
    WRITE(out_unit,*) "Matlambda    = "//TO_string(Matlambda)
    WRITE(out_unit,*) "Cav1lambda   = "//TO_string(Cav1lambda)
    WRITE(out_unit,*) "Cav2lambda   = "//TO_string(Cav2lambda)
    WRITE(out_unit,*) "CoeffDipMomt = "//TO_string(CoeffDipMomt)
    WRITE(out_unit,*) "MatNb        = "//TO_string(MatNb)
    WRITE(out_unit,*) "Cav1Nb       = "//TO_string(Cav1Nb)
    WRITE(out_unit,*) "Cav2Nb       = "//TO_string(Cav2Nb)
    WRITE(out_unit,*) "NB           = "//TO_string(NB)
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
    CALL Action(Op_phi, DipMomtxPosxI, Phi, Verbose=Verbose, Debug=.FALSE.)
    TotH(:,J) = TotH(:,J) + Op_phi * Cav1w * Matlambda * Cav1lambda * CoeffDipMomt
    
    Op_phi = ZERO
    CALL Action(Op_phi, DipMomtxIxPos, Phi, Verbose=Verbose, Debug=.FALSE.)
    TotH(:,J) = TotH(:,J) + Op_phi * Cav2w * Matlambda * Cav2lambda * CoeffDipMomt
  END DO
  IF (Debug) WRITE(out_unit,*)
  IF (Debug) CALL Write_Mat(TotH, out_unit, NB, info="TotH")

  CALL diagonalization(TotH, Eigenenergies, Eigenstates)
  IF (Debug) WRITE(out_unit,*)
  IF (Debug) CALL Write_Vec(Eigenenergies(1:4), out_unit, 1, info="EigenEnergies(1:4)")
  
  CALL Construct_Ana_Normal_modes(Ana_Normal_modes, Matw, Matm, Cav1w, Cav2w, Matlambda, Cav1lambda, Cav2lambda, CoeffDipMomt, De&
  &bug_opt=.TRUE.)
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
  DO i_3 = 0, Cav2Nb-1
    DO i_2 = 0, Cav1Nb-1
      DO i_1 = 0, MatNb-1
        CALL Compute_Ana_Eigenenergie(Ana_Eigenenergies(J), i_1, i_2, i_3, Ana_Normal_modes(0), Ana_Normal_modes(1)&
        &, Ana_Normal_modes(2), Debug_opt=.FALSE.)
        J = J + 1
      END DO
    END DO
  END DO
  IF (Debug) WRITE(out_unit,*)
  IF (Debug) CALL Write_Vec(Ana_Eigenenergies, out_unit, 1, info="Ana_Eigenenergies")

  DO J = 1, 4
    min_index = MINLOC(Ana_Eigenenergies, dim=1) - 1 ! "-1" because MINLOC does not take into account that the indexes were renamed (?)
      ! WRITE(out_unit,*) "MININDEX = "//TO_string(min_index)
    CALL Equal_tensor(error_opnd, EigenEnergies(J), Ana_Eigenenergies(min_index))
    CALL Logical_Test(test_opnd, error_opnd, test2=.FALSE., info="State "//TO_string(J-1))
    IF (Debug .OR. error_opnd) THEN
      WRITE(out_unit,*) "J, EigenEnergies(J), Ana_Eigenenergies(min_index) = "//TO_string(J)//", "//TO_string(EigenEnergies(J))//&
      &", "//TO_string(Ana_Eigenenergies(min_index))
    END IF
    Ana_Eigenenergies(min_index) = HUGE(1)
  END DO


  !----------------------------Testing the writing---------------------------
  CALL Write(DipMomtxPosxI)

  
  !----------------------------Testing the deallocation---------------------------
  CALL Dealloc(DipMomtxPosxI, Verbose=Verbose, Debug=Debug)
  CALL Logical_Test(test_opnd, ALLOCATED(DipMomtxPosxI%tab_indexes_mat_op), test2=.FALSE., info="DipMomtxPosxI%tab_mat deallocated")
  CALL Logical_Test(test_opnd, ALLOCATED(DipMomtxPosxI%tab_indexes_cav_op), test2=.FALSE., info="DipMomtxPosxI%tab_cav deallocated")
  
  CALL Finalize_Test(test_opnd)


  CONTAINS


  SUBROUTINE Construct_Ana_Normal_modes(N_modes, Matw_loc, Matm_loc, Cav1w_loc, Cav2w_loc, Matlambda_loc, Cav1lambda_loc, &
    &Cav2lambda_loc, CoeffDipMomt_loc, Debug_opt)
    USE QDUtil_m
    IMPLICIT NONE 

    real(kind=Rkind),    intent(inout) :: N_modes(3)
    real(kind=Rkind),    intent(in)    :: Matw_loc
    real(kind=Rkind),    intent(in)    :: Matm_loc
    real(kind=Rkind),    intent(in)    :: Cav1w_loc
    real(kind=Rkind),    intent(in)    :: Cav2w_loc
    real(kind=Rkind),    intent(in)    :: Matlambda_loc
    real(kind=Rkind),    intent(in)    :: Cav1lambda_loc
    real(kind=Rkind),    intent(in)    :: Cav2lambda_loc
    real(kind=Rkind),    intent(in)    :: CoeffDipMomt_loc
    logical, optional,   intent(in)    :: Debug_opt

    real(kind=Rkind)                   :: L1, L2 ! analytical coupling terms
    real(kind=Rkind)                   :: MWH(3,3)
    real(kind=Rkind)                   :: N_coos(3,3)
    logical                            :: Debug_local = .FALSE.


    IF (PRESENT(Debug_opt)) THEN; Debug_local = Debug_opt
    ELSE; Debug_local = .FALSE.; END IF 


    L1 = Cav1lambda_loc * Matlambda_loc * CoeffDipMomt_loc * Cav1w_loc / SQRT(Matm_loc)
    L2 = Cav2lambda_loc * Matlambda_loc * CoeffDipMomt_loc * Cav2w_loc / SQRT(Matm_loc)

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "--- System paramaters given to Construct_Ana_Normal_modes"
      WRITE(out_unit,*) "Matw_loc         = "//TO_string(Matw_loc)
      WRITE(out_unit,*) "Matm_loc         = "//TO_string(Matm_loc)
      WRITE(out_unit,*) "Cav1w_loc        = "//TO_string(Cav1w_loc)
      WRITE(out_unit,*) "Cav2w_loc        = "//TO_string(Cav2w_loc)
      WRITE(out_unit,*) "Matlambda_loc    = "//TO_string(Matlambda_loc)
      WRITE(out_unit,*) "Cav1lambda_loc   = "//TO_string(Cav1lambda_loc)
      WRITE(out_unit,*) "Cav2lambda_loc   = "//TO_string(Cav2lambda_loc)
      WRITE(out_unit,*) "CoeffDipMomt_loc = "//TO_string(CoeffDipMomt_loc)
      WRITE(out_unit,*) "L1               = "//TO_string(L1)
      WRITE(out_unit,*) "L2               = "//TO_string(L2)
    END IF 

    MWH      = ZERO
    MWH(1,1) = Matw_loc **2
    MWH(2,2) = Cav1w_loc**2
    MWH(3,3) = Cav2w_loc**2
    MWH(1,2) = L1
    MWH(2,1) = L1
    MWH(1,3) = L2
    MWH(3,1) = L2

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
