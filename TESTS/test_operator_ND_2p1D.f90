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
  logical             :: Debug   = .FALSE.

  logical             :: Dense   = .FALSE.
  TYPE(Operator_ND_t) :: HxIxI
  TYPE(Operator_ND_t) :: IxHxI
  TYPE(Operator_ND_t) :: IxIxH
  TYPE(Operator_ND_t) :: DipMomtxIxPos
  TYPE(Operator_ND_t) :: IxDipMomtxPos
  real(kind=Rkind)    :: Matw
  real(kind=Rkind)    :: Matm
  real(kind=Rkind)    :: Cavw
  real(kind=Rkind)    :: Cavlambda
  real(kind=Rkind)    :: CoeffDipMomt

  real(kind=Rkind)    ::    Phi(12)                                                                             ! an any vector representing the excitation state/wavefunction = a linear combination of the basis functions from the canonical basis set over \mathbb{R} /!\ Not normalized yet !
  real(kind=Rkind)    :: Op_phi(12)                                                                             ! an any vector representing the excitation state/wavefunction = a linear combination of the basis functions from the canonical basis set over \mathbb{R} /!\ Not normalized yet !
  real(kind=Rkind)    :: TotH(12,12)                                                                             ! an any vector representing the excitation state/wavefunction = a linear combination of the basis functions from the canonical basis set over \mathbb{R} /!\ Not normalized yet !
  real(kind=Rkind)    :: REigval(12)                                                                             ! an any vector representing the excitation state/wavefunction = a linear combination of the basis functions from the canonical basis set over \mathbb{R} /!\ Not normalized yet !
  real(kind=Rkind)    :: REigvec(12,12)                                                                             ! an any vector representing the excitation state/wavefunction = a linear combination of the basis functions from the canonical basis set over \mathbb{R} /!\ Not normalized yet !

  TYPE(test_t)        :: test_opnd
  logical             :: error_opnd = .FALSE.

  integer             :: J


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
  CALL Get(Matw, "w", "Matter", 1)
  CALL Get(Matm, "m", "Matter", 1)
  CALL Get(Cavw, "w", "Cavity", 1)
  CALL Get(Cavlambda, "lambda", "Cavity", 1)
  CoeffDipMomt = ONE ! assumed linear here
  IF (Debug) THEN
    WRITE(out_unit,*)
    WRITE(out_unit,*) "--- System paramaters"
    WRITE(out_unit,*) "Matw         = "//TO_string(Matw)
    WRITE(out_unit,*) "Matm         = "//TO_string(Matm)
    WRITE(out_unit,*) "Cavw         = "//TO_string(Cavw)
    WRITE(out_unit,*) "Cavlambda    = "//TO_string(Cavlambda)
    WRITE(out_unit,*) "CoeffDipMomt = "//TO_string(CoeffDipMomt)
  END IF 

  TotH = ZERO
  DO J = 1, 12
    Phi = ZERO
    Phi(J) = ONE

    Op_phi = ZERO
    CALL Action(Op_phi, HxIxI, Phi, Verbose=Verbose, Debug=Debug)
    TotH(:,J) = TotH(:,J) + Op_phi

    Op_phi = ZERO
    CALL Action(Op_phi, IxHxI, Phi, Verbose=Verbose, Debug=Debug)
    TotH(:,J) = TotH(:,J) + Op_phi
    
    Op_phi = ZERO
    CALL Action(Op_phi, IxIxH, Phi, Verbose=Verbose, Debug=Debug)
    TotH(:,J) = TotH(:,J) + Op_phi
    
    Op_phi = ZERO
    CALL Action(Op_phi, DipMomtxIxPos, Phi, Verbose=Verbose, Debug=Debug)
    TotH(:,J) = TotH(:,J) + Op_phi
    
    Op_phi = ZERO
    CALL Action(Op_phi, IxDipMomtxPos, Phi, Verbose=Verbose, Debug=Debug)
    TotH(:,J) = TotH(:,J) + Op_phi
  END DO

  CALL Write_Mat(TotH, out_unit, 12, info="TotH")

  !----------------------------Testing the writing---------------------------
  CALL Write(IxDipMomtxPos)

  
  !----------------------------Testing the deallocation---------------------------
  CALL Dealloc(IxDipMomtxPos, Verbose=Verbose, Debug=Debug)
  CALL Logical_Test(test_opnd, ALLOCATED(IxDipMomtxPos%tab_indexes_mat_op), test2=.FALSE., info="IxDipMomtxPos%tab_mat deallocated")
  CALL Logical_Test(test_opnd, ALLOCATED(IxDipMomtxPos%tab_indexes_cav_op), test2=.FALSE., info="IxDipMomtxPos%tab_cav deallocated")
  
  CALL Finalize_Test(test_opnd)


  CONTAINS


  SUBROUTINE Construct_ana_TotH_psi_real(Ana_TotH_psi_R1_real, Matw_loc, Matm_loc, Cavw_loc, Cavlambda_loc, CoeffDipMomt_loc, Psi&
    &, Debug_opt)
    USE QDUtil_m
    IMPLICIT NONE 
!----------------------------------------------------------------------------------!
!  /!\            /!\            /!\            /!\            /!\            /!\  !
!   The following analytical form of the WF resulting from the action of the       !
!   H_{tot}^{1p1D} on the R1 \Psi^{1p1D} is done assuming the Psi is written on    !
!   the same basis set than inside the code i.e. the eigenvectors of the           !
!   H_{tot}^{1p1D} are set in the same order as used in the Action procedure :     !
!   \bigl\{\ket{00}, \ket{10}, \ket{20}, \ket{01}, \ket{11}, \ket{21}\bigl\}.      !
!   This has to be taken into account as writing the analytical hamiltonian matrix !
!   and as applying it to the WF vector (pay attention to the coefficients order)  !
!  /!\            /!\            /!\            /!\            /!\            /!\  !
!----------------------------------------------------------------------------------!

    real(kind=Rkind),    intent(inout) :: Ana_TotH_psi_R1_real(12)
    real(kind=Rkind),    intent(in)    :: Matw_loc
    real(kind=Rkind),    intent(in)    :: Matm_loc
    real(kind=Rkind),    intent(in)    :: Cavw_loc
    real(kind=Rkind),    intent(in)    :: Cavlambda_loc
    real(kind=Rkind),    intent(in)    :: CoeffDipMomt_loc
    real(kind=Rkind),    intent(in)    :: Psi(12)
    logical, optional,   intent(in)    :: Debug_opt

    real(kind=Rkind)                   :: A                     ! analytical coupling term
    logical                            :: Debug_local = .FALSE.

    IF (PRESENT(Debug_opt)) THEN; Debug_local = Debug_opt
    ELSE; Debug_local = .FALSE.; END IF 

    A = Cavlambda_loc * CoeffDipMomt_loc * SQRT(Cavw_loc / (Matw_loc*Matm_loc))

    Ana_TotH_psi_R1_real(1) =   (  Matw_loc +   Cavw_loc)*Psi(1) + A * Psi(5)
    Ana_TotH_psi_R1_real(2) =   (3*Matw_loc +   Cavw_loc)*Psi(2) + A * ( Psi(4) + SQRT(TWO)*Psi(12) )
    Ana_TotH_psi_R1_real(3) =   (5*Matw_loc +   Cavw_loc)*Psi(3) + A * SQRT(TWO) * Psi(5)
    Ana_TotH_psi_R1_real(4) =   (  Matw_loc + 3*Cavw_loc)*Psi(4) + A * Psi(2)
    Ana_TotH_psi_R1_real(5) = 3*(  Matw_loc +   Cavw_loc)*Psi(5) + A * ( Psi(1) + SQRT(TWO)*Psi(3) )
    Ana_TotH_psi_R1_real(12) =   (5*Matw_loc + 3*Cavw_loc)*Psi(12) + A * SQRT(TWO) * Psi(2)
    
    Ana_TotH_psi_R1_real = Ana_TotH_psi_R1_real / 2

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "Resulting matrix of the TotH action over the operand :"
      CALL Write_Vec(Ana_TotH_psi_R1_real, out_unit, Size(Ana_TotH_psi_R1_real), info="Ana_TotH_psi_R1_real")
    END IF

  END SUBROUTINE


  SUBROUTINE Construct_ana_TotH_psi_complex(Ana_TotH_psi_R1_complex, Matw_loc, Matm_loc, Cavw_loc, Cavlambda_loc, CoeffDipMomt_lo&
    &c, Psi, Debug_opt)
    USE QDUtil_m
    IMPLICIT NONE 
!----------------------------------------------------------------------------------!
!  /!\            /!\            /!\            /!\            /!\            /!\  !
!   The following analytical form of the WF resulting from the action of the       !
!   H_{tot}^{1p1D} on the R1 \Psi^{1p1D} is done assuming the Psi is written on    !
!   the same basis set than inside the code i.e. the eigenvectors of the           !
!   H_{tot}^{1p1D} are set in the same order as used in the Action procedure :     !
!   \bigl\{\ket{00}, \ket{10}, \ket{20}, \ket{01}, \ket{11}, \ket{21}\bigl\}.      !
!   This has to be taken into account as writing the analytical hamiltonian matrix !
!   and as applying it to the WF vector (pay attention to the coefficients order)  !
!  /!\            /!\            /!\            /!\            /!\            /!\  !
!----------------------------------------------------------------------------------!

    complex(kind=Rkind), intent(inout) :: Ana_TotH_psi_R1_complex(12)
    real(kind=Rkind),    intent(in)    :: Matw_loc
    real(kind=Rkind),    intent(in)    :: Matm_loc
    real(kind=Rkind),    intent(in)    :: Cavw_loc
    real(kind=Rkind),    intent(in)    :: Cavlambda_loc
    real(kind=Rkind),    intent(in)    :: CoeffDipMomt_loc
    complex(kind=Rkind), intent(in)    :: Psi(12)
    logical, optional,   intent(in)    :: Debug_opt

    real(kind=Rkind)                   :: A ! analytical coupling term
    logical                            :: Debug_local = .FALSE.

    IF (PRESENT(Debug_opt)) THEN; Debug_local = Debug_opt
    ELSE; Debug_local = .FALSE.; END IF 

    A = Cavlambda_loc * CoeffDipMomt_loc * SQRT(Cavw_loc / (Matw_loc*Matm_loc))
    IF (Debug_local) WRITE(out_unit,*) "A = "//TO_string(A)

    Ana_TotH_psi_R1_complex(1) =   (  Matw_loc +   Cavw_loc)*Psi(1) + A * Psi(5)
    Ana_TotH_psi_R1_complex(2) =   (3*Matw_loc +   Cavw_loc)*Psi(2) + A * ( Psi(4) + SQRT(TWO)*Psi(12) )
    Ana_TotH_psi_R1_complex(3) =   (5*Matw_loc +   Cavw_loc)*Psi(3) + A * SQRT(TWO) * Psi(5)
    Ana_TotH_psi_R1_complex(4) =   (  Matw_loc + 3*Cavw_loc)*Psi(4) + A * Psi(2)
    Ana_TotH_psi_R1_complex(5) = 3*(  Matw_loc +   Cavw_loc)*Psi(5) + A * ( Psi(1) + SQRT(TWO)*Psi(3) )
    Ana_TotH_psi_R1_complex(12) =   (5*Matw_loc + 3*Cavw_loc)*Psi(12) + A * SQRT(TWO) * Psi(2)
    
    Ana_TotH_psi_R1_complex = Ana_TotH_psi_R1_complex / 2

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "Resulting matrix of the TotH action over the operand :"
      CALL Write_Vec(Ana_TotH_psi_R1_complex, out_unit, Size(Ana_TotH_psi_R1_complex), info="Ana_TotH_psi_R1_complex")
    END IF

  END SUBROUTINE


END PROGRAM
