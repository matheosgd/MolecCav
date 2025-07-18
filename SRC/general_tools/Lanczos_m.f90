!==================================================================================================
!==================================================================================================
!
! This file is part of MolecCav.
!
!==================================================================================================
!
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
!
!==================================================================================================
!
! README : HAVE NOT BEEN TESTED IN THIS PROGRAM.
! Module that contains all the procedures related to the implementation of the Lanczos iterative d-
! iagonalization methode in MolecCav. In a nutshell : how to construct the Krylov basis, how to au-
! gment it, how to construct the triband Hamiltonian matrix on this basis, and how to change the b-
! asis on which a set a vector is expressed from one to an other (which is not really related to L-
! anczos and may actually be moved to Algebra).
! This module is in the "-6" level.
!
! Initialize : Interface for the MolecCav_Initialize_krylov_basis procedure.
!
! Augment : Interface for the MolecCav_Augment_krylov_basis procedure.
!
! Construct_TribandeH : Interface for the MolecCav_Construct_TribandeH procedure.
!
! BasisChange : Interface for the MolecCav_BasisChange_2TO1_R2 procedure.
!
! MolecCav_Initialize_krylov_basis :
!
! MolecCav_Augment_krylov_basis :
!
! MolecCav_Construct_TribandeH :
!
! MolecCav_BasisChange_2TO1_R2 :
!
!==================================================================================================
!==================================================================================================
MODULE Lanczos_m
  !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
  USE QDUtil_m
  USE Algebra_m
  IMPLICIT NONE
  

  PRIVATE
  
  PUBLIC Initialize, Augment, Construct_TribandeH, BasisChange


  INTERFACE Initialize
    MODULE PROCEDURE MolecCav_Initialize_krylov_basis
  END INTERFACE
  INTERFACE Augment
    MODULE PROCEDURE MolecCav_Augment_krylov_basis
  END INTERFACE
  INTERFACE Construct_TribandeH
    MODULE PROCEDURE MolecCav_Construct_TribandeH
  END INTERFACE
  INTERFACE BasisChange
    MODULE PROCEDURE MolecCav_BasisChange_2TO1_R2
  END INTERFACE


  CONTAINS


  SUBROUTINE MolecCav_Initialize_krylov_basis(KBasis, TotH, Psi, Nb_krylov, Debug)
    USE QDUtil_m
    USE Algebra_m
    USE Sum_of_products_m
    IMPLICIT NONE
  
    real(kind=Rkind), allocatable, intent(inout) :: KBasis(:,:) ! the basis of the krylov basis set expressed in the wavefunction/hamiltonian/complete/original basis set
    TYPE(Sum_of_products_t),       intent(in)    :: TotH        ! the hamiltonian in the original basis set
    real(kind=Rkind),              intent(in)    :: Psi(:)      ! the wavefunction in the original basis set
    integer,                       intent(in)    :: Nb_krylov   ! the SIZE of the krylov basis
    logical, optional,             intent(in)    :: Debug

    real(kind=Rkind), allocatable                :: KBasis_non_ortho(:,:) ! the intermediary KBasis before it is orthonormalized
    real(kind=Rkind), allocatable                :: Intermediary(:,:)     ! for the orthonormality check
    real(kind=Rkind), allocatable                :: Identity(:,:)         ! also
    integer                                      :: i
    logical                                      :: Debug_local
  
    !--- Debugging options --------------------------------
    IF (PRESENT(Debug)) THEN; Debug_local = Debug
    ELSE; Debug_local = .FALSE.; END IF

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "o Arguments of MolecCav_Initialize_krylov_basis :"
      WRITE(out_unit,*) "The <<KBasis>> argument    :" 
      CALL Write_Mat(KBasis, out_unit, SIZE(KBasis, dim=2), info="KBasis")
      WRITE(out_unit,*) "The <<TotH>> argument      :" 
      CALL Write(TotH)
      WRITE(out_unit,*) "The <<Psi>> argument       :" 
      CALL Write_Vec(Psi, out_unit, 1, info="Psi")
      WRITE(out_unit,*) "The <<Nb_krylov>> argument :"//TO_string(Nb_krylov)
      FLUSH(out_unit)
    END IF

    !--- Checking dimensions ------------------------------
    ! IF (.NOT. ALLOCATED(KBasis)) THEN
    !   WRITE(out_unit,*) "### KBasis has to have been allocated BEFORE calling MolecCav_Initialize_krylov_basis. Please check init&
    !   &ialization."
    !   STOP "### KBasis has to have been allocated BEFORE calling MolecCav_Initialize_krylov_basis. Please check initialization."
    ! END IF 

    ! IF (.NOT. ALLOCATED(Psi)) THEN
    !   WRITE(out_unit,*) "### KBasis cannot be constructed in MolecCav_Initialize_krylov_basis if Psi is not allocated. Please che&
    !   &ck initialization."
    !   STOP "### KBasis cannot be constructed in MolecCav_Initialize_krylov_basis if Psi is not allocated. Please check initializa&
    !   &tion."
    ! END IF

    IF (Nb_krylov /= SIZE(KBasis, dim=2)) THEN
      WRITE(out_unit,*) "### The number of vectors KBasis is allocated to does not match Nb_krylov. Please check initialization."
      WRITE(out_unit,*) "    SIZE(KBasis, dim=2) = "//TO_string(Size(KBasis, dim=2))//"; Nb_krylov = "//TO_string(Nb_krylov)
      STOP "### The number of vectors KBasis is allocated to does not match Nb_krylov. Please check initialization."
    END IF 

    IF (SIZE(Psi, dim=1) /= Size(KBasis, dim=1)) THEN
      WRITE(out_unit,*) "### The sizes of the KBasis's vectors do not match Psi's vector SIZE i.e. the original basis set SIZE. P&
      &lease check initialization."
      STOP "### The sizes of the KBasis's vectors do not match Psi's vector SIZE i.e. the original basis set SIZE. Please check i&
      &nitialization."
    END IF 

    !--- Creation of the Krylov basis ---------------------
    ALLOCATE(KBasis_non_ortho(SIZE(Psi), Nb_krylov))

    KBasis_non_ortho(:,:) = ZERO
    KBasis_non_ortho(:,1) = Psi(:)
  
    DO i = 2, Nb_krylov
      CALL Action(KBasis_non_ortho(:,i), TotH, KBasis_non_ortho(:,i-1), Verbose=0, Debug=Debug_local)
    END DO
  
    !--- Orthonormalization -------------------------------
    ALLOCATE(KBasis(SIZE(Psi), Nb_krylov))
    CALL Gram_schmidt(KBasis, KBasis_non_ortho)
    CALL Gram_schmidt(KBasis, KBasis)
   
    !--- Check ortho --------------------------------------
    ALLOCATE(Identity(Nb_krylov,Nb_krylov))
    Identity(:,:) = ZERO

    DO i = 1, Nb_krylov
      Identity(i,i) = 1
    END DO

    ALLOCATE(Intermediary(Nb_krylov,Nb_krylov))
    Intermediary = MATMUL(TRANSPOSE(KBasis), KBasis)
    ! Intermediary = matmul(conjg(transpose(KBasis)), KBasis)

    IF (MAXVAL(ABS(Intermediary-Identity)) > 1E-10_Rkind) THEN
      WRITE(out_unit,*) "################### WARNING ################## WARNING ################## WARNING ##################"
      WRITE(out_unit,*) "                             The KBasis is not orthonormal up to 10E-10 "
      WRITE(out_unit,*) "################### WARNING ################## WARNING ################## WARNING ##################"
    END IF
    
    IF (Debug_local) WRITE(out_unit,*) "--- The constructed KBasis :"
    IF (Debug_local) CALL Write_Mat(KBasis, out_unit, SIZE(KBasis, dim=2), info="KBasis")

  END SUBROUTINE MolecCav_Initialize_krylov_basis

 
  SUBROUTINE MolecCav_Augment_krylov_basis(KBasis, TotH, Debug)
    USE QDUtil_m
    USE Sum_of_products_m
    IMPLICIT NONE
  
    real(kind=Rkind), allocatable, intent(inout) :: KBasis(:,:)
    TYPE(Sum_of_products_t),       intent(in)    :: TotH
    logical, optional,             intent(in)    :: Debug

    real(kind=Rkind), allocatable                :: TemporaryKB(:,:)
    real(kind=Rkind), allocatable                :: V(:)
    real(kind=Rkind), allocatable                :: Intermediary(:,:)
    real(kind=Rkind), allocatable                :: Identity(:,:)
    integer                                      :: Nb_krylov
    integer                                      :: NB ! SIZE of the original basis set
    integer                                      :: i
    logical                                      :: Debug_local

    !--- Debugging options --------------------------------
    IF (PRESENT(Debug)) THEN; Debug_local = Debug
    ELSE; Debug_local = .FALSE.; END IF

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "o Arguments of MolecCav_Augment_krylov_basis :"
      WRITE(out_unit,*) "The <<KBasis>> argument :" 
      IF (ALLOCATED(KBasis)) THEN; CALL Write_Mat(KBasis, out_unit, SIZE(KBasis, dim=2), info="KBasis")
      ELSE; WRITE(out_unit,*) "is NOT allocated !"; END IF ! but why do we need to keep it allocatable ?
      WRITE(out_unit,*) "The <<TotH>> argument :" 
      CALL Write(TotH)
      FLUSH(out_unit)
    END IF

    !--- Initialization -----------------------------------
    NB          = SIZE(KBasis, dim = 1)
    Nb_krylov   = SIZE(KBasis, dim = 2)
    TemporaryKB = KBasis
    ALLOCATE(V(NB))

    DEALLOCATE(KBasis)
    ALLOCATE(KBasis(NB, Nb_krylov+1))

    !--- Construction -----------------------------------
    CALL Action(V, TotH, TemporaryKB(:,Nb_krylov), Debug=Debug) 
    
    KBasis(:,Nb_krylov+1) = V(:)

    DO i = 1, Nb_krylov
      KBasis(:,Nb_krylov+1) = KBasis(:,Nb_krylov+1) - DOT_PRODUCT(TemporaryKB(:,i), V(:))*TemporaryKB(:,i)
    END DO
    KBasis(:,Nb_krylov+1) = KBasis(:,Nb_krylov+1) / SQRT(DOT_PRODUCT(KBasis(:,Nb_krylov+1),KBasis(:,Nb_krylov+1)))

    KBasis(:,1:Nb_krylov) = TemporaryKB
    DEALLOCATE(TemporaryKB, V)

    !--- Check ortho --------------------------------------
    ALLOCATE(Identity(Nb_krylov+1,Nb_krylov+1))
    Identity(:,:) = ZERO
    DO i = 1, Nb_krylov+1
      Identity(i,i) = 1
    END DO

    ALLOCATE(Intermediary(Nb_krylov+1,Nb_krylov+1))
    Intermediary = MATMUL(TRANSPOSE(KBasis), KBasis)
    ! Intermediary = MATMUL(CONJG(TRANSPOSE(KBasis)), KBasis)

    IF (MAXVAL(ABS(Intermediary-Identity)) > 1E-10_Rkind) THEN
      WRITE(out_unit,*) "################### WARNING ################## WARNING ################## WARNING ##################"
      WRITE(out_unit,*) "                             The KBasis is not orthonormal up to 10E-10 "
      WRITE(out_unit,*) "################### WARNING ################## WARNING ################## WARNING ##################"
    END IF
  
    IF (Debug_local) WRITE(out_unit,*) "--- The augmented KBasis :"
    IF (Debug_local) CALL Write_Mat(KBasis, out_unit, SIZE(KBasis, dim=2), info="KBasis")

  END SUBROUTINE MolecCav_Augment_krylov_basis


  SUBROUTINE MolecCav_Construct_TribandeH(TribandH, KBasis, TotH, Debug)
    USE QDUtil_m
    USE Sum_of_products_m
    IMPLICIT NONE

    real(kind=Rkind),        intent(inout) :: TribandH(:,:)
    real(kind=Rkind),        intent(in)    :: KBasis(:,:)
    TYPE(Sum_of_products_t), intent(in)    :: TotH
    logical, optional,       intent(in)    :: Debug

    integer                                :: i, Nb_krylov
    logical                                :: Debug_local

    !--- Debugging options --------------------------------
    IF (PRESENT(Debug)) THEN; Debug_local = Debug
    ELSE; Debug_local = .FALSE.; END IF

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "o Arguments of MolecCav_Construct_TribandeH :"
      WRITE(out_unit,*) "The <<TribandH>> argument :" 
      CALL Write_Mat(TribandH, out_unit, SIZE(TribandH, dim=2), info="TribandH")
      WRITE(out_unit,*) "The <<KBasis>> argument :" 
      CALL Write_Mat(KBasis, out_unit, SIZE(KBasis, dim=2), info="KBasis")
      WRITE(out_unit,*) "The <<TotH>> argument :" 
      CALL Write(TotH)
      FLUSH(out_unit)
    END IF

    !--- Checking dimensions ------------------------------
    ! IF (.NOT. ALLOCATED(TribandH)) THEN
    !   WRITE(out_unit,*) "### TribandH has to have been allocated BEFORE calling MolecCav_Construct_TribandeH. Please check initia&
    !   &lization."
    !   STOP "### TribandH has to have been allocated BEFORE calling MolecCav_Construct_TribandeH. Please check initialization."
    ! END IF 

    ! IF (.NOT. ALLOCATED(KBasis)) THEN
    !   WRITE(out_unit,*) "### TribandH cannot be constructed in MolecCav_Construct_TribandeH if KBasis is not allocated. Please ch&
    !   &eck initialization."
    !   STOP "### TribandH cannot be constructed in MolecCav_Construct_TribandeH if KBasis is not allocated. Please check initializ&
    !   &ation."
    ! END IF

    IF (SIZE(KBasis, dim=1) /= SIZE(TribandH, dim=1)) THEN
      WRITE(out_unit,*) "### The sizes of the KBasis's vectors do not match TribandH's vector size i.e. the original basis set si&
      &ze. Please check initialization."
      STOP "### The sizes of the KBasis's vectors do not match TribandH's vector size i.e. the original basis set size. Please ch&
      &eck initialization."
    END IF 

    IF (SIZE(KBasis, dim=2) /= SIZE(TribandH, dim=2)) THEN
      WRITE(out_unit,*) "### The number of vectors KBasis is allocated to does not match TribandH. Please check initialization."
      WRITE(out_unit,*) "    size(KBasis, dim=2) = "//TO_string(SIZE(KBasis, dim=2))//"; size(TribandH, dim=2) = "//TO_string(SIZ&
      &E(TribandH, dim=2))
      STOP "### The number of vectors KBasis is allocated to does not match TribandH. Please check initialization."
    END IF 

    !--- Allocation ---------------------------------------
    Nb_krylov = SIZE(KBasis, dim = 2)    
   
    !--- Construction -------------------------------------
    DO i = 1, Nb_krylov
      CALL Action(TribandH(:,i), TotH, KBasis(:,i), Debug=Debug_local)
    END DO 

    IF (Debug_local) WRITE(out_unit,*) "--- The constructed TribandH :"
    IF (Debug_local) CALL Write_Mat(TribandH, out_unit, SIZE(TribandH, dim=2), info="TribandH")

  END SUBROUTINE MolecCav_Construct_TribandeH


  SUBROUTINE MolecCav_BasisChange_2TO1_R2(Psi_1, Psi_2, ChangeM_1TO2, Debug)
    USE QDUtil_m
    IMPLICIT NONE

    real(kind=Rkind),  intent(inout) :: Psi_1(:,:)        ! vectors defined on the basis B_1
    real(kind=Rkind),  intent(in)    :: Psi_2(:,:)        ! vectors defined on the basis B_2
    real(kind=Rkind),  intent(in)    :: ChangeM_1TO2(:,:) ! change-of-basis matrixfrom Basis B_1 to B_2 i.e. B_2 expressed in B_1
    logical, optional, intent(in)    :: Debug

    integer                          :: i, j, Nb_2
    logical                          :: Debug_local

    !--- Debugging options --------------------------------
    IF (PRESENT(Debug)) THEN; Debug_local = Debug
    ELSE; Debug_local = .FALSE.; END IF

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "o Arguments of MolecCav_BasisChange_2TO1_R2 :"
      WRITE(out_unit,*) "The <<Psi_1>> argument :" 
      CALL Write_Mat(Psi_1, out_unit, SIZE(Psi_1, dim=2), info="Psi_1")
      WRITE(out_unit,*) "The <<Psi_2>> argument :" 
      CALL Write_Mat(Psi_2, out_unit, SIZE(Psi_2, dim=2), info="Psi_2")
      WRITE(out_unit,*) "The <<ChangeM_1TO2>> argument :" 
      CALL Write_Mat(ChangeM_1TO2, out_unit, SIZE(ChangeM_1TO2, dim=2), info="ChangeM_1TO2")
      FLUSH(out_unit)
    END IF

    !--- Checking dimensions ------------------------------
    ! IF (.NOT. ALLOCATED(Psi_1)) THEN
    !   WRITE(out_unit,*) "### Psi_1 has to have been allocated BEFORE calling MolecCav_BasisChange_2TO1_R2. Please check initializ&
    !   &ation."
    !   STOP "### Psi_1 has to have been allocated BEFORE calling MolecCav_BasisChange_2TO1_R2. Please check initialization."
    ! END IF 

    ! IF (.NOT. ALLOCATED(Psi_2)) THEN
    !   WRITE(out_unit,*) "### Psi_1 cannot be constructed in MolecCav_BasisChange_2TO1_R2 if Psi_2 is not allocated. Please check &
    !   &initialization."
    !   STOP "### Psi_1 cannot be constructed in MolecCav_BasisChange_2TO1_R2 if Psi_2 is not allocated. Please check initializatio&
    !   &n."
    ! END IF

    ! IF (.NOT. ALLOCATED(ChangeM_1TO2)) THEN
    !   WRITE(out_unit,*) "### Psi_1 cannot be constructed in MolecCav_BasisChange_2TO1_R2 if ChangeM_1TO2 is not allocated. Please&
    !   & check initialization."
    !   STOP "### Psi_1 cannot be constructed in MolecCav_BasisChange_2TO1_R2 if ChangeM_1TO2 is not allocated. Please check initia&
    !   &lization."
    ! END IF

    IF (SIZE(Psi_1, dim=1) /= SIZE(ChangeM_1TO2, dim=1)) THEN
      WRITE(out_unit,*) "### The sizes of the Psi_1's vectors do not match ChangeM_1TO2's vector size. Please check initialization."
      STOP "### The sizes of the Psi_1's vectors do not match Psi_2's vector size. Please check initialization."
    END IF 

    IF (SIZE(Psi_1, dim=2) /= SIZE(Psi_2, dim=2)) THEN 
      WRITE(out_unit,*) "### The number of vectors held in Psi_1 does not match the Psi_2 one. Please check initialization."
      WRITE(out_unit,*) "    SIZE(Psi_1, dim=2) = "//TO_string(SIZE(Psi_1, dim=2))//"; SIZE(Psi_2, dim=2) = "//TO_string(SIZ&
      &E(Psi_2, dim=2))
      STOP "### The number of vectors held in Psi_1 does not match the Psi_2 one. Please check initialization."
    END IF 

    IF (SIZE(Psi_2, dim=1) /= SIZE(ChangeM_1TO2, dim=2)) THEN
      WRITE(out_unit,*) "### The sizes of the Psi_1's vectors do not match ChangeM_1TO2's number of vectors. Please check initial&
      &ization."
      STOP "### The sizes of the Psi_1's vectors do not match ChangeM_1TO2's number of vectors. Please check initialization."
    END IF 

    !--- Change of basis ----------------------------------
    Psi_1(:,:) = ZERO
    Psi_1 = MATMUL(ChangeM_1TO2, Psi_2)

    IF (Debug_local) WRITE(out_unit,*) "--- The constructed Psi_1 :"
    IF (Debug_local) CALL Write_Mat(Psi_1, out_unit, SIZE(Psi_1, dim=2), info="Psi_1")

  END SUBROUTINE MolecCav_BasisChange_2TO1_R2


END MODULE