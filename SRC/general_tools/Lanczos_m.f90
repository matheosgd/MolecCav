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
! README :
! to be written soon
!==================================================================================================
!==================================================================================================
MODULE Lanczos_m
  !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
  USE QDUtil_m
  IMPLICIT NONE
  

  PRIVATE
  
  PUBLIC Initialize, Increment, Gram_schmidt


  INTERFACE Initialize
    MODULE PROCEDURE MolecCav_Initialize_krylov_basis
  END INTERFACE
  INTERFACE Augment
    MODULE PROCEDURE MolecCav_Augment_krylov_basis
  END INTERFACE
  INTERFACE Gram_schmidt
    MODULE PROCEDURE MolecCav_Gram_schmidt
  END INTERFACE
  INTERFACE MolecCav_Construct_H_tribande
    MODULE PROCEDURE MolecCav_Construct_H_tribande
  END INTERFACE
  INTERFACE KrylovTOBasis
    MODULE PROCEDURE MolecCav_KrylovTOBasis
  END INTERFACE
  INTERFACE MolecCav_Construct_psi_approx
    MODULE PROCEDURE MolecCav_Construct_psi_approx
  END INTERFACE
  

  CONTAINS


  SUBROUTINE MolecCav_Initialize_krylov_basis(KBasis, TotH, Psi, Nb_krylov, Debug)
    USE QDUtil_m
    USE Sum_of_products_m
    IMPLICIT NONE
  
    real(kind=Rkind),        intent(inout) :: KBasis(:,:) ! the basis of the krylov basis set expressed in the wavefunction/hamiltonian/complete/original basis set
    TYPE(Sum_of_products_m), intent(in)    :: TotH        ! the hamiltonian in the original basis set
    real(kind=Rkind),        intent(in)    :: Psi(:)      ! the wavefunction in the original basis set
    integer,                 intent(in)    :: Nb_krylov   ! the size of the krylov basis
    logical, optional,       intent(in)    :: Debug

    logical                                :: Debug_local
    real(kind=Rkind), allocatable          :: KBasis_non_ortho(:,:) ! the intermediary KBasis before it is orthonormalized
    real(kind=Rkind), allocatable          :: Intermediary(:,:)     ! for the orthonormality check
    real(kind=Rkind), allocatable          :: Identity(:,:)         ! also
    integer                                :: i
  
    !------------------------------------------------------Debugging options-----------------------------------------------------
    IF (PRESENT(Debug)) THEN; Debug_local = Debug
    ELSE; Debug_local = .FALSE.; END IF

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "o Arguments of MolecCav_Initialize_krylov_basis :"
      WRITE(out_unit,*) "The <<KBasis>> argument :" 
      CALL Write_Mat(KBasis, out_unit, SIZE(KBasis, dim=2), info="KBasis")
      WRITE(out_unit,*) "The <<TotH>> argument :" 
      CALL Write(TotH)
      WRITE(out_unit,*) "The <<Psi>> argument :" 
      CALL Write_Vec(Psi, out_unit, 1, info="Psi")
      WRITE(out_unit,*) "The <<Nb_krylov>> argument :"//TO_string(Nb_krylov)
      FLUSH(out_unit)
    END IF

    !-----------------------------------------------------Checking dimensions----------------------------------------------------
    IF (.NOT. ALLOCATED(KBasis)) THEN
      WRITE(out_unit,*) "### KBasis has to have been allocated BEFORE calling MolecCav_Initialize_krylov_basis. Please check init&
      &ialization."
      STOP "### KBasis has to have been allocated BEFORE calling MolecCav_Initialize_krylov_basis. Please check initialization."
    END IF 

    IF (.NOT. ALLOCATED(Psi)) THEN
      WRITE(out_unit,*) "### KBasis cannot be constructed in MolecCav_Initialize_krylov_basis if Psi is not allocated. Please che&
      &ck initialization."
      STOP "### KBasis cannot be constructed in MolecCav_Initialize_krylov_basis if Psi is not allocated. Please check initializa&
      &tion."
    END IF

    IF (Nb_krylov /= SIZE(KBasis, dim=2)) THEN
      WRITE(out_unit,*) "### The number of vectors KBasis is allocated to does not match Nb_krylov. Please check initialization."
      WRITE(out_unit,*) "    SIZE(KBasis, dim=2) = "//TO_string(Size(KBasis, dim=2))//"; Nb_krylov = "//TO_string(Nb_krylov)
      STOP "### The number of vectors KBasis is allocated to does not match Nb_krylov. Please check initialization."
    END IF 

    IF (SIZE(Psi, dim=1) /= Size(KBasis, dim=1)) THEN
      WRITE(out_unit,*) "### The sizes of the KBasis's vectors do not match Psi's vector size i.e. the original basis set size. P&
      &lease check initialization."
      STOP "### The sizes of the KBasis's vectors do not match Psi's vector size i.e. the original basis set size. Please check i&
      &nitialization."
    END IF 

    !-----------------------Creation of the Krylov basis-----------------------
    ALLOCATE(KBasis_non_ortho(SIZE(Psi), Nb_krylov))

    KBasis_non_ortho(:,:) = ZERO
    KBasis_non_ortho(:,1) = Psi(:)
  
    DO i = 2, Nb_krylov
      CALL Action(KBasis_non_ortho(:,i), TotH, KBasis_non_ortho(:,i-1), Verbose=0, Debug=Debug_local)
    END DO
  
    !----------------------------Orthonormalization----------------------------
    CALL Gram_schmidt(KBasis, KBasis_non_ortho)
    CALL Gram_schmidt(KBasis, KBasis)
   
    !--------------------------------Check ortho-------------------------------
    ALLOCATE(Identity(Nb_krylov,Nb_krylov))
    Identity(:,:) = ZERO

    DO i = 1, Nb_krylov
      Identity(i,i) = 1
    END DO

    ALLOCATE(S(Nb_krylov,Nb_krylov))
    Intermediary = matmul(transpose(KBasis), KBasis)
    ! Intermediary = matmul(conjg(transpose(KBasis)), KBasis)

    IF (MAXVAL(ABS(Intermediary)-Identity) > 10E-1O) THEN
      WRITE(out_unit,*) "################### WARNING ################## WARNING ################## WARNING ##################"
      WRITE(out_unit,*) "                             The KBasis is not orthonormal up to 10E-10 "
      WRITE(out_unit,*) "################### WARNING ################## WARNING ################## WARNING ##################"
    END IF
    
  END SUBROUTINE MolecCav_Initialize_krylov_basis

 
  SUBROUTINE MolecCav_Augment_krylov_basis(KBasis, TotH)
    USE QDUtil_m
    USE Sum_of_products_m
    IMPLICIT NONE
  
    real(kind=Rkind), allocatable, intent(inout) :: KBasis(:,:)
    TYPE(Sum_of_products_t),       intent(in)    :: TotH

    real(kind=Rkind), allocatable                :: TemporaryKB(:,:)
    real(kind=Rkind), allocatable                :: V(:)
    real(kind=Rkind), allocatable                :: Intermediary(:,:)
    real(kind=Rkind), allocatable                :: Identity(:,:)
    integer                                      :: Nb_krylov
    integer                                      :: NB ! size of the original basis set
    integer                                      :: i

    !------------------------------Initialization------------------------------
    NB          = SIZE(KBasis, dim = 1)
    Nb_krylov   = SIZE(KBasis, dim = 2)
    TemporaryKB = KBasis
    ALLOCATE(V(NB))

    DEALLOCATE(KBasis)
    ALLOCATE(KBasis(NB, Nb_krylov+1))

    !-------------------------------Construction-------------------------------
    CALL Action(V, TotH, TemporaryKB(:,Nb_krylov), Verbose=Verbose, Debug=Debug) 
    
    KBasis(:,Nb_krylov+1) = V(:)

    DO i = 1, Nb_krylov
      KBasis(:,Nb_krylov+1) = KBasis(:,Nb_krylov+1) - dot_product(TemporaryKB(:,i), V(:))*TemporaryKB(:,i)
    END DO
    K(:,m+1) = K(:,m+1) / sqrt(dot_product(K(:,m+1),K(:,m+1)))

    K(:,1:m) = TemporaryKB
    DEALLOCATE(TemporaryKB, V)

    !--------------------------------Check ortho-------------------------------
    ALLOCATE(Identity(Nb_krylov+1,Nb_krylov+1))
    Identity(:,:) = ZERO
    DO i = 1, m+1
      Identity(i,i) = 1
    END DO
    ALLOCATE(S(m+1,m+1))

    S = matmul(conjg(transpose(K)), K)
    WRITE(out_unit, *) 'Smax = ', maxval(abs(S(:,:))-Identity(:,:))
    
  END SUBROUTINE MolecCav_Augment_krylov_basis


  SUBROUTINE MolecCav_Gram_schmidt(K, Q)
    USE QDUtil_m
    IMPLICIT NONE
    
    real(kind=Rkind), allocatable, intent(inout) :: K(:,:)
    real(kind=Rkind),              intent(in)    :: Q(:,:)

    real(kind=Rkind), allocatable                :: v(:)
    integer                                         :: im, jm, n, m

    n = size(Q, dim=1)
    m = size(Q, dim=2)    
    ALLOCATE(K(n,m))
    ALLOCATE(v(n))
    K(:,:) = ZERO
    v(:) = ZERO
   
    K(:,1) = Q(:,1) / sqrt(dot_product(Q(:,1),Q(:,1)))
    
    DO im = 2, m
      v(:) = Q(:, im)
      DO jm = 1, im-1
        v(:) = v(:) - dot_product(Q(:,im), K(:,jm)) * K(:,jm)
      END DO
      K(:, im) = v(:) / sqrt(dot_product(v(:),v(:)))
    END DO
            
  END SUBROUTINE MolecCav_Gram_schmidt
  
  
  SUBROUTINE MolecCav_Construct_H_tribande(triband_H, K, H)
    USE QDUtil_m
    IMPLICIT NONE

    real(kind=Rkind), intent(inout) :: triband_H(:,:)
    real(kind=Rkind), intent(in)    :: K(:,:)
    real(kind=Rkind),    intent(in)    :: H(:,:)

    integer                            :: ib, m

    !--------------------------------Allocation--------------------------------
    m = size(K, dim = 2)    
   
    !-------------------------------Construction-------------------------------
    triband_H = matmul(conjg(transpose(K)), matmul(H, K))
    !WRITE(out_unit, *) '------------------------------------'
    !CALL WRITE_Mat(triband_H,out_unit,5,info='triband_H') !ça prenait trop de place dans results
      
  END SUBROUTINE MolecCav_Construct_H_tribande


  SUBROUTINE MolecCav_KrylovTOBasis(Vec_B, Vec_K, K)
    USE QDUtil_m
    IMPLICIT NONE

    real(kind=Rkind), intent(inout) :: Vec_B(:,:)                                ! mat des vecteurs colonnes def sur a base des sinus
    real(kind=Rkind), intent(in)    :: Vec_K(:,:)                                ! mat des vecteurs colonnes def sur a base de Krylov
    real(kind=Rkind), intent(in)    :: K(:,:)

    integer                            :: i, j, m

    m = size(K, dim = 2)
    Vec_B(:,:) = ZERO
    Vec_B = matmul(K, Vec_K)

  END SUBROUTINE MolecCav_KrylovTOBasis


  SUBROUTINE MolecCav_Construct_psi_approx(Psidt, psi, Vec_B, Valp, dt)
    USE QDUtil_m
    IMPLICIT NONE

    real(kind=Rkind), intent(inout) :: Psidt(:)
    real(kind=Rkind), intent(in)    :: Vec_B(:,:), Valp(:), psi(:)
    real(kind=Rkind),    intent(in)    :: dt

    real(kind=Rkind), allocatable   :: c1(:), c2(:)
    integer                            :: m, i
   
    m = size(Vec_B, dim = 2)
    ALLOCATE(c1(m), c2(m))

    c1 = matmul(conjg(transpose(Vec_B)), Psi)
    c2(:) = c1(:)*exp(-EYE*dt*Valp(:))
  
    WRITE(out_unit, *) 'coeff m-1 = ', c2(m-1), 'coeff m = ', c2(m)

    Psidt = matmul(Vec_B, c2)

  END SUBROUTINE MolecCav_Construct_psi_approx
  

END MODULE