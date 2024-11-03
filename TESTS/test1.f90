PROGRAM test1
  USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : real64 
  IMPLICIT NONE

  TYPE :: MC_harmo_basis_t    ! MC = MolecCav
    integer           :: D    ! Number of the basis/OHL/mode/dimension
    integer           :: Nb_D ! number of basis vectors
    real(kind=real64) :: w_D  ! Eigenpulsation associated with the OHL
  END TYPE

  TYPE(MC_harmo_basis_t) :: Basis_mode_1
 
  WRITE(*,*) 'This is a test.'

END PROGRAM
