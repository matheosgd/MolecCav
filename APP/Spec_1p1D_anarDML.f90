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
PROGRAM Spec_1p1D_anar
  !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
  USE QDUtil_m
  USE Operator_ND_m
  USE Sum_of_products_m
  USE Transition_spectrum_m
  IMPLICIT NONE


  integer                       :: Verbose = 20!40
  logical                       :: Debug   = .FALSE.
  integer                       :: nioint, niospec, nioHanar, err_io

  logical                       :: Dense   = .TRUE.
  TYPE(Sum_of_products_t)       :: TotH
  TYPE(Sum_of_products_t)       :: DipMomt

  real(kind=Rkind)              :: Matm, Matw, Cavw, Matlambda, Cavlambda, lambda
  integer                       :: Nb_1, Nb_2, NB, J

  real(kind=Rkind), allocatable :: H_anar_in(:,:), H_anar(:,:)
  real(kind=Rkind), allocatable :: Mu_anar_in(:,:), Mu_anar(:,:)

  real(kind=Rkind), allocatable :: AnarREigvec(:,:)
  real(kind=Rkind), allocatable :: AnarREigval(:)

  real(kind=Rkind), allocatable :: Phi(:)
  real(kind=Rkind), allocatable :: TotH_matrix(:,:)
  real(kind=Rkind), allocatable :: TotREigvec(:,:)
  real(kind=Rkind), allocatable :: TotREigval(:)

  real(kind=Rkind)              :: E_threshold
  integer                       :: N_states
  TYPE(Transition_spectrum_t)   :: TranSpec

  real(kind=Rkind)              :: Gamma, Start_plot, Stop_plot, Step_plot, Conversion, Energy, Intensity
  integer                       :: I
  real(kind=Rkind)              :: Nmodes(2), Ncoos(2,2)


  !-------------------------Sum_of_products operators initialization-------------------------
  CALL Initialize_totH(TotH, in_unit, Dense=Dense, Verbose=Verbose, Debug=Debug)
  IF (Debug) THEN
    WRITE(out_unit,*)
    WRITE(out_unit,*) "--------------TotH object constructed by MolecCav_Initialize_total_hamiltonian--------------"
    CALL Write(TotH)
    WRITE(out_unit,*) "------------End TotH object constructed by MolecCav_Initialize_total_hamiltonian------------"
  END IF

  CALL Initialize_dipmomt(DipMomt, in_unit, Dense=Dense, Verbose=Verbose, Debug=Debug)
  IF (Debug) THEN
    WRITE(out_unit,*)
    WRITE(out_unit,*) "--------------DipMomt object constructed by MolecCav_Initialize_dipmomt--------------"
    CALL Write(DipMomt)
    WRITE(out_unit,*) "------------End DipMomt object constructed by MolecCav_Initialize_dipmomt------------"
  END IF

  !-------------------------System initialization-------------------------
  CALL Get(Matm, "m", "Matter", 1)
  CALL Get(Matw, "w", "Matter", 1)
  CALL Get(Cavw, "w", "Cavity", 1)
  CALL Get(Matlambda, "lambda", "Matter", 1)
  CALL Get(Cavlambda, "lambda", "Cavity", 1)
  lambda = Matlambda * Cavlambda
  CALL Get(Nb_1, "Nb", "Matter", 1)
  CALL Get(Nb_2, "Nb", "Cavity", 1)
  NB = Nb_1 * Nb_2

  WRITE(out_unit,*)
  WRITE(out_unit,*) "--- System parameters :"
  WRITE(out_unit,*) "Matm      = "//TO_string(Matm)
  WRITE(out_unit,*) "Matw      = "//TO_string(Matw)
  WRITE(out_unit,*) "Cavw      = "//TO_string(Cavw)
  WRITE(out_unit,*) "Matlambda = "//TO_string(Matlambda)
  WRITE(out_unit,*) "Cavlambda = "//TO_string(Cavlambda)
  WRITE(out_unit,*) "lambda    = "//TO_string(lambda)
  WRITE(out_unit,*) "MatNb     = "//TO_string(Nb_1)
  WRITE(out_unit,*) "CavNb     = "//TO_string(Nb_2)

  !============================================================================
  !!!! read anaharmonic Hamiltonian matrix + diagonalization (test)
  OPEN(NEWUNIT = nioHanar, FILE = 'DATA/H30x30', FORM = 'formatted', ACTION = 'read', IOSTAT=err_io)
  WRITE(out_unit,*) "### err_io = "//TO_string(err_io)

  ALLOCATE(H_anar_in(30,30))
  CALL Read_Mat(H_anar_in, nioHanar, 5, err_io)
  WRITE(out_unit,*) "### err_io = "//TO_string(err_io)
  CLOSE(nioHanar)

  ALLOCATE(H_anar(Nb_1,Nb_1))
  H_anar = H_anar_in(1:Nb_1, 1:Nb_1)
  CALL Write_Mat(H_anar, out_unit, SIZE(H_anar, dim=2), info="H_anar")

  ALLOCATE(AnarREigval(Nb_1))
  ALLOCATE(AnarREigvec(Nb_1,Nb_1))
  CALL diagonalization(H_anar, AnarREigval, AnarREigvec)
  WRITE(out_unit,*)
  CALL Write_Vec(AnarREigval, out_unit, 1, info="MatH_anar EigenEnergies [Ha]")
  !============================================================================
  !============================================================================
  !!!! read anaharmonic Dipole moment (1 component) matrix
  OPEN(NEWUNIT = nioHanar, FILE = 'DATA/Mu30x30', FORM = 'formatted', ACTION = 'read', IOSTAT=err_io)
  WRITE(out_unit,*) "### err_io = "//TO_string(err_io)

  ALLOCATE(Mu_anar_in(30,30))
  CALL Read_Mat(Mu_anar_in, nioHanar, 5, err_io)
  WRITE(out_unit,*) "### err_io = "//TO_string(err_io)
  CLOSE(nioHanar)

  ALLOCATE(Mu_anar(Nb_1,Nb_1))
  Mu_anar = Mu_anar_in(1:Nb_1, 1:Nb_1)
  CALL Write_Mat(Mu_anar, out_unit, SIZE(Mu_anar, dim=2), info="Mu_anar")
  !============================================================================

  !============================================================================
  WRITE(out_unit,*) "--- HARMONIC tab_mat_ops"
  DO I=1, SIZE(tab_mat_ops)
    CALL Write(tab_mat_ops(I))
  END DO
  tab_mat_ops(1)%Tab_op(1)%Dense_val = H_anar !Hamiltonian
  tab_mat_ops(1)%Tab_op(4)%Dense_val = Mu_anar ! dipmomt (1 component)
  WRITE(out_unit,*) "--- ANARHARMONIC tab_mat_ops"
  DO I=1, SIZE(tab_mat_ops)
    CALL Write(tab_mat_ops(I),More=.TRUE.)
  END DO
  !============================================================================
  ALLOCATE(TotH_matrix(NB, NB))
  ALLOCATE(Phi(NB))
  TotH_matrix = ZERO

  DO J = 1, NB
    Phi = ZERO
    Phi(J) = ONE
    CALL Action(TotH_matrix(:,J), TotH, Phi, Verbose=0, Debug=.FALSE.)
  END DO

  WRITE(out_unit,*)
  WRITE(out_unit,*) "*** RESULTING MATRIX OF TOTH"
  CALL Write_Mat(TotH_matrix, out_unit, NB, info="TotH_matrix")

  ALLOCATE(TotREigval(NB))
  ALLOCATE(TotREigvec(NB,NB))

  CALL diagonalization(TotH_matrix, TotREigval, TotREigvec)
  WRITE(out_unit,*)
  CALL Write_Vec(TotREigval, out_unit, 1, info="TotH EigenEnergies")
  

  !----------------------------Computing of the spectra---------------------------
  ! CALL Initialize(TranSpec, REigvec, DipMomt, E_threshold=E_threshold, REigval=REigval, Nb_states=N_states, Verbose=Verbose, Deb&
  ! &ug=Debug)
  ! CALL Initialize(TranSpec, REigvec, DipMomt, E_threshold=E_threshold, REigval=REigval, Nb_states=N_states, Verbose=Verbose, Deb&
  ! &ug=Debug)
  CALL Initialize(TranSpec, TotREigvec, DipMomt, TotREigval, Nb_states=10, Verbose=0, Debug=Debug)

  WRITE(out_unit,*) "Spectrum information"
  CALL Write_Vec(TranSpec%tab_energies, out_unit, SIZE(TranSpec%tab_energies), info="Transition Energies")
  CALL Write_Vec(TranSpec%tab_ints,     out_unit, SIZE(TranSpec%tab_ints),     info="Transition Intensities")


  OPEN(NEWUNIT = nioint,  FILE = 'OUT/TransInts_ANAR_1p1D_wmat'//TO_string(REAL(Matw,kind=RkS))//'_wcav'//TO_string(REAL(Cavw,kin&
  &d=RkS))//'_lamb'//TO_string(REAL(lambda,kind=RkS))//'.txt',  FORM = 'formatted', ACTION = 'write', POSITION = 'rewind')

  OPEN(NEWUNIT = niospec, FILE = 'OUT/Spectrum_ANAR_1p1D_wmat'//TO_string(REAL(Matw,kind=RkS))//'_wcav'//TO_string(REAL(Cavw,kind&
  &=RkS))//'_lamb'//TO_string(REAL(lambda,kind=RkS))//'.txt',  FORM = 'formatted', ACTION = 'write', POSITION = 'rewind')

  !----------------------------computing the spectra---------------------------
  WRITE(nioint, *) "Transition Energy -------- Transition intensity"
  DO J = 1, TranSpec%N_trstns
    IF (TranSpec%tab_ints(J)>1E-10) WRITE(nioint, *) TranSpec%tab_energies(J), TranSpec%tab_ints(J)
    IF (TranSpec%tab_ints(J)<1E-10) WRITE(nioint, *) TranSpec%tab_energies(J), 0
  END DO


  Conversion = 219474.6                        ! 1 Ha = Conversion.cm-1
  Gamma      = 1.0                             ! => 10cm-1
  ! Start_plot = 5.8661769977076117E-003 - 1E-4  ! in Ha 
  ! Stop_plot  = 2.9871171804932286E-002 + 1E-4  ! in Ha
  Start_plot = 1.8030078508707852E-002 - 1E-4  ! in Ha 
  Stop_plot  = 1.8165407767964686E-002 + 1E-4  ! in Ha
  !--- for the 2 photons experiments ----------------------
  ! Start_plot = 3.6083210204059454E-002 - 1E-4  ! in Ha 
  ! Stop_plot  = 3.6185678195890060E-002 + 1E-4  ! in Ha
  Step_plot  = (Stop_plot - Start_plot) / 9000 ! divide by desired number of points

  WRITE(niospec, *) "Energy -------- Transition intensity"
  DO I = 1, 9000
    Energy    = (Start_plot + I*Step_plot)*Conversion
    Intensity = ZERO
    DO J = 1, TranSpec%N_trstns
      IF (TranSpec%tab_ints(J) > ONETENTH**10) THEN
        Intensity = Intensity + TranSpec%tab_ints(J)*Lorentzian(Energy, TranSpec%tab_energies(J)*Conversion, Gamma)
      END IF
    END DO  
    WRITE(niospec, *) Energy, Intensity
  END DO


END PROGRAM
