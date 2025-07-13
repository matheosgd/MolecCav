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
  USE Algebra_m
  USE Operator_ND_m
  USE Sum_of_products_m
  USE Transition_spectrum_m
  IMPLICIT NONE


  integer                       :: Verbose = 40
  logical                       :: Debug   = .TRUE.
  integer                       :: nioint, niospec, nioHanar, err_io

  logical                       :: Dense   = .TRUE.
  TYPE(Sum_of_products_t)       :: TotH
  TYPE(Sum_of_products_t)       :: DipMomt

  real(kind=Rkind)              :: Matm, Matw, Cavw, Matlambda, Cavlambda, lambda
  integer                       :: Nb_1, Nb_2, NB, J

  real(kind=Rkind), allocatable :: H_anar_in(:,:), H_anar(:,:)
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

  OPEN(NEWUNIT = nioHanar, FILE = '/home/segaud/MolecCav/DATA/H30x30', FORM = 'formatted', ACTION = 'read', POSITION = 'rewind', &
  &IOSTAT=err_io)
  WRITE(out_unit,*) "### err_io = "//TO_string(err_io)

  ALLOCATE(H_anar_in(30,30))
  CALL Read_Mat(H_anar_in, nioHanar, 5, err_io)
  WRITE(out_unit,*) "### err_io = "//TO_string(err_io)

  ALLOCATE(H_anar(Nb_1,Nb_1))
  H_anar = H_anar_in(1:Nb_1, 1:Nb_1)
  CALL Write_Mat(H_anar, out_unit, SIZE(H_anar, dim=2), info="H_anar")

  ALLOCATE(AnarREigval(Nb_1))
  ALLOCATE(AnarREigvec(Nb_1,Nb_1))
  CALL diagonalization(H_anar, AnarREigval, AnarREigvec)
  WRITE(out_unit,*)
  CALL Write_Vec(AnarREigval, out_unit, 1, info="MatH_anar EigenEnergies [Ha]")

  WRITE(out_unit,*)
  WRITE(out_unit,*) "--- HARMONIC tab_mat_ops"
  DO I=1, SIZE(tab_mat_ops)
    CALL Write(tab_mat_ops(I))
  END DO
  tab_mat_ops(1)%Tab_op(1)%Dense_val = H_anar
  WRITE(out_unit,*) "--- ANARHARMONIC tab_mat_ops"
  DO I=1, SIZE(tab_mat_ops)
    CALL Write(tab_mat_ops(I))
  END DO

  ALLOCATE(TotH_matrix(NB, NB))
  ALLOCATE(Phi(NB))
  TotH_matrix = ZERO

  DO J = 1, NB
    Phi = ZERO
    Phi(J) = ONE
    CALL Action(TotH_matrix(:,J), TotH, Phi, Verbose=Verbose, Debug=.FALSE.)
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
  CALL Initialize(TranSpec, TotREigvec, DipMomt, TotREigval, Nb_states=10, Verbose=Verbose, Debug=Debug)

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
!######### stretching #############
    ! Intensity = TranSpec%tab_ints(1)*1E4*Lorentzian(Energy, TranSpec%tab_energies(1)*Conversion, Gamma) &
    ! &         + TranSpec%tab_ints(4)*1E0*Lorentzian(Energy, TranSpec%tab_energies(4)*Conversion, Gamma) &
    ! &         + TranSpec%tab_ints(6)*1E5*Lorentzian(Energy, TranSpec%tab_energies(6)*Conversion, Gamma) &
    ! &         + TranSpec%tab_ints(8)*1E2*Lorentzian(Energy, TranSpec%tab_energies(8)*Conversion, Gamma)
!##################################
    Intensity = TranSpec%tab_ints(1)*Lorentzian(Energy, TranSpec%tab_energies(1)*Conversion, Gamma) &
    &         + TranSpec%tab_ints(2)*Lorentzian(Energy, TranSpec%tab_energies(2)*Conversion, Gamma) &
    &         + TranSpec%tab_ints(3)*Lorentzian(Energy, TranSpec%tab_energies(3)*Conversion, Gamma) &
    &         + TranSpec%tab_ints(4)*Lorentzian(Energy, TranSpec%tab_energies(4)*Conversion, Gamma) &
    &         + TranSpec%tab_ints(5)*Lorentzian(Energy, TranSpec%tab_energies(5)*Conversion, Gamma) &
    &         + TranSpec%tab_ints(6)*Lorentzian(Energy, TranSpec%tab_energies(6)*Conversion, Gamma) &
    &         + TranSpec%tab_ints(7)*Lorentzian(Energy, TranSpec%tab_energies(7)*Conversion, Gamma) &
    &         + TranSpec%tab_ints(8)*Lorentzian(Energy, TranSpec%tab_energies(8)*Conversion, Gamma) &
    &         + TranSpec%tab_ints(9)*Lorentzian(Energy, TranSpec%tab_energies(9)*Conversion, Gamma)
    !======================================================
    !
    ! For the 2 photons experiment (just remove 2 first tr-
    ! ansitions) : if you try you are not supposed to see 
    ! anything on the spectra, transition intensities are
    ! really null.
    !
    !======================================================
    ! Intensity = TranSpec%tab_ints(4)*Lorentzian(Energy, TranSpec%tab_energies(4)*Conversion, Gamma) &
    ! &         + TranSpec%tab_ints(5)*Lorentzian(Energy, TranSpec%tab_energies(5)*Conversion, Gamma) &
    ! &         + TranSpec%tab_ints(6)*Lorentzian(Energy, TranSpec%tab_energies(6)*Conversion, Gamma) &
    ! &         + TranSpec%tab_ints(7)*Lorentzian(Energy, TranSpec%tab_energies(7)*Conversion, Gamma) &
    ! &         + TranSpec%tab_ints(8)*Lorentzian(Energy, TranSpec%tab_energies(8)*Conversion, Gamma) &
    ! &         + TranSpec%tab_ints(9)*Lorentzian(Energy, TranSpec%tab_energies(9)*Conversion, Gamma)
    WRITE(niospec, *) Energy, Intensity
  END DO

  !----------------------------computing the Nmodes---------------------------
  CALL Compute_normal_modes(Nmodes, Ncoos, Matm, Matw, Cavw, Matlambda, Cavlambda, ONE, .TRUE.)

  
  CONTAINS


  SUBROUTINE Compute_normal_modes(Nmodes_l, Ncoos_l, Matm_l, Matw_l, Cavw_l, Matlambda_l, Cavlambda_l, CoeffDipMomt_l, Debug_l)
    USE QDUtil_m
    IMPLICIT NONE 

    real(kind=Rkind),  intent(inout) :: Nmodes_l(2)                             ! VP of the MWH
    real(kind=Rkind),  intent(inout) :: Ncoos_l(2,2)                     ! \overrightarrow{VP} of the MWH
    real(kind=Rkind),  intent(in)    :: Matm_l
    real(kind=Rkind),  intent(in)    :: Matw_l
    real(kind=Rkind),  intent(in)    :: Cavw_l
    real(kind=Rkind),  intent(in)    :: Matlambda_l 
    real(kind=Rkind),  intent(in)    :: Cavlambda_l
    real(kind=Rkind),  intent(in)    :: CoeffDipMomt_l
    logical, optional, intent(in)    :: Debug_l

    real(kind=Rkind)                 :: Cross
    real(kind=Rkind)                 :: MWH(2,2)
    logical                          :: Debug_local = .TRUE.


    !------------------------------------------------------Debugging options-----------------------------------------------------
    IF (PRESENT(Debug_l))   THEN; Debug_local = Debug_l
    ELSE; Debug_local = .FALSE.; END IF

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "--- Parameters for the mass-weightened Hessian :"
      WRITE(out_unit,*) "Matm      = "//TO_string(Matm_l)
      WRITE(out_unit,*) "Matw      = "//TO_string(Matw_l)
      WRITE(out_unit,*) "Cavw      = "//TO_string(Cavw_l)
      WRITE(out_unit,*) "Matlambda = "//TO_string(Matlambda_l)
      WRITE(out_unit,*) "Cavlambda = "//TO_string(Cavlambda_l)
      WRITE(out_unit,*) "lambda    = "//TO_string(Matlambda_l*Cavlambda_l)
      WRITE(out_unit,*) "--- End MWH parameters"
    END IF 


    !------------------------------------------------------Computing-----------------------------------------------------
    Cross    = Cavlambda_l*Matlambda_l*Cavw_l*CoeffDipMomt_l / SQRT(Matm_l)
    MWH      = ZERO
    MWH(1,1) = Matw**2
    MWH(2,2) = Cavw**2
    MWH(1,2) = Cross
    MWH(2,1) = MWH(1,2)

    IF (Debug_local) THEN
      CALL Write_Mat(MWH, out_unit, Size(MWH, dim=2), info="MWH")
    END IF

    CALL diagonalization(MWH, Nmodes_l, Ncoos_l)
    IF (Debug_local) CALL Write_Vec(Nmodes_l,       out_unit, Size(Nmodes_l),              info="Normal modes")
    IF (Debug_local) CALL Write_Mat(Ncoos_l, out_unit, Size(Ncoos_l, dim=2), info="Normal coordniates")

    WRITE(out_unit,*)
    WRITE(out_unit,*) "--- Normal modes checking :"
    DO i = 1, Size(Nmodes_l)
      IF (Nmodes_l(i) >= 0) THEN
        WRITE(out_unit,*) TO_string(i)//"^{th} Normal coordinate has positive squared frequency, lead&
                         &ing to w_"//TO_string(i)//" = "//TO_string(SQRT(Nmodes_l(i)))
      ELSE
        WRITE(out_unit,*) TO_string(i)//"^{th} Normal coordinate has NEGATIVE squared frequency, lead&
        &ing to w_"//TO_string(i)//" = "//TO_string(EYE*SQRT(-Nmodes_l(i)))
      END IF
    END DO
    WRITE(out_unit,*) "Expected ZPE by half-sum of the total system eigenpulsations (from MWH): "//TO&
                    &_string( ( SQRT(Nmodes_l(1))+SQRT(Nmodes_l(2)) )/2 )

  END SUBROUTINE Compute_normal_modes


END PROGRAM
