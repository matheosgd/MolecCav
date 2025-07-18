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
PROGRAM MolecCav
  !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
  USE QDUtil_m
  USE Sum_of_products_m
  USE Transition_spectrum_m
  IMPLICIT NONE


  integer                       :: Verbose = 40
  logical                       :: Debug   = .FALSE.
  integer                       :: nioint, niospec

  logical                       :: Dense   = .FALSE.
  TYPE(Sum_of_products_t)       :: TotH
  TYPE(Sum_of_products_t)       :: DipMomt

  real(kind=Rkind)              :: Matm, Matw, Cavw, Matlambda, Cavlambda
  integer                       :: Nb_mode, NB, J
  real(kind=Rkind), allocatable :: Phi(:)
  real(kind=Rkind), allocatable :: TotH_matrix(:,:)
  real(kind=Rkind), allocatable :: REigvec(:,:)
  real(kind=Rkind), allocatable :: REigval(:)

  real(kind=Rkind)              :: E_threshold
  integer                       :: N_states
  TYPE(Transition_spectrum_t)   :: TranSpec

  real(kind=Rkind)              :: Gamma, Start_plot, Stop_plot, Step_plot, Conversion, Energy, Intensity
  integer                       :: I


  !-------------------------Sum_of_products operators initialization-------------------------
  CALL Initialize_totH(TotH, in_unit, Dense=Dense, Verbose=Verbose, Debug=Debug)
  IF (Debug) THEN
    WRITE(out_unit,*)
    WRITE(out_unit,*) "--- TotH object constructed by MolecCav_Initialize_total_hamiltonian"
    CALL Write(TotH)
  END IF

  CALL Initialize_dipmomt(DipMomt, in_unit, Dense=Dense, Verbose=Verbose, Debug=Debug)
  IF (Debug) THEN
    WRITE(out_unit,*)
    WRITE(out_unit,*) "--- DipMomt object constructed by MolecCav_Initialize_dipmomt"
    CALL Write(DipMomt)
  END IF

  !-------------------------System initialization-------------------------
  NB = 1
  DO I = 1, TotH%N_mat
    WRITE(out_unit,*)
    WRITE(out_unit,*) "--- System parameters :"
    CALL Get(Matm,  "m", "Matter", I)
    WRITE(out_unit,*) "Mat"//TO_string(I)//"m      = "//TO_string(Matm)
    CALL Get(Matw, "w", "Matter", I)
    WRITE(out_unit,*) "Mat"//TO_string(I)//"w     = "//TO_string(Matw)
    CALL Get(Matlambda, "lambda", "Matter", I)
    WRITE(out_unit,*) "Mat"//TO_string(I)//"lambda = "//TO_string(Matlambda)
    CALL Get(Nb_mode, "Nb", "Matter", I)
    NB = NB*Nb_mode
    WRITE(out_unit,*) "Mat"//TO_string(I)//"Nb    = "//TO_string(Nb_mode)
  END DO
  DO I = 1, TotH%N_cav
    CALL Get(Cavw,  "w", "Cavity", I)
    WRITE(out_unit,*) "Cav"//TO_string(I)//"w     = "//TO_string(Cavw)
    CALL Get(Cavlambda, "lambda", "Cavity", I)
    WRITE(out_unit,*) "Cav"//TO_string(I)//"lambda = "//TO_string(Cavlambda)
    CALL Get(Nb_mode, "Nb", "Cavity", I)
    NB = NB*Nb_mode
    WRITE(out_unit,*) "Cav"//TO_string(I)//"Nb    = "//TO_string(Nb_mode)
  END DO 

  ALLOCATE(TotH_matrix(NB, NB))
  ALLOCATE(Phi(NB))
  TotH_matrix = ZERO

  DO J = 1, NB
    Phi = ZERO
    Phi(J) = ONE
    CALL Action(TotH_matrix(:,J), TotH, Phi, Verbose=Verbose, Debug=Debug)
  END DO

  WRITE(out_unit,*)
  WRITE(out_unit,*) "*** RESULTING MATRIX OF TOTH"
  !CALL Write_Mat(TotH_matrix, out_unit, NB, info="TotH_matrix")

  ALLOCATE(REigval(NB))
  ALLOCATE(REigvec(NB,NB))

  CALL diagonalization(TotH_matrix, REigval, REigvec)
  WRITE(out_unit,*)
  !CALL Write_Vec(REigval, out_unit, 1, info="EigenEnergies")
  

  !----------------------------Computing of the spectra---------------------------
  ! CALL Initialize(TranSpec, REigvec, DipMomt, E_threshold=E_threshold, REigval=REigval, Nb_states=N_states, Verbose=Verbose, Deb&
  ! &ug=Debug)
  ! CALL Initialize(TranSpec, REigvec, DipMomt, E_threshold=E_threshold, REigval=REigval, Nb_states=N_states, Verbose=Verbose, Deb&
  ! &ug=Debug)
  CALL Initialize(TranSpec, REigvec, DipMomt, REigval, Nb_states=10, Verbose=Verbose, Debug=Debug)

  WRITE(out_unit,*) "Spectrum information"
  CALL Write_Vec(TranSpec%tab_energies, out_unit, SIZE(TranSpec%tab_energies), info="Transition Energies")
  CALL Write_Vec(TranSpec%tab_ints,     out_unit, SIZE(TranSpec%tab_ints),     info="Transition Intensities")


  OPEN(NEWUNIT = nioint, FILE = 'OUT/TransInts_MolecCav_wmat'//TO_string(REAL(Mat2w,kind=RkS))//'_DTpm'//TO_string(REAL(Mat2w-Mat1w,ki&
  &nd=RkS))//'_wcav'//TO_string(REAL(Cavw,kind=RkS))//'_lamb'//TO_string(REAL(lambda,kind=RkS))//'.txt',  FORM = 'formatted', ACT&
  &ION = 'write', POSITION = 'rewind')

  OPEN(NEWUNIT = niospec, FILE = 'OUT/Spectrum_MolecCav_wmat'//TO_string(REAL(Mat2w,kind=RkS))//'_DTpm'//TO_string(REAL(Mat2w-Mat1w,ki&
  &nd=RkS))//'_wcav'//TO_string(REAL(Cavw,kind=RkS))//'_lamb'//TO_string(REAL(lambda,kind=RkS))//'.txt',  FORM = 'formatted', ACT&
  &ION = 'write', POSITION = 'rewind')

  !----------------------------computing of the spectra---------------------------
  WRITE(nioint, *) "Transition Energy -------- Transition intensity"
  DO J = 1, TranSpec%N_trstns
    IF (TranSpec%tab_ints(J)>1E-10) WRITE(nioint, *) TranSpec%tab_energies(J), TranSpec%tab_ints(J)
    IF (TranSpec%tab_ints(J)<1E-10) WRITE(nioint, *) TranSpec%tab_energies(J), 0
  END DO


  Conversion = 219474.6                       ! 1 Ha = Conversion.cm-1
  Gamma      = 1.0                            ! => 10cm-1
  Start_plot = 1.8794295214982763E-002 - 3E-4 ! in Ha 
  Stop_plot  = 1.8986432304458810E-002 + 3E-4 ! in Ha
  Step_plot  = (Stop_plot - Start_plot) / 900 ! divide by desired number of points

  WRITE(niospec, *) "Energy -------- Transition intensity"
  DO I = 1, 900
    Energy    = (Start_plot + I*Step_plot)*Conversion
    Intensity = TranSpec%tab_ints(1)*Lorentzian(Energy, TranSpec%tab_energies(1)*Conversion, Gamma)
    DO I = 2, 9
    Intensity = Intensity TranSpec%tab_ints(I)*Lorentzian(Energy, TranSpec%tab_energies(I)*Conversion, Gamma)
    END DO
    WRITE(niospec, *) Energy, Intensity, 2.0
  END DO

  
  CONTAINS


  SUBROUTINE Compute_normal_modes(Nmodes_l, Ncoos_l, Matm_l, Mat1w_l, Mat2w_l, Mat3w_l, Cavw_l, Matlambda_l, Cavlambda_l, &
    &CoeffDipMomt_l, Debug_l)
    USE QDUtil_m
    IMPLICIT NONE 

    real(kind=Rkind),  intent(inout) :: Nmodes_l(4)                      ! VP of the MWH
    real(kind=Rkind),  intent(inout) :: Ncoos_l(4,4)                     ! \overrightarrow{VP} of the MWH
    real(kind=Rkind),  intent(in)    :: Matm_l
    real(kind=Rkind),  intent(in)    :: Mat1w_l
    real(kind=Rkind),  intent(in)    :: Mat2w_l
    real(kind=Rkind),  intent(in)    :: Mat3w_l
    real(kind=Rkind),  intent(in)    :: Cavw_l
    real(kind=Rkind),  intent(in)    :: Matlambda_l 
    real(kind=Rkind),  intent(in)    :: Cavlambda_l
    real(kind=Rkind),  intent(in)    :: CoeffDipMomt_l
    logical, optional, intent(in)    :: Debug_l

    real(kind=Rkind)                 :: Cross_1
    real(kind=Rkind)                 :: Cross_2
    real(kind=Rkind)                 :: Cross_3
    real(kind=Rkind)                 :: MWH(4,4)
    logical                          :: Debug_local = .TRUE.


    !------------------------------------------------------Debugging options-----------------------------------------------------
    IF (PRESENT(Debug_l))   THEN; Debug_local = Debug_l
    ELSE; Debug_local = .FALSE.; END IF

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "--- Parameters for the mass-weightened Hessian :"
      WRITE(out_unit,*) "Matm      = "//TO_string(Matm_l)
      WRITE(out_unit,*) "Mat1w     = "//TO_string(Mat1w_l)
      WRITE(out_unit,*) "Mat2w     = "//TO_string(Mat2w_l)
      WRITE(out_unit,*) "Mat3w     = "//TO_string(Mat3w_l)
      WRITE(out_unit,*) "Cavw      = "//TO_string(Cavw_l)
      WRITE(out_unit,*) "Matlambda = "//TO_string(Matlambda_l)
      WRITE(out_unit,*) "Cavlambda = "//TO_string(Cavlambda_l)
      WRITE(out_unit,*) "lambda    = "//TO_string(Matlambda_l*Cavlambda_l)
      WRITE(out_unit,*) "--- End MWH parameters"
    END IF 


    !------------------------------------------------------Computing-----------------------------------------------------
    Cross_1  = Cavlambda_l*Matlambda_l*Cavw_l*CoeffDipMomt_l / SQRT(Matm_l)
    Cross_2  = Cavlambda_l*Matlambda_l*Cavw_l*CoeffDipMomt_l / SQRT(Matm_l)
    Cross_3  = Cavlambda_l*Matlambda_l*Cavw_l*CoeffDipMomt_l / SQRT(Matm_l)
    MWH      = ZERO
    MWH(1,1) = Mat1w**2
    MWH(2,2) = Mat2w**2
    MWH(3,3) = Mat3w**2
    MWH(4,4) = Cavw**2
    MWH(1,4) = Cross_1
    MWH(4,1) = MWH(1,4)
    MWH(2,4) = Cross_2
    MWH(4,2) = MWH(2,4)
    MWH(3,4) = Cross_3
    MWH(4,3) = MWH(3,4)

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
