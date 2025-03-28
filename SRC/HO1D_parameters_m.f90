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
! README :
! The module to initialize the HO by reading its parameters from the namelist.    
!==================================================================================================
!==================================================================================================
MODULE HO1D_parameters_m
  !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
  USE QDUtil_m                                                                                                                   ! gives Rkind=real64; out_unit=OUTPUT_UNIT; INPUT_UNIT=in_unit; EYE=i and other numbers; TO_LOWERCASE; TO_UPPERCASE; Write_Mat and Write_Vec; TO_string;... We thereby use ZERO instead of 0.0_real64
  IMPLICIT NONE


  TYPE :: HO1D_parameters_t                                                                                                      ! NB: everything is initialized at values that are not supposed to make it possible of the cavity mode lecture/creation have successfully been executed
    integer          :: Nb     = 0                                                                                               ! the size of the basis set used hereafter to describe the HO. The basis functions are the eigenvectors/functions of this HO, so the Hamiltonian matrix is diagonal in this basis.
    real(kind=Rkind) :: w      = ZERO                                                                                            ! the eigenpulsation associated with the HO
    real(kind=Rkind) :: m      = ZERO                                                                                            ! the mass associated with the HO
    real(kind=Rkind) :: eq_pos = -ONE                                                                                            ! the equilibrium position of the HO
  END TYPE


  INTERFACE Read_HO1D_parameters
    MODULE PROCEDURE MolecCav_Read_HO1D_parameters
  END INTERFACE
  INTERFACE Write_HO1D_parameters
    MODULE PROCEDURE MolecCav_Write_HO1D_parameters
  END INTERFACE
    

  CONTAINS


  SUBROUTINE MolecCav_Read_HO1D_parameters(HO1D_para, nio, Verbose, Debug)                                                                       ! nio is the label of the file from which the values have to be drawn.
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    IMPLICIT NONE
    
    TYPE(HO1D_parameters_t), intent(inout) :: HO1D_para   
    integer,                 intent(in)    :: nio
    integer,                 intent(in)    :: Verbose
    logical,                 intent(in)    :: Debug

    integer                                :: Nb, err_io                                                                             ! the number of basis vectors of the HO, and an error control variable
    real(kind=Rkind)                       :: w, m, eq_pos                                                                           ! eigenpulsation, mass, molecule-coupling strength, and equilibrium position associated with this HO
    integer                                :: Verbose_local = 0                                                                      ! controls the level of information printing (result and calculations, not intermediary results which are related to Debug). From 20 to 24 at this layer
    logical                                :: Debug_local   = .FALSE.                                                                ! controls the printing of intermediary results for checking and debugging

    IF (PRESENT(Verbose)) Verbose_local = Verbose
    IF (PRESENT(Debug))   Debug_local   = Debug

    IF (Verbose_local > 20) WRITE(out_unit,*) 
    IF (Verbose_local > 20) WRITE(out_unit,*) "-------------------------------------------------INITIALIZING THE HO1D PARAMETERS-&
                                              &------------------------------------------------"; FLUSH(out_unit)

    NAMELIST /HO_1/ Nb, w, m, eq_pos                                                                                             ! declare the nml HO_1 and specify the parameter's list to be found within

    !----------------------------------------------Initialization to default values----------------------------------------------
    Nb     = 1
    w      = ZERO
    m      = ZERO
    eq_pos = -ONE
 
    !-----------------------------------------------------Reading of the nml-----------------------------------------------------
    IF (Verbose_local > 22) WRITE(out_unit,*) 
    IF (Verbose_local > 22) WRITE(out_unit,*) "--------------------------------------------Reading the namelist of the HO1D param&
                                              &eters-------------------------------------------"
    
    READ(nio, nml = HO_1, iostat = err_io)                                                                                       ! assign the values read in the nml to the declared list of parameters

    IF (Debug_local) THEN
      WRITE(out_unit,*) "--- The namelist parameters are read as :"
      WRITE(out_unit, nml = HO_1)
      WRITE(out_unit,*) "--- End of the namelist parameters"
    END IF
    
    !-----------------------------------------------------Check reading error----------------------------------------------------
    IF(err_io /= 0) THEN
      WRITE(out_unit,*); WRITE(out_unit,*) "### Error reading the nml in Read_HO1D_parameters (err_io/=0)"
      WRITE(out_unit,*) "    err_io = "//TO_string(err_io)
      STOP "### Error Read_HO1D_parameters. Please check basis data"
    END IF

    IF (Nb == 0) THEN
      WRITE(out_unit,*) "The number of basis vector associated to any HO CANNOT be 0 (what are are you going to study if there is&
                       & no system ???). Please check the data file '.nml'"
      STOP "### The number of basis vector associated to any HO CANNOT be 0 (what are are you going to study if there is&
                       & no system ???). Please check the data file '.nml'"
    END IF
    
    !----------------------------------------Construction of the Cavity_mode_t type object---------------------------------------
    IF (Verbose_local > 22) WRITE(out_unit,*) 
    IF (Verbose_local > 22) WRITE(out_unit,*) "--------------------------------------------Constructing the HO1D_parameters_t----&
                                              &---------------------------------------"
    HO1D_para%Nb     = Nb
    HO1D_para%w      = w
    HO1D_para%m      = m
    HO1D_para%eq_pos = eq_pos

    IF (Verbose_local > 20) WRITE(out_unit,*) 
    IF (Verbose_local > 20) WRITE(out_unit,*) "-------------------------------------------------HO1D PARAMETERS INITIALIZED------&
                                              &-------------------------------------------"; FLUSH(out_unit)

    IF (Verbose_local > 21) THEN
      IF (Verbose_local < 23) WRITE(out_unit,*)
      WRITE(out_unit,*) "--- Cavity mode constructed by MolecCav_Read_HO1D_parameters :"
      CALL Write_cavity_mode(HO1D_para)
      WRITE(out_unit,*) "--- End Cavity mode constructed by MolecCav_Read_HO1D_parameters"
    END IF

  END SUBROUTINE MolecCav_Read_HO1D_parameters

  
  SUBROUTINE MolecCav_Write_HO1D_parameters(HO1D_para)
    TYPE(HO1D_parameters_t), intent(in) :: HO1D_para

    WRITE(out_unit,*) "_____________________________________The parameters of the 1D HO____________________________________"
    WRITE(out_unit,*) "|Basis set size of the HO (HO1D_para%Nb)                                     | ", HO1D_para%Nb
    WRITE(out_unit,*) "|____________________________________________________________________________|______________________"
    WRITE(out_unit,*) "|Eigenpulsation of the HO (HO1D_para%w)                                      | ", HO1D_para%w
    WRITE(out_unit,*) "|____________________________________________________________________________|______________________"
    WRITE(out_unit,*) "|Mass associated with the HO (HO1D_para%m)                                   | ", HO1D_para%m
    WRITE(out_unit,*) "|____________________________________________________________________________|______________________"
    WRITE(out_unit,*) "|Equilibrium position of the HO (HO1D_para%eq_pos)                           | ", HO1D_para%eq_pos
    WRITE(out_unit,*) "|____________________________________________________________________________|______________________"
    FLUSH(out_unit)

  END SUBROUTINE MolecCav_Write_HO1D_parameters


END MODULE