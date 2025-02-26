!==================================================================================================
!==================================================================================================
!This file is part of MolecCav.
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
MODULE Cavity_mode_m
  !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
  USE QDUtil_m                                                                 ! gives Rkind=real64; out_unit=OUTPUT_UNIT; INPUT_UNIT=in_unit; EYE=i and other numbers; TO_LOWERCASE; TO_UPPERCASE;... We thereby use ZERO instead of 0.0_real64
  IMPLICIT NONE


  TYPE :: Cavity_mode_t                                                        ! MC = MolecCav NB: everything is initialized at values that are not supposed to make it possible of the cavity mode lecture/creation have successfully been executed
    integer          :: D      = 0                                             ! label of the HO/mode/dimension/associated basis set
    integer          :: Nb     = 0                                             ! number of basis vectors associated with the HO D
    real(kind=Rkind) :: w      = ZERO                                          ! eigenpulsation associated with the HO D
    real(kind=Rkind) :: m      = ZERO                                          ! mass associated with the HO D
    real(kind=Rkind) :: lambda = -ONE                                          ! strength parameter of the coupling between the mode D and the molecule
    real(kind=Rkind) :: eq_pos = -ONE                                          ! equilibrium position of the HO
  END TYPE


  INTERFACE Read_cavity_mode
    MODULE PROCEDURE MolecCav_Read_cavity_mode
  END INTERFACE
  INTERFACE Write_cavity_mode
    MODULE PROCEDURE MolecCav_Write_cavity_mode
  END INTERFACE
    

  CONTAINS


  SUBROUTINE MolecCav_Read_cavity_mode(Mode, nio)                              ! nio is the label of the file from which the values have to be drawn.
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    IMPLICIT NONE
    
    TYPE(Cavity_mode_t),    intent(inout) :: Mode   
    integer,                intent(in)    :: nio

    integer                               :: D, Nb, err_io                     ! label of the basis/HO/mode/dimension, its number of basis vectors, and an error control variable
    real(kind=Rkind)                      :: w, m, lambda, eq_pos              ! eigenpulsation, mass, molecule-coupling strength, and equilibrium position associated with this HO
    logical, parameter                    :: Debug = .FALSE.

    NAMELIST /HO_1/ D, Nb, w, m, lambda, eq_pos                                ! declare the nml HO_1 and specify the parameter's list to be found within

    !----------------------Initialization to default values--------------------
    D      = 0
    Nb     = 1
    w      = ZERO
    m      = ZERO
    lambda = -ONE
    eq_pos = -ONE
 
    !------------------------------Reading of the nml--------------------------
    WRITE(out_unit,*) 
    WRITE(out_unit,*) '********************************************************************************'
    WRITE(out_unit,*) '**************************** READING BASIS OF THE HO ***************************'
    WRITE(out_unit,*) '********************************************************************************'
    
    READ(nio, nml = HO_1, iostat = err_io)                                     ! assign the values read in the nml to the declared list of parameters

    IF (Debug) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "-----------------------The namelist parameters are read as----------------------"
      WRITE(out_unit, nml = HO_1)
      WRITE(out_unit,*) "-------------------------End of the namelist parameters-------------------------"
    END IF
    
    !------------------------------Check reading error-------------------------
    IF(err_io < 0) THEN
      WRITE(out_unit,*) ''
      WRITE(out_unit,*) '#######################################################'
      WRITE(out_unit,*) '######## Error in Read_cavity_mode (err_io<0) #########'
      WRITE(out_unit,*) '#######################################################'
      WRITE(out_unit,*) '################ err_io = ', err_io, '################'
      STOP '################# Check basis data ################'

    ELSE IF( err_io > 0) THEN
      WRITE(out_unit,*) ''
      WRITE(out_unit,*) '#######################################################'
      WRITE(out_unit,*) '######## Error in Read_cavity_mode (err_io>0) #########'
      WRITE(out_unit,*) '#######################################################'
      WRITE(out_unit,*) '################ err_io = ', err_io, '################'
      STOP '################# Check basis data ################'

    END IF

    IF (Nb == 0) THEN
      WRITE(out_unit,*) "The number of basis vector associated to any HO CANNOT be 0 (what are are you going to study if there is&
                       & no system ???). Please check the data file '.nml'"
      STOP "The number of basis vector associated to any HO CANNOT be 0 (what are are you going to study if there is&
                       & no system ???). Please check the data file '.nml'"
    END IF
    
    !---------------Construction of the Cavity_mode_t type object-----------
    Mode%D      = D
    Mode%Nb     = Nb
    Mode%w      = w
    Mode%m      = m
    Mode%lambda = lambda
    Mode%eq_pos = eq_pos

    WRITE(out_unit,*) 
    WRITE(out_unit,*) '********************************************************************************'
    WRITE(out_unit,*) '************************** BASIS OF THE HO CONSTRUCTED *************************'
    WRITE(out_unit,*) '********************************************************************************'

    IF (Debug) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "--------------Cavity mode constructed by MolecCav_Read_cavity_mode--------------"
      CALL Write_cavity_mode(Mode)
      WRITE(out_unit,*) "------------End Cavity mode constructed by MolecCav_Read_cavity_mode------------"
    END IF

  END SUBROUTINE MolecCav_Read_cavity_mode

  
  SUBROUTINE MolecCav_Write_cavity_mode(Mode)
    TYPE(Cavity_mode_t), intent(in) :: Mode

    WRITE(out_unit,*) "____________________________________The associated HO cavity mode___________________________________"
    WRITE(out_unit,*) "|Index of the cavity mode Mode%D                                             | ", Mode%D
    WRITE(out_unit,*) "|____________________________________________________________________________|______________________"
    WRITE(out_unit,*) "|Basis set size of the cavity mode Mode%Nb                                   | ", Mode%Nb
    WRITE(out_unit,*) "|____________________________________________________________________________|______________________"
    WRITE(out_unit,*) "|Eigenpulsation of the cavity mode Mode%w                                    | ", Mode%w
    WRITE(out_unit,*) "|____________________________________________________________________________|______________________"
    WRITE(out_unit,*) "|Mass of the associated HO Mode%m                                            | ", Mode%m
    WRITE(out_unit,*) "|____________________________________________________________________________|______________________"
    WRITE(out_unit,*) "|Coupling strength between the matter and this cavity mode Mode%lambda       | ", Mode%lambda
    WRITE(out_unit,*) "|____________________________________________________________________________|______________________"
    WRITE(out_unit,*) "|Equilibrium position of this cavity mode's HO Mode%eq_pos                   | ", Mode%eq_pos
    WRITE(out_unit,*) "|____________________________________________________________________________|______________________"
    FLUSH(out_unit)

  END SUBROUTINE MolecCav_Write_cavity_mode


END MODULE